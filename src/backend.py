"""IQM Resonance backend connection and topology-aware qubit selection."""

from __future__ import annotations

import os

import numpy as np
import rustworkx as rx
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

try:
    from iqm.qiskit_iqm import IQMProvider
    IQM_AVAILABLE = True
except ImportError:
    IQM_AVAILABLE = False

RESONANCE_URL = "https://resonance.meetiqm.com"


def get_backend(token: str | None = None, device: str = "emerald"):
    """Connect to IQM Resonance and return a Qiskit backend.

    If `token` is None and IQM_TOKEN env var is unset, falls back to AerSimulator.
    The IQM client errors if both env var and arg are set, so we pass arg only
    when the env var is absent.
    """
    env_token = os.environ.get("IQM_TOKEN")
    effective = token or env_token
    if effective and IQM_AVAILABLE:
        if env_token:
            provider = IQMProvider(RESONANCE_URL, quantum_computer=device)
        else:
            provider = IQMProvider(RESONANCE_URL, quantum_computer=device, token=token)
        return provider.get_backend()
    print("No IQM token — using Aer simulator.")
    return AerSimulator()


def is_simulator(backend) -> bool:
    return isinstance(backend, AerSimulator)


def get_qubit_metrics(backend) -> tuple[dict[int, dict], dict[tuple[int, int], float]]:
    """Pull per-qubit + per-CZ-pair calibration metrics from the IQM device.

    Returns (qubit_metrics, cz_fidelities) where:
      qubit_metrics = {idx: {readout_fidelity, error_0_to_1, error_1_to_0,
                              one_q_fidelity, t1, t2}}
      cz_fidelities = {(idx_a, idx_b): fidelity}  (sorted tuple)
    """
    if is_simulator(backend):
        return {}, {}
    qms = backend.client.get_quality_metric_set()
    name_to_idx = {f"QB{i+1}": i for i in range(backend.num_qubits)}
    qubits: dict[int, dict] = {i: {} for i in range(backend.num_qubits)}
    pairs: dict[tuple[int, int], float] = {}
    for o in qms.observations:
        f = o.dut_field
        v = o.value
        if f.startswith("metrics.ssro.measure.constant."):
            parts = f.split(".")
            idx = name_to_idx.get(parts[4])
            if idx is None: continue
            metric = ".".join(parts[5:])
            if metric == "fidelity": qubits[idx]["readout_fidelity"] = v
            elif metric == "error_0_to_1": qubits[idx]["error_0_to_1"] = v
            elif metric == "error_1_to_0": qubits[idx]["error_1_to_0"] = v
        elif f.startswith("metrics.rb.prx.drag_crf_sx."):
            parts = f.split(".")
            idx = name_to_idx.get(parts[4])
            if idx is None: continue
            if parts[5].startswith("fidelity"):
                qubits[idx]["one_q_fidelity"] = v
        elif f.startswith("metrics.rb.clifford.uz_cz."):
            parts = f.split(".")
            pair = parts[4]
            if "__" not in pair: continue
            qa, qb = pair.split("__")
            ia, ib = name_to_idx.get(qa), name_to_idx.get(qb)
            if ia is None or ib is None: continue
            if parts[5].startswith("fidelity"):
                pairs[(min(ia, ib), max(ia, ib))] = v
        elif f.startswith("characterization.model."):
            parts = f.split(".")
            idx = name_to_idx.get(parts[2])
            if idx is None: continue
            if parts[3] == "t1_time": qubits[idx]["t1"] = v
            elif parts[3] == "t2_echo_time": qubits[idx]["t2"] = v
    return qubits, pairs


DEFAULT_THRESHOLDS = dict(
    min_readout_fidelity=0.90,
    min_one_q_fidelity=0.99,
    min_t1=10e-6,
    min_t2=5e-6,
    min_cz_fidelity=0.90,
)


def filter_qubits_by_metrics(
    qubit_metrics: dict[int, dict],
    cz_fidelities: dict[tuple[int, int], float] | None = None,
    min_readout_fidelity: float = 0.90,
    min_one_q_fidelity: float = 0.99,
    min_t1: float = 10e-6,
    min_t2: float = 5e-6,
    min_cz_fidelity: float = 0.90,
) -> set[int]:
    """Return qubit indices passing per-qubit thresholds and (optionally) having
    at least one CZ partner above min_cz_fidelity.
    """
    good: set[int] = set()
    for idx, m in qubit_metrics.items():
        if not (m.get("readout_fidelity", 0) >= min_readout_fidelity
                and m.get("one_q_fidelity", 0) >= min_one_q_fidelity
                and m.get("t1", 0) >= min_t1
                and m.get("t2", 0) >= min_t2):
            continue
        good.add(idx)
    # CZ-pair filter: a qubit is excluded if ALL its CZ pairs are below threshold
    if cz_fidelities:
        bad_cz: set[int] = set()
        per_qubit_pairs: dict[int, list[float]] = {}
        for (a, b), f in cz_fidelities.items():
            per_qubit_pairs.setdefault(a, []).append(f)
            per_qubit_pairs.setdefault(b, []).append(f)
        for idx in list(good):
            fids = per_qubit_pairs.get(idx, [])
            if not fids or max(fids) < min_cz_fidelity:
                bad_cz.add(idx)
        good -= bad_cz
    return good


def _enumerate_subgrids(backend, rows: int, cols: int,
                        max_grids: int = 500,
                        excluded_qubits: set[int] | None = None) -> list[list[int]]:
    """Find all contiguous rows×cols rectangular sub-grids on the device topology.

    Uses VF2 subgraph isomorphism. Returns layouts in row-major order matching
    the cluster_2d circuit's logical qubit numbering. Up to `max_grids` results.
    """
    target = rx.PyGraph()
    target.add_nodes_from(range(rows * cols))
    for r in range(rows):
        for c in range(cols - 1):
            target.add_edge(r * cols + c, r * cols + c + 1, None)
    for r in range(rows - 1):
        for c in range(cols):
            target.add_edge(r * cols + c, (r + 1) * cols + c, None)

    excluded = excluded_qubits or set()
    device = rx.PyGraph()
    device.add_nodes_from(range(backend.num_qubits))
    seen_edges: set[tuple[int, int]] = set()
    for a, b in backend.coupling_map:
        if a in excluded or b in excluded:
            continue
        e = (min(a, b), max(a, b))
        if e in seen_edges:
            continue
        seen_edges.add(e)
        device.add_edge(a, b, None)

    grids: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for m in rx.vf2_mapping(device, target, subgraph=True, induced=False):
        inv = {tgt: dev for dev, tgt in m.items()}
        layout = [inv[i] for i in range(rows * cols)]
        key = tuple(layout)
        if key in seen:
            continue
        seen.add(key)
        grids.append(layout)
        if len(grids) >= max_grids:
            break
    return grids


CZ_GATE_TIME_S    = 100e-9   # IQM Emerald CZ gate duration  ~100 ns
MEASUREMENT_TIME_S = 1.5e-6  # IQM measurement duration       ~1.5 µs (typical)


def _compute_edge_weights(
    qubit_metrics: dict[int, dict],
    cz_fidelities: dict[tuple[int, int], float],
    readout_weight: float = 0.0,
    t1_weight: float = 1.0,
    t2_weight: float = 1.0,
    one_q_weight: float = 1.0,
    measurement_decoherence: bool = False,
    gate_time_s: float = CZ_GATE_TIME_S,
    measurement_time_s: float = MEASUREMENT_TIME_S,
) -> dict[tuple[int, int], float]:
    """Comprehensive edge cost capturing all dominant noise channels.

      w(u,v) = (1 - F_CZ(u,v))                                 CZ gate error
             + readout_weight · ((1-F_RO(u)) + (1-F_RO(v)))    (default 0; QREM)
             + t1_weight · t_eff · (1/T1(u) + 1/T1(v))         amplitude damping
             + t2_weight · t_eff · (1/T2(u) + 1/T2(v))         dephasing
             + one_q_weight · ((1-F_1Q(u)) + (1-F_1Q(v)))      1Q gate error

    where t_eff = gate_time + (measurement_time / 2 if measurement_decoherence
    else 0). Including measurement decoherence accounts for the fact that
    qubits decay during the ~1.5 µs readout, comparable to several CZ gates.

    First-order infidelity for amplitude damping over time t is t/T1, and
    for dephasing it is t/T2 — coefficients match these textbook formulae
    (no spurious 1/2 factor).
    """
    t_eff = gate_time_s
    if measurement_decoherence:
        t_eff += 0.5 * measurement_time_s   # /2: measurement contributes half
                                            # because not every CZ is followed
                                            # immediately by readout
    out: dict[tuple[int, int], float] = {}
    for (a, b), fid in cz_fidelities.items():
        m_a = qubit_metrics.get(a, {})
        m_b = qubit_metrics.get(b, {})
        w = 1.0 - fid

        if readout_weight:
            ro_a = m_a.get("readout_fidelity", 1.0)
            ro_b = m_b.get("readout_fidelity", 1.0)
            w += readout_weight * ((1.0 - ro_a) + (1.0 - ro_b))

        if t1_weight:
            t1_a = m_a.get("t1") or 1e9
            t1_b = m_b.get("t1") or 1e9
            w += t1_weight * t_eff * (1.0 / t1_a + 1.0 / t1_b)

        if t2_weight:
            t2_a = m_a.get("t2") or 1e9
            t2_b = m_b.get("t2") or 1e9
            w += t2_weight * t_eff * (1.0 / t2_a + 1.0 / t2_b)

        if one_q_weight:
            f1_a = m_a.get("one_q_fidelity", 1.0)
            f1_b = m_b.get("one_q_fidelity", 1.0)
            w += one_q_weight * ((1.0 - f1_a) + (1.0 - f1_b))

        out[(min(a, b), max(a, b))] = w
    return out


def predict_stabilizer_fidelities(
    qubits: list[int],
    logical_edges: list[tuple[int, int]],
    qubit_metrics: dict[int, dict],
    cz_fidelities: dict[tuple[int, int], float],
    gate_time_s: float = CZ_GATE_TIME_S,
    measurement_time_s: float = MEASUREMENT_TIME_S,
    measurement_decoherence: bool = False,
    include_readout: bool = False,
) -> dict[int, float]:
    """Predict ⟨g_i⟩ for each logical qubit i.

    Default: ignores readout (assume parity-QREM is applied post-hoc).
    Set ``include_readout=True`` for the raw (un-mitigated) prediction.

    Heuristic model: each stabilizer g_i = X_i ⊗ ⊗_(j∈N(i)) Z_j is degraded by:
      - CZ errors on edges incident to {i} ∪ N(i) (depolarizing channel
        approximation: each error reduces ⟨g_i⟩ by factor F_CZ)
      - Per-qubit 1Q gate error (1 H per qubit prep; +1 H if X-measured)
      - Per-qubit decoherence over total time t_eff = gate + 0.5·readout
      - (Optional) Per-qubit readout error

    First-order coefficients: amplitude-damping infidelity = t/T1, dephasing = t/T2.
    """
    from src.circuits.graph_state import neighbours_from_edges

    n = len(qubits)
    nbrs = neighbours_from_edges(n, logical_edges)
    t_eff = gate_time_s + (0.5 * measurement_time_s if measurement_decoherence else 0.0)
    out: dict[int, float] = {}
    for i in range(n):
        involved_logical = [i] + nbrs[i]
        involved_phys = [qubits[lq] for lq in involved_logical]
        f = 1.0

        # CZ fidelities on edges incident to any qubit in the i-neighborhood
        for (a, b) in logical_edges:
            if a in involved_logical or b in involved_logical:
                pa, pb = qubits[a], qubits[b]
                cz = cz_fidelities.get((min(pa, pb), max(pa, pb)), 0.99)
                f *= cz

        for pq in involved_phys:
            mp = qubit_metrics.get(pq, {})
            # 1Q error: 1 H during prep + 1 more H if measured in X basis.
            # We don't know the coloring here so use 1 average; ~0.05% effect.
            f *= mp.get("one_q_fidelity", 1.0)
            # First-order decoherence: t/T1 + t/T2 (no spurious 1/2)
            t1 = mp.get("t1") or 1e9
            t2 = mp.get("t2") or 1e9
            f *= max(0.0, 1.0 - t_eff / t1 - t_eff / t2)
            if include_readout:
                f *= mp.get("readout_fidelity", 1.0)
        out[i] = f
    return out


def _local_search_swap(
    qubits_in_tree: set[int],
    tree_edges: list[tuple[int, int]],
    edge_weight: dict[tuple[int, int], float],
    nbrs: dict[int, set[int]],
    max_iterations: int = 50,
) -> tuple[set[int], list[tuple[int, int]], float]:
    """Local search: try replacing each tree edge with a non-tree edge that keeps
    the tree connected and reduces total weight. Repeat until no improvement.
    """
    in_tree = set(qubits_in_tree)
    edges = [tuple(sorted(e)) for e in tree_edges]
    total_w = sum(edge_weight.get(e, 1.0) for e in edges)

    for _ in range(max_iterations):
        improved = False
        for i, e_remove in enumerate(list(edges)):
            # Removing e_remove splits the tree into two components.
            # Find them via BFS in tree_edges \ {e_remove}.
            adj_t: dict[int, list[int]] = {q: [] for q in in_tree}
            for ee in edges:
                if ee == e_remove:
                    continue
                adj_t[ee[0]].append(ee[1])
                adj_t[ee[1]].append(ee[0])
            # BFS from one endpoint of removed edge
            visited = {e_remove[0]}
            stack = [e_remove[0]]
            while stack:
                u = stack.pop()
                for v in adj_t[u]:
                    if v not in visited:
                        visited.add(v)
                        stack.append(v)
            comp_a, comp_b = visited, in_tree - visited

            # Look for a replacement edge crossing the cut, with smaller weight
            best_replacement = None
            best_w = edge_weight.get(e_remove, 1.0)
            for u in comp_a:
                for v in nbrs.get(u, set()):
                    if v not in comp_b:
                        continue
                    e_new = (min(u, v), max(u, v))
                    if e_new == e_remove:
                        continue
                    w_new = edge_weight.get(e_new, 1.0)
                    if w_new < best_w:
                        best_w = w_new
                        best_replacement = e_new

            if best_replacement is not None:
                w_old = edge_weight.get(e_remove, 1.0)
                edges[i] = best_replacement
                total_w += (best_w - w_old)
                improved = True
                break  # restart loop after any improvement

        if not improved:
            break
    return in_tree, edges, total_w


def select_best_tree(
    backend,
    target_n: int,
    apply_thresholds: bool = True,
    min_readout_fidelity: float = 0.90,
    min_one_q_fidelity: float = 0.99,
    min_t1: float = 10e-6,
    min_t2: float = 5e-6,
    min_cz_fidelity: float = 0.90,
    readout_weight: float = 0.0,
    t1_weight: float = 1.0,
    t2_weight: float = 1.0,
    one_q_weight: float = 1.0,
    measurement_decoherence: bool = False,
    use_local_search: bool = True,
) -> dict:
    """Find the best target_n-qubit connected subgraph of the device.

    Cost minimised per tree edge:
      (1 - F_CZ)                                              CZ gate infidelity
      + t1_weight * t_eff * (1/T1_u + 1/T1_v)                amplitude damping
      + t2_weight * t_eff * (1/T2_u + 1/T2_v)                dephasing
      + one_q_weight * ((1-F_1Q_u) + (1-F_1Q_v))             1Q gate error
      (readout_weight * RO penalty — off by default since QREM corrects RO)

    where t_eff = gate_time + 0.5 * measurement_time when
    `measurement_decoherence=True`. Default OFF: measurement-time decoherence
    is partly captured by readout fidelity (which QREM mitigates), so
    including it again here over-counts. Set to True for a conservative
    lower bound on prep fidelity.

    Quality threshold filter (default ON) excludes qubits with marginal
    calibration that empirically drift more than reported metrics suggest.

    Returns dict with: qubits, edges, logical_edges, weight, coloring,
    predicted_W, predicted_stab_fidelities, n_starts_tried.
    """
    if is_simulator(backend):
        # Synthetic line graph; predict ideal W since no noise on simulator
        edges = [(i, i + 1) for i in range(target_n - 1)]
        coloring = [i % 2 for i in range(target_n)]
        return {
            "qubits": list(range(target_n)),
            "edges": [(i, i + 1) for i in range(target_n - 1)],
            "weight": 0.0,
            "logical_edges": edges,
            "coloring": coloring,
            "n_starts_tried": 1,
            "predicted_W": float(target_n),
            "predicted_stab_fidelities": {i: 1.0 for i in range(target_n)},
        }

    qm, cz = get_qubit_metrics(backend)

    if apply_thresholds:
        good = filter_qubits_by_metrics(
            qm, cz, min_readout_fidelity, min_one_q_fidelity, min_t1, min_t2,
            min_cz_fidelity)
    else:
        # Only exclude qubits with completely zero CZ fidelity (truly broken).
        good = set(range(backend.num_qubits))
        bad = {a for (a, b), f in cz.items() if f <= 0.0} | \
              {b for (a, b), f in cz.items() if f <= 0.0}
        good = good - bad

    if len(good) < target_n:
        raise ValueError(f"Only {len(good)} qubits available — need {target_n}.")

    # Comprehensive edge cost: CZ + decoherence (T1, T2) + 1Q gate error.
    # Readout deliberately excluded (QREM mitigates it post-hoc).
    edge_weight = _compute_edge_weights(
        qm, cz,
        readout_weight=readout_weight,
        t1_weight=t1_weight,
        t2_weight=t2_weight,
        one_q_weight=one_q_weight,
        measurement_decoherence=measurement_decoherence,
    )

    # Restricted weighted graph
    nbrs_unweighted: dict[int, set[int]] = {q: set() for q in good}
    nbrs: dict[int, list[tuple[int, float]]] = {q: [] for q in good}
    for (a, b), fid in cz.items():
        if a not in good or b not in good:
            continue
        if fid < min_cz_fidelity:
            continue
        w = edge_weight.get((min(a, b), max(a, b)), 1.0 - fid)
        nbrs[a].append((b, w))
        nbrs[b].append((a, w))
        nbrs_unweighted[a].add(b)
        nbrs_unweighted[b].add(a)

    # Prim's algorithm from each seed; track lightest tree of exactly target_n nodes.
    import heapq
    best = None
    for seed in good:
        if len(nbrs[seed]) == 0:
            continue
        in_tree = {seed}
        tree_edges: list[tuple[int, int, float]] = []
        heap: list[tuple[float, int, int]] = []
        for nb, w in nbrs[seed]:
            heapq.heappush(heap, (w, seed, nb))
        total_w = 0.0
        while heap and len(in_tree) < target_n:
            w, u, v = heapq.heappop(heap)
            if v in in_tree:
                continue
            in_tree.add(v)
            tree_edges.append((u, v, w))
            total_w += w
            for nb, w2 in nbrs[v]:
                if nb not in in_tree:
                    heapq.heappush(heap, (w2, v, nb))
        if len(in_tree) != target_n:
            continue
        if best is None or total_w < best["weight"]:
            best = {
                "qubits": sorted(in_tree),  # physical indices
                "phys_edges": [(a, b) for a, b, _ in tree_edges],
                "weight": total_w,
            }

    if best is None:
        raise ValueError(f"No connected {target_n}-qubit subtree on filtered graph.")

    # Local search refinement: try edge swaps that reduce total weight while
    # keeping the tree connected and the same set of nodes.
    if use_local_search:
        in_tree_refined, edges_refined, weight_refined = _local_search_swap(
            set(best["qubits"]),
            [(a, b) for a, b in best["phys_edges"]],
            edge_weight,
            nbrs_unweighted,
        )
        if weight_refined < best["weight"]:
            best = {
                "qubits": sorted(in_tree_refined),
                "phys_edges": edges_refined,
                "weight": weight_refined,
            }

    # Build a stable mapping physical→logical (row-major over `qubits` order).
    qubits = best["qubits"]
    phys_to_log = {p: i for i, p in enumerate(qubits)}
    logical_edges = [(phys_to_log[a], phys_to_log[b]) for a, b in best["phys_edges"]]

    from src.circuits.graph_state import two_coloring
    coloring = two_coloring(target_n, logical_edges)
    if coloring is None:
        raise RuntimeError("Spanning tree should always be 2-colorable.")

    # Predict per-stabilizer fidelity and W for the chosen tree
    pred = predict_stabilizer_fidelities(qubits, logical_edges, qm, cz)
    predicted_W = sum(pred.values())

    return {
        "qubits": qubits,
        "edges": best["phys_edges"],
        "logical_edges": logical_edges,
        "weight": best["weight"],
        "coloring": coloring,
        "n_starts_tried": len(good),
        "predicted_stab_fidelities": pred,
        "predicted_W": predicted_W,
    }


def select_best_tree_empirical(
    backend,
    target_n: int,
    edge_fidelities: dict[tuple[int, int], float],
    *,
    min_edge_F: float = 0.5,
    use_local_search: bool = True,
) -> dict:
    """Same multi-start Prim's MST search as `select_best_tree`, but using a
    *measured* per-edge graph-state fidelity F_ij as the cost.

    Edge weight  w(u,v) = 1 - F_ij  for any (u,v) in `edge_fidelities` with
    F_ij >= `min_edge_F`. Edges below the threshold (and not natively present
    in `edge_fidelities`) are excluded from the search.

    `edge_fidelities` is the output of
    `src.diagnostics.run_edge_map` post-processed via `fidelity_map`.

    Returns the same dict shape as `select_best_tree` (qubits, edges,
    logical_edges, weight, coloring, n_starts_tried). `predicted_W` here is
    the sum of measured F_ij over tree edges — a direct empirical proxy for
    GME witness quality (no calibration model involved).
    """
    if is_simulator(backend):
        edges = [(i, i + 1) for i in range(target_n - 1)]
        coloring = [i % 2 for i in range(target_n)]
        return {
            "qubits": list(range(target_n)),
            "edges": edges,
            "weight": 0.0,
            "logical_edges": edges,
            "coloring": coloring,
            "n_starts_tried": 1,
            "predicted_W": float(target_n),
            "predicted_stab_fidelities": {i: 1.0 for i in range(target_n)},
            "measured_edge_F": {e: 1.0 for e in edges},
            "min_edge_F_used": 1.0,
            "mean_edge_F_used": 1.0,
        }

    norm = lambda a, b: (min(a, b), max(a, b))
    fids = {norm(a, b): F for (a, b), F in edge_fidelities.items()
            if F >= min_edge_F}
    if not fids:
        raise ValueError("No edges above min_edge_F threshold.")

    edge_weight = {e: 1.0 - F for e, F in fids.items()}
    nbrs_unweighted: dict[int, set[int]] = {}
    nbrs: dict[int, list[tuple[int, float]]] = {}
    qubits_present: set[int] = set()
    for (a, b), w in edge_weight.items():
        qubits_present.update((a, b))
        nbrs_unweighted.setdefault(a, set()).add(b)
        nbrs_unweighted.setdefault(b, set()).add(a)
        nbrs.setdefault(a, []).append((b, w))
        nbrs.setdefault(b, []).append((a, w))

    if len(qubits_present) < target_n:
        raise ValueError(f"Only {len(qubits_present)} qubits in measured map "
                          f"above F={min_edge_F}; need {target_n}.")

    import heapq
    best = None
    for seed in qubits_present:
        if not nbrs.get(seed):
            continue
        in_tree = {seed}
        tree_edges: list[tuple[int, int, float]] = []
        heap: list[tuple[float, int, int]] = []
        for nb, w in nbrs[seed]:
            heapq.heappush(heap, (w, seed, nb))
        total_w = 0.0
        while heap and len(in_tree) < target_n:
            w, u, v = heapq.heappop(heap)
            if v in in_tree:
                continue
            in_tree.add(v)
            tree_edges.append((u, v, w))
            total_w += w
            for nb, w2 in nbrs[v]:
                if nb not in in_tree:
                    heapq.heappush(heap, (w2, v, nb))
        if len(in_tree) != target_n:
            continue
        if best is None or total_w < best["weight"]:
            best = {
                "qubits": sorted(in_tree),
                "phys_edges": [(a, b) for a, b, _ in tree_edges],
                "weight": total_w,
            }

    if best is None:
        raise ValueError(f"No connected {target_n}-qubit subtree above "
                          f"F={min_edge_F}.")

    if use_local_search:
        in_tree_r, edges_r, w_r = _local_search_swap(
            set(best["qubits"]),
            [(a, b) for a, b in best["phys_edges"]],
            edge_weight,
            nbrs_unweighted,
        )
        if w_r < best["weight"]:
            best = {"qubits": sorted(in_tree_r), "phys_edges": edges_r,
                    "weight": w_r}

    qubits = best["qubits"]
    phys_to_log = {p: i for i, p in enumerate(qubits)}
    logical_edges = [(phys_to_log[a], phys_to_log[b]) for a, b in best["phys_edges"]]

    from src.circuits.graph_state import two_coloring
    coloring = two_coloring(target_n, logical_edges)
    if coloring is None:
        raise RuntimeError("Spanning tree should always be 2-colorable.")

    measured_F = {(min(a, b), max(a, b)): fids[(min(a, b), max(a, b))]
                  for (a, b) in best["phys_edges"]}
    # Empirical proxy for stabilizer fidelity: each qubit's local F is the
    # geometric mean of incident measured edge fidelities. Sum -> predicted_W.
    pred: dict[int, float] = {}
    for i, p in enumerate(qubits):
        incident = [F for e, F in measured_F.items() if p in e]
        pred[i] = float(np.prod(incident) ** (1.0 / max(1, len(incident)))) if incident else 1.0
    predicted_W = sum(pred.values())

    return {
        "qubits": qubits,
        "edges": best["phys_edges"],
        "logical_edges": logical_edges,
        "weight": best["weight"],
        "coloring": coloring,
        "n_starts_tried": len(qubits_present),
        "predicted_stab_fidelities": pred,
        "predicted_W": predicted_W,
        "measured_edge_F": measured_F,
        "min_edge_F_used": min(measured_F.values()),
        "mean_edge_F_used": float(np.mean(list(measured_F.values()))),
    }


def select_best_subgrid(backend, rows: int, cols: int,
                         readout_mode: str = "fidelity",
                         excluded_qubits: set[int] | None = None,
                         apply_thresholds: bool = True,
                         min_readout_fidelity: float = 0.95,
                         min_one_q_fidelity: float = 0.998,
                         min_t1: float = 30e-6,
                         min_t2: float = 20e-6,
                         min_cz_fidelity: float = 0.95,
                         ) -> tuple[list[int], float, int]:
    """Pick the highest-fidelity contiguous rows×cols sub-grid.

    1. Optionally filter out bad qubits via thresholds on calibration metrics.
    2. Enumerate every contiguous rectangular sub-grid over the surviving qubits.
    3. Score each via IQM's CostEvaluator (live CZ + readout calibration).
    4. Return (layout, cost, n_candidates).

    Layout matches the cluster_2d circuit's row-major qubit ordering, so passing
    `initial_layout=layout` to transpile gives zero SWAP overhead.
    """
    if is_simulator(backend):
        return list(range(rows * cols)), 0.0, 1

    excluded = set(excluded_qubits or set())
    if apply_thresholds:
        qm, cz = get_qubit_metrics(backend)
        good = filter_qubits_by_metrics(
            qm, cz, min_readout_fidelity, min_one_q_fidelity, min_t1, min_t2,
            min_cz_fidelity)
        all_qubits = set(range(backend.num_qubits))
        excluded = excluded | (all_qubits - good)

    candidates = _enumerate_subgrids(backend, rows, cols, excluded_qubits=excluded)
    if not candidates:
        raise ValueError(
            f"No {rows}×{cols} contiguous sub-grid on the {backend.num_qubits - len(excluded)} "
            f"qubits passing thresholds. Lower thresholds or pick a smaller grid.")

    try:
        from iqm.qubit_selector.qubit_selector import (
            CostEvaluator, CostFunction, ReadoutMode,
        )
    except ImportError:
        return candidates[0], 0.0, len(candidates)

    rm = {"none": ReadoutMode.NONE, "fidelity": ReadoutMode.FIDELITY,
          "qndness": ReadoutMode.QNDNESS}.get(readout_mode, ReadoutMode.FIDELITY)

    n = rows * cols
    canonical = QuantumCircuit(n)
    canonical.h(range(n))
    for r in range(rows):
        for c in range(cols - 1):
            canonical.cz(r * cols + c, r * cols + c + 1)
    for r in range(rows - 1):
        for c in range(cols):
            canonical.cz(r * cols + c, (r + 1) * cols + c)

    evaluator = CostEvaluator(
        backend=backend, quantum_circuit=canonical,
        cost_function=CostFunction.GATE_COST_CZ,
        readoutmode=rm, layouts=candidates,
    )
    top_layouts, scores = evaluator.get_top_layouts(num_layouts=1)
    return top_layouts[0], scores[0], len(candidates)
