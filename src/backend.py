"""IQM Resonance backend connection and topology-aware qubit selection."""

from __future__ import annotations

import os

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


def filter_qubits_by_metrics(
    qubit_metrics: dict[int, dict],
    cz_fidelities: dict[tuple[int, int], float] | None = None,
    min_readout_fidelity: float = 0.95,
    min_one_q_fidelity: float = 0.998,
    min_t1: float = 30e-6,
    min_t2: float = 20e-6,
    min_cz_fidelity: float = 0.95,
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


def select_best_tree(
    backend,
    target_n: int,
    apply_thresholds: bool = True,
    min_readout_fidelity: float = 0.90,
    min_one_q_fidelity: float = 0.99,
    min_t1: float = 10e-6,
    min_t2: float = 5e-6,
    min_cz_fidelity: float = 0.90,
) -> dict:
    """Find the best target_n-qubit connected subgraph of the device.

    Strategy: filter qubits by per-qubit thresholds, restrict device coupling
    to the survivors with edge weight = 1 − CZ_fidelity (use 1.0 for unknown
    edges so they are unattractive), then for every starting qubit run Prim's
    algorithm to grow a tree of size target_n. Return the lightest such tree.

    Returns a dict with:
      qubits   : list of physical qubit indices (size target_n)
      edges    : list of (a, b) physical edges in the spanning tree (size target_n - 1)
      weight   : sum of edge weights (Σ 1 − CZ_fidelity over tree edges)
      logical_edges : same edges, relabeled to 0..target_n-1 in row-major order
      coloring : 2-coloring of the tree (0/1 per logical qubit)
      n_starts_tried : number of seeds explored
    """
    if is_simulator(backend):
        # Synthetic: line graph 0–1–2–…–(n-1)
        edges = [(i, i + 1) for i in range(target_n - 1)]
        coloring = [i % 2 for i in range(target_n)]
        return {
            "qubits": list(range(target_n)),
            "edges": [(i, i + 1) for i in range(target_n - 1)],
            "weight": 0.0,
            "logical_edges": edges,
            "coloring": coloring,
            "n_starts_tried": 1,
        }

    qm, cz = get_qubit_metrics(backend)

    if apply_thresholds:
        good = filter_qubits_by_metrics(
            qm, cz, min_readout_fidelity, min_one_q_fidelity, min_t1, min_t2,
            min_cz_fidelity)
    else:
        good = set(range(backend.num_qubits))

    if len(good) < target_n:
        raise ValueError(f"Only {len(good)} qubits pass thresholds — need {target_n}.")

    # Restricted weighted graph: keep edges with both endpoints in `good` AND
    # CZ fidelity above min_cz_fidelity.
    nbrs: dict[int, list[tuple[int, float]]] = {q: [] for q in good}
    for (a, b), fid in cz.items():
        if a not in good or b not in good:
            continue
        if fid < min_cz_fidelity:
            continue
        w = 1.0 - fid  # smaller = better
        nbrs[a].append((b, w))
        nbrs[b].append((a, w))

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

    # Build a stable mapping physical→logical (row-major over `qubits` order).
    qubits = best["qubits"]
    phys_to_log = {p: i for i, p in enumerate(qubits)}
    logical_edges = [(phys_to_log[a], phys_to_log[b]) for a, b in best["phys_edges"]]

    from src.circuits.graph_state import two_coloring
    coloring = two_coloring(target_n, logical_edges)
    if coloring is None:
        raise RuntimeError("Spanning tree should always be 2-colorable.")

    return {
        "qubits": qubits,
        "edges": best["phys_edges"],
        "logical_edges": logical_edges,
        "weight": best["weight"],
        "coloring": coloring,
        "n_starts_tried": len(good),
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
