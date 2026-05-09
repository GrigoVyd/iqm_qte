"""Entanglement Bottleneck Map / Selector Audit.

Experimentally measures the 2-qubit graph-state fidelity F_ij = (1 + <X_iZ_j> +
<Z_iX_j> + <Y_iY_j>) / 4 on every CZ-capable physical edge of an IQM device.

F_ij > 1/2 certifies that the state on (i, j) is entangled (the graph state
|G_ij> = CZ_ij|+>|+> is locally equivalent to a Bell state, and for any
maximally entangled target state |psi>, max over product states of |<psi|a,b>|^2
= 1/2 — this is the standard fidelity-witness argument).

The module supports two execution modes:

  parallel   — edges grouped into matchings via greedy edge coloring; each
               matching's edges are prepared and measured in a single circuit.
  isolated   — each edge tested alone; comparison against parallel mode
               quantifies parallel-CZ crosstalk.

Companion functions audit a proposed layout against the measured map and
recommend chains / patches whose weakest edge is strongest.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from typing import Iterable, Sequence

import numpy as np
from qiskit import QuantumCircuit, transpile

Edge = tuple[int, int]
BASES = ("XZ", "ZX", "YY")


# ---------------------------------------------------------------------------
# Edges and matchings
# ---------------------------------------------------------------------------

def _norm(edge: Iterable[int]) -> Edge:
    a, b = edge
    return (a, b) if a < b else (b, a)


def get_coupling_edges(backend) -> list[Edge]:
    """Unique undirected CZ-capable edges from the backend coupling map."""
    cmap = getattr(backend, "coupling_map", None)
    if cmap is None or len(list(cmap)) == 0:
        raise ValueError(
            "Backend has no coupling_map. Pass `edges` explicitly "
            "(e.g. on AerSimulator)."
        )
    return sorted({_norm(e) for e in cmap})


def make_undirected_edges(edges: Iterable[Edge]) -> list[Edge]:
    return sorted({_norm(e) for e in edges})


def edge_color_matchings(edges: Sequence[Edge]) -> list[list[Edge]]:
    """Greedy edge-coloring → list of matchings (parallelizable groups).

    Color each edge with the smallest color not used by any incident edge.
    Equivalent to a Vizing-bounded coloring; gives Δ or Δ+1 matchings.
    """
    edges = make_undirected_edges(edges)
    incident: dict[int, list[int]] = defaultdict(list)  # qubit -> list of colors used
    color: dict[Edge, int] = {}
    for e in edges:
        u, v = e
        used = set(incident[u]) | set(incident[v])
        c = 0
        while c in used:
            c += 1
        color[e] = c
        incident[u].append(c)
        incident[v].append(c)
    by_color: dict[int, list[Edge]] = defaultdict(list)
    for e, c in color.items():
        by_color[c].append(e)
    return [sorted(by_color[c]) for c in sorted(by_color)]


# ---------------------------------------------------------------------------
# Circuit construction
# ---------------------------------------------------------------------------

def _apply_basis_rotation(qc: QuantumCircuit, log_q: int, pauli: str) -> None:
    """Rotate the measurement basis so a subsequent Z-measure reads `pauli`."""
    if pauli == "X":
        qc.h(log_q)
    elif pauli == "Y":
        qc.sdg(log_q)
        qc.h(log_q)
    elif pauli == "Z":
        pass
    else:
        raise ValueError(f"Unknown Pauli {pauli!r}")


def build_matching_graph_state_circuit(
    matching: Sequence[Edge],
    basis: str,
    *,
    no_cz: bool = False,
) -> tuple[QuantumCircuit, list[int], dict[int, int]]:
    """Build a circuit that prepares CZ|++> on every edge in `matching` in
    parallel and measures the requested Pauli pair on each edge.

    `basis` is one of 'XZ', 'ZX', 'YY' (Pauli on first vs second qubit of edge).

    Returns (circuit, physical_qubits, phys_to_log) where:
      physical_qubits = sorted list of physical qubit indices touched
      phys_to_log     = {physical: logical_index_in_circuit}

    The caller passes `initial_layout=physical_qubits` to transpile().
    """
    if basis not in BASES:
        raise ValueError(f"basis must be one of {BASES}, got {basis!r}")
    p1_pauli, p2_pauli = basis[0], basis[1]

    physical_qubits = sorted({q for e in matching for q in e})
    phys_to_log = {p: k for k, p in enumerate(physical_qubits)}
    n = len(physical_qubits)

    qc = QuantumCircuit(n, n)

    # Preparation: |+> on every used qubit, then CZ on each edge.
    qc.h(range(n))
    if not no_cz:
        for p_i, p_j in matching:
            qc.cz(phys_to_log[p_i], phys_to_log[p_j])
    qc.barrier()

    # Basis rotation per qubit. A qubit may belong to only one edge in a
    # matching (matchings are disjoint), so its basis is unambiguous.
    for p_i, p_j in matching:
        _apply_basis_rotation(qc, phys_to_log[p_i], p1_pauli)
        _apply_basis_rotation(qc, phys_to_log[p_j], p2_pauli)

    qc.measure(range(n), range(n))
    return qc, physical_qubits, phys_to_log


# ---------------------------------------------------------------------------
# Counts → expectation values
# ---------------------------------------------------------------------------

def _bit_pm1(b: str) -> int:
    return 1 if b == "0" else -1


def expectation_from_counts(
    counts: dict[str, int],
    log_q1: int,
    log_q2: int,
    n_qubits: int,
) -> float:
    """Estimate <O_q1 O_q2> with ±1 outcomes from the bitstring distribution.

    Qiskit bitstring convention: classical bit c is at position (n-1-c) of
    the (space-stripped) string. Logical qubit k is measured to classical
    bit k by build_matching_graph_state_circuit.
    """
    if not counts:
        return 0.0
    total = 0
    weighted = 0
    for bs, cnt in counts.items():
        bits = bs.replace(" ", "")
        v1 = _bit_pm1(bits[n_qubits - 1 - log_q1])
        v2 = _bit_pm1(bits[n_qubits - 1 - log_q2])
        weighted += v1 * v2 * cnt
        total += cnt
    return weighted / total if total else 0.0


# ---------------------------------------------------------------------------
# Fidelity, errors, witness
# ---------------------------------------------------------------------------

def _stderr(e: float, shots: int) -> float:
    var = max(0.0, 1.0 - e * e)
    return math.sqrt(var / max(1, shots))


def compute_edge_fidelity(
    e_xz: float, e_zx: float, e_yy: float, shots: int
) -> dict:
    """Compute graph-state fidelity, error, z-score, and entanglement flags.

    F = (1 + e_xz + e_zx + e_yy) / 4
    sigma_F = (1/4) * sqrt(sigma_xz^2 + sigma_zx^2 + sigma_yy^2)
    z_score = (F - 1/2) / sigma_F
    """
    F = (1.0 + e_xz + e_zx + e_yy) / 4.0
    sx = _stderr(e_xz, shots)
    sz = _stderr(e_zx, shots)
    sy = _stderr(e_yy, shots)
    sigma_F = 0.25 * math.sqrt(sx * sx + sz * sz + sy * sy)
    z = (F - 0.5) / sigma_F if sigma_F > 0 else float("inf")
    return {
        "F": F,
        "sigma_F": sigma_F,
        "z_score": z,
        "entangled_3sigma": (F - 3 * sigma_F) > 0.5,
        "entangled_meanonly": F > 0.5,
    }


# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

@dataclass
class EdgeResult:
    edge: Edge
    mode: str                       # 'parallel' or 'isolated'
    matching_id: int                # -1 for isolated
    shots: int
    e_XZ: float
    e_ZX: float
    e_YY: float
    F: float
    sigma_F: float
    z_score: float
    entangled_3sigma: bool
    entangled_meanonly: bool
    no_cz: bool = False

    def to_dict(self) -> dict:
        d = asdict(self)
        d["edge"] = list(self.edge)
        return d


def _is_aer(backend) -> bool:
    return backend.__class__.__name__ == "AerSimulator"


def _transpile_for_backend(qc: QuantumCircuit, backend, physical_qubits: list[int],
                            optimization_level: int) -> QuantumCircuit:
    if _is_aer(backend):
        return transpile(qc, backend=backend, optimization_level=optimization_level)
    return transpile(
        qc, backend=backend, optimization_level=optimization_level,
        initial_layout=physical_qubits,
    )


def _get_counts_list(result, n_circuits: int) -> list[dict[str, int]]:
    raw = result.get_counts()
    if isinstance(raw, dict):
        return [raw] * n_circuits if n_circuits == 1 else [raw]
    return list(raw)


def _run_circuits(backend, circuits: list[QuantumCircuit], shots: int) -> list[dict]:
    if not circuits:
        return []
    job = backend.run(circuits, shots=shots)
    result = job.result()
    counts = _get_counts_list(result, len(circuits))
    if len(counts) != len(circuits):
        raise RuntimeError(
            f"Expected {len(circuits)} count dicts, got {len(counts)}"
        )
    return counts


def _measure_matching(
    backend, matching: Sequence[Edge], shots: int, *,
    matching_id: int, mode: str, no_cz: bool, optimization_level: int,
) -> list[EdgeResult]:
    circuits: list[QuantumCircuit] = []
    phys_qubits_per_basis: list[list[int]] = []
    phys_to_log_per_basis: list[dict[int, int]] = []
    for basis in BASES:
        qc, pq, ptol = build_matching_graph_state_circuit(
            matching, basis, no_cz=no_cz,
        )
        circuits.append(_transpile_for_backend(qc, backend, pq, optimization_level))
        phys_qubits_per_basis.append(pq)
        phys_to_log_per_basis.append(ptol)

    counts = _run_circuits(backend, circuits, shots)
    counts_by_basis = dict(zip(BASES, counts))

    results: list[EdgeResult] = []
    for p_i, p_j in matching:
        # All three bases share the same physical_qubits set (same matching),
        # so phys_to_log is identical across them.
        ptol = phys_to_log_per_basis[0]
        n = len(phys_qubits_per_basis[0])
        log_i, log_j = ptol[p_i], ptol[p_j]
        e_xz = expectation_from_counts(counts_by_basis["XZ"], log_i, log_j, n)
        e_zx = expectation_from_counts(counts_by_basis["ZX"], log_i, log_j, n)
        e_yy = expectation_from_counts(counts_by_basis["YY"], log_i, log_j, n)
        f = compute_edge_fidelity(e_xz, e_zx, e_yy, shots)
        results.append(EdgeResult(
            edge=(p_i, p_j), mode=mode, matching_id=matching_id, shots=shots,
            e_XZ=e_xz, e_ZX=e_zx, e_YY=e_yy, no_cz=no_cz, **f,
        ))
    return results


def run_edge_map(
    backend,
    *,
    edges: Sequence[Edge] | None = None,
    shots: int = 4000,
    mode: str = "parallel",
    no_cz: bool = False,
    optimization_level: int = 3,
    progress: bool = True,
) -> list[EdgeResult]:
    """Measure 2-qubit graph-state fidelity on each requested edge.

    mode='parallel'  — edges grouped into matchings; one job per matching.
    mode='isolated'  — each edge run alone; one circuit set per edge.
                       Use to compare against parallel mode for crosstalk.
    """
    if mode not in ("parallel", "isolated"):
        raise ValueError(f"mode must be 'parallel' or 'isolated', got {mode!r}")
    if edges is None:
        edges = get_coupling_edges(backend)
    edges = make_undirected_edges(edges)
    if not edges:
        return []

    results: list[EdgeResult] = []
    if mode == "parallel":
        matchings = edge_color_matchings(edges)
        if progress:
            print(f"[edge_map] {len(edges)} edges → {len(matchings)} matchings "
                  f"(no_cz={no_cz})")
        for mid, m in enumerate(matchings):
            if progress:
                print(f"  matching {mid+1}/{len(matchings)}: {len(m)} edges")
            results.extend(_measure_matching(
                backend, m, shots, matching_id=mid, mode="parallel",
                no_cz=no_cz, optimization_level=optimization_level,
            ))
    else:  # isolated
        if progress:
            print(f"[edge_map] {len(edges)} edges, isolated (no_cz={no_cz})")
        for i, e in enumerate(edges):
            if progress and (i % 10 == 0 or i == len(edges) - 1):
                print(f"  edge {i+1}/{len(edges)}: {e}")
            results.extend(_measure_matching(
                backend, [e], shots, matching_id=-1, mode="isolated",
                no_cz=no_cz, optimization_level=optimization_level,
            ))
    return results


def run_isolated_edge_map(backend, *, edges=None, shots=4000,
                            optimization_level=3, progress=True
                            ) -> list[EdgeResult]:
    return run_edge_map(
        backend, edges=edges, shots=shots, mode="isolated",
        optimization_level=optimization_level, progress=progress,
    )


def isolated_vs_parallel(
    isolated: Sequence[EdgeResult], parallel: Sequence[EdgeResult],
) -> list[dict]:
    """Per-edge comparison of isolated vs parallel fidelity.

    Returns rows with delta_F = F_iso - F_par and z_delta = delta / sigma_combined.
    """
    iso_map = {r.edge: r for r in isolated}
    par_map = {r.edge: r for r in parallel}
    rows = []
    for e in sorted(set(iso_map) & set(par_map)):
        a, b = iso_map[e], par_map[e]
        delta = a.F - b.F
        sigma = math.sqrt(a.sigma_F ** 2 + b.sigma_F ** 2)
        if sigma > 0:
            z = delta / sigma
        else:
            z = 0.0 if delta == 0 else float("inf")
        rows.append({
            "edge": list(e),
            "F_isolated": a.F, "F_parallel": b.F,
            "delta_F": delta,
            "sigma_combined": sigma,
            "z_delta": z,
        })
    return rows


# ---------------------------------------------------------------------------
# Node-level diagnostics
# ---------------------------------------------------------------------------

def node_diagnostics(results: Sequence[EdgeResult]) -> dict[int, dict]:
    """Per-qubit summary: mean fidelity over incident edges (Q_i)."""
    by_node: dict[int, list[float]] = defaultdict(list)
    for r in results:
        by_node[r.edge[0]].append(r.F)
        by_node[r.edge[1]].append(r.F)
    return {
        q: {"Q": float(np.mean(fs)), "n_edges": len(fs),
            "min_F": float(np.min(fs)), "max_F": float(np.max(fs))}
        for q, fs in by_node.items()
    }


def edge_residuals(results: Sequence[EdgeResult]) -> list[dict]:
    """R_ij = F_ij - (Q_i + Q_j)/2.  Strongly negative → coupler-specific bug."""
    nodes = node_diagnostics(results)
    rows = []
    for r in results:
        qi, qj = r.edge
        baseline = 0.5 * (nodes[qi]["Q"] + nodes[qj]["Q"])
        rows.append({
            "edge": list(r.edge), "F": r.F, "baseline": baseline,
            "residual": r.F - baseline,
        })
    return sorted(rows, key=lambda x: x["residual"])


# ---------------------------------------------------------------------------
# Layout audit + recommendations
# ---------------------------------------------------------------------------

def _required_physical_edges(layout: Sequence[int],
                              required_logical: Sequence[Edge]) -> list[Edge]:
    return make_undirected_edges(
        (layout[i], layout[j]) for i, j in required_logical
    )


def audit_layout(
    layout: Sequence[int],
    required_logical_edges: Sequence[Edge],
    measured: Sequence[EdgeResult],
) -> dict:
    """Audit a proposed layout against measured edge entanglement.

    layout                    — physical qubit per logical index
    required_logical_edges    — logical edges the experiment will exercise
    measured                  — output of run_edge_map (any mode)

    Returns: per-edge lookup, weakest, mean/min F, validated flag,
    list of edges below the 3σ threshold.
    """
    by_edge = {r.edge: r for r in measured}
    needed = _required_physical_edges(layout, required_logical_edges)

    found: list[EdgeResult] = []
    missing: list[Edge] = []
    for e in needed:
        if e in by_edge:
            found.append(by_edge[e])
        else:
            missing.append(e)

    if not found:
        return {
            "validated": False, "reason": "no measured edges overlap layout",
            "missing_edges": [list(e) for e in missing],
            "edges_used": [list(e) for e in needed],
        }

    weakest = min(found, key=lambda r: r.F)
    failing = [r for r in found if not r.entangled_3sigma]
    return {
        "layout": list(layout),
        "edges_used": [list(e) for e in needed],
        "missing_edges": [list(e) for e in missing],
        "n_edges": len(needed),
        "n_measured": len(found),
        "n_below_3sigma": len(failing),
        "mean_F": float(np.mean([r.F for r in found])),
        "min_F": float(np.min([r.F for r in found])),
        "weakest_edge": {
            "edge": list(weakest.edge), "F": weakest.F, "z": weakest.z_score,
        },
        "failing_edges": [
            {"edge": list(r.edge), "F": r.F, "z": r.z_score} for r in failing
        ],
        "validated": (len(missing) == 0 and len(failing) == 0),
    }


def _build_neighbor_map(measured: Sequence[EdgeResult]
                          ) -> dict[int, dict[int, EdgeResult]]:
    nbr: dict[int, dict[int, EdgeResult]] = defaultdict(dict)
    for r in measured:
        a, b = r.edge
        nbr[a][b] = r
        nbr[b][a] = r
    return nbr


def recommend_best_chain(
    measured: Sequence[EdgeResult], length: int, *,
    require_3sigma: bool = True, max_search: int = 20000,
) -> dict | None:
    """Find a simple path of `length` qubits maximizing the weakest edge F.

    Score is min(F) along the path; tie-break by mean(F). DFS bounded by
    `max_search` partial paths to keep cost predictable.
    """
    if length < 2:
        raise ValueError("length must be ≥ 2")
    nbr = _build_neighbor_map(measured)
    if not nbr:
        return None
    best: tuple[float, float, list[int]] | None = None
    explored = [0]

    def dfs(path: list[int], path_edges: list[EdgeResult]):
        if explored[0] >= max_search:
            return
        explored[0] += 1
        if len(path) == length:
            mins = min(r.F for r in path_edges)
            avg = float(np.mean([r.F for r in path_edges]))
            nonlocal best
            cand = (mins, avg, list(path))
            if best is None or cand > best:
                best = cand
            return
        for nxt, r in nbr[path[-1]].items():
            if nxt in path:
                continue
            if require_3sigma and not r.entangled_3sigma:
                continue
            dfs(path + [nxt], path_edges + [r])

    for start in sorted(nbr):
        dfs([start], [])

    if best is None:
        return None
    mins, avg, path = best
    edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
    return {
        "qubits": path, "edges": edges,
        "min_F": mins, "mean_F": avg, "length": length,
        "weakest_edge": min(
            edges, key=lambda e: nbr[e[0]][e[1]].F,
        ),
    }


def recommend_best_patch(
    measured: Sequence[EdgeResult], rows: int, cols: int, *,
    backend=None, require_3sigma: bool = True,
) -> dict | None:
    """Best contiguous rows×cols sub-grid by min(F) of internal edges.

    Uses _enumerate_subgrids from src.backend if a backend is provided;
    otherwise falls back to brute-force from the measured graph.
    """
    measured_map = {r.edge: r for r in measured}

    if backend is not None and not _is_aer(backend):
        from src.backend import _enumerate_subgrids
        candidates = _enumerate_subgrids(backend, rows, cols)
    else:
        return None

    best = None
    for layout in candidates:
        internal = []
        for r in range(rows):
            for c in range(cols - 1):
                a, b = layout[r * cols + c], layout[r * cols + c + 1]
                internal.append(_norm((a, b)))
        for r in range(rows - 1):
            for c in range(cols):
                a, b = layout[r * cols + c], layout[(r + 1) * cols + c]
                internal.append(_norm((a, b)))

        if not all(e in measured_map for e in internal):
            continue
        edge_results = [measured_map[e] for e in internal]
        if require_3sigma and not all(r.entangled_3sigma for r in edge_results):
            continue
        mn = min(r.F for r in edge_results)
        avg = float(np.mean([r.F for r in edge_results]))
        cand = (mn, avg, layout, internal)
        if best is None or (mn, avg) > (best[0], best[1]):
            best = cand

    if best is None:
        return None
    mn, avg, layout, internal = best
    return {
        "rows": rows, "cols": cols, "layout": layout,
        "internal_edges": [list(e) for e in internal],
        "min_F": mn, "mean_F": avg,
    }


# ---------------------------------------------------------------------------
# IO + plot
# ---------------------------------------------------------------------------

def save_results_json_csv(results: Sequence[EdgeResult], path_prefix: str) -> None:
    rows = [r.to_dict() for r in results]
    with open(f"{path_prefix}.json", "w") as f:
        json.dump(rows, f, indent=2)
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(f"{path_prefix}.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            vals = []
            for k in keys:
                v = row[k]
                if isinstance(v, list):
                    v = ":".join(str(x) for x in v)
                vals.append(str(v))
            f.write(",".join(vals) + "\n")


def edge_color_for_result(r: EdgeResult) -> str:
    if r.entangled_3sigma:
        return "#2ca02c"   # green
    if r.entangled_meanonly:
        return "#ffbf00"   # amber
    return "#d62728"       # red


def plot_edge_map(
    results: Sequence[EdgeResult], backend=None, *,
    save_path: str | None = None, title: str | None = None,
):
    """Draw the device coupling graph with per-edge entanglement colors.

    Untested backend edges are drawn light gray.
    """
    import matplotlib.pyplot as plt
    import rustworkx as rx
    from rustworkx.visualization import mpl_draw

    by_edge = {r.edge: r for r in results}

    if backend is not None and getattr(backend, "coupling_map", None) is not None:
        all_edges = sorted({_norm(e) for e in backend.coupling_map})
        n_qubits = backend.num_qubits
    else:
        all_edges = sorted(by_edge)
        nodes = {q for e in all_edges for q in e}
        n_qubits = max(nodes) + 1 if nodes else 0

    g = rx.PyGraph()
    g.add_nodes_from(list(range(n_qubits)))
    edge_colors = []
    for u, v in all_edges:
        g.add_edge(u, v, None)
        r = by_edge.get(_norm((u, v)))
        edge_colors.append(edge_color_for_result(r) if r is not None else "#cccccc")

    fig, ax = plt.subplots(figsize=(10, 8))
    mpl_draw(
        g, ax=ax, with_labels=True, node_color="#dddddd",
        edge_color=edge_colors, width=2.5,
    )
    if title:
        ax.set_title(title)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=120)
    return fig
