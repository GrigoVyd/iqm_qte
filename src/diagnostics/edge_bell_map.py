"""Per-edge graph-state fidelity map.

Measures F_ij = (1 + <X_i Z_j> + <Z_i X_j> + <Y_i Y_j>) / 4 on every native
CZ edge of the device. F_ij > 1/2 certifies entanglement on that edge
(|G_ij> = CZ |+>|+> is locally equivalent to a Bell state, so the same
fidelity-witness argument as for Bell pairs applies).

Edges are grouped into matchings via greedy edge coloring so each matching
runs in a single circuit (parallel CZ layer + parallel measurement).
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from typing import Iterable, Sequence

import numpy as np
from qiskit import QuantumCircuit, transpile

Edge = tuple[int, int]
BASES = ("XZ", "ZX", "YY")


def _norm(edge: Iterable[int]) -> Edge:
    a, b = edge
    return (a, b) if a < b else (b, a)


def get_coupling_edges(backend) -> list[Edge]:
    cmap = getattr(backend, "coupling_map", None)
    if cmap is None or len(list(cmap)) == 0:
        raise ValueError("Backend has no coupling_map; pass `edges` explicitly.")
    return sorted({_norm(e) for e in cmap})


def edge_color_matchings(edges: Sequence[Edge]) -> list[list[Edge]]:
    """Greedy edge coloring -> list of disjoint matchings."""
    edges = sorted({_norm(e) for e in edges})
    incident: dict[int, list[int]] = defaultdict(list)
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


def _basis_rotation(qc: QuantumCircuit, q: int, pauli: str) -> None:
    if pauli == "X":
        qc.h(q)
    elif pauli == "Y":
        qc.sdg(q); qc.h(q)
    elif pauli != "Z":
        raise ValueError(pauli)


def _build_matching_circuit(matching: Sequence[Edge], basis: str
                            ) -> tuple[QuantumCircuit, list[int], dict[int, int]]:
    p1, p2 = basis[0], basis[1]
    phys = sorted({q for e in matching for q in e})
    p2l = {p: k for k, p in enumerate(phys)}
    n = len(phys)
    qc = QuantumCircuit(n, n)
    qc.h(range(n))
    for a, b in matching:
        qc.cz(p2l[a], p2l[b])
    qc.barrier()
    for a, b in matching:
        _basis_rotation(qc, p2l[a], p1)
        _basis_rotation(qc, p2l[b], p2)
    qc.measure(range(n), range(n))
    return qc, phys, p2l


def _bit_pm1(b: str) -> int:
    return 1 if b == "0" else -1


def expectation_from_counts(counts: dict[str, int], li: int, lj: int, n: int) -> float:
    if not counts:
        return 0.0
    total = 0
    weighted = 0
    for bs, cnt in counts.items():
        bits = bs.replace(" ", "")
        v1 = _bit_pm1(bits[n - 1 - li])
        v2 = _bit_pm1(bits[n - 1 - lj])
        weighted += v1 * v2 * cnt
        total += cnt
    return weighted / total if total else 0.0


def _stderr(e: float, shots: int) -> float:
    return math.sqrt(max(0.0, 1.0 - e * e) / max(1, shots))


def compute_edge_fidelity(e_xz: float, e_zx: float, e_yy: float, shots: int) -> dict:
    F = (1.0 + e_xz + e_zx + e_yy) / 4.0
    sx, sz, sy = (_stderr(x, shots) for x in (e_xz, e_zx, e_yy))
    sigma = 0.25 * math.sqrt(sx * sx + sz * sz + sy * sy)
    z = (F - 0.5) / sigma if sigma > 0 else float("inf")
    return {"F": F, "sigma_F": sigma, "z_score": z,
            "entangled_3sigma": (F - 3 * sigma) > 0.5,
            "entangled_meanonly": F > 0.5}


@dataclass
class EdgeResult:
    edge: Edge
    matching_id: int
    shots: int
    e_XZ: float
    e_ZX: float
    e_YY: float
    F: float
    sigma_F: float
    z_score: float
    entangled_3sigma: bool
    entangled_meanonly: bool

    def to_dict(self) -> dict:
        d = asdict(self)
        d["edge"] = list(self.edge)
        return d


def _is_aer(backend) -> bool:
    return backend.__class__.__name__ == "AerSimulator"


def _transpile(qc, backend, phys, opt):
    if _is_aer(backend):
        return transpile(qc, backend=backend, optimization_level=opt)
    return transpile(qc, backend=backend, optimization_level=opt, initial_layout=phys)


def _run_matching(backend, matching, shots, mid, opt) -> list[EdgeResult]:
    circuits, ptols, ns = [], [], []
    for basis in BASES:
        qc, phys, p2l = _build_matching_circuit(matching, basis)
        circuits.append(_transpile(qc, backend, phys, opt))
        ptols.append(p2l)
        ns.append(len(phys))
    job = backend.run(circuits, shots=shots)
    raw = job.result().get_counts()
    if isinstance(raw, dict):
        raw = [raw]
    counts_by_basis = dict(zip(BASES, raw))
    out = []
    for a, b in matching:
        p2l, n = ptols[0], ns[0]
        li, lj = p2l[a], p2l[b]
        e_xz = expectation_from_counts(counts_by_basis["XZ"], li, lj, n)
        e_zx = expectation_from_counts(counts_by_basis["ZX"], li, lj, n)
        e_yy = expectation_from_counts(counts_by_basis["YY"], li, lj, n)
        f = compute_edge_fidelity(e_xz, e_zx, e_yy, shots)
        out.append(EdgeResult(edge=(a, b), matching_id=mid, shots=shots,
                              e_XZ=e_xz, e_ZX=e_zx, e_YY=e_yy, **f))
    return out


def run_edge_map(backend, *, edges: Sequence[Edge] | None = None,
                  shots: int = 4000, optimization_level: int = 3,
                  progress: bool = True) -> list[EdgeResult]:
    """Measure F_ij on every requested edge in parallel matchings."""
    if edges is None:
        edges = get_coupling_edges(backend)
    edges = sorted({_norm(e) for e in edges})
    if not edges:
        return []
    matchings = edge_color_matchings(edges)
    if progress:
        print(f"[edge_map] {len(edges)} edges -> {len(matchings)} matchings")
    out: list[EdgeResult] = []
    for mid, m in enumerate(matchings):
        if progress:
            print(f"  matching {mid+1}/{len(matchings)}: {len(m)} edges")
        out.extend(_run_matching(backend, m, shots, mid, optimization_level))
    return out


def fidelity_map(results: Sequence[EdgeResult]) -> dict[Edge, float]:
    """{(a,b) -> F_ij} dict for downstream consumers (e.g. tree selection)."""
    return {_norm(r.edge): r.F for r in results}


def plot_edge_map(results: Sequence[EdgeResult], backend=None, *,
                   save_path: str | None = None, title: str | None = None,
                   layout: dict | None = None):
    """Color-coded coupling-graph view: green=3sigma entangled, amber=mean only,
    red=not entangled, gray=untested."""
    import matplotlib.pyplot as plt
    import rustworkx as rx
    from rustworkx.visualization import mpl_draw

    by_edge = {_norm(r.edge): r for r in results}
    if backend is not None and getattr(backend, "coupling_map", None) is not None:
        all_edges = sorted({_norm(e) for e in backend.coupling_map})
        nq = backend.num_qubits
    else:
        all_edges = sorted(by_edge)
        nq = (max({q for e in all_edges for q in e}) + 1) if all_edges else 0

    g = rx.PyGraph()
    g.add_nodes_from(list(range(nq)))
    colors = []
    for e in all_edges:
        g.add_edge(e[0], e[1], None)
        r = by_edge.get(e)
        if r is None:
            colors.append("#cccccc")
        elif r.entangled_3sigma:
            colors.append("#2ca02c")
        elif r.entangled_meanonly:
            colors.append("#ffbf00")
        else:
            colors.append("#d62728")

    fig, ax = plt.subplots(figsize=(11, 8))
    kw = {"ax": ax, "with_labels": True, "node_color": "#dddddd",
          "edge_color": colors, "width": 2.5}
    if layout is not None:
        kw["pos"] = layout
    mpl_draw(g, **kw)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=120, bbox_inches="tight")
    return fig


def fidelity_heatmap(results: Sequence[EdgeResult], backend, *,
                      save_path: str | None = None, title: str | None = None):
    """N x N symmetric F-matrix heatmap (NaN where edge is not native)."""
    import matplotlib.pyplot as plt
    nq = backend.num_qubits
    M = np.full((nq, nq), np.nan)
    for r in results:
        a, b = r.edge
        M[a, b] = r.F
        M[b, a] = r.F
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(M, cmap="RdYlGn", vmin=0.25, vmax=1.0)
    ax.set_xlabel("Qubit j"); ax.set_ylabel("Qubit i")
    ax.set_title(title or "Measured pairwise graph-state fidelity F_ij")
    plt.colorbar(im, ax=ax, label="F_ij  (>0.5 = entangled)")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=120, bbox_inches="tight")
    return fig
