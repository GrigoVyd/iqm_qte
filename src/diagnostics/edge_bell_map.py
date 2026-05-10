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


def plot_bottleneck_map(
    results: Sequence[EdgeResult],
    backend,
    *,
    save_path: str | None = None,
    title: str | None = None,
    excluded_qubits: Iterable[int] | None = None,
    excluded_edges: Iterable[Edge] | None = None,
    figsize: tuple[float, float] = (13, 13),
    annotate_qubits: bool = True,
    show_weakest: bool = True,
):
    """Polished single-chip "entanglement bottleneck" view.

    Draws the device's exact diamond-grid layout, with every measured CZ pair
    coloured by its graph-state fidelity F_ij (RdYlGn from 0.5 to 1.0 — the
    bottom of the colorbar matches the entanglement threshold).

    Highlights:
      - the weakest measured edge (dashed black overlay);
      - the weakest qubit by mean incident F (black double ring);
      - excluded couplers / qubits (dashed light-gray + open circles).

    Title gets a one-line statistical summary; bottom legend explains the
    three accent styles.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize
    from src.visualization import _device_layout

    pos = _device_layout(backend)
    nq = backend.num_qubits
    by_edge = {_norm(r.edge): r for r in results}
    measured_edges = set(by_edge.keys())
    all_edges = {_norm(e) for e in backend.coupling_map}
    expl_excl_edges = {_norm(e) for e in (excluded_edges or [])}
    excluded_set_q = set(excluded_qubits or ())

    # Edges that aren't in today's measurement set count as excluded
    not_measured = (all_edges - measured_edges) | expl_excl_edges

    # Per-qubit Q = mean F over incident *measured* edges.
    incidence: dict[int, list[float]] = {}
    for (a, b), r in by_edge.items():
        incidence.setdefault(a, []).append(r.F)
        incidence.setdefault(b, []).append(r.F)
    q_score = {q: float(np.mean(fs)) for q, fs in incidence.items()}

    # Discover weakest edge & qubit
    weakest_edge = min(by_edge.values(), key=lambda r: r.F) if by_edge else None
    if q_score:
        weakest_qubit = min(q_score, key=lambda k: q_score[k])
    else:
        weakest_qubit = None

    # Qubits that have *no* measured incident edge → effectively excluded
    auto_excluded_q = {q for q in range(nq) if q not in incidence and q not in excluded_set_q}
    excluded_q_total = excluded_set_q | auto_excluded_q

    fig, ax = plt.subplots(figsize=figsize)
    cmap = plt.cm.RdYlGn
    norm = Normalize(vmin=0.5, vmax=1.0)

    # --- background: excluded couplers as dashed light-gray ---
    for e in not_measured:
        if e[0] not in pos or e[1] not in pos:
            continue
        x1, y1 = pos[e[0]]; x2, y2 = pos[e[1]]
        ax.plot([x1, x2], [y1, y2], color="#bdbdbd", linewidth=1.6,
                 linestyle="--", alpha=0.85, zorder=1)

    # --- measured edges, coloured by F ---
    for e, r in by_edge.items():
        x1, y1 = pos[e[0]]; x2, y2 = pos[e[1]]
        ax.plot([x1, x2], [y1, y2], color=cmap(norm(r.F)),
                 linewidth=4.5, alpha=0.95, zorder=2,
                 solid_capstyle="round")

    # --- weakest-edge overlay (dashed black) ---
    if show_weakest and weakest_edge is not None:
        a, b = weakest_edge.edge
        x1, y1 = pos[a]; x2, y2 = pos[b]
        ax.plot([x1, x2], [y1, y2], color="#1a1a1a", linewidth=4.0,
                 linestyle=(0, (4, 3)), zorder=4, solid_capstyle="round")

    # --- nodes ---
    R_MEAS = 0.34
    R_EXCL = 0.32
    for q in range(nq):
        if q not in pos:
            continue
        x, y = pos[q]
        if q in excluded_q_total:
            ax.add_patch(patches.Circle((x, y), radius=R_EXCL,
                                          facecolor="white",
                                          edgecolor="#bdbdbd",
                                          linewidth=1.2, linestyle="--",
                                          zorder=5))
            if annotate_qubits:
                ax.text(x, y, f"QB{q+1}", ha="center", va="center",
                         fontsize=7.5, color="#9e9e9e", zorder=6)
            continue
        Q = q_score.get(q, float("nan"))
        col = cmap(norm(Q)) if not math.isnan(Q) else (0.85, 0.85, 0.85, 1.0)
        ax.add_patch(patches.Circle((x, y), radius=R_MEAS,
                                      facecolor=col,
                                      edgecolor="white",
                                      linewidth=1.5, zorder=5))
        if annotate_qubits:
            txt_col = "white" if (Q < 0.83 or Q > 0.93) else "black"
            ax.text(x, y, f"QB{q+1}", ha="center", va="center",
                     fontsize=8.5, color=txt_col, fontweight="bold", zorder=6)

    # --- weakest-qubit overlay (black double ring) ---
    if show_weakest and weakest_qubit is not None and weakest_qubit in pos:
        x, y = pos[weakest_qubit]
        ax.add_patch(patches.Circle((x, y), radius=R_MEAS + 0.18,
                                      facecolor="none", edgecolor="black",
                                      linewidth=2.0, zorder=7))

    # --- frame ---
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_facecolor("#fbfbfd")
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    pad = 1.2
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)

    # --- title + subtitle ---
    n_meas = len(by_edge); n_total = len(all_edges)
    if results:
        Fs = [r.F for r in results]
        f_min, f_mean = min(Fs), float(np.mean(Fs))
        weakest = weakest_edge
        we_str = (f"min F = {f_min:.3f} on QB{weakest.edge[0]+1}-QB{weakest.edge[1]+1}"
                  if weakest else f"min F = {f_min:.3f}")
        wq_str = (f"weakest qubit QB{weakest_qubit+1} (Q={q_score[weakest_qubit]:.3f})"
                  if weakest_qubit is not None else "")
        shots_str = (f" · {results[0].shots} shots/circuit"
                     if getattr(results[0], "shots", None) else "")
        subtitle = (f"{n_meas}/{n_total} edges measured{shots_str} · "
                    f"mean F = {f_mean:.3f}, {we_str}, {wq_str}")
    else:
        subtitle = "no measurements"
    full_title = title or "Entanglement bottleneck map"
    ax.set_title(f"{full_title}\n{subtitle}", fontsize=13,
                  fontweight="bold", pad=14)

    # --- colorbar (right) ---
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.02, shrink=0.85)
    cb.set_label("graph-state fidelity F  (separable bound = 0.5)",
                  fontsize=11)
    cb.ax.tick_params(labelsize=10)

    # --- bottom legend ---
    legend_handles = []
    if show_weakest and weakest_edge is not None:
        we = weakest_edge
        legend_handles.append(Line2D([0], [0], color="#1a1a1a", linewidth=4.0,
                                       linestyle=(0, (4, 3)),
                                       label=f"weakest edge (QB{we.edge[0]+1}–QB{we.edge[1]+1})"))
    if show_weakest and weakest_qubit is not None:
        legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                       markerfacecolor="white",
                                       markeredgecolor="black", markeredgewidth=2,
                                       markersize=14,
                                       label=f"weakest qubit (QB{weakest_qubit+1})"))
    legend_handles.append(Line2D([0], [0], color="#bdbdbd", linewidth=1.8,
                                   linestyle="--",
                                   label="excluded coupler (not in today's calibration)"))
    ax.legend(handles=legend_handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.02), ncol=len(legend_handles),
              frameon=False, fontsize=10)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=130, bbox_inches="tight")
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
