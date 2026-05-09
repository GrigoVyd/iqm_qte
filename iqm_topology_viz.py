"""IQM-Resonance-style topology visualization.

Two main entry points:
  - plot_device_topology(...)       static plot, optionally with a tree highlighted
  - animate_tree_growth(...)        GIF showing how the chosen tree scales with n
"""

from __future__ import annotations

from collections import deque

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import numpy as np


# ---------- layout ----------

# Exact IQM Emerald (Aphrodite) layout, transcribed from the ASCII art in
# iqm.cirq_iqm.devices.aphrodite.Aphrodite docstring. Coordinates are in the
# diamond-grid convention used on the Resonance dashboard: each unit step in
# x or y corresponds to one diagonal step on the chip. Qubit i (0-indexed)
# maps to physical name QB{i+1}.
EMERALD_POSITIONS: dict[int, tuple[float, float]] = {
    # row 10 (top)
    53: (1, 10), 50: (3, 10), 45: (5, 10), 38: (7, 10),
    # row 9
    52: (0, 9), 49: (2, 9), 44: (4, 9), 37: (6, 9), 30: (8, 9),
    # row 8
    51: (-1, 8), 48: (1, 8), 43: (3, 8), 36: (5, 8), 29: (7, 8),
    # row 7
    47: (0, 7), 42: (2, 7), 35: (4, 7), 28: (6, 7), 21: (8, 7),
    # row 6
    46: (-1, 6), 41: (1, 6), 34: (3, 6), 27: (5, 6), 20: (7, 6), 13: (9, 6),
    # row 5
    40: (0, 5), 33: (2, 5), 26: (4, 5), 19: (6, 5), 12: (8, 5),
    # row 4
    39: (-1, 4), 32: (1, 4), 25: (3, 4), 18: (5, 4), 11: (7, 4), 6: (9, 4),
    # row 3
    31: (0, 3), 24: (2, 3), 17: (4, 3), 10: (6, 3), 5: (8, 3),
    # row 2
    23: (1, 2), 16: (3, 2), 9: (5, 2), 4: (7, 2), 1: (9, 2),
    # row 1
    22: (0, 1), 15: (2, 1), 8: (4, 1), 3: (6, 1), 0: (8, 1),
    # row 0 (bottom)
    14: (1, 0), 7: (3, 0), 2: (5, 0),
}


def _device_layout(backend) -> dict[int, tuple[float, float]]:
    """Hardcoded IQM Emerald layout if the backend has 54 qubits.
    Falls back to BFS-diamond walk for unknown topologies.
    """
    if backend.num_qubits == 54:
        return {i: (float(x), float(y)) for i, (x, y) in EMERALD_POSITIONS.items()}
    # Fallback: BFS-diamond walk
    coupling: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for a, b in backend.coupling_map:
        e = (min(a, b), max(a, b))
        if e in seen:
            continue
        seen.add(e)
        coupling.append(e)

    adj: dict[int, set[int]] = {}
    for a, b in coupling:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    if not adj:
        return {i: (float(i), 0.0) for i in range(backend.num_qubits)}

    DIRS = [(1, 1), (-1, 1), (-1, -1), (1, -1)]

    start = max(adj, key=lambda x: len(adj[x]))
    pos: dict[int, tuple[int, int]] = {start: (0, 0)}
    used: set[tuple[int, int]] = {(0, 0)}
    q = deque([start])

    while q:
        u = q.popleft()
        ux, uy = pos[u]
        for v in sorted(adj[u], key=lambda x: -len(adj.get(x, []))):
            if v in pos:
                continue
            for dx, dy in DIRS:
                cand = (ux + dx, uy + dy)
                if cand in used:
                    continue
                pos[v] = cand
                used.add(cand)
                q.append(v)
                break
            else:
                for r in range(2, 8):
                    placed = False
                    for dx, dy in [(r, 0), (-r, 0), (0, r), (0, -r),
                                    (r, r), (-r, r), (r, -r), (-r, -r)]:
                        cand = (ux + dx, uy + dy)
                        if cand not in used:
                            pos[v] = cand
                            used.add(cand)
                            q.append(v)
                            placed = True
                            break
                    if placed:
                        break

    placed_set = set(pos.keys())
    isolated = [i for i in range(backend.num_qubits) if i not in placed_set]
    if placed_set:
        max_x = max(p[0] for p in pos.values())
        for i, q_node in enumerate(isolated):
            pos[q_node] = (max_x + 4 + (i % 2), -2 + 2 * (i // 2))
    else:
        for i, q_node in enumerate(isolated):
            pos[q_node] = (i, 0)

    return {k: (float(v[0]), float(v[1])) for k, v in pos.items()}


# ---------- color helpers ----------

def _bipartite_coloring(backend) -> dict[int, int]:
    """2-color the device coupling graph (BFS). Returns {qubit_idx: 0 or 1}."""
    adj: dict[int, set[int]] = {}
    for a, b in backend.coupling_map:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    color: dict[int, int] = {}
    for start in range(backend.num_qubits):
        if start in color:
            continue
        if start not in adj:
            color[start] = 0
            continue
        color[start] = 0
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj.get(u, []):
                if v not in color:
                    color[v] = 1 - color[u]
                    q.append(v)
    return color


def _color_for_metric(value: float | None, metric: str) -> tuple:
    """Map a metric value to a viridis color. Returns gray if missing."""
    if value is None:
        return (0.8, 0.8, 0.8, 1.0)
    if metric == "readout_fidelity":
        norm = (value - 0.85) / 0.14
    elif metric == "one_q_fidelity":
        norm = (value - 0.99) / 0.01
    elif metric == "t1":
        norm = min((value * 1e6 - 5) / 95, 1.0)
    elif metric == "t2":
        norm = min((value * 1e6 - 1) / 49, 1.0)
    else:
        norm = float(value)
    norm = max(0.0, min(1.0, norm))
    return plt.cm.viridis(norm)


def _cz_color(fid: float | None) -> tuple:
    if fid is None:
        return (0.78, 0.78, 0.78, 1.0)
    return plt.cm.RdYlGn(max(0.0, min(1.0, (fid - 0.85) / 0.15)))


# ---------- main static plot ----------

def plot_device_topology(
    backend,
    metrics: dict[int, dict] | None = None,
    cz_fidelities: dict[tuple[int, int], float] | None = None,
    color_by: str = "readout_fidelity",
    highlight_qubits: list[int] | None = None,
    highlight_edges: list[tuple[int, int]] | None = None,
    title: str | None = None,
    ax=None,
    show_labels: bool = True,
    pos: dict[int, tuple[float, float]] | None = None,
    spotlight: bool = True,
):
    """Render the device topology with an optional sub-tree highlighted.

    spotlight: if True, dim non-selected qubits/edges so the tree pops.
    """
    if pos is None:
        pos = _device_layout(backend)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
    else:
        fig = ax.figure

    n = backend.num_qubits
    hi_qubits = set(highlight_qubits or [])
    hi_edges = {(min(a, b), max(a, b)) for a, b in (highlight_edges or [])}
    has_highlights = bool(hi_qubits or hi_edges)

    bg_alpha = 0.22 if (spotlight and has_highlights) else 0.65

    bp_color = _bipartite_coloring(backend) if color_by == "bipartite" else {}
    BP_COLORS = {0: '#3aa56b', 1: '#7b3aa5'}

    # --- Layer 1: dim background edges ---
    for a, b in backend.coupling_map:
        if a >= b:
            continue
        x1, y1 = pos[a]
        x2, y2 = pos[b]
        e = (a, b)
        if e in hi_edges:
            continue
        col = _cz_color(cz_fidelities.get(e) if cz_fidelities else None)
        ax.plot([x1, x2], [y1, y2], color=col, linewidth=1.5,
                alpha=bg_alpha, zorder=1)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        d = patches.RegularPolygon((mx, my), 4, radius=0.16,
                                    orientation=np.pi/4,
                                    facecolor=col, edgecolor='none',
                                    alpha=bg_alpha, zorder=2)
        ax.add_patch(d)

    # --- Layer 2: dim background nodes ---
    for i in range(n):
        if i in hi_qubits:
            continue
        x, y = pos[i]
        if color_by == "bipartite":
            col = BP_COLORS[bp_color.get(i, 0)]
        else:
            col = _color_for_metric(
                metrics.get(i, {}).get(color_by) if metrics else None, color_by)
        circle = patches.Circle((x, y), radius=0.32, facecolor=col,
                                 edgecolor='white', linewidth=1.0,
                                 alpha=bg_alpha, zorder=3)
        ax.add_patch(circle)
        if show_labels:
            ax.text(x, y, f"QB{i+1}", ha='center', va='center',
                    fontsize=6.5, color='white', alpha=bg_alpha + 0.2,
                    fontweight='bold', zorder=4)

    # --- Layer 3: highlighted edges (route) ---
    for (a, b) in hi_edges:
        if a not in pos or b not in pos:
            continue
        x1, y1 = pos[a]
        x2, y2 = pos[b]
        ax.plot([x1, x2], [y1, y2], color='#fff4e0', linewidth=10,
                alpha=0.9, zorder=10, solid_capstyle='round')
        ax.plot([x1, x2], [y1, y2], color='#ff3b3b', linewidth=5.5,
                zorder=11, solid_capstyle='round')

    # --- Layer 4: highlighted nodes ---
    for i in hi_qubits:
        x, y = pos[i]
        if color_by == "bipartite":
            col = BP_COLORS[bp_color.get(i, 0)]
        else:
            col = _color_for_metric(
                metrics.get(i, {}).get(color_by) if metrics else None, color_by)
        ax.add_patch(patches.Circle((x, y), radius=0.55,
                                     facecolor='none', edgecolor='#ff3b3b',
                                     linewidth=4, zorder=12))
        ax.add_patch(patches.Circle((x, y), radius=0.42,
                                     facecolor=col, edgecolor='white',
                                     linewidth=1.8, zorder=13))
        if show_labels:
            ax.text(x, y, f"QB{i+1}", ha='center', va='center',
                    fontsize=8, color='white', fontweight='bold', zorder=14)

    # --- Frame ---
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ['top', 'right', 'left', 'bottom']:
        ax.spines[s].set_visible(False)
    ax.set_facecolor('#f5f5f7')

    if pos:
        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        pad = 1.2
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)

    if title:
        ax.set_title(title, fontsize=13, fontweight='bold', pad=12)

    if has_highlights:
        legend_elems = [
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor='#888', markeredgecolor='#ff3b3b',
                   markeredgewidth=2.5, markersize=14, label=f'Selected qubit ({len(hi_qubits)})'),
            Line2D([0], [0], color='#ff3b3b', linewidth=4,
                   label=f'Route edge ({len(hi_edges)})'),
        ]
        ax.legend(handles=legend_elems, loc='upper right', fontsize=9,
                  frameon=True, facecolor='white', framealpha=0.95)

    return fig, ax


# ---------- animation ----------

def animate_prim_steps(
    backend,
    target_n: int,
    metrics: dict[int, dict] | None = None,
    cz_fidelities: dict[tuple[int, int], float] | None = None,
    seed_qubit: int | None = None,
    excluded_qubits: set[int] | None = None,
    out_path: str = "gme_prim_animation.gif",
    color_by: str = "readout_fidelity",
    interval_ms: int = 700,
):
    """Animate Prim's algorithm building the spanning tree edge-by-edge."""
    import heapq
    import matplotlib.animation as animation

    pos = _device_layout(backend)
    excluded = excluded_qubits or set()

    nbrs: dict[int, list[tuple[int, float]]] = {}
    for (a, b), fid in (cz_fidelities or {}).items():
        if a in excluded or b in excluded:
            continue
        w = 1.0 - fid
        nbrs.setdefault(a, []).append((b, w))
        nbrs.setdefault(b, []).append((a, w))

    if seed_qubit is None:
        if metrics:
            seed_qubit = max(
                (i for i in nbrs.keys()),
                key=lambda i: metrics.get(i, {}).get(color_by, 0))
        else:
            seed_qubit = next(iter(nbrs))

    in_tree = {seed_qubit}
    tree_edges: list[tuple[int, int]] = []
    heap: list[tuple[float, int, int]] = []
    for nb, w in nbrs.get(seed_qubit, []):
        heapq.heappush(heap, (w, seed_qubit, nb))

    frames: list[tuple[set[int], list[tuple[int, int]], list[tuple[int, int]]]] = [
        (set(in_tree), [], [(u, v) for _, u, v in heap])
    ]

    while heap and len(in_tree) < target_n:
        w, u, v = heapq.heappop(heap)
        if v in in_tree:
            continue
        in_tree.add(v)
        tree_edges.append((u, v))
        for nb, w2 in nbrs.get(v, []):
            if nb not in in_tree:
                heapq.heappush(heap, (w2, v, nb))
        cands = [(uu, vv) for _, uu, vv in heap if vv not in in_tree]
        frames.append((set(in_tree), list(tree_edges), cands))

    fig, ax = plt.subplots(figsize=(11, 11))
    n_frames = len(frames)

    def draw_frame(i: int):
        ax.clear()
        in_tree_now, tree_edges_now, cand_edges = frames[i]
        plot_device_topology(
            backend, metrics=metrics, cz_fidelities=cz_fidelities,
            color_by=color_by,
            highlight_qubits=list(in_tree_now),
            highlight_edges=tree_edges_now,
            title=f"Prim's algorithm — step {i}/{n_frames-1}  "
                  f"(|tree| = {len(in_tree_now)} qubits, {len(tree_edges_now)} edges)",
            ax=ax, pos=pos, spotlight=True,
        )
        for u, v in cand_edges:
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            ax.plot([x1, x2], [y1, y2], color='#ffa726',
                    linewidth=2.5, linestyle='--', alpha=0.8, zorder=8)

    anim = animation.FuncAnimation(
        fig, draw_frame, frames=n_frames, interval=interval_ms, repeat=True,
    )
    anim.save(out_path, writer=animation.PillowWriter(
        fps=max(1, 1000 // interval_ms)))
    plt.close(fig)
    return out_path


def animate_tree_growth(
    backend,
    trees_by_n: dict[int, dict],
    metrics: dict[int, dict] | None = None,
    cz_fidelities: dict[tuple[int, int], float] | None = None,
    out_path: str = "gme_tree_growth.gif",
    color_by: str = "readout_fidelity",
    interval_ms: int = 1100,
    title_prefix: str = "IQM Emerald — optimal spanning tree",
):
    """Build a GIF showing the chosen tree at each n in ascending order."""
    import matplotlib.animation as animation

    pos = _device_layout(backend)
    sorted_ns = sorted(trees_by_n.keys())

    fig, ax = plt.subplots(figsize=(11, 11))

    def draw_frame(n_value: int):
        ax.clear()
        t = trees_by_n[n_value]
        plot_device_topology(
            backend, metrics=metrics, cz_fidelities=cz_fidelities,
            color_by=color_by,
            highlight_qubits=t['qubits'],
            highlight_edges=t['edges'],
            title=f'{title_prefix} — n = {n_value} qubits  (tree weight {t["weight"]:.3f})',
            ax=ax, pos=pos, show_labels=True,
            spotlight=True,
        )

    anim = animation.FuncAnimation(
        fig, draw_frame, frames=sorted_ns, interval=interval_ms, repeat=True,
    )
    anim.save(out_path, writer=animation.PillowWriter(fps=max(1, 1000 // interval_ms)))
    plt.close(fig)
    return out_path
