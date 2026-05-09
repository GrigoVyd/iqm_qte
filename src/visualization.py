"""IQM-Resonance-style topology visualization.

Colors qubits by a chosen metric (T1, T2, readout fidelity, 1Q fidelity).
Highlights a chosen subset (the tree we picked) and the edges of the chosen
graph state.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import numpy as np
import rustworkx as rx


# Hardcoded IQM Emerald qubit positions (54 qubits) — derived from a BFS over
# the device coupling graph using diamond-grid offsets, with collision-aware
# backtracking. Yields the rotated-square-lattice layout you see on the IQM
# Resonance dashboard. Coordinates are in arbitrary units.
def _device_layout(backend) -> dict[int, tuple[float, float]]:
    """2D coordinates per qubit, computed once via a smart BFS that maps the
    bipartite topology to a diamond grid. Each placed qubit's neighbours are
    tried in 4 directions (NE, NW, SE, SW); we pick the direction that
    minimises future placement conflicts (most-degrees-first heuristic).
    """
    from collections import deque
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

    DIRS = [(1, 1), (-1, 1), (-1, -1), (1, -1)]   # NE, NW, SW, SE

    # Start from a high-degree central qubit so the layout grows in all directions.
    start = max(adj, key=lambda x: len(adj[x]))
    pos: dict[int, tuple[int, int]] = {start: (0, 0)}
    used: set[tuple[int, int]] = {(0, 0)}
    q = deque([start])

    while q:
        u = q.popleft()
        ux, uy = pos[u]
        # Place each unplaced neighbour: pick first free diamond offset.
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
                # All 4 diamond slots occupied — extend further out.
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

    # Place any qubits not in the BFS-reached component (isolated/dangling).
    placed_set = set(pos.keys())
    isolated = [i for i in range(backend.num_qubits) if i not in placed_set]
    # Find the bounding box of placed qubits
    if placed_set:
        max_x = max(p[0] for p in pos.values())
        for i, q_node in enumerate(isolated):
            pos[q_node] = (max_x + 4 + (i % 2), -2 + 2 * (i // 2))
    else:
        for i, q_node in enumerate(isolated):
            pos[q_node] = (i, 0)

    # Convert to floats
    return {k: (float(v[0]), float(v[1])) for k, v in pos.items()}


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
):
    """Plot device topology in the IQM dashboard style.

    metrics: per-qubit metrics from get_qubit_metrics; coloring uses metrics[i][color_by].
    cz_fidelities: per-pair fidelity, used to color edge diamonds.
    highlight_qubits: bigger marker + bold border.
    highlight_edges: (a, b) physical edges drawn thicker and in red.
    """
    if pos is None:
        pos = _device_layout(backend)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
    else:
        fig = ax.figure

    n = backend.num_qubits
    highlight_qubits = set(highlight_qubits or [])
    highlight_edge_set = {(min(a, b), max(a, b)) for a, b in (highlight_edges or [])}

    # --- edges ---
    for a, b in backend.coupling_map:
        if a >= b:  # undirected, plot once
            continue
        x1, y1 = pos[a]
        x2, y2 = pos[b]
        e = (a, b)

        # Edge color based on CZ fidelity (or default gray)
        if cz_fidelities and e in cz_fidelities:
            cz_f = cz_fidelities[e]
            edge_color = plt.cm.RdYlGn((cz_f - 0.85) / 0.15)  # 0.85→red, 1.0→green
        else:
            edge_color = "#cccccc"

        is_highlight = e in highlight_edge_set
        ax.plot([x1, x2], [y1, y2],
                color="#e84545" if is_highlight else edge_color,
                linewidth=4.5 if is_highlight else 1.3,
                zorder=2 if is_highlight else 1,
                alpha=1.0 if is_highlight else 0.55)

        # Diamond mid-marker (IQM dashboard style)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        diamond_color = edge_color if not is_highlight else "#e84545"
        diamond = patches.RegularPolygon(
            (mx, my), 4, radius=0.18, orientation=np.pi/4,
            facecolor=diamond_color,
            edgecolor="black" if is_highlight else "none",
            linewidth=1.2 if is_highlight else 0,
            zorder=3 if is_highlight else 2)
        ax.add_patch(diamond)

    # --- nodes ---
    for i in range(n):
        x, y = pos[i]
        # Color by chosen metric
        color = "#cccccc"
        if metrics and i in metrics:
            v = metrics[i].get(color_by)
            if v is not None:
                # Map color_by-specific range to [0,1]
                if color_by == "readout_fidelity":
                    norm = (v - 0.85) / 0.14   # 0.85 red → 0.99 green
                elif color_by == "one_q_fidelity":
                    norm = (v - 0.99) / 0.01
                elif color_by == "t1":
                    norm = min((v * 1e6 - 5) / 95, 1.0)   # 5µs red → 100µs green
                elif color_by == "t2":
                    norm = min((v * 1e6 - 1) / 49, 1.0)   # 1µs red → 50µs green
                else:
                    norm = min(max(float(v), 0), 1)
                norm = max(0.0, min(1.0, norm))
                color = plt.cm.viridis(norm)

        is_highlight = i in highlight_qubits
        size = 0.42 if is_highlight else 0.32
        circle = patches.Circle(
            (x, y), radius=size,
            facecolor=color,
            edgecolor="#e84545" if is_highlight else "white",
            linewidth=3 if is_highlight else 1.5,
            zorder=5)
        ax.add_patch(circle)

        if show_labels:
            ax.text(x, y, f"QB{i+1}", ha="center", va="center",
                    fontsize=7 if is_highlight else 6.5,
                    color="white", fontweight="bold", zorder=6)

    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ["top", "right", "left", "bottom"]:
        ax.spines[s].set_visible(False)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")

    # Legend
    if highlight_qubits or highlight_edge_set:
        legend_elems = []
        if highlight_qubits:
            legend_elems.append(Line2D([0], [0], marker="o", color="w",
                markerfacecolor="#888", markeredgecolor="#e84545",
                markeredgewidth=2.5, markersize=12, label="Selected qubit"))
        if highlight_edge_set:
            legend_elems.append(Line2D([0], [0], color="#e84545", linewidth=3,
                                        label="Tree edge (CZ gate)"))
        ax.legend(handles=legend_elems, loc="upper right", fontsize=9, frameon=True)

    return fig, ax
