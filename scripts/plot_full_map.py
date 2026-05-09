"""Render the full Emerald entanglement bottleneck map.

Loads results/full_map.json (no hardware needed) and produces a publication-ish
heat-map of the chip:

  - layout via Kamada-Kawai on the coupling graph (preserves grid topology)
  - edges colored continuously by F (viridis from F=0.5 to F=1.0)
  - node size and outline color encode Q_i (mean F over incident edges)
  - bottleneck qubit/edge highlighted

Companion bar chart sorts edges by F, marks the 0.5 separability threshold and
the 3σ-certification line.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

JSON_PATH = Path("results/full_map.json")
OUT_PNG = Path("results/full_map_chip.png")
OUT_BAR_PNG = Path("results/full_map_edges_sorted.png")


def load_results():
    rows = json.loads(JSON_PATH.read_text())
    edges = {tuple(d["edge"]): d for d in rows}
    nodes = sorted({q for e in edges for q in e})
    # Q_i = mean F over incident edges
    incident = defaultdict(list)
    for (a, b), d in edges.items():
        incident[a].append(d["F"])
        incident[b].append(d["F"])
    Q = {q: float(np.mean(incident[q])) for q in nodes}
    return edges, nodes, Q


def make_layout(nodes, edges):
    """Kamada-Kawai on the coupling graph: preserves graph distances, so a
    square lattice comes out looking like a square lattice."""
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for (a, b) in edges:
        g.add_edge(a, b)
    # KK is deterministic up to global rotation/reflection — fine for a static
    # plot. Use BFS-based initial guess for better grid alignment.
    pos = nx.kamada_kawai_layout(g)
    return pos


def plot_chip(edges, nodes, Q, pos):
    fig, ax = plt.subplots(figsize=(13, 11))

    # Continuous edge colormap from F=0.5 (red) to F=1.0 (deep green)
    f_min_cmap, f_max_cmap = 0.5, 1.0
    cmap = mpl.colormaps["RdYlGn"]
    norm = mpl.colors.Normalize(vmin=f_min_cmap, vmax=f_max_cmap)

    # Draw edges
    for (a, b), d in edges.items():
        x0, y0 = pos[a]
        x1, y1 = pos[b]
        F = d["F"]
        col = cmap(norm(max(f_min_cmap, F)))
        # Width: thicker for higher F (encodes z indirectly through F too)
        w = 1.5 + 4.0 * max(0.0, (F - 0.5) / 0.5)
        ax.plot([x0, x1], [y0, y1], color=col, linewidth=w,
                solid_capstyle="round", zorder=2)

    # Identify weakest qubit and edge for annotation
    worst_q = min(Q.items(), key=lambda kv: kv[1])
    worst_e = min(edges.items(), key=lambda kv: kv[1]["F"])

    # Draw nodes — color by Q, size proportional to Q (worse qubits stand out)
    qs = np.array([Q[n] for n in nodes])
    xs = np.array([pos[n][0] for n in nodes])
    ys = np.array([pos[n][1] for n in nodes])
    sizes = 200 + 700 * (qs - qs.min()) / max(1e-9, qs.max() - qs.min())
    # Use the same colormap for nodes for visual consistency
    node_colors = cmap(norm(qs))
    sc = ax.scatter(xs, ys, s=sizes, c=node_colors, edgecolors="black",
                     linewidths=1.0, zorder=3)
    for n in nodes:
        x, y = pos[n]
        is_worst = (n == worst_q[0])
        ax.text(x, y, f"{n}", ha="center", va="center",
                fontsize=8 if not is_worst else 10,
                fontweight="bold" if is_worst else "normal",
                zorder=4, color="black" if Q[n] > 0.75 else "white")

    # Highlight worst edge with dashed black outline
    a, b = worst_e[0]
    x0, y0 = pos[a]; x1, y1 = pos[b]
    ax.plot([x0, x1], [y0, y1], color="black", linewidth=8,
            linestyle=(0, (4, 2)), alpha=0.55, zorder=1)

    # Highlight worst qubit with red ring
    wx, wy = pos[worst_q[0]]
    ax.scatter([wx], [wy], s=2000, facecolors="none",
                edgecolors="black", linewidths=2.5, zorder=5)

    # Colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("graph-state fidelity F  (separable bound = 0.5)",
                    fontsize=10)
    cbar.ax.axhline(0.5, color="black", linewidth=1.5)

    # Stats summary in title
    Fs = np.array([d["F"] for d in edges.values()])
    ax.set_title(
        f"Entanglement bottleneck map — IQM Emerald  "
        f"(81 edges · 1000 shots/circuit)\n"
        f"mean F = {Fs.mean():.3f}    "
        f"min F = {Fs.min():.3f} on edge {worst_e[0]}    "
        f"weakest qubit Q{worst_q[0]} (Q = {worst_q[1]:.3f})",
        fontsize=12,
    )
    ax.set_aspect("equal")
    ax.axis("off")

    # Legend
    handles = [
        plt.Line2D([0], [0], color="black", lw=2.5, linestyle=(0, (4, 2)),
                   label=f"weakest edge ({worst_e[0]})"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="none", markeredgecolor="black",
                   markeredgewidth=2.5, markersize=18,
                   label=f"weakest qubit (Q{worst_q[0]})"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="gray", markersize=8,
                   label="node = qubit  (size, color = mean F over incident edges)"),
    ]
    ax.legend(handles=handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.05), ncol=3, fontsize=9,
              frameon=False)

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_PNG}")


def plot_edges_sorted(edges):
    rows = sorted(edges.values(), key=lambda d: d["F"])
    Fs = np.array([d["F"] for d in rows])
    sigmas = np.array([d["sigma_F"] for d in rows])
    labels = [f"{d['edge'][0]}-{d['edge'][1]}" for d in rows]
    n = len(rows)

    fig, ax = plt.subplots(figsize=(14, 5))
    cmap = mpl.colormaps["RdYlGn"]
    norm = mpl.colors.Normalize(vmin=0.5, vmax=1.0)
    colors = cmap(norm(np.clip(Fs, 0.5, 1.0)))
    ax.bar(range(n), Fs, yerr=3 * sigmas, color=colors,
           edgecolor="0.2", linewidth=0.4, error_kw=dict(lw=0.6, capsize=1.5))
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1,
               label="separability bound (F = 0.5)")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_ylabel("F  (error bars = 3σ)")
    ax.set_xlabel("edge (sorted by F)")
    ax.set_title(
        "All 81 edges of IQM Emerald, sorted by graph-state fidelity   "
        f"·   mean = {Fs.mean():.3f}  ·  min = {Fs.min():.3f}  ·  max = {Fs.max():.3f}",
        fontsize=11,
    )
    ax.set_ylim(0.4, 1.02)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_BAR_PNG, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_BAR_PNG}")


def main() -> int:
    if not JSON_PATH.exists():
        sys.exit(f"missing {JSON_PATH} — run with shots first")
    edges, nodes, Q = load_results()
    pos = make_layout(nodes, edges)
    plot_chip(edges, nodes, Q, pos)
    plot_edges_sorted(edges)
    return 0


if __name__ == "__main__":
    sys.exit(main())
