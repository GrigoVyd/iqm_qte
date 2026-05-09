"""Render the full Emerald entanglement bottleneck map using the actual chip
coordinates from IQM's published topology image (1:1 visual match).

Differences vs `plot_full_map.py`:
  - positions taken from the IQM Emerald layout (read off the published figure),
    not Kamada-Kawai
  - qubit labels use IQM's 1-indexed QBn naming
  - excluded qubits (no entry in today's coupling map) are still drawn as
    light gray hollow circles, with dashed couplers, for context

The mapping QB1..QB54  →  qiskit index 0..53 is just (n - 1).
Coordinates are in the rotated-grid system: qubit at lattice (i, j) is at
(i, j) screen position; this is already what gives the chip its diamond shape.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

JSON_PATH = Path("results/full_map.json")
OUT_PNG = Path("results/full_map_chiplayout.png")

# (x, y) for QB1 .. QB54 — read off the IQM Resonance Emerald figure.
# Rotated square lattice, neighbours are diagonal (Δx, Δy) ∈ {±(1, ±1)}.
QB_POS: dict[int, tuple[int, int]] = {
    1:  (10, 1),  2:  (11, 2),  3:  (7, 0),   4:  (8, 1),   5:  (9, 2),
    6:  (10, 3),  7:  (11, 4),  8:  (5, 0),   9:  (6, 1),   10: (7, 2),
    11: (8, 3),   12: (9, 4),   13: (10, 5),  14: (11, 6),  15: (3, 0),
    16: (4, 1),   17: (5, 2),   18: (6, 3),   19: (7, 4),   20: (8, 5),
    21: (9, 6),   22: (10, 7),  23: (2, 1),   24: (3, 2),   25: (4, 3),
    26: (5, 4),   27: (6, 5),   28: (7, 6),   29: (8, 7),   30: (9, 8),
    31: (10, 9),  32: (2, 3),   33: (3, 4),   34: (4, 5),   35: (5, 6),
    36: (6, 7),   37: (7, 8),   38: (8, 9),   39: (9, 10),  40: (1, 4),
    41: (2, 5),   42: (3, 6),   43: (4, 7),   44: (5, 8),   45: (6, 9),
    46: (7, 10),  47: (1, 6),   48: (2, 7),   49: (3, 8),   50: (4, 9),
    51: (5, 10),  52: (1, 8),   53: (2, 9),   54: (3, 10),
}


def load_results():
    rows = json.loads(JSON_PATH.read_text())
    edges = {tuple(sorted(d["edge"])): d for d in rows}
    nodes_with_data = sorted({q for e in edges for q in e})
    incident = defaultdict(list)
    for (a, b), d in edges.items():
        incident[a].append(d["F"])
        incident[b].append(d["F"])
    Q = {q: float(np.mean(incident[q])) for q in nodes_with_data}
    return edges, nodes_with_data, Q


def all_chip_edges() -> set[tuple[int, int]]:
    """Every nearest-neighbour pair on the rotated-square lattice that has
    both endpoints in the QB_POS dict — i.e. the *physical* coupling map."""
    pos_to_qb = {pos: qb for qb, pos in QB_POS.items()}
    edges = set()
    for qb, (x, y) in QB_POS.items():
        for dx, dy in [(1, 1), (1, -1)]:  # half the directions = each edge once
            nbr = pos_to_qb.get((x + dx, y + dy))
            if nbr is not None:
                a, b = qb - 1, nbr - 1                # qiskit indices
                edges.add((min(a, b), max(a, b)))
    return edges


def main() -> int:
    if not JSON_PATH.exists():
        sys.exit(f"missing {JSON_PATH}")
    measured, nodes_with_data, Q = load_results()
    chip_edges = all_chip_edges()

    fig, ax = plt.subplots(figsize=(13, 11))

    f_min, f_max = 0.5, 1.0
    cmap = mpl.colormaps["RdYlGn"]
    norm = mpl.colors.Normalize(vmin=f_min, vmax=f_max)

    # 1. Draw every physical edge — measured ones colored by F, missing ones
    #    dashed gray (those couplers are not in today's coupling map).
    for e in chip_edges:
        a, b = e
        qb_a, qb_b = a + 1, b + 1
        x0, y0 = QB_POS[qb_a]
        x1, y1 = QB_POS[qb_b]
        if e in measured:
            F = measured[e]["F"]
            col = cmap(norm(max(f_min, F)))
            w = 1.5 + 5.0 * max(0.0, (F - 0.5) / 0.5)
            ax.plot([x0, x1], [y0, y1], color=col, linewidth=w,
                    solid_capstyle="round", zorder=2)
        else:
            ax.plot([x0, x1], [y0, y1], color="0.7", linewidth=1.5,
                    linestyle=(0, (3, 3)), zorder=1)

    # 2. Identify worst spot
    worst_q = min(Q.items(), key=lambda kv: kv[1])
    worst_e = min(measured.items(), key=lambda kv: kv[1]["F"])

    # 3. Draw qubit nodes
    for qb, (x, y) in QB_POS.items():
        idx = qb - 1
        if idx in Q:
            q_val = Q[idx]
            face = cmap(norm(q_val))
            edge = "black"
            text_col = "black" if q_val > 0.78 else "white"
            ax.scatter([x], [y], s=820, c=[face], edgecolors=edge,
                       linewidths=1.0, zorder=3)
            ax.text(x, y, f"QB{qb}", ha="center", va="center",
                    fontsize=8.5, color=text_col, zorder=4,
                    fontweight="bold" if idx == worst_q[0] else "normal")
        else:
            ax.scatter([x], [y], s=820, facecolors="white",
                       edgecolors="0.6", linewidths=1.0,
                       linestyle="--", zorder=3)
            ax.text(x, y, f"QB{qb}", ha="center", va="center",
                    fontsize=8.5, color="0.5", zorder=4)

    # 4. Highlights
    a, b = worst_e[0]
    x0, y0 = QB_POS[a + 1]; x1, y1 = QB_POS[b + 1]
    ax.plot([x0, x1], [y0, y1], color="black", linewidth=8,
            linestyle=(0, (4, 2)), alpha=0.55, zorder=1.5)
    wx, wy = QB_POS[worst_q[0] + 1]
    ax.scatter([wx], [wy], s=2200, facecolors="none",
                edgecolors="black", linewidths=2.5, zorder=5)

    # 5. Colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("graph-state fidelity F  (separable bound = 0.5)",
                    fontsize=10)

    # 6. Title + legend
    Fs = np.array([d["F"] for d in measured.values()])
    n_meas = len(measured)
    n_total = len(chip_edges)
    ax.set_title(
        f"Entanglement bottleneck map — IQM Emerald  (chip layout)\n"
        f"{n_meas}/{n_total} edges measured · 1000 shots/circuit · "
        f"mean F = {Fs.mean():.3f}, min F = {Fs.min():.3f} on "
        f"QB{worst_e[0][0]+1}-QB{worst_e[0][1]+1}, "
        f"weakest qubit QB{worst_q[0]+1} (Q={worst_q[1]:.3f})",
        fontsize=12,
    )
    ax.set_aspect("equal")
    ax.axis("off")

    handles = [
        plt.Line2D([0], [0], color="black", lw=2.5, linestyle=(0, (4, 2)),
                   label=f"weakest edge (QB{worst_e[0][0]+1}-QB{worst_e[0][1]+1})"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="none", markeredgecolor="black",
                   markeredgewidth=2.5, markersize=18,
                   label=f"weakest qubit (QB{worst_q[0]+1})"),
        plt.Line2D([0], [0], color="0.7", lw=2, linestyle=(0, (3, 3)),
                   label="excluded coupler (not in today's calibration)"),
    ]
    ax.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.04), ncol=3, fontsize=9, frameon=False)

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    print(f"Saved: {OUT_PNG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
