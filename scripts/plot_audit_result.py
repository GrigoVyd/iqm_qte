"""Render the 2x3 audit result as a labelled subgraph plot.

Loads results/audit_gme_2x3.json (no hardware needed) and produces:
  - results/audit_gme_2x3_subgraph.png  — focused 6-qubit layout, edges
                                          colored + labelled with F (3σ flag)
  - per-edge legend with the entanglement-classification color code
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import rustworkx as rx
from rustworkx.visualization import mpl_draw

JSON_PATH = Path("results/audit_gme_2x3.json")
OUT_PNG = Path("results/audit_gme_2x3_subgraph.png")


def color_for(F: float, ent_3sigma: bool) -> str:
    if ent_3sigma:
        return "#2ca02c"
    if F > 0.5:
        return "#ffbf00"
    return "#d62728"


def main() -> int:
    if not JSON_PATH.exists():
        sys.exit(f"missing {JSON_PATH} — run audit_gme_2x3_with_crosstalk.py first")
    data = json.loads(JSON_PATH.read_text())

    layout = data["layout"]
    edges_par = data["edges_parallel"]
    audit = data["audit_parallel"]
    rows, cols = data["rows"], data["cols"]

    # Position each PHYSICAL qubit at its LOGICAL grid coordinate so the
    # plot reflects how the cluster circuit actually uses the chip.
    pos_phys: dict[int, tuple[float, float]] = {}
    for log_idx, phys in enumerate(layout):
        r, c = log_idx // cols, log_idx % cols
        pos_phys[phys] = (c, -r)

    nodes = list(pos_phys.keys())
    node_idx = {q: k for k, q in enumerate(nodes)}
    g = rx.PyGraph()
    g.add_nodes_from([str(q) for q in nodes])

    edge_colors = []
    for e in edges_par:
        a, b = e["edge"]
        g.add_edge(node_idx[a], node_idx[b], None)
        edge_colors.append(color_for(e["F"], e["entangled_3sigma"]))

    fig, ax = plt.subplots(figsize=(9, 6))
    mpl_draw(
        g, ax=ax,
        pos={node_idx[q]: pos_phys[q] for q in nodes},
        with_labels=True, labels=lambda s: f"Q{s}",
        node_color="#dddddd", node_size=1100, font_size=11,
        edge_color=edge_colors, width=4.0,
    )

    # Manual edge labels at midpoints (rustworkx edge_labels rendering is poor)
    for e in edges_par:
        a, b = e["edge"]
        x = (pos_phys[a][0] + pos_phys[b][0]) / 2
        y = (pos_phys[a][1] + pos_phys[b][1]) / 2
        ax.annotate(
            f"{e['F']:.3f}",
            xy=(x, y), xytext=(0, 6), textcoords="offset points",
            ha="center", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.25", fc="white",
                       ec="0.7", alpha=0.9),
        )

    title = (
        f"Selector audit on {data['device']} — layout {layout}\n"
        f"validated={audit['validated']}  "
        f"mean F={audit['mean_F']:.3f}  "
        f"min F={audit['min_F']:.3f}  "
        f"({audit['n_edges'] - audit['n_below_3sigma']}/{audit['n_edges']} "
        f"edges 3σ-entangling)"
    )
    ax.set_title(title, fontsize=11)

    # Legend
    legend_elems = [
        plt.Line2D([0], [0], color="#2ca02c", lw=4,
                    label="F − 3σ > 0.5  (certified entangling)"),
        plt.Line2D([0], [0], color="#ffbf00", lw=4,
                    label="F > 0.5 but not 3σ-significant"),
        plt.Line2D([0], [0], color="#d62728", lw=4,
                    label="F ≤ 0.5  (not entangled)"),
    ]
    ax.legend(handles=legend_elems, loc="lower center",
               bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=9, frameon=False)

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_PNG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
