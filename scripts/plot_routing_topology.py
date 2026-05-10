"""Topology plot for the cluster-routing sweep on Garnet.

Loads results/cluster_routing_garnet_sweep.json (the routes + selector layout)
and visualises:

  - all Garnet qubits + couplers (full chip topology, Kamada-Kawai layout)
  - the 2x5 logical sub-grid the selector picked, drawn as light blue edges
  - each route's path overlaid in colour:
        green if both Bell + CHSH passed at 3σ
        amber if Bell passed, CHSH failed
        red   if Bell witness failed
  - endpoints A and B labelled per route
  - F and |S| values per route in the legend

The first run reads the coupling map directly from the backend (no shots, free)
and caches it to results/garnet_topology.json. Subsequent runs skip the backend
call.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

JSON_RESULTS = Path("results/cluster_routing_garnet_sweep.json")
TOPO_CACHE = Path("results/garnet_topology.json")
OUT_PNG = Path("results/cluster_routing_topology.png")


def load_or_fetch_topology(device: str) -> tuple[list[tuple[int, int]], int]:
    if TOPO_CACHE.exists():
        d = json.loads(TOPO_CACHE.read_text())
        if d.get("device") == device:
            return [tuple(e) for e in d["edges"]], d["num_qubits"]
    # Fetch from backend
    from src.backend import get_backend, is_simulator
    backend = get_backend(token=os.environ.get("IQM_TOKEN"), device=device)
    if is_simulator(backend):
        sys.exit(f"need IQM token to fetch {device} topology")
    edges = sorted({(min(a, b), max(a, b)) for a, b in backend.coupling_map})
    n = backend.num_qubits
    TOPO_CACHE.write_text(json.dumps({
        "device": device, "num_qubits": n,
        "edges": [list(e) for e in edges],
    }, indent=2))
    print(f"  cached topology to {TOPO_CACHE}")
    return edges, n


def status_color(res: dict) -> str:
    bell_ok = res["bell"]["entangled_3sigma"]
    chsh_ok = res.get("chsh", {}).get("pass_3sigma", False)
    if bell_ok and chsh_ok:
        return "#2ca02c"
    if bell_ok and not chsh_ok:
        return "#ffbf00"
    return "#d62728"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--device", default="garnet")
    args = p.parse_args()

    if not JSON_RESULTS.exists():
        sys.exit(f"missing {JSON_RESULTS}")
    routes = json.loads(JSON_RESULTS.read_text())
    routes = [r for r in routes if not r.get("no_cz", False)]
    routes.sort(key=lambda r: r["L"])

    edges, n_qubits = load_or_fetch_topology(args.device)

    # We need the selector layout. It is implicit in the result: each route's
    # "path" is in *logical* indices of the 2x5 sub-grid. The mapping
    # logical→physical was the same for every route in this sweep (one
    # `select_best_subgrid` call). We don't store it directly in the JSON, so
    # we re-derive it by re-running the selector (free, no shots).
    from src.backend import get_backend, select_best_subgrid
    backend = get_backend(token=os.environ.get("IQM_TOKEN"),
                            device=args.device)
    rows, cols = routes[0]["rows"], routes[0]["cols"]
    layout, _, _ = select_best_subgrid(backend, rows, cols)
    print(f"  layout: {layout}")

    # Build the chip graph
    g = nx.Graph()
    g.add_nodes_from(range(n_qubits))
    g.add_edges_from(edges)
    pos = nx.kamada_kawai_layout(g)

    fig, ax = plt.subplots(figsize=(15, 11))

    # Draw all chip edges in light gray
    for (u, v) in edges:
        x0, y0 = pos[u]; x1, y1 = pos[v]
        ax.plot([x0, x1], [y0, y1], color="#cccccc",
                linewidth=1.5, zorder=1)

    # Draw the selector's 2x5 logical sub-grid edges in light blue
    sub_edges_logical: set[tuple[int, int]] = set()
    for r in range(rows):
        for c in range(cols - 1):
            sub_edges_logical.add((r * cols + c, r * cols + c + 1))
    for r in range(rows - 1):
        for c in range(cols):
            sub_edges_logical.add((r * cols + c, (r + 1) * cols + c))

    for (a, b) in sub_edges_logical:
        u, v = layout[a], layout[b]
        x0, y0 = pos[u]; x1, y1 = pos[v]
        ax.plot([x0, x1], [y0, y1], color="#a8c8e8",
                linewidth=3, zorder=2)

    # Per logical edge: which is the shortest route that introduced it?
    # That is the route whose status determines this edge's color, telling the
    # story "passed up to here, broke at this hop".
    edge_origin: dict[tuple[int, int], dict] = {}
    for res in sorted(routes, key=lambda r: r["L"]):
        for a, b in zip(res["path"], res["path"][1:]):
            key = (min(a, b), max(a, b))
            edge_origin.setdefault(key, res)

    for (a, b), origin in edge_origin.items():
        u, v = layout[a], layout[b]
        x0, y0 = pos[u]; x1, y1 = pos[v]
        ax.plot([x0, x1], [y0, y1], color=status_color(origin),
                linewidth=5.5, zorder=3, alpha=0.9,
                solid_capstyle="round")

    # Draw all qubits
    xs = [pos[q][0] for q in range(n_qubits)]
    ys = [pos[q][1] for q in range(n_qubits)]
    in_layout = set(layout)
    node_colors = ["#dddddd" if q not in in_layout else "white"
                   for q in range(n_qubits)]
    edge_colors = ["#888888" if q not in in_layout else "black"
                   for q in range(n_qubits)]
    ax.scatter(xs, ys, s=600, c=node_colors,
                edgecolors=edge_colors, linewidths=1.4, zorder=4)
    for q in range(n_qubits):
        ax.text(*pos[q], f"QB{q+1}", ha="center", va="center",
                fontsize=8, zorder=5,
                color="#444" if q not in in_layout else "black",
                fontweight="bold" if q in in_layout else "normal")

    # Mark endpoints A, B of the longest route
    longest = max(routes, key=lambda r: r["L"])
    A_phys, B_phys = layout[longest["path"][0]], layout[longest["path"][-1]]
    for label, q in [("A", A_phys), ("B", B_phys)]:
        ax.scatter(*pos[q], s=1500, facecolors="none",
                    edgecolors="black", linewidths=2.5, zorder=6)
        ax.annotate(label, pos[q], xytext=(15, 15),
                     textcoords="offset points",
                     fontsize=14, fontweight="bold")

    # Title + legend
    ax.set_title(
        "Cluster routing on IQM Garnet — 2×5 logical sub-grid\n"
        f"selector layout (logical 0..9) → physical {layout}",
        fontsize=12,
    )
    ax.set_aspect("equal")
    ax.axis("off")

    handles = [
        plt.Line2D([0], [0], color="#cccccc", lw=2,
                   label="device coupler (not used)"),
        plt.Line2D([0], [0], color="#a8c8e8", lw=3,
                   label="2×5 sub-grid edges (cluster CZs run here)"),
    ]
    for res in routes:
        L = res["L"]
        F = res["bell"]["F"]
        S = res.get("chsh", {}).get("abs_S", float("nan"))
        handles.append(plt.Line2D(
            [0], [0], color=status_color(res), lw=4,
            label=f"L={L:>2}  route   F={F:.3f}, |S|={S:.2f}",
        ))
    ax.legend(handles=handles, loc="upper left",
               bbox_to_anchor=(1.0, 1.0), fontsize=10, frameon=False)

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_PNG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
