"""Coupling map helpers and topology visualisation."""

from __future__ import annotations

import matplotlib.pyplot as plt
import rustworkx as rx
from rustworkx.visualization import mpl_draw


def coupling_map_to_graph(coupling_map: list[list[int]]) -> rx.PyGraph:
    """Convert a Qiskit coupling map (list of [a,b] edges) to a rustworkx graph."""
    g = rx.PyGraph()
    nodes: set[int] = set()
    for a, b in coupling_map:
        nodes.add(a)
        nodes.add(b)
    max_node = max(nodes) if nodes else 0
    # Add nodes with their index as label
    for _ in range(max_node + 1):
        g.add_node(None)
    for a, b in coupling_map:
        if not g.has_edge(a, b):
            g.add_edge(a, b, None)
    return g


def plot_topology(backend, highlight_qubits: list[int] | None = None,
                  title: str = "Qubit topology") -> None:
    """Plot the device coupling map, optionally highlighting selected qubits."""
    coupling = list(backend.coupling_map)
    n = backend.num_qubits

    g = rx.PyGraph()
    for _ in range(n):
        g.add_node(None)
    for a, b in coupling:
        if not g.has_edge(a, b):
            g.add_edge(a, b, None)

    colors = []
    for i in range(n):
        if highlight_qubits and i in highlight_qubits:
            colors.append("#e84545")
        else:
            colors.append("#32a8a4")

    fig, ax = plt.subplots(figsize=(10, 8))
    mpl_draw(g, ax=ax, with_labels=True, node_color=colors,
             font_color="white", node_size=400)
    ax.set_title(title)
    plt.tight_layout()
    plt.show()


def get_reduced_coupling_map(full_coupling_map: list[list[int]],
                             qubit_indices: list[int]) -> list[list[int]]:
    """Return only edges where both endpoints are in qubit_indices."""
    qubit_set = set(qubit_indices)
    return [list(e) for e in full_coupling_map if set(e).issubset(qubit_set)]


def relabel_to_contiguous(qubit_indices: list[int],
                          coupling_map: list[list[int]]) -> tuple[list[list[int]], dict[int, int]]:
    """Relabel qubit_indices to 0..n-1 and update coupling_map accordingly.

    Returns (new_coupling_map, old_to_new mapping).
    """
    old_to_new = {old: new for new, old in enumerate(qubit_indices)}
    new_cm = [[old_to_new[a], old_to_new[b]] for a, b in coupling_map]
    return new_cm, old_to_new


def coloring_2d_grid(rows: int, cols: int) -> list[int]:
    """Return checkerboard 2-coloring (0/1) for a rows×cols grid.

    Index mapping: qubit (r, c) → r*cols + c.
    Black qubits (color=0): (r+c) even. White (color=1): (r+c) odd.
    """
    return [(r + c) % 2 for r in range(rows) for c in range(cols)]
