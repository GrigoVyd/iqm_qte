"""2D cluster state circuit builder for square lattice topology.

Circuit depth = 3 layers regardless of qubit count:
  Layer 1: H on all qubits
  Layer 2: CZ on all horizontal (row) edges
  Layer 3: CZ on all vertical (column) edges

CZ is a native gate on IQM hardware — no decomposition needed.
"""

from __future__ import annotations

from qiskit import QuantumCircuit

from src.backend import QubitLayout


def build_cluster_2d(layout: QubitLayout, rows: int, cols: int) -> QuantumCircuit:
    """Build a 2D cluster state on the given grid layout.

    layout must come from get_best_grid(backend, rows, cols).
    Qubit mapping: layout.qubit_indices[r*cols + c] is grid position (r,c).

    Returns a QuantumCircuit on rows*cols qubits (local indices 0..n-1).
    """
    n = rows * cols
    qc = QuantumCircuit(n, name=f"Cluster-{rows}x{cols}")

    # Layer 1: Hadamard all
    qc.h(range(n))

    # Layer 2: CZ along horizontal edges (same row, adjacent columns)
    qc.barrier()
    for r in range(rows):
        for c in range(cols - 1):
            q_left = r * cols + c
            q_right = r * cols + c + 1
            qc.cz(q_left, q_right)

    # Layer 3: CZ along vertical edges (adjacent rows, same column)
    qc.barrier()
    for r in range(rows - 1):
        for c in range(cols):
            q_top = r * cols + c
            q_bot = (r + 1) * cols + c
            qc.cz(q_top, q_bot)

    return qc


def build_cluster_2d_no_measure(rows: int, cols: int) -> QuantumCircuit:
    """Cluster state without measurement — for use as a sub-circuit."""
    n = rows * cols
    qc = QuantumCircuit(n, name=f"Cluster-{rows}x{cols}-state")
    qc.h(range(n))
    qc.barrier()
    for r in range(rows):
        for c in range(cols - 1):
            qc.cz(r * cols + c, r * cols + c + 1)
    qc.barrier()
    for r in range(rows - 1):
        for c in range(cols):
            qc.cz(r * cols + c, (r + 1) * cols + c)
    return qc


def get_grid_edges(rows: int, cols: int) -> tuple[list[tuple[int,int]], list[tuple[int,int]]]:
    """Return (horizontal_edges, vertical_edges) for a rows×cols grid.

    Qubit index convention: (r, c) → r*cols + c.
    """
    h_edges = [(r*cols+c, r*cols+c+1) for r in range(rows) for c in range(cols-1)]
    v_edges = [(r*cols+c, (r+1)*cols+c) for r in range(rows-1) for c in range(cols)]
    return h_edges, v_edges


def checkerboard_coloring(rows: int, cols: int) -> list[int]:
    """Checkerboard 2-coloring: 0 (black) if (r+c) even, 1 (white) if odd.

    Used by the GME witness to split qubits into two measurement groups.
    """
    return [(r + c) % 2 for r in range(rows) for c in range(cols)]
