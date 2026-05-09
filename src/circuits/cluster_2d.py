"""2D cluster state circuit on a square lattice.

Depth = 3 layers regardless of qubit count:
  Layer 1: H on all qubits
  Layer 2: CZ on all horizontal (row) edges
  Layer 3: CZ on all vertical (column) edges

H and CZ are both native on IQM hardware — no decomposition.
"""

from __future__ import annotations

from qiskit import QuantumCircuit


def build_cluster_2d_no_measure(rows: int, cols: int) -> QuantumCircuit:
    """Build a depth-3 cluster state on a rows×cols logical grid.

    Logical qubit (r, c) → index r*cols + c (row-major).
    """
    n = rows * cols
    qc = QuantumCircuit(n, name=f"Cluster-{rows}x{cols}")
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


def checkerboard_coloring(rows: int, cols: int) -> list[int]:
    """Checkerboard 2-coloring: 0 (black) if (r+c) even, 1 (white) if odd."""
    return [(r + c) % 2 for r in range(rows) for c in range(cols)]
