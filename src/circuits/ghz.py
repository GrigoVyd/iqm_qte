"""Topology-aware GHZ state circuit builder."""

from __future__ import annotations

from qiskit import QuantumCircuit

from src.backend import QubitLayout


def build_ghz(layout: QubitLayout) -> QuantumCircuit:
    """Build an n-qubit GHZ circuit respecting the given qubit chain.

    The layout must be a linear chain (from get_best_chain). The circuit
    applies H to the first qubit, then CNOT along the chain — no SWAPs needed
    as long as consecutive qubits in layout.qubit_indices are coupled.

    Returns a QuantumCircuit on n qubits (local indices 0..n-1).
    The caller must transpile with the appropriate coupling_map.
    """
    n = len(layout.qubit_indices)
    qc = QuantumCircuit(n, name=f"GHZ-{n}")
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i) if _star_topology(layout) else qc.cx(i - 1, i)
    qc.measure_all()
    return qc


def build_ghz_chain(n: int) -> QuantumCircuit:
    """Build a GHZ circuit as a linear CNOT chain (qubits 0..n-1, no layout needed).

    Use when you have already selected a linear chain of coupled qubits.
    """
    qc = QuantumCircuit(n, name=f"GHZ-chain-{n}")
    qc.h(0)
    for i in range(1, n):
        qc.cx(i - 1, i)
    qc.measure_all()
    return qc


def build_ghz_no_measure(n: int) -> QuantumCircuit:
    """GHZ state without measurement — used as a sub-circuit for witnesses."""
    qc = QuantumCircuit(n, name=f"GHZ-{n}-state")
    qc.h(0)
    for i in range(1, n):
        qc.cx(i - 1, i)
    return qc


def _star_topology(layout: QubitLayout) -> bool:
    """True if the coupling map forms a star (one central qubit connected to all others)."""
    n = len(layout.qubit_indices)
    if n <= 2:
        return False
    degree: dict[int, int] = {}
    for a, b in layout.coupling_map:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    return max(degree.values()) == n - 1
