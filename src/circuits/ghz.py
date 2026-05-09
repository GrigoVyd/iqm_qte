"""GHZ state circuit (linear CNOT chain). Used by the Mermin variety layer."""

from __future__ import annotations

from qiskit import QuantumCircuit


def build_ghz_chain(n: int) -> QuantumCircuit:
    """GHZ circuit on qubits 0..n-1 as a linear CNOT chain, with measurement."""
    qc = QuantumCircuit(n, name=f"GHZ-chain-{n}")
    qc.h(0)
    for i in range(1, n):
        qc.cx(i - 1, i)
    qc.measure_all()
    return qc


def build_ghz_no_measure(n: int) -> QuantumCircuit:
    """GHZ state preparation only — used as a sub-circuit for witnesses."""
    qc = QuantumCircuit(n, name=f"GHZ-{n}-state")
    qc.h(0)
    for i in range(1, n):
        qc.cx(i - 1, i)
    return qc
