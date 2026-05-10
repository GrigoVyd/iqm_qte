"""F-gate (Diker 2016) cascade preparation of |W_n>.

|W_n> = (1/sqrt(n)) sum_{k=0..n-1} |0...1_k...0>

The F-gate F_k = R_y(-theta_k) . CZ . R_y(theta_k) with
theta_k = arccos(sqrt(1/(n - k + 1))) implements an amplitude-splitter
along a chain of nearest-neighbour qubits. A final CNOT ladder undoes the
phase pattern. Total: ~3n CZ-equivalent gates, depth O(n), no SWAPs when
mapped onto a Hamiltonian path.
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit


def _f_gate(qc: QuantumCircuit, ctrl: int, tgt: int, theta: float) -> None:
    qc.ry(-theta, tgt)
    qc.cz(ctrl, tgt)
    qc.ry(theta, tgt)


def build_w_state(n: int) -> QuantumCircuit:
    """Return the n-qubit W-state preparation circuit on a linear chain."""
    qc = QuantumCircuit(n, name=f"W{n}")
    qc.x(0)
    for k in range(1, n):
        theta = float(np.arccos(np.sqrt(1.0 / (n - k + 1))))
        _f_gate(qc, k - 1, k, theta)
    for k in range(n - 1):
        qc.cx(k + 1, k)
    return qc


def w_state_z_circuit(n: int) -> QuantumCircuit:
    qc = build_w_state(n)
    qc.measure_all()
    return qc


def w_state_x_circuit(n: int) -> QuantumCircuit:
    qc = build_w_state(n)
    qc.h(range(n))
    qc.measure_all()
    return qc
