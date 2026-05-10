"""Genuine Multipartite Entanglement (GME) witness for 2D cluster states.

Based on:
  - G. Tóth, O. Gühne, "Detecting Genuine Multipartite Entanglement with Two
    Local Measurements", Phys. Rev. Lett. 94, 060501 (2005).
    arXiv:quant-ph/0405165
  - G. Tóth, O. Gühne, "Entanglement detection in the stabilizer formalism",
    Phys. Rev. A 72, 022340 (2005). arXiv:quant-ph/0501020
  - Y. Zhou, Q. Zhao, X. Yuan, X. Ma, "Detecting multipartite entanglement
    structure with minimal resources", npj Quantum Information 5, 83 (2019).
    arXiv:1904.05001

For a 2-colorable graph state (e.g. 2D cluster on square lattice), GME can be
certified using only 2 measurement settings:

  Setting A: for each BLACK qubit i, measure X_i ⊗ ⊗_{j∈N(i)} Z_j
  Setting B: for each WHITE qubit i, measure X_i ⊗ ⊗_{j∈N(i)} Z_j

The witness value is:
  W = Σ_i ⟨g_i⟩

where g_i is the stabilizer of qubit i: g_i = X_i ⊗ ⊗_{j∈N(i)} Z_j.

For the ideal cluster state: ⟨g_i⟩ = +1 for all i, so W_ideal = n.
For any fully biseparable state: W ≤ n - 1 (classical bound).

W > n - 1  ⟹  state is GENUINELY MULTIPARTITELY ENTANGLED.

The checkerboard coloring ensures:
  - Setting A measures all stabilizers of black qubits (neighbours are white → Z)
  - Setting B measures all stabilizers of white qubits (neighbours are black → Z)
Each setting requires measuring X on one colour class and Z on the other —
these commute within each setting, so both can be read from a single shot.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit


def build_gme_circuits(
    state_circuit: QuantumCircuit,
    rows: int,
    cols: int,
) -> tuple[QuantumCircuit, QuantumCircuit]:
    """Return two measurement circuits (setting A, setting B) for the GME witness.

    state_circuit: cluster state preparation (without measurement).
    rows, cols: grid dimensions. Qubit (r,c) → r*cols+c.

    Setting A: black qubits (r+c even) → X basis; white qubits → Z basis.
    Setting B: white qubits (r+c odd)  → X basis; black qubits → Z basis.

    Returns (circuit_A, circuit_B).
    """
    n = rows * cols

    def is_black(r, c):
        return (r + c) % 2 == 0

    def make_circuit(measure_black_in_x: bool) -> QuantumCircuit:
        qc = state_circuit.copy()
        for r in range(rows):
            for c in range(cols):
                q = r * cols + c
                if is_black(r, c) == measure_black_in_x:
                    qc.h(q)   # X basis: H before measure
                # else: Z basis, no rotation needed
        qc.measure_all()
        return qc

    circuit_a = make_circuit(measure_black_in_x=True)   # black→X, white→Z
    circuit_b = make_circuit(measure_black_in_x=False)  # white→X, black→Z
    return circuit_a, circuit_b


def compute_gme_witness(
    counts_a: dict[str, int],
    counts_b: dict[str, int],
    rows: int,
    cols: int,
) -> dict:
    """Compute the GME witness value W = Σ_i ⟨g_i⟩.

    Each stabilizer g_i = X_i ⊗ Z_{neighbours} contributes ±1.

    For setting A (black qubits measured in X):
      ⟨g_i⟩ for black qubit i = average of (x_i * Π_{j∈N(i)} z_j)
      where x_i = ±1 from qubit i's X measurement, z_j = ±1 from Z measurement.

    For setting B (white qubits measured in X):
      ⟨g_i⟩ for white qubit i similarly.
    """
    n = rows * cols

    def bit_to_pm1(bit: str) -> int:
        return 1 if bit == "0" else -1

    def neighbours(r: int, c: int) -> list[int]:
        nbrs = []
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                nbrs.append(nr * cols + nc)
        return nbrs

    def is_black(r, c):
        return (r + c) % 2 == 0

    def stabilizer_expval(counts: dict[str, int], qubit: int,
                           nbr_qubits: list[int]) -> float:
        total = sum(counts.values())
        exp = 0.0
        for bitstring, cnt in counts.items():
            bits = bitstring.replace(" ", "")
            # Qiskit bitstrings: bit index i = position (n-1-i) from left
            def get_bit(q):
                return bits[n - 1 - q]

            val = bit_to_pm1(get_bit(qubit))
            for nb in nbr_qubits:
                val *= bit_to_pm1(get_bit(nb))
            exp += val * cnt / total
        return exp

    stabilizer_values = {}

    for r in range(rows):
        for c in range(cols):
            q = r * cols + c
            nbrs = neighbours(r, c)
            if is_black(r, c):
                # stabilizer g_q measured in setting A
                ev = stabilizer_expval(counts_a, q, nbrs)
            else:
                # stabilizer g_q measured in setting B
                ev = stabilizer_expval(counts_b, q, nbrs)
            stabilizer_values[q] = ev

    W = sum(stabilizer_values.values())

    # Classical (biseparable) bound = n - 1
    bisep_bound = n - 1
    violation = W - bisep_bound
    ideal = float(n)

    return {
        "W": W,
        "W_ideal": ideal,
        "biseparable_bound": bisep_bound,
        "violation": violation,
        "violation_fraction": violation / ideal,
        "is_gme": W > bisep_bound,
        "n_qubits": n,
        "stabilizer_values": stabilizer_values,
    }


def gme_significance(W: float, n: int, shots: int) -> float:
    """Estimate z-score of GME violation.

    Each stabilizer expectation value has std ≈ 1/√shots.
    W = sum of n stabilizers, so std(W) ≈ √n / √shots.
    Violation = W - (n-1).
    """
    std_W = np.sqrt(n) / np.sqrt(shots)
    bisep_bound = n - 1
    return (W - bisep_bound) / std_W if std_W > 0 else 0.0


def run_gme(backend, state_circuit: QuantumCircuit, rows: int, cols: int,
            shots: int = 10000) -> dict:
    """Run GME witness on backend and return full results dict."""
    from qiskit import transpile as qk_transpile

    circ_a, circ_b = build_gme_circuits(state_circuit, rows, cols)
    transpiled = qk_transpile([circ_a, circ_b], backend=backend, optimization_level=3)

    job = backend.run(transpiled, shots=shots)
    raw = job.result().get_counts()
    if isinstance(raw, dict):
        counts_a = counts_b = raw
    else:
        counts_a, counts_b = raw[0], raw[1]

    result = compute_gme_witness(counts_a, counts_b, rows, cols)
    result["significance_sigma"] = gme_significance(result["W"], rows * cols, shots)
    return result
