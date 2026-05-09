"""Readout error mitigation utilities.

Two strategies:
1. Full QREM: calibrate a 2^n × 2^n transition matrix (practical for n ≤ 10)
2. Parity-QREM: calibrate only 2-qubit parity subspaces (scales to 50+ qubits)

For the GME witness we only need parity products of stabilizers,
so parity-QREM is the right choice for large systems.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit


# ---------------------------------------------------------------------------
# Full QREM (small systems, n ≤ 10)
# ---------------------------------------------------------------------------

def build_calibration_circuits(n: int) -> list[QuantumCircuit]:
    """Build 2^n calibration circuits (one per computational basis state)."""
    circuits = []
    for state in range(2 ** n):
        qc = QuantumCircuit(n, name=f"cal_{state:0{n}b}")
        for q in range(n):
            if (state >> q) & 1:
                qc.x(q)
        qc.measure_all()
        circuits.append(qc)
    return circuits


def compute_calibration_matrix(cal_results: list[dict[str, int]], n: int) -> np.ndarray:
    """Compute the n-qubit transition matrix A where A[i,j] = P(measure i | prepared j).

    cal_results[j] = counts when state j was prepared.
    """
    dim = 2 ** n
    A = np.zeros((dim, dim))
    for j, counts in enumerate(cal_results):
        total = sum(counts.values())
        for bitstring, cnt in counts.items():
            bits = bitstring.replace(" ", "")
            i = int(bits, 2)
            A[i, j] = cnt / total
    return A


def apply_qrem(counts: dict[str, int], A: np.ndarray) -> dict[str, int]:
    """Apply full QREM: invert calibration matrix and correct counts vector."""
    n = int(np.log2(A.shape[0]))
    dim = 2 ** n
    total = sum(counts.values())

    prob_vec = np.zeros(dim)
    for bitstring, cnt in counts.items():
        bits = bitstring.replace(" ", "")
        i = int(bits, 2)
        prob_vec[i] = cnt / total

    A_inv = np.linalg.pinv(A)
    corrected = A_inv @ prob_vec
    corrected = np.clip(corrected, 0, None)
    corrected /= corrected.sum()

    return {f"{i:0{n}b}": int(round(corrected[i] * total)) for i in range(dim)
            if corrected[i] > 0}


# ---------------------------------------------------------------------------
# Parity-QREM (large systems)
# ---------------------------------------------------------------------------

def build_parity_calibration_circuits(qubits: list[int], n_total: int) -> list[QuantumCircuit]:
    """Build 4 circuits per qubit pair for parity calibration.

    For each pair (i, j), prepare |00⟩, |01⟩, |10⟩, |11⟩ and measure.
    Returns 4 * C(len(qubits), 2) circuits.
    """
    circuits = []
    pairs = [(qubits[i], qubits[j]) for i in range(len(qubits)) for j in range(i+1, len(qubits))]
    for qi, qj in pairs:
        for state in range(4):
            qc = QuantumCircuit(n_total, name=f"pcal_{qi}_{qj}_{state:02b}")
            if state & 1:
                qc.x(qi)
            if state & 2:
                qc.x(qj)
            qc.measure_all()
            circuits.append((qi, qj, state, qc))
    return circuits


def parity_corrected_expval(counts: dict[str, int], qubit_i: int, qubit_j: int,
                             n_total: int, cal_matrix_2q: np.ndarray) -> float:
    """Compute readout-mitigated ZZ parity ⟨ZiZj⟩ from counts.

    cal_matrix_2q: 4×4 transition matrix for the (i,j) qubit pair.
    """
    total = sum(counts.values())
    prob = np.zeros(4)
    for bitstring, cnt in counts.items():
        bits = bitstring.replace(" ", "")
        bi = int(bits[n_total - 1 - qubit_i])
        bj = int(bits[n_total - 1 - qubit_j])
        idx = bi + 2 * bj
        prob[idx] += cnt / total

    A_inv = np.linalg.pinv(cal_matrix_2q)
    corrected = np.clip(A_inv @ prob, 0, None)
    corrected /= corrected.sum() if corrected.sum() > 0 else 1

    # ZZ eigenvalue: +1 if same parity, -1 if different
    zz = corrected[0] - corrected[1] - corrected[2] + corrected[3]
    return float(zz)


def run_full_qrem(backend, circuits: list[QuantumCircuit], n: int,
                  shots: int = 1000) -> tuple[list[dict[str,int]], np.ndarray]:
    """Run calibration + mitigate a list of n-qubit circuits.

    Returns (mitigated_counts_list, calibration_matrix).
    Only practical for n ≤ 10.
    """
    from qiskit import transpile as qk_transpile

    cal_circuits = build_calibration_circuits(n)
    all_circuits = cal_circuits + circuits
    transpiled = qk_transpile(all_circuits, backend=backend, optimization_level=1)
    job = backend.run(transpiled, shots=shots)
    raw = job.result().get_counts()
    if isinstance(raw, dict):
        raw = [raw]

    cal_results = raw[:2**n]
    exp_results = raw[2**n:]

    A = compute_calibration_matrix(cal_results, n)
    mitigated = [apply_qrem(c, A) for c in exp_results]
    return mitigated, A
