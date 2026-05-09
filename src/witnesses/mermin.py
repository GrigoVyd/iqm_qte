"""Mermin n-qubit inequality witness.

Uses the operator M_n = Re((X+iY)^⊗n), which expands to:
  M_n = Σ_{k=0,2,4,...} (-1)^{k/2} * Σ_{|S|=k} (⊗_{i∈S} Y_i ⊗_{j∉S} X_j)

i.e. sum of ALL n-fold Pauli products with an EVEN number of Y's,
with alternating signs: (+1) for k=0,4,8,...  and (-1) for k=2,6,10,...

This operator is maximally violated by the standard n-qubit GHZ state:
  |GHZ_n⟩ = (|0⟩^n + |1⟩^n) / √2

Bounds:
  Classical (LHV): |⟨M_n⟩| ≤ 2^(n/2)          [from |Re(z)| ≤ |z| = (√2)^n]
  Quantum GHZ:      ⟨M_n⟩  = 2^(n-1)

  Violation ratio = 2^(n-1) / 2^(n/2) = 2^(n/2 - 1)  (grows exponentially with n)

For n=3: classical ≤ 2√2 ≈ 2.83,  GHZ gives 4
For n=5: classical ≤ 4√2 ≈ 5.66,  GHZ gives 16
For n=7: classical ≤ 8√2 ≈ 11.3,  GHZ gives 64

Number of measurement circuits = 2^(n-1)  (all even-k subsets), which is
feasible for n ≤ 7 (64 circuits) within IQM's 200-circuit batch limit.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from qiskit import QuantumCircuit


# ---------------------------------------------------------------------------
# Operator enumeration
# ---------------------------------------------------------------------------

def mermin_terms(n: int) -> list[tuple[tuple[str, ...], float]]:
    """Return all terms of M_n = Re((X+iY)^⊗n).

    Returns list of (pauli_string, coefficient) where pauli_string is a
    tuple of length n with entries 'X' or 'Y', coefficient is ±1.

    Only even-k terms appear (k = number of Y's). Coefficient = (-1)^(k/2).
    """
    terms = []
    for k in range(0, n + 1, 2):          # k = 0, 2, 4, ...
        coeff = (-1) ** (k // 2)
        for positions in combinations(range(n), k):
            ps = tuple("Y" if i in positions else "X" for i in range(n))
            terms.append((ps, float(coeff)))
    return terms


def classical_bound(n: int) -> float:
    """LHV upper bound: |⟨M_n⟩| ≤ 2^(n/2)."""
    return 2.0 ** (n / 2)


def quantum_maximum(n: int) -> float:
    """GHZ quantum maximum: ⟨M_n⟩ = 2^(n-1)."""
    return 2.0 ** (n - 1)


def num_circuits(n: int) -> int:
    """Number of distinct measurement circuits needed."""
    from math import comb
    return sum(comb(n, k) for k in range(0, n + 1, 2))


# ---------------------------------------------------------------------------
# Circuit construction
# ---------------------------------------------------------------------------

def _basis_rotation(qc: QuantumCircuit, qubit: int, basis: str) -> None:
    """Rotate qubit into measurement basis before computational-basis measurement."""
    if basis == "X":
        qc.h(qubit)
    elif basis == "Y":
        qc.sdg(qubit)
        qc.h(qubit)


def build_mermin_circuits(
    n: int, state_circuit: QuantumCircuit
) -> dict[tuple[str, ...], QuantumCircuit]:
    """Build one measurement circuit per unique Pauli string in M_n.

    state_circuit: prepares the GHZ-n state (without measurement).
    Returns dict mapping pauli_string → QuantumCircuit.
    """
    terms = mermin_terms(n)
    unique_paulis = {ps for ps, _ in terms}

    circuits = {}
    for ps in unique_paulis:
        qc = state_circuit.copy()
        for i, basis in enumerate(ps):
            _basis_rotation(qc, i, basis)
        qc.measure_all()
        circuits[ps] = qc
    return circuits


# ---------------------------------------------------------------------------
# Result computation
# ---------------------------------------------------------------------------

def compute_mermin(
    results: dict[tuple[str, ...], dict[str, int]], n: int
) -> dict:
    """Compute ⟨M_n⟩ from measurement counts.

    results: dict mapping pauli_string → counts dict.
    """
    terms = mermin_terms(n)
    total_qubits = n

    expectation = 0.0
    for ps, coeff in terms:
        counts = results[ps]
        total = sum(counts.values())
        exp_val = 0.0
        for bitstring, cnt in counts.items():
            bits = bitstring.replace(" ", "")
            # Qiskit little-endian: bit i = bits[n-1-i]
            eigenvalue = 1.0
            for bit in bits:
                eigenvalue *= (1 if bit == "0" else -1)
            exp_val += eigenvalue * cnt / total
        expectation += coeff * exp_val

    cb = classical_bound(n)
    qm = quantum_maximum(n)
    violation = abs(expectation) - cb

    return {
        "M_n": expectation,
        "classical_bound": cb,
        "quantum_maximum": qm,
        "violation": violation,
        "violation_ratio": abs(expectation) / cb if cb > 0 else 0.0,
        "n_qubits": n,
        "n_circuits": num_circuits(n),
    }


def run_mermin(
    backend, n: int, state_circuit: QuantumCircuit, shots: int = 4000
) -> dict:
    """Run Mermin-n witness on backend and return results.

    Note: for n ≥ 8, num_circuits(n) = 2^(n-1) > 128, split into batches.
    """
    from qiskit import transpile as qk_transpile

    meas_circuits = build_mermin_circuits(n, state_circuit)
    keys = list(meas_circuits.keys())
    circuit_list = list(meas_circuits.values())

    # Batch in groups of 200 (IQM limit)
    BATCH = 200
    all_counts: list[dict] = []
    for start in range(0, len(circuit_list), BATCH):
        batch = circuit_list[start:start + BATCH]
        transpiled = qk_transpile(batch, backend=backend, optimization_level=3)
        job = backend.run(transpiled, shots=shots)
        raw = job.result().get_counts()
        if isinstance(raw, dict):
            raw = [raw]
        all_counts.extend(raw)

    results = {keys[i]: all_counts[i] for i in range(len(keys))}
    return compute_mermin(results, n)
