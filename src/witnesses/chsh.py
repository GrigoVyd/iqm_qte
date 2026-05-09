"""CHSH inequality witness.

Prepares the singlet state |ψ⟩ = (|01⟩ - |10⟩)/√2 and measures
the CHSH operator S = A1(B1-B2) + A2(B1+B2).

Classical bound: |⟨S⟩| ≤ 2
Quantum maximum: |⟨S⟩| = 2√2 ≈ 2.828

Alice's bases: A1=Z, A2=X
Bob's bases:   B1 = (Z+X)/√2  (45° rotation), B2 = (Z-X)/√2
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit


def build_singlet() -> QuantumCircuit:
    """Prepare the singlet state (|01⟩ - |10⟩)/√2 on 2 qubits."""
    qc = QuantumCircuit(2, name="singlet")
    qc.h(0)
    qc.cx(0, 1)
    qc.x(1)
    qc.z(0)
    return qc


def build_chsh_circuits() -> dict[str, QuantumCircuit]:
    """Return the 4 measurement circuits for CHSH (A1B1, A1B2, A2B1, A2B2).

    Each circuit prepares the singlet, rotates into the measurement basis,
    and measures both qubits.

    Measurement conventions (eigenvalue ±1 from bit 0/1):
      A1 (Z basis): no rotation on qubit 0
      A2 (X basis): H on qubit 0
      B1 (45° in ZX): Ry(-π/4) on qubit 1
      B2 (135° in ZX): Ry(-π/4) + H on qubit 1
    """
    circuits = {}
    for a in (1, 2):
        for b in (1, 2):
            qc = build_singlet()
            if a == 2:
                qc.h(0)                     # Alice: X basis
            qc.ry(-np.pi / 4, 1)            # Bob: rotate 45°
            if b == 2:
                qc.h(1)                     # Bob: second basis
            qc.measure_all()
            circuits[f"A{a}B{b}"] = qc
    return circuits


def corr_from_counts(counts: dict[str, int]) -> float:
    """Compute ZZ correlation ⟨AB⟩ = P(same) - P(different) from counts."""
    total = sum(counts.values())
    same = counts.get("00", 0) + counts.get("11", 0)
    diff = counts.get("01", 0) + counts.get("10", 0)
    return (same - diff) / total


def compute_s_value(results: dict[str, dict[str, int]]) -> float:
    """Compute CHSH S value from counts for all 4 measurement settings.

    S = ⟨A1B1⟩ - ⟨A1B2⟩ + ⟨A2B1⟩ + ⟨A2B2⟩

    Classical bound: |S| ≤ 2. Quantum max: |S| = 2√2 ≈ 2.828.
    """
    c11 = corr_from_counts(results["A1B1"])
    c12 = corr_from_counts(results["A1B2"])
    c21 = corr_from_counts(results["A2B1"])
    c22 = corr_from_counts(results["A2B2"])
    return c11 - c12 + c21 + c22


def violation_significance(s: float, shots_per_setting: int) -> float:
    """Estimate statistical significance (z-score) of CHSH violation.

    Each correlation estimator has std ≈ 1/√shots. S sums 4 correlators,
    so std(S) ≈ 2/√shots.
    """
    std_s = 2.0 / np.sqrt(shots_per_setting)
    violation = abs(s) - 2.0
    return violation / std_s if std_s > 0 else 0.0


def run_chsh(backend, shots: int = 4000, transpile_fn=None) -> dict:
    """Run CHSH on backend, return full results dict.

    transpile_fn: optional callable(circuits, backend) -> transpiled circuits.
    """
    from qiskit import transpile as qk_transpile

    circuits = build_chsh_circuits()
    circuit_list = list(circuits.values())
    keys = list(circuits.keys())

    if transpile_fn:
        circuit_list = transpile_fn(circuit_list, backend)
    else:
        circuit_list = qk_transpile(circuit_list, backend=backend, optimization_level=3)

    job = backend.run(circuit_list, shots=shots)
    raw = job.result().get_counts()
    if isinstance(raw, dict):
        raw = [raw]

    results = {key: raw[i] for i, key in enumerate(keys)}
    s = compute_s_value(results)
    sig = violation_significance(s, shots)

    return {
        "S": s,
        "S_ideal": -2 * np.sqrt(2),
        "classical_bound": 2.0,
        "violation": abs(s) - 2.0,
        "significance_sigma": sig,
        "counts": results,
    }
