"""Classical shadows for entanglement witnessing.

Protocol (Huang, Kueng, Preskill 2020 + locally-entangled extensions 2024):
  1. Prepare state ρ
  2. For each shadow: apply random single-qubit Clifford U (X, Y, or Z basis)
     then measure in computational basis → get bitstring b
  3. Classical shadow snapshot: σ = U† |b⟩⟨b| U  (local inverse)
  4. Estimate ⟨O⟩ = median-of-means over snapshots

This gives O(log M / ε²) samples for M Pauli observables simultaneously,
independent of system size n.

Implemented here:
  - Random Pauli shadow circuits (X/Y/Z basis per qubit)
  - Estimating arbitrary Pauli expectation values
  - Rényi-2 entropy estimation via purity Tr(ρ_A²)
  - Entanglement detection via purity: S_2 > 0 iff Tr(ρ_A²) < 1
"""

from __future__ import annotations

import random
from collections import defaultdict

import numpy as np
from qiskit import QuantumCircuit


PAULI_BASES = ("X", "Y", "Z")


# ---------------------------------------------------------------------------
# Shadow circuit generation
# ---------------------------------------------------------------------------

def build_shadow_circuits(
    state_circuit: QuantumCircuit,
    num_shadows: int = 1000,
    seed: int | None = None,
) -> tuple[list[QuantumCircuit], list[tuple[str, ...]]]:
    """Generate num_shadows measurement circuits with random Pauli bases.

    Returns (circuits, bases) where bases[i] is a tuple of 'X'/'Y'/'Z' per qubit.
    Each circuit = state_circuit + random single-qubit rotations + measure_all.
    """
    rng = random.Random(seed)
    n = state_circuit.num_qubits
    circuits = []
    bases = []

    for _ in range(num_shadows):
        basis = tuple(rng.choice(PAULI_BASES) for _ in range(n))
        qc = state_circuit.copy()
        for q, b in enumerate(basis):
            if b == "X":
                qc.h(q)
            elif b == "Y":
                qc.sdg(q)
                qc.h(q)
            # Z: no rotation
        qc.measure_all()
        circuits.append(qc)
        bases.append(basis)

    return circuits, bases


def _bit_to_pm1(bit: str) -> int:
    return 1 if bit == "0" else -1


# ---------------------------------------------------------------------------
# Shadow post-processing
# ---------------------------------------------------------------------------

def process_shadow_results(
    results: list[dict[str, int]],
    bases: list[tuple[str, ...]],
) -> list[dict]:
    """Convert raw counts + bases into shadow snapshot records.

    Each snapshot record: {"basis": basis, "outcome": most_common_bitstring, "prob": ...}
    For efficiency, we store the most-likely outcome (or sample from distribution).
    """
    snapshots = []
    for counts, basis in zip(results, bases):
        total = sum(counts.values())
        # Sample one outcome proportionally (expectation over shots)
        # For median-of-means we store all outcomes weighted
        snapshots.append({"basis": basis, "counts": counts, "total": total})
    return snapshots


def estimate_pauli_expval(
    snapshots: list[dict],
    pauli: tuple[str, ...],
) -> float:
    """Estimate ⟨P⟩ for Pauli P = ⊗_i σ_{pauli[i]} using classical shadows.

    Only snapshots where each qubit's basis matches the Pauli contribute.
    (Snapshots with mismatched basis are discarded for that Pauli.)

    For matched snapshots: the estimator is 3^k * ∏_i eigenvalue_i
    where k = number of non-identity Pauli operators.
    """
    n = len(pauli)
    # qubits where Pauli is non-identity
    active = [i for i, p in enumerate(pauli) if p != "I"]
    k = len(active)
    factor = 3 ** k  # correction for random Pauli sampling

    estimates = []
    for snap in snapshots:
        basis = snap["basis"]
        counts = snap["counts"]
        total = snap["total"]
        # Check if basis matches Pauli on all active qubits
        if any(basis[i] != pauli[i] for i in active):
            continue
        # Average eigenvalue product over shots
        exp_val = 0.0
        for bitstring, cnt in counts.items():
            bits = bitstring.replace(" ", "")
            ev = factor
            for i in active:
                ev *= _bit_to_pm1(bits[n - 1 - i])  # Qiskit little-endian
            exp_val += ev * cnt / total
        estimates.append(exp_val)

    if not estimates:
        return 0.0
    return float(np.median(estimates))


def estimate_purity(snapshots: list[dict], subsystem: list[int] | None = None) -> float:
    """Estimate Tr(ρ_A²) for subsystem A via classical shadows.

    Uses the unbiased estimator:
      Tr(ρ_A²) ≈ (3^|A|) * E[(-1)^{b·b'} * 1_{basis match}]
    over pairs of independent snapshots with matching bases on A.

    For a pure state: Tr(ρ_A²) = 1.
    For a maximally mixed subsystem of size m: Tr(ρ_A²) = 1/2^m.
    Rényi-2 entropy: S_2 = -log2(Tr(ρ_A²))
    """
    n = len(snapshots[0]["basis"])
    if subsystem is None:
        subsystem = list(range(n))
    m = len(subsystem)
    factor = 3 ** m

    # Pair up independent snapshots
    half = len(snapshots) // 2
    group1 = snapshots[:half]
    group2 = snapshots[half:2*half]

    purity_estimates = []
    for s1, s2 in zip(group1, group2):
        b1, b2 = s1["basis"], s2["basis"]
        # Bases must match on all subsystem qubits
        if any(b1[i] != b2[i] for i in subsystem):
            continue

        # For each pair of outcomes, compute (-3)^{hamming distance on subsystem}
        c1, t1 = s1["counts"], s1["total"]
        c2, t2 = s2["counts"], s2["total"]

        pair_val = 0.0
        for bits1, cnt1 in c1.items():
            for bits2, cnt2 in c2.items():
                bits1 = bits1.replace(" ", "")
                bits2 = bits2.replace(" ", "")
                # Product of (-1)^{bit_diff} on subsystem
                sign = 1.0
                for i in subsystem:
                    if bits1[n-1-i] != bits2[n-1-i]:
                        sign *= -1
                pair_val += sign * (cnt1/t1) * (cnt2/t2)

        purity_estimates.append(factor * pair_val)

    if not purity_estimates:
        return 1.0
    return float(np.mean(purity_estimates))


def renyi2_entropy(purity: float) -> float:
    """S_2 = -log2(Tr(ρ²)). Positive value proves entanglement across the cut."""
    purity = max(purity, 1e-10)
    return -np.log2(purity)


# ---------------------------------------------------------------------------
# High-level runner
# ---------------------------------------------------------------------------

def run_shadows(
    backend,
    state_circuit: QuantumCircuit,
    num_shadows: int = 1000,
    shots_per_shadow: int = 10,
    subsystems: list[list[int]] | None = None,
    seed: int = 42,
) -> dict:
    """Run classical shadows protocol on backend.

    Due to IQM's 200-circuit batch limit, shadows are split into batches.
    shots_per_shadow: shots per random basis circuit.
    subsystems: list of qubit subsets for purity estimation.
    """
    from qiskit import transpile as qk_transpile

    circuits, bases = build_shadow_circuits(state_circuit, num_shadows, seed=seed)

    # Batch into groups of 200 (IQM limit)
    BATCH = 200
    all_counts: list[dict[str, int]] = []
    for start in range(0, len(circuits), BATCH):
        batch = circuits[start:start+BATCH]
        transpiled = qk_transpile(batch, backend=backend, optimization_level=1)
        job = backend.run(transpiled, shots=shots_per_shadow)
        raw = job.result().get_counts()
        if isinstance(raw, dict):
            raw = [raw]
        all_counts.extend(raw)

    snapshots = process_shadow_results(all_counts, bases)

    results = {"num_shadows": num_shadows, "purities": {}, "renyi2": {}}

    if subsystems is None:
        n = state_circuit.num_qubits
        # Default: half-system cut
        subsystems = [list(range(n // 2))]

    for sub in subsystems:
        key = str(sub)
        p = estimate_purity(snapshots, sub)
        results["purities"][key] = p
        results["renyi2"][key] = renyi2_entropy(p)

    results["snapshots"] = snapshots  # keep for custom Pauli estimation
    return results
