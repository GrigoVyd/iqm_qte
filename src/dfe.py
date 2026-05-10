"""Direct Fidelity Estimation (DFE) for stabilizer / graph states.

For a stabilizer state |G⟩ with stabilizer group S = ⟨g_1, …, g_n⟩,

    |G⟩⟨G| = (1/2^n) Σ_{S ∈ S} S

Measure ⟨S⟩ for K randomly-sampled S ∈ S, average. Variance is bounded
for stabilizer states (Flammia & Liu 2011), so K ≈ O(1/ε²) suffices for
precision ε.

Reference: Flammia & Liu, "Direct Fidelity Estimation from Few Pauli
Measurements", Phys. Rev. Lett. 106, 230501 (2011).
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Pauli

from src.circuits.graph_state import neighbours_from_edges


def stabilizer_pauli(
    subset: list[int],
    n: int,
    edges: list[tuple[int, int]],
) -> Pauli:
    """Compute the Pauli string for the stabilizer ∏_{i ∈ subset} g_i.

    g_i = X_i ⊗ ⊗_{j ∈ N(i)} Z_j. Multiplying via Qiskit's Pauli class
    handles the symplectic algebra and phase tracking for us.
    """
    nbrs = neighbours_from_edges(n, edges)
    # Start with identity
    label = ['I'] * n
    p = Pauli(''.join(label))
    for i in subset:
        # Build g_i: X on i, Z on each neighbor
        gi_label = ['I'] * n
        # Qiskit Pauli uses big-endian: leftmost char is qubit n-1
        gi_label[n - 1 - i] = 'X'
        for j in nbrs[i]:
            gi_label[n - 1 - j] = 'Z'
        gi = Pauli(''.join(gi_label))
        p = p.compose(gi)
    return p


def sample_random_stabilizer(
    n: int, edges: list[tuple[int, int]], rng: np.random.Generator
) -> tuple[Pauli, list[int]]:
    """Sample a uniformly-random element of the stabilizer group.

    Returns (Pauli string, subset of generator indices used).
    """
    # Uniform over 2^n stabilizers: each generator independently in/out
    subset = [i for i in range(n) if rng.integers(0, 2)]
    return stabilizer_pauli(subset, n, edges), subset


def build_pauli_measurement_circuit(
    state_circuit: QuantumCircuit, pauli: Pauli
) -> QuantumCircuit:
    """Append basis-rotation gates and measurement to measure ⟨pauli⟩.

    Per qubit:
      X factor → apply H before Z-basis measurement
      Y factor → apply S† H
      Z factor → measure as-is
      I factor → measure (any basis OK; we'll skip in product calculation)
    """
    n = state_circuit.num_qubits
    qc = state_circuit.copy()
    # Pauli.to_label() returns big-endian string; qubit q corresponds to char n-1-q
    label = str(pauli)
    if label.startswith('-i') or label.startswith('+i') or label.startswith('-') or label.startswith('+'):
        # Strip phase prefix; we record it separately
        pass
    label_clean = pauli.to_label()
    # to_label may include phase prefix; strip it
    if label_clean and label_clean[0] in '+-':
        if len(label_clean) > 1 and label_clean[1] == 'i':
            label_clean = label_clean[2:]
        else:
            label_clean = label_clean[1:]
    # Now label_clean has length n; char at position k corresponds to qubit n-1-k
    for k, ch in enumerate(label_clean):
        q = n - 1 - k
        if ch == 'X':
            qc.h(q)
        elif ch == 'Y':
            qc.sdg(q)
            qc.h(q)
        # Z and I: no rotation
    qc.measure_all()
    return qc


def pauli_expval_from_counts(
    counts: dict[str, int], pauli: Pauli, n: int,
) -> float:
    """Estimate ⟨pauli⟩ from measurement counts.

    For each shot, eigenvalue = ∏_(non-I qubits) (-1)^bit × overall_phase.
    Average over shots.
    """
    label_clean = pauli.to_label()
    # Strip phase prefix and detect overall phase
    overall_phase = +1
    if label_clean and label_clean[0] in '+-':
        sign = -1 if label_clean[0] == '-' else +1
        if len(label_clean) > 1 and label_clean[1] == 'i':
            # ±i phase: stabilizer should not have this; warn
            return float('nan')
        else:
            overall_phase = sign
            label_clean = label_clean[1:]

    non_i_qubits = [n - 1 - k for k, ch in enumerate(label_clean) if ch != 'I']

    total = sum(counts.values())
    if total == 0:
        return 0.0
    accum = 0.0
    for bs, cnt in counts.items():
        bs = bs.replace(' ', '')
        v = 1
        for q in non_i_qubits:
            bit = bs[n - 1 - q]
            v *= 1 if bit == '0' else -1
        accum += v * cnt
    return overall_phase * accum / total


def direct_fidelity_estimation(
    backend,
    state_circuit: QuantumCircuit,
    n: int,
    edges: list[tuple[int, int]],
    initial_layout: list[int] | None = None,
    n_samples: int = 200,
    shots_per_sample: int = 100,
    seed: int = 42,
) -> dict:
    """Estimate F(|G⟩, ρ) by sampling K random stabilizers, measuring each,
    and averaging.

    Standard error ~ 1 / √K (for stabilizer states, the per-sample variance
    is ≤ 1, giving Hoeffding-type concentration).

    Returns dict with: F_estimate, F_std_err, n_samples, n_shots_total,
    samples (list of per-sample expvals).
    """
    rng = np.random.default_rng(seed)
    samples = []
    all_circuits = []
    used_paulis = []
    for _ in range(n_samples):
        pauli, subset = sample_random_stabilizer(n, edges, rng)
        used_paulis.append(pauli)
        circ = build_pauli_measurement_circuit(state_circuit, pauli)
        all_circuits.append(circ)

    transpiled = transpile(
        all_circuits, backend=backend,
        initial_layout=initial_layout, optimization_level=3,
    )
    job = backend.run(transpiled, shots=shots_per_sample)
    counts_list = job.result().get_counts()
    if isinstance(counts_list, dict):
        counts_list = [counts_list]

    for cnt, pauli in zip(counts_list, used_paulis):
        expval = pauli_expval_from_counts(cnt, pauli, n)
        samples.append(expval)

    samples_arr = np.array(samples)
    F_est = float(samples_arr.mean())
    F_std_err = float(samples_arr.std(ddof=1) / np.sqrt(len(samples_arr)))
    return {
        "F_estimate": F_est,
        "F_std_err": F_std_err,
        "n_samples": n_samples,
        "shots_per_sample": shots_per_sample,
        "n_shots_total": n_samples * shots_per_sample,
        "samples": samples,
        "job_id": getattr(job, 'job_id', lambda: None)(),
    }
