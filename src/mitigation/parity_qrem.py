"""Parity readout error mitigation (parity QREM).

For an observable that depends only on the parity of a subset of qubits — like
a stabilizer expectation ⟨g_i⟩ = ⟨X_i ⊗ Z_{N(i)}⟩ — the per-qubit readout error
factorizes. Each qubit i contributes a multiplicative correction:

    c_i = 1 / (P(0|0)_i + P(1|1)_i − 1)

where P(b|b)_i is the probability of correctly reading bit b on qubit i.
The corrected stabilizer expectation is:

    ⟨g_i⟩_mitigated = ⟨g_i⟩_raw × ∏_{j ∈ {i} ∪ N(i)} c_j

References:
- Bravyi, Sheldon et al., "Mitigating measurement errors in multiqubit experiments",
  Phys. Rev. A 103, 042605 (2021).

Calibration data: pulled from the IQM device's quality metric set
(error_0_to_1, error_1_to_0 per qubit) — no extra calibration circuits needed.
"""

from __future__ import annotations


def correction_factors_from_metrics(qubit_metrics: dict[int, dict]) -> dict[int, float]:
    """Per-qubit parity-QREM correction factors c_i from device metrics.

    qubit_metrics: {qubit_idx: {error_0_to_1, error_1_to_0, ...}}.
    Returns {qubit_idx: c_i}; missing qubits get c_i = 1.0 (no correction).
    """
    factors: dict[int, float] = {}
    for idx, m in qubit_metrics.items():
        e01 = m.get("error_0_to_1")
        e10 = m.get("error_1_to_0")
        if e01 is None or e10 is None:
            factors[idx] = 1.0
            continue
        denom = (1 - e01) + (1 - e10) - 1  # = 1 - e01 - e10
        if denom <= 0:
            factors[idx] = 1.0
        else:
            factors[idx] = 1.0 / denom
    return factors


def apply_qrem_to_stabilizers(
    stabilizer_values: dict[int, float],
    rows: int,
    cols: int,
    layout: list[int],
    correction_factors: dict[int, float],
) -> dict[int, float]:
    """Apply parity-QREM to per-qubit stabilizer expectations.

    stabilizer_values: {logical_qubit: ⟨g_i⟩} from compute_gme_witness.
    layout: physical qubit per logical index (logical i → layout[i]).
    correction_factors: {physical_qubit: c}.
    """
    def neighbours(r: int, c: int) -> list[int]:
        nbrs = []
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                nbrs.append(nr * cols + nc)
        return nbrs

    out: dict[int, float] = {}
    for r in range(rows):
        for c in range(cols):
            q_logical = r * cols + c
            involved_logical = [q_logical] + neighbours(r, c)
            involved_physical = [layout[i] for i in involved_logical]
            factor = 1.0
            for p in involved_physical:
                factor *= correction_factors.get(p, 1.0)
            out[q_logical] = stabilizer_values[q_logical] * factor
    return out
