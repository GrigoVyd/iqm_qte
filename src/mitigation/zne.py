"""Zero-Noise Extrapolation (ZNE) for the GME witness.

Idea: run the same logical circuit at multiple noise levels α, observe
W(α), then extrapolate to α = 0.

We scale noise by **CZ gate folding**: replace each CZ U with U·U·U (α=3),
or U·U·U·U·U (α=5). Since CZ² = I, an odd number of repetitions has the
same logical effect as one CZ, but the gate-error budget multiplies by α.

Reference: Temme, Bravyi, Gambetta, "Error mitigation for short-depth
quantum circuits", PRL 119, 180509 (2017).
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit


def fold_cz_gates(circuit: QuantumCircuit, scale_factor: int) -> QuantumCircuit:
    """Return a new circuit where each CZ is repeated `scale_factor` times.

    scale_factor must be a positive odd integer (1, 3, 5, …) so the
    logical effect is unchanged.

    A `barrier` is inserted between consecutive folded copies so a
    later transpilation pass with optimization_level=3 does not
    cancel CZ·CZ pairs back to identity.
    """
    if scale_factor < 1 or scale_factor % 2 == 0:
        raise ValueError("scale_factor must be a positive odd integer")
    new_qc = QuantumCircuit(*circuit.qregs, *circuit.cregs, name=f"{circuit.name}-x{scale_factor}")
    for instr in circuit.data:
        if instr.operation.name == 'cz':
            for k in range(scale_factor):
                if k > 0:
                    new_qc.barrier(*instr.qubits)
                new_qc.append(instr.operation, instr.qubits, instr.clbits)
        else:
            new_qc.append(instr.operation, instr.qubits, instr.clbits)
    return new_qc


def extrapolate_to_zero_noise(
    scales: list[float],
    values: list[float],
    method: str = "linear",
) -> tuple[float, dict]:
    """Fit values vs scales and extrapolate to scale = 0.

    Methods:
      'linear': fit W(α) = a + b·α, return a
      'quadratic': fit W(α) = a + b·α + c·α², return a
      'exponential': fit W(α) = A·exp(-γ·α), return A

    Returns (W_extrapolated, fit_info).
    """
    s = np.array(scales, dtype=float)
    v = np.array(values, dtype=float)
    if method == "linear":
        coeffs = np.polyfit(s, v, 1)         # [b, a]
        a, b = coeffs[1], coeffs[0]
        return float(a), {"a": a, "b": b, "method": "linear"}
    elif method == "quadratic":
        if len(s) < 3:
            raise ValueError("quadratic fit needs ≥3 points")
        coeffs = np.polyfit(s, v, 2)         # [c, b, a]
        a, b, c = coeffs[2], coeffs[1], coeffs[0]
        return float(a), {"a": a, "b": b, "c": c, "method": "quadratic"}
    elif method == "exponential":
        # Fit log(v) = log(A) - γ·α
        if any(x <= 0 for x in v):
            raise ValueError("exponential fit requires positive values")
        log_v = np.log(v)
        coeffs = np.polyfit(s, log_v, 1)     # [-γ, log A]
        gamma = -coeffs[0]
        A = float(np.exp(coeffs[1]))
        return A, {"A": A, "gamma": gamma, "method": "exponential"}
    else:
        raise ValueError(f"unknown method '{method}'")
