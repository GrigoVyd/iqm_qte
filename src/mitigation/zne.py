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


def _fit_once(s: np.ndarray, v: np.ndarray, method: str) -> float:
    """One fit, return W(0). Internal helper used by both the public
    extrapolate function and the bootstrap variance estimator.
    """
    if method == "linear":
        coeffs = np.polyfit(s, v, 1)
        return float(coeffs[1])
    elif method == "quadratic":
        coeffs = np.polyfit(s, v, 2)
        return float(coeffs[2])
    elif method == "exponential":
        if np.any(v <= 0):
            raise ValueError("exponential fit requires positive values")
        coeffs = np.polyfit(s, np.log(v), 1)
        return float(np.exp(coeffs[1]))
    else:
        raise ValueError(f"unknown method '{method}'")


def extrapolate_to_zero_noise(
    scales: list[float],
    values: list[float],
    method: str = "linear",
) -> tuple[float, dict]:
    """Fit values vs scales and extrapolate to scale = 0.

    Methods:
      'linear': fit W(α) = a + b·α, return a
      'quadratic': fit W(α) = a + b·α + c·α², return a (needs ≥3 points)
      'exponential': fit W(α) = A·exp(-γ·α), return A
    """
    s = np.array(scales, dtype=float)
    v = np.array(values, dtype=float)
    if method == "quadratic" and len(s) < 3:
        raise ValueError("quadratic fit needs ≥3 points")
    a = _fit_once(s, v, method)
    info = {"method": method, "a": a}
    if method == "linear":
        b = float(np.polyfit(s, v, 1)[0])
        info["b"] = b
    elif method == "quadratic":
        coeffs = np.polyfit(s, v, 2)
        info.update({"b": float(coeffs[1]), "c": float(coeffs[0])})
    elif method == "exponential":
        coeffs = np.polyfit(s, np.log(v), 1)
        info.update({"A": a, "gamma": float(-coeffs[0])})
    return a, info


def extrapolate_with_bootstrap(
    scales: list[float],
    raw_counts_a: list[dict[str, int]],
    raw_counts_b: list[dict[str, int]],
    n: int,
    edges: list[tuple[int, int]],
    coloring: list[int],
    correction_factors: dict[int, float] | None = None,
    layout: list[int] | None = None,
    method: str = "linear",
    n_bootstrap: int = 200,
    seed: int = 42,
) -> dict:
    """ZNE extrapolation with proper uncertainty via shot-level bootstrap.

    Re-samples each scale's measurement counts with replacement `n_bootstrap`
    times, recomputes W per scale, refits, and reports the mean and
    standard deviation of the extrapolated W.

    Returns dict with: W_zne_mean, W_zne_std, W_zne_per_bootstrap,
    sigma_above_bound, method.

    The σ-above-bound here properly propagates extrapolation noise — much
    more honest than the single-measurement formula.
    """
    from src.witnesses.gme_graph import compute_gme_witness_graph

    rng = np.random.default_rng(seed)
    s = np.array(scales, dtype=float)

    # Pre-compute the bitstring/count arrays per scale for resampling
    def to_arrays(counts_dict):
        bitstrings, weights = [], []
        for bs, c in counts_dict.items():
            bitstrings.append(bs)
            weights.append(c)
        return np.array(bitstrings), np.array(weights, dtype=float)

    arrs_a = [to_arrays(c) for c in raw_counts_a]
    arrs_b = [to_arrays(c) for c in raw_counts_b]

    W_extrapolated_samples = []
    for _ in range(n_bootstrap):
        W_per_scale = []
        for i, _ in enumerate(scales):
            # resample with replacement at each scale
            bsA, wA = arrs_a[i]
            bsB, wB = arrs_b[i]
            total_A = int(wA.sum())
            total_B = int(wB.sum())
            sampled_A_idx = rng.choice(len(bsA), size=total_A, replace=True, p=wA/wA.sum())
            sampled_B_idx = rng.choice(len(bsB), size=total_B, replace=True, p=wB/wB.sum())
            cA = {}
            for k in sampled_A_idx:
                cA[bsA[k]] = cA.get(bsA[k], 0) + 1
            cB = {}
            for k in sampled_B_idx:
                cB[bsB[k]] = cB.get(bsB[k], 0) + 1
            res = compute_gme_witness_graph(cA, cB, n, edges, coloring)
            stab = res['stabilizer_values']
            if correction_factors is not None and layout is not None:
                from src.mitigation.parity_qrem import apply_qrem_to_stabilizers_graph
                stab = apply_qrem_to_stabilizers_graph(
                    stab, n, edges, layout, correction_factors)
            W_per_scale.append(sum(stab.values()))

        try:
            W_zero = _fit_once(s, np.array(W_per_scale), method)
            W_extrapolated_samples.append(W_zero)
        except Exception:
            continue

    arr = np.array(W_extrapolated_samples)
    bound = n - 1
    return {
        "method": method,
        "W_zne_mean": float(arr.mean()),
        "W_zne_std":  float(arr.std()),
        "W_zne_per_bootstrap": arr.tolist(),
        "sigma_above_bound": float((arr.mean() - bound) / arr.std()) if arr.std() > 0 else float('inf'),
        "n_bootstrap": len(arr),
    }
