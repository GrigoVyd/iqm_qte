"""W-state observables: Z-basis fidelity to |W_n>, X-basis non-linear witness.

Z-basis: P(single excitation, summed over all n positions) = sum_k P(|0..1_k..0>).
For an ideal W-state this equals 1; for the classical mixture of single
excitations it also equals 1, so this is a *necessary* but not sufficient
fidelity proxy. Combined with the X-witness below it pins down quantum
behaviour.

X-basis: <X_i X_j> for i != j is +2/n on |W_n>, 0 on the classical
single-excitation mixture. The mean over all pairs is the witness; values
significantly above 0 prove non-classical coherence.
"""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np


def _strip(bs: str) -> str:
    return bs.replace(" ", "")


def parse_counts(counts: dict[str, int], n: int) -> dict[str, int]:
    """Normalise IQM bitstrings: split on spaces, take last register, trim to n bits."""
    out: dict[str, int] = {}
    for bs, c in counts.items():
        parts = bs.strip().split()
        b = parts[-1][-n:]
        out[b] = out.get(b, 0) + c
    return out


def z_fidelity(counts: dict[str, int], n: int) -> float:
    """Sum of probabilities of single-excitation strings (Hamming weight 1)."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    s = 0
    for bs, c in counts.items():
        b = _strip(bs)
        if len(b) != n:
            continue
        if b.count("1") == 1:
            s += c
    return s / total


def pair_correlators(counts: dict[str, int], n: int) -> list[float]:
    """Return the n*(n-1)/2 ordered-pair <X_i X_j> values from raw counts.

    Counts are assumed to be in the X basis (i.e. a Hadamard layer was
    applied before the Z-measurement that produced these bits).
    """
    total = sum(counts.values())
    if total == 0:
        return [0.0] * (n * (n - 1) // 2)
    out: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            v = 0
            for bs, c in counts.items():
                b = _strip(bs)
                ei = 1 - 2 * int(b[-(i + 1)])
                ej = 1 - 2 * int(b[-(j + 1)])
                v += ei * ej * c
            out.append(v / total)
    return out


def x_witness(counts: dict[str, int], n: int,
               return_pairs: bool = False
               ) -> float | tuple[float, list[float]]:
    """Mean pairwise X-correlator. >0 (above shot noise) ⇒ non-classical.

    If `return_pairs=True` also returns the per-pair list, useful for
    computing a conservative variance estimate.
    """
    pairs = pair_correlators(counts, n)
    mean = float(np.mean(pairs)) if pairs else 0.0
    if return_pairs:
        return mean, pairs
    return mean


def witness_significance(counts: dict[str, int], n: int) -> dict:
    """Return witness mean + conservative sigma + sigma above 0.

    The conservative sigma is the sample-std across pairs / sqrt(n_pairs);
    this absorbs the inter-pair correlations from sharing qubits.
    """
    mean, pairs = x_witness(counts, n, return_pairs=True)
    if not pairs:
        return {"mean": 0.0, "sigma": float("inf"), "z": 0.0,
                "ideal_2_over_n": 2.0 / n, "n_pairs": 0}
    sigma = float(np.std(pairs)) / math.sqrt(len(pairs))
    z = mean / sigma if sigma > 0 else float("inf")
    return {"mean": mean, "sigma": sigma, "z": z,
            "ideal_2_over_n": 2.0 / n, "n_pairs": len(pairs)}
