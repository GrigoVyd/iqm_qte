"""Genuine Multipartite Entanglement (GME) witness for arbitrary connected
bipartite graph states (generalisation of the 2D cluster witness).

For a connected 2-colorable graph G = (V, E) with stabilizers
    g_i = X_i ⊗ ⊗_{j ∈ N(i)} Z_j

the witness is
    W = Σ_i ⟨g_i⟩.

Bounds:
    W ≤ n − 1   for any biseparable state (Tóth & Gühne 2005)
    W = n       for the ideal graph state.

So W > n − 1 ⟹ Genuine Multipartite Entanglement.

Two measurement settings suffice (the 2-coloring trick):
    Setting A: black qubits in X basis, white qubits in Z basis
    Setting B: white qubits in X basis, black qubits in Z basis
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, transpile

from src.circuits.graph_state import neighbours_from_edges


def build_gme_circuits_graph(
    state_circuit: QuantumCircuit,
    coloring: list[int],
) -> tuple[QuantumCircuit, QuantumCircuit]:
    """Two measurement circuits (setting A, setting B) for the graph-state GME witness.

    coloring: list of {0, 1} per qubit (0 = "black", 1 = "white").
    """
    def make(measure_color_in_x: int) -> QuantumCircuit:
        qc = state_circuit.copy()
        for q, c in enumerate(coloring):
            if c == measure_color_in_x:
                qc.h(q)
        qc.measure_all()
        return qc

    return make(0), make(1)


def compute_gme_witness_graph(
    counts_a: dict[str, int],
    counts_b: dict[str, int],
    n: int,
    edges: list[tuple[int, int]],
    coloring: list[int],
) -> dict:
    """Compute W = Σ_i ⟨g_i⟩ from raw measurement counts on an arbitrary graph state."""
    nbrs = neighbours_from_edges(n, edges)

    def bit(bs: str, q: int) -> int:
        return 1 if bs[n - 1 - q] == "0" else -1

    def stab_expval(counts: dict[str, int], qubit: int, neighbours: list[int]) -> float:
        total = sum(counts.values())
        exp = 0.0
        for bs, cnt in counts.items():
            bs = bs.replace(" ", "")
            v = bit(bs, qubit)
            for nb in neighbours:
                v *= bit(bs, nb)
            exp += v * cnt / total
        return exp

    stabilizer_values: dict[int, float] = {}
    for q in range(n):
        if coloring[q] == 0:
            ev = stab_expval(counts_a, q, nbrs[q])
        else:
            ev = stab_expval(counts_b, q, nbrs[q])
        stabilizer_values[q] = ev

    W = sum(stabilizer_values.values())
    bisep_bound = n - 1
    ideal = float(n)

    return {
        "W": W,
        "W_ideal": ideal,
        "biseparable_bound": bisep_bound,
        "violation": W - bisep_bound,
        "violation_fraction": (W - bisep_bound) / ideal,
        "is_gme": W > bisep_bound,
        "n_qubits": n,
        "n_edges": len(set((min(a, b), max(a, b)) for a, b in edges if a != b)),
        "stabilizer_values": stabilizer_values,
    }


def gme_significance_graph(W: float, n: int, shots: int) -> float:
    """Same significance estimator as the rectangular witness: σ ≈ √n / √shots."""
    std_W = np.sqrt(n) / np.sqrt(shots)
    return (W - (n - 1)) / std_W if std_W > 0 else 0.0


def run_gme_graph(backend, state_circuit: QuantumCircuit, n: int,
                  edges: list[tuple[int, int]], coloring: list[int],
                  shots: int = 4000, initial_layout: list[int] | None = None) -> dict:
    """Run the generic GME witness end-to-end (transpile, submit, post-process)."""
    circ_a, circ_b = build_gme_circuits_graph(state_circuit, coloring)
    kwargs = {"backend": backend, "optimization_level": 3}
    if initial_layout is not None:
        kwargs["initial_layout"] = initial_layout
    transpiled = transpile([circ_a, circ_b], **kwargs)
    job = backend.run(transpiled, shots=shots)
    raw = job.result().get_counts()
    counts_a, counts_b = (raw, raw) if isinstance(raw, dict) else (raw[0], raw[1])
    res = compute_gme_witness_graph(counts_a, counts_b, n, edges, coloring)
    res["significance_sigma"] = gme_significance_graph(res["W"], n, shots)
    return res
