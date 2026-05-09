"""Aer-simulator smoke test for the Entanglement Bottleneck Map module.

Checks:
  1. With CZ on every edge, F ≈ 1 on Aer (every edge certified entangling).
  2. With no CZ (control), F ≈ 0.25 → no edge falsely flagged entangled.
  3. Edge coloring produces matchings with no shared qubit.
  4. Bitstring ordering: a hand-crafted state gives the expected sign.
  5. Audit + chain/patch recommenders run without error on the simulator map.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.backend import get_backend
from src.diagnostics import (
    audit_layout,
    edge_color_matchings,
    expectation_from_counts,
    isolated_vs_parallel,
    node_diagnostics,
    recommend_best_chain,
    run_edge_map,
)


def test_edge_coloring_matchings_disjoint():
    edges = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3)]
    matchings = edge_color_matchings(edges)
    for m in matchings:
        seen = set()
        for a, b in m:
            assert a not in seen and b not in seen, f"Matching not disjoint: {m}"
            seen.update((a, b))
    flat = sorted(e for m in matchings for e in m)
    assert flat == sorted(edges), "Edge coloring lost or duplicated edges"
    print(f"  matchings: {matchings}")


def test_bit_ordering_simple():
    # 2 logical qubits; counts where qubit 0 = 1 and qubit 1 = 0 → bitstring "01".
    counts = {"01": 100}
    e = expectation_from_counts(counts, 0, 1, 2)
    # qubit 0 (rightmost char '1') → -1, qubit 1 ('0') → +1, product = -1.
    assert abs(e + 1) < 1e-9, f"Expected -1, got {e}"
    counts = {"00": 50, "11": 50}
    e = expectation_from_counts(counts, 0, 1, 2)
    # both 0 → +1*+1 = +1; both 1 → -1*-1 = +1. Mean = +1.
    assert abs(e - 1) < 1e-9, f"Expected +1, got {e}"


def test_full_map_aer_with_cz():
    backend = get_backend()
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (1, 4)]
    results = run_edge_map(backend, edges=edges, shots=4000,
                            mode="parallel", progress=False)
    assert len(results) == len(edges)
    for r in results:
        assert r.F > 0.95, f"Edge {r.edge} F={r.F:.3f} too low on noiseless Aer"
        assert r.entangled_3sigma, f"Edge {r.edge} not certified on Aer"
    print(f"  with CZ:    F = {[round(r.F, 3) for r in results]}  → all entangled")


def test_no_cz_control_aer():
    backend = get_backend()
    edges = [(0, 1), (1, 2), (2, 3)]
    results = run_edge_map(backend, edges=edges, shots=4000,
                            mode="parallel", no_cz=True, progress=False)
    for r in results:
        assert r.F < 0.4, (
            f"No-CZ edge {r.edge} F={r.F:.3f} suspiciously high (≤0.5 expected)"
        )
        assert not r.entangled_3sigma, f"No-CZ edge {r.edge} falsely entangled"
    print(f"  no CZ:      F = {[round(r.F, 3) for r in results]}  → none entangled")


def test_isolated_vs_parallel_aer_consistency():
    backend = get_backend()
    edges = [(0, 1), (2, 3)]
    par = run_edge_map(backend, edges=edges, shots=4000,
                        mode="parallel", progress=False)
    iso = run_edge_map(backend, edges=edges, shots=4000,
                        mode="isolated", progress=False)
    cmp = isolated_vs_parallel(iso, par)
    for row in cmp:
        assert abs(row["delta_F"]) < 0.05, (
            f"Aer iso vs par delta too large: {row}"
        )
    print(f"  iso vs par delta_F: {[round(r['delta_F'], 3) for r in cmp]}")


def test_audit_and_recommend_aer():
    backend = get_backend()
    edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
    results = run_edge_map(backend, edges=edges, shots=4000,
                            mode="parallel", progress=False)
    nodes = node_diagnostics(results)
    assert set(nodes) == {0, 1, 2, 3, 4}
    audit = audit_layout(
        layout=[0, 1, 2, 3, 4],
        required_logical_edges=[(0, 1), (1, 2), (2, 3), (3, 4)],
        measured=results,
    )
    assert audit["validated"], audit
    chain = recommend_best_chain(results, length=4, require_3sigma=True)
    assert chain is not None and chain["min_F"] > 0.95
    print(f"  recommended best 4-chain qubits: {chain['qubits']}  "
          f"min_F={chain['min_F']:.3f}")


if __name__ == "__main__":
    print("=== edge coloring → disjoint matchings ===")
    test_edge_coloring_matchings_disjoint()
    print("\n=== bitstring ordering sanity ===")
    test_bit_ordering_simple()
    print("\n=== full map on Aer (with CZ) ===")
    test_full_map_aer_with_cz()
    print("\n=== no-CZ control on Aer ===")
    test_no_cz_control_aer()
    print("\n=== isolated vs parallel (Aer should agree) ===")
    test_isolated_vs_parallel_aer_consistency()
    print("\n=== audit + chain recommendation ===")
    test_audit_and_recommend_aer()
    print("\n=== ALL DIAGNOSTIC TESTS PASSED ===")
