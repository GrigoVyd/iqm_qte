"""Aer simulator smoke test for all kept modules."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("=== backend ===")
from src.backend import get_backend
backend = get_backend()
print(f"Backend: {backend}")

print("\n=== cluster 2D circuit ===")
from src.circuits.cluster_2d import build_cluster_2d_no_measure
qc = build_cluster_2d_no_measure(3, 3)
print(f"3x3 cluster: {qc.num_qubits} qubits, depth {qc.depth()}")

print("\n=== rectangular GME witness ===")
from src.witnesses.gme_cluster import run_gme
res = run_gme(backend, qc, 3, 3, shots=4000)
print(f"W = {res['W']:.4f}  bound={res['biseparable_bound']}  "
      f"GME={res['is_gme']}  sig={res['significance_sigma']:.1f}")
assert res['is_gme'], "GME should be certified on simulator"

print("\n=== arbitrary graph state + GME witness ===")
from src.circuits.graph_state import build_graph_state, two_coloring
from src.witnesses.gme_graph import (build_gme_circuits_graph,
                                       compute_gme_witness_graph,
                                       fidelity_lower_bound)
from qiskit import transpile
n = 6
edges = [(i, i+1) for i in range(n-1)]
coloring = two_coloring(n, edges)
state = build_graph_state(n, edges)
circ_a, circ_b = build_gme_circuits_graph(state, coloring)
tr = transpile([circ_a, circ_b], backend, optimization_level=0)
res2 = backend.run(tr, shots=4000).result()
out = compute_gme_witness_graph(res2.get_counts(0), res2.get_counts(1), n, edges, coloring)
fb = fidelity_lower_bound(res2.get_counts(0), res2.get_counts(1), n, edges, coloring)
print(f"6-qubit path graph: W = {out['W']:.4f}  GME={out['is_gme']}  F_lb={fb['F_lower_bound']:.3f}")
assert out['is_gme']

print("\n=== ZNE folding ===")
from src.mitigation.zne import fold_cz_gates, extrapolate_to_zero_noise
folded = fold_cz_gates(state, 3)
n_cz_folded = folded.count_ops().get('cz', 0)
print(f"Folded x3: {n_cz_folded} CZs (originally {state.count_ops().get('cz', 0)})")
assert n_cz_folded == 3 * state.count_ops().get('cz', 0)

print("\n=== DFE on simulator ===")
from src.dfe import direct_fidelity_estimation
dfe = direct_fidelity_estimation(backend, state, n, edges, n_samples=100, shots_per_sample=100)
print(f"DFE F = {dfe['F_estimate']:.4f} +/- {dfe['F_std_err']:.4f}")
assert dfe['F_estimate'] > 0.95, "Noiseless simulator should give F~1"

print("\n=== ALL TESTS PASSED ===")
