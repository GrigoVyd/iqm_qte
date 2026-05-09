"""Aer-simulator smoke test for all kept modules."""
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

print("\n=== GME witness ===")
from src.witnesses.gme_cluster import run_gme
res = run_gme(backend, qc, 3, 3, shots=4000)
print(f"W = {res['W']:.4f}  bound={res['biseparable_bound']}  ideal={res['W_ideal']}  "
      f"GME={res['is_gme']}  sig={res['significance_sigma']:.1f}")
assert res['is_gme'], "GME should be certified on simulator"

print("\n=== CHSH (variety layer) ===")
from src.witnesses.chsh import run_chsh
res = run_chsh(backend, shots=4000)
print(f"|S| = {abs(res['S']):.4f}  bound=2  violation={res['violation']:.4f}")
assert res['violation'] > 0

print("\n=== Mermin-3 (variety layer) ===")
from src.witnesses.mermin import run_mermin
from src.circuits.ghz import build_ghz_no_measure
res = run_mermin(backend, 3, build_ghz_no_measure(3), shots=4000)
print(f"M_3 = {res['M_n']:.4f}  classical={res['classical_bound']:.4f}  "
      f"quantum={res['quantum_maximum']:.4f}")
assert res['violation'] > 0

print("\n=== ALL TESTS PASSED ===")
