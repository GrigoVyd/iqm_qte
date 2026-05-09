import sys
sys.path.insert(0, r"C:\Users\david\Documents\iqm_hackathon")

print("=== Testing backend ===")
from src.backend import get_backend
backend = get_backend()
print(f"Backend: {backend}")

print("\n=== Testing GHZ circuit ===")
from src.circuits.ghz import build_ghz_chain, build_ghz_no_measure
qc = build_ghz_chain(5)
print(f"GHZ-5: {qc.num_qubits} qubits, depth {qc.depth()}")

print("\n=== Testing Cluster 2D circuit ===")
from src.circuits.cluster_2d import build_cluster_2d_no_measure
qc_cluster = build_cluster_2d_no_measure(3, 3)
print(f"Cluster 3x3: {qc_cluster.num_qubits} qubits, depth {qc_cluster.depth()} (expect ~3 ignoring barriers)")

print("\n=== Testing CHSH on simulator ===")
from src.witnesses.chsh import run_chsh
result = run_chsh(backend, shots=8000)
print(f"S = {result['S']:.4f}  (ideal: {result['S_ideal']:.4f}, bound: +/-2)")
print(f"Violation: {result['violation']:.4f}  Significance: {result['significance_sigma']:.1f}sigma")
assert result['violation'] > 0, "CHSH should violate classical bound on simulator!"

print("\n=== Testing Mermin-3 on simulator ===")
from src.witnesses.mermin import run_mermin
ghz3 = build_ghz_no_measure(3)
res = run_mermin(backend, 3, ghz3, shots=8000)
print(f"M_3 = {res['M_n']:.4f}  classical_bound={res['classical_bound']:.4f}  quantum_max={res['quantum_maximum']:.4f}")
print(f"Violation: {res['violation']:.4f}")
assert res['violation'] > 0, "Mermin-3 should violate on simulator!"

print("\n=== Testing GME witness on simulator ===")
from src.witnesses.gme_cluster import run_gme
rows, cols = 3, 3
cluster = build_cluster_2d_no_measure(rows, cols)
res = run_gme(backend, cluster, rows, cols, shots=8000)
print(f"W = {res['W']:.4f}  bisep_bound={res['biseparable_bound']}  ideal={res['W_ideal']}")
print(f"GME certified: {res['is_gme']}  significance: {res['significance_sigma']:.1f}sigma")
assert res['is_gme'], "GME witness should certify cluster state on simulator!"

print("\n=== Testing Classical Shadows on simulator ===")
from src.witnesses.classical_shadows import run_shadows
ghz4 = build_ghz_no_measure(4)
res = run_shadows(backend, ghz4, num_shadows=200, shots_per_shadow=50,
                  subsystems=[[0, 1]], seed=42)
purity = list(res["purities"].values())[0]
s2 = list(res["renyi2"].values())[0]
print(f"Purity Tr(rho_A^2) = {purity:.4f}  Renyi-2 entropy = {s2:.4f}")

print("\n=== ALL TESTS PASSED ===")
