"""Sanity check ZNE on n=8 (already-certified) using new bootstrap sigma."""
from __future__ import annotations
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from qiskit import transpile
from src.backend import get_backend, select_best_tree, get_qubit_metrics
from src.circuits.graph_state import build_graph_state
from src.witnesses.gme_graph import build_gme_circuits_graph
from src.mitigation.parity_qrem import correction_factors_from_metrics
from src.mitigation.zne import fold_cz_gates, extrapolate_with_bootstrap

N = 8
SCALES = [1, 3, 5]
SHOTS = 1500

backend = get_backend(device='emerald')
qm, _ = get_qubit_metrics(backend)
cfs = correction_factors_from_metrics(qm)

# Reuse the previous job result instead of running again
import sys as _sys
JOB_ID = '019e0df5-53e1-78f1-a248-bcf89fb40565'   # n=8 ZNE job from previous run
_reuse = '--rerun' not in _sys.argv

tree = select_best_tree(backend, N)
print(f'Tree (n={N}): qubits={tree["qubits"]}')
print(f'predicted W: {tree["predicted_W"]:.3f} / {N}')

state = build_graph_state(N, tree['logical_edges'])
circ_a, circ_b = build_gme_circuits_graph(state, tree['coloring'])

all_circs = []
for s in SCALES:
    fa = fold_cz_gates(circ_a, s)
    fb = fold_cz_gates(circ_b, s)
    transpiled = transpile([fa, fb], backend=backend,
                            initial_layout=tree['qubits'], optimization_level=3)
    all_circs.extend(transpiled)
    print(f'  scale {s}: {transpiled[0].count_ops().get("cz", 0)} CZs, depth {transpiled[0].depth()}')

if _reuse:
    print(f'\nReusing job {JOB_ID} (no fresh shots).')
    job = backend.retrieve_job(JOB_ID)
    counts = job.result().get_counts()
else:
    print(f'\nSubmitting {len(all_circs)} circuits ...')
    t0 = time.time()
    job = backend.run(all_circs, shots=SHOTS)
    print(f'  Job ID: {job.job_id()}')
    counts = job.result().get_counts()
    print(f'  Done in {time.time()-t0:.1f}s\n')

counts_a = [counts[2*i]   for i in range(len(SCALES))]
counts_b = [counts[2*i+1] for i in range(len(SCALES))]

# Bootstrap ZNE
res = extrapolate_with_bootstrap(
    SCALES, counts_a, counts_b, N,
    tree['logical_edges'], tree['coloring'],
    correction_factors=cfs, layout=tree['qubits'],
    method='linear', n_bootstrap=300,
)
print(f"Bootstrap ZNE (linear, +QREM):")
print(f"  W mean   = {res['W_zne_mean']:.3f}")
print(f"  W std    = {res['W_zne_std']:.3f}")
print(f"  Bound    = {N-1}")
print(f"  sigma above  = {res['sigma_above_bound']:+.2f}")
print(f"  GME      = {'YES' if res['W_zne_mean'] > N-1 else 'no'}")

# Compare to single-measurement sigma at scale=1
import numpy as np
from src.witnesses.gme_graph import compute_gme_witness_graph
from src.mitigation.parity_qrem import apply_qrem_to_stabilizers_graph
res1 = compute_gme_witness_graph(counts_a[0], counts_b[0], N, tree['logical_edges'], tree['coloring'])
stab_mit = apply_qrem_to_stabilizers_graph(res1['stabilizer_values'], N,
                                             tree['logical_edges'], tree['qubits'], cfs)
W_mit = sum(stab_mit.values())
sigma_naive = (W_mit - (N-1)) / (np.sqrt(N) / np.sqrt(SHOTS))
print(f"\nSingle-shot scale=1, +QREM:")
print(f"  W = {W_mit:.3f}  (naive sigma = {sigma_naive:+.2f})")
