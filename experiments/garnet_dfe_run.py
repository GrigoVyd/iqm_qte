"""DFE for all sizes on Garnet, batched into a single IQM job.

5 sizes x 100 stabilizer samples x 100 shots = 50K total shots, ~18 tokens.
"""
from __future__ import annotations
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

import numpy as np
from qiskit import transpile
from src.backend import get_backend, select_best_tree
from src.circuits.graph_state import build_graph_state
from src.dfe import sample_random_stabilizer, build_pauli_measurement_circuit, pauli_expval_from_counts

NS = [6, 8, 12, 16, 20]
N_SAMPLES = 100
SHOTS = 100
SEED = 42

backend = get_backend(device='garnet')
print(f'Backend: {backend}, {backend.num_qubits} qubits')

# Build all circuits across all sizes -- one job
trees = {}
all_circuits = []
descriptors = []   # (n, sample_idx, pauli)
rng = np.random.default_rng(SEED)
for n in NS:
    tree = select_best_tree(backend, n)
    trees[n] = tree
    state = build_graph_state(n, tree['logical_edges'])
    for k in range(N_SAMPLES):
        pauli, _ = sample_random_stabilizer(n, tree['logical_edges'], rng)
        circ = build_pauli_measurement_circuit(state, pauli)
        all_circuits.append(circ)
        descriptors.append((n, k, pauli))

print(f'Building {len(all_circuits)} measurement circuits across n in {NS} ...')
# Transpile in chunks per size (different initial_layouts)
transpiled = []
idx = 0
for n in NS:
    chunk = all_circuits[idx:idx + N_SAMPLES]
    tr = transpile(chunk, backend=backend, initial_layout=trees[n]['qubits'],
                   optimization_level=3)
    transpiled.extend(tr)
    idx += N_SAMPLES

print(f'Submitting {len(transpiled)} circuits in one job, {SHOTS} shots/circuit ...')
t0 = time.time()
job = backend.run(transpiled, shots=SHOTS)
print(f'  Job ID: {job.job_id()}')
counts_list = job.result().get_counts()
elapsed = time.time() - t0
print(f'  Done in {elapsed:.1f}s\n')

# Process: bin samples per n, average the per-stabilizer expectations
results = []
for ni, n in enumerate(NS):
    expvals = []
    for k in range(N_SAMPLES):
        i = ni * N_SAMPLES + k
        n2, k2, pauli = descriptors[i]
        assert n2 == n and k2 == k
        ev = pauli_expval_from_counts(counts_list[i], pauli, n)
        if not np.isnan(ev):
            expvals.append(ev)
    arr = np.array(expvals)
    F_est = float(arr.mean())
    F_std_err = float(arr.std(ddof=1) / np.sqrt(len(arr)))
    results.append({
        'n': n,
        'F_estimate': F_est, 'F_std_err': F_std_err,
        'n_samples_used': len(arr),
        'samples': expvals,
    })
    print(f'n={n:>2}: F = {F_est:.4f} +/- {F_std_err:.4f}  '
          f'({len(arr)}/{N_SAMPLES} samples)')

print('\n=== Garnet DFE summary ===')
for r in results:
    print(f'  n={r["n"]:>2}: F = {r["F_estimate"]:.4f} +/- {r["F_std_err"]:.4f}')

with open('experiments/garnet_dfe_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print('\nSaved experiments/garnet_dfe_results.json')
