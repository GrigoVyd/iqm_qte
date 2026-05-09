"""Comprehensive showcase on IQM Garnet — same pipeline as 05_showcase.ipynb.

Runs the full mitigation stack (raw, +QREM, +ZNE, +QREM+ZNE) for several n
values in a single batched job. Saves results to garnet_showcase_results.json
for direct comparison with the Emerald showcase.
"""
from __future__ import annotations
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

import numpy as np
from qiskit import transpile
from src.backend import (get_backend, get_qubit_metrics,
                          filter_qubits_by_metrics, select_best_tree)
from src.circuits.graph_state import build_graph_state, neighbours_from_edges
from src.witnesses.gme_graph import (build_gme_circuits_graph,
                                       compute_gme_witness_graph,
                                       fidelity_lower_bound)
from src.mitigation.parity_qrem import (correction_factors_from_metrics,
                                          apply_qrem_to_stabilizers_graph)
from src.mitigation.zne import fold_cz_gates, extrapolate_with_bootstrap

NS = [6, 8, 12, 16, 20]
SCALES = [1, 3, 5]
SHOTS = 1500

backend = get_backend(device='garnet')
qm, cz = get_qubit_metrics(backend)
cfs = correction_factors_from_metrics(qm)
print(f'Backend: {backend}, {backend.num_qubits} qubits, {len(cz)} CZ pairs')
good = filter_qubits_by_metrics(qm, cz)
print(f'Qubits passing thresholds: {len(good)}/{backend.num_qubits}')

# Pick trees
trees = {}
for n in NS:
    t = select_best_tree(backend, n)
    trees[n] = t
    print(f'  n={n}: cost={t["weight"]:.4f} predW={t["predicted_W"]:.3f} '
          f'qubits={t["qubits"]}')

# Build all circuits and submit as ONE job
all_circs = []
descriptors = []
for n in NS:
    t = trees[n]
    state = build_graph_state(n, t['logical_edges'])
    circ_a, circ_b = build_gme_circuits_graph(state, t['coloring'])
    for s in SCALES:
        fa = fold_cz_gates(circ_a, s); fb = fold_cz_gates(circ_b, s)
        tr = transpile([fa, fb], backend=backend, initial_layout=t['qubits'],
                       optimization_level=3)
        all_circs.extend(tr)
        descriptors.append((n, s, 'a')); descriptors.append((n, s, 'b'))

print(f'\nSubmitting {len(all_circs)} circuits to Garnet ...')
t0 = time.time()
job = backend.run(all_circs, shots=SHOTS)
print(f'  Job ID: {job.job_id()}')
counts_list = job.result().get_counts()
print(f'  Done in {time.time()-t0:.1f}s\n')

# Group counts by (n, scale, ab)
counts_by = {}
for i, (n, s, ab) in enumerate(descriptors):
    counts_by.setdefault((n, s), {})[ab] = counts_list[i]

# Process every n × every mitigation
results = []
for n in NS:
    t = trees[n]; edges = t['logical_edges']; coloring = t['coloring']
    qubits = t['qubits']
    cA_list = [counts_by[(n, s)]['a'] for s in SCALES]
    cB_list = [counts_by[(n, s)]['b'] for s in SCALES]

    # Single-scale W (raw and +QREM)
    res_s1 = compute_gme_witness_graph(cA_list[0], cB_list[0], n, edges, coloring)
    W_raw = res_s1['W']
    stab_q = apply_qrem_to_stabilizers_graph(
        res_s1['stabilizer_values'], n, edges, qubits, cfs)
    W_qrem = sum(stab_q.values())

    # ZNE: bootstrap with and without QREM
    bs_zne = extrapolate_with_bootstrap(
        SCALES, cA_list, cB_list, n, edges, coloring,
        correction_factors=None, layout=None, method='linear', n_bootstrap=200)
    bs_qz = extrapolate_with_bootstrap(
        SCALES, cA_list, cB_list, n, edges, coloring,
        correction_factors=cfs, layout=qubits, method='linear', n_bootstrap=200)

    # Fidelity lower bound
    fb = fidelity_lower_bound(cA_list[0], cB_list[0], n, edges, coloring)

    bound = n - 1
    sn_std = np.sqrt(n) / np.sqrt(SHOTS)
    results.append({
        'n': n, 'qubits': qubits, 'logical_edges': edges,
        'predicted_W': t['predicted_W'], 'tree_cost': t['weight'],
        'bound': bound,
        'W_raw': W_raw, 'sigma_raw': (W_raw-bound)/sn_std,
        'W_qrem': W_qrem, 'sigma_qrem': (W_qrem-bound)/sn_std,
        'W_zne': bs_zne['W_zne_mean'], 'sigma_zne': bs_zne['sigma_above_bound'],
        'W_zne_std': bs_zne['W_zne_std'],
        'W_qz': bs_qz['W_zne_mean'],   'sigma_qz': bs_qz['sigma_above_bound'],
        'W_qz_std': bs_qz['W_zne_std'],
        'is_gme_raw': W_raw>bound, 'is_gme_qrem': W_qrem>bound,
        'is_gme_zne': bs_zne['W_zne_mean']>bound,
        'is_gme_qz':  bs_qz['W_zne_mean']>bound,
        'F_lower_bound': fb['F_lower_bound'], 'P_A': fb['P_A'], 'P_B': fb['P_B'],
    })
    r = results[-1]
    print(f'\nn={n:>2}  bound={bound}')
    print(f'  raw      : W={W_raw:6.3f}  sigma={r["sigma_raw"]:+6.2f}  GME={"YES" if r["is_gme_raw"] else "no"}')
    print(f'  +QREM    : W={W_qrem:6.3f}  sigma={r["sigma_qrem"]:+6.2f}  GME={"YES" if r["is_gme_qrem"] else "no"}')
    print(f'  +ZNE     : W={r["W_zne"]:6.3f} +/-{r["W_zne_std"]:.2f}  sigma={r["sigma_zne"]:+6.2f}  GME={"YES" if r["is_gme_zne"] else "no"}')
    print(f'  +QREM+ZNE: W={r["W_qz"]:6.3f} +/-{r["W_qz_std"]:.2f}  sigma={r["sigma_qz"]:+6.2f}  GME={"YES" if r["is_gme_qz"] else "no"}')
    print(f'  F lower  : {r["F_lower_bound"]:.3f}  (P_A={r["P_A"]:.3f}, P_B={r["P_B"]:.3f})')

print('\n=== SUMMARY (Garnet) ===')
for k, lbl in [('is_gme_raw', 'raw      '), ('is_gme_qrem', '+QREM    '),
               ('is_gme_zne', '+ZNE     '), ('is_gme_qz',  '+QREM+ZNE')]:
    m = max((r['n'] for r in results if r[k]), default=0)
    print(f'  Max n certified, {lbl}: {m}')

with open('experiments/garnet_showcase_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print('\nSaved experiments/garnet_showcase_results.json')
