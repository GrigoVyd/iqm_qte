"""Compare 'old' (threshold + CZ-only) vs 'new' (no threshold + quality weight)
tree-selection strategies in the SAME hardware job, so calibration drift is
shared and the comparison is fair.

Submits 8 circuits in a single job (4 sizes × 2 strategies × 2 settings/size)
divided into 2 jobs (one per strategy) but back-to-back.
"""
from __future__ import annotations
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

import numpy as np
from qiskit import transpile

from src.backend import get_backend, select_best_tree, get_qubit_metrics
from src.circuits.graph_state import build_graph_state, neighbours_from_edges
from src.witnesses.gme_graph import build_gme_circuits_graph, compute_gme_witness_graph, gme_significance_graph
from src.mitigation.parity_qrem import correction_factors_from_metrics

NS = [12, 16, 20]
SHOTS = 1000   # smaller to keep cost down

OLD_KW = dict(apply_thresholds=True,  min_readout_fidelity=0.90, min_one_q_fidelity=0.99,
              min_t1=10e-6, min_t2=5e-6, min_cz_fidelity=0.90,
              readout_weight=0.0, use_local_search=False)
NEW_KW = dict(apply_thresholds=False, readout_weight=0.5, use_local_search=True)

backend = get_backend(device='emerald')
qm, _ = get_qubit_metrics(backend)
cfs = correction_factors_from_metrics(qm)

# Pick trees with both strategies
trees = {}
for n in NS:
    trees[n] = {'old': select_best_tree(backend, n, **OLD_KW),
                'new': select_best_tree(backend, n, **NEW_KW)}
    print(f'n={n}: old qubits={trees[n]["old"]["qubits"]}')
    print(f'      new qubits={trees[n]["new"]["qubits"]}')
    print(f'      same? {sorted(trees[n]["old"]["qubits"]) == sorted(trees[n]["new"]["qubits"])}')

# Build all circuits
all_circuits, descriptors = [], []
for n in NS:
    for strat in ['old', 'new']:
        t = trees[n][strat]
        state = build_graph_state(n, t['logical_edges'])
        circ_a, circ_b = build_gme_circuits_graph(state, t['coloring'])
        tr = transpile([circ_a, circ_b], backend=backend,
                       initial_layout=t['qubits'], optimization_level=3)
        all_circuits.extend(tr)
        descriptors.append((n, strat, t))
        descriptors.append((n, strat, t))

print(f'\nSubmitting {len(all_circuits)} circuits in one job ...')
t0 = time.time()
job = backend.run(all_circuits, shots=SHOTS)
print(f'  Job ID: {job.job_id()}')
counts_list = job.result().get_counts()
print(f'  Done in {time.time()-t0:.1f}s\n')

# Group counts by tree (every 2 are A/B for one tree)
results = []
for i in range(0, len(all_circuits), 2):
    n, strat, t = descriptors[i]
    cA, cB = counts_list[i], counts_list[i+1]
    res = compute_gme_witness_graph(cA, cB, n, t['logical_edges'], t['coloring'])
    sigma_raw = gme_significance_graph(res['W'], n, SHOTS)

    # QREM
    nbrs = neighbours_from_edges(n, t['logical_edges'])
    stab_mit = {}
    for q in range(n):
        involved = [q] + nbrs[q]
        factor = 1.0
        for lq in involved:
            factor *= cfs.get(t['qubits'][lq], 1.0)
        stab_mit[q] = res['stabilizer_values'][q] * factor
    W_mit = sum(stab_mit.values())
    sigma_mit = gme_significance_graph(W_mit, n, SHOTS)

    raw_arr = np.array(list(res['stabilizer_values'].values()))
    mit_arr = np.array(list(stab_mit.values()))
    results.append({
        'n': n, 'strategy': strat,
        'W_raw': res['W'], 'W_mit': W_mit,
        'sigma_raw': sigma_raw, 'sigma_mit': sigma_mit,
        'g_raw_mean': float(raw_arr.mean()),
        'g_mit_mean': float(mit_arr.mean()),
        'is_gme_raw': res['is_gme'], 'is_gme_mit': W_mit > n - 1,
    })

print(f"{'n':>4} {'strat':>6} {'W raw':>7} {'sig raw':>8} {'g raw':>7} {'W mit':>7} {'sig mit':>8} {'g mit':>7} {'GME':>5}")
print('-' * 70)
for r in results:
    gme = 'YES' if r['is_gme_mit'] else 'no '
    print(f"{r['n']:>4} {r['strategy']:>6} {r['W_raw']:>7.2f} {r['sigma_raw']:>8.2f} "
          f"{r['g_raw_mean']:>7.3f} {r['W_mit']:>7.2f} {r['sigma_mit']:>8.2f} "
          f"{r['g_mit_mean']:>7.3f} {gme:>5}")

with open('strategy_comparison.json', 'w') as f:
    json.dump(results, f, indent=2)
print('\nSaved strategy_comparison.json')
