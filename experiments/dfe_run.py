"""Run Direct Fidelity Estimation on hardware for our spanning trees.

Costs ~ K * shots_per_sample shots per size. Default: 200 * 100 = 20K shots.
"""
from __future__ import annotations
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from src.backend import get_backend, select_best_tree, get_qubit_metrics
from src.circuits.graph_state import build_graph_state
from src.witnesses.gme_graph import compute_gme_witness_graph, fidelity_lower_bound, build_gme_circuits_graph
from src.mitigation.parity_qrem import correction_factors_from_metrics, apply_qrem_to_stabilizers_graph
from src.dfe import direct_fidelity_estimation

NS = [8, 20]                  # sizes to estimate fidelity for
N_SAMPLES = 200               # DFE samples per size
SHOTS_PER_SAMPLE = 100

backend = get_backend(device='emerald')
qm, _ = get_qubit_metrics(backend)
cfs = correction_factors_from_metrics(qm)

results = []
for n in NS:
    print(f'\n=== n = {n} ===', flush=True)
    tree = select_best_tree(backend, n)
    state = build_graph_state(n, tree['logical_edges'])

    t0 = time.time()
    dfe = direct_fidelity_estimation(
        backend, state, n, tree['logical_edges'],
        initial_layout=tree['qubits'],
        n_samples=N_SAMPLES, shots_per_sample=SHOTS_PER_SAMPLE,
    )
    elapsed = time.time() - t0

    # Also compute fidelity lower bound from a separate 2-circuit measurement
    # (the graph-state witness measurement). This shares no shots with DFE.
    circ_a, circ_b = build_gme_circuits_graph(state, tree['coloring'])
    from qiskit import transpile
    tr = transpile([circ_a, circ_b], backend=backend,
                    initial_layout=tree['qubits'], optimization_level=3)
    job_lb = backend.run(tr, shots=1500)
    raw = job_lb.result().get_counts()
    cA, cB = raw if not isinstance(raw, list) else (raw[0], raw[1])
    lb = fidelity_lower_bound(cA, cB, n, tree['logical_edges'], tree['coloring'])

    results.append({
        'n': n,
        'F_estimate': dfe['F_estimate'],
        'F_std_err': dfe['F_std_err'],
        'F_lower_bound': lb['F_lower_bound'],
        'P_A': lb['P_A'], 'P_B': lb['P_B'],
        'n_samples': dfe['n_samples'],
        'shots_per_sample': dfe['shots_per_sample'],
        'elapsed_s': elapsed,
    })
    print(f'  DFE: F = {dfe["F_estimate"]:.4f} +/- {dfe["F_std_err"]:.4f}'
          f'  ({dfe["n_samples"]} samples, {dfe["shots_per_sample"]} shots each)')
    print(f'  Lower bound (Toth-Guhne): F >= {lb["F_lower_bound"]:.4f}'
          f'  (P_A={lb["P_A"]:.3f}, P_B={lb["P_B"]:.3f})')
    print(f'  DFE elapsed: {elapsed:.1f}s')

print('\n' + '=' * 60)
print('Summary')
print('=' * 60)
for r in results:
    print(f"  n={r['n']:>2}: F = {r['F_estimate']:.4f} +/- {r['F_std_err']:.4f}"
          f"  (lower bound F >= {r['F_lower_bound']:.3f})")

with open('experiments/dfe_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print('\nSaved experiments/dfe_results.json')
