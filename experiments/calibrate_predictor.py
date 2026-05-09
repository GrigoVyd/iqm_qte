"""Calibrate cost-function scale factors against measured per-stabilizer
fidelities from the showcase run.

Model: predicted_⟨g_i⟩ = ∏ (CZ_fid)^α_CZ · ∏ (1 − α_T1 · t/T1) · ∏ (1 − α_T2 · t/T2) · ∏ (F_1Q)^α_1Q

We fit the four α scale factors by minimizing squared error vs measured
⟨g_i⟩ (with QREM applied) on the existing showcase data.
"""
from __future__ import annotations
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

import numpy as np
from scipy.optimize import minimize

from src.backend import (get_backend, get_qubit_metrics, select_best_tree,
                          CZ_GATE_TIME_S)
from src.circuits.graph_state import build_graph_state, neighbours_from_edges
from src.witnesses.gme_graph import (build_gme_circuits_graph,
                                       compute_gme_witness_graph)
from src.mitigation.parity_qrem import (correction_factors_from_metrics,
                                          apply_qrem_to_stabilizers_graph)
from qiskit import transpile

NS = [6, 8, 12, 16, 20]
SHOTS = 1500


def predicted_g_with_alphas(
    qubits, logical_edges, qubit_metrics, cz_fidelities,
    a_cz=1.0, a_t1=1.0, a_t2=1.0, a_1q=1.0,
    t_g=CZ_GATE_TIME_S,
):
    """Same shape as predict_stabilizer_fidelities but with per-channel α."""
    n = len(qubits)
    nbrs = neighbours_from_edges(n, logical_edges)
    out = {}
    for i in range(n):
        involved = [i] + nbrs[i]
        f = 1.0
        for (a, b) in logical_edges:
            if a in involved or b in involved:
                pa, pb = qubits[a], qubits[b]
                cz = cz_fidelities.get((min(pa, pb), max(pa, pb)), 0.99)
                f *= cz**a_cz   # raise CZ fidelity to α_CZ
        for lq in involved:
            pq = qubits[lq]
            mp = qubit_metrics.get(pq, {})
            f *= mp.get("one_q_fidelity", 1.0)**a_1q
            t1 = mp.get("t1") or 1e9
            t2 = mp.get("t2") or 1e9
            f *= max(0.0, 1.0 - a_t1 * t_g / t1 - a_t2 * t_g / t2)
        out[i] = f
    return out


def main():
    # Load showcase results, which has per-qubit measurements
    with open('experiments/showcase_results.json') as f:
        showcase = json.load(f)

    # Re-run on Aer simulator wouldn't help; we need actual measured ⟨g_i⟩.
    # The showcase run JSON only stored per-stabilizer dicts indirectly.
    # For now we re-fetch the original job by ID and recompute per-stab.

    backend = get_backend(device='emerald')
    qm, cz = get_qubit_metrics(backend)
    cfs = correction_factors_from_metrics(qm)

    # Reuse the showcase job's counts via job ID stored in the notebook
    # (hardcoded here from the latest run)
    SHOWCASE_JOB_ID = '019e0dff-bd29-7580-933e-27a975a1c733'
    print(f'Reusing showcase job {SHOWCASE_JOB_ID} ...')
    job = backend.retrieve_job(SHOWCASE_JOB_ID)
    counts_list = job.result().get_counts()

    # Reconstruct (n, scale) → counts mapping (3 scales, 2 settings each, 5 sizes)
    # Order matches the showcase notebook: for n in TARGET_NS, for s in SCALES,
    #   add transpiled[0] (a) and transpiled[1] (b)
    SCALES = [1, 3, 5]
    counts_by = {}
    idx = 0
    for n in NS:
        for s in SCALES:
            counts_by[(n, s, 'a')] = counts_list[idx]; idx += 1
            counts_by[(n, s, 'b')] = counts_list[idx]; idx += 1

    # Build training data: (qubits, logical_edges, measured ⟨g_i⟩+QREM) per stabilizer
    measured = []   # (qubits, logical_edges, qubit_metrics, cz_fidelities, dict_of_g)
    for n in NS:
        tree = select_best_tree(backend, n)
        cA = counts_by[(n, 1, 'a')]; cB = counts_by[(n, 1, 'b')]
        res = compute_gme_witness_graph(cA, cB, n, tree['logical_edges'], tree['coloring'])
        stab_q = apply_qrem_to_stabilizers_graph(
            res['stabilizer_values'], n, tree['logical_edges'], tree['qubits'], cfs)
        measured.append((tree['qubits'], tree['logical_edges'], stab_q))

    n_total_stabs = sum(len(m[2]) for m in measured)
    print(f'Total stabilizer measurements available: {n_total_stabs}')

    def loss(alphas, weight_decay=0.01):
        a_cz, a_t1, a_t2, a_1q = alphas
        total = 0.0
        n_pts = 0
        for (qubits, edges, meas_g) in measured:
            pred = predicted_g_with_alphas(
                qubits, edges, qm, cz,
                a_cz=a_cz, a_t1=a_t1, a_t2=a_t2, a_1q=a_1q)
            for i, p in pred.items():
                m = meas_g[i]
                total += (p - m)**2
                n_pts += 1
        # L2 regularization toward α=1 (no surprises)
        total += weight_decay * sum((a - 1.0)**2 for a in alphas)
        return total / n_pts

    # Initial guess
    x0 = [1.0, 1.0, 1.0, 1.0]
    res = minimize(loss, x0, method='Nelder-Mead',
                    options={'xatol': 1e-4, 'fatol': 1e-6, 'disp': True})
    a_cz, a_t1, a_t2, a_1q = res.x

    print()
    print(f'Optimal scale factors:')
    print(f'  alpha_CZ = {a_cz:.4f}   (1.0 = trust calibration; >1 = errors are larger than reported)')
    print(f'  alpha_T1 = {a_t1:.4f}')
    print(f'  alpha_T2 = {a_t2:.4f}')
    print(f'  alpha_1Q = {a_1q:.4f}')
    print(f'Final RMS error: {np.sqrt(res.fun):.4f}')

    # Compute baseline (α=1) error
    baseline_err = np.sqrt(loss([1, 1, 1, 1], weight_decay=0))
    fitted_err = np.sqrt(loss(res.x, weight_decay=0))
    print(f'\nBaseline (α=1) RMS error  : {baseline_err:.4f}')
    print(f'After calibration RMS    : {fitted_err:.4f}')
    print(f'Improvement              : {(1 - fitted_err/baseline_err)*100:.1f}%')

    # Save calibrated weights
    out = {
        'alpha_CZ': float(a_cz), 'alpha_T1': float(a_t1),
        'alpha_T2': float(a_t2), 'alpha_1Q': float(a_1q),
        'baseline_rms': float(baseline_err),
        'calibrated_rms': float(fitted_err),
        'n_data_points': n_total_stabs,
    }
    with open('experiments/calibrated_alphas.json', 'w') as f:
        json.dump(out, f, indent=2)
    print('\nSaved experiments/calibrated_alphas.json')


if __name__ == '__main__':
    main()
