"""Minimal first hardware test: GME witness on a 2x3 cluster state (6 qubits).

Run with:
    python hardware_test_gme.py

Reads IQM_TOKEN from environment. Falls back to Aer simulator if absent.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.backend import get_backend, is_simulator, select_best_subgrid
from src.circuits.cluster_2d import build_cluster_2d_no_measure
from src.witnesses.gme_cluster import build_gme_circuits, compute_gme_witness, gme_significance

ROWS, COLS = 2, 3
SHOTS = 4000
DEVICE = "emerald"


def select_qubits(backend, rows: int, cols: int) -> list[int] | None:
    if is_simulator(backend):
        return list(range(rows * cols))
    try:
        layout, cost, n = select_best_subgrid(backend, rows, cols, readout_mode='fidelity')
        print(f"  Selected ({n} candidates): {layout}  cost={cost:.4f}")
        return layout
    except Exception as e:
        print(f"  Qubit selection failed ({e}); letting Qiskit pick layout.")
        return None


def main() -> int:
    token = os.environ.get("IQM_TOKEN")
    print("=" * 60)
    print(f"GME hardware test: {ROWS}x{COLS} = {ROWS*COLS} qubits, {SHOTS} shots")
    print("=" * 60)

    backend = get_backend(token=token, device=DEVICE)
    print(f"Backend: {backend}")
    print(f"Hardware: {not is_simulator(backend)}")
    print()

    print("[1/4] Selecting qubits ...")
    initial_layout = select_qubits(backend, ROWS, COLS)
    print()

    print("[2/4] Building cluster state circuit (depth 3) ...")
    state_circuit = build_cluster_2d_no_measure(ROWS, COLS)
    print(f"  Qubits: {state_circuit.num_qubits}")
    print(f"  Gates: {dict(state_circuit.count_ops())}")
    print()

    print("[3/4] Building 2 GME measurement circuits ...")
    circ_a, circ_b = build_gme_circuits(state_circuit, ROWS, COLS)
    print(f"  Setting A depth: {circ_a.depth()}")
    print(f"  Setting B depth: {circ_b.depth()}")
    print()

    print("[4/4] Transpiling and running on hardware ...")
    from qiskit import transpile
    kwargs = {"backend": backend, "optimization_level": 3}
    if initial_layout is not None and not is_simulator(backend):
        kwargs["initial_layout"] = initial_layout
    transpiled = transpile([circ_a, circ_b], **kwargs)
    print(f"  Transpiled depths: A={transpiled[0].depth()}, B={transpiled[1].depth()}")
    print(f"  Submitting job ...")
    t0 = time.time()
    job = backend.run(transpiled, shots=SHOTS)
    try:
        job_id = job.job_id()
        print(f"  Job ID: {job_id}")
    except Exception:
        pass
    result = job.result()
    elapsed = time.time() - t0
    print(f"  Job completed in {elapsed:.1f} s")
    print()

    counts = result.get_counts()
    if isinstance(counts, dict):
        counts_a = counts_b = counts
    else:
        counts_a, counts_b = counts[0], counts[1]

    res = compute_gme_witness(counts_a, counts_b, ROWS, COLS)
    n = ROWS * COLS
    sigma = gme_significance(res["W"], n, SHOTS)

    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"  W (measured)         = {res['W']:.4f}")
    print(f"  W ideal              = {res['W_ideal']:.0f}  (perfect cluster state)")
    print(f"  Biseparable bound    = {res['biseparable_bound']}  (W must exceed this)")
    print(f"  Violation            = {res['violation']:.4f}")
    print(f"  Fraction of ideal    = {res['violation_fraction']*100:.1f}%")
    print(f"  Statistical sigma    = {sigma:.2f}")
    print(f"  GME CERTIFIED?       = {'YES' if res['is_gme'] else 'NO'}")
    print()
    print("Per-qubit stabilizer expectation values:")
    for q, ev in res["stabilizer_values"].items():
        print(f"  qubit {q}: <g_{q}> = {ev:+.4f}")
    print()
    return 0 if res["is_gme"] else 1


if __name__ == "__main__":
    sys.exit(main())
