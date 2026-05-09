import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

from qiskit import QuantumCircuit, transpile
from iqm.qiskit_iqm import IQMProvider
from iqm.qubit_selector.qubit_selector import (
    CalibrationDataManager,
    CostEvaluator,
    ReadoutMode,
    CostFunction,
)

# ── Configuration ─────────────────────────────────────────────────────────────
N_QUBITS = 30
SHOTS    = 200

# ── Connect to IQM Resonance ──────────────────────────────────────────────────
api_token = os.environ.get("IQM_TOKEN")
if not api_token:
    for p in [Path(__file__).parent / ".secrets" / "iqm_api_key",
              Path(__file__).parent.parent / ".secrets" / "iqm_api_key"]:
        if p.exists():
            api_token = p.read_text().strip()
            print(f"Using token from {p}")
            break
if not api_token:
    raise RuntimeError("Set IQM_TOKEN env var or place token in .secrets/iqm_api_key")
os.environ["IQM_TOKEN"] = api_token

provider = IQMProvider("https://resonance.meetiqm.com", quantum_computer="emerald")
backend  = provider.get_backend()
print(f"Connected to: {backend.name}  ({backend.num_qubits} qubits)")


# ── Circuit preparation (W-state, F-gate / Diker method) ──────────────────────
def _f_gate(qc: QuantumCircuit, ctrl: int, tgt: int, theta: float) -> None:
    qc.ry(-theta, tgt)
    qc.cz(ctrl, tgt)
    qc.ry(theta, tgt)


def build_w_state_n(n: int) -> QuantumCircuit:
    qc = QuantumCircuit(n, name="W{}".format(n))
    qc.x(0)
    for k in range(1, n):
        theta = np.arccos(np.sqrt(1.0 / (n - k + 1)))
        _f_gate(qc, k - 1, k, theta)
    for k in range(n - 1):
        qc.cx(k + 1, k)
    return qc


qc_w = build_w_state_n(N_QUBITS)
print(f"W{N_QUBITS}: depth={qc_w.depth()}, ops={dict(qc_w.count_ops())}")

# ── Layout selection ──────────────────────────────────────────────────────────
qc_for_selector = build_w_state_n(N_QUBITS)
qc_for_selector.measure_all()

calibration_data = CalibrationDataManager().get_calibration_fidelities(backend)

layouts, costs = CostEvaluator(
    backend=backend,
    quantum_circuit=qc_for_selector,
    cost_function=CostFunction.GATE_COST_CZ,
    readoutmode=ReadoutMode.QNDNESS,
    num_trials=2000,
).get_top_layouts(num_layouts=20)

print(f"\nTop 5 layouts (lowest cost = best fidelity):")
for rank, (layout, cost) in enumerate(zip(layouts[:5], costs[:5]), 1):
    names = [backend.index_to_qubit_name(q) for q in layout]
    print(f"  Rank {rank}: {names}  cost={cost:.4f}")

best_layout = list(layouts[0])
best_names  = [backend.index_to_qubit_name(q) for q in best_layout]
print(f"\nSelected qubits: {best_names}")

# ── Build Z-basis and X-basis circuits ────────────────────────────────────────
qc_z = build_w_state_n(N_QUBITS)
qc_z.measure_all()

qc_x = build_w_state_n(N_QUBITS)
qc_x.h(range(N_QUBITS))
qc_x.measure_all()

# ── Transpile onto best layout ────────────────────────────────────────────────
def transpile_on_best(qc: QuantumCircuit) -> QuantumCircuit:
    return transpile(qc, backend=backend, initial_layout=best_layout, optimization_level=3)

qc_z_t = transpile_on_best(qc_z)
qc_x_t = transpile_on_best(qc_x)
print(f"\nZ-basis depth after transpilation: {qc_z_t.depth()}")
print(f"X-basis depth after transpilation: {qc_x_t.depth()}")

# ── Run on hardware ───────────────────────────────────────────────────────────
print(f"\nSubmitting {SHOTS} shots × 2 circuits  (N={N_QUBITS})...")
job_z = backend.run(qc_z_t, shots=SHOTS)
job_x = backend.run(qc_x_t, shots=SHOTS)
print(f"Z-basis job: {job_z.job_id()}")
print(f"X-basis job: {job_x.job_id()}")

counts_z = job_z.result().get_counts()
counts_x = job_x.result().get_counts()
print("Z-basis counts:", counts_z)
print("X-basis counts:", counts_x)

# ── Z-basis analysis ──────────────────────────────────────────────────────────
total_z  = sum(counts_z.values())
probs_z  = {k: v / total_z for k, v in counts_z.items()}
w_bitstrings = {format(1 << i, "0{}b".format(N_QUBITS)) for i in range(N_QUBITS)}
fidelity_z   = sum(probs_z.get(s, 0) for s in w_bitstrings)

print(f"\nZ-basis fidelity: {fidelity_z:.3f}  (ideal: 1.000)")
print("Marginal P(q_i=1)  (ideal: {:.3f}):".format(1 / N_QUBITS))
for qi in range(N_QUBITS):
    p1 = sum(v for k, v in probs_z.items() if k[-(qi + 1)] == "1")
    print(f"  q{qi}: {p1:.3f}")

# ── X-basis analysis ──────────────────────────────────────────────────────────
total_x = sum(counts_x.values())
probs_x = {k: v / total_x for k, v in counts_x.items()}


def expectation_xi_xj(probs, i, j):
    val = 0.0
    for bitstr, prob in probs.items():
        ei = 1 - 2 * int(bitstr[-(i + 1)])
        ej = 1 - 2 * int(bitstr[-(j + 1)])
        val += ei * ej * prob
    return val


corr_matrix = np.zeros((N_QUBITS, N_QUBITS))
for i in range(N_QUBITS):
    for j in range(i + 1, N_QUBITS):
        c = expectation_xi_xj(probs_x, i, j)
        corr_matrix[i, j] = c
        corr_matrix[j, i] = c

ideal_corr = 2 / N_QUBITS
avg_corr   = corr_matrix[np.triu_indices(N_QUBITS, k=1)].mean()
shot_stderr = 1 / np.sqrt(SHOTS)

print(f"\nEntanglement witness W for W{N_QUBITS}:")
print(f"  Classical mixture (separable): 0.0000")
print(f"  Ideal W state (entangled):     {ideal_corr:.4f}")
print(f"  Hardware measurement:          {avg_corr:.4f}")
print(f"  Shot noise std-error:          ±{shot_stderr:.4f}")

# ── Plots ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
vmax = max(ideal_corr * 1.2, 0.01)
im = axes[0].imshow(corr_matrix, cmap="RdYlGn", vmin=-vmax, vmax=vmax, aspect="auto")
axes[0].set_title("Measured <Xi Xj>  (hardware, {} shots)".format(SHOTS))
axes[0].set_xlabel("Qubit j"); axes[0].set_ylabel("Qubit i")
plt.colorbar(im, ax=axes[0], label="Correlator value")

ideal_matrix = np.full((N_QUBITS, N_QUBITS), ideal_corr)
np.fill_diagonal(ideal_matrix, 0)
im2 = axes[1].imshow(ideal_matrix, cmap="RdYlGn", vmin=-vmax, vmax=vmax, aspect="auto")
axes[1].set_title("Ideal <Xi Xj> = 2/N = {:.3f}  (W{})".format(ideal_corr, N_QUBITS))
axes[1].set_xlabel("Qubit j"); axes[1].set_ylabel("Qubit i")
plt.colorbar(im2, ax=axes[1], label="Correlator value")

plt.suptitle(
    "X-basis pairwise correlators — W{} state\n"
    "Classical: 0 everywhere; Quantum: {:.3f} off-diagonal".format(N_QUBITS, ideal_corr),
    fontsize=11, y=1.02,
)
plt.tight_layout()
plt.savefig("w{}_correlators.png".format(N_QUBITS), dpi=100, bbox_inches="tight")
plt.show()
print("Saved w{}_correlators.png".format(N_QUBITS))
