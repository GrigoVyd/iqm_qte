# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Challenge Overview

**ETH Quantum Hackathon 2026 — "Witness My Entanglement!"** (IQM Challenge)

Prove quantum entanglement on IQM hardware in the most compelling, scalable, and flexible way possible. Three goals:
1. Prove entanglement exists for at least one quantum state (rigorously — no classical explanation)
2. Demonstrate entanglement across many qualitatively distinct states
3. Entangle as many qubits as possible

## Scoring Weights

| Criterion | Weight |
|-----------|--------|
| Theoretical correctness | 20% |
| Implementation sophistication (hardware optimization) | 30% |
| Number of qubits entangled | 10–40% (team chooses) |
| Flexibility / variety of states | 10–40% (team chooses) |

Criteria 3+4 must sum to 40% and each be ≥10%. **Bonus**: scalability to 100s of qubits, identifying bottlenecks, loophole mitigation.

## Hardware

**IQM Resonance** cloud platform. Available devices:
- **Emerald** — 54 qubits, square lattice topology
- **Garnet** — square lattice topology
- **Sirius** — star topology, NOT supported by CUDA-Q

Qubit naming: IQM labels qubits `QB1`, `QB2`, … (1-indexed). Qiskit maps these to 0-indexed integers (`QB1` → index `0`, `QB5` → index `4`). Always use `backend.qubit_name_to_index(name)` to convert.

Check queue lengths and per-device fidelity metrics on the Resonance website before submitting jobs. Credits are limited.

## Environment Setup

**Python 3.11 or 3.12 required** (not latest Python).

```bash
pip install --upgrade pip
pip install "iqm-client[qiskit]"   # Qiskit + IQM
pip install "qrisp[iqm]"           # Qrisp + IQM (optional)
pip install numpy==2.0 matplotlib pylatexenc
```

## Connecting to Hardware

**Qiskit (primary framework):**
```python
from iqm.qiskit_iqm import IQMProvider
provider = IQMProvider("https://resonance.meetiqm.com",
                       quantum_computer="emerald",  # or "garnet"
                       token="YOUR_API_TOKEN")
backend = provider.get_backend()
```

**Qrisp (alternative):**
```python
from qrisp.interface import IQMBackend
quantum_computer = IQMBackend(api_token="YOUR_API_TOKEN", device_instance="garnet")
```

API token: Resonance website → initials (top right) → API token → refresh icon. Token is per-user and can be regenerated if lost (old one deactivates).

## Key Patterns

**Transpile for device topology — always use optimization_level=3:**
```python
from qiskit import transpile
qc_transpiled = transpile(qc, backend=backend, optimization_level=3)
job = backend.run(qc_transpiled, shots=10000)
counts = job.result().get_counts()
```

**Target specific qubits (select best-performing ones from Resonance metrics):**
```python
qubits = [backend.qubit_name_to_index(name) for name in ["QB5", "QB6", "QB10"]]
reduced_coupling_map = [list(edge) for edge in backend.coupling_map
                        if set(edge).issubset(set(qubits))]
qc_transpiled = transpile(qc, backend, coupling_map=reduced_coupling_map, optimization_level=3)
```

**Batch multiple circuits in one job (saves credits and queue time):**
```python
job = backend.run(circuit_list, shots=1000)
# Retrieve later by job ID:
job_id = job.job_id()
retrieved = backend.retrieve_job(job_id)
```

**Visualize coupling map:**
```python
from rustworkx.visualization import mpl_draw
mpl_draw(backend.coupling_map.graph, arrows=True, with_labels=True, node_color='#32a8a4')
```

## Entanglement Theory

**CHSH inequality** — two-qubit entanglement witness:
- Classical bound: |⟨S⟩| ≤ 2
- Quantum maximum: |⟨S⟩| = 2√2 ≈ 2.828
- Violation proves non-classical correlations (entanglement)
- Uses singlet state |ψ⟩ = (1/√2)(|01⟩ − |10⟩) and four measurement basis combinations

**GHZ states** — multipartite entanglement:
- |GHZ_n⟩ = (1/√2)(|00…0⟩ + |11…1⟩)
- Cannot factor out ANY qubit's state → true multipartite entanglement
- Circuit: H on first qubit, then CNOT chain to all others
- Must respect coupling map — not all qubits directly connected

**Other valid methods**: Mermin inequalities (generalize CHSH to n qubits), entanglement witnesses, quantum games, W states.

**Key distinction**: GHZ-like states from `fake_GHZ = (1/2)(|000⟩+|001⟩+|110⟩+|111⟩)` only have pairwise entanglement — must prove genuine multipartite entanglement.

## Project Structure

All challenge materials are in `ETHQHack2026-main/`:
- `README.ipynb` — full challenge description, CHSH example code, scoring criteria
- `connecting_to_Resonance.ipynb` — connection snippets for all 4 frameworks (Qrisp, CUDA-Q, Qiskit, Cirq)
- `Routing_to_specific_qubits.ipynb` — GHZ state with topology-aware qubit selection

New code should go in the project root or a `src/` directory. Use Jupyter notebooks for hardware experiments (easier job ID retrieval); Python scripts for reusable utilities.

## Submission

Google Form submission (link in README.ipynb). Deliver: code + key results. Multiple submissions allowed.
