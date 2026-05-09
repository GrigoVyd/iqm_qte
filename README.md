# Witness My Entanglement — ETH Quantum Hackathon 2026

**IQM Challenge: Prove quantum entanglement on IQM hardware in the most compelling, scalable, and flexible way possible.**

## Headline result

**Genuine Multipartite Entanglement certified on IQM Emerald** for a 6-qubit 2D cluster state:

| | |
|---|---|
| Witness W | **5.48** / 6 (ideal) |
| Biseparable bound | 5 (W must exceed) |
| Significance | **12.3σ** |
| W / W_ideal | **91.3%** |
| Circuit depth | 3 (constant in n) |
| Measurement settings | 2 (constant in n) |
| Hardware job time | 7.5 s |

This is the smallest demonstration. The same circuit and witness scale to 50+ qubits with **no algorithmic change** — only one job per size — making this approach a direct path to 1000-qubit GME on next-generation hardware.

---

## Table of Contents

1. [Strategy: why 2D cluster + GME witness](#strategy-why-2d-cluster--gme-witness)
2. [Hardware](#hardware)
3. [Theory](#theory)
4. [Procedure walkthrough](#procedure-walkthrough)
5. [Supporting layers](#supporting-layers)
6. [Project structure](#project-structure)
7. [Running it yourself](#running-it-yourself)
8. [References](#references)

---

## Strategy: why 2D cluster + GME witness

### The challenge

Three goals:
1. **Prove** entanglement rigorously (no classical loophole)
2. Demonstrate **variety** of states/witnesses
3. Maximise **qubit count**

Scoring weights chosen: **20% qubits / 20% variety / 30% implementation / 20% theory**.

### What the prior art looks like

- **IQM Garnet 20-qubit GHZ (2024)**: F > 0.5, published. The "obvious" approach is already done.
- **IBM Eagle 127-qubit GHZ (2023)**: Linear chain, requires O(n) circuit depth.
- **MIT iQuHACK 2026 winner (Topological Ducks)**: 16 qubits with CHSH + graph states + qubit routing optimization (OR-Tools). Won on routing, not witness theory.

### Our angle

IQM Emerald is a **square lattice**. We use the *native* entanglement structure of that hardware: 2D cluster states.

| Property | GHZ chain | 2D cluster |
|---|---|---|
| Circuit depth | O(n) | **3** (constant) |
| Native gates needed | CNOT chain | H + CZ (both native on IQM) |
| Topology required | Linear path | Square grid (matches Emerald exactly) |
| Measurement settings for GME | O(n) | **2** (constant) |
| Robustness | Loss of one qubit kills it | Local correlations survive failures |

A 6-qubit cluster state and a 54-qubit cluster state run **the same circuit** with the same depth. This is what unlocks scalability to 1000+ qubits on future hardware.

---

## Hardware

[IQM Resonance](https://resonance.meetiqm.com) — IQM Emerald (54 qubits, square lattice).

| Spec | Value |
|---|---|
| Native gates | PRX (Phased-X), CZ |
| 1Q fidelity | 99.93% |
| 2Q fidelity | 99.5% |
| T₁ | 0.96 ms |
| Topology | Square lattice |
| Max shots/job | 20,000 |
| Max circuits/batch | 200 |

CZ along grid edges is native — our cluster circuit transpiles with **zero overhead**.

---

## Theory

### The cluster state

For a graph G = (V, E):

$$|C_G\rangle = \prod_{(i,j) \in E} \mathrm{CZ}_{ij} \cdot H^{\otimes n} |0\rangle^{\otimes n}$$

For a square-lattice graph, every edge corresponds to a native CZ on IQM hardware. The circuit has 3 layers regardless of n:

```
Layer 1:  H on every qubit              (n parallel H gates)
Layer 2:  CZ on every horizontal edge   (parallel)
Layer 3:  CZ on every vertical edge     (parallel)
```

### The stabilizers

The cluster state is the unique +1 eigenstate of n stabilizer operators:

$$g_i = X_i \otimes \bigotimes_{j \in N(i)} Z_j$$

where N(i) are i's graph neighbors. By construction, $\langle g_i \rangle = +1$ for all i in the perfect cluster state.

### The witness

$$\boxed{W = \sum_{i=1}^n \langle g_i \rangle}$$

**Theorem (Tóth & Gühne 2005)**: For any biseparable state, W ≤ n − 1. Therefore:

$$W > n - 1 \;\;\Longrightarrow\;\; \text{state is Genuinely Multipartitely Entangled}$$

A biseparable state is one that splits as $|\psi_A\rangle \otimes |\psi_B\rangle$ across some bipartition. GME means **no such splitting exists** — all n qubits are entangled together as one indivisible whole.

### Why only 2 measurement settings — the 2-coloring trick

The square lattice is bipartite (2-colorable like a chessboard):

```
B W B W B W
W B W B W B
B W B W B W
```

Every stabilizer $g_i$ measures one qubit in X (qubit i itself) and its neighbors in Z. Because every neighbor is the *opposite* color, two settings cover all stabilizers:

- **Setting A**: black qubits in X basis, white qubits in Z basis → reads ⟨g_i⟩ for every black qubit
- **Setting B**: white qubits in X basis, black qubits in Z basis → reads ⟨g_i⟩ for every white qubit

That's the entire witness. Two circuits, regardless of n. This is the central insight from the Tóth-Gühne 2005 papers and Zhou et al. 2019.

### Statistical significance

$$\sigma = \frac{W - (n-1)}{\sqrt{n}/\sqrt{N_{\text{shots}}}}$$

assuming variance ≤ 1 per stabilizer (true for ±1 eigenvalues). Standard threshold for "certified": σ > 3.

---

## Procedure walkthrough

The full step-by-step procedure with code, circuit diagrams, and hardware execution is in [`experiments/00_gme_walkthrough.ipynb`](experiments/00_gme_walkthrough.ipynb).

Summary of the steps (for a 2×3 = 6-qubit grid):

1. Connect to IQM Emerald via API token
2. Build the depth-3 cluster state circuit
3. Apply the checkerboard coloring → 2 measurement circuits
4. Submit one job (2 circuits × 4000 shots)
5. From bitstrings, compute each ⟨g_i⟩ = mean of $(-1)^{b_i + \sum_{j∈N(i)} b_j}$
6. Sum: W = Σ⟨g_i⟩
7. Check W > n−1

Per-qubit results from our hardware run:

| Qubit | ⟨g_i⟩ |
|---|---|
| 0 | +0.892 |
| 1 | +0.901 |
| 2 | +0.922 |
| 3 | +0.924 |
| 4 | +0.895 |
| 5 | +0.944 |
| **Sum (W)** | **5.477** |

All stabilizers within 0.89–0.94 — clean, consistent, no broken qubit.

---

## Supporting variety layer

While the cluster GME witness is our main result, **CHSH + Mermin** gives a separate mathematical framework (Bell inequalities) for the variety score:

- **CHSH** on a singlet: |S| > 2 proves 2-qubit entanglement (Tsirelson bound 2√2)
- **Mermin-n** on GHZ-n: |M_n| > 2^{n/2} classical bound; quantum value 2^{n−1} grows exponentially

Code: [`src/witnesses/chsh.py`](src/witnesses/chsh.py), [`src/witnesses/mermin.py`](src/witnesses/mermin.py).
Notebook: [`experiments/03_chsh_mermin.ipynb`](experiments/03_chsh_mermin.ipynb).

---

## Project structure

```
iqm_hackathon/
├── README.md                          ← this file
├── CLAUDE.md                          ← assistant instructions
├── hardware_test_gme.py               ← minimal hardware smoke test
├── test_smoke.py                      ← all modules on Aer simulator
│
├── src/
│   ├── backend.py                     ← IQM connection + sub-grid selector
│   ├── circuits/
│   │   ├── cluster_2d.py              ← 2D cluster state (depth 3)
│   │   └── ghz.py                     ← GHZ (variety layer)
│   └── witnesses/
│       ├── gme_cluster.py             ← GME witness (main)
│       ├── chsh.py                    ← CHSH on singlet (variety)
│       └── mermin.py                  ← Mermin M_n (variety)
│
└── experiments/
    ├── 00_gme_walkthrough.ipynb       ← step-by-step 2×3 walkthrough (start here)
    ├── 02_gme_scaling.ipynb           ← scaling sweep on hardware
    └── 03_chsh_mermin.ipynb           ← CHSH + Mermin variety layer
```

---

## Running it yourself

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows; on Linux/Mac: source .venv/bin/activate
pip install "iqm-client[qiskit]" numpy matplotlib scipy jupyter
```

**Without hardware** (Aer simulator):
```bash
python test_smoke.py
```

**With hardware** (IQM Resonance):
```bash
$env:IQM_TOKEN = "your_token_here"     # Windows PowerShell
python hardware_test_gme.py            # smallest demo
jupyter lab experiments/00_gme_walkthrough.ipynb
```

The smoke test on the simulator gives W = n exactly. The hardware test gives W slightly less than n, with σ depending on noise.

---

## References

**The witness**

1. G. Tóth, O. Gühne, *"Detecting Genuine Multipartite Entanglement with Two Local Measurements"*, **Phys. Rev. Lett. 94, 060501 (2005)**. [arXiv:quant-ph/0405165](https://arxiv.org/abs/quant-ph/0405165) — foundational 2-setting witness for cluster/GHZ states.

2. G. Tóth, O. Gühne, *"Entanglement detection in the stabilizer formalism"*, **Phys. Rev. A 72, 022340 (2005)**. [arXiv:quant-ph/0501020](https://arxiv.org/abs/quant-ph/0501020) — the W = Σ⟨g_i⟩ ≤ n−1 biseparable bound.

3. Y. Zhou, Q. Zhao, X. Yuan, X. Ma, *"Detecting multipartite entanglement structure with minimal resources"*, **npj Quantum Information 5, 83 (2019)**. [arXiv:1904.05001](https://arxiv.org/abs/1904.05001) — explicit 2-setting result for cluster states via 2-colorability.

4. N. K. H. Li, X. Dai, M. H. Muñoz-Arias, K. Reuer, M. Huber, N. Friis, *"Detecting genuine multipartite entanglement in multi-qubit devices with restricted measurements"*, **arXiv:2504.21076 (2025)** — recent extension to k-inseparability and SDP-optimized witnesses for restricted-measurement devices.

**Cluster states**

5. H. J. Briegel, R. Raussendorf, *"Persistent entanglement in arrays of interacting particles"*, **Phys. Rev. Lett. 86, 910 (2001)** — defining paper for cluster states.

**Bell inequalities (Layer 1)**

6. J. F. Clauser, M. A. Horne, A. Shimony, R. A. Holt, *"Proposed experiment to test local hidden-variable theories"*, **Phys. Rev. Lett. 23, 880 (1969)**.

7. N. D. Mermin, *"Extreme quantum entanglement in a superposition of macroscopically distinct states"*, **Phys. Rev. Lett. 65, 1838 (1990)**.

---

*ETH Quantum Hackathon 2026 — Team submission, branch `david`.*
