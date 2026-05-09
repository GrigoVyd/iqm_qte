# Witness My Entanglement — ETH Quantum Hackathon 2026

**IQM Challenge: Prove quantum entanglement on IQM hardware in the most compelling, scalable, and flexible way possible.**

## Headline results

We deploy **two qualitatively distinct entanglement witnesses on three different families of multipartite states**, all on real IQM hardware, and reach **GME on 20 qubits with the entire Garnet chip**.

### 1. Graph-state GME witness (main result)

Optimal spanning-tree graph state + parity-QREM + zero-noise extrapolation. Garnet, full chip:

| n | bound | W +QREM | σ +QREM | W +QREM+ZNE | σ +QREM+ZNE | GME |
|---|---|---|---|---|---|---|
| 6 | 5 | 5.99 | +15.6 | 6.13 | +23.2 | ✓ |
| 8 | 7 | 8.04 | +14.2 | 8.32 | +20.8 | ✓ |
| 12 | 11 | 11.70 | +7.8 | 12.25 | +15.1 | ✓ |
| 16 | 15 | 15.53 | +5.2 | 16.40 | +13.9 | ✓ |
| **20** | **19** | **19.18** | **+1.6** | **20.38** | **+13.1** | **✓** |

Theory: Tóth & Gühne 2005 PRL+PRA. Witness = $\sum_i \langle g_i\rangle \le n-1$ for any biseparable state. We generalised from rectangular clusters to **arbitrary 2-colorable graph states**, picking minimum-weight spanning trees on the live device topology — a tree on $n$ qubits has only $n-1$ CZ gates (vs $\sim 2n$ for a grid), so prep fidelity stays high.

### 2. GHZ vs W state fidelity comparison (variety)

Same Diker / F-gate W-state preparation, different N. Direct hardware fidelity comparison on Emerald & Garnet:

| n | F (W state) | F (GHZ state) | ΔF (W − GHZ) |
|---|---|---|---|
| 3 | 0.94 | 0.96 | −0.02 |
| 5 | 0.89 | 0.88 | +0.02 |
| 7 | 0.50 | 0.38 | +0.13 |
| 10 | 0.35 | 0.24 | +0.11 |
| 15 | 0.25 | 0.14 | +0.12 |
| 19 | 0.19 | 0.09 | **+0.10** |

**W states maintain higher fidelity than GHZ states as n grows** — empirical demonstration of GHZ fragility (single-qubit loss collapses entanglement) vs W robustness (loss leaves residual entanglement).

### 3. Direct Fidelity Estimation (DFE)

Flammia & Liu 2011 — sample random stabilizers, average. **Quantitative state fidelity, not just witness violation.** Garnet at every n:

| n | F (DFE) |
|---|---|
| 6 | 0.832 ± 0.007 |
| 8 | 0.712 ± 0.010 |
| 12 | 0.626 ± 0.012 |
| 16 | 0.502 ± 0.010 |
| 20 | 0.334 ± 0.013 |

### 4. W-state non-linear entanglement witness

Z-basis statistics alone cannot distinguish a W state from a classical mixture of single-excitation strings. The X-basis pairwise correlator $\langle X_i X_j\rangle$ does: quantum gives $2/N$, classical gives $0$. Detected on hardware via the F-gate Diker preparation.

---

## Table of Contents

1. [Strategy and why we win](#strategy-and-why-we-win)
2. [Hardware](#hardware)
3. [Theory I — Graph-state GME witness](#theory-i--graph-state-gme-witness)
4. [Theory II — W states + non-linear witness](#theory-ii--w-states--non-linear-witness)
5. [Pipeline: from device metrics to certified W](#pipeline-from-device-metrics-to-certified-w)
6. [Variety table](#variety-table)
7. [Project structure](#project-structure)
8. [Running it yourself](#running-it-yourself)
9. [References](#references)

---

## Strategy and why we win

### The challenge

Three goals:
1. **Prove** entanglement rigorously (no classical loophole)
2. Demonstrate **variety** of states/witnesses
3. Maximise **qubit count**

Scoring weights chosen: **20% qubits / 20% variety / 30% implementation / 20% theory**.

### What the prior art looks like

- **IQM Garnet 20-qubit GHZ (2024)**: F > 0.5, published. Standard GHZ approach.
- **IBM Eagle 127-qubit GHZ (2023)**: Linear chain, requires O(n) circuit depth.
- **MIT iQuHACK 2026 winner (Topological Ducks)**: 16 qubits with CHSH + graph states + OR-Tools routing optimization. Won on routing, not witness theory.

### Our approach in five prongs

| Prong | What it does | Effect |
|---|---|---|
| **2-colorable graph states** | Use any connected bipartite subgraph instead of rectangular grids | n−1 CZ gates (tree) vs ~2n (grid) — 30% fewer errors |
| **Threshold filter + multi-start Prim's** | Drop weak qubits, route around bad CZ pairs via min-weight spanning tree on live calibration data | Truly optimal sub-tree per n |
| **Parity QREM** | Per-qubit readout correction from device-reported error rates | +0.10 lift to ⟨g_i⟩ — no extra calibration circuits |
| **Zero-noise extrapolation** | CZ folding (with barriers to prevent transpiler cancellation) → linear fit → extrapolate to α=0, with bootstrap σ | Pushes n=20 from σ=+1.6 to σ=+13 on Garnet |
| **Variety: W vs GHZ vs graph states** | Three qualitatively different state families, two different witnesses (sum-of-stabilizers + non-linear pairwise correlator) | Honest 4 distinct hardware demonstrations |

Each prong is independently sound and stacks with the others.

---

## Hardware

[IQM Resonance](https://resonance.meetiqm.com) — both Emerald (54q) and Garnet (20q).

| Spec | Emerald | Garnet |
|---|---|---|
| Qubits | 54 | 20 |
| Native gates | PRX (Phased-X), CZ | PRX, CZ |
| 1Q fidelity | ~99.93% | ~99.93% |
| 2Q fidelity (CZ) | 92–99% (varies per pair) | 95–99% |
| T₁ | up to 0.96 ms (varies) | up to 100 µs |
| Topology | Bipartite — 81 native CZ pairs | Bipartite — 30 native CZ pairs |
| Qubits passing our default thresholds | 50 / 54 | **20 / 20** |

**Key empirical finding**: Garnet has uniformly cleaner calibration. Every Garnet qubit passes our default thresholds; on Emerald 4 are excluded. Garnet is the better demonstration platform for our witness.

---

## Theory I — Graph-state GME witness

### The graph state

For any graph G = (V, E):

$$|G\rangle = \prod_{(i,j) \in E} \mathrm{CZ}_{ij} \cdot H^{\otimes n} |0\rangle^{\otimes n}$$

H on every qubit, then CZ on every edge. Constant-time per layer.

### The stabilizer

For each qubit i:

$$g_i = X_i \otimes \bigotimes_{j \in N(i)} Z_j$$

By construction, $\langle g_i \rangle = +1$ for the perfect graph state.

### The witness (Tóth & Gühne 2005)

$$W = \sum_{i=1}^n \langle g_i \rangle$$

**Theorem**: For any biseparable state, $W \leq n - 1$. Therefore:

$$W > n - 1 \implies \text{Genuine Multipartite Entanglement}$$

This bound holds for **any connected 2-colorable graph state**.

### Why a spanning tree wins

For a connected graph on n nodes:
- Any spanning tree has exactly n − 1 edges (= n − 1 CZ gates) — minimum possible
- A 2D rectangular grid has ~2n edges
- Trees are always 2-colorable

So a tree gives **fewer CZ gates**, hence higher prep fidelity, hence higher $\langle g_i\rangle$, hence W can clear n−1 at larger sizes. We empirically verified this in `01_emerald_showcase.ipynb` — a complete bipartite graph on 8 qubits fails to certify GME under noise that easily certifies a tree.

### Why only 2 measurement settings

Same trick as for cluster states: BFS-2-color the chosen graph (always possible — it's bipartite). Then:
- **Setting A**: 0-colored qubits in X basis, 1-colored in Z basis → reads ⟨g_i⟩ for every 0-colored qubit
- **Setting B**: swap → reads ⟨g_i⟩ for every 1-colored qubit

Two circuits, regardless of n.

---

## Theory II — W states + non-linear witness

### W state preparation (Diker / F-gate method)

$$|W_N\rangle = \frac{1}{\sqrt{N}}\sum_{k=0}^{N-1}|0\cdots 1_k \cdots 0\rangle$$

A single excitation distributed coherently across all $N$ qubits. The F-gate cascade prepares it in $O(N)$ two-qubit gates:

```
|10...0⟩ → F-gate chain → CNOT correction → |W_N⟩
```

Where each F-gate $F_k = R_y(-\theta_k) \cdot CZ \cdot R_y(\theta_k)$ with $\theta_k = \arccos\sqrt{1/(N-k+1)}$.

### Non-linear entanglement witness

The Z-basis statistics — single qubit excited with probability $1/N$, joint excitation forbidden — are **also reproduced by a classical mixture** $\rho_{\text{classical}} = \frac{1}{N}\sum_k |0\cdots 1_k \cdots 0\rangle\langle\cdots|$.

The X-basis breaks the tie:

| State | $\langle X_i X_j\rangle$ for $i \neq j$ |
|---|---|
| Quantum $\|W_N\rangle$ | $2/N$ |
| Classical mixture | $0$ |

A measured $\langle X_i X_j\rangle > 0$ at confidence above shot noise witnesses **genuine quantum coherence** — a different framework than the stabilizer-sum witness.

---

## Pipeline: from device metrics to certified W

```
┌─────────────────────────────────────────────────────────┐
│ 1. Pull live calibration from IQM API                   │
│    - Per-qubit: T1, T2, readout fidelity, 1Q gate fid   │
│    - Per-CZ-pair: gate fidelity                         │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 2. Threshold filter qubits                              │
│    Drop if RO < 0.90, 1Q < 0.99, T1 < 10µs, T2 < 5µs    │
│    Or if best CZ partner has fidelity < 0.90            │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 3. Min-weight spanning tree (Prim's) + local search     │
│    Edge weight = (1 − F_CZ) + T1/T2/1Q penalties         │
│    Try every survivor as seed, keep lightest tree       │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 4. Predict ⟨g_i⟩ for each i a-priori                    │
│    (CZ × T1/T2 × 1Q channel multiplication)             │
│    If predicted W > n-1: feasible. Else: warn user.     │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 5. Build folded circuits at scales α = 1, 3, 5          │
│    Replace each CZ with α copies separated by barriers  │
│    Same logical effect (CZ² = I), α× the gate noise     │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 6. Submit ALL circuits in one batched job               │
│    Drift-fair comparison across n × scales              │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 7. Apply parity-QREM per stabilizer                     │
│    c_i = 1/(P(0|0)_i + P(1|1)_i − 1)                    │
│    Multiply ⟨g_i⟩ by ∏ c_j over involved qubits          │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 8. ZNE extrapolation with shot-level bootstrap σ        │
│    Fit W(α=1,3,5) linearly, extrapolate to α=0          │
│    Resample 200 times for proper variance propagation    │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ 9. Verdict: W > n−1 ⟹ GME at significance σ              │
│    Plus DFE for direct fidelity estimate                 │
└─────────────────────────────────────────────────────────┘
```

The entire pipeline runs end-to-end in [`experiments/01_emerald_showcase.ipynb`](experiments/01_emerald_showcase.ipynb) and [`experiments/02_garnet_showcase.ipynb`](experiments/02_garnet_showcase.ipynb).

---

## Variety table

We deploy **3 distinct multipartite state families** with **2 distinct entanglement witnesses** on real hardware:

| State | Preparation | Witness type | n range | Notebook |
|---|---|---|---|---|
| Graph state (spanning tree) | H + CZ on tree edges | Stabilizer sum (Tóth-Gühne) | 6 – 20 | `01`, `02` |
| GHZ state | H + CNOT ladder | Probability fidelity | 3 – 20 | `03` |
| W state | F-gate (Diker) cascade | Probability fidelity + non-linear $\langle X_i X_j\rangle$ | 3 – 30 | `03`, `04` |

Plus Direct Fidelity Estimation (Flammia-Liu 2011) for state quality, complementing all three.

---

## Project structure

```
iqm_hackathon/
├── README.md                       ← this file
├── CLAUDE.md                       ← assistant instructions
├── test_smoke.py                   ← Aer smoke test
│
├── src/
│   ├── backend.py                  ← IQM connection, threshold filter,
│   │                                  cost function, predictor, sub-tree finder
│   ├── circuits/
│   │   ├── cluster_2d.py           ← rectangular 2D cluster state
│   │   └── graph_state.py          ← arbitrary graph state + 2-coloring
│   ├── witnesses/
│   │   ├── gme_cluster.py          ← rectangular cluster GME witness
│   │   └── gme_graph.py            ← MAIN: arbitrary graph GME witness +
│   │                                  Tóth-Gühne fidelity lower bound
│   ├── mitigation/
│   │   ├── parity_qrem.py          ← readout correction
│   │   └── zne.py                  ← CZ folding + bootstrap extrapolation
│   ├── dfe.py                      ← Direct Fidelity Estimation
│   └── visualization.py            ← exact IQM dashboard layouts (Emerald + Garnet),
│                                     spotlight + Prim's animation
│
└── experiments/
    ├── 00_gme_walkthrough.ipynb         ← step-by-step 2×3 cluster intro
    ├── 01_emerald_showcase.ipynb        ← MAIN: comprehensive Emerald run
    │                                       (raw / +QREM / +ZNE / +QREM+ZNE + DFE)
    ├── 02_garnet_showcase.ipynb         ← MAIN: same on Garnet (cleaner chip)
    ├── 03_ghz_vs_w_comparison.ipynb     ← VARIETY: GHZ vs W fidelity scaling
    ├── 04_w_state_entanglement.ipynb    ← VARIETY: W-state non-linear witness
    └── run_w_state.py                   ← W-state runner script
```

---

## Running it yourself

```bash
python -m venv .venv
.venv\Scripts\activate              # Windows
pip install "iqm-client[qiskit]" iqm-qubit-selector \
            numpy matplotlib scipy jupyter rustworkx networkx
```

**Without hardware** (Aer simulator):
```bash
python test_smoke.py
```

**With hardware** (IQM Resonance):
```bash
$env:IQM_TOKEN = "your_token_here"
jupyter lab experiments/02_garnet_showcase.ipynb
```

Hardware budget: full Emerald + Garnet sweeps (witnesses + ZNE + DFE) ≈ 50 IQM tokens.

---

## References

**Graph-state GME witness**

1. G. Tóth, O. Gühne, *"Detecting Genuine Multipartite Entanglement with Two Local Measurements"*, **Phys. Rev. Lett. 94, 060501 (2005)**. [arXiv:quant-ph/0405165](https://arxiv.org/abs/quant-ph/0405165)
2. G. Tóth, O. Gühne, *"Entanglement detection in the stabilizer formalism"*, **Phys. Rev. A 72, 022340 (2005)**. [arXiv:quant-ph/0501020](https://arxiv.org/abs/quant-ph/0501020)
3. Y. Zhou et al., *"Detecting multipartite entanglement structure with minimal resources"*, **npj Quantum Information 5, 83 (2019)**. [arXiv:1904.05001](https://arxiv.org/abs/1904.05001)

**W-state preparation**

4. F. Diker, *"Deterministic construction of arbitrary W states with quadratically increasing number of two-qubit gates"*, **arXiv:1606.09290 (2016)**.

**Mitigation**

5. K. Temme, S. Bravyi, J. M. Gambetta, *"Error mitigation for short-depth quantum circuits"*, **PRL 119, 180509 (2017)** — ZNE.
6. S. Bravyi et al., *"Mitigating measurement errors in multiqubit experiments"*, **Phys. Rev. A 103, 042605 (2021)** — QREM.
7. S. Flammia, Y.-K. Liu, *"Direct Fidelity Estimation from Few Pauli Measurements"*, **PRL 106, 230501 (2011)** — DFE.

**Graph states**

8. H. J. Briegel, R. Raussendorf, *"Persistent entanglement in arrays of interacting particles"*, **PRL 86, 910 (2001)**.

---

*ETH Quantum Hackathon 2026 — Team submission, branch `submission`. 20-qubit GME on Garnet certified.*
