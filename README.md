# Witness My Entanglement — ETH Quantum Hackathon 2026

**Q.te team · IQM Challenge** — prove quantum entanglement on IQM hardware (Emerald 54 q + Garnet 20 q) in the most compelling, scalable, and flexible way possible.

## ⭐ Start here

The submission is a single self-contained notebook:

> **[experiments/witness_my_entanglement.ipynb](experiments/witness_my_entanglement.ipynb)**
>
> Three threads — W states, graph-state GME, routed Bell pair — plus
> a routing-coefficient outlook (gri's branch). Every figure and table
> in this README is produced from this notebook. A single
> `RERUN_HW = False` flag at the top reproduces every plot offline
> from the saved JSONs in `experiments/consolidated_results/`; flip
> to `True` to re-collect data on hardware (~12 batched jobs).

There is also a **scrolling presentation site** at
[`presentation/index.html`](presentation/index.html) — open in any
browser, walk through with the arrow keys, press `F` for fullscreen,
click any plot to zoom.

## Headline results

We deploy **four qualitatively distinct entanglement witnesses on three different families of multipartite states**, all on real IQM hardware, and reach **GME on 20 qubits with the entire Garnet chip**.

### 1. Graph-state GME witness (main result)

Optimal spanning-tree graph state + parity-QREM + zero-noise extrapolation. Both devices, sweeps from `experiments/witness_my_entanglement.ipynb`. The physical maximum is $W=n$ (perfect state); the GME bound is $W=n-1$.

**Garnet (full 20-qubit chip):**

| n | bound (n−1) | W +QREM | σ +QREM | W +QREM+ZNE (raw / clipped) | σ +QREM+ZNE | GME |
|---|---|---|---|---|---|---|
| 6 | 5 | 6.04 | +19.1 | 6.17 / **6.00** | +28.2 | ✓ |
| 12 | 11 | 11.98 | +12.6 | 12.52 / **12.00** | +24.5 | ✓ |
| **20** | **19** | **19.25** | **+2.5** | 20.45 / **20.00** | **+15.0** | **✓** |

**Emerald (54q):**

| n | bound (n−1) | W +QREM | σ +QREM | W +QREM+ZNE (raw / clipped) | σ +QREM+ZNE | GME |
|---|---|---|---|---|---|---|
| 6 | 5 | 6.06 | +19.4 | 6.34 / **6.00** | +35.7 | ✓ |
| 12 | 11 | 11.88 | +11.3 | 12.84 / **12.00** | +28.2 | ✓ |
| 20 | 19 | 18.11 | −8.9 | 19.75 / **19.75** | **+7.6** | ✓ (after ZNE) |

The certified result is **+QREM** (physical at every $n$, above the GME bound on Garnet at every n including the full chip). On Emerald n=20 the calibration-cost tree picks a region where +QREM alone is below the bound; ZNE then lifts the witness above it at +7.6σ. +QREM+ZNE values that exceed the physical ceiling $W=n$ are clipped — the 1–2% overshoot is the linear-fit systematic. Empirical-edge-fidelity routing (§2.4) closes most of the n=20 Emerald gap, see below.

Theory: Tóth & Gühne 2005 PRL+PRA. Witness = $\sum_i \langle g_i\rangle \le n-1$ for any biseparable state. We generalised from rectangular clusters to **arbitrary 2-colorable graph states**, picking minimum-weight spanning trees on the live device topology — a tree on $n$ qubits has only $n-1$ CZ gates (vs $\sim 2n$ for a grid), so prep fidelity stays high.

### 2. GHZ vs W state fidelity comparison (legacy variety result)

Earlier iteration compared the F-gate W preparation against an H+CNOT-ladder GHZ at matched n on both devices. W stayed higher than GHZ at every n past ~5 qubits (e.g. n=19 on Emerald: F_W=0.19 vs F_GHZ=0.09) — empirical demonstration of GHZ fragility (single-qubit loss collapses entanglement) vs W robustness. Kept in [legacy/03_ghz_vs_w_comparison.ipynb](experiments/legacy/03_ghz_vs_w_comparison.ipynb); not part of the consolidated submission run.

### 3. Free fidelity lower bound (Tóth-Gühne)

Same two-circuit data that gives W also gives a **free** lower bound on the state fidelity to the target graph state: $F \geq \langle P_A\rangle + \langle P_B\rangle - 1$ where $\langle P_X\rangle$ is the projection onto the +1 eigenspace of all stabilizers in setting $X$. From the consolidated run:

| Device | n | ⟨P_A⟩ | ⟨P_B⟩ | F_lb |
|--------|---|---|---|---|
| Garnet  | 6  | 0.872 | 0.903 | **0.776** |
| Garnet  | 12 | 0.733 | 0.753 | **0.486** |
| Garnet  | 20 | 0.511 | 0.497 | **0.008** |
| Emerald | 6  | 0.871 | 0.904 | **0.775** |
| Emerald | 12 | 0.702 | 0.795 | **0.497** |
| Emerald | 20 | 0.358 | 0.434 | −0.208 (improves to **+0.11** with empirical-edge tree, §6) |

DFE (Flammia-Liu 2011) was used in earlier iterations (legacy notebooks 01/02) for an independent quantitative fidelity check; it agrees with the bound above within shot noise. Implementation kept in [src/dfe.py](src/dfe.py).

### 4. W-state non-linear entanglement witness

Z-basis statistics alone cannot distinguish a W state from a classical mixture of single-excitation strings. The X-basis pairwise correlator $\langle X_i X_j\rangle$ does: quantum gives $2/N$, classical gives $0$. The full table per device is in §5; the headline is that on **both** Emerald and Garnet at n ∈ {5, 10, 15, 19} the measured witness is +25σ to +60σ above the classical bound 0.

### 5. W-state scaling — both devices

Same F-gate W preparation across n ∈ {5, 10, 15, 19} on **both** Emerald and Garnet, with beam-search Hamiltonian-path routing. n=19 is the maximum chain length on Garnet (the 20-qubit Apollo lattice has no Hamiltonian path of length 20). Conservative σ = std across pairs / √n_pairs:

| Device | n | Z fidelity F_z | ⟨W_x⟩ | quantum ideal 2/n | σ_avg | σ above classical 0 |
|--------|---|---|---|---|---|---|
| Garnet  | 5  | 0.913 | 0.384 | 0.400 | 0.0064 | **+59.8** |
| Garnet  | 10 | 0.800 | 0.170 | 0.200 | 0.0058 | **+29.1** |
| Garnet  | 15 | 0.630 | 0.114 | 0.133 | 0.0041 | **+27.9** |
| Garnet  | 19 | 0.449 | 0.082 | 0.105 | 0.0032 | **+25.6** |
| Emerald | 5  | 0.885 | 0.388 | 0.400 | 0.0088 | **+43.8** |
| Emerald | 10 | 0.748 | 0.174 | 0.200 | 0.0068 | **+25.7** |
| Emerald | 15 | 0.587 | 0.114 | 0.133 | 0.0033 | **+35.0** |
| Emerald | 19 | 0.345 | 0.093 | 0.105 | 0.0030 | **+30.7** |

The non-linear witness rules out the classical single-excitation mixture decisively (≥ +25σ) at every n on both chips, even where Z-basis fidelity to ideal-|W⟩ has dropped below 0.5. Z-fidelity tracks chain quality and decays roughly exponentially in n, as expected from the F-gate cascade depth.

### 6. Routed Bell pair via measurement-based teleportation (Garnet)

Pick a path A–q₁–…–B through the cluster, measure all internals in X
(and any 2D off-path neighbours in Z). What remains on (A, B) is locally
equivalent to a Bell pair, with a Pauli byproduct determined by the
internal-X outcomes. We correct the byproduct in postprocessing and verify
both Bell fidelity and the CHSH inequality at the endpoints.

| Path L | F (with CZ) | F (no-CZ control) | \|S\| (with) | \|S\| (no-CZ) | Bell 3σ | CHSH 3σ |
|--------|-------------|-------------------|------------|------------|---------|---------|
| 3 | **0.895** | 0.494 | **2.500** | 1.384 | ✓ | ✓ |
| 5 | **0.819** | 0.493 | **2.222** | 1.413 | ✓ | ✓ |
| 7 | **0.743** | 0.495 | 1.918 | 1.330 | ✓ | ✗ |

Bell entanglement is certified at every L; CHSH violates the classical
bound (\|S\| > 2) at L=3 and L=5 and lands just below at L=7. The no-CZ
control flatlines at the separable bound F=0.5 and the maximum-mixture
\|S\|=√2 — exactly the expected behaviour of a state with no cluster
entanglement to extract from. Implementation: [src/witnesses/routed_bell.py](src/witnesses/routed_bell.py).

### 7. Empirical-fidelity tree selection (Anna's edge map → Prim's)

Measure $F_{ij} = (1+\langle X_iZ_j\rangle+\langle Z_iX_j\rangle+\langle Y_iY_j\rangle)/4$ on every native CZ pair (greedy edge-coloring → matchings: 4 matchings cover all 30 Garnet edges and all 81 Emerald edges). Plug $w(i,j)=1-F_{ij}^{\text{measured}}$ into the same multi-start Prim's MST as the calibration tree. Head-to-head on hardware:

| Device | n | Method | W_raw | W +QREM | F_lower_bound |
|--------|---|--------|-------|---------|---------------|
| Garnet | 12 | calibration | 10.50 | 11.88 | 0.46 |
| Garnet | 12 | **empirical** | **10.59** | **12.07** | **0.50** |
| Emerald | 20 | calibration | 15.96 | 18.12 | **−0.20** |
| Emerald | 20 | **empirical** | **17.32** | **18.69** | **+0.11** |

The empirical tree wins on every metric on both devices. On Emerald n=20 it picks an entirely different qubit set (qubits 0–13/19/23/24/32/33/40, vs the calibration tree's 22-53 region), and the change **flips F_lb from negative to positive** — a real qualitative win, not just a numerical nudge. Saved at [empirical_vs_calibration.png](experiments/consolidated_results/empirical_vs_calibration.png).

### 8. Routing strategies head-to-head

Same n=15 W-state circuit, two routing approaches in one batched Garnet job for direct comparison (from the consolidated run):

| Metric | IQM Qubit Selector | Beam-search chain |
|---|---|---|
| Z-basis fidelity F_z | 0.313 | **0.679** |
| ⟨W_x⟩ (X-basis witness) | 0.085 | **0.114** |
| σ above classical 0 | +14.3 | **+31.2** |

Custom beam-search chain routing for linear-interaction circuits dramatically outperforms general-purpose layout selection on this specific circuit family. **The IQM Qubit Selector is excellent for arbitrary circuits, but trees / chains beat it when the circuit's interaction graph is itself a tree / chain.**

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
| GHZ state | H + CNOT ladder | Probability fidelity | 3 – 19 | `03` |
| W state | F-gate (Diker) cascade | Probability fidelity + non-linear $\langle X_i X_j\rangle$ | 3 – 30 | `03`, `04`, `05` |

Plus Direct Fidelity Estimation (Flammia-Liu 2011) for state quality, complementing all three.

**Implementation sophistication** is also showcased explicitly: notebook `06` runs the same W-state circuit through both the IQM Qubit Selector (CostEvaluator) and a custom beam-search chain router that exploits linear-interaction structure for swap-free routing. Direct head-to-head hardware comparison.

---

## Project structure

```
iqm_hackathon/
├── README.md                       ← this file
├── CLAUDE.md                       ← assistant instructions
├── test_smoke.py                   ← Aer smoke test
│
├── src/
│   ├── backend.py                  ← IQM connection, threshold filter, MST tree
│   │                                  selectors (calibration + empirical-F variants)
│   ├── circuits/
│   │   ├── cluster_2d.py           ← rectangular 2D cluster state
│   │   ├── graph_state.py          ← arbitrary graph state + 2-coloring
│   │   └── w_state.py              ← W-state F-gate (Diker) cascade
│   ├── routing/beam_chain.py       ← beam-search Hamiltonian-path router
│   ├── witnesses/
│   │   ├── gme_graph.py            ← Tóth-Gühne stabilizer-sum witness
│   │   ├── gme_cluster.py          ← rectangular cluster variant
│   │   ├── w_witness.py            ← Z-fid + non-linear pairwise X-witness
│   │   └── routed_bell.py          ← cluster-MBQC routed Bell pair + CHSH
│   ├── mitigation/
│   │   ├── parity_qrem.py          ← readout correction (no extra cal shots)
│   │   └── zne.py                  ← CZ folding (with barriers) + bootstrap
│   ├── diagnostics/edge_bell_map.py ← per-CZ-pair fidelity map ("Anna's" map)
│   ├── dfe.py                      ← Direct Fidelity Estimation
│   └── visualization.py            ← exact IQM dashboard layouts + Prim's animation
│
├── experiments/
│   ├── witness_my_entanglement.ipynb    ← ⭐ MAIN: full story end-to-end on both devices
│   │                                       Part 1 (W-states) + Part 2 (graph GME +
│   │                                       empirical-F tree) + Part 3 (routed Bell)
│   │                                       + §1.4 routing-coefficient optimization
│   ├── consolidated_results/            ← JSON + PNG outputs of the consolidated run
│   └── legacy/                          ← per-topic notebooks from the iterative phase
│
├── presentation/                   ← scroll-snapping deck (open index.html in browser)
└── qte_brand_kit/                  ← logo, palette, slide title backdrop
```

---

## Running it yourself

```bash
python -m venv .venv
.venv\Scripts\activate              # Windows
pip install "iqm-client[qiskit]" iqm-qubit-selector qiskit-aer \
            numpy matplotlib scipy jupyter rustworkx networkx
```

**Reproduce every plot offline (no hardware, no IQM token)** — open
[`experiments/witness_my_entanglement.ipynb`](experiments/witness_my_entanglement.ipynb)
in JupyterLab and run all. The first cell sets `RERUN_HW = False`
by default, so every section reads its saved JSON from
`experiments/consolidated_results/` and reproduces the figures
exactly. With **no token at all**, set `QTE_OFFLINE=1` to skip the
`.secrets/iqm_api_key` lookup and force the Aer-simulator shim
(Garnet / Emerald topology preserved):

```bash
QTE_OFFLINE=1 jupyter nbconvert --to notebook --execute \
    --inplace experiments/witness_my_entanglement.ipynb
```

**Quick sanity check** (Aer-only smoke test of every src/ helper):
```bash
python test_smoke.py
```

**With hardware** (IQM Resonance):
```bash
$env:IQM_TOKEN = "your_token_here"
jupyter lab experiments/witness_my_entanglement.ipynb
```

The consolidated notebook reads `RERUN_HW` at the top: set `True` to submit
fresh jobs (one batched job per device per experiment block, ~8 jobs total),
or `False` to reproduce every figure from the saved JSONs in
`experiments/consolidated_results/`.

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
