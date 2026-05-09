"""IQM Resonance backend connection and topology-aware qubit selection."""

from __future__ import annotations

import os

import rustworkx as rx
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

try:
    from iqm.qiskit_iqm import IQMProvider
    IQM_AVAILABLE = True
except ImportError:
    IQM_AVAILABLE = False

RESONANCE_URL = "https://resonance.meetiqm.com"


def get_backend(token: str | None = None, device: str = "emerald"):
    """Connect to IQM Resonance and return a Qiskit backend.

    If `token` is None and IQM_TOKEN env var is unset, falls back to AerSimulator.
    The IQM client errors if both env var and arg are set, so we pass arg only
    when the env var is absent.
    """
    env_token = os.environ.get("IQM_TOKEN")
    effective = token or env_token
    if effective and IQM_AVAILABLE:
        if env_token:
            provider = IQMProvider(RESONANCE_URL, quantum_computer=device)
        else:
            provider = IQMProvider(RESONANCE_URL, quantum_computer=device, token=token)
        return provider.get_backend()
    print("No IQM token — using Aer simulator.")
    return AerSimulator()


def is_simulator(backend) -> bool:
    return isinstance(backend, AerSimulator)


def _enumerate_subgrids(backend, rows: int, cols: int,
                        max_grids: int = 500) -> list[list[int]]:
    """Find all contiguous rows×cols rectangular sub-grids on the device topology.

    Uses VF2 subgraph isomorphism. Returns layouts in row-major order matching
    the cluster_2d circuit's logical qubit numbering. Up to `max_grids` results.
    """
    target = rx.PyGraph()
    target.add_nodes_from(range(rows * cols))
    for r in range(rows):
        for c in range(cols - 1):
            target.add_edge(r * cols + c, r * cols + c + 1, None)
    for r in range(rows - 1):
        for c in range(cols):
            target.add_edge(r * cols + c, (r + 1) * cols + c, None)

    device = rx.PyGraph()
    device.add_nodes_from(range(backend.num_qubits))
    seen_edges: set[tuple[int, int]] = set()
    for a, b in backend.coupling_map:
        e = (min(a, b), max(a, b))
        if e in seen_edges:
            continue
        seen_edges.add(e)
        device.add_edge(a, b, None)

    grids: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for m in rx.vf2_mapping(device, target, subgraph=True, induced=False):
        inv = {tgt: dev for dev, tgt in m.items()}
        layout = [inv[i] for i in range(rows * cols)]
        key = tuple(layout)
        if key in seen:
            continue
        seen.add(key)
        grids.append(layout)
        if len(grids) >= max_grids:
            break
    return grids


def select_best_subgrid(backend, rows: int, cols: int,
                         readout_mode: str = "fidelity") -> tuple[list[int], float, int]:
    """Pick the highest-fidelity contiguous rows×cols sub-grid.

    1. Enumerate every contiguous rectangular sub-grid on the device.
    2. Score each via IQM's CostEvaluator (live CZ + readout calibration).
    3. Return (layout, cost, n_candidates).

    Layout matches the cluster_2d circuit's row-major qubit ordering, so passing
    `initial_layout=layout` to transpile gives zero SWAP overhead.
    """
    if is_simulator(backend):
        return list(range(rows * cols)), 0.0, 1

    candidates = _enumerate_subgrids(backend, rows, cols)
    if not candidates:
        raise ValueError(f"No {rows}×{cols} contiguous sub-grid found on device.")

    try:
        from iqm.qubit_selector.qubit_selector import (
            CostEvaluator, CostFunction, ReadoutMode,
        )
    except ImportError:
        return candidates[0], 0.0, len(candidates)

    rm = {"none": ReadoutMode.NONE, "fidelity": ReadoutMode.FIDELITY,
          "qndness": ReadoutMode.QNDNESS}.get(readout_mode, ReadoutMode.FIDELITY)

    n = rows * cols
    canonical = QuantumCircuit(n)
    canonical.h(range(n))
    for r in range(rows):
        for c in range(cols - 1):
            canonical.cz(r * cols + c, r * cols + c + 1)
    for r in range(rows - 1):
        for c in range(cols):
            canonical.cz(r * cols + c, (r + 1) * cols + c)

    evaluator = CostEvaluator(
        backend=backend, quantum_circuit=canonical,
        cost_function=CostFunction.GATE_COST_CZ,
        readoutmode=rm, layouts=candidates,
    )
    top_layouts, scores = evaluator.get_top_layouts(num_layouts=1)
    return top_layouts[0], scores[0], len(candidates)
