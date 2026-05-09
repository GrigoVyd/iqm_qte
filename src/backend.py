"""IQM Resonance backend setup and qubit selection utilities."""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator

try:
    from iqm.qiskit_iqm import IQMProvider
    IQM_AVAILABLE = True
except ImportError:
    IQM_AVAILABLE = False

RESONANCE_URL = "https://resonance.meetiqm.com"


def get_backend(token: str | None = None, device: str = "emerald"):
    """Connect to IQM Resonance and return a Qiskit backend.

    Falls back to Aer simulator if no token is provided.
    Token can also be set via IQM_TOKEN environment variable.
    """
    token = token or os.environ.get("IQM_TOKEN")
    if token and IQM_AVAILABLE:
        provider = IQMProvider(RESONANCE_URL, quantum_computer=device, token=token)
        return provider.get_backend()
    print("No token provided — using Aer simulator.")
    return AerSimulator()


def is_simulator(backend) -> bool:
    return isinstance(backend, AerSimulator)


# ---------------------------------------------------------------------------
# Qubit selection
# ---------------------------------------------------------------------------

@dataclass
class QubitLayout:
    """Selected qubit subset with its coupling map."""
    qubit_names: list[str]          # IQM names, e.g. ["QB5", "QB6", "QB10"]
    qubit_indices: list[int]        # Qiskit 0-based indices
    coupling_map: list[list[int]]   # edges within this subset (Qiskit indices)
    # fidelity scores for reference
    scores: dict[str, float] = field(default_factory=dict)


def get_calibration_scores(backend) -> dict[str, float]:
    """Return a per-qubit score (higher = better) from backend properties.

    Uses average gate error across all gates on that qubit.
    On a simulator, returns uniform scores.
    """
    if is_simulator(backend):
        n = backend.configuration().n_qubits if hasattr(backend, "configuration") else 30
        return {f"QB{i+1}": 1.0 for i in range(n)}

    props = backend.properties()
    scores: dict[str, float] = {}
    for q_idx, qubit_name in enumerate(backend.qubit_name_to_index.__self__.qubit_name_to_index if False else []):
        pass  # placeholder — filled below

    # Build index->name map from backend
    try:
        # iqm-client exposes qubit_name_to_index as a method
        name_to_idx: dict[str, int] = {
            name: backend.qubit_name_to_index(name)
            for name in backend.qubit_names
        }
        idx_to_name = {v: k for k, v in name_to_idx.items()}
    except AttributeError:
        # Fallback: name qubits QB1..QBn
        n = backend.num_qubits
        idx_to_name = {i: f"QB{i+1}" for i in range(n)}

    for idx, name in idx_to_name.items():
        try:
            # Average readout error + gate errors as a penalty
            readout_err = props.readout_error(idx)
            gate_errs = [props.gate_error(g.gate, [idx])
                         for g in props.gates if idx in g.qubits and len(g.qubits) == 1
                         if props.gate_error(g.gate, [idx]) is not None]
            avg_gate_err = float(np.mean(gate_errs)) if gate_errs else 0.01
            scores[name] = 1.0 - readout_err - avg_gate_err
        except Exception:
            scores[name] = 0.5
    return scores


def _idx_to_name_map(backend) -> dict[int, str]:
    try:
        return {backend.qubit_name_to_index(n): n for n in backend.qubit_names}
    except AttributeError:
        return {i: f"QB{i+1}" for i in range(backend.num_qubits)}


def _name_to_idx_map(backend) -> dict[str, int]:
    return {v: k for k, v in _idx_to_name_map(backend).items()}


def get_best_chain(backend, n: int, scores: dict[str, float] | None = None) -> QubitLayout:
    """Find a length-n path in the coupling graph maximising qubit fidelity.

    Uses a greedy best-first search: start from highest-scoring qubit,
    always extend to the highest-scoring unvisited neighbour.
    """
    scores = scores or get_calibration_scores(backend)
    coupling = list(backend.coupling_map) if not is_simulator(backend) else [
        [i, i+1] for i in range(29)
    ]
    idx_to_name = _idx_to_name_map(backend)

    # Build adjacency
    adj: dict[int, set[int]] = {}
    for a, b in coupling:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    def score(idx: int) -> float:
        return scores.get(idx_to_name.get(idx, ""), 0.0)

    best_chain: list[int] = []
    best_score = -1.0

    for start in adj:
        chain = [start]
        visited = {start}
        while len(chain) < n:
            last = chain[-1]
            neighbours = [nb for nb in adj.get(last, []) if nb not in visited]
            if not neighbours:
                break
            nxt = max(neighbours, key=score)
            chain.append(nxt)
            visited.add(nxt)

        if len(chain) == n:
            total = sum(score(q) for q in chain)
            if total > best_score:
                best_score = total
                best_chain = chain

    if not best_chain:
        raise ValueError(f"Could not find a chain of length {n} in coupling map.")

    names = [idx_to_name[i] for i in best_chain]
    edges = [[best_chain[i], best_chain[i+1]] for i in range(len(best_chain)-1)]
    return QubitLayout(qubit_names=names, qubit_indices=best_chain,
                       coupling_map=edges, scores={n: scores[n] for n in names})


def get_best_grid(backend, rows: int, cols: int,
                  scores: dict[str, float] | None = None) -> QubitLayout:
    """Find a rows×cols rectangular subgraph in the coupling map.

    IQM Emerald is a square lattice so a perfect grid always exists.
    Picks the subgraph with the highest average qubit score.
    """
    scores = scores or get_calibration_scores(backend)
    idx_to_name = _idx_to_name_map(backend)
    name_to_idx = {v: k for k, v in idx_to_name.items()}

    if is_simulator(backend):
        # Build a synthetic 6×9 grid for testing
        total = rows * cols
        indices = list(range(total))
        names = [f"QB{i+1}" for i in indices]
        h_edges = [[r*cols+c, r*cols+c+1] for r in range(rows) for c in range(cols-1)]
        v_edges = [[r*cols+c, (r+1)*cols+c] for r in range(rows-1) for c in range(cols)]
        return QubitLayout(qubit_names=names, qubit_indices=indices,
                           coupling_map=h_edges+v_edges,
                           scores={n: 1.0 for n in names})

    coupling = list(backend.coupling_map)
    adj: dict[int, set[int]] = {}
    for a, b in coupling:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    def score(idx: int) -> float:
        return scores.get(idx_to_name.get(idx, ""), 0.0)

    # Try to find a valid grid via BFS/enumeration
    # For a square lattice, we look for a contiguous rows×cols rectangle
    best_layout: QubitLayout | None = None
    best_avg = -1.0

    all_qubits = sorted(adj.keys())

    for seed in all_qubits:
        # BFS to find a rows×cols grid starting from seed
        grid = _try_build_grid(seed, rows, cols, adj)
        if grid is None:
            continue
        flat = [grid[r][c] for r in range(rows) for c in range(cols)]
        avg = sum(score(q) for q in flat) / len(flat)
        if avg > best_avg:
            best_avg = avg
            h_edges = [[grid[r][c], grid[r][c+1]]
                       for r in range(rows) for c in range(cols-1)]
            v_edges = [[grid[r][c], grid[r+1][c]]
                       for r in range(rows-1) for c in range(cols)]
            names = [idx_to_name[q] for q in flat]
            best_layout = QubitLayout(
                qubit_names=names, qubit_indices=flat,
                coupling_map=h_edges+v_edges,
                scores={idx_to_name[q]: scores.get(idx_to_name[q], 0.5) for q in flat}
            )

    if best_layout is None:
        raise ValueError(f"Could not find a {rows}×{cols} grid in coupling map. "
                         "Try a smaller grid or use get_best_chain().")
    return best_layout


def _try_build_grid(seed: int, rows: int, cols: int,
                    adj: dict[int, set[int]]) -> list[list[int]] | None:
    """Attempt to build a rows×cols grid of qubits starting from seed.

    Returns a 2D list grid[row][col] of qubit indices, or None if impossible.
    Uses BFS along horizontal and vertical directions.
    """
    # First find a horizontal chain of length cols from seed
    row0 = _greedy_path(seed, cols, adj, forbidden=set())
    if row0 is None:
        return None

    grid = [row0]
    used = set(row0)

    for r in range(1, rows):
        prev_row = grid[r-1]
        new_row = []
        for c, parent in enumerate(prev_row):
            neighbours = adj.get(parent, set()) - used
            # Pick neighbour not already in new_row
            candidates = [nb for nb in neighbours if nb not in new_row]
            if not candidates:
                return None
            # Prefer neighbour that also connects to previous column's down-neighbour
            chosen = candidates[0]
            new_row.append(chosen)
            used.add(chosen)

        # Verify horizontal connectivity within new_row
        for c in range(len(new_row)-1):
            if new_row[c+1] not in adj.get(new_row[c], set()):
                return None
        grid.append(new_row)

    return grid


def _greedy_path(start: int, length: int, adj: dict[int, set[int]],
                 forbidden: set[int]) -> list[int] | None:
    path = [start]
    visited = {start} | forbidden
    while len(path) < length:
        last = path[-1]
        candidates = adj.get(last, set()) - visited
        if not candidates:
            return None
        path.append(next(iter(candidates)))
        visited.add(path[-1])
    return path


def transpile_for_backend(circuits, backend, coupling_map=None, optimization_level=3):
    """Transpile one or more circuits for the given backend."""
    return transpile(circuits, backend=backend,
                     coupling_map=coupling_map,
                     optimization_level=optimization_level)
