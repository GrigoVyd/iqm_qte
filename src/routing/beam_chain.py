"""Beam-search Hamiltonian-path routing for chain-shaped circuits (W-states).

Given the device coupling map and live calibration, find the length-n
hardware path that maximises the multiplicative log-fidelity score
(CZ x readout x T1 x T2 surviving the chain duration). Greedy beam-search
keeps the best `beam_width` partial paths at each step.

This wins over the IQM Qubit Selector when the circuit's interaction graph
is itself a Hamiltonian path: the selector treats the layout as an
optimisation over arbitrary subgraphs, while beam-search exploits the
chain structure directly and never inserts SWAPs.
"""
from __future__ import annotations

import math
from typing import Callable, Iterable


def beam_search_chain(
    n: int,
    adj: dict[int, set[int]],
    cz_fid: dict[frozenset, float],
    ro_fid: dict[int, float],
    t1_fn: Callable[[int, float], float],
    t2_fn: Callable[[int, float], float],
    gate_ns: float = 80.0,
    beam_width: int = 100,
    default_cz: float = 0.90,
    default_ro: float = 0.95,
) -> tuple[list[int], float]:
    """Return (chain, multiplicative_score) for the best length-n simple path.

    `t1_fn(phys, wait_ns) -> survival_probability` is the per-qubit decay
    factor over the wait time `wait_ns`. Same for `t2_fn`. The score
    accumulates log-survival per qubit at the *remaining* idle time it must
    survive after being prepared.
    """
    def ls(x: float) -> float:
        return math.log(max(x, 1e-12))

    # Initial beam: each qubit as a single-vertex path. Score = readout +
    # decoherence over the full circuit duration (qubit must survive until
    # the end).
    full_wait = (n - 1) * 3 * gate_ns
    beam = [
        (
            ls(ro_fid.get(q, default_ro))
            + ls(t1_fn(q, full_wait))
            + ls(t2_fn(q, full_wait)),
            [q],
        )
        for q in adj
    ]
    beam.sort(reverse=True)
    beam = beam[:beam_width]

    for step in range(n - 1):
        pos = step + 1
        wait_ns = (n - 1 - pos) * 3 * gate_ns
        candidates = []
        for log_score, path in beam:
            last = path[-1]
            used = set(path)
            for nbr in adj.get(last, ()):
                if nbr in used:
                    continue
                new_log = (
                    log_score
                    + ls(cz_fid.get(frozenset((last, nbr)), default_cz))
                    + ls(ro_fid.get(nbr, default_ro))
                    + ls(t1_fn(nbr, wait_ns))
                    + ls(t2_fn(nbr, wait_ns))
                )
                candidates.append((new_log, path + [nbr]))
        if not candidates:
            break
        candidates.sort(reverse=True)
        beam = candidates[:beam_width]

    complete = [(s, p) for s, p in beam if len(p) == n]
    if not complete:
        raise RuntimeError(
            f"No hardware path of length {n} found "
            f"(reduce n or increase beam_width)."
        )
    best_log, best_path = complete[0]
    return best_path, math.exp(best_log)


def calibration_to_routing_inputs(backend, cal: dict) -> dict:
    """Convert IQM CalibrationDataManager output into the dicts expected by
    `beam_search_chain` (CZ keyed by frozenset of qubit indices etc.).

    `cal` is the dict returned by
    `iqm.qubit_selector.qubit_selector.CalibrationDataManager().get_calibration_fidelities(backend)`.
    """
    def _parse_pair(d):
        out = {}
        for k, v in d.items():
            try:
                names = eval(k)
                idxs = tuple(backend.qubit_name_to_index(nm) for nm in names)
                out[frozenset(idxs)] = float(v)
            except Exception:
                pass
        return out

    def _parse_qubit(d):
        out = {}
        for nm, v in d.items():
            try:
                out[backend.qubit_name_to_index(nm)] = float(v)
            except Exception:
                pass
        return out

    return {
        "cz_fid": _parse_pair(cal.get("CZ", {})),
        "ro_fid": _parse_qubit(cal.get("readout", {})),
        "t1_us": _parse_qubit(cal.get("t1", {})),
        "t2_us": _parse_qubit(cal.get("t2", {})),
    }


def make_decoherence_fns(t1_us: dict[int, float], t2_us: dict[int, float],
                          default_t1_us: float = 200.0,
                          default_t2_us: float = 100.0):
    """Closures: phys, wait_ns -> survival probability."""
    import math as _math

    def t1_fn(phys: int, wait_ns: float) -> float:
        return _math.exp(-wait_ns / (t1_us.get(phys, default_t1_us) * 1e3))

    def t2_fn(phys: int, wait_ns: float) -> float:
        return _math.exp(-wait_ns / (t2_us.get(phys, default_t2_us) * 1e3))

    return t1_fn, t2_fn


def adjacency_from_backend(backend) -> dict[int, set[int]]:
    adj: dict[int, set[int]] = {}
    for u, v in backend.coupling_map.get_edges():
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    return adj


def mean_cz_duration_ns(backend, default_ns: float = 80.0) -> float:
    durations = []
    target = getattr(backend, "target", None)
    if target is not None and "cz" in target:
        for _, props in target["cz"].items():
            if props and getattr(props, "duration", None):
                durations.append(props.duration)
    return (sum(durations) / len(durations) * 1e9) if durations else default_ns
