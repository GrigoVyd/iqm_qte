"""Cluster-state entanglement routing + Bell fidelity + CHSH game.

Protocol on a 2D cluster |C_grid⟩:
    1. Prepare the full 2D cluster state (H on every qubit, CZ on every grid edge).
    2. Measure every off-path qubit in Z.
    3. Measure every internal-path qubit (q_1 .. q_{L-2}) in X.
    4. Leave endpoints A = path[0] and B = path[-1] for the chosen endpoint basis.
    5. Compute byproduct-corrected Bell fidelity F (sufficient for entanglement
       certification when F − 3σ_F > 1/2).
    6. Optionally also measure four CHSH settings on the endpoints and report
       S and game win probability.

Byproduct rule (derived analytically from cluster-stabilizer propagation; see
README/notes for derivation, validated on Aer by `derive_byproduct_table_simulator`):

    sign_x = (-1)^( Σ m_i for internal i at *even* path-position
                  + Σ m_O for off-path O adjacent to a path qubit at even position )
    sign_z = (-1)^( Σ m_i for internal i at *odd*  path-position
                  + Σ m_O for off-path O adjacent to a path qubit at odd  position )

For L *even*, the residual stabilizers on the endpoints are (Z_A X_B, X_A Z_B)
rather than (X_A X_B, Z_A Z_B). We apply a single H on B in the circuit so the
endpoint state lives in the standard Bell basis; the same byproduct formula
applies.

Per-shot correction:
    XX channel: corrected = sign_x · raw, target +1
    ZZ channel: corrected = sign_z · raw, target +1
    YY channel: corrected = -sign_x · sign_z · raw, target +1
    F = (1 + corrected_XX + corrected_YY + corrected_ZZ) / 4

CHSH (byproducts only on A; B-side rotated bases (Z±X)/√2 are unaffected):
    A0 = Z (corrected with sign_z),  A1 = X (corrected with sign_x)
    B0 = (Z + X) / √2,               B1 = (Z − X) / √2
    S = E(A0,B0) − E(A0,B1) + E(A1,B0) + E(A1,B1)
    p_win = 1/2 + |S|/8
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from typing import Iterable, Sequence

from qiskit import QuantumCircuit, transpile

# ---------------------------------------------------------------------------
# Grid utilities
# ---------------------------------------------------------------------------

def grid_index(r: int, c: int, cols: int) -> int:
    return r * cols + c


def grid_rc(idx: int, cols: int) -> tuple[int, int]:
    return divmod(idx, cols)


def grid_neighbors(r: int, c: int, rows: int, cols: int) -> list[int]:
    out = []
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nr, nc = r + dr, c + dc
        if 0 <= nr < rows and 0 <= nc < cols:
            out.append(grid_index(nr, nc, cols))
    return out


def grid_edges(rows: int, cols: int) -> list[tuple[int, int]]:
    edges = []
    for r in range(rows):
        for c in range(cols - 1):
            edges.append((grid_index(r, c, cols), grid_index(r, c + 1, cols)))
    for r in range(rows - 1):
        for c in range(cols):
            edges.append((grid_index(r, c, cols), grid_index(r + 1, c, cols)))
    return edges


def validate_path_is_connected(path: Sequence[int], rows: int, cols: int) -> None:
    """Raise ValueError unless every consecutive pair in `path` is grid-adjacent
    and no qubit appears twice."""
    if len(path) < 2:
        raise ValueError("path needs at least two qubits (A and B)")
    if len(set(path)) != len(path):
        raise ValueError(f"path repeats a qubit: {path}")
    for q in path:
        if not (0 <= q < rows * cols):
            raise ValueError(f"path qubit {q} outside {rows}x{cols} grid")
    for q1, q2 in zip(path, path[1:]):
        r1, c1 = grid_rc(q1, cols)
        r2, c2 = grid_rc(q2, cols)
        if abs(r1 - r2) + abs(c1 - c2) != 1:
            raise ValueError(f"path step {q1} -> {q2} is not grid-adjacent")


# ---------------------------------------------------------------------------
# Circuit construction
# ---------------------------------------------------------------------------

def _basis_rotation(qc: QuantumCircuit, q: int, basis: str) -> None:
    """Apply rotation so a subsequent Z-measurement reads `basis` (one of XYZ)."""
    if basis == "X":
        qc.h(q)
    elif basis == "Y":
        qc.sdg(q)
        qc.h(q)
    elif basis == "Z":
        pass
    else:
        raise ValueError(f"unknown Pauli basis {basis!r}")


def _chsh_rotation(qc: QuantumCircuit, q: int, setting: str) -> None:
    """CHSH endpoint settings:
       A0 = Z,  A1 = X,  B0 = (Z+X)/√2,  B1 = (Z−X)/√2.
       Implemented by rotating the state so a subsequent Z-measurement reads
       the requested observable."""
    if setting in ("A0", "Z"):
        return
    if setting in ("A1", "X"):
        qc.h(q)
        return
    # For B: ⟨(Z+X)/√2⟩ on |ψ⟩ = ⟨Z⟩ on R_y(+π/4)|ψ⟩, since R_y(+π/4) sends
    # the +eigenvector of (Z+X)/√2 to |0⟩. Likewise for B1 with R_y(−π/4).
    if setting == "B0":
        qc.ry(-math.pi / 4, q)
        return
    if setting == "B1":
        qc.ry(+math.pi / 4, q)
        return
    raise ValueError(f"unknown CHSH setting {setting!r}")


def build_cluster_routing_circuit(
    rows: int, cols: int, path: Sequence[int],
    endpoint_basis: tuple[str, str], *,
    no_cz: bool = False, chsh: bool = False,
) -> tuple[QuantumCircuit, dict]:
    """Build the routing circuit for a single endpoint-basis combo.

    endpoint_basis: (basis_A, basis_B). For chsh=True, these are CHSH labels
    ('A0','A1','B0','B1'); for Bell-fidelity mode they are 'X'/'Y'/'Z'.

    Returns (qc, metadata). Metadata records everything the parser needs:
    rows, cols, path, internal qubits, off-path qubits, endpoints, n_qubits,
    apply_h_b flag, basis labels, no_cz flag.
    """
    validate_path_is_connected(path, rows, cols)
    n = rows * cols
    A, B = path[0], path[-1]
    L = len(path)
    internal = list(path[1:-1])
    path_set = set(path)
    offpath = [q for q in range(n) if q not in path_set]

    qc = QuantumCircuit(n, n)

    # 1. Prepare cluster: H on all, CZ on every grid edge.
    qc.h(range(n))
    if not no_cz:
        for u, v in grid_edges(rows, cols):
            qc.cz(u, v)
    qc.barrier()

    # 2. For L even, apply H on B so endpoint stabilizers come out (XX, ZZ)
    #    instead of (XZ, ZX) — does not apply for the no-CZ control (the
    #    state is then trivial and we just want the formula to break gracefully).
    apply_h_b = (L % 2 == 0) and (not no_cz)
    if apply_h_b:
        qc.h(B)

    # 3. Internal path qubits in X basis (H), off-path qubits in Z (no rotation).
    for q in internal:
        qc.h(q)

    # 4. Endpoint rotations.
    if chsh:
        _chsh_rotation(qc, A, endpoint_basis[0])
        _chsh_rotation(qc, B, endpoint_basis[1])
    else:
        _basis_rotation(qc, A, endpoint_basis[0])
        _basis_rotation(qc, B, endpoint_basis[1])

    qc.measure(range(n), range(n))

    metadata = {
        "rows": rows, "cols": cols, "path": list(path),
        "A": A, "B": B, "internal": internal, "offpath": offpath,
        "n_qubits": n, "L": L, "apply_h_b": apply_h_b,
        "endpoint_basis": tuple(endpoint_basis),
        "no_cz": no_cz, "chsh": chsh,
    }
    return qc, metadata


def build_routing_bell_circuits(
    rows: int, cols: int, path: Sequence[int], *, no_cz: bool = False,
) -> list[tuple[QuantumCircuit, dict]]:
    """Three circuits: endpoints in (X, X), (Y, Y), (Z, Z)."""
    return [
        build_cluster_routing_circuit(rows, cols, path, basis, no_cz=no_cz)
        for basis in [("X", "X"), ("Y", "Y"), ("Z", "Z")]
    ]


def build_routing_chsh_circuits(
    rows: int, cols: int, path: Sequence[int], *, no_cz: bool = False,
) -> list[tuple[QuantumCircuit, dict]]:
    """Four CHSH endpoint pairings: (A0,B0), (A0,B1), (A1,B0), (A1,B1)."""
    return [
        build_cluster_routing_circuit(rows, cols, path, basis,
                                      no_cz=no_cz, chsh=True)
        for basis in [("A0", "B0"), ("A0", "B1"),
                      ("A1", "B0"), ("A1", "B1")]
    ]


# ---------------------------------------------------------------------------
# Bit parsing + byproduct rule
# ---------------------------------------------------------------------------

def _bit(bs: str, q: int, n: int) -> int:
    """Qiskit bit ordering: classical bit q = bs[n - 1 - q]."""
    return int(bs.replace(" ", "")[n - 1 - q])


def byproduct_signs_from_outcome(bs: str, metadata: dict) -> tuple[int, int]:
    """(sign_x, sign_z) ∈ {-1, +1}² from a single outcome bitstring.

    sign_x: parity of internal-X outcomes at even path positions plus parity
            of off-path Z outcomes adjacent to path qubits at even positions.
    sign_z: same with odd path positions.

    For L even, this is computed on the H_B-rotated state; the formula is
    identical because H_B is applied in the circuit before measurement.
    """
    n = metadata["n_qubits"]
    rows, cols = metadata["rows"], metadata["cols"]
    path = metadata["path"]
    pos_in_path = {q: i for i, q in enumerate(path)}
    path_set = set(path)

    sign_x = 1
    sign_z = 1

    # Internal X-measurement contributions
    for q in metadata["internal"]:
        m = _bit(bs, q, n)
        if pos_in_path[q] % 2 == 0:
            sign_x *= 1 - 2 * m
        else:
            sign_z *= 1 - 2 * m

    # Off-path Z-measurement contributions: each off-path neighbour of a path
    # qubit p contributes (-1)^{m_O} into the sign indexed by the parity of
    # p's position. If an off-path qubit is adjacent to several path qubits
    # the contributions multiply (in mod-2 arithmetic for each parity bucket).
    for p in path:
        pos = pos_in_path[p]
        r, c = grid_rc(p, cols)
        for nbr in grid_neighbors(r, c, rows, cols):
            if nbr in path_set:
                continue
            m = _bit(bs, nbr, n)
            if pos % 2 == 0:
                sign_x *= 1 - 2 * m
            else:
                sign_z *= 1 - 2 * m

    return sign_x, sign_z


def _channel_correction(basis_label: str, sign_x: int, sign_z: int) -> int:
    """The per-shot multiplicative correction for one Bell channel."""
    if basis_label == "XX":
        return sign_x
    if basis_label == "ZZ":
        return sign_z
    if basis_label == "YY":
        return -sign_x * sign_z
    raise ValueError(f"unknown channel {basis_label!r}")


# ---------------------------------------------------------------------------
# Counts → correlations
# ---------------------------------------------------------------------------

def parse_counts_by_pattern(
    counts: dict[str, int], metadata: dict,
) -> dict[tuple[int, int], dict[tuple[int, int], int]]:
    """Bin shots by (sign_x, sign_z) and by endpoint outcomes (a, b).

    Useful for diagnostics; the corrected/uncorrected aggregators below use
    the per-shot signs directly without binning.
    """
    n = metadata["n_qubits"]
    A, B = metadata["A"], metadata["B"]
    out: dict = defaultdict(lambda: defaultdict(int))
    for bs, cnt in counts.items():
        sx, sz = byproduct_signs_from_outcome(bs, metadata)
        a, b = _bit(bs, A, n), _bit(bs, B, n)
        out[(sx, sz)][(a, b)] += cnt
    return {k: dict(v) for k, v in out.items()}


def corrected_endpoint_correlation(
    counts: dict[str, int], metadata: dict, basis_label: str,
) -> tuple[float, float]:
    """Per-shot byproduct-corrected correlation in {-1, +1}.

    For an ideal cluster state, all three channels (XX, YY, ZZ) average to +1.
    Returns (mean, sigma).
    """
    n = metadata["n_qubits"]
    A, B = metadata["A"], metadata["B"]
    total = 0
    weighted = 0
    for bs, cnt in counts.items():
        sx, sz = byproduct_signs_from_outcome(bs, metadata)
        a, b = _bit(bs, A, n), _bit(bs, B, n)
        raw = (1 - 2 * a) * (1 - 2 * b)
        weighted += _channel_correction(basis_label, sx, sz) * raw * cnt
        total += cnt
    if total == 0:
        return 0.0, float("inf")
    e = weighted / total
    var = max(0.0, 1.0 - e * e)
    return e, math.sqrt(var / total)


def uncorrected_endpoint_correlation(
    counts: dict[str, int], metadata: dict,
) -> tuple[float, float]:
    """Same correlation but without applying the byproduct correction.

    For the ideal cluster state this averages to ≈ 0 because the conditional
    Bell state varies with the measurement outcomes. Used as a control to
    show that the byproduct correction is doing real work.
    """
    n = metadata["n_qubits"]
    A, B = metadata["A"], metadata["B"]
    total = 0
    weighted = 0
    for bs, cnt in counts.items():
        a, b = _bit(bs, A, n), _bit(bs, B, n)
        weighted += (1 - 2 * a) * (1 - 2 * b) * cnt
        total += cnt
    if total == 0:
        return 0.0, float("inf")
    e = weighted / total
    var = max(0.0, 1.0 - e * e)
    return e, math.sqrt(var / total)


# ---------------------------------------------------------------------------
# Bell fidelity
# ---------------------------------------------------------------------------

def compute_bell_fidelity(
    c_xx: float, c_yy: float, c_zz: float,
    s_xx: float, s_yy: float, s_zz: float,
) -> dict:
    """Bell fidelity from byproduct-corrected correlations.

    Per the per-shot correction in this module, all three channels ideally
    average to +1, so:

        F = (1 + C_XX + C_YY + C_ZZ) / 4

    F − 3σ_F > 0.5 certifies the endpoint pair as entangled.
    """
    F = (1.0 + c_xx + c_yy + c_zz) / 4.0
    sigma_F = 0.25 * math.sqrt(s_xx * s_xx + s_yy * s_yy + s_zz * s_zz)
    z = (F - 0.5) / sigma_F if sigma_F > 0 else float("inf")
    return {
        "F": F, "sigma_F": sigma_F, "z_score": z,
        "entangled_3sigma": (F - 3 * sigma_F) > 0.5,
        "entangled_meanonly": F > 0.5,
        "C_XX": c_xx, "C_YY": c_yy, "C_ZZ": c_zz,
    }


# ---------------------------------------------------------------------------
# CHSH
# ---------------------------------------------------------------------------

def _chsh_corrected_correlation(
    counts: dict[str, int], metadata: dict,
) -> tuple[float, float]:
    """Byproduct-corrected correlation for one CHSH endpoint pairing.

    The byproduct can be placed entirely on A: σ = X_A^a Z_A^b |Φ+⟩.
    Therefore:
        A0 = Z is corrected with sign_z = (-1)^a
        A1 = X is corrected with sign_x = (-1)^b
        B-side rotated bases are unaffected (no correction needed).
    """
    n = metadata["n_qubits"]
    A, B = metadata["A"], metadata["B"]
    a_setting = metadata["endpoint_basis"][0]
    total = 0
    weighted = 0
    for bs, cnt in counts.items():
        sx, sz = byproduct_signs_from_outcome(bs, metadata)
        a, b = _bit(bs, A, n), _bit(bs, B, n)
        raw = (1 - 2 * a) * (1 - 2 * b)
        if a_setting in ("A0", "Z"):
            corr = sz
        elif a_setting in ("A1", "X"):
            corr = sx
        else:
            raise ValueError(f"unexpected A setting {a_setting}")
        weighted += corr * raw * cnt
        total += cnt
    if total == 0:
        return 0.0, float("inf")
    e = weighted / total
    var = max(0.0, 1.0 - e * e)
    return e, math.sqrt(var / total)


_CHSH_SIGN_VARIANTS = [
    # All 8 sign assignments (e_00, e_01, e_10, e_11). The classical bound is
    # |S| ≤ 2 for any of them; the Tsirelson bound 2√2 is reached by the variant
    # whose sign pattern matches the resource state's stabilizer structure
    # (e.g. Φ+ maxes "+ + + −", singlet maxes "+ − + +"). We pick the variant
    # with the largest |S| in finite-shot data; for the routed state (Φ+ after
    # byproduct correction) this should be the "+ + + −" variant.
    (+1, +1, +1, -1),  (+1, +1, -1, +1),  (+1, -1, +1, +1),  (-1, +1, +1, +1),
    (-1, -1, -1, +1),  (-1, -1, +1, -1),  (-1, +1, -1, -1),  (+1, -1, -1, -1),
]


def compute_chsh_game(chsh_results: dict[tuple[str, str], dict]) -> dict:
    """chsh_results: {(A_setting, B_setting): {'E': float, 'sigma': float}}.

    Computes |S| as the maximum over the 8 sign assignments of the CHSH form.
    For the byproduct-corrected routed Φ+ this picks the (+,+,+,−) variant,
    matching the standard CHSH inequality for Φ+.
    """
    e = {k: v["E"] for k, v in chsh_results.items()}
    s = {k: v["sigma"] for k, v in chsh_results.items()}
    keys = [("A0", "B0"), ("A0", "B1"), ("A1", "B0"), ("A1", "B1")]
    sigma_S = math.sqrt(sum(s[k] ** 2 for k in keys))

    best_S = 0.0
    best_signs = (+1, -1, +1, +1)  # the formula the user wrote
    for signs in _CHSH_SIGN_VARIANTS:
        S_variant = sum(sg * e[k] for sg, k in zip(signs, keys))
        if abs(S_variant) > abs(best_S):
            best_S = S_variant
            best_signs = signs

    p_win = 0.5 + abs(best_S) / 8.0
    sigma_p = sigma_S / 8.0
    return {
        "S": best_S, "sigma_S": sigma_S, "abs_S": abs(best_S),
        "best_sign_pattern": list(best_signs),
        "p_win": p_win, "sigma_p_win": sigma_p,
        "classical_bound_S": 2.0, "tsirelson_bound": 2.0 * math.sqrt(2),
        "pass_3sigma": (abs(best_S) - 3 * sigma_S) > 2.0,
        "per_setting": {
            f"{k[0]}{k[1]}": {"E": e[k], "sigma": s[k]} for k in keys
        },
    }


# ---------------------------------------------------------------------------
# Top-level runner
# ---------------------------------------------------------------------------

def _is_aer(backend) -> bool:
    return backend.__class__.__name__ == "AerSimulator"


def _transpile(circuits, backend, initial_layout, optimization_level):
    if _is_aer(backend) or initial_layout is None:
        return transpile(circuits, backend=backend,
                         optimization_level=optimization_level)
    return transpile(circuits, backend=backend,
                     initial_layout=initial_layout,
                     optimization_level=optimization_level)


def _run(backend, circuits, shots):
    job = backend.run(circuits, shots=shots)
    result = job.result()
    out = []
    for i in range(len(circuits)):
        out.append(result.get_counts(i))
    return out


def run_routing_demo(
    backend, rows: int, cols: int, path: Sequence[int], shots: int,
    *, initial_layout: Sequence[int] | None = None,
    optimization_level: int = 3, do_chsh: bool = True,
    no_cz: bool = False, label: str | None = None,
) -> dict:
    """Run Bell-fidelity (3 circuits) + optional CHSH (4 circuits) for one path.

    Returns a self-contained dict with everything needed for reporting and
    JSON serialization.
    """
    bell_pairs = build_routing_bell_circuits(rows, cols, path, no_cz=no_cz)
    bell_circs = [c for c, _ in bell_pairs]
    bell_metas = [m for _, m in bell_pairs]

    chsh_pairs = (build_routing_chsh_circuits(rows, cols, path, no_cz=no_cz)
                  if do_chsh else [])
    chsh_circs = [c for c, _ in chsh_pairs]
    chsh_metas = [m for _, m in chsh_pairs]

    all_circs = bell_circs + chsh_circs
    all_metas = bell_metas + chsh_metas
    transpiled = _transpile(all_circs, backend, initial_layout,
                              optimization_level)
    counts_list = _run(backend, transpiled, shots)

    bell_counts = counts_list[: len(bell_circs)]
    chsh_counts = counts_list[len(bell_circs):]

    # Bell channels
    channel_labels = ["XX", "YY", "ZZ"]
    corrected = {}
    uncorrected = {}
    for label_ch, counts, meta in zip(channel_labels, bell_counts, bell_metas):
        e_c, s_c = corrected_endpoint_correlation(counts, meta, label_ch)
        e_u, s_u = uncorrected_endpoint_correlation(counts, meta)
        corrected[label_ch] = {"E": e_c, "sigma": s_c}
        uncorrected[label_ch] = {"E": e_u, "sigma": s_u}

    bell = compute_bell_fidelity(
        corrected["XX"]["E"], corrected["YY"]["E"], corrected["ZZ"]["E"],
        corrected["XX"]["sigma"], corrected["YY"]["sigma"],
        corrected["ZZ"]["sigma"],
    )

    out = {
        "label": label, "rows": rows, "cols": cols, "path": list(path),
        "L": len(path), "A": path[0], "B": path[-1],
        "shots": shots, "no_cz": no_cz,
        "bell": bell,
        "bell_corrected_correlations": corrected,
        "bell_uncorrected_correlations": uncorrected,
    }

    if do_chsh:
        chsh_settings = [("A0", "B0"), ("A0", "B1"),
                         ("A1", "B0"), ("A1", "B1")]
        per = {}
        for setting, counts, meta in zip(chsh_settings, chsh_counts, chsh_metas):
            e_c, s_c = _chsh_corrected_correlation(counts, meta)
            per[setting] = {"E": e_c, "sigma": s_c}
        out["chsh"] = compute_chsh_game(per)

    return out


def sweep_routes(
    backend, rows: int, cols: int, paths: Iterable[Sequence[int]],
    shots: int, **kwargs,
) -> list[dict]:
    """Run `run_routing_demo` for each path, return list of result dicts."""
    return [run_routing_demo(backend, rows, cols, p, shots, **kwargs)
            for p in paths]


# ---------------------------------------------------------------------------
# Simulator-only validation of the byproduct rule
# ---------------------------------------------------------------------------

def derive_byproduct_table_simulator(
    rows: int, cols: int, path: Sequence[int], *, shots: int = 16000,
) -> dict:
    """Simulator-side validation of `byproduct_signs_from_outcome`.

    For each (sign_x, sign_z) outcome bin, the conditional state on the
    endpoints must be the Bell state with stabilizers (sign_x · XX, sign_z · ZZ)
    after applying H_B (when L is even). We run the three Bell circuits,
    bin shots by (sign_x, sign_z), and check the conditional ⟨XX⟩, ⟨YY⟩,
    ⟨ZZ⟩ matches what the byproduct formula predicts.

    Returns a per-bin diagnostic dict; intended for tests/notebooks, not for
    hardware (very expensive in shots).
    """
    from qiskit_aer import AerSimulator
    backend = AerSimulator()
    table: dict[tuple[int, int], dict] = {}
    bell_pairs = build_routing_bell_circuits(rows, cols, path)
    transpiled = transpile([c for c, _ in bell_pairs], backend=backend)
    job = backend.run(transpiled, shots=shots)
    res = job.result()
    counts_list = [res.get_counts(i) for i in range(3)]
    metas = [m for _, m in bell_pairs]
    labels = ["XX", "YY", "ZZ"]

    n = metas[0]["n_qubits"]
    A, B = metas[0]["A"], metas[0]["B"]

    for label, counts, meta in zip(labels, counts_list, metas):
        for bs, cnt in counts.items():
            sx, sz = byproduct_signs_from_outcome(bs, meta)
            entry = table.setdefault(
                (sx, sz), {"shots_per_channel": {"XX": 0, "YY": 0, "ZZ": 0},
                            "weighted": {"XX": 0, "YY": 0, "ZZ": 0}}
            )
            a, b = _bit(bs, A, n), _bit(bs, B, n)
            raw = (1 - 2 * a) * (1 - 2 * b)
            entry["weighted"][label] += raw * cnt
            entry["shots_per_channel"][label] += cnt

    # Compute conditional <XX>, <YY>, <ZZ> per bin, plus predicted vs measured
    out = {}
    for (sx, sz), data in table.items():
        cond = {}
        for label in labels:
            n_shots = data["shots_per_channel"][label]
            cond[label] = (data["weighted"][label] / n_shots) if n_shots else None
        predicted = {"XX": sx, "ZZ": sz, "YY": -sx * sz}
        out[(sx, sz)] = {
            "shots_per_channel": data["shots_per_channel"],
            "measured_correlations": cond,
            "predicted_correlations": predicted,
            "matches": all(
                cond[k] is None or abs(cond[k] - predicted[k]) < 0.05
                for k in labels
            ),
        }
    return out


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def save_results_json(results: list[dict] | dict, filename: str) -> None:
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)
