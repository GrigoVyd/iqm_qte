"""Routed Bell pair via measurement-based teleportation on a cluster state.

Protocol (chain version, 2D cluster handled as the path subset):

    A -- q_1 -- q_2 -- ... -- q_{L-2} -- B

1. Prepare the linear cluster on the path:  H on every path qubit, CZ on
   adjacent path pairs.
2. (For 2D clusters) measure all *off-path* qubits in Z to disentangle them
   from the path.
3. Measure all internal path qubits q_1, ..., q_{L-2} in X.
4. Measure the endpoints A, B in the chosen Pauli basis (XX, YY or ZZ for
   Bell fidelity; or rotated bases for CHSH).
5. The endpoint state is locally equivalent to a Bell pair, modulo a Pauli
   byproduct that depends on the internal X-outcomes. We correct the
   correlators with sign factors derived from a noiseless statevector
   simulation of the same circuit.

Conditional sign correction: for every shot, look up s_O(m_internals) — the
ideal correlator value under the noiseless protocol with the same
internal-X bits. Multiply the observed ±1 endpoint outcome by s, and
average. The result is the *byproduct-corrected* correlator on the Bell
target |Φ+⟩ (which we engineer to be the noiseless target by an extra
H_B at the end of cluster prep).

Bell fidelity for |Φ+⟩:

    F = (1 + ⟨X_A X_B⟩_c - ⟨Y_A Y_B⟩_c + ⟨Z_A Z_B⟩_c) / 4

with shot-noise

    σ(C) = sqrt((1 - C^2) / N),
    σ_F  = (1/4) sqrt(σ_XX^2 + σ_YY^2 + σ_ZZ^2).

Entanglement is certified when F − 3 σ_F > 1/2 (the separable bound).

Optional CHSH game: measure endpoints in
    A_0 = Z, A_1 = X, B_0 = (Z+X)/√2, B_1 = (Z-X)/√2
and compute
    S = E(A0,B0) − E(A0,B1) + E(A1,B0) + E(A1,B1).
Pass: |S| − 3 σ_S > 2.
"""

from __future__ import annotations

import math
from itertools import product
from typing import Iterable

import numpy as np
from qiskit import QuantumCircuit, transpile

# CHSH endpoint Pauli combinations: (A_axis, B_axis) where A_axis is "Z" or
# "X" and B_axis is one of the rotated bases. We measure these as physical
# rotations applied before a Z-measurement.
_CHSH_BASES: dict[tuple[int, int], tuple[str, float]] = {
    # (a, b) -> (A_pauli_label, B_rotation_angle) where the angle is theta in
    #          R_y(-2θ) applied to B before a Z-measurement.
    # B0 = (Z+X)/√2 = R_y(π/4) Z R_y(-π/4) → measure in basis rotated by π/8
    # B1 = (Z-X)/√2 = R_y(-π/4) Z R_y(π/4) → rotated by -π/8
    (0, 0): ("Z", math.pi / 8),
    (0, 1): ("Z", -math.pi / 8),
    (1, 0): ("X", math.pi / 8),
    (1, 1): ("X", -math.pi / 8),
}


# ---------------------------------------------------------------------------
# Circuit construction
# ---------------------------------------------------------------------------

def build_path_cluster_circuit(
    path: list[int],
    *,
    endpoint_basis: str = "ZZ",
    chsh_setting: tuple[int, int] | None = None,
    cluster_qubits: Iterable[int] | None = None,
    n_total_qubits: int | None = None,
    skip_cz: bool = False,
) -> tuple[QuantumCircuit, list[int]]:
    """Build a cluster + measurement circuit for routing entanglement along
    `path` from endpoint A=path[0] to endpoint B=path[-1].

    `cluster_qubits` is the *full* set of cluster qubits to prepare (path ∪
    off-path neighbours for a 2D cluster); off-path qubits will be measured in
    Z. Defaults to `path` itself (pure 1D path cluster).

    `endpoint_basis`: "XX", "YY", "ZZ" — Pauli basis for endpoint Bell-fidelity
    measurement. Ignored when `chsh_setting` is given.

    `chsh_setting`: optional (a, b) ∈ {0,1}^2 selecting one of the four CHSH
    bases (A∈{Z,X} on endpoint A; B rotated to (Z±X)/√2 on endpoint B).

    `skip_cz`: control variant — same circuit minus the cluster CZ layer.
    Bell fidelity should not exceed 0.5 in this case.

    Returns (circuit, clbit_layout) where clbit_layout[i] is the qubit
    measured into classical bit i. The endpoints are clbit 0 and clbit
    len(path)-1; internals are clbits 1..len(path)-2; off-path qubits go
    after that.
    """
    if endpoint_basis not in {"XX", "YY", "ZZ"} and chsh_setting is None:
        raise ValueError(f"endpoint_basis must be XX/YY/ZZ, got {endpoint_basis!r}")
    cluster_qubits = list(cluster_qubits) if cluster_qubits is not None else list(path)
    if path[0] not in cluster_qubits or path[-1] not in cluster_qubits:
        raise ValueError("path endpoints must be in cluster_qubits")
    for q in path:
        if q not in cluster_qubits:
            raise ValueError(f"path qubit {q} not in cluster_qubits")

    if n_total_qubits is None:
        n_total_qubits = max(cluster_qubits) + 1
    n_meas = len(cluster_qubits)
    qc = QuantumCircuit(n_total_qubits, n_meas)

    # 1. Prepare cluster: H on every cluster qubit, CZ on every cluster edge.
    for q in cluster_qubits:
        qc.h(q)
    if not skip_cz:
        # Build the cluster's edge set: pairs in `cluster_qubits` that are
        # adjacent on the path or otherwise in cluster_edges (only path-adjacent
        # for the 1D protocol; 2D users pass a richer cluster_qubits set and
        # supply edges via the wrapper that builds them).
        for i in range(len(path) - 1):
            qc.cz(path[i], path[i + 1])
    qc.barrier()

    # 2. Off-path qubits → Z basis (no rotation needed; just measure later).
    off_path = [q for q in cluster_qubits if q not in path]

    # 3. Internal path qubits → X basis.
    internals = path[1:-1]
    for q in internals:
        qc.h(q)

    # 3b. Parity fix-up. The 1D-cluster wire applies H^{L-2} to the propagated
    # state on B (independent of m, mod Pauli byproducts). For L odd the
    # residual is H, taking the cluster's CZ|++⟩ into |Φ+⟩ already. For L even
    # the residual is I and we need an extra H on B to land on |Φ+⟩.
    if (len(path) % 2) == 0:
        qc.h(path[-1])

    # 4. Endpoint rotations.
    A, B = path[0], path[-1]
    if chsh_setting is None:
        for q, p in ((A, endpoint_basis[0]), (B, endpoint_basis[1])):
            if p == "X":
                qc.h(q)
            elif p == "Y":
                qc.sdg(q); qc.h(q)
            elif p == "Z":
                pass
    else:
        a_label, theta = _CHSH_BASES[chsh_setting]
        if a_label == "X":
            qc.h(A)
        # B rotation: measure in basis (cosθ Z + sinθ X) ≡ rotate by -2θ around y, then Z-measure.
        qc.ry(-2 * theta, B)

    # 5. Measure: clbit ordering = [endpoints A, internals..., endpoint B,
    # off-path...]. We fix it by measuring path[0]..path[-1] in order, then
    # off-path. This places A at clbit 0, internals at clbits 1..L-2, B at
    # clbit L-1.
    clbit_layout: list[int] = []
    for clbit, q in enumerate(path):
        qc.measure(q, clbit)
        clbit_layout.append(q)
    for clbit, q in enumerate(off_path, start=len(path)):
        qc.measure(q, clbit)
        clbit_layout.append(q)

    return qc, clbit_layout


# ---------------------------------------------------------------------------
# Sign table from noiseless simulation
# ---------------------------------------------------------------------------

def _statevec_after_cluster(path_len: int, *, skip_cz: bool = False) -> np.ndarray:
    """Statevector after cluster prep + internal-H rotation, no measurements.

    Returned in Qiskit order: index = sum_q (bit_q << q), i.e. q_0 = LSB.
    """
    qc = QuantumCircuit(path_len)
    for q in range(path_len):
        qc.h(q)
    if not skip_cz:
        for i in range(path_len - 1):
            qc.cz(i, i + 1)
    for k in range(1, path_len - 1):
        qc.h(k)  # X-basis rotation on internals
    # Same parity fix-up as build_path_cluster_circuit.
    if (path_len % 2) == 0:
        qc.h(path_len - 1)
    from qiskit.quantum_info import Statevector
    return Statevector.from_instruction(qc).data


_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _project_endpoints(sv: np.ndarray, path_len: int,
                        m_tuple: tuple[int, ...]) -> np.ndarray:
    """Slice the statevector to the 4 amplitudes of (q_0, q_{path_len-1})
    consistent with internal outcomes m_tuple. Returns un-normalised 4-vector
    indexed by 2*j_B + j_A."""
    out = np.zeros(4, dtype=complex)
    for ja, jb in product((0, 1), repeat=2):
        bits = [ja] + list(m_tuple) + [jb]
        idx = sum(b << i for i, b in enumerate(bits))
        out[2 * jb + ja] = sv[idx]
    return out


def _correlator(state4: np.ndarray, op_a: str, op_b: str) -> float:
    """⟨ψ| O_A ⊗ O_B |ψ⟩ on a 2-qubit state with A=q_0 (LSB), B=q_1."""
    norm = float(np.real(state4.conj() @ state4))
    if norm < 1e-12:
        return 0.0
    psi = state4 / np.sqrt(norm)
    O = np.kron(_PAULI[op_b], _PAULI[op_a])
    return float(np.real(psi.conj() @ O @ psi))


def derive_sign_table(path_len: int, settings=("XX", "YY", "ZZ"),
                       skip_cz: bool = False) -> dict[tuple, dict[str, int]]:
    """For each internal-outcome tuple m, run noiseless statevector sim and
    return the *ideal* correlator sign for each setting.

    sign = +1 / -1 if |⟨O⟩| ≈ 1 noiselessly; 0 if the correlator is not
    pinned (the Bell sector for that m doesn't determine that Pauli).
    """
    n_internal = max(0, path_len - 2)
    sv = _statevec_after_cluster(path_len, skip_cz=skip_cz)
    table: dict[tuple, dict[str, int]] = {}
    if n_internal == 0:
        m_iter = [()]
    else:
        m_iter = list(product((0, 1), repeat=n_internal))
    for m_tuple in m_iter:
        ψ = _project_endpoints(sv, path_len, m_tuple)
        if float(np.real(ψ.conj() @ ψ)) < 1e-12:
            continue
        signs: dict[str, int] = {}
        for s in settings:
            c = _correlator(ψ, s[0], s[1])
            if abs(c) < 1e-6:
                signs[s] = 0
            else:
                signs[s] = 1 if c > 0 else -1
        table[m_tuple] = signs
    return table


def derive_chsh_sign_table(path_len: int) -> dict[tuple, dict[tuple[int, int], float]]:
    """For each m_tuple, return the ideal CHSH correlator value E(a,b) for each
    of the four settings on the noiseless protocol. Used for byproduct
    correction in the CHSH game.
    """
    n_internal = max(0, path_len - 2)
    sv = _statevec_after_cluster(path_len)
    table: dict[tuple, dict[tuple[int, int], float]] = {}
    m_iter = [()] if n_internal == 0 else list(product((0, 1), repeat=n_internal))
    for m_tuple in m_iter:
        ψ = _project_endpoints(sv, path_len, m_tuple)
        if float(np.real(ψ.conj() @ ψ)) < 1e-12:
            continue
        E: dict[tuple[int, int], float] = {}
        # CHSH ideal correlators on a Bell pair: E(A,B) = ±1/√2 depending on
        # the basis combination. We compute by direct simulation.
        for ab in ((0, 0), (0, 1), (1, 0), (1, 1)):
            a_label, theta = _CHSH_BASES[ab]
            # A axis: Z (a=0) or X (a=1). B axis: cosθ Z + sinθ X.
            A = _PAULI[a_label]
            B = math.cos(2 * theta) * _PAULI["Z"] + math.sin(2 * theta) * _PAULI["X"]
            O = np.kron(B, A)
            norm = float(np.real(ψ.conj() @ ψ))
            if norm < 1e-12:
                E[ab] = 0.0
            else:
                psi = ψ / math.sqrt(norm)
                E[ab] = float(np.real(psi.conj() @ O @ psi))
        table[m_tuple] = E
    return table


# ---------------------------------------------------------------------------
# Counts → corrected correlators
# ---------------------------------------------------------------------------

def _strip(bs: str) -> str:
    return bs.replace(" ", "")


def _read_outcome(bitstring: str, clbit: int) -> int:
    """Bit at classical-bit index clbit from a Qiskit-style bitstring (LSB first
    after stripping whitespace; Qiskit prints MSB first)."""
    return int(bitstring[len(bitstring) - 1 - clbit])


def corrected_correlator(counts: dict[str, int], path_len: int,
                          setting: str, sign_table: dict[tuple, dict[str, int]]
                          ) -> tuple[float, float, int]:
    """⟨O_A O_B⟩ corrected for the byproduct sector.

    For each shot:
      m = tuple of internal X outcomes (clbits 1..L-2)
      raw_o = (-1)^(end_A + end_B) where end_X are clbits 0 and L-1
      corrected = sign_table[m][setting] * raw_o
    Average over shots that have a defined sign.

    Returns (corrected_mean, sigma, n_used).
    """
    total_w = 0.0
    total_n = 0
    discarded = 0
    for bs, c in counts.items():
        b = _strip(bs)
        if len(b) < path_len:
            continue
        m = tuple(_read_outcome(b, k) for k in range(1, path_len - 1))
        end_A = _read_outcome(b, 0)
        end_B = _read_outcome(b, path_len - 1)
        raw = 1 - 2 * (end_A ^ end_B)
        signs = sign_table.get(m)
        if signs is None or signs.get(setting, 0) == 0:
            discarded += c
            continue
        s = signs[setting]
        total_w += s * raw * c
        total_n += c
    if total_n == 0:
        return 0.0, float("inf"), 0
    mean = total_w / total_n
    sigma = math.sqrt(max(0.0, 1.0 - mean * mean) / total_n)
    return mean, sigma, total_n


def raw_correlator(counts: dict[str, int], path_len: int) -> tuple[float, float, int]:
    """Uncorrected ⟨O_A O_B⟩ — average raw endpoint parity, no byproduct sign
    correction. Should be biased / noisy compared to corrected; useful as a
    control to demonstrate the value of correction."""
    total_w = 0.0
    total_n = 0
    for bs, c in counts.items():
        b = _strip(bs)
        if len(b) < path_len:
            continue
        end_A = _read_outcome(b, 0)
        end_B = _read_outcome(b, path_len - 1)
        raw = 1 - 2 * (end_A ^ end_B)
        total_w += raw * c
        total_n += c
    if total_n == 0:
        return 0.0, float("inf"), 0
    mean = total_w / total_n
    sigma = math.sqrt(max(0.0, 1.0 - mean * mean) / total_n)
    return mean, sigma, total_n


def chsh_corrected_correlator(counts: dict[str, int], path_len: int,
                                ab_setting: tuple[int, int],
                                ideal_table: dict[tuple, dict[tuple[int, int], float]]
                                ) -> tuple[float, float, int]:
    """Byproduct-corrected CHSH correlator E(a,b) from a single setting's counts.

    For each shot's internal-X outcomes m we know the noiseless ideal E_m
    (sign × ~1/√2). The byproduct flips the sign of E_m relative to the
    canonical |Φ+⟩ value E_{|Φ+⟩}(a,b). We correct each shot by

        multiplier(m) = sign(E_m) × sign(E_{|Φ+⟩}(a,b)),

    so the corrected mean recovers E_{|Φ+⟩}(a,b) (with its native sign), and
    S = E00 + E01 + E10 - E11 → 2√2 noiselessly.
    """
    n_internal = max(0, path_len - 2)
    zero_m = tuple([0] * n_internal)
    ref = ideal_table.get(zero_m, {}).get(ab_setting, 0.0)
    target_sign = 1 if ref > 0 else (-1 if ref < 0 else 0)
    if target_sign == 0:
        return 0.0, float("inf"), 0

    total_w = 0.0
    total_n = 0
    for bs, c in counts.items():
        b = _strip(bs)
        if len(b) < path_len:
            continue
        m = tuple(_read_outcome(b, k) for k in range(1, path_len - 1))
        end_A = _read_outcome(b, 0)
        end_B = _read_outcome(b, path_len - 1)
        raw = 1 - 2 * (end_A ^ end_B)
        ideal_m = ideal_table.get(m, {}).get(ab_setting, 0.0)
        if abs(ideal_m) < 1e-6:
            continue
        s_m = 1 if ideal_m > 0 else -1
        # multiplier × raw → s_m · target_sign · raw (signed estimate of E_target)
        total_w += target_sign * s_m * raw * c
        total_n += c
    if total_n == 0:
        return 0.0, float("inf"), 0
    mean = total_w / total_n
    sigma = math.sqrt(max(0.0, 1.0 - mean * mean) / total_n)
    return mean, sigma, total_n


# ---------------------------------------------------------------------------
# Bell fidelity + CHSH
# ---------------------------------------------------------------------------

def bell_fidelity(C_XX: float, C_YY: float, C_ZZ: float,
                   sigma_XX: float, sigma_YY: float, sigma_ZZ: float
                   ) -> dict:
    """Bell fidelity to the *protocol's* Bell-equivalent target state.

    With byproduct-absorbing sign correction (`corrected_correlator` here),
    the conditioned correlators all approach +1 noiselessly regardless of
    which Bell state the protocol's noiseless target is locally equivalent to.
    The Bell fidelity is therefore

        F = (1 + ⟨XX⟩_c + ⟨YY⟩_c + ⟨ZZ⟩_c) / 4 → 1   (noiseless)

    Equivalent (after absorbing sign(⟨O⟩_{|Φ+⟩}) into the multiplier) to
    F = (1 + ⟨XX⟩ - ⟨YY⟩ + ⟨ZZ⟩)/4 against |Φ+⟩.
    """
    F = (1 + C_XX + C_YY + C_ZZ) / 4
    sigma_F = 0.25 * math.sqrt(sigma_XX**2 + sigma_YY**2 + sigma_ZZ**2)
    return {
        "F": F,
        "sigma_F": sigma_F,
        "F_minus_half": F - 0.5,
        "z_above_half": (F - 0.5) / sigma_F if sigma_F > 0 else float("inf"),
        "entangled_3sigma": (F - 3 * sigma_F) > 0.5,
    }


def chsh_S(E00: float, E01: float, E10: float, E11: float,
            s00: float, s01: float, s10: float, s11: float) -> dict:
    """CHSH inequality value from the four corrected correlators.

    With our angle convention (A0=Z, A1=X; B0=(Z+X)/√2, B1=(Z-X)/√2) the
    Tsirelson-saturating combination is

        S = E(A0,B0) + E(A0,B1) + E(A1,B0) - E(A1,B1).

    Classical bound |S| ≤ 2; Tsirelson 2√2 ≈ 2.828.
    """
    S = E00 + E01 + E10 - E11
    sigma_S = math.sqrt(s00**2 + s01**2 + s10**2 + s11**2)
    return {
        "S": S,
        "abs_S": abs(S),
        "sigma_S": sigma_S,
        "p_win": 0.5 + abs(S) / 8,
        "sigma_p": sigma_S / 8,
        "passes_3sigma": (abs(S) - 3 * sigma_S) > 2,
        "z_above_2": (abs(S) - 2) / sigma_S if sigma_S > 0 else float("inf"),
    }


# ---------------------------------------------------------------------------
# High-level runner (simulator + hardware via the same interface)
# ---------------------------------------------------------------------------

def run_routed_bell(
    backend,
    path: list[int],
    *,
    shots: int = 4000,
    do_chsh: bool = False,
    skip_cz: bool = False,
    n_total_qubits: int | None = None,
    optimization_level: int = 3,
    progress: bool = True,
) -> dict:
    """End-to-end runner. Submits one job per measurement setting (3 for Bell,
    +4 for CHSH if requested), batched. Returns a dict with raw counts +
    corrected correlators + Bell fidelity (+ CHSH if requested).
    """
    if n_total_qubits is None:
        n_total_qubits = (max(path) + 1) if hasattr(backend, "num_qubits") else len(path)
        try:
            n_total_qubits = backend.num_qubits
        except Exception:
            pass
    L = len(path)

    # Bell circuits: XX, YY, ZZ
    bell_circuits = []
    for basis in ("XX", "YY", "ZZ"):
        qc, _ = build_path_cluster_circuit(
            path, endpoint_basis=basis, skip_cz=skip_cz,
            n_total_qubits=n_total_qubits,
        )
        bell_circuits.append(qc)

    chsh_circuits = []
    if do_chsh:
        for ab in ((0, 0), (0, 1), (1, 0), (1, 1)):
            qc, _ = build_path_cluster_circuit(
                path, chsh_setting=ab, skip_cz=skip_cz,
                n_total_qubits=n_total_qubits,
            )
            chsh_circuits.append(qc)

    all_circuits = bell_circuits + chsh_circuits

    if progress:
        print(f"[routed_bell] path={path} L={L} shots={shots} "
              f"circuits={len(all_circuits)}{' (with CHSH)' if do_chsh else ''}")

    transpiled = transpile(all_circuits, backend=backend,
                            optimization_level=optimization_level)
    job = backend.run(transpiled, shots=shots)
    try:
        if progress:
            print(f"  job id: {job.job_id()}")
    except Exception:
        pass
    counts_list = job.result().get_counts()
    if isinstance(counts_list, dict):
        counts_list = [counts_list]

    bell_counts = dict(zip(("XX", "YY", "ZZ"), counts_list[:3]))
    chsh_counts = {}
    if do_chsh:
        for ab, c in zip(((0, 0), (0, 1), (1, 0), (1, 1)), counts_list[3:7]):
            chsh_counts[ab] = c

    # Sign tables
    sign_table = derive_sign_table(L, skip_cz=skip_cz)

    out: dict = {
        "path": list(path), "shots": shots, "L": L,
        "skip_cz": skip_cz,
        "counts": {b: dict(bell_counts[b]) for b in bell_counts},
    }

    # Corrected + raw correlators for each Bell setting
    corr = {}
    raw = {}
    for b in ("XX", "YY", "ZZ"):
        c, s, _ = corrected_correlator(bell_counts[b], L, b, sign_table)
        rc, rs, _ = raw_correlator(bell_counts[b], L)
        corr[b] = {"value": c, "sigma": s}
        raw[b] = {"value": rc, "sigma": rs}
    out["corrected_correlators"] = corr
    out["raw_correlators"] = raw

    # Bell fidelity (corrected)
    bf = bell_fidelity(corr["XX"]["value"], corr["YY"]["value"], corr["ZZ"]["value"],
                        corr["XX"]["sigma"], corr["YY"]["sigma"], corr["ZZ"]["sigma"])
    out["bell_fidelity"] = bf

    # Raw (uncorrected) Bell fidelity for control comparison
    raw_bf = bell_fidelity(raw["XX"]["value"], raw["YY"]["value"], raw["ZZ"]["value"],
                            raw["XX"]["sigma"], raw["YY"]["sigma"], raw["ZZ"]["sigma"])
    out["bell_fidelity_uncorrected"] = raw_bf

    if do_chsh:
        chsh_table = derive_chsh_sign_table(L)
        Es = {}
        for ab in ((0, 0), (0, 1), (1, 0), (1, 1)):
            c, s, _ = chsh_corrected_correlator(chsh_counts[ab], L, ab, chsh_table)
            Es[ab] = {"value": c, "sigma": s}
        out["chsh_correlators"] = Es
        out["chsh"] = chsh_S(Es[(0, 0)]["value"], Es[(0, 1)]["value"],
                              Es[(1, 0)]["value"], Es[(1, 1)]["value"],
                              Es[(0, 0)]["sigma"], Es[(0, 1)]["sigma"],
                              Es[(1, 0)]["sigma"], Es[(1, 1)]["sigma"])
        out["counts"].update({f"chsh_{ab}": dict(chsh_counts[ab]) for ab in chsh_counts})

    return out


def print_summary(result: dict) -> None:
    """Pretty-print one routed_bell run."""
    L = result["L"]; path = result["path"]
    bf = result["bell_fidelity"]
    bfu = result["bell_fidelity_uncorrected"]
    print(f"path L={L}: A=q{path[0]} → B=q{path[-1]}  ({path})")
    print(f"  shots={result['shots']}  skip_cz={result['skip_cz']}")
    print("  Corrected:")
    for b in ("XX", "YY", "ZZ"):
        c = result["corrected_correlators"][b]
        print(f"    C_{b} = {c['value']:+.4f} ± {c['sigma']:.4f}")
    print(f"    F = {bf['F']:.4f} ± {bf['sigma_F']:.4f}   "
          f"F-0.5 = {bf['F_minus_half']:+.4f}   "
          f"3σ entangled: {bf['entangled_3sigma']}")
    print("  Uncorrected (control):")
    print(f"    F_uncorr = {bfu['F']:.4f} ± {bfu['sigma_F']:.4f}")
    if "chsh" in result:
        ch = result["chsh"]
        print(f"  CHSH: |S| = {ch['abs_S']:.4f} ± {ch['sigma_S']:.4f}   "
              f"p_win = {ch['p_win']:.4f}   "
              f"3σ pass: {ch['passes_3sigma']}")
