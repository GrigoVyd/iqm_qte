"""Aer-only smoke tests for the cluster-routing demo.

Validates, in this order:
  1. Pure 1D wires of length 2, 3, 4, 5: corrected Bell fidelity ≈ 1,
     uncorrected aggregation ≈ 0 (showing byproduct correction is necessary).
  2. The byproduct rule itself, by binning shots per (sign_x, sign_z) and
     comparing to the analytical prediction.
  3. No-CZ control: corrected Bell fidelity ≤ 0.5, certainly not entangling.
  4. 2D cluster paths (with off-path Z measurements) of length 3 and 5:
     corrected Bell fidelity ≈ 1.
  5. CHSH on the routed endpoints: |S| ≈ 2√2 (Tsirelson).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qiskit_aer import AerSimulator

from src.mbqc import (
    derive_byproduct_table_simulator,
    run_routing_demo,
)


def _check(label: str, cond: bool, info: str = "") -> None:
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {label}  {info}")
    if not cond:
        raise AssertionError(label)


def test_1d_wires_lengths_2_to_5():
    backend = AerSimulator()
    print("=== 1D wires (no off-path qubits): length 2..5 ===")
    for L in (2, 3, 4, 5):
        path = list(range(L))
        res = run_routing_demo(backend, rows=1, cols=L, path=path,
                                shots=8000, do_chsh=False)
        F = res["bell"]["F"]
        c = res["bell_corrected_correlations"]
        u = res["bell_uncorrected_correlations"]
        _check(f"L={L} corrected F ≈ 1", F > 0.97, f"F={F:.4f}")
        # For L >= 3 there is at least one internal X-measurement whose outcome
        # randomizes the byproduct, so the uncorrected aggregation should be
        # near 0 in at least one channel. For L = 2 there are no internal
        # measurements at all so corrected == uncorrected == 1 trivially.
        if L >= 3:
            # Depending on L's parity, one channel may still be deterministic
            # (e.g. for L=3, uncorrected ⟨XX⟩ = +1 because sign_x has no
            # contributions). The byproduct rule is non-trivially needed if
            # AT LEAST ONE channel's uncorrected aggregation is near 0.
            min_uncorr = min(abs(u[k]["E"]) for k in u)
            _check(f"L={L} uncorrected ≪ corrected (some channel ≈ 0)",
                    min_uncorr < 0.1,
                    f"min|uncorrected E|={min_uncorr:.3f}, "
                    f"all={[round(u[k]['E'],3) for k in ('XX','YY','ZZ')]}")
        print(f"      corrected XX,YY,ZZ = "
              f"{c['XX']['E']:+.3f} {c['YY']['E']:+.3f} {c['ZZ']['E']:+.3f}")


def test_byproduct_rule_against_simulator():
    print("\n=== byproduct rule self-consistency on 1D L=4 path ===")
    table = derive_byproduct_table_simulator(rows=1, cols=4,
                                              path=[0, 1, 2, 3], shots=20000)
    # All non-empty bins must match analytical prediction
    for (sx, sz), entry in table.items():
        m = entry["measured_correlations"]
        p = entry["predicted_correlations"]
        for ch in ("XX", "YY", "ZZ"):
            if m[ch] is None:
                continue
            _check(f"bin (sx={sx:+d},sz={sz:+d}) {ch}",
                    abs(m[ch] - p[ch]) < 0.08,
                    f"meas={m[ch]:+.3f}, pred={p[ch]:+d}")


def test_no_cz_control():
    backend = AerSimulator()
    print("\n=== no-CZ control on L=4 wire (entanglement should vanish) ===")
    res = run_routing_demo(backend, rows=1, cols=4, path=[0, 1, 2, 3],
                            shots=4000, do_chsh=False, no_cz=True)
    F = res["bell"]["F"]
    _check("no-CZ F not entangling", F <= 0.5 + 0.05,
            f"F={F:.4f} (should not be near 1)")


def test_2d_cluster_with_offpath():
    backend = AerSimulator()
    print("\n=== 2D cluster routing (off-path Z measurements present) ===")
    # 2x3 grid — path along the top row, off-path qubits below
    path = [0, 1, 2]   # logical-grid indices for the first row
    res = run_routing_demo(backend, rows=2, cols=3, path=path,
                            shots=8000, do_chsh=False)
    F = res["bell"]["F"]
    _check("2x3 top-row L=3 F ≈ 1", F > 0.95, f"F={F:.4f}")

    # 3x3 grid — path snaking through with off-path qubits adjacent
    path = [0, 3, 4, 5, 8]   # 0->3->4->5->8 (down, right, right, down)
    res = run_routing_demo(backend, rows=3, cols=3, path=path,
                            shots=8000, do_chsh=False)
    F = res["bell"]["F"]
    _check("3x3 snake L=5 F ≈ 1", F > 0.95, f"F={F:.4f}")


def test_chsh_on_routed_endpoints():
    backend = AerSimulator()
    print("\n=== CHSH on routed endpoints ===")
    # 1D L=4 wire: should give |S| ≈ 2√2
    res = run_routing_demo(backend, rows=1, cols=4, path=[0, 1, 2, 3],
                            shots=8000, do_chsh=True)
    chsh = res["chsh"]
    _check("L=4 |S| > 2 (classical violation)",
            chsh["abs_S"] > 2.5, f"|S|={chsh['abs_S']:.3f}")
    _check("L=4 |S| approaches 2√2",
            abs(chsh["abs_S"] - 2.0 * 2 ** 0.5) < 0.15,
            f"|S|={chsh['abs_S']:.3f} (target 2√2 ≈ 2.828)")
    print(f"      p_win={chsh['p_win']:.3f} (classical bound 0.75)")

    # 2D 2x3 path with off-path qubits
    res = run_routing_demo(backend, rows=2, cols=3, path=[0, 1, 2],
                            shots=8000, do_chsh=True)
    chsh = res["chsh"]
    _check("2x3 routed |S| > 2", chsh["abs_S"] > 2.5,
            f"|S|={chsh['abs_S']:.3f}")


if __name__ == "__main__":
    test_1d_wires_lengths_2_to_5()
    test_byproduct_rule_against_simulator()
    test_no_cz_control()
    test_2d_cluster_with_offpath()
    test_chsh_on_routed_endpoints()
    print("\n=== ALL CLUSTER-ROUTING TESTS PASSED ===")
