"""CLI: Cluster Entanglement Routing + CHSH Game.

Builds a 2D cluster on a rows×cols grid, then routes entanglement onto two
distant endpoint qubits A and B by measuring all off-path qubits in Z and all
internal-path qubits in X. The endpoint pair is verified with a Bell-fidelity
witness (F − 3σ_F > 0.5) and a CHSH game (|S| − 3σ_S > 2).

Examples:
    # Aer simulator demo, default 1D wires of several lengths:
    python scripts/run_cluster_routing_demo.py --rows 1 --cols 5

    # Single explicit path on 2x3 grid:
    python scripts/run_cluster_routing_demo.py --rows 2 --cols 3 --path 0,1,2

    # Sweep multiple paths on 3x3 with no-CZ control:
    python scripts/run_cluster_routing_demo.py --rows 3 --cols 3 \
        --path 0,1,2 --path 0,3,4,5,8 --no-cz-control

    # Hardware run (uses select_best_subgrid for the layout):
    python scripts/run_cluster_routing_demo.py --rows 2 --cols 3 \
        --path 0,1,2 --shots 4000

Reads IQM_TOKEN from environment. Falls back to AerSimulator if absent.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.backend import get_backend, is_simulator, select_best_subgrid
from src.mbqc import (
    run_routing_demo, save_results_json, validate_path_is_connected,
)


def parse_path(s: str) -> list[int]:
    return [int(x) for x in s.split(",")]


def print_result(res: dict) -> None:
    bell = res["bell"]
    print(f"  route       : {res['path']}  (L={res['L']})")
    print(f"  endpoints   : (A, B) = ({res['A']}, {res['B']})")
    print(f"  shots       : {res['shots']}    no_cz={res['no_cz']}")
    cc = res["bell_corrected_correlations"]
    print(f"  C_XX,YY,ZZ  : {cc['XX']['E']:+.4f}  {cc['YY']['E']:+.4f}  "
          f"{cc['ZZ']['E']:+.4f}    (corrected, all should be ≈ +1)")
    uc = res["bell_uncorrected_correlations"]
    print(f"  uncorrected : {uc['XX']['E']:+.4f}  {uc['YY']['E']:+.4f}  "
          f"{uc['ZZ']['E']:+.4f}    (control: shows byproduct correction works)")
    print(f"  F           : {bell['F']:.4f}    σ_F = {bell['sigma_F']:.4f}")
    print(f"  F − 0.5     : {bell['F'] - 0.5:+.4f}    z = {bell['z_score']:.2f}")
    print(f"  entangled (3σ): {'YES' if bell['entangled_3sigma'] else 'no'}    "
          f"entangled (mean): {'YES' if bell['entangled_meanonly'] else 'no'}")
    if "chsh" in res:
        ch = res["chsh"]
        print(f"  CHSH |S|    : {ch['abs_S']:.4f}    σ_S = {ch['sigma_S']:.4f}    "
              f"(classical ≤ 2, Tsirelson 2√2 ≈ 2.828)")
        print(f"  p_win       : {ch['p_win']:.4f}    "
              f"(classical ≤ 0.75; CHSH pass: "
              f"{'YES' if ch['pass_3sigma'] else 'no'})")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--device", default="emerald")
    p.add_argument("--rows", type=int, default=1)
    p.add_argument("--cols", type=int, default=5)
    p.add_argument("--path", action="append",
                   help="comma-separated qubit indices, can be repeated")
    p.add_argument("--shots", type=int, default=4000)
    p.add_argument("--no-chsh", action="store_true",
                   help="skip CHSH circuits (Bell-fidelity only)")
    p.add_argument("--no-cz-control", action="store_true",
                   help="also run a no-CZ control on each path")
    p.add_argument("--out", default="results/cluster_routing.json")
    p.add_argument("--auto-layout", action="store_true",
                   help="on hardware, pick the best rows×cols sub-grid via "
                        "select_best_subgrid and remap the path")
    args = p.parse_args()

    backend = get_backend(token=os.environ.get("IQM_TOKEN"), device=args.device)
    print(f"Backend: {backend}    hardware={not is_simulator(backend)}")

    paths = [parse_path(s) for s in (args.path or [])]
    if not paths:
        # Default sweep: progressively longer paths along the first row
        paths = [list(range(L)) for L in (2, 3, 4, 5)
                 if L <= args.rows * args.cols]
        print(f"(no --path given) sweeping default 1D wires: {paths}")

    for path in paths:
        validate_path_is_connected(path, args.rows, args.cols)

    initial_layout = None
    if args.auto_layout and not is_simulator(backend):
        layout, cost, n_cand = select_best_subgrid(backend, args.rows,
                                                     args.cols)
        initial_layout = layout
        print(f"\nselector layout: {layout}  (cost={cost:.4f}, "
              f"candidates={n_cand})")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    all_results: list[dict] = []

    t0 = time.time()
    for i, path in enumerate(paths):
        label = f"L={len(path)}_{'-'.join(map(str, path))}"
        print(f"\n--- route {i+1}/{len(paths)}: {label} ---")
        res = run_routing_demo(
            backend, args.rows, args.cols, path, shots=args.shots,
            initial_layout=initial_layout, do_chsh=not args.no_chsh,
            no_cz=False, label=label,
        )
        print_result(res)
        all_results.append(res)

        if args.no_cz_control:
            print(f"  ... no-CZ control:")
            ctrl = run_routing_demo(
                backend, args.rows, args.cols, path, shots=args.shots,
                initial_layout=initial_layout, do_chsh=not args.no_chsh,
                no_cz=True, label=f"{label}_noCZ",
            )
            ctrl_F = ctrl["bell"]["F"]
            print(f"      no-CZ F = {ctrl_F:.4f}    "
                  f"(must NOT exceed 0.5 + 3σ to validate the witness)")
            if "chsh" in ctrl:
                print(f"      no-CZ |S| = {ctrl['chsh']['abs_S']:.4f}    "
                      f"(must NOT exceed 2)")
            all_results.append(ctrl)

    elapsed = time.time() - t0
    print(f"\n{len(paths)} routes done in {elapsed:.1f} s.")

    save_results_json(all_results, args.out)
    print(f"Saved: {args.out}")

    # Summary table
    print("\n=== SUMMARY ===")
    print(f"  {'label':<32} {'L':>3} {'F':>8} {'σ_F':>8} "
          f"{'F-0.5':>8} {'3σ?':>5} {'|S|':>8} {'CHSH?':>6}")
    for r in all_results:
        b = r["bell"]
        s = r.get("chsh", {})
        s_val = f"{s['abs_S']:.3f}" if s else "—"
        s_pass = ("YES" if s.get("pass_3sigma") else "no") if s else "—"
        print(f"  {r['label']:<32} {r['L']:>3} {b['F']:>8.4f} "
              f"{b['sigma_F']:>8.4f} {b['F'] - 0.5:>+8.4f} "
              f"{'YES' if b['entangled_3sigma'] else 'no':>5} "
              f"{s_val:>8} {s_pass:>6}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
