"""CLI: run the Entanglement Bottleneck Map on a backend.

Examples:
    # Mode A — full hardware map (parallel matchings)
    python scripts/run_entanglement_bottleneck_map.py --shots 4000

    # Add a no-CZ control to the same run
    python scripts/run_entanglement_bottleneck_map.py --shots 4000 --no-cz-control

    # Isolated mode (one edge per circuit) and isolated/parallel comparison
    python scripts/run_entanglement_bottleneck_map.py --shots 4000 --isolated --compare

    # Mode B — audit a layout produced by the existing selector
    python scripts/run_entanglement_bottleneck_map.py --shots 4000 \
        --audit-grid 2x3

Reads IQM_TOKEN from environment. Falls back to AerSimulator if absent (in that
case --edges must be provided since the simulator has no coupling map).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.backend import get_backend, is_simulator, select_best_subgrid
from src.diagnostics import (
    audit_layout,
    edge_residuals,
    get_coupling_edges,
    isolated_vs_parallel,
    node_diagnostics,
    plot_edge_map,
    run_edge_map,
    save_results_json_csv,
)


def _parse_grid(s: str) -> tuple[int, int]:
    r, c = s.lower().split("x")
    return int(r), int(c)


def _parse_edges(s: str) -> list[tuple[int, int]]:
    out = []
    for part in s.split(","):
        a, b = part.split(":")
        out.append((int(a), int(b)))
    return out


def _grid_logical_edges(rows: int, cols: int) -> list[tuple[int, int]]:
    es = []
    for r in range(rows):
        for c in range(cols - 1):
            es.append((r * cols + c, r * cols + c + 1))
    for r in range(rows - 1):
        for c in range(cols):
            es.append((r * cols + c, (r + 1) * cols + c))
    return es


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--device", default="emerald")
    p.add_argument("--shots", type=int, default=4000)
    p.add_argument("--isolated", action="store_true",
                   help="run isolated mode (one edge per circuit)")
    p.add_argument("--compare", action="store_true",
                   help="run BOTH parallel and isolated, then compare")
    p.add_argument("--no-cz-control", action="store_true",
                   help="also run the no-CZ control on the same edges")
    p.add_argument("--edges",
                   help="explicit edge list as 'a:b,c:d,...' (required on Aer)")
    p.add_argument("--audit-grid", help="audit best rows x cols grid layout, "
                                          "e.g. 2x3")
    p.add_argument("--out", default="results/edge_map",
                   help="output prefix for .json/.csv (no extension)")
    p.add_argument("--plot", default=None, help="save edge-map PNG to this path")
    args = p.parse_args()

    token = os.environ.get("IQM_TOKEN")
    backend = get_backend(token=token, device=args.device)
    print(f"Backend: {backend}  hardware={not is_simulator(backend)}")

    if args.edges:
        edges = _parse_edges(args.edges)
    elif is_simulator(backend):
        sys.exit("On AerSimulator, --edges is required (e.g. '0:1,1:2,2:3').")
    else:
        edges = get_coupling_edges(backend)
    print(f"Edges to test: {len(edges)}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    if args.compare:
        print("\n--- parallel mode ---")
        par = run_edge_map(backend, edges=edges, shots=args.shots, mode="parallel")
        print("\n--- isolated mode ---")
        iso = run_edge_map(backend, edges=edges, shots=args.shots, mode="isolated")
        save_results_json_csv(par, args.out + "_parallel")
        save_results_json_csv(iso, args.out + "_isolated")
        cmp = isolated_vs_parallel(iso, par)
        with open(args.out + "_compare.json", "w") as f:
            json.dump(cmp, f, indent=2)
        primary = par
    else:
        mode = "isolated" if args.isolated else "parallel"
        primary = run_edge_map(backend, edges=edges, shots=args.shots, mode=mode)
        save_results_json_csv(primary, args.out)

    if args.no_cz_control:
        print("\n--- no-CZ control (same edges) ---")
        ctrl = run_edge_map(backend, edges=edges, shots=args.shots,
                            mode="parallel", no_cz=True)
        save_results_json_csv(ctrl, args.out + "_nocz")
        bad = [r for r in ctrl if r.F > 0.5 + 3 * r.sigma_F]
        print(f"  no-CZ edges falsely above 0.5+3σ: {len(bad)} "
              f"(should be ≈0; non-zero indicates a bug)")

    elapsed = time.time() - t0
    print(f"\nMap done in {elapsed:.1f} s. {len(primary)} edges measured.")
    n_ent = sum(1 for r in primary if r.entangled_3sigma)
    print(f"  certified entangling (3σ): {n_ent}/{len(primary)}")
    if primary:
        Fs = [r.F for r in primary]
        print(f"  F: min={min(Fs):.3f}  mean={sum(Fs)/len(Fs):.3f}  "
              f"max={max(Fs):.3f}")
    if not is_simulator(backend):
        worst = sorted(primary, key=lambda r: r.F)[:5]
        print("  weakest 5 edges:")
        for r in worst:
            print(f"    {r.edge}: F={r.F:.3f}  z={r.z_score:.1f}  "
                  f"3σ-entangled={r.entangled_3sigma}")

    nodes = node_diagnostics(primary)
    if nodes:
        worst_q = sorted(nodes.items(), key=lambda kv: kv[1]["Q"])[:5]
        print("  weakest 5 qubits (mean F over incident edges):")
        for q, info in worst_q:
            print(f"    Q{q}: Q={info['Q']:.3f}  ({info['n_edges']} edges)")

    resid = edge_residuals(primary)[:5]
    if resid:
        print("  most coupler-specific weak edges (most negative residual):")
        for row in resid:
            print(f"    {row['edge']}: F={row['F']:.3f}  "
                  f"baseline={row['baseline']:.3f}  R={row['residual']:+.3f}")

    if args.audit_grid:
        rows, cols = _parse_grid(args.audit_grid)
        if is_simulator(backend):
            layout = list(range(rows * cols))
            print("(simulator) using trivial layout for audit demo")
        else:
            layout, cost, n_cand = select_best_subgrid(backend, rows, cols)
            print(f"\nselector picked layout {layout} (cost={cost:.4f}, "
                  f"candidates={n_cand})")
        req = _grid_logical_edges(rows, cols)
        audit = audit_layout(layout, req, primary)
        with open(args.out + "_audit.json", "w") as f:
            json.dump(audit, f, indent=2)
        print(f"  validated={audit.get('validated')}  "
              f"weakest_edge={audit.get('weakest_edge')}  "
              f"failing_edges={len(audit.get('failing_edges', []))}")

    if args.plot:
        plot_edge_map(primary, backend=backend, save_path=args.plot,
                      title="Entanglement bottleneck map")
        print(f"  plot: {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
