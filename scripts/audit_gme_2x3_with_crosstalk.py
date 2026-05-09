"""Targeted hardware test for the Entanglement Bottleneck Map.

Audits the SAME 2x3 layout that select_best_subgrid picks for the GME demo,
in BOTH parallel and isolated modes. Three results from one job:

  1. Audit       — every internal edge of the layout, F + 3σ-classification
  2. Bottleneck  — the weakest edge (and whether it explains the W deficit)
  3. Crosstalk   — ΔF = F_isolated − F_parallel per edge

Run:
    python scripts/audit_gme_2x3_with_crosstalk.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.backend import get_backend, is_simulator, select_best_subgrid
from src.diagnostics import (
    audit_layout, isolated_vs_parallel, node_diagnostics, run_edge_map,
)

ROWS, COLS = 2, 3
SHOTS = 4000
DEVICE = "emerald"
OUT_PREFIX = "results/audit_gme_2x3"


def grid_logical_edges(rows: int, cols: int) -> list[tuple[int, int]]:
    es = []
    for r in range(rows):
        for c in range(cols - 1):
            es.append((r * cols + c, r * cols + c + 1))
    for r in range(rows - 1):
        for c in range(cols):
            es.append((r * cols + c, (r + 1) * cols + c))
    return es


def main() -> int:
    Path(OUT_PREFIX).parent.mkdir(parents=True, exist_ok=True)
    print(f"=== Audit + crosstalk on best {ROWS}x{COLS} layout ({DEVICE}) ===\n")

    backend = get_backend(token=os.environ.get("IQM_TOKEN"), device=DEVICE)
    print(f"Backend: {backend}  hardware={not is_simulator(backend)}\n")

    if is_simulator(backend):
        layout = list(range(ROWS * COLS))
        print(f"(simulator) trivial layout: {layout}")
    else:
        layout, cost, n_cand = select_best_subgrid(backend, ROWS, COLS)
        print(f"Selector layout: {layout}")
        print(f"  cost={cost:.4f}  ({n_cand} candidates considered)\n")

    logical_edges = grid_logical_edges(ROWS, COLS)
    physical_edges = sorted({
        (min(layout[i], layout[j]), max(layout[i], layout[j]))
        for i, j in logical_edges
    })
    print(f"Internal edges to test ({len(physical_edges)}):")
    for e in physical_edges:
        print(f"  {e}")
    print()

    t0 = time.time()
    print(f"--- Parallel mode ({SHOTS} shots) ---")
    par = run_edge_map(backend, edges=physical_edges, shots=SHOTS,
                       mode="parallel", progress=True)
    print(f"\n--- Isolated mode ({SHOTS} shots) ---")
    iso = run_edge_map(backend, edges=physical_edges, shots=SHOTS,
                       mode="isolated", progress=True)
    elapsed = time.time() - t0
    print(f"\nBoth modes done in {elapsed:.1f} s.\n")

    # 1. Audit
    audit = audit_layout(layout, logical_edges, par)
    print("=" * 60)
    print(f"1. AUDIT (parallel mode, the one the GME demo actually uses)")
    print("=" * 60)
    print(f"  validated (all edges 3σ-entangling): {audit['validated']}")
    print(f"  mean F = {audit['mean_F']:.3f}    "
          f"min F = {audit['min_F']:.3f}    "
          f"failing (3σ): {audit['n_below_3sigma']}/{audit['n_edges']}")
    print(f"  weakest edge: {audit['weakest_edge']}")
    if audit["failing_edges"]:
        print("  edges below 3σ threshold:")
        for fe in audit["failing_edges"]:
            print(f"    {fe}")

    # 2. Per-edge breakdown
    print()
    print("=" * 60)
    print("2. PER-EDGE BREAKDOWN")
    print("=" * 60)
    print(f"  {'edge':<14}{'F_par':>8}{'F_iso':>8}{'ΔF':>9}"
          f"{'z_iso':>8}{'3σ-par':>9}")
    par_map = {r.edge: r for r in par}
    iso_map = {r.edge: r for r in iso}
    for e in physical_edges:
        rp, ri = par_map[e], iso_map[e]
        delta = ri.F - rp.F
        flag = "yes" if rp.entangled_3sigma else "NO"
        print(f"  {str(e):<14}{rp.F:>8.3f}{ri.F:>8.3f}{delta:>+9.3f}"
              f"{ri.z_score:>8.1f}{flag:>9}")

    # 3. Crosstalk summary
    cmp = isolated_vs_parallel(iso, par)
    print()
    print("=" * 60)
    print("3. CROSSTALK (isolated − parallel)")
    print("=" * 60)
    deltas = [c["delta_F"] for c in cmp]
    z_deltas = [c["z_delta"] for c in cmp]
    print(f"  mean ΔF = {sum(deltas)/len(deltas):+.3f}")
    print(f"  max  ΔF = {max(deltas):+.3f}")
    sig = [c for c in cmp if abs(c["z_delta"]) > 3]
    print(f"  edges with |z_ΔF| > 3 (significant crosstalk): {len(sig)}/{len(cmp)}")
    if sig:
        for c in sorted(sig, key=lambda x: x["delta_F"], reverse=True):
            print(f"    {c['edge']}: ΔF={c['delta_F']:+.3f}  "
                  f"z={c['z_delta']:+.1f}")

    # Node-level
    print()
    print("=" * 60)
    print("4. NODE-LEVEL (mean F over incident edges, parallel mode)")
    print("=" * 60)
    nodes = node_diagnostics(par)
    for q in sorted(nodes):
        info = nodes[q]
        print(f"  Q{q}: Q={info['Q']:.3f}  ({info['n_edges']} edges)")

    # Save everything
    out = {
        "device": DEVICE, "rows": ROWS, "cols": COLS, "shots": SHOTS,
        "layout": list(layout),
        "physical_edges": [list(e) for e in physical_edges],
        "audit_parallel": audit,
        "edges_parallel": [r.to_dict() for r in par],
        "edges_isolated": [r.to_dict() for r in iso],
        "crosstalk": cmp,
        "node_diagnostics_parallel": {str(k): v for k, v in nodes.items()},
    }
    with open(OUT_PREFIX + ".json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUT_PREFIX}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
