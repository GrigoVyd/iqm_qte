"""Two-panel plot for the cluster-routing scaling sweep.

Reads results/cluster_routing_garnet_sweep.json (or any sweep JSON via --in)
and produces:
  - left  panel: Bell fidelity F vs path length L (with 3σ error bars and
                 the 0.5 separability bound + 3σ-cert threshold marked)
  - right panel: CHSH |S| vs L (with classical bound 2 and Tsirelson 2√2)

Both panels share the same x axis. Routes that pass the witness are colored
green, those that fail amber.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="in_path",
                   default="results/cluster_routing_garnet_sweep.json")
    p.add_argument("--out", default="results/cluster_routing_scaling.png")
    p.add_argument("--device", default="Garnet")
    args = p.parse_args()

    src = Path(args.in_path)
    if not src.exists():
        sys.exit(f"missing {src}")
    rows = json.loads(src.read_text())

    # Drop no-CZ controls if any; keep only data routes
    rows = [r for r in rows if not r.get("no_cz", False)]
    rows.sort(key=lambda r: r["L"])

    Ls = np.array([r["L"] for r in rows])
    Fs = np.array([r["bell"]["F"] for r in rows])
    sFs = np.array([r["bell"]["sigma_F"] for r in rows])
    pass_F = np.array([r["bell"]["entangled_3sigma"] for r in rows])

    have_chsh = all("chsh" in r for r in rows)
    if have_chsh:
        Ss = np.array([r["chsh"]["abs_S"] for r in rows])
        sSs = np.array([r["chsh"]["sigma_S"] for r in rows])
        pass_S = np.array([r["chsh"]["pass_3sigma"] for r in rows])

    fig, axes = plt.subplots(1, 2 if have_chsh else 1, figsize=(13, 5))
    ax_F = axes[0] if have_chsh else axes
    ax_S = axes[1] if have_chsh else None

    # Bell-fidelity panel
    colors_F = ["#2ca02c" if p else "#ff7f0e" for p in pass_F]
    ax_F.errorbar(Ls, Fs, yerr=3 * sFs, fmt="o", capsize=4,
                   markersize=10, ecolor="0.3", elinewidth=1.2,
                   markerfacecolor="white", zorder=3)
    for L, F, sF, c in zip(Ls, Fs, sFs, colors_F):
        ax_F.plot([L], [F], "o", markersize=10, color=c, zorder=4)
    ax_F.axhline(0.5, color="black", linestyle="--", linewidth=1,
                  label="separability bound (F = 0.5)")
    ax_F.axhline(1.0, color="0.7", linestyle=":", linewidth=1,
                  label="ideal F = 1")
    ax_F.set_xticks(Ls)
    ax_F.set_xlabel("path length L  (qubits in route)")
    ax_F.set_ylabel("Bell fidelity F  (3σ error bars)")
    ax_F.axhline(0.25, color="0.7", linestyle=":", linewidth=1,
                  label="F = 1/4 (separable / random)")
    ax_F.set_ylim(0.0, 1.05)
    ax_F.set_title(
        f"Routed endpoint Bell fidelity — {args.device}\n"
        f"endpoint pair (path[0], path[-1]) verified via byproduct-corrected XX/YY/ZZ"
    )
    ax_F.legend(loc="lower left", fontsize=9)
    ax_F.grid(alpha=0.25)
    for L, F in zip(Ls, Fs):
        ax_F.annotate(f"{F:.3f}", (L, F), xytext=(8, 0),
                       textcoords="offset points", fontsize=9,
                       va="center")

    if have_chsh:
        colors_S = ["#2ca02c" if p else "#ff7f0e" for p in pass_S]
        ax_S.errorbar(Ls, Ss, yerr=3 * sSs, fmt="o", capsize=4,
                       markersize=10, ecolor="0.3", elinewidth=1.2,
                       markerfacecolor="white", zorder=3)
        for L, S, sS, c in zip(Ls, Ss, sSs, colors_S):
            ax_S.plot([L], [S], "o", markersize=10, color=c, zorder=4)
        ax_S.axhline(2.0, color="black", linestyle="--", linewidth=1,
                      label="classical bound (|S| = 2)")
        ax_S.axhline(2.0 * math.sqrt(2), color="0.7", linestyle=":",
                      linewidth=1, label=f"Tsirelson 2√2 ≈ 2.828")
        ax_S.set_xticks(Ls)
        ax_S.set_xlabel("path length L  (qubits in route)")
        ax_S.set_ylabel("CHSH |S|  (3σ error bars)")
        ax_S.set_ylim(0, 3.0)
        ax_S.set_title(
            f"CHSH game on routed endpoints — {args.device}\n"
            f"|S| > 2 ⇒ classical bound violated; p_win > 0.75"
        )
        ax_S.legend(loc="lower left", fontsize=9)
        ax_S.grid(alpha=0.25)
        for L, S in zip(Ls, Ss):
            ax_S.annotate(f"{S:.3f}", (L, S), xytext=(8, 0),
                           textcoords="offset points", fontsize=9,
                           va="center")

    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
