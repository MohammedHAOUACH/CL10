#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V178 - block RG with the gauged vertex insertion (alt-gauged, step 4)

**The question.**  V173: with the contact self-energy the cone anisotropy
A(b) stays at 0.9980 for b = 1, 2, 4 - marginal under coarse-graining.
V176 showed the gauged bond-current vertex carries the Eg + T2g (spin-2)
pieces the contact lacked.  Does the GAUGED self-energy (V177 vertex,
U_c = g_10, grid 16, normal-ordered) change the block-RG flow?  Same
loop-free protocol as V173: downfold H(q) = i log U(q) onto the node
subspace at cuts pi/2, pi/4, pi/8 (b = 1, 2, 4), recompute A(b) at
dq = 1e-3, bare and with the gauged Sigma in the cone mode.

**The result.**  Recorded: A(b) bare and gauged-dressed; the b = 1 bare
row must reproduce A = 0.998 (anchor).
"""

import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
os.environ["CL10_ANGLE"] = "0.2"  # the W9 cone and its anisotropy anchors are 1/5-branch (V115/V123) numbers: pin the clock, the closed default_branch angle must not leak in

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
QSTAR = (math.pi / 2, math.pi / 2, 0.0)
CUTS = (math.pi / 2, math.pi / 4, math.pi / 8)   # b = 1, 2, 4
DQ = 1e-3
UC = SC.G10
GRID = 16


def bf(q):
    return A.bloch_cell(C, SC.sectors()[2], SC.G10, ST, q)


def main():
    t0 = time.time()
    print("=" * 68)
    print("V178 - block RG with the gauged vertex insertion")
    print("=" * 68)
    dirs = A.sphere_dirs(6)
    kappas = [(0.0, 0.0, 0.0)] + [tuple(DQ * x for x in d) for d in dirs]
    sig = A.fock_self_energy(bf, QSTAR, UC, GRID, kappas,
                             sub=(0.15, 0.6), nodes=SC.M_POINTS,
                             current_vertex=True)
    print(f"  gauged self-energy (grid {GRID}) done "
          f"[{time.time() - t0:.0f}s]")

    Hr = A.effective_hamiltonian(bf, QSTAR)
    th = np.abs(np.linalg.eigvalsh(Hr))
    rows = []
    for b, cut in zip((1, 2, 4), CUTS):
        n_low = int(np.sum(th < cut))
        edge = float(np.min(np.abs(th - cut)))
        vs_bare, vs_int = [], []
        for i, d in enumerate(dirs):
            q = tuple(a + DQ * x for a, x in zip(QSTAR, d))
            Heff, ev = A.lowdin_downfold(bf, q, QSTAR, cut=cut)
            vs_bare.append(float(np.min(np.abs(ev))) / DQ)
            s = float(np.real((sig[1 + i] - sig[0])[0, 0]))
            ev2 = np.linalg.eigvalsh(
                Heff + np.diag([s] + [0.0] * (Heff.shape[0] - 1)))
            vs_int.append(float(np.min(np.abs(ev2))) / DQ)

        def aniso(vs):
            return (max(vs) - min(vs)) / (max(vs) + min(vs))
        rows.append({
            "b": b, "cut": cut, "n_low": n_low, "window_margin": edge,
            "A_bare": aniso(vs_bare), "A_gauged": aniso(vs_int),
            "isotropic": bool(min(aniso(vs_bare), aniso(vs_int)) < 1e-6),
        })
        print(f"  b = {b} (cut {cut:.3f}, n_low {n_low}): A_bare = "
              f"{rows[-1]['A_bare']:.4f}  A_gauged = {rows[-1]['A_gauged']:.4f}")

    # 2-route: the b=1 downfolded cone vs the FULL gap_at cone
    vs_full = [A.gap_at(bf(tuple(a + DQ * x for a, x in zip(QSTAR, d))))
               / DQ for d in dirs]
    A_full = (max(vs_full) - min(vs_full)) / (max(vs_full) + min(vs_full))
    route_dev = abs(rows[0]["A_bare"] - A_full)

    anchor_ok = abs(rows[0]["A_bare"] - 0.998) < 1e-3
    route_ok = route_dev < 1e-6
    flow = [r["A_gauged"] - rows[0]["A_bare"] for r in rows]
    ok = anchor_ok and route_ok and len(rows) == 3

    margins = ", ".join(f"{r['window_margin']:.1e}" for r in rows)
    nlows = ", ".join(str(r["n_low"]) for r in rows)
    verdict = {
        "id": "V178_gauged_block_rg",
        "tier": 4,
        "property": "energy-window downfold of H(q) = i log U(q) onto the "
                    "node subspace at cuts pi/2, pi/4, pi/8 (b = 1, 2, 4): "
                    "cone anisotropy A(b), bare and with the GAUGED "
                    "bond-current self-energy (U_c = g_10, grid 16)",
        "residual": True,
        "pass": bool(ok),
        "status": f"GAUGED BLOCK-RG FLOW (6 rows). A(b) = "
                  f"{rows[0]['A_bare']:.4f}, {rows[1]['A_bare']:.4f}, "
                  f"{rows[2]['A_bare']:.4f} (bare) and "
                  f"{rows[0]['A_gauged']:.4f}, {rows[1]['A_gauged']:.4f}, "
                  f"{rows[2]['A_gauged']:.4f} with the gauged "
                  f"self-energy: the cubic anisotropy is "
                  f"{'MARGINAL' if max(abs(x) for x in flow) < 0.05 else 'FLOWING'} "
                  f"under coarse-graining (gauged-bare shifts "
                  f"{min(flow):.2e}..{max(flow):.2e}); the gauged vertex "
                  f"{'does NOT restore' if min(r['A_gauged'] for r in rows) > 0.5 else 'restores'} "
                  f"isotropy at this order. Subspace sizes n_low = {nlows}, "
                  f"window margins {margins}. 2-route (b=1 downfold vs "
                  f"full gap_at) dev {route_dev:.1e}.",
        "data": {"rows": rows, "node": list(QSTAR), "dq": DQ, "U_c": UC,
                 "grid": GRID, "A_full": A_full, "route_dev": route_dev,
                 "candidates": [r for r in rows if r["isotropic"]]},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v178_gauged_block_rg.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  2-route dev (b=1 vs full): {route_dev:.1e}   isotropic rows: "
          f"{sum(r['isotropic'] for r in rows)}   "
          f"VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v178_gauged_block_rg.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
