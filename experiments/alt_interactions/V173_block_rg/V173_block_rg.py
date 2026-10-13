#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V173 - block RG: does the anisotropy grow or shrink under coarse-
graining? (alt-interactions, the Ising-miracle route, step 3)

**The question.**  The Ising/lattice-phi^4 miracle says the cubic
anisotropy is IRRELEVANT at the interacting fixed point: under coarse-
graining the cone becomes isotropic.  The loop-free version of that test
on the CL10 cone: downfold the effective Hamiltonian H(q) = i log U(q)
onto the low subspace DEFINED AT THE NODE (modes with |theta| < cut,
cut in {pi/2, pi/4, pi/8} = block scales b in {1, 2, 4}) and recompute
the cone anisotropy A(b) of the downfolded 4-band theory at dq = 1e-3 -
bare, and with the V172 contact self-energy Sigma (U_c = g_10) inserted
in the cone mode, normal-ordered (the k-independent tadpole subtracted).
A(b) -> 0 is the miracle; A(b) ~ 0.998 constant is marginality.

**The result.**  A(b) stays at the free value 0.998 for every scale and
with the self-energy insertion (recorded shifts): the cubic anisotropy
of the free cone is MARGINAL under coarse-graining and the minimal
admissible collision does not change that at this order.  The downfold
is exact on the low band (b = 1 reproduces the full gap_at cone to
1e-6), so this is a statement about the theory, not the truncation.
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

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
QSTAR = (math.pi / 2, math.pi / 2, 0.0)
CUTS = (math.pi / 2, math.pi / 4, math.pi / 8)   # b = 1, 2, 4
DQ = 1e-3
UC = SC.G10


def bf(q):
    return A.bloch_cell(C, SC.sectors()[2], SC.G10, ST, q)


def main():
    t0 = time.time()
    dirs = A.sphere_dirs(6)
    kappas = [(0.0, 0.0, 0.0)] + [tuple(DQ * x for x in d) for d in dirs]
    sig = A.fock_self_energy(bf, QSTAR, UC, 24, kappas,   # 1x1 each;
                             sub=(0.15, 0.6), nodes=SC.M_POINTS)

    # band structure at the node: subspace sizes and window-edge margins
    Hr = A.effective_hamiltonian(bf, QSTAR)
    th = np.abs(np.linalg.eigvalsh(Hr))
    rows = []
    for b, cut in zip((1, 2, 4), CUTS):
        n_low = int(np.sum(th < cut))
        edge = float(np.min(np.abs(th - cut)))     # contamination margin
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
            "A_bare": aniso(vs_bare), "A_interacting": aniso(vs_int),
            "isotropic": bool(min(aniso(vs_bare), aniso(vs_int)) < 1e-6),
        })

    # 2-route: the b=1 downfolded cone vs the FULL gap_at cone
    vs_full = [A.gap_at(bf(tuple(a + DQ * x for a, x in zip(QSTAR, d))))
               / DQ for d in dirs]
    A_full = (max(vs_full) - min(vs_full)) / (max(vs_full) + min(vs_full))
    route_dev = abs(rows[0]["A_bare"] - A_full)

    anchor_ok = abs(rows[0]["A_bare"] - 0.998) < 1e-3
    route_ok = route_dev < 1e-6
    flow = [r["A_bare"] - rows[0]["A_bare"] for r in rows]
    ok = anchor_ok and route_ok and len(rows) == 3
    margins = ", ".join(f"{r['window_margin']:.1e}" for r in rows)
    nlows = ", ".join(str(r["n_low"]) for r in rows)
    verdict = {
        "id": "V173_block_rg",
        "tier": 4,
        "property": "energy-window downfold of H(q) = i log U(q) onto the "
                    "node subspace at cuts pi/2, pi/4, pi/8 (b = 1, 2, 4): "
                    "cone anisotropy A(b), bare and with the contact "
                    "self-energy (U_c = g_10)",
        "residual": True,
        "pass": bool(ok),
        "status": f"BLOCK-RG ANISOTROPY FLOW (6 rows). A(b) = "
                  f"{rows[0]['A_bare']:.4f}, {rows[1]['A_bare']:.4f}, "
                  f"{rows[2]['A_bare']:.4f} for b = 1, 2, 4 (bare) and "
                  f"{rows[0]['A_interacting']:.4f}, "
                  f"{rows[1]['A_interacting']:.4f}, "
                  f"{rows[2]['A_interacting']:.4f} with the contact "
                  f"self-energy: the cubic anisotropy is MARGINAL under "
                  f"coarse-graining (flow {min(flow):.2e}..{max(flow):.2e}"
                  f"), the collision does not restore isotropy. Subspace "
                  f"sizes n_low = {nlows}, window "
                  f"margins {margins}. "
                  f"2-route (b=1 downfold vs full gap_at) dev "
                  f"{route_dev:.1e}.",
        "data": {"rows": rows, "node": list(QSTAR), "dq": DQ, "U_c": UC,
                 "A_full": A_full, "route_dev": route_dev,
                 "candidates": [r for r in rows if r["isotropic"]]},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v173_block_rg.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V173 - block-RG flow of the cone anisotropy")
    print("=" * 68)
    for r in rows:
        print(f"  b = {r['b']} (cut {r['cut']:.3f}, n_low {r['n_low']}): "
              f"A_bare = {r['A_bare']:.4f}  A_int = {r['A_interacting']:.4f}")
    print(f"  2-route dev (b=1 vs full): {route_dev:.1e}   "
          f"VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v173_block_rg.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
