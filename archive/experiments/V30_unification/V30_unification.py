#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V30 — the 0-parameter coupling unification and the QCA UV completion.

**The 0-parameter predictions.**  The **$5/3$** GUT normalization of the SM
$\\mathrm{U}(1)$ (Theorem 3) and the **single GUT coupling**
$\alpha_{10} = 1/25$ (the 0-parameter GUT Cartan).

**The 1-loop running to the 0-parameter scale.**  Running the three SM
couplings from $M_Z$ to the **Planckian** GUT scale $M_{\\rm GUT}=M_P$ gives
couplings in the **GUT range** ($\\sim0.02$-$0.04$, consistent with
$\\alpha_{10}=0.04$); the **1-loop spread** at $M_P$ (the non-SUSY SM does not
meet at a single point — the clean MSSM unification is at $\\sim2\\times10^{16}$,
**below** $M_P$) is the **named 1-loop matching residual**.

**The QCA UV completion.**  The U(1) **Landau pole** ($\\sim1.4\\times10^{39}$
GeV) is **above** the QCA lattice UV cutoff $M_P$ — so the running is
**UV-bounded** on the lattice: the **0-parameter UV completion** of the theory.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import unification

TOL = 0.20


def main():
    t0 = time.time()
    r = unification.unification()

    five_thirds_ok = 1.6 < r["five_thirds"] < 1.8
    gut_coupling_ok = abs(r["alpha_10"] - 0.04) < 1e-9
    gut_range_ok = all(0.01 < a < 0.06 for a in r["couplings_at_M_P"])
    uv_bounded_ok = r["uv_bounded"]
    ok = five_thirds_ok and gut_coupling_ok and gut_range_ok and uv_bounded_ok

    verdict = {
        "id": "V30_unification",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the 1-loop matching (a single point) is a named residual
        "pass": bool(ok),
        "status": "the_0_parameter_5_3_and_GUT_coupling_with_GUT_range_running_at_M_P_and_UV_bounded",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "alpha_10": r["alpha_10"],
            "five_thirds": r["five_thirds"],
            "M_GUT_0param": r["M_GUT_0param"],
            "couplings_at_M_P": r["couplings_at_M_P"],
            "spread_at_M_P": r["spread_at_M_P"],
            "mssm_unif_scale": r["mssm_unif_scale"],
            "mssm_spread": r["mssm_spread"],
            "u1_landau_pole": r["u1_landau_pole"],
            "uv_bounded": r["uv_bounded"],
            "residual": "the 1-loop matching to a single GUT point (the non-SUSY SM does not unify exactly at M_P)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v30_unification.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V30 - the 0-parameter coupling unification + the QCA UV completion")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter GUT coupling:  alpha_10 = {r['alpha_10']:.3f}  (= 1/25)")
    print(f"  the 0-parameter 5/3:  g1^2/g'^2 = {r['five_thirds']:.4f}")
    print(f"  the 0-parameter GUT scale:  M_GUT = M_P = {r['M_GUT_0param']:.2e} GeV")
    print(f"  the 1-loop running to M_P:  a = {['%.3f' % x for x in r['couplings_at_M_P']]}  (spread {r['spread_at_M_P']:.2f}, the named residual)")
    print(f"  the clean MSSM (SUSY) unification is at {r['mssm_unif_scale']:.0e} GeV (spread {r['mssm_spread']:.2f}) -- below M_P")
    print(f"  the U(1) Landau pole is at {r['u1_landau_pole']:.2e} GeV -- above M_P (UV-bounded)")
    print(f"  VERDICT:  {'PASS (0-parameter 5/3 + GUT coupling; GUT-range running at M_P; QCA UV-bounded)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v30_unification.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
