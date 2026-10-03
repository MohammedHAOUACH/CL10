#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V55 — the baryon asymmetry as a **product of 0-parameter numbers** (W5).

The W5 residual (Section 13): the **CP-violating factor** $\eta_{CP}$ in the
heavy-$N$ decay (V40 takes $\eta_{CP} = 10^{-3}$, the observed range
$10^{-3}\text{–}10^{-2}$).  This experiment presents the baryon asymmetry
as a **product of 0-parameter numbers** (the physical GUT-leptogenesis
structure):

    \Delta B \;=\; \frac{28}{79}\;\varepsilon_{CP}\;D

with $\frac{28}{79}$ = the **SM sphaleron**, $\varepsilon_{CP} = g_0^4/16\pi$
= the **1-loop CP-odd** asymmetry ($g_0 = 1/8 = 1/2^3$, **3** = the spatial
dimension), and $D = (1/16\pi)^2$ = the **2-loop dilution**.  **Every number
is 0-parameter** (root / QFT / SM); the loop counts are $1$ (CP) $+$
$2$ (dilution) $= 3$.  The result is $\Delta B = 6.8e-10$, a **ratio 1.12**
to the observed $6.1e-10$ (a **12%** match).  The **named residual**: the
CP **phase** (which makes $\varepsilon_{CP}$ nonzero) and the 2-digit ratio
1.12.
"""
import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import eta_cp

TOL = 1e-12


def main():
    t0 = time.time()
    r = eta_cp.verify()
    d = r["r"]
    st = d["structure"]
    ok = r["ok"]

    verdict = {
        "id": "V55_eta_cp",
        "tier": 6,
        "property": "P5,P6",
        "residual": True,   # the exact 2-digit CP magnitude and the exact Delta_B are the named residual
        "pass": bool(ok),
        "status": ("the_baryon_asymmetry_Delta_B=(28/79)x_eps_CP_x_D_is_a_product_of_0-parameter_numbers_"
                   "(28/79_the_SM_sphaleron; eps_CP=g0^4/16pi_the_1-loop_CP-odd_asymmetry_with_g0=1/2^3; "
                   "D=(1/16pi)^2_the_2-loop_dilution): Delta_B=6.8e-10_ratio_1.12_to_the_observed_6.1e-10; "
                   "the_CP_phase_and_the_2-digit_ratio_are_the_named_residual"),
        "tol": TOL,
        "data": {
            "backend": d["backend"],
            "loop_factor": d["loop_factor"],
            "eta_CP_0param": d["eta_CP_0param"],
            "eta_CP_observed": d["eta_CP_observed"],
            "eta_CP_over_observed": d["eta_CP_over_observed"],
            "g0": d["g0"],
            "cp_odd_loop_factor": d["cp_odd_loop_factor"],
            "sphaleron": d["sphaleron"],
            "delta_B_0param": d["delta_B_0param"],
            "delta_B_observed": d["delta_B_observed"],
            "delta_B_ratio": d["delta_B_ratio"],
            "delta_B_v40": d["delta_B_v40"],
            "delta_B_ratio_v40": d["delta_B_ratio_v40"],
            "structure": {
                "sphaleron": st["sphaleron"],
                "eps_CP": st["eps_CP"],
                "D_dilution": st["D_dilution"],
                "g0": st["g0"],
                "g0_is_1_over_2_cubed": st["g0_is_1_over_2_cubed"],
                "loop_count": st["loop_count"],
            },
            "residual": ("the CP phase (which makes eps_CP nonzero, an O(1) "
                         "input) and the 2-digit ratio 1.12 — the 28/79, "
                         "g_0 = 1/2^3, 16pi, and loop counts are 0-parameter"),
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v55_eta_cp.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V55 - the baryon asymmetry as a product of 0-parameter numbers (W5)")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print("  Delta_B = (28/79) x eps_CP x D, every factor 0-parameter:")
    print(f"    - sphaleron   = 28/79 = {st['sphaleron']:.4f}   (the SM)")
    print(f"    - eps_CP      = g_0^4/16pi = {st['eps_CP']:.3e}   (1-loop; g_0 = 1/2^3, 3 = spatial)")
    print(f"    - D (dilution)= (1/16pi)^2 = {st['D_dilution']:.3e}   (2-loop)")
    print(f"    - loop counts = 1 (CP) + 2 (dilution) = {st['loop_count']}")
    print(f"  0-parameter Delta_B = {d['delta_B_0param']:.3e}   (observed {d['delta_B_observed']:.1e};  ratio {d['delta_B_ratio']:.3f})")
    print(f"  V40 comparison:  Delta_B = {d['delta_B_v40']:.3e}  (ratio {d['delta_B_ratio_v40']:.3f})")
    print(f"  named residual (W5):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the baryon asymmetry is a product of 0-parameter numbers; 12% match; the CP phase is the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v55_eta_cp.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
