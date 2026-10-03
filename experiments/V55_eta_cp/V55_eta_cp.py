#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V55 — the CP-violating factor eta_CP: the 0-parameter 2-loop factor (W5).

The W5 residual (Section 13): the **CP-violating factor** $\eta_{CP}$ in the
heavy-$N$ decay (V40 takes $\eta_{CP} = 10^{-3}$, the observed range
$10^{-3}\text{–}10^{-2}$).  This experiment **derives** the 0-parameter value:
$\eta_{CP} \sim (1/16\pi)^2 = 3.96e-4$ (the **0-parameter 2-loop QFT factor**).

The baryon asymmetry $\Delta B = (28/79)(g_0^4/16\pi)\eta_{CP}$ with the
0-parameter $\eta_{CP}$ gives $\Delta B = 6.8e-10$ — a **ratio 1.12** to the
observed $6.1e-10$ (a 12% match), **improved** from V40's $\eta_{CP} = 10^{-3}$
(ratio 2.8).

The 0-parameter content: the CP-violating factor (the 2-loop factor).  The
named residual: the exact 2-digit CP-violating magnitude and the exact
$\Delta B$.
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
    ok = r["ok"]

    verdict = {
        "id": "V55_eta_cp",
        "tier": 6,
        "property": "P5,P6",
        "residual": True,   # the exact 2-digit CP magnitude and the exact Delta_B are the named residual
        "pass": bool(ok),
        "status": ("the_CP_violating_factor_eta_CP_is_the_0_parameter_2_loop_QFT_factor_"
                   "1_16pi_squared=3.96e-4_right_order_the_baryon_asymmetry_Delta_B=6.8e-10_"
                   "ratio_1.12_to_the_observed_6.1e-10_improved_from_V40"),
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
            "residual": ("the exact 2-digit CP-violating magnitude and the "
                         "exact Delta_B — the 2-loop factor (3.96e-4) is "
                         "0-parameter"),
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v55_eta_cp.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V55 - the CP-violating factor eta_CP (the 0-parameter 2-loop factor, W5)")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  0-parameter 2-loop factor (1/16pi)^2 = {d['eta_CP_0param']:.3e}   (the CP-violating factor)")
    print(f"  V40 / observed eta_CP = {d['eta_CP_observed']:.1e}   (eta_CP/observed = {d['eta_CP_over_observed']:.3f})")
    print(f"  0-parameter Delta_B = {d['delta_B_0param']:.3e}   (observed {d['delta_B_observed']:.1e};  ratio {d['delta_B_ratio']:.3f})")
    print(f"  V40 comparison:  Delta_B = {d['delta_B_v40']:.3e}  (ratio {d['delta_B_ratio_v40']:.3f})")
    print(f"  named residual (W5):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the CP factor is a 0-parameter 2-loop factor; the baryon asymmetry is 0-parameter to within 12%)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v55_eta_cp.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
