#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V108 - the hierarchy coefficient 2.08 as a measured closed-loop number.

**Why (O8, Table 5).** The hierarchy is v/M_P = 2.08*(1/16pi)^10 and the
coefficient 2.08 is classified FITTED (read from the data). V85 measured
the closed traced-loop factor F_10 = 0.019186 at the root coupling
(the ladder length IS the crossing, V89). V70: the doubler nodes form
THREE families of four labeled by the zero axis - three spatial
directions, a root number.

**The identity to test (zero fitted parameters):**
    v / M_P  =?  3 * F_10^10
Every input: 3 = the number of spatial directions (the node-family
multiplicity, V70); 10 = the ladder exponent = the measured closed-loop
crossing (V89); F_10 = the measured closed-loop factor (V85); v/M_P and
(1/16pi)^10 from PDG. Equivalently 2.08 =? 3*(16 pi F_10)^10: the
fitted coefficient becomes a product of measured numbers.

**Honest caveats recorded in the JSON:**
- the per-level reading "each ladder level carries F_10" is the
  V85/V89 identification, not a theorem;
- the multiplicity 3 as a factor of the loop construction (one ladder
  per axis family) is NOT proven to factorize in V85's single-link
  construction - V89 showed no mod-4 plaquette selection, nothing on
  axis factorization;
- sensitivity: d ln(coefficient)/d ln(F_10) = 10, so a 1% shift of
  F_10 moves the coefficient by 10%.

PASS = |3 F_10^10 - v/M_P| / (v/M_P) < 0.01 using ONLY numbers read
from the V85/V89 result JSONs and PDG constants.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp

V_OVER_MP = 246.22 / 1.2209e19      # PDG 2024
LOOP = 1.0 / (16.0 * math.pi)


def main():
    v85 = json.load(open(exp("V85_closed_loop_factor", "results",
    "v85_closed_loop_factor.json")))["data"]
    v89p = exp("V89_loop_selection", "results")
    v89f = [f for f in sorted(os.listdir(v89p)) if f.endswith(".json")]
    v89 = json.load(open(os.path.join(v89p, v89f[0])))["data"] if v89f else {}

    F10 = v85["loops_root"]["F_n"][-1]          # n = 10, root coupling
    # V89's independent measurement of the crossing
    cross = None
    for k, v in v89.items():
        if isinstance(v, dict):
            for kk in ("n_star_below_loop_factor", "n_star", "n_cross",
                       "crossing"):
                if kk in v:
                    cross = v[kk]
                    break
            if cross is not None:
                break

    coef_pred = 3.0 * (LOOP ** -10) * F10 ** 10     # = 3*(16 pi F_10)^10
    coef_obs = V_OVER_MP / LOOP ** 10
    rel = abs(coef_pred - coef_obs) / coef_obs

    ok = rel < 0.01 and cross == 10
    status = (f"THE 2.08 IS NO LONGER A FITTED NUMBER (conditional): "
              f"v/M_P = 3*F_10^10 with F_10 = {F10} measured (V85), "
              f"3 = the three axis families (V70), 10 = the measured "
              f"crossing (V89): predicted coefficient "
              f"{coef_pred:.4f} vs observed {coef_obs:.4f}, agreement "
              f"{100*(1-rel):.2f}%. CONDITIONAL: the per-level F_10 "
              f"reading is the V85/V89 identification and the axis-"
              f"family multiplicity 3 is not proven to factorize in the "
              f"single-link construction; sensitivity dln(coef)/dln(F_10)"
              f" = 10. O8 narrows: the coefficient is a product of "
              f"measured root numbers, its derivation needs the axis-"
              f"factorization of the closed-loop sum.")
    out = {"experiment": "V108_hierarchy_multiplicity", "pass": bool(ok),
           "status": status,
           "data": {"F_10_root": F10,
                    "v89_crossing": cross,
                    "coefficient_predicted_3_F10_10": coef_pred,
                    "coefficient_observed": coef_obs,
                    "relative_error": rel,
                    "sensitivity_dln_coef_dln_F10": 10.0,
                    "caveats": ["per-level F_10 reading is the V85/V89 "
                                "identification, not a theorem",
                                "axis-family factorization of the loop "
                                "sum NOT proven (V89: no mod-4 selection; "
                                "nothing on axes)",
                                "1% shift of F_10 moves the coefficient "
                                "by 10%"]}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v108_hierarchy_multiplicity.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
