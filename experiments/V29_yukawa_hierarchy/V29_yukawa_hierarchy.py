#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V29 — the 0-parameter Yukawa structure: the democratic base, the
tribimaximal mixing, which diagonalizes it: only sin^2(12) = 1/3 is derived
0-parameter (the unique heavy eigenvector); sin^2(23) and sin^2(13) sweep [0, 2/3]
under the O(2) freedom of the 2-fold light subspace), and the mass-hierarchy
residual.

**The 0-parameter result.**  The three families are the three two-planes of the
$(1,3)$ 3-space.  The GUT structure gives a **democratic** (degenerate) family
base — the three families are **equal** at the 0-parameter level (a verified
consequence of the root).  Diagonalizing the **democratic** matrix gives the
**tribimaximal** mixing
$$\\sin^2\\theta_{12}=1/3,\\ \\sin^2\\theta_{23}=1/2,\\ \\sin^2\\theta_{13}=0,$$
close to the observed PMNS ($0.31, 0.58, 0.02$) (mean deviation $\\sim 0.04$).

**Honest status of each angle** (measured here, not asserted).  The democratic
base has spectrum {3,0,0}: the heavy eigenvector (1,1,1)/sqrt(3) is **unique**,
but the light eigenspace is **2-fold degenerate**, so every diagonalizer is an
O(2) rotation in that subspace.  Sweeping the O(2) freedom:

- sin^2(theta_12) = 1/3 is **truly derived 0-parameter**: it is the squared
  overlap of any family basis vector with the **unique** heavy eigenvector, and
  it does not move under the O(2) freedom.  Observed 0.307 — an **8%** match,
  and the one genuine 0-parameter mixing prediction of the model.
- sin^2(theta_23) and sin^2(theta_13) **sweep [0, 2/3]**: they are **not**
  derived.  The TBM values (1/2, 0) are the member of the family **adapted to
  the transposition (23)** of S3.  Since (12), (13) and (23) are **equivalent**
  symmetries of |1><1|, nothing in the democratic base singles out (23): which
  transposition breaks is an **additional named input**.

**The named residual.**  The **mass hierarchy** (the three observed levels,
$m_t/m_b\\sim40$, $m_\\tau/m_e\\sim3500$) is **not** 0-parameter: the
democratic base is **degenerate**, so the three-level split requires the
**10-Higgs vev** breaking (a GUT-scale input) — the **exact** ratios are the
named residual.  The **structure** (democratic base $+$ 10-Higgs breaking $+$
tribimaximal mixing) **is** 0-parameter **as a structure** — but only
sin^2(theta_12) = 1/3 is 0-parameter **as a number**; theta_23 and theta_13
are an O(2) convention.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import yukawa

TOL = 1e-12


def main():
    t0 = time.time()
    y = yukawa.yukawa_matrix()
    mix = yukawa.tribimaximal_mixing()
    h = yukawa.hierarchy_residual()
    obs = mix["obs"]
    dev = (abs(mix["sin2_12"] - obs[0]) + abs(mix["sin2_23"] - obs[1])
           + abs(mix["sin2_13"] - obs[2])) / 3
    ok = (y["rank"] == 1 and y["n_distinct_levels"] == 2
          and y["off_diagonal_Dirac"]
          and dev < 0.1 and h["degenerate"] and y["n_gen"] == 3)

    verdict = {
        "id": "V29_yukawa_hierarchy",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the mass hierarchy is a named residual
        "pass": bool(ok),
        "status": "the_democratic_base_derives_sin2_12=1/3_0_parameter_(8%_from_observed);_sin2_23_and_sin2_13_sweep_[0,2/3]_under_the_O2_freedom_(the_TBM_is_the_23-adapted_choice)",
        "tol": TOL,
        "data": {
            "backend": mix["backend"],
            "n_gen": y["n_gen"],
            "mass_matrix_rank": y["rank"],
            "mass_matrix_spectrum": y["spectrum"],
            "n_distinct_levels": y["n_distinct_levels"],
            "democratic": y["democratic"],
            "singular": y["singular"],
            "off_diagonal_Dirac": y["off_diagonal_Dirac"],
            "tribimaximal": {"sin2_12": mix["sin2_12"], "sin2_23": mix["sin2_23"], "sin2_13": mix["sin2_13"]},
            "sin2_12_is_0param_derived": mix["sin2_12_is_0param_derived"],
            "sin2_23_is_0param_derived": mix["sin2_23_is_0param_derived"],
            "sin2_13_is_0param_derived": mix["sin2_13_is_0param_derived"],
            "sin2_23_O2_range": mix["sin2_23_O2_range"],
            "sin2_13_O2_range": mix["sin2_13_O2_range"],
            "diagonalizes_democratic": mix["diagonalizes_democratic"],
            "observed_PMNS": {"sin2_12": obs[0], "sin2_23": obs[1], "sin2_13": obs[2]},
            "mean_mixing_deviation": dev,
            "base_spectrum": h["base_spectrum"],
            "degenerate_base": h["degenerate"],
            "residual": h["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v29_yukawa_hierarchy.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V29 - the 0-parameter Yukawa structure (democratic base + tribimaximal)")
    print("=" * 68)
    print(f"  backend (GPU):  {mix['backend']}")
    print(f"  [V28] the {y['n_gen']}x{y['n_gen']} democratic mass matrix: rank {y['rank']} "
          f"(spectrum {['%.3g' % x for x in y['spectrum']]}: 1 heavy + 2-fold light), Dirac off-diagonal")
    print("  [V29] the democratic base -> the TRIBIMAXIMAL mixing (it diagonalizes |1><1|: "
          f"{mix['diagonalizes_democratic']}):")
    print(f"    sin^2(12)={mix['sin2_12']:.3f} (obs {obs[0]:.3f}), sin^2(23)={mix['sin2_23']:.3f} (obs {obs[1]:.3f}), sin^2(13)={mix['sin2_13']:.3f} (obs {obs[2]:.3f})")
    print(f"    HONEST: sin^2(12)=1/3 is DERIVED 0-parameter (the unique heavy eigenvector, "
          f"8% from observed);")
    print(f"            sin^2(23) sweeps {tuple(round(x,3) for x in mix['sin2_23_O2_range'])} and "
          f"sin^2(13) sweeps {tuple(round(x,3) for x in mix['sin2_13_O2_range'])} (the O(2) freedom) —")
    print(f"            the TBM values are the (23)-adapted CHOICE, not a derivation.")
    print(f"  the mean 0-parameter vs observed PMNS deviation:  {dev:.3f}")
    print(f"  the 0-parameter base spectrum (degenerate):  {['%.3g' % x for x in h['base_spectrum']]}")
    print(f"  named residual:  {h['residual']}")
    print(f"  VERDICT:  {'PASS (sin^2(12)=1/3 derived 0-parameter ~ observed; sin^2(23), sin^2(13) and the 3-level hierarchy are named residuals)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v29_yukawa_hierarchy.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
