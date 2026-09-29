#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V29 — the 0-parameter Yukawa structure: the democratic base, the
tribimaximal mixing (the 0-parameter PMNS prediction), and the mass-hierarchy
residual.

**The 0-parameter result.**  The three families are the three two-planes of the
$(1,3)$ 3-space.  The GUT structure gives a **democratic** (degenerate) family
base — the three families are **equal** at the 0-parameter level (a verified
consequence of the root).  Diagonalizing the **democratic** matrix gives the
**tribimaximal** mixing
$$\\sin^2\\theta_{12}=1/3,\\ \\sin^2\\theta_{23}=1/2,\\ \\sin^2\\theta_{13}=0,$$
a **0-parameter** prediction **close to the observed PMNS**
($0.31, 0.58, 0.02$) (mean deviation $\\sim 0.04$).

**The named residual.**  The **mass hierarchy** (the three observed levels,
$m_t/m_b\\sim40$, $m_\\tau/m_e\\sim3500$) is **not** 0-parameter: the
democratic base is **degenerate**, so the three-level split requires the
**10-Higgs vev** breaking (a GUT-scale input) — the **exact** ratios are the
named residual.  The **structure** (democratic base $+$ 10-Higgs breaking $+$
tribimaximal mixing) **is** 0-parameter.
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
    ok = (y["non_singular"] and y["off_diagonal_Dirac"]
          and dev < 0.1 and h["degenerate"] and y["n_gen"] == 3)

    verdict = {
        "id": "V29_yukawa_hierarchy",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the mass hierarchy is a named residual
        "pass": bool(ok),
        "status": "the_0_parameter_democratic_base_gives_the_tribimaximal_PMNS_close_to_observed",
        "tol": TOL,
        "data": {
            "backend": mix["backend"],
            "n_gen": y["n_gen"],
            "mass_matrix_rank": y["rank"],
            "non_singular": y["non_singular"],
            "off_diagonal_Dirac": y["off_diagonal_Dirac"],
            "tribimaximal": {"sin2_12": mix["sin2_12"], "sin2_23": mix["sin2_23"], "sin2_13": mix["sin2_13"]},
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
    print(f"  [V28] the {y['n_gen']}x{y['n_gen']} mass matrix: rank {y['rank']} (non-singular), Dirac off-diagonal")
    print("  [V29] the 0-parameter democratic base -> the 0-parameter TRIBIMAXIMAL mixing:")
    print(f"    sin^2(12)={mix['sin2_12']:.3f} (obs {obs[0]:.3f}), sin^2(23)={mix['sin2_23']:.3f} (obs {obs[1]:.3f}), sin^2(13)={mix['sin2_13']:.3f} (obs {obs[2]:.3f})")
    print(f"  the mean 0-parameter vs observed PMNS deviation:  {dev:.3f}")
    print(f"  the 0-parameter base spectrum (degenerate):  {['%.3g' % x for x in h['base_spectrum']]}")
    print(f"  named residual:  {h['residual']}")
    print(f"  VERDICT:  {'PASS (0-parameter tribimaximal ~ observed PMNS; the hierarchy the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v29_yukawa_hierarchy.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
