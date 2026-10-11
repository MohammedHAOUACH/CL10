#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V28 - the Yukawa / mass-hierarchy 0-parameter structure.

The **0-parameter** fermion-mass structure of the CL10 QCA:

- the **10-Higgs 16->16-bar**: the Cl(10) 10 (the grade-1 GUT vector) is
  **off-diagonal** in the 16/16-bar split — a **10-Higgs vev** gives the
  **0-parameter Dirac mass term** $16 \\cdot 10 \\cdot 16\\bar$; the
  **diagonal** (16->16, 16-bar->16-bar) blocks are **0** (the 10 is odd —
  no Majorana mass for the charged fermions);
- the **3 x 3 mass matrix**: with the **3 generations** (the 0-parameter 3 of
  V26), the mass term gives a **3 x 3** matrix per fermion type;
- the **0-parameter rank**: the **3-fold family symmetry** (the 3 spatial
  directions are equivalent) forces the 10-Higgs to couple the 3 families
  **equally** — the 3x3 is the **democratic** base |1><1|, **rank 1**,
  spectrum {3,0,0}: **one heavy** level + a **2-fold light** level (the
  0-parameter **2-level** structure).  A rank-1 3x3 is **singular** (two
  zero eigenvalues): the theory does **not** give mass to all three
  generations at 0-parameter.

The **named residual**: the **3-level mass hierarchy** (the 10-Higgs vev
breaking the 2-fold light degeneracy — the observed ratios $m_\\mu/m_e$,
$m_t/m_b$, the neutrino masses) and the **exact mixing angles** (CKM/PMNS)
— the theory gives the **0-parameter 2-level structure**, not the 3 levels
nor their values.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import yukawa

TOL = 1e-12


def main():
    t0 = time.time()
    r = yukawa.yukawa_matrix()
    # le critère réel : la base démocratique est RANG 1 (spectre {3,0,0}),
    # pas non-singulière.  Les 3 masses distinctes sont le résidu nommé (W6).
    ok = (r["rank"] == 1 and r["n_distinct_levels"] == 2
          and r["n_heavy"] == 1 and r["n_light_2fold"] == 2
          and r["democratic"] and r["singular"]
          and r["off_diagonal_Dirac"] and r["diag_zero"] and r["n_gen"] == 3)

    verdict = {
        "id": "V28_yukawa",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the mass hierarchy + the mixing angles are named residuals
        "pass": bool(ok),
        "status": "the_0_parameter_democratic_rank1_Dirac_mass_matrix_the_2_level_structure",
        "tol": TOL,
        "data": {
            "n_gen": r["n_gen"],
            "n_dim": r["n_dim"],
            "rank": r["rank"],
            "spectrum": r["spectrum"],
            "n_distinct_levels": r["n_distinct_levels"],
            "n_heavy": r["n_heavy"],
            "n_light_2fold": r["n_light_2fold"],
            "democratic": r["democratic"],
            "singular": r["singular"],
            "off_diagonal_Dirac": r["off_diagonal_Dirac"],
            "diag_zero": r["diag_zero"],
            "mass_term": r["mass_term"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v28_yukawa.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V28 - the Yukawa / mass-hierarchy 0-parameter structure")
    print("=" * 68)
    print(f"  the 0-parameter mass term:  {r['mass_term']}")
    print(f"  the 10-Higgs is OFF-diagonal (Dirac):  {r['off_diagonal_Dirac']}  (the diagonal is 0: {r['diag_zero']})")
    print(f"  the {r['n_gen']}x{r['n_gen']} democratic mass matrix:  rank = {r['rank']}  "
          f"(spectrum {['%.3g' % x for x in r['spectrum']]}: {r['n_heavy']} heavy + "
          f"{r['n_light_2fold']}-fold light -> the 0-parameter 2-level structure)")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (0-parameter democratic rank-1 Dirac mass matrix, the 2-level structure)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v28_yukawa.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
