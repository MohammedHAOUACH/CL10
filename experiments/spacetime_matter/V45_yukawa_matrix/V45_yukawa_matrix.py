#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V45 — the 3-generation Yukawa mass matrix: the 10-to-1 projection (O6, C4).

The **10-to-1 reduction** (Section 12): the ten residuals O1-O10 are
projections of **one** 0-parameter object — the GUT QCA step ``U`` (the
32-level spectrum at ``g_10 = 1/5``).  This experiment is the **O6
projection**: the **3-generation Yukawa mass matrix** and the **mixing**,
read off from that one object, with **no fit**:

1. **The one object -> the 3 families.**  The 32-D Weyl module is **one**
   generation; the **3** families are the **3 two-planes** of the $(1,3)$
   3-space, so the 3-generation step is the family-symmetric
   ``U_96 = I_3 (x) U_32`` and its spectrum is **3 x (32 levels)** = 96
   levels — the one object, three times.

2. **The 0-parameter 3x3 mass matrix (the democratic base).**  At the
   symmetric (0-parameter) level the 3 families are **equal**: ``M =
   |1><1|`` (all-ones, rank 1), spectrum ``{3, 0, 0}`` — **one heavy**
   level + a **2-fold light** level (the 0-parameter **2-level structure**).

3. **The mixing (the TBM), with its honest split.**  The TBM **diagonalizes**
   the democratic base.  Of its angles, only $\\sin^2 12 = 1/3$ is **derived
   0-parameter**: it is the squared overlap with the **unique** heavy eigenvector
   $(1,1,1)/\\sqrt3$, and it does not move under the residual symmetry.  The
   light eigenspace is **2-fold degenerate**, so every diagonalizer is an $O(2)$
   rotation in it: $\\sin^2 23$ and $\\sin^2 13$ **sweep $[0, 2/3]$** (measured
   here).  The TBM values $1/2, 0$ are the member **adapted to the transposition
   (23)** of $S_3$ — and $(12),(13),(23)$ are equivalent symmetries of $|1\\rangle\\langle1|$,
   so the base does not single out $(23)$.  Close to the observed PMNS
   (0.31, 0.58, 0.02).

4. **The 3-level hierarchy (named).**  The **3-level** mass hierarchy
   (3 distinct masses $m_1 < m_2 < m_3$) is the **10-Higgs vev** breaking
   (it lifts the 2-fold light degeneracy); the **exact ratios** are the
   **named residual**.

5. **The 10-Higgs Dirac structure (0-parameter, V28).**  The grade-1 GUT
   vector (the 10) is **off-diagonal** in the 16/16-bar split, so the
   0-parameter mass term $16\\,10\\,16\\bar$ is the **Dirac** term.

**The named residual (O6).**  The **exact** 3-gen hierarchy (the 10-Higgs
vev breaking, the ratios $m_1:m_2:m_3$), the angles $\\sin^2 23$ and
$\\sin^2 13$ (the $O(2)$ choice of which $S_3$ transposition breaks), and the
**non-trivial CKM** (with $H_u = H_d$ both democratic the computed CKM is the
**identity**) — the one object does not fix them.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import yukawa_3x3

TOL = 1e-12


def main():
    t0 = time.time()
    r = yukawa_3x3.w6_verdict()
    one = r["one_object"]
    mm = r["mass_matrix"]
    mix = r["mixing"]
    hier = r["hierarchy"]
    fs = r["family_symmetry"]
    ok = r["ok"]

    verdict = {
        "id": "V45_yukawa_matrix",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the exact 3-gen hierarchy + the exact CKM/PMNS (O6) are the named residual
        "pass": bool(ok),
        "status": "the_0_parameter_3x3_Yukawa_mass_matrix_democratic_base_the_TBM_mixing_the_2_level_structure_the_exact_hierarchy_is_the_named_residual",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "g10": r["g10"],
            "n_levels_1gen": one["n_levels_1gen"],
            "n_levels_3gen": one["n_levels_3gen"],
            "family_triplication": one["family_triplication"],
            "zero_point_1gen": one["zero_point_1gen"],
            "det_U_minus_1": one["det_U_minus_1"],
            "mass_matrix": {
                "eigenvalues": mm["eigenvalues"],
                "rank": mm["rank"],
                "n_heavy": mm["n_heavy"],
                "n_light_2fold": mm["n_light_2fold"],
                "n_distinct_levels": mm["n_distinct_levels"],
                "democratic": mm["democratic"],
            },
            "family_symmetry": {
                "n_space": fs["n_space"],
                "n_families_C32": fs["n_families_C32"],
                "S_a2_err": fs["S_a2_err"],
                "anticommutator_Sa_Sb": fs["anticommutator_Sa_Sb"],
                "spatial_norms": fs["spatial_norms"],
                "norm_spread": fs["norm_spread"],
                "equivalent_3fold": fs["equivalent_3fold"],
                "democratic_spectrum": fs["democratic_spectrum"],
                "two_level": fs["two_level"],
            },
            "mixing": {
                "sin2_12": mix["sin2_12"],
                "sin2_23": mix["sin2_23"],
                "sin2_13": mix["sin2_13"],
                "obs_pmns": mix["obs_pmns"],
                "obs_ckm": mix["obs_ckm"],
                "dev_pmns": mix["dev_pmns"],
                "dev_ckm": mix["dev_ckm"],
                "diagonalizes_democratic": mix["diagonalizes_democratic"],
                "sin2_12_is_0param_derived": mix["sin2_12_is_0param_derived"],
                "sin2_23_is_0param_derived": mix["sin2_23_is_0param_derived"],
                "sin2_13_is_0param_derived": mix["sin2_13_is_0param_derived"],
                "sin2_23_O2_range": mix["sin2_23_O2_range"],
                "sin2_13_O2_range": mix["sin2_13_O2_range"],
                "mixing_status": mix["mixing_status"],
            },
            "hierarchy": {
                "n_levels_0param": hier["n_levels_0param"],
                "n_levels_observed": hier["n_levels_observed"],
                "structure": hier["structure"],
                "degenerate_light": hier["degenerate_light"],
            },
            "higgs_dirac": r["higgs_dirac"],
            "s2_unitary_err": r["s2_unitary_err"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v45_yukawa_matrix.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V45 - the 3-gen Yukawa mass matrix: the 10-to-1 projection (O6, C4)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the ONLY input:  g_10 = {r['g10']:.4f} (= 1/5)")
    print(f"  [1] the one object -> 3 families:  {one['n_levels_1gen']} levels (1 gen) x 3 = {one['n_levels_3gen']} levels (3 gen)")
    print(f"  [2] 3-fold family symmetry (the DERIVATION of the democratic 3x3):  the 3 spatial directions are "
          f"{'EQUIVALENT' if fs['equivalent_3fold'] else 'NOT'} (each S_a^2=-I, anticommuting, equal norm) -> the 10-Higgs couples the 3 families equally")
    print(f"  [3] 0-parameter 3x3 mass matrix (the democratic base):  eigenvalues "
          f"{['%.3g' % x for x in mm['eigenvalues']]}  ({mm['n_distinct_levels']}-level: {mm['n_heavy']} heavy + {mm['n_light_2fold']}-fold light, rank {mm['rank']})")
    obs = mix["obs_pmns"]
    hh = r["mixing_honest"]
    print(f"  [4] mixing (the TBM, which diagonalizes |1><1|: {mix['diagonalizes_democratic']}):  "
          f"sin^2(12)={mix['sin2_12']:.3f} (obs {obs[0]:.3f}), "
          f"sin^2(23)={mix['sin2_23']:.3f} (obs {obs[1]:.3f}), "
          f"sin^2(13)={mix['sin2_13']:.3f} (obs {obs[2]:.3f})")
    print("      HONEST: sin^2(12) = 1/3 is DERIVED 0-parameter (the unique heavy eigenvector);")
    print(f"              sin^2(23) sweeps [{hh['sin2_23_range'][0]:.3f}, {hh['sin2_23_range'][1]:.3f}] and")
    print(f"              sin^2(13) sweeps [{hh['sin2_13_range'][0]:.3f}, {hh['sin2_13_range'][1]:.3f}] under the O(2) freedom")
    print("              of the 2-fold light subspace — the TBM values are the (23)-adapted CHOICE.")
    print(f"  [5] 0-parameter hierarchy:  {hier['structure']}  (the {hier['n_levels_observed']}-level is the named residual);  "
          f"mean 0-parameter vs observed PMNS deviation:  {mix['dev_pmns']:.3f}")
    print(f"  [6] 10-Higgs Dirac term (the 0-parameter mass):  rank {r['higgs_dirac']['rank']} (16->16-bar, off-diagonal)")
    print(f"  [7] named residual (O6):  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the 0-parameter 3x3 matrix + TBM + 2-level structure; the 3-level/exact mixing is the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v45_yukawa_matrix.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
