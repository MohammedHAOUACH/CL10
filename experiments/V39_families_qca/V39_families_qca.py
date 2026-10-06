#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V39 — the three families **inside** the QCA dynamics (C4): the
3-generation 0-parameter unitary QCA step on the 96-D module.

**The 0-parameter result.**  The 3-generation internal space is
$H_{3f} = \\bigoplus_{k=1}^{3} H_{32}$ (a **96-D** module: **48 Weyl** +
**48 CPT-conjugate**): the 3 families are the **3 two-planes** $\\binom{3}{2}$
of the $(1,3)$ 3-space (the V26 reading), each family one **32-D** copy.
The **3-generation QCA step** is the **family-symmetric** 0-parameter step

$$U_{96} = \\mathrm{blockdiag}(U_{32},\\, U_{32},\\, U_{32}),\\qquad
U_{32} = \\prod_{a<b} [\\cos(g_{10})\\,I - i\\,\\sin(g_{10})\\, T_{ab}],$$

i.e. the **45** Spin(10) plane rotations at $g_{10} = 0.2$ act **identically
on each family** — the 0-parameter statement that the 45 do **not**
distinguish the 3 two-planes (the **family symmetry**).  The 0-parameter
**family mixing** of the GUT is the **tribimaximal** 3x3 unitary (the V29
angles $\\sin^2\\theta_{12}=1/3$, $\\sin^2\\theta_{23}=1/2$, $\\theta_{13}=0$), a
property of the neutrino/Yukawa sector, not of the time-evolution step.  Of
these only $\\sin^2\\theta_{12}=1/3$ is 0-parameter **as a number**; $\\theta_{23}$
and $\\theta_{13}$ are an $O(2)$ convention (Theorem 13, V60).

**The checks.**
1. the 96x96 step is **exactly unitary** ($\\|U U^{\\dagger} - I\\| \\sim
   10^{-15}$);
2. the **family symmetry**: $U_{96}$ is blockdiag of **3 equal** 32x32
   blocks (the 45 do not distinguish the 3 planes);
3. the step **preserves the per-family 16/16-bar split**: the 45 commute
   with the per-family pseudoscalar (the V22/V33 logic), and the 96-D
   step commutes with the 96-D pseudoscalar;
4. the **tribimaximal mixing** (V29) is unitary **and diagonalizes the
   democratic base** (the check the earlier code lacked); of its angles only
   $\\sin^2\\theta_{12}=1/3$ is 0-parameter **as a number** (the unique heavy
   eigenvector) — $\\theta_{23}, \\theta_{13}$ are an $O(2)$ convention;
5. the **dimension**: 48 Weyl $= 16 \\times 3$ and 48 CPT-conjugate
   (the antimatter).

**The named residual.**  the why-3 (W3, as before); the family **mass**
hierarchy (the 10-Higgs breaking, as before).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import family

TOL = 1e-12


def main():
    t0 = time.time()
    r = family.family_qca()
    ok = (r["unitary_err"] < TOL
          and r["family_symmetry_off_block"] < TOL
          and r["family_symmetry_block_diff"] < TOL
          and r["max_45_pseudoscalar_commutator"] < 1e-9
          and r["step_32_pseudoscalar_commutator"] < 1e-9
          and r["step_96_pseudoscalar_commutator"] < 1e-9
          and r["tribimaximal_unitary_err"] < TOL
          and r["tribimaximal_angles_ok"]
          # the regression guard: the TBM must DIAGONALIZE the democratic base.
          # The old code checked only the angles (the moduli), which a bimaximal
          # matrix with the wrong sign also satisfies.
          and r["tribimaximal_diagonalizes_democratic"]
          and r["n_weyl_3gen"] == r["n_weyl_1gen"] * 3 == 48
          and r["n_cpt_3gen"] == 48
          and r["n_module_3gen"] == 96)

    verdict = {
        "id": "V39_families_qca",
        "tier": 1,
        "property": "P1",
        "residual": True,   # the why-3 (W3) and the family mass hierarchy are named
        "pass": bool(ok),
        "status": "the_three_families_inside_the_qca_dynamics_the_96d_step_is_0_parameter_unitary_and_family_symmetric",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "g10": r["g10"],
            "n_planes": r["n_planes"],
            "n_module_1gen": r["n_module_1gen"],
            "n_module_3gen": r["n_module_3gen"],
            "n_weyl_1gen": r["n_weyl_1gen"],
            "n_weyl_3gen": r["n_weyl_3gen"],
            "n_cpt_3gen": r["n_cpt_3gen"],
            "unitary_err": r["unitary_err"],
            "family_symmetry_off_block": r["family_symmetry_off_block"],
            "family_symmetry_block_diff": r["family_symmetry_block_diff"],
            "max_45_pseudoscalar_commutator": r["max_45_pseudoscalar_commutator"],
            "step_32_pseudoscalar_commutator": r["step_32_pseudoscalar_commutator"],
            "step_96_pseudoscalar_commutator": r["step_96_pseudoscalar_commutator"],
            "tribimaximal_unitary_err": r["tribimaximal_unitary_err"],
            "tribimaximal_sin2": r["tribimaximal_sin2"],
            "tribimaximal_diagonalizes_democratic": r["tribimaximal_diagonalizes_democratic"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v39_families_qca.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V39 - the 3 families inside the QCA dynamics (the 3-generation 0-parameter step)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  [check 1] the {r['n_module_3gen']}x{r['n_module_3gen']} step is EXACTLY unitary:")
    print(f"    ||U U^dag - I||_max = {r['unitary_err']:.3e}")
    print("  [check 2] the family symmetry (the 45 do not distinguish the 3 planes):")
    print(f"    the off-block entries:  max = {r['family_symmetry_off_block']:.3e}  (0: blockdiag)")
    print(f"    the 3 blocks are EQUAL: max|block_i - block_1| = {r['family_symmetry_block_diff']:.3e}")
    print("  [check 3] the per-family 16/16-bar split is preserved:")
    print(f"    the 45 vs the per-family pseudoscalar:  max comm = {r['max_45_pseudoscalar_commutator']:.3e}")
    print(f"    the 32-D GUT step vs its pseudoscalar:  max comm = {r['step_32_pseudoscalar_commutator']:.3e}")
    print(f"    the 96-D step vs the 96-D pseudoscalar: max comm = {r['step_96_pseudoscalar_commutator']:.3e}")
    s = r["tribimaximal_sin2"]
    print(f"  [check 4] the tribimaximal mixing (V29):  sin^2(12) = {s['s12']:.4f} (= 1/3), "
          f"sin^2(23) = {s['s23']:.4f} (= 1/2), sin^2(13) = {s['s13']:.1e} (= 0)")
    print(f"    ||U_TBM U_TBM^T - I||_max = {r['tribimaximal_unitary_err']:.3e}  (unitary)")
    print(f"    U_TBM diagonalizes the democratic base |1><1|: {r['tribimaximal_diagonalizes_democratic']}")
    print(f"    (only sin^2(12)=1/3 is 0-parameter as a number; sin^2(23), sin^2(13) are an O(2) choice)")
    print(f"  [check 5] dimension:  {r['n_weyl_3gen']} Weyl = {r['n_weyl_1gen']} x 3;  {r['n_cpt_3gen']} CPT-conjugate")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the 3 families are inside the QCA dynamics; the 96x96 step is 0-parameter unitary and family-symmetric)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v39_families_qca.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
