#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V59 — the **unified** solution to the four limits of the theory.

The four "limits" (``doc/remarques.txt``) are not four independent gaps:
they are the **four projections of one 0-parameter invariant**,
``10 = 2 x 5`` (the ``2`` = the chiral 16/16-bar split; the ``5`` = the
qubit count of the 32-D module), **plus one irreducible CP phase**.  This
experiment computes the unification (``cl10.four_limits``) and reduces the
residual set to its minimum:

- **L1 (conventionality):** the 7 algebraic invariants are ``g_10``-invariant
  (the 16/16-bar CPT mirror holds for **every** angle, residual ``~1e-15``)
  -- **0 free parameters in the algebra**; the dynamics is the one-parameter
  family ``U(g_10)`` whose selected value ``1/5`` is a 0-parameter
  **selection rule** (the *reading* is the audit's mapping convention).
- **L2 (hierarchy) + L4 (DM):** the **same** loop count ``n = 10`` is the
  **unique** ``n`` in ``1..20`` within **one decade** of **both** closures:
  the hierarchy (``2.08 (1/16\pi)^10``) and the DM WIMP mass
  (``1.18 x 100 GeV = 118`` GeV) -- **two statements of one scale** (the
  "9-to-1"), **not** two independent constraints; the strict factor-2 test
  of both is **empty**.  The 2-digit coefficients are **named residuals**
  (``observed / 10-loop factor``).
- **L3 (baryogenesis):** the rank-1 (1 heavy ``N``) Yukawa has **0** CP
  phases (no CP violation, no baryon); the 2-step (2 heavy ``N``) seesaw
  has **1** irreducible CP phase (``Im det(Y) = |sin delta|``).  The baryon
  asymmetry needs **exactly one** CP phase; the 0-parameter **factors**
  (the 1-loop ``g_0^4/16\pi`` rate at the 0-parameter coupling
  ``g_0 = 1/8`` (a mapping convention), the ``(1/16\pi)^2`` dilution, the
  ``28/79`` sphaleron) give the **right order** (``6.8e-10``, ratio
  ``1.11``); the 0-parameter QCA alone gives ``epsilon_CP = 0`` (the
  CP-odd theorem).

The **reduced residual set**: **0** free parameters in the algebra +
**1** 't Hooft angle of the dynamics (selected value = a 0-parameter
selection rule) + **2** mapping conventions (``g_0``, ``M_GUT = M_P``) +
**2** two-digit coefficients (``2.08``, ``1.18`` -- **named**, read from
the data) + **1** irreducible CP phase (``delta``).  Nothing else.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import four_limits

TOL = 1e-9


def main():
    t0 = time.time()
    r = four_limits.verify(TOL)
    d = r["r"]
    ok = bool(r["ok"])

    l1 = d["L1_conventionality"]
    l24 = d["L2_hierarchy_plus_L4_DM"]
    l3 = d["L3_baryogenesis"]
    cp = d["coupling_provenance"]

    verdict = {
        "id": "V59_four_limits",
        "tier": 6,
        "property": "P1",
        "residual": True,   # 0 free in the algebra; 3 named (2 two-digit coeffs + 1 CP phase) + 3 mapping conventions
        "pass": ok,
        "status": ("the_4_limits_reduce_to_0_free_parameters_IN_THE_ALGEBRA_(the_7_invariants_hold_for_every_g10)"
                   "plus_1_t_Hooft_angle_of_the_dynamics_whose_value_g10=1/5_is_a_mapping_convention_not_a_selection_rule;"
                   "the_4_limits_are_the_4_projections_of_one_0-parameter_invariant_10=2x5_plus_one_CP_phase;"
                   "reduced_to_3_named_(2_two-digit_coeffs_2.08_1.18_read_from_data_and_1_CP_phase)_plus_3_mapping_conventions;"
                   "unique_n=10_within_one_decade_of_both_statements_of_one_scale_(the_strict_factor2_test_is_empty);"
                   "baryogenesis_needs_one_CP_phase_(the_0-param_QCA_alone_gives_zero_by_the_CP-odd_theorem);"
                   "delta_is_INDEPENDENT_of_g10_(the_QCA_is_real/CP-even)"),
        "tol": TOL,
        "data": {
            "backend": d["backend"],
            "unifying_invariant": d["unifying_invariant"],
            "L1_g10_invariant": bool(l1["g10_invariant"]),
            "L1_max_mirror_residual": float(l1["max_mirror_residual_over_grid"]),
            "L1_grid": l1["grid"],
            "L2L4_unique_n": int(l24["n_best"]),
            "L2L4_n_is_2x5": bool(l24["n_gen_is_2x5"]),
            "L2L4_unique_of_both": bool(l24["unique_n_within_1decade_of_both"]),
            "L2L4_coeff_hierarchy": float(l24["coeff_hierarchy"]),
            "L2L4_coeff_dm": float(l24["coeff_dm"]),
            "L2L4_candidates": l24["n_candidates_within_1decade"],
            "L2L4_factor2_test_empty": bool(l24["factor2_test_of_both_empty"]),
            "L3_rank1_cp_phases": int(l3["rank1_cp_phases"]),
            "L3_twoN_cp_phases": int(l3["minimal_cp_phases_required"]),
            "L3_cpodd_at_0": float(l3["cpodd_at_delta_0"]),
            "L3_cpodd_at_pi2": float(l3["cpodd_at_delta_pi2"]),
            "L3_delta_B_0param": float(l3["delta_B_0param"]),
            "L3_delta_B_observed": float(l3["delta_B_observed"]),
            "L3_delta_B_ratio": float(l3["delta_B_ratio"]),
            "free_parameters_in_algebra": int(d["residual_set"]["free_parameters_in_algebra"]),
            "free_parameters_in_dynamics": int(d["residual_set"]["free_parameters_in_dynamics"]),
            "g10_is_root_reading": bool(cp["g10_is_root_reading"]),
            "g0_is_root_reading": bool(cp["g0_is_root_reading"]),
            "g10_root_reading": float(cp["g10_root_reading"]),
            "g0_root_reading": float(cp["g0_root_reading"]),
            "qca_cp_even": bool(cp["qca_cp_even"]),
            "delta_independent_of_g10": bool(cp["delta_independent_of_g10"]),
            "mapping_conventions": cp["mapping_conventions"],
            "residual_set": d["residual_set"],
            "total_named_residuals": int(d["total_named_residuals"]),
            "verdict": d["verdict"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v59_four_limits.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V59 - the unified solution to the four limits  [tier 6, P1]")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  unifying invariant:  {d['unifying_invariant']}")
    print(f"  [L1]  g_10-invariant CPT mirror:  {l1['g10_invariant']}  "
          f"(max residual {l1['max_mirror_residual_over_grid']:.2e} over g_10=0.1..0.5)")
    print(f"  [L2+L4]  unique loop count n = {l24['n_best']} (= 2 x 5):  unique-within-1-decade-of-both = {l24['unique_n_within_1decade_of_both']}  (strict factor-2 test: {l24['factor2_test_of_both_empty'] and 'empty'})")
    print(f"            named residuals:  hierarchy = {l24['coeff_hierarchy']:.3f},  DM = {l24['coeff_dm']:.3f}")
    print(f"  [L3]  rank-1 CP phases = {l3['rank1_cp_phases']} (none);  2-step (2N) = {l3['minimal_cp_phases_required']} (one)  "
          f"|Delta B| = {l3['delta_B_0param']:.2e} (ratio {l3['delta_B_ratio']:.2f})")
    print("  [couplings]  g_10 = 1/5 (selection rule);  g_0 = 1/8 (selection rule);  -> 0 free in the algebra, 1 angle in the dynamics")
    print(f"  reduced residual set:  {d['residual_set']}  (total {d['total_named_residuals']} named)")
    print(f"  VERDICT:  {'PASS (0 free in algebra; 3 named = 2 two-digit coeffs + 1 CP phase)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v59_four_limits.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
