#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V52 - the honest parameter audit: the solution to the four limits of the
theory.

**The four limits** — the honest note "what a grand physicist should also
flag" — are all facets of **one** meta-problem: the label
**"0-parameter"** is applied, undifferentiated, to numbers of **four**
different epistemic statuses.  This experiment is the **audit**: it
**computes** (on the GPU, the CuPy engine) the four resolutions and
**cross-checks** the root by a **second route** (the symbolic ``Cl(10,0)``
of the ``kingdon`` Clifford-algebra library):

**Point 1 — "0-parameter is partly conventional".**
Every "0-parameter" number is **reclassified** into its true provenance
(``provenance_table``), the audit's **$8+0+3+1+5+6+1$** split over 24 rows: **algebraic
0-param** (a property of the ``Cl(10)`` root: the 16/16-bar, the 45, the
5/3, the anomalies, the 32-level spectrum, ``det U = 1``, the **16/16-bar
spectral mirror**, the **closed clock angle** — the closure axiom (V131)
forces the order-8 root ``g_10 = pi/4``), **free 1-param** (empty:
``g10_status`` shows the UV-finiteness / 2-loop Landau-pole condition fixes
**nothing** — the pole lies below ``M_P`` for **every** ``g_10``; the
selection is algebraic), **convention** (``g_0 = 1/8``, the
``g_10 = 1/n_qubits = 1/5`` reading with its ``alpha_10 = g^2`` coupling
anchor, and the ladder base ``1/(16 pi)`` —
mapping / scale-identification choices, not predictions), **coincidence**
(the ``v/E_P = 2.02e-17`` 0.18% match), **named residual** (``eta_CP``,
``Lambda_obs``, the sub-Planck DM scale (the rung-10 ladder value, carried
by the cell singlet B), the DM portal strength ``lambda_B``, the exact
``m_nu``), **fitted coefficient** (the twelve numbers of the Table 6 fitted
class: the loop clock ``g_obs = 0.1999980`` — booked fitted since the V233
arbitration, its non-derivability from the root a theorem of the reserved
classes (V309, V326) — the desert thresholds, the six ``kappa_s``,
``n_2step = 3``, the DM-mass coefficient ``1.18``, the pair portal
``lambda_*``), and **conditional** (the hierarchy coefficient
``2.710 = 3(16 pi F_208)^10``, conditional on the node-set choice).

**Point 2 — "the hierarchy is fixed, not explained".**
The **smallest** level of the 32-level QCA spectrum at the root clock is
``pi/4 = 0.785`` (a mass ``0.785 * M_P ~ 9.6e18 GeV``) — the **only**
sub-Planck scale in the single-step 0-parameter structure.  It is **not**
``2.02e-17``: the
**exponent** ``-17`` is **not** a 0-parameter consequence of the single
step.  The 0-parameter content is the **existence** of the hierarchy (two
widely separated scales); the **exponent** is the named residual (O8).  The
0.18% match is re-labelled a **2-digit coincidence**.

**Point 3 — "baryogenesis is added, not produced" (the weakest point).**
The **0-parameter CP-odd theorem**: **every** CP-odd (T-odd) invariant of
the 0-parameter QCA **vanishes** — ``Im Tr`` of every 2- and 3-rotation
product (full), a large sample of 4-rotation products, **every**
B-violating ``h_i R_ab h_j R_cd`` sandwich (full), and the CP-odd
invariant of the **rank-1, real** democratic 0-parameter Yukawa.  The
0-parameter QCA is **CP-conserving**, **B-conserving** (the 45 commute with
the pseudoscalar), and **in equilibrium** (unitary) — it satisfies **none**
of the three **Sakharov** conditions.  **Therefore it cannot generate the
baryon asymmetry — a 0-parameter theorem** (not an assertion).  The
asymmetry **requires** the **named** CP-violating phase ``eta_CP``; the
0-parameter **structure that enables** it is the **2-step seesaw with at
least two inequivalent heavy N** (a single rank-1 Yukawa has **no**
CP-odd invariant).  The GUT-leptogenesis **kinematics** (0-parameter) then
gives the **right order** of ``Delta B``.

**Point 4 — "the 0-parameter dark matter does not match the observed".**
The 0-parameter abundance gap is re-computed (``Omega_0param/Omega_obs
~ 1e-52``), **and** the key fact added: the **only** sub-Planck scale in
the 0-parameter structure is the smallest spectral gap
(``0.0115 * M_P ~ 1e17 GeV``) — **~14 orders of magnitude above** the
observed DM (``~1e2-1e3 GeV``).  The observed DM **requires a sub-Planck
scale not in the 0-parameter structure** — a **quantified**, honest failure
(O7).

**The unifying 10-to-1 statement.**  The four "failures" are **honest,
quantified residuals** of the one 0-parameter object; the **genuinely
0-parameter** content (the algebraic invariants) **holds** and is computed
here; the **exact values** are the **named** projections.  The audit turns
four hand-wavy "honest limits" into **four computed theorems + a provenance
table** — the grand-physicist contribution.

The **verdict**: ``pass`` when the **genuinely 0-parameter invariants hold**
(the CP-odd theorem, the 16/16-bar spectral mirror, ``det U = 1``, the
2-route root) **and** the four limits are **quantified** (the hierarchy
exponent not in the spectrum, the clock closed by the root, the DM gap, the
provenance table complete).  ``residual`` is **True** (it names the exact
values the theory does not fix).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import parameter_audit

TOL = 1e-12


def main():
    t0 = time.time()
    r = parameter_audit.verify()
    kr = r["two_route"]
    sp = r["spectral_mirror"]
    cp = r["cp_odd_theorem"]
    g10 = r["g10_status"]
    hi = r["hierarchy_status"]
    dm = r["dm_status"]
    tb = r["provenance_table"]
    ok = bool(r["ok"])
    verdict = {
        "id": "V52_parameter_audit",
        "tier": 6,
        "property": "P5,P6",
        "residual": True,
        "pass": ok,
        "status": ("the_honest_parameter_audit: the '0-parameter' label is reclassified "
                  "(the 8+0+3+1+5+6+1 split over 24 rows: algebraic 0-param / free 1-param / "
                  "convention / coincidence / named residual / fitted coefficient / "
                  "conditional; the loop clock g_obs = 0.1999980 is booked fitted - no "
                  "closed clock passes either empirical contact (V233) and its "
                  "non-derivability from the root is a theorem of the reserved classes "
                  "(V309, V326)); the 0-parameter CP-odd invariants all "
                  "vanish (the QCA satisfies none of the 3 Sakharov conditions and CANNOT "
                  "generate the baryon asymmetry - a theorem); the 16/16-bar spectra are "
                  "exact CPT mirrors (antimatter = time-reversed matter; holds for ALL "
                  "g_10); the clock is fixed by closure (order-8 root) on the internal "
                  "sector, but the loop-hierarchy reading clock is a FIT (V233); the 1/5 "
                  "reading is a "
                  "convention and M_GUT=M_P is abandoned by V64; the hierarchy "
                  "exponent and the DM scale are NOT in the 0-parameter "
                  "single-step spectrum (quantified residuals O8/O7)"),
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "two_route_root": {
                "kingdon_available": kr["kingdon_available"],
                "kingdon_sq_minus_1": kr["kingdon_sq_minus_1"],
                "kingdon_anticommute": kr["kingdon_anticommute"],
                "kingdon_w2_scalar": kr["kingdon_w2_scalar"],
                "pauli_sq_minus_I": kr["pauli_sq_minus_I"],
                "pauli_anticommute": kr["pauli_anticommute"],
                "pauli_w2_plus_I": kr["pauli_w2_plus_I"],
                "two_route_agree": kr["two_route_agree"],
            },
            "spectral_mirror": {
                "g10": sp["g10"],
                "w_pure_imaginary": sp["w_pure_imaginary"],
                "commute_Uw": sp["commute_Uw"],
                "preserves_chiral_split": sp["preserves_chiral_split"],
                "mirror_resid": sp["mirror_resid"],
                "mirror_ok": sp["mirror_ok"],
                "mirror_resid_all_g10": sp["mirror_resid_all_g10"],
                "mirror_holds_for_all_g10": sp["mirror_holds_for_all_g10"],
                "g10_grid": sp["g10_grid"],
                "zp_16": sp["zp_16"],
                "zp_16bar": sp["zp_16bar"],
                "det_U_minus_1": sp["det_U_minus_1"],
                "sl32": sp["sl32"],
                "omega_min": sp["omega_min"],
                "omega_min_mass_GeV": sp["omega_min_mass_GeV"],
            },
            "cp_odd_theorem": {
                "max_cpodd_2rot": cp["max_cpodd_2rot"],
                "max_cpodd_3rot": cp["max_cpodd_3rot"],
                "max_cpodd_4rot_sample": cp["max_cpodd_4rot_sample"],
                "max_cpodd_b_sandwich": cp["max_cpodd_b_sandwich"],
                "cpodd_democratic_yukawa": cp["cpodd_democratic_yukawa"],
                "n4_sample": cp["n4_sample"],
                "all_cp_odd_zero": cp["all_cp_odd_zero"],
                "sakharov": cp["sakharov"],
            },
            "g10_status": {
                "g10_0param": g10["g10_0param"],
                "uv_finiteness_fixes_g10": g10["uv_finiteness_fixes_g10"],
                "pole_2loop_below_MP_for_all_g10": g10["pole_2loop_below_MP_for_all_g10"],
                "g10_is_free_t_hooft_angle": g10["g10_is_free_t_hooft_angle"],
            },
            "hierarchy_status": {
                "v_over_EP_observed": hi["v_over_EP_observed"],
                "v_over_EP_claimed_0param": hi["v_over_EP_claimed_0param"],
                "match_rel_err": hi["match_rel_err"],
                "match_is_2digit_coincidence": hi["match_is_2digit_coincidence"],
                "omega_min": hi["omega_min"],
                "omega_min_mass_GeV": hi["omega_min_mass_GeV"],
                "spectrum_exponent_matches_observed": hi["spectrum_exponent_matches_observed"],
                "hierarchy_exists_0param": hi["hierarchy_exists_0param"],
            },
            "dm_status": {
                "soliton_mass_GeV": dm["soliton_mass_GeV"],
                "abundance_gap_0param_over_obs": dm["abundance_gap_0param_over_obs"],
                "omega_min_mass_GeV": dm["omega_min_mass_GeV"],
                "dm_observed_range_GeV": dm["dm_observed_range_GeV"],
                "smallest_sub_planck_scale_above_observed_DM_orders":
                    dm["smallest_sub_planck_scale_above_observed_DM_orders"],
                "observed_DM_not_0param": dm["observed_DM_not_0param"],
            },
            "provenance_table": tb,
            "invariants_ok": r["invariants_ok"],
            "limits_ok": r["limits_ok"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v52_parameter_audit.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V52 - the honest parameter audit (the solution to the four limits of the theory)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  2-route root (kingdon symbolic vs Pauli GPU):  {kr['two_route_agree']}")
    print(f"  16/16-bar spectral mirror (CPT = spectral conjugation):  "
          f"resid {sp['mirror_resid']:.2e}  (all g_10: {sp['mirror_resid_all_g10']:.2e}, "
          f"ok {sp['mirror_holds_for_all_g10']})")
    print(f"  all CP-odd invariants vanish:  {cp['all_cp_odd_zero']}  "
          f"(the QCA satisfies none of the 3 Sakharov conditions)")
    print(f"  UV-finiteness does not fix the clock (closure does):  {not g10['uv_finiteness_fixes_g10']}")
    print(f"  hierarchy exponent not in the spectrum:  {not hi['spectrum_exponent_matches_observed']}")
    print(f"  observed DM not the 0-parameter (O7, quantified):  {dm['observed_DM_not_0param']}")
    print(f"  provenance classes:  {tb['counts_by_provenance']}")
    print(f"  0-parameter invariants hold:  {r['invariants_ok']}   "
          f"the four limits quantified:  {r['limits_ok']}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  "
          f"(the audit: the 0-parameter invariants hold; the four limits are quantified)")
    print("  JSON -> " + os.path.join(outdir, "v52_parameter_audit.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
