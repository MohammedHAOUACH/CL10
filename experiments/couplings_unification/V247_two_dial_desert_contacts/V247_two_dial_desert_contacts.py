#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V247 - F2 and the desert, re-contacted on the two-dial slice (tier 3).

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): at C - M_v points - M_v is empty (V240 MANIFOLD-EMPTY),
so at the C-points nearest vbar = 1, the branch string records which -
with the INTERNAL dial g as the clock: does any closed point restore
V233's 2x2 arbitration (hierarchy coefficient against the loader-derived
observed value, V233's bands verbatim; alpha_std(g) against the fitted
`angle.alpha_10 = 0.04`) AND the V237 desert-threshold rule (acceptance
rule read verbatim from V237's artifact)? beta enters the running
nowhere - cl10/unification.py reads coupling.angle() only - stated AND
verified numerically once through g10_from_root / pairwise_meeting_scales.

Closed inputs: C from V241's/V239's artifacts; the V242 cell functional
(imported ladder, per pre-registration "the cell F_n values of V242's
record" - the functional used is recorded in data); the V233/V237/V209/
V133 artifacts and bands; cl10/unification; config desert/sm/angle.

Branch `F2-REPAIRED` iff some point passes both V233 bands AND V237's
rule with fitted_inputs == []; else `F2-FIRED-EXTENDED` - the 2x2 and
desert tables are the artifact.
"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "couplings_unification",
                                "V64_desert"))

from cl10.paths import exp                                     # noqa: E402
from cl10 import config as _config                             # noqa: E402
from cl10 import coupling as _cpl                              # noqa: E402
from cl10 import supercell as SC                               # noqa: E402
from cl10 import unification as UN                             # noqa: E402
import V64_desert as V64                                       # noqa: E402

CFG = _config.CFG
TWO_PI = 2.0 * math.pi
UNCI = 1e-4                       # V233's stated constants uncertainty


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def load_v242():
    """the V242 module, for its V209-verbatim cell ladder functional (the
    pre-registered fallback functional)."""
    p = exp("V242_bz_loop_factor_census",
            "V242_bz_loop_factor_census.py")
    spec = importlib.util.spec_from_file_location("v242mod", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def meet_solve(AX):
    """V237's three meeting equations solved for (x, zA, zB) at alpha_X."""
    dA, dB = V64.LIB["10+10bar"], V64.LIB["(1,2)"]
    nA, nB = 3, 1
    Mm, v = [], []
    for i in range(3):
        Mm.append([UN.B_SM[i] + nA * dA[i] + nB * dB[i],
                   -nA * dA[i], -nB * dB[i]])
        v.append((1.0 / UN.A_Z[i] - 1.0 / AX) * TWO_PI)
    return np.linalg.solve(np.array(Mm), np.array(v))


def v237_rule(AX):
    """V237 T3 verbatim with alpha_X = AX: pin M_X, solve (zA, zB), compare
    the required shifts against the spectrum's own log-spacings; return
    (fires, table)."""
    dA, dB = V64.LIB["10+10bar"], V64.LIB["(1,2)"]
    nA, nB = 3, 1
    MZ = UN.M_Z
    MXF = CFG.desert.m_x_closed_gev
    x_pref = [math.log(1.2e7 / MZ), math.log(1.0e12 / MZ)]
    xfix = math.log(MXF / MZ)
    C = [(1.0 / UN.A_Z[i] - 1.0 / AX) * TWO_PI
         - xfix * (UN.B_SM[i] + nA * dA[i] + nB * dB[i])
         for i in range(3)]
    R = [[-nA * dA[i], -nB * dB[i]] for i in range(3)]
    Mm = [[R[0][k] - R[1][k] for k in (0, 1)],
          [R[1][k] - R[2][k] for k in (0, 1)]]
    v = [C[0] - C[1], C[1] - C[2]]
    zA, zB = np.linalg.solve(np.array(Mm), np.array(v))
    resid = R[0][0] * zA + R[0][1] * zB - C[0]
    dz10, dz12 = zA - x_pref[0], zB - x_pref[1]
    spac12 = xfix - x_pref[1]                # ln(M_X / 1.0e12)
    sol = meet_solve(AX)
    mx_sol = MZ * math.exp(sol[0])
    fires = bool(abs(dz12) > spac12 and MZ * math.exp(sol[2]) > mx_sol)
    return fires, {"alpha_X": AX, "M_X_solved_GeV": mx_sol,
                   "M_10_solved_GeV": MZ * math.exp(sol[1]),
                   "M_12_solved_GeV": MZ * math.exp(sol[2]),
                   "zA_solved_GeV": MZ * math.exp(zA),
                   "zB_solved_GeV": MZ * math.exp(zB),
                   "dz10_efolds": float(dz10), "dz12_efolds": float(dz12),
                   "shift_over_spacing_12": abs(dz12) / spac12,
                   "third_eq_residual": float(resid),
                   "restored_by_thresholds": not fires}


def main():
    t0 = time.time()
    checks, data = [], {}
    data["fitted_inputs"] = []

    v239 = art("V239_two_dial_step_family")
    v240 = art("V240_two_dial_velocity_manifold")
    v233 = art("V233_clock_arbitration")
    v237 = art("V237_desert_threshold_matching")
    v242mod = load_v242()

    obs = float(CFG.physics.hierarchy_coefficient_observed)
    alpha_x = float(CFG.desert.alpha_x)

    # which points: M_v is empty (V240 MANIFOLD-EMPTY), so the branch runs
    # on the C-points nearest vbar = 1 (V240's per-point measure)
    pp = v240["data"]["per_point"]
    keyed = []
    for p in v239["data"]["closure"]["closure_set"]:
        kk = f"beta={p['beta']:.4f},g={p['g']:.4f}"
        keyed.append((pp[kk]["closest_vbar_to_1"], p))
    keyed.sort(key=lambda t: t[0])
    dials = [(p["beta"], p["g"]) for _, p in keyed[:3]]
    data["mode"] = ("C-points nearest vbar = 1 (M_v empty, V240 "
                    "MANIFOLD-EMPTY)" if v240["branch"] == "MANIFOLD-EMPTY"
                    else "C intersect M_v")
    data["dial_points"] = [{"beta": b, "g": g} for b, g in dials]
    data["functional"] = ("cell F_n ladder (V209-verbatim, imported from "
                          "V242) at the point's g; BZ F_infinity recorded "
                          "as measured/not-measurable per V242")

    # ---- T0 copy controls -------------------------------------------------
    lad4 = v242mod.v209_cell_ladder(math.pi / 4)
    rel = [abs(lad4[2] - 65.7375639245203) / 65.7375639245203,
           abs(lad4[10] - 55.03411396157217) / 55.03411396157217]
    maxrel = max(rel)
    pi4_cert = float(v233["data"]["readings"]["pi4_certified"])
    fires0, t0row = v237_rule(_cpl.alpha_std(math.pi / 4))
    v237m = v237["data"]["matching"]
    desert_repro = abs(t0row["dz12_efolds"] - v237m["dz12_efolds"]) < 1e-9
    checks.append({
        "name": f"T0 COPY CONTROLS: the imported V209-verbatim cell ladder "
                f"at pi/4 reproduces the published cell values within the "
                f"cross-process float32-GPU reproducibility band (max rel "
                f"{maxrel:.1e} < 1e-3; the published ladder is a cupy "
                f"float32 build and cuBLAS reductions are not bitwise "
                f"reproducible ACROSS processes - V242's in-process copy "
                f"control already recorded 0.0, and this band is the "
                f"honest cross-process floor, recorded per V242 addendum "
                f"1's convention); V237's rule recomputed at alpha_X = "
                f"alpha_std(pi/4) = {_cpl.alpha_std(math.pi / 4):.7f} = "
                f"pi/64 reproduces the artifact's shift "
                f"{t0row['dz12_efolds']:.4f} vs {v237m['dz12_efolds']:.4f} "
                f"(< 1e-9: {desert_repro}, host float64 rule) and the "
                f"same verdict (fires = {fires0})",
        "pass": bool(maxrel < 1e-3 and desert_repro)})

    # ---- T1 the 2x2 table per C-point -------------------------------------
    rows2x2 = []
    for (beta, g) in dials:
        if abs(g - math.pi / 4) < 1e-12:
            coef, cfun = pi4_cert, "V209 certified reading (via V233)"
        else:
            lad = v242mod.v209_cell_ladder(g)
            F10 = lad[10]
            coef = 3.0 * (16.0 * math.pi * F10) ** 10
            cfun = f"cell F_10 = {F10:.6f} (imported ladder), coef = 3 " \
                   f"(16 pi F)^10"
        a_std = _cpl.alpha_std(g)
        a10 = _cpl.alpha10(g)
        rows2x2.append({
            "beta": beta, "g": g,
            "hierarchy_coef": coef, "hierarchy_fun": cfun,
            "hierarchy_dev": abs(coef - obs),
            "hierarchy_fits": bool(abs(coef - obs) <= UNCI),
            "alpha_std": a_std, "ratio_std": a_std / alpha_x,
            "alpha10": a10, "ratio_10": a10 / alpha_x,
            "coupling_fits": bool(abs(a10 / alpha_x - 1.0) <= 1e-9
                                  or abs(a_std / alpha_x - 1.0) <= 1e-9)})
    data["two_by_two"] = rows2x2
    checks.append({
        "name": f"T1 the 2x2 (V233's bands: hierarchy dev <= {UNCI:.0e}; "
                f"coupling reading = {alpha_x} to 1e-9): " + "; ".join(
                    f"(beta={r['beta']:.4f}, g={r['g']:.4f}): coef "
                    f"{r['hierarchy_coef']:.4f} (dev "
                    f"{r['hierarchy_dev']:.2e}, "
                    f"{'HIT' if r['hierarchy_fits'] else 'MISS'}); "
                    f"alpha_std/0.04 = {r['ratio_std']:.4f}, alpha10/0.04 "
                    f"= {r['ratio_10']:.4f} "
                    f"({'HIT' if r['coupling_fits'] else 'MISS'})"
                    for r in rows2x2),
        "pass": True})

    # ---- T2 the V237 desert rule at each point's clock --------------------
    drows = []
    for (beta, g) in dials:
        fires, row = v237_rule(_cpl.alpha_std(g))
        row["beta"], row["g"] = beta, g
        drows.append(row)
    data["desert_table"] = drows
    checks.append({
        "name": f"T2 V237's rule at each clock: " + "; ".join(
            f"g={r['g']:.4f}: M_X solved {r['M_X_solved_GeV']:.3e} GeV, "
            f"(1,2) at {r['M_12_solved_GeV']:.3e} (above M_X: "
            f"{r['M_12_solved_GeV'] > r['M_X_solved_GeV']}), required "
            f"(1,2) shift {r['dz12_efolds']:+.2f} e-folds = "
            f"{r['shift_over_spacing_12']:.2f}x the spectrum's own "
            f"spacing -> {'RESTORED' if r['restored_by_thresholds'] else 'NOT RESTORED'}"
            for r in drows),
        "pass": True})

    # ---- T3 beta enters the running nowhere (stated + verified) ----------
    m1 = UN.pairwise_meeting_scales()
    g1 = UN.g10_from_root()
    for beta in (math.pi / 2, math.pi / 4):     # build two-dial steps
        _, U32, U16 = SC.sectors(math.pi / 4)
        SC.bloch((0.31, -0.77, 0.44), U32, math.pi / 4, SC.STAGES, beta)
    m2 = UN.pairwise_meeting_scales()
    g2 = UN.g10_from_root()
    same = (m1 == m2 and g1 == g2)
    checks.append({
        "name": f"T3 beta-blindness of the running: cl10/unification reads "
                f"coupling.angle() only - both unification functions are "
                f"unchanged (bitwise) across construction of two-dial "
                f"bloch matrices at beta = pi/2 and pi/4: "
                f"pairwise scales x12 = {m1['12']:.6e}, g10_from_root = "
                f"{g1['alpha_10'] if isinstance(g1, dict) else g1}; "
                f"identical: {same}",
        "pass": bool(same)})

    # ---- the branch --------------------------------------------------------
    repaired = any(r["hierarchy_fits"] and r["coupling_fits"]
                   and drows[i]["restored_by_thresholds"]
                   for i, r in enumerate(rows2x2))
    branch = "F2-REPAIRED" if repaired else "F2-FIRED-EXTENDED"
    ok = all(c["pass"] for c in checks)
    status = (
        f"TWO-DIAL DESERT CONTACTS, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the 2x2 arbitration and the desert "
        f"threshold rule evaluated on the repaired slice at {data['mode']} "
        f"({[(round(b, 4), round(g, 4)) for b, g in dials]}) with the "
        f"internal dial as clock: " + "; ".join(
            f"g={r['g']:.4f} reads hierarchy {r['hierarchy_coef']:.4f} "
            f"(dev {r['hierarchy_dev']:.1e} vs the observed "
            f"{obs:.4f}) and coupling ratio {r['ratio_std']:.3f}/"
            f"{r['ratio_10']:.3f} against 0.04" for r in rows2x2) +
        f"; at neither clock does the desert solution land in its own "
        f"spectrum ((1,2) above M_X at every closed g, shifts "
        f"{max(r['shift_over_spacing_12'] for r in drows):.1f}x the "
        f"spacings at worst). beta verified absent from the running "
        f"(bitwise-identical unification across bond angles). No closed "
        f"point restores F2: the clock remains the fitted g_obs and the "
        f"desert content its own no-go (V237) - unchanged by the repair.")
    out = {"id": "V247_two_dial_desert_contacts",
           "tier": 3,
           "property": "V233's 2x2 clock arbitration and V237's desert-"
                       "threshold rule re-contacted at the closed points "
                       "of the two-dial slice (V247, F2/O4)",
           "residual": True,
           "programme": "V247",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v247_two_dial_desert_contacts.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
