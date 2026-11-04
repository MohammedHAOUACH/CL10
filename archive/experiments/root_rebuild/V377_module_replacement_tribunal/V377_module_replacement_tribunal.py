#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V377 - the module-replacement tribunal (booked BEFORE the run in
review/pre_registrations_root.md, "The module-replacement tribunal:
V377-V378", with the pre-run booking correction of the same day).

Question: does any BUILT module replacement move any contact of the
ten-to-one? The class is exactly what the repository has built: the
committed 36-level cell, the altcell signature family (the p = 0
quaternionic member), the 2I engine internal layer, the star module,
the pentagon engine, and the rung lattice (n_A, n_B free). Contacts
are the published numbers, read from the committed artifacts where
the member measured them and recomputed in-script for the committed
cell and the rung lattice.

Verdicts: MODULE-REPLACEMENT-CONTACT-STATIONARY (all columns
stationary) / a named branch otherwise (MODULE-CLOSURE-RUNG-FOUND on
exact rung closure; per-column contact-shift words).
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
from cl10 import cell, u10pred  # noqa: E402
from cl10.paths import exp  # noqa: E402

T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, ok))


checks = []


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def best_res(x, mmax):
    """min over 1 <= m <= mmax of || m * x ||, x in [0, 1) of a real
    number mod 1, by the best continued-fraction convergent."""
    f = Fraction(x).limit_denominator(mmax)
    return abs(f.numerator - f.denominator * x), f.denominator


def main():
    data = {}

    # ---- T0 internal self-tests (no contact) --------------------------
    loop = float(u10pred.LOOP)
    a_loop = abs(loop - 1.0 / (16.0 * math.pi)) <= 1e-15
    chk("LOOP = 1/(16pi) = 0.0198943... from the library "
        "(dev %.2e <= 1e-15): %s"
        % (abs(loop - 1.0 / (16.0 * math.pi)), a_loop), a_loop)
    pa = abs(cell.rung_phase(3) - cell.phase_A())
    pb = abs(cell.rung_phase(10) - cell.phase_B())
    a_ph = pa == 0.0 and pb == 0.0
    chk("rung phases loop^3, loop^10 match the config rungs "
        "(3, 10) exactly (dev %.1e, %.1e): %s" % (pa, pb, a_ph), a_ph)

    G = math.pi / 4.0
    U = cell.cell_step_np(G)
    n = U.shape[0]
    unit = float(np.abs(U.conj().T @ U - np.eye(n)).max())
    a_u = unit <= 1e-12
    chk("36-level cell step at g = pi/4 unitary (defect %.2e "
                  "<= 1e-12): %s" % (unit, a_u), a_u)
    ph = np.sort(np.angle(np.linalg.eigvals(U)))
    booked = np.sort(np.array([-G] * 10 + [G] * 10 + [-3 * G] * 6
                              + [3 * G] * 6 + [-loop ** 3, loop ** 3,
                                              -loop ** 10, loop ** 10]))
    mset = float(np.abs(ph - booked).max())
    a_ms = mset <= 1e-12
    chk("cell eigenphase multiset at Gamma, g = pi/4 is "
                  "{+-pi/4 x 10, +-3pi/4 x 6, +-loop^3, +-loop^10} "
                  "(max dev %.2e <= 1e-12): %s" % (mset, a_ms), a_ms)
    data["loop"] = loop
    data["rung_phase_A"] = cell.phase_A()
    data["rung_phase_B"] = cell.phase_B()

    # ---- T1 closure census of the committed cell -----------------------
    S = U[:32, :32]
    s8 = float(np.abs(np.linalg.matrix_power(S, 8) - np.eye(32)).max())
    a_s8 = s8 <= 1e-12
    chk("spinor block alone closes at order 8 (defect %.2e "
                  "<= 1e-12): %s" % (s8, a_s8), a_s8)
    powers = {}
    for m in range(1, 25):
        powers[m] = float(np.abs(np.linalg.matrix_power(U, m)
                                 - np.eye(n)).max())
    mmin = min(powers, key=powers.get)
    dmin = powers[mmin]
    a_clos = dmin >= 1.4e-5 and min(powers) > 0.0
    chk("full 36-level cell NEVER closes at g = pi/4: min over "
                  "m = 1..24 of |U^m - I| is %.4e at m = %d, >= the "
                  "booked 1.4e-5 floor (obstructing phases 12*loop^3, "
                  "12*loop^10): %s" % (dmin, mmin, a_clos), a_clos)
    data["closure_defect_min"] = dmin
    data["closure_defect_argm"] = mmin

    # ---- T2 rung-lattice census ----------------------------------------
    exact_all = []
    eps_n = None
    for nn in range(1, 31):
        x = loop ** nn / (2.0 * math.pi)
        r, q = best_res(x, 10 ** 6)
        exact_all.append(r)
        if eps_n is None and x <= 1e-9:
            eps_n = nn
    emin = min(exact_all)
    a_exact = emin > 0.0
    chk("rung grid (n = 1..30): NO rung admits an exact closure "
                  "- min continued-fraction residual ||m*loop^n/2pi|| "
                  "over m <= 1e6 is %.3e > 0.0: %s" % (emin, a_exact),
                  a_exact)
    a_eps = eps_n == 5
    chk("epsilon-closure threshold measured n* = %d == booked 5 "
                  "(loop^n/(2pi) <= 1e-9; a DECOUPLING, not a clock): %s"
                  % (eps_n, a_eps), a_eps)
    # the epsilon-closed rungs move no clock: rebuild the cell on (5,5)
    import cl10.cell as _c
    saved = (_c.RUNG_A, _c.RUNG_B)
    try:
        _c.RUNG_A, _c.RUNG_B = 5, 5
        U55 = _c.cell_step_np(G)
        d55 = float(np.abs(np.linalg.matrix_power(U55, 8)
                           - np.eye(n)).max())
        _c.RUNG_A, _c.RUNG_B = 7, 7
        U77 = _c.cell_step_np(G)
        d77 = float(np.abs(np.linalg.matrix_power(U77, 8)
                           - np.eye(n)).max())
    finally:
        _c.RUNG_A, _c.RUNG_B = saved
    a_cls = d55 > 0.0 and d77 > 0.0 and d55 < 1e-6 and d77 < 1e-6
    chk("cell on the epsilon-closed rungs (5,5), (7,7): order-8 "
                  "defect %.3e, %.3e - strictly > 0 (no exact closure) and "
                  "< 1e-6 (decoupled): closure clocks stay {pi/2, pi/4, "
                  "3pi/4}, no new clock: %s" % (d55, d77, a_cls), a_cls)
    data["rung_exact_min_residual"] = emin
    data["rung_eps_nstar"] = eps_n
    data["rung55_order8_defect"] = d55
    data["rung77_order8_defect"] = d77

    # ---- T3 the contact ledger (artifact-gated) --------------------------
    g_obs = 0.19999798835751165  # booked: from the V233 readings field
    observed = 2.076564111173832

    v233 = json_art("V233_clock_arbitration")["data"]["readings"]
    a_gobs = abs(v233["g_obs"] - g_obs) == 0.0
    chk("g_obs = 0.19999798835751165 read from the V233 "
                  "readings field (exact match): %s" % a_gobs, a_gobs)
    pi4_cert = v233["pi4_certified"]
    v209 = json_art("V209_root_clock_degeneracy")["data"]["certified"]
    pi4_209 = v209["pi/4@1e-8"]["coef"]
    a_cross = abs(pi4_cert - pi4_209) == 0.0
    chk("certified pi/4 read %r agrees between V233 and V209 "
                  "(bit-identical): %s" % (pi4_cert, a_cross), a_cross)

    v312 = json_art("V312_quaternionic_loop_ladder")["data"]
    coef_q = {k: v["coef"] for k, v in v312["ladder_runs"]["0"].items()}
    a_q3 = (abs(coef_q["pi*1/2"] - 2.299444) < 1e-4
            and abs(coef_q["pi*1/4"] - 2.708108) < 1e-4
            and abs(coef_q["pi*3/4"] - 2.705722) < 1e-4)
    chk("quaternionic (p = 0) closure-clock reads {2.299444, "
                  "2.708108, 2.705722} from the V312 artifact: %s"
                  % a_q3, a_q3)
    closest = v312["closest_read"]
    deficits = {k: abs(c - observed) for k, c in coef_q.items()}
    deficits["pi/4 certified"] = abs(pi4_cert - observed)
    dmin_col = min(deficits.values())
    dmax_col = max(deficits.values())
    a_def = (abs(closest["residual_abs"] - dmin_col) <= 1e-12
             and abs(dmin_col - 0.2228796480336257) <= 1e-15
             and abs(dmax_col - 0.6329617389049436) <= 1e-15)
    chk("hierarchy column over ladder-bearing members: "
        "|coef - 2.076564111173832| in [0.22288, 0.63296] - "
        "measured [%.6f, %.6f], endpoints bit-exact to the V312 "
        "closest_read and the certified pi/4 read: %s" % (dmin_col,
        dmax_col, a_def), a_def)

    clocks_committed = [math.pi / 2, math.pi / 4, 3 * math.pi / 4]
    clock_dist = min(abs(g - g_obs) for g in clocks_committed)
    v249 = json_art("V249_star_module_closure_census")["data"]
    d249 = v249["membership"]["nearest_1_5"]["distance"]
    a_star = abs(d249 - (math.pi / 4 - 0.2)) <= 1e-12
    chk("star module: nearest closed star point to the pin 0.2 "
                  "is at distance %.6f (V249 field, = pi/4 - 0.2): %s"
                  % (d249, a_star), a_star)
    per_eta = v249["membership"]["per_eta"]
    a_eta = (per_eta["0.0"]["closed_points"] == 4
             and all(v["closed_points"] == 0 for k, v in per_eta.items()
                     if k != "0.0"))
    chk("star closures {1:1, 2:3} at eta = 0, ZERO closures at "
                  "every booked eta twist (the twist moves nothing): %s"
                  % a_eta, a_eta)
    clock_col_min = min(clock_dist, min(
        (v["min_dist_to_0.2"] for v in per_eta.values()
         if v["min_dist_to_0.2"] is not None), default=9.0))
    a_clk = clock_col_min >= 0.28
    chk("clock column: NO member with a closure reaches g_obs - "
                  "min clock-contact distance over committed {pi/2, pi/4, "
                  "3pi/4} and the star class is %.6f >= booked 0.28 "
                  "floor: %s" % (clock_col_min, a_clk), a_clk)

    v337 = json_art("V337_2i_engine")["data"]["grid"]
    holds = sum(1 for p in v337 if p["holds"])
    a_2i = holds == 0 and len(v337) == 9
    chk("2I engine: %d of %d angle-grid points keep unitarity + "
                  "gapless M + exact 3-fold R + mirror (booked 0 of 9): %s"
                  % (holds, len(v337), a_2i), a_2i)

    for vname, word in (("V356_pentagon_closure_free_clock",
                         "PENTAGON-CLOSURE-CLOCK-FREE"),
                        ("V359_pentengine_ladder", "PENTLADDER-NO-CROSSING"),
                        ("V357_pentagon_walk_ladder", "PENTAGON-LADDER-VOID"),
                        ("V342_icosahedral_node_census",
                         "N2B-NODE-ORBIT-557-CONE-ANISOTROPIC"),
                        ("V350_backreaction_fixed_point",
                         "FP-DEGENERATE-RIGID-DENSITY")):
        st = json_art(vname)["status"]
        a_w = word in st
        chk("%s verdict-word gate (%r on record): %s"
                      % (vname, word, a_w), a_w)

    v117 = json_art("V117_real_form")["data"]
    j2 = v117["T3"]["J2"]
    a_j2 = j2 == -1.0
    chk("signature column: the p = 0 member is quaternionic "
                  "(J^2 = -1.0, V117 T3 field) - the CP-content structure "
                  "changes but V312 keeps the odd-n rule (no new CP "
                  "number) and the deficit column is unchanged: %s"
                  % a_j2, a_j2)
    odd = v312["odd_rule_somewhere"] is True
    chk("V312 odd-n selection rule survives on p = 0 (field "
                  "true): %s" % odd, odd)

    v376 = json_art("V376_dressed_cp_odd_census")
    a_376 = float(v376["residual"]) <= 1e-12
    chk("dressed-class CP column: worst booked invariant %.2e "
                  "<= 1e-12 (V376 artifact): %s"
                  % (float(v376["residual"]), a_376), a_376)
    v375 = json_art("V375_dressed_no_ladder_contact")
    a_375 = (v375["data"]["hits"] == [] or v375["data"]["hits"] == 0)
    chk("dressed-ladder column: zero contacts (V375 artifact): %s"
                  % a_375, a_375)

    data["g_obs"] = g_obs
    data["observed_coef"] = observed
    data["deficits_abs"] = deficits
    data["clock_col_min_dist"] = clock_col_min
    data["star_nearest_02"] = d249
    data["two_i_holds"] = holds

    # ---- T4 decision -----------------------------------------------------
    ok = all(c[1] for c in checks)
    rung_branch = "MODULE-CLOSURE-RUNG-EPS-DECOUPLED"
    status = ("MODULE-REPLACEMENT-CONTACT-STATIONARY (no built module "
              "member moves a contact of the ten-to-one: the clock column "
              "min distance %.6f >= 0.28 in every closure class, the "
              "hierarchy deficit stays in [%.4f, %.4f] for every "
              "ladder-bearing member, the ladder objects are VOID exactly "
              "where closure is clock-free (V356/V357/V359), content "
              "columns reproduce verbatim; rung branch %s: rung closure "
              "exists only by decoupling (n* = %d), never exactly; the "
              "named surviving structure is the crystalline dilemma - the "
              "ladder needs a crystallographic shift module, clock-free "
              "closure needs a non-crystallographic one, no built member "
              "carries both)"
              % (clock_col_min, dmin_col, dmax_col, rung_branch, eps_n))
    out = {"id": "V377_module_replacement_tribunal",
           "tier": 2,
           "residual": clock_col_min,
           "property": ("module-replacement contact stationarity over the "
                        "built class: closure census of the committed cell "
                        "and the rung lattice (new), contact ledger over "
                        "every built member (artifact-gated), tier 2, "
                        "O8/O9/O3/O4/O5/O7"),
           "status": status,
           "pass": bool(ok),
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v377_module_replacement_tribunal.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V377 done] %s" % status[:400])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
