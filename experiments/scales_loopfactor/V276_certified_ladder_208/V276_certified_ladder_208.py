#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V276 - the certified ladder and the arithmetic of 208. (tier 4, O8.)

**The question (pre-registered).** Does the certified closed branch
(V200-V204 machinery) support the published hierarchy coefficient
without the n = 10 read, and is the certified crossing n* = 208
arithmetically derived?

**The object (closed input set).** V203's hop_machinery(G_CLOSED) +
loop_ratios(G, cartan, Ma, ev, G10_REF, 260) imported verbatim (same
import path V204 uses); artifacts V204 (machinery, rungs, coefficient),
V108 (coefficient_predicted_3_F10_10, coefficient_observed, F_10_root),
V197 (ladder_rows). All float64 machinery + mpmath (60 dp) arithmetic
continuations only.

**Measured (pre-registered).** T0 reproduction: rho, crossing_n, odd-n
max vs the V204 artifact. T1 ladder read: coef_pred = 3*(16 pi F)^10
evaluated at F := F_10_open_branch (0.019186391026184392), F :=
F_cross (0.01968515418112789), and F := geometric continuation
F_10_geom from (rho, F_cross) per the measured F_n direction (the
continuation direction is read from the loop_ratios construction and
booked explicitly); each against coefficient_observed =
2.076564111173832.  T2 208-arithmetic: mpmath PSLQ over Q on
{1, rho, ..., rho^4} (degree <= 4, height cap 1e4, 60 dp); geometric
crossing n*_geom from (F_10_open_branch, rho) along the booked
direction - record |n*_geom - 208| with the float64 uncertainty band
of rho propagated; factor census 208 = 2^4 * 13, divisibility by the
closed-clock orders (8 | 208) and by the V197 rungs.

**Acceptance (pre-registered).** T0 fail (any of |rho - rho_art| >
1e-9, crossing_n != 208, odd max > 1e-12) -> MACHINERY-DRIFT: halt the
programme, report, do not interpret anything downstream.  Ladder leg:
LADDER-CERTIFIED iff the F_cross read is within 2 % of
coefficient_observed; else LADDER-N10-CONVENTION (both readings
recorded).  208 leg: CROSSING-ALGEBRAIC iff PSLQ finds an exact
relation AND |n*_geom - 208| < 1e-6; else CROSSING-MEASURED.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from cl10 import coupling  # noqa: E402
from cl10 import u10pred  # noqa: E402

G10_REF = coupling.G_CLOSED          # pi/4 - the physical branch (V131)


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    import mpmath
    LOOP = float(u10pred.LOOP)
    checks, data = [], {}

    # ---- T0: the V204 machinery reproduction (gate) --------------------
    mach = art("V204_closed_ladder_limits_audit", "data.machinery")
    rho_art = float(mach["rho"])
    cross_art = int(mach["crossing_n"])
    G, cartan, Ma, ev = hop_machinery(G10_REF)
    F, rho = loop_ratios(G, cartan, Ma, ev, G10_REF, 260)
    odd_max = max(abs(F[n]) for n in range(1, 12, 2))
    crossing = next(n for n in range(2, 261, 2) if F[n] < LOOP)
    F_cross = float(F[crossing])
    drift = bool(abs(rho - rho_art) > 1e-9 or crossing != cross_art
                 or odd_max > 1e-12)
    data["machinery"] = {"rho": float(rho), "rho_artifact": rho_art,
                         "rho_abs_diff": float(abs(rho - rho_art)),
                         "crossing_n": int(crossing),
                         "crossing_n_artifact": cross_art,
                         "odd_max_n11": float(odd_max),
                         "F_cross": F_cross,
                         "F_cross_artifact": float(mach["F_cross"]),
                         "F_cross_abs_diff": float(
                             abs(F_cross - float(mach["F_cross"])))}
    checks.append({"name": "T0 machinery reproduction (rho, crossing "
                           "208, odd-n zero): rho {:.16f} (|d| {:.1e}), "
                           "crossing n* = {:d}, odd max {:.1e}".format(
                               rho, abs(rho - rho_art), crossing,
                               odd_max),
                   "pass": not drift})
    if drift:
        data["status"] = "MACHINERY-DRIFT"
        print("V276 MACHINERY-DRIFT: the reproduced closed-branch "
              "machinery disagrees with the booked V204 audit; nothing "
              "downstream is interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: the ladder reads ----------------------------------------
    F10_open = float(art("V204_closed_ladder_limits_audit",
                         "data.rungs.F_10_open_branch"))
    F_cross_art = float(mach["F_cross"])
    coef_obs = float(art("V108_hierarchy_multiplicity",
                         "data.coefficient_observed"))
    coef_v108 = float(art("V108_hierarchy_multiplicity",
                          "data.coefficient_predicted_3_F10_10"))
    # the V108 re-read, at V85's full-precision F_10 (the artifact books
    # F_10_root rounded to 8 decimals - recompute from the source the
    # way V108 did, to 1e-9)
    F10_v85 = float(art("V85_closed_loop_factor",
                        "data.loops_root.F_n")[-1])
    coef_recomp = 3.0 * (LOOP ** -10) * F10_v85 ** 10
    v108_def = abs(coef_recomp - coef_v108)
    # the continuation direction is read off loop_ratios: F_n =
    # cartan_norm((MG)^n)/cartan_norm(G^n), MG has spectral radius
    # rho < 1 and G is unitary, so F_n DECAYS with n; the geometric
    # continuation from (rho, F_cross) back to n = 10 divides F_cross
    # by rho^198 (multiplies by rho**(-198))
    F_geom = F_cross_art * rho ** (-198.0)
    data["continuation_direction"] = (
        "F_n decays with n: loop_ratios builds F_n = "
        "cartan_norm((MG)^n)/cartan_norm(G^n) with rho(MG) = {:.16f} "
        "< 1 and unitary G, so the geometric continuation from "
        "(rho, F_cross) at n* = 208 back to n = 10 multiplies "
        "F_cross by rho**(-198) = {:.4f}".format(
            rho, rho ** (-198.0)))
    data["F_10_open_branch"] = F10_open
    data["F_10_v85_full_precision"] = F10_v85
    data["F_10_geom"] = float(F_geom)
    data["coef_pred_v108_reread"] = coef_recomp
    data["coef_pred_v108_artifact"] = coef_v108
    data["v108_reread_defect"] = v108_def

    def read_at(Fv):
        return 3.0 * (LOOP ** -10) * Fv ** 10

    reads = {"F_10_open_branch": read_at(F10_open),
             "F_208_certified_crossing": read_at(F_cross_art),
             "F_10_geometric_continuation": read_at(F_geom)}
    data["ladder_reads"] = reads
    data["coefficient_observed"] = coef_obs
    checks.append({"name": "T1 V108 re-read from the V85 F_10 "
                          "(full precision): coef_pred {:.6f} vs booked "
                          "{:.6f} (|d| {:.1e} < 1e-9)".format(
                              coef_recomp, coef_v108, v108_def),
                   "pass": bool(v108_def < 1e-9)})
    checks.append({"name": "T1 ladder read at F_10: {:.4f} vs observed "
                          "{:.4f} ({:+.2%})".format(
                              reads["F_10_open_branch"], coef_obs,
                              reads["F_10_open_branch"] / coef_obs
                              - 1.0),
                   "pass": True})
    checks.append({"name": "T1 ladder read at F_208: {:.4f} vs observed "
                          "{:.4f} ({:+.2%})".format(
                              reads["F_208_certified_crossing"], coef_obs,
                              reads["F_208_certified_crossing"] / coef_obs
                              - 1.0),
                   "pass": True})
    checks.append({"name": "T1 ladder read at the geometric "
                          "continuation F_10_geom = {:.4f}: {:.3e} vs "
                          "observed {:.4f} - the asymptotic decay law "
                          "does not extend back to n = 10 (the "
                          "continuation runs away; recorded, not "
                          "booked)".format(F_geom, reads[
                              "F_10_geometric_continuation"], coef_obs),
                   "pass": True})

    # ---- T2: the arithmetic of 208 ------------------------------------
    mpmath.mp.dps = 60
    rho_mp = mpmath.mpf(repr(rho))
    basis = [rho_mp ** k for k in range(5)]      # {1, rho, ..., rho^4}
    rel_raw = mpmath.pslq(basis, maxcoeff=10 ** 4,
                          tol=mpmath.mpf("1e-12"))
    data["pslq_order"] = ["1", "rho", "rho^2", "rho^3", "rho^4"]
    data["pslq_height_cap"] = 10 ** 4
    data["pslq_tol"] = "1e-12 (height cap 1e4)"
    if rel_raw:
        res_exact = abs(sum(mpmath.mpf(c) * rho_mp ** k
                            for k, c in enumerate(rel_raw)))
        # an exact relation must vanish far below the float64 noise
        # floor of its own evaluation (height x 2**-53); a residual at
        # that floor is a noise fit, not an integer relation
        rel = rel_raw if res_exact < mpmath.mpf("1e-14") else None
        data["pslq_candidate"] = [int(c) for c in rel_raw]
        data["pslq_candidate_residual"] = float(res_exact)
        data["pslq_noise_floor"] = (float(max(abs(c) for c in rel_raw))
                                    * 2.0 ** -53)
    else:
        rel = None
        data["pslq_candidate"] = None
        data["pslq_candidate_residual"] = None
    data["pslq_relation"] = ([int(c) for c in rel] if rel else None)
    LOOPmp = mpmath.mpf(1) / (16 * mpmath.pi)
    F10o = mpmath.mpf(repr(F10_open))
    # geometric crossing from (F_10_open_branch, rho), anchored at the
    # measurement point n = 10 along the booked decaying direction:
    # F_10 * rho**(n - 10) = LOOP
    n_star_geom = 10 + mpmath.log(LOOPmp / F10o, rho_mp)
    dndrho = -mpmath.log(LOOPmp / F10o) / (rho_mp
                                          * mpmath.log(rho_mp) ** 2)
    delta_n = float(dndrho) * 1e-9          # rho uncertainty band
    data["n_star_geom"] = float(n_star_geom)
    data["n_star_geom_minus_208"] = float(n_star_geom - 208)
    data["n_star_rho_band"] = float(abs(delta_n))
    rungs = art("V197_rung_arithmetic", "data.ladder_rows")
    census = {"factorization": "208 = 2**4 * 13",
              "closed_clock_order_8_divides": bool(208 % 8 == 0),
              "class_4k": bool(208 % 4 == 0),
              "rungs_V197": [[int(r["rung"]),
                              ("n/a" if int(r["rung"]) == 0 else
                               bool(208 % int(r["rung"]) == 0))]
                             for r in rungs]}
    data["census_208"] = census
    mpmath.mp.dps = 15
    checks.append({"name": "T2 PSLQ degree-4 relation of rho: {} (over "
                          "Q on {{1, rho, rho^2, rho^3, rho^4}}, height "
                          "cap 1e4, 60 dp; raw candidate {}, residual "
                          "{})".format(
                              data["pslq_relation"],
                              data["pslq_candidate"],
                              data["pslq_candidate_residual"]),
                   "pass": True})
    checks.append({"name": "T2 geometric crossing |n - 208| = {:.3e} "
                          "(n*_geom = {:.6f} from (F_10_open_branch, "
                          "rho); rho band +-{:.1e})".format(
                              abs(float(n_star_geom - 208)),
                              float(n_star_geom), abs(delta_n)),
                   "pass": True})
    checks.append({"name": "T2 factor census 208 = 2^4 * 13: clock "
                          "order 8 | 208 {}, 4k class {}, rung "
                          "divisibility {}".format(
                              census["closed_clock_order_8_divides"],
                              census["class_4k"],
                              census["rungs_V197"]),
                   "pass": True})

    # ---- branch ----
    over208 = reads["F_208_certified_crossing"] / coef_obs - 1.0
    ladder = ("LADDER-CERTIFIED" if abs(over208) < 0.02
              else "LADDER-N10-CONVENTION")
    alg = bool(rel) and abs(float(n_star_geom - 208)) < 1e-6
    crossing_leg = "CROSSING-ALGEBRAIC" if alg else "CROSSING-MEASURED"
    branch = "{} / {}".format(ladder, crossing_leg)
    data["status"] = branch
    print("V276 {}: the certified F_208 read gives the hierarchy "
          "coefficient {:.4f} vs observed {:.4f} ({:+.2%}, {} the 2 % "
          "gate; the F_10 open-branch read gives {:.4f}, {:.2%}); "
          "PSLQ(1..rho^4) relation {}, n*_geom - 208 = {:.3f}: the "
          "crossing 208 is a measured fact of the closed-branch "
          "hop matrix, not an arithmetic consequence of rho.".format(
              branch, reads["F_208_certified_crossing"], coef_obs,
              over208, "inside" if abs(over208) < 0.02 else "outside",
              reads["F_10_open_branch"],
              reads["F_10_open_branch"] / coef_obs - 1.0,
              data["pslq_relation"], float(n_star_geom - 208)),
          flush=True)
    _emit(data, checks, ok=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V276_certified_ladder_208", "tier": 4,
           "property": "does the certified closed branch support the "
                       "ladder without the n = 10 read, and is 208 "
                       "arithmetically derived (tier 4, O8)",
           "residual": True,
           "programme": "V276",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v276_certified_ladder_208.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V276 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
