#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V324 - the closest-theory ledger, version 2 (tier 1, hard,
theory-hunt; the arbiter extended with the two supplier slots).

**The question (pre-registered).** With the joint-centre rung-17
multipliers (V323 band hits, density supplier-grade 0.464 as booked)
and the V319 zero-fitted rung-10 weights (density-DISQUALIFIED
1.45/1.46 as booked), does any zero-fitted row of the V318
five-contact ledger pass all five contacts, and what is the new
closest-row distance?

**Rows.** All 48 V318 ledger rows (stored coef, S17/3, provenance,
zero_fitted) x chi-slot {1, 255/247, 117/109, 90/91} x W-slot
{1, 1+LOOP^3, exp(+LOOP^3 n/17)} = 576 rows.

**Conventions (fixed pre-run).** coef' = coef_row x W^10; rung-17
leg' = (S17row/3) x W(17) x chi; C2 recomputed as membership of
M_P LOOP^10 W(10) in physics.wimp_window_gev (W(10) deviates by
< 1e-5 - membership is stable, booked by the recomputation); C1, C5
carried. Bands V318's exactly.

**Density qualifiers (pre-fixed).** A row passing all five is a
CONSISTENT candidate formally (LEDGER2-CONSISTENT-CANDIDATE); its
SUCCESS TIER (the only paper-relevant one) additionally requires
every supplied pass to be supplier-grade. The chi slot is
supplier-grade (V323 note); the W slot is density-disqualified
(V319 note). Since no V318 row passed all five, any candidate's C4
pass is W-supplied hence density-disqualified: booked as
"consistent with density-disqualified coefficient supply" - a
bound, not a discovery, no paper provenance change. The §17 named
distance token updates to the closest-row leaf either way.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402

from cl10 import config as _config                        # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
TOL_SCALE = float(CFG.pillar_surgery.ledger_scale_rel_tol)
M_ATM = float(CFG.neutrino.m_atm_ev)
WIMP = CFG.physics.wimp_window_gev


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    t0 = time.time()
    checks, data = [], {}
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    DM_BASE = float(u10pred.M_PLANCK) * LOOP ** 10

    v318 = json.load(open(os.path.join(
        exp("V318_closest_theory_ledger"), "results",
        "v318_closest_theory_ledger.json")))
    ledger = v318["data"]["ledger"]
    if len(ledger) != 48:
        print("V324 MACHINERY-DRIFT: V318 ledger has %d rows, expected 48"
              % len(ledger), flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- W slot (V319 algebraic weights; re-derived from LOOP) ---------
    import math
    Ws = [("identity", 1.0, 1.0),
          ("1+LOOP^3 (V319 family A, density-disqualified)",
           1.0 + LOOP ** 3, 1.0 + LOOP ** 3),
          ("exp(+LOOP^3 n/17) (V319 family B, density-disqualified)",
           math.exp(LOOP ** 3 * 10.0 / 17.0), math.exp(LOOP ** 3))]
    wA_leaf = art("V319_rung_weight_double_target", "data.members")
    # T0b leaf: recompute V319's C4 norm for family A on the best row and
    # compare to the booked leaf
    best318 = v318["data"]["best_zero_fitted_row"]
    leaf_members = [m for m in wA_leaf if m["family"] == "A"
                    and m["rule"].startswith("1+1*LOOP^3")]
    leaf_norm10 = leaf_members[0]["norm_rung10"]
    leaf_W10 = leaf_members[0]["W10"]

    # ---- chi slot (V323 band hits, supplier-grade density 0.464) -------
    chis = [("identity", 1.0),
            ("255/247 (3 t2/sum t)", 255.0 / 247.0),
            ("117/109 (t3/mean_mw)", 117.0 / 109.0),
            ("90/91 (m-weighted gap)", 90.0 / 91.0)]

    rows = []
    for r0 in ledger:
        for (wn, w10, w17) in Ws:
            for (chin, chi) in chis:
                coef = r0["coef"] * w10 ** 10
                leg17 = r0["S17_over_3_eV"] * w17 * chi
                n1 = 0.0
                n2 = 0.0 if WIMP[0] <= DM_BASE * w10 <= WIMP[1] else 1e9
                n3 = abs(leg17 - M_ATM) / M_ATM / TOL_SCALE
                n4 = abs(coef - COEF_OBS) / BAND_C
                n5 = r0["normalized"][4]
                norms = [n1, n2, n3, n4, n5]
                rows.append({
                    "clock": r0["clock"], "weight": r0["weight"],
                    "selection": r0["selection"],
                    "W": wn, "W10": w10, "W17": w17,
                    "chi": chin, "coef": coef, "leg17_eV": leg17,
                    "normalized": norms,
                    "worst_normalized": max(norms),
                    "zero_fitted": r0["zero_fitted"],
                    "provenance": r0["provenance"]
                                  + ["algebraic (W slot, V319)",
                                     "algebraic (chi slot, V323)"],
                })
    rows.sort(key=lambda r: (r["worst_normalized"], r["clock"], r["W"],
                             r["chi"]))
    data["n_rows"] = len(rows)

    # ---- T0a: identity re-score reproduces V318 -------------------------
    ident = [r for r in rows if r["W"].startswith("identity")
             and r["chi"] == "identity"]
    bzf_leaf = best318["normalized"]
    mine = next(r for r in ident if r["clock"] == best318["clock"]
                and r["weight"] == best318["weight"])
    t0a = max(abs(a - b) for a, b in zip(mine["normalized"], bzf_leaf))
    checks.append({"name": "T0a identity re-score of the V318 best "
                   "zero-fitted row reproduces the booked normalized "
                   "quintuple to {:.1e} < 1e-9".format(t0a),
                   "pass": bool(t0a < 1e-9)})

    # ---- T0b: W = 1+LOOP^3 reproduces V319 leaves -----------------------
    t0bW = abs(Ws[1][1] - float(leaf_W10)) / float(leaf_W10)
    famA = next(r for r in rows if r["W"].startswith("1+LOOP^3")
                and r["chi"] == "identity"
                and r["clock"] == best318["clock"]
                and r["weight"] == best318["weight"])
    t0bN = abs(famA["normalized"][3] - leaf_norm10) / leaf_norm10
    checks.append({"name": "T0b family-A W10 {:.13f} vs V319 leaf {:.13f} "
                   "(rel {:.1e}) and C4 norm {:.9f} vs leaf {:.9f} "
                   "(rel {:.1e} < 1e-6)".format(Ws[1][1], float(leaf_W10),
                                                t0bW, famA["normalized"][3],
                                                leaf_norm10, t0bN),
                   "pass": bool(t0bW < 1e-12 and t0bN < 1e-6)})

    if not all(c["pass"] for c in checks):
        data["status"] = "MACHINERY-DRIFT"
        print("V324 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- branch ----------------------------------------------------------
    cands = [r for r in rows if r["zero_fitted"]
             and r["worst_normalized"] <= 1.0]
    closest = rows[0]
    for r in rows:
        r["all_inside"] = bool(max(r["normalized"]) <= 1.0)
    qualified = [r for r in cands if r["W"].startswith("identity")]
    branch = ("LEDGER2-CONSISTENT-CANDIDATE" if cands
              else "LEDGER2-CLOSEST-DECLARED")
    data["status"] = branch
    data["closest_row"] = closest
    data["n_candidates"] = len(cands)
    data["n_candidates_identity_W"] = len(qualified)
    data["qualifier"] = ("candidates are consistent with "
                         "density-DISQUALIFIED coefficient supply (W slot, "
                         "V319 rule): a bound, not a discovery; no paper "
                         "re-provenance" if cands and not qualified else
                         "no density qualifier needed" if qualified else
                         "no candidate")
    checks.append({"name": "branch {}: {} zero-fitted all-pass rows "
                   "({} without W); closest row worst_normalized "
                   "{:.7f} (clock {}, W {}, chi {}, norms {}); "
                   "qualifier: {}".format(
                       branch, len(cands), len(qualified),
                       closest["worst_normalized"], closest["clock"][:48],
                       closest["W"][:40], closest["chi"],
                       ["{:.4f}".format(v) for v in closest["normalized"]],
                       data["qualifier"]), "pass": True})
    data["top20"] = rows[:20]
    data["wall_time"] = time.time() - t0
    _emit(data, checks, ok=all(c["pass"] for c in checks))
    print("[V324] {} closest {:.7f} candidates {} ({})".format(
        branch, closest["worst_normalized"], len(cands),
        data["qualifier"]), flush=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V324_closest_theory_ledger_v2", "tier": 1,
           "property": "five-contact closest-theory ledger re-scored with "
                       "the joint-centre rung-17 multipliers and the V319 "
                       "rung-10 weights, density qualifiers pre-fixed "
                       "(tier 1, hard, O8/O10)",
           "residual": True, "programme": "theory-hunt V324",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v324_closest_theory_ledger_v2.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V324 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
