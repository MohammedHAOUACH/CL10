#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V327 - is the portal coupling lambda_B an algebraic number?
(tier 3, O7, theory-hunt closure programme.)

**The question (pre-registered).** Is the portal strength
lambda_B = 4.0623088619285984e-14 (V195's booked production leaf,
data.production.lambda_fi - the single named number of the dark-matter
sector, O7) a number of the booked closed class - does a supplied rule
generate it?

**Closed class (booked pre-run, verbatim).**
(i) integer polynomials: mpmath findpoly, degrees 1..6, height <= 20, at
    x = lambda_B at 50 dp (V309's PSLQ convention); a candidate is booked
    only with substitution residual < 1e-14 at 50 dp; candidates with
    c_0 = 0 are rejected (a polynomial divisible by x vanishes near
    x ~ 4e-14 by construction, not by rule - booked guard).
(ii) rational census p/q, q <= 5000, against lambda_B/LOOP^k for integer
    k = 1..90, band 1e-6 on the ratio; a hit is booked only if the rule
    lambda_B = (p/q) LOOP^k survives re-verification at 50 dp with
    LOOP = 1/(16 pi) exact and lambda_B the artifact leaf (the residual
    booked with its size relative to lambda_B). AMENDED by the new
    pre-registration of 2026-11-01 (first-run tautology diagnostic
    booked there): band RELATIVE on the ratio, and the census runs with
    the symmetric height bound max(p, q) <= 5000; the unbounded-hit
    count is booked as the amendment's mechanism leaf. FOUND from this
    class needs the re-verified residual relative to the leaf below
    1e-12 (amendment booking 2); 1e-6-band contacts at the generic
    continued-fraction floor 1/q^2 are booked as contacts, not rules.
(iii) the V326 closed dictionary class: every member c of
    V326 data.dictionary_class.members times LOOP^k, k = 1..90, against
    lambda_B at the 5% band (ledger_scale_rel_tol).

**T0 (pre-registered).** V195 data.production.lambda_fi and
V193_freeze_in data.variants[0].lambda_fi_solution equal bit-for-bit
(booked abs diff, gate == 0.0). Failure: MACHINERY-DRIFT, nothing
interpreted.

**Branch (pre-registered).** LAMBDA-B-ALGEBRAIC-FOUND (any class hits:
the poly / (p/q, k) / (c, k) plus the verified residual booked) /
LAMBDA-B-NO-RELATION-DECLARED (all three classes miss: the closest
contact of each class and its band multiple booked).
"""
import json
import os
import sys

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402

from cl10 import config as _config                        # noqa: E402

CFG = _config.CFG
LOOP = float(CFG.physics.loop)       # 1/(16 pi), the booked loop factor
DEN_MAX = 5000                       # the booked rational census range


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def poly_string(coeffs):
    """An integer-coefficient list [c_0..c_d] as a readable polynomial."""
    terms = []
    for d, c in enumerate(coeffs):
        if c == 0:
            continue
        x = "" if d == 0 else ("x" if d == 1 else "x^%d" % d)
        terms.append(("%+d%s" % (c, x)) if terms else ("%d%s" % (c, x)))
    return " ".join(terms)


def main():
    import mpmath

    dmax = int(CFG.pillar_surgery.pslq_max_degree)
    hcap = int(CFG.pillar_surgery.pslq_max_coeff)
    band = float(CFG.pillar_surgery.clock_band_abs)      # 1e-6 on the ratio
    tol_dict = float(CFG.pillar_surgery.ledger_scale_rel_tol)   # 5% band
    checks, data = [], {}

    lam = float(art("V195_cell_redefinition", "data.production.lambda_fi"))

    # ---- T0 gate: the two sibling freeze-in/production leaves agree ----
    # (the V193 variants are stored as a list; the booked path names the
    #  first variant - the dotted-key art() cannot index, read directly.)
    with open(exp("V193_freeze_in", "results",
                  "v193_freeze_in.json")) as fh:
        v193 = json.load(fh)
    lam193 = float(v193["data"]["variants"][0]["lambda_fi_solution"])
    t0_eq = abs(lam - lam193)
    data["t0"] = {"v195_lambda_fi": lam, "v193_variant0_lambda_fi": lam193,
                  "abs_diff": t0_eq}
    checks.append({"name": "T0 the booked portal value: V195 "
                          "data.production.lambda_fi = {:.17g} and V193 "
                          "data.variants[0].lambda_fi_solution equal "
                          "bit-for-bit (|d| {:.1e}, gate == 0.0)".format(
                              lam, t0_eq),
                   "pass": bool(t0_eq == 0.0)})
    if t0_eq != 0.0:
        data["status"] = "MACHINERY-DRIFT"
        print("V327 MACHINERY-DRIFT: the V195/V193 portal anchors are not "
              "the same number; nothing is interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- class (i): findpoly degrees 1..dmax at 50 dp ------------------
    mpmath.mp.dps = 50
    x = mpmath.mpf(repr(lam))
    noise = mpmath.mpf("1e-14")      # substitution floor (V309 convention)
    found_poly = None
    per_degree = {}
    for d in range(1, dmax + 1):
        try:
            rel = mpmath.findpoly(x, d, maxcoeff=hcap)
        except ValueError:
            # pslq raises when x^d falls below the working tolerance; with
            # x ~ 4e-14 any height-bounded relation with c_0 != 0 would
            # need |c_0| = |sum_{k>=1} c_k x^k| < 1e-11, impossible for a
            # nonzero integer - no relation can hide above such a degree.
            rel = None
        if rel is None:
            per_degree[d] = None
            continue
        coeffs = [int(c) for c in rel]
        res = abs(sum(mpmath.mpf(c) * x ** k for k, c in enumerate(rel)))
        height = max(abs(c) for c in coeffs)
        trivial = bool(coeffs[0] == 0)
        per_degree[d] = {"coeffs": coeffs, "poly": poly_string(coeffs),
                         "residual_50dp": float(res), "height": height,
                         "c0_zero_rejected": trivial}
        if trivial:
            # x^m * Q(x) evaluates ~0 near x ~ 4e-14 by construction:
            # the booked guard rejects it (V309 never needed it, its x
            # was O(1/5)).
            continue
        if res < noise:
            found_poly = {"degree": d, "height": height, "coeffs": coeffs,
                          "poly": poly_string(coeffs),
                          "residual_50dp": float(res)}
    data["findpoly"] = {"degrees": "1..%d" % dmax, "height_cap": hcap,
                        "substitution_noise_floor": float(noise),
                        "per_degree": {str(d): per_degree[d]
                                       for d in per_degree}}
    mpmath.mp.dps = 15
    checks.append({"name": "class (i) findpoly degrees 1..%d, height %d: %s"
                   % (dmax, hcap,
                      "no nontrivial relation" if found_poly is None
                      else "verified {}".format(found_poly)),
                   "pass": True})

    # ---- class (ii): rational census against lambda_B/LOOP^k ----------
    # Amended class (new pre-registration 2026-11-01, first-run
    # diagnostic booked there): band RELATIVE 1e-6 on the ratio, and the
    # SYMMETRIC height bound max(p, q) <= 5000 - an unbounded numerator
    # makes the ratio test a float identity at large loop powers
    # (342,334-hit first-run diagnostic). The unbounded count is booked
    # as the amendment's mechanism leaf.
    qs = np.arange(1, DEN_MAX + 1, dtype=float)
    per_k = []
    hits = []
    n_unbounded_diag = 0
    for k in range(1, 91):
        xk = lam / LOOP ** k
        ps = np.round(xk * qs)
        sel = ps >= 1.0
        if not bool(np.any(sel)):
            continue
        errs = np.abs(ps[sel] / qs[sel] - xk)
        n_unbounded_diag += int(np.count_nonzero(errs < band))
        keep = ps[sel] <= DEN_MAX
        if not bool(np.any(keep)):
            continue
        qv, pv, ev = qs[sel][keep], ps[sel][keep], errs[keep]
        rel = ev / xk
        inside = np.nonzero(rel < band)[0]
        jm = int(np.argmin(rel))
        per_k.append({"k": k, "x_k": xk, "count_inside": int(inside.size),
                      "closest_p": int(pv[jm]), "closest_q": int(qv[jm]),
                      "closest_rel_err": float(rel[jm])})
        for i in inside:
            hits.append({"k": k, "p": int(pv[i]), "q": int(qv[i]),
                         "err_rel": float(rel[i])})
    found_rat = None
    verified_best_contact = None
    if hits:
        # 50-dp re-verification: LOOP = 1/(16 pi) exact, lam the leaf.
        # FOUND only if the re-verified residual RELATIVE to the leaf is
        # < 1e-12 (amendment booking 2: the 1e-6-band contacts sit at
        # the generic continued-fraction floor 1/q^2 and are booked as
        # contacts, not rules).
        mpmath.mp.dps = 50
        best = min(hits, key=lambda h: h["err_rel"])
        lam_mp = mpmath.mpf(repr(lam))
        LOOPM = mpmath.mpf(1) / (16 * mpmath.pi)
        rule = (mpmath.mpf(best["p"]) / mpmath.mpf(best["q"])
                * LOOPM ** best["k"])
        resid = abs(lam_mp - rule)
        relv = float(resid / lam_mp)
        cfg_loop_dev = float(abs(LOOPM - mpmath.mpf(repr(LOOP))) / LOOPM)
        verified_best_contact = dict(best, residual_abs_50dp=float(resid),
                                     residual_rel=relv,
                                     diophantine_floor_1_over_q2=
                                         1.0 / best["q"] ** 2,
                                     cfg_loop_rel_dev_vs_exact=
                                         cfg_loop_dev)
        if relv < 1e-12:
            found_rat = verified_best_contact
        mpmath.mp.dps = 15
    data["rational_census"] = {
        "den_max": DEN_MAX, "band_rel_on_ratio": band,
        "symmetric_height_bound": DEN_MAX,
        "found_gate_residual_rel": 1e-12,
        "amendment_unbounded_p_hit_count_first_run_mechanism":
            n_unbounded_diag,
        "ks": "1..90", "n_hits": len(hits),
        "hits_first5": hits[:5],
        "per_k_closest": per_k,
        "verified_best_contact": verified_best_contact,
        "verified_rule": found_rat}
    checks.append({"name": "class (ii, amended) rationals max(p,q) <= "
                           "{:d} inside relative {:.0e} of lambda_B/LOOP^k "
                           "over k = 1..90: {:d} contacts ({} unbounded-p "
                           "hits are the booked amendment mechanism); best "
                           "verified contact {} -> FOUND only under the "
                           "1e-12 residual gate{} "
                   .format(DEN_MAX, band, len(hits), n_unbounded_diag,
                           verified_best_contact,
                           " (rule supplied)" if found_rat else " (residual "
                           "{:.2e} sits at the 1/q^2 floor - contact, not "
                           "rule)".format(verified_best_contact[
                               "residual_rel"])
                           if (verified_best_contact and not found_rat)
                           else ""),
                   "pass": True})

    # ---- class (iii): the V326 dictionary class x LOOP^k ---------------
    members = [float(s.split("=")[1]) for s in
               art("V326_phase_dictionary_bridge",
                   "data.dictionary_class.members")]
    best3 = None
    contacts3 = []
    for c in members:
        for k in range(1, 91):
            val = c * LOOP ** k
            rel = abs(val / lam - 1.0)
            if best3 is None or rel < best3["rel_to_lambda_B"]:
                best3 = {"c": c, "k": k, "value": val,
                         "rel_to_lambda_B": rel,
                         "band_multiples": rel / tol_dict}
            if rel <= tol_dict:
                contacts3.append({"c": c, "k": k, "value": val,
                                  "rel_to_lambda_B": rel})
    data["dictionary_class"] = {"members": members, "ks": "1..90",
                                "band_rel": tol_dict,
                                "contacts": contacts3,
                                "closest": best3}
    checks.append({"name": "class (iii) {} dictionary members x LOOP^k, "
                           "k = 1..90, {:.0%} band: {} contacts; closest "
                           "c = {} at k = {}: {:.4e} rel = {:.1f} band-"
                           "multiples".format(len(members), tol_dict,
                                              len(contacts3),
                                              best3["c"], best3["k"],
                                              best3["rel_to_lambda_B"],
                                              best3["band_multiples"]),
                   "pass": True})

    # ---- branch ---------------------------------------------------------
    rule = None
    if found_poly is not None:
        rule = {"class": "i", "kind": "integer polynomial",
                **found_poly}
    elif found_rat is not None:
        rule = {"class": "ii", "kind": "rational x loop power",
                "rule": "lambda_B = ({:d}/{:d}) * LOOP^{:d}".format(
                    found_rat["p"], found_rat["q"], found_rat["k"]),
                **found_rat}
    elif contacts3:
        rule = {"class": "iii", "kind": "dictionary member x loop power",
                **contacts3[0]}
    branch = ("LAMBDA-B-ALGEBRAIC-FOUND" if rule is not None
              else "LAMBDA-B-NO-RELATION-DECLARED")
    data["status"] = branch
    data["relation"] = rule
    print("V327 {}: lambda_B = {:.17g}; findpoly 1..{} cap {} -> {}; "
          "rational hits: {}; dictionary closest c = {} k = {} at {:.3e} "
          "rel ({:.1f} band-multiples).".format(
              branch, lam, dmax, hcap,
              found_poly["poly"] if found_poly else "none",
              len(hits), best3["c"], best3["k"], best3["rel_to_lambda_B"],
              best3["band_multiples"]), flush=True)
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V327_algebraicity_lambda_B", "tier": 3,
           "property": "is the portal coupling lambda_B a number of the "
                       "booked algebraic/rational-loop/dictionary closed "
                       "class - could a root rule generate it "
                       "(tier 3, O7)",
           "residual": True,
           "programme": "theory-hunt V327",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v327_algebraicity_lambda_B.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V327 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
