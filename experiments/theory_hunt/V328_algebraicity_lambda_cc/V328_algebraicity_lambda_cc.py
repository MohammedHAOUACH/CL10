#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V328 - is the trace-phase breaking an algebraic number?
(tier 3, O1, theory-hunt closure programme.)

**The question (pre-registered).** Is the trace-phase breaking size
eps_crit = 1.5791367041742974e-120 (V78_lambda_protection booked leaf,
data.eps_crit - the value equivalent to the observed 1e-122 M_P^4, the
named bound of O1) a number of the booked closed class - does a
supplied rule generate it?

**Closed class (booked pre-run, verbatim).**
(i) integer polynomials: mpmath findpoly, degrees 1..6, height <= 20,
    at x = eps_crit/LOOP^70 (rescaled: x is O(0.2); float-safe at
    1e-120), 50 dp, substitution noise floor 1e-14, and the V327
    c_0 = 0 guard verbatim.
(ii) rational census p/q against eps_crit/LOOP^k for integer k = 1..200,
    band 1e-6 RELATIVE on the ratio, symmetric height bound
    max(p, q) <= 5000 (both per the amendments booked in
    review/pre_registrations.md - the V309 DEN_MAX range with the
    numerator bound; ratios beyond the bound contribute no rows);
    FOUND only under the 1e-12 relative-residual gate (amendment 2),
    50-dp re-verified with LOOP = 1/(16 pi) exact and eps_crit the
    artifact leaf.
(iii) the V326 closed dictionary class members c times LOOP^k,
    k = 1..200, against eps_crit at the 5% band (ledger_scale_rel_tol).

**T0 (pre-registered).** V78 data.eps_crit = 1.5791367041742974e-120:
(a) bare loop exponent re-derived, log(eps_crit)/log(LOOP) against the
paper contact 70.42, gate |d| < 0.01; (b) the paper's rounded leaf:
|eps_crit/1.6e-120 - 1| = 0.01304 booked, gate < 0.02. Failure:
MACHINERY-DRIFT, nothing interpreted.

**Branch (pre-registered).** LAMBDA-CC-ALGEBRAIC-FOUND (any class hits
under its booked gate: the poly / (p/q, k) / (c, k) plus the verified
residual booked) / LAMBDA-CC-NO-RELATION-DECLARED (all miss: closest k
and band multiple booked per class).
"""
import json
import math
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
    band = float(CFG.pillar_surgery.clock_band_abs)      # 1e-6 rel on ratio
    tol_dict = float(CFG.pillar_surgery.ledger_scale_rel_tol)   # 5% band
    checks, data = [], {}

    eps = float(art("V78_lambda_protection", "data.eps_crit"))

    # ---- T0 gate: the booked leaf and the paper's loop exponent --------
    exp_loop = math.log(eps) / math.log(LOOP)
    t0a = abs(exp_loop - 70.42)
    t0b = abs(eps / 1.6e-120 - 1.0)
    data["t0"] = {"eps_crit": eps, "bare_loop_exponent": exp_loop,
                  "paper_contact_70_42_abs_diff": t0a,
                  "paper_rounded_1p6e120_rel_diff": t0b}
    checks.append({"name": "T0 the booked trace-phase leaf {:.17g}: bare "
                          "loop exponent {:.5f} reproduces the paper "
                          "contact 70.42 (|d| {:.3e} < 0.01); |eps/"
                          "1.6e-120 - 1| = {:.5f} < 0.02 (the paper "
                          "rounds the leaf to two significant "
                          "figures)".format(eps, exp_loop, t0a, t0b),
                   "pass": bool(t0a < 0.01 and t0b < 0.02)})
    if not (t0a < 0.01 and t0b < 0.02):
        data["status"] = "MACHINERY-DRIFT"
        print("V328 MACHINERY-DRIFT: the V78 trace-phase anchors do not "
              "re-derive; nothing is interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- class (i): findpoly at x = eps/LOOP^70, 50 dp -----------------
    x0 = eps / LOOP ** 70                # rescale: x0 is O(0.2)
    mpmath.mp.dps = 50
    x = mpmath.mpf(repr(x0))
    noise = mpmath.mpf("1e-14")
    found_poly = None
    per_degree = {}
    for d in range(1, dmax + 1):
        try:
            rel = mpmath.findpoly(x, d, maxcoeff=hcap)
        except ValueError:               # pslq nonzero-vector guard (V327)
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
            continue
        if res < noise:
            found_poly = {"degree": d, "height": height, "coeffs": coeffs,
                          "poly": poly_string(coeffs),
                          "residual_50dp": float(res)}
    data["findpoly"] = {"rescaled_at": "eps_crit / LOOP^70", "x": x0,
                        "degrees": "1..%d" % dmax, "height_cap": hcap,
                        "substitution_noise_floor": float(noise),
                        "per_degree": {str(d): per_degree[d]
                                       for d in per_degree}}
    mpmath.mp.dps = 15
    checks.append({"name": "class (i) findpoly of eps_crit/LOOP^70 = "
                           "{:.9g}, degrees 1..%d, height %d: %s"
                   % (dmax, hcap,
                      "no nontrivial relation" if found_poly is None
                      else "verified {}".format(found_poly)),
                   "pass": True})

    # ---- class (ii): rational census against eps/LOOP^k, k = 1..200 ----
    # Amended class (booked in review/pre_registrations.md): relative
    # band on the ratio, symmetric height bound max(p, q) <= 5000,
    # FOUND only under the 1e-12 relative-residual gate.
    qs = np.arange(1, DEN_MAX + 1, dtype=float)
    log_eps = math.log(eps)
    log_loop = math.log(LOOP)
    per_k = []
    hits = []
    n_unbounded_diag = 0
    for k in range(1, 201):
        # float LOOP**k underflows for k >~ 79: evaluate in log space
        # (relative evaluation noise ~1e-13 << the 1e-6 band).
        xk = math.exp(log_eps - k * log_loop)
        if xk > 1e150:
            continue
        ps = np.round(xk * qs)
        sel = ps >= 1.0
        if not bool(np.any(sel)):
            continue
        errs = np.abs(ps[sel] / qs[sel] - xk)
        n_unbounded_diag += int(np.count_nonzero(errs < band * xk))
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
        mpmath.mp.dps = 50
        best = min(hits, key=lambda h: h["err_rel"])
        eps_mp = mpmath.mpf(repr(eps))
        LOOPM = mpmath.mpf(1) / (16 * mpmath.pi)
        rule = (mpmath.mpf(best["p"]) / mpmath.mpf(best["q"])
                * LOOPM ** best["k"])
        resid = abs(eps_mp - rule)
        relv = float(resid / eps_mp)
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
        "amendment_unbounded_p_hit_count_mechanism": n_unbounded_diag,
        "ks": "1..200", "n_hits": len(hits),
        "hits_first5": hits[:5],
        "per_k_closest": per_k,
        "verified_best_contact": verified_best_contact,
        "verified_rule": found_rat}
    checks.append({"name": "class (ii, amended) rationals max(p,q) <= "
                           "{:d} inside relative {:.0e} of eps_crit/"
                           "LOOP^k over k = 1..200: {:d} contacts at k "
                           "with a representable bounded row ({} unbounded-"
                           "band hits are the booked mechanism); best "
                           "verified contact {} under the 1e-12 gate{}"
                   .format(DEN_MAX, band, len(hits), n_unbounded_diag,
                           verified_best_contact,
                           " (rule supplied)" if found_rat else ""),
                   "pass": True})

    # ---- class (iii): the V326 dictionary class x LOOP^k ---------------
    members = [float(s.split("=")[1]) for s in
               art("V326_phase_dictionary_bridge",
                   "data.dictionary_class.members")]
    best3 = None
    contacts3 = []
    for c in members:
        for k in range(1, 201):
            # log-space relative distance (both values ~1e-120-scale):
            delta = math.log(c) + k * log_loop - log_eps
            if best3 is None or abs(delta) < best3["abs_log_distance"]:
                best3 = {"c": c, "k": k, "abs_log_distance": abs(delta),
                         "value": c * math.exp(k * log_loop)
                         if k * log_loop > -700 else
                         math.exp(k * log_loop) * c}
            if abs(math.exp(delta) - 1.0) <= tol_dict and abs(delta) < 0.5:
                contacts3.append({"c": c, "k": k, "delta_log": delta,
                                  "rel_to_eps_crit":
                                      abs(math.exp(delta) - 1.0)})
    closest_rel3 = abs(math.exp(best3["abs_log_distance"]) - 1.0) \
        if best3["abs_log_distance"] < 50 else None
    data["dictionary_class"] = {"members": members, "ks": "1..200",
                                "band_rel": tol_dict,
                                "contacts": contacts3,
                                "closest": best3,
                                "closest_rel_to_eps_crit": closest_rel3,
                                "closest_band_multiples":
                                    closest_rel3 / tol_dict
                                    if closest_rel3 is not None else None}
    checks.append({"name": "class (iii) {} dictionary members x LOOP^k, "
                           "k = 1..200, {:.0%} band: {} contacts; closest "
                           "c = {} at k = {}: |dlog| = {:.4f} rel {:.3e} "
                           "= {:.1f} band-multiples".format(
                               len(members), tol_dict, len(contacts3),
                               best3["c"], best3["k"],
                               best3["abs_log_distance"], closest_rel3,
                               closest_rel3 / tol_dict),
                   "pass": True})

    # ---- branch ---------------------------------------------------------
    rule = None
    if found_poly is not None:
        rule = {"class": "i", "kind": "integer polynomial", **found_poly}
    elif found_rat is not None:
        rule = {"class": "ii", "kind": "rational x loop power",
                "rule": "eps_crit = ({:d}/{:d}) * LOOP^{:d}".format(
                    found_rat["p"], found_rat["q"], found_rat["k"]),
                **found_rat}
    elif contacts3:
        rule = {"class": "iii", "kind": "dictionary member x loop power",
                **contacts3[0]}
    branch = ("LAMBDA-CC-ALGEBRAIC-FOUND" if rule is not None
              else "LAMBDA-CC-NO-RELATION-DECLARED")
    data["status"] = branch
    data["relation"] = rule
    print("V328 {}: eps_crit = {:.17g} (bare loop exponent {:.4f}); "
          "findpoly of eps/LOOP^70 -> {}; rational contacts: {}; "
          "dictionary closest c = {} k = {} at {:.3e} rel.".format(
              branch, eps, exp_loop,
              found_poly["poly"] if found_poly else "none",
              len(hits), best3["c"], best3["k"], closest_rel3), flush=True)
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V328_algebraicity_lambda_cc", "tier": 3,
           "property": "is the trace-phase breaking eps_crit a number of "
                       "the booked algebraic/rational-loop/dictionary "
                       "closed class - could a rule supply the 1.6e-120 "
                       "(tier 3, O1)",
           "residual": True,
           "programme": "theory-hunt V328",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v328_algebraicity_lambda_cc.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V328 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
