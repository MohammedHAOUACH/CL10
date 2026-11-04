#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V316 - the weight double-contact (tier 2, hard, O8 + neutrino,
theory-hunt).

**The question (pre-registered).** Does a booked-class reweighting of the
traced-word sum move BOTH rung contacts simultaneously inside their bands -
rung-10: |coef(g=1/5) - COEF_OBS| <= pillar_surgery.coef_band_abs;
rung-17: |S17(variant)/3 - neutrino.m_atm_ev|/m_atm <=
pillar_surgery.ledger_scale_rel_tol - WITHOUT touching the mechanism or
the clock (g = 1/5 fixed)?

**The object (closed class, mechanism and clock fixed; pre-registered).**
W1 eta-law: V200's F10_of lifted verbatim (the eta-law functional
cnorm((M(eta)G)^n)/cnorm(G^n)) at the booked eta set
{angle.branches.closed = {pi/4}, V202's data.selection.nearest_crossing}.
W2 phase-class renormalization: the traced words reweighted by V202's
phase-class index k (the even sum of half-eigenvalues of the word) with
multipliers {1, cos(k g), cos(2 k g)} at the clock g; the class matrix
vectors P_k^(n) are built exactly by V202's prefix construction as a
Laurent polynomial in the class index (identical objects, n = 1..20); the
m = 1 member is gated against loop_ratios (< 1e-9 - V202's "the phase
classes ARE the census" equality).
W3 node-set channel factors: V205's block-3 k=0 propagator-trace filter
lifted verbatim (w = SC.phases(SC.bloch(q, U16(g), g10 = g)),
f_s = |sum_j exp(i Re w_j)|/128) at each of the 12 SC.M_POINTS nodes at
g = 1/5; all C(12,3) = 220 triples; the coefficient is read with the
triple's combined factor per the booked required-sum-factor relation
(V203/V204/V205): coef(T) = (sum_{s in T} f_s/f_bar) x (F_10/LOOP)^10,
f_bar = the 12-factor mean at the clock - the normalization that
reproduces the booked sum-factor-3 plain reading when factors are equal
(stated pre-run; the node set is combinatorial, never fitted).
Rung-17 leg per variant (pre-registered): S17^var = nv.loop_ladder(17) x
nv.EV_PER_GEV x |F_17^var / F_17^plain| - the plain member has ratio
exactly 1 and reproduces the booked V314 S17/3 leaf; W3 leaves the F
machinery untouched (ratio 1). Ladders run nmax = 20 (F_n for n <= 12 is
bit-identical to the nmax=12 leg: loop_ratios is incremental).

**Measured (pre-registered).** T0: plain rung-10 read reproduces V206's
coef_at_fifth to < 1e-9 AND the plain rung-17 leg reproduces V314's
booked S17/3 to < 1e-9. The W3 leg runs exactly 220 triples (counted in
the artifact). Every variant's two normalized residuals and their max.

**Acceptance (pre-registered).** WEIGHT-DOUBLE-CONTACT iff some variant
passes both bands (its full reweight rule booked); WEIGHT-NO-DOUBLE
otherwise, with the full table ranked by the max of the two normalized
residuals, top 20 rows and the per-family best booked.
"""
import itertools
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                               # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))
sys.path.insert(0, exp("V85_closed_loop_factor"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from V85_closed_loop_factor import cartan_norm            # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
G5 = float(CFG.angle.branches.implemented)    # the rung-10 clock 1/5
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
TOL_SCALE = float(CFG.pillar_surgery.ledger_scale_rel_tol)
M_ATM = float(CFG.neutrino.m_atm_ev)
NMAX = 20                          # rung-17 leg needs F_17 (pre-reg)


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def read_coef(Fv, LOOP):
    """The published channel reading 3 (16 pi F)^10 = 3 (F/LOOP)^10
    (V308's read_coef, copied)."""
    return 3.0 * (Fv / LOOP) ** 10


def eta_ladder(G, Ma, ev, cartan, Gnorm, eta, n):
    """V200's F10_of eta-law lifted to length n:
    cnorm((M(eta) G)^n)/cnorm(G^n)."""
    M = sum(np.exp(1j * (eta / 2.0) * e) * m
            for e, m in zip(ev, Ma))
    MG = M @ G
    return cartan_norm(np.linalg.matrix_power(MG, n), cartan) / Gnorm[n]


def class_polynomial(G, Ma, half_ev, nmax):
    """V202's phase-class decomposition as a matrix Laurent polynomial:
    {n: {k: P_k^(n)}}, (M(g) G)^n = sum_k e^{i k g} P_k^(n). V202's
    prefix expansion W_new = W @ B0[lam] carried at the matrix level
    (its C_k are exactly the Cartan vectors of these P_k)."""
    B0 = [m @ G for m in Ma]
    P = {int(h): B.copy() for h, B in zip(half_ev, B0)}
    out = {}
    for n in range(1, nmax + 1):
        out[n] = {k: v.copy() for k, v in P.items()}
        if n == nmax:
            break
        Q = {}
        for k, Wm in P.items():
            for h, B in zip(half_ev, B0):
                kk = k + int(h)
                Q[kk] = Q.get(kk, np.zeros_like(Wm)) + Wm @ B
        P = Q
    return out


def main():
    t0 = time.time()
    checks, data = [], {}
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)

    # ---- the fixed mechanism at the fixed clock ------------------------
    G, cartan, Ma, ev = hop_machinery(G5)
    ev = [float(e) for e in ev]
    Gnorm = {n: cartan_norm(np.linalg.matrix_power(G, n), cartan)
             for n in range(1, NMAX + 1)}
    F, _rho = loop_ratios(G, cartan, Ma, ev, G5, NMAX)

    # ---- T0a: plain rung-10 read reproduces V206 coef_at_fifth --------
    coef5_art = float(art("V206_clock_tension",
                          "data.data_clock.coef_at_fifth"))
    coef_plain = read_coef(F[10], LOOP)
    t0a = abs(coef_plain - coef5_art)
    data["t0a_coef"] = {"read": coef_plain, "artifact": coef5_art,
                        "abs_diff": t0a}
    checks.append({"name": "T0a plain rung-10 read at g = 1/5 reproduces "
                          "V206's booked coef_at_fifth {:.15f} (|d| "
                          "{:.1e} < 1e-9)".format(coef5_art, t0a),
                   "pass": bool(t0a < 1e-9)})

    # ---- T0b: plain rung-17 leg reproduces V314's booked S17/3 --------
    S17 = float(nv.loop_ladder(17)) * float(nv.EV_PER_GEV)
    s17_art = float(art("V314_four_scale_ledger",
                        "data.rederived.S17_over_3_eV"))
    t0b = abs(S17 / 3.0 - s17_art)
    data["t0b_S17"] = {"read_eV": S17 / 3.0, "artifact": s17_art,
                       "abs_diff": t0b}
    checks.append({"name": "T0b plain rung-17 leg reproduces V314's booked "
                          "S17/3 {:.9f} eV (|d| {:.1e} < 1e-9)".format(
                              s17_art, t0b),
                   "pass": bool(t0b < 1e-9)})
    if not (t0a < 1e-9 and t0b < 1e-9):
        data["status"] = "MACHINERY-DRIFT"
        print("V316 MACHINERY-DRIFT: the rung-10/rung-17 machinery no "
              "longer reproduces the booked leaves; nothing interpreted",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- W2 machinery: the class polynomial + its census gate ---------
    # the class index must be the exact integer half-eigenvalue (V202's
    # own guard: even eigenvalues to 1e-9; int() truncation of a noisy
    # float would mis-merge classes - the bug the first run caught)
    ev_even = all(abs(e - round(e)) < 1e-9 and int(round(e)) % 2 == 0
                  for e in ev)
    if not ev_even:
        data["status"] = "MACHINERY-DRIFT"
        print("V316 MACHINERY-DRIFT: the link Casimir eigenvalues are "
              "not the booked even integers -32/-12/+4", flush=True)
        _emit(data, checks, ok=False)
        return 1
    half_ev = [int(round(e / 2.0)) for e in ev]
    Pk = class_polynomial(G, Ma, half_ev, NMAX)
    recon_err = 0.0
    for n in range(1, NMAX + 1):
        Bn = sum(np.exp(1j * k * G5) * Mm for k, Mm in Pk[n].items())
        Fn_rec = cartan_norm(Bn, cartan) / Gnorm[n]
        recon_err = max(recon_err, abs(Fn_rec - abs(F[n]))
                        / max(abs(F[n]), 1e-30))
    data["w2_census_gate"] = {"max_rel_recon_err_vs_loop_ratios":
                              recon_err, "n_range": [1, NMAX]}
    checks.append({"name": "T0c (W2 class gate) the m = 1 phase-class "
                          "polynomial reproduces loop_ratios' F_n at the "
                          "clock for every n = 1..{} to {:.1e} < 1e-9 "
                          "(V202's phase-classes-are-the-census equality "
                          "carried at the matrix level)".format(
                              NMAX, recon_err),
                   "pass": bool(recon_err < 1e-9)})
    if recon_err >= 1e-9:
        data["status"] = "MACHINERY-DRIFT"
        print("V316 MACHINERY-DRIFT: the phase-class ladder does not "
              "reconstruct the booked loop ladder", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- W3 machinery: the 12 M-node channel factors at the clock ----
    _B, _U32, U16c = SC.sectors(G5)
    f_nodes = {}
    for i, q in enumerate(SC.M_POINTS):
        w = SC.phases(SC.bloch(np.array(q, float), U16c, g10=G5))
        f_nodes[str(i)] = {"q": [float(v) for v in q],
                           "factor": float(abs(
                               sum(np.exp(1j * float(np.real(x)))
                                   for x in w)) / 128.0)}
    fvals = [v["factor"] for v in f_nodes.values()]
    fbar = sum(fvals) / len(fvals)
    data["w3_node_factors"] = {"clock": G5, "f_bar": fbar,
                               "factors": f_nodes,
                               "spread_max_over_fbar":
                                   max(abs(v - fbar) for v in fvals)
                                   / fbar}
    triples = list(itertools.combinations(range(12), 3))
    checks.append({"name": "W3 node-set census runs exactly {} triples "
                          "over the 12 M nodes at g = 1/5; node factor "
                          "spread about the mean max {:.3e} (if all "
                          "factors are equal the 220 readings coincide "
                          "with the plain sum-factor-3 reading - the "
                          "booked required-sum-factor composition)".format(
                              len(triples),
                              max(abs(v - fbar) for v in fvals) / fbar),
                   "pass": bool(len(triples) == 220)})

    # ---- the variant table --------------------------------------------
    variants = []

    def rung_legs(Fn10, Fn17):
        coef = read_coef(Fn10, LOOP)
        # the plain F_17 at g = 1/5 is not exactly zero (the exact odd-n
        # cancellation is booked only at the order-8 clock class, V206),
        # but the ratio is guarded so a vanishing denominator degrades
        # to an out-of-band reading, never a crash
        s17v = (S17 * (abs(Fn17) / abs(F[17]))
                if abs(F[17]) > 1e-18 else float("inf"))
        n10 = abs(coef - COEF_OBS) / BAND_C
        n17 = abs(s17v / 3.0 - M_ATM) / M_ATM / TOL_SCALE
        return coef, s17v, n10, n17

    # plain (the reference row)
    c0, s0, n010, n017 = rung_legs(F[10], F[17])
    variants.append({"family": "plain", "rule":
                     "plain rung-10 read 3(F_10/LOOP)^10 at g = 1/5",
                     "prov": "algebraic machinery, config clock",
                     "coef": c0, "S17_over_3_eV": s0 / 3.0,
                     "norm_rung10": n010, "norm_rung17": n017,
                     "worst_normalized": max(n010, n017)})
    # W1 eta-law at the booked eta set
    eta_art = float(art("V202_loop_angle_selection",
                        "data.selection.nearest_crossing"))
    eta_set = {"pi/4 (angle.branches.closed)":
               float(CFG.angle.branches.closed),
               "V202 eta* (data.selection.nearest_crossing)": eta_art}
    for label, eta in eta_set.items():
        F10e = eta_ladder(G, Ma, ev, cartan, Gnorm, eta, 10)
        F17e = eta_ladder(G, Ma, ev, cartan, Gnorm, eta, 17)
        c, s, n10, n17 = rung_legs(F10e, F17e)
        variants.append({"family": "W1-eta", "rule":
                         "V200 eta-law F_of(eta) at eta = {:.9f} [{}]"
                         .format(eta, label),
                         "eta": eta, "prov":
                         ("config-declared (closed branch angle)"
                          if "closed" in label else "measured (V202 "
                          "crossing)"),
                         "F10": float(F10e), "coef": c,
                         "S17_over_3_eV": s / 3.0,
                         "norm_rung10": n10, "norm_rung17": n17,
                         "worst_normalized": max(n10, n17)})
    # W2 phase-class multipliers {cos(k g), cos(2 k g)} (m = 1 is the
    # plain member, gated as T0c and not duplicated as a variant)
    for mname, mfun in (("cos(k*g)", lambda k: math.cos(k * G5)),
                        ("cos(2*k*g)", lambda k: math.cos(2 * k * G5))):
        F10w = (cartan_norm(sum(mfun(k) * np.exp(1j * k * G5) * Mm
                                for k, Mm in Pk[10].items()),
                            cartan) / Gnorm[10])
        F17w = (cartan_norm(sum(mfun(k) * np.exp(1j * k * G5) * Mm
                                for k, Mm in Pk[17].items()),
                            cartan) / Gnorm[17])
        c, s, n10, n17 = rung_legs(F10w, F17w)
        variants.append({"family": "W2-phaseclass",
                         "rule": "class-k multiplier {} on the traced "
                                 "words at clock g = 1/5".format(mname),
                         "prov": "combinatorial (V202 phase classes)",
                         "F10": float(F10w), "coef": c,
                         "S17_over_3_eV": s / 3.0,
                         "norm_rung10": n10, "norm_rung17": n17,
                         "worst_normalized": max(n10, n17)})
    # W3 all 220 node triples
    base_c = (F[10] / LOOP) ** 10
    for T in triples:
        comb = sum(f_nodes[str(s)]["factor"] for s in T) / fbar
        coef = comb * base_c
        s17v = S17
        n10 = abs(coef - COEF_OBS) / BAND_C
        n17 = abs(s17v / 3.0 - M_ATM) / M_ATM / TOL_SCALE
        variants.append({"family": "W3-triples",
                         "rule": "node triple {} (M_POINTS order), "
                                 "combined sum factor {:.6f}".format(
                                     list(T), comb),
                         "triple": list(T), "sum_factor": float(comb),
                         "prov": "combinatorial (node set)",
                         "coef": coef, "S17_over_3_eV": s17v / 3.0,
                         "norm_rung10": n10, "norm_rung17": n17,
                         "worst_normalized": max(n10, n17)})

    variants.sort(key=lambda v: v["worst_normalized"])
    double = [v for v in variants
              if v["norm_rung10"] <= 1.0 and v["norm_rung17"] <= 1.0]
    data["coefficient_observed"] = COEF_OBS
    data["band_coef"] = BAND_C
    data["band_scale_rel"] = TOL_SCALE
    data["m_atm_eV"] = M_ATM
    data["variants"] = variants
    data["n_variants"] = len(variants)
    fam_best = {}
    for v in variants:
        if v["family"] not in fam_best:
            fam_best[v["family"]] = v
    data["per_family_best"] = fam_best
    data["top20"] = variants[:20]

    branch = "WEIGHT-DOUBLE-CONTACT" if double else "WEIGHT-NO-DOUBLE"
    data["status"] = branch
    if double:
        data["double_contact_rule"] = double[0]["rule"]
    checks.append({"name": "branch {}: {} variants; double band (|d| <= "
                          "{:.0e} on coef and rel <= {:.0e} on S17/3) "
                          "hit by {}; best row [{}] worst-normalized "
                          "{:.4f} (coef {:.7f}, S17/3 {:.9f} eV); "
                          "per-family bests: {}".format(
                              branch, len(variants), BAND_C, TOL_SCALE,
                              len(double), "none" if not variants
                              else variants[0]["family"],
                              variants[0]["worst_normalized"],
                              variants[0]["coef"],
                              variants[0]["S17_over_3_eV"],
                              {k: "{:.4f}".format(v["worst_normalized"])
                               for k, v in fam_best.items()}),
                   "pass": True})
    print("V316 {}: {} variants (plain 1, W1 2, W2 2, W3 220); best "
          "worst-normalized {:.4f} [{}]; double-band hits {}.".format(
              branch, len(variants), variants[0]["worst_normalized"],
              variants[0]["family"], len(double)), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V316_weight_double_contact", "tier": 2,
           "property": "booked reweightings of the traced-word sum "
                       "against BOTH rung contacts (rung-10 coefficient "
                       "and rung-17 neutrino) at the fixed clock and "
                       "mechanism (tier 2, hard, O8+neutrino)",
           "residual": True, "programme": "theory-hunt V316",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v316_weight_double_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V316 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
