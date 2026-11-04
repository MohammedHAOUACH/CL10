#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V318 - the closest-theory ledger, five contacts (tier 1, hard,
theory-hunt; the arbiter and the deliverable).

**The question (pre-registered).** Among ALL compound theories = (clock
from the booked set) x (weight variant) x (selection), ranked on FIVE
contacts, which comes closest to the measured values, and does any
zero-fitted row pass them all?

**Rows (pre-registered).** Clocks: 1/5 (angle.branches.implemented,
fitted), the V309+V310 corrected clock (V314's row-(d) construction:
g_cross_rung10 minus the booked V310 supplied-hit value), every V315
closure clock carrying its odd-rule status, and the V311 p=0 closure
clocks pi/2, pi/4, 3pi/4 (deduplicated against the V315 set). Weights:
plain, the V316 WEIGHT-DOUBLE-CONTACT winners if any else the top-3 by
the double ranking, and the W1 eta values. Selection: none, the V317 S3
unique triple if a cell fires, the V317 S1 member if a cell fires.
Cross-product cap 200; if over, keep top-3 per dimension and book the
pruning rule (mandatory clocks {1/5, corrected, pi/2, pi/4, 3pi/4}; the
remaining V315 closure clocks enter only if the product fits, else only
the top-3 by plain coefficient contact - booked in the check string).
EVERY weight functional is read at the ROW's clock (V314's per-row clock
convention); rung-17 legs use the variant F-ratio at the row clock as
V316 books; the W2 multipliers are class functions of the clock,
evaluated at the row clock.

**Contacts (V314's scoring block copied verbatim, normalized-residual
convention unchanged).** C1 seesaw M_P*LOOP^3 (algebraic, normalized 0);
C2 DM window physics.wimp_window_gev (band ledger_scale_rel_tol); C3
rung-17 neutrino S17^var/3 vs neutrino.m_atm_ev (same band); C4 rung-10
coefficient vs COEF_OBS (coef_band_abs); C5 the new desert-coupling
contact: V247's booked criterion - v237_rule at alpha_std(row clock),
residual |dz12_efolds - V237-booked|, band 1e-9 (V247's T0 desert
reproduction rule, host float64), the row's shift_over_spacing_12 booked
against the booked 1.1171513611688386.

**Provenance (pre-registered).** every clock/weight/selection tagged
{algebraic, config, measured, combinatorial, fitted}; zero-fitted = no
fitted input anywhere in the row.

**T0 (pre-registered).** V314's 11 rows re-run inside V318 reproduce the
booked best_zero_fitted_row.worst_normalized leaf to < 1e-9 (the current
gated artifact leaf wins over the plan's superseded literal - pinned in
the pre-registration).

**Acceptance (pre-registered).** LEDGER-CONSISTENT-THEORY-CANDIDATE iff
some zero-fitted row passes all five contacts; LEDGER-CLOSEST-THEORY-
DECLARED otherwise: the full ranked table booked and THE closest row
named explicitly with its five normalized residuals, its provenance and
its distance (worst normalized). Closest is not true; the booking keeps
that distinction.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))
sys.path.insert(0, exp("V85_closed_loop_factor"))
sys.path.insert(0, exp("V247_two_dial_desert_contacts"))

import numpy as np                                        # noqa: E402
from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from V85_closed_loop_factor import cartan_norm            # noqa: E402
from V247_two_dial_desert_contacts import v237_rule       # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import coupling as _cpl                         # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
NMAX = 260                       # V314's plain-read range (verbatim)
NMAX_W = 20                      # V316's weight-ladder range
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
TOL_SCALE = float(CFG.pillar_surgery.ledger_scale_rel_tol)
BAND_C5 = 1e-9                   # V247's booked desert reproduction band


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def artj(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def plain_coef_at(g):
    """Re-derived rung-10 / crossing read of the standard ladder (V314
    verbatim)."""
    G, cartan, Ma, ev = hop_machinery(g)
    F, _rho = loop_ratios(G, cartan, Ma, ev, g, NMAX)
    n_star = next((n for n in range(2, NMAX + 1, 2)
                   if abs(F[n]) < LOOP), None)
    rung10 = 3.0 * (abs(F[10]) / LOOP) ** 10
    cross = (None if n_star is None
             else 3.0 * (abs(F[n_star]) / LOOP) ** 10)
    return {"rung10": rung10, "n_star": n_star, "crossing": cross}


LOOP = float(u10pred.LOOP)
MP = float(u10pred.M_PLANCK)
COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)


def weight_read(g, wdesc):
    """The weight functional at the ROW clock g: returns (coef,
    rung-17 F-ratio). 'plain' = V314's per-row rung-10 read with ratio 1;
    W1 = V200 eta-law at the eta parameter (V316); W2 = class-polynomial
    multiplier at the clock (V316, multipliers are class functions of the
    row clock); W3 = V205 node-factor sum composition at the clock with
    the triple (F machinery untouched: ratio 1)."""
    fam = wdesc["family"]
    if fam == "plain":
        return plain_coef_at(g)["rung10"], 1.0
    G, cartan, Ma, ev = hop_machinery(g)
    Gnorm = {n: cartan_norm(np.linalg.matrix_power(G, n), cartan)
             for n in (10, 17)}
    ev = [float(e) for e in ev]
    Fp, _ = loop_ratios(G, cartan, Ma, ev, g, NMAX_W)
    if fam == "W1-eta":
        eta = wdesc["eta"]

        def Fn_eta(n):
            M = sum(np.exp(1j * (eta / 2.0) * e) * m
                    for e, m in zip(ev, Ma))
            return (cartan_norm(np.linalg.matrix_power(M @ G, n),
                                cartan) / Gnorm[n])
        F10w, F17w = Fn_eta(10), Fn_eta(17)
        ratio = (abs(F17w) / abs(Fp[17]) if abs(Fp[17]) > 1e-18
                 else float("inf"))
        return 3.0 * (F10w / LOOP) ** 10, ratio
    if fam == "W2-phaseclass":
        # class polynomial at the row clock (V316's exact construction;
        # exact integer half-eigenvalues per V202's even-eigenvalue guard)
        half_ev = [int(round(e / 2.0)) for e in ev]
        B0 = [m @ G for m in Ma]
        P = {int(h): B.copy() for h, B in zip(half_ev, B0)}
        Pk = {}
        for n in range(1, NMAX_W + 1):
            Pk[n] = {k: v.copy() for k, v in P.items()}
            if n == NMAX_W:
                break
            Q = {}
            for k, Wm in P.items():
                for h, B in zip(half_ev, B0):
                    kk = k + int(h)
                    Q[kk] = Q.get(kk, np.zeros_like(Wm)) + Wm @ B
            P = Q
        mfun = (lambda k: math.cos(k * g)) if wdesc["m"] == "cos(k*g)" \
            else (lambda k: math.cos(2 * k * g))

        def Fw(n):
            Bn = sum(mfun(k) * np.exp(1j * k * g) * Mm
                     for k, Mm in Pk[n].items())
            return cartan_norm(Bn, cartan) / Gnorm[n]
        F10w, F17w = Fw(10), Fw(17)
        ratio = (abs(F17w) / abs(Fp[17]) if abs(Fp[17]) > 1e-18
                 else float("inf"))
        return 3.0 * (F10w / LOOP) ** 10, ratio
    # W3-triples: node factors of V205 block 3 at the row clock
    _B, _U32, U16c = SC.sectors(g)
    fvals = []
    for q in SC.M_POINTS:
        w = SC.phases(SC.bloch(np.array(q, float), U16c, g10=g))
        fvals.append(float(abs(sum(np.exp(1j * float(np.real(x)))
                                   for x in w)) / 128.0))
    fbar = sum(fvals) / len(fvals)
    comb = sum(fvals[i] for i in wdesc["triple"]) / fbar
    base_c = (abs(Fp[10]) / LOOP) ** 10
    return comb * base_c, 1.0


def main():
    t0 = time.time()
    checks, data = [], {}

    # ---- contacts 1-3 re-derived (V314 verbatim) -----------------------
    rows197 = art("V197_rung_arithmetic", "data.ladder_rows")
    book = {int(r["rung"]): float(r["value_gev"]) for r in rows197}
    scale_derived = {n: MP * LOOP ** n for n in (3, 10, 17)}
    gate = {n: abs(scale_derived[n] - book[n]) / abs(book[n])
            for n in (3, 10, 17)}
    S17 = float(nv.loop_ladder(17)) * float(nv.EV_PER_GEV)
    m_atm = float(CFG.neutrino.m_atm_ev)
    wimp_lo, wimp_hi = [float(x) for x in CFG.physics.wimp_window_gev]
    gate["S17"] = abs(S17 - float(art("V303_junction_energy_r",
                                      "data.contact.S17_eV"))) \
        / float(art("V303_junction_energy_r", "data.contact.S17_eV"))
    t0ok = all(v < TOL_SCALE for v in gate.values())
    data["scale_gates"] = {str(k): v for k, v in gate.items()}
    checks.append({"name": "T0a re-derivation gates (rel < {:.2f}): {}"
                   .format(TOL_SCALE, {str(k): "{:.1e}".format(v)
                                       for k, v in gate.items()}),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V318 MACHINERY-DRIFT: ladder scales do not re-derive",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    c1 = {"residual": 0.0, "band": None, "normalized": 0.0,
          "detail": "seesaw M_P*LOOP^3 = {:.6e} GeV (algebraic, rung A "
                    "= {})".format(scale_derived[3],
                                   int(CFG.cell.rung_A))}
    inside = wimp_lo <= scale_derived[10] <= wimp_hi
    edge = (0.0 if inside else
            min(abs(scale_derived[10] - wimp_lo),
                abs(scale_derived[10] - wimp_hi)) / wimp_lo)
    c2 = {"residual": edge, "band": TOL_SCALE,
          "normalized": edge / TOL_SCALE,
          "detail": "DM M_P*LOOP^10 = {:.4f} GeV vs window [{:.0f}, "
                    "{:.0f}]: {}".format(scale_derived[10], wimp_lo,
                                         wimp_hi,
                                         "IN WINDOW" if inside else
                                         "outside")}

    # ---- T0b: V314's 11 rows re-run inside (verbatim construction) ----
    g5 = float(CFG.angle.branches.implemented)
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    hits310 = art("V310_coefficient_expression_census", "data.supplied")
    unitary_p = [int(p) for p in art("V311_quaternionic_root_step",
                                     "data.unitary_p")]
    sig311 = art("V311_quaternionic_root_step", "data.signatures")
    v312 = artj("V312_quaternionic_loop_ladder")
    ladder312 = v312["data"].get("ladder_runs", {})
    scan = [float(g) for g in CFG.pillar_surgery.scan_clocks]

    v314_rows = []
    rd = plain_coef_at(g5)
    v314_rows.append({"label": "1/5", "g": g5, "coef": rd["rung10"],
                      "fitted": True})
    rd_o = plain_coef_at(g_obs)
    v314_rows.append({"label": "g_obs", "g": g_obs, "coef": rd_o["rung10"],
                      "fitted": True})
    for p in unitary_p:
        for key, r in ladder312.get(str(p), {}).items():
            if r["coef"] is None:
                continue
            v314_rows.append({"label": "p={:d} closure {}".format(p, key),
                              "g": r["clock"], "coef": r["coef"],
                              "fitted": False})
    gap_hits = [h for h in hits310 if "clock_gap" in h["target"]]
    if gap_hits:
        h = gap_hits[0]
        g_cand = g_cross - h["value"]
        rd_c = plain_coef_at(g_cand)
        v314_rows.append({"label": "V310 corrected", "g": g_cand,
                          "coef": rd_c["rung10"], "fitted": False})
    for g in scan:
        rd_s = plain_coef_at(g)
        if rd_s["crossing"] is None:
            continue
        v314_rows.append({"label": "scan {}".format(g), "g": g,
                          "coef": rd_s["crossing"], "fitted": False})
    nu_rel = abs(S17 / 3.0 - m_atm) / m_atm
    scored314 = []
    for r in v314_rows:
        c4r = abs(r["coef"] - COEF_OBS)
        worst = max(0.0, edge / TOL_SCALE, nu_rel / TOL_SCALE,
                    c4r / BAND_C)
        scored314.append({"row": r, "worst": worst,
                          "zero_fitted": not r["fitted"]})
    scored314.sort(key=lambda s: (s["worst"], s["row"]["label"]))
    best314 = next((s for s in scored314 if s["zero_fitted"]), None)
    leaf314 = float(art("V314_four_scale_ledger",
                        "data.best_zero_fitted_row.worst_normalized"))
    t0b = abs(best314["worst"] - leaf314) if best314 else float("inf")
    data["t0b_v314_rows"] = {"n_rows": len(v314_rows),
                             "best_zf_worst_rederived":
                             (best314["worst"] if best314 else None),
                             "artifact_leaf": leaf314,
                             "abs_diff": t0b}
    checks.append({"name": "T0b V314's {} rows re-run inside V318 "
                          "reproduce the booked best-zero-fitted "
                          "worst-normalized {:.16f} (|d| {:.1e} < "
                          "1e-9)".format(len(v314_rows), leaf314, t0b),
                   "pass": bool(t0b < 1e-9)})
    if t0b >= 1e-9:
        data["status"] = "MACHINERY-DRIFT"
        print("V318 MACHINERY-DRIFT: the V314 rows do not reproduce the "
              "booked leaf; nothing extended", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- C5 booked criterion (V237/V247) ------------------------------
    a_dz12 = float(art("V237_desert_threshold_matching",
                       "data.matching.dz12_efolds"))
    a_shift = float(art("V237_desert_threshold_matching",
                        "data.matching.shift_over_spacing_12"))

    # ---- the row dimensions -------------------------------------------
    v315 = artj("V315_closure_quantizer_census")
    v316 = artj("V316_weight_double_contact")
    v317 = artj("V317_selection_conjunction")

    clock_dims = []   # (label, g, prov, extra)
    clock_dims.append({"label": "1/5 (angle.branches.implemented)",
                       "g": g5, "prov": "fitted (1/n_qubits reference "
                       "spectrum)", "odd_exact": None})
    g_cand = g_cross - gap_hits[0]["value"]
    clock_dims.append({"label": "g_cross_rung10 - ({}) the V310 "
                              "clock-gap hit".format(
                                  gap_hits[0]["expression"]),
                       "g": g_cand,
                       "prov": "algebraic/conventional/measured (V310 "
                               "CLASS-PRE census hit)", "odd_exact": None})
    seen_g = {round(g5, 12), round(g_cand, 12)}
    closure_clocks = {}
    for v in v315["data"]["variants"]:
        for k in v["closures"]:
            closure_clocks[k] = v315["data"]["clock_leg"][k]
    for k, leg in sorted(closure_clocks.items(),
                         key=lambda kv: kv[1]["residual"]):
        if round(leg["g"], 12) in seen_g:
            continue
        seen_g.add(round(leg["g"], 12))
        clock_dims.append({
            "label": "V315 closure {} (odd-rule {})".format(
                k, "exact" if leg["odd_rule_exact"] else "broken"),
            "g": leg["g"],
            "prov": "algebraic/combinatorial (V315 closure census)",
            "odd_exact": bool(leg["odd_rule_exact"]),
            "plain_residual": leg["residual"]})
    # V311 p=0 clocks are in the V315 set by construction (identity is a
    # class member); any missing ones would enter here:
    for key, row in sorted(sig311.get(str(unitary_p[0]), {})
                           .get("closure_clocks", {}).items()):
        gk = math.pi * int(key.split("*")[1].split("/")[0]) \
            / int(key.split("/")[1])
        if round(gk, 12) in seen_g:
            continue
        seen_g.add(round(gk, 12))
        clock_dims.append({"label": "V311 p=0 closure {} (k={:d})".format(
                              key, row["k"]), "g": gk,
                           "prov": "algebraic (V311 closure census)",
                           "odd_exact": None})

    weight_dims = [{"family": "plain", "rule": "plain rung-10 read at "
                    "the row clock (V314 convention)",
                    "prov": "algebraic machinery",
                    "norm316": None}]
    if v316["branch"] == "WEIGHT-DOUBLE-CONTACT":
        wsel = [v for v in v316["data"]["variants"]
                if v["norm_rung10"] <= 1.0 and v["norm_rung17"] <= 1.0]
    else:
        wsel = v316["data"]["variants"][1:4]      # top-3 by the double
        # ranking (the plain row is already the 'plain' dimension member)
    for v in wsel:
        wd = {"family": v["family"], "rule": v["rule"],
              "prov": v["prov"], "norm316": v["worst_normalized"]}
        if v["family"] == "W1-eta":
            wd["eta"] = float(v["eta"])
        elif v["family"] == "W2-phaseclass":
            wd["m"] = ("cos(k*g)" if "cos(k*g)" in v["rule"]
                       else "cos(2*k*g)")
        elif v["family"] == "W3-triples":
            wd["triple"] = [int(i) for i in v["triple"]]
        weight_dims.append(wd)
    for v in v316["data"]["variants"]:
        if v["family"] == "W1-eta" and not any(
                w.get("eta") == float(v["eta"]) for w in weight_dims):
            weight_dims.append({"family": "W1-eta",
                                "rule": v["rule"] + " (W1 eta set)",
                                "prov": v["prov"], "eta": float(v["eta"]),
                                "norm316": v["worst_normalized"]})

    sel_dims = [{"label": "none", "prov": "trivial"}]
    s3 = v317["data"].get("F3S3", {})
    s3f = v317["data"].get("F2S3", {})
    s3p = v317["data"].get("F4S3", {})
    if (s3.get("unique") or s3f.get("unique") or s3p.get("unique")):
        sel_dims.append({"label": "V317 S3 unique triple {}".format(
                             s3.get("best_triple") or
                             s3f.get("best_triple") or
                             s3p.get("best_triple")),
                         "prov": "combinatorial (V317 selected node "
                                 "set)"})
    if v317["data"].get("F1S1", {}).get("unique"):
        sel_dims.append({"label": "V317 S1 member (T* minimizer)",
                         "prov": "measured (V313/V317 selection)"})

    n_prod = len(clock_dims) * len(weight_dims) * len(sel_dims)
    pruned = False
    if n_prod > 200:
        pruned = True
        # pre-registered pruning, booked verbatim: when the cross product
        # exceeds 200, mandatory clocks stay ({1/5, corrected, and the
        # V311 p=0 clocks pi/2, pi/4, 3pi/4 named by the row list}) and
        # the remaining V315 closure clocks enter only as the top-3 by
        # plain rung-10 coefficient contact at the clock (the V315
        # clock_leg residual); the weight dimension keeps its booked
        # composition (plain + top-3 + W1 eta values); selections all.
        mand = [g5, g_cand, math.pi / 2, math.pi / 4, 3 * math.pi / 4]
        keep = [c for c in clock_dims
                if any(abs(c["g"] - m) < 1e-12 for m in mand)]
        rest = sorted((c for c in clock_dims
                       if not any(abs(c["g"] - m) < 1e-12
                                  for m in mand)),
                      key=lambda c: c.get("plain_residual", 1e30))[:3]
        clock_dims = keep + rest
        n_prod = (len(clock_dims) * len(weight_dims) * len(sel_dims))

    # ---- score every compound row --------------------------------------
    rows = []
    for cd in clock_dims:
        for wd in weight_dims:
            coef, ratio = weight_read(cd["g"], wd)
            s17v = S17 * ratio
            nu_rel_v = abs(s17v / 3.0 - m_atm) / m_atm
            for sd in sel_dims:
                c3 = {"residual": nu_rel_v, "band": TOL_SCALE,
                      "normalized": nu_rel_v / TOL_SCALE,
                      "detail": "S17/3 = {:.9f} eV vs m_atm {:.4f} eV: "
                                "{:+.4%} (rung-17 F-ratio {})".format(
                                    s17v / 3.0, m_atm,
                                    s17v / 3.0 / m_atm - 1.0,
                                    "{:.6g}".format(ratio))}
                c4r = abs(coef - COEF_OBS)
                c4 = {"residual": c4r, "band": BAND_C,
                      "normalized": c4r / BAND_C,
                      "detail": "coef {:.7f} vs observed {:.7f} (|d| "
                                "{:.2e}, band {:.0e})".format(
                                    coef, COEF_OBS, c4r, BAND_C)}
                fires5, row5 = v237_rule(_cpl.alpha_std(cd["g"]))
                c5r = abs(row5["dz12_efolds"] - a_dz12)
                c5 = {"residual": c5r, "band": BAND_C5,
                      "normalized": c5r / BAND_C5,
                      "detail": "desert (1,2)-shift {:.15g} e-folds vs "
                                "booked {:.15g} (|d| {:.2e}, band "
                                "{:.0e}); shift/spacing {:.9f} vs booked "
                                "{:.9f}; rule fires: {}".format(
                                    row5["dz12_efolds"], a_dz12, c5r,
                                    BAND_C5, row5["shift_over_spacing_12"],
                                    a_shift, fires5)}
                contacts = [c1, c2, c3, c4, c5]
                worst = max(c["normalized"] for c in contacts)
                prov = [cd["prov"], wd["prov"], sd["prov"]]
                fitted = any("fitted" in p for p in prov)
                rows.append({
                    "clock": cd["label"], "clock_g": cd["g"],
                    "weight": wd["rule"], "selection": sd["label"],
                    "coef": coef, "S17_over_3_eV": s17v / 3.0,
                    "shift_over_spacing_12":
                        row5["shift_over_spacing_12"],
                    "provenance": prov,
                    "normalized": [c["normalized"] for c in contacts],
                    "contacts": [c["detail"] for c in contacts],
                    "worst_normalized": worst,
                    "zero_fitted": not fitted,
                    "all_inside": all(
                        c["normalized"] <= 1.0 for c in contacts)})
    rows.sort(key=lambda r: (r["worst_normalized"], r["clock"]))
    data["n_rows"] = len(rows)
    data["ledger"] = rows
    data["pruned"] = pruned
    data["dimensions"] = {"clocks": [c["label"] for c in clock_dims],
                          "weights": [w["rule"] for w in weight_dims],
                          "selections": [s["label"] for s in sel_dims]}

    # ---- branch ---------------------------------------------------------
    winner = next((r for r in rows if r["all_inside"]
                   and r["zero_fitted"]), None)
    branch = ("LEDGER-CONSISTENT-THEORY-CANDIDATE" if winner
              else "LEDGER-CLOSEST-THEORY-DECLARED")
    data["status"] = branch
    closest = rows[0]
    data["closest_row"] = closest
    if winner:
        data["winner"] = winner
    best_zf = next((r for r in rows if r["zero_fitted"]), None)
    data["best_zero_fitted_row"] = best_zf
    checks.append({"name": "branch {}: {} compound rows (cap 200, pruned "
                          "{}); zero-fitted rows {}; closest row [clock "
                          "| weight | selection] worst-normalized "
                          "{:.6f}; C5 passes {} of {} rows (band {:.0e} "
                          "on the booked dz12)".format(
                              branch, len(rows), pruned,
                              sum(1 for r in rows if r["zero_fitted"]),
                              closest["worst_normalized"],
                              sum(1 for r in rows if r["normalized"][4]
                                  <= 1.0), len(rows), BAND_C5),
                   "pass": True})
    print("V318 {}: {} rows; closest: worst-normalized {:.6f} | clock {} "
          "| weight {} | selection {} | five normalized {} | provenance "
          "{} | zero-fitted {}".format(
              branch, len(rows), closest["worst_normalized"],
              closest["clock"], closest["weight"], closest["selection"],
              ["{:.3e}".format(x) for x in closest["normalized"]],
              closest["provenance"], closest["zero_fitted"]), flush=True)
    for line in rows[:20]:
        print("  worst-norm {:>9.3e}  zf {:>5}  {} || {}".format(
            line["worst_normalized"], line["zero_fitted"],
            line["clock"][:58], line["weight"][:48]), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V318_closest_theory_ledger", "tier": 1,
           "property": "the five-contact closest-theory ledger over the "
                       "booked clock x weight x selection cross product, "
                       "with provenance tags and the zero-fitted success "
                       "tier (tier 1, hard, the arbiter)",
           "residual": True, "programme": "theory-hunt V318",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v318_closest_theory_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V318 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
