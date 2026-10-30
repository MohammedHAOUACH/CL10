#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V319 - the rung-weight census, closed class double-target (tier 2,
hard, O8, theory-hunt).

**The question (pre-registered).** V318's closest zero-fitted row fails
exactly one contact, C4 (coefficient at the V310 corrected clock,
worst-normalized 1.205). Does a rung-dependent weight built ONLY of named
constants (LOOP - the booked ladder base; alpha_10 - config conventional),
no dial tuned to a target, land BOTH rung contacts inside their bands at
the zero-fitted corrected clock - and is the census then distinguishable
from density (the V306 discipline, quantitative rule in the
pre-registration)?

**Class (closed, fixed pre-run; V306 rule: targets never inside W).**
A: W(n) = 1 +/- LOOP^k, k = 1..6 (12).
B: W(n) = exp(+/- LOOP^k * n/rung17), k = 1..6, rung17 read from the V197
   ladder rows (12).
C: W(n) = 1 +/- alpha_10^a * LOOP^b, a,b in {1,2}x{1,2} (8).
D: fitted reference (EXCLUDED from the success tier): lambda solved to
   land rung-10 exactly on COEF_OBS.
32 parameter-free members.

**Reading conventions (booked).** coef' = C * W(10)^10; rung-17 leg' =
(S17/3) * W(17) (odd rungs carry ratio 1 exactly - certified odd rule;
V316's leg convention); C1 algebraic-normalized; C2 = rung-10 mass x
W(10) vs the wimp window; C5 = V247 v237_rule at alpha_std(clock),
weight-free.

**T0 (pre-registered).** V206 coef_at_fifth and the V318 closest-row coef
leaf reproduced by the plain reads to < 1e-9; V197/V303 scale legs rel <
ledger_scale_rel_tol.

**Coincidence discipline (fixed pre-run).** Per family: raw hits/total and
the log-uniform expected random hits over the rung-10 window; a member is
SUPPLIER-GRADE only if its family's expected random hits < 0.5, else every
in-band hit is booked as supplied-by-density (the V306 lesson).
Consistency leg (not branch-deciding): the certified F_208 crossing
deficit recomputed under each in-band W.

**Acceptance (pre-registered words).** LEDGER-CONSISTENT-THEORY-CANDIDATE
iff >= 1 supplier-grade parameter-free member's corrected-clock row passes
all five contacts zero-fitted; RUNG-WEIGHT-NO-SUPPLIER otherwise, with the
full ranked table, every in-band hit booked with its family density, and
the certified-crossing consistency.
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
sys.path.insert(0, exp("V247_two_dial_desert_contacts"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from V247_two_dial_desert_contacts import v237_rule       # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import coupling as _cpl                         # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
NMAX = 260                       # V314/V318 plain-read range (verbatim)
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
TOL_SCALE = float(CFG.pillar_surgery.ledger_scale_rel_tol)
BAND_C5 = 1e-9                   # V247's booked desert reproduction band

LOOP = float(u10pred.LOOP)
MP = float(u10pred.M_PLANCK)
COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def plain_coef_at(g):
    """Re-derived rung-10 read of the standard ladder (V314/V318
    verbatim)."""
    G, cartan, Ma, ev = hop_machinery(g)
    F, _rho = loop_ratios(G, cartan, Ma, ev, g, NMAX)
    rung10 = 3.0 * (abs(F[10]) / LOOP) ** 10
    return {"rung10": rung10, "F": {n: abs(F[n])
                                    for n in (10, 208)}}


def _emit(data, checks, ok):
    out = {"id": "V319_rung_weight_double_target", "tier": 2,
           "property": "closed-class rung-weight census against the "
                       "double rung target at the zero-fitted corrected "
                       "clock, with the V306 density discipline (tier 2, "
                       "hard, O8)",
           "residual": True, "programme": "theory-hunt V319",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v319_rung_weight_double_target.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V319 done] {}".format(out["branch"]))


def main():
    t0 = time.time()
    checks, data = [], {}

    # ---- T0a: scale legs re-derived (V314/V318 convention) -----------
    rows197 = art("V197_rung_arithmetic", "data.ladder_rows")
    book = {int(r["rung"]): float(r["value_gev"]) for r in rows197}
    rung17 = max(book)
    scale_derived = {n: MP * LOOP ** n for n in (3, 10, 17)}
    gate = {n: abs(scale_derived[n] - book[n]) / abs(book[n])
            for n in (3, 10, 17)}
    S17 = float(nv.loop_ladder(17)) * float(nv.EV_PER_GEV)
    s17_leaf = float(art("V303_junction_energy_r", "data.contact.S17_eV"))
    gate["S17"] = abs(S17 - s17_leaf) / s17_leaf
    t0a = all(v < TOL_SCALE for v in gate.values())
    data["scale_gates"] = {str(k): v for k, v in gate.items()}
    checks.append({"name": "T0a re-derivation gates (rel < {:.2f}): {}"
                   .format(TOL_SCALE, {str(k): "{:.1e}".format(v)
                                       for k, v in gate.items()}),
                   "pass": bool(t0a)})
    if not t0a:
        data["status"] = "MACHINERY-DRIFT"
        print("V319 MACHINERY-DRIFT: ladder scales do not re-derive",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the two plain reads (T0b) ------------------------------------
    g5 = float(CFG.angle.branches.implemented)
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    hits310 = art("V310_coefficient_expression_census", "data.supplied")
    gap_hits = [h for h in hits310 if "clock_gap" in h["target"]]
    g_cand = g_cross - float(gap_hits[0]["value"])
    C_15 = plain_coef_at(g5)["rung10"]
    leaf15 = float(art("V206_clock_tension", "data.data_clock.coef_at_fifth"))
    C_CORR = plain_coef_at(g_cand)["rung10"]
    leafc = float(art("V318_closest_theory_ledger", "data.closest_row.coef"))
    d15, dc = abs(C_15 - leaf15), abs(C_CORR - leafc)
    t0b = (d15 < 1e-9 and dc < 1e-9)
    checks.append({"name": "T0b plain reads reproduce V206 coef_at_fifth "
                   "{:.15g} (|d| {:.1e}) and the V318 closest-row coef "
                   "{:.15g} (|d| {:.1e}), both < 1e-9".format(
                       leaf15, d15, leafc, dc), "pass": bool(t0b)})
    if not t0b:
        data["status"] = "MACHINERY-DRIFT"
        print("V319 MACHINERY-DRIFT: plain reads do not reproduce the "
              "booked leaves", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- fixed per-clock legs at the corrected clock ------------------
    m_atm = float(CFG.neutrino.m_atm_ev)
    wimp_lo, wimp_hi = [float(x) for x in CFG.physics.wimp_window_gev]
    a_dz12 = float(art("V237_desert_threshold_matching",
                       "data.matching.dz12_efolds"))
    a_shift = float(art("V237_desert_threshold_matching",
                        "data.matching.shift_over_spacing_12"))
    fires5, row5 = v237_rule(_cpl.alpha_std(g_cand))
    c5r = abs(row5["dz12_efolds"] - a_dz12)
    c5 = {"residual": c5r, "band": BAND_C5, "normalized": c5r / BAND_C5,
          "detail": "desert (1,2)-shift {:.15g} vs booked {:.15g} (|d| "
                    "{:.2e}, band {:.0e}); shift/spacing {:.9f} vs "
                    "booked {:.9f}; weight-free".format(
                        row5["dz12_efolds"], a_dz12, c5r, BAND_C5,
                        row5["shift_over_spacing_12"], a_shift)}
    c1 = {"residual": 0.0, "band": None, "normalized": 0.0,
          "detail": "seesaw rung identity, algebraic-normalized (booked)"}
    nu_rel = abs(S17 / 3.0 - m_atm) / m_atm

    # required rung-10 multiplier window (for the density leg)
    wlo = math.log((COEF_OBS - BAND_C) / C_CORR) / 10.0
    whi = math.log((COEF_OBS + BAND_C) / C_CORR) / 10.0

    # certified F_208 leg (consistency, not branch-deciding): the
    # certified crossing lives at the ROOT clock (V209/V276); at g_cand
    # the 208-rung ratio underflows float64 - booked as fact, not used
    g_cert = float(CFG.angle.branches.closed)
    F208 = plain_coef_at(g_cert)["F"][208]
    F208_cand = plain_coef_at(g_cand)["F"][208]
    data["F208_root"] = F208
    data["F208_at_g_cand_float64"] = F208_cand

    a10 = float(CFG.angle.alpha_10)

    def members():
        for s in (1.0, -1.0):
            for k in range(1, 7):
                yield ("A", "1{:+d}*LOOP^{}".format(int(s), k),
                       lambda n, s=s, k=k: 1.0 + s * LOOP ** k,
                       "algebraic (ladder-base power, V276 named input)")
        for s in (1.0, -1.0):
            for k in range(1, 7):
                yield ("B", "exp({:+d}*LOOP^{}*n/{:d})".format(int(s), k,
                                                               rung17),
                       lambda n, s=s, k=k: math.exp(
                           s * LOOP ** k * n / rung17),
                       "algebraic (rung-linear drift, ladder-base power)")
        for s in (1.0, -1.0):
            for a in (1, 2):
                for b in (1, 2):
                    yield ("C", "1{:+d}*a10^{}*LOOP^{}".format(int(s), a,
                                                               b),
                           lambda n, s=s, a=a, b=b: 1.0
                           + s * a10 ** a * LOOP ** b,
                           "algebraic/config (alpha_10 conventional, V52)")

    lam = math.log(COEF_OBS / C_CORR) / 100.0   # W(n)=exp(lam*n) lands
    # rung-10 exactly: (W(10))^10 = exp(lam*100) = COEF/C_CORR. (The
    # pre-registered anchor form exp(lam*(n-10)) is degenerate at n = 10,
    # W(10) = 1 identically - fixed as a booked definition correction,
    # tier-irrelevant: D is excluded from the success tier by rule.)

    scored, families = [], {}
    for fam, rule, W, prov in list(members()) + [
            ("D", "exp(lam*n) fitted lam",
             lambda n: math.exp(lam * n), "fitted (solved)")]:
        w10, w17, w208 = W(10.0), W(float(rung17)), W(208.0)
        coef = C_CORR * w10 ** 10
        n10 = abs(coef - COEF_OBS) / BAND_C
        n17 = abs((S17 / 3.0) * w17 - m_atm) / m_atm / TOL_SCALE
        # C2 with the weight applied to the rung-10 mass
        mdm = scale_derived[10] * w10
        inside = wimp_lo <= mdm <= wimp_hi
        edge2 = (0.0 if inside else
                 min(abs(mdm - wimp_lo), abs(mdm - wimp_hi)) / wimp_lo)
        n2 = edge2 / TOL_SCALE
        cert = 1.0 - F208 * w208 / LOOP
        row = {"family": fam, "rule": rule, "prov": prov,
               "W10": w10, "W17": w17, "coef": coef,
               "norm_rung10": n10, "norm_rung17": n17, "norm_rung10_DM": n2,
               "in_both_bands": bool(n10 <= 1.0 and n17 <= 1.0
                                     and n2 <= 1.0),
               "certified_F208_deficit": cert,
               "worst_normalized": max(n10, n17, n2)}
        scored.append(row)
        families.setdefault(fam, []).append(row)

    # density leg (V306 discipline, quantitative, fixed pre-run)
    density = {}
    for fam, rows in families.items():
        mags = [abs(math.log(r["W10"])) for r in rows if r["W10"] != 1.0]
        if not mags:
            density[fam] = {"members": len(rows),
                            "hits": sum(1 for r in rows
                                        if r["in_both_bands"]),
                            "expected_random_hits": 0.0,
                            "window_lnW10": [wlo, whi],
                            "supplier_grade": False}
            continue
        lo, hi = min(mags), max(mags)
        span = math.log(hi / lo) if lo > 0 and hi > lo else float("inf")
        # log-uniform placement of candidate magnitudes over the family
        # log-range; hit iff inside the (positive) window magnitudes
        p = (math.log(whi / wlo) / span
             if span not in (0.0, float("inf")) and wlo > 0 else 1.0)
        exp_hits = len(mags) * min(p, 1.0)
        hits = sum(1 for r in rows if r["in_both_bands"])
        density[fam] = {"members": len(rows), "hits": hits,
                        "expected_random_hits": exp_hits,
                        "window_lnW10": [wlo, whi],
                        "supplier_grade": bool(exp_hits < 0.5)}

    suppliers = [r for r in scored
                 if r["in_both_bands"] and r["family"] != "D"
                 and density[r["family"]]["supplier_grade"]]

    # full five-contact ledger for supplier rows at the corrected clock
    ledger = []
    for r in suppliers:
        c2 = {"residual": 0.0, "band": TOL_SCALE, "normalized": 0.0,
              "detail": "DM rung-10 x W(10) = {:.4f} GeV vs window: {}"
                        .format(scale_derived[10] * r["W10"],
                                "IN WINDOW" if r["norm_rung10_DM"] == 0.0
                                else "edge")}
        c3 = {"residual": nu_rel * r["W17"], "band": TOL_SCALE,
              "normalized": r["norm_rung17"],
              "detail": "rung-17 leg x W(17) = {:.9f} eV vs {:.4f} (rel "
                        "{:.4%}, band {:.0%})".format(
                            (S17 / 3.0) * r["W17"], m_atm,
                            abs((S17 / 3.0) * r["W17"] - m_atm) / m_atm,
                            TOL_SCALE)}
        c4 = {"residual": abs(r["coef"] - COEF_OBS), "band": BAND_C,
              "normalized": r["norm_rung10"],
              "detail": "coef {:.7f} vs observed {:.7f} (|d| {:.3e}, "
                        "band {:.0e})".format(r["coef"], COEF_OBS,
                                              abs(r["coef"] - COEF_OBS),
                                              BAND_C)}
        contacts = [c1, c2, c3, c4, c5]
        provs = ["algebraic/conventional/measured (V310 census hit)",
                 r["prov"], "trivial"]
        ledger.append({"clock": "V310 corrected g = {:.15g}".format(g_cand),
                       "weight": r["rule"], "selection": "none",
                       "provenance": provs,
                       "zero_fitted": not any("fitted" in p
                                              for p in provs),
                       "normalized": [c["normalized"] for c in contacts],
                       "worst_normalized": max(c["normalized"]
                                               for c in contacts),
                       "all_pass": all(c["normalized"] <= 1.0
                                       for c in contacts)})
    zf_pass = [L for L in ledger if L["all_pass"] and L["zero_fitted"]]

    scored.sort(key=lambda r: (r["worst_normalized"], r["rule"]))
    data["clock_g_cand"] = g_cand
    data["C_CORR"] = C_CORR
    data["C_15"] = C_15
    data["S17_over_3_eV"] = S17 / 3.0
    data["members"] = scored
    data["density"] = density
    data["suppliers"] = [r["rule"] for r in suppliers]
    data["ledger"] = ledger
    data["fitted_reference"] = {"rule": "exp(lam*n)",
                                "lam": lam,
                                "worst_normalized": max(
                                    r["worst_normalized"] for r in scored
                                    if r["family"] == "D")}
    hits_all = sum(1 for r in scored
                   if r["in_both_bands"] and r["family"] != "D")

    branch = ("LEDGER-CONSISTENT-THEORY-CANDIDATE" if zf_pass
              else "RUNG-WEIGHT-NO-SUPPLIER")
    data["status"] = branch
    checks.append({"name": "branch {}: 32 parameter-free members + 1 "
                   "fitted reference; in-band hits (parameter-free) {}; "
                   "supplier-grade {}; zero-fitted full-ledger pass {}; "
                   "density booked {}"
                   .format(branch, hits_all, len(suppliers), len(zf_pass),
                           {f: "{}/{} (E={:.2f})".format(
                               density[f]["hits"], density[f]["members"],
                               density[f]["expected_random_hits"])
                               for f in sorted(density)}),
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    top = scored[0]
    print("  best member {}: {} worst {:.4f}".format(
        top["family"], top["rule"], top["worst_normalized"]), flush=True)
    for r in scored:
        if r["in_both_bands"] and r["family"] != "D":
            print("  IN-BAND {} {} norm10 {:.3f} norm17 {:.3f} "
                  "certified F208 deficit {:.3e} (plain {:.3e})"
                  .format(r["family"], r["rule"], r["norm_rung10"],
                          r["norm_rung17"], r["certified_F208_deficit"],
                          1.0 - F208 / LOOP), flush=True)
    for L in ledger:
        print("  LEDGER {} | {} | {} | norm {} | worst {:.6f} | zf {}"
              .format(L["clock"][:40], L["weight"], L["selection"],
                      ["{:.3e}".format(x) for x in L["normalized"]],
                      L["worst_normalized"], L["zero_fitted"]), flush=True)
    print("  density: {}".format(
        {f: "{}/{} hits, E={:.2f}, grade {}".format(
            density[f]["hits"], density[f]["members"],
            density[f]["expected_random_hits"],
            density[f]["supplier_grade"]) for f in sorted(density)}),
        flush=True)
    data["wall_time"] = time.time() - t0
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
