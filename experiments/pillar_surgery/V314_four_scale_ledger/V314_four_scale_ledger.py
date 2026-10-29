#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V314 - the four-scale ledger (tier 2, hard, O8/O10, pillar-surgery).

**The question (pre-registered).** assemble every candidate supplied by
V308-V313 into compound scale-theories and score them against the four
measured scale contacts.

**The four contacts (all re-derived in-experiment from config paths and
artifacts, re-derivation gate pillar_surgery.ledger_scale_rel_tol):**
1. seesaw M_P LOOP^3 = 9.6e13 GeV (algebraic, no band issue - the
   V197 rung-3 row re-derived),
2. DM M_P LOOP^10 = 118 GeV vs the observed window
   physics.wimp_window_gev (window membership, V197 rung-10 row),
3. neutrino ladder M_P LOOP^17 and the candidate S17/3 vs
   neutrino.m_atm_ev (the -2.7 % booked gap of V303's contact),
4. the hard one: the rung-10 coefficient vs the loader-derived
   hierarchy_coefficient_observed at the row's clock, band
   pillar_surgery.coef_band_abs.

**Compound theories (rows).** clock rule x coefficient variant, the
candidate lists read from the V308-V313 artifacts: clocks = the fitted
1/5 (angle.branches.implemented), the fitted V206 inversion g_obs, the
V311 p-unitary closure clocks (algebraic, from the V311 artifact), the
V310 census-hit corrected clock g_cross_rung10 - alpha_10*g_impl*
alpha1_mz (algebraic/conventional/measured inputs), and the V313 scan
clocks (conventional).  Coefficient variants = the plain rung-10 ladder
read re-derived at the row clock (V203 machinery), the V308 variant
list, the V312 quaternionic-ladder reads at the p-closure clocks
(artifact), and the standard-ladder crossing read at scan clocks.
Cross product < 50 rows; every cell carries its provenance.

**Measured (pre-registered).** every row scored on all four contacts
(normalized residual: residual / booked band; the algebraic seesaw and
the in-window DM contact score 0), the full ledger table printed sorted
by worst-contact residual.

**Acceptance (pre-registered).** LEDGER-CONSISTENT-THEORY-CANDIDATE iff
some row: all four contacts inside their booked bands AND every input
provenance algebraic/conventional/measured - zero fitted; that row is
then the true formula and the roll-up reclassifies the paper.
LEDGER-BEST-DEFICIT-DECLARED otherwise: the closest row's worst
residual booked with the full table - the current programme's bound on
its best compound theory.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from cl10 import config as _config                         # noqa: E402
from cl10 import neutrino as nv                            # noqa: E402
from cl10 import u10pred                                   # noqa: E402

CFG = _config.CFG
NMAX = 260


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
    t0 = time.time()
    checks, data = [], {}

    LOOP = float(u10pred.LOOP)
    MP = float(u10pred.M_PLANCK)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    BAND_C = float(CFG.pillar_surgery.coef_band_abs)
    TOL_SCALE = float(CFG.pillar_surgery.ledger_scale_rel_tol)
    m_atm = float(CFG.neutrino.m_atm_ev)
    wimp_lo, wimp_hi = [float(x) for x in CFG.physics.wimp_window_gev]

    # ---- contacts 1-3 re-derived, gated against the booked artifacts ----
    rows197 = art("V197_rung_arithmetic", "data.ladder_rows")
    book = {int(r["rung"]): float(r["value_gev"]) for r in rows197}
    scale_derived = {n: MP * LOOP ** n for n in (3, 10, 17)}
    gate = {n: abs(scale_derived[n] - book[n]) / abs(book[n])
            for n in (3, 10, 17)}
    S17 = float(nv.loop_ladder(17)) * float(nv.EV_PER_GEV)
    S3_book = float(art("V303_junction_energy_r", "data.contact.S17_eV"))
    m3_book = float(art("V303_junction_energy_r",
                        "data.contact.m3_onset_eV"))
    gate["S17"] = abs(S17 - S3_book) / S3_book
    gate["m3"] = abs(S17 / 3.0 - m3_book) / m3_book
    t0ok = all(v < TOL_SCALE for v in gate.values())
    data["scale_gates"] = {str(k): v for k, v in gate.items()}
    data["rederived"] = {"seesaw_gev": scale_derived[3],
                         "dm_gev": scale_derived[10],
                         "nu_ladder_gev": scale_derived[17],
                         "S17_eV": S17, "S17_over_3_eV": S17 / 3.0,
                         "wimp_window_gev": [wimp_lo, wimp_hi]}
    checks.append({"name": "T0 re-derivation gates (rel < {:.2f}): {} "
                           "- DM {:.4f} GeV in window [{:.0f}, {:.0f}]".format(
                               TOL_SCALE,
                               {str(k): "{:.1e}".format(v)
                                for k, v in gate.items()},
                               scale_derived[10], wimp_lo, wimp_hi),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V314 MACHINERY-DRIFT: the ladder scales do not re-derive "
              "inside the ledger gate; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- contacts 1-3 values (identical across all rows) -----------------
    c1 = {"residual": 0.0, "band": None, "normalized": 0.0,
          "detail": "seesaw M_P*LOOP^3 = {:.6e} GeV (algebraic, rung A "
                    "= {})".format(scale_derived[3], int(CFG.cell.rung_A))}
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
    nu_rel = abs(S17 / 3.0 - m_atm) / m_atm
    c3 = {"residual": nu_rel, "band": TOL_SCALE,
          "normalized": nu_rel / TOL_SCALE,
          "detail": "S17/3 = {:.9f} eV vs m_atm {:.4f} eV: {:+.4%} "
                    "(the booked V303 gap); bare ladder M_P*LOOP^17 = "
                    "{:.6e} eV rel {:+.2%}".format(
                        S17 / 3.0, m_atm, S17 / 3.0 / m_atm - 1.0,
                        scale_derived[17] * 1e9,
                        scale_derived[17] * 1e9 / m_atm - 1.0)}

    # ---- candidate clocks and coefficient variants from the artifacts ---
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    v310 = json.load(open(os.path.join(
        exp("V310_coefficient_expression_census"), "results",
        "v310_coefficient_expression_census.json")))
    hits310 = v310["data"]["supplied"]
    g5 = float(CFG.angle.branches.implemented)
    unitary_p = [int(p) for p in art("V311_quaternionic_root_step",
                                     "data.unitary_p")]
    sig311 = art("V311_quaternionic_root_step", "data.signatures")
    v312 = json.load(open(os.path.join(
        exp("V312_quaternionic_loop_ladder"), "results",
        "v312_quaternionic_loop_ladder.json")))
    ladder312 = v312["data"].get("ladder_runs", {})
    scan = [float(g) for g in CFG.pillar_surgery.scan_clocks]

    def plain_coef_at(g):
        """Re-derived rung-10 / crossing read of the standard ladder."""
        G, cartan, Ma, ev = hop_machinery(g)
        F, _rho = loop_ratios(G, cartan, Ma, ev, g, NMAX)
        n_star = next((n for n in range(2, NMAX + 1, 2)
                       if abs(F[n]) < LOOP), None)
        rung10 = 3.0 * (abs(F[10]) / LOOP) ** 10
        cross = (None if n_star is None
                 else 3.0 * (abs(F[n_star]) / LOOP) ** 10)
        return {"rung10": rung10, "n_star": n_star, "crossing": cross}

    rows = []
    # (a) fitted fifth clock + V308 plain read (re-derived)
    rd = plain_coef_at(g5)
    rows.append({"clock": "1/5 (angle.branches.implemented)",
                 "clock_g": g5, "clock_prov": "fitted (1/n_qubits "
                                              "reference spectrum)",
                 "coef_expr": "plain rung-10 read 3(F_10/LOOP)^10 at "
                              "the clock (V308 variant a)",
                 "coef_prov": "algebraic machinery, fitted clock",
                 "coef": rd["rung10"]})
    # (b) fitted V206 inversion clock + plain read
    rd_o = plain_coef_at(g_obs)
    rows.append({"clock": "g_obs (V206 data-clock inversion)",
                 "clock_g": g_obs, "clock_prov": "fitted (inverted "
                                                 "from the coefficient)",
                 "coef_expr": "plain rung-10 read at g_obs",
                 "coef_prov": "algebraic machinery, fitted clock",
                 "coef": rd_o["rung10"]})
    # (c) V311 unitary-signature closure clocks + V312 ladder reads
    for p in unitary_p:
        for key, r in ladder312.get(str(p), {}).items():
            if r["coef"] is None:
                continue
            rows.append({"clock": "Cl(0,10)-family p={:d} closure {} "
                                  "(k={:d})".format(
                                      p, key,
                                      sig311[str(p)]["closure_clocks"][key]
                                      ["k"]),
                         "clock_g": r["clock"],
                         "clock_prov": "algebraic (V311 closure census)",
                         "coef_expr": "quaternionic-root ladder crossing "
                                      "read at {} (V312)".format(key),
                         "coef_prov": "algebraic (signature ladder)",
                         "coef": r["coef"]})
    # (d) V310 clock-gap corrected clock + rung-10 read there
    gap_hits = [h for h in hits310
                if "clock_gap" in h["target"]]
    if gap_hits:
        h = gap_hits[0]
        g_cand = g_cross - h["value"]
        rd_c = plain_coef_at(g_cand)
        rows.append({"clock": "g_cross_rung10 - ({}), the V310 "
                              "clock-gap hit".format(h["expression"]),
                     "clock_g": g_cand,
                     "clock_prov": "algebraic/conventional/measured "
                                   "(V310 CLASS-PRE census hit)",
                     "coef_expr": "plain rung-10 read at the corrected "
                                  "clock",
                     "coef_prov": "algebraic machinery",
                     "coef": rd_c["rung10"],
                     "gap_expr": h["expression"],
                     "gap_hit_value": h["value"],
                     "gap_hit_residual": h["abs_residual"]})
    # (e) V313 scan clocks (dictionary candidates) + standard crossing read
    for g in scan:
        rd_s = plain_coef_at(g)
        if rd_s["crossing"] is None:
            continue
        rows.append({"clock": "scan clock g = {:.9g} (V313 "
                              "dictionary scan)".format(g),
                     "clock_g": g,
                     "clock_prov": "conventional (scan set)",
                     "coef_expr": "standard-ladder crossing read at "
                                  "the clock",
                     "coef_prov": "algebraic machinery, scan clock",
                     "coef": rd_s["crossing"]})

    # ---- score every row -------------------------------------------------
    scored = []
    for r in rows:
        c4_res = abs(r["coef"] - COEF_OBS)
        c4 = {"residual": c4_res, "band": BAND_C,
              "normalized": c4_res / BAND_C,
              "detail": "coef {} vs observed {:.7f} (|d| {:.2e}, band "
                        "{:.0e})".format("{:.7f}".format(r["coef"]),
                                         COEF_OBS, c4_res, BAND_C)}
        norm = max(c1["normalized"], c2["normalized"], c3["normalized"],
                   c4["normalized"])
        fitted = "fitted" in r["clock_prov"] or "fitted" in r["coef_prov"]
        scored.append({"row": r, "contacts": [c1, c2, c3, c4],
                       "worst_normalized": norm,
                       "zero_fitted": not fitted})
    scored.sort(key=lambda s: (s["worst_normalized"], s["row"]["clock"]))
    data["coefficient_observed"] = COEF_OBS
    data["band_coef"] = BAND_C
    data["neutrino_contact_rel"] = nu_rel
    ledger = []
    for s in scored:
        ledger.append({
            "clock": s["row"]["clock"],
            "clock_g": s["row"]["clock_g"],
            "clock_prov": s["row"]["clock_prov"],
            "coef_expr": s["row"]["coef_expr"],
            "coef": s["row"]["coef"],
            "contacts": [c["detail"] for c in s["contacts"]],
            "normalized": [c["normalized"] for c in s["contacts"]],
            "worst_normalized": s["worst_normalized"],
            "zero_fitted": s["zero_fitted"]})
    data["ledger"] = ledger
    for line in ledger:
        print("  {:>62s}  coef {:>12.7f}  worst-norm {:>9.3e}  "
              "fitted: {}".format(line["clock"][-62:], line["coef"],
                                  line["worst_normalized"],
                                  not line["zero_fitted"]), flush=True)

    # ---- branch ----------------------------------------------------------
    winner = None
    for s in scored:
        all_inside = (s["contacts"][0]["residual"] == 0.0
                      and s["contacts"][1]["normalized"] < 1.0
                      and s["contacts"][2]["normalized"] < 1.0
                      and s["contacts"][3]["residual"] < BAND_C)
        if all_inside and s["zero_fitted"]:
            winner = s
            break
    branch = ("LEDGER-CONSISTENT-THEORY-CANDIDATE" if winner
              else "LEDGER-BEST-DEFICIT-DECLARED")
    data["status"] = branch
    if winner:
        data["winner"] = {"clock": winner["row"]["clock"],
                          "coef_expr": winner["row"]["coef_expr"],
                          "worst_normalized": winner["worst_normalized"]}
    else:
        best_zf = next((s for s in scored if s["zero_fitted"]), None)
        best_any = scored[0]
        data["best_zero_fitted_row"] = {
            "clock": best_zf["row"]["clock"],
            "coef": best_zf["row"]["coef"],
            "worst_normalized": best_zf["worst_normalized"],
            "binding_contact": max(
                range(4), key=lambda i: best_zf["contacts"][i]
                ["normalized"]),
            "contacts": [c["detail"] for c in best_zf["contacts"]]} \
            if best_zf else None
        data["best_any_row"] = {
            "clock": best_any["row"]["clock"],
            "coef": best_any["row"]["coef"],
            "worst_normalized": best_any["worst_normalized"]}
    checks.append({"name": "branch {}: {} rows scored, best zero-fitted "
                           "worst-normalized {} (binding: neutrino rel "
                           "gap {} normalized {:.3f})".format(
                               branch, len(scored),
                               "n/a" if not data.get(
                                   "best_zero_fitted_row") else
                               "{:.3e}".format(
                                   data["best_zero_fitted_row"]
                                   ["worst_normalized"]),
                               "{:.4%}".format(nu_rel),
                               c3["normalized"]),
                   "pass": True})
    print("V314 {}: ledger over {} compound rows; the neutrino contact "
          "({:+.2%}) alone gives every row a normalized {:.3f}; best "
          "zero-fitted row {}.".format(
              branch, len(scored), S17 / 3.0 / m_atm - 1.0,
              c3["normalized"],
              "none" if not data.get("best_zero_fitted_row") else
              data["best_zero_fitted_row"]["clock"]), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V314_four_scale_ledger", "tier": 2,
           "property": "compound four-scale contact ledger over the "
                       "candidate clocks and coefficient variants "
                       "supplied by V308-V313 (tier 2, hard, O8/O10)",
           "residual": True, "programme": "pillar-surgery V314",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v314_four_scale_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V314 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
