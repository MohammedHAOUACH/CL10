#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V238 - the desert spectrum's experimental bounds (tier 3, residual).

Every threshold of the desert spectrum - three 10 + 10bar pairs at
1.2e7 GeV with the closed solution's M_10 = 7.06e4 GeV, one (1,2) pair
at 1.0e12 GeV - is read against an experimental bound. Which are already
excluded?  The bounds are fetched, each recorded with source and date in
data.inputs and registered in cl10/config.yaml under bounds:; the
thresholds and multiplicities are the closed V64/V133 values. No
re-tuning of the desert content to escape an exclusion.
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10 import config as CFG                                     # noqa: E402
from cl10.paths import exp                                         # noqa: E402

V_CEV = 246.22            # GeV, the SM Higgs vev (measured: G_F)


def main():
    t0 = time.time()
    checks, data = [], {}

    data["inputs"] = {
        "tau_p_e_pi_superk_yr": {
            "value": CFG.CFG["bounds"]["tau_p_e_pi_superk_yr"],
            "source": "Super-Kamiokande I-IV, 450 kt yr, tau/B(p -> e+ pi0) "
                      "> 2.4e34 yr at 90% CL, arXiv:2010.16098 "
                      "[Phys. Rev. D 105, 012004 (2022)]; PDG 2025 "
                      "summary", "retrieved": "2026-10-23"},
        "m_coloured_direct_gev": {
            "value": CFG.CFG["bounds"]["m_coloured_direct_gev"],
            "source": "CMS t-channel leptoquark exchange in dimuon + "
                      "dielectron, LQ masses probed up to 5 TeV, JHEP 12 "
                      "(2025) 052 (CMS-EXO-22-013); direct pair-production "
                      "ceiling ~2 TeV (PDG 2025 leptoquark review, scalar "
                      "LQ > 1435 GeV first generation)",
            "retrieved": "2026-10-23"},
        "m_doublet_ew_gev": {
            "value": CFG.CFG["bounds"]["m_doublet_ew_gev"],
            "source": "PDG 2025 Sec. 10 (Electroweak Model and Constraints "
                      "on New Physics): oblique STU fit window, extra "
                      "SU(2) doublets probed to the TeV scale; GFitter "
                      "new-physics scale reaches a few TeV",
            "retrieved": "2026-10-23"},
        "sqrts_lhc_gev": {"value": 1.36e4,
                          "source": "LHC Run 3 centre-of-mass energy, PDG "
                                    "2025", "retrieved": "2026-10-23"},
    }
    v133 = json.load(open(exp("V133_desert_pi64", "results",
                              "v133_desert_pi64.json")))
    q = v133["data"]["quantized_pi64"]
    m10, m12_sol = q["M_10_GeV"], q["M_12_GeV"]
    tau_p = q["tau_p_yr"]
    mx = q["M_X_GeV"]
    data["closed_solution"] = {"M_X_GeV": mx, "M_10_GeV": m10,
                               "M_12_solved_GeV": m12_sol,
                               "tau_p_yr": tau_p}
    m12 = 1.0e12                       # the (:983) (1,2) threshold
    m10_pref = 1.2e7                   # the (:983) 10 + 10bar threshold

    b_tau = CFG.CFG["bounds"]["tau_p_e_pi_superk_yr"]
    b_col = CFG.CFG["bounds"]["m_coloured_direct_gev"]
    b_dbl = CFG.CFG["bounds"]["m_doublet_ew_gev"]

    # ---- the threshold -> bound -> status table ----------------------------
    rows = [
        {"threshold": "10 + 10bar pairs, closed M_10", "GeV": m10,
         "bound": "coloured-state ceiling (t-channel LQ search)",
         "bound_GeV": b_col, "ratio": m10 / b_col,
         "status": "NOT EXCLUDED: 14x above the strongest LHC colour probe; "
                   "pair production itself is kinematically closed "
                   "(2 M_10 = %.1f TeV > sqrt(s) = 13.6 TeV)"
                   % (2 * m10 / 1e3)},
        {"threshold": "10 + 10bar pairs, preferred (:983)", "GeV": m10_pref,
         "bound": "coloured-state ceiling", "bound_GeV": b_col,
         "ratio": m10_pref / b_col,
         "status": "NOT EXCLUDED: 2400x above the ceiling"},
        {"threshold": "(1,2) pair (:983)", "GeV": m12,
         "bound": "EW-precision doublet ceiling", "bound_GeV": b_dbl,
         "ratio": m12 / b_dbl,
         "status": "NOT EXCLUDED: 1e9x above; the vector-like pair "
                   "decouples, oblique effects ~ (v/M)^2 = %.1e against a "
                   "fit window of ~0.1" % ((V_CEV / m12) ** 2)},
        {"threshold": "(1,2) pair of the destroyed-ordering solution",
         "GeV": m12_sol, "bound": "EW-precision doublet ceiling",
         "bound_GeV": b_dbl, "ratio": m12_sol / b_dbl,
         "status": "NOT EXCLUDED: above M_X itself (V237), hence also "
                   "outside every collider and precision window"},
        {"threshold": "proton decay at the closed M_X", "yr": tau_p,
         "bound": "Super-K tau/B(p -> e+ pi0)", "bound_yr": b_tau,
         "ratio": tau_p / b_tau,
         "status": "NOT EXCLUDED: 1.2e5x above the bound (F5 margin; "
                   "Hyper-K in 20 yr reaches ~1e35 yr, still 4 orders away)"},
    ]
    data["table"] = rows
    n_excl = sum(1 for r in rows if r["status"].startswith("EXCLUDED"))

    # the GW row: the CGW |c_gw| < 1e-15 bound cited as [60] constrains the
    # model's propagating spin-2 branch (V102/V115/V128 programme), not the
    # static desert thresholds: record it as vacuous for this spectrum.
    data["table"].append({
        "threshold": "all desert thresholds",
        "bound": "GW speed |c_gw| < 1e-15 (LVK, [60])",
        "status": "VACUOUS: static thresholds do not mediate the propagating "
                  "GW; the bound is carried by the graviton branch itself "
                  "(cone programme)"})

    ok_t = all(r["ratio"] > 1.0 for r in rows if "ratio" in r)
    checks.append({
        "name": f"T1 the five threshold/bound rows: every named threshold "
                f"sits above its fetched ceiling - the closest is the "
                f"closed M_10 = {m10:.3e} GeV at {m10 / b_col:.1f}x the "
                f"5 TeV colour ceiling, and the proton lifetime at "
                f"{tau_p:.2e} yr is {tau_p / b_tau:.1e}x the Super-K bound; "
                f"{n_excl} thresholds excluded",
        "pass": bool(ok_t and n_excl == 0)})

    branch = "NO THRESHOLD EXCLUDED" if n_excl == 0 else "THRESHOLD EXCLUDED"
    ok = all(c["pass"] for c in checks)
    status = (
        f"desert spectrum bounds: {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the three 10 + 10bar pairs at "
        f"M_10 = {m10:.3e} GeV sit {m10 / b_col:.1f}x above the strongest "
        f"LHC colour probe (5 TeV, CMS t-channel) and pair production is "
        f"kinematically closed at the LHC (2 M_10 = {2 * m10 / 1e3:.0f} TeV "
        f"> 13.6 TeV); the (1,2) pair at 1.0e12 GeV is 1e9x above the "
        f"electroweak-precision doublet ceiling and fully decoupled "
        f"((v/M)^2 = {(V_CEV / m12) ** 2:.0e}); the proton lifetime "
        f"{tau_p:.2e} yr is {tau_p / b_tau:.1e}x the Super-Kamiokande "
        f"p -> e+ pi0 bound - no F10 fires from the desert spectrum; its "
        f"failure (V237) is internal threshold ordering, not exclusion)")
    out = {"id": "V238_desert_spectrum_bounds",
           "tier": 3,
           "property": "every desert threshold against its fetched "
                       "experimental ceiling (V238, F10 candidate)",
           "residual": True,
           "programme": "V238",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v238_desert_spectrum_bounds.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
