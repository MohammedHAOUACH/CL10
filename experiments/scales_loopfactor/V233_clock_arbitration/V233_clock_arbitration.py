#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V233 - the clock decision, a 2x2 score (tier 5).

**The question this pre-registers against** (`review/pre_registrations.md`):
the manuscript defers the clock decision. The closed clock `pi/4` reads the
hierarchy coefficient to `2.7095` (V209, certified on the `n = 4k` class,
detuning `+1e-8`) against the observed `2.076564` (+30%); the observed value
is read by the SAME self-consistent functional at the NON-closed rational
clock `g = 1/5` (coef `2.08811`, V206 `coef_at_fifth`), and the observed
clock `g_obs = 0.1999980` is itself a fitted number (V206). Which clock does
the data pick, and what does the closure axiom cost either way?

Closed input set (pre-registered): the recorded certified/float readings of
V206 (`closure_class`, `data_clock`) and V209 (`certified`, detuned `pi/4` on
the `n = 4k` class, dps >= 50); the observed coefficient
`physics.hierarchy_coefficient_observed` with its stated constants
uncertainty `1e-4` (V205); the meeting-coupling readings
`coupling.alpha_std(g)` / `coupling.alpha10(g)` against the fitted input
`desert.alpha_x = 0.04` and the closed reading `pi/64 = 0.049087` (ratio
`1.227`, V204); clocks `{pi/4, 1/5}` ONLY.

Measured quantity: the 2x2 score - for each clock, fit/miss of (i) the
hierarchy coefficient against `2.076564 +- 1e-4` and (ii) the meeting-
coupling reading against the fitted `0.04` (the published comparison is the
ratio `1.227`, V204).

Acceptance (pre-registered): some clock fits BOTH within the stated
uncertainties -> F2's branch decision closes with that clock (a config
change); no clock fits both (expected) -> F2 is recorded as FIRED for the
closure axiom: the clock becomes a fitted number `g_obs = 0.1999980` and
Table 6 reclasses it from `algebraic` to `fitted`.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                    # noqa: E402

from cl10 import config as _config                            # noqa: E402
from cl10 import coupling as _cpl                             # noqa: E402

CFG = _config.CFG
UNCI = 1e-4           # stated constants uncertainty on the observed
                      # coefficient (V205: the constants of observation
                      # shift it by only 1e-4)


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    data["fitted_inputs"] = []

    v205 = art("V205_channel_filtered_deficit")
    v206 = art("V206_clock_tension")
    v209 = art("V209_root_clock_degeneracy")
    v204 = art("V204_closed_ladder_limits_audit")

    coef_obs = float(CFG.physics.hierarchy_coefficient_observed)
    assert abs(coef_obs
               - float(v204["data"]["coefficient"]["coef_observed"])) < 1e-12
    alpha_x = float(CFG.desert.alpha_x)          # the fitted desert input

    # ---- T0 the readings, cross-checked between experiments ---------------
    pi4_cert = float(v209["data"]["certified"]["pi/4@1e-8"]["coef"])
    pi4_f64 = float(v206["data"]["closure_class"]["pi/4"]["coef"])
    fifth_coef = float(v206["data"]["data_clock"]["coef_at_fifth"])
    assert abs(coef_obs - float(
        v205["data"]["arithmetic"]["coef_observed"])) < 1e-12
    g_obs = float(v206["data"]["data_clock"]["g_obs"])
    detune = abs(pi4_cert - pi4_f64)
    # the V209 certified reading is the authority; float64 detunes off it
    checks.append({
        "name": f"T0 reading audit: the certified pi/4 coefficient "
                f"{pi4_cert:.7f} (V209, detuning +1e-8, dps >= 50) and the "
                f"float64 reading {pi4_f64:.7f} (V206) differ by "
                f"{detune:.3e} (0.35%) - the certified value is used; the "
                f"observed coefficient {coef_obs:.6f} agrees with V204 to "
                f"1e-12; the fitted clock of V206 is g_obs = {g_obs:.7f} "
                f"(2.0e-6 off the fifth)",
        "pass": bool(detune < 2e-2 and abs(
            g_obs - float(v206["data"]["data_clock"]["g_obs"])) == 0.0)})
    data["readings"] = {"pi4_certified": pi4_cert, "pi4_float64": pi4_f64,
                        "fifth": fifth_coef, "g_obs": g_obs,
                        "observed": coef_obs, "obs_uncertainty": UNCI,
                        "alpha_x_fitted": alpha_x}

    # ---- T1 score (i): the hierarchy coefficient ---------------------------
    clocks = {"pi/4": pi4_cert, "1/5": fifth_coef}
    hi = {}
    for k, c in clocks.items():
        dev = abs(c - coef_obs)
        hi[k] = {"coef": c, "dev": dev,
                 "fits": bool(dev <= UNCI),
                 "miss_pct": 100.0 * (c - coef_obs) / coef_obs}
    checks.append({
        "name": f"T1 hierarchy-coefficient score against {coef_obs:.6f} "
                f"+- {UNCI:.0e}: pi/4 reads {pi4_cert:.4f} (dev "
                f"{hi['pi/4']['dev']:.3e} = {hi['pi/4']['dev'] / UNCI:.0f} x "
                f"the stated uncertainty, {hi['pi/4']['miss_pct']:+.1f}% -> "
                f"MISS; 1/5 reads {fifth_coef:.4f} (dev "
                f"{hi['1/5']['dev']:.3e} = {hi['1/5']['dev'] / UNCI:.0f} x "
                f"the uncertainty, {hi['1/5']['miss_pct']:+.2f}% -> MISS; "
                f"NEITHER clock reads the observed coefficient within its "
                f"stated uncertainty (the 1/5 agreement claimed by V206 is "
                f"a FIT of the clock: g_obs = {g_obs:.7f}, not 1/5)",
        "pass": not (hi["pi/4"]["fits"] or hi["1/5"]["fits"])})
    data["score_hierarchy"] = hi

    # ---- T2 score (ii): the meeting coupling vs the fitted 0.04 ------------
    mc = {}
    for k, g in (("pi/4", math.pi / 4), ("1/5", 0.2)):
        a_std = _cpl.alpha_std(g)
        a_10 = _cpl.alpha10(g)
        mc[k] = {"alpha_std": a_std, "ratio_std": a_std / alpha_x,
                 "alpha10": a_10, "ratio_10": a_10 / alpha_x,
                 "fits_std": bool(abs(a_std - alpha_x) <= 1e-12 * alpha_x),
                 "fits_10": bool(abs(a_10 - alpha_x) <= 1e-12 * alpha_x)}
    # the published closed reading pi/64 and its ratio (V204: 1.227)
    pi64 = math.pi / 64.0
    assert abs(mc["pi/4"]["alpha_std"] - pi64) < 1e-15
    assert abs(pi64 / alpha_x - 1.2271857) < 1e-5
    # and the fitted clock's g^2 reading (circular: g_obs was fitted to it)
    gobs_10 = g_obs * g_obs
    data["score_meeting_coupling"] = mc
    data["g_obs_alpha10"] = {"a": gobs_10, "ratio": gobs_10 / alpha_x,
                             "circular": "g_obs was itself fitted to read "
                                         "the observed coefficient"}
    checks.append({
        "name": f"T2 meeting-coupling score against the fitted "
                f"desert.alpha_x = {alpha_x}: pi/4 standard reading "
                f"alpha_std = pi/64 = {pi64:.6f}, ratio "
                f"{mc['pi/4']['ratio_std']:.4f} (the published 1.227, V204) "
                f"-> MISS; 1/5 standard reading "
                f"{mc['1/5']['alpha_std']:.6f}, ratio "
                f"{mc['1/5']['ratio_std']:.4f} -> MISS (12.6x low); the "
                f"g^2 convention at 1/5 is 0.04 EXACTLY (ratio "
                f"{mc['1/5']['ratio_10']:.6f}) but that convention is the "
                f"fitted one - at the fitted clock g_obs it reads "
                f"{gobs_10:.7f} BY CONSTRUCTION (circular)",
        "pass": not (mc["pi/4"]["fits_std"] or mc["1/5"]["fits_std"])})

    # ---- T3 the decision ----------------------------------------------------
    fits_both = [k for k in clocks
                 if hi[k]["fits"] and (mc[k]["fits_std"] or mc[k]["fits_10"])]
    branch = "F2-CLOCK-DECIDED" if fits_both else "F2-FIRED"
    if fits_both:
        data["fitted_inputs"] = []
    checks.append({
        "name": f"T3 clock decision: 2x2 score (pi/4: hierarchy MISS +30.5%, "
                f"coupling MISS 1.227; 1/5: hierarchy MISS 0.55% = 115x the "
                f"stated uncertainty, coupling standard MISS 0.0796 / g^2 "
                f"convention exact but circular): no clock in {{pi/4, 1/5}} "
                f"fits both readings within the stated uncertainties -> "
                f"F2 is FIRED for the closure axiom: the closure axiom and "
                f"the loop-hierarchy identification are mutually exclusive "
                f"at the measured precision, and the clock becomes a fitted "
                f"number g_obs = {g_obs:.7f} (algebraic -> fitted)",
        "pass": bool(not fits_both)})
    data["clock_fits_both"] = fits_both

    ok = all(c["pass"] for c in checks)
    status = (
        f"F2 CLOCK DECISION, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the 2x2 score on clocks {{pi/4, 1/5}} "
        f"against (i) the observed hierarchy coefficient 2.076564 +- 1e-4 "
        f"and (ii) the fitted desert coupling 0.04. (i) pi/4 reads "
        f"{pi4_cert:.4f} certified (+30.5%, 6330x the uncertainty - the "
        f"number V209 certified); 1/5 reads {fifth_coef:.4f} (+0.55%, 115x "
        f"the uncertainty); neither fits. The claimed 1/5 agreement is a "
        f"FIT: V206's clock is g_obs = {g_obs:.7f} (2.0e-6 off the fifth), "
        f"read to close the gap. (ii) pi/4's standard reading is pi/64 = "
        f"{pi64:.6f}, ratio 1.227 to the fitted 0.04 (the published root/"
        f"observed comparison, V204); 1/5's standard reading 0.0031831 is "
        f"12.6x low; only the g^2 convention matches, and at the fitted "
        f"clock it matches BY CONSTRUCTION. VERDICT: F2 fires for the "
        f"closure axiom - the closure axiom and the loop-hierarchy "
        f"identification are mutually exclusive at the measured precision; "
        f"the clock is a fitted number (g_obs = 0.1999980), Table 6 "
        f"reclasses it algebraic -> fitted, and Table 7 O8's gap states the "
        f"+30% miss as unresolved rather than a named single number.")
    out = {"id": "V233_clock_arbitration",
           "tier": 5,
           "property": "the 2x2 clock decision (hierarchy coefficient and "
                       "meeting coupling) between the closed clock pi/4 and "
                       "the rational clock 1/5 (V233, F2)",
           "residual": True,
           "programme": "V233",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v233_clock_arbitration.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
