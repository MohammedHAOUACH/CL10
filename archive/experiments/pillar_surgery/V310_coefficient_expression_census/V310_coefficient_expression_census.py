#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V310 - the coefficient expression census over booked constants
(tier 2, hard, O8, pillar-surgery).

**The question (pre-registered).** does a bounded expression over BOOKED
constants supply the -0.55% coefficient correction (multiplicative
target: expr = observed/fifth-clock ratio form AND the direct form
= 2.0765641 +- coef_band_abs) or the clock gap 1.35e-4 booked by V309?

**The object.** the V306 generator pattern (itertools.permutations over
a booked dict, pairs with {+,-,x,/}), extended to depth
pillar_surgery.expr_depth_max = 3: pairs, then one more booked constant.
The class is booked in two layers: CLASS-PRE, the pre-registered ops
{+,-,x,/ both orders} plus the pow^{+-1} reading p^-1 = 1/p, decides
the branch word; CLASS-POW, the superset additionally letting the
exponent be a booked constant (p**c, p**(-c), c**p, guarded), is scored
and booked as annotation only - its power coincidences never decide the
word.  Base constants, every one read from its config path or artifact
(no retyped decimals): LOOP (u10pred.LOOP),
angle.alpha_10, the three angle.branches values, physics.m_planck_gev,
physics.v_electroweak_gev, the loader-derived
hierarchy_coefficient_observed, sm.m_z_gev, sm.alpha1_mz/alpha2_mz/
alpha3_mz, sm.m_unification_obs_gev, sm.delta_b_obs, neutrino
.m_atm_ev/m_sol_ev/sum_mnu_max_ev, physics.omega_dm_h2, V206 artifact
g_obs/g_cross_rung10/coef_at_fifth, V276 artifact rho/F_cross, V206
closure-class coef(pi/4).  A candidate may not restate its own target:
the constants that DEFINE a target (coef_obs for the direct form,
coef_obs/coef_fifth for the ratio form, g_obs/g_cross_rung10 for the
gap) are removed from that target's bases - V306's self-pair exclusion
lifted to the depth-3 class, booked in the artifact.

**Measured (pre-registered).** the five smallest residuals per target
with their expression strings; the targets are re-derived in-experiment:
ratio = coefficient_observed / coef_at_fifth with band coef_band_abs/
coef_at_fifth (the ratio form of the booked coefficient band), direct =
coefficient_observed with band coef_band_abs, gap = g_cross_rung10 -
g_obs (V309's booking recomputed from the V206 artifact) with band
clock_band_abs.

**Acceptance (pre-registered).** COEF-EXPRESSION-SUPPLIED iff some
expression lands inside its target's booked band; COEF-CENSUS-NEGATIVE
otherwise - a bounded no-go booked with the class size and the five
closest misses per target.
"""
import itertools
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402

from cl10 import config as _config  # noqa: E402
from cl10 import u10pred  # noqa: E402

CFG = _config.CFG


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
    L = float(u10pred.LOOP)
    coef_obs = float((u10pred.V_EW / u10pred.M_PLANCK) / L ** 10)
    band_c = float(CFG.pillar_surgery.coef_band_abs)
    band_g = float(CFG.pillar_surgery.clock_band_abs)
    depth = int(CFG.pillar_surgery.expr_depth_max)
    checks, data = [], {}

    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    coef_fifth = float(art("V206_clock_tension",
                           "data.data_clock.coef_at_fifth"))
    coef_pi4 = float(art("V206_clock_tension",
                         "data.closure_class.pi/4.coef"))
    rho_276 = float(art("V276_certified_ladder_208", "data.machinery.rho"))
    f_cross = float(art("V276_certified_ladder_208",
                        "data.machinery.F_cross"))

    B = {
        "LOOP": L,
        "alpha_10": float(CFG.angle.alpha_10),
        "g_impl": float(CFG.angle.branches.implemented),
        "g_closed": float(CFG.angle.branches.closed),
        "g_half": float(CFG.angle.branches.half),
        "M_P": float(CFG.physics.m_planck_gev),
        "v_EW": float(CFG.physics.v_electroweak_gev),
        "coef_obs": coef_obs,
        "m_Z": float(CFG.sm.m_z_gev),
        "alpha1_mz": float(CFG.sm.alpha1_mz),
        "alpha2_mz": float(CFG.sm.alpha2_mz),
        "alpha3_mz": float(CFG.sm.alpha3_mz),
        "M_unif_obs": float(CFG.sm.m_unification_obs_gev),
        "delta_b": float(CFG.sm.delta_b_obs),
        "m_atm": float(CFG.neutrino.m_atm_ev),
        "m_sol": float(CFG.neutrino.m_sol_ev),
        "sum_mnu_max": float(CFG.neutrino.sum_mnu_max_ev),
        "omega_dm_h2": float(CFG.physics.omega_dm_h2),
        "g_obs": g_obs,
        "g_cross10": g_cross,
        "coef_fifth": coef_fifth,
        "coef_pi4": coef_pi4,
        "rho_V276": rho_276,
        "F_cross_V276": f_cross,
    }

    # ---- the three booked targets ---------------------------------------
    gap = g_cross - g_obs
    targets = {
        "coef_ratio(=coef_obs/coef_fifth)": (coef_obs / coef_fifth,
                                             band_c / coef_fifth),
        "coef_direct(=coef_obs)": (coef_obs, band_c),
        "clock_gap(=g_cross10-g_obs)": (gap, band_g),
    }
    # a census may not restate its own target (V306's self-pair exclusion
    # lifted to the depth-3 class): the constants that DEFINE a target are
    # removed from the bases for that target's candidate list
    drop_for = {"coef_ratio(=coef_obs/coef_fifth)": {"coef_obs",
                                                     "coef_fifth"},
                "coef_direct(=coef_obs)": {"coef_obs"},
                "clock_gap(=g_cross10-g_obs)": {"g_obs", "g_cross10"}}
    levels = {}  # per target: (base entries, depth-2 candidates from them)
    for tname in targets:
        keep = {n: v for n, v in B.items() if n not in drop_for[tname]}
        c2 = [(v, n) for n, v in keep.items()]
        for (n1, a), (n2, b) in itertools.permutations(keep.items(), 2):
            c2 += [(a + b, f"({n1}+{n2})"), (a - b, f"({n1}-{n2})"),
                   (a * b, f"({n1}*{n2})"), (a / b, f"({n1}/{n2})")]
        levels[tname] = (keep, c2)

    # ---- depth 3 per target: one more booked constant.  The class is
    # booked in TWO layers, both scored, neither rebalanced after the
    # fact:  CLASS-PRE = the pre-registered ops {+,-,x,/ both orders}
    # plus the pow^{+-1} reading p^-1 = 1/p (p^1 = p reproduces the
    # depth-2 list);  CLASS-POW = the superset additionally allowing the
    # exponent to be a booked constant, p**c, p**(-c), c**p (guarded).
    # The BRANCH word is decided on CLASS-PRE; CLASS-POW hits are booked
    # as annotations only.
    def gen3(keep, c2, strict):
        out = []
        for pv, pn in c2:
            for n2, b in keep.items():
                terms = [(pv + b, f"({pn}+{n2})"), (pv - b, f"({pn}-{n2})"),
                         (b - pv, f"({n2}-{pn})"),
                         (pv * b, f"({pn}*{n2})"),
                         (pv / b, f"({pn}/{n2})") if b != 0.0 else None,
                         (b / pv, f"({n2}/{pn})") if pv != 0.0 else None]
                if pv != 0.0:
                    terms.append((1.0 / pv, f"({pn})^-1"))
                if not strict:
                    for ex, tag in ((b, "^"), (-b, "^-")):
                        try:
                            w = pv ** ex
                        except (OverflowError, ZeroDivisionError,
                                ValueError):
                            w = None
                        if isinstance(w, float) and w != 0.0 \
                                and math.isfinite(w) \
                                and abs(math.log10(abs(w))) < 300:
                            terms.append((w, f"({pn}){tag}{n2}"))
                    if 0.0 < pv:
                        try:
                            w = b ** pv
                        except (OverflowError, ZeroDivisionError,
                                ValueError):
                            w = None
                        if isinstance(w, float) and w != 0.0 \
                                and math.isfinite(w) \
                                and abs(math.log10(abs(w))) < 300:
                            terms.append((w, f"({n2})^({pn})"))
                for e in terms:
                    if e is None:
                        continue
                    v, nm = e
                    if v is not None and math.isfinite(v):
                        out.append((v, nm))
        return out

    scored, scored_pow_only = [], []
    n_pre3 = n_pow3 = 0
    for tname, (tv, tb) in targets.items():
        keep, c2 = levels[tname]
        pre = gen3(keep, c2, strict=True)
        wide = gen3(keep, c2, strict=False)
        n_pre3 += len(pre)
        n_pow3 += len(wide) - len(pre)
        for val, name in list(pre) + list(c2):
            scored.append((abs(val - tv), tname, name, float(val), tb))
        pre_names = {pn for _p, pn in pre}
        for val, name in wide:
            if name not in pre_names:
                scored_pow_only.append((abs(val - tv), tname, name,
                                        float(val), tb))
    scored.sort(key=lambda x: x[0])
    scored_pow_only.sort(key=lambda x: x[0])

    data["targets"] = {t: {"value": tv, "band_abs": tb,
                           "bases_exclude": sorted(drop_for[t])}
                       for t, (tv, tb) in targets.items()}
    data["class"] = {"base_constants": list(B.keys()), "n_base": len(B),
                     "n_candidates_per_target_depth2": {
                         t: len(levels[t][1]) for t in targets},
                     "n_scored_class_pre": len(scored),
                     "n_depth3_pre": n_pre3,
                     "n_depth3_pow_superset_extra": n_pow3,
                     "depth_max": depth,
                     "class_pre_ops": "+,-,x,/ both orders + 1/p "
                                      "(the pow^+-1 reading)",
                     "class_pow_extra": "p**c, p**(-c), c**p (guarded) - "
                                        "annotation only, NOT branch-"
                                        "deciding",
                     "tautology_rule": "the target-defining constants are "
                                       "removed from that target's bases "
                                       "(V306's self-pair exclusion "
                                       "lifted to the class)"}
    per_target = {}
    for t in targets:
        hits = [s for s in scored if s[1] == t][:5]
        per_target[t] = [{"abs_residual": r, "expression": n, "value": v,
                          "band_abs": tb} for r, _, n, v, tb in hits]
    data["closest_per_target"] = per_target
    supplied = [s for s in scored if s[0] <= s[4]]
    data["supplied"] = [{"target": t, "expression": n, "value": v,
                         "abs_residual": r, "band_abs": tb}
                        for r, t, n, v, tb in supplied[:10]]
    supplied_pow = [s for s in scored_pow_only if s[0] <= s[4]]
    data["closest_pow_superset_per_target"] = {
        t: [{"abs_residual": r, "expression": n, "value": v,
             "band_abs": tb} for r, _, n, v, tb
            in [s for s in scored_pow_only if s[1] == t][:5]]
        for t in targets}
    data["supplied_pow_superset"] = [
        {"target": t, "expression": n, "value": v, "abs_residual": r,
         "band_abs": tb} for r, t, n, v, tb in supplied_pow[:10]]
    branch = ("COEF-EXPRESSION-SUPPLIED" if supplied
              else "COEF-CENSUS-NEGATIVE")
    data["status"] = branch

    checks.append({"name": "T0 targets re-derived from booked sources: "
                           "ratio {:.9f}, direct {:.7f}, gap {:.9e}; "
                           "class: {} base constants, {} scored (depth "
                           "{}, {} CLASS-PRE depth-3 terms, +{} CLASS-POW "
                           "annotation terms), no target restated".format(
                               coef_obs / coef_fifth, coef_obs,
                               gap, len(B), len(scored), depth,
                               n_pre3, n_pow3),
                   "pass": True})
    tnames = list(targets)
    checks.append({"name": "T1 minimum residual per target (CLASS-PRE): "
                           "{} -> closest {} at {:.3e}; {} -> {} at "
                           "{:.3e}; {} -> {} at {:.3e}".format(
                               tnames[0],
                               per_target[tnames[0]][0]["expression"],
                               per_target[tnames[0]][0]["abs_residual"],
                               tnames[1],
                               per_target[tnames[1]][0]["expression"],
                               per_target[tnames[1]][0]["abs_residual"],
                               tnames[2],
                               per_target[tnames[2]][0]["expression"],
                               per_target[tnames[2]][0]["abs_residual"]),
                   "pass": True})
    checks.append({"name": "branch {}: {} CLASS-PRE expression(s) inside "
                           "their booked bands; CLASS-POW superset adds "
                           "{} annotated hits (power coincidences, not "
                           "branch-deciding)".format(branch, len(supplied),
                                                     len(supplied_pow)),
                   "pass": True})

    print("V310 {}: over the booked class ({} bases, {} CLASS-PRE "
          "candidates + {} CLASS-POW annotations scored against 3 "
          "targets) the CLASS-PRE closest misses are {} -> ratio "
          "(res {:.2e}, band {:.1e}), {} -> direct (res {:.2e}, band "
          "{:.1e}), {} -> gap (res {:.2e}, band {:.1e}); CLASS-POW "
          "annotated hits in band: {}.".format(
              branch, len(B), len(scored), len(scored_pow_only),
              per_target["coef_ratio(=coef_obs/coef_fifth)"][0]["expression"],
              per_target["coef_ratio(=coef_obs/coef_fifth)"][0]["abs_residual"],
              targets["coef_ratio(=coef_obs/coef_fifth)"][1],
              per_target["coef_direct(=coef_obs)"][0]["expression"],
              per_target["coef_direct(=coef_obs)"][0]["abs_residual"],
              targets["coef_direct(=coef_obs)"][1],
              per_target["clock_gap(=g_cross10-g_obs)"][0]["expression"],
              per_target["clock_gap(=g_cross10-g_obs)"][0]["abs_residual"],
              targets["clock_gap(=g_cross10-g_obs)"][1],
              len(supplied_pow)),
          flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V310_coefficient_expression_census", "tier": 2,
           "property": "does a bounded depth-3 expression over booked "
                       "constants supply the -0.55% coefficient "
                       "correction or the 1.35e-4 clock gap "
                       "(tier 2, hard, O8)",
           "residual": True, "programme": "pillar-surgery V310",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v310_coefficient_expression_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V310 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
