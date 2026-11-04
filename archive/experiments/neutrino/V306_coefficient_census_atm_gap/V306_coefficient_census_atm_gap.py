#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V306 - does any booked coefficient supply the atmospheric gap?
(tier 2, hard, O2.)
**The question (pre-registered).** The energy-selected point sits
-2.70% under the atmospheric scale (V303): m_3(2) = S17/3 = 0.0487 eV
against M_ATM_EV = 0.0501 eV, i.e. r_atm = 1.9191154 against the onset
r* = 2. Does ANY booked constant, or a bounded combination of two,
supply the gap coefficient to 1e-3 absolute? The candidate class was
fixed before the run: base constants read from their artifacts/config
(no retyped decimals), every a, a +/- b, a*b, a/b for distinct a,b;
both target forms scored (the multiplicative rho_req and the additive
Delta-r = 2 - r_atm in r units).
**Measured.** The five smallest residuals with their expressions; the
exact minimum; the class size. Branch: R-COEFFICIENT-SUPPLIED (min
residual <= 1e-3 - the supplying expression booked, and it owes a
derivation test before any reclassification: a two-number coincidence
is not a mechanism) / R-COEFFICIENT-NOT-SUPPLIED (bounded no-go: class
size and closest misses booked).
"""
import glob
import itertools
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")



HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402
from cl10 import u10pred                                  # noqa: E402
from cl10 import config as _cfg                           # noqa: E402


def art(name):
    hits = glob.glob(os.path.join(exp(name, "results"), "*.json"))
    assert len(hits) == 1, (name, hits)
    with open(hits[0]) as fh:
        return json.load(fh)


def find_key(obj, key):
    """Booked value located by key name in a verdict tree (provenance by
    name, not by position)."""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            r = find_key(v, key)
            if r is not None:
                return r
    elif isinstance(obj, list):
        for v in obj:
            r = find_key(v, key)
            if r is not None:
                return r
    return None


def main():
    t0 = time.time()
    checks, data = [], {}
    S17 = nv.loop_ladder(17) * nv.EV_PER_GEV
    m_atm = nv.M_ATM_EV
    r_atm = S17 / m_atm - 1.0
    d_r = 2.0 - r_atm
    rho_req = m_atm / (S17 / 3.0) - 1.0

    v276 = art("V276_certified_ladder_208")
    v206 = art("V206_clock_tension")
    rho = float(v276["data"]["machinery"]["rho"])
    f_cross = float(v276["data"]["machinery"]["F_cross"])
    f_open = float(v276["data"]["F_10_open_branch"])
    g_obs = float(find_key(v206["data"], "g_obs"))
    from cl10 import coupling
    a10 = float(coupling.ALPHA_10)
    assert abs(a10 - 0.04) < 1e-12, a10

    L = float(u10pred.LOOP)
    B = {
        "LOOP": L, "LOOP^2": L * L, "LOOP^3": L ** 3, "LOOP^4": L ** 4,
        "alpha_10": a10, "pi/64": math.pi / 64.0,
        "g_obs-1/5": g_obs - 0.2,
        "1-rho_V276": 1.0 - rho,
        "F_cross/F_10_open-1": f_cross / f_open - 1.0,
        "F_10_open/LOOP-1": f_open / L - 1.0,
        "1-F_cross/LOOP": 1.0 - f_cross / L,
        "LOOP-F_10_open": L - f_open,
        "1/208": 1.0 / 208.0, "208": 208.0,
        "1/2": 0.5, "1/3": 1.0 / 3.0, "1/4": 0.25, "1/5": 0.2, "1/8": 0.125,
    }
    targets = {"rho_req": rho_req, "Delta_r": d_r}
    cands = []
    for name, a in B.items():
        cands.append((a, name))
    for (n1, a), (n2, b) in itertools.permutations(B.items(), 2):
        cands += [(a + b, f"{n1}+{n2}"), (a - b, f"{n1}-{n2}"),
                  (a * b, f"{n1}*{n2}"), (a / b, f"{n1}/{n2}")]
    scored = []
    for tname, tv in targets.items():
        for val, name in cands:
            scored.append((abs(val - tv), tname, name, float(val)))
    scored.sort(key=lambda x: x[0])
    data["targets"] = {"S17_eV": float(S17), "M_ATM_EV": float(m_atm),
                       "r_atm": float(r_atm), "Delta_r": float(d_r),
                       "rho_req": float(rho_req)}
    data["class"] = {"base_constants": list(B.keys()), "n_base": len(B),
                     "n_candidates": len(cands),
                     "n_scored": len(scored), "tolerance": 1e-3}
    top = [{"abs_residual": r, "target": t, "expression": n, "value": v}
           for r, t, n, v in scored[:5]]
    data["closest"] = top
    supplied = [c for c in scored if c[0] <= 1e-3]
    data["supplied"] = [{"target": t, "expression": n, "value": v,
                         "abs_residual": r} for r, t, n, v in supplied]
    branch = ("R-COEFFICIENT-SUPPLIED" if supplied
              else "R-COEFFICIENT-NOT-SUPPLIED")
    checks.append({"name": "T0 targets re-derived from booked sources: "
                          "r_atm = %.7f, Delta_r = %.7f, rho_req = %.6f; "
                          "class: %d base constants, %d expressions, %d "
                          "scored" % (r_atm, d_r, rho_req, len(B), len(cands),
                                      len(scored)), "pass": True})
    checks.append({"name": "T1 min residual %.2e (%s for %s): %s"
                          % (scored[0][0], scored[0][2], scored[0][1],
                             branch),
                   "pass": True})
    data["status"] = branch
    print("V306 {}: over the booked class ({} base constants, {} "
          "expressions x 2 targets) the closest misses are {} -> "
          "Delta_r = {:.7f} (residual {:.2e}) and {} -> rho_req = "
          "{:.6f} (residual {:.2e}); no booked coefficient pair supplies "
          "the 2.7 % atmospheric gap to 1e-3{}".format(
              branch, len(B), len(cands),
              top[0]["expression"], d_r, top[0]["abs_residual"],
              top[1]["expression"], rho_req, top[1]["abs_residual"],
              "" if not supplied else
              " EXCEPT: " + "; ".join(
                  "{} vs {} (res {:.1e})".format(s[2], s[1], s[0])
                  for s in supplied[:4])), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V306_coefficient_census_atm_gap", "tier": 2,
           "property": "booked-coefficient census of the -2.70 % "
                       "atmospheric gap at the energy-selected point",
           "residual": True, "programme": "V306",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v306_coefficient_census_atm_gap.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V306 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
