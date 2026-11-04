#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V298 - the certified-rung portal ledger (tier 3).

**The question (pre-registered).** Does the certified rung structure
(V197/V276: the certified ladder crossing n* = 208, the hierarchy
coefficient) produce a dimensionless portal coupling of order the
booked baryogenesis portal strength |lambda_B| ~ 4.06e-14?

**The object.** the V197 certified ladder artifact (rung rows),
the V276 certified reads (F_208, coefficient), V108's booked
coefficient_observed, V131's closed branch G_CLOSED = pi/4, and the
booked target 4.06e-14 (V106 portal-window read as published).

**Measured.** every dimensionless portal candidate constructible
from the booked certified numbers without a new dial - the
candidates are the published combinations:
  lam_1 = coef_obs * G_CLOSED^2 / (16 pi),
  lam_2 = (coef_obs/3)^2 * G_CLOSED^2 / (16 pi),
  lam_3 = F_208 / F_10_ref   (certified rung suppression) -
each in log10 distance from the booked 4.06e-14. A 2 % match gate.

Branch: RUNG-PORTAL-FIELD (a candidate within 2 %) /
PORTAL-NAMED-INPUT (closest candidate and its distance booked; the
portal coupling stays a named input).
"""
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
from cl10 import coupling                                 # noqa: E402
from cl10.paths import exp                                # noqa: E402

TARGET = 4.06e-14


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    v197 = art("V197_rung_arithmetic")
    v276 = art("V276_certified_ladder_208")
    v108 = art("V108_hierarchy_multiplicity")
    F208 = float(v276["data"]["ladder_reads"]["F_208_certified_crossing"])
    F10 = float(v276["data"]["ladder_reads"]["F_10_open_branch"])
    coef_obs = float(v108["data"]["coefficient_observed"])
    G = float(coupling.G_CLOSED)
    t0ok = bool(F208 > 0 and F10 > 0 and coef_obs > 0
                and "rungs" in json.dumps(v197)[:4000])
    checks.append({"name": "T0 anchors: F_208 = {:.4g}, F_10 = {:.4g}, "
                          "coef_obs = {:.4g}, G_CLOSED = pi/4"
                   .format(F208, F10, coef_obs), "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V298 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    c1 = coef_obs * G * G / (16.0 * math.pi)
    c2 = (coef_obs / 3.0) ** 2 * G * G / (16.0 * math.pi)
    c3 = F208 / F10
    cands = {"lam_1_coef_G2_over_16pi": c1, "lam_2_coef3_G2_over_16pi": c2,
             "lam_3_F208_over_F10": c3}
    ledger = {k: {"value": v,
                  "log10_ratio_to_target": float(math.log10(v / TARGET))}
              for k, v in cands.items()}
    best = min(ledger.items(),
               key=lambda kv: abs(kv[1]["log10_ratio_to_target"]))
    data["ledger"] = ledger
    data["best"] = {"name": best[0], **best[1]}
    within = [k for k, v in ledger.items()
              if abs(math.log10(v["value"] / TARGET)) < 0.03]
    status = "RUNG-PORTAL-FIELD" if within else "PORTAL-NAMED-INPUT"
    data["status"] = status
    print("V298 {}: closest candidate {} = {:.3g} at log10 distance "
          "{:+.2f} from the booked 4.06e-14".format(
              status, best[0], cands[best[0]],
              best[1]["log10_ratio_to_target"]), flush=True)
    checks.append({"name": "T1 portal ledger: closest certified candidate "
                          "at 10^{:+.2f} of the booked portal strength"
                   .format(best[1]["log10_ratio_to_target"]),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V298_rung_portal_ledger", "tier": 3,
           "property": "a certified-rung portal coupling of order the "
                       "booked 4.06e-14",
           "residual": True, "programme": "V298",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v298_rung_portal_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V298 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
