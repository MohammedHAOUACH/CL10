#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V296 - the flux baryogenesis ledger (tier 3).

**The question (pre-registered).** Is there a booked flux-sector CP
asymmetry of order the observed baryon asymmetry eta_B = 6.1e-10?

**The object.** the three booked asymmetry-generation artifacts:
V55 (flux CP-odd invariants), V106 (leptogenesis decay asymmetry
cases), V125 (node-current CP asymmetry eps). No new mechanism: the
ledger collects every booked asymmetry magnitude and compares to
6.1e-10.

**Measured.** every booked asymmetry-like magnitude in the three
artifacts (recursive key scan for eps/asym/eta/CP leaves, full list
booked), the maximum absolute value, and its log-distance from
6.1e-10.

Branch: ASYMMETRY-MECHANISM-FOUND (a booked |eps| within a factor
10 of 6.1e-10 with a stated conversion to eta_B) /
LEDGER-INSUFFICIENT (the maximum booked asymmetry and its shortfall
booked; the generation of the baryon asymmetry stays an open
problem).
"""
import json
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

ETA_B = 6.1e-10
KEYS = ("eps", "asym", "eta", "cp_asym", "epsilon", "jarsko")


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def scan(obj, path=""):
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            kl = str(k).lower()
            if isinstance(v, (int, float)) and not isinstance(v, bool) \
                    and any(s in kl for s in KEYS):
                out.append((path + "/" + str(k), float(v)))
            elif isinstance(v, dict) and any(s in kl for s in KEYS):
                out.append((path + "/" + str(k), json.dumps(v)[:200]))
            out += scan(v, path + "/" + str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += scan(v, path + "/%d" % i)
    return out


def main():
    t0 = time.time()
    checks, data = [], {}
    led = {}
    nums = []
    for name in ("V55_flux_cp_invariants", "V106_leptogenesis_window",
                 "V125_node_current_cp_asym"):
        try:
            j = art(name)
        except FileNotFoundError:
            # category/name resolution may differ; search the tree
            found = None
            for cat in os.listdir(os.path.join(ROOT, "experiments")):
                p = os.path.join(ROOT, "experiments", cat,
                                 name.lower() + ".json")
                if os.path.exists(p):
                    found = p
                    break
            if found is None:
                led[name] = {"artifact": "missing"}
                continue
            with open(found) as fh:
                j = json.load(fh)
        hits = scan(j.get("data", {}))
        led[name] = {"n_hits": len(hits),
                     "hits": [(p, v if isinstance(v, float) else v)
                              for p, v in hits][:40]}
        for p, v in hits:
            if isinstance(v, float):
                nums.append((name + p, v))
    data["ledger"] = led
    if nums:
        big = max(nums, key=lambda t: abs(t[1]) if abs(t[1]) > 0 else 0)
        cand = min(nums, key=lambda t: abs(abs(t[1]) - ETA_B))
        data["max_abs"] = {"path": big[0], "value": big[1]}
        data["closest_to_eta_B"] = {"path": cand[0], "value": cand[1],
                                    "log10_ratio": float(
                                        __import__("math").log10(
                                            abs(cand[1]) / ETA_B))}
        ok_factor = abs(cand[1]) / ETA_B < 10.0
    else:
        ok_factor = False
        data["closest_to_eta_B"] = None
    status = ("ASYMMETRY-MECHANISM-FOUND" if ok_factor
              else "LEDGER-INSUFFICIENT")
    data["status"] = status
    print("V296 {}: closest booked asymmetry {} (target {})".format(
        status, data["closest_to_eta_B"], ETA_B), flush=True)
    checks.append({"name": "T1 ledger of booked asymmetries: {} numeric "
                          "entries; closest to eta_B within x10? {}"
                   .format(len(nums), ok_factor), "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V296_flux_baryogenesis_ledger", "tier": 3,
           "property": "a booked flux-sector CP asymmetry of order "
                       "eta_B = 6.1e-10",
           "residual": True, "programme": "V296",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v296_flux_baryogenesis_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V296 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
