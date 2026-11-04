#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V297 - the induced-action ledger (tier 3).

**The question (pre-registered).** Does any booked certified-ladder
loop factor convert through the certified normalization into an
induced gravitational action coefficient within 2 % of the booked
Planck normalization?

**The object.** the V283 metric-trace ladder rows (the same
closed-branch hop machinery, n, N grid, F values as booked in the
V283 artifact), the V276 certified crossing read F_208 and
coefficient_observed (V108).

**Measured.** for every booked ladder row with a stable crossing:
the ratio F_row / F_208 (the dimensionless induced-action
coefficient in the certified normalization), its deviation from the
booked Planck-side ratio (coefficient_observed / 3: the published
hierarchy normalization carried as the Planck-side number); the
minimum over rows; the 2 % gate.

Branch: TRACE-ACTION-FOUND (a row within 2 %) /
TRACE-ACTION-ABSENT (closest row and its deviation booked; the
graviton-to-Einstein step stays a named open gap).
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


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def find_rows(obj, out):
    """ladder-shaped rows: dicts with an n and an N key."""
    if isinstance(obj, dict):
        if "n" in obj and "N" in obj:
            out.append(obj)
        for v in obj.values():
            find_rows(v, out)
    elif isinstance(obj, list):
        for v in obj:
            find_rows(v, out)


def main():
    t0 = time.time()
    checks, data = [], {}
    v283 = art("V283_metric_trace_loop")
    v276 = art("V276_certified_ladder_208")
    v108 = art("V108_hierarchy_multiplicity")
    F208 = float(v276["data"]["ladder_reads"]["F_208_certified_crossing"])
    coef_obs = float(v108["data"]["coefficient_observed"])
    planck_ratio = coef_obs / 3.0
    rows = []
    find_rows(v283["data"], rows)
    # the V283 ladder is booked as {Nkey: {n: F}}
    lad = v283["data"].get("ladder", {})
    if isinstance(lad, dict):
        for Nk, nv in lad.items():
            if not isinstance(nv, dict):
                continue
            Nlab = int("".join(ch for ch in str(Nk) if ch.isdigit())
                       or Nk)
            for nk, F in nv.items():
                if isinstance(F, (int, float)):
                    rows.append({"n": int(nk), "N": int(Nlab),
                                 "F": float(F)})
    rows = [r for r in rows if isinstance(r.get("F"), (int, float))
            or isinstance(r.get("F_n"), (int, float))
            or isinstance(r.get("Fval"), (int, float))]

    def Fof(r):
        for k in ("F", "F_n", "Fval", "F_cross"):
            if isinstance(r.get(k), (int, float)):
                return float(r[k])
        return None

    rows = [r for r in rows if Fof(r) is not None]
    t0ok = bool(rows) and F208 > 0
    checks.append({"name": "T0 anchors: {} ladder rows in V283, "
                          "F_208 = {:.3e}".format(len(rows), F208),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V297 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    ledger = []
    for r in rows:
        F = Fof(r)
        ratio = F / F208
        dev = abs(ratio - planck_ratio) / planck_ratio
        ledger.append({"n": int(r["n"]), "N": int(r["N"]),
                       "F": F, "ratio": ratio, "dev": dev})
    ledger.sort(key=lambda r: r["dev"])
    best = ledger[0]
    data["ledger_top"] = ledger[:10]
    data["planck_ratio"] = planck_ratio
    data["best"] = best
    status = ("TRACE-ACTION-FOUND" if best["dev"] < 0.02
              else "TRACE-ACTION-ABSENT")
    data["status"] = status
    print("V297 {}: best row n={} N={} ratio {:.4g} vs booked {:.4g} "
          "({:.1%})".format(status, best["n"], best["N"], best["ratio"],
                            planck_ratio, best["dev"]), flush=True)
    checks.append({"name": "T1 induced-action ledger: closest certified "
                          "row {:.1%} from the Planck-side ratio"
                   .format(best["dev"]), "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V297_trace_action_ledger", "tier": 3,
           "property": "a certified ladder factor converting to the "
                       "booked induced-action (Planck) normalization",
           "residual": True, "programme": "V297",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v297_trace_action_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V297 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
