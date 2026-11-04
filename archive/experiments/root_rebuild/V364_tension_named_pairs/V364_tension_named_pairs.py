#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V364 - the tension ledger: is the elected dial beta_elect named by
a pair of machine-named quantities? (booked BEFORE the run in
pre_registrations_root.md, programme B). Target read from the V355
artifact; pool harvested mechanically from 8 artifacts + the 24
provenance rows + 4 audit constants; exact Decimal arithmetic at 40
digits; ops direct, +, -, - (swapped), *, /; relative tol 1e-12
(absolute 1e-15 below 1e-6); chance booked <= 4e-17 expected false
hits. Verdicts: BETA-NAMED-DIRECT / BETA-NAMED-PAIR /
BETA-CONVENTION-CHAIN (hits internal to the build's own chain) /
BETA-NO-NAMED-PAIR (honest negative, pass True) / MACHINERY-DRIFT.
"""
import json
import os
import sys
from decimal import Decimal, getcontext

getcontext().prec = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)

FILES = [
    "archive/experiments/scales_loopfactor/V206_clock_tension/results/v206_clock_tension.json",
    "archive/experiments/couplings_unification/V247_two_dial_desert_contacts/results/v247_two_dial_desert_contacts.json",
    "archive/experiments/theory_hunt/V335_named_clock_roots/results/v335_named_clock_roots.json",
    "archive/experiments/root_rebuild/V350_backreaction_fixed_point/results/v350_backreaction_fixed_point.json",
    "archive/experiments/root_rebuild/V352_oneloop_external_backreaction/results/v352_oneloop_external_backreaction.json",
    "archive/experiments/root_rebuild/V353_inverse_external_sector/results/v353_inverse_external_sector.json",
    "archive/experiments/root_rebuild/V354_elect_read_root_degeneracy/results/v354_elect_read_root_degeneracy.json",
    "archive/experiments/root_rebuild/V355_elect_one_dial/results/v355_elect_one_dial.json",
]


def leaves(obj, path, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            leaves(v, path + "/" + str(k), out)
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            leaves(v, path + "[%d]" % i, out)
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)):
        x = float(obj)
        if x == x and abs(x) < 1e6:
            out.append((path, Decimal(repr(x))))


def main():
    checks, data = [], {}
    pool = []
    missing = [f for f in FILES if not os.path.exists(os.path.join(ROOT, f))]
    for f in FILES:
        p = os.path.join(ROOT, f)
        if not os.path.exists(p):
            continue
        tag = f.split("/")[2]
        d = json.load(open(p))
        if "data" in d:
            leaves(d["data"], tag, pool)
    from cl10 import parameter_audit as pa
    for row in pa.provenance_table()["rows"]:
        try:
            v = Decimal(repr(float(row["value"])))
        except (TypeError, ValueError):
            continue
        pool.append(("provenance:" + row["quantity"], v))
    for cname in ("E", "G10", "M_P", "V_EW"):
        try:
            pool.append(("audit:" + cname, Decimal(repr(float(getattr(pa, cname))))))
        except (TypeError, ValueError):
            pass
    if missing:
        checks.append(("T1 pool harvest: missing artifacts %s" % missing, False))
        _emit(data, checks, False, "MACHINERY-DRIFT")
        return 1

    # target + clone removal
    v355 = json.load(open(os.path.join(ROOT, FILES[-1])))
    target = Decimal(repr(float(v355["data"]["beta_elect"])))
    v354 = json.load(open(os.path.join(ROOT, FILES[-2])))
    anchor = Decimal(repr(float(
        v354["data"]["census"]["J1-identity"]["beta"])))
    anchor_ok = abs(anchor - target) <= Decimal("1e-15") * max(
        Decimal(1), abs(target))
    checks.append(("T0 anchor: beta_elect (V355) == census J1 beta "
                   "(V354) to 1e-15 rel: %s" % anchor_ok, anchor_ok))
    def near_target(x):
        return abs(x - target) <= Decimal("1e-14") * max(Decimal(1),
                                                         abs(target))
    pool = [(p, v) for (p, v) in pool if not near_target(v)]
    n = len(pool)
    pool_ok = (n <= 2000)
    checks.append(("T1 pool census: N = %d leaves from %d artifacts + "
                   "table + constants (booked <= 2000): %s"
                   % (n, len(FILES), pool_ok), pool_ok))
    if not pool_ok:
        _emit(data, checks, False, "MACHINERY-DRIFT")
        return 1

    def hit(opval):
        if abs(opval) < Decimal("1e-6"):
            return abs(opval - target) <= Decimal("1e-15")
        return abs(opval - target) <= Decimal("1e-12") * abs(opval)

    hits = []
    for i, (px, x) in enumerate(pool):
        if hit(x):
            hits.append(("direct", px, "", x))
    for i in range(n):
        px, x = pool[i]
        for j in range(n):
            if i == j:
                continue
            py, y = pool[j]
            for op, val in (("+", x + y), ("-", x - y), ("*", x * y)):
                if hit(val):
                    hits.append((op, px, py, val))
            if y != 0:
                q = x / y
                if hit(q):
                    hits.append(("/", px, py, q))

    def internal(a, b):
        return (("V354" in a or "V355" in a) and
                ("V354" in b or "V355" in b))
    named = [h for h in hits
             if not (h[0] == "direct" or internal(h[1], h[2]))]
    chains = [h for h in hits if h[0] != "direct" and internal(h[1], h[2])]
    direct = [h for h in hits if h[0] == "direct"]
    data["pool_size"] = n
    data["n_tests_bound"] = n * (5 * (n - 1) + 1)
    data["hits"] = [{"op": h[0], "x": h[1], "y": h[2],
                     "value": str(h[3])} for h in hits[:40]]
    data["target"] = str(target)
    if direct:
        status = "BETA-NAMED-DIRECT-" + direct[0][1]
        ok, word_ok = True, True
    elif named:
        h = named[0]
        status = "BETA-NAMED-PAIR-%s:%s:%s" % (h[0], h[1], h[2])
        ok, word_ok = True, True
    elif chains:
        h = chains[0]
        status = ("BETA-CONVENTION-CHAIN-%s:%s:%s (naming negative, "
                  "chain internal to the build)" % (h[0], h[1], h[2]))
        ok, word_ok = True, True
    else:
        status = "BETA-NO-NAMED-PAIR"
        ok, word_ok = True, True
    checks.append(("T2 scan over %d tests: %d hits (%d external-named, "
                   "%d internal-chain, %d direct): %s"
                   % (data["n_tests_bound"], len(hits), len(named),
                      len(chains), len(direct), word_ok), word_ok))
    ok = ok and anchor_ok and pool_ok
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:500])
                    for c in checks), flush=True)
    print("[V364 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0


def _emit(data, checks, ok, status):
    out = {"id": "V364_tension_named_pairs", "tier": 2,
           "residual": True,
           "property": "is the one dial of the elected theory named by "
                       "a pair of machine-named published quantities",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v364_tension_named_pairs.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
