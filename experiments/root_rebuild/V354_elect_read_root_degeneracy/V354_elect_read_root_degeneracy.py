#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V354 - the elected theory's decisive build test at the exact
fifth (booked before the run in review/pre_registrations_root.md).

ONE external dressing, TWO conditions per channel J:
  READ: 3 (|F~_10(1/5; alpha J)| / LOOP)^10 = COEF_OBS;
  ROOT: rootB(1/5, alpha J) = 1/5 (self-consistent born clock).
Degeneracy of the two amplitudes + published-arithmetic census
gates the booked completion V355. GPU 1 via V352's identity-certified
batched scans (T0e 8.8e-16).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V350_backreaction_fixed_point"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))
from cl10.paths import exp  # noqa: E402
import V350_backreaction_fixed_point as V350  # noqa: E402
import V352_oneloop_external_backreaction as V352  # noqa: E402

FIFTH = 0.2


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def readF(alpha, Ju):
    F10 = abs(V352.Fn_corr(FIFTH, alpha * Ju)(10))
    return 3.0 * (F10 / V352.LOOP) ** 10


def main():
    checks, data = [], {}
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    c15_art = float(art("V206_clock_tension",
                        "data.data_clock.coef_at_fifth"))
    zeroJ = np.zeros((32, 32), dtype=np.complex128)

    # ---- T0 tree guards --------------------------------------------------
    r0 = V352.rootB(FIFTH, zeroJ)
    dev_tree = abs(r0 - g_obs) if r0 else 9.9
    c15 = readF(0.0, zeroJ)
    dev_read = abs(c15 - c15_art)
    data["tree_offset_g_obs_minus_fifth"] = g_obs - FIFTH
    checks.append(("T0 tree guards: rootB(1/5, 0) = %.10f lands on the "
                   "published g_obs to %.1e (<= 1e-8: the paper's own "
                   "%+.3e offset re-derived); reading at 1/5 = %.6f vs "
                   "the V206 coef_at_fifth artifact %.6f (|d| %.1e <= "
                   "1e-9)" % (r0 if r0 else -1, dev_tree,
                              g_obs - FIFTH, c15, c15_art, dev_read),
                   dev_tree <= 1e-8 and dev_read <= 1e-9))

    channels = (("J1-identity", np.eye(32, dtype=np.complex128)),
                ("J2-cartan2", sum(C @ C for C in V350.CARTAN)),
                ("J3-step(1/5)", V350.step(FIFTH)))
    res = {}
    for tag, Ju in channels:
        def hR(a, Ju=Ju):
            return readF(a, Ju) - V352.COEF_OBS

        def hT(a, Ju=Ju):
            r = V352.rootB(FIFTH, a * Ju)
            return None if r is None else r - FIFTH

        xs = np.linspace(-2.0, 2.0, 41)
        # READ amplitude(s)
        vs = [hR(x) for x in xs]
        a_read = []
        for i in range(len(xs) - 1):
            if vs[i] * vs[i + 1] < 0.0:
                lo, hi = xs[i], xs[i + 1]
                plo = vs[i]
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    vm = hR(mid)
                    if vm * plo > 0.0:
                        lo, plo = mid, vm
                    else:
                        hi = mid
                a_read.append(0.5 * (lo + hi))
        # ROOT amplitude(s)
        vt = [hT(x) for x in xs]
        a_root = []
        for i in range(len(xs) - 1):
            if vt[i] is None or vt[i + 1] is None:
                continue
            if vt[i] * vt[i + 1] < 0.0:
                lo, hi = xs[i], xs[i + 1]
                for _ in range(30):
                    mid = 0.5 * (lo + hi)
                    vm = hT(mid)
                    vlo = hT(lo)
                    if vm is None or vlo is None or vm * vlo > 0.0:
                        lo = mid
                    else:
                        hi = mid
                a_root.append(0.5 * (lo + hi))
        entry = {"alpha_read": a_read, "alpha_root": a_root}
        # degeneracy pairs
        deg = []
        for ar in a_read:
            for ax in a_root:
                if abs(ar - ax) <= 1e-8:
                    deg.append({"alpha": 0.5 * (ar + ax),
                                "gap": abs(ar - ax)})
        entry["degeneracies"] = deg
        entry["closest_gap"] = min((abs(ar - ax) for ar in a_read
                                    for ax in a_root), default=None)
        res[tag] = entry

    # ---- published-arithmetic census of any degenerate alpha -----------
    eng_ints = [3, 8, 11, 16, 32, 45, 64]
    census = {}
    for tag, e in res.items():
        for d in e["degeneracies"]:
            a = d["alpha"]
            b = a / V350.LOOP
            hits = []
            for q in range(1, 61):
                m = round(a * q)
                if abs(a - m / q) <= 1e-8 and abs(m) <= 60:
                    hits.append("alpha=%d/%d" % (m, q))
                    break
            for t in eng_ints:
                rr = a / t
                if abs(rr - round(rr)) <= 1e-8 and round(rr) != 0:
                    hits.append("alpha=%d x %d" % (round(rr), t))
            census[tag] = {"alpha": a, "beta": b, "hits": hits}
    data["elect"] = res
    data["census"] = census

    checks.append(("V354 per channel: READ amplitudes %s / ROOT "
                   "amplitudes %s / closest gap %s / degeneracies %s"
                   % ({t: [round(x, 8) for x in e["alpha_read"]]
                      for t, e in res.items()},
                     {t: [round(x, 8) for x in e["alpha_root"]]
                      for t, e in res.items()},
                     {t: e["closest_gap"] for t, e in res.items()},
                     {t: e["degeneracies"] for t, e in res.items()}),
                   True))

    deg_hit = any(c.get("hits") for c in census.values())
    deg_any = any(e["degeneracies"] for e in res.values())
    solved_both = all(e["alpha_read"] and e["alpha_root"]
                      for e in res.values())
    if not all(c[1] for c in checks):
        status, ok = "MACHINERY-DRIFT", False
    elif deg_hit:
        status, ok = "ELECT-READ-ROOT-DEGENERATE", True
    elif deg_any:
        status, ok = "ELECT-READ-ROOT-DEGENERATE-NO-ARITH", True
    elif solved_both:
        status, ok = "ELECT-SPLIT", True
    elif deg_any:
        status, ok = "ELECT-READ-ROOT-DEGENERATE-NO-ARITH", True
    else:
        status, ok = "ELECT-NO-READ-OR-NO-ROOT", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for tag, e in res.items():
        print("   %-13s read=%s root=%s closest-gap=%s"
              % (tag, ["%+.8f" % x for x in e["alpha_read"]],
                 ["%+.8f" % x for x in e["alpha_root"]],
                 ("%.3e" % e["closest_gap"])
                 if e["closest_gap"] is not None else "-"), flush=True)
    print("[V354 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V354_elect_read_root_degeneracy", "tier": 1,
           "property": "one dressing, two conditions: can the fifth "
                       "be simultaneously read and born in the "
                       "published arithmetic",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v354_elect_read_root_degeneracy.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
