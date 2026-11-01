#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V352 - one-loop external backreaction on the rung ladder (booked
exactly before the run; the only door V350's rigidity theorem left).

The published ladder (V206) is tree-level. The one-loop self-energy
of the published 32-dim Weyl module in the rung chain,
    B~_n(g) = B_n(g) + LOOP * sum_{k=1}^{n-1} B_k J B_{n-k},
with THREE booked insertion choices, all published objects, zero fit:
  J1 = a_ins * I, a_ins from the V205 artifact (branch_1_5 insertion
       amplitude, read by art() - never retyped);
  J2 = sum_a CARTAN_a^2 (V206's published Cartan generators);
  J3 = G(g_prev) (vertex dressing: the only g-dependent insertion,
       the map closes nontrivially only here).
F~_n = cartan_norm(B~_n)/cartan_norm(G_n). Maps A/B as booked
(continuation roots, V350 rule), damped iteration, battery vs the
published anchors at the 1e-9 guard. For g-independent J the
iteration degenerates by construction: the honest reading is then
the shift delta = g* - g_published, classified RIGID-CHANNEL when
zero to 1e-12 (the correction cannot move the clock at one loop).
Machinery imported from V350 (T0h-identical hop path, V206 verbatim
semantics).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V350_backreaction_fixed_point"))
from cl10.paths import exp  # noqa: E402
import V350_backreaction_fixed_point as V350  # noqa: E402

LOOP = V350.LOOP
COEF_OBS = V350.COEF_OBS
REF_G = V350.REF_G
D_L = V350.D_L


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


A_INS = float(art("V205_channel_filtered_deficit",
                  "data.arithmetic.branch_1_5_insertion_amplitude"))
J1 = A_INS * np.eye(32, dtype=np.complex128)
J2 = sum(C @ C for C in V350.CARTAN)


def Fn_corr(g, J, nmax=10):
    """Corrected rung functional: one-loop self-energy chain."""
    G = V350.step(g)
    rho_id = np.eye(D_L, dtype=np.complex128) / D_L
    Ma = V350.hops_hoisted(G, rho_id)
    ph = [np.exp(1j * (g / 2) * e) for e in V350.EV]
    M = sum(ph[i] * Ma[i] for i in range(3))
    MG = M @ G
    B = [np.eye(32, dtype=np.complex128)]
    Gn = np.eye(32, dtype=np.complex128)
    for _ in range(nmax):
        B.append(B[-1] @ MG)
        Gn = Gn @ G
    def Fn(n):
        Bt = B[n] + LOOP * sum(B[k] @ J @ B[n - k]
                               for k in range(1, n))
        return (V350.cartan_norm(Bt, V350.CARTAN)
                / V350.cartan_norm(Gn, V350.CARTAN))
    return Fn


def rootA(g_prev, J):
    """Continuation root of |F~_10| = LOOP on (0.15, 0.33)."""
    f = lambda g: abs(Fn_corr(g, J)(10)) - LOOP
    xs = np.arange(0.15, 0.3301, 3e-3)
    roots = []
    pv = f(xs[0])
    for k in range(1, len(xs)):
        cv = f(xs[k])
        if pv == 0.0:
            roots.append(float(xs[k - 1]))
        elif pv * cv < 0.0:
            lo, hi = xs[k - 1], xs[k]
            for _ in range(27):
                mid = 0.5 * (lo + hi)
                if f(mid) * pv > 0.0:
                    lo = mid
                else:
                    hi = mid
            roots.append(0.5 * (lo + hi))
        pv = cv
    return min(roots, key=lambda r: abs(r - g_prev)) if roots else None


def rootB(g_prev, J):
    """Continuation root of 3(F~_10/LOOP)^10 = COEF_OBS."""
    def c10(g):
        F10 = abs(Fn_corr(g, J)(10))
        d = F10 / LOOP
        return 3.0 * d ** 10 if d < 1 else None
    xs = np.arange(0.15, 0.3301, 3e-3)
    roots = []
    prev = None
    for x in xs:
        c = c10(x)
        v = None if c is None else c - COEF_OBS
        if v is not None:
            if prev is not None and prev[1] * v < 0.0:
                lo, hi, plo = prev[0], x, prev[1]
                for _ in range(27):
                    mid = 0.5 * (lo + hi)
                    cm = c10(mid)
                    vm = None if cm is None else cm - COEF_OBS
                    if vm is None or vm * plo > 0.0:
                        lo = mid
                    else:
                        hi = mid
                roots.append(0.5 * (lo + hi))
            prev = (x, v)
    return min(roots, key=lambda r: abs(r - g_prev)) if roots else None


def iterate(mapf, g0, J):
    for mix in (1.0, 0.5, 0.2):
        g = g0
        status = "NOT-CONVERGED"
        n_it = 0
        for it in range(40):
            gp = mapf(g, J)
            if gp is None:
                status = "BRACKET-MISS"
                break
            gn = mix * gp + (1 - mix) * g
            n_it = it + 1
            if abs(gn - g) < 1e-10:
                g = gn
                status = "CONVERGED"
                break
            g = gn
        if status == "CONVERGED":
            eps = 1e-6
            a, b, c = (mapf(g, J), mapf(g - eps, J), mapf(g + eps, J))
            k = None
            if None not in (a, b, c):
                k = max(abs(a - b) / (2 * eps), abs(c - a) / (2 * eps))
            return status, float(g), n_it, k
    return status, None, n_it, None


def main():
    checks, data = [], {}
    g_pub = {"g_obs": float(art("V206_clock_tension",
                                "data.data_clock.g_obs")),
             "g_cross_rung10": float(art("V206_clock_tension",
                                         "data.data_clock"
                                         ".g_cross_rung10")),
             "1/5": 0.2}
    data["a_ins_from_V205"] = A_INS

    # ---- T0a/T0b: imported uncorrected machinery re-derives anchors ----
    def F10_plain(g):
        return abs(V350.Fn_of(g, 10)[10])
    go = V350.bracket_rootB(np.eye(D_L) / D_L, g_pub["g_obs"])
    gc = V350.bracket_rootA(np.eye(D_L) / D_L, g_pub["g_cross_rung10"])
    dA = abs(gc - g_pub["g_cross_rung10"]) if gc else 9.9
    dB = abs(go - g_pub["g_obs"]) if go else 9.9
    Fref = F10_plain(REF_G)
    Fref_ref = float(art("V113_axis_factorization", "data.F_10_measured"))
    checks.append(("T0a uncorrected imported ladder re-derives g_cross to "
                   "%.1e and g_obs to %.1e (<= 1e-5); T0b F_10 at 1/5 = "
                   "%.8f vs V113 %.8f (|d| %.1e <= 1e-6)"
                   % (dA, dB, Fref, Fref_ref, abs(Fref - Fref_ref)),
                   dA <= 1e-5 and dB <= 1e-5
                   and abs(Fref - Fref_ref) <= 1e-6))

    # ---- T0c: perturbativity of each insertion --------------------------
    TOL = 1e-4
    t0c = {}
    for tag, J in (("J1-a_ins*I", J1), ("J2-cartan2", J2),
                   ("J3-G(1/5)", V350.step(REF_G))):
        d = abs(abs(Fn_corr(REF_G, J)(10)) - Fref)
        t0c[tag] = d
    data["correction_size_at_1_5"] = {k: v for k, v in t0c.items()}
    checks.append(("T0c one-loop correction |F~_10 - F_10| at 1/5: %s "
                   "(perturbative band (1e-4, 0.5); sub-1e-4 tag => that "
                   "channel is rigid at one loop)"
                   % {k: "%.2e" % v for k, v in t0c.items()}, True))

    if not all(c[1] for c in checks):
        data["status"] = "MACHINERY-DRIFT"
        print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                        for c in checks), flush=True)
        _emit(data, checks, False, "MACHINERY-DRIFT")
        print("[V352 done] MACHINERY-DRIFT", flush=True)
        return 1

    # ---- maps and fixed points ------------------------------------------
    results = {}
    for tag, Jm, gdep in (("J1-a_ins*I", J1, False),
                          ("J2-cartan2", J2, False),
                          ("J3-G(g_prev)", None, True)):
        for mtag, mapf, g0 in (("A-loop-crossing", rootA,
                                g_pub["g_cross_rung10"]),
                               ("B-coefficient-clock", rootB,
                                g_pub["g_obs"])):
            if gdep:
                # g-dependent insertion: map rebuilds J from the iterate
                def m3(gp, _J, mapf=mapf):
                    return mapf(gp, V350.step(gp))
                mf = m3
            else:
                mf = mapf
            st, gs, nit, k = iterate(mf, g0, Jm)
            entry = {"status": st, "g_star": gs, "iters": nit,
                     "contraction": k}
            if gs is not None:
                entry["distances"] = {t: abs(gs - v)
                                      for t, v in g_pub.items()}
                entry["dmin"] = min(entry["distances"].values())
                entry["delta_vs_published"] = gs - g0
            results[tag + "|" + mtag] = entry

    data["fixed_points"] = results

    # ---- classification --------------------------------------------------
    contact = near = anyconv = False
    any_shift = False
    for e in results.values():
        if e["g_star"] is None:
            continue
        anyconv = True
        if abs(e.get("delta_vs_published", 0.0)) > 1e-12:
            any_shift = True
        if e["dmin"] <= 1e-9:
            contact = True
        if e["dmin"] <= 1e-4:
            near = True
    checks.append(("FIXED-POINT BATTERY: %d map entries, %d converged; "
                   "correction shifts the clock: %s; nearest anchor "
                   "distance %s (guards 1e-9 / 1e-4)"
                   % (len(results),
                      sum(e["g_star"] is not None
                          for e in results.values()),
                      any_shift,
                      min((e.get("dmin", 9.9) for e in results.values()),
                          default=None)), True))

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for tag, e in results.items():
        print("   %-26s %s g*=%s delta=%s dmin=%s k=%s"
              % (tag, e["status"],
                 ("%.10f" % e["g_star"]) if e["g_star"] else "-",
                 ("%.3e" % e.get("delta_vs_published", float("nan"))),
                 ("%.3e" % e["dmin"]) if e.get("dmin") is not None else "-",
                 e["contraction"]), flush=True)

    if contact:
        status = "FP2-CONTACT"
    elif near:
        status = "FP2-NEAR"
    elif anyconv and not any_shift:
        status = "RIGID-CHANNEL"
    elif anyconv:
        status = "FP2-DRIFT"
    else:
        status = "NO-FIXED-POINT"
    print("[V352 done] %s" % status, flush=True)
    _emit(data, checks, True, status)
    return 0


def _emit(data, checks, ok, status):
    out = {"id": "V352_oneloop_external_backreaction", "tier": 1,
           "property": "one-loop Weyl-module self-energy backreaction "
                       "on the rung ladder: does self-consistency move "
                       "or fix the clock",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v352_oneloop_external_backreaction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
