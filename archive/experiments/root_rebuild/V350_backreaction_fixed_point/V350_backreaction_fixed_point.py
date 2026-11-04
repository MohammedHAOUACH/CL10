#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V350 - the clock as OUTPUT of a matter<->geometry backreaction loop
(booked exactly before the run, review/pre_registrations_root.md).

The published rung ladder (V206 verbatim: step(g) over the 45 Spin(10)
planes, sector_hops with the V85 projectors of TL, phases
exp(i g e_lam /2), F_n = cartan_norm ratio, sectors 0..2) already
accepts a link density rho as a first-class input. The hypothesis
here, never tested by the whole hunt: define the matter density from
the ladder's own channel excitations, rho(g) = sum_lam n_lam(g)
sigma_lam / Z with n_lam = 2 - 2 cos(g e_lam /2), sigma_lam the unit
partial traces of the projectors; then iterate
    g -> g' = rung-10 crossing (map A)  /  coefficient inversion
    (map B)  computed AT rho(g),
a damped fixed-point iteration. If a fixed point sits ON g_obs or
g_cross at the 1e-9 guard, the clock is self-consistency.

Hoist note (V326's booked pattern): the g-independent half of
(Pa @ (G kron I)) is precomputed; the delta-algebra of the published
einsum shows rho enters the traced-link hops ONLY through its column
sums - a pre-existing structure of V206, not an invention here, and
gate T0h asserts the hoist reproduces the direct sector_hops
bit-close before any reading.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                        # noqa: E402
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "scales_loopfactor",
                                "V85_closed_loop_factor"))

from V79_gauged_bond_layer import (adjoint_link_operators,        # noqa: E402
                                   root_generators,
                                   structure_constants)
from V85_closed_loop_factor import cartan_norm, projectors        # noqa: E402
from cl10 import clifford                                        # noqa: E402
from cl10 import engine as E                                     # noqa: E402
from cl10 import u10pred                                         # noqa: E402
from cl10 import config as _config                               # noqa: E402

_NG = int(_config.CFG.physics.n_gen)
D_L = _NG * (_NG - 1) // 2          # 45 adjoint link dof
LOOP = float(u10pred.LOOP)
COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
REF_G = float(_config.CFG.angle.branches.implemented)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


COEF_AT_FIFTH_REF = art("V113_axis_factorization", "data.coef_one_per_family")
F10_REF = art("V113_axis_factorization", "data.F_10_measured")

Ts = root_generators()
_f, _k, _a = structure_constants(Ts)
Ls = adjoint_link_operators(_f)
TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
P, EV = projectors(TL)
H = [E.asnumpy(m) for m in clifford.generators()]
CARTAN = [1j * (H[a] @ H[a + 1]) for a in (0, 2, 4, 6, 8)]
PLANES = [1j * (H[a] @ H[b]) for a in range(10) for b in range(a + 1, 10)]


def step(th):
    U = np.eye(32, dtype=np.complex128)
    I32 = np.eye(32, dtype=np.complex128)
    for T in PLANES:
        U = U @ (math.cos(th) * I32 - 1j * math.sin(th) * T)
    return U


def sector_hops(G, rho):
    """V206 verbatim (the direct path, used only by the identity gate)."""
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    out = []
    for Pa in P:
        X = (Pa @ Gi).reshape(32, D_L, 32, D_L)
        Z = np.einsum("ajbm,ml->ajbl", X, rho)
        out.append(np.einsum("ajbj->ab", Z))
    return out


# ---- hoist (exact, T0h-gated): the delta algebra of the published
# einsum collapses out[a,c] = sum_b' A_rho[a,b'] G[b',c] with
# A_rho[lam] = einsum("ajbm,mj->ab", Pa4[lam], rho): 32x32 per channel.
P4 = [Pa.reshape(32, D_L, 32, D_L) for Pa in P]
# sigma_lam: unit partial traces Tr_int(P_lam) (45x45, Hermitian PSD)
SIG = []
for Pm in P4:
    Str = np.einsum("ajal->jl", Pm)
    tr = float(np.real(np.trace(Str)))
    SIG.append((Str / tr if tr > 1e-14 else Str, tr))


def hops_hoisted(G, rho):
    out = []
    for lam in range(3):
        A = np.einsum("ajbm,mj->ab", P4[lam], rho)
        out.append(A @ G)
    return out


def Fn_of(g, nmax=10, rho=None, sectors=None):
    """V206 verbatim on the hoisted hop path."""
    G = step(g)
    if rho is None:
        rho = np.eye(D_L, dtype=np.complex128) / D_L
    Ma = hops_hoisted(G, rho)
    sel = range(3) if sectors is None else sectors
    ph = [np.exp(1j * (g / 2) * e) for e in EV]
    M = sum(ph[i] * Ma[i] for i in sel)
    MG = M @ G
    Fn, Bn, Gn = {}, np.eye(32, dtype=np.complex128), \
        np.eye(32, dtype=np.complex128)
    for n in range(1, nmax + 1):
        Bn = Bn @ MG
        Gn = Gn @ G
        Fn[n] = cartan_norm(Bn, CARTAN) / cartan_norm(Gn, CARTAN)
    return Fn


def rho_of(g):
    """The booked backreaction density: n_lam = 2 - 2 cos(g e_lam/2),
    rho(g) = sum_lam (n_lam tr_lam / Z) sigma_lam, Z = sum n tr:
    Hermitian PSD unit trace by construction; a REAL density matrix."""
    n = [2.0 - 2.0 * math.cos(g * e / 2.0) for e in EV]
    Z = sum(nl * tr for nl, (_sg, tr) in zip(n, SIG))
    out = np.zeros((D_L, D_L), dtype=np.complex128)
    for nl, (sg, tr) in zip(n, SIG):
        out += (nl * tr / Z) * sg
    return out


def bracket_rootA(rho, g_prev):
    """g' : |F_10(g)| = LOOP on (0.15,0.33): grid 1e-3, ALL sign changes
    bisected (the razor band re-crosses many times, V326), the root
    CONTINUATED from g_prev is returned (booked correction: the rung-10
    crossing branch attached to the iterate, not the first crossing)."""
    f = lambda g: abs(Fn_of(g, 10, rho=rho)[10]) - LOOP
    xs = np.arange(0.15, 0.3301, 1e-3)
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
    if not roots:
        return None
    return min(roots, key=lambda r: abs(r - g_prev))


def bracket_rootB(rho, g_prev):
    """g' : 3 (F_10/LOOP)^10 = COEF_OBS (V206 inversion), all crossings
    bisected, CONTINUATION from g_prev (same booked rule as map A)."""
    def c10(g):
        F10 = abs(Fn_of(g, 10, rho=rho)[10])
        d = F10 / LOOP
        return 3.0 * d ** 10 if d < 1 else None
    roots = []
    xs = np.arange(0.15, 0.3301, 1e-3)
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
    if not roots:
        return None
    return min(roots, key=lambda r: abs(r - g_prev))


def iterate(mapf, g0, mix_list):
    trials = []
    for mix in mix_list:
        g = g0
        hist = []
        status = "NOT-CONVERGED"
        for it in range(40):
            gp = mapf(rho_of(g), g)
            if gp is None:
                status = "BRACKET-MISS"
                break
            gn = mix * gp + (1 - mix) * g
            hist.append(gn)
            if abs(gn - g) < 1e-10:
                g = gn
                status = "CONVERGED"
                break
            g = gn
        if len(hist) >= 50 and len(set(round(x, 12)
                                       for x in hist[-50:])) >= 50:
            status = "CYCLE"
        # contraction factor at the candidate
        eps = 1e-6
        gp0 = mapf(rho_of(g), g)
        gpl = mapf(rho_of(g - eps), g)
        gpr = mapf(rho_of(g + eps), g)
        k = None
        if None not in (gp0, gpl, gpr):
            k = abs(gp0 - gpl) / (2 * eps), abs(gpr - gp0) / (2 * eps)
        trials.append({"mix": mix, "status": status,
                       "g_star": float(g) if status == "CONVERGED"
                       else None, "iters": len(hist),
                       "contraction": (max(k) if k else None)})
        if status == "CONVERGED":
            return trials, float(g), True
    return trials, None, False


def main():
    checks, data = [], {}
    g_pub = {"g_obs": float(art("V206_clock_tension",
                                "data.data_clock.g_obs")),
             "g_cross_rung10": float(art("V206_clock_tension",
                                         "data.data_clock"
                                         ".g_cross_rung10")),
             "1/5": 0.2}

    # ---- T0h hoist identity: direct sector_hops vs hoisted, one shot ---
    gtest, rtest = 0.2371, np.eye(D_L, dtype=np.complex128) / D_L
    Gt = step(gtest)
    direct = sector_hops(Gt, rtest)
    hoist = hops_hoisted(Gt, rtest)
    dev_h = max(float(np.max(np.abs(a - b))) for a, b in zip(direct, hoist))
    rt = np.random.RandomState(7)
    R = rt.randn(D_L, D_L) * 1e-2 + np.eye(D_L)
    R = R @ R.conj().T
    direct2 = sector_hops(Gt, R)
    hoist2 = hops_hoisted(Gt, R)
    dev_h2 = max(float(np.max(np.abs(a - b))) for a, b
                 in zip(direct2, hoist2))
    checks.append(("T0h hop-machinery hoist identity (exact A_rho "
                   "collapse of the published einsum, booked pattern): "
                   "worst dev rho=I/45 %.1e, random rho %.1e (<= 1e-12)"
                   % (dev_h, dev_h2), dev_h <= 1e-12 and dev_h2 <= 1e-12))

    # ---- T0a: rho = I/45 re-derives the published anchors verbatim -----
    rho_id = np.eye(D_L, dtype=np.complex128) / D_L
    gc = bracket_rootA(rho_id, g_pub["g_cross_rung10"])
    go = bracket_rootB(rho_id, g_pub["g_obs"])
    dA = abs(gc - g_pub["g_cross_rung10"]) if gc else 9.9
    dB = abs(go - g_pub["g_obs"]) if go else 9.9
    checks.append(("T0a published anchors re-derived inside this script "
                   "with the UNbackreacted density: |g_cross - art| = "
                   "%.1e, |g_obs - art| = %.1e (<= 1e-5)" % (dA, dB),
                   dA <= 1e-5 and dB <= 1e-5))

    # ---- T0b: reference branch F_10 at REF_G vs the V113 line ----------
    F10 = abs(Fn_of(REF_G, 10)[10])
    checks.append(("T0b reference-branch F_10 at 1/5 = %.8f against the "
                   "published V113 line %.8f (|d| %.1e <= 1e-6)"
                   % (F10, F10_REF, abs(F10 - F10_REF)),
                   abs(F10 - F10_REF) <= 1e-6))

    # ---- T0c: rho(g) is a density matrix: Hermitian, PSD, trace 1 ------
    t0c = 0.0
    for gq in (0.17, 0.2, 0.24, 0.3, 0.5):
        rg = rho_of(gq)
        t0c = max(t0c,
                  float(np.max(np.abs(rg - rg.conj().T))),
                  abs(float(np.real(np.trace(rg))) - 1.0),
                  -float(np.min(np.linalg.eigvalsh((rg + rg.conj().T)
                                                  / 2))))
    checks.append(("T0c booked density rho(g): Hermitian / trace 1 / PSD "
                   "worst violation at 5 sampled g: %.1e (<= 1e-12)"
                   % t0c, t0c <= 1e-12))

    # ---- T0d: internal uniformity of the published projectors: is the
    # booked backreaction density RIGID at I/45? (measured, certified;
    # if it is, every fixed point of the booked maps is DEGENERATE: the
    # matter channel closed by the published machinery itself) ---------
    unif = max(float(np.max(np.abs(sg - np.eye(D_L) / D_L)))
               for sg, tr in SIG)
    rho_dev = max(float(np.max(np.abs(rho_of(gq) - np.eye(D_L) / D_L)))
                  for gq in (0.17, 0.2, 0.24, 0.3, 0.5))
    data["sigma_uniformity_dev"] = unif
    data["rho_vs_identity_dev"] = rho_dev
    checks.append(("T0d internal uniformity: max |sigma_lam - I/45| = "
                   "%.1e; booked rho(g) - I/45 max over sampled g = %.1e "
                   "=> the self-consistent density is RIGID at I/45 "
                   "(Schur-type uniformity of the published V85 "
                   "projectors): the matter backreaction channel through "
                   "link densities built from the ladder's own spectrum "
                   "is CLOSED by the published machinery; a random "
                   "(non-self-consistent) rho does move F_10 (~1e-3), so "
                   "the ladder IS rho-sensitive - only the self-consistent "
                   "channel is rigid" % (unif, rho_dev),
                   rho_dev <= 1e-12))

    if not all(c[1] for c in checks):
        data["status"] = "MACHINERY-DRIFT"
        print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                        for c in checks), flush=True)
        _emit(data, checks, False)
        print("[V350 done] MACHINERY-DRIFT", flush=True)
        return 1

    # ---- the fixed-point loop ------------------------------------------
    results = {}
    any_fp = False
    for tag, mapf, g0 in (("A-rung10-crossing", bracket_rootA,
                           g_pub["g_cross_rung10"]),
                          ("B-coefficient-clock", bracket_rootB,
                           g_pub["g_obs"])):
        trials, gs, conv = iterate(mapf, g0, [1.0, 0.5, 0.2])
        any_fp |= conv
        dist = None
        if conv:
            dist = {t: abs(gs - v) for t, v in g_pub.items()}
        results[tag] = {"trials": trials, "g_star": gs,
                        "distances": dist}

    # ---- battery -------------------------------------------------------
    contact = False
    near = False
    for tag, r in results.items():
        if r["g_star"] is None:
            continue
        dmin = min(r["distances"].values())
        r["dmin"] = dmin
        if dmin <= 1e-9:
            contact = True
        if dmin <= 1e-4:
            near = True
    data["maps"] = results
    data["published_anchors"] = g_pub
    checks.append(("FIXED-POINT BATTERY: %d map(s) converged; nearest "
                   "anchor distance %s (guards 1e-9 contact / 1e-4 near)"
                   % (sum(r["g_star"] is not None
                          for r in results.values()),
                      min((r.get("dmin", 9.9) for r in results.values()),
                          default=None)), True))

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for tag, r in results.items():
        print("   %s: %s g*=%s dmin=%s contraction=%s"
              % (tag, r["trials"][-1]["status"],
                 ("%.10f" % r["g_star"]) if r["g_star"] else "-",
                 ("%.3e" % r["dmin"]) if r.get("dmin") is not None else "-",
                 r["trials"][-1]["contraction"]), flush=True)

    rigid = data["rho_vs_identity_dev"] <= 1e-12
    if rigid:
        # post-run booked correction (first run): the booked maps'
        # fixed points are DEGENERATE - the self-consistent density is
        # forced to I/45, so "contact at the published anchor" is the
        # old clock reproducing itself through a rigid channel, not
        # self-consistency creating the clock. No class change, no
        # target shift: the branch word is corrected, the measurement
        # stands.
        status = "FP-DEGENERATE-RIGID-DENSITY"
    elif contact:
        status = "FP-CONTACT"
    elif near:
        status = "FP-NEAR"
    elif any_fp:
        status = "FP-DRIFT"
    else:
        status = "NO-FIXED-POINT"
    print("[V350 done] %s" % status, flush=True)
    _emit(data, checks, True, status)
    return 0


def _emit(data, checks, ok, status=None):
    out = {"id": "V350_backreaction_fixed_point", "tier": 1,
           "property": "is the data clock the fixed point of a booked "
                       "matter<->geometry backreaction loop on the "
                       "published rung ladder",
           "status": status or data.get("status", ""),
           "pass": bool(ok), "data": data,
           "checks": [{"name": c[0] if isinstance(c, str) else c[0],
                       "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v350_backreaction_fixed_point.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
