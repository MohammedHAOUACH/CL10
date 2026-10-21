#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V206 - the clock tension: closure, the data clock, and the dead zone (O8/O4).

**Why (O8, Section 12.1; and O4).** V203 measured the self-consistent
closed branch: step and traced-link phase both at pi/4, first crossing of
1/(16 pi) at n* = 208 with deficit 0.9895, coefficient 2.700 against the
observed 2.0765 (+30%); V205 showed the deficit is channel-blind and that
the required insertion amplitude equals the 1/5-branch deficit to 0.055%
- the fingerprint that the data reads another clock. This experiment closes
the question MEASURED, with zero free parameters, over the WHOLE closure
class (V131: the internal step is exactly periodic only at pi/2 and pi/4):

1. the closure-class readings of the hierarchy coefficient: pi/4 gives
   2.7002 (+30.0%); pi/2 - the degenerate branch with no mass generation
   (V134) - gives 2.2994 (+10.7%, crossing n* = 102); 3*pi/4 gives 2.7084.
   The closest reading over the whole closure class is 10.7% from the
   observed value: the closure axiom CANNOT reproduce the hierarchy
   coefficient (certificate under the V131 theorem, which exhausts the
   closed angles).
2. the data-clock inversion: the observed coefficient 2.0765 is read by
   the SAME self-consistent functional at rung 10 (n* = 10) at the angle
   g_obs = 0.1999980 - the rung-10 rational clock 1/5 to relative 1e-5.
   The rung-10 crossing of 1/(16 pi) itself sits at g_cross = 0.20013,
   1.3e-4 from 1/5. 1/5 is exactly the O4 clock: alpha_10 = g^2 = 1/25,
   the data-preferred meeting coupling. Two independent observables - the
   desert meeting and the hierarchy coefficient - read the SAME non-closed
   clock.
3. the dead zone and the exactness boundary: at pi/3, pi/5, pi/6, 1/2,
   1/3 the first crossing collapses to n = 3; the exact odd-length
   cancellation holds ONLY on the order-8 closed clock family (F_3, F_5
   below 2e-15 at pi/4 and 3pi/4), blows up at pi/2 (1e17 - the order-4
   branch has no exact ladder) and is tiny but nonzero at 1/5 (0.24).
   The exact ladder exists only on the order-8 closed clock, yet the data
   reads a non-closed clock: this is the named tension.
4. no principled object bridges it: the Cartan-5 link mixture deepens the
   pi/4 miss (0.9915 against 0.9895), roots-only links collapse the ladder
   (n* <= 4), and no single Casimir sector crosses anywhere near (best
   0.047 of the observed).

VERDICT: CLOCK TENSION. The closure axiom and the loop-hierarchy
identification cannot both hold: the closure class misses the observed
coefficient by at least 10.7% (the physical closed branch by 30.1%),
while the observed value is read by the rung-10 clock 1/5 to 1e-5 - the
same 1/5 the O4 meeting coupling prefers. The named residual is the clock
itself: O4 and O8 are ONE tension, not two.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "scales_loopfactor",
                                "V85_closed_loop_factor"))

from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)
from V85_closed_loop_factor import cartan_norm, projectors  # noqa: E402
from cl10 import clifford  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import u10pred  # noqa: E402

D_L = 45
LOOP = 1.0 / (16.0 * math.pi)
COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
PAIRS = [(a, b) for a in range(10) for b in range(a + 1, 10)]
IDX = {p: i for i, p in enumerate(PAIRS)}
CARTAN_PAIRS = [(0, 2), (2, 4), (4, 6), (6, 8), (0, 8)]
IC = [IDX[tuple(sorted(p))] for p in CARTAN_PAIRS]

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
    """Ma_lambda(rho) = Tr_L[(P_lambda (G x I)) (I x rho)]."""
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    out = []
    for Pa in P:
        X = (Pa @ Gi).reshape(32, D_L, 32, D_L)
        Z = np.einsum("ajbm,ml->ajbl", X, rho)
        out.append(np.einsum("ajbj->ab", Z))
    return out


def Fn_of(g, nmax, rho=None, sectors=None):
    G = step(g)
    if rho is None:
        rho = np.eye(D_L, dtype=np.complex128) / D_L
    Ma = sector_hops(G, rho)
    sel = range(3) if sectors is None else sectors
    ph = [np.exp(1j * (g / 2) * e) for e in EV]
    M = sum(ph[i] * Ma[i] for i in sel)
    MG = M @ G
    Fn, Bn, Gn = {}, np.eye(32, dtype=np.complex128), np.eye(32, dtype=np.complex128)
    for n in range(1, nmax + 1):
        Bn = Bn @ MG
        Gn = Gn @ G
        Fn[n] = cartan_norm(Bn, CARTAN) / cartan_norm(Gn, CARTAN)
    return Fn


def first_crossing(Fn, nmax):
    for n in range(2, nmax + 1):
        if 1e-12 < abs(Fn[n]) < LOOP:
            return n, abs(Fn[n])
    return None, None


def reading(Fn, nmax):
    n, F = first_crossing(Fn, nmax)
    d = F / LOOP if F else None
    return n, d, (3.0 * d ** 10 if d else None)


def main():
    checks, data = [], {}

    # ---- T1 the closure-class readings ---------------------------------
    closed = {}
    for tag, g in [("pi/2", math.pi / 2), ("pi/4", math.pi / 4),
                   ("3pi/4", 3 * math.pi / 4)]:
        Fn = Fn_of(g, 600)
        n, d, coef = reading(Fn, 600)
        closed[tag] = {"n_star": n, "deficit": d, "coef": coef,
                       "miss_pct": 100 * (coef / COEF_OBS - 1)}
    pin4, pin2 = closed["pi/4"], closed["pi/2"]
    checks.append({
        "name": f"T1 closure-class readings of the hierarchy coefficient"
                f" (self-consistent, step = phase): pi/4 crosses at"
                f" n* = {pin4['n_star']} deficit {pin4['deficit']:.6f}"
                f" coef {pin4['coef']:.4f} ({pin4['miss_pct']:+.2f}%);"
                f" pi/2 (the V134 degenerate branch: no mass generation)"
                f" crosses at n* = {pin2['n_star']} deficit"
                f" {pin2['deficit']:.6f} coef {pin2['coef']:.4f}"
                f" ({pin2['miss_pct']:+.2f}%); 3pi/4 coef"
                f" {closed['3pi/4']['coef']:.4f}. The closest reading over"
                f" the WHOLE closure class (V131 theorem: pi/2 and pi/4"
                f" only) misses the observed {COEF_OBS:.4f} by"
                f" {abs(pin2['miss_pct']):.2f}%: the closure axiom cannot"
                f" reproduce the coefficient",
        "pass": bool(pin4["n_star"] == 208
                     and abs(pin4["deficit"] - 0.989585) < 1e-4
                     and pin2["n_star"] == 102
                     and abs(pin2["deficit"] - 0.973758) < 1e-4
                     and min(abs(c["miss_pct"]) for c in closed.values()) > 10.0)})
    data["closure_class"] = closed

    # ---- T2 the data-clock inversion -----------------------------------
    def coef10(g):
        F10 = abs(Fn_of(g, 10)[10])
        return 3.0 * (F10 / LOOP) ** 10 if F10 / LOOP < 1 else None

    lo, hi = 0.1985, 0.2
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        c = coef10(mid)
        if c is None or c > COEF_OBS:
            hi = mid
        else:
            lo = mid
    g_obs = 0.5 * (lo + hi)
    c15 = coef10(0.2)

    def f10(g):
        return abs(Fn_of(g, 10)[10])

    lo2, hi2 = 0.2, 0.21
    for _ in range(60):
        mid = 0.5 * (lo2 + hi2)
        if f10(mid) < LOOP:
            lo2 = mid
        else:
            hi2 = mid
    g_cross = 0.5 * (lo2 + hi2)

    checks.append({
        "name": f"T2 data-clock inversion: the observed coefficient is read"
                f" by the SAME functional at rung 10 at g_obs = {g_obs:.7f}"
                f" (bracket width {hi - lo:.1e}); |g_obs - 1/5| ="
                f" {abs(g_obs - 0.2):.1e} - the rung-10 rational clock to"
                f" relative {abs(g_obs - 0.2) / 0.2:.0e}. The rung-10"
                f" crossing of 1/(16 pi) sits at g_cross = {g_cross:.6f},"
                f" {abs(g_cross - 0.2):.1e} above 1/5: the data-preferred"
                f" rational clock is AT the rung-10 crossing. At exactly"
                f" 1/5 the self-consistent reading is coef = {c15:.4f}"
                f" ({100 * (c15 / COEF_OBS - 1):+.2f}%). g_obs is the O4"
                f" clock: alpha_10 = g^2 = 0.0400 = 1/25, the data-"
                f"preferred meeting coupling: TWO independent observables"
                f" read the same non-closed clock",
        "pass": bool(abs(g_obs - 0.2) < 1e-5 and abs(c15 - 2.0881) < 1e-3
                     and abs(g_cross - 0.2) < 1e-3)})
    data["data_clock"] = {"g_obs": g_obs, "g_obs_minus_fifth": g_obs - 0.2,
                          "g_cross_rung10": g_cross, "coef_at_fifth": c15}

    # ---- T3 the dead zone and the exactness boundary --------------------
    dead = {}
    for tag, g in [("pi/3", math.pi / 3), ("pi/5", math.pi / 5),
                   ("pi/6", math.pi / 6), ("1/2", 0.5), ("1/3", 1 / 3)]:
        Fn = Fn_of(g, 30)
        n, d, coef = reading(Fn, 30)
        dead[tag] = {"n_star": n, "deficit": d}
    odd_closed = {}
    for tag, g in [("pi/2", math.pi / 2), ("pi/4", math.pi / 4),
                   ("3pi/4", 3 * math.pi / 4), ("1/5", 0.2)]:
        Fn = Fn_of(g, 6)
        odd_closed[tag] = max(abs(Fn[3]), abs(Fn[5]))
    checks.append({
        "name": f"T3 the dead zone: at {', '.join(dead)} the first crossing"
                f" is n = {sorted(set(v['n_star'] for v in dead.values()))}"
                f" (the loop ladder collapses at length 3). The exact"
                f" odd-length cancellation holds ONLY on the order-8 closed"
                f" clock family: max(|F_3|,|F_5|) = "
                f"{', '.join(f'{t}:{v:.1e}' for t, v in odd_closed.items())}"
                f" - it holds at pi/4 and 3pi/4, blows up at pi/2 (the"
                f" order-4 branch has no exact ladder at all) and is tiny"
                f" but nonzero at 1/5. The exact ladder exists only on the"
                f" order-8 closed clock, yet the data reads a non-closed"
                f" clock",
        "pass": bool(all(v["n_star"] == 3 for v in dead.values())
                     and odd_closed["pi/4"] < 1e-12
                     and odd_closed["3pi/4"] < 1e-12
                     and odd_closed["pi/2"] > 1e6
                     and odd_closed["1/5"] > 1e-6)})
    data["dead_zone"] = dead
    data["odd_selection"] = odd_closed

    # ---- T4 no principled object bridges the gap ------------------------
    rhoC = np.zeros((D_L, D_L), dtype=np.complex128)
    for i in IC:
        rhoC[i, i] = 1.0 / 5
    rhoR = np.eye(D_L, dtype=np.complex128) / D_L
    for i in IC:
        rhoR[i, i] = 0.0
    rhoR = rhoR / 40.0
    variants = {}
    for tag, g in [("pi/4", math.pi / 4), ("1/5", 0.2)]:
        for rtag, rho in [("cartan5", rhoC), ("roots40", rhoR)]:
            Fn = Fn_of(g, 600, rho=rho)
            n, d, coef = reading(Fn, 600)
            variants[f"{tag}:{rtag}"] = {"n_star": n, "deficit": d,
                                         "coef": coef}
    sector_best = 0.0
    for mu in range(3):
        Fn = Fn_of(math.pi / 4, 600, sectors=(mu,))
        _, _, coef = reading(Fn, 600)
        sector_best = max(sector_best, coef or 0.0)
    checks.append({
        "name": f"T4 no principled object bridges it: Cartan-5 link at pi/4"
                f" deepens the miss (deficit"
                f" {variants['pi/4:cartan5']['deficit']:.6f} vs baseline"
                f" 0.989585); roots-only links collapse the ladder"
                f" (n* = {variants['pi/4:roots40']['n_star']}); the best"
                f" single Casimir sector at pi/4 reaches only"
                f" {sector_best / COEF_OBS:.3f} of the observed coefficient",
        "pass": bool(variants["pi/4:cartan5"]["deficit"] > 0.989585
                     and variants["pi/4:roots40"]["n_star"] <= 4
                     and sector_best < 0.1 * COEF_OBS)})
    data["variants"] = variants
    data["sector_best_over_obs"] = sector_best / COEF_OBS

    status = (f"CLOCK TENSION, MEASURED (O8/O4, zero free parameters). "
              f"(i) The closure class ({', '.join(closed)}) reads the "
              f"hierarchy coefficient as {pin4['coef']:.4f} (physical "
              f"branch, {pin4['miss_pct']:+.1f}%) and {pin2['coef']:.4f} "
              f"(degenerate pi/2 branch, {pin2['miss_pct']:+.1f}%): the "
              f"closest reading over the whole closed set misses by "
              f"{abs(pin2['miss_pct']):.1f}% - the closure axiom cannot "
              f"produce the coefficient. (ii) The observed value is read by "
              f"the SAME functional at rung 10 at g_obs = {g_obs:.7f}, the "
              f"rung-10 rational clock 1/5 to relative "
              f"{abs(g_obs - 0.2) / 0.2:.0e}, and the rung-10 crossing sits "
              f"at g_cross = {g_cross:.6f} - {abs(g_cross - 0.2):.1e} from "
              f"1/5: nature's clock IS the rational rung-10 clock, the same "
              f"clock O4's meeting coupling alpha_10 = g^2 = 1/25 prefers. "
              f"(iii) The exact ladder (odd-length cancellation below "
              f"2e-15) exists only on the order-8 closed clock family; at "
              f"pi/2 the odd loops blow up (1e17) and at 1/5 they are tiny "
              f"but nonzero - nature's clock supports no exact ladder. "
              f"(iv) No principled link mixture or sector reading bridges "
              f"the gap. VERDICT: the hierarchy and the desert meeting are "
              f"ONE tension - the closed clock pi/4 is not the clock nature "
              f"reads at the loop level; the closure axiom and the "
              f"loop-hierarchy identification cannot both hold, and the "
              f"named residual is the clock itself.")
    out = {"id": "V206_clock_tension",
           "tier": 4,
           "property": "closure-class readings of the hierarchy coefficient "
                       "vs the data clock: measured one tension (V206)",
           "residual": True,
           "programme": "V206",
           "branch": "W-CLOCK",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": all(c["pass"] for c in checks)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v206_clock_tension.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
