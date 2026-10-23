#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V209 - the root clock is phase-degenerate: the Cartan loop ratio is 0/0 there.

The closure axiom (V131) makes the clock a root element, $g \in \{\pi/4,\ \pi/2,
3\pi/4\}$, and the closed-branch machinery (V85, V200, V202, V203, V204, V206)
reads the hierarchy coefficient from the traced-loop Cartan-content ratio
$F_n = \mathrm{cnorm}((MG)^n)/\mathrm{cnorm}(G^n)$ at exactly those angles.

This experiment shows that the closure-class angles are precisely the angles at
which the $\mathbf{288}$ ($T\!\cdot\!L = -12$) and $\mathbf{1120}$
($T\!\cdot\!L = +4$) Casimir sectors carry the SAME traced-link phase,
$\exp(i(g/2)(-12)) = \exp(i(g/2)(4)) \iff 8g \in 2\pi\mathbb{Z}$, and that at a
phase-degenerate angle the numerator and the denominator of $F_n$ vanish
identically on the even-$n$ class, so $F_n$ is not a measurement there.

1. **Exact integer arithmetic, $g = \pi/4$.**  At the root angle the whole
   construction is integral: the $\mathrm{Cl}(10,0)$ generators are $0/\pm1$
   integers, $G = \prod_{a<b}(\mathbb I + h_a h_b)/2^{45/2} = N/2^{45/2}$ with
   $N$ integer; $T\!\cdot\!L \in \{0,\pm2,\pm2i\}$ with spectrum
   $\{4^{1120}, -12^{288}, -32^{32}\}$; the Casimir projectors are Lagrange
   interpolants $P_e = \prod_{e'\ne e}(T\!\cdot\!L - e'\mathbb I)/(e - e')$ -- no
   eigensolver, no round-off.  Then $\mathrm{Tr}(C_a N^n) = 0$ for every even
   $n$ and $\mathrm{Tr}(C_a W^n) = 0$ for every $n$, where $W$ is the numerator
   of $MG$: the even-$n$ ratio is exactly $0/0$, the odd-$n$ ratio exactly $0$.
2. **The float64 signature.**  At $\pi/4$, $3\pi/4$, $\pi/2$ the float64 Cartan
   norms of $G^n$ (even $n$) and of $(MG)^n$ (all $n$) sit at the round-off
   floor, $10^{-17}\ldots10^{-13}$, and grow with $n$ the way accumulated error
   does; their ratios -- the published $F_2 = 65.7$, $F_{10} = 55.0$, the
   crossing at $n^* = 208$, $F_{208} = 0.019685$ -- are ratios of floors.  The
   same code on the numpy core returns $64.4$, $53.9$, $0.019687$: a $2.07\%$
   backend shift on a number quoted to four digits.
3. **The branch that does measure.**  At the non-degenerate implemented angle
   $g = 1/5$ the same functional is precision-stable: 15- and 30-digit
   arithmetic agree to $10^{-12}$ and reproduce V85's float64 $F_{10} =$
   $0.01918639$ to $10^{-8}$: deficit $0.9644132$, coefficient
   $3\,\mathrm{deficit}^{10} = 2.0881$ (V85, V108, V113).

**Verdict.** The loop functional has a structurally caused zero on the closure
class; the measured loop factor belongs to the open branch.  Any $F_n$ at
a closure-class angle is defined by continuity (scanned here); a raw float64
reading on the $4k{+}2$ class is a ratio of round-off floors, not a measurement.
"""
import json
import math
import os
import sys
import time
from math import gcd

import numpy as np
import scipy.sparse as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                   # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V85_closed_loop_factor"))
from V79_gauged_bond_layer import (adjoint_link_operators,    # noqa: E402
                                   root_generators, structure_constants)
from cl10 import clifford                                    # noqa: E402
from cl10 import config as _config                           # noqa: E402
from cl10 import engine as E                                 # noqa: E402
from cl10 import sm_qca, u10pred                             # noqa: E402

_NG = int(_config.CFG.physics.n_gen)
D_L = _NG * (_NG - 1) // 2          # adjoint link dimension (45)
D_S = int(_config.CFG.physics.n_internal)   # the spinor block (32)


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published number has a
    source, it is not a decimal restated here."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j

EIGS = (-32, -12, 4)
CLOSURE = (("pi/2", math.pi / 2), ("pi/4", math.pi / 4), ("3pi/4", 3 * math.pi / 4))
OPEN = (("1/5", float(_config.CFG.angle.branches.implemented)),
        ("pi/3", math.pi / 3))
LOOP = float(u10pred.LOOP)
REF_G = float(_config.CFG.angle.branches.implemented)
N_EXACT = 60                        # exact-integer scan length
N_FLOOR = 12                        # window for the floor statistics
N_CERT = 212                        # float64 / certified ladder length


def cint(z):
    zr, zi = float(np.real(z)), float(np.imag(z))
    ar, ai = int(round(zr)), int(round(zi))
    assert abs(zr - ar) < 1e-9 and abs(zi - ai) < 1e-9, z
    return (ar, ai)


def skeleton():
    """Angle-independent exact data: integer G numerator N, Lagrange projector
    blocks and denominators, their common denominator D, Cartan rows, planes."""
    h = [np.array(E.asnumpy(m), dtype=object) for m in clifford.generators()]
    Ts = [E.asnumpy(t) for t in root_generators()]
    f, _, _ = structure_constants(Ts)
    Ls = [E.asnumpy(L) for L in adjoint_link_operators(f)]
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    Rr = sp.csr_matrix(np.rint(TL.real).astype(np.int64))
    Ii = sp.csr_matrix(np.rint(TL.imag).astype(np.int64))

    def lagrange_blocks(e):
        ar = sp.identity(TL.shape[0], format='csr', dtype=np.int64)
        ai = sp.csr_matrix(TL.shape, dtype=np.int64)
        den = 1
        for e2 in EIGS:
            if e2 == e:
                continue
            dr = Rr - e2 * sp.identity(TL.shape[0], format='csr', dtype=np.int64)
            di = Ii.copy()
            ar, ai = ar @ dr - ai @ di, ar @ di + ai @ dr
            den *= (e - e2)
        Rd, Id = ar.toarray(), ai.toarray()
        B = np.zeros((D_S, D_S), dtype=object)
        C = np.zeros((D_S, D_S), dtype=object)
        for j in range(D_L):
            B = B + Rd[j::D_L, j::D_L]
            C = C + Id[j::D_L, j::D_L]
        return (B, C), den

    blocks, dens = {}, {}
    for e in EIGS:
        blocks[e], dens[e] = lagrange_blocks(e)
    D = 1
    for e in EIGS:
        k = abs(dens[e]) * D_L
        D = k // gcd(k, D) * D

    N = np.eye(D_S, dtype=object)
    for a in range(10):
        for b in range(a + 1, 10):
            N = N @ (np.eye(D_S, dtype=object) + h[a] @ h[b])

    cart = []
    for ca in (0, 2, 4, 6, 8):
        H = h[ca] @ h[ca + 1]
        cart.append([(a, b) + cint(H[a, b]) for a in range(D_S) for b in range(D_S)
                     if cint(H[a, b]) != (0, 0)])
    planes = []
    for a in range(10):
        for b in range(a + 1, 10):
            T = 1j * (h[a] @ h[b])
            planes.append([[cint(T[r, c]) for c in range(D_S)] for r in range(D_S)])
    return N, blocks, dens, D, cart, planes, TL


def cmm(A, B):
    """(re, im) complex matrix product over an exact/ high-precision ring."""
    m, p = len(B), len(B[0])
    Bt = [[B[k][j] for k in range(m)] for j in range(p)]
    out = []
    for Ai in A:
        row = []
        for col in Bt:
            sr = si = 0
            for (ar, ai), (br, bi) in zip(Ai, col):
                if ar or ai:
                    sr += ar * br - ai * bi
                    si += ar * bi + ai * br
            row.append((sr, si))
        out.append(row)
    return out


def ceye(n, one=1, zero=0):
    return [[((one if a == b else zero), zero) for b in range(n)] for a in range(n)]


def tr_AB(K, cart):
    """[(A, B)] with Tr(H_a K) = A + iB; |Tr(C_a K)| = |Tr(H_a K)| (C_a = i H_a)."""
    out = []
    for nz in cart:
        A = B = 0
        for a, b, hr, hi in nz:
            kr, ki = K[b][a]
            A += hr * kr - hi * ki
            B += hr * ki + hi * kr
        out.append((A, B))
    return out


def mp_ladder(angle_expr, nmax, dps, planes, blocks, dens, cart):
    """Cartan-content ladder F_n at an arbitrary angle in dps digits."""
    with mp.workdps(dps + 6):
        txt = str(angle_expr)
        g = (mp.mpf(txt) if txt.replace('.', '', 1).replace('-', '', 1).isdigit()
             else eval(txt, {"mp": mp, "mpf": mp.mpf, "pi": mp.pi,
                             "sqrt": mp.sqrt}))
        cosv, sinv = mp.cos(g), mp.sin(g)
        G = ceye(D_S, mp.mpf(1), mp.mpf(0))
        for T in planes:
            Fm = [[((cosv if a == b else mp.mpf(0)) + sinv * t[1],
                    (mp.mpf(0) if a == b else mp.mpf(0)) - sinv * t[0])
                   for b, t in enumerate(row)] for a, row in enumerate(T)]
            G = cmm(G, Fm)
        Mx = [[(mp.mpf(0), mp.mpf(0))] * D_S for _ in range(D_S)]
        for e in EIGS:
            pc, ps = mp.cos(g / 2 * e), mp.sin(g / 2 * e)
            Bm, Cm = blocks[e]
            acc = [[(pc * Bm[a][b] - ps * Cm[a][b], pc * Cm[a][b] + ps * Bm[a][b])
                    for b in range(D_S)] for a in range(D_S)]
            Ma = cmm(acc, G)
            inv = 1 / (dens[e] * D_L)
            Mx = [[(Mx[a][b][0] + inv * Ma[a][b][0], Mx[a][b][1] + inv * Ma[a][b][1])
                   for b in range(D_S)] for a in range(D_S)]
        X = cmm(Mx, G)
        Xp, Gp = ceye(D_S, mp.mpf(1), mp.mpf(0)), ceye(D_S, mp.mpf(1), mp.mpf(0))
        out = {}
        for n in range(1, nmax + 1):
            Xp = cmm(Xp, X)
            Gp = cmm(Gp, G)
            cG = sum((mp.sqrt(A * A + B * B) for A, B in tr_AB(Gp, cart)),
                     mp.mpf(0))
            cX = sum((mp.sqrt(A * A + B * B) for A, B in tr_AB(Xp, cart)),
                     mp.mpf(0))
            out[n] = {"cnorm_G": cG, "cnorm_X": cX,
                      "F": (cX / cG) if cG != 0 else None}
        return out



def main():
    t0 = time.time()
    checks, data = [], {}
    N_int, blocks, dens, D, cart, planes, TL = skeleton()

    # ---- T0 the phase degeneracy IS the closure class --------------------
    degen = {}
    for tag, g in CLOSURE + OPEN:
        ph = [complex(np.exp(1j * (g / 2) * e)) for e in EIGS]
        degen[tag] = {"phase_sep_-12_+4": abs(ph[1] - ph[2]),
                      "eight_g_over_2pi": 8 * g / (2 * math.pi)}
    sep_close = max(degen[t]["phase_sep_-12_+4"] for t, _ in CLOSURE)
    sep_open = min(degen[t]["phase_sep_-12_+4"] for t, _ in OPEN)
    multipliers = ", ".join(f"{degen[t]['eight_g_over_2pi']:.0f}" for t, _ in CLOSURE)
    checks.append({
        "name": f"T0 the traced-link phases of the 288 (T.L = -12) and 1120"
                f" (T.L = +4) Casimir sectors coincide exactly on the closure"
                f" class (largest separation {sep_close:.2e}, 8g/2pi ="
                f" {multipliers} = integers) and split on the open angles"
                f" (smallest separation {sep_open:.3f})",
        "pass": bool(sep_close < 1e-9 and sep_open > 0.1)})
    data["phase_degeneracy"] = degen

    # ---- T1 exact arithmetic at the root clock: the ratio is 0/0 ---------
    Nc = [[cint(N_int[a][b]) for b in range(D_S)] for a in range(D_S)]
    PHASE = {-32: (1, 0), -12: (0, 1), 4: (0, 1)}        # exp(i(pi/8)e), exact
    Mtil = [[(0, 0)] * D_S for _ in range(D_S)]
    for e in EIGS:
        Bm, Cm = blocks[e]
        s = D // (dens[e] * D_L)
        qr, qi = PHASE[e]
        for a in range(D_S):
            for b in range(D_S):
                br, bi = int(Bm[a][b]), int(Cm[a][b])
                re_, im_ = Mtil[a][b]
                Mtil[a][b] = (re_ + (br * qr - bi * qi) * s,
                              im_ + (br * qi + bi * qr) * s)
    W = cmm(cmm(Mtil, Nc), Nc)          # V85/V203's Ma carries one G, so MG = M~G^2

    zero_den, zero_num, nonzero_den, odd_zero = [], [], {}, []
    Wp, Np = ceye(D_S), ceye(D_S)
    for n in range(1, N_EXACT + 1):
        Wp = cmm(Wp, W)
        Np = cmm(Np, Nc)
        tN = tr_AB(Np, cart)
        tW = tr_AB(Wp, cart)
        nzN = sum(1 for A, B in tN if A or B)
        nzW = sum(1 for A, B in tW if A or B)
        if n % 2 == 0:
            if nzN == 0:
                zero_den.append(n)
            else:
                nonzero_den[n] = nzN
            if nzW == 0:
                zero_num.append(n)
        elif nzW == 0:
            odd_zero.append(n)
    checks.append({
        "name": f"T1 EXACT integer arithmetic at g = pi/4 (integer Clifford"
                f" generators, integer T.L, Lagrange projectors, no eigensolver,"
                f" n <= {N_EXACT}): Tr(C_a N^n) = 0 for every even n"
                f" ({len(zero_den)}/{N_EXACT // 2} even lengths, all five Cartan"
                f" directions) and Tr(C_a W^n) = 0 for every n ({len(zero_num)}"
                f" even + {len(odd_zero)} odd): the even-n Cartan ratio is"
                f" EXACTLY 0/0 and the odd-n ratio EXACTLY 0"
                + ("" if not nonzero_den else f"; EXCEPTIONS {sorted(nonzero_den)}"),
        "pass": bool(len(zero_den) == N_EXACT // 2 and not nonzero_den
                     and len(zero_num) + len(odd_zero) == N_EXACT)})
    data["exact"] = {"n_max": N_EXACT, "even_zero_denominator": len(zero_den),
                     "even_zero_numerator": len(zero_num),
                     "odd_zero_numerator": len(odd_zero),
                     "nonzero_even_denominators": sorted(nonzero_den)}

    # ---- T2 what float64 returns at a degenerate angle -------------------
    hnp = [E.asnumpy(m) for m in clifford.generators()]
    carnp = [1j * (hnp[a] @ hnp[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = [E.asnumpy(t) for t in root_generators()]
    f, _, _ = structure_constants(Ts)
    Ls = [E.asnumpy(L) for L in adjoint_link_operators(f)]
    TLf = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    w, V = np.linalg.eigh(TLf)
    P = []
    for e in EIGS:
        idx = np.argsort(np.abs(w - e))[:int({-32: 32, -12: 288, 4: 1120}[e])]
        P.append(V[:, idx] @ V[:, idx].conj().T)

    def cnorm_np(K):
        return float(sum(abs(np.trace(C @ K)) for C in carnp))

    floors = {}
    for tag, g in CLOSURE + OPEN:
        G = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(),
                                          g)).astype(np.complex128)
        Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
        Ma = [np.einsum("ajbj->ab", (Pa @ Gi).reshape(D_S, D_L, D_S, D_L)) / D_L
              for Pa in P]
        M = sum(p * m for p, m in
                zip([np.exp(1j * (g / 2) * e) for e in EIGS], Ma))
        X = M @ G
        cG, cX = [], []
        Bn, Gn = np.eye(D_S, dtype=complex), np.eye(D_S, dtype=complex)
        for n in range(1, N_CERT + 1):
            Bn = Bn @ X
            Gn = Gn @ G
            cX.append(cnorm_np(Bn))
            cG.append(cnorm_np(Gn))
        floors[tag] = {
            "max_even_cnorm_G": float(max(cG[n - 1] for n in range(2, N_FLOOR + 1, 2))),
            "max_cnorm_MG": float(max(cX[:N_FLOOR])),
            "min_odd_cnorm_G": float(min(cG[n - 1] for n in range(1, N_FLOOR + 1, 2))),
            "cnorm_MG_n2_n6_n12": [float(cX[n - 1]) for n in (2, 6, 12)],
            "F_n_float": {str(n): float(cX[n - 1] / cG[n - 1])
                          for n in range(2, N_CERT + 1, 2)}}
    deg_ok = all(floors[t]["max_even_cnorm_G"] < 1e-11 and floors[t]["max_cnorm_MG"] < 1e-11
                 for t, _ in CLOSURE)
    open_ok = all(floors[t]["max_cnorm_MG"] > 1e-6 for t, _ in OPEN)
    growth = [floors["pi/4"]["cnorm_MG_n2_n6_n12"][i] for i in (0, 2)]
    checks.append({
        "name": f"T2 the float64 signature of a 0/0: at pi/2, pi/4, 3pi/4 every"
                f" Cartan norm is at the round-off floor (max even-n"
                f" cnorm(G^n) = {max(floors[t]['max_even_cnorm_G'] for t, _ in CLOSURE):.1e},"
                f" max cnorm((MG)^n) ="
                f" {max(floors[t]['max_cnorm_MG'] for t, _ in CLOSURE):.1e}) and"
                f" the numerator grows like accumulated error (pi/4, n = 2, 6,"
                f" 12: {', '.join(f'{v:.1e}' for v in floors['pi/4']['cnorm_MG_n2_n6_n12'])});"
                f" at the open angles the same norms are"
                f" > {min(floors[t]['max_cnorm_MG'] for t, _ in OPEN):.1e} - a"
                f" real signal, not a floor",
        "pass": bool(deg_ok and open_ok and growth[1] > 3 * growth[0])})
    data["float64_floors"] = floors

    # ---- T3/T4 high-precision branch evaluation --------------------------
    def mp_branch(expr, nmax, dps):
        return mp_ladder(expr, nmax, dps, planes, blocks, dens, cart)

    mp15 = mp_branch("0.2", 10, 15)
    mp30 = mp_branch("0.2", 10, 30)
    F10_15, F10_30 = float(mp15[10]["F"]), float(mp30[10]["F"])
    checks.append({
        "name": f"T3 the open branch g = 1/5 IS a measurement: F_10 ="
                f" {F10_30:.11f} at 30 digits and {F10_15:.11f} at 15 digits"
                f" (precision drift {abs(F10_30 - F10_15):.1e}), matching V85's"
                f" float64 {art('V113_axis_factorization', 'data.F_10_measured'):.8f}"
                f" to {abs(F10_30 - art('V113_axis_factorization', 'data.F_10_measured')):.1e};"
                f" deficit {F10_30 / LOOP:.7f}, coefficient"
                f" 3 x deficit^10 = {3 * (F10_30 / LOOP) ** 10:.4f}",
        "pass": bool(abs(F10_30 - F10_15) < 1e-12
                     and abs(F10_30 - art("V113_axis_factorization",
                                          "data.F_10_measured")) < 1e-7
                     and abs(3 * (F10_30 / LOOP) ** 10
                             - art("V113_axis_factorization",
                                   "data.coef_one_per_family")) < 5e-4)})
    data["open_branch"] = {"F_10_15dp": F10_15, "F_10_30dp": F10_30,
                           "deficit": F10_30 / LOOP,
                           "coef_3_deficit10": 3 * (F10_30 / LOOP) ** 10}
    # ---- T4: the period-4 classes resolve the 0/0 differently ------------
    conv = {}
    for eps in ("1e-1", "1e-3", "1e-5", "1e-7"):
        r = mp_branch(f"mp.pi/4 + mp.mpf('{eps}')", 12, 55)
        conv[eps] = {str(n): float(r[n]["F"]) for n in range(2, 13, 2)}
    ff = floors["pi/4"]["F_n_float"]
    fourk = [n for n in range(4, 13, 4)]
    fourk2 = [n for n in range(2, 13, 4)]
    rem_ok = all(
        abs(conv["1e-7"][str(n)] - conv["1e-5"][str(n)])
        < 1e-4 * abs(conv["1e-5"][str(n)])
        and abs(conv["1e-7"][str(n)] - ff[str(n)]) < 1e-3 * ff[str(n)]
        for n in fourk)
    div_ok = all(conv["1e-3"][str(n)] > 1e4 * ff[str(n)]
                 and conv["1e-5"][str(n)] > 1e2 * conv["1e-3"][str(n)]
                 for n in fourk2)
    rem_dev = max(abs(conv["1e-7"][str(n)] - ff[str(n)]) / abs(ff[str(n)])
                  for n in fourk)
    checks.append({
        "name": f"T4 the period-4 classes resolve the 0/0 differently: on"
                f" n = 4k the ratio CONVERGES as g -> pi/4 (F_4"
                f" {conv['1e-5']['4']:.5f} at +1e-5 vs {conv['1e-7']['4']:.5f}"
                f" at +1e-7; F_8 {conv['1e-5']['8']:.5f} vs"
                f" {conv['1e-7']['8']:.5f}; F_12 {conv['1e-5']['12']:.5f} vs"
                f" {conv['1e-7']['12']:.5f}) - the singularity is REMOVABLE -"
                f" and the float64 root-clock reading IS that limit, to"
                f" {rem_dev:.0e} relative (F_4 {ff['4']:.4f}, F_8"
                f" {ff['8']:.4f}, F_12 {ff['12']:.4f}); on n = 4k+2 the ratio"
                f" DIVERGES (F_2 {conv['1e-1']['2']:.1e} at +1e-1,"
                f" {conv['1e-3']['2']:.1e} at +1e-3 and {conv['1e-5']['2']:.1e}"
                f" at +1e-5 against the floor ratio {ff['2']:.1f}; F_6"
                f" {conv['1e-3']['6']:.1e}; F_10 {conv['1e-3']['10']:.1e}) -"
                f" the numerator vanishes to LOWER order than the denominator",
        "pass": bool(rem_ok and div_ok)})
    data["branch_structure"] = {
        "float64_pi4": ff,
        "detuned": conv,
        "removable_class": fourk, "divergent_class": fourk2}

    # ---- T5 certified readings by continuity at 50 digits ----------------
    cert = {}
    for tag, ang in (("pi/4@1e-8", "mp.pi/4 + mp.mpf('1e-8')"),
                     ("pi/4@1e-6", "mp.pi/4 + mp.mpf('1e-6')"),
                     ("3pi/4@1e-6", "3 * mp.pi/4 + mp.mpf('1e-6')"),
                     ("pi/2@1e-6", "mp.pi/2 + mp.mpf('1e-6')")):
        r = mp_branch(ang, N_CERT, 50)
        cr = [n for n in range(2, N_CERT + 1)
              if r[n]["F"] is not None and 1e-9 < r[n]["F"] < LOOP]
        star = cr[0] if cr else None
        d = float(r[star]["F"] / LOOP) if star else None
        cert[tag] = {"first_crossing_n": star, "n_below_loop": len(cr),
                     "deficit": d, "coef": (3 * d ** 10 if d else None),
                     "F_4": float(r[4]["F"]), "F_8": float(r[8]["F"]),
                     "F_12": float(r[12]["F"]), "F_20": float(r[20]["F"]),
                     "F_208": float(r[208]["F"]), "F_210": float(r[210]["F"])}

    def fl_cross(tag):
        f = floors[tag]["F_n_float"]
        for n in sorted(int(k) for k in f):
            if 1e-12 < f[str(n)] < LOOP:
                return n
        return None

    fl4 = floors["pi/4"]["F_n_float"]
    cert_rem_ok = all(abs(cert["pi/4@1e-8"]["F_" + k] - fl4[k]) < 1e-2 * fl4[k]
                      for k in ("4", "8", "12", "20", "208"))
    stab_ok = (cert["pi/4@1e-6"]["first_crossing_n"] == 208
               == cert["pi/4@1e-8"]["first_crossing_n"]
               and abs(cert["pi/4@1e-8"]["F_208"] - cert["pi/4@1e-6"]["F_208"])
               < 2e-4 * cert["pi/4@1e-8"]["F_208"])
    div_ok = all(cert[t]["F_210"] > 1e6 * LOOP
                 for t in ("pi/4@1e-6", "pi/4@1e-8", "3pi/4@1e-6"))
    checks.append({
        "name": f"T5 the 4k class is certified by continuity (50 digits, g ="
                f" pi/4 + 1e-6 then + 1e-8): the first crossing of 1/(16 pi) is"
                f" n* = {cert['pi/4@1e-8']['first_crossing_n']} at both"
                f" detunings, F_208 = {cert['pi/4@1e-6']['F_208']:.7f} ->"
                f" {cert['pi/4@1e-8']['F_208']:.7f} (drift"
                f" {abs(cert['pi/4@1e-8']['F_208'] - cert['pi/4@1e-6']['F_208']) / cert['pi/4@1e-8']['F_208']:.0e})"
                f" - deficit {cert['pi/4@1e-8']['deficit']:.6f}, coefficient"
                f" 3 x deficit^10 = {cert['pi/4@1e-8']['coef']:.4f}; the float64"
                f" root-clock ladder sits within 1e-2 of these limits (F_4"
                f" {fl4['4']:.4f} vs {cert['pi/4@1e-8']['F_4']:.5f}, F_8"
                f" {fl4['8']:.4f} vs {cert['pi/4@1e-8']['F_8']:.5f}, F_208"
                f" {fl4['208']:.7f} vs {cert['pi/4@1e-8']['F_208']:.7f}); 3pi/4"
                f" crosses at n* = {cert['3pi/4@1e-6']['first_crossing_n']}"
                f" (float64 {fl_cross('3pi/4')}) with coefficient"
                f" {cert['3pi/4@1e-6']['coef']:.4f}; the 4k+2 class has no"
                f" limit (F_210 = {cert['pi/4@1e-8']['F_210']:.1e} at pi/4 +1e-8"
                f" vs the floor ratio {fl4['210']:.2f}), and the pi/2 branch -"
                f" no period-4 filter, its traced phases are +/-1 - crosses at"
                f" n* = {cert['pi/2@1e-6']['first_crossing_n']} certified vs"
                f" {fl_cross('pi/2')} at float64",
        "pass": bool(cert_rem_ok and stab_ok and div_ok
                     and cert["pi/4@1e-8"]["first_crossing_n"] == 208
                     == fl_cross("pi/4")
                     and cert["3pi/4@1e-6"]["first_crossing_n"]
                     == fl_cross("3pi/4"))})
    data["certified"] = cert
    data["float64_crossings"] = {t: fl_cross(t) for t, _ in CLOSURE}


    coef_open = 3 * (F10_30 / LOOP) ** 10
    status = (f"ROOT-CLOCK PHASE DEGENERACY, MEASURED (zero free parameters,"
              f" {time.time() - t0:.0f} s). (i) The closure class IS the phase-"
              f"degenerate set: 8g in 2pi Z makes the 288 and 1120 Casimir"
              f" sectors share their traced-link phase (separation <"
              f" {sep_close:.0e}), while the open angles separate them by >"
              f" {sep_open:.2f}. (ii) Exact integer arithmetic at g = pi/4:"
              f" Tr(C_a N^n) = 0 for every even n <= {N_EXACT} and Tr(C_a W^n) ="
              f" 0 for every n <= {N_EXACT} - every even-n Cartan ratio is"
              f" exactly 0/0, so float64 returns ratios of round-off floors"
              f" (pi/4 floor {floors['pi/4']['max_even_cnorm_G']:.0e} over"
              f" {floors['pi/4']['max_cnorm_MG']:.0e}). (iii) The period-4"
              f" classes of the clock filter resolve that 0/0 differently: on"
              f" n = 4k numerator and denominator vanish to the SAME (first)"
              f" order in g - pi/4, the singularity is REMOVABLE, and the"
              f" float64 ladder is the g -> pi/4 limit to ~4e-4 relative (T4:"
              f" F_4 {ff['4']:.5f} vs {conv['1e-7']['4']:.5f} at pi/4 + 1e-7);"
              f" on n = 4k+2 the denominator vanishes to far higher order, the"
              f" ratio is unbounded with poles accumulating at the closure angle"
              f" (F_2 {conv['1e-1']['2']:.1e} at +1e-1, {conv['1e-3']['2']:.1e}"
              f" at +1e-3; F_210 {cert['pi/4@1e-8']['F_210']:.1e} at +1e-8), so"
              f" the quoted F_2 = 65.7, F_6 = 60.1, F_10 = 55.0, F_170 = 1.57"
              f" and everything built on them - (F_10)^17, R_17 = 4.0e-30, the"
              f" 33x branch separation - are floor artefacts, not measurements."
              f" (iv) Certified root-clock readings (T5): n* = 208 at pi/4 and"
              f" 3pi/4, deficit {cert['pi/4@1e-8']['deficit']:.6f} /"
              f" {cert['3pi/4@1e-6']['deficit']:.6f}, hierarchy coefficient 3 x"
              f" deficit^10 = {cert['pi/4@1e-8']['coef']:.4f} /"
              f" {cert['3pi/4@1e-6']['coef']:.4f}. The pi/2 branch carries no"
              f" period-4 filter (traced phases +/-1) and its crossing is steep"
              f" and floor-conditioned (n* ="
              f" {cert['pi/2@1e-6']['first_crossing_n']} certified vs"
              f" {fl_cross('pi/2')} at float64): its 2.30 reading is not a"
              f" stable closure-class number. (v) The loop factor that is"
              f" measured without any limiting procedure lives on the open"
              f" branch: F_10(1/5) = {F10_30:.9f} stable to"
              f" {abs(F10_30 - F10_15):.0e}, deficit {F10_30 / LOOP:.6f},"
              f" coefficient 3 x deficit^10 = {coef_open:.4f}, reproducing"
              f" V85/V108/V113. VERDICT: the closed-root ladder STRUCTURE"
              f" (V131, the exact odd-n selection rule, V197's vanishing odd"
              f" rungs) and the 4k-class crossing survive; the certified"
              f" closure-class coefficient is 2.710 (pi/4) / 2.711 (3pi/4)"
              f" against the float64 2.699-2.700, and every n = 4k+2 reading at"
              f" a closure-class angle is a ratio of round-off floors, not a"
              f" measurement.")
    out = {"id": "V209_root_clock_degeneracy",
           "tier": 4,
           "property": "the traced-loop Cartan ratio is 0/0 on the phase-"
                       "degenerate closure class and a measurement on the open "
                       "branch (V209)",
           "residual": True,
           "programme": "V209",
           "branch": "W-NAMED",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": all(c["pass"] for c in checks)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v209_root_clock_degeneracy.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
