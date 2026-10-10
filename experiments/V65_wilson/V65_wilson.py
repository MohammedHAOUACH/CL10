#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V65 - the Wilson mechanism: does the QCA REMOVE the half-zone nodes (W9), and
does the removal PREDICT the desert spectrum (W4)?

**The question.**  V63 proved the doubling is present: the 16-sector one-step
matrix is gapless at the 12 half-zone momenta (+-pi/2, +-pi/2, 0).  V64 proved the measured
couplings REQUIRE vectorlike matter at two thresholds (M_10 = 1.2e7,
M_(1,2) = 1.0e12 GeV).  If the SAME mechanism - a Wilson-type term of the
step - both gaps the nodes and fixes their masses, then W9's pathology
becomes W4's prediction.

**Part A - the mechanism (W9).**  The Cartan of Spin(10) in the root:
Sigma_k = i h_{2k-1} h_{2k}, commuting Hermitian involutions; in the 16 basis
they are simultaneously diagonal (the 16 weights, +-1 sign vectors of odd
parity).  The hypercharge Y and the SU(5)-preserving X are Cartan
combinations - IDENTIFIED here by matching eigenvalue multisets to the
0-parameter SM charges (both the Georgi-Glashow and the flipped SU(5) are
found).  The Wilson factor per step,
    W(q) = exp(-i w (1-cos qx + 1-cos qy + 1-cos qz) (a X + b Y)),
is internal (the 16 sector stays closed), exactly unitary, and W(0) = I: the
zone-centre gap is untouched.  At the 12 M momenta sum(1-cos q_j) = 2.
The generator A = w sum_j (1 - cos q_j) D is strictly range 1, but exp(-iA)
is only QUASI-local (Bessel tails): the strict light cone of Theorem 5
becomes a Lieb-Robinson cone - MEASURED here on the real-space step.  The
amplitude in shells of one block beyond the bond cone decays exponentially
(worst ratio 0.25 per block at w = 0.1): 5.6e-06 in the first shell, 3.1e-09
at the Lieb-Robinson line r = 1.2 n + 6 after 24 stages.

**Part B - the thresholds (W4).**  The node mass M_f = w |a X_f + b Y_f|
M_P / 3 (the /3: one step is 6 Planck times, the node gap is
2w|D_f| step-phase).  Uniform w splits ALL 12 nodes identically: 12 copies
of the full split 16 -> Delta b = (32, 32, 32) against the required ~
(2.4, 4, 3.9): OVERPRODUCTION, reported not hidden.  The selective reading
(3 of 12 nodes split, 9 stay at the cutoff) is SOLVED against the meeting
equations with alpha_10 = 1/25 fixed: one angle theta = (a : b) and one
scale w predict six fragment thresholds; the scan finds whether a physical
solution exists and what it predicts.
"""

import itertools
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import clifford
from cl10 import engine as E
from cl10 import supercell as SC
from cl10 import unification as UN

M_P, M_Z, A_Z, B_SM = UN.M_P, UN.M_Z, UN.A_Z, UN.B_SM
TWO_PI = 2.0 * math.pi
INV_A10 = 25.0
Z_MIN = math.log(1e6 / M_Z)
X_MAX = math.log(M_P / M_Z)
G10 = SC.G10


# ---------------------------------------------------------------------------
# the Cartan of the root and the 16 weights
# ---------------------------------------------------------------------------
def cartan_weights():
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    B, _U32, _U16 = SC.sectors()
    Dk = [B.conj().T @ C @ B for C in cartan]
    # the five commute; algebraically independent coefficients -> no
    # degeneracy of the combination, so eigh gives a COMMON eigenbasis
    M = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(M)
    lam = np.array([[float(np.real(np.vdot(V[:, s], Dk[k] @ V[:, s])))
                     for k in range(5)] for s in range(16)])
    off = max(float(np.abs(V.conj().T @ Dk[k] @ V - np.diagflat(
        [np.vdot(V[:, s], Dk[k] @ V[:, s]) for s in range(16)])).max())
        for k in range(5))
    return lam, off


def identify_XY(lam):
    """All Cartan combinations X whose eigenvalue multiset is {1x10, -3x5,
    5x1} (both SU(5) embeddings), and Y with the SM hypercharge multiset."""
    y_target = np.sort(np.array([1 / 6] * 6 + [-2 / 3] * 3 + [1 / 3] * 3
                                 + [-1 / 2] * 2 + [1.0, 0.0]))
    x_target = np.sort(np.array([1.0] * 10 + [-3.0] * 5 + [5.0]))

    def hits(target, grid, chunk=200_000):
        out = []
        it = itertools.product(grid, repeat=5)
        while True:
            C = np.array(list(itertools.islice(it, chunk)))
            if len(C) == 0:
                return out
            v = np.sort(C @ lam.T, axis=1)
            for h in np.where(np.all(np.abs(v - target) < 1e-9, axis=1))[0]:
                out.append(C[h])

    xs = hits(x_target, [float(k) for k in range(-4, 7)])
    # dedupe up to overall sign
    ux = []
    for c in xs:
        if not any(np.allclose(c, s * u) for u in ux for s in (1, -1)):
            ux.append(c)
    ys = hits(y_target, [k / 12.0 for k in range(-12, 13)])[:1]
    return ux, ys


# ---------------------------------------------------------------------------
# Part A - the Wilson factor on the one-step matrix
# ---------------------------------------------------------------------------
def step_wilson(q, U16, w, ddiag):
    M = SC.bloch(q, U16)
    ph = w * sum(1.0 - math.cos(qj) for qj in q)
    W = np.diag(np.exp(-1j * ph * ddiag))
    return np.kron(np.eye(8, dtype=complex), W) @ M


def gap_of(M):
    return float(np.min(np.abs(np.sort(np.angle(np.linalg.eigvals(M))))))


def wilson_kernel_decay(w, cmax, L=32):
    """Position-space kernel of exp(-i w c sum_j (1-cos q_j)) for the largest
    charge cmax: Bessel tails.  Returns (tail at r>=2 / centre, tail r>=4)."""
    qs = 2.0 * np.pi * np.arange(L) / L
    Wq = np.exp(-1j * w * cmax * (3.0 - (np.cos(qs)[:, None, None]
                                        + np.cos(qs)[None, :, None]
                                        + np.cos(qs)[None, None, :])))
    k = np.fft.ifftn(Wq)
    rad = (np.minimum(np.arange(L), L - np.arange(L))[:, None, None]
           + np.minimum(np.arange(L), L - np.arange(L))[None, :, None]
           + np.minimum(np.arange(L), L - np.arange(L))[None, None, :])
    a = np.abs(k)
    return float(a[rad >= 2].max() / a[0, 0, 0]), float(a[rad >= 4].max()
                                                        / a[0, 0, 0])


def bessel_J(w, kmax=8, n=200001):
    """J_k(w) by the standard integral J_k = (1/pi) int_0^pi cos(k t - w sin t) dt.

    The block-lattice kernel of exp(-i w c (1 - cos q)) is e^{-iwc} (-i)^m
    J_m(wc) in the block displacement m, so |J_m(wc)| is the predicted tail."""
    t = np.linspace(0.0, np.pi, n)
    st = np.sin(t)
    out = []
    for k in range(kmax + 1):
        f = np.cos(k * t - w * st)
        out.append(float(np.trapezoid(f, t) / np.pi))
    return out



def light_cone_wilson(Lb, n_stages, w, ddiag, U16, tol=1e-10, seed=0):
    """Real-space 16-sector evolution on the block lattice: the six bond
    stages of the step plus the Wilson factor once per step.  The bond
    stages are strictly causal (one Planck cell per stage, Theorem 5); the
    Wilson factor is a block-momentum multiplication, so it is only
    QUASI-local: its kernel in block displacement m is |J_m(w c)|.  The
    test measures the amplitude in shells of block distance beyond the
    bond cone and compares it with that Bessel prediction."""
    rng = np.random.default_rng(seed)
    f = np.zeros((Lb, Lb, Lb, SC.NSUB, 16), dtype=complex)
    v = rng.standard_normal(16) + 1j * rng.standard_normal(16)
    x0 = (Lb // 2, Lb // 2, Lb // 2)
    f[x0 + (0,)] = v / np.linalg.norm(v)
    qs = 2.0 * np.pi * np.arange(Lb) / Lb
    s1 = 3.0 - (np.cos(qs)[:, None, None] + np.cos(qs)[None, :, None]
                + np.cos(qs)[None, None, :])
    Wfac = np.exp(-1j * w * s1[..., None, None]
                  * ddiag[None, None, None, None, :])
    R = np.arange(Lb)
    sub = np.arange(SC.NSUB)
    ex, ey, ez = sub & 1, (sub >> 1) & 1, (sub >> 2) & 1
    dd = np.abs(R - x0[0])
    dd = np.minimum(dd, Lb - dd)
    dist = ((2 * dd[:, None, None, None] + ex)
            + (2 * dd[None, :, None, None] + ey)
            + (2 * dd[None, None, :, None] + ez))
    lr = 1.2 * n_stages + 6.0   # the reported reference line (1 cell per
                                # stage + one Wilson cell per step)
    radii = []
    outside = 0.0
    shells = [0.0] * 7          # max amplitude k blocks beyond the bond cone
    for t in range(n_stages):
        f = SC.stage(f, U16, t)
        if t % 6 == 5:
            f = np.fft.ifftn(np.fft.fftn(f, axes=(0, 1, 2)) * Wfac,
                             axes=(0, 1, 2))
        nrm = np.sqrt(np.sum(np.abs(f) ** 2, axis=-1))
        if t == n_stages - 1:
            outside = float(nrm[dist > lr].max())
            cone = float(n_stages)      # the bond light cone, Planck cells
            for k in range(len(shells)):
                sel = (dist > cone + 2 * k) & (dist <= cone + 2 * (k + 1))
                if sel.any():
                    shells[k] = float(nrm[sel].max())
        idx = np.argwhere(nrm > tol)
        rad = max(int(dist[tuple(p)]) for p in idx) if len(idx) else 0
        radii.append(rad)
    return radii, outside, shells


# ---------------------------------------------------------------------------
# Part B - thresholds vs the meeting equations
# ---------------------------------------------------------------------------
def solve_theta(theta, frs, n_nodes=3):
    """Given the charge direction (a, b) = (cos, sin), solve the two meeting
    equations for (z0, x); return (residual_F3, z0, x, zs)."""
    ca = math.cos(theta)
    cb = math.sin(theta)
    cf = [abs(ca * f["X"] + cb * f["Y"]) for f in frs]
    if min(cf) < 1e-9:
        return None

    def eqs(z0, x):
        zs = [z0 + math.log(c) for c in cf]
        F = []
        for i in range(3):
            s = INV_A10 + B_SM[i] * x / TWO_PI
            for k, f in enumerate(frs):
                act = max(0.0, x - zs[k])
                if act > 0:
                    s += n_nodes * f["db"][i] * act / TWO_PI
            F.append(s - 1.0 / A_Z[i])
        return F, zs

    best = None
    for z0s, xs in ((20.0, 36.0), (15.0, 34.0), (25.0, 38.0), (10.0, 30.0)):
        z0, x = z0s, xs
        for _ in range(120):
            F, zs = eqs(z0, x)
            h = 1e-5
            J = np.zeros((2, 2))
            for j in range(2):
                Fp, _ = eqs(z0 + (h if j == 0 else 0), x + (h if j == 1 else 0))
                J[0, j] = (Fp[0] - F[0]) / h
                J[1, j] = (Fp[1] - F[1]) / h
            try:
                d = np.linalg.solve(J, [-F[0], -F[1]])
            except np.linalg.LinAlgError:
                break
            z0 = max(-10.0, min(60.0, z0 + d[0]))
            x = max(5.0, min(50.0, x + d[1]))
            if abs(d[0]) < 1e-10 and abs(d[1]) < 1e-10:
                break
        F, zs = eqs(z0, x)
        if abs(F[0]) < 1e-6 and abs(F[1]) < 1e-6:
            best = (F[2], z0, x, zs)
            break
    return best


def scan_roots(frs, n_nodes=3, n_theta=720):
    sols = []
    for k in range(n_theta):
        th = math.pi * k / n_theta
        r = solve_theta(th, frs, n_nodes)
        if r is not None:
            sols.append((th, *r))
    roots = []
    for a, b in zip(sols, sols[1:]):
        if a[1] * b[1] < 0:
            lo, hi = a, b
            for _ in range(60):
                mid = (lo[0] + hi[0]) / 2.0
                rm = solve_theta(mid, frs, n_nodes)
                if rm is None:
                    break
                if lo[1] * rm[0] <= 0:
                    hi = (mid, *rm)
                else:
                    lo = (mid, *rm)
            th = (lo[0] + hi[0]) / 2.0
            r = solve_theta(th, frs, n_nodes)
            if r and abs(r[0]) < 1e-6:
                res, z0, x, zs = r
                w = math.exp(z0) * 3.0 * M_Z / M_P
                Ms = [math.exp(z) * M_Z for z in zs]
                phys = (Z_MIN <= x <= X_MAX and w > 0
                        and all(Z_MIN <= z <= X_MAX for z in zs))
                roots.append({
                    "theta_rad": th, "w": w, "M_X_GeV": math.exp(x) * M_Z,
                    "thresholds_GeV": {f["label"]: Ms[i]
                                       for i, f in enumerate(frs)},
                    "physical": phys,
                })
    uniq = []
    for r in roots:
        if not any(abs(r["theta_rad"] - u["theta_rad"]) < 1e-3 for u in uniq):
            uniq.append(r)
    return uniq


def fragment_table(Xs, Ys):
    frs = {}
    for s in range(16):
        key = (round(Xs[s], 6), round(Ys[s], 6))
        frs.setdefault(key, {"X": key[0], "Y": key[1], "n": 0})["n"] += 1
    LABEL = {(1.0, 1 / 6): "Q(3,2)", (1.0, -2 / 3): "u^c", (1.0, 1.0): "e^c",
             (-3.0, 1 / 3): "d^c", (-3.0, -1 / 2): "L(1,2)", (5.0, 0.0): "N",
             (-3.0, -2 / 3): "u^c", (1.0, 1 / 3): "d^c", (1.0, 0.0): "N",
             (5.0, 1.0): "e^c", (-3.0, -1 / 2): "L(1,2)"}
    out = []
    for key, f in frs.items():
        n, X, Y = f["n"], f["X"], f["Y"]
        if n == 6:
            T2, T3 = 3 / 2, 1.0
        elif n == 3:
            T2, T3 = 0.0, 1 / 2
        elif n == 2:
            T2, T3 = 1 / 2, 0.0
        else:
            T2, T3 = 0.0, 0.0
        T1 = (3 / 5) * n * Y * Y
        f["label"] = LABEL.get(key, f"({X:g},{Y:g})")
        f["db"] = ((2 / 3) * 2 * T1, (2 / 3) * 2 * T2, (2 / 3) * 2 * T3)
        out.append(f)
    return sorted(out, key=lambda f: f["label"])

PHYS_GG = {(1.0, round(1 / 6, 6)), (1.0, round(-2 / 3, 6)), (1.0, 1.0),
           (-3.0, round(1 / 3, 6)), (-3.0, round(-1 / 2, 6)), (5.0, 0.0)}
PHYS_FL = {(1.0, round(1 / 6, 6)), (1.0, round(1 / 3, 6)), (1.0, 0.0),
           (-3.0, round(-2 / 3, 6)), (-3.0, round(-1 / 2, 6)), (5.0, 1.0)}


def physical_su5(Xs, Ys):
    pairs = {(round(Xs[s], 6), round(Ys[s], 6)) for s in range(16)}
    if pairs == PHYS_GG:
        return "Georgi-Glashow"
    if pairs == PHYS_FL:
        return "flipped"
    return None

def main():
    t0 = time.time()

    lam, offdiag = cartan_weights()
    x_coefs, y_coefs = identify_XY(lam)
    assert x_coefs and y_coefs, "X/Y not found"
    y_coef = y_coefs[0]
    Ys = lam @ y_coef
    phys = [(xc, physical_su5(lam @ xc, Ys)) for xc in x_coefs]
    phys = [(xc, kind) for xc, kind in phys if kind]
    assert phys, "no physical SU(5) X found"

    # ---- Part A: the Wilson factor on the implemented Floquet --------------
    _B, _U32, U16 = SC.sectors()
    checks = {}
    M0 = SC.bloch((0.0, 0.0, 0.0), U16)
    gap0_plain = gap_of(M0)
    # generic direction: the Georgi-Glashow X + Y
    xc_gg = next(xc for xc, kind in phys if kind == "Georgi-Glashow")
    ddiag = lam @ xc_gg + Ys
    checks["gap0"] = gap0_plain
    checks["zone_centre_unchanged"] = abs(
        gap0_plain - gap_of(step_wilson((0.0, 0.0, 0.0), U16, 0.3,
                                           ddiag))) < 1e-12
    U_F = step_wilson((0.7, -1.1, 2.0), U16, 0.3, ddiag)
    checks["unitarity"] = float(np.abs(
        U_F @ U_F.conj().T - np.eye(128, dtype=complex)).max())
    tail2, tail4 = wilson_kernel_decay(0.1, float(np.abs(ddiag).max()))
    checks["kernel_tail_r2_w01"] = tail2
    checks["kernel_tail_r4_w01"] = tail4
    _t2p, tail4_phys = wilson_kernel_decay(1e-8, float(np.abs(ddiag).max()))
    checks["kernel_tail_r4_wphys"] = tail4_phys
    node_gaps = {}
    for w in (0.01, 0.05, 0.1, 0.2):
        g = min(gap_of(step_wilson(q, U16, w, ddiag)) for q in SC.M_POINTS)
        node_gaps[w] = g
    checks["doubler_gaps_open"] = all(g > 1e-6 for g in node_gaps.values())
    checks["doubler_gaps"] = node_gaps
    # real-space light cone with the Wilson factor: the bond stages keep the
    # strict cone (1 Planck cell per stage, Theorem 5); the Wilson factor is
    # a block-momentum multiplication, so beyond the bond cone quasi-locality
    # means EXPONENTIAL decay of the amplitude in block distance.  The test
    # measures the amplitude in shells of one block and requires the decay
    # ratio per block to stay below 1/2; the single-factor kernel |J_m(w c)|
    # of the largest Wilson charge is reported as the reference scale.
    radii, outside, shells = light_cone_wilson(32, 24, 0.1, ddiag, U16)
    checks["front_radius_at_24"] = int(radii[-1])
    checks["amplitude_outside_LR_cone"] = outside
    checks["tail_shells_beyond_cone"] = shells
    ratios = [shells[k + 1] / shells[k] for k in range(len(shells) - 1)
              if shells[k] > 0.0 and shells[k + 1] > 0.0]
    checks["tail_decay_ratio_max"] = max(ratios)
    cmax = float(np.abs(ddiag).max())
    checks["wilson_charge_max"] = cmax
    checks["bessel_kernel_Jm"] = [abs(v) for v
                                  in bessel_J(0.1 * cmax, kmax=len(shells))][1:]

    # ---- Part B -------------------------------------------------------------
    results_b = {}
    for i, (xc, kind) in enumerate(phys):
        Xs = lam @ xc
        fr_list = fragment_table(Xs, Ys)
        db_uniform = [sum(12 * f["db"][j] for f in fr_list) for j in range(3)]
        uniq = scan_roots(fr_list)
        results_b[f"X_{i}"] = {
            "X_coef": xc.tolist(),
            "su5": kind,
            "fragments": [{"label": f["label"], "X": f["X"], "Y": f["Y"],
                           "n": f["n"],
                           "db": [round(v, 4) for v in f["db"]]}
                          for f in fr_list],
            "uniform_overproduction_delta_b": [round(v, 2)
                                               for v in db_uniform],
            "selective_roots": uniq,
        }

    ok = (checks["zone_centre_unchanged"] and checks["doubler_gaps_open"]
          and checks["unitarity"] < 1e-10
          and checks["kernel_tail_r4_wphys"] < 1e-15
          and checks["amplitude_outside_LR_cone"] < 1e-6
          and checks["tail_decay_ratio_max"] < 0.5
          and checks["tail_shells_beyond_cone"][0] < 1e-3)
    n_phys = sum(1 for r in results_b.values()
                 for root in r["selective_roots"] if root["physical"])
    verdict = {
        "experiment": "V65_wilson",
        "pass": ok,
        "status": ("W9 MECHANISM FOUND: an internal Wilson factor of the "
                   "SAME step, exactly unitary with W(0)=I, gaps all 12 "
                   "doubler nodes; strictly causal bond stages plus a "
                   "QUASI-local Wilson factor - beyond the bond cone the "
                   "amplitude decays exponentially in block distance "
                   f"(worst ratio {checks['tail_decay_ratio_max']:.2f} per "
                   f"block, first shell {checks['tail_shells_beyond_cone'][0]:.1e}; "
                   f"front radius {checks['front_radius_at_24']} Planck "
                   "cells after 24 stages, "
                   f"{checks['amplitude_outside_LR_cone']:.1e} past the "
                   "Lieb-Robinson line); X and Y identified as Cartan "
                   "combinations of the root (both SU(5) embeddings found); "
                   "W4: uniform w overproduces (Delta b = 32), the "
                   f"selective 3-node reading has {n_phys} physical roots",
                   ),
        "data": {
            "cartan_offdiag": offdiag,
            "Y_coef": y_coef.tolist(),
            "part_A_checks": checks,
            "part_B": results_b,
        },
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v65_wilson.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V65 - the Wilson mechanism (W9) and the desert thresholds (W4)")
    print("=" * 68)
    print(f"  Cartan offdiag = {offdiag:.2e};  Y_coef = {y_coef}")
    print(f"  physical SU(5) X candidates: {len(phys)}")
    print(f"  [A] zone-centre gap {gap0_plain:.4f} unchanged: "
          f"{checks['zone_centre_unchanged']}; unitarity "
          f"{checks['unitarity']:.2e}")
    print(f"      kernel tail r>=2: {tail2:.2e}, r>=4: {tail4:.2e} "
          f"(quasi-local, Bessel)")
    print(f"      doubler gaps vs w: "
          + ", ".join(f"w={w}: {g:.4f}" for w, g in node_gaps.items()))
    print(f"      light cone at 24 steps: radius {radii[-1]}, amplitude "
          f"outside Lieb-Robinson cone {outside:.2e}")
    print("      tail beyond the bond cone (block shells): "
          + ", ".join(f"m={m + 1}: {s:.1e}"
                      for m, s in enumerate(shells)))
    print(f"      worst decay ratio {checks['tail_decay_ratio_max']:.2f} per "
          f"block; single-factor kernel |J_m(w cmax)| for cmax = {cmax:.3f}: "
          + ", ".join(f"{p:.1e}" for p in checks["bessel_kernel_Jm"]))
    for name, r in results_b.items():
        print(f"  [B {name} ({r['su5']})] uniform Delta b = "
              f"{r['uniform_overproduction_delta_b']} -> OVERPRODUCED")
        print(f"      selective 3-node roots: {len(r['selective_roots'])}")
        for root in r["selective_roots"]:
            print(f"        theta = {root['theta_rad']:.3f}  "
                  f"w = {root['w']:.2e}  M_X = {root['M_X_GeV']:.2e}  "
                  f"physical = {root['physical']}")
            for lab, M in root["thresholds_GeV"].items():
                print(f"          {lab:8s} M = {M:.2e} GeV")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - {verdict['status'][0]}")
    print(f"  JSON -> " + os.path.join(outdir, "v65_wilson.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
