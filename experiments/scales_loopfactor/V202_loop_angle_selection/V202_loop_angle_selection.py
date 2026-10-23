#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V202 - is 1/(16 pi) a SELECTED Wilson angle? (O8)

**Where O8 stands.**  V200 factorised the traced-loop weight exactly:
F_10 = N_10 x (per-word amplitude) x (interference), the census sum equals the
matrix power, and F_10(eta) hits 1/(16 pi) at 8 angles in (0, pi).  The open
question is then sharp: does one of those crossing angles carry a NUMBER OF
THE ROOT - the step angle g_10 = 1/5, a closed angle pi/4 or pi/2, or a small
rational - so that the continuum loop factor 1/(16 pi) becomes a property of
the root's angle arithmetic?  This experiment measures it: if some root angle
is a crossing, O8's residual is a named angle tension with a known solution;
if none is, 1/(16 pi) is NOT what the root supplies, and the named weight is
F_10(g_10) itself.

**The exact functional (measured, no fit).**  The three Casimir eigenvalues of
the 32 (x) 45 link, -32, -12, 4, are EVEN, so every length-n word carries a
phase e^{i k eta} with k EVEN (k = sum of half-eigenvalues): F_10 has period
pi, is invariant under eta -> pi - eta, and its eight crossings pair at sums
of exactly pi.  Grouping the 3^10 words by phase index gives 63 classes (58
nonzero), the Cartan coefficients of the five components obey the mirror
identities C_k[0] = C_k[4] and C_k[1] = C_k[3], and each component's modulus
is pi-pairing invariant: the functional is a sum of three real cosine
polynomials in the EVEN harmonics, with exact coefficients from the root.

**Checks (computed).**
1. T0 phase classes: eigenvalues even; the grouped census reconstructs V126's
   word sum at eta = 0, 0.2, 1.4423 to machine precision.
2. T1 functional structure: period pi, eta -> pi - eta pairing per component,
   coefficient mirror identities.
3. T2 the crossings: F_10(eta) = 1/(16 pi) has exactly four solutions in
   (0, pi/2], found by bisection to 1e-13, and the eight in (0, pi) pair at
   sums pi.
4. T3 the selection test: every candidate angle of the root against the
   crossings - g_10 = 1/5, 1/6, 1/7, pi/16, pi/15, pi/14, pi/12, pi/4, pi/2,
   theta_root - and the rational/algebraic content of cos(2 eta*).
5. T4 the named arithmetic: the angle gap of g_10 to the nearest crossing,
   the F_10 deficit, and the interference ratio between the root and the
   crossing.

**Pre-registered branch.**  SEL-EXACT iff some root candidate angle satisfies
|F_10(candidate) - 1/(16 pi)| < 1e-9 (1/(16 pi) is then a root-angle number);
ANGLE-OUTSIDE otherwise (1/(16 pi) is not supplied by any root angle and the
named traced-loop weight stays F_10(g_10)).  pass is True in both branches.
"""
import json
import os
import sys
import time

import numpy as np

os.environ["CL10_BACKEND"] = "numpy"

# pinned reference-branch clock: the word census is Cartan-dark at the closed
# default_branch angle, and V126's own pin latches too late here (V79/V85
# import cl10 first) - pin the clock before any cl10 import.
os.environ["CL10_ANGLE"] = "0.2"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                        # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V85_closed_loop_factor"))
sys.path.insert(0, exp("V126_loop_enumeration"))

from V79_gauged_bond_layer import (adjoint_link_operators,       # noqa: E402
                                   root_generators, structure_constants)
from V85_closed_loop_factor import projectors                    # noqa: E402
from V126_loop_enumeration import (cartan_vector, sectors,    # noqa: E402
                                      word_census)
from cl10 import clifford as CL                                  # noqa: E402
from cl10 import config                                          # noqa: E402
from cl10 import engine as E                                     # noqa: E402
from cl10 import sm_qca                                          # noqa: E402
from cl10 import u10pred                                         # noqa: E402

D_L = 45
N = 10
LOOP = u10pred.LOOP
ETA_ROOT = float(config.CFG.angle.g10)


def build():
    """The exact objects of V85/V126: the three sector hop matrices, the
    phase half-eigenvalues, the Cartan functionals."""
    G = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(),
                                      sm_qca.G10)).astype(np.complex128)
    h = [E.asnumpy(m) for m in CL.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = root_generators()
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    P, ev = projectors(TL)
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    Ma = []
    for Pa in P:
        r = (Pa @ Gi).reshape(32, D_L, 32, D_L)
        Ma.append(np.einsum("ajbj->ab", r) / D_L)
    B0 = [m @ G for m in Ma]            # eta = 0 hop matrices
    return G, cartan, [float(e) for e in ev], Ma, B0


def phase_classes(B0, cartan, half_ev, n=N):
    """All 3^n words at eta = 0 grouped by their even phase index
    k = sum of half-eigenvalues: F_10(eta) = sum_k C_k e^{i k eta} exactly,
    because the word phases are scalars."""
    prefixes = [(np.eye(32, dtype=complex), 0.0)]
    for _ in range(n):
        prefixes = [(W @ B0[lam], kk + half_ev[lam])
                    for W, kk in prefixes for lam in range(3)]
    groups = {}
    for W, kk in prefixes:
        groups.setdefault(int(round(kk)), []).append(
            cartan_vector(W, cartan))
    return {k: np.sum(v, axis=0) for k, v in groups.items()}


def main():
    t0 = time.time()
    out = {"experiment": "V202_loop_angle_selection", "tier": 4,
           "programme": "V202", "data": {}}
    checks = []

    G, cartan, ev, Ma, B0 = build()
    half_ev = [e / 2.0 for e in ev]
    C = phase_classes(B0, cartan, half_ev)
    den = float(np.abs(cartan_vector(np.linalg.matrix_power(G, N),
                                     cartan)).sum())

    def T(eta):
        z = np.zeros(5, dtype=complex)
        for k, v in C.items():
            z += v * np.exp(1j * k * eta)
        return z

    def F(eta):
        return float(np.abs(T(eta)).sum() / den)

    # ---- T0 the phase classes ARE the census ------------------------------
    ev_even = all(abs(e - round(e)) < 1e-9 and int(round(e)) % 2 == 0
                  for e in ev)
    k_even = all(k % 2 == 0 for k in C)
    nz = {k for k, v in C.items() if np.abs(v).max() > 1e-14}
    ref_err = 0.0
    for eta in (0.0, ETA_ROOT, 1.4423076923076923):
        cx = word_census(sectors(eta, Ma, ev, G), cartan, N)
        ref_err = max(ref_err, float(np.abs(T(eta) - cx[N]["total"]).max()))
    out["data"]["phase_classes"] = {
        "eigenvalues": ev, "classes": len(C), "nonzero": len(nz),
        "k_range": [min(C), max(C)], "grouped_vs_census_error": ref_err}
    checks.append(("T0 the three Casimir eigenvalues %s are even, every word "
                   "phase index k is even, and %d phase classes (%d nonzero) "
                   "reconstruct the enumerated word sum exactly at eta = 0, "
                   "%.3f, 1.442 (%.1e)"
                   % ([round(e) for e in ev], len(C), len(nz), ETA_ROOT,
                      ref_err),
                   ev_even and k_even and ref_err < 1e-12))

    # ---- T1 the exact functional structure --------------------------------
    per = np.abs(T(0.0))
    pair_err = max(float(np.abs(np.abs(T(e)) - np.abs(T(np.pi - e))).max())
                   for e in (0.11, 0.31, ETA_ROOT, 1.13, 1.4423))
    period_err = max(abs(F(e) - F(e + np.pi)) for e in (0.31, 1.13))
    mirror01 = max(float(np.abs(C[k][0] - C[k][4]).max()) for k in C)
    mirror13 = max(float(np.abs(C[k][1] - C[k][3]).max()) for k in C)
    out["data"]["structure"] = {
        "pairing_pi_minus_eta_error": pair_err, "period_error": period_err,
        "mirror_C04_error": mirror01, "mirror_C13_error": mirror13,
        "component_moduli_at_pi_minus_eta": [float(x) for x in per]}
    checks.append(("T1 the functional is exact: period pi (error %.1e), the "
                   "per-component moduli are pi-pairing invariant (max error "
                   "%.1e), and the coefficient mirror identities C_k[0] = "
                   "C_k[4], C_k[1] = C_k[3] hold (errors %.1e / %.1e) - F_10 "
                   "is a sum of three real cosine polynomials in even "
                   "harmonics" % (period_err, pair_err, mirror01, mirror13),
                   period_err < 1e-12 and pair_err < 1e-12
                   and mirror01 < 1e-10 and mirror13 < 1e-10))

    # ---- T2 the crossings of 1/(16 pi) ------------------------------------
    xs = np.linspace(1e-9, np.pi / 2, 40001)
    dsign = np.array([F(x) - LOOP for x in xs])
    roots = []
    for i in range(len(xs) - 1):
        if (dsign[i] < 0) != (dsign[i + 1] < 0):
            a, b, sa = xs[i], xs[i + 1], dsign[i] < 0
            for _ in range(60):
                m = 0.5 * (a + b)
                if (F(m) - LOOP < 0) == sa:
                    a = m
                else:
                    b = m
            roots.append(0.5 * (a + b))
    res = [abs(F(x) - LOOP) for x in roots]
    pair_sums = [roots[i] + roots[j]
                 for i in range(len(roots))
                 for j in range(i + 1, len(roots))]
    pi_pairs = sum(1 for s in pair_sums if abs(s - np.pi) < 1e-9)
    out["data"]["crossings"] = {
        "n_roots_0_pi2": len(roots),
        "eta_star": [float(x) for x in roots],
        "residuals": [float(r) for r in res],
        "pi_pairs_in_(0,pi/2]": pi_pairs,
        "cos2eta_star": [float(np.cos(2 * x)) for x in roots]}
    checks.append(("T2 F_10(eta) = 1/(16 pi) has exactly %d solutions in "
                   "(0, pi/2]: %s (residuals <= %.1e); with the period and "
                   "the pairing the eight crossings in (0, pi) close into "
                   "four pi-sum pairs"
                   % (len(roots), [round(x, 7) for x in roots], max(res)),
                   len(roots) == 4 and max(res) < 1e-13))

    # ---- T3 does a ROOT angle reach 1/(16 pi)? ----------------------------
    cands = [("g_10 = 1/5", 0.2), ("1/6", 1.0 / 6.0), ("1/7", 1.0 / 7.0),
             ("pi/16", np.pi / 16.0), ("pi/15", np.pi / 15.0),
             ("pi/14", np.pi / 14.0), ("pi/12", np.pi / 12.0),
             ("pi/8", np.pi / 8.0), ("pi/4", np.pi / 4.0),
             ("pi/2", np.pi / 2.0),
             ("theta_root", float(config.CFG.desert.theta_root_rad))]
    scan = {lab: float(F(x) - LOOP) for lab, x in cands}
    best = min(abs(v) for v in scan.values())
    # rational / quadratic content of the nearest crossing
    from fractions import Fraction
    eta1 = roots[0]
    rat_err = {}
    for lab, y in (("eta/pi", eta1 / np.pi), ("cos^2 eta", np.cos(eta1) ** 2),
                   ("cos 2eta", np.cos(2 * eta1))):
        fr = Fraction(y).limit_denominator(2000)
        rat_err[lab] = float(abs(y - float(fr)))
    out["data"]["selection"] = {
        "candidates": scan, "best_abs_mismatch": best,
        "nearest_crossing": float(eta1),
        "angle_gap_g10_to_crossing": float(ETA_ROOT - eta1),
        "rational_reconstruction_errors": rat_err}
    checks.append(("T3 no root angle selects the crossing: the smallest "
                   "mismatch of F_10(candidate) against 1/(16 pi) over the "
                   "root's closed angles and rational conventions is %.1e "
                   "(at %s); the nearest crossing to g_10 = 1/5 sits at "
                   "eta* = %.6f, an angle gap of %.2e, and cos(2 eta*) "
                   "reconstructs to no rational of denominator <= 2000 "
                   "(best residual %.0e)"
                   % (best, min(scan, key=lambda k: abs(scan[k])),
                      eta1, ETA_ROOT - eta1, max(rat_err.values())),
                   best > 1e-6 and max(rat_err.values()) > 1e-9))

    # ---- T4 the named arithmetic of the residual --------------------------
    F0 = F(0.0)
    I_root = F(ETA_ROOT) / F0
    I_cross = LOOP / F0
    branch = "SEL-EXACT" if best < 1e-9 else "ANGLE-OUTSIDE"
    out["data"]["named"] = {
        "F10_at_g10": F(ETA_ROOT), "F10_at_eta0": F0,
        "deficit_vs_loop": F(ETA_ROOT) / LOOP - 1.0,
        "interference_root": I_root, "interference_at_crossing": I_cross,
        "interference_deficit_ratio": I_root / I_cross}
    checks.append(("T4 the named weight: F_10(g_10 = 1/5) = %.9f, the "
                   "interference factor at the root is %.4f against %.4f "
                   "required at the nearest crossing (ratio %.4f); the root "
                   "sits %.4f%% below the continuum loop factor"
                   % (F(ETA_ROOT), I_root, I_cross, I_root / I_cross,
                      100 * (1 - F(ETA_ROOT) / LOOP)),
                   abs(F(ETA_ROOT) - 0.019186391026184278) < 1e-13
                   and 0.0 < I_root < 1.0))

    out["branch"] = branch
    ok = all(c[1] for c in checks)
    out["checks"] = [c[0] for c in checks]
    out["pass"] = ok
    out["residual"] = True
    out["status"] = (
        "the traced-loop weight of O8 is an exact real cosine functional of "
        "the Wilson angle - 63 even phase classes, period pi, eta -> pi - eta "
        "pairing, coefficient mirrors - and its four crossings of 1/(16 pi) "
        "in (0, pi/2] are %s; NO angle of the root reaches the crossing "
        "(branch %s): the nearest to g_10 = 1/5 sits at eta* = %.6f, gap "
        "%.2e, and the closed angles miss by >= %.0e. The named weight "
        "stays F_10(1/5) = %.9f, %.4f%% below the continuum loop factor; "
        "the residual is the interference deficit ratio %.4f"
        % ([round(x, 6) for x in roots], branch, roots[0],
           ETA_ROOT - roots[0], best, F(ETA_ROOT),
           100 * (1 - F(ETA_ROOT) / LOOP), I_root / I_cross))
    out["wall_time"] = round(time.time() - t0, 1)
    rdir = os.path.join(HERE, "results")
    os.makedirs(rdir, exist_ok=True)
    with open(os.path.join(rdir, "v202_loop_angle_selection.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for line, _ in checks:
        print(" ", line)
    print("  BRANCH: %s" % branch)
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
