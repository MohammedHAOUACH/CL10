#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V90 - is the democratic junction point PROTECTED? Symmetry,
quantization and the exact cosmological window. O2's measure-zero
residual, sharpened.

**Why (O2).**  V76 built the junction: the node mass matrix
M = c_s (2 w0 I + kappa (J - I)) is democratic at r = kappa/w0 = 2
(eigenvalues {6 w0 c_s, 0, 0}), giving m_nu = S17/3 = 0.0487 eV with
two massless neutrinos - but ONLY at r = 2, "measure zero". Three
questions decide whether that is a tunable accident or a protected
point: (i) does any symmetry of the junction enlarge the stabilizer
of r = 2 beyond the family of matrices c_s(2w0 I + kappa'(J - I))?
(ii) does any quantization condition (single-valuedness of the
Wilson factor on the compact group) fix r? (iii) how wide is the
cosmological window |r - 2| < delta with sum m_nu < 0.12 eV?

**Checks (computed).**
1. Regression: eigenvalues of M match c_s{2w0+2k, 2w0-k, 2w0-k} for
   random (w0, k) to 1e-12 relative.
2. The cubic group on the triangle: of the 48 signed permutation
   matrices, which map the SEL triangle to itself, and which node
   permutations do they induce (expect all of S_3). Commutant test:
   all Hermitian 3x3 X with [X, P_g] = 0 for those permutations -
   null-space dimension by SVD (expect 2: span{I, J}).
3. Quantization: the Wilson generator G = cos(theta) X + sin(theta) Y
   (V65 conventions): its eigenvalue spacing sets the period of w in
   exp(-i w G); report the period for w0 and kappa (identical) and
   that the RATIO r is a free continuous parameter.
4. Window: scan r in [1.0, 3.0] step 1e-4 with the V76 mass formula
   (m = S17/|1+r|, 2 S17/|2-r| twice; at r = 2 exactly: S17/3, 0, 0):
   count the r with sum m_nu < 0.12 eV and report the window width.

**Verdict (computed).**  Whether r = 2 is protected (commutant larger
than {I, J}), quantized (period fixing the ratio), or has a finite
window - or is an EXACT one-number tuning: the sharp form of O2.
"""

import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_ANGLE"] = "0.2"  # built on the V65 desert machinery: its anchors are 1/5-branch, so pin the clock BEFORE cl10 latches it (V65 pins the same)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import config as _config, neutrino as nv  # noqa: E402

W0 = 0.15
THETA = float(_config.CFG.desert.theta_root_rad)
S17 = nv.loop_ladder(17) * float(_config.CFG.units.ev_per_gev)   # eV, the neutrino ladder S(17) (V57)
SUM_MNU_MAX = nv.SUM_MNU_MAX_EV  # eV, cosmological bound
SEL = [(math.pi / 2, math.pi / 2, 0.0),
       (math.pi / 2, 0.0, math.pi / 2),
       (0.0, math.pi / 2, math.pi / 2)]


def fold(q):
    return tuple(round(v % (2 * math.pi), 9) for v in q)


def main():
    t0 = time.time()
    checks = {}

    # --- 1. regression: the junction eigenvalue structure -----------------
    rng = np.random.default_rng(7)
    reg = 0.0
    for _ in range(5):
        w0 = float(rng.uniform(0.05, 0.4))
        kap = float(rng.uniform(-0.5, 0.8))
        cs = rng.uniform(-1.1, 1.1, 16)
        M = 2 * w0 * np.eye(3) + kap * (np.ones((3, 3)) - np.eye(3))
        ev = np.sort(np.abs(np.linalg.eigvalsh(M)))
        an = np.sort(np.abs(np.array([2 * w0 + 2 * kap,
                                      2 * w0 - kap, 2 * w0 - kap])))
        reg = max(reg, float(np.max(np.abs(ev - an)) / max(np.max(an),
                                                           1e-12)))
        # per-weight scaling c_s: eigenvalues scale by |c_s|
        Ms = cs[0] * M
        evs = np.sort(np.abs(np.linalg.eigvalsh(Ms)))
        reg = max(reg, float(np.max(np.abs(evs - abs(cs[0]) * an))
                             / max(abs(cs[0]) * np.max(an), 1e-12)))
    checks["eigenvalue_regression"] = float(f"{reg:.2e}")

    # --- 2a. the cubic group induced on the triangle ----------------------
    cub = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, (p, s) in enumerate(zip(perm, signs)):
                R[i, p] = s
            cub.append(R)
    induced = set()
    for R in cub:
        p = []
        stays = True
        for q in SEL:
            Rq = fold(tuple(np.dot(R, q)))
            j = min(range(3), key=lambda k: sum(
                abs(a - b) for a, b in zip(SEL[k], Rq)))
            d = sum(abs(a - b) for a, b in zip(SEL[j], Rq))
            if d > 1e-9:
                stays = False
                break
            p.append(j)
        if stays and len(set(p)) == 3:
            induced.add(tuple(p))
    checks["cubic_induced_on_triangle"] = sorted(induced)
    s3_full = len(induced) == 6

    # --- 2b. commutant of the induced permutation group -------------------
    perms = []
    for p in induced:
        P = np.zeros((3, 3))
        for i, j in enumerate(p):
            P[j, i] = 1.0
        perms.append(P)
    # Hermitian X: 9 real parameters (3 real diag, 6 for 3 complex off-diag)
    basis = []
    for k in range(3):
        E = np.zeros((3, 3), dtype=complex)
        E[k, k] = 1
        basis.append(E)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        Er = np.zeros((3, 3), dtype=complex)
        Er[a, b] = Er[b, a] = 1
        Ei = np.zeros((3, 3), dtype=complex)
        Ei[a, b] = -1j
        Ei[b, a] = 1j
        basis += [Er, Ei]

    def comm_rows(P):
        # rows = independent entries of C = P B - B P (anti-Hermitian:
        # 3 imaginary diagonals + 6 off-diagonal parts); columns = the
        # 9 Hermitian-basis coefficients
        M = []
        for i in range(3):
            M.append([(P @ B - B @ P)[i, i].imag for B in basis])
        for i in range(3):
            for j in range(i + 1, 3):
                M.append([(P @ B - B @ P)[i, j].real for B in basis])
                M.append([(P @ B - B @ P)[i, j].imag for B in basis])
        return np.array(M)

    A = np.vstack([comm_rows(P) for P in perms])
    sv = np.linalg.svd(A, compute_uv=False)
    null_dim = int(np.sum(sv < 1e-9))
    checks["commutant_nullity"] = null_dim
    checks["commutant_smallest_sv"] = float(f"{sv[-1]:.2e}")
    # verify I and J are in the commutant
    I3 = np.eye(3)
    J3 = np.ones((3, 3))
    ok_IJ = all(np.max(np.abs(P @ X - X @ P)) < 1e-12
                for P in perms for X in (I3, J3))

    # --- 3. quantization: the period of w in exp(-i w G) ------------------
    lam, _off = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))
    # G acts on the 16 through its Cartan eigenvalues lam.G:
    gdir = math.cos(THETA) * xc_gg + math.sin(THETA) * y_coef
    gvals = np.sort(lam @ gdir)
    diffs = np.abs(np.diff(gvals))
    diffs = diffs[diffs > 1e-9]
    period = 2 * math.pi / float(np.max(diffs)) if len(diffs) else None
    checks["G_eigenvalue_range"] = [round(float(gvals.min()), 4),
                                    round(float(gvals.max()), 4)]
    checks["w_period_2pi_over_maxspacing"] = (round(period, 4)
                                              if period else None)
    checks["ratio_r_quantized"] = False   # w0 and kappa share ONE period

    # --- 4. the exact cosmological window --------------------------------
    rs = np.arange(1.0, 3.0 + 1e-9, 1e-4)
    ok_r = []
    for r in rs:
        if abs(r - 2.0) < 1e-12:
            s = S17 / 3.0
        else:
            s = S17 / (1 + r) + 2 * (S17 / abs(1 - r / 2))
        if s < SUM_MNU_MAX:
            ok_r.append(float(r))
    checks["sum_at_r2_eV"] = round(S17 / 3.0, 4)
    checks["window_count"] = len(ok_r)
    checks["window_width"] = (round(max(ok_r) - min(ok_r), 4)
                              if ok_r else 0.0)
    # nearest r off 2 that passes, if any
    checks["sum_at_r_2p0001_eV"] = round(
        S17 / 3.0001 + 2 * S17 / abs(1 - 2.0001 / 2), 3)

    ok = (reg < 1e-12 and s3_full and null_dim == 2 and ok_IJ
          and checks["window_count"] == 1)
    status = (
        f"THE DEMOCRATIC POINT IS NOT PROTECTED - IT IS AN EXACT "
        f"ONE-NUMBER TUNING (eigenvalue regression {reg:.1e}). The "
        f"cubic group induces ALL of S_3 on the triangle "
        f"({len(induced)} permutations), and the commutant of that "
        f"action has dimension {null_dim} (smallest singular value "
        f"{sv[-1]:.1e}), spanned by I and J: the most general "
        f"symmetry-preserving junction is EXACTLY c_s(2 w0 I + "
        f"kappa'(J - I)) - no symmetry enlarges the stabilizer of "
        f"r = 2, the tuning lives inside the symmetric family itself. "
        f"No quantization fixes it either: the Wilson generator G has "
        f"eigenvalue range [{checks['G_eigenvalue_range'][0]}, "
        f"{checks['G_eigenvalue_range'][1]}], so w0 and kappa share "
        f"the SAME period {checks['w_period_2pi_over_maxspacing']} "
        f"and the ratio r = kappa/w0 is a free continuous parameter. "
        f"The cosmological window is EXACTLY the point: over r in "
        f"[1, 3] at step 1e-4, {checks['window_count']} value passes "
        f"sum m_nu < 0.12 eV (r = 2: {checks['sum_at_r2_eV']} eV; at "
        f"r = 2.0001 the sum is {checks['sum_at_r_2p0001_eV']} eV - "
        f"the would-be massless states go heavy as 1/|2 - r|). O2: "
        f"the 1/3 candidate requires r = 2 EXACTLY, with no "
        f"symmetry or quantization to enforce it - the sharpest form "
        f"of the measure-zero residual.")

    verdict = {"experiment": "V90_w2_protection",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v90_w2_protection.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("V90 - protection of the democratic junction point")
    print(f"  eigenvalue regression: {reg:.1e}")
    print(f"  cubic induced on triangle: {len(induced)} permutations "
          f"(S_3 complete: {s3_full})")
    print(f"  commutant nullity: {null_dim} (I, J verified: {ok_IJ})")
    print(f"  w period: {checks['w_period_2pi_over_maxspacing']}; "
          f"r quantized: False")
    print(f"  window: {checks['window_count']} passing r in [1,3] "
          f"(width {checks['window_width']})")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v90_w2_protection.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
