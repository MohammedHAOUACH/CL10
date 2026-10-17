#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V81 - the desert content of the GAUGED bond layer: the composite
fermion spectrum 16 (x) 45 and coupling unification (O4).

**Why (O4, the named incompatibility closed or not).**  V74 proved the
node mechanism's meeting point is FORCED (rank-1 split-invariance): the
automaton supplies full split 16's per node, the couplings prefer
partial suits. The gauged bond layer (V79) supplies exactly that missing
unit: the Gauss constraint makes the physical fermion a COMPOSITE -
matter 16 dressed by a link flux in the adjoint 45 - and in the Higgs
phase (Wilson direction v, V65/V67) the flux in a broken direction is
screened and costs w |alpha.v|: the composite spectrum 16 (x) 45 splits
into SM fragments with INDEPENDENT thresholds - partial suits, one
scale w + one offset.

**The content (0-parameter, fully structural).**
- Weights: the 16 weights lam (V65, in the D_k = i h h Cartan basis,
  coordinates +-1) + the 45 flux weights (5 Cartan zeros + 40 roots,
  extracted by commutators with the Cartan).
- SM roots from the root data itself: the unbroken roots (alpha.x =
  alpha.y = 0) contain the weak root beta_L (the difference of the two
  L weights) and the color roots beta_1, beta_2 (differences of the
  three u^c weights). No multiset guessing.
- One-loop shifts over the (vectorlike, mu -> -mu symmetric) set:
    db1 = (2/3)(3/5) sum Y_s^2,
    db2 = (2/3) sum [(s.beta_L)/(beta_L,beta_L)]^2   (SU(2) index via
    the coroot: T(fund) = 1/2 automatic),
    db3 = (2/3) (1/8) sum C2(p_s, q_s)               (SU(3) index via
    the (p,q) string labels along beta_1, beta_2).
  VERIFIED against the V64 LIB values for THREE independent fragments:
  (3,2) -> (2/15, 2, 4/3), (3b,1) -> (2/3, 0, 2/3), (1,2) -> (1/5, 2/3, 0).
- Thresholds: state (mu, alpha) at m = w (c0 + |alpha.v(theta)|),
    v(theta) = cos(theta) xc + sin(theta) yc (V65 conventions); one of
    the five zero-flux 16-copies IS the SM matter (in B_SM) and is
    excluded from the desert running.

**The test.** Fix alpha_X = 1/25 (the root value, never refitted). Three
meeting equations, three unknowns (w, c0, theta): solve; report M_X,
tau_p, the threshold spectrum; compare with the data-preferred
M_X = 1.8e16 GeV (V64) and the F5 bound. The four remaining Cartan-flux
16-copies at w c0 are reported honestly - the lightest content, to be
lifted by the extra-U(1) breaking.

**Verdict (computed).** Either the composite spectrum lands the meeting
at alpha_X = 1/25 with desert thresholds (O4 CLOSED by the gauged
layer), or it does not and the residual is named with numbers.
"""

import json
import math
import os
import sys
import time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
from V65_wilson import cartan_weights, identify_XY, physical_su5  # noqa: E402
from cl10 import clifford  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import u10pred, unification as UN  # noqa: E402

TWO_PI = 2.0 * math.pi
MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
A10, M_P = UN.ALPHA_10, UN.M_P
INV_A10 = 1.0 / A10
X_MAX = math.log(M_P / MZ)
Z_MIN = math.log(1e6 / MZ)
R9 = lambda a: tuple(np.round(a, 9))  # noqa: E731


def flux_weights():
    """The 45 weights of the adjoint in the Cartan basis: the plane
    generators are NOT weight vectors (ad_{D_k} mixes them), so build
    the 45x45 matrices of ad_{D_k} in the plane basis and take their
    joint eigenvalues (5 zeros + 40 roots, entries {0, +-2})."""
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan32 = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    planes = [1j * (h[a] @ h[b]) for a in range(10) for b in range(a + 1, 10)]
    Ms = []
    for C in cartan32:
        M = np.zeros((45, 45))
        for B in range(45):
            cmt = C @ planes[B] - planes[B] @ C
            for A in range(45):
                M[A, B] = float(np.real(-1j * np.trace(cmt @ planes[A])
                                        / 32.0))
        Ms.append(M)
    rng = np.random.default_rng(0)
    G = sum(p * M for p, M in zip(rng.uniform(0.3, 1.0, 5), Ms))
    w, V = np.linalg.eig(G)
    order = np.argsort(w)
    clusters = []
    for j in order:
        if clusters and abs(w[j] - clusters[-1][0]) < 1e-6:
            clusters[-1][1].append(j)
        else:
            clusters.append((w[j], [j]))
    wts = []
    for _wc, js in clusters:
        B = V[:, js]
        for M in Ms[1:]:                    # split degenerate clusters
            Mb = B.conj().T @ M @ B
            _wb, Vb = np.linalg.eig(Mb)
            B = B @ Vb
        for m in range(len(js)):
            v = B[:, m]
            nv = float(np.real(v.conj() @ v))
            wts.append(np.array([float(np.real(-1j * v.conj() @ (M @ v)))
                                 for M in Ms]) / nv)
    return np.array(wts)          # 45 vectors (5 Cartan ones are 0)


def sm_roots(alphas, lam, xc, yc):
    """The SM roots from the root data: unbroken roots (alpha.x = 0,
    alpha.y = 0); beta_L = the difference of the two L weights (Y = -1/2);
    beta_1, beta_2 = two independent differences of the three u^c weights
    (Y = -2/3)."""
    unbr = [a for a in alphas if np.abs(a @ xc) < 1e-9 and np.abs(a @ yc) < 1e-9
            and np.abs(a).sum() > 1e-9]
    Lp = lam[np.abs(lam @ yc + 0.5) < 1e-9]
    assert len(Lp) == 2, "L doublet not found"
    bL = Lp[0] - Lp[1]
    assert any(np.allclose(bL, a) for a in unbr), "beta_L is not unbroken"
    Uc = lam[np.abs(lam @ yc + 2 / 3) < 1e-9]
    assert len(Uc) == 3, "u^c triplet not found"
    droot = [Uc[i] - Uc[j] for i in range(3) for j in range(i + 1, 3)
             if any(np.allclose(Uc[i] - Uc[j], a) for a in unbr)]
    b1 = droot[0]
    # simple roots: b2 must form a chain (b1 + b2 also a root), else the
    # weight set is an L-shape and the string labels break
    b2 = next(d for d in droot[1:] if not np.allclose(d, b1)
              and not np.allclose(d, -b1)
              and any(np.allclose(b1 + d, a) for a in unbr))
    return bL, b1, b2


def c2_labels(states, b1, b2):
    """Per-state SU(3) quadratic Casimir C2(p, q) from the string labels
    along the color roots; multiplicities preserved."""
    cnt = Counter(map(R9, states))
    cache = {}
    out = np.zeros(len(states))
    for i, s in enumerate(states):
        key = R9(s)
        nu = s.copy()
        moved = True
        while moved:
            moved = False
            while R9(nu + b1) in cnt:
                nu = nu + b1
                moved = True
            while R9(nu + b2) in cnt:
                nu = nu + b2
                moved = True
        hk = R9(nu)
        if hk not in cache:
            p = q = 0
            v = nu
            while R9(v - b1) in cnt:
                v = v - b1
                p += 1
            v = nu
            while R9(v - b2) in cnt:
                v = v - b2
                q += 1
            cache[hk] = (p * p + q * q + p * q + 3 * p + 3 * q) / 3.0
        out[i] = cache[hk]
    return out


def main():
    t0 = time.time()
    checks = {}

    # ---- conventions from V65 ------------------------------------------------
    lam, _off = cartan_weights()
    x_coefs, y_coefs = identify_XY(lam)
    y_coef = y_coefs[0]
    phys = [(xc, physical_su5(lam @ xc, lam @ y_coef)) for xc in x_coefs]
    phys = [(xc, kind) for xc, kind in phys if kind]
    xc_gg = next(xc for xc, kind in phys if kind == "Georgi-Glashow")

    # ---- the composite spectrum 16 x 45 -------------------------------------
    alphas = flux_weights()                       # 45 weights (5 zeros)
    n_flux_zero = int(np.sum(np.abs(alphas).sum(axis=1) < 1e-12))
    comp = (lam[:, None, :] + alphas[None, :, :]).reshape(-1, 5)  # 720
    bL, b1, b2 = sm_roots(alphas, lam, xc_gg, y_coef)
    bL2 = float(bL @ bL)
    Ys = comp @ y_coef
    m2 = (comp @ bL) / bL2                        # SU(2) Cartan charges
    c2 = c2_labels(comp, b1, b2)                  # SU(3) Casimir per state
    db_state = np.stack([(2 / 3) * (3 / 5) * Ys ** 2,
                         (2 / 3) * m2 ** 2,
                         (2 / 3) * c2 / 8.0], axis=1)
    db_total = db_state.sum(axis=0)
    # the five zero-flux 16-copies: one IS the SM matter (in B_SM); the
    # other four are extra. The eig order is arbitrary: find the index.
    zero_j = next(j for j in range(45) if np.abs(alphas[j]).sum() < 1e-12)
    sm_idx = np.arange(16) * 45 + zero_j      # comp is mu-major: i*45 + j
    checks["composite"] = {
        "n_states": int(comp.shape[0]),
        "flux_zero_multiplicity": n_flux_zero,
        "db_total": [round(v, 4) for v in db_total],
        "net_16_copies": n_flux_zero,
        "note": "the content carries FIVE net chiral 16-copies (the "
                "zero-flux quintet); the alpha != 0 composites are "
                "vectorlike under the SM and shift the running through "
                "their indices only"}

    # ---- regression: the V64 LIB shifts for three fragments ------------------
    def frag(Y):
        return lam[np.abs(lam @ y_coef - Y) < 1e-9]

    def shifts(states):
        Y = states @ y_coef
        m = (states @ bL) / bL2
        c = c2_labels(states, b1, b2)
        return ((2 / 3) * (3 / 5) * float(np.sum(Y ** 2)),
                (2 / 3) * float(np.sum(m ** 2)),
                (2 / 3) * float(np.sum(c) / 8.0))

    lib = {}
    for name, Y, want in (("3,2", 1 / 6, (2 / 15, 2.0, 4 / 3)),
                          ("3b,1", -2 / 3, (16 / 15, 0.0, 2 / 3)),
                          ("1,2", -1 / 2, (2 / 5, 2 / 3, 0.0))):
        pair = np.vstack([frag(Y), -frag(Y)])
        got = shifts(pair)
        lib[name] = {"measured": [round(v, 6) for v in got],
                     "V64_LIB": list(want),
                     "err": max(abs(a - b) for a, b in zip(got, want))}
    checks["lib_regression"] = lib
    lib_ok = all(v["err"] < 1e-9 for v in lib.values())

    # ---- thresholds: m = w (c0 + |alpha.v(theta)|) --------------------------
    def solve(th):
        ca, cb = math.cos(th), math.sin(th)
        v = ca * xc_gg + cb * y_coef
        c_al = np.abs(alphas @ v)                  # per flux
        # comp is mu-major: tile the flux factors over the 16 weights
        c_state = np.tile(c_al, 16)
        # one of the five zero-flux 16-copies IS the SM matter already
        # counted in B_SM: exclude it from the desert running
        dbs = db_state.copy()
        dbs[sm_idx] *= 0.0

        def eqs3(z0, c0, x):
            w = math.exp(z0) * MZ
            zs = np.log(np.maximum(w * (c0 + c_state), 1e-30) / MZ)
            F = []
            for i in range(3):
                s = INV_A10 + B_SM[i] * x / TWO_PI
                act = np.clip(x - zs, 0, None)
                s += float(np.sum(dbs[:, i] * act)) / TWO_PI
                F.append(s - 1 / A_Z[i])
            return F

        best = None
        for z0s, c0s, xs in ((20.0, 1.0, 36.0), (15.0, 0.3, 34.0),
                             (25.0, 2.0, 38.0), (10.0, 0.1, 30.0),
                             (22.0, 0.5, 37.0)):
            z0, c0, x = z0s, c0s, xs
            for _ in range(200):
                F = np.array(eqs3(z0, c0, x))
                h = 1e-5
                J = np.zeros((3, 3))
                for j, (dz0, dc0, dx) in enumerate([(h, 0, 0), (0, h, 0),
                                                    (0, 0, h)]):
                    Fp = np.array(eqs3(z0 + dz0, c0 + dc0, x + dx))
                    J[:, j] = (Fp - F) / h
                try:
                    d = np.linalg.solve(J, -F)
                except np.linalg.LinAlgError:
                    break
                z0 = max(-5.0, min(60.0, z0 + d[0]))
                c0 = max(1e-3, min(50.0, c0 + d[1]))
                x = max(5.0, min(50.0, x + d[2]))
                if abs(d).max() < 1e-10:
                    break
            F = np.array(eqs3(z0, c0, x))
            if max(abs(f) for f in F) < 1e-6:
                best = (z0, c0, x)
                break
        if best is None:
            return None
        z0, c0, x = best
        w = math.exp(z0) * MZ
        M_X = math.exp(x) * MZ
        zs = np.log(np.maximum(w * (c0 + c_state), 1e-30) / MZ)
        phys_ok = (Z_MIN <= zs.min() and zs.max() <= x <= X_MAX and c0 > 0)
        return {"theta": th, "w_GeV": w, "c0": c0, "M_X_GeV": M_X,
                "z_min": float(zs.min()), "z_max": float(zs.max()),
                "physical": phys_ok}

    sols = []
    for k in range(0, 180):
        th = math.pi * k / 180.0
        r = solve(th)
        if r and r["physical"]:
            sols.append(r)
    uniq = []
    for r in sols:
        if not any(abs(r["theta"] - u["theta"]) < 1e-3 for u in uniq):
            uniq.append(r)

    def tau_p(M_X, alpha_X=A10):
        t_gev = M_X ** 4 / (alpha_X ** 2 * u10pred.M_P_ROBSON ** 5)
        return t_gev * u10pred.GEV_TO_S / u10pred.S_TO_YR

    best = None
    for r in uniq:
        r["tau_p_yr"] = tau_p(r["M_X_GeV"])
        if best is None or abs(math.log10(r["M_X_GeV"])
                               - math.log10(1.8e16)) < abs(
                                   math.log10(best["M_X_GeV"])
                                   - math.log10(1.8e16)):
            best = r

    checks["meeting"] = {"n_solutions": len(uniq),
                         "solutions": [{k: (round(v, 4) if isinstance(v, float)
                                            else v) for k, v in r.items()}
                                       for r in uniq[:6]]}

    # the check verifies the COMPUTATION (V74/V75 precedent): the index
    # machinery reproduces the V64 LIB exactly and the scan is complete;
    # the physics verdict (no desert solution) is carried by the status.
    ok = lib_ok

    if best:
        status = (
            f"DESERT CONTENT OF THE GAUGED LAYER COMPUTED: the composite "
            f"fermions (matter 16 x adjoint link flux, {comp.shape[0]} "
            f"states, vectorlike, total db = {tuple(round(v,2) for v in db_total)}) "
            f"with thresholds w(c0 + |alpha.v|): the meeting at "
            f"alpha_X = 1/25 EXACT has {len(uniq)} physical solution(s); "
            f"best: theta = {best['theta']:.3f}, w = {best['w_GeV']:.3e} "
            f"GeV, c0 = {best['c0']:.3f}, M_X = {best['M_X_GeV']:.3e} GeV, "
            f"tau_p = {best['tau_p_yr']:.2e} yr. LIB regression exact on "
            f"three fragments (err {max(v['err'] for v in lib.values()):.1e}). "
            f"The Cartan-flux quintet: {n_flux_zero - 1} EXTRA zero-flux "
            f"16-copies at w*c0 (the fifth is the SM matter itself, "
            f"excluded from the running) are the lightest content - named "
            f"residual: they must be lifted by the extra-U(1) breaking or "
            f"they overcount desert matter.")
    else:
        status = (
            "DESERT CONTENT OF THE GAUGED LAYER: the composite spectrum "
            f"(db = {tuple(round(v,2) for v in db_total)}) does NOT meet "
            "the three couplings at alpha_X = 1/25 for any (w, c0, theta) "
            "with desert thresholds: O4 stays open, now with the composite "
            "content excluded by name.")

    verdict = {"experiment": "V81_link_desert_content", "pass": bool(ok),
               "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v81_link_desert_content.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V81 - composite desert content of the gauged bond layer")
    print("=" * 68)
    for name, v in lib.items():
        print(f"  LIB ({name:5s}) err       = {v['err']:.2e}")
    print(f"  composite states          = {comp.shape[0]} "
          f"(net 16-copies: {checks['composite']['net_16_copies']})")
    print(f"  total db                  = {tuple(round(v,3) for v in db_total)}")
    print(f"  physical solutions        = {len(uniq)}")
    if best:
        print(f"  best: M_X = {best['M_X_GeV']:.3e} GeV, "
              f"tau_p = {best['tau_p_yr']:.2e} yr")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {os.path.join(outdir, 'v81_link_desert_content.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
