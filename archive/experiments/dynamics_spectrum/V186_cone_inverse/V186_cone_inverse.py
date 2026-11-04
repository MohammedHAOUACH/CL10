#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V186 - the cone inverse problem: search the entangling-bond space for a
Clifford-type tangent algebra at a gapless cubic multiplet.

**Why.**  V185 proved the separation lemma: for the implemented step every
tangent operator is (sub-site) x (internal identity), so no internal choice
can move the cone - the bond must entangle the two layers. V120 (one Peierls
phase) and V121 (one su(2) triplet, one strength per axis) tried two POINTS
of that space and closed nothing. This experiment treats the whole
axis-covariant entangling family as an inverse problem:

    hop along axis j carries  exp(+/- i phi_j D_j),
    D_j any even Hermitian involution of the 16 sector,
    (D_x, D_y, D_z) any triple from a structured pool,
    phi_j any strengths, with (theta, w) of the Wilson factor retuned per bond.

and asks for a point where, at the cubic corner R = (pi/2,pi/2,pi/2), the
exactly 3-fold multiplet is (i) gapless, (ii) carries a non-commuting tangent
algebra (comm_ratio > 0.5: a cone needs the S_j to close a Clifford/spin
algebra, V118), and (iii) has direction-independent speeds (anisotropy < 0.1).

**The correction this experiment also establishes.**  V121 reported that the
three-axis bond "breaks the degree-1 Fourier support that proves Theorem 5
(leakage 0.27)". Re-measured here with the Wilson factor OFF, the
three-axis bond is EXACTLY radius-1 (leakage 8.7e-17): the 0.27 leakage is the
Wilson Bessel tail, present identically for the implemented bond. The
entangling bond does NOT cost the causal radius - the price V121 named was
misattributed to the bond.

**Speed (GPU).**  M(q) = W(q) . sum_sigma c_sigma(q) P_sigma with 27 constant
128x128 matrices per bond (V121's Fourier form). The (theta, w) retune grid
is one batched tensordot + batched eigvals on the engine backend (cupy when
the card is free); the tangent metrics run on host (128x128, cheap).

Verdict: PASS = the protocol is valid (structure, regression, unitarity,
radius-1 of the bond, baseline reproduces V118's three planes) AND the scan
reports its best cone candidate either way; the physics verdict (cone found
/ no-go floor) is in the status with numbers.
"""
import itertools
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

# The cone programme runs on the CPU core ALWAYS, whatever the caller exported.
# Measured on this host: the scan is a batch of independent 128x128 Floquet
# eigendecompositions, ~40x faster on 12 OpenBLAS-single-threaded CPU workers
# than on a 3090 (cupy has no batched eigvals: the "batch" is a serial loop of
# ~50 ms GPU round-trips). The GPU route also needs a CUDA context per pool
# worker, and the scan then ran for hours without finishing (measured twice).
# Pinned BEFORE any cl10 import so the engine never selects a GPU backend here,
# and the pin is inherited by the pool workers.
os.environ["CL10_BACKEND"] = "numpy"

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# reference spectrum: the g = 1/5 Wilson-family member (paper 8.4):
# the entangling-bond search space (Propositions 4-5: gapless region
# Abelian to 1e-11, non-commuting region gapped by >= 8e-3) is the
# dispersive spectrum; the clock pinned pair is searched by
# V192/V203 instead.
os.environ["CL10_ANGLE"] = "0.2"
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))

import V65_wilson as V65                      # noqa: E402
from tqdm import tqdm                         # noqa: E402
from scipy.optimize import minimize           # noqa: E402
from cl10 import clifford as CL               # noqa: E402
from cl10 import config as cl10_config        # noqa: E402
from cl10 import supercell as SC              # noqa: E402
from cl10 import engine as E                  # noqa: E402
from cl10 import parpool                      # noqa: E402

_B16, _U32, _U16 = SC.sectors()
B16, U16 = E.asnumpy(_B16), E.asnumpy(_U16)
# raw backend module for batched linear algebra (the V102 pattern): cupy
# when the card is live, numpy otherwise.
_XP = E._E if E._BACKEND == "cupy" else np
I16 = np.eye(16, dtype=complex)
G10, ORDER = SC.G10, SC.STAGES
C = (math.pi / 2,) * 3                        # the cubic corner R
MNODE = (math.pi / 2, math.pi / 2, 0.0)       # an M doubler

_lam, _off = V65.cartan_weights()
_v65 = json.load(open(exp("V65_wilson",
                          "results", "v65_wilson.json")))
_gg = next(r for r in _v65["data"]["part_B"].values()
           if r["su5"] == "Georgi-Glashow")
XCOEF = np.array(_gg["X_coef"])
YCOEF = np.array(_v65["data"]["Y_coef"])
DX16 = _lam @ XCOEF                           # Wilson Cartan weights (16)
DY16 = _lam @ YCOEF

_h = [E.asnumpy(m) for m in CL.generators()]
T16 = {(a, b): B16.conj().T @ (1j * (_h[a] @ _h[b])) @ B16
       for a in range(10) for b in range(a + 1, 10)}


def Tg(a, b):
    """T_ab on the 16, any order (T_ba = -T_ab)."""
    return T16[(a, b)] if a < b else -T16[(b, a)]


def Rrot(phi, sign, D):
    """exp(i sign phi D) - exact because D^2 = I, D Hermitian."""
    return math.cos(phi) * I16 + 1j * sign * math.sin(phi) * D


def to_full(A):
    return A.transpose(0, 2, 1, 3).reshape(128, 128)


class Step:
    """M(q) = W(q) . sum_sigma c_sigma(q) P_sigma for one entangling bond.

    bond = ((D_x, phi_x), (D_y, phi_y), (D_z, phi_z)): the hop along axis j
    carries exp(+/- i phi_j D_j) with the sign of the hop sense. phi = 0 is
    exactly the implemented step."""

    def __init__(self, bond):
        Ui = to_full(np.array([[U16 if e == f else np.zeros((16, 16))
                                for f in range(8)] for e in range(8)])
                     .reshape(8, 8, 16, 16))
        stages = []
        for even, j in ORDER:
            D, phi = bond[j]
            c, s = math.cos(G10), math.sin(G10)
            diag = np.zeros((8, 8, 16, 16), dtype=complex)
            plus = np.zeros((8, 8, 16, 16), dtype=complex)
            minus = np.zeros((8, 8, 16, 16), dtype=complex)
            bit = 1 << j
            for e in range(8):
                diag[e, e] = c * I16
            for e in range(8):
                if (e & bit) == (0 if even else bit):
                    o = e ^ bit
                    sign = 1 if (e & bit) == 0 else -1
                    hop = -1j * s * Rrot(phi, sign, D)
                    hopr = -1j * s * Rrot(phi, -sign, D)
                    if even:
                        diag[e, o] = hop
                        diag[o, e] = hopr
                    else:
                        plus[e, o] = hop
                        minus[o, e] = hopr
            if even:
                stages.append([(None, 0, to_full(diag) @ Ui)])
            else:
                stages.append([(None, 0, to_full(diag) @ Ui),
                               (j, 1, to_full(plus) @ Ui),
                               (j, -1, to_full(minus) @ Ui)])
        terms = [({}, np.eye(128, dtype=complex))]
        for st in stages:
            new = []
            for t in terms:
                for axis, sgn, mat in st:
                    sig = dict(t[0])
                    if axis is not None:
                        sig[axis] = sig.get(axis, 0) + sgn
                    new.append((sig, mat @ t[1]))
            terms = new
        self.sig = np.array([[t[0].get(0, 0), t[0].get(1, 0), t[0].get(2, 0)]
                             for t in terms])
        self.P = np.array([t[1] for t in terms])
        self.Px = _XP.asarray(self.P)         # raw-backend copy for batches

    def M(self, q, w, th):
        ph = w * sum(1.0 - math.cos(qj) for qj in q)
        d = np.cos(th) * DX16 + np.sin(th) * DY16
        wd = np.exp(-1j * ph * np.tile(d, 8))
        coeffs = np.exp(1j * (self.sig @ np.asarray(q)))
        X = np.tensordot(coeffs, self.P, axes=([0], [0]))
        return wd[:, None] * X

    def batch_topgap(self, ws, ths, q=C, m=3):
        """Batched top-of-multiplet gap over the (w, theta) grid, on the
        raw backend (GPU when cupy is live)."""
        qv = np.asarray(q)
        # coefficients are tiny (one per Fourier term): compute on the host,
        # then move to the backend (cupy refuses numpy operands in @).
        coeffs = _XP.asarray(np.exp(1j * (self.sig @ qv)))
        X = _XP.tensordot(coeffs, self.Px, axes=([0], [0]))
        ws = np.asarray(ws, dtype=float)
        ths = np.asarray(ths, dtype=float)
        d = (np.cos(ths)[:, None] * DX16[None, :]
             + np.sin(ths)[:, None] * DY16[None, :])       # (B,16)
        ph = ws * sum(1.0 - math.cos(qj) for qj in q)      # (B,)
        wd = _XP.asarray(np.exp(-1j * ph[:, None]
                                * np.tile(d, (1, 8))))     # (B,128), tile
        # (NOT repeat: the internal index is the FAST one, matching
        #  kron(I_8, diag(d)) in cl10.supercell and Step.M above)
        Xb = X[None, :, :] * wd[:, :, None]
        try:
            ev = _XP.linalg.eigvals(Xb)
        except Exception:                                  # stacked fallback
            ev = _XP.stack([_XP.linalg.eigvals(Xb[i])
                            for i in range(Xb.shape[0])])
        a = _XP.sort(_XP.abs(_XP.angle(ev)), axis=1)
        return E.asnumpy(a[:, m - 1])


def tangent(st, q, m, h=1e-5, w=0.0, th=0.0):
    """S_j = i P^dag (d_j M) P on the m-dim near-gapless subspace, plus the
    V118/V121 metrics: comm_ratio and the Clifford defect."""
    ev, V = np.linalg.eig(st.M(q, w, th))
    sel = np.argsort(np.abs(np.angle(ev)))[:m]
    Q, _ = np.linalg.qr(V[:, sel])
    S = []
    for j in range(3):
        qp, qm = list(q), list(q)
        qp[j] += h
        qm[j] -= h
        D = (st.M(tuple(qp), w, th) - st.M(tuple(qm), w, th)) / (2 * h)
        S.append(1j * (Q.conj().T @ D @ Q))
    nrm = [float(np.linalg.norm(s)) for s in S]
    mean2 = sum(n * n for n in nrm) / 3.0
    comm = max(float(np.linalg.norm(S[a] @ S[b] - S[b] @ S[a]))
               for a, b, _ in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
    comm_ratio = comm / mean2 if mean2 > 0 else 0.0
    c2 = sum(float(np.linalg.norm(s @ s)) for s in S) / 3.0
    defect = 0.0
    if c2 > 0:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            defect = max(defect, float(np.linalg.norm(
                S[a] @ S[b] + S[b] @ S[a])) / (2 * c2))
        for s in S:
            defect = max(defect, abs(float(np.linalg.norm(s @ s)) - c2) / c2)
    return S, comm_ratio, defect


def sphere(n=6):
    out = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                r = math.sqrt(i * i + j * j + k * k)
                if r > 0:
                    out.append((i / r, j / r, k / r))
    return out


def speeds(S, dirs):
    spd = []
    for kd in dirs:
        Ax = kd[0] * S[0] + kd[1] * S[1] + kd[2] * S[2]
        spd.extend(sorted(abs(x) for x in np.linalg.eigvalsh(Ax)))
    spd = np.array(spd)
    nz = spd[spd > 1e-6]
    if not nz.size:
        return None, 0.0, 0.0, int((spd < 1e-6).sum()), len(spd)
    return (float((nz.max() - nz.min()) / (nz.max() + nz.min())),
            float(nz.min()), float(nz.max()),
            int((spd < 1e-6).sum()), len(spd))


def fourier_support(bond, w=0.0):
    """Degree-1 Fourier support of M(q) along each axis (5-point DFT): the
    causal-radius-1 statement. Tested with the Wilson factor OFF (w=0) and
    ON, so the bond's locality and the Wilson tail are separated."""
    st = Step(bond)
    err = 0.0
    for j in range(3):
        mats = []
        for k in range(5):
            q = [0.61, -1.13, 2.27]
            q[j] = 2 * math.pi * k / 5
            mats.append(st.M(tuple(q), w, 0.4))
        for n in (2, 3):
            F = sum(mats[k] * np.exp(-2j * math.pi * n * k / 5)
                    for k in range(5)) / 5
            err = max(err, float(np.abs(F).max()))
    return err


def tune(st, q, m, grid=18, wmax=3.0, rounds=6, nmax=300):
    """Golden-section refinement of (w, theta) minimising the top level of
    the m-dimensional multiplet at fixed bond (V121's tune)."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0

    def top(wv, thv):
        return float(st.batch_topgap([wv], [thv], q, m)[0])

    best = (9.0, None)
    for it in range(grid):
        th = 2 * math.pi * it / grid
        for w in np.linspace(0.0, wmax, 19):
            v = top(float(w), th)
            if v < best[0]:
                best = (v, (float(w), th))
    w, th = best[1]
    for _ in range(rounds):
        for var in ("w", "th"):
            lo, hi = (w - 0.1, w + 0.1) if var == "w" else (th - 0.2, th + 0.2)
            c1, c2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
            f1 = top(c1, th) if var == "w" else top(w, c1)
            f2 = top(c2, th) if var == "w" else top(w, c2)
            for _ in range(25):
                if f1 < f2:
                    hi, c2, f2 = c2, c1, f1
                    c1 = hi - gr * (hi - lo)
                    f1 = top(c1, th) if var == "w" else top(w, c1)
                else:
                    lo, c1, f1 = c1, c2, f2
                    c2 = lo + gr * (hi - lo)
                    f2 = top(c2, th) if var == "w" else top(w, c2)
            if var == "w":
                w = 0.5 * (lo + hi)
            else:
                th = 0.5 * (lo + hi)
    # Nelder-Mead polish on (w, theta): the gapless corner family is a curve
    # in the plane and the golden-section windows are 1-D - V118's 1e-10
    # points need the 2-D polish.
    res = minimize(lambda x: top(float(x[0]), float(x[1])), [w, th],
                   method="Nelder-Mead",
                   options={"maxfev": nmax, "fatol": 1e-13, "xatol": 1e-10})
    if res.fun < top(w, th):
        w, th = float(res.x[0]), float(res.x[1])
    return w, th, top(w, th)


def measure(st, q, m, dirs, w=None, th=None):
    if w is None or th is None:
        w2, th2, gtop = tune(st, q, m)
    else:
        w2, th2 = float(w), float(th)
        gtop = float(st.batch_topgap([w2], [th2], q, m)[0])
    S, cr, defect = tangent(st, q, m, 1e-5, w2, th2)
    an, vmin, vmax, nz0, nd = speeds(S, dirs)
    a = np.sort(np.abs(np.angle(np.linalg.eigvals(st.M(q, w2, th2)))))
    return {"w": w2, "theta": th2, "gap": float(a[0]), "top": gtop,
            "spread": float(a[m - 1] - a[0]), "comm_ratio": cr,
            "clifford_defect": defect, "anisotropy": an,
            "v_min": vmin, "v_max": vmax}


# ---------------------------------------------------------------------------
# the bond pool
# ---------------------------------------------------------------------------

def bond_pool(seed=186, n_random=40, n_disjoint=60):
    """Structured triples (D_x, D_y, D_z) of even Hermitian involutions.

    (a) su(2) triples: for each index triple (a,b,c), the cyclic planes
        (T_ab, T_bc, T_ca) - [Sigma_a,Sigma_b] = -2i eps Sigma_c exactly.
        V121's SIG is the (0,1,2) member.
    (b) disjoint commuting triples (T_ab, T_cd, T_ef), sampled.
    (c) random even involutions: D = sign(H), H a random real combination
        of the 45 T_ab - the general family, Monte-Carlo."""
    rng = np.random.default_rng(seed)
    pool = []
    for (a, b, c) in itertools.combinations(range(10), 3):
        pool.append((("su2", (a, b, c)),
                     (Tg(a, b), Tg(b, c), Tg(c, a))))
    idx6 = list(itertools.combinations(range(10), 6))
    for six in rng.permutation(idx6)[:n_disjoint]:
        (a, b, c, d, e, f) = six
        pool.append((("disj", six),
                     (T16[(a, b)], T16[(c, d)], T16[(e, f)])))
    keys = list(T16.keys())
    for _ in range(n_random):
        Ds = []
        for _ax in range(3):
            H = np.zeros((16, 16), dtype=complex)
            for coef in rng.normal(0, 1, 12):
                H += coef * T16[keys[rng.integers(len(keys))]]
            H = 0.5 * (H + H.conj().T)
            ev, V = np.linalg.eigh(H)
            Ds.append(V @ np.diag(np.sign(ev)) @ V.conj().T)
        pool.append((("rand", int(rng.integers(1 << 30))), tuple(Ds)))
    return pool


def _scan_worker(item):
    """One bond of the (w,theta) grid scan (module-level for pickling)."""
    (tag, Ds, pat, sgn, ws, ths) = item
    bond = tuple((Ds[j], pat * sgn[j]) for j in range(3))
    st = Step(bond)
    g = st.batch_topgap(np.repeat(ws, len(ths)),
                        np.tile(ths, len(ws)), C, 3)
    k = int(np.argmin(g))
    return (float(g[k]), tag, bond, float(ws[k // len(ths)]),
            float(ths[k % len(ths)]))


def _refine_worker(item):
    """Golden-section + Nelder-Mead retune and cone-metric measurement of
    one top-scan candidate (module-level for pickling)."""
    (tag, bond, dirs) = item
    st = Step(bond)
    m = measure(st, C, 3, dirs)
    return (m, tag, bond)


def main():
    t0 = time.time()
    out = {"experiment": "V186_cone_inverse", "data": {}}
    checks = []
    dirs = sphere(6)
    ZERO = ((np.eye(16, dtype=complex), 0.0),) * 3

    # ---- T0: structure, regression, unitarity, and the RADIUS-1
    #          correction: the bond is exactly radius-1 with the Wilson
    #          factor OFF; the 0.27 leakage of V121 is the Wilson tail. ----
    st0 = Step(ZERO)
    reg = float(np.abs(st0.M(C, 0.0, 0.0) - SC.bloch(C, U16)).max())
    rng = np.random.default_rng(186)
    unit = 0.0
    pool_probe = bond_pool(seed=7, n_random=3, n_disjoint=1)
    for _ in range(8):
        tag, Ds = pool_probe[rng.integers(len(pool_probe))]
        phi = tuple(float(x) for x in rng.uniform(-1.5, 1.5, 3))
        st = Step(tuple(zip(Ds, phi)))
        q = tuple(float(x) for x in rng.uniform(-math.pi, math.pi, 3))
        A = st.M(q, 0.0, 0.0)
        unit = max(unit, float(np.abs(A @ A.conj().T - np.eye(128)).max()))
    lr_bond = fourier_support(((Tg(0, 1), 0.5), (Tg(1, 2), 0.9),
                               (Tg(2, 0), 1.3)), w=0.0)
    lr_bond_w = fourier_support(((Tg(0, 1), 0.5), (Tg(1, 2), 0.9),
                                 (Tg(2, 0), 1.3)), w=0.7)
    lr_impl_w = fourier_support(ZERO, w=0.7)
    checks.append(("T0 phi=0 reproduces the implemented step (%.1e); the "
                   "entangling step is exactly unitary (%.1e); and the "
                   "three-axis bond is EXACTLY radius-1 with the Wilson "
                   "factor OFF (degree-1 leakage %.1e) - the 0.27 leakage "
                   "V121 attributed to the bond is the Wilson Bessel tail, "
                   "measured %.1e for the bond and %.1e for the IMPLEMENTED "
                   "bond at the same w: the entangling bond does NOT cost "
                   "the causal radius"
                   % (reg, unit, lr_bond, lr_bond_w, lr_impl_w),
                   reg < 1e-12 and unit < 1e-12 and lr_bond < 1e-12
                   and lr_bond_w > 1e-3 and abs(lr_bond_w - lr_impl_w) < 0.05))
    out["data"]["structure"] = {
        "phi0_regression": reg, "unitarity": unit,
        "bond_radius1_wilson_off": lr_bond,
        "bond_radius1_wilson_on": lr_bond_w,
        "implemented_radius1_wilson_on": lr_impl_w}

    # ---- T1: baseline - reproduce V118 at the corner: gapless triplet,
    #          commuting tangent, three nodal planes. ----
    st_imp = Step(ZERO)
    base = measure(st_imp, C, 3, dirs)
    checks.append(("T1 baseline at R: the implemented step's corner triplet "
                   "is driven gapless by the Wilson factor (top %.1e, "
                   "spread %.1e) with a COMMUTING tangent algebra "
                   "(comm_ratio %.3f < 0.5) and anisotropy %s - V118's three "
                   "nodal planes reproduced"
                   % (base["top"], base["spread"], base["comm_ratio"],
                      ("%.4f" % base["anisotropy"])
                      if base["anisotropy"] is not None else "n/a"),
                   base["top"] < 1e-8 and base["spread"] < 1e-8
                   and base["comm_ratio"] < 0.5))
    out["data"]["baseline_R"] = base

    # ---- T2: the pool scan. For every bond triple and strength pattern,
    #          the (theta,w) retune grid is one batched GPU eigval call. --
    smoke = cl10_config.smoke()               # CL10_SMOKE, via cl10/config.py
    pool = bond_pool(n_random=(2 if smoke else 40),
                     n_disjoint=(3 if smoke else 60))
    if smoke:
        pool = [p for p in pool if p[0][0] == "su2"][:8]
    patterns = [(p, s) for p in ((0.5,) if smoke else (0.3, 0.9))
                for s in (((1, 1, 1),) if smoke else
                          ((1, 1, 1), (1, 1, -1)))]
    ws = np.linspace(0.0, 3.0, 9 if smoke else 19)
    ths = np.linspace(0.0, 2 * math.pi, 8 if smoke else 24)
    nscan = len(pool) * len(patterns)
    print("  scan: %d bonds (grid %dx%d, %d workers)"
          % (nscan, len(ws), len(ths), parpool.njobs_default()), flush=True)
    bar = tqdm(total=nscan, desc="V186 scan", unit="bond",
               mininterval=30, file=sys.stdout)
    items = [(tag, Ds, pat, sgn, ws, ths)
             for tag, Ds in pool for pat, sgn in patterns]
    scan = parpool.pmap(_scan_worker, items,
                        progress=lambda d, n: bar.update(1))
    bar.close()
    scan.sort(key=lambda r: r[0])
    out["data"]["scan_n"] = len(scan)
    out["data"]["scan_best10"] = [
        {"gap_top": r[0], "family": r[1][0], "idx": str(r[1][1]),
         "phi": [b[1] for b in r[2]]} for r in scan[:10]]
    print("  scan: %d bonds, best top-gap %.1e (%s %s phi=%s)"
          % (len(scan), scan[0][0], scan[0][1][0], scan[0][1][1],
             [round(b[1], 2) for b in scan[0][2]]), flush=True)

    # ---- T3: refine the top candidates (golden-section (w,theta)), then
    #          measure the cone metrics. ----
    top = scan[:3 if smoke else 20]
    bar = tqdm(total=len(top), desc="V186 refine", unit="bond",
               mininterval=30, file=sys.stdout)
    refined = parpool.pmap(_refine_worker,
                            [(tag, bond, dirs)
                             for _, tag, bond, _, _ in top],
                            progress=lambda d, n: bar.update(1))
    bar.close()
    refined.sort(key=lambda r: r[0]["top"] + r[0]["spread"])
    out["data"]["refined"] = [
        {"family": t[0], "idx": str(t[1]), "phi": [b[1] for b in bd], **mm}
        for mm, t, bd in refined[:10]]

    # ---- T4: continuous optimisation of the best candidate: Nelder-Mead
    #          over (phi_x, phi_y, phi_z, w, theta) minimising top+spread. --
    from scipy.optimize import minimize
    mm, tag, bond = refined[0]
    Ds = tuple(b[0] for b in bond)

    def loss(x):
        st = Step(tuple((Ds[j], float(x[j])) for j in range(3)))
        return float(st.batch_topgap([float(x[3])], [float(x[4])], C, 3)[0])

    x0 = np.array([bond[0][1], bond[1][1], bond[2][1], mm["w"], mm["theta"]])
    bar = tqdm(total=80 if smoke else 400, desc="V186 optimise",
               unit="eval", mininterval=30, file=sys.stdout)

    def _cb(*_a):
        bar.update(1)

    res = minimize(loss, x0, method="Nelder-Mead", callback=_cb,
                   options={"maxfev": 80 if smoke else 400,
                            "xatol": 1e-6, "fatol": 1e-10})
    bar.close()
    st_opt = Step(tuple((Ds[j], float(res.x[j])) for j in range(3)))
    opt = measure(st_opt, C, 3, dirs, float(res.x[3]), float(res.x[4]))
    opt["phi"] = [float(res.x[j]) for j in range(3)]
    opt["family"] = tag[0]
    opt["idx"] = str(tag[1])
    out["data"]["optimised"] = opt
    print("  optimised: top %.1e spread %.1e comm %.3f A %s"
          % (opt["top"], opt["spread"], opt["comm_ratio"],
             ("%.4f" % opt["anisotropy"]) if opt["anisotropy"] is not None
             else "n/a"), flush=True)

    # ---- T5: the cone criterion over everything measured. ----
    allm = [base] + [r[0] for r in refined] + [opt]
    cone = [m for m in allm
            if m["top"] < 1e-8 and m["comm_ratio"] > 0.5
            and m["anisotropy"] is not None and m["anisotropy"] < 0.1]
    found = bool(cone)
    gmin = min(m["top"] + m["spread"] for m in allm)
    crmax = max(m["comm_ratio"] for m in allm)
    amin = min((m["anisotropy"] for m in allm if m["anisotropy"] is not None),
               default=1.0)
    checks.append(("T2-T4 the cone criterion (gapless triplet top < 1e-8 AND "
                   "comm_ratio > 0.5 AND anisotropy < 0.1) is %s over %d "
                   "scanned bonds + %d refined + 1 optimised: best gap+spread "
                   "%.1e, max comm_ratio %.3f, min anisotropy %.4f"
                   % ("MET - a cone exists" if found else "NOT met",
                      len(scan), len(refined), gmin, crmax, amin),
                   True))     # physics verdict either way; protocol validity
    out["data"]["cone_found"] = found

    # ---- T6: does the best entangling bond remove the M-node doubling? --
    mrows = []
    for lbl, bd in (("implemented", ZERO),
                    ("best", tuple((Ds[j], opt["phi"][j]) for j in range(3)))):
        st = Step(bd)
        w, th, gtop = tune(st, MNODE, 1)
        S, cr, defect = tangent(st, MNODE, 1, 1e-5, w, th)
        an, vmin, vmax, _, _ = speeds(S, dirs)
        mrows.append({"bond": lbl, "gap": gtop, "comm_ratio": cr,
                      "anisotropy": an})
    out["data"]["M_node"] = mrows
    checks.append(("T6 the M-node doubler survives the best entangling bond "
                   "(gap %.1e vs %.1e implemented): the bond programme does "
                   "not remove the doubling either"
                   % (mrows[1]["gap"], mrows[0]["gap"]),
                   mrows[1]["gap"] < 1e-4))

    ok = all(c[1] for c in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "CONE INVERSE PROBLEM: %s. The "
        "axis-covariant entangling family - hop along axis j carries "
        "exp(+/- i phi_j D_j), D_j from a pool of %d triples (all 120 su(2) "
        "plane triples, %d disjoint, %d random even involutions) x strength "
        "patterns, with (theta,w) retuned per bond - was scanned (%d bonds, "
        "batched GPU eigvals) and the best was continuously optimised. "
        "STRUCTURE CORRECTION: the entangling bond is exactly radius-1 with "
        "the Wilson factor OFF (leakage %.1e); V121's 0.27 leakage is the "
        "Wilson tail (implemented bond: %.1e at the same w) - the bond does "
        "NOT cost the causal radius. BASELINE: the implemented corner "
        "triplet is gapless (top %.1e) with commuting tangent (comm %.3f), "
        "V118 reproduced. RESULT: %s - best gap+spread %.1e, max "
        "comm_ratio %.3f (a cone needs > 0.5), min anisotropy %.4f (a cone "
        "needs < 0.1). The M-node doubling survives the best bond (gap "
        "%.1e). O9 %s."
        % ("CONE FOUND" if found else "NO CONE IN THE FAMILY",
           len(pool), 60, 40, len(scan), lr_bond, lr_impl_w,
           base["top"], base["comm_ratio"],
           "a Clifford-type tangent algebra was reached" if found
           else "no bond in the family reaches a Clifford-type tangent "
                "algebra at a gapless cubic multiplet",
           gmin, crmax, amin, mrows[1]["gap"],
           "CLOSED by construction" if found else "stays open: the "
           "obstruction is not the bond strength or the generator choice - "
           "it is the separation lemma (V185) plus the gaplessness rigidity "
           "measured here"))

    res_dir = os.path.join(HERE, "results")
    os.makedirs(res_dir, exist_ok=True)
    out["wall_time"] = time.time() - t0
    with open(os.path.join(res_dir, "v186_cone_inverse.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V186 - the cone inverse problem")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
