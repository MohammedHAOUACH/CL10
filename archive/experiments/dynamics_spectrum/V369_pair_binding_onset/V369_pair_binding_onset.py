#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V369 - pair-binding onset at the certified layer interaction
(booked BEFORE the run; interaction programme first cut). J.L identity
(T0), kernel shift -77 (T1), channel tower + V366/V368 anchors (T2),
18-eta Casimir slope (T3); channel projectors from the k=400 SM run;
true channel internal bands via compressed P_c^dag (U32 (x) U32) P_c;
eta scan of the composite-vs-nearest-channels gap over 12 booked K
plus K=0 and the anchor worst-K. Verdicts: BINDING-ONSET-<eta*> /
NO-ONSET-IN-SCAN / CASIMIR-SLOPE-DRIFT / ANCHOR-DRIFT /
TOWER-TRUNCATED.
"""
import json
import os
import sys

import cupy as cp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)
from cl10 import supercell as SC  # noqa: E402
from cl10.paths import exp  # noqa: E402

GRID_N = 6


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def modpi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def circmin(c, T):
    i = np.searchsorted(T, c)
    n = len(T)
    best = np.full(len(c), np.inf)
    for sh in (i - 1, i, i + 1):
        d = np.abs(modpi(c - T[sh % n]))
        best = np.minimum(best, d)
    return float(best.min())


def main():
    checks, data = [], {}

    def log(m):
        print(m, flush=True)

    Ts = root_generators()
    f, _k, _a = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    _Bm, U32, _U16 = SC.sectors(SC.G10)
    dim = 32 * 32 * 45

    Qsp = [(sp.kron(sp.kron(sp.csr_matrix(Ts[a]), sp.eye(32,
                                                         format="csr")),
                    sp.eye(45, format="csr"))
           + sp.kron(sp.kron(sp.eye(32, format="csr"),
                             sp.csr_matrix(Ts[a])),
                     sp.eye(45, format="csr"))
           + sp.kron(sp.kron(sp.eye(32, format="csr"),
                             sp.eye(32, format="csr")),
                     sp.csr_matrix(Ls[a]))).tocsr()
           for a in range(45)]

    def JLv(x):                       # J.L = sum_a Q_a Q_a  (same stack)
        return sum(q @ (q @ x) for q in Qsp)

    Tg = [cp.asarray(t) for t in Ts]
    Lg = [cp.asarray(L) for L in Ls]

    def c2_matvec(x):
        h = cp.asarray(x).reshape(32, 32, 45)
        y = cp.zeros_like(h)
        for a in range(45):
            z = (cp.einsum("ab,bgm->agm", Tg[a], h)
                 + cp.einsum("ab,gbm->gam", Tg[a], h)
                 + cp.einsum("ab,gmb->gma", Lg[a], h))
            y += (cp.einsum("ab,bgm->agm", Tg[a], z)
                  + cp.einsum("ab,gbm->gam", Tg[a], z)
                  + cp.einsum("ab,gmb->gma", Lg[a], z))
        return cp.asnumpy(y.ravel())

    # ------------------------------------------------- T0 Casimir identity
    # the identity compares 1/2 (C2tot - C21 - C22 - C2L) against the
    # explicit current-link contraction sum_a J_a L_a with
    # J_a = T(1)_a + T(2)_a:
    t0 = 0.0
    rng = np.random.default_rng(0x99)
    for _ in range(6):
        v = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)

        def T1a(a, vec):
            return np.einsum("ab,bgm->agm", Ts[a], vec)

        def T2a(a, vec):
            return np.einsum("ab,gbm->gam", Ts[a], vec)

        def La(a, vec):
            return np.einsum("ab,gmb->gma", Ls[a], vec)

        hv = v.reshape(32, 32, 45)
        jlv = np.zeros_like(hv)
        for a in range(45):
            j = T1a(a, hv) + T2a(a, hv)               # current J_a
            jlv += La(a, j)                           # J_a L_a
        c2v = JLv(v).reshape(32, 32, 45)
        b1 = np.zeros_like(hv)
        b2 = np.zeros_like(hv)
        b3 = np.zeros_like(hv)
        for a in range(45):
            b1 += T1a(a, T1a(a, hv))
            b2 += T2a(a, T2a(a, hv))
            b3 += La(a, La(a, hv))
        rhs = 0.5 * (c2v - b1 - b2 - b3)
        t0 = max(t0, float(np.max(np.abs(jlv - rhs))
                           / max(np.abs(jlv).max(), 1e-12)))
    checks.append(("T0 Casimir identity ||J.L - (C2tot-C21-C22-C2L)/2|| "
                   "rel <= 1e-9 on 6 random states: %.1e: %s"
                   % (t0, t0 < 1e-9), t0 < 1e-9))

    # ------------------------------------------------ kernel + T1 (-77)
    op = spla.LinearOperator((dim, dim), matvec=c2_matvec,
                             dtype=np.complex128)
    cache = os.path.join(HERE, "results", "v369_sm400.npz")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache):
        z = np.load(cache)
        w, V = z["w"], z["V"]
    else:
        w, v = spla.eigsh(op, k=400, which="SM")
        idx = np.argsort(w)
        w, V = w[idx], np.asarray(v)[:, idx]
        np.savez(cache, w=w, V=V)
    ker = np.abs(w) < 1e-9
    W = V[:, ker]
    Gm = W.conj().T @ W
    gv, GU = np.linalg.eigh(Gm)
    W = W @ (GU * gv ** -0.5) @ GU.conj().T
    nk = W.shape[1]
    JLW = np.column_stack([JLv(W[:, d]) for d in range(nk)])
    t1 = float(np.max(np.abs(JLW + 77.0 * W)))
    checks.append(("T1 kernel: dim %d (booked 2), ||(J.L + 77) W|| <= "
                   "1e-12: %.1e: %s" % (nk, t1, nk == 2 and t1 < 1e-12),
                   nk == 2 and t1 < 1e-12))

    # ------------------------------------------- T2 tower + V368 anchor
    nz = w[~ker]
    groups = []
    cur = [nz[0]]
    for x in nz[1:]:
        if x - cur[-1] < 1e-5:
            cur.append(x)
        else:
            groups.append(cur)
            cur = [x]
    groups.append(cur)
    chan = [(float(np.mean(g)), len(g)) for g in groups[:6]]
    log("channel tower (c, mult in first 400): %s"
        % [(round(c, 6), m) for c, m in chan])
    edge = float(w[-1])
    trunc = any(abs(edge - c) < 1e-5 for c, _ in chan[:3])
    mono = all(chan[i][0] < chan[i + 1][0] for i in range(len(chan) - 1))
    t2a = abs(chan[0][0] - 36.0) < 1e-6 and chan[0][1] >= 2

    v368 = json_art("V368_torus_rigidity_window")
    anchor_val = float(v368["data"]["gaps_min"])
    U6 = np.linalg.matrix_power(U32, 6)
    iota = np.angle(np.linalg.eigvals(U6))
    grid = [(2 * np.pi * i / GRID_N, 2 * np.pi * j / GRID_N,
             2 * np.pi * l / GRID_N)
            for i in range(GRID_N) for j in range(GRID_N)
            for l in range(GRID_N)]
    betas = np.array([np.angle(np.linalg.eigvals(SC.bond_product(q)))
                      for q in grid])
    v251 = json_art("V251_gauged_two_particle_sector")
    rows = {kk: float(vv) for kk, vv in
            v251["data"]["sensitivity"]["omega_rows"].items()}
    phys = [vv for vv in rows.values() if abs(vv) > 1e-9]
    tri = np.array([(a, b, c) for a in range(GRID_N)
                    for b in range(GRID_N) for c in range(GRID_N)])
    idxs = np.arange(216)
    lin = np.array([36, 6, 1])
    rngA = np.random.default_rng(0x51)
    g0 = []
    worstK = 0
    for K in range(216):
        k2 = ((tri - tri[K]) % GRID_N) @ lin
        comp = modpi((betas[idxs][:, :, None]
                      + betas[k2][:, None, :]).ravel()[:, None]
                     + np.array(phys)[None, :]).ravel()
        samp = rngA.choice(216, size=24, replace=False)
        s2 = ((tri[samp] - tri[K]) % GRID_N) @ lin
        sb = (betas[samp][:, :, None]
              + betas[s2][:, None, :]).ravel()
        cont = np.sort(modpi(sb[:, None]
                             + (iota[:, None]
                                + iota[None, :]).ravel()[None, :]
                             ).ravel())
        gg = circmin(comp, cont)
        g0.append(gg)
        if gg < g0[worstK]:
            worstK = K
    g0 = np.array(g0)
    anchor_ok = abs(g0.min() - anchor_val) < 1e-9
    checks.append(("T2 anchor: V368 eta=0 min gap %.6e reproduced to "
                   "1e-9: %.6e: %s" % (anchor_val, g0.min(), anchor_ok),
                   anchor_ok))
    checks.append(("T2 tower: head channel 36 mult %d >= 2, monotone %s, "
                   "not truncated %s: %s"
                   % (chan[0][1], mono, not trunc,
                      bool(t2a and mono and not trunc)),
                   bool(t2a and mono and not trunc)))

    # ---------------------------------- channel projectors + true bands
    chan3 = [c for c, _ in chan[:3]]
    Ph = {}
    bandgates = []
    U2r = np.kron(U32, U32).reshape(32, 32, 32, 32)
    for cval in chan3:
        sel = np.abs(w - cval) < 1e-5
        P = V[:, sel]
        # booked hygiene (V366 precedent): Lanczos returns NON-ORTHOGONAL
        # bases on degenerate clusters; Gram inverse-sqrt within the
        # cluster preserves the C2 eigenspace.
        Gc = P.conj().T @ P
        gv2, GU2 = np.linalg.eigh(Gc)
        assert gv2.min() > 1e-8, "channel Gram rank-deficient"
        P = P @ (GU2 * gv2 ** -0.5) @ GU2.conj().T
        opP = np.column_stack([op @ P[:, d] for d in range(P.shape[1])])
        resid = float(np.max(np.abs(opP - cval * P)))
        orth = float(np.max(np.abs(P.conj().T @ P
                                   - np.eye(P.shape[1]))))
        Pm = P.reshape(32, 32, 45, P.shape[1])
        UPp = np.einsum("abAB,ABmd->abmd", U2r, Pm)
        M = P.conj().T @ UPp.reshape(dim, P.shape[1])
        ph = np.sort(np.angle(np.linalg.eigvals(M)))
        Ph[cval] = ph
        bandgates.append(resid < 1e-8 and orth < 1e-10)
        log("channel %.4f: dim %d, C2 resid %.1e, orth %.1e, %d internal "
            "phases" % (cval, P.shape[1], resid, orth, len(ph)))
    checks.append(("channel projectors: C2 resid <= 1e-8 and P^dag P = I "
                   "<= 1e-10 on all three: %s" % all(bandgates),
                   all(bandgates)))

    # -------------------------------------------------- eta scan
    Ks = [int(x) for x in np.random.default_rng(0x51).choice(
        216, size=12, replace=False)]
    Ks = sorted(set(Ks + [0, worstK]))
    rngK = {K: np.random.default_rng(0x51 + K).choice(
        216, size=24, replace=False) for K in Ks}

    def gap_at(eta):
        gmin = np.inf
        for K in Ks:
            k2 = ((tri - tri[K]) % GRID_N) @ lin
            comp = modpi((betas[idxs][:, :, None]
                          + betas[k2][:, None, :]).ravel()[:, None]
                         + np.array(phys)[None, :]
                         - 77.0 * eta).ravel()
            sidx = rngK[K]
            s2 = ((tri[sidx] - tri[K]) % GRID_N) @ lin
            sb = (betas[sidx][:, :, None]
                  + betas[s2][:, None, :]).ravel()
            cont = np.sort(modpi(np.concatenate(
                [(sb[:, None] + Ph[cval][None, :]
                  + 0.5 * (cval - 154.0) * eta).ravel()
                 for cval in chan3])))
            gmin = min(gmin, circmin(comp, cont))
        return gmin

    etas = [0.05 * i for i in range(1, 21)]
    gaps = {e: gap_at(e) for e in etas}
    for e in sorted(gaps):
        log("eta %.2f: gap %.6f" % (e, gaps[e]))
    onset = [e for e in etas if gaps[e] > 1e-3]
    eta_star = None
    if onset:
        lo, hi = onset[0] - 0.05, onset[0]
        for _ in range(10):
            mid = 0.5 * (lo + hi)
            if gap_at(mid) > 1e-3:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-4:
                break
        eta_star = hi
    es = np.array([e for e in etas if e >= 0.55])
    gs = np.array([gaps[e] for e in es])
    slope = float(np.polyfit(es, gs, 1)[0]) if len(es) >= 4 else 0.0
    slope_ok = abs(slope - 18.0) < 0.5
    checks.append(("T3 18-eta Casimir slope: fit %.3f in 18+-0.5: %s"
                   % (slope, slope_ok), slope_ok))

    v101 = json_art("V101_bond_current_ground_state")
    e0r = float(v101["data"]["E0"])
    eta_anchor_ok = abs(e0r + 5802.398188196636) < 1e-6
    checks.append(("E-eta anchor: V101 published E0 = %.6f eta re-read: "
                   "%s" % (e0r, eta_anchor_ok), eta_anchor_ok))

    data["tower"] = [(c, m) for c, m in chan]
    data["gaps"] = {str(e): gaps[e] for e in etas}
    data["eta_star"] = eta_star
    data["slope"] = slope
    data["anchor_g0_min"] = float(g0.min())
    data["Ks"] = Ks

    gates_ok = (t0 < 1e-9 and nk == 2 and t1 < 1e-12
                and anchor_ok and t2a and mono and not trunc
                and all(bandgates) and eta_anchor_ok)
    if not gates_ok:
        status = "TOWER-TRUNCATED" if trunc else "V369-ANCHOR-DRIFT"
        ok = False
    elif not slope_ok:
        status = "CASIMIR-SLOPE-DRIFT-%.3f" % slope
        ok = False
    elif eta_star is not None:
        status = ("BINDING-ONSET-%.4f (composite band exits the true "
                  "nearest-channel continuum band above eta*; Casimir "
                  "slope %.3f vs booked 18; eta declared-not-read in "
                  "the published theory - conditional sentence only, "
                  "O-card open)" % (eta_star, slope))
        ok = True
    else:
        status = ("NO-ONSET-IN-SCAN (gap stays under 1e-3 through "
                  "eta = 1.00 - honest red for the first interaction "
                  "cut)")
        ok = True
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V369 done] %s" % status, flush=True)
    out = {"id": "V369_pair_binding_onset", "tier": 2, "residual": True,
           "property": "does the certified bond-current interaction "
                       "open a binding window for the singlet pair "
                       "below the nearest continuum channel",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v369_pair_binding_onset.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
