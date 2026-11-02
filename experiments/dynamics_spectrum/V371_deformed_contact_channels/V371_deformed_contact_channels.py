#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V371 - deformed contact channels (formally booked BEFORE the run
in review/pre_registrations_root.md; executable first cut of the V371
shape; the full torus lattice-sum spectral question stays booked).

Under the published bond-current deformation e^{i eta Sum_b J_b.L_b},
do the certified singlet channels acquire channel-dependent phases
that RE-OPEN a gap of the pair against the two-quanta continuum at
the contact cut (door n°2's opening condition at the scope V372
licenses: pair at contact bond, one flux line of length <= 2).

Machine (booked):
 T0 kernel counts: ell=1 fiber (46080) ker = 2 (V366/V369
    machinery rebuilt); ell=2 singlet basis built inside the shared
    isotypic subspace (2+90+420) x (1+45+210) = 131072 with ker
    C2 = 6 (V372 count re-measured in-script); every assembled fiber
    vector carries ||C2 v|| <= 1e-8 and Gram defect <= 1e-12.
 A  deformation spectra: 2x2 ell=1 block must reproduce -64 I_2
    (V370 in-repo anchor); 6x6 D2 eigenvalues d_1..d_6 booked as
    the non-phase fact; J.L spectrum [jmin, jmax] on the full E4
    fiber (the continuum's deformation range, measured).
 C  window scan at the V368 kinematics (verbatim machinery, seed
    0x51): free-step min gap must reproduce V368's 2.017083e-05 to
    1e-9; eta grid [0,1] step 0.02; continuum points smeared by the
    measured deformation band eta*(jmax-jmin); channel k composite
    set shifted by eta*d_k. GATE BOOKED: re-opened iff
    max_eta window >= 0.05 rad/step.
Verdicts: THRESHOLD-REOPENED-CONTACT-CUT / NO-THRESHOLD-CONTACT-CUT
/ MACHINERY-DRIFT.
"""
import json
import os
import sys

import numpy as np
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

KTOL = 1e-8
GRID_N = 6


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def modpi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def circmin(c, T):
    c = np.asarray(c)
    j = np.searchsorted(T, c)
    best = np.full(c.shape, np.inf)
    for jj in (j - 1, j):
        valid = (jj >= 0) & (jj < len(T))
        d = np.where(valid, np.abs(c - T[np.clip(jj, 0, len(T) - 1)]),
                     np.inf)
        best = np.minimum(best, d)
    best = np.minimum(best, np.abs(c - T[0] - 2 * np.pi))
    best = np.minimum(best, np.abs(c - T[-1] + 2 * np.pi))
    return best


def _get_cp():
    try:
        import cupy
        cupy.cuda.Stream  # noqa: B018
        return cupy
    except Exception:      # noqa: BLE01
        return None


def main():
    checks, data = [], {}
    cp = _get_cp()

    def log(m):
        print(m, flush=True)

    Ts = [np.asarray(t, dtype=np.complex128) for t in root_generators()]
    f, _k, _a = structure_constants(Ts)
    Ls = [np.asarray(l, dtype=np.complex128)
          for l in adjoint_link_operators(f)]

    # ---------- ell=1: E4 fiber kernel + the -64 anchor ----------
    def e4_apply(v, opJL=False):
        h = v.reshape(32, 32, 45)
        if opJL:
            y = np.zeros_like(h)
            for a in range(45):
                t = np.einsum("Ff,abf->abF", Ls[a], h)
                y += np.einsum("Aa,abF->AbF", Ts[a], t)
            return y.ravel()
        y = np.zeros_like(h)
        for a in range(45):
            q = (np.einsum("Aa,abf->Abf", Ts[a], h)
                 + np.einsum("Bb,abf->aBf", Ts[a], h)
                 + np.einsum("Ff,abf->abF", Ls[a], h))
            y += (np.einsum("Aa,abf->Abf", Ts[a], q)
                  + np.einsum("Bb,abf->aBf", Ts[a], q)
                  + np.einsum("Ff,abf->abF", Ls[a], q))
        return y.ravel()

    op1 = spla.LinearOperator((46080, 46080), matvec=e4_apply,
                              dtype=np.complex128)
    w1, v1 = spla.eigsh(op1, k=10, which="SM", return_eigenvectors=True)
    idx = np.argsort(np.abs(w1))
    w1, v1 = w1[idx], np.asarray(v1)[:, idx]
    W1 = v1[:, np.abs(w1) < KTOL]
    ker1 = W1.shape[1]
    G = W1.conj().T @ W1
    gv, GU = np.linalg.eigh(G)
    assert gv.min() > 1e-8, "ell=1 kernel Gram rank-deficient"
    W1 = W1 @ (GU * gv ** -0.5) @ GU.conj().T
    res1 = float(np.max(np.abs(np.stack(
        [e4_apply(W1[:, p]) for p in range(ker1)], axis=1))))
    # J.L on the ell=1 kernel: must be -64 I_2 (V370)
    JW1 = np.stack([e4_apply(W1[:, p], opJL=True) for p in range(ker1)],
                   axis=1)
    D1 = W1.conj().T @ JW1
    a1 = (ker1 == 2 and res1 <= KTOL
          and np.max(np.abs(D1 + 64.0 * np.eye(2))) <= 1e-6)
    checks.append(("T0 ell=1: ker %d (2), ||C2 W|| %.1e (<=1e-8), "
                   "J.L block diag %s (V370 anchor -64 I_2): %s"
                   % (ker1, res1, np.round(np.diag(D1), 6), a1), a1))

    # J.L spectrum on the whole E4 fiber: continuum smear band
    opJL = spla.LinearOperator((46080, 46080),
                               matvec=lambda x: e4_apply(x, opJL=True),
                               dtype=np.complex128)
    hi = spla.eigsh(opJL, k=3, which="LA", return_eigenvectors=False)
    lo = spla.eigsh(opJL, k=3, which="SA", return_eigenvectors=False)
    jmax, jmin = float(np.max(hi)), float(np.min(lo))
    data["JL_fiber_band"] = [jmin, jmax]
    log("J.L band on E4 fiber: [%.3f, %.3f]" % (jmin, jmax))

    # ---------- ell=2: exact singlet basis in the shared sector ----------
    I32, I45 = np.eye(32), np.eye(45)
    QMa = [np.kron(t, I32) + np.kron(I32, t) for t in Ts]      # 1024
    QFa = [np.kron(l, I45) + np.kron(I45, l) for l in Ls]      # 2025
    CM = np.zeros((1024, 1024), dtype=np.complex128)
    for q in QMa:
        CM += q @ q
    CF = np.zeros((2025, 2025), dtype=np.complex128)
    for q in QFa:
        CF += q @ q
    vm, Um = np.linalg.eigh((CM + CM.conj().T) / 2.0)
    vf, Uf = np.linalg.eigh((CF + CF.conj().T) / 2.0)

    def block(vals, vecs, target, tol=1e-6):
        sel = np.abs(vals - target) <= tol * max(1.0, abs(target))
        return vecs[:, sel]

    MB = np.concatenate([block(vm, Um, 0.0), block(vm, Um, 64.0),
                         block(vm, Um, 96.0)], axis=1)
    FB = np.concatenate([block(vf, Uf, 0.0), block(vf, Uf, 64.0),
                         block(vf, Uf, 96.0)], axis=1)
    nm, nf = MB.shape[1], FB.shape[1]
    shape_ok = (nm, nf) == (512, 256)
    checks.append(("ell=2 shared-sector bases: (%d,%d) = (512,256)"
                   " (2+90+420 x 1+45+210): %s" % (nm, nf, shape_ok),
                   shape_ok))

    Mred = [MB.conj().T @ q @ MB for q in QMa]     # 512x512
    Fred = [FB.conj().T @ q @ FB for q in QFa]     # 256x256
    ndim = nm * nf

    def red_matvec(x):
        X = x.reshape(nm, nf)
        Y = np.zeros_like(X)
        for a in range(45):
            Y += Mred[a] @ X + X @ Fred[a].T
        return Y.ravel()
    Mredc = [cp.asarray(m) for m in Mred] if cp else None
    Fredc = [cp.asarray(m) for m in Fred] if cp else None
    if cp:
        def red_matvec_gpu(x):
            X = cp.asarray(x).reshape(nm, nf)
            Y = cp.zeros_like(X)
            for a in range(45):
                Y += Mredc[a] @ X + X @ Fredc[a].T
            return cp.asnumpy(Y.ravel())
        op2 = spla.LinearOperator((ndim, ndim), matvec=red_matvec_gpu,
                                  dtype=np.complex128)
    else:
        op2 = spla.LinearOperator((ndim, ndim), matvec=red_matvec,
                                  dtype=np.complex128)
    w2, v2 = spla.eigsh(op2, k=14, which="SM", return_eigenvectors=True)
    idx = np.argsort(np.abs(w2))
    w2, v2 = w2[idx], np.asarray(v2)[:, idx]
    Xk = v2[:, np.abs(w2) < KTOL]
    ker2 = Xk.shape[1]
    G2 = Xk.conj().T @ Xk
    g2v, g2u = np.linalg.eigh(G2)
    assert g2v.min() > 1e-8, "ell=2 sector Gram rank-deficient"
    Xk = Xk @ (g2u * g2v ** -0.5) @ g2u.conj().T
    # lift to the full fiber and gate residuals there
    Wf = np.stack([((MB @ Xk[:, p]) @ FB.T).ravel()
                   for p in range(ker2)], axis=1)
    def c2_full(v):
        h = v.reshape(32, 32, 45, 45)
        y = np.zeros_like(h)
        for a in range(45):
            q = (np.einsum("Aa,abfg->Abfg", Ts[a], h)
                 + np.einsum("Bb,abfg->aBfg", Ts[a], h)
                 + np.einsum("Pf,abfg->abPg", Ls[a], h)
                 + np.einsum("Qg,abfg->abfQ", Ls[a], h))
            y += (np.einsum("Aa,abfg->Abfg", Ts[a], q)
                  + np.einsum("Bb,abfg->aBfg", Ts[a], q)
                  + np.einsum("Pf,abfg->abPg", Ls[a], q)
                  + np.einsum("Qg,abfg->abfQ", Ls[a], q))
        return y
    res2 = float(max(np.max(np.abs(c2_full(Wf[:, p]).reshape(
        32, 32, 45, 45))) for p in range(ker2)))
    gw = float(np.max(np.abs(Wf.conj().T @ Wf - np.eye(ker2))))
    a2 = (ker2 == 6 and res2 <= 1e-8 and gw <= 1e-12)
    checks.append(("T0 ell=2: lifted ker %d (V372's 6), max ||C2 v||"
                   " %.1e (<=1e-8), Gram defect %.1e: %s"
                   % (ker2, res2, gw, a2), a2))

    # ---------- B: deformation 6x6 on the lifted basis ----------
    def jl2(v):
        h = v.reshape(32, 32, 45, 45)
        y = np.zeros_like(h)
        for a in range(45):
            t1 = np.einsum("Pf,abfg->abPg", Ls[a], h)
            y += np.einsum("Aa,abPg->AbPg", Ts[a], t1)
            t2 = np.einsum("Qg,abfg->abfQ", Ls[a], h)
            y += np.einsum("Bb,abfQ->aBfQ", Ts[a], t2)
        return y
    JW2 = np.stack([jl2(Wf[:, p]) for p in range(ker2)], axis=1)
    D2 = Wf.conj().T @ JW2.reshape(2073600, ker2)
    D2 = (D2 + D2.conj().T) / 2.0
    dk = np.linalg.eigvalsh(D2)
    data["D2_eigs"] = [float(x) for x in dk]
    data["D2_nonphase"] = float(np.max(np.abs(dk + 64.0)))
    log("ell=2 deformation eigs d_k: %s (max |d+64| %.3f)"
        % (np.round(dk, 3), data["D2_nonphase"]))

    # ---------- C: window scan at the V368 kinematics ----------
    _Bm, U32, _U16 = SC.sectors(SC.G10)
    U6 = np.linalg.matrix_power(U32, 6)
    iota = np.angle(np.linalg.eigvals(U6))
    grid = [(2 * np.pi * i / GRID_N, 2 * np.pi * j / GRID_N,
             2 * np.pi * l / GRID_N)
            for i in range(GRID_N) for j in range(GRID_N)
            for l in range(GRID_N)]
    betas = np.array([np.angle(np.linalg.eigvals(SC.bond_product(q)))
                      for q in grid])
    v251 = json_art("V251_gauged_two_particle_sector")
    rows = {k: float(v) for k, v in
            v251["data"]["sensitivity"]["omega_rows"].items()}
    phys = [v for v in rows.values() if abs(v) > 1e-9]
    rng = np.random.default_rng(0x51)
    idxg = np.arange(GRID_N ** 3)
    tri = np.array([(a, b, c) for a in range(GRID_N)
                    for b in range(GRID_N) for c in range(GRID_N)])
    cont_sets, comp_sets = [], []
    for K in range(GRID_N ** 3):
        k1 = idxg
        k2 = ((tri - tri[K]) % GRID_N) @ np.array([36, 6, 1])
        bA, bB = betas[k1], betas[k2]
        comp = modpi((bA[:, :, None] + bB[:, None, :]).ravel()[:, None]
                     + np.array(phys)[None, :]).ravel()
        samp = rng.choice(216, size=24, replace=False)
        s1, s2 = samp, ((tri[samp] - tri[K]) % GRID_N) @ np.array([36,
                                                                  6, 1])
        sb = (betas[s1][:, :, None] + betas[s2][:, None, :]).ravel()
        cont = np.sort(modpi(sb[:, None] + (iota[:, None]
                                            + iota[None, :]).ravel()
                             [None, :]).ravel())
        cont_sets.append(cont)
        comp_sets.append(comp)
    gaps_free = [float(circmin(comp_sets[K], cont_sets[K]).min())
                 for K in range(GRID_N ** 3)]
    free_min = min(gaps_free)
    v368 = json_art("V368_torus_rigidity_window")
    v368_min = float(v368["data"]["gaps_min"])
    anchor368 = abs(free_min - v368_min) < 1e-9
    checks.append(("C anchor: free-step min gap %.9f reproduces V368"
                   " %.9f to 1e-9: %s" % (free_min, v368_min, anchor368),
                   anchor368))

    etas = np.round(np.arange(0.0, 1.001, 0.02), 3)
    all_d = np.concatenate([np.diag(D1), dk])
    smear_unit = jmax - jmin
    windows = []
    for eta in etas:
        sm = eta * smear_unit
        wmin = np.inf
        if sm >= 2 * np.pi:
            windows.append(0.0)
            continue
        for k, d in enumerate(all_d):
            for K in range(GRID_N ** 3):
                gk = circmin(modpi(comp_sets[K] + eta * d),
                             cont_sets[K]).min() - sm / 2.0
                wmin = min(wmin, max(0.0, gk))
        windows.append(float(wmin))
    windows = np.array(windows)
    best_eta = int(np.argmax(windows))
    data["window_scan"] = {"etas": [float(e) for e in etas],
                           "windows": [float(w) for w in windows],
                           "best_eta": float(etas[best_eta]),
                           "best_window": float(windows[best_eta])}
    log("window scan: best %.4f at eta %.2f (smear unit %.1f;"
        " eta*smear >= 2pi zeroes every window)"
        % (windows[best_eta], etas[best_eta], smear_unit))

    okm = a1 and a2 and shape_ok and anchor368
    reopened = windows[best_eta] >= 0.05 and best_eta > 0
    if not okm:
        status, ok = "V371-MACHINERY-DRIFT", False
    elif reopened:
        status = ("THRESHOLD-REOPENED-CONTACT-CUT (eta* = %.2f: min-"
                  "over-K channel window %.4f >= 0.05 rad/step - the "
                  "bond-current deformation re-isolates a singlet "
                  "channel against the smeared continuum at the "
                  "contact cut)" % (etas[best_eta],
                                    windows[best_eta]))
        ok = True
    else:
        status = ("NO-THRESHOLD-CONTACT-CUT (the deformation splits "
                  "the ell=2 channels - max |d_k + 64| = %.1f, so it "
                  "is NOT a phase - but no eta clears the smeared "
                  "two-quanta continuum at the contact cut: best "
                  "window %.4f at eta %.2f, and every eta with "
                  "eta*(jmax-jmin) >= 2pi (>= %.3f) zeroes every "
                  "window by continuum smear alone)"
                  % (data["D2_nonphase"], windows[best_eta],
                     etas[best_eta], 2 * np.pi / smear_unit))
        ok = True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V371 done] %s" % status, flush=True)
    out = {"id": "V371_deformed_contact_channels", "tier": 2,
           "residual": True,
           "property": "does the bond-current deformation re-open a "
                       "pair window at the contact cut (door n°2 "
                       "first decidable cut)",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v371_deformed_contact_channels.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
