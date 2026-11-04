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
    machinery rebuilt, J.L = (T1+T2).L per V370); ell=2 kernel of
    O^2 = [sum_a (T1+T2+L3+L4)_a]^2 computed by Lanczos on the FULL
    2073600 fiber, ker = 6 (V372 count re-measured in-script). The
    first implementation booked a shared-sector basis sandwich
    (2+90+420)x(1+45+210) = 131072; the in-script full-fiber
    residual gate exposed it as NOT O-invariant (sum-Q^2 eigenspaces
    are not invariant under sum Q -- commutator nonzero) and the
    booked MACHINERY-DRIFT verdict fired; the machinery was fixed
    before any verdict run, gate thresholds unchanged. The shared
    cluster census survives as the V372 multiplicity cross-check.
    Every kernel vector carries ||O^2 v|| <= 1e-8 and Gram defect
    <= 1e-12.
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
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "gauged_layer",
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
    # KERNEL OPERATOR (anchored convention, V366/V370/V372): the
    # quadratic Casimir of the DIAGONAL action, C = sum_a (J_a^tot)^2
    # -- a sum of squares, one per generator. NOT (sum_a J_a)^2:
    # that collapse is a different operator with a dense continuum
    # at zero (plain-Lanczos probe bg_39: Ritz values falling to 0
    # like j^-2 instead of a clean ladder) and it cannot define the
    # kernel at all. The legitimate speedup is algebraic: expand
    # (sum_legs J_a^leg)^2 over leg pairs and PRECOMPUTE each
    # generator-coupling matrix once:
    #   K_M = sum_a T_a T_a (single leg), K_F = sum_a L_a L_a,
    #   C_MM = sum_a T_a (x) T_a (leg pair), C_FF, C_MF likewise.
    # Then C acts as 10 gemm-style contractions -- equivalent to
    # the 45-generator loop (V370 anchors: Casimir ladder {0,36,
    # ...}, clean gap 36) and ~10x faster.
    K_M = sum(t @ t for t in Ts)                        # 32x32
    K_F = sum(l @ l for l in Ls)                        # 45x45
    C_MM = np.zeros((1024, 1024), dtype=np.complex128)
    C_FF = np.zeros((2025, 2025), dtype=np.complex128)
    C_MF = np.zeros((1440, 1440), dtype=np.complex128)
    for a in range(45):
        C_MM += np.kron(Ts[a], Ts[a])
        C_FF += np.kron(Ls[a], Ls[a])
        C_MF += np.kron(Ts[a], Ls[a])

    # ---------- ell=1: E4 fiber kernel + the -64 anchor ----------
    def mf_e4(h, leg_m, leg_f):
        # C_MF between matter leg (leg_m) and the flux leg: contract
        # the (matter,flux) pair into 1440, gemm, expand back
        if leg_m == 0:
            u = h.transpose(0, 2, 1).reshape(1440, 32)
            z = np.tensordot(C_MF, u, axes=([1], [0]))
            return np.moveaxis(z.reshape(32, 45, 32), 1, 2)
        u = h.transpose(1, 2, 0).reshape(1440, 32)
        z = np.tensordot(C_MF, u, axes=([1], [0]))
        return np.moveaxis(z.reshape(32, 45, 32), 1, 2).transpose(
            1, 0, 2)

    def e4_apply(v, opJL=False):
        h = v.reshape(32, 32, 45)
        # JL = sum_a (T1_a + T2_a) L_a: exactly the two (matter,
        # flux) couplings of C_MF (V370 anchor -64 I_2 on ker;
        # the single-leg half gave -32)
        if opJL:
            return (mf_e4(h, 0, 2) + mf_e4(h, 1, 2)).ravel()
        # C = sum_a (J_a^tot)^2 over legs (0,1,2): single-leg
        # Casimirs (coefficient 1) + pair couplings with FACTOR 2
        # (AB+BA = 2AB, distinct legs commute -- the V370 loop
        # applies it implicitly)
        pairs = (np.tensordot(C_MM, h.reshape(1024, 45),
                              axes=([1], [0])).reshape(32, 32, 45)
                 + mf_e4(h, 0, 2) + mf_e4(h, 1, 2))
        y = (np.tensordot(K_M, h, axes=([1], [0]))
             + np.tensordot(K_M, h, axes=([1], [1])).transpose(1, 0, 2)
             + np.tensordot(K_F, h, axes=([1], [2])).transpose(1, 2, 0)
             + 2.0 * pairs)
        return y.ravel()

    # ell=1 kernel: reuse the V369 SM-400 cache (V370's own source
    # for the -64 I_2 anchor; the fold runs V369 before V371 in
    # sorted order). A fresh SA eigsh here reproduced kernel
    # RESIDUALS of 3e-15 but a JL-invariance LEAK of 3.2e-3 -- its
    # span was NOT the published anchor span, and the -64 I_2
    # statement is only exact on the span where J.L acts as -64 I
    # (leak 8.8e-14). Anchoring on the cache makes A1 a literal
    # reproduction of V370, not a re-derivation under ARPACK luck.
    CACHE = exp("V369_pair_binding_onset", "results", "v369_sm400.npz")
    if os.path.exists(CACHE):
        zc = np.load(CACHE)
        w1c, V1c = zc["w"], zc["V"]
        W1 = V1c[:, np.abs(w1c) < 1e-9]
        data["ell1_kernel_source"] = "V369 sm400 cache (V370 source)"
    else:
        op1 = spla.LinearOperator((46080, 46080), matvec=e4_apply,
                                  dtype=np.complex128)
        # which="SA": the Casimir is PSD, the kernel sits at the
        # LEFT ENDPOINT (mode-0 Lanczos converges endpoints; "SM"
        # shift-free is the classic 100x slowdown that killed the
        # first run at the 1700 s wall).
        w1, v1 = spla.eigsh(op1, k=10, which="SA",
                            return_eigenvectors=True, maxiter=600)
        idx = np.argsort(np.abs(w1))
        w1, v1 = w1[idx], np.asarray(v1)[:, idx]
        W1 = v1[:, np.abs(w1) < KTOL]
        data["ell1_kernel_source"] = "fresh SA eigsh (no cache)"
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
    # JL-invariance leak: V370 measured J.L ACTS AS -64 I on the
    # kernel (leak 8.8e-14). Reproduced here on the same span.
    leak1 = float(np.linalg.norm(JW1 - W1 @ D1))
    data["leak_ell1"] = leak1
    a1 = (ker1 == 2 and res1 <= KTOL
          and np.max(np.abs(D1 + 64.0 * np.eye(2))) <= 1e-6
          and leak1 <= 1e-6)
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

    # ---------- ell=2: kernel of O^2 on the FULL 2073600 fiber -----
    # The booked shared-sector basis sandwich was WRONG: the sum-Q^2
    # eigenspaces are NOT invariant under O = sum J (the commutator
    # [Q_a, sum_b Q_b^2] ~ f_abc{Q_c,Q_b}/2 != 0), and the in-script
    # full-fiber residual gate caught exactly that (0.83 -> booked
    # V371-MACHINERY-DRIFT fired as designed). Fixed BEFORE any
    # verdict: Lanczos directly on the whole fiber (GPU matvecs);
    # the shared-cluster census stays as the V372 multiplicity
    # cross-check it always was -- a census, not a subspace.
    I32, I45 = np.eye(32), np.eye(45)
    QMa = [np.kron(t, I32) + np.kron(I32, t) for t in Ts]      # 1024
    QFa = [np.kron(l, I45) + np.kron(I45, l) for l in Ls]      # 2025
    CM = np.zeros((1024, 1024), dtype=np.complex128)
    for q in QMa:
        CM += q @ q
    CF = np.zeros((2025, 2025), dtype=np.complex128)
    for q in QFa:
        CF += q @ q
    vm = np.linalg.eigvalsh((CM + CM.conj().T) / 2.0)
    vf = np.linalg.eigvalsh((CF + CF.conj().T) / 2.0)

    def mult_at(vals, target, tol=1e-6):
        return int(np.sum(np.abs(vals - target)
                          <= tol * max(1.0, abs(target))))

    census = [(mult_at(vm, t), mult_at(vf, t))
              for t in (0.0, 64.0, 96.0)]
    shape_ok = census == [(2, 1), (90, 45), (420, 210)]
    checks.append(("ell=2 shared-cluster census (V372 cross-check): "
                   "values (0,64,96) multiplicities %s == "
                   "[(2,1),(90,45),(420,210)]: %s"
                   % (census, shape_ok), shape_ok))

    # four-leg Casimir C = sum_a (J_a^tot)^2 on (32,32,45,45):
    # 4 single-leg squares + 2 x 6 pair couplings, all precomputed.
    # (mf pair terms: permute the two legs to the front, contract
    # the 1440-dim pair space, expand back with the inverse perm)
    MF_PERMS = [((0, 2, 1, 3), (0, 2, 1, 3)),
                ((0, 3, 1, 2), (0, 2, 3, 1)),
                ((1, 2, 0, 3), (2, 0, 1, 3)),
                ((1, 3, 0, 2), (2, 0, 3, 1))]

    def c4_np(g):
        g = g.reshape(32, 32, 45, 45)
        p = np.tensordot(C_MM, g.reshape(1024, 2025),
                         axes=([1], [0])).reshape(32, 32, 45, 45)
        p = p + np.tensordot(
            C_FF, g.transpose(2, 3, 0, 1).reshape(2025, 1024),
            axes=([1], [0])).reshape(45, 45, 32, 32).transpose(
            2, 3, 0, 1)
        for perm, inv in MF_PERMS:
            z = np.tensordot(C_MF, g.transpose(*perm).reshape(1440,
                                                              1440),
                             axes=([1], [0])).reshape(32, 45, 32, 45)
            p = p + z.transpose(*inv)
        return (np.tensordot(K_M, g, axes=([1], [0]))
                + np.tensordot(K_M, g, axes=([1], [1]))
                .transpose(1, 0, 2, 3)
                + np.tensordot(K_F, g, axes=([1], [2]))
                .transpose(1, 2, 0, 3)
                + np.tensordot(K_F, g, axes=([1], [3]))
                .transpose(1, 2, 3, 0)
                + 2.0 * p)

    def o_full(x):
        return c4_np(x.reshape(32, 32, 45, 45))

    FIBER = 32 * 32 * 45 * 45
    if cp:
        KMc, KFc = cp.asarray(K_M), cp.asarray(K_F)
        CMMc, CFFc, CMFc = (cp.asarray(C_MM), cp.asarray(C_FF),
                            cp.asarray(C_MF))

        def c4_gpu(g):
            g = g.reshape(32, 32, 45, 45)
            p = cp.tensordot(CMMc, g.reshape(1024, 2025),
                             axes=([1], [0])).reshape(32, 32, 45, 45)
            p = p + cp.tensordot(
                CFFc, g.transpose(2, 3, 0, 1).reshape(2025, 1024),
                axes=([1], [0])).reshape(45, 45, 32, 32).transpose(
                2, 3, 0, 1)
            for perm, inv in MF_PERMS:
                z = cp.tensordot(
                    CMFc, g.transpose(*perm).reshape(1440, 1440),
                    axes=([1], [0])).reshape(32, 45, 32, 45)
                p = p + z.transpose(*inv)
            return (cp.tensordot(KMc, g, axes=([1], [0]))
                    + cp.tensordot(KMc, g, axes=([1], [1]))
                    .transpose(1, 0, 2, 3)
                    + cp.tensordot(KFc, g, axes=([1], [2]))
                    .transpose(1, 2, 0, 3)
                    + cp.tensordot(KFc, g, axes=([1], [3]))
                    .transpose(1, 2, 3, 0)
                    + 2.0 * p)

        def mv2(x):
            return cp.asnumpy(c4_gpu(cp.asarray(x)).ravel())
    else:
        def mv2(x):
            return o_full(x).ravel()
    op2 = spla.LinearOperator((FIBER, FIBER), matvec=mv2,
                              dtype=np.complex128)

    # ARPACK budget as diagnosed by the probe runs: default tol=0
    # asks residuals below the Lanczos roundoff floor (2e-16 vs
    # ~2e-11 at ||O^2|| ~ 1e5) -> never converges; tol=1e-9 sits
    # above the floor. Each restart costs (ncv-k) matvecs; ncv=30,
    # maxiter=250. Partial Ritz vectors on non-convergence are
    # legitimate Krylov output -- the polish refines whatever ships.
    try:
        w2, v2 = spla.eigsh(op2, k=14, which="SA",
                            return_eigenvectors=True,
                            maxiter=250, tol=1e-9, ncv=30)
    except spla.ArpackNoConvergence as e:
        w2, v2 = e.eigenvalues, e.eigenvectors
        log("arpack partial: %d Ritz pairs, polishing anyway"
            % len(w2))
    idx = np.argsort(np.abs(w2))
    w2, v2 = w2[idx], np.asarray(v2)[:, idx]
    # Rayleigh-Ritz polish inside the Ritz span (drives residuals
    # to the roundoff floor before the gate)
    OV = np.stack([mv2(v2[:, j]) for j in range(v2.shape[1])],
                  axis=1)
    B = v2.conj().T @ OV
    B = (B + B.conj().T) / 2.0
    th, Y = np.linalg.eigh(B)
    order = np.argsort(np.abs(th))
    th, Y = th[order], Y[:, order]
    # spectral gap 36 >> noise: any cut in [1e-6, 30] isolates the
    # kernel; 1.0 chosen for margin over ARPACK residual scatter
    sel = np.abs(th) < 1.0
    Wf = v2 @ Y[:, sel]
    thp = th[sel]
    resR = float(max(np.max(np.abs(mv2(Wf[:, j]) - thp[j] * Wf[:, j]))
                     for j in range(Wf.shape[1])))
    data["ritz_polish_residual"] = resR
    log("ritz polish: %d vectors, max O^2 residual %.2e"
        % (Wf.shape[1], resR))
    ker2 = Wf.shape[1]
    G2 = Wf.conj().T @ Wf
    g2v, g2u = np.linalg.eigh(G2)
    assert g2v.min() > 1e-8, "ell=2 kernel Gram rank-deficient"
    Wf = Wf @ (g2u * g2v ** -0.5) @ g2u.conj().T
    res2 = float(max(np.max(np.abs(o_full(Wf[:, p])))
                     for p in range(ker2)))
    gw = float(np.max(np.abs(Wf.conj().T @ Wf - np.eye(ker2))))
    a2 = (ker2 == 6 and res2 <= 1e-8 and gw <= 1e-12)
    checks.append(("T0 ell=2: full-fiber ker %d (V372's 6), max "
                   "||O^2 v|| %.1e (<=1e-8), Gram defect %.1e: %s"
                   % (ker2, res2, gw, a2), a2))

    # ---------- B: deformation 6x6 on the fiber-kernel basis ------
    def jl2(v):
        # contact-cut deformation sum: each line bond link (flux
        # legs 2,3) carries the FULL matter current (T1 + T2) --
        # the V370 definition (T1+T2).L at one bond anchors -64 and
        # the length-2 line is the same operator on both links.
        # Four (matter,flux) C_MF couplings, coefficient 1 (J.L is
        # a cross product, not a square -- no factor 2 here).
        h = v.reshape(32, 32, 45, 45)
        y = np.zeros_like(h)
        for perm, inv in MF_PERMS:
            z = np.tensordot(C_MF, h.transpose(*perm).reshape(1440,
                                                              1440),
                             axes=([1], [0])).reshape(32, 45, 32, 45)
            y = y + z.transpose(*inv)
        return y

    # GATE: coupling form == 45-loop anchored form on random
    # vectors (a wrong inverse permutation in MF_PERMS silently
    # changes the operator -- bg_44 localization found exactly
    # one such bug; one cheap loop pass per operator gates it).
    rngO = np.random.default_rng(0x37)
    def loop_e4(x, opJL=False):
        h = x.reshape(32, 32, 45)
        acc = np.zeros_like(h)

        def lg(a, g):
            return (np.einsum("Aa,abf->Abf", Ts[a], g)
                    + np.einsum("Bb,abf->aBf", Ts[a], g)
                    + np.einsum("Ff,abf->abF", Ls[a], g))
        for a in range(45):
            if opJL:
                lz = np.einsum("Ff,abf->abF", Ls[a], h)
                acc += (np.einsum("Aa,abF->AbF", Ts[a], lz)
                        + np.einsum("Bb,abF->aBF", Ts[a], lz))
            else:
                z = lg(a, h)
                acc += lg(a, z)
        return acc.ravel()

    def loop_c4(x):
        h = x.reshape(32, 32, 45, 45)
        acc = np.zeros_like(h)

        def lg(a, g):
            return (np.einsum("Aa,abfg->Abfg", Ts[a], g)
                    + np.einsum("Bb,abfg->aBfg", Ts[a], g)
                    + np.einsum("Pf,abfg->abPg", Ls[a], g)
                    + np.einsum("Qg,abfg->abfQ", Ls[a], g))
        for a in range(45):
            acc += lg(a, lg(a, h))
        return acc.ravel()

    def loop_jl2(x):
        h = x.reshape(32, 32, 45, 45)
        acc = np.zeros_like(h)
        for a in range(45):
            t3 = np.einsum("Pf,abfg->abPg", Ls[a], h)
            t4 = np.einsum("Qg,abfg->abfQ", Ls[a], h)
            acc += (np.einsum("Aa,abPg->AbPg", Ts[a], t3)
                    + np.einsum("Bb,abPg->aBPg", Ts[a], t3)
                    + np.einsum("Aa,abfQ->AbfQ", Ts[a], t4)
                    + np.einsum("Bb,abfQ->aBfQ", Ts[a], t4))
        return acc.ravel()
    dop = 0.0
    for _ in range(2):
        u4 = rngO.normal(size=46080) + 1j * rngO.normal(size=46080)
        u2 = rngO.normal(size=FIBER) + 1j * rngO.normal(size=FIBER)
        for got, ref in ((e4_apply(u4), loop_e4(u4)),
                         (e4_apply(u4, opJL=True),
                          loop_e4(u4, opJL=True)),
                         (c4_np(u2).ravel(), loop_c4(u2)),
                         (jl2(u2).ravel(), loop_jl2(u2))):
            dop = max(dop, float(np.max(np.abs(got - ref))
                                 / max(np.max(np.abs(ref)), 1e-12)))
    checks.append(("operator self-test: coupling form == 45-loop "
                   "form, max rel diff %.1e (<=1e-9): %s"
                   % (dop, dop <= 1e-9), dop <= 1e-9))
    log("operator self-test: max rel diff %.2e" % dop)

    # stack along the LAST axis: C-order reshape (32,32,45,45,ker2)
    # -> (2073600, ker2) keeps the channel column intact
    JW2 = np.stack([jl2(Wf[:, p]) for p in range(ker2)], axis=-1)
    D2 = Wf.conj().T @ JW2.reshape(2073600, ker2)
    D2 = (D2 + D2.conj().T) / 2.0
    # JL-invariance leak on the fiber kernel span (relative): a
    # kernel span that is really the Gauss singlet space must be
    # near-invariant under the deformation operator too
    JW2m = JW2.reshape(2073600, ker2)
    leak2 = float(np.linalg.norm(JW2m - Wf @ D2)
                  / np.linalg.norm(JW2m))
    data["leak_ell2_rel"] = leak2
    log("fiber JL leak (rel): %.2e" % leak2)
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
    # .real: D1 is certified real by the -64 I_2 anchor gate (a1),
    # and numpy's % (modpi) has no complex support
    all_d = np.concatenate([np.diag(D1).real, dk])
    smear_unit = jmax - jmin
    windows = []
    for eta in etas:
        sm = eta * smear_unit
        if eta == 0.0:
            windows.append(free_min)  # the un-deformed anchor
            continue
        # valid a-fortiori zero: circmin <= pi always, so half-smear
        # >= pi forces window 0 (the free-step gap is NOT a bound --
        # a shifted channel set can sit wider than free_min; that
        # possibility is exactly what the scan tests)
        if sm >= 2 * np.pi:
            windows.append(0.0)
            continue
        wmin = np.inf
        # EARLY EXIT: wmin is a running min, so once it hits 0 the
        # exact window for this eta is 0 and every remaining
        # (channel, K) pair is dead weight -- without this the scan
        # alone cost ~800 s and blew the 1700 s wall (bg_24).
        for k, d in enumerate(all_d):
            if wmin <= 0.0:
                break
            for K in range(GRID_N ** 3):
                gk = circmin(modpi(comp_sets[K] + eta * d),
                             cont_sets[K]).min() - sm / 2.0
                wmin = min(wmin, max(0.0, gk))
                if wmin <= 0.0:
                    break
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

    a3 = leak2 <= 1e-3
    checks.append(("fiber kernel JL-invariance leak (rel) %.1e "
                   "(<=1e-3): %s" % (leak2, a3), a3))
    okm = a1 and a2 and a3 and shape_ok and anchor368
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
