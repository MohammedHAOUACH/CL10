#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V373 - the full-lattice pair+string machinery (booked BEFORE any
run in review/pre_registrations_root.md; beyond_step programme, first
cut; tier 1; suite 286).

Property: does the gauge-invariant pair+flux-line Hilbert space of the
Lambda^3 torus, at fixed total momentum K, pair separation r over ALL
grid separations, one straight geodesic flux line of length
ell = d_geo(r), EXIST as working machinery - gated against everything
the reality-interface phase already measured?

Machine (booked, verbatim from the registry):
 - Imported verbatim: V79 adjoint_link_operators / root_generators /
   structure_constants (the V371 import); V101 zero_modes, the
   intertwiner kernel solver (J_sub W + W L^T), as the model for the
   chain solver - generalized to the MEASURED certified block
   structure (link-pair channel blocks, tensored with a single 45
   link for ell = 3), same operator and ARPACK class, applied by
   batched matmul/einsum (no 45^3 generator stack is ever
   materialized).
 - The kernel dimensions of the chain problem are MEASURED per chain
   length, never assumed: the fusion semiring of the 45-algebra
   beyond r in {1, 45, 210} is a measurement OUTPUT. The counts come
   from deflation rounds of the matrix-free normal operator: a
   CAPPED ARPACK pass (k = 4, ncv = 12, 300 matvec cap) on the
   shifted operator M - 32I (half the measured Casimir gap of 64),
   so the kernel is the spectral extreme ARPACK must find -
   ARPACK's relative stopping criterion can never converge a zero
   eigenvalue in place, and unshifted "SA" was measured to return
   only the 64+ junk. Acceptance is the conjunction (Rayleigh
   quotient <= 1e-4 AND residual <= 1e-5) on the TRUE M, both far
   below the gap; the cap bounds the no-kernel case and the
   deflation makes the count exact.
 - Effective Hermitian generator (booked form):
     H = H_hop (x) I_flux
       + I_matter (x) sum_{links in string} sum_a L_a L_a
       + eta * (bond-current deformation in the certified channel
                 basis).
   The per-link string Casimir coefficient is MEASURED from the
   adjoint representation (K_F spectrum), not assumed. The pair
   hopping H_hop is the V368 free two-quanta kinematics (seed 0x51)
   carried by a real-space one-particle bond propagation whose Bloch
   image is verified against SC.bond_product to 1e-12.
 Gates:
  A1 N=2 contact reduction reproduces V371 verbatim: ker counts 2
     (ell=1, RE-DERIVED in-script - no V369 cache) and 6 (ell=2),
     the deformation set d_k = {-96,-96,-64,-64,0,0}, and
     J.L = -64 I_2 on the ell=1 kernel with leak <= 1e-12.
  A2 eta = 0 free min gap reproduces V368's 2.017083e-05 to 1e-9.
  A3 chain-kernel structure (a length-ell line carries a chain of
     ell + 1 links, booked) anchors on V100/V101: pair-channel census
     (2 singlets, 90 of 45, 420 of 210); link-pair census I_2 = 1,
     single 45, 210 channel with the channel generator blocks gated
     against the measured Casimir per block; Hom(pair, 45) block
     counts == the T1 contact-fiber Lanczos kernel; chain ell=1
     kernel 6 with block pattern 2+2+2 over r = 1, 45, 210 AND the
     glued identity 6 x m_45^2 == the published 24-state V100/V101
     kernel (gram ranks 24/24, V101 top degeneracy 6); measured
     I_3 = 1; chain ell=2 kernel recorded (output); chain ell=3
     measured in the SINGLET sector (I_4 measured against V100's 6;
     non-singlet 4-link blocks exceed the honest runtime budget and
     are DECLARED UNMEASURED, never silent); the 2-link kernel
     tensors spanning the SAME 6-dim fiber kernel V371 measured
     (overlap gate); the deformation set reproduced in the kernel
     basis; and a second solver route (materialized link-pair
     generators) cross-checking the singlet 2- and 3-link counts.
Verdicts: FULL-LATTICE-MACHINERY-BUILT / V373-MACHINERY-DRIFT.
Scope booked: one flux line, straight and geodesic; the deformation
convention is the V370/V371 one (each line link carries the full
matter current (T1+T2).L); channel bases measured, never assumed;
the 4-link non-singlet sector declared UNMEASURED at runtime.
"""
import hashlib
import importlib.util
import json
import math
import os
import shutil
import sys
import time

# host courtesy: this box runs other services; cap BLAS and let the
# scheduler prefer them (the 15h swap-storm run blocked the machine).
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

# NOT numpy-pinnable (measured 2026-11-05, twice): the big chain blocks
# (s ~ 1.3e7) go through the host ARPACK fallback, whose workd is
# 3*ncv*s fp64 - ~10 GB - and whose matvecs take minutes each. The
# pinned run died at the cgroup cap after 2716 s; the unpinned run on a
# FULL card OOMed the card and the host (59/62 GiB). The published route
# is the GPU big-block path of hom_null (fp32 Lanczos + fp64 polish);
# it needs a GPU CARD WITH ROOM. The fold runner therefore pins it to a
# free card from the outside (CUDA_VISIBLE_DEVICES) and runs the whole
# sweep under a cgroup memory cap, so a runaway experiment is killed
# (visible FAILED) instead of freezing the box.

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

_spec = importlib.util.spec_from_file_location(
    "v101", exp("V101_bond_current_ground_state",
                "V101_bond_current_ground_state.py"))
V101 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V101)

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
    t0 = time.time()
    checks, data = [], {}
    cp = _get_cp()
    gpu = cp is not None

    def log(m):
        print("[%6.1fs] %s" % (time.time() - t0, m), flush=True)

    # ------------------------------------------------------------------
    # couplings, V371 verbatim
    # ------------------------------------------------------------------
    Ts = [np.asarray(t, dtype=np.complex128) for t in root_generators()]
    f, _k, _a = structure_constants(Ts)
    Ls = [np.asarray(l, dtype=np.complex128)
          for l in adjoint_link_operators(f)]
    K_M = sum(t @ t for t in Ts)                        # 32x32
    K_F = sum(l @ l for l in Ls)                        # 45x45
    C_MM = np.zeros((1024, 1024), dtype=np.complex128)
    C_FF = np.zeros((2025, 2025), dtype=np.complex128)
    C_MF = np.zeros((1440, 1440), dtype=np.complex128)
    for a in range(45):
        C_MM += np.kron(Ts[a], Ts[a])
        C_FF += np.kron(Ls[a], Ls[a])
        C_MF += np.kron(Ts[a], Ls[a])
    log("couplings precomputed (45 generators, adjoint links)")

    # ------------------------------------------------------------------
    # T1: ell=1 E4 fiber kernel, RE-DERIVED in-script (no cache); the
    # V371 lesson: Rayleigh-Ritz polish inside the Ritz span, then the
    # JL-leak gate certifies the RIGHT span.
    # ------------------------------------------------------------------
    def mf_e4(h, leg_m, leg_f):
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
        if opJL:
            return (mf_e4(h, 0, 2) + mf_e4(h, 1, 2)).ravel()
        pairs = (np.tensordot(C_MM, h.reshape(1024, 45),
                              axes=([1], [0])).reshape(32, 32, 45)
                 + mf_e4(h, 0, 2) + mf_e4(h, 1, 2))
        y = (np.tensordot(K_M, h, axes=([1], [0]))
             + np.tensordot(K_M, h, axes=([1], [1])).transpose(1, 0, 2)
             + np.tensordot(K_F, h, axes=([1], [2])).transpose(1, 2, 0)
             + 2.0 * pairs)
        return y.ravel()

    op1 = spla.LinearOperator((46080, 46080), matvec=e4_apply,
                              dtype=np.complex128)
    w1, v1 = spla.eigsh(op1, k=10, which="SA", return_eigenvectors=True,
                        maxiter=600, tol=1e-9, ncv=30)
    OV1 = np.stack([e4_apply(v1[:, j]) for j in range(v1.shape[1])],
                   axis=1)
    B1 = v1.conj().T @ OV1
    B1 = (B1 + B1.conj().T) / 2.0
    th1, Y1 = np.linalg.eigh(B1)
    o1 = np.argsort(np.abs(th1))
    th1, Y1 = th1[o1], Y1[:, o1]
    sel1 = np.abs(th1) < 1.0            # Casimir gap 36 >> noise
    W1 = v1 @ Y1[:, sel1]
    G1 = W1.conj().T @ W1
    g1v, g1u = np.linalg.eigh(G1)
    assert g1v.min() > 1e-8, "ell=1 kernel Gram rank-deficient"
    W1 = W1 @ (g1u * g1v ** -0.5) @ g1u.conj().T
    ker1 = W1.shape[1]
    res1 = float(np.max(np.abs(np.stack(
        [e4_apply(W1[:, p]) for p in range(ker1)], axis=1))))
    JW1 = np.stack([e4_apply(W1[:, p], opJL=True) for p in range(ker1)],
                   axis=1)
    D1 = W1.conj().T @ JW1
    D1 = (D1 + D1.conj().T) / 2.0
    leak1 = float(np.linalg.norm(JW1 - W1 @ D1))
    d1e = np.linalg.eigvalsh(D1)
    a1_jl = (ker1 == 2 and res1 <= KTOL
             and np.max(np.abs(d1e + 64.0)) <= 1e-12
             and leak1 <= 1e-12)
    checks.append(("A1a ell=1 in-script (no cache): ker %d (2), "
                   "||C2 W|| %.1e, J.L eig %s (-64,-64) <=1e-12, "
                   "span leak %.1e (<=1e-12): %s"
                   % (ker1, res1, np.round(d1e, 9), leak1, a1_jl),
                   a1_jl))
    log("ell=1: ker %d, J.L eig %s, leak %.1e"
        % (ker1, np.round(d1e, 9), leak1))

    # ------------------------------------------------------------------
    # T2: ell=2, the V371 c4 machinery verbatim (full 2073600 fiber,
    # GPU matvecs where cupy is present).
    # ------------------------------------------------------------------
    MF_PERMS = [((0, 2, 1, 3), (0, 2, 1, 3)),
                ((0, 3, 1, 2), (0, 2, 3, 1)),
                ((1, 2, 0, 3), (2, 0, 1, 3)),
                ((1, 3, 0, 2), (2, 0, 3, 1))]

    def c4_np(g):
        g = np.reshape(g, (32, 32, 45, 45))
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

    FIBER2 = 32 * 32 * 45 * 45

    if gpu:
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
            return c4_np(x).ravel()

    op2 = spla.LinearOperator((FIBER2, FIBER2), matvec=mv2,
                              dtype=np.complex128)
    try:
        w2, v2 = spla.eigsh(op2, k=14, which="SA",
                            return_eigenvectors=True,
                            maxiter=250, tol=1e-9, ncv=30)
    except spla.ArpackNoConvergence as e:
        w2, v2 = e.eigenvalues, e.eigenvectors
        log("arpack partial ell=2: %d Ritz pairs, polishing anyway"
            % len(w2))
    idx = np.argsort(np.abs(w2))
    w2, v2 = w2[idx], np.asarray(v2)[:, idx]
    OV2 = np.stack([mv2(v2[:, j]) for j in range(v2.shape[1])], axis=1)
    B2 = v2.conj().T @ OV2
    B2 = (B2 + B2.conj().T) / 2.0
    th2, Y2 = np.linalg.eigh(B2)
    o2 = np.argsort(np.abs(th2))
    th2, Y2 = th2[o2], Y2[:, o2]
    Wf = v2 @ Y2[:, np.abs(th2) < 1.0]
    res2 = float(max(np.max(np.abs(mv2(Wf[:, p])))
                     for p in range(Wf.shape[1])))
    ker2 = Wf.shape[1]
    G2 = Wf.conj().T @ Wf
    g2v, g2u = np.linalg.eigh(G2)
    assert g2v.min() > 1e-8, "ell=2 kernel Gram rank-deficient"
    Wf = Wf @ (g2u * g2v ** -0.5) @ g2u.conj().T
    gw2 = float(np.max(np.abs(Wf.conj().T @ Wf - np.eye(ker2))))
    a1_k2 = (ker2 == 6 and res2 <= 1e-8 and gw2 <= 1e-12)
    checks.append(("A1b ell=2 full-fiber: ker %d (6), ||C4 v|| %.1e, "
                   "Gram defect %.1e: %s" % (ker2, res2, gw2, a1_k2),
                   a1_k2))
    log("ell=2: ker %d, max residual %.1e" % (ker2, res2))

    def jl2(v):
        # the V371 deformation on the length-2 line: every line link
        # carries the full matter current (T1 + T2).L (V370 anchor).
        h = v.reshape(32, 32, 45, 45)
        y = np.zeros_like(h)
        for perm, inv in MF_PERMS:
            z = np.tensordot(C_MF, h.transpose(*perm).reshape(1440,
                                                              1440),
                             axes=([1], [0])).reshape(32, 45, 32, 45)
            y = y + z.transpose(*inv)
        return y

    JW2 = np.stack([jl2(Wf[:, p]) for p in range(ker2)], axis=-1)
    JW2m = JW2.reshape(FIBER2, ker2)
    D2 = Wf.conj().T @ JW2m
    D2 = (D2 + D2.conj().T) / 2.0
    leak2 = float(np.linalg.norm(JW2m - Wf @ D2) / np.linalg.norm(JW2m))
    dk = np.linalg.eigvalsh(D2)
    booked_dk = np.array([-96.0, -96.0, -64.0, -64.0, 0.0, 0.0])
    a1_dk = np.max(np.abs(dk - booked_dk)) <= 1e-9
    checks.append(("A1c deformation set d_k == {-96,-96,-64,-64,0,0} "
                   "exactly: max dev %.2e (<=1e-9); fiber JL leak "
                   "%.1e (<=1e-3): %s"
                   % (np.max(np.abs(dk - booked_dk)), leak2,
                      a1_dk and leak2 <= 1e-3),
                   a1_dk and leak2 <= 1e-3))
    log("ell=2 deformation d_k: %s" % np.round(dk, 6))

    # ------------------------------------------------------------------
    # T3: operator self-tests, fast coupling form == 45-loop anchored
    # form (the V371 permanent guard, tightened to 1e-14).
    # ------------------------------------------------------------------
    rngO = np.random.default_rng(0x37)

    def loop_e4(x, opJL=False):
        h = x.reshape(32, 32, 45)
        acc = np.zeros_like(h)
        for a in range(45):
            if opJL:
                lz = np.einsum("Ff,abf->abF", Ls[a], h)
                acc += (np.einsum("Aa,abF->AbF", Ts[a], lz)
                        + np.einsum("Bb,abF->aBF", Ts[a], lz))
            else:
                def lg(g):
                    return (np.einsum("Aa,abf->Abf", Ts[a], g)
                            + np.einsum("Bb,abf->aBf", Ts[a], g)
                            + np.einsum("Ff,abf->abF", Ls[a], g))
                acc += lg(lg(h))
        return acc.ravel()

    def loop_c4(x):
        h = x.reshape(32, 32, 45, 45)
        acc = np.zeros_like(h)
        for a in range(45):
            def lg(g):
                return (np.einsum("Aa,abfg->Abfg", Ts[a], g)
                        + np.einsum("Bb,abfg->aBfg", Ts[a], g)
                        + np.einsum("Pf,abfg->abPg", Ls[a], g)
                        + np.einsum("Qg,abfg->abfQ", Ls[a], g))
            acc += lg(lg(h))
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

    def relerr(a, b):
        a = np.asarray(a).ravel()
        b = np.asarray(b).ravel()
        return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)),
                                                 1e-12))

    dop = 0.0
    u4 = rngO.normal(size=46080) + 1j * rngO.normal(size=46080)
    dop = max(dop, relerr(e4_apply(u4), loop_e4(u4)),
              relerr(e4_apply(u4, opJL=True), loop_e4(u4, opJL=True)))
    u2 = rngO.normal(size=FIBER2) + 1j * rngO.normal(size=FIBER2)
    dop = max(dop, relerr(c4_np(u2), loop_c4(u2)),
              relerr(jl2(u2), loop_jl2(u2)))
    a3op = dop <= 1e-14
    checks.append(("operator self-test: V371 fast coupling form == "
                   "45-loop anchored form on e4/C4/JL2, max rel diff "
                   "%.2e (<=1e-14): %s" % (dop, a3op), a3op))
    data["operator_selftest"] = dop
    log("operator self-test: max rel diff %.2e" % dop)
    del u4, u2

    # ------------------------------------------------------------------
    # T4: pair-channel and link-pair-channel spectra, MEASURED (the
    # V100/V101 anchors); the fusion data enters only as OUTPUT.
    # ------------------------------------------------------------------
    I32, I45 = np.eye(32), np.eye(45)
    Jpair = [np.kron(t, I32) + np.kron(I32, t) for t in Ts]     # 1024
    C2pair = sum(J @ J for J in Jpair)
    evp, evpv = np.linalg.eigh(C2pair)
    del C2pair

    def chan_split(vals):
        out = {}
        for i, v in enumerate(vals):
            v = 0.0 if v < 1e-6 else float(v)
            hit = None
            for kk in out:
                if abs(kk - v) <= 1e-5 * max(1.0, abs(v)):
                    hit = kk
                    break
            if hit is None:
                out[v] = [i]
            else:
                out[hit].append(i)
        return out

    pair_census = {}
    pair_ch = {}
    for cv0, idxs in chan_split(evp).items():
        cv = round(cv0, 6)
        E = evpv[:, idxs]
        pair_ch[cv] = (E, [E.conj().T @ J @ E for J in Jpair],
                       len(idxs))
        pair_census[cv] = len(idxs)
    del Jpair

    # link-pair Casimir spectrum via dense eigh on 2025 (small); the
    # fusion census AND the channel generator blocks (the measured
    # certified blocks the T5 solvers act on) come out of this one
    # diagonalization; no 2025^2 generator stack is materialized.
    C2ll = np.zeros((2025, 2025), dtype=np.complex128)
    for a in range(45):
        Qa = np.kron(Ls[a], I45) + np.kron(I45, Ls[a])
        C2ll += Qa @ Qa
    evl, evlv = np.linalg.eigh(C2ll)
    del C2ll
    ll_census = {}
    ll_ch = {}
    for cv0, idxs in chan_split(evl).items():
        cv = round(cv0, 6)
        E = evlv[:, idxs]
        n = len(idxs)
        Eb = E.reshape(45, 45, n)              # (link1, link2, chan)
        gens = []
        for a in range(45):
            La = Ls[a]
            gL1 = np.einsum("ac,cbn->abn", La, Eb).reshape(2025, n)
            # (I (x) L_a) E: out[alpha,beta,n] = sum_g Eb[alpha,g,n]
            # L_a[beta,g]  (L is antisymmetric: index order matters)
            gL2 = np.einsum("agn,bg->abn", Eb, La).reshape(2025, n)
            t = (np.tensordot(E.conj(), gL1, axes=([0], [0]))
                 + np.tensordot(E.conj(), gL2, axes=([0], [0])))
            gens.append(t)
        ll_ch[cv] = (E, gens, n)
        ll_census[cv] = n
    # gate the certified block generators against the measured
    # Casimir: on each block sum gL^2 must equal cv * I (a silently
    # zeroed or transposed block action would pass every downstream
    # count vacuously otherwise - the V371 machinery lesson again).
    gdefect = 0.0
    for cv, (E, gens, n) in ll_ch.items():
        Cg = sum(G @ G for G in gens)
        gdefect = max(gdefect, float(np.max(np.abs(Cg - cv * np.eye(n)))))
    a_g = gdefect <= 1e-8
    checks.append(("link-pair channel generator blocks: sum gL^2 == "
                   "measured Casimir value per block, max defect "
                   "%.2e (<=1e-8): %s" % (gdefect, a_g), a_g))
    data["ll_genblock_defect"] = gdefect
    pc_ok = (pair_census.get(0.0) == 2 and pair_census.get(64.0) == 90
             and pair_census.get(96.0) == 420)
    checks.append(("A3a pair-channel census (measured): Casimir dims "
                   "%s; (0:2, 64:90, 96:420) V371-census content: %s"
                   % (pair_census, pc_ok), pc_ok))
    a3_ll = (ll_census.get(0.0) == 1 and ll_census.get(64.0) == 45
             and ll_census.get(96.0) == 210)
    checks.append(("A3b link-pair fusion (measured): I_2 = %d (V100: "
                   "1), single 45 copy %s, 210 channel %s: %s"
                   % (ll_census.get(0.0), ll_census.get(64.0) == 45,
                      ll_census.get(96.0) == 210, a3_ll), a3_ll))
    data["pair_census"] = {str(k): int(v) for k, v
                           in sorted(pair_census.items())}
    data["link_pair_census"] = {str(k): int(v) for k, v
                                in sorted(ll_census.items())}
    log("pair channels %s" % sorted(pair_census.items()))
    log("link-pair channels %s" % sorted(ll_census.items()))

    # ------------------------------------------------------------------
    # T5: chain-kernel dimensions MEASURED per chain length by exact
    # intertwiner counting: the V101 zero_modes operator
    # (J W + W Lambda^T = 0). The right action is taken in the
    # certified block structure MEASURED by T4 (link-pair channel
    # blocks, tensor a single 45 link for ell = 3), so every solve
    # runs on a small block; Hom is additive over the blocks, nothing
    # is assumed about the multiplicities and no 45^3 generator
    # stack is ever materialized.
    # ------------------------------------------------------------------
    def hom_mv(Jgens, right):
        """matrix-free normal operator of the intertwiner equation:
        W -> sum_a (J_a W + W Lambda_a^T, then J_a^dag Q + Q Lambda_a^*);
        right = ("leg", gens) single-block action or ("prod", gensA,
        gensB) tensor of two block actions applied by einsum; GPU
        matvec on the large block spaces. Returns (mv, mvX, mvX32,
        s, n, rdims): mv is the host-boundary fp64 op, mvX the
        GPU-resident fp64 op, mvX32 the complex64 GPU op (candidate
        generation on the big blocks; the RTX 3090 fp64 GEMM rate
        makes fp64 matvecs ~10 s at s ~ 1.8e7)."""
        Jnp = np.stack(Jgens)
        n = Jnp.shape[1]
        if right[0] == "leg":
            Lnp = np.stack(right[1])
            rdims = (Lnp.shape[1],)
        else:
            GL = np.stack(right[1])
            GB = np.stack(right[2])
            rdims = (GL.shape[1], GB.shape[1])
        D = int(np.prod(rdims))
        s = n * D
        big = gpu and s >= 60000

        if not big:
            JX = Jnp
            JTX = np.ascontiguousarray(Jnp.conj().transpose(0, 2, 1))
            if right[0] == "leg":
                G1 = Lnp
                G1T = np.ascontiguousarray(Lnp.transpose(0, 2, 1))
                G1c = np.ascontiguousarray(Lnp.conj())
            else:
                G1 = GL
                G2 = GB
                G1c = np.ascontiguousarray(GL.conj())
                G2c = np.ascontiguousarray(GB.conj())
            if right[0] == "leg":
                def mvX(v):
                    w = v.reshape(n, D)
                    q = np.matmul(JX, w[None]) + np.matmul(w[None], G1T)
                    res = np.matmul(JTX, q) + np.matmul(q, G1c)
                    return res.sum(0).reshape(s)
            else:
                def mvX(v):
                    w = v.reshape(n, rdims[0], rdims[1])
                    res = np.zeros_like(w)
                    for a in range(45):
                        q = np.matmul(JX[a], w.reshape(n, D)).reshape(
                            (n, rdims[0], rdims[1]))
                        q += np.einsum("jcb,ac->jab", w, G1[a])
                        q += np.einsum("jad,bd->jab", w, G2[a])
                        res += np.matmul(JTX[a], q.reshape(n, D)).reshape(
                            (n, rdims[0], rdims[1]))
                        res += np.einsum("jcb,ca->jab", q, G1c[a])
                        res += np.einsum("jad,db->jab", q, G2c[a])
                    return res.reshape(s)
            mv = mvX
            mvX32 = None
            return mv, mvX, mvX32, s, n, rdims

        # ---- GPU-resident big path (GEMM form; the cupy einsum
        # dispatch was measured ~20x below ZGEMM rate on these
        # shapes) ----
        JX = cp.asarray(Jnp)
        JTX = cp.ascontiguousarray(JX.conj().transpose(0, 2, 1))
        if right[0] == "leg":
            G1 = cp.asarray(Lnp)
            G1T = cp.ascontiguousarray(G1.transpose(0, 2, 1))
            G1c = cp.ascontiguousarray(G1.conj())

            def mvX(v):
                w = v.reshape(n, D)
                q = cp.matmul(JX, w[None]) + cp.matmul(w[None], G1T)
                res = cp.matmul(JTX, q) + cp.matmul(q, G1c)
                return res.sum(0).reshape(s)
        else:
            G1 = cp.asarray(GL)
            G2 = cp.asarray(GB)
            G1c = cp.ascontiguousarray(G1.conj())
            G2c = cp.ascontiguousarray(G2.conj())
            A, B = rdims
            nA = n * A
            G2T = cp.ascontiguousarray(G2.transpose(0, 2, 1))
            G1cT = cp.ascontiguousarray(G1c.transpose(0, 2, 1))

            def mvX(v):
                w = v.reshape(n, A, B)
                res = cp.zeros_like(w)
                wf = w.reshape(nA, B)
                wt = cp.ascontiguousarray(
                    w.transpose(1, 0, 2)).reshape(A, n * B)
                for a in range(45):
                    q = cp.matmul(JX[a], w.reshape(n, D))
                    q = q.reshape(n, A, B) + cp.matmul(
                        G1[a], wt).reshape(A, n, B).transpose(1, 0, 2)
                    q = q + cp.matmul(wf, G2T[a]).reshape(n, A, B)
                    res += cp.matmul(JTX[a],
                                     q.reshape(n, D)).reshape(n, A, B)
                    qt = cp.ascontiguousarray(
                        q.transpose(1, 0, 2)).reshape(A, n * B)
                    res += cp.matmul(G1cT[a], qt).reshape(
                        A, n, B).transpose(1, 0, 2)
                    res += cp.matmul(q.reshape(nA, B),
                                     G2c[a]).reshape(n, A, B)
                return res.reshape(s)

        def mv(v):
            return cp.asnumpy(mvX(cp.asarray(v)))

        # complex64 twin (candidate generation): same code path,
        # cast operator stack.
        JX32 = JX.astype(cp.complex64)
        JTX32 = JTX.astype(cp.complex64)
        if right[0] == "leg":
            G132 = G1.astype(cp.complex64)
            G1T32 = G1T.astype(cp.complex64)
            G1c32 = G1c.astype(cp.complex64)

            def mvX32(v):
                w = v.reshape(n, D)
                q = cp.matmul(JX32, w[None]) + cp.matmul(w[None], G1T32)
                res = cp.matmul(JTX32, q) + cp.matmul(q, G1c32)
                return res.sum(0).reshape(s)
        else:
            G132 = G1.astype(cp.complex64)
            G232 = G2.astype(cp.complex64)
            G1c32 = G1c.astype(cp.complex64)
            G2c32 = G2c.astype(cp.complex64)
            G2T32 = G232.transpose(0, 2, 1)
            G1cT32 = G1c32.transpose(0, 2, 1)

            def mvX32(v):
                w = v.reshape(n, A, B)
                res = cp.zeros_like(w)
                wf = w.reshape(nA, B)
                wt = cp.ascontiguousarray(
                    w.transpose(1, 0, 2)).reshape(A, n * B)
                for a in range(45):
                    q = cp.matmul(JX32[a], w.reshape(n, D))
                    q = q.reshape(n, A, B) + cp.matmul(
                        G132[a], wt).reshape(A, n, B).transpose(1, 0, 2)
                    q = q + cp.matmul(wf, G2T32[a]).reshape(n, A, B)
                    res += cp.matmul(JTX32[a],
                                     q.reshape(n, D)).reshape(n, A, B)
                    qt = cp.ascontiguousarray(
                        q.transpose(1, 0, 2)).reshape(A, n * B)
                    res += cp.matmul(G1cT32[a], qt).reshape(
                        A, n, B).transpose(1, 0, 2)
                    res += cp.matmul(q.reshape(nA, B),
                                     G2c32[a]).reshape(n, A, B)
                return res.reshape(s)
        return mv, mvX, mvX32, s, n, rdims

    def hom_null(Jgens, right, maxvec=40):
        """the COMPLETE kernel by deflation rounds. Small blocks
        (s < 6e4): a CAPPED ARPACK pass (k = 4, ncv = 12, 300 matvec
        cap) on the shifted operator M - 32I over the orthogonal
        complement of the modes found so far, so the intertwiner
        kernel is the spectral extreme (the shift is half the
        measured Casimir gap of 64; ARPACK's relative criterion can
        never converge lambda = 0, and unshifted "SA" fills the
        slots with junk and misses the kernel - measured). Big
        blocks: complex64 Lanczos candidates on the same shifted
        operator, polished exactly in fp64 by a Rayleigh-Krylov
        step (memory + fp64-rate bound on the 3090s; see the branch
        comment). Both routes accept on the TRUE M in fp64 only:
        Rayleigh quotient |w.Mw| <= 1e-4 and residual ||M w|| <=
        1e-5, both orders of magnitude below the gap, so no nonzero
        mode can fake acceptance. A round accepting nothing exhausts
        the kernel; deflation + the filters make the count exact,
        not best effort, and the cap bounds the no-kernel case."""
        mv, mvX64, mvX32, s, n, rdims = hom_mv(Jgens, right)
        found = []
        if s >= 60000 and gpu:
            # release the previous block's cached segments (the pool
            # holds reserved memory across shapes; a fresh 3.4 GB
            # allocation on top of 24 GB reserved is the measured
            # OOM mode)
            cp.get_default_memory_pool().free_all_blocks()
        while len(found) < maxvec:
            if s >= 60000 and gpu:
                # big-block route (memory AND fp64-rate bound): the
                # host scipy path needs workd = 3*ncv*s complex128
                # (~10 GB at s ~ 1.8e7; the box has ~9 GB free while
                # other services hold the rest - the first run
                # swapped and died after 15 h), and the 3090 fp64
                # GEMM rate makes fp64 matvecs ~10 s at that size.
                # Candidates therefore come from restarted Lanczos
                # (window 12, full reorthogonalization) in
                # COMPLEX64 on the shifted M - 32I (the RTX 3090 runs
                # this ~20x faster; the shifted kernel cluster at
                # theta = -32 is resolved far below the fp32 noise
                # floor of ~1e-4 for this spectrum, and junk sits at
                # +32), then each candidate batch is polished
                # EXACTLY in fp64 by one Rayleigh-Krylov step on
                # B = orth([W, (M - 32I) W]) - images come from
                # {M W, M^2 W}, so polishing costs 2k matvecs.
                # Acceptance is unchanged and computed on the TRUE
                # M in fp64: Rayleigh |w.Mw| <= 1e-4 AND residual
                # ||M w|| <= 1e-5. A converged shifted junk minimum
                # (theta > -16) exhausts the block.
                Fx64 = [cp.asarray(f, dtype=cp.complex128)
                        for f in found]
                Fx32 = [u.astype(cp.complex64) for u in Fx64]

                def deflX(v, _Fx=Fx32):
                    for u in _Fx:
                        v -= cp.vdot(u, v) * u
                    return v

                def newq():
                    q = (cp.random.standard_normal(s, dtype=cp.float32)
                         + 1j * cp.random.standard_normal(
                             s, dtype=cp.float32))
                    q = deflX(q)
                    nrm = float(cp.linalg.norm(q))
                    return q / nrm if nrm > 1e-12 else None

                q = newq()
                if q is None:
                    break
                steps = 0
                cands = []
                junkstreak = 0
                phantom = 0
                got = 0
                exhausted = False

                def polish(_s=s):
                    """fp64 Rayleigh-Krylov polish of the fp32
                    candidate batch: Wc, MV = M Wc, MMV = M^2 Wc;
                    B = [Wc, (M - 32I) Wc] with MB = [MV, (M - 32I)
                    MV]; accepted only if the TRUE-M fp64 filters
                    pass. Batches of 4 keep the fp64 stack under
                    ~7 GB at s ~ 1.8e7 (k = 12 was a measured OOM
                    on the 24 GB card)."""
                    nonlocal got, phantom
                    k = len(cands)
                    Wc = cp.linalg.qr(cp.column_stack(cands)
                                      .astype(cp.complex128))[0]
                    MV = cp.column_stack([mvX64(Wc[:, j])
                                          for j in range(k)])
                    MMV = cp.column_stack([mvX64(MV[:, j])
                                           for j in range(k)])
                    B = cp.hstack([Wc, MV - 32.0 * Wc])
                    MB = cp.hstack([MV, MMV - 32.0 * MV])
                    del Wc, MV, MMV
                    Gb = B.conj().T @ MB
                    _, R = cp.linalg.qr(B)
                    Ri = np.linalg.inv(cp.asnumpy(R))
                    H = Ri.conj().T @ cp.asnumpy(Gb) @ Ri
                    del Gb
                    H = (H + H.conj().T) / 2.0
                    _, U = np.linalg.eigh(H)
                    acc = 0
                    for p in range(2 * k):
                        c = Ri @ U[:, p]
                        w = B @ cp.asarray(c)
                        nw = float(cp.linalg.norm(w))
                        if nw < 1e-8:
                            continue
                        rw = MB @ cp.asarray(c)
                        w = w / nw
                        rw = rw / nw
                        if float(cp.linalg.norm(rw)) > 1e-5:
                            continue
                        if float(abs(np.real(np.vdot(
                                cp.asnumpy(w), cp.asnumpy(rw))))) > 1e-4:
                            continue
                        for u in Fx64:
                            w -= u * cp.vdot(u, w)
                        if float(cp.linalg.norm(w)) < 1e-6:
                            continue
                        w = w / cp.linalg.norm(w)
                        found.append(cp.asnumpy(w))
                        Fx64.append(w)
                        Fx32.append(w.astype(cp.complex64))
                        got += 1
                        acc += 1
                        if len(found) >= maxvec:
                            break
                    del B, MB
                    cands.clear()
                    phantom = 0 if acc else phantom + 1
                    log("big block s=%d: polish batch k=%d accepted %d"
                        " (found %d, steps %d)"
                        % (_s, k, acc, len(found), steps))

                # exhaustion rules: (i) 3 consecutive windows whose
                # lowest shifted mode sits ABOVE -16 (the fp64
                # residual-estimate test cannot be transplanted to
                # fp32: the estimate stalls at the fp32 noise floor);
                # (ii) two consecutive candidate batches whose fp64
                # filters reject everything = fp32 phantoms, the
                # kernel is spent. (i) is what happens on junk
                # blocks; (ii) caps the phantom hunt after the
                # kernel was already deflated (measured: a full 240
                # step round otherwise). The step cap itself is
                # 240 for a first pass and 72 once something has
                # been accepted: true kernel vectors converge from
                # a random start in well under 40 steps for this
                # spectral gap, and the fp64 polish matvecs (~50 s
                # at s ~ 1.8e7) make over-long re-rounds expensive.
                cap_r = 160 if not found else 72
                while steps < cap_r:
                    V, al, be = [], [], []
                    v, beta = q, cp.float32(0.0)
                    while len(V) < 12 and steps < cap_r:
                        u = mvX32(v) - cp.float32(32.0) * v
                        steps += 1
                        if V:
                            u -= beta * V[-1]
                        for vv in V:
                            u -= cp.vdot(vv, u) * vv
                        a = float(cp.real(cp.vdot(v, u)))
                        u -= a * v
                        for vv in V:
                            u -= cp.vdot(vv, u) * vv
                        b = float(cp.linalg.norm(u))
                        al.append(a)
                        be.append(b)
                        V.append(v)
                        if b < 1e-6:
                            break
                        beta = b
                        v = u / b
                    mm = len(V)
                    Tm = np.diag(np.asarray(al))
                    if mm > 1:
                        bb = np.asarray(be[:-1])
                        Tm += np.diag(bb, 1) + np.diag(bb, -1)
                    th, Y = np.linalg.eigh(Tm)
                    y = Y[:, 0]
                    w = complex(y[0]) * V[0]
                    for i in range(1, mm):
                        w = w + complex(y[i]) * V[i]
                    w = w / cp.linalg.norm(w)
                    if th[0] < -16.0:
                        cands.append(w)
                        junkstreak = 0
                        if len(cands) >= 4:
                            polish()
                            if phantom >= 2:
                                break
                    else:
                        junkstreak += 1
                        if junkstreak >= 3:
                            exhausted = True
                            break
                    q = deflX(w)
                    nrm = float(cp.linalg.norm(q))
                    if nrm < 1e-8:
                        break
                    q = q / nrm
                del V, q
                if cands:
                    polish()
                if not got:
                    if not exhausted:
                        log("big block s=%d: solver cap reached with "
                            "no acceptance (declared complete on the "
                            "exhaustion logic; recheck if counts "
                            "disagree with the second route)" % s)
                    break
                continue
            if found:
                Qm = np.column_stack(found)

                def mv_d(v, Qm=Qm):
                    v = v - Qm @ (Qm.conj().T @ v)
                    out = mv(v)
                    return out - Qm @ (Qm.conj().T @ out)
            else:
                mv_d = mv
            if s <= 8:
                Mm = np.stack([mv_d(np.eye(s, dtype=np.complex128)[:, j])
                               for j in range(s)], axis=1)
                vals, vecs = np.linalg.eigh(Mm)
            else:
                # M - 32I: the kernel cluster becomes the spectral
                # EXTREME (the shift is half the measured Casimir gap
                # of 64). ARPACK's relative stopping criterion can
                # never declare lambda = 0 converged, and unshifted
                # "SA" was measured to fill all k slots with the
                # 64+ junk, missing the kernel entirely. The shift
                # controls only which Ritz pairs ARPACK returns; the
                # acceptance filter below still tests the TRUE M.
                def mv_s(v, _f=mv_d):
                    return _f(v) - 32.0 * v
                try:
                    vals, vecs = spla.eigsh(
                        spla.LinearOperator((s, s), matvec=mv_s,
                                            dtype=np.complex128),
                        k=min(4, s - 2), ncv=min(12, s - 1),
                        which="SA", tol=1e-6, maxiter=300)
                except spla.ArpackNoConvergence as anc:
                    vals = np.real(anc.eigenvalues)
                    vecs = anc.eigenvectors
            got = 0
            for j in range(len(vals)):
                w = np.asarray(vecs[:, j]).astype(np.complex128)
                rw = mv(w)
                if np.linalg.norm(rw) > 1e-5:
                    continue
                if abs(float(np.real(np.vdot(w, rw)))) > 1e-4:
                    continue
                for u in found:
                    w -= u * np.vdot(u, w)
                nw = np.linalg.norm(w)
                if nw < 1e-8:
                    continue
                found.append(w / nw)
                got += 1
            if not got:
                break
        return [f.reshape((n,) + rdims) for f in found]

    Jr0 = pair_ch[0.0][1]

    # per-block kill-safety: every measured chain block is cached to
    # its own npz (dir results/v373_ck, script-sha guarded - stale
    # dirs from an edited script are wiped). A killed run resumes
    # the T5 census at block granularity: the 3 h lost twice to
    # coarse-grained restarts cannot happen again.
    CKD = os.path.join(HERE, "results", "v373_ck")
    mysha = hashlib.sha256(
        open(os.path.abspath(__file__), "rb").read()).hexdigest()
    if os.path.isdir(CKD):
        shf = os.path.join(CKD, ".sha")
        if (not os.path.exists(shf)
                or open(shf).read().strip() != mysha):
            shutil.rmtree(CKD)
    os.makedirs(CKD, exist_ok=True)
    with open(os.path.join(CKD, ".sha"), "w") as fh:
        fh.write(mysha)

    def _fmt(x):
        return "%.6g" % x

    def _ck_sig(Jgens, right):
        """content signature of the operator stack a block was solved with.
        The fusion channels carry degenerate eigenspaces (a 20-fold at
        cv 36), whose eigenbasis LAPACK/cuSOLVER may rotate freely; two
        processes build the same channels from last-bit-different coupling
        bytes (measured 2026-11-05: the fold's subprocess cache-hit the
        previous process's vectors, the span check collapsed to 0.005).
        A cache hit is only valid if the generator bytes reproduce
        EXACTLY; otherwise the block is re-solved in this process."""
        h = hashlib.sha256()
        h.update(np.ascontiguousarray(np.stack(Jgens),
                                      dtype=np.complex128).tobytes())
        h.update(np.ascontiguousarray(np.stack(right[1]),
                                      dtype=np.complex128).tobytes())
        if right[0] == "prod":
            h.update(np.ascontiguousarray(np.stack(right[2]),
                                          dtype=np.complex128).tobytes())
        return h.hexdigest()[:16]

    def _ck_get(stage, key, sig):
        fn = os.path.join(CKD, "%s_%s.npz" % (stage, key))
        if not os.path.exists(fn):
            return None
        with np.load(fn) as z:
            if "sig" not in z.files or str(z["sig"]) != sig:
                return None
            return [z["v%d" % i] for i in range(len(z.files) - 1)]

    def _ck_put(stage, key, vecs, sig):
        np.savez(os.path.join(CKD, "%s_%s.npz" % (stage, key)),
                 sig=sig,
                 **{"v%d" % i: np.asarray(v, dtype=np.complex128)
                    for i, v in enumerate(vecs)})

    # booked chain convention: a flux line of length ell carries a
    # chain of ell + 1 links; the chain-kernel counts
    # dim Hom(pair, 45^(ell+1)) are MEASURED per chain length over
    # the measured channel blocks (Hom is additive over blocks). The
    # 1-link Hom is the contact-fiber cross-check against T1, not a
    # booked chain length.
    Nl1, Uf = {}, {}
    for cv, (E, Jr, nR) in pair_ch.items():
        sig = _ck_sig(Jr, ("leg", Ls))
        v = _ck_get("n1", _fmt(cv), sig)
        if v is None:
            v = hom_null(Jr, ("leg", Ls))
            _ck_put("n1", _fmt(cv), v, sig)
        Nl1[cv] = len(v)
        if v:
            Uf[cv] = v
    ker1_hom = sum(Nl1.values())
    a3_k1 = ker1_hom == ker1
    checks.append(("A3c Hom(pair, 45) block counts: kernel %d == T1"
                   " full-fiber Lanczos ker %d (contact fiber): %s"
                   % (ker1_hom, ker1, a3_k1), a3_k1))

    # chain ell=1 (two links, the V99/V101 three-site chain):
    Nc1, U2 = {}, {}
    c1pat = {}
    for cv, (E, Jr, nR) in pair_ch.items():
        t = 0
        for lkv, (El, gL, nL) in ll_ch.items():
            sig = _ck_sig(Jr, ("leg", gL))
            v = _ck_get("e1", "%s|%s" % (_fmt(cv), _fmt(lkv)), sig)
            if v is None:
                v = hom_null(Jr, ("leg", gL))
                _ck_put("e1", "%s|%s" % (_fmt(cv), _fmt(lkv)), v, sig)
            if v:
                U2[(cv, lkv)] = v
                c1pat["%s|%s" % (round(cv, 4),
                                 round(lkv, 4))] = len(v)
            t += len(v)
        Nc1[cv] = t
    ker_c1 = sum(Nc1.values())
    pat_ok = (len(c1pat) == 3 and sorted(c1pat.values()) == [2, 2, 2])
    a3_k2 = (ker_c1 == ker2 and pat_ok)
    checks.append(("A3d chain ell=1 (45x45): measured kernel %d == "
                   "V371-fiber Lanczos ker 6 AND the V100 shared-"
                   "channel content r = 1, 45, 210 with pair-side "
                   "multiplicities 2 x (single ll copies): block "
                   "pattern %s: %s" % (ker_c1,
                                       sorted(c1pat.items()),
                                       a3_k2), a3_k2))
    # the published V99/V100/V101 glue identity: 24 = m_45^2 x 6
    # (V99: dim ker = m45^2 * Nsum), with the 24-state glued
    # kernel gram-rank 24 and the V101 top degeneracy 6.
    v100d = json_art("V100_vertex_correction")
    v101d = json_art("V101_bond_current_ground_state")
    m45 = pair_census.get(64.0, 0) // 45
    a_v101 = (ker_c1 * m45 ** 2 == int(v100d["data"]["n_states"])
              and int(v100d["data"]["gram_rank"]) == 24
              and int(v101d["data"]["gram_rank"]) == 24
              and int(v101d["data"]["degeneracy"]) == 6)
    checks.append(("A3d2 published-kernel anchor: measured Hom %d "
                   "x measured m_45^2 %d == V100 glued kernel %s "
                   "states, gram ranks 24/24 (V100/V101), V101 "
                   "top degeneracy 6: %s"
                   % (ker_c1, m45 ** 2, v100d["data"]["n_states"],
                      a_v101), a_v101))

    # chain ell=2 (three links): blocks (lkv x single 45 link)
    Nc2, U3 = {}, {}
    for cv, (E, Jr, nR) in pair_ch.items():
        t = 0
        for lkv, (El, gL, nL) in ll_ch.items():
            log("chain ell=2 block cv %s lkv %s" % (round(cv, 3),
                                                    round(lkv, 3)))
            sig = _ck_sig(Jr, ("prod", gL, Ls))
            v = _ck_get("n2", "%s|%s" % (_fmt(cv), _fmt(lkv)), sig)
            if v is None:
                v = hom_null(Jr, ("prod", gL, Ls))
                _ck_put("n2", "%s|%s" % (_fmt(cv), _fmt(lkv)), v, sig)
            t += len(v)
            if v:
                U3[(cv, lkv)] = v
        Nc2[cv] = t
    ker_c2 = sum(Nc2.values())
    # I_3 = dim(45^3)^G measured from the singlet pair channels:
    # Hom(singlet channels, 45^3) = (2 copies) x I_3; V100 books 1.
    n_sing3 = Nc2.get(0.0, 0)
    nsing = pair_census.get(0.0, 1)
    a3_i3 = n_sing3 == nsing
    checks.append(("A3e I_3 = dim(45^3)^G: singlet-channel Hom "
                   "count %d == %d copies x V100's I_3 = 1: %s"
                   % (n_sing3, nsing, a3_i3), a3_i3))

    # chain ell=3 (four links, the booked prediction output): the
    # SINGLET pair channel is measured over all (lkv x lkv') blocks;
    # non-singlet sectors reach block sizes n x 945 x 945 beyond the
    # honest runtime budget and are DECLARED UNMEASURED (booked
    # scope limit, recorded here and in run notes, never silent).
    Nc3s = 0
    U4 = {}
    for lkvA, (EA, gA, nA) in ll_ch.items():
        for lkvB, (EB, gB, nB) in ll_ch.items():
            log("chain ell=3 singlet block %s x %s"
                % (round(lkvA, 3), round(lkvB, 3)))
            sig = _ck_sig(Jr0, ("prod", gA, gB))
            v = _ck_get("n3s", "%s|%s" % (_fmt(lkvA), _fmt(lkvB)), sig)
            if v is None:
                v = hom_null(Jr0, ("prod", gA, gB))
                _ck_put("n3s", "%s|%s" % (_fmt(lkvA), _fmt(lkvB)), v, sig)
            Nc3s += len(v)
            if v:
                U4[(lkvA, lkvB)] = v
    i4_meas = Nc3s // nsing if nsing else -1
    v100_i4 = int(round(float(v100d["data"]["I"]["4"][0])))
    a_i4 = i4_meas == v100_i4
    checks.append(("A3f I_4 = dim(45^4)^G singlet-sector "
                   "measurement: %d Hom copies / %d singlet "
                   "copies = I_4 = %d == V100's published 6 "
                   "(block sum over all lkv x lkv' channels of "
                   "45x45; non-singlet pair sectors DECLARED "
                   "UNMEASURED): %s"
                   % (Nc3s, nsing, i4_meas, a_i4), a_i4))

    data["chain_kernel_dims"] = {"links1_contact": int(ker1_hom),
                                 "chain_ell1_links2": int(ker_c1),
                                 "chain_ell2_links3": int(ker_c2),
                                 "chain_ell3_links4_singlet":
                                     int(Nc3s)}
    data["chain_hom_counts"] = {
        "contact": {str(k): int(v) for k, v in sorted(Nl1.items())},
        "chain_ell1": {str(k): int(v) for k, v
                       in sorted(Nc1.items())},
        "chain_ell2": {str(k): int(v) for k, v
                       in sorted(Nc2.items())}}
    data["chain_ell1_block_pattern"] = c1pat
    data["I3_measured"] = n_sing3 // nsing if nsing else -1
    data["I4_measured_singlet_sector"] = i4_meas
    log("measured chain kernels: contact 1-link %d, ell=1 %d, "
        "ell=2 %d, ell=3 singlet-sector %d (I_4 = %d vs V100 %d)"
        % (ker1_hom, ker_c1, ker_c2, Nc3s, i4_meas, v100_i4))

    # ------------------------------------------------------------------
    # T6: string Casimir coefficient, MEASURED from the adjoint rep
    # (the booked "not assumed" term of the generator).
    # ------------------------------------------------------------------
    kv = np.linalg.eigvalsh((K_F + K_F.conj().T) / 2.0)
    c_spread = float(kv.max() - kv.min())
    c_F = float(np.mean(kv))
    a_cf = c_spread <= 1e-8
    checks.append(("string Casimir per link: sum_a L_a L_a "
                   "spectrum flat at %.6f (spread %.1e <= 1e-8): %s"
                   % (c_F, c_spread, a_cf), a_cf))
    data["string_casimir_per_link"] = c_F


    # ------------------------------------------------------------------
    # T7: A2 -- the V368 torus kinematics verbatim (seed 0x51): the
    # free-step min gap must reproduce 2.017083e-05 to 1e-9.
    # ------------------------------------------------------------------
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
    gaps_free = []
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
        gaps_free.append(float(circmin(comp, cont).min()))
    free_min = min(gaps_free)
    v368 = json_art("V368_torus_rigidity_window")
    v368_min = float(v368["data"]["gaps_min"])
    a2 = abs(free_min - v368_min) < 1e-9
    checks.append(("A2 free min gap %.9f reproduces V368 %.9f to 1e-9:"
                   " %s" % (free_min, v368_min, a2), a2))
    data["free_min_repro"] = free_min
    gsort = np.argsort(gaps_free)
    data["worstK_index"] = int(gsort[0])
    data["worstK_tri"] = [int(x) for x in tri[int(gsort[0])]]
    data["gap_order"] = sorted(int(i) for i in gsort[:12])

    # ------------------------------------------------------------------
    # T8: real-space one-particle bond propagation P1 on the N^3
    # torus (the pair-hopping machinery of V374): its Bloch image
    # must equal SC.bond_product(q) to 1e-12 at sample momenta.
    # ------------------------------------------------------------------
    NSUB = SC.NSUB

    def one_particle_prop(N):
        cells = N ** 3
        dim = NSUB * cells
        S = np.zeros((dim, dim), dtype=np.complex128)
        for c in range(cells):
            for e in range(NSUB):
                S[NSUB * c + e, NSUB * c + e] = 1.0

        def cellv(ci):
            return ((ci[0] % N) * N + (ci[1] % N)) * N + (ci[2] % N)

        for (even, axis) in SC.STAGES:
            cc, ss = math.cos(SC.G10), math.sin(SC.G10)
            M = np.zeros((dim, dim), dtype=np.complex128)
            bit = 1 << axis
            for c in range(cells):
                i, r = divmod(c, N * N)
                j, l = divmod(r, N)
                ci = [i, j, l]
                for e in range(NSUB):
                    if (e & bit) != (0 if even else bit):
                        continue
                    o = e ^ bit
                    M[NSUB * c + e, NSUB * c + e] += cc
                    M[NSUB * c + o, NSUB * c + o] += cc
                    if even:
                        M[NSUB * c + e, NSUB * c + o] += -1j * ss
                        M[NSUB * c + o, NSUB * c + e] += -1j * ss
                    else:
                        # stage_matrix: for (e & bit) == 0 the entry
                        # S[e,o] carries e^{-i q_axis}: real space
                        # psi[o] read from cell c + e_axis; the S[o,e]
                        # entry carries e^{+i q_axis}: psi[e] from
                        # cell c - e_axis.
                        cp1 = list(ci)
                        cp1[axis] += 1
                        cm1 = list(ci)
                        cm1[axis] -= 1
                        M[NSUB * c + e, NSUB * cellv(cp1) + o] \
                            += -1j * ss
                        M[NSUB * c + o, NSUB * cellv(cm1) + e] \
                            += -1j * ss
            S = M @ S
        return S

    P1 = one_particle_prop(GRID_N)
    un = float(np.max(np.abs(P1.conj().T @ P1 - np.eye(len(P1)))))

    def bloch_of(P, q):
        N = GRID_N
        B = np.zeros((NSUB, NSUB), dtype=np.complex128)
        for c in range(N ** 3):
            i, r = divmod(c, N * N)
            j, l = divmod(r, N)
            ph = np.exp(1j * (q[0] * i + q[1] * j + q[2] * l))
            B += ph * P[0:NSUB, NSUB * c:NSUB * c + NSUB]
        return B

    bdef = 0.0
    # sample momenta must lie on the torus reciprocal lattice
    # (multiples of 2 pi / N) or the wrapped Bloch sums alias
    for q in [(0.0, 0.0, 0.0), (2 * np.pi * 3 / GRID_N, 0.0, 0.0),
              (2 * np.pi * 1 / GRID_N, 2 * np.pi * 2 / GRID_N, 0.0),
              (2 * np.pi * 3 / GRID_N, 0.0, 2 * np.pi * 1 / GRID_N)]:
        got = bloch_of(P1, q)
        ref = SC.bond_product(q)
        bdef = max(bdef, float(np.max(np.abs(got - ref))))
    a_p1 = un <= 1e-12 and bdef <= 1e-12
    checks.append(("P1 real-space propagation: unitarity defect "
                   "%.1e, Bloch vs SC.bond_product %.1e (both <= "
                   "1e-12): %s" % (un, bdef, a_p1), a_p1))
    log("P1: unitarity %.1e, bloch defect %.1e" % (un, bdef))

    def dgeo(tri_v, N):
        return sum(min(v, N - v) for v in tri_v)

    tri2 = np.array([(a, b, c) for a in range(2) for b in range(2)
                     for c in range(2) if (a, b, c) != (0, 0, 0)])
    shell2 = {}
    for t in tri2:
        d = dgeo(tuple(int(x) for x in t), 2)
        shell2[d] = shell2.get(d, 0) + 1
    data["n2_shell_counts"] = shell2
    log("N=2 separation shells (d_geo): %s" % shell2)

    # ------------------------------------------------------------------
    # T9: the kernel basis carried forward - the T5 solver tensors
    # W (pair channel x 45^ell), labelled by pair channel (the V101
    # build_kernel pattern: the channel eigenspaces are the
    # co-intertwiners); deformation spectra per ell in this basis
    # must reproduce -64 (ell=1) and the {-96,-64,0} set
    # (ell=2); the ell=2 span must equal the V371 fiber
    # kernel (overlap gate); a second solver route on the singlet
    # channel confirms the chain counts.
    # ------------------------------------------------------------------
    basis1 = [(cv, w) for cv, vs in Uf.items() for w in vs]
    basis2 = []
    for (cv, lkv), vm2 in U2.items():
        El, gL, nL = ll_ch[lkv]
        Ell = El.reshape(45, 45, nL)
        ElT = np.ascontiguousarray(Ell.reshape(2025, nL).T)
        for u in vm2:
            basis2.append((cv, (u @ ElT).reshape(-1, 45, 45)))
    basis3 = []
    for (cv, lkv), vm3 in U3.items():
        El, gL, nL = ll_ch[lkv]
        Ell = El.reshape(45, 45, nL)
        EllT = np.ascontiguousarray(Ell.reshape(2025, nL).T)
        for u in vm3:
            # M[j,a,b,g] = sum_k W[j,k,g] E[a,b,k]: the third link
            # stays a separate fiber axis (a label collision with a
            # block-basis index silently dropped it and broke the
            # 3-link deformation spectrum)
            M = (np.ascontiguousarray(u.transpose(0, 2, 1))
                 .reshape(-1, nL) @ EllT)
            basis3.append((cv, M.reshape(u.shape[0], 45, 45, 45)
                           .transpose(0, 2, 3, 1)))
    basis4 = []
    for (lkvA, lkvB), vm4 in U4.items():
        EA, gA, nA = ll_ch[lkvA]
        EB, gB, nB = ll_ch[lkvB]
        EAb = EA.reshape(45, 45, nA)
        EBb = EB.reshape(45, 45, nB)
        for u in vm4:
            # GEMM-form contraction of "jkl,abk,cdn->jabcd" (the
            # unoptimized einsum dispatch on the (945,945) blocks
            # was measured to stall for over an hour per tensor)
            basis4.append((0.0, np.tensordot(
                np.tensordot(u, EAb, axes=([1], [2])),
                EBb, axes=([1], [2]))))
    log("kernel-tensor counts: contact %d, ell1 %d, ell2 %d, "
        "ell3-singlet %d"
        % (len(basis1), len(basis2), len(basis3), len(basis4)))
    bc_ok = (len(basis1) == ker1_hom and len(basis2) == ker_c1
             and len(basis3) == ker_c2 and len(basis4) == Nc3s)
    checks.append(("chain-kernel tensors carried forward match the "
                   "measured dims (%d,%d,%d,%d): %s"
                   % (len(basis1), len(basis2), len(basis3),
                      len(basis4), bc_ok), bc_ok))

    def jl_on_factored(cv, M):
        Jr = pair_ch[cv][1]
        links = M.ndim - 1
        acc = np.zeros_like(M)
        for a in range(45):
            t = np.tensordot(Jr[a], M, axes=([1], [0]))
            for ln in range(links):
                tm = np.moveaxis(t, ln + 1, 0)
                acc += np.moveaxis(
                    np.tensordot(Ls[a], tm.reshape(45, -1),
                                 axes=([1], [0])).reshape(tm.shape),
                    0, ln + 1)
        return acc

    def decomp(basis):
        # the deformation operator is block-diagonal over pair
        # channels (each sector carries its own Jr on its own W
        # dimension): decompose per sector and concatenate the
        # eigenvalue multisets.
        groups = {}
        for cv, M in basis:
            groups.setdefault(cv, []).append(M)
        out = []
        for cv, Ms in groups.items():
            n = len(Ms)
            flat = [np.asarray(M).ravel() for M in Ms]
            Bm = np.zeros((n, n), dtype=np.complex128)
            Am = np.zeros((n, n), dtype=np.complex128)
            ys = [jl_on_factored(cv, M) for M in Ms]
            for p in range(n):
                for q in range(n):
                    Bm[p, q] = np.vdot(flat[p], flat[q])
                    Am[p, q] = np.vdot(flat[q], ys[p])
            Bm = (Bm + Bm.conj().T) / 2.0
            Am = (Am + Am.conj().T) / 2.0
            bv, bu = np.linalg.eigh(Bm)
            sel = bv > 1e-9 * bv.max()
            S = bu[:, sel] * bv[sel] ** -0.5
            Dg = S.conj().T @ Am @ S
            Dg = (Dg + Dg.conj().T) / 2.0
            out.append(np.linalg.eigvalsh(Dg))
        return np.concatenate(out) if out else np.array([])

    dk1 = decomp(basis1)
    dk2 = np.sort(decomp(basis2))
    dk3 = np.sort(decomp(basis3))
    dk4 = np.sort(decomp(basis4)) if basis4 else np.array([])
    log("kernel-basis deformation spectra: contact %s chain-ell1 %s "
        "chain-ell2 %s chain-ell3-singlet %s"
        % (np.round(dk1, 6), np.round(dk2, 6), np.round(dk3, 6),
           np.round(dk4, 6)))
    dk1_ok = np.max(np.abs(dk1 + 64.0)) <= 1e-6
    checks.append(("kernel-basis deformation on the contact fiber "
                   "reproduces the V370 -64 anchor: max dev %.2e "
                   "(<=1e-6): %s" % (np.max(np.abs(dk1 + 64.0)),
                                     dk1_ok), dk1_ok))
    dk2_ok = np.max(np.abs(dk2 - booked_dk)) <= 1e-6
    checks.append(("kernel-basis deformation on the 2-link chain "
                   "reproduces the fiber set {-96,-96,-64,-64,0,0}: "
                   "max dev %.2e (<=1e-6): %s"
                   % (np.max(np.abs(dk2 - booked_dk)), dk2_ok),
                   dk2_ok))

    # 2-link chain span consistency: the block-built basis must span
    # the SAME kernel V371's full-fiber Lanczos measured.
    ov = 0.0
    if basis2 and ker2:
        Fs = np.stack([np.tensordot(pair_ch[c][0], M,
                                    axes=([1], [0])).ravel()
                       for (c, M) in basis2], axis=1)
        fv, fu = np.linalg.eigh(Fs.conj().T @ Fs)
        Fs = Fs @ (fu * np.maximum(fv, 1e-12) ** -0.5) @ fu.conj().T
        O = Wf.conj().T @ Fs
        sv = np.linalg.svd(O, compute_uv=False)
        ov = float(np.min(sv))
    a_ov = ov >= 1.0 - 1e-8
    checks.append(("2-link chain block-span == V371 fiber kernel "
                   "span (min overlap singular value %.10f >= "
                   "1-1e-8): %s" % (ov, a_ov), a_ov))

    # independent cross-check of the chain counts with a second
    # solver route on the singlet pair channel: a dense right action
    # (materialized 2025^2 link-pair generators) versus the block
    # matvecs used above.
    Lnp2 = [np.kron(Ls[a], I45) + np.kron(I45, Ls[a])
            for a in range(45)]
    vmf2 = hom_null(Jr0, ("leg", Lnp2))
    vmf3 = hom_null(Jr0, ("prod", Lnp2, Ls))
    a_fs = (len(vmf2) == Nc1.get(0.0, -1)
            and len(vmf3) == Nc2.get(0.0, -1))
    checks.append(("second solver route (materialized link-pair gens)"
                   ", singlet channel: 2-link %d == block %d, 3-link "
                   "%d == block %d: %s"
                   % (len(vmf2), Nc1.get(0.0, -1), len(vmf3),
                      Nc2.get(0.0, -1), a_fs), a_fs))
    data["singlet_direct_counts"] = {"links2": int(len(vmf2)),
                                     "links3": int(len(vmf3))}
    data["dk_factored"] = {"contact": [float(x) for x in dk1],
                           "chain_ell1": [float(x) for x in dk2],
                           "chain_ell2": [float(x) for x in dk3],
                           "chain_ell3_singlet":
                               [float(x) for x in dk4]}
    data["dk_fiber_ell2"] = [float(x) for x in dk]

    # cache the channel machinery for V374
    for ell, basis in ((1, basis1), (2, basis2), (3, basis3),
                       (4, basis4)):
        arrs = {}
        metas = []
        for i, (cv, M) in enumerate(basis):
            arrs["m%d" % i] = np.asarray(M, dtype=np.complex128)
            metas.append(cv)
        np.savez(os.path.join(HERE, "results",
                              "v373_channels_ell%d.npz" % ell),
                 metas=np.asarray(metas, dtype=np.float64), **arrs)
    np.savez(os.path.join(HERE, "results", "v373_channels.npz"),
             P1=P1,
             chain_kernel_dims=np.asarray([ker1_hom, ker_c1, ker_c2,
                                           Nc3s]),
             string_casimir=np.asarray([c_F]),
             dk1=np.asarray(dk1), dk2=np.asarray(dk),
             dk3=np.asarray(dk3), dk4=np.asarray(dk4))

    # ------------------------------------------------------------------
    ok = all(c[1] for c in checks)
    if ok:
        status = ("FULL-LATTICE-MACHINERY-BUILT (gauge-invariant "
                  "pair+straight-string machinery measured per chain "
                  "length: chain-kernel dims ell=1 %d, ell=2 %d, "
                  "ell=3 singlet-sector %d (I_4 = %d, booked "
                  "prediction output); V371 contact reproduction ker "
                  "2/6 + d_k set + -64 I_2 leak %.1e; V368 free-gap "
                  "reproduction %.1e; string Casimir %.6f/link "
                  "measured; I_2=1, I_3=1 and the V100/V101 24-state"
                  " glued-kernel anchor reproduced from measured "
                  "fusion data)"
                  % (ker_c1, ker_c2, Nc3s, i4_meas, leak1,
                     abs(free_min - v368_min), c_F))
    else:
        status = ("V373-MACHINERY-DRIFT (a booked gate failed; no "
                  "spectral claim is made)")

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V373 done] %s" % status, flush=True)
    out = {"id": "V373_full_lattice_pair_string", "tier": 1,
           "residual": False,
           "property": "the full-lattice pair+straight-flux-string "
                       "machinery, gated on every reality-interface "
                       "anchor (the V371 full lattice-sum machine, "
                       "built at last)",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v373_full_lattice_pair_string.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    log("wall %.1f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
