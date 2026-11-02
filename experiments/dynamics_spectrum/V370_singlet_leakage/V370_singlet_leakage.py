#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V370 - singlet leakage under the certified bond-current operator
(booked BEFORE the run; interaction programme second cut). T0: the
corrected V369 decomposition as machine gate:
1/2(C2tot - C2_1 - C2_2 - C2_L) = J_1.J_2 + J.L. Gauss on the kernel
gives J.L W = -(90 + 2 J_1.J_2) W, so leakage is controlled by the
current-current contraction. Uses the V369 ARPACK cache (= kernel +
four lowest channels, exactly 400 columns). Verdicts: LEAKAGE-ZERO /
LEAKS-<l>-DOM-<c> / GATE-FAIL. Tier 2, zero fitted.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "dynamics_spectrum",
                                "V369_pair_binding_onset"))
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)

dim = 32 * 32 * 45
CACHE = os.path.join(ROOT, "experiments", "dynamics_spectrum",
                     "V369_pair_binding_onset", "results", "v369_sm400.npz")


def main():
    checks, data = [], {}

    def log(m):
        print(m, flush=True)

    Ts = root_generators()
    f, _k, _a = structure_constants(Ts)
    Ls = adjoint_link_operators(f)

    def T1a(a, t):
        return np.einsum("ab,bgm->agm", Ts[a], t)

    def T2a(a, t):
        return np.einsum("ab,gbm->gam", Ts[a], t)

    def La(a, t):
        return np.einsum("ab,gmb->gma", Ls[a], t)

    def C2tot(x):
        h = x.reshape(32, 32, 45)
        acc = np.zeros_like(h)
        for a in range(45):
            z = T1a(a, h) + T2a(a, h) + La(a, h)
            acc += T1a(a, z) + T2a(a, z) + La(a, z)
        return acc.ravel()

    def JL(x):                        # sum_a (T1_a + T2_a) L_a
        h = x.reshape(32, 32, 45)
        acc = np.zeros_like(h)
        for a in range(45):
            acc += T1a(a, La(a, h)) + T2a(a, La(a, h))
        return acc.ravel()

    # ---------------- T0 corrected decomposition (the V369 error, made
    # a checkable identity): 1/2(C2tot - C21 - C22 - C2L) = J1.J2 + J.L
    rng = np.random.default_rng(0x99)
    t0 = 0.0
    for _ in range(6):
        v = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        h = v.reshape(32, 32, 45)
        c2 = C2tot(v).reshape(32, 32, 45)
        b1 = np.zeros_like(h)
        b2 = np.zeros_like(h)
        b3 = np.zeros_like(h)
        j12 = np.zeros_like(h)
        for a in range(45):
            b1 += T1a(a, T1a(a, h))
            b2 += T2a(a, T2a(a, h))
            b3 += La(a, La(a, h))
            j12 += T1a(a, T2a(a, h))
        jl = JL(v).reshape(32, 32, 45)
        lhs = 0.5 * (c2 - b1 - b2 - b3)
        rhs = j12 + jl
        t0 = max(t0, float(np.max(np.abs(lhs - rhs))
                           / max(np.abs(lhs).max(), 1e-12)))
    checks.append(("T0 corrected identity 1/2(C2tot-blocks) = "
                   "1/2 J1.J2 sym + J.L (rel <= 1e-9): %.1e: %s"
                   % (t0, t0 < 1e-9), t0 < 1e-9))

    if os.path.exists(CACHE):
        z = np.load(CACHE)
        w, V = z["w"], z["V"]
    else:
        # fallback on a fresh checkout (V369 normally builds this
        # cache first; the gate runs experiments in sorted order):
        # same k=400 SM run, cupy matvec, V369 verbatim machinery.
        import cupy as cp
        import scipy.sparse.linalg as spla
        Tgc = [cp.asarray(t) for t in Ts]
        Lgc = [cp.asarray(L) for L in Ls]

        def c2_matvec(x):
            h = cp.asarray(x).reshape(32, 32, 45)
            y = cp.zeros_like(h)
            for a in range(45):
                z2 = (cp.einsum("ab,bgm->agm", Tgc[a], h)
                      + cp.einsum("ab,gbm->gam", Tgc[a], h)
                      + cp.einsum("ab,gmb->gma", Lgc[a], h))
                y += (cp.einsum("ab,bgm->agm", Tgc[a], z2)
                      + cp.einsum("ab,gbm->gam", Tgc[a], z2)
                      + cp.einsum("ab,gmb->gma", Lgc[a], z2))
            return cp.asnumpy(y.ravel())
        op = spla.LinearOperator((dim, dim), matvec=c2_matvec,
                                 dtype=np.complex128)
        w0, v0 = spla.eigsh(op, k=400, which="SM")
        order0 = np.argsort(w0)
        w, V = w0[order0], np.asarray(v0)[:, order0]
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        np.savez(CACHE, w=w, V=V)
    op = None
    ker = np.abs(w) < 1e-9
    W = V[:, ker]
    G = W.conj().T @ W
    gv, GU = np.linalg.eigh(G)
    W = W @ (GU * gv ** -0.5) @ GU.conj().T
    nk = W.shape[1]
    c2res = float(np.max(np.abs(np.column_stack([C2tot(W[:, d])
                                                 for d in range(nk)]))))
    gate_k = nk == 2 and c2res < 1e-9
    checks.append(("kernel gate: dim %d, ||C2 W|| %.1e: %s"
                   % (nk, c2res, gate_k), gate_k))

    # ------------------------------------ M1 intra-kernel block of J.L
    JLW = np.column_stack([JL(W[:, d]) for d in range(nk)])
    B = W.conj().T @ JLW
    B = 0.5 * (B + B.conj().T)
    be = np.linalg.eigvalsh(B)
    log("M1 B = W^dag (J.L) W eigenvalues: %.6f %.6f" % (be[0], be[1]))

    # --------------------------------------- M2 channel leakage
    ltot = float(np.linalg.norm(JLW - W @ B))
    nz = np.abs(w) >= 1e-9
    vals = w[nz]
    order = np.argsort(vals)
    groups, cur = [], [order[0]]
    for o in order[1:]:
        if vals[o] - vals[cur[-1]] < 1e-5:
            cur.append(o)
        else:
            groups.append(cur)
            cur = [o]
    groups.append(cur)
    # booked hygiene correction (V366 precedent, cross-cluster class):
    # the projector-free total leakage bounds EVERY channel weight by
    # ||eps||, so per-channel weights are meaningful only in the
    # LEAKING branch (ltot > 1e-12), where they would need a properly
    # cross-orthogonalized channel basis (the SM cache carries O(1)
    # CROSS-cluster overlap at k=400 - first pass measured l_36 =
    # 90.51 = 64*sqrt(2) = ||J.L W|| against a total of 7.3e-14: pure
    # Lanczos contamination).
    lchan, residue = {}, 0.0
    if ltot > 1e-12:
        lchan, dom, lsum2 = {}, None, 0.0
        prev = [W]
        for g in groups:
            cval = float(np.mean(vals[g]))
            P = V[:, np.array(g)]
            for Q in prev:
                P = P - Q @ (Q.conj().T @ P)
            Gc = P.conj().T @ P
            gv2, GU2 = np.linalg.eigh(Gc)
            assert gv2.min() > 1e-8, "cluster Gram rank-deficient: " \
                "cached SM basis insufficient for channel resolution"
            P = P @ (GU2 * gv2 ** -0.5) @ GU2.conj().T
            lc = float(np.linalg.norm(P.conj().T @ JLW))
            lchan[round(cval, 4)] = lc
            lsum2 += lc * lc
            if dom is None or lc > lchan[dom]:
                dom = round(cval, 4)
            prev.append(P)
        residue = max(0.0, ltot ** 2 - lsum2) ** 0.5
    log("M2 leakage total l = %.6e; channel weights %s (moot when the "
        "total is machine-zero); residue %.1e"
        % (ltot, lchan if lchan else "not-resolved (zero branch)",
           residue))

    data["B_eigs"] = [float(be[0]), float(be[1])]
    data["l_total"] = ltot
    data["l_channels"] = lchan
    data["residue"] = residue
    gates_ok = t0 < 1e-9 and gate_k
    if not gates_ok:
        status = "V370-GATE-FAIL"
        ok = False
    elif ltot <= 1e-12:
        status = ("LEAKAGE-ZERO (the certified bond-current interaction "
                  "preserves gauge neutrality of the singlet pair at "
                  "first order; intra-kernel shifts %.6f / %.6f; binding "
                  "still needs the channel spectral comparison - V371 "
                  "question)" % (be[0], be[1]))
        ok = True
    else:
        dc = max(lchan, key=lchan.get)
        status = ("LEAKS-%.3e-DOM-%s (the published bond-current "
                  "interaction mixes the neutral pair into charged "
                  "channels at first order, dominant channel %s weight "
                  "%.3e, residue outside the four lowest channels "
                  "%.1e; the naive interaction route to a physical "
                  "sector closes at this cut)"
                  % (ltot, dc, dc, lchan[dc], residue))
        ok = True
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V370 done] %s" % status, flush=True)
    out = {"id": "V370_singlet_leakage", "tier": 2, "residual": True,
           "property": "does the certified bond-current interaction "
                       "preserve or mix the neutral pair subspace",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "v370_singlet_leakage.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
