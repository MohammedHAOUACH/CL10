#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V366 - singlet-extension census (reality-interface first cut,
booked BEFORE the run in pre_registrations_root.md). Kernel of
C2 = sum_a Q^a Q^a (diagonal certified Spin(10) gauge charge) on five
extended fibers built from the VERBATIM published generators; anchor
E0 reproduces V259 (published fiber ker = 0, min nonzero C2 =
45.000000). Group-theory anchors booked: ker(E1)=1, ker(E2)=1,
ker(E3)=0, ker(E4)=2, ker(E5)=9. Verdicts: SINGLET-ADMIT-* /
EXT-EMPTY / EXT-GROUP-MISMATCH / MACHINERY-DRIFT.
"""
import json
import os
import sys

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
from cl10 import altcell as AC  # noqa: E402

KTOL = 1e-9
PREDS = {"E1": 1, "E2": 1, "E3": 0, "E4": 2, "E5": 9}


def _ext_pair(TA, TB, Ls):
    """(A (x) B (x) link45) diagonal charges, 45 sparse matrices."""
    dA, dB = TA[0].shape[0], TB[0].shape[0]
    IA, IB = sp.eye(dA, format="csr"), sp.eye(dB, format="csr")
    I45 = sp.eye(45, format="csr")
    Q = []
    for a in range(45):
        M = (sp.kron(sp.kron(sp.csr_matrix(TA[a]), IB), I45)
             + sp.kron(sp.kron(IA, sp.csr_matrix(TB[a])), I45)
             + sp.kron(sp.kron(IA, IB), sp.csr_matrix(Ls[a])))
        Q.append(M.tocsr())
    return Q


def main():
    checks, data = [], {"fibers": {}}

    def log(m):
        print(m, flush=True)

    Ts = root_generators()                    # 45 on the 32
    f, _k, _a = structure_constants(Ts)
    Ls = adjoint_link_operators(f)            # 45 adjoint on 45
    _Bm, T16, _f16 = AC.weyl_generators()     # 45 on the 16
    T16c = [-t.T.copy() for t in T16]         # conjugate 16-bar: -T^T

    # E0 published fiber (anchor), V259's construction verbatim
    QP = [sp.kron(sp.eye(8, format="csr"), sp.csr_matrix(t)) for t in T16]
    QQ = [sp.kron(sp.csr_matrix(T), sp.eye(45, format="csr"))
          + sp.kron(sp.eye(32, format="csr"), sp.csr_matrix(L))
          for T, L in zip(Ts, Ls)]
    Q0 = []
    for a in range(45):
        Z = sp.lil_matrix((1568, 1568), dtype=np.complex128)
        Z[:128, :128] = QP[a].tolil()
        Z[128:, 128:] = QQ[a].tolil()
        Q0.append(Z.tocsr())
    C0 = sum(Q @ Q for Q in Q0)
    w0 = spla.eigsh(C0, k=4, which="SM", return_eigenvectors=False)
    w0.sort()
    ker0 = int(np.sum(np.abs(w0) < KTOL))
    e0_ok = (ker0 == 0 and abs(w0[0] - 45.0) < 1e-6)
    checks.append(("E0 anchor (published fiber, V259 verbatim): ker %d, "
                   "min C2 %.6f (booked 0 / 45.000000): %s"
                   % (ker0, w0[0], e0_ok), e0_ok))
    log("E0 anchor: ker %d, min C2 %.6f" % (ker0, w0[0]))

    I16a = sp.eye(16, format="csr")
    fibers = {
        "E1": [(sp.kron(sp.csr_matrix(T16[a]), I16a)
                + sp.kron(I16a, sp.csr_matrix(T16c[a]))).tocsr()
               for a in range(45)],                          # no link
        "E2": _ext_pair(T16, T16c, Ls),
        "E3": _ext_pair(T16, T16, Ls),
        "E4": _ext_pair(Ts, Ts, Ls),
    }
    I168, I16 = sp.eye(1568, format="csr"), sp.eye(16, format="csr")
    fibers["E5"] = [(sp.kron(Q0[a], I16)
                     + sp.kron(I168, sp.csr_matrix(T16c[a]))).tocsr()
                    for a in range(45)]

    ker_dims = {}
    for name in ("E1", "E2", "E3", "E4", "E5"):
        Q = fibers[name]
        dim = Q[0].shape[0]
        if dim <= 16384:
            C2 = sum(a @ a for a in Q)
            op = C2
        else:
            op = spla.LinearOperator(
                (dim, dim), matvec=lambda x, Q=Q: sum(a @ (a @ x) for a in Q),
                dtype=np.complex128)
        k = min(dim - 2, PREDS[name] + 8)
        w, v = spla.eigsh(op, k=k, which="SM", return_eigenvectors=True)
        W = np.asarray(v)
        idx = np.argsort(np.abs(w))
        w, W = w[idx], W[:, idx]
        ker = int(np.sum(np.abs(w) < KTOL))
        resid = 0.0
        if ker:
            Wk = W[:, :ker]
            G = Wk.conj().T @ Wk
            gv, GU = np.linalg.eigh(G)
            assert gv.min() > 1e-8, "%s kernel Gram rank-deficient" % name
            Wk = Wk @ (GU * gv ** -0.5) @ GU.conj().T
            C2W = sum(a @ (a @ Wk) for a in Q)
            resid = float(np.max(np.abs(C2W)))
            wnorm = float(np.max(np.abs(Wk.conj().T @ Wk - np.eye(ker))))
        else:
            wnorm = 0.0
        ker_dims[name] = ker
        data["fibers"][name] = {
            "dim": dim, "ker": ker,
            "min_eigs": [float(z) for z in w[:6]],
            "ker_resid": resid, "ker_nonunitarity": wnorm}
        log("%s dim %d: ker %d (pred %d), min eigs %s, resid %.1e, "
            "nonunit %.1e" % (name, dim, ker, PREDS[name],
                              ["%.4f" % z for z in w[:4]], resid, wnorm))
        ok = (resid < 1e-6 and wnorm < 1e-12) if ker else True
        checks.append(("%s kernel basis: dim %d, ||C2 W|| %.1e, W^dag W "
                       "defect %.1e: %s" % (name, ker, resid, wnorm, ok),
                       ok))

    gm = all(ker_dims[n] == PREDS[n] for n in PREDS)
    checks.append(("group-theory anchors ker = %s (booked %s): %s"
                   % (ker_dims, PREDS, gm), gm))

    ok = e0_ok and gm
    if gm and any(v > 0 for v in ker_dims.values()):
        status = ("SINGLET-ADMIT-" + "-".join("%s:%d" % (n, ker_dims[n])
                                              for n in PREDS)
                  + " (emptiness is a content theorem, not an action "
                  "theorem; V367 cascade live: dispersion + F-legs at "
                  "the singlet sectors)")
    elif gm and all(v == 0 for v in ker_dims.values()):
        status = "EXT-EMPTY"
    else:
        status = "EXT-GROUP-MISMATCH"
        ok = False
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V366 done] %s" % status, flush=True)
    out = {"id": "V366_singlet_extension_census", "tier": 2,
           "residual": True,
           "property": "does the certified Spin(10) gauge action admit "
                       "gauge singlets in ANY extended matter content",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v366_singlet_extension_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
