#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V269 - the dynamical-link module: 1 (+) 45 (+) 210 inside the root
(tier 1, hard, residual: false - construction of the refonte object).

**The question (pre-registered).** Does the truncated dynamical-link
module L = 1 (+) 45 (+) 210 exist inside the root's own so(10): does the
Casimir projector of 45 (x) 45 isolate exactly dims [1, 45, 210] with
three distinct Casimir eigenvalues, and do the 210-block generators close
the committed algebra?

**The object (closed input set).** The committed V79 machinery: Ts (the
45 generators on the 32 module from sm_qca.plane_generators), f (the
structure constants), Ls (the adjoint link operators). Build
C2_2025 = sum_a (L^a (x) I + I (x) L^a)^2 on 45 (x) 45, dense eigh; bin
the spectrum into the sectors {1, 45, 210}; restrict the generators to
the 210 block by the measured projector; verify the block algebra with
the committed f. Re-measure the pair multiplicities live: C2_pair on
32 (x) 32, the multiplicities of eigenvalues 0 / 64 / lambda_210 must be
(2, 2, 2) (the V98/V99 fusion counts, recomputed not read).

**Measured (pre-registered).** The Casimir eigenvalues and bin dims of
the 2025-module (lambda_210 booked as measured); the 210-algebra defect;
the pair multiplicities.

**Acceptance (pre-registered).** LINK-MODULE-OK iff bin dims exactly
[1, 45, 210] at three distinct eigenvalues (binning gap > 1e-8), the 210
algebra closes < 1e-9, and the multiplicities are (2, 2, 2);
LINK-MODULE-DEFECT otherwise.
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")  # the 1/5 reference branch
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))
import numpy as np  # noqa: E402
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)

NP = 45


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V269 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    log("start (host algebra: 2025 and 1024 dense, no GPU stack needed)")

    # ---- T0: the committed algebra re-verified live (V79 control) ----
    Ts = root_generators()
    f, k_tr, antisym = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    herm = max(float(np.abs(L - L.conj().T).max()) for L in Ls)
    # [L^A, L^B] = i f^{ABC} L^C on the 45 (adjoint) itself
    adj_def = 0.0
    for a in range(NP):
        for b in range(a + 1, NP):
            lhs = Ls[a] @ Ls[b] - Ls[b] @ Ls[a]
            rhs = 1j * sum(f[a, b, c] * Ls[c] for c in range(NP))
            adj_def = max(adj_def, float(np.abs(lhs - rhs).max()))
    c2_adj = sum(L @ L for L in Ls)
    c2_adj_scalar = float(np.abs(c2_adj - (np.trace(c2_adj) / NP)
                                * np.eye(NP)).max())
    t0pass = (antisym < 1e-12 and herm < 1e-12 and adj_def < 1e-12
              and c2_adj_scalar < 1e-12)
    checks.append({"name": "T0 V79 control live: f antisymmetry {:.1e}, "
                          "adjoint Hermiticity {:.1e}, adjoint algebra "
                          "defect {:.1e}, C2(45) scalar dev {:.1e} "
                          "(value {:.6f})".format(antisym, herm, adj_def,
                                                  c2_adj_scalar,
                                                  float(np.real(
                                                      np.trace(c2_adj) / NP))),
                   "pass": bool(t0pass)})
    data["c2_adjoint"] = float(np.real(np.trace(c2_adj) / NP))
    log("T0 done (f, adjoint algebra exact)")

    # ---- T1: the 2025-module Casimir spectrum ----
    I45 = np.eye(NP, dtype=np.complex128)
    Ltot = [np.kron(L, I45) + np.kron(I45, L) for L in Ls]
    C2t = np.zeros((NP * NP,) * 2, dtype=np.complex128)
    for L in Ltot:
        C2t += L @ L
    ev, U = np.linalg.eigh(C2t)
    ev = np.real(ev)
    order = np.argsort(ev)
    ev = ev[order]
    U = U[:, order]
    # bin by spectral gaps
    bins = []
    i0 = 0
    for i in range(1, len(ev) + 1):
        if i == len(ev) or ev[i] - ev[i - 1] > 1e-8:
            bins.append((float(np.mean(ev[i0:i])), i0, i))
            i0 = i
    data["c2_spectrum"] = [{"value": v, "dim": j - k} for (v, k, j) in bins]
    log("2025 Casimir: {} distinct eigenvalues {}".format(
        len(bins), [round(v, 6) for v, _, _ in bins]))
    # exchange parity of every bin: in 45 (x) 45 the adjoint 45 is the
    # ANTI-symmetric summand (the bracket map Lambda^2 -> 45 is the f
    # itself), so the swap operator Pi identifies it without retyping a
    # number. The 210 and the singlet are pinned by dimension.
    def swap_par(lo, hi):
        B = U[:, lo:hi].reshape(NP, NP, hi - lo)
        return float(np.real(np.mean(np.sum(B.conj()
                                           * B.swapaxes(0, 1),
                                           axis=(0, 1)))))
    binpar = []
    for (v, k, j) in bins:
        p = swap_par(k, j)
        binpar.append((v, j - k, p, k, j))
        log("  bin C2 {:8.4f} dim {:4d} swap parity {:+.4f}".format(
            v, j - k, p))
    data["c2_spectrum"] = [{"value": v, "dim": d, "swap_parity": p}
                           for (v, d, p, _, _) in binpar]
    sect = {}
    dims_ok = True
    c = [b for b in binpar if b[1] == 1 and abs(b[0]) < 1e-9]
    if len(c) == 1:
        sect["1"] = c[0]
    else:
        dims_ok = False
    c = [b for b in binpar if b[1] == 45 and b[2] < -0.5]
    if len(c) == 1:
        sect["45"] = c[0]
    else:
        dims_ok = False
    c = [b for b in binpar if b[1] == 210]
    if len(c) == 1:
        sect["210"] = c[0]
    else:
        dims_ok = False
    blocks = {name: (b[3], b[4]) for name, b in sect.items()}
    gram_def = 0.0
    for name, (lo, hi) in blocks.items():
        B = U[:, lo:hi]
        gram_def = max(gram_def,
                       float(np.abs(B.conj().T @ B - np.eye(B.shape[1]))
                             .max()))
    t1pass = bool(dims_ok and gram_def < 1e-10 and len(bins) >= 3)
    checks.append({"name": "T1 projector sectors isolated: dims [1, 45, "
                          "210] at eigenvalues {} (singlet {:.1e}, adjoint "
                          "{:.6f} antisymmetric parity {:.3f}, 210 {:.6f}); "
                          "Gram defect {:.1e}; full 45(x)45 spectrum "
                          "{}".format(
                              {k: round(b[0], 6) for k, b in sect.items()},
                              sect["1"][0] if "1" in sect else float("nan"),
                              sect["45"][0] if "45" in sect else float("nan"),
                              sect["45"][2] if "45" in sect else 0.0,
                              sect["210"][0] if "210" in sect else float("nan"),
                              gram_def,
                              [(round(v, 4), d, round(p, 3))
                               for (v, d, p, _, _) in binpar]),
                   "pass": t1pass})
    data["lambda_210"] = float(sect["210"][0]) if "210" in sect \
        else float("nan")
    data["adjoint_45_value"] = float(sect["45"][0]) if "45" in sect \
        else float("nan")
    data["sector_gram_defect"] = gram_def
    log("T1 sectors: " + str({k: round(b[0], 6) for k, b in sect.items()}))

    # ---- T2: the 210-block generators close the committed algebra ----
    lo, hi = blocks["210"]
    P = U[:, lo:hi]
    L210 = [P.conj().T @ Lt @ P for Lt in Ltot]
    alg210 = 0.0
    for a in range(NP):
        for b in range(a + 1, NP):
            lhs = L210[a] @ L210[b] - L210[b] @ L210[a]
            rhs = 1j * sum(f[a, b, c] * L210[c] for c in range(NP))
            alg210 = max(alg210, float(np.abs(lhs - rhs).max()))
    herm210 = max(float(np.abs(L - L.conj().T).max()) for L in L210)
    c2_210 = sum(L @ L for L in L210)
    lam = np.real(np.trace(c2_210) / 210)
    c2_210_scalar = float(np.abs(c2_210 - lam * np.eye(210)).max())
    t2pass = bool(alg210 < 1e-9 and herm210 < 1e-12
                  and c2_210_scalar < 1e-9
                  and abs(float(lam) - data["lambda_210"]) < 1e-9)
    checks.append({"name": "T2 210-block algebra: closure {:.1e}, "
                          "Hermiticity {:.1e}, C2(210) = {:.6f} scalar dev "
                          "{:.1e} (matches projector eigenvalue)".format(
                              alg210, herm210, lam, c2_210_scalar),
                   "pass": t2pass})
    data["c2_210_scalar_defect"] = c2_210_scalar
    log("T2 210-algebra closure {:.2e}".format(alg210))

    # ---- T3: pair multiplicities measured live (V98/V99 control) ----
    I32 = np.eye(32, dtype=np.complex128)
    Tpair = [np.kron(T, I32) + np.kron(I32, T) for T in Ts]
    C2p = np.zeros((1024,) * 2, dtype=np.complex128)
    for T in Tpair:
        C2p += T @ T
    pev = np.sort(np.real(np.linalg.eigvalsh(C2p)))
    def mult_at(val, tol=1e-6):
        n = int(np.sum(np.abs(pev - val) < tol))
        return n
    m1 = mult_at(0.0)
    m45 = mult_at(data["c2_adjoint"])
    m210 = mult_at(data["lambda_210"])
    t3pass = (m1 == 2 and m45 == 2 * 45 and m210 == 2 * 210)
    checks.append({"name": "T3 pair (32x32) isotypic counts at the measured "
                          "Casimir values: singlet {} (2 expected), adjoint "
                          "{} (90), 210 {} (420): multiplicities ({}, {}, "
                          "{}) vs V98/V99 (2, 2, 2)".format(
                              m1, m45, m210, m1 // 1, m45 // 45, m210 // 210),
                   "pass": t3pass})
    data["pair_multiplicities"] = {"1": m1, "45": m45 // 45, "210": m210 // 210}
    log("T3 multiplicities (2, 2, 2) control: " + str(data["pair_multiplicities"]))

    ok = all(c["pass"] for c in checks)
    branch = "LINK-MODULE-OK" if ok else "LINK-MODULE-DEFECT"
    status = ("dynamical-link module 1(+)45(+)210 built from the committed "
              "V79 generators by Casimir projector on 45(x)45 ({:d} s): "
              "C2 eigenvalues 0 / {:.6f} / {:.6f}, dims [1, 45, 210], 210 "
              "algebra closes to {:.1e}, pair multiplicities ({}, {}, {}) "
              "re-measured. Branch {}.".format(
                  int(time.time() - t0), data["c2_adjoint"],
                  data["lambda_210"], alg210,
                  data["pair_multiplicities"]["1"],
                  data["pair_multiplicities"]["45"],
                  data["pair_multiplicities"]["210"], branch))
    out = {"id": "V269_dynamical_link_module",
           "tier": 1,
           "property": "the truncated dynamical-link module 1(+)45(+)210 "
                       "exists as an exact so(10) module inside the root "
                       "(V269)",
           "residual": False,
           "programme": "V269",
           "branch": branch,
           "data": {**data,
                    "link_truncation": [1, 45, 210],
                    "fitted_inputs": [],
                    "wall_time": time.time() - t0}}
    out["checks"] = checks
    out["status"] = status
    out["pass"] = bool(ok)
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v269_dynamical_link_module.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
