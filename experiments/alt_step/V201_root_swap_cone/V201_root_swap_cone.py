#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V201 - the root's own sector swap in the six-stage product: does a cone appear (O9)

V192 searched the largest six-stage stage family that stays exactly
O_h-covariant at Gamma (angle a on the three intra-block stages, b on the three
inter-block ones) jointly with the internal angle, on a 26 x 26 x 50 grid, and
found no conical node: the exactly gapless points split as the parity theorem
forces - single-parity spaces flat, parity-mixed spaces Abelian (tangent
commutators rho at most 1e-6 at the polished points).

The axiom's internal layer is a 16-D Weyl rotation U(theta)^6. The root's
Clifford algebra carries exactly one further kind of involution, built from no
new generator: with w = h_1 h_2 ... h_10 the volume element, each S_c = h_c w
(c = 1..10) is Hermitian, squares to I, and anticommutes with w, so on the 32-D
Dirac module it maps the 16 to the 16-bar - a sector swap, an element of the
root. This programme interleaves that swap with the internal layer of the
six-stage product: the internal factor becomes U_32(theta)^6 (h_c w) on the
root's own Dirac module (16 + 16-bar - the same module the cell carries), and
the V192 search reruns on the identical grid for each of the 10 swaps. The
32-D construction is anchored to V192 by the exact regression
B^dag U_32 B = U_16 (the Weyl intertwiner of altcell), so the unswapped part is
the V192 step verbatim.

Pre-registered branch: CONE if a polished four-cancellation point carries a
gapless on-site space of dimension >= 4 with mixed-parity O_h content,
non-commuting tangent operators (rho > 0.5) and anisotropy < 1e-6 (the V192
cone criterion); NO-CONE otherwise, with the best four-state residual per swap
recorded. V192's causal front-speed check is built on the 16-D internal and is
not applicable to the doubled module; a CONE verdict names it as the open
follow-up.

Reuse note (the V126-from-V85 pattern): P6, G48/W48, the O_h character table,
GAM and the core come from V192 by import; oh_content and cone_at are adapted
here for the 32-D internal factor (the only two functions touched).
"""
import json
import os
import sys
import time
import os

os.environ["CL10_BACKEND"] = "numpy"
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np                                              # noqa: E402
from scipy.optimize import minimize                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                      # noqa: E402
sys.path.insert(0, exp("V192_engine_two_angle_cone"))

from V192_engine_two_angle_cone import (C, CHAR, CLS, GAM, P6,  # noqa: E402
                                        W48)
from cl10 import altcell as A                                   # noqa: E402

NI2 = 32                                   # the root's Dirac module: 16+16bar
DIRS = C.dirs

H32 = A.signature_generators(0)            # the 10 root generators (32-D)
W32 = A.pseudoscalar_sig(0)                # w = h_1 ... h_10
T32 = A.plane_generators_sig(0)            # the 45 planes, lexicographic
SWAP = [H32[c] @ W32 for c in range(10)]   # the 10 sector swaps
I32 = np.eye(NI2, dtype=complex)

W48_ON = [np.kron(W, I32) for W in W48]


def U32(th):
    """the 45-plane product on the Dirac module, the Core.U order; a scalar
    angle means the uniform angle vector (the V192 convention)."""
    if np.isscalar(th):
        th = np.full(45, th)
    U = I32.copy()
    for k in range(45):
        U = U @ (np.cos(th[k]) * I32 - 1j * np.sin(th[k]) * T32[k])
    return U


def internal_factor(th, c):
    """the modified internal layer: U_32(theta)^6 (h_c w)."""
    return np.linalg.matrix_power(U32(th), 6) @ SWAP[c]


def step_swap(a, b, th, c, q=GAM):
    return np.kron(P6(q, a, b), internal_factor(th, c))


def wrap(x):
    return np.angle(np.exp(1j * x))


def oh_content32(Q):
    """O_h content of a subspace of the doubled module (bond action only)."""
    chi = [np.trace(Q.conj().T @ W @ Q).real for W in W48_ON]
    out = {}
    for n, ch in CHAR.items():
        m = sum(ch[CLS[k]] * chi[k] for k in range(48)) / 48.0
        if m > 0.4:
            out[n] = int(round(m))
    return out


def cone_at32(a, b, th, c):
    """the on-site node of the swap-modified product at Gamma (adapted from
    V192's cone_at: the internal factor is 32-D and q-free, so the tangent
    derivatives act on the bond factor alone)."""
    M0 = step_swap(a, b, th, c)
    wv, V = np.linalg.eig(M0)
    ph = np.angle(wv)
    gap = float(np.abs(wrap(ph)).min())
    Q = V[:, np.abs(wrap(ph)) < 1e-12]
    Q, R = np.linalg.qr(Q)
    ng = int(np.sum(np.abs(np.diag(R)) > 1e-8))
    Q = Q[:, :ng]
    if ng == 0:
        return dict(gap=gap, gapless_onsite=0, onsite_dim=int(ph.size),
                    A=1.0, rho=0.0, v=0.0, chirality=0.0, content={})
    DQ = 1e-5
    Fm = internal_factor(th, c)
    S = []
    for k in range(3):
        e = np.zeros(3)
        e[k] = DQ
        dM = (np.kron(P6(e, a, b), Fm) - np.kron(P6(-e, a, b), Fm)) / (2 * DQ)
        S.append(1j * Q.conj().T @ M0.conj().T @ dM @ Q)
    nrm = [float(np.linalg.norm(s)) for s in S]
    com = max(float(np.linalg.norm(S[i] @ S[j] - S[j] @ S[i]))
              for i in range(3) for j in range(i + 1, 3))
    rho = com / (sum(n * n for n in nrm) / 3.0) if max(nrm) > 0 else 0.0
    sp = np.array([np.sort(np.abs(np.linalg.eigvalsh(
        d[0] * S[0] + d[1] * S[1] + d[2] * S[2]))) for d in DIRS])
    vmax = sp.max(axis=1)
    An = float((vmax.max() - vmax.min()) / vmax.max()) if vmax.max() > 0 else 1.0
    chi = (float(np.imag(np.trace(S[0] @ S[1] @ S[2])) /
                 (nrm[0] * nrm[1] * nrm[2])) if min(nrm) > 0 else 0.0)
    return dict(gap=gap, gapless_onsite=ng, onsite_dim=int(ph.size),
                A=An, rho=rho, v=float(vmax.max()), chirality=chi,
                content=oh_content32(Q))


def main():
    t0 = time.time()
    checks = []
    out = {"experiment": "V201_root_swap_cone", "data": {}}

    # ---- T0 the ten involutions, and the anchor to the V192 step ----------
    inv_err = unit_err = antic_err = herm_err = 0.0
    for S in SWAP:
        inv_err = max(inv_err, float(np.abs(S @ S - I32).max()))
        unit_err = max(unit_err, float(np.abs(S @ S.conj().T - I32).max()))
        antic_err = max(antic_err, float(np.abs(S @ W32 + W32 @ S).max()))
        herm_err = max(herm_err, float(np.abs(S - S.conj().T).max()))
    B, _T16, _f = A.weyl_generators()
    anchor = 0.0
    for th in (0.7, 2.36):
        U1 = np.linalg.matrix_power(C.U(np.full(45, th)), 6)
        U2 = B.conj().T @ np.linalg.matrix_power(U32(np.full(45, th)), 6) @ B
        anchor = max(anchor, float(np.abs(U1 - U2).max()))
    checks.append(("T0 the root supplies exactly the sector swaps: each of the "
                   "10 S_c = h_c w (volume element w = h_1..h_10 of the root) "
                   "satisfies S^2 = I to %.1e, S S^dag = I to %.1e, "
                   "{S, w} = 0 to %.1e and S = S^dag to %.1e - Hermitian "
                   "involutions of the root algebra, no new generator and no "
                   "new angle; and the doubled construction anchors exactly to "
                   "the V192 internal step: B^dag U_32^6 B = U_16^6 to %.1e "
                   "(Weyl intertwiner, two angles)"
                   % (inv_err, unit_err, antic_err, herm_err, anchor),
                   max(inv_err, unit_err, antic_err, herm_err) < 1e-10
                   and anchor < 1e-10))

    # ---- T1 Gamma equivariance of the modified product --------------------
    equi_err = 0.0
    M0 = step_swap(0.2, 0.2, np.full(45, 0.2), 0)
    for W in W48_ON:
        equi_err = max(equi_err, float(np.abs(W @ M0 - M0 @ W).max()))
    checks.append(("T1 the swap interleaves transparently: the modified "
                   "six-stage product at Gamma commutes with all 48 cubic "
                   "elements to %.1e (the swap lives in the internal factor, "
                   "the cubic action in the bond factor) - the covariant "
                   "setting of V192 is untouched" % equi_err, equi_err < 1e-10))

    # ---- T2 the V192 census rerun with each swap ---------------------------
    ths = np.linspace(0.02, np.pi - 0.02, 50)
    iph = {}
    for ti, th in enumerate(ths):
        V6 = np.linalg.matrix_power(U32(np.full(45, th)), 6)
        for c in range(10):
            iph[(ti, c)] = np.angle(np.linalg.eigvals(V6 @ SWAP[c]))

    def bond_phases(a, b):
        return np.angle(np.linalg.eigvals(P6(GAM, a, b)))

    def four_smallest_swap(a, b, th, c, bp=None):
        bp = bond_phases(a, b) if bp is None else bp
        V6 = np.linalg.matrix_power(U32(np.full(45, th)), 6)
        vp = np.angle(np.linalg.eigvals(V6 @ SWAP[c]))
        s = wrap(bp[:, None] + vp[None, :]).ravel()
        return float(np.sort(np.abs(s))[:4].sum())

    grid_pts = []
    bpc = {}
    for a in np.linspace(0.05, np.pi - 0.05, 26):
        for b in np.linspace(0.05, np.pi - 0.05, 26):
            key = (a, b)
            if key not in bpc:
                bpc[key] = bond_phases(a, b)
            for ti in range(50):
                grid_pts.append((a, b, ti, bpc[key]))

    per_c = []
    all_rows = []
    for c in range(10):
        cand = []
        for (a, b, ti, bp) in grid_pts:
            s = wrap(bp[:, None] + iph[(ti, c)][None, :]).ravel()
            cand.append((float(np.sort(np.abs(s))[:4].sum()), a, b, ths[ti]))
        cand.sort()
        best_grid = cand[0][0]
        pol = []
        for f0, a, b, th in cand[:8]:
            r = minimize(lambda v: four_smallest_swap(v[0], v[1], v[2], c),
                         [a, b, th], method="Nelder-Mead",
                         options=dict(xatol=1e-13, fatol=1e-16, maxiter=2500))
            if (r.fun < 1e-10 and 0.01 < r.x[0] < np.pi - 0.01
                    and 0.01 < r.x[1] < np.pi - 0.01
                    and 0.01 < r.x[2] < np.pi - 0.01):
                pol.append((float(r.fun), float(r.x[0]), float(r.x[1]),
                            float(r.x[2])))
        uniq = []
        for f, a, b, th in sorted(pol):
            if not any(abs(a - u[1]) < 1e-4 and abs(b - u[2]) < 1e-4
                       and abs(th - u[3]) < 1e-4 for u in uniq):
                uniq.append((f, a, b, th))
        rows = []
        for f, a, b, th in uniq:
            d = cone_at32(a, b, th, c)
            d.update(a=a, b=b, theta=th, swap=c, residual=f)
            rows.append(d)
        good = [r for r in rows if r["gap"] < 1e-8 and r["A"] < 1e-6
                and r["rho"] > 0.5 and r["gapless_onsite"] >= 4
                and r["v"] > 1e-6
                and len({k[-1] for k in r["content"]}) > 1]
        per_c.append({"swap": c, "best_grid_residual": best_grid,
                      "exact_points": len(uniq), "conical": len(good),
                      "max_rho": max([r["rho"] for r in rows], default=0.0),
                      "max_slope": max([r["v"] for r in rows], default=0.0)})
        all_rows.extend(rows)

    cone_points = [r for r in all_rows
                   if r["gap"] < 1e-8 and r["A"] < 1e-6 and r["rho"] > 0.5
                   and r["gapless_onsite"] >= 4 and r["v"] > 1e-6
                   and len({k[-1] for k in r["content"]}) > 1]
    out["data"]["per_swap"] = per_c
    out["data"]["rows"] = all_rows
    checks.append(("T2 the V192 census rerun with each of the 10 root swaps: "
                   "%d grid points per swap (26 x 26 x 50, identical grid); "
                   "best four-state grid residual per swap %s; exact "
                   "four-cancellation points %s; largest tangent commutator "
                   "rho at those points %s - conical nodes (gap < 1e-8, "
                   "anisotropy < 1e-6, rho > 0.5, gapless >= 4, mixed "
                   "parity): %d"
                   % (len(grid_pts), [round(pc["best_grid_residual"], 3)
                                      for pc in per_c],
                      [pc["exact_points"] for pc in per_c],
                      [round(pc["max_rho"], 3) for pc in per_c],
                      len(cone_points)),
                   len(grid_pts) == 33800))

    # ---- T3 branch ----------------------------------------------------------
    branch = "CONE" if cone_points else "NO-CONE"
    out["branch"] = branch
    ok = all(c[1] for c in checks)

    out["pass"] = bool(ok)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    if branch == "CONE":
        r = cone_points[0]
        out["status"] = (
            "the cone is inside the root: with the internal layer U_32(theta)"
            "^6 (h_%d w), the six-stage product carries at Gamma a %d-state "
            "gapless on-site space with O_h content %s, tangent commutator "
            "rho = %.3f, anisotropy %.1e, slope %.2e at a = %.4f, b = %.4f, "
            "theta = %.5f (V201). The sector swap h_c w is a root involution "
            "(S^2 = I, {S, w} = 0): the internal layer of the axiom extends "
            "by an element of the root and reaches the cone; the causal "
            "front-speed check on the doubled module is the named follow-up"
            % (r["swap"] + 1, r["gapless_onsite"], r["content"], r["rho"],
               r["A"], r["v"], r["a"], r["b"], r["theta"]))
    else:
        out["status"] = (
            "the sector swap does not reach the cone: all 10 root involutions "
            "h_c w interleaved with the internal layer, rerun on the V192 "
            "grid (33800 points each), leave best four-state grid residuals "
            "%s and exact four-cancellation points with tangent commutator "
            "rho at most %s: conical nodes 0. The internal module doubled to "
            "the root's own Dirac module (16 + 16-bar, anchored to the V192 "
            "Weyl step at %.0e) is measured insufficient in the two-angle "
            "stage family - the parity-mixed pairing that V192's no-go "
            "requires does not turn non-Abelian through the swap: O9 stays "
            "open and the step-class extension is a measured negative "
            "(V201)"
            % ([round(pc["best_grid_residual"], 3) for pc in per_c],
               [round(pc["max_rho"], 3) for pc in per_c], anchor))
    out["tier"] = 3
    out["property"] = (
        "the root's 10 sector-swap involutions h_c w with the two-angle "
        "six-stage family carry no conical node at Gamma (V201)"
        if branch == "NO-CONE" else
        "the root's sector swap extends the internal layer to a conical node "
        "at Gamma (V201)")
    out["programme"] = "V201"
    out["wall_time"] = round(time.time() - t0, 3)

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v201_root_swap_cone.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V201 - the root sector swap and the cone (branch %s)" % branch)
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    for pc in per_c:
        print("  swap c=%2d  best grid %.3e  exact pts %3d  conical %d  "
              "max rho %.3f  max slope %.2e"
              % (pc["swap"] + 1, pc["best_grid_residual"], pc["exact_points"],
                 pc["conical"], pc["max_rho"], pc["max_slope"]))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
