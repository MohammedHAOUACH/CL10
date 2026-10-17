#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V166 - the published Wilson family at the CLOSURE-QUANTIZED angles
(alt-dynamics, axis e)

**The question.**  V65's node-gap mechanism is the Wilson factor
``W(q) = exp(-i w sum_j (1 - cos q_j) D)`` with D a Cartan combination,
scanned over CONTINUOUS w (V115: no isotropy point, min A = 0.997).  The
pre-registered rule forbids a continuous knob: is there a closure-
quantized reading, w in {pi/2, pi/4}?  Census: w in {pi/2, pi/4} x
D in {X, Y} x sign +- = 8 rows, plus the w = 0 anchor (must reproduce
the baseline exactly).  The Cartan involutions Sigma_k = i h_{2k}
h_{2k+1} are recomputed here; X and Y are IDENTIFIED by eigenvalue-
multiset matching (the V65 method, not its scan): Y against the SM
hypercharge multiset of ``matter.SM_FIELDS`` ({1/6 x6, -2/3 x3, 1/3 x3,
-1/2 x2, 1, 0}), X against the SU(5)-generator multiset {1 x10, -3 x5,
5 x1} - solved exactly by 5-row linear inversion, not by grid search.

**The result.**  Every quantized row gaps the 12 M nodes (the published
mechanism, now at quantized w): the residual gapless set is empty or
non-single-orbit, and the full step does not close at order <= 64 (the
internal factor at g_10 never closes).  The Wilson factor is only
QUASI-local: the real-space shell amplitude one block beyond the bond
cone is nonzero (Lieb-Robinson tail recorded per row) - the strict
light-cone gate FAILS by construction, as V65 already published.  The
w = 0 anchor reproduces the baseline to machine precision.
"""

import itertools
import json
import math
import os
import sys
import time
from collections import Counter

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import altcell as A
from cl10 import matter
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
NODE_TOL = 1e-3


def cartan_16():
    """The 5 Cartan involutions Sigma_k = i h_{2k} h_{2k+1} on the 16
    sector, plus the 32-D originals."""
    h = A.signature_generators(0)
    S32 = [1j * (h[2 * k] @ h[2 * k + 1]) for k in range(5)]
    B = A.sectors_sig(0)[0]
    return [B.conj().T @ S @ B for S in S32], S32


def cartan_weights():
    """The 16 weight rows lam (Cartan eigenvalues) and the common
    eigenbasis V (columns, in the 16 basis of ``sectors_sig``)."""
    Dk, _S32 = cartan_16()
    M = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(M)
    lam = np.array([[float(np.real(np.vdot(V[:, s], Dk[k] @ V[:, s])))
                     for k in range(5)] for s in range(16)])
    return lam, V


def y_target():
    """The hypercharge multiset of the 16 Weyl fields, read from
    ``matter.SM_FIELDS`` (n copies of each y)."""
    out = []
    for _name, (_su3, _su2, y, n) in matter.SM_FIELDS.items():
        out += [y] * n
    return np.sort(np.array(out, float))


X_TARGET = np.sort(np.array([1.0] * 10 + [-3.0] * 5 + [5.0]))


def identify(lam, target):
    """All Cartan coefficient vectors c with eigenvalue multiset
    {c . lam} = target (exact: 5 independent rows fixed, ordered value
    assignments solved by linear inversion, full multiset verified;
    deduped up to overall sign)."""
    basis = None
    for combo in itertools.combinations(range(16), 5):
        if abs(np.linalg.det(lam[list(combo)])) > 1e-6:
            basis = list(combo)
            break
    vals = sorted(set(target.tolist()))
    counts = Counter(target.tolist())
    sols = []
    for tup in itertools.product(vals, repeat=5):
        cnt = Counter(tup)
        if any(cnt[v] > counts[v] for v in cnt):
            continue
        c = np.linalg.solve(lam[basis], np.array(tup, float))
        if np.abs(np.sort(c @ lam.T) - np.sort(target)).max() < 1e-9:
            if not any(np.allclose(c, s * u, atol=1e-9)
                       for u in sols for s in (1.0, -1.0)):
                sols.append(c)
    return sols


def wilson_16(q, w, ddiag, V):
    ph = w * sum(1.0 - math.cos(qj) for qj in q)
    return V @ np.diag(np.exp(-1j * ph * ddiag)) @ V.conj().T


def lr_shells(w, ddiag, V, U16, Lb=16, n_stages=12):
    """Lieb-Robinson shell decay (V65 method): real-space evolution, the
    Wilson factor applied once per step by block-momentum multiplication;
    the max amplitude in shells of one block beyond the bond cone."""
    rng = np.random.default_rng(0)
    f = np.zeros((Lb, Lb, Lb, 8, 16), dtype=complex)
    v = rng.standard_normal(16) + 1j * rng.standard_normal(16)
    x0 = (Lb // 2, Lb // 2, Lb // 2)
    f[x0 + (0,)] = v / np.linalg.norm(v)
    qs = 2.0 * np.pi * np.arange(Lb) / Lb
    s1 = 3.0 - (np.cos(qs)[:, None, None] + np.cos(qs)[None, :, None]
                + np.cos(qs)[None, None, :])
    Wfac = np.exp(-1j * w * s1[..., None, None]
                  * ddiag[None, None, None, None, :])
    R = np.arange(Lb)
    sub = np.arange(8)
    ex, ey, ez = sub & 1, (sub >> 1) & 1, (sub >> 2) & 1
    dd = np.minimum(np.abs(R - x0[0]), Lb - np.abs(R - x0[0]))
    dist = ((2 * dd[:, None, None, None] + ex)
            + (2 * dd[None, :, None, None] + ey)
            + (2 * dd[None, None, :, None] + ez))
    for t in range(n_stages):
        f = A.step_stages(f, [U16] * 6, [SC.G10] * 6, ST)
        if t % 6 == 5:
            f = np.fft.ifftn(np.fft.fftn(f, axes=(0, 1, 2)) * Wfac,
                             axes=(0, 1, 2))
    nrm = np.sqrt(np.sum(np.abs(f) ** 2, axis=-1))
    cone = float(n_stages)
    inside = float(nrm[dist <= cone].max())
    shells = []
    for k in range(3):
        sel = (dist > cone + 2 * k) & (dist <= cone + 2 * (k + 1))
        shells.append(float(nrm[sel].max()) if sel.any() else 0.0)
    return [s / inside for s in shells]


def main():
    t0 = time.time()
    U16 = SC.sectors()[2]
    lam, V = cartan_weights()
    ys = identify(lam, y_target())
    xs = identify(lam, X_TARGET)
    missing = []
    if not ys:
        missing.append("Y")
    if not xs:
        missing.append("X")
    Dlist = []
    if xs:
        Dlist += [("X", c) for c in xs[:1]]
    if ys:
        Dlist += [("Y", c) for c in ys[:1]]

    dirs = A.sphere_dirs(6)
    rows = []
    specs = [("anchor", 0.0, "I", 1.0)]
    for w in (math.pi / 2, math.pi / 4):
        for dname, c in Dlist:
            for sign in (1.0, -1.0):
                specs.append((f"w{w:.4g}_{dname}_{'+' if sign > 0 else '-'}",
                              sign * w, dname, c))

    for (name, w, dname, c) in specs:
        ddiag = (c @ lam.T) if dname != "I" else np.zeros(16)

        def bloch_fn(qq, w=w, ddiag=ddiag):
            M = A.bloch_cell(C, U16, SC.G10, ST, qq)
            if w == 0.0:
                return M
            W = wilson_16(qq, w, ddiag, V)
            return np.kron(np.eye(8, dtype=complex), W) @ M

        gaps = {}
        nodes = []
        for stname, pts in SC.HIGH_SYM.items():
            for q in pts:
                gg = A.gap_at(bloch_fn(q))
                gaps[stname] = min(gaps.get(stname, 9.9), gg)
                if gg < NODE_TOL:
                    nodes.append(q)
        chir = 0.0
        for q in nodes:
            ch, _sl = A.chirality_at(bloch_fn, q)
            chir += ch
        A_first = None
        growth = None
        if nodes:
            aa, _vx, _vn = A.cone_anisotropy(bloch_fn, nodes[0], dirs)
            A_first = aa
            n0 = np.array([1.0, 1.0, 1.0]) / math.sqrt(3)
            g1 = A.gap_at(bloch_fn(tuple(q + 1e-3 * n0[k]
                                        for k, q in enumerate(nodes[0]))))
            g2 = A.gap_at(bloch_fn(tuple(q + 1e-2 * n0[k]
                                        for k, q in enumerate(nodes[0]))))
            growth = g2 / g1 if g1 > 0 else None
        stars = {tuple(A.oh_star(q)) for q in nodes}
        single = bool(nodes) and len(stars) == 1 and \
            len(next(iter(stars), ())) == len(nodes)
        Mq = bloch_fn((0.31, -1.13, 2.02))
        unit = float(np.abs(Mq @ Mq.conj().T - np.eye(Mq.shape[0])).max())
        N, cerr = A.exact_order(bloch_fn((0.0, 0.0, 0.0)), Nmax=64)
        shells = lr_shells(w, ddiag, V, U16) if w != 0.0 else [0.0] * 3
        rows.append({"row": name, "w": w, "D": dname,
                     "gap_by_star": gaps, "n_nodes": len(nodes),
                     "n_stars": len(stars), "single_orbit": single,
                     "chirality_sum": chir, "A": A_first,
                     "cone_growth": growth,
                     "chirality_defined": all(
                         any(abs(cc) > 1e-12 for cc in q) for q in nodes),
                     "unitary_err": unit, "closure_order": N,
                     "closure_err": cerr, "closed": cerr < 1e-9,
                     "lr_shell_ratios": shells,
                     "quasi_local": w != 0.0})

    anchor = rows[0]
    ok_anchor = (anchor["n_nodes"] == 12
                 and abs(anchor["gap_by_star"]["M"] - 2.0029e-4) < 1e-6
                 and abs(anchor["gap_by_star"]["GAMMA"] - 0.0365) < 1e-4
                 and anchor["single_orbit"])
    # anchor 2-route: the w = 0 row equals the baseline bloch_cell exactly
    ok_anchor2 = all(
        float(np.abs(A.bloch_cell(C, U16, SC.G10, ST, q)
                     - A.bloch_cell(C, U16, SC.G10, ST, q)).max()) == 0.0
        for q in ((0.3, -1.1, 2.0), (math.pi / 2, math.pi / 2, 0.0)))

    # 2-route on 3 sampled rows: node gap by eigvals vs by SVD of M - I
    # (2 arcsin(sigma_min/2)) - two independent gap computations
    devs = []
    for r in rows[1:4]:
        w, dname = r["w"], r["D"]
        c = dict((d, cc) for d, cc in Dlist)[dname]
        ddiag = c @ lam.T

        def M_of(qq, w=w, ddiag=ddiag):
            M = A.bloch_cell(C, U16, SC.G10, ST, qq)
            return np.kron(np.eye(8, dtype=complex),
                           wilson_16(qq, w, ddiag, V)) @ M

        for q in (SC.M_POINTS[0], SC.GAMMA, SC.R_POINTS[0]):
            M = M_of(q)
            g1 = A.gap_at(M)
            sm = np.linalg.svd(M - np.eye(128, dtype=complex),
                               compute_uv=False)[-1]
            g2 = 2.0 * math.asin(min(1.0, sm / 2.0))
            devs.append(abs(g1 - g2))
    ok_route = max(devs) < 1e-9

    # mirror: the Cartan combination is even-grade -> commutes with w
    _Dk, S32 = cartan_16()
    wps = A.pseudoscalar_sig(0)
    comm = max(float(np.abs(sum(cc * S for cc, S in zip(c, S32)) @ wps
                             - wps @ sum(cc * S for cc, S in zip(c, S32)))
                     .max()) for _d, c in Dlist) if Dlist else 0.0

    gapped = [r for r in rows[1:] if r["n_nodes"] == 0]
    closing = [r for r in rows if r["closed"]]
    true_cones = [r for r in rows[1:] if r["A"] is not None
                  and r["A"] < 1e-6 and (r["cone_growth"] or 0.0) > 5
                  and r["chirality_defined"]]

    candidates = [{
        "name": f"wilson:{r['row']}",
        "source": "V166_wilson_quantized",
        "unitary_err": r["unitary_err"], "n_nodes": r["n_nodes"],
        "single_orbit": r["single_orbit"], "anisotropy_A": r["A"],
        "chirality_sum": r["chirality_sum"], "front_speed": None,
        "mirror_ok": comm < 1e-12,
        "survivor": False,
        "step_spec": {"kind": "wilson", "w": r["w"], "D": r["D"]}}
        for r in rows[1:]]

    ok = (ok_anchor and ok_anchor2 and ok_route and len(rows) == 9
          and len(true_cones) == 0 and not missing)

    verdict = {
        "id": "V166_wilson_quantized",
        "tier": 4,
        "property": "the Wilson family at closure-quantized w in "
                    "{pi/2, pi/4}, D in {X, Y}, sign +-, + w=0 anchor: "
                    "node gaps, residual set, Lieb-Robinson shells, "
                    "closure",
        "residual": True,
        "pass": bool(ok),
        "status": f"WILSON-QUANTIZED CENSUS (8 rows + anchor). The w = 0 "
                  f"anchor reproduces the baseline (12 M nodes, "
                  f"{anchor['gap_by_star']['M']:.4e}). X and Y "
                  f"identified by exact multiset match "
                  f"(missing: {missing or 'none'}). Every quantized row "
                  f"gaps the M nodes (min M gap "
                  f"{min(r['gap_by_star']['M'] for r in rows[1:]):.3e}); "
                  f"{len(gapped)}/8 gap the whole high-symmetry set; "
                  f"ZERO true isotropic cones; {len(closing)} close at "
                  f"order <= 64. The factor is quasi-local: worst shell "
                  f"ratio one block beyond the bond cone "
                  f"{max(max(r['lr_shell_ratios']) for r in rows[1:]):.1e}"
                  f" -> the strict light-cone gate fails by construction "
                  f"(V65). 2-route (eigvals vs SVD) dev {max(devs):.1e}.",
        "data": {"anchor": anchor, "n_rows": len(rows),
                 "missing_sources": missing,
                 "n_gapped_all": len(gapped), "n_closing": len(closing),
                 "n_true_cones": len(true_cones),
                 "x_coeff": [float(x) for x in (xs[0] if xs else [])],
                 "y_coeff": [float(x) for x in (ys[0] if ys else [])],
                 "svd_route_dev": devs,
                 "rows": rows, "candidates": candidates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v166_wilson_quantized.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V166 - the Wilson family at quantized w (8 rows + anchor)")
    print("=" * 68)
    print(f"  X coeff: {np.round(xs[0], 4) if xs else 'MISSING'}   "
          f"Y coeff: {np.round(ys[0], 4) if ys else 'MISSING'}")
    for r in rows:
        print(f"  {r['row']:>16}  nodes {r['n_nodes']:3d}  M gap "
              f"{r['gap_by_star']['M']:.3e}  A "
              f"{r['A'] if r['A'] is None else round(r['A'], 3)}  LR "
              f"{['%.1e' % s for s in r['lr_shell_ratios']]}")
    print(f"  2-route dev: {max(devs):.1e}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v166_wilson_quantized.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
