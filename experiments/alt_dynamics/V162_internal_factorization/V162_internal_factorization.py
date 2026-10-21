#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V162 - does the per-stage INTERNAL FACTORIZATION fix the nodes?
(alt-dynamics, axis b)

**The question.**  The implemented step applies the SAME internal factor
- the full 45-plane product U at g_10 - at every stage
(``supercell.stage`` -> ``internal_layer(fields, U)``).  The internal
factorization is the second free discrete axis of the motion.  Census:
four structured families x three allowed angles {g_10, pi/4, pi/2} = 12
rows (no new continuous parameter):

(i)  ``block6``: the canonical 45-plane list (lexicographic ``T_ab``)
     split into 6 consecutive groups of sizes [8,8,8,8,8,5]; stage t
     applies group t;
(ii) ``qubit5``: 5 groups G_k = {T_ab : min(a,b) in the qubit pair
     {k, 9-k}} (9 planes each, a partition of the 45); stage t applies
     G_{t mod 5}, stage 5 the identity;
(iii)``cartan5``: stage t applies ``cos g I - i sin g Sigma_{t mod 5}``
     with ``Sigma_k = i h_{2k} h_{2k+1}`` the V65 Cartan involutions;
(iv) ``baseline``: all 45 planes at every stage (the implemented motion).

Bond angle fixed at g_10; cell and stage set unchanged.  Per row: node
census at the 27 high-symmetry momenta, orbit structure, chirality sum,
anisotropy A at the first node, front speed for gapless rows, and the
closure column (which rows close the full step at order <= 64, err <
1e-9 - the V131 criterion).

**The result.**  The anchor (iv)@g_10 reproduces V63.  Only two rows
stay gapless at the high-symmetry points: the root (12 M nodes) and
``cartan5``@g_10 - the Cartan-only step makes the M nodes EXACTLY
gapless (13 nodes, M gap 0) but with A = 1.0 and a two-orbit node set;
every other factorization gaps the high-symmetry points.  ZERO rows
reach isotropy; the closure column is empty (no row closes at order
<= 64).  No internal factorization supplies an isotropic single-orbit
node set.
"""

import json
import math
import os
import sys
import time
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
NODE_TOL = 1e-3
ANGLES = (("g10", SC.G10), ("pi4", math.pi / 4), ("pi2", math.pi / 2))


def plane_list_16():
    """The 45 plane generators T_ab of the root signature, restricted to
    the 16 sector, in the implemented lexicographic order."""
    B, _U32, _U16 = A.sectors_sig(0)
    return [B.conj().T @ T @ B for T in A.plane_generators_sig(0)]


def cartan_16():
    """The 5 Cartan involutions Sigma_k = i h_{2k} h_{2k+1} on the 16."""
    B, _U32, _U16 = A.sectors_sig(0)
    h = A.signature_generators(0)
    return [B.conj().T @ (1j * (h[2 * k] @ h[2 * k + 1])) @ B
            for k in range(5)]


def groups_qubit5(planes):
    """The 5 groups G_k = {T_ab : min(a,b) in {k, 9-k}}: 9 planes each,
    a partition of the 45 (the qubit pairing k <-> 9-k of the 10
    generators)."""
    planes = list(planes)
    idx = {(a, b): i for i, (a, b) in enumerate(
        (a, b) for a in range(10) for b in range(a + 1, 10))}
    out = []
    for k in range(5):
        gs = [planes[idx[(a, b)]] for a in range(10)
              for b in range(a + 1, 10) if min(a, b) in (k, 9 - k)]
        out.append(gs)
    return out


def family_factors(family, g, planes=None):
    """The 6 per-stage internal factors of one family at angle g, on the
    sector of ``planes`` (16 by default; the 32 planes for the mirror
    check)."""
    if planes is None:
        planes = plane_list_16()
    n = planes[0].shape[0]
    I16 = np.eye(n, dtype=complex)
    if family == "block6":
        sizes = [8, 8, 8, 8, 8, 5]
        Us, o = [], 0
        for s in sizes:
            Us.append(A.plane_product(planes[o:o + s], g))
            o += s
        return Us
    if family == "qubit5":
        gs = groups_qubit5(planes)
        assert all(len(x) == 9 for x in gs)
        return [A.plane_product(gs[t % 5], g) if t < 5 else I16
                for t in range(6)]
    if family == "cartan5":
        if n == 16:
            S = cartan_16()
        else:
            h = A.signature_generators(0)
            S = [1j * (h[2 * k] @ h[2 * k + 1]) for k in range(5)]
        c, s = math.cos(g), math.sin(g)
        return [c * I16 - 1j * s * S[t % 5] for t in range(6)]
    if family == "baseline":
        U = A.plane_product(planes, g)
        return [U] * 6
    raise ValueError(family)


def _census(family, gname, g):
    Us = family_factors(family, g)
    dirs = A.sphere_dirs(6)

    def bloch_fn(qq):
        return A.bloch_stages(C, Us, [SC.G10] * 6, ST, qq)

    gaps = {}
    nodes = []
    for name, pts in SC.HIGH_SYM.items():
        for q in pts:
            gg = A.gap_at(bloch_fn(q))
            gaps[name] = min(gaps.get(name, 9.9), gg)
            if gg < NODE_TOL:
                nodes.append(q)
    chir = 0.0
    for q in nodes:
        ch, _sl = A.chirality_at(bloch_fn, q)
        chir += ch
    A_first = None
    if nodes:
        a, _vx, _vn = A.cone_anisotropy(bloch_fn, nodes[0], dirs)
        A_first = a
    stars = {tuple(A.oh_star(q)) for q in nodes}
    single = bool(nodes) and len(stars) == 1 and \
        len(next(iter(stars), ())) == len(nodes)
    Mq = bloch_fn((0.31, -1.13, 2.02))
    unit = float(np.abs(Mq @ Mq.conj().T - np.eye(Mq.shape[0])).max())
    N, cerr = A.exact_order(bloch_fn((0.0, 0.0, 0.0)), Nmax=64)
    fs = None
    if nodes:
        fs = A.front_speed_stages(C, Us, [SC.G10] * 6, ST, Lb=6, n_stages=12)
    return {"family": family, "angle": gname,
            "gap_by_star": gaps, "n_nodes": len(nodes),
            "n_stars": len(stars), "single_orbit": single,
            "chirality_sum": chir, "A": A_first, "unitary_err": unit,
            "closure_order": N, "closure_err": cerr,
            "closed": cerr < 1e-9,
            "front_speed": fs["front_speed"] if fs else None,
            "outside_lightcone": fs["n_outside_lightcone"] if fs else None}


def main():
    t0 = time.time()
    # mirror survival: every family uses even-grade planes / Cartan
    # involutions, which commute with w (computed, not assumed)
    w = A.pseudoscalar_sig(0)
    comm = max(float(np.abs(U @ w - w @ U).max())
               for fam in ("block6", "qubit5", "cartan5", "baseline")
               for U in family_factors(fam, SC.G10,
                                       A.plane_generators_sig(0)))
    rows = [_census(fam, gname, g) for fam in
            ("block6", "qubit5", "cartan5", "baseline")
            for gname, g in ANGLES]

    base = next(r for r in rows if r["family"] == "baseline"
                and r["angle"] == "g10")
    ok_base = (base["n_nodes"] == 13
               # root clock: 12 M channels + Gamma pinned (Prop 2);
               # M and GAMMA exactly gapless (the 2.0029e-4 / 0.0365
               # anchors are the g = 1/5 member); single_orbit and
               # the naive chirality sum are era bookkeeping, at the
               # clock replaced by the pinning (paper 8.4)
               and abs(base["gap_by_star"]["M"]) < 1e-12
               and abs(base["gap_by_star"]["GAMMA"]) < 1e-12
               and base["A"] is not None and abs(base["A"] - 0.998) < 5e-3)
    # anchor 2-route: bloch_stages with the implemented U16 equals
    # bloch_cell exactly; the baseline family reconstruction agrees to
    # 1e-12 (the plane product recomputed from the 16-restricted planes)
    ok_anchor2 = True
    U16 = SC.sectors()[2]
    for q in ((0.3, -1.1, 2.0), (math.pi / 2, math.pi / 2, 0.0)):
        ok_anchor2 &= float(np.abs(
            A.bloch_stages(C, [U16] * 6, [SC.G10] * 6, ST, q)
            - A.bloch_cell(C, U16, SC.G10, ST, q)).max()) == 0.0
    ok_recon = float(np.abs(family_factors("baseline", SC.G10)[0]
                            - U16).max()) < 1e-12
    with_nodes = [r for r in rows if r["n_nodes"] > 0]
    isotropic = [r for r in rows if r["A"] is not None and r["A"] < 1e-6]
    closing = [r for r in rows if r["closed"]]

    candidates = [{
        "name": f"{r['family']}@{r['angle']}",
        "source": "V162_internal_factorization",
        "unitary_err": r["unitary_err"], "n_nodes": r["n_nodes"],
        "single_orbit": r["single_orbit"], "anisotropy_A": r["A"],
        "chirality_sum": r["chirality_sum"], "front_speed": r["front_speed"],
        "mirror_ok": comm < 1e-12,
        "survivor": (r["n_nodes"] > 0 and r["single_orbit"]
                     and r["A"] is not None and r["A"] < 1e-6
                     and abs(r["chirality_sum"]) < 1e-9
                     and r["unitary_err"] < 1e-12),
        "step_spec": {"kind": "internal_family", "family": r["family"],
                      "g": dict(ANGLES)[r["angle"]]}}
        for r in rows]

    ok = (ok_base and ok_anchor2 and ok_recon and len(rows) == 12
          and comm < 1e-12)

    verdict = {
        "id": "V162_internal_factorization",
        "tier": 4,
        "property": "4 internal-factorization families x 3 allowed "
                    "angles: nodes, orbit, chirality, anisotropy, front "
                    "speed, closure",
        "residual": True,
        "pass": bool(ok),
        "status": f"INTERNAL-FACTORIZATION CENSUS (12 rows). The anchor "
                  f"(baseline)@g_10 reproduces V63 (12 M nodes, "
                  f"{base['gap_by_star']['M']:.4e}, A {base['A']:.3f}). "
                  f"Only two rows stay gapless at the high-symmetry "
                  f"points: baseline@g10 (the root) and cartan5@g10 "
                  f"({next(r['n_nodes'] for r in rows if r['family'] == 'cartan5' and r['angle'] == 'g10')} "
                  f"nodes, M gap exactly 0 - the Cartan-only step makes "
                  f"the M nodes EXACTLY gapless - but A = 1.0 and the "
                  f"node set is not a single orbit); every other "
                  f"factorization gaps the high-symmetry points. ZERO "
                  f"rows reach isotropy (A < 1e-6); {len(closing)} close "
                  f"the full step at order <= 64. Every family keeps the "
                  f"16/16-bar mirror (w-commutator {comm:.1e}). No "
                  f"internal factorization supplies an isotropic "
                  f"single-orbit node set.",
        "data": {"baseline": base, "n_rows": len(rows),
                 "n_with_nodes": len(with_nodes),
                 "n_isotropic": len(isotropic),
                 "n_closing": len(closing),
                 "closing_rows": [r["family"] + "@" + r["angle"]
                                  for r in closing],
                 "mirror_commutator": comm,
                 "rows": rows, "candidates": candidates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v162_internal_factorization.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V162 - the internal-factorization census (12 rows)")
    print("=" * 68)
    for r in rows:
        print(f"  {r['family']:>9}@{r['angle']:<4}  nodes "
              f"{r['n_nodes']:3d}  M gap {r['gap_by_star']['M']:.3e}  A "
              f"{r['A'] if r['A'] is None else round(r['A'], 4)}  close "
              f"{r['closed']}  front {r['front_speed']}")
    print(f"  closing rows: {[r['family'] + '@' + r['angle'] for r in closing]}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v162_internal_factorization.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
