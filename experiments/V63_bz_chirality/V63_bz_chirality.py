#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V63 - the band structure of the automaton cell: is the 16 a **single**
chiral sector, or does the Nielsen-Ninomiya doubling survive the dynamics (W9)?

**What had never been done.**  The publication claimed "no doublers" as an
argument (the module is on-site, so the Nielsen-Ninomiya theorem does not
apply *as stated*) and named the Brillouin-zone analysis as the open item.
This experiment **runs that analysis** on the cell of the automaton: the exact
one-step Bloch matrix of ``cl10.supercell``, validated against the real-space
step on a 2^3 block torus, scanned over the block zone.

**The six results.**
1. **Chirality is a dynamical symmetry.**  The bond stages act on the
   *sub-site* index and the internal layer is even, so the step commutes with
   the pseudoscalar: the 16 and the 16-bar are invariant blocks of the
   256-dimensional one-step matrix (leakage ~1e-15).
2. **The Bloch decomposition is exact.**  The step is one operator, so the
   one-step matrix on the block torus is the direct sum of its Bloch blocks:
   the 1024 eigenphases of the real-space 2^3 step equal the union of the 8
   128-level spectra at q in {0, pi}^3.
3. **The zone center is gapped.**  min|omega(GAMMA)| = 0.0365 per step: the
   would-be massless matter is NOT at q = 0.
4. **The gapless points are the M points.**  12 gapless nodes at
   (+-pi/2, +-pi/2, 0) and permutations - the 12 edge-midpoint momenta of the
   block zone, one star under O_h (little group C2v) - with linear dispersion
   |v| ~ 0.14.
5. **Nielsen-Ninomiya is satisfied**: the 12 node chiralities sum to 0.  The
   doubling is present and the theory has no mechanism (Wilson term, overlap,
   domain wall) that removes it.
6. **The R corners are one star**: the 8 signed corners (+-pi/2, +-pi/2, +-pi/2)
   form a single O_h orbit (little group C3v, order 6), and at every one of
   them the lowest multiplet is exactly 3-fold degenerate - the only exact
   degeneracy of the zone and the only exactly massless multiplet of the step.

**The named residual (W9).**  Removing the half-zone nodes - or showing that the
gapped zone-center sector, not the gapless node momenta, is the physical
fermion sector.
"""

import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np

from cl10 import clifford
from cl10 import engine as E
from cl10 import supercell as SC

TOL = 1e-9


def sector_closure():
    """The 16 <-> 16-bar off-diagonal block of the one-step matrix, and the
    16-sector block itself (must equal the 128-level Bloch matrix)."""
    B, U32, U16 = SC.sectors()
    q = (0.31, -0.77, 0.44)
    M = SC.bloch(q, U32)
    w = np.asarray(clifford.pseudoscalar(), dtype=np.complex128)
    ev, evec = np.linalg.eig(w)
    Bb = np.linalg.qr(evec[:, np.argsort(np.imag(ev))[:16]])[0]
    V16 = np.kron(np.eye(SC.NSUB), B)
    V16b = np.kron(np.eye(SC.NSUB), Bb)
    leak = float(np.abs(V16.conj().T @ M @ V16b).max())
    block_err = float(np.abs(V16.conj().T @ M @ V16 - SC.bloch(q, U16)).max())
    return leak, block_err


def torus_validation():
    """The real-space one-step operator on the 2^3 block torus vs the union of
    its 8 Bloch blocks (the step is one operator, so the equality is exact)."""
    _, _, U16 = SC.sectors()
    Lb = 2
    dim = Lb ** 3 * SC.NSECTOR_ONSITE
    F = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        f = np.zeros((Lb, Lb, Lb, SC.NSUB, SC.NSECTOR), dtype=complex)
        f.flat[k] = 1.0
        F[:, k] = SC.step(f, U16).ravel()
    a = np.sort(np.angle(np.linalg.eigvals(F)))
    b = []
    for qx in (0.0, math.pi):
        for qy in (0.0, math.pi):
            for qz in (0.0, math.pi):
                b.append(SC.phases(SC.bloch((qx, qy, qz), U16)))
    b = np.sort(np.concatenate(b))
    return float(np.abs(a - b).max())


def nodes():
    """The 12 half-zone representatives: gap, axis slopes, chirality."""
    _, _, U16 = SC.sectors()
    out = []
    for q in SC.M_POINTS:
        g = SC.gap(q, U16)
        sl = []
        for ax in range(3):
            if abs(q[ax]) < 1e-12:
                continue                       # the flat direction of the node
            dq = 1e-4
            qa = list(q); qa[ax] -= dq
            qb = list(q); qb[ax] += dq
            pa = np.angle(np.linalg.eigvals(SC.bloch(tuple(qa), U16)))
            pb = np.angle(np.linalg.eigvals(SC.bloch(tuple(qb), U16)))
            sl.append(float((pb[np.argmin(np.abs(pb))]
                             - pa[np.argmin(np.abs(pa))]) / (2 * dq)))
        out.append({"q": [round(v, 6) for v in q], "gap": g,
                    "slope": [round(v, 5) for v in sl],
                    "chirality": float(np.sign(sl[0] * sl[1]))})
    return out


def main():
    t0 = time.time()
    _, _, U16 = SC.sectors()
    leak, block_err = sector_closure()
    dev = torus_validation()
    gap0 = SC.gap(SC.GAMMA, U16)
    _, U32, _ = SC.sectors()
    gap0_32 = SC.gap(SC.GAMMA, U32)
    nd = nodes()
    chirality_sum = float(sum(n["chirality"] for n in nd))
    v_typ = float(np.mean([abs(s) for n in nd for s in n["slope"]]))
    max_M_gap = float(max(n["gap"] for n in nd))

    scan = []
    for g10 in (0.2, 0.15, 0.1):
        _, U32g, U16g = SC.sectors(g10)
        scan.append({"g10": g10,
                     "gap_q0": SC.gap(SC.GAMMA, U16g, g10),
                     "gap_node": SC.gap((math.pi / 2, 0.0, -math.pi / 2), U16g, g10)})

    # the R corner: one star, the full O_h, an exactly 3-fold lowest multiplet
    r_mult = sorted({SC.degeneracy(q, U16)[0][1] for q in SC.R_POINTS})
    r_low = [SC.degeneracy(q, U16)[0][0] for q in SC.R_POINTS]
    r_gap_spread = float(max(r_low) - min(r_low))

    ok = (leak < TOL and block_err < TOL and dev < 1e-9
          and gap0 > 1e-3 and max_M_gap < 1e-3
          and abs(chirality_sum) < 1e-9 and len(nd) == 12
          and r_mult == [3] and SC.stabilizer_order(SC.R_POINTS[0]) == 6
          and SC.stabilizer_order(SC.M_POINTS[0]) == 4)

    verdict = {
        "id": "V63_bz_chirality",
        "tier": 2,
        "property": "the 16-sector one-step spectrum over the block Brillouin zone",
        "residual": True,          # the half-zone nodes are present: W9 is NOT closed
        "pass": bool(ok),
        "status": "chirality_is_a_dynamical_symmetry_but_the_half_zone_nodes_are_"
                  "present_12_gapless_node_momenta_chirality_sum_0_and_the_R_"
                  "points_carry_exactly_3fold_cubic_multiplets",
        "tol": TOL,
        "data": {
            "backend": E._BACKEND,
            "g10": float(SC.G10),
            "stages_per_step": len(SC.STAGES),
            "sector_leakage_16_to_16bar": leak,
            "sector_block_err": block_err,
            "bloch_vs_torus_dev": dev,
            "gap_q0_16": gap0,
            "gap_q0_32": gap0_32,
            "n_doubler_nodes": len(nd),
            "max_M_node_gap": max_M_gap,
            "v_typical": v_typ,
            "chirality_sum": chirality_sum,
            "nodes": nd,
            "g10_scan": scan,
            "M_star": [list(q) for q in SC.M_STAR],
            "R_star": [list(q) for q in SC.R_STAR],
            "little_group_M": SC.stabilizer_order(SC.M_POINTS[0]),
            "little_group_R": SC.stabilizer_order(SC.R_POINTS[0]),
            "R_lowest_multiplicity": r_mult,
            "R_gap_spread_over_the_star": r_gap_spread,
            "residual": "the half-zone nodes are not removed (no Wilson/overlap/"
                        "domain-wall term in the theory); the gap at q=0 is not "
                        "protected in g_10",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v63_bz_chirality.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V63 - the band structure of the automaton cell (the W9 doubling check)")
    print("=" * 68)
    print(f"  backend:  {E._BACKEND}   g_10 = {SC.G10:.4f}   "
          f"{len(SC.STAGES)} stages per step")
    print(f"  [1] 16-sector closure of the step:  16->16bar leakage = {leak:.2e}, "
          f"16 block = the 128-level Bloch matrix to {block_err:.2e}")
    print(f"  [2] the Bloch blocks vs the real-space 2^3 block torus:  "
          f"max deviation = {dev:.2e}")
    print(f"  [3] the zone center:  min|omega(GAMMA)| = {gap0:.5f} (16), "
          f"{gap0_32:.5f} (32)  ->  GAPPED: the matter is not at q=0")
    print(f"  [4] the M points (+-pi/2, +-pi/2, 0):  {len(nd)} nodes, "
          f"max gap = {max_M_gap:.2e}, |v| = {v_typ:.4f}, "
          f"sum chirality = {chirality_sum:.0f}")
    print("  [5] g_10 scan (the gap is not protected):  "
          + "  ".join(f"g10={s['g10']}: {s['gap_q0']:.4f}/{s['gap_node']:.4f}"
                      for s in scan))
    print(f"  [6] the R corners:  star {len(SC.R_STAR)}, little group "
          f"{SC.stabilizer_order(SC.R_POINTS[0])}, lowest multiplet "
          f"{r_mult[0]}-fold at every corner;  M: star {len(SC.M_STAR)}, "
          f"little group {SC.stabilizer_order(SC.M_POINTS[0])}")
    print(f"  VERDICT:  {'PASS (the doubling is present at the M points; W9 stays open)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v63_bz_chirality.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
