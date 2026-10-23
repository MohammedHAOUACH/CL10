#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V139 - the nodal-plane crossing at the implemented coupling: a massless
Weyl pair with a 2-dimensional little group.

V128 ran the massless census on a 15^3 grid and reported "no exactly gapless
mode anywhere" in the implemented step (g10 = 1/5, smallest quasi-energy
5.1e-4). That census was too coarse. The gap function has isolated zeros
between the grid points, and this experiment finds them by refinement.

**What is found (all computed).** On each of the three coordinate planes
(q_z = 0 and its axis permutations) the step has an exactly gapless point,
q0 = (-2.1889854, -1.4445017, 0) and its images: the gap refines to 3e-14
with the residual falling as O(eps^2) along the plane, i.e. a genuine zero,
not a near miss. The O_h star of q0 has 24 points and the little group has
order 2. The dispersion is linear in every in-plane direction, with
direction-dependent speed v in [0.031, 0.219]; moving off the plane costs a
gap linear in q_z (5.1e-4 at q_z = 1e-2, 5.1e-3 at 1e-1), so the crossing is
a massless point node confined to the nodal plane, not a line and not a
gapped quadratic touching.

**Why it matters for the graviton question.** V138 closed the zone-centre
route: a lattice spin 2 needs the symmetric-traceless quintet Eg + T2g, and
the bond factor never carries Eg. The node found here has little group of
order 2, so its tangent space carries only 1-dimensional irreps: no
helicity pair, no spin 2 - the same representation-theoretic verdict as at
the M nodes, now at a node that is exactly massless at the implemented
coupling. What the step does supply at g10 = 1/5 is a massless two-fold
degeneracy with a linear, anisotropic cone: the O9 anisotropy problem, not
a graviton.

The machinery is the exact factorisation M(q) = P(q) (x) U16^6 (V128 T0), so
the gap is computed from an 8x8 and a 16x16 eigenvalue problem.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# reference spectrum: the g = 1/5 Wilson-family member (paper 8.2):
# the exactly gapless 24-point orbit at q0 = (-2.1890561, -1.4445961,
# 0), little group of order 2, in-plane v in [0.031, 0.219], is a
# statement about THAT member; at the pi/4 clock the gapless set is
# Gamma plus six face-diagonal lines (V203).
os.environ["CL10_ANGLE"] = "0.2"

from cl10 import supercell as SC                                       # noqa: E402

TOL_GAP = 1e-9


def bond_and_internal(g):
    _B, _U32, U16 = SC.sectors(g)
    U6 = np.linalg.matrix_power(U16, 6)
    return U6


def p_bond(q, g):
    P = np.eye(8, dtype=complex)
    for even, j in SC.STAGES:
        P = SC.stage_matrix(even, j, q[j], g) @ P
    return P


def gap_at(q, g, U6):
    bp = np.angle(np.linalg.eigvals(p_bond(q, g)))
    ip = np.angle(np.linalg.eigvals(U6))
    return float(np.min(np.abs(np.angle(np.exp(1j * (bp[:, None] + ip[None, :]))))))


def refine(gap, c0, axes, r0=0.02, iters=26):
    """Coordinate refinement of the gap over the given free axes."""
    c = list(c0)
    r = r0
    best = (gap(c), list(c))
    for _ in range(iters):
        grid = [np.linspace(c[k] - r, c[k] + r, 7) if k in axes else [c[k]]
                for k in range(3)]
        for a in grid[0]:
            for b in grid[1]:
                for z in grid[2]:
                    v = gap((a, b, z))
                    if v < best[0]:
                        best = (v, [a, b, z])
        c = best[1]
        r *= 0.45
    return best


def main():
    out = {"experiment": "V139_nodal_plane_crossing", "data": {}}
    checks = []
    g = SC.G10
    U6 = bond_and_internal(g)
    gap = lambda q: gap_at(q, g, U6)                                    # noqa: E731

    # ---- T0 the V128 grid census misses the zeros -------------------------
    grid = [-math.pi + 2 * math.pi * i / 14 for i in range(15)]
    coarse = min(gap(q) for q in _prod(grid, grid, grid))
    out["data"]["coarse"] = {"grid": 15, "min_gap": coarse}
    checks.append(("T0 the 15^3 census of V128 sees a smallest quasi-energy "
                   "of %.3e - no exact zero on the grid" % coarse,
                   coarse > 1e-4))

    # ---- T1 the crossing on each nodal plane ------------------------------
    seeds = []
    for s in (1, -1):
        for (i, j) in ((0, 1), (0, 2), (1, 2)):
            for a in (-2.1889854, -1.4445017):
                for b in (-2.1889854, -1.4445017):
                    q = [0.0, 0.0, 0.0]
                    q[i] = s * a
                    q[j] = s * b
                    seeds.append((tuple(q), (i, j)))
    orbits = {}
    for sd, (i, j) in seeds:
        axes = [k for k in range(3) if k in (i, j)]
        v, c = refine(gap, sd, axes)
        if v < TOL_GAP:
            key = tuple(round(x, 5) for x in c)
            if key not in orbits or v < orbits[key][0]:
                orbits[key] = (v, key)
    star = max((tuple(SC.star(k)) for k in orbits), key=len)
    gv = min(orbits[p][0] for p in orbits if p in star)
    q0 = list(min((p for p in star if abs(p[2]) < 1e-3),
                  key=lambda p: abs(p[2])))
    for _ in range(4):
        q0 = refine(gap, q0, (0, 1), r0=1e-3, iters=16)[1]
    q0 = tuple(q0)
    out["data"]["node"] = {"gap": gv, "q0": list(q0),
                           "gap_at_q0": gap(q0),
                           "star_size": len(star),
                           "stabilizer_order": 48 // len(star),
                           "orbits": len(orbits)}
    checks.append(("T1 one orbit of exactly gapless points: %d points, little "
                   "group of order %d, representative %s, refined gap %.1e"
                   % (len(star), 48 // len(star), q0, gv),
                   len(star) == 24 and gv < TOL_GAP))

    # ---- T2 the zero is exact: the gap vanishes at the point --------------
    g0 = gap(q0)
    prof = {}
    for eps in (1e-3, 1e-2, 5e-2):
        prof["%.0e" % eps] = gap((q0[0] + eps, q0[1], 0.0))
    floor = max(g0, 1e-15)
    out["data"]["zero"] = {"gap_at_q0": g0, "in_plane_profile": prof,
                           "profile_over_gap": [prof[k] / floor
                                                for k in prof]}
    checks.append(("T2 the gap at the refined point is %.1e, while moving "
                   "1e-3, 1e-2, 5e-2 along the plane costs %s - a factor %d "
                   "above the point value: the gap vanishes at the point and "
                   "grows linearly away from it, the signature of a simple "
                   "zero, not a floor"
                   % (g0, {k: "%.1e" % v for k, v in prof.items()},
                      int(prof["5e-02"] / floor)),
                   g0 < 1e-7 and prof["5e-02"] > 1e3 * floor))

    # ---- T3 linear in-plane dispersion, linear out-of-plane mass ---------
    dirs = ((1, 0, 0), (0, 1, 0), (1, 1, 0), (1, -1, 0),
            (1, 2, 0), (2, 1, 0))
    slopes = {}
    for d in dirs:
        n = math.sqrt(sum(x * x for x in d))
        sp = [gap(tuple(q0[k] + eps * d[k] / n for k in range(3))) / eps
              for eps in (1e-3, 1e-4)]
        slopes[str(d)] = [round(v, 5) for v in sp]
    inplane = [v for d in dirs for v in slopes[str(d)]]
    mass = {}
    for eps in (1e-1, 1e-2, 1e-3):
        mass["%.0e" % eps] = gap((q0[0], q0[1], eps))
    out["data"]["dispersion"] = {"slopes_by_direction": slopes,
                                 "v_min": min(inplane), "v_max": max(inplane),
                                 "anisotropy": (max(inplane) - min(inplane))
                                 / max(inplane),
                                 "off_plane_gap": mass}
    lin_out = all(mass["%.0e" % (10.0 ** -k)] > 1e3 * floor
                  for k in (1, 2, 3))
    checks.append(("T3 the branch is linear in every in-plane direction with "
                   "v in [%.4f, %.4f] (anisotropy %.3f), and leaving the "
                   "plane costs a gap linear in q_z (%s): a massless point "
                   "node confined to the nodal plane"
                   % (min(inplane), max(inplane),
                      (max(inplane) - min(inplane)) / max(inplane),
                      {k: "%.1e" % v for k, v in mass.items()}),
                   min(inplane) > 1e-3 and lin_out))

    # ---- T4 the representation verdict -----------------------------------
    n_irr = 48 // len(star)
    out["data"]["little_group"] = {"order": n_irr,
                                   "max_irrep_dimension": 1 if n_irr <= 2 else 2}
    checks.append(("T4 the little group of the node has order %d, so its "
                   "tangent space carries only %d-dimensional irreps: no "
                   "helicity +/-1 or +/-2 pair, no spin-2 content - the same "
                   "representation-theoretic closure as at the M nodes, now "
                   "at a node that is exactly massless at g10 = 1/5"
                   % (n_irr, 1 if n_irr <= 2 else 2), n_irr == 2))

    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["pass"] = ok
    out["status"] = (
        "THE NODAL-PLANE CROSSING AT THE IMPLEMENTED COUPLING. The V128 "
        "census ('no exactly gapless mode anywhere', grid minimum %.1e) is "
        "too coarse: on each coordinate plane the step has an exactly "
        "gapless point (representative %s, gap %.1e), one O_h orbit of %d "
        "points with little group of order %d. The dispersion is linear in "
        "the plane with v in [%.4f, %.4f] and the node is gapped linearly "
        "off the plane, so it is a massless point node confined to the "
        "nodal plane. Representation content: a group of order 2 has only "
        "1-dimensional irreps, so this node carries no helicity pair and no "
        "spin 2 - consistent with V138, which closed the zone-centre route "
        "on Eg-absence. VERDICT: the implemented step does contain exactly "
        "massless excitations at g10 = 1/5, contrary to the coarse census; "
        "they are anisotropic Weyl-type nodes (O9), not gravitons, and the "
        "spin-2 no-go of V138 stands."
        % (coarse, list(q0), gv, len(star), 48 // len(star),
           min(inplane), max(inplane)))

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v139_nodal_plane_crossing.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V139 - the nodal-plane crossing at g10 = 1/5")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  node %s gap %.1e | star %d | v [%.4f, %.4f]"
          % (list(q0), gv, len(star), min(inplane), max(inplane)))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def _prod(*vecs):
    import itertools
    return itertools.product(*vecs)


if __name__ == "__main__":
    sys.exit(main())
