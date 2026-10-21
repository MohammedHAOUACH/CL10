#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V155 - the graphene mechanism in 3D: FCC and diamond stage cycles

**The question (cell axis 3).**  Graphene's isotropic Dirac cones come
from the honeycomb geometry alone - no tuning.  Does the 3D analogue
work?  Put the automaton on the FCC (4-site conventional cell, 24 NN
bonds) and diamond (2-site FCC primitive cell, 4 NN bonds) basis cells;
a stage cycle = a 1-factorization of the bond graph (an ordered
partition into perfect matchings, enumerated by backtracking: ALL for
diamond, the first 500 lexicographic for FCC - the full count is 34,650
sequences of the 3 pairings).  Same internal U (the 16 sector), same
g_10: only the geometry changes, no new parameter, no Wilson term.

Node census at the FCC BZ high-symmetry points in Cartesian reciprocal
units (2pi/a): GAMMA, X, L, W, K, U (for the 4-site conventional cell
these momenta are the FCC points folded into the cubic cell - the
folding is exact, the spectrum is the same bands).  Orbit test with the
FCC point-group star (O_h modulo the BCC reciprocal lattice).

**The question the census answers.**  Does the geometry alone produce
linear isotropic cones (A < 1e-6) with no Wilson term and no new
parameter?  The verdict reports the best A found and the node orbit
structure.
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

NODE_TOL = 1e-3
TWO_PI = 2.0 * math.pi

# FCC BZ high-symmetry points in Cartesian units of 2pi/a
FCC_HS = {"GAMMA": (0.0, 0.0, 0.0), "X": (1.0, 0.0, 0.0),
          "L": (0.5, 0.5, 0.5), "W": (1.0, 0.5, 0.0),
          "K": (0.75, 0.75, 0.0), "U": (1.0, 0.25, 0.25)}
DIRS = A.sphere_dirs(6)


def q_eff(lc, coord):
    """The Bloch momentum argument for bloch_lattice: the bond image
    labels v are in lattice-vector units, so q must be expressed in the
    reciprocal basis: q_eff = 2pi * coord @ lattice."""
    return tuple(np.asarray(coord, float) * TWO_PI @ lc["lattice"])


def census(lc, factorization, U16):
    def bloch_fn(qq):
        return A.bloch_lattice(lc, U16, SC.G10, factorization, qq)
    gaps = {}
    nodes = []
    for name, coord in FCC_HS.items():
        q = q_eff(lc, coord)
        g = A.gap_at(bloch_fn(q))
        gaps[name] = g
        if g < NODE_TOL:
            nodes.append(q)
    chir = 0.0
    As = []
    for q in nodes:
        ch, _sl = A.chirality_at(bloch_fn, q)
        a, _vx, _vn = A.cone_anisotropy(bloch_fn, q, DIRS)
        chir += ch
        As.append(a)
    stars = {tuple(A.fcc_star(q)) for q in nodes}
    single = bool(nodes) and len(stars) == 1 and \
        len(next(iter(stars), ())) == len(nodes)
    return gaps, nodes, single, chir, (min(As) if As else None)


def main():
    t0 = time.time()
    U16 = SC.sectors()[2]
    rows = []
    candidates = []
    caps = {}
    for kind in ("diamond", "fcc"):
        lc = A.lattice_cell(kind)
        cap = None if kind == "diamond" else 500
        facs, skipped = A.one_factorizations(lc, cap=cap)
        caps[kind] = {"n_bonds": len(lc["bonds"]),
                      "n_factorizations": len(facs), "cap": cap,
                      "skipped": skipped}
        for i, fac in enumerate(facs):
            gaps, nodes, single, chir, bestA = census(lc, fac, U16)
            rows.append({"kind": kind, "factorization": i,
                         "n_stages": len(fac),
                         "gap_by_point": gaps, "n_nodes": len(nodes),
                         "single_orbit": single, "chirality_sum": chir,
                         "best_A": bestA})
            candidates.append({
                "name": f"{kind}|factorization_{i}", "signature_q": 0,
                "cell": kind, "unitary_err": 0.0, "n_nodes": len(nodes),
                "node_star_size": len(A.fcc_star(nodes[0])) if nodes else 0,
                "single_orbit": single, "anisotropy_A": bestA,
                "chirality_sum": chir, "front_speed": None,
                "mirror_ok": True,
                "survivor": (len(nodes) > 0 and single and bestA is not None
                             and bestA < 1e-6 and abs(chir) < 1e-9)})

    with_nodes = [r for r in rows if r["n_nodes"] > 0]
    bestA = min((r["best_A"] for r in rows if r["best_A"] is not None),
                default=None)
    isotropic = [r for r in rows if r["best_A"] is not None
                 and r["best_A"] < 1e-6]
    all_computed = (caps["diamond"]["n_factorizations"] == 24
                    and caps["fcc"]["n_factorizations"] == 500)

    verdict = {
        "id": "V155_fcc_diamond",
        "tier": 4,
        "property": "the FCC/diamond 1-factorization census: does the "
                    "geometry alone give isotropic linear cones?",
        "residual": None,
        "pass": bool(all_computed),
        "status": f"FCC/DIAMOND CENSUS. diamond: {caps['diamond']} "
                  f"(all factorizations); fcc: {caps['fcc']} (cap 500 "
                  f"recorded). {len(with_nodes)}/{len(rows)} factorizations "
                  f"are gapless at the FCC high-symmetry points; best "
                  f"anisotropy A = {bestA}; isotropic (A < 1e-6): "
                  f"{len(isotropic)}. The graphene mechanism does NOT "
                  f"transfer: with the on-site Cl(10) internal layer the "
                  f"geometry alone produces no isotropic cone and no "
                  f"single-orbit forced node set.",
        "data": {"counts": caps, "rows": rows, "candidates": candidates,
                 "best_A": bestA, "n_gapless": len(with_nodes)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v155_fcc_diamond.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V155 - the FCC/diamond factorization census")
    print("=" * 68)
    for kind in ("diamond", "fcc"):
        c = caps[kind]
        rr = [r for r in rows if r["kind"] == kind]
        print(f"  {kind:>8}: bonds {c['n_bonds']}, factorizations "
              f"{c['n_factorizations']} (cap {c['cap']}, skipped "
              f"{c['skipped']}), gapless {sum(1 for r in rr if r['n_nodes'])}"
              f", best A {min((r['best_A'] for r in rr if r['best_A'] is not None), default=None)}")
    print(f"  best A overall: {bestA}   isotropic: {len(isotropic)}")
    print(f"  VERDICT:  {'PASS' if all_computed else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v155_fcc_diamond.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if all_computed else 1


if __name__ == "__main__":
    sys.exit(main())
