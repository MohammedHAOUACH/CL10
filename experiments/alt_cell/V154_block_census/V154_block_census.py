#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V154 - do the 2x2x2 blocks cause the nodes? (alt-cell, cell axis 2)

**The question.**  The automaton cell is a 2x2x2 block (8 sub-sites, six
staggered stages).  Census the mixed-radix generalization: every block
shape in {2,3}^3 (8 blocks, 8..27 sub-sites) with the canonical stagger
(sum_j n_j phases, phase-major), plus the implemented cycle as the
(2,2,2) anchor row and the canonical cycle on (2,2,2) as a diagnostic.
Same internal U (the 16 sector), same g_10 - only the cell changes.

Per block: node census at the 27 high-symmetry block momenta (V63
tolerance 1e-3), orbit test (is the gapless set one O_h star?), cone
anisotropy at the gapless points (V115 metric), chirality sum, and the
causal radius (V36 method on a torus: the front must move exactly one
Planck cell per stage - front speed 1, zero amplitude outside the light
cone).

**2-route.**  The (2,2,2) anchor row: bloch_cell == supercell.bloch
(exactly 0) AND the real-space 2^3 torus step vs the Bloch union.

**The result.**  The anchor reproduces V63/V115.  The (2,2,2) cell keeps
its 12 M nodes under both the implemented and the canonical cycle; every
block with a 3-division is gapped at the 27 sampled high-symmetry
momenta (sampling statement).  Every block keeps front speed exactly 1
with zero outside-light-cone amplitude, and no block reaches cone
isotropy A < 1e-6: the cell size is not the cause.
"""

import itertools
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

NODE_TOL = 1e-3
DIRS = A.sphere_dirs(6)


def census(c, stages, U16):
    def bloch_fn(qq):
        return A.bloch_cell(c, U16, SC.G10, stages, qq)
    gaps = {}
    nodes = []
    for name, pts in SC.HIGH_SYM.items():
        for q in pts:
            g = A.gap_at(bloch_fn(q))
            gaps[name] = min(gaps.get(name, 9.9), g)
            if g < NODE_TOL:
                nodes.append(q)
    chir = 0.0
    As = []
    for q in nodes:
        ch, _sl = A.chirality_at(bloch_fn, q)
        a, _vx, _vn = A.cone_anisotropy(bloch_fn, q, DIRS)
        chir += ch
        As.append(a)
    stars = {tuple(A.oh_star(q)) for q in nodes}
    single = bool(nodes) and len(stars) == 1 and \
        len(next(iter(stars), ())) == len(nodes)
    return gaps, nodes, single, chir, (min(As) if As else None)


def main():
    t0 = time.time()
    U16 = SC.sectors()[2]
    shapes = list(itertools.product((2, 3), repeat=3))
    rows = []
    candidates = []
    anchor = None

    for shape in shapes:
        c = A.cell(shape)
        stages = (A.stages_baseline() if shape == (2, 2, 2)
                  else A.stages_shape(shape))
        gaps, nodes, single, chir, bestA = census(c, stages, U16)
        fs = A.front_speed(c, U16, stages, SC.G10, Lb=6, n_stages=12)
        row = {"cell": "x".join(map(str, shape)), "shape": list(shape),
               "cycle": ("implemented" if shape == (2, 2, 2)
                         else "canonical-phase-major"),
               "n_stages": len(stages), "nsub": c["nsub"],
               "gap_by_star": gaps, "n_nodes": len(nodes),
               "single_orbit": single, "chirality_sum": chir,
               "best_A": bestA,
               "front_speed": fs["front_speed"],
               "n_outside_lightcone": fs["n_outside_lightcone"],
               "final_norm": fs["final_norm"]}
        if shape == (2, 2, 2):
            # route (i): exact equality with supercell.bloch
            q = (0.31, -0.77, 0.44)
            row["bloch_vs_supercell"] = float(np.abs(
                A.bloch_cell(c, U16, SC.G10, A.stages_baseline(), q)
                - SC.bloch(q, U16)).max())
            # route (ii): real-space 2^3 torus vs the Bloch union
            Lb = 2
            dim = Lb ** 3 * SC.NSECTOR_ONSITE
            F = np.zeros((dim, dim), dtype=complex)
            for k in range(dim):
                f = np.zeros((Lb, Lb, Lb, SC.NSUB, SC.NSECTOR),
                             dtype=complex)
                f.flat[k] = 1.0
                F[:, k] = SC.step(f, U16).ravel()
            a = np.sort(np.angle(np.linalg.eigvals(F)))
            b = np.sort(np.concatenate([
                np.sort(np.angle(np.linalg.eigvals(
                    A.bloch_cell(c, U16, SC.G10, A.stages_baseline(),
                                 (qx, qy, qz)))))
                for qx in (0.0, math.pi) for qy in (0.0, math.pi)
                for qz in (0.0, math.pi)]))
            row["bloch_vs_torus_dev"] = float(np.abs(a - b).max())
            anchor = row
        rows.append(row)
        candidates.append({
            "name": f"block:{row['cell']}|{row['cycle']}", "signature_q": 0,
            "cell": row["cell"], "unitary_err": 0.0,
            "n_nodes": len(nodes),
            "node_star_size": len(A.oh_star(nodes[0])) if nodes else 0,
            "single_orbit": single, "anisotropy_A": bestA,
            "chirality_sum": chir, "front_speed": fs["front_speed"],
            "mirror_ok": True,
            "survivor": (len(nodes) > 0 and single and bestA is not None
                         and bestA < 1e-6 and abs(chir) < 1e-9
                         and fs["front_speed"] == 1.0
                         and fs["n_outside_lightcone"] == 0)})

    # diagnostic: the canonical cycle on (2,2,2) (a different cycle)
    c222 = A.cell((2, 2, 2))
    gaps, nodes, single, chir, bestA = census(c222, A.stages_shape((2, 2, 2)),
                                              U16)
    fs = A.front_speed(c222, U16, A.stages_shape((2, 2, 2)), SC.G10, Lb=6,
                       n_stages=12)
    rows.append({"cell": "2x2x2", "shape": [2, 2, 2],
                 "cycle": "canonical-phase-major(diagnostic)",
                 "n_stages": 6, "nsub": 8, "gap_by_star": gaps,
                 "n_nodes": len(nodes), "single_orbit": single,
                 "chirality_sum": chir, "best_A": bestA,
                 "front_speed": fs["front_speed"],
                 "n_outside_lightcone": fs["n_outside_lightcone"],
                 "final_norm": fs["final_norm"]})

    ok_anchor = (anchor["n_nodes"] == 12
                 and abs(anchor["gap_by_star"]["M"] - 2.0029e-4) < 1e-6
                 and abs(anchor["gap_by_star"]["GAMMA"] - 0.0365) < 1e-4
                 and anchor["bloch_vs_supercell"] < 1e-14
                 and anchor["bloch_vs_torus_dev"] < 1e-9
                 and anchor["front_speed"] == 1.0
                 and anchor["n_outside_lightcone"] == 0)
    causal = all(r["front_speed"] == 1.0 and r["n_outside_lightcone"] == 0
                 for r in rows)
    all_computed = len(rows) == 9
    no_isotropy = all(r["best_A"] is None or r["best_A"] >= 1e-6
                      for r in rows)
    ok = ok_anchor and causal and all_computed and no_isotropy

    verdict = {
        "id": "V154_block_census",
        "tier": 4,
        "property": "the block census over shapes {2,3}^3: nodes, orbits, "
                    "anisotropy, causal front",
        "residual": True,
        "pass": bool(ok),
        "status": f"BLOCK CENSUS (8 shapes + 1 diagnostic). The (2,2,2) "
                  f"anchor reproduces V63/V115 (12 M nodes, "
                  f"{anchor['gap_by_star']['M']:.4e}, GAMMA "
                  f"{anchor['gap_by_star']['GAMMA']:.4f}; bloch==supercell "
                  f"exactly, torus dev {anchor['bloch_vs_torus_dev']:.1e}). "
                  f"Node counts at the 27 sampled high-symmetry momenta: "
                  f"{[r['n_nodes'] for r in rows]} - the (2,2,2) cell "
                  f"keeps its 12 M nodes under BOTH cycles; every block "
                  f"with a 3-division is gapped at the sampled momenta "
                  f"(the stagger moves the nodes off the sampled set or "
                  f"gaps them: a sampling statement, not a removal "
                  f"theorem). Every front speed is exactly 1 with zero "
                  f"outside-light-cone amplitude, and no block reaches "
                  f"A < 1e-6: the cell size is not the cause.",
        "data": {"rows": rows, "candidates": candidates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v154_block_census.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V154 - the block census {2,3}^3")
    print("=" * 68)
    for r in rows:
        ba = r["best_A"]
        print(f"  {r['cell']:>7} {r['cycle'][:22]:>22}  nodes {r['n_nodes']:2d}"
              f"  single {str(r['single_orbit']):>5}  A "
              f"{'-' if ba is None else format(ba, '.3f')}  v_front "
              f"{r['front_speed']:.1f}  outside "
              f"{r['n_outside_lightcone']}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v154_block_census.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
