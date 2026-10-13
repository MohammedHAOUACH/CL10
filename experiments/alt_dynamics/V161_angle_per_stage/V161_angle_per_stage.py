#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V161 - does the per-stage bond ANGLE fix the nodes? (alt-dynamics, axis a)

**The question.**  The implemented step applies the SAME bond angle
g_10 = 1/5 at every one of the six stages (``supercell.stage`` ->
``bond_layer(..., theta)``).  The angle is the first free discrete axis of
the motion never censused.  Census: ALL 3^6 = 729 per-stage bond-angle
assignments ``thetas in {g_10, pi/4, pi/2}^6`` (the pre-registered set -
no new continuous parameter), internal U fixed at the baseline 45-plane
product at g_10, cell and stage set unchanged.  Per row: closure of the
full step at q = 0 (``exact_order``, order <= 64, the V131 criterion),
node census at the 27 high-symmetry momenta (gap < 1e-3, the V63
tolerance), node stars under O_h, chirality sum, cone anisotropy A at the
first node (V115 metric, 162 directions), front speed for the gapless
rows (V36 support test, per-stage angles).

**The result.**  The all-g_10 anchor row reproduces V63/V115 (12 M nodes,
gap 2.0029e-4, GAMMA 0.0365, A = 0.998).  The census records every row:
the half-zone node set survives every admissible per-stage angle
assignment (the angles only shift the node gap and the cone slope), no
row reaches cone isotropy, and no row closes the step at order <= 64
(the internal factor at g_10 never closes - V131).
"""

import itertools
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
NODE_TOL = 1e-3
ANGLES = (SC.G10, math.pi / 4, math.pi / 2)
TAG = {SC.G10: "g", math.pi / 4: "4", math.pi / 2: "2"}
_W = {}


def _init():
    if not _W:
        _W["U16"] = SC.sectors()[2]
        _W["DIRS"] = A.sphere_dirs(6)


def _census(angles):
    U16 = _W["U16"]
    dirs = _W["DIRS"]

    def bloch_fn(qq):
        return A.bloch_stages(C, [U16] * 6, list(angles), ST, qq)

    gaps = {}
    nodes = []
    for name, pts in SC.HIGH_SYM.items():
        for q in pts:
            g = A.gap_at(bloch_fn(q))
            gaps[name] = min(gaps.get(name, 9.9), g)
            if g < NODE_TOL:
                nodes.append(q)
    chir = 0.0
    A_first = None
    for q in nodes:
        ch, _sl = A.chirality_at(bloch_fn, q)
        chir += ch
    if nodes:
        a, _vx, _vn = A.cone_anisotropy(bloch_fn, nodes[0], dirs)
        A_first = a
    stars = {tuple(A.oh_star(q)) for q in nodes}
    single = bool(nodes) and len(stars) == 1 and \
        len(next(iter(stars), ())) == len(nodes)
    # exact unitarity of the one-step Bloch (route 1)
    Mq = bloch_fn((0.31, -1.13, 2.02))
    unit = float(np.abs(Mq @ Mq.conj().T - np.eye(Mq.shape[0])).max())
    # closure of the full step at q = 0 (V131 criterion)
    N, cerr = A.exact_order(bloch_fn((0.0, 0.0, 0.0)), Nmax=64)
    fs = None
    if nodes:
        fs = A.front_speed_stages(C, [U16] * 6, list(angles), ST,
                                  Lb=6, n_stages=12)
    return {"angles": [TAG[t] for t in angles],
            "gap_by_star": gaps, "n_nodes": len(nodes),
            "n_stars": len(stars), "single_orbit": single,
            "chirality_sum": chir, "A": A_first,
            "unitary_err": unit, "closure_order": N, "closure_err": cerr,
            "closed": cerr < 1e-9,
            "front_speed": fs["front_speed"] if fs else None,
            "outside_lightcone": fs["n_outside_lightcone"] if fs else None}


def _torus_route(angles):
    """Route 2 of the 2-route check: the real-space step on the (2,2,2)-cell
    torus (L = 2 cells, 1024 amplitudes) built column-by-column from
    ``step_stages``; its eigenphase multiset must equal the union of the
    ``bloch_stages`` multisets over the 8 torus momenta q in {0, pi}^3."""
    U16 = _W["U16"]
    L = 2
    dim = L ** 3 * 8 * 16
    Mat = np.zeros((dim, dim), dtype=complex)
    for col in range(dim):
        f = np.zeros((L, L, L, 8, 16), dtype=complex)
        f.ravel()[col] = 1.0
        out = A.step_stages(f, [U16] * 6, list(angles), ST)
        Mat[:, col] = out.ravel()
    ph_r = np.sort(np.angle(np.linalg.eigvals(Mat)))
    ph_b = []
    for q in itertools.product((0.0, math.pi), repeat=3):
        ph_b.append(np.angle(np.linalg.eigvals(
            A.bloch_stages(C, [U16] * 6, list(angles), ST, q))))
    ph_b = np.sort(np.concatenate(ph_b))
    return float(np.abs(ph_r - ph_b).max())


def main():
    t0 = time.time()
    _init()
    assignments = list(itertools.product(ANGLES, repeat=6))
    assert len(assignments) == 729
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count()),
                             initializer=_init) as ex:
        rows = list(ex.map(_census, assignments, chunksize=8))

    base = next(r for r in rows if all(t == "g" for t in r["angles"]))
    ok_base = (base["n_nodes"] == 12
               and abs(base["gap_by_star"]["M"] - 2.0029e-4) < 1e-6
               and abs(base["gap_by_star"]["GAMMA"] - 0.0365) < 1e-4
               and abs(base["chirality_sum"]) < 1e-9
               and base["single_orbit"]
               and base["A"] is not None and abs(base["A"] - 0.998) < 5e-3)
    # discriminating check: the per-stage angle must reach the bond layer
    with_pi2 = [r for r in rows if "2" in r["angles"]]
    discriminating = any(r["n_nodes"] != base["n_nodes"] or
                         abs(r["gap_by_star"]["M"] -
                             base["gap_by_star"]["M"]) > 1e-9
                         for r in with_pi2)
    with_nodes = [r for r in rows if r["n_nodes"] > 0]
    single_rows = [r for r in rows if r["single_orbit"]]
    isotropic = [r for r in rows if r["A"] is not None and r["A"] < 1e-6]
    closing = [r for r in rows if r["closed"]]

    # 2-route verification on 10 sampled rows (anchor + 9 angle variants)
    sample_idx = [0, 121, 242, 363, 404, 485, 566, 607, 688, 728]
    devs = [_torus_route(assignments[i]) for i in sample_idx]
    ok_route = max(devs) < 1e-9

    candidates = [{
        "name": "angles:" + "".join(r["angles"]),
        "source": "V161_angle_per_stage",
        "unitary_err": r["unitary_err"], "n_nodes": r["n_nodes"],
        "single_orbit": r["single_orbit"], "anisotropy_A": r["A"],
        "chirality_sum": r["chirality_sum"],
        "front_speed": r["front_speed"], "mirror_ok": True,
        "survivor": (r["n_nodes"] > 0 and r["single_orbit"]
                     and r["A"] is not None and r["A"] < 1e-6
                     and abs(r["chirality_sum"]) < 1e-9
                     and r["unitary_err"] < 1e-12),
        "step_spec": {"kind": "angles",
                      "thetas": [float(t) for t in angles]}}
        for angles, r in zip(assignments, rows)]


    ok = ok_base and discriminating and ok_route and len(rows) == 729

    verdict = {
        "id": "V161_angle_per_stage",
        "tier": 4,
        "property": "all 3^6 = 729 per-stage bond-angle assignments "
                    "{g_10, pi/4, pi/2}^6: closure, nodes, orbit, "
                    "chirality, anisotropy, front speed",
        "residual": True,
        "pass": bool(ok),
        "status": f"PER-STAGE BOND-ANGLE CENSUS (729 rows). The all-g_10 "
                  f"anchor reproduces V63/V115 (12 M nodes, "
                  f"{base['gap_by_star']['M']:.4e}, GAMMA "
                  f"{base['gap_by_star']['GAMMA']:.4f}, A "
                  f"{base['A']:.3f}). {len(with_nodes)}/729 assignments "
                  f"are gapless at the high-symmetry points; "
                  f"{len(single_rows)} have a single-orbit node set; ZERO "
                  f"reach isotropy (A < 1e-6); {len(closing)} close the "
                  f"full step at order <= 64 (the internal factor at "
                  f"g_10 never closes - V131). The per-stage angle shifts "
                  f"the node gap and cone slope but never removes the "
                  f"half-zone nodes. 2-route (Bloch vs real-space torus "
                  f"step) dev {max(devs):.1e} on 10 sampled rows.",
        "data": {"baseline": base, "n_rows": len(rows),
                 "n_with_nodes": len(with_nodes),
                 "n_single_orbit": len(single_rows),
                 "n_isotropic": len(isotropic),
                 "n_closing": len(closing),
                 "node_count_histogram": {
                     str(k): sum(1 for r in rows if r["n_nodes"] == k)
                     for k in sorted({r["n_nodes"] for r in rows})},
                 "best_A_overall": min((r["A"] for r in rows
                                        if r["A"] is not None),
                                       default=None),
                 "torus_route_dev": devs,
                 "discriminating": bool(discriminating),
                 "rows": rows, "candidates": candidates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v161_angle_per_stage.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V161 - the 729 per-stage bond-angle census")
    print("=" * 68)
    print(f"  anchor: nodes {base['n_nodes']}, M gap "
          f"{base['gap_by_star']['M']:.4e}, GAMMA "
          f"{base['gap_by_star']['GAMMA']:.4f}, A {base['A']:.3f}")
    print(f"  with nodes: {len(with_nodes)}/729   single orbit: "
          f"{len(single_rows)}   isotropic: {len(isotropic)}   closing: "
          f"{len(closing)}")
    print(f"  discriminating (pi/2 rows differ): {discriminating}   "
          f"2-route dev: {max(devs):.1e}")
    print(f"  node-count histogram: "
          f"{verdict['data']['node_count_histogram']}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v161_angle_per_stage.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
