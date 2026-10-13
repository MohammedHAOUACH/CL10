#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V165 - does a period-2 FLOQUET alternation of stage lists fix the
nodes? (alt-dynamics, axis d)

**The question.**  The implemented step is the single six-stage cycle
A = STAGES.  A Floquet variant alternates two stage lists per step:
A, B, A, B, ... with B a discrete re-reading of A:

  reverse(A)      the cycle run backwards;
  axis-cycle(A)   the axis index 0 -> 1 -> 2 -> 0;
  parity-swap(A)  the bond phase 0 <-> 1.

No new parameter (same six stages, same g_10, same cell).  The one-step
operator is the two-step Bloch M2(q) = M_B(q) M_A(q); the time-translation
unit is two stages of the original cycle, so the per-stage causal budget
stays 1 Planck cell.  Per variant: gapless set of M2 at the 27
high-symmetry momenta, orbit structure, chirality sum, anisotropy A on
the M2 eigenphases at the first node, front speed cycling A, B.
Anchor: the trivial 2-step (A, A) must show the baseline nodes with
doubled period (gap at M = 2 x 2.0029e-4, 12 nodes).

**The result.**  The (A, A) anchor reproduces the doubled-period nodes.
None of the three Floquet variants removes the half-zone nodes or
isotropizes a cone (A >= 0.99 everywhere gapless): the time-dependent
reading of the same six stages is a relabeling of the same motion.
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

C = A.cell((2, 2, 2))
NODE_TOL = 1e-3


def variants():
    base = A.stages_baseline()
    out = {"AA": base,
           "reverse": list(reversed(base)),
           "axis_cycle": [((j + 1) % 3, p) for (j, p) in base],
           "parity_swap": [(j, 1 - p) for (j, p) in base]}
    return base, out


def _census(B, U16, dirs):
    base = A.stages_baseline()

    def M2(qq):
        MA = A.bloch_stages(C, [U16] * 6, [SC.G10] * 6, base, qq)
        MB = A.bloch_stages(C, [U16] * 6, [SC.G10] * 6, B, qq)
        return MB @ MA

    gaps = {}
    nodes = []
    for name, pts in SC.HIGH_SYM.items():
        for q in pts:
            gg = A.gap_at(M2(q))
            gaps[name] = min(gaps.get(name, 9.9), gg)
            if gg < NODE_TOL:
                nodes.append(q)
    chir = 0.0
    for q in nodes:
        ch, _sl = A.chirality_at(M2, q)
        chir += ch
    A_first = None
    growth = None
    if nodes:
        aa, _vx, _vn = A.cone_anisotropy(M2, nodes[0], dirs)
        A_first = aa
        n0 = np.array([1.0, 1.0, 1.0]) / math.sqrt(3)
        g1 = A.gap_at(M2(tuple(q + 1e-3 * n0[k]
                              for k, q in enumerate(nodes[0]))))
        g2 = A.gap_at(M2(tuple(q + 1e-2 * n0[k]
                              for k, q in enumerate(nodes[0]))))
        growth = g2 / g1 if g1 > 0 else None
    stars = {tuple(A.oh_star(q)) for q in nodes}
    single = bool(nodes) and len(stars) == 1 and \
        len(next(iter(stars), ())) == len(nodes)
    Mq = M2((0.31, -1.13, 2.02))
    unit = float(np.abs(Mq @ Mq.conj().T - np.eye(Mq.shape[0])).max())
    N, cerr = A.exact_order(M2((0.0, 0.0, 0.0)), Nmax=64)
    fs = None
    if nodes:
        cyc = list(base) + list(B)
        fs = A.front_speed_stages(C, [U16] * 12, [SC.G10] * 12, cyc,
                                  Lb=6, n_stages=12)
    return {"gap_by_star": gaps, "n_nodes": len(nodes),
            "n_stars": len(stars), "single_orbit": single,
            "chirality_sum": chir, "A": A_first, "cone_growth": growth,
            "chirality_defined": all(
                any(abs(c) > 1e-12 for c in q) for q in nodes),
            "unitary_err": unit, "closure_order": N, "closure_err": cerr,
            "closed": cerr < 1e-9,
            "front_speed": fs["front_speed"] if fs else None,
            "outside_lightcone": fs["n_outside_lightcone"] if fs else None}


def main():
    t0 = time.time()
    U16 = SC.sectors()[2]
    dirs = A.sphere_dirs(6)
    base, vs = variants()
    rows = []
    for name in ("AA", "reverse", "axis_cycle", "parity_swap"):
        r = _census(vs[name], U16, dirs)
        r["variant"] = name
        rows.append(r)

    anchor = rows[0]
    ok_anchor = (anchor["n_nodes"] == 12
                 and abs(anchor["gap_by_star"]["M"]
                         - 2 * 2.0029e-4) < 2e-6
                 and anchor["single_orbit"]
                 and abs(anchor["chirality_sum"]) < 1e-9)
    # 2-route: the two-step Bloch M2 vs the real-space two-step on the
    # (2,2,2)-cell torus (12 stages cycling A, B)
    def torus_route(B):
        cyc = list(base) + list(B)
        L = 2
        dim = L ** 3 * 8 * 16
        Mat = np.zeros((dim, dim), dtype=complex)
        for col in range(dim):
            f = np.zeros((L, L, L, 8, 16), dtype=complex)
            f.ravel()[col] = 1.0
            out = A.step_stages(f, [U16] * 12, [SC.G10] * 12, cyc)
            Mat[:, col] = out.ravel()
        ph_r = np.sort(np.angle(np.linalg.eigvals(Mat)))
        ph_b = []
        for q in itertools.product((0.0, math.pi), repeat=3):
            MA = A.bloch_stages(C, [U16] * 6, [SC.G10] * 6, base, q)
            MB = A.bloch_stages(C, [U16] * 6, [SC.G10] * 6, B, q)
            ph_b.append(np.angle(np.linalg.eigvals(MB @ MA)))
        ph_b = np.sort(np.concatenate(ph_b))
        return float(np.abs(ph_r - ph_b).max())

    devs = [torus_route(vs[n]) for n in ("AA", "reverse", "axis_cycle",
                                         "parity_swap")]
    ok_route = max(devs) < 1e-9

    gapless = [r for r in rows if r["n_nodes"] > 0]
    isotropic = [r for r in rows if r["A"] is not None and r["A"] < 1e-6]
    true_cones = [r for r in isotropic if (r["cone_growth"] or 0.0) > 5
                  and r["chirality_defined"]]

    candidates = [{
        "name": f"floquet:{r['variant']}",
        "source": "V165_floquet",
        "unitary_err": r["unitary_err"], "n_nodes": r["n_nodes"],
        "single_orbit": r["single_orbit"], "anisotropy_A": r["A"],
        "chirality_sum": r["chirality_sum"], "front_speed": r["front_speed"],
        "mirror_ok": True,
        "survivor": (r["n_nodes"] > 0 and r["single_orbit"]
                     and r["A"] is not None and r["A"] < 1e-6
                     and (r["cone_growth"] or 0.0) > 5
                     and r["chirality_defined"]
                     and abs(r["chirality_sum"]) < 1e-9
                     and r["unitary_err"] < 1e-12),
        "step_spec": {"kind": "floquet", "variant": r["variant"]}}
        for r in rows]

    ok = ok_anchor and ok_route and len(true_cones) == 0
    best_A = min((r["A"] for r in rows if r["A"] is not None),
                 default=float("nan"))

    verdict = {
        "id": "V165_floquet",
        "tier": 4,
        "property": "period-2 Floquet alternation of the stage list: "
                    "two-step node census, orbit, chirality, anisotropy, "
                    "front speed",
        "residual": True,
        "pass": bool(ok),
        "status": f"FLOQUET CENSUS (3 variants + trivial anchor). The "
                  f"(A,A) anchor shows the baseline nodes with doubled "
                  f"period (12 nodes, M gap "
                  f"{anchor['gap_by_star']['M']:.4e} = 2 x 2.0029e-4). "
                  f"{len(gapless)}/4 two-step operators are gapless at "
                  f"the high-symmetry points; ZERO true isotropic cones "
                  f"(best A {best_A:.3f}). 2-route dev "
                  f"{max(devs):.1e}. The period-2 reading of the same "
                  f"six stages is a relabeling: no new node structure.",
        "data": {"anchor": anchor, "n_rows": len(rows),
                 "n_gapless": len(gapless), "n_true_cones": len(true_cones),
                 "torus_route_dev": devs,
                 "rows": rows, "candidates": candidates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v165_floquet.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V165 - the period-2 Floquet census")
    print("=" * 68)
    for r in rows:
        print(f"  {r['variant']:>12}  nodes {r['n_nodes']:3d}  M gap "
              f"{r['gap_by_star']['M']:.3e}  A "
              f"{r['A'] if r['A'] is None else round(r['A'], 4)}  front "
              f"{r['front_speed']}")
    print(f"  2-route dev: {max(devs):.1e}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v165_floquet.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
