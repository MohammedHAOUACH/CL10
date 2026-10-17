#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V163 - does a q-INDEPENDENT mass stage gap the nodes? (alt-dynamics, axis c)

**The question.**  The published node-gap mechanism (V65) is the
q-DEPENDENT Wilson factor.  The alternative: one extra internal factor per
step, ``U_step = M_ab(g) * U^6`` with ``M_ab(g) = cos g I - i sin g T_ab``
a single-plane rotation at an allowed angle - q-independent, so it keeps
strict causality.  Census: all 45 planes x g in {pi/2, pi/4} = 90 rows.
Per row: node gaps at the 27 high-symmetry momenta (does a single plane
rotation gap the 12 M nodes WITHOUT the q-dependent Wilson form?), the
remaining gapless set, anisotropy A at the first node, mirror leakage
(the 32-sector commutator with w), closure of the full step at q = 0.

**The result.**  81/90 rows gap EVERY high-symmetry point (the
mechanism works - but it gaps everything, so the fermions are all
massive).  9 rows keep points under the 1e-3 node tolerance; the two
with A < 1e-6 are massive zone-centre points (cone growth ~0.9, no
linear dispersion, chirality undefined), not Dirac cones.  Mirror
leakage is exactly 0 (even-grade planes) and no row closes the full
step at order <= 64.  A q-independent mass stage is a mass generator,
not a node selector: it cannot be the mechanism the open problems need.
"""

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

import itertools
import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
NODE_TOL = 1e-3
ANGLES = (("pi2", math.pi / 2), ("pi4", math.pi / 4))


def main():
    t0 = time.time()
    B, U32, U16 = A.sectors_sig(0)
    planes32 = A.plane_generators_sig(0)
    planes16 = [B.conj().T @ T @ B for T in planes32]
    w = A.pseudoscalar_sig(0)
    dirs = A.sphere_dirs(6)
    I16 = np.eye(16, dtype=complex)
    pairs = [(a, b) for a in range(10) for b in range(a + 1, 10)]
    assert len(pairs) == 45

    rows = []
    for (a, b), (gname, g) in [(p, ag) for p in pairs for ag in ANGLES]:
        Mab16 = math.cos(g) * I16 - 1j * math.sin(g) * planes16[
            pairs.index((a, b))]
        Mab32 = math.cos(g) * np.eye(32, dtype=complex) - 1j * math.sin(
            g) * planes32[pairs.index((a, b))]
        leak = float(np.abs(Mab32 @ w - w @ Mab32).max())

        def bloch_fn(qq, Mab=Mab16):
            M = A.bloch_cell(C, U16, SC.G10, ST, qq)
            return np.kron(np.eye(8, dtype=complex), Mab) @ M

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
        growth = None
        if nodes:
            aa, _vx, _vn = A.cone_anisotropy(bloch_fn, nodes[0], dirs)
            A_first = aa
            # a true cone: the gap at dq = 1e-2 grows ~10x over dq = 1e-3;
            # a massive point saturates (ratio ~1)
            n0 = np.array([1.0, 1.0, 1.0]) / math.sqrt(3)
            g1 = A.gap_at(bloch_fn(tuple(q + 1e-3 * n0[k] for k, q
                                        in enumerate(nodes[0]))))
            g2 = A.gap_at(bloch_fn(tuple(q + 1e-2 * n0[k] for k, q
                                         in enumerate(nodes[0]))))
            growth = g2 / g1 if g1 > 0 else None
        stars = {tuple(A.oh_star(q)) for q in nodes}
        single = bool(nodes) and len(stars) == 1 and \
            len(next(iter(stars), ())) == len(nodes)
        Mq = bloch_fn((0.31, -1.13, 2.02))
        unit = float(np.abs(Mq @ Mq.conj().T - np.eye(Mq.shape[0])).max())
        N, cerr = A.exact_order(bloch_fn((0.0, 0.0, 0.0)), Nmax=64)
        fs = None
        if nodes:
            Us = [U16] * 6 + [Mab16]
            th = [SC.G10] * 6 + [0.0]
            stx = list(ST) + [[]]
            fs = A.front_speed_stages(C, Us, th, stx, Lb=6, n_stages=12)
        rows.append({"plane": (a, b), "angle": gname,
                     "gap_by_star": gaps, "n_nodes": len(nodes),
                     "n_stars": len(stars), "single_orbit": single,
                     "chirality_sum": chir, "A": A_first,
                     "cone_growth": growth,
                     "chirality_defined": all(
                         any(abs(c) > 1e-12 for c in q) for q in nodes),
                     "mirror_leakage": leak, "unitary_err": unit,
                     "closure_order": N, "closure_err": cerr,
                     "closed": cerr < 1e-9,
                     "front_speed": fs["front_speed"] if fs else None,
                     "outside_lightcone":
                         fs["n_outside_lightcone"] if fs else None})

    gapped_all = [r for r in rows if r["n_nodes"] == 0]
    with_nodes = [r for r in rows if r["n_nodes"] > 0]
    isotropic = [r for r in rows if r["A"] is not None and r["A"] < 1e-6]
    true_cones = [r for r in isotropic if (r["cone_growth"] or 0.0) > 5
                  and r["chirality_defined"]]
    closing = [r for r in rows if r["closed"]]
    ok_leak = max(r["mirror_leakage"] for r in rows) < 1e-12

    # 2-route on 3 sampled rows: the real-space torus step (6 bond
    # stages + internal layers + the mass layer, built column-by-column
    # from step_stages / internal_layer) must reproduce the union of the
    # Bloch eigenphase multisets over the 8 torus momenta q in {0, pi}^3
    def torus_route(a, b, g):
        Mab16 = math.cos(g) * I16 - 1j * math.sin(g) * planes16[
            pairs.index((a, b))]
        L = 2
        dim = L ** 3 * 8 * 16
        Mat = np.zeros((dim, dim), dtype=complex)
        for col in range(dim):
            f = np.zeros((L, L, L, 8, 16), dtype=complex)
            f.ravel()[col] = 1.0
            f = A.step_stages(f, [U16] * 6, [SC.G10] * 6, ST)
            f = SC.internal_layer(f, Mab16)
            Mat[:, col] = f.ravel()
        ph_r = np.sort(np.angle(np.linalg.eigvals(Mat)))
        ph_b = []
        for q in itertools.product((0.0, math.pi), repeat=3):
            Mab = math.cos(g) * I16 - 1j * math.sin(g) * planes16[
                pairs.index((a, b))]
            Mq = np.kron(np.eye(8, dtype=complex), Mab) \
                @ A.bloch_cell(C, U16, SC.G10, ST, q)
            ph_b.append(np.angle(np.linalg.eigvals(Mq)))
        ph_b = np.sort(np.concatenate(ph_b))
        return float(np.abs(ph_r - ph_b).max())

    sample = [(0, 1, math.pi / 2), (4, 9, math.pi / 4),
              (0, 6, math.pi / 4)]
    devs = [torus_route(a, b, g) for (a, b, g) in sample]
    ok_route = max(devs) < 1e-9

    candidates = [{
        "name": f"mass:({r['plane'][0]},{r['plane'][1]})@{r['angle']}",
        "source": "V163_mass_stage",
        "unitary_err": r["unitary_err"], "n_nodes": r["n_nodes"],
        "single_orbit": r["single_orbit"], "anisotropy_A": r["A"],
        "chirality_sum": r["chirality_sum"], "front_speed": r["front_speed"],
        "mirror_ok": r["mirror_leakage"] < 1e-12,
        "survivor": (r["n_nodes"] > 0 and r["single_orbit"]
                     and r["A"] is not None and r["A"] < 1e-6
                     and (r["cone_growth"] or 0.0) > 5
                     and r["chirality_defined"]
                     and abs(r["chirality_sum"]) < 1e-9
                     and r["unitary_err"] < 1e-12),
        "step_spec": {"kind": "mass", "plane": list(r["plane"]),
                      "g": dict(ANGLES)[r["angle"]]}}
        for r in rows]

    ok = (len(rows) == 90 and ok_leak and ok_route
          and len(true_cones) == 0)

    verdict = {
        "id": "V163_mass_stage",
        "tier": 4,
        "property": "one extra q-independent internal plane rotation per "
                    "step: 45 planes x {pi/2, pi/4}, node gaps, residual "
                    "gapless set, mirror leakage, closure",
        "residual": True,
        "pass": bool(ok),
        "status": f"MASS-STAGE CENSUS (90 rows). {len(gapped_all)}/90 rows "
                  f"gap EVERY high-symmetry point: a single-plane "
                  f"q-independent rotation gaps the 12 M nodes but gaps "
                  f"the whole spectrum. {len(with_nodes)} rows keep "
                  f"points under the 1e-3 node tolerance; the two with "
                  f"A < 1e-6 are MASSIVE zone-centre points (cone growth "
                  f"{[round(r['cone_growth'], 2) for r in isotropic]}, "
                  f"chirality undefined), not Dirac cones: ZERO rows "
                  f"supply an isotropic single-orbit gapless set "
                  f"({len(true_cones)} true cones). Mirror leakage "
                  f"exactly {max(r['mirror_leakage'] for r in rows):.1e} "
                  f"(even-grade planes); {len(closing)} rows close at "
                  f"order <= 64. 2-route dev {max(devs):.1e}. A "
                  f"q-independent mass stage is a mass generator, not a "
                  f"node selector.",
        "data": {"n_rows": len(rows), "n_gapped_all": len(gapped_all),
                 "n_with_nodes": len(with_nodes),
                 "n_isotropic": len(isotropic),
                 "n_true_cones": len(true_cones),
                 "n_closing": len(closing),
                 "max_mirror_leakage": max(r["mirror_leakage"]
                                           for r in rows),
                 "torus_route_dev": devs,
                 "min_M_gap": min(r["gap_by_star"]["M"] for r in rows),
                 "rows": rows, "candidates": candidates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v163_mass_stage.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V163 - the q-independent mass-stage census (90 rows)")
    print("=" * 68)
    print(f"  gapped everywhere: {len(gapped_all)}/90   with nodes: "
          f"{len(with_nodes)}   isotropic: {len(isotropic)}   closing: "
          f"{len(closing)}")
    print(f"  min M gap over rows: "
          f"{min(r['gap_by_star']['M'] for r in rows):.4e}")
    print(f"  max mirror leakage: "
          f"{max(r['mirror_leakage'] for r in rows):.1e}   2-route dev: "
          f"{max(devs):.1e}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v163_mass_stage.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
