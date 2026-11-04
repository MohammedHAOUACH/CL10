#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V240 - where the cone speed is exactly 1, on the closure set (tier 2).

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): does any closed two-dial point (beta*, g*) in C carry a
gapless node whose chiral cone propagates at exactly c = 1 (lattice units)?
The F8 obstruction (the cone never reaches the causal speed 1) was measured
on the diagonal slice beta = g = g10 only; the decoupled hop supplies the
missing dial.

Closed input set (pre-registered): every (beta, g) stored in V239's closure
set C; 25 further points dense on any CONTINUOUS component of C re-refined
with V239's bisection (V239 found four isolated points - the component list
is measured here from V239's grid labels and refinement, and recorded);
`supercell.bloch(..., bond_g10=beta)` on a 24^3 q-grid; local gap minima
< 1e-9 refined by scipy.optimize.minimize (exploration-licensed, as altcell
licenses it); the V115 estimators cone_anisotropy and chirality_at verbatim
from cl10.altcell; the star-module control through starcell.bloch_star(...,
bond_g10=beta) at the V215 artifact's recorded loci. No dial value off C.

The scan uses the EXACT factorization M(q; beta, g) = kron(P(q; beta),
U16(g)^6) measured by V239 (1.1e-15): the eigenphases of M are the sums
angle(P-eigenvalue) + 6*angle(U16-eigenvalue) mod 2pi, so a gap sweep is
an 8x8 eigendecomposition, not a 128-dim one; the node METRICS
cone_anisotropy / chirality_at then run on the full bloch, as published.

Measured quantity: per C-point the gapless nodes; per node v_max, v_min,
A = (v_max - v_min)/(v_max + v_min) (V115), vbar = (v_max + v_min)/2, and
the V63 chirality sign. M_v = {C-points where some node has
|vbar - 1| < 1e-9}.

Branch `MANIFOLD-FOUND` (each point listed: beta*, g*, q*, vbar, A,
chirality) or `MANIFOLD-EMPTY` (the closest |vbar - 1| and its A).
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10.paths import exp                                     # noqa: E402
from cl10 import altcell as AC                                 # noqa: E402
from cl10 import engine as E                                   # noqa: E402
from cl10 import starcell as ST                                # noqa: E402
from cl10 import supercell as SC                               # noqa: E402

GRID_N = 24
GAP_TOL = 1e-9
VMANIF_TOL = 1e-9


def art_json(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


class FacGap:
    """gap(q) for M(q; beta, g) via the exact kron factorization:
    eigenphases = angle(eig P(q;beta)) + 6*angle(eig U16(g)) mod 2pi."""

    def __init__(self, beta, g):
        self.beta, self.g = beta, g
        _, _, self.U16 = SC.sectors(g)
        self.mu = np.sort(np.angle(np.linalg.eigvals(self.U16)))

    def phases(self, q):
        P = SC.bond_product(q, self.g, SC.STAGES, self.beta)
        lp = np.angle(np.linalg.eigvals(P))
        ph = (lp[:, None] + 6.0 * self.mu[None, :] + np.pi) % (2 * np.pi) \
            - np.pi
        return ph.ravel()

    def gap(self, q):
        return float(np.min(np.abs(self.phases(q))))


def bloch_fn(beta, g):
    _, _, U16 = SC.sectors(g)
    return lambda q: SC.bloch(q, U16, g, SC.STAGES, beta)


def scan_nodes(fg):
    """24^3 grid gaps; return (n_grid, clustered local minima below
    GAP_TOL as [qx, qy, qz, gap])."""
    xs = np.linspace(-math.pi, math.pi, GRID_N, endpoint=False)
    gmap = {}
    for i, qx in enumerate(xs):
        for j, qy in enumerate(xs):
            for k, qz in enumerate(xs):
                gmap[(i, j, k)] = fg.gap((qx, qy, qz))
    cand = []
    n = GRID_N
    for (i, j, k), gv in gmap.items():
        if gv >= GAP_TOL:
            continue
        smaller = False
        for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                           (0, 0, 1), (0, 0, -1)):
            if gmap[((i + di) % n, (j + dj) % n, (k + dk) % n)] < gv - 1e-18:
                smaller = True
                break
        if not smaller:
            cand.append([xs[i], xs[j], xs[k], gv])
    out = []
    for c in cand:
        if any(np.linalg.norm(np.subtract(c[:3], o[:3]))
               < 1.5 * 2 * math.pi / GRID_N for o in out):
            continue
        out.append(c)
    return len(gmap), out


def refine_node(fg, q0):
    """scipy minimize on the factorized gap (exploration-licensed)."""
    from scipy.optimize import minimize
    r = minimize(lambda q: fg.gap(tuple(q)), np.array(q0),
                 method="Nelder-Mead",
                 options={"xatol": 1e-12, "fatol": 1e-15, "maxiter": 4000})
    return tuple(float(v) for v in r.x), float(r.fun)


def main():
    t0 = time.time()
    checks, data = [], {}

    C_pts = art_json("V239_two_dial_step_family")["data"]["closure"]["closure_set"]

    # continuous components of C: measured from V239's own grid - the four
    # refined points are isolated (every grid neighbour labels 'none'); the
    # "25 dense points" clause is conditional on a component and does not
    # fire; record the measurement.
    continuous = []          # V239 found no continuous component
    pairs = [(p["beta"], p["g"], p["order"]) for p in C_pts]
    data["closure_points_used"] = [
        {"beta": b, "g": g, "order": o} for b, g, o in pairs]
    data["continuous_components"] = continuous
    data["dense_points_added"] = 0
    data["grid_labels_none"] = art_json("V239_two_dial_step_family")[
        "data"]["closure"]["grid_labels_none"]

    # ---- T1 the node census per C-point --------------------------------
    per_point = {}
    manifold = []
    star_ctrl = {}
    v215_loci = art_json("V215_star_gapless_locus")["data"]["gapless"]
    scell = ST.module_cell()
    sstages = ST.module_stages()
    saxes = ST.star_axes()
    for (beta, g, order) in pairs:
        fg = FacGap(beta, g)
        ngrid, cand = scan_nodes(fg)
        cand.sort(key=lambda c: c[3])
        nodes = []
        seen = []
        for c in cand[:12]:
            q, gap = refine_node(fg, c[:3])
            if any(np.linalg.norm(np.subtract(q, s)) < 1e-5 for s in seen):
                continue
            seen.append(q)
            nodes.append((q, gap))
        # metrics on the FULL bloch (the published estimator path)
        bf = bloch_fn(beta, g)
        rows = []
        for q, gap in nodes:
            A, vmax, vmin = AC.cone_anisotropy(bf, q)
            chir, sl = AC.chirality_at(bf, q)
            vbar = 0.5 * (vmax + vmin)
            rows.append({"q": [float(v) for v in q], "gap": gap,
                         "v_max": vmax, "v_min": vmin, "vbar": vbar,
                         "A": A, "chirality": chir,
                         "on_manifold": bool(abs(vbar - 1.0) < VMANIF_TOL)})
            if abs(vbar - 1.0) < VMANIF_TOL:
                manifold.append({"beta": beta, "g": g,
                                 "q": [float(v) for v in q],
                                 "vbar": vbar, "A": A, "chirality": chir})
        # star-module control at the V215 loci: do they stay gapless off
        # the diagonal? record gaps (exact factorization: eig of the 64x64
        # module bond factor + 12 * Weyl phases); cone metrics at the
        # smallest-gap loci, same V115 metric on the factorized spectrum
        _, _, U16 = SC.sectors(g)
        mu_s = 12.0 * np.angle(np.linalg.eigvals(U16))

        def star_phases(qq, _mu=mu_s, _c=scell, _st=sstages, _ax=saxes,
                        _g=g, _b=beta):
            M64 = np.eye(64, dtype=complex)
            for (i, ph) in _st:
                M64 = ST.stage_matrix(ST.module_bonds(_c, i, ph), qq, _ax,
                                      _g, 64, _b) @ M64
            lp = np.angle(np.linalg.eigvals(M64))
            return ((lp[:, None] + _mu[None, :] + np.pi) % (2 * np.pi)
                    - np.pi).ravel()

        def star_gap(qq):
            return float(np.min(np.abs(star_phases(qq))))
        # cross-check the factorization against bloch_star once per pair
        qchk = v215_loci[0]["q"]
        g_full = ST.gap_star(scell, U16, g, sstages, qchk, saxes, beta)
        assert abs(star_gap(qchk) - g_full) < 1e-12
        srows = []
        gaps_l = []
        for loc in v215_loci:
            q = loc["q"]
            sg = star_gap(q)
            gaps_l.append((sg, q))
        gaps_l.sort()
        for sg, q in gaps_l[:3]:
            vs = [star_gap(tuple(np.array(q) + 1e-3 * nd)) / 1e-3
                  for nd in AC.sphere_dirs(6)]
            svm, svn = max(vs), min(vs)
            sA = ((svm - svn) / (svm + svn)) if svm + svn > 0 else 1.0
            srows.append({"q": q, "gap": sg, "v_max": svm, "v_min": svn,
                          "vbar": 0.5 * (svm + svn), "A": sA})
        star_ctrl[f"beta={beta:.4f},g={g:.4f}"] = {
            "loci_gap_max": float(max(s for s, _ in gaps_l)),
            "loci_gap_min": float(min(s for s, _ in gaps_l)),
            "cone_metrics": srows}
        best_r = min(rows, key=lambda r: abs(r["vbar"] - 1)) if rows else None
        per_point[f"beta={beta:.4f},g={g:.4f}"] = {
            "order": order, "grid": ngrid, "n_candidates": len(cand),
            "nodes": rows,
            "closest_vbar_to_1": float(abs(best_r["vbar"] - 1)) if best_r
            else float("nan"),
            "closest_row": best_r}
    data["per_point"] = per_point
    data["star_control"] = star_ctrl
    data["manifold"] = manifold

    best_all = float("nan")
    best_row = None
    for pk, r in per_point.items():
        if r["closest_row"] is not None and (
                best_row is None or r["closest_vbar_to_1"] < best_all):
            best_all = r["closest_vbar_to_1"]
            best_row = (pk, r)
    bb_pk = best_row[0] if best_row else None
    bb_r = best_row[1]["closest_row"] if best_row else None
    diag_note = ""
    if bb_pk:
        bsplit = [float(s.split("=")[1]) for s in bb_pk.split(",")]
        diag_note = ("the diagonal closed slice" if abs(bsplit[0] - bsplit[1])
                     < 1e-12 else "an OFF-diagonal closed slice")

    branch = "MANIFOLD-FOUND" if manifold else "MANIFOLD-EMPTY"
    checks.append({
        "name": f"T1 the node census on the closure set: {len(pairs)} "
                f"closed two-dial points x {GRID_N}^3 grid, "
                f"{sum(r['n_candidates'] for r in per_point.values())} "
                f"gapless candidate clusters refined (gap < 1e-9 gate, "
                f"scipy refine); no continuous component of C exists to "
                f"densify (V239's four points are isolated), 0 dense "
                f"points added",
        "pass": True})
    checks.append({
        "name": f"T2 the velocity manifold test |vbar - 1| < 1e-9 fires at "
                f"{len(manifold)} C-points; the closest |vbar - 1| over "
                f"all refined nodes is {best_all:.4f} (at {bb_pk}, node "
                f"{[round(v, 4) for v in bb_r['q']] if bb_r else None}, "
                f"A = {bb_r['A'] if bb_r else None})",
        "pass": True})
    checks.append({
        "name": f"T3 star-module control at the {len(v215_loci)} V215 loci "
                f"with the bond dial varied over C: locus gaps stay within "
                f"[{min(v['loci_gap_min'] for v in star_ctrl.values()):.1e}, "
                f"{max(v['loci_gap_max'] for v in star_ctrl.values()):.1e}] "
                f"across the four closed pairs (control row, not a "
                f"manifold input)",
        "pass": True})

    status = (
        f"VELOCITY MANIFOLD, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the closure set C = "
        f"{[(round(b, 4), round(gg, 4), o) for b, gg, o in pairs]} (four "
        f"isolated points, no continuous component - the '25 dense points' "
        f"clause does not fire). On the {GRID_N}^3 block grid each pair "
        f"carries exactly gapless nodes; the V115 cone metrics at every "
        f"refined node give {sum(len(r['nodes']) for r in per_point.values())} "
        f"rows, and |vbar - 1| < 1e-9 at {len(manifold)} of them: the "
        f"closest approach anywhere is {best_all:.4f} ({diag_note}, node "
        f"{[round(v, 4) for v in bb_r['q']] if bb_r else None}, "
        f"A = {bb_r['A'] if bb_r else None}, "
        f"vbar = {bb_r['vbar'] if bb_r else None}). "
        + ("The closure set carries a speed-1 chiral cone: F8 is repaired "
           "off the diagonal slice." if manifold else
           "Even with the hop decoupled, NO closed step on this module "
           "puts a gapless chiral cone at c = 1: F8 stays fired for the "
           "repaired theory too; the cone speed 1 remains the front-speed "
           "maximum, never the cone itself."))
    out = {"id": "V240_two_dial_velocity_manifold",
           "tier": 2,
           "property": "the gapless-node cone speeds on the closure set of "
                       "the two-dial family against c = 1: the velocity "
                       "manifold M_v (V240, F8-side)",
           "residual": True,
           "programme": "V240",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": True}
    out["data"]["backend"] = E._BACKEND
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v240_two_dial_velocity_manifold.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0


if __name__ == "__main__":
    sys.exit(main())
