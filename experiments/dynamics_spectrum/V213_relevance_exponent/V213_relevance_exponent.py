#!/usr/bin/env python3
"""V213 - the relevance exponent of the anisotropic tangent term at the
measured node: does coarse-graining remove the anisotropy or not?

The continuum question of Section 13.1 is not "is the cone isotropic at one
resolution" - it is whether the anisotropic terms in the tangent dispersion
are RELEVANT, MARGINAL or IRRELEVANT at the gapless point.  The earlier reading
of this (V173, and its gauged twin V178) reported the anisotropy of the node
subspace at block scales b = 1, 2, 4, found it unchanged at 0.998, and called
the cubic anisotropy marginal.  That verdict is not available from the data it
describes, for two independent reasons, and V213 measures the slope with both
removed.

(i) V173 measured the anisotropy at the SAME probe radius (dq = 1e-3 in lattice
    units) at every b, so its three numbers differ only in the dimension of the
    downfolded subspace: holding the probe fixed in lattice units while the
    block grows is holding the probe at a momentum that grows by b each step,
    the opposite of the Wilsonian step.  V213 rescales the window pi/(2b) AND
    the probe b r_0 - it follows fixed physical momentum - and it recomputes the
    fixed-probe convention alongside, so the flaw is visible in the data rather
    than asserted.
(ii) V173's value 0.998 is essentially the saturated value of its own measure.
    The gapless set through the node is a line (the GAMMA-M family,
    V132/V203/V212): probed over a sphere of directions, the directions inside
    the nodal set have EXACTLY zero velocity, so (v_max - v_min)/(v_max +
    v_min) = 1 at every radius and carries no exponent.  V213 separates the two
    pieces: the fraction of exactly-nodal directions (the codimension of the
    gapless set, measured at every resolution) and the spread of the NONZERO
    transverse velocities, whose log-log slope is the relevance exponent p:
    p > 0 irrelevant, |p| within continuum.relevance_marginal_band marginal,
    p < 0 relevant.

T1 the momentum ladder.  At the exactly gapless node of the clock in force (the
   M representative, machine-zero gap at the pi/4 root clock, V132/V203,
   re-verified here) the lowest eigenphase E(n-hat, r) of the implemented
   one-step matrix is measured on the V115 direction grid at every radius in
   cl10/config.yaml: continuum.radii.  The leading exponent alpha is the log-log
   slope of the energy (1 for a linear touching); the nodal fraction and the
   transverse anisotropy A(r) come from the same data.

T2 the block ladder.  The Lowdin downfold of H(q) = i log U(q) onto the node
   subspace |theta| < pi/(2b) for b = 1, 2, 4 (continuum.rg_cuts), probe at
   b r_0: the slope of log A against log b is the same exponent read from a
   different ladder.

T3 the invariant-theory prediction, from the measured little group of the node
   (cl10/continuum.py): with d the degree of the first anisotropic invariant of
   the node's little group and alpha the measured leading exponent, the leading
   relative anisotropy is expected to scale as r^(d - alpha) when d > alpha,
   and not at all (p = 0, exactly marginal) when the little group already
   allows an anisotropic term at the leading order itself.

T4 the accuracy controls: the b = 1 downfold against the full one-step matrix at
   the same radius, the window-edge margins of every subspace, and the number of
   usable points in each fit.

The verdict is the decision the continuum route turns on: a positive p outside
the marginal band would mean the implemented O_h module flows to an isotropic
cone and the module can be kept; a marginal or relevant p - or a gapless set
whose codimension does not change with resolution - means no coarse-graining
inside this module produces rotational invariance, and the module itself has to
change.
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
os.environ["CL10_BACKEND"] = "numpy"   # the alt programmes are numpy-native

import numpy as np

from cl10 import altcell as A
from cl10 import config as _config
from cl10 import continuum as K
from cl10 import supercell as SC

CFG = _config.CFG
CELL = A.cell((2, 2, 2))
STAGES = A.stages_baseline()
RADII = [float(r) for r in CFG.continuum.radii]
CUTS = [float(c) for c in CFG.continuum.rg_cuts]     # pi/(2b): b = 1, 2, 4
BLOCKS = [2 ** i for i in range(len(CUTS))]
GRID = int(CFG.continuum.direction_grid)
BAND = float(CFG.continuum.relevance_marginal_band)
CENSUS_DEGREES = int(CFG.continuum.census_degrees)


def bloch(q):
    return A.bloch_cell(CELL, SC.sectors()[2], SC.G10, STAGES, q)


def lowest(q):
    return float(np.min(np.abs(np.angle(np.linalg.eigvals(bloch(q))))))


def power_ladder(xs, ys):
    """log-log slope and prefactor over the usable (positive) points."""
    good = [(x, y) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(good) < 2:
        return None, None, len(good)
    x = np.array([math.log(a) for a, _ in good])
    y = np.array([math.log(b) for _, b in good])
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(math.exp(intercept)), len(good)


def velocity_spread(vs):
    """(v_max - v_min)/(v_max + v_min): the V115/V173 anisotropy measure."""
    return (max(vs) - min(vs)) / (max(vs) + min(vs)) if len(vs) > 1 else 0.0


def nodal_fraction(vs, tol):
    """The fraction of probed directions that stay exactly gapless."""
    return sum(1 for v in vs if v <= tol) / len(vs)


def tangent_field(q0, r, dirs, tol):
    """The tangent dispersion at distance r from the node over the direction
    grid: velocities v = E/r, split into the exactly-nodal directions (the
    ones that stay on the gapless set) and the transverse ones, with both
    anisotropy measures - saturated over all directions, and over the
    transverse velocities only."""
    vs = [lowest(tuple(np.array(q0, float) + r * d)) / r for d in dirs]
    trans = [v for v in vs if v > tol]
    return {"r": r, "v_max": max(vs), "v_mean": float(np.mean(vs)),
            "nodal_directions": len(vs) - len(trans),
            "nodal_fraction": nodal_fraction(vs, tol),
            "v_trans_min": min(trans) if trans else 0.0,
            "v_trans_max": max(trans) if trans else 0.0,
            "A_saturated": velocity_spread(vs),
            "A_transverse": velocity_spread(trans),
            "E_mean": float(np.mean(vs)) * r}


def main():
    t0 = time.time()
    tol = float(CFG.tolerances.standard)
    clock = float(CFG.angle.g10)
    branch = str(CFG.angle.branch)
    node = tuple(float(x) for x in SC.HIGH_SYM["M"][0])
    dirs = A.sphere_dirs(GRID)

    # ------------------------------------------------------- T1: radius ladder
    node_gap = lowest(node)
    ladder = [tangent_field(node, r, dirs, tol) for r in RADII]
    p_radius, cA, nA = power_ladder([row["r"] for row in ladder],
                                    [row["A_transverse"] for row in ladder])
    alpha, cv, nv = power_ladder([row["r"] for row in ladder],
                                 [row["E_mean"] for row in ladder])

    # --------------------------------- T2: block ladder, rescaled probe/window
    r_ref = RADII[0]
    blocks = []
    for b, cut in zip(BLOCKS, CUTS):
        r = r_ref * b
        vs_df, vs_full, vs_fixed = [], [], []
        for d in dirs:
            q = tuple(np.array(node, float) + r * d)
            Heff, ev = A.lowdin_downfold(bloch, q, node, cut=cut)
            vs_df.append(float(np.min(np.abs(ev))) / r)
            vs_full.append(lowest(q) / r)
        if b == BLOCKS[0]:                 # probe b r_0 == r_0: same downfold
            vs_fixed = list(vs_df)
        else:
            for d in dirs:                 # the V173 convention, for contrast
                q = tuple(np.array(node, float) + r_ref * d)
                Heff, ev = A.lowdin_downfold(bloch, q, node, cut=cut)
                vs_fixed.append(float(np.min(np.abs(ev))) / r_ref)
        th = np.abs(np.linalg.eigvalsh(A.effective_hamiltonian(bloch, node)))
        blocks.append({"b": b, "cut": cut, "probe_r": r,
                       "n_low": int(np.sum(th < cut)),
                       "window_margin": float(np.min(np.abs(th - cut))),
                       "A_rescaled": velocity_spread(vs_df),
                       "nodal_fraction_rescaled": nodal_fraction(vs_df, tol),
                       "A_fixed_probe": velocity_spread(vs_fixed),
                       "nodal_fraction_fixed_probe": nodal_fraction(vs_fixed, tol),
                       "A_full_matrix": velocity_spread(vs_full)})
    p_block, _, nB = power_ladder([row["b"] for row in blocks],
                                  [row["A_rescaled"] for row in blocks])
    p_fixed, _, _ = power_ladder([row["b"] for row in blocks],
                                 [row["A_fixed_probe"] for row in blocks])

    # ------------------------------------ T3: the little group of the node
    lg = K.locus_census([("M", node)], dmax=CENSUS_DEGREES)[0]
    d_aniso = lg["first_anisotropic_degree"]
    leading_order = round(alpha) if alpha is not None else None
    p_pred = 0.0 if (d_aniso is None or leading_order is None
                     or d_aniso <= leading_order) else float(d_aniso - leading_order)

    # ------------------------------------------------ T4: accuracy controls
    route_dev = abs(blocks[0]["A_rescaled"] - blocks[0]["A_full_matrix"])
    cut_ok = all(abs(cut - math.pi / (2 * b)) < tol
                 for b, cut in zip(BLOCKS, CUTS))
    nodal_fracs = [row["nodal_fraction"] for row in ladder]

    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    chk("T1 the node is exactly gapless at the clock in force",
        node_gap < tol,
        f"node {list(node)} gap {node_gap:.2e} at the {branch} clock "
        f"g = {clock:.6f}")
    chk("T1 the tangent dispersion is linear in the probe radius",
        alpha is not None and abs(alpha - 1.0) < 0.1,
        f"leading exponent alpha = {alpha:.4f} over {nv} radii")
    chk("T1 the gapless set through the node keeps its codimension at every "
        "resolution",
        all(f > 0 for f in nodal_fracs)
        and max(nodal_fracs) - min(nodal_fracs) < 1e-9,
        f"{ladder[0]['nodal_directions']} of {len(dirs)} directions are exactly "
        f"nodal at every radius: fractions "
        + ", ".join("%.4f" % f for f in nodal_fracs))
    chk("T1 the transverse-velocity anisotropy has enough usable radii",
        nA >= 4, f"{nA} of {len(RADII)} radii usable")
    chk("T2 the windows are the block scales they claim to be",
        cut_ok, "cuts " + ", ".join(f"{c:.4f}" for c in CUTS)
        + " = pi/(2b) for b = " + ", ".join(str(b) for b in BLOCKS))
    chk("T2 the two ladders read the same exponent",
        p_radius is not None and p_block is not None
        and abs(p_radius - p_block) < 0.35,
        f"radius ladder p = {p_radius:.4f} ({nA} points); block ladder with "
        f"rescaled window and probe p = {p_block:.4f} ({nB} points)")
    chk("T2 the V173 fixed-probe convention carries no scaling information",
        p_fixed is not None and abs(p_fixed) < 0.05,
        f"fixed-probe slope {p_fixed:.4f} against rescaled-probe slope "
        f"{p_block:.4f}: the earlier marginality verdict measured the subspace "
        f"dimension, not the renormalisation")
    chk("T4 the b = 1 downfold is exact on the low band",
        route_dev < 1e-6,
        f"|A(downfold, b = 1) - A(full matrix)| = {route_dev:.2e}")
    chk("T4 every block window is clean of edge contamination",
        all(row["window_margin"] > tol for row in blocks),
        "; ".join(f"b {row['b']}: n_low {row['n_low']}, margin "
                  f"{row['window_margin']:.2e}" for row in blocks))
    chk("T3 the measured exponent matches the little-group prediction in sign "
        "class",
        (p_pred == 0.0 and (p_radius is not None and p_radius <= BAND))
        or (p_pred > 0.0 and p_radius is not None and p_radius > BAND),
        f"little group of {list(node)}: order {lg['little_group_order']}, "
        f"quadratic invariants {lg['quadratic_invariants']}, first anisotropic "
        f"degree {d_aniso}, leading exponent alpha {alpha:.3f}: predicted "
        f"p = {p_pred:.3f}, measured {p_radius:.3f}")

    p = p_radius
    irrelevant = p is not None and p > BAND
    relevant = p is not None and p < -BAND
    label = "IRRELEVANT" if irrelevant else ("RELEVANT" if relevant else "MARGINAL")
    flow = ("coarse-graining inside the implemented module drives the cone "
            "isotropic" if irrelevant else
            "no coarse-graining inside the implemented module removes the "
            "anisotropy: the module itself has to change")
    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "relevance exponent of the anisotropic tangent term at the "
                    "exactly gapless node: the transverse-velocity anisotropy "
                    "A(r) over the config radius ladder, and A(b) of the "
                    "Lowdin-downfolded node subspace at b = 1, 2, 4 with both "
                    "the window pi/(2b) and the probe b r_0 rescaled, against "
                    "the fixed-probe convention of V173/V178, with the "
                    "codimension of the gapless set measured at every resolution",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"RELEVANCE EXPONENT at the {branch} clock g = {clock:.6f} at the "
            f"node {list(node)} (gap {node_gap:.1e}). (T1) tangent dispersion "
            f"over {len(dirs)} directions: leading exponent alpha = "
            f"{alpha:.4f} (linear); the gapless set through the node holds its "
            f"codimension at every resolution - {ladder[0]['nodal_directions']} "
            f"of {len(dirs)} directions are exactly nodal (velocity <="
            f"{tol:.0e}) at every radius of the ladder, so the saturated "
            f"anisotropy (v_max - v_min)/(v_max + v_min) is pinned at "
            f"{ladder[0]['A_saturated']:.4f} and carries no exponent, which is "
            f"what V173's 0.998 was. Transverse anisotropy: "
            + ", ".join(f"{row['A_transverse']:.4f} at r = {row['r']:.0e}"
                        for row in ladder)
            + f"; log-log slope p = {p_radius:.4f} over {nA} radii. "
              f"(T2) block ladder with window pi/(2b) and probe b r_0: A = "
              + ", ".join(f"{row['A_rescaled']:.4f} (b {row['b']})"
                          for row in blocks)
            + f", slope {p_block:.4f} - the same exponent from an independent "
              f"ladder (agreement {abs(p_radius - p_block):.3f}). The "
              f"V173/V178 fixed-probe convention gives "
              + ", ".join(f"{row['A_fixed_probe']:.4f} (b {row['b']})"
                          for row in blocks)
            + f", slope {p_fixed:.4f}: it measures the downfolded subspace, "
              f"not the flow. (T3) the node's little group (order "
              f"{lg['little_group_order']}, {lg['quadratic_invariants']} "
              f"quadratic invariants) already allows an anisotropic invariant "
              f"at degree {d_aniso} <= the leading order {leading_order}, predicting p = {p_pred:.3f} "
              f"against the measured {p_radius:.3f}. (T4) the b = 1 downfold "
              f"reproduces the full matrix to {route_dev:.1e}; window margins "
              + ", ".join(f"{row['window_margin']:.1e}" for row in blocks)
            + "; subspace sizes "
            + ", ".join(str(row["n_low"]) for row in blocks)
            + f". VERDICT: the anisotropic term is {label} at the measured node "
              f"(p = {p if p is not None else float('nan'):.3f}, marginal band "
              f"+/-{BAND:.0e}) and the gapless set never loses its "
              f"codimension, so {flow}."),
        "data": {"clock_branch": branch, "g10": clock, "node": list(node),
                 "node_gap": node_gap, "directions": len(dirs),
                 "radius_ladder": ladder, "p_radius": p_radius,
                 "p_radius_prefactor": cA, "n_radius_points": nA,
                 "leading_order": leading_order,
                 "leading_exponent": alpha, "leading_coefficient": cv,
                 "blocks": blocks, "p_block": p_block,
                 "p_fixed_probe": p_fixed, "n_block_points": nB,
                 "blocks_used": BLOCKS, "little_group": lg,
                 "predicted_p": p_pred, "route_dev": route_dev,
                 "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v213_relevance_exponent.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V213 - the relevance exponent of the node anisotropy")
    print("=" * 68)
    for row in ladder:
        print(f"  r = {row['r']:.1e}: v_max {row['v_max']:.4f}  nodal "
              f"{row['nodal_directions']:3d}/{len(dirs)}  A_sat "
              f"{row['A_saturated']:.4f}  A_trans {row['A_transverse']:.4f}")
    print(f"  radius ladder: alpha = {alpha:.4f}   p = {p_radius:.4f} "
          f"({nA} points)")
    for row in blocks:
        print(f"  b = {row['b']} (cut {row['cut']:.3f}, r {row['probe_r']:.1e}, "
              f"n_low {row['n_low']}): A rescaled {row['A_rescaled']:.4f} | "
              f"fixed probe {row['A_fixed_probe']:.4f} | full "
              f"{row['A_full_matrix']:.4f}")
    print(f"  block ladder p = {p_block:.4f}   fixed-probe slope = {p_fixed:.4f}")
    print(f"  little group order {lg['little_group_order']}, first anisotropic "
          f"degree {d_aniso}, alpha {alpha:.3f}, predicted p = {p_pred:.3f}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
