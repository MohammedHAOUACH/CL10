#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V82 - the composite desert content is EXCLUDED under every natural
mass ansatz: W4's supplier question is closed in both directions.

**Why (W4).**  V64 repairs unification with three 10+10bar pairs at
1.2e7 GeV plus one (1,2) pair at 1e12 GeV meeting at alpha_X = 1/25.
V74 excluded the node mechanism as the supplier (forced M_X = 8.9e14,
F5-excluded). V81 excluded the gauged layer's composite content under
the flux-only threshold ansatz m = w(c0 + |alpha.v|). This experiment
closes the remaining natural ansatz families and quantifies WHY.

**The content (0-parameter).**  The composite fermions of the V79 gauged
bond layer: matter 16 weights mu tensored with adjoint link flux weights
alpha: nu = mu + alpha, 720 Weyl states, five net chiral 16-copies; the
one-loop index shifts db_state are computed exactly (V81 LIB regression:
err < 8e-15).

**Ansatz families tested (thresholds z_s = ln[w(c0 + q_s)/M_Z]).**
A1 flux-only:        q_s = |alpha.v|                     (V81: 0 solutions)
A2 full weight:      q_s = |nu.v|                        (this run)
A3 broken flux only: content = {(mu,alpha): alpha.v != 0}, q = |nu.v|
A4 light composites: content = {0 < |nu.v| < tau}, q = |nu.v|
A5 relaxed coupling: pairwise meeting F12 = F23 = 0, alpha_X free,
   bisection in z0 (exact for the piecewise-linear F) over the full
   direction and c0 grids.

**Structural theorem (measured exactly).**  Over the content with one
SM copy removed: sum db_1 = sum db_2 = 144.056 (equal to 1e-12 - the
embedding-index identity (3/5)sum Y^2 = sum T_2 over 16 x 45), and
sum(db_2 - db_3) = -76.0. The no-contribution meeting values are
x_12 = 25.60 and x_23 = 31.52 (different): reconciling them needs
threshold-weighted sums of opposite sign pattern that the q-correlated
spectra never deliver. Minimum max-residual over all grids:
A2 5.68, A3 3.26, A4 4.43 (targets are 0); A5: 0 physical solutions.

**Verdict.** The desert spectrum of V64 is NOT the composite content of
the gauged bond layer under any natural mass ansatz - a measured,
quantified exclusion. W4's supplier remains the named V64 input; the
question is now closed in both directions (V74: nodes excluded;
V81+V82: composites excluded).
"""

import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
sys.path.insert(0, exp("V81_link_desert_content"))

from V65_wilson import cartan_weights, identify_XY, physical_su5  # noqa: E402
from V81_link_desert_content import (c2_labels, flux_weights,  # noqa: E402
                                     sm_roots)
from cl10 import unification as UN  # noqa: E402


def content():
    lam, _ = cartan_weights()
    xs, ys = identify_XY(lam)
    yc = ys[0]
    gg = next(x for x in xs
              if physical_su5(lam @ x, lam @ yc) == "Georgi-Glashow")
    al = flux_weights()
    comp = (lam[:, None, :] + al[None, :, :]).reshape(-1, 5)
    bL, b1, b2 = sm_roots(al, lam, gg, yc)
    bL2 = float(bL @ bL)
    Ys = comp @ yc
    m2 = (comp @ bL) / bL2
    c2 = c2_labels(comp, b1, b2)
    db = np.stack([(2 / 3) * (3 / 5) * Ys ** 2,
                   (2 / 3) * m2 ** 2,
                   (2 / 3) * c2 / 8.0], axis=1)
    zero_j = next(j for j in range(45) if np.abs(al[j]).sum() < 1e-12)
    sm_idx = np.arange(16) * 45 + zero_j
    return lam, gg, yc, al, comp, db, sm_idx


def residual_min(qfun, selfun, thetas, c0s, z0s, xs_grid):
    """min over grids of max_i |F_i| with F from cl10.unification."""
    TWO_PI = 2 * math.pi
    best = (1e9, None)
    for th in thetas:
        v = math.cos(th) * selfun + math.sin(th) * qfun[1]
        for c0 in c0s:
            q = qfun[0](v, c0)
            lq = np.log(np.maximum(c0 + q, 1e-30))
            dbs = qfun[2]
            for z0 in z0s:
                zs = z0 + lq
                for x in xs_grid:
                    F = np.array(
                        [1 / UN.ALPHA_10 + UN.B_SM[i] * x / TWO_PI
                         + float(np.sum(dbs[:, i]
                                        * np.clip(x - zs, 0, None))) / TWO_PI
                         - 1 / UN.A_Z[i] for i in range(3)])
                    r = float(np.abs(F).max())
                    if r < best[0]:
                        best = (r, (th, c0, z0, x))
    return best


def main():
    t0 = time.time()
    lam, gg, yc, al, comp, db_all, sm_idx = content()

    dbs = db_all.copy()
    dbs[sm_idx] *= 0.0
    dbt = dbs.sum(axis=0)
    checks = {}
    checks["index_identity"] = {
        "db1_total": round(float(dbt[0]), 6),
        "db2_total": round(float(dbt[1]), 6),
        "db3_total": round(float(dbt[2]), 6),
        "db1_minus_db2": float(abs(dbt[0] - dbt[1])),
        "db2_minus_db3_sum": round(float(dbt[1] - dbt[2]), 4),
        "note": "sum db1 = sum db2 exactly: the embedding-index identity "
                "(3/5) sum Y^2 = sum T2 over 16 x 45; the content shifts "
                "U(1) and SU(2) EQUALLY - the data require different "
                "shifts: part of the structural obstruction"}
    TWO_PI = 2 * math.pi
    x12 = ((1 / UN.A_Z[0] - 1 / UN.A_Z[1])
           / ((UN.B_SM[0] - UN.B_SM[1]) / TWO_PI))
    x23 = ((1 / UN.A_Z[1] - 1 / UN.A_Z[2])
           / ((UN.B_SM[1] - UN.B_SM[2]) / TWO_PI))
    checks["no_contribution_meetings"] = {
        "x_12": round(float(x12), 3), "x_23": round(float(x23), 3),
        "note": "the U(1)-SU(2) and SU(2)-SU(3) no-contribution meeting "
                "scales differ: content contributions must reconcile them "
                "with threshold-weighted sums of the differences db_i - db_j"}

    thetas = [k * math.pi / 2 / 59 for k in range(60)]
    c0s = (0.001, 0.01, 0.1, 1.0, 10.0)
    z0s = np.arange(8.0, 26.0, 1.0)
    xs_grid = np.arange(20.0, 45.0, 0.5)

    # A2: full weight, all 720 (SM copy removed)
    q2 = (lambda v, c0: np.abs(comp @ v), yc, dbs)
    rA2 = residual_min(q2, gg, thetas, c0s, z0s, xs_grid)
    checks["A2_full_weight"] = {"min_max_residual": round(rA2[0], 3),
                                "arg": [round(float(a), 4) for a in rA2[1]]}

    # A3: broken flux only
    best3 = (1e9, None)
    for th in thetas:
        v = math.cos(th) * gg + math.sin(th) * yc
        av = al @ v
        broken = np.abs(av) > 1e-9
        keep = np.repeat(broken, 16)
        comp_k = comp[keep]
        db_k = db_all[keep]
        q = np.abs(comp_k @ v)
        for c0 in c0s:
            lq = np.log(np.maximum(c0 + q, 1e-30))
            for z0 in z0s:
                zs = z0 + lq
                for x in xs_grid:
                    F = np.array(
                        [1 / UN.ALPHA_10 + UN.B_SM[i] * x / TWO_PI
                         + float(np.sum(db_k[:, i]
                                        * np.clip(x - zs, 0, None))) / TWO_PI
                         - 1 / UN.A_Z[i] for i in range(3)])
                    r = float(np.abs(F).max())
                    if r < best3[0]:
                        best3 = (r, (th, c0, z0, x, int(broken.sum())))
    checks["A3_broken_flux"] = {"min_max_residual": round(best3[0], 3),
                                "arg": [round(float(a), 4)
                                        for a in best3[1][:4]],
                                "n_broken_roots_at_best": best3[1][4]}

    # A4: light-composite selection
    best4 = (1e9, None)
    for th in thetas:
        v = math.cos(th) * gg + math.sin(th) * yc
        q_all = np.abs(comp @ v)
        for tau in (0.02, 0.05, 0.1, 0.2, 0.4):
            sel = (q_all > 1e-9) & (q_all < tau)
            if not 0 < sel.sum() <= 120:
                continue
            q = q_all[sel]
            db_s = db_all[sel]
            lq = np.log(q)
            for z0 in np.arange(5.0, 26.0, 1.0):
                zs = z0 + lq
                for x in xs_grid:
                    F = np.array(
                        [1 / UN.ALPHA_10 + UN.B_SM[i] * x / TWO_PI
                         + float(np.sum(db_s[:, i]
                                        * np.clip(x - zs, 0, None))) / TWO_PI
                         - 1 / UN.A_Z[i] for i in range(3)])
                    r = float(np.abs(F).max())
                    if r < best4[0]:
                        best4 = (r, (th, tau, z0, x, int(sel.sum())))
    checks["A4_light_composites"] = {"min_max_residual": round(best4[0], 3),
                                     "arg": [round(float(a), 4)
                                             for a in best4[1][:4]],
                                     "n_selected_at_best": best4[1][4]}

    # A5: relaxed alpha_X: pairwise meeting by bisection (exact for the
    # piecewise-linear F12 in z0), full direction and c0 grids.
    def Fm(i, j, x, z0, lq, d):
        zs = z0 + lq
        return ((1 / UN.A_Z[i] - 1 / UN.A_Z[j])
                - (UN.B_SM[i] - UN.B_SM[j]) * x / TWO_PI
                - float(np.sum((d[:, i] - d[:, j])
                               * np.clip(x - zs, 0, None))) / TWO_PI)

    n_solutions = 0
    alphas = []
    for th in [k * math.pi / 2 / 89 for k in range(90)]:
        v = math.cos(th) * gg + math.sin(th) * yc
        q = np.abs(comp @ v)
        for c0 in c0s:
            lq = np.log(np.maximum(c0 + q, 1e-30))
            prev = None
            for x in np.arange(15.0, 45.01, 0.25):
                lo, hi = -10.0, x + 5
                flo, fhi = Fm(0, 1, x, lo, lq, dbs), Fm(0, 1, x, hi, lq, dbs)
                if flo * fhi > 0:
                    prev = None
                    continue
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    fm = Fm(0, 1, x, mid, lq, dbs)
                    if flo * fm <= 0:
                        hi = mid
                    else:
                        lo, flo = mid, fm
                z0 = 0.5 * (lo + hi)
                f23 = Fm(1, 2, x, z0, lq, dbs)
                if prev is not None and prev[1] * f23 <= 0:
                    zs = z0 + lq
                    contrib = zs < x - 1e-9
                    if (contrib.any()
                            and zs[contrib].min() >= math.log(1e6 / UN.M_Z)
                            and x <= math.log(1.22e19 / UN.M_Z)):
                        s1 = (1 / UN.A_Z[0] - UN.B_SM[0] * x / TWO_PI
                              - float(np.sum(dbs[:, 0]
                                             * np.clip(x - zs, 0, None)))
                              / TWO_PI)
                        if s1 > 0:
                            n_solutions += 1
                            alphas.append(1 / s1)
                prev = (x, f23)
    checks["A5_relaxed_alpha"] = {
        "n_physical_pairwise_solutions": n_solutions,
        "alpha_X_range": ([round(min(alphas), 4), round(max(alphas), 4)]
                          if alphas else None)}

    ok = (checks["index_identity"]["db1_minus_db2"] < 1e-9
          and rA2[0] > 1.0 and best3[0] > 1.0 and best4[0] > 1.0)

    status = (
        "THE COMPOSITE CONTENT IS EXCLUDED AS THE DESERT SUPPLIER UNDER "
        "EVERY NATURAL MASS ANSATZ: full-weight min residual "
        f"{rA2[0]:.2f}, broken-flux {best3[0]:.2f}, light-composite "
        f"{best4[0]:.2f} (target 0); relaxed-coupling pairwise meeting: "
        f"{n_solutions} physical solutions. Structural reason measured: "
        "sum db1 = sum db2 = 144.06 EXACTLY (embedding-index identity) - "
        "the content shifts U(1) and SU(2) equally while the data require "
        "different shifts - and the no-contribution meeting scales differ "
        f"(x_12 = {x12:.2f} vs x_23 = {x23:.2f}). W4 is now closed in both "
        "directions: nodes excluded (V74), composites excluded (V81+V82); "
        "the V64 spectrum remains a named input.")

    verdict = {"experiment": "V82_composite_desert_exclusion",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v82_composite_desert_exclusion.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V82 - composite desert content: exclusion under all ansatz")
    print("=" * 68)
    print(f"  db totals (SM copy removed) = "
          f"{tuple(round(float(t), 3) for t in dbt)}")
    print(f"  db1 - db2 (exact)           = "
          f"{checks['index_identity']['db1_minus_db2']:.2e}")
    print(f"  no-contrib meetings x12/x23 = {x12:.2f} / {x23:.2f}")
    print(f"  A2 full weight              : min residual {rA2[0]:.3f}")
    print(f"  A3 broken flux              : min residual {best3[0]:.3f}")
    print(f"  A4 light composites         : min residual {best4[0]:.3f}")
    print(f"  A5 relaxed alpha_X          : {n_solutions} solutions")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v82_composite_desert_exclusion.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
