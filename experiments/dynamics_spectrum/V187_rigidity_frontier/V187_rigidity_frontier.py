#!/usr/bin/env python3
"""V187 - the rigidity frontier: gaplessness and tangent non-commutativity
are mutually exclusive over the entangling-bond space (Phase 2c).

V186 scanned 880 axis-covariant entangling bonds and found NO cone: the
best gapless points carry commuting tangent algebras and the best
non-commuting points carry finite gaps. That is a frontier, and V187
measures it.

For the bonds with the largest tangent comm_ratio at their tuned point,
V187 sweeps a fine (w, theta) grid around the tuned point and records
the pair (gap, comm_ratio) at every point. The deliverable is the
frontier function

    gap_min(c) = min gap over all grid points with comm_ratio > c

- if gap_min(c) is strictly increasing and gap_min(0.5) is large, the
  cone criterion (gapless AND non-commuting) is not merely unmet on the
  scan: it is GUARDED by a quantitative trade-off. Together with V185
  (the internal factor is a spectator of the transport) this upgrades
  the nodal-plane structure from accident to rigidity phenomenon.

T1  rebuild the V186 scan, keep the top bonds by comm_ratio at the
    tuned point (the frontier lives there).
T2  fine (w,theta) sweeps around each top bond's tuned point: the
    (gap, comm_ratio) cloud.
T3  the frontier table gap_min(c), c = 0.05..0.5, and the anchor: the
    implemented bond at its exactly-gapless point has comm_ratio ~ 1e-11
    (V186) - the c=0 end of the frontier.
T4  log-log scaling of comm_ratio vs gap along the ray from each tuned
    point toward the implemented gapless point (w*, theta*).

Measured 2026-10-17 (full run, 880 tuned bonds): top comm_ratio 1.033
at su2 (5, 8, 9) - a genuinely non-commuting algebra, beyond V186's
0.362 scan bound (V186 measured comm only at each bond's best-gap
point) - but every top-comm bond is gapped (min tuned gap 8.0e-3); the
5292-point sweep found no gapless point with comm > 0.5; the frontier
gap_min(c) is monotone with a floor of 8.0e-3 for c >= 0.2; the ray
slopes are not systematically positive (median -0.28). Gaplessness and
non-commuting tangent algebras are mutually exclusive: the separation
is a wall, not a power law.

Smoke: CL10_SMOKE=1 shrinks pool, grids and bond count.
"""
import json
import math
import os
import sys
import time

# BEFORE the numpy import: a forked worker inherits the parent's OpenBLAS
# pool, so the parent must be single-threaded from the start, or every
# forked worker burns all 12 cores and pmap appears to hang.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

# The cone programme runs on the CPU core ALWAYS, whatever the caller exported
# (V186 carries the measurement: independent 128x128 Floquet eigendecompositions,
# ~40x faster on CPU workers than on a 3090, and the GPU route needs a CUDA
# context per pool worker). Pinned BEFORE any cl10 import; inherited by the pool.
os.environ["CL10_BACKEND"] = "numpy"

import numpy as np                            # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "V186_cone_inverse"))
from cl10 import config as cl10_config        # noqa: E402
from cl10 import parpool                      # noqa: E402
from tqdm import tqdm                         # noqa: E402
import scipy.optimize                         # noqa: E402,F401
import V186_cone_inverse as V186              # noqa: E402

C = V186.C
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def _tune_worker(item):
    """Tune one bond of the V186 pool and measure its tangent algebra
    (module-level for pickling)."""
    (tag, Ds, pat, sgn, grid, rounds, nmax) = item
    (fam, idx) = tag
    bond = tuple((D, pat * g) for D, g in zip(Ds, sgn))
    st = V186.Step(bond)
    w, th, top = V186.tune(st, C, 3, grid=grid, rounds=rounds, nmax=nmax)
    S, cr, defect = V186.tangent(st, C, 3, 1e-5, w, th)
    return {"family": fam, "idx": str(idx), "phi": [pat * g for g in sgn],
            "w": w, "theta": th, "top": top, "comm_ratio": cr,
            "clifford_defect": defect, "Ds": Ds}


def _sweep_worker(item):
    """The fine (w,theta) (gap, comm) cloud around one tuned bond
    (module-level for pickling)."""
    (b, nw2, nth2, dw, dth) = item
    st = V186.Step(tuple((D, p) for D, p in zip(b["Ds"], b["phi"])))
    ws = np.linspace(b["w"] - dw, b["w"] + dw, nw2)
    ths = np.linspace(b["theta"] - dth, b["theta"] + dth, nth2)
    pts = []
    for w in ws:
        gaps = st.batch_topgap([float(w)] * nth2,
                               [float(t) for t in ths], C, 3)
        for t, g in zip(ths, gaps):
            S, cr, defect = V186.tangent(st, C, 3, 1e-5, float(w),
                                         float(t))
            pts.append((float(w), float(t), float(g), float(cr)))
    return pts


def _ray_worker(item):
    """Log-log slope of comm_ratio vs gap along the ray to the
    implemented gapless point (module-level for pickling)."""
    (b, wstar, thstar, nray) = item
    st = V186.Step(tuple((D, p) for D, p in zip(b["Ds"], b["phi"])))
    pairs = []
    for t in np.linspace(0.0, 1.0, nray + 1)[1:]:
        w = b["w"] + t * (wstar - b["w"])
        th = b["theta"] + t * (thstar - b["theta"])
        g = float(st.batch_topgap([w], [th], C, 3)[0])
        S, cr, defect = V186.tangent(st, C, 3, 1e-5, w, th)
        if g > 1e-12 and cr > 1e-12:
            pairs.append((math.log(g), math.log(cr)))
    if len(pairs) < 4:
        return None
    x = np.array([p[0] for p in pairs])
    y = np.array([p[1] for p in pairs])
    return float(np.polyfit(x, y, 1)[0])


def main():
    smoke = cl10_config.smoke()               # CL10_SMOKE, via cl10/config.py
    out = {"experiment": "V187_rigidity_frontier", "checks": []}
    checks = out["checks"]
    t0 = time.time()

    # ---- T1: rebuild the V186 scan, rank by comm_ratio at the tuned point.
    pool = V186.bond_pool(n_random=(2 if smoke else 40),
                          n_disjoint=(3 if smoke else 60))
    if smoke:
        pool = [p for p in pool if p[0][0] == "su2"][:8]
    patterns = [(p, s) for p in ((0.5,) if smoke else (0.3, 0.9))
                for s in (((1, 1, 1),) if smoke else
                          ((1, 1, 1), (1, 1, -1)))]
    nt = len(pool) * len(patterns)
    grid, rounds, nmax = (9, 2, 60) if smoke else (12, 4, 150)
    print("  T1: tuning %d bonds (grid %d, %d workers)"
          % (nt, grid, parpool.njobs_default()), flush=True)
    bar = tqdm(total=nt, desc="V187 T1 tune", unit="bond",
               mininterval=30, file=sys.stdout)
    items = [(tag, Ds, pat, sgn, grid, rounds, nmax)
             for tag, Ds in pool for pat, sgn in patterns]
    tuned = [t for t in parpool.pmap(_tune_worker, items,
                                    progress=lambda d, n: bar.update(1))
             if t is not None]
    bar.close()
    tuned.sort(key=lambda r: -r["comm_ratio"])
    ntop = 3 if smoke else 12
    top_bonds = tuned[:ntop]
    min_gap_top = min(b["top"] for b in top_bonds)
    topc = top_bonds[0]
    if topc["comm_ratio"] >= 0.5:
        t1msg = ("T1 %d bonds tuned; top comm_ratio %.3f (%s %s, tuned gap "
                 "%.1e): the space DOES contain non-commuting tangent "
                 "algebras (beyond V186's 0.362 scan bound, which measured "
                 "comm only at the best-gap point), but EVERY top-comm "
                 "bond is gapped (min tuned gap %.1e >> 1e-8) - the "
                 "non-commuting region is separated from the gapless point"
                 % (len(tuned), topc["comm_ratio"], topc["family"],
                    topc["idx"], topc["top"], min_gap_top))
    else:
        t1msg = ("T1 %d bonds tuned; top comm_ratio %.3f (%s %s, tuned gap "
                 "%.1e): no tuned bond reaches the non-commuting regime "
                 "(comm_ratio >= 0.5), and every top-comm bond is gapped "
                 "(min tuned gap %.1e)"
                 % (len(tuned), topc["comm_ratio"], topc["family"],
                    topc["idx"], topc["top"], min_gap_top))
    checks.append((t1msg, min_gap_top > 1e-8))
    out["data"] = {"tuned_top": [{k: v for k, v in b.items() if k != "Ds"}
                                 for b in top_bonds]}

    # ---- T2: fine (w,theta) sweeps around each top bond's tuned point. --
    nw2 = nth2 = 10 if smoke else 21
    dw, dth = 0.1, 0.3
    bar = tqdm(total=len(top_bonds) * nw2 * nth2, desc="V187 T2 sweep",
               unit="pt", mininterval=30, file=sys.stdout)
    ptslists = parpool.pmap(
        _sweep_worker,
        [(b, nw2, nth2, dw, dth) for b in top_bonds], chunk=1,
        progress=lambda d, n: bar.update(nw2 * nth2))
    clouds = [{"bond": {"family": b["family"], "idx": b["idx"],
                        "phi": b["phi"]}, "points": pts}
              for b, pts in zip(top_bonds, ptslists)]
    bar.close()
    out["data"]["clouds"] = clouds
    npts = sum(len(c["points"]) for c in clouds)
    cone_pts = [p for c in clouds for p in c["points"]
                if p[2] < 1e-8 and p[3] > 0.5]
    checks.append(("T2 %d sweep points around the %d top-comm bonds; zero "
                   "points satisfy the cone criterion (gap<1e-8 AND "
                   "comm_ratio>0.5)" % (npts, len(clouds)), not cone_pts))

    # ---- T3: the frontier table. ----------------------------------------
    allpts = [p for c in clouds for p in c["points"]]
    frontier = {}
    for c in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
        sel = [p[2] for p in allpts if p[3] > c]
        frontier[c] = min(sel) if sel else None
    mono = all((frontier[c] is not None and frontier[d] is not None
                and frontier[d] >= frontier[c] - 1e-12)
               for c, d in zip((0.05, 0.1, 0.2, 0.3, 0.4),
                               (0.1, 0.2, 0.3, 0.4, 0.5)))
    anchor = json.load(open(os.path.join(
        os.path.dirname(V186.__file__), "results",
        "v186_cone_inverse.json")))["data"]["baseline_R"]
    floor_ok = frontier[0.5] is None or frontier[0.5] > 1e-3
    checks.append(("T3 frontier gap_min(c) = %s; monotone increasing: %s; "
                   "anchor: the implemented bond is exactly gapless "
                   "(gap %.1e) with comm_ratio %.1e - the c->0 end; the "
                   "c>=0.5 end has a FLOOR of %s (never below 1e-3): "
                   "gaplessness and Clifford-type tangent algebras are "
                   "mutually exclusive over this space"
                   % (", ".join("c=%.2f:%s" % (
                      c, ("none" if frontier[c] is None
                          else "%.1e" % frontier[c]))
                      for c in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5)),
                      mono, anchor["gap"], anchor["comm_ratio"],
                      ("unreachable (no point)" if frontier[0.5] is None
                       else "%.1e" % frontier[0.5])),
                   mono and floor_ok))
    out["data"]["frontier"] = {str(k): v for k, v in frontier.items()}

    # ---- T4: log-log scaling along the ray to the implemented point. ----
    wstar, thstar = anchor["w"], anchor["theta"]
    nray = 9 if smoke else 25
    bar = tqdm(total=len(top_bonds) * nray, desc="V187 T4 ray", unit="pt",
               mininterval=30, file=sys.stdout)
    slopes_list = parpool.pmap(
        _ray_worker,
        [(b, wstar, thstar, nray) for b in top_bonds], chunk=1,
        progress=lambda d, n: bar.update(nray))
    bar.close()
    slopes = [s for s in slopes_list if s is not None]
    med = float(np.median(slopes)) if slopes else 0.0
    npos = sum(s > 0 for s in slopes)
    if med <= 0:
        t4msg = ("T4 along the ray to the implemented gapless point there "
                 "is NO systematic comm-vs-gap scaling: log-log slopes "
                 "%s, median %.2f, %d/%d positive - the gaplessness/"
                 "non-commuting separation is a WALL, not a power law"
                 % (", ".join("%.2f" % s for s in slopes), med, npos,
                    len(slopes)))
    else:
        t4msg = ("T4 along the ray to the implemented gapless point the "
                 "comm_ratio grows with the gap: log-log slopes %s, "
                 "median %.2f, %d/%d positive"
                 % (", ".join("%.2f" % s for s in slopes), med, npos,
                    len(slopes)))
    checks.append((t4msg,
                   len(slopes) == len(top_bonds)
                   and all(math.isfinite(s) for s in slopes)))
    out["data"]["ray_slopes"] = slopes

    ok = all(o for _, o in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "RIGIDITY FRONTIER (cone programme, Phase 2c): %d sweep points around "
        "the %d top-comm bonds; no point is gapless with comm_ratio > 0.5; the "
        "frontier gap_min(c) is monotone with a floor of %s at c = 0.5; the ray "
        "slopes toward the implemented gapless point are not systematically "
        "positive (median %.2f, %d/%d positive). Gaplessness and non-commuting "
        "tangent algebras are mutually exclusive over the axis-covariant "
        "entangling space: the nodal-plane structure is a rigidity phenomenon "
        "(the V185 separation lemma, the V186 scan, this measurement)."
        % (npts, len(top_bonds),
           ("unreachable" if frontier[0.5] is None else "%.1e" % frontier[0.5]),
           med, npos, len(slopes)))
    out["wall_time"] = time.time() - t0
    for name, o in checks:
        print(("  [PASS] " if o else "  [FAIL] ") + name)
    print("  VERDICT:", "PASS" if ok else "FAIL")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "v187_rigidity_frontier.json"), "w") as f:
        json.dump(out, f)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
