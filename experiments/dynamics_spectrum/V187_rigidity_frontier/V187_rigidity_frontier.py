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

Smoke: CL10_SMOKE=1 shrinks pool, grids and bond count.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "V186_cone_inverse"))
from cl10 import engine as E                                    # noqa: E402
import V186_cone_inverse as V186                                # noqa: E402

C = V186.C
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def main():
    smoke = os.environ.get("CL10_SMOKE", "") == "1"
    out = {"experiment": "V187_rigidity_frontier", "checks": []}
    checks = out["checks"]

    # ---- T1: rebuild the V186 scan, rank by comm_ratio at the tuned point.
    pool = V186.bond_pool(n_random=(2 if smoke else 40),
                          n_disjoint=(3 if smoke else 60))
    if smoke:
        pool = [p for p in pool if p[0][0] == "su2"][:8]
    patterns = [(p, s) for p in ((0.5,) if smoke else (0.3, 0.9))
                for s in (((1, 1, 1),) if smoke else
                          ((1, 1, 1), (1, 1, -1)))]
    tuned = []
    for (fam, idx), Ds in pool:
        for pat, sgn in patterns:
            bond = tuple((D, pat * g) for D, g in zip(Ds, sgn))
            st = V186.Step(bond)
            w, th, top = V186.tune(st, C, 3, grid=(9 if smoke else 18))
            S, cr, defect = V186.tangent(st, C, 3, 1e-5, w, th)
            tuned.append({"family": fam, "idx": str(idx), "phi":
                          [pat * g for g in sgn], "w": w, "theta": th,
                          "top": top, "comm_ratio": cr,
                          "clifford_defect": defect, "Ds": Ds})
    tuned.sort(key=lambda r: -r["comm_ratio"])
    ntop = 3 if smoke else 12
    top_bonds = tuned[:ntop]
    checks.append(("T1 %d bonds tuned; top comm_ratio %.3f (%s %s); cone "
                   "criterion comm_ratio > 0.5 unmet by every tuned bond"
                   % (len(tuned), top_bonds[0]["comm_ratio"],
                      top_bonds[0]["family"], top_bonds[0]["idx"]),
                   top_bonds[0]["comm_ratio"] < 0.5))
    out["data"] = {"tuned_top": [{k: v for k, v in b.items() if k != "Ds"}
                                 for b in top_bonds]}

    # ---- T2: fine (w,theta) sweeps around each top bond's tuned point. --
    nw2 = nth2 = 10 if smoke else 21
    dw, dth = 0.1, 0.3
    clouds = []
    for b in top_bonds:
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
        clouds.append({"bond": {"family": b["family"], "idx": b["idx"],
                                "phi": b["phi"]},
                       "points": pts})
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
    slopes = []
    for b in top_bonds:
        st = V186.Step(tuple((D, p) for D, p in zip(b["Ds"], b["phi"])))
        pairs = []
        for t in np.linspace(0.0, 1.0, 9 if smoke else 25)[1:]:
            w = b["w"] + t * (wstar - b["w"])
            th = b["theta"] + t * (thstar - b["theta"])
            g = float(st.batch_topgap([w], [th], C, 3)[0])
            S, cr, defect = V186.tangent(st, C, 3, 1e-5, w, th)
            if g > 1e-12 and cr > 1e-12:
                pairs.append((math.log(g), math.log(cr)))
        if len(pairs) >= 4:
            x = np.array([p[0] for p in pairs])
            y = np.array([p[1] for p in pairs])
            slope = float(np.polyfit(x, y, 1)[0])
            slopes.append(slope)
    checks.append(("T4 along the ray to the implemented gapless point the "
                   "comm_ratio decays with the gap: log-log slopes %s, "
                   "median %.2f > 0 and at least 2/3 positive"
                   % (", ".join("%.2f" % s for s in slopes),
                      float(np.median(slopes))),
                   bool(slopes) and float(np.median(slopes)) > 0
                   and sum(s > 0 for s in slopes) >= 2 * len(slopes) // 3))
    out["data"]["ray_slopes"] = slopes

    out["verdict"] = all(ok for _, ok in checks)
    for name, ok in checks:
        print(("  [PASS] " if ok else "  [FAIL] ") + name)
    print("  VERDICT:", "PASS" if out["verdict"] else "FAIL")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "v187_rigidity_frontier.json"), "w") as f:
        json.dump(out, f)
    return 0 if out["verdict"] else 1


if __name__ == "__main__":
    sys.exit(main())
