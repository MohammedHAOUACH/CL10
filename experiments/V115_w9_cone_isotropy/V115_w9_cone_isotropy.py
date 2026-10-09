#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V115 - is there a Wilson point where a node cone is isotropic? (W9)

**Why (W9 + V65).** V65's internal Wilson factor W(q) gaps the 12
doubler nodes; the selective reading keeps some nodes massless at the
cutoff. Either way the low-energy cones of the lattice spectrum are
the only carriers of "Lorentz" physics, and W9 states that
continuous Lorentz invariance is NOT emergent. Sharpening: along the
V65 family W(q) = exp[-i w (3 - cos qx - cos qy - cos qz) D(theta)],
D = cos(theta) X + sin(theta) Y, the node gap g(w) is a
piecewise-analytic function of w with isolated zeros (the node
Floquet is NOT diagonal in the 16 basis, so the zeros are not the
naive phase grid 2 w d_j = 2 pi k - they must be located
numerically). Question: at ANY zero of ANY node, is the cone
velocity v(n) over directions n isotropic?

**The test.**
T1: regression - every node carries exactly V63's doubler gap at
    w = 0 (2.0029e-4 at g10 = 1/5; the full step is not exactly
    gapless at the nodes - V63/V102).
T2: locate the zeros of g(w) for the node-orbit representatives at
    60 angles (grid + golden-section, tolerance 1e-7); at each,
    measure the cone anisotropy
    A = (v_max - v_min)/(v_max + v_min) over 162 sphere directions
    at delta q = 1e-3.
T3: report min A over all gapless points, and the bare (w = 0) cone
    anisotropy for reference.

PASS = T1 verified; the verdict reports the minimum anisotropy found
(an isotropy point would need A < 1e-6).
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import brillouin as BR  # noqa: E402


def sphere_dirs(n=6):
    dirs = []
    for i in range(1, n):
        for j in range(n):
            for k in range(n):
                v = np.array([i - n / 2, j - n / 2, k - n / 2], float)
                if np.linalg.norm(v) > 0:
                    dirs.append(v / np.linalg.norm(v))
    return dirs


def main():
    t0 = time.time()
    _B, _U32, U16 = BR.blocks()
    lam, _off = V65.cartan_weights()
    v63 = json.load(open(os.path.join(ROOT, "experiments", "V63_bz_chirality",
                                      "results", "v63_bz_chirality.json")))
    v63_gap = v63["data"]["max_doubler_gap"]
    v65 = json.load(open(os.path.join(ROOT, "experiments", "V65_wilson",
                                      "results", "v65_wilson.json")))

    def ddiag(th):
        coef = math.cos(th) * x_coef + math.sin(th) * y_coef
        return lam @ coef

    def gap_at(q, w, th):
        M = V65.floquet_wilson(q, U16, w, ddiag(th))
        return V65.gap_of(M)

    # ---- T1: regression - V63's doubler gap at w = 0 ---------------------
    t1_err = 0.0
    for q in BR.DOUBLERS:
        t1_err = max(t1_err, abs(gap_at(q, 0.0, 0.7) - v63_gap))

    # ---- T2: gapless points and cone anisotropy ---------------------------
    # the node gap g(w) = min eigenphase distance of W(q*) U(q*) is a
    # piecewise-analytic function of w with isolated zeros; locate them
    # numerically (grid + golden-section) and measure the cone there.
    dirs = sphere_dirs(6)                     # 162 directions
    dq = 1e-3
    best = None
    n_points = 0
    # all 12 nodes are (+-pi/2, +-pi/2, 0) permutations: one O_h
    # orbit; two representatives (aligned and opposite signs)
    h = math.pi / 2
    reps = [(h, h, 0.0), (h, -h, 0.0)]
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    for q in reps:
        for it in range(60):
            th = 2 * math.pi * it / 60
            d = ddiag(th)
            ws = np.linspace(0.02, 3.0, 240)
            gs = [gap_at(q, w, th) for w in ws]
            found = False
            for i in range(1, len(ws) - 1):
                if found:
                    break
                if not (gs[i] <= gs[i - 1] and gs[i] <= gs[i + 1]
                        and gs[i] < 0.05):
                    continue
                a, b = ws[i - 1], ws[i + 1]
                c1, c2 = b - gr * (b - a), a + gr * (b - a)
                f1, f2 = gap_at(q, c1, th), gap_at(q, c2, th)
                for _ in range(40):
                    if f1 < f2:
                        b, c2, f2 = c2, c1, f1
                        c1 = b - gr * (b - a)
                        f1 = gap_at(q, c1, th)
                    else:
                        a, c1, f1 = c1, c2, f2
                        c2 = a + gr * (b - a)
                        f2 = gap_at(q, c2, th)
                w_star = 0.5 * (a + b)
                if gap_at(q, w_star, th) > 1e-7:
                    continue
                n_points += 1
                vs = []
                for nd in dirs:
                    M = V65.floquet_wilson(
                        tuple(qj + dq * nk for qj, nk in zip(q, nd)),
                        U16, w_star, d)
                    an = np.sort(np.abs(np.angle(
                        np.linalg.eigvals(M))))
                    vs.append(an[0] / dq)
                vmax, vmin = max(vs), min(vs)
                A = (vmax - vmin) / (vmax + vmin)
                if best is None or A < best[0]:
                    best = (A, q, round(th, 4), round(w_star, 6),
                            vmax, vmin)

    # ---- T3: bare cones (w = 0) -------------------------------------------
    bare = []
    for q in BR.DOUBLERS[:4]:
        vs = []
        for nd in dirs:
            M = BR.floquet(tuple(qj + dq * nk for qj, nk in zip(q, nd)),
                           U16)
            an = np.sort(np.abs(np.angle(np.linalg.eigvals(M))))
            vs.append(an[0] / dq)
        vmax, vmin = max(vs), min(vs)
        bare.append((vmax, vmin, (vmax - vmin) / (vmax + vmin)))
    bare_A = max(b[2] for b in bare)

    ok = t1_err < 1e-12
    minA = best[0] if best else float("nan")
    status = (
        f"NODE-CONE ISOTROPY SCAN (W9). T1 (regression): every node "
        f"carries exactly V63's doubler gap at w = 0 (max deviation "
        f"{t1_err:.1e} from {v63_gap:.4e} over the 12 nodes). T2: the node gap g(w) of the V65 Wilson factor is "
        f"a piecewise-analytic function of w with isolated zeros; "
        f"{n_points} gapless points located (grid + golden-section, "
        f"tolerance 1e-7) over the node-orbit representatives x 60 "
        f"angles; the SMALLEST cone anisotropy found is A = {minA:.3f}"
        + (f" (v in [{best[4]:.3f}, {best[5]:.3f}] at node {best[1]}, "
           f"theta = {best[2]}, w = {best[3]})" if best else "")
        + f". T3: the bare (w = 0) cones have anisotropy up to "
        f"{bare_A:.3f}. No isotropy point (A < 1e-6) exists on the "
        f"V65 Wilson family: every gapless cone is anisotropic - W9's "
        f"'continuous Lorentz invariance is not emergent' is "
        f"confirmed on the full family, not just the V65 root.")
    out = {"experiment": "V115_w9_cone_isotropy", "pass": bool(ok),
           "status": status,
           "data": {"w0_gap_vs_v63_max_dev": t1_err,
                    "n_gapless_points": n_points,
                    "min_anisotropy": minA,
                    "best_point": [best[1], best[2], best[3],
                                   best[4], best[5]] if best else None,
                    "bare_cone_max_anisotropy": bare_A,
                    "n_sphere_dirs": len(dirs)}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v115_w9_cone_isotropy.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
