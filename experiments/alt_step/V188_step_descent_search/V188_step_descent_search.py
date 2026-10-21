#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V188 - Adam descent on the 45-plane step angles, factorized route
(alt-step, 2/4)

**The search.**  V184 defined the cone loss L(theta) = A + 100 * omega(q*)^2
on the (2,2,2) cell with the baseline stages, anchored to the implemented step
at theta = g_10 * 1.  This experiment runs the pre-registered alternating
scheme from 25 seeds: 10 cycles of [refine the node q* (Nelder-Mead on the
gap), then 30 Adam steps on theta at that node (lr 0.02)], keeping the best
iterate of the trajectory (Adam overshoots on this rugged landscape).  Seeds:
the baseline g_10 * 1; 16 Uniform(0, 2 pi) draws (rng 185); 8 draws from
{pi/4, pi/2}^45 (rng 186).

**No gradient.**  The step factorizes exactly,

    M(q, theta) = P(q) (x) U(theta)^ns

- the bond layer P(q) carries no step angle at all - so the loss is a max over
the (i, j) phase pairs of cos(phi_i(q) + psi_j(theta)): piecewise smooth,
kinked wherever the arg-max pair changes.  There is no gradient to
differentiate.  The descent therefore takes central finite differences on the
45 angles, 90 forwards per Adam step, each forward a 16x16 eigendecomposition
plus a 215 x 128 pair max - cheaper than the autodiff it replaces, and the
whole search runs on the CPU core in minutes instead of an hour on the card.
V184's own anchor (d) established that its eigenvector subgradient and finite
differences agree to 2.7e-8, so this is the same subgradient.  The factorized
core is V189's `Core` object, imported from V189: the descent and the census
evaluate one and the same function.

**Candidacy (pre-registered, unchanged).**  A < 1e-6 AND omega(q*) < 1e-6: an
isotropic gapless cone reached by descent.  The hard tier (unitarity - exact
by construction here, closure order <= 64, front speed 1, chirality) is
applied to candidates in V190.  A negative result is a valid deliverable -
gates are never relaxed.

**Anchors.**  (a) the baseline seed row reproduces V184's A(theta0) to 1e-3
(A is a knife-edge ratio, dA/domega_min ~ 5e3, and the node landing jitters);
(b) at one node of the baseline step, the factorized and the dense (128^3)
routes agree: lam_max to 1e-12, and the cone metric to 1e-3 - the residual
there is the arccos amplification of a 1e-15 eigenvalue error at a gapless
point, not a difference of spectra.
"""

import importlib.util
import json
import os
import sys
import time
import os

os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import supercell as SC

_spec = importlib.util.spec_from_file_location(
    "V189_cell_step_joint_search",
    os.path.join(ROOT, "experiments", "alt_step",
                 "V189_cell_step_joint_search",
                 "V189_cell_step_joint_search.py"))
V189 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V189)

G10, DQ = SC.G10, V189.DQ
N_CYCLES = 10
N_TH = 30
LR_TH = 0.02
FD = 1e-5


def metrics(ph_st, ph_v):
    """The V115 cone metric and the gap at one node, from the phase match."""
    s = ph_st[:, :, None] + ph_v[None, None, :]
    w = np.abs(np.angle(np.exp(1j * s))).min(axis=(1, 2)) / DQ
    vmax, vmin = w.max(), w.min()
    return float((vmax - vmin) / (vmax + vmin)) if vmax + vmin > 0 else 1.0


def main():
    t0 = time.time()
    print("=" * 68)
    print("V188 - Adam descent over the 45-plane step angles "
          "(factorized step, numpy CPU, finite-difference subgradient)")
    print("=" * 68)
    C = V189.Core((2, 2, 2), "baseline")

    def loss(th, q, ph_st):
        ph_v = C.v_phases(th)
        om = C.gap(q, ph_v)
        return metrics(ph_st, ph_v) + 100.0 * om * om

    def grad(th, q, ph_st):
        """Central differences on the 45 angles: the loss is a max over phase
        pairs, so this is a subgradient of the same object V184 differentiated
        (the two agree to 2.7e-8, V184 anchor (d))."""
        g = np.zeros(45)
        for k in range(45):
            tp, tm = th.copy(), th.copy()
            tp[k] += FD
            tm[k] -= FD
            g[k] = (loss(tp, q, ph_st) - loss(tm, q, ph_st)) / (2.0 * FD)
        return g

    def adam(th, q, ph_st):
        m = v = np.zeros(45)
        for t in range(1, N_TH + 1):
            g = grad(th, q, ph_st)
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g * g
            th = th - LR_TH * (m / (1.0 - 0.9 ** t)) / (
                np.sqrt(v / (1.0 - 0.999 ** t)) + 1e-8)
        return th

    seeds = [("baseline g_10*1", np.full(45, G10))]
    rng = np.random.default_rng(185)
    for s in range(16):
        seeds.append((f"uniform-{s}", rng.uniform(0.0, 2.0 * np.pi, 45)))
    rng2 = np.random.default_rng(186)
    for s in range(8):
        seeds.append((f"quarter-{s}",
                      rng2.choice([np.pi / 4, np.pi / 2],
                                  45).astype(np.float64)))

    rows = []
    for name, th in seeds:
        th_in = th.copy()
        ph_v = C.v_phases(th)
        q, _ = C.find_node(ph_v)
        ph_st = C._p_phases_many(q[None, :] + DQ * C.dirs)
        a0, o0 = metrics(ph_st, ph_v), C.gap(q, ph_v)
        best = (a0 + 100.0 * o0 * o0, th.copy())
        for _ in range(N_CYCLES):
            th = adam(th, q, ph_st)
            ph_v = C.v_phases(th)
            q, _ = C.find_node(ph_v)
            ph_st = C._p_phases_many(q[None, :] + DQ * C.dirs)
            ll = loss(th, q, ph_st)
            if ll < best[0]:
                best = (ll, th.copy())
        th = best[1]
        ph_v = C.v_phases(th)
        q, _ = C.find_node(ph_v)
        ph_st = C._p_phases_many(q[None, :] + DQ * C.dirs)
        aa, om = metrics(ph_st, ph_v), C.gap(q, ph_v)
        rows.append({"seed": name, "L": aa + 100.0 * om * om, "A": aa,
                     "A0": a0, "omega_node": om,
                     "dtheta": float(np.linalg.norm(th - th_in)),
                     "q_star": [float(x) for x in q],
                     "theta_star": [float(x) for x in th]})
        print(f"  {name:>16}: A {aa:.4f}  omega(q*) {om:.2e}  "
              f"L {rows[-1]['L']:.4f}  ||dtheta|| {rows[-1]['dtheta']:.3f} "
              f"[{time.time() - t0:.0f}s]")

    # ---- anchor (a): the theta0 row reproduces V184's baseline A.  Gate 1e-3:
    # A is a knife-edge ratio (dA/domega_min ~ 5e3) and the Nelder-Mead node
    # landing jitters ~1e-5 in q*, which moves A by ~1e-4. ----
    p184 = os.path.join(ROOT, "experiments", "alt_step",
                        "V184_step_loss_landscape", "results",
                        "v184_step_loss_landscape.json")
    with open(p184) as fh:
        A184 = json.load(fh)["data"]["rows"][0]["A"]
    devA = abs(rows[0]["A0"] - A184)

    # ---- anchor (b): the factorized route against V184's dense (128^3) route,
    # at one node of the baseline step.  lam_max is the well-conditioned
    # comparison; the cone metric carries the arccos amplification. ----
    R = V189.V184.build_core()
    th0 = np.full(45, G10)
    q_ref = np.array(R["find_node"](th0))
    lam_ref = float(R["lam_max"](q_ref, th0))
    a_ref = float(R["aniso"](q_ref, th0))
    ph0 = C.v_phases(th0)
    lam_here = float(np.cos(C.gap(q_ref, ph0)))
    ph_st = C._p_phases_many(q_ref[None, :] + DQ * C.dirs)
    dev_lam, dev_A = abs(lam_ref - lam_here), abs(a_ref - metrics(ph_st, ph0))

    cands = [r for r in rows
             if r["A"] < 1e-6 and r["omega_node"] < 1e-6]
    best = min(rows, key=lambda r: r["A"])

    ok = devA < 1e-3 and dev_lam < 1e-12 and dev_A < 1e-3 and len(rows) == 25

    verdict = {
        "id": "V188_step_descent_search",
        "tier": 4,
        "property": "Adam descent on the V184 cone loss over the 45 per-plane "
                    "angles of the (2,2,2)-cell step, 25 pre-registered seeds, "
                    "alternating node/step scheme (10 cycles x [Nelder-Mead "
                    "node + 30 theta-steps]); the step is evaluated through "
                    "its exact factorization P(q) (x) U(theta)^ns and the "
                    "descent takes finite-difference subgradients",
        "residual": None,
        "pass": bool(ok),
        "status": f"DESCENT RUN (25 seeds, factorized step, no gradient). "
                  f"Baseline row reproduces V184 A to {devA:.1e}; the "
                  f"factorized route reproduces V184's dense route at one "
                  f"node to lam_max {dev_lam:.1e} and A {dev_A:.1e}. Best A "
                  f"over all seeds: {best['A']:.4f} (seed {best['seed']}); "
                  f"candidates (A < 1e-6 and gapless): {len(cands)}. Baseline "
                  f"final A = {rows[0]['A']:.4f} (||dtheta|| "
                  f"{rows[0]['dtheta']:.3f}); random-seed omega(q*) min "
                  f"{min(r['omega_node'] for r in rows[1:]):.2e}. Route: "
                  f"numpy CPU.",
        "data": {"rows": rows, "n_cycles": N_CYCLES, "lr_th": LR_TH,
                 "fd_step": FD, "anchor_dev_A": devA,
                 "anchor_lam_dev": dev_lam, "anchor_A_dev": dev_A,
                 "n_candidates": len(cands), "best": best},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v188_step_descent_search.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v188_step_descent_search.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
