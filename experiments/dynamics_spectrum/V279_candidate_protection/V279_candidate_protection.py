#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V279 - is the candidate protected? (tier 2, hard, O10)

**The question (pre-registered).** Do on-site perturbations that
respect the Cesaro commutant of the step leave the node gapless? If
V277 and V278 produced no candidate, the test runs as a CONTROL on
the star-module two-fold line (V219 anchor, V220/V227 machinery) and
its verdict carries a "control:" prefix and gates nothing.

**The object (closed input set).** The candidate step by value from
the V277/V278 artifacts (import of their modules rebuilds the
construction exactly: the artifacts carry the 18/30 best parameters);
or the published star-module step at the closed clock at the V219
booked line point q0 = (0.825816499798531, 0, 0) with m = 2.
Machinery: V220 `commutant` (Cesaro, CESARO = 240), `cluster_frame`,
`circle_holonomy`; V227 `phase_gap`, `refine`, `exponent_away`
(their modules imported, not restated). Perturbation: the Cesaro
average C of a random Hermitian (seed 279), scaled to ||C|| = 1,
switched on as exp(-i theta C) on the left of the step; gap versus
theta at thetas (1e-3, 1e-2) (V220's THETAS).

**Measured.** mode; commutant dimension proxy = the cluster-
separation and the obstruction ||U C U^-1 - C|| booked by V220's
commutant; the node gap; induced gaps at both thetas; the opening
(gap(1e-2) - gap(0)) and the slope (gap(1e-2) - gap(1e-3))/9e-3;
the obstruction norm; the circle holonomy charge of the node; for
the control, the transverse exponent of V227 (None = line).

**Acceptance.** pass = the test ran and the defect is booked; status
CANDIDATE-PROTECTED iff gap(1e-2) < gap(0) + 1e-6, else
CANDIDATE-MASSIVE with the opening slope and the invariant-projector
norm booked.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("CL10_BACKEND", "numpy")  # launcher may force cupy

import numpy as np                                        # noqa: E402
from scipy.linalg import expm                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
sys.path.insert(0, exp("V277_bond_decoration_cone_scan"))
sys.path.insert(0, exp("V278_range2_parity_cone"))
sys.path.insert(0, exp("V220_nodal_line_charge"))
sys.path.insert(0, exp("V227_locus_doubling_charges"))

import V277_bond_decoration_cone_scan as V277             # noqa: E402
import V278_range2_parity_cone as V278                    # noqa: E402
import V220_nodal_line_charge as V220                     # noqa: E402
import V227_locus_doubling_charges as V227                # noqa: E402
from cl10 import starcell as X                            # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import engine as E                              # noqa: E402

GCLOSED = math.pi / 4.0


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    t0 = time.time()
    checks, data = [], {}

    b277 = art("V277_bond_decoration_cone_scan", "branch")
    b278 = art("V278_range2_parity_cone", "branch")

    # ---- assemble the node and its step ----------------------------------
    if b277.endswith("FOUND") and not b278.endswith("FOUND"):
        mode, level = "candidate:V277", 128
        x = art("V277_bond_decoration_cone_scan", "data.search.best.x")
        amps = [(x[3 * t], x[3 * t + 1], x[3 * t + 2]) for t in range(6)]
        q_star = (0.0, 0.0, 0.0)

        def step(q):
            return V277.bloch_dec(q, amps)
    elif b278.endswith("FOUND"):
        mode, level = "candidate:V278", 128
        p = np.array(art("V278_range2_parity_cone", "data.search.best.x"))
        q_star = (0.0, 0.0, 0.0)

        def step(q):
            return V278.bloch30(q, p)
    else:
        mode, level = "control:star-module two-fold line", 1024
        q0 = list(art("V219_nodal_line_anisotropy", "data.q0"))
        CELL = X.module_cell()
        STAGES = X.module_stages(order="phase-major")
        AXES = np.asarray(X.star_axes(), float)
        B16, _, U16 = SC.sectors(GCLOSED)

        def step(q):
            return E.asnumpy(X.bloch_star(CELL, U16, GCLOSED, STAGES, q,
                                          AXES, bond_g10=GCLOSED))
        # T0 (control): the V219 booked node gap at q0
        gap_anchor = float(art("V219_nodal_line_anisotropy",
                               "data.cone_star.gap"))
        gq = V227.phase_gap(np.asarray(q0), step)
        t0_ok = bool(abs(gq - gap_anchor) < 1e-9)
        data["t0_control_node"] = {"q0": q0, "gap_measured": gq,
                                   "gap_V219_artifact": gap_anchor}
        checks.append({"name": "T0 control node: V219 line point "
                              "q0 = {} re-derived gap {:.2e} vs booked "
                              "{:.2e} (tol 1e-9)".format(q0, gq,
                                                         gap_anchor),
                       "pass": t0_ok})
        if not t0_ok:
            data["status"] = "MACHINERY-DRIFT"
            print("V279 MACHINERY-DRIFT: the V219 control node does not "
                  "reproduce; nothing downstream interpreted", flush=True)
            _emit(data, checks, ok=False)
            return 1
        q_star, _ = V227.refine(np.asarray(q0), step)
        data["refined_q"] = [float(v) for v in q_star]

    # ---- T0: the node itself ---------------------------------------------
    ph = np.sort(np.abs(np.angle(np.linalg.eigvals(step(np.asarray(
        q_star, float))))))
    m = int((ph < V220.NODE_TOL).sum())
    gap0 = float(ph.min())
    data["node"] = {"mode": mode, "level": level, "m": m,
                    "gap": gap0, "m_next": float(ph[m]) if m < len(ph) else None}

    # ---- the commutant and the switched-on invariant ---------------------
    U0 = np.asarray(step(np.asarray(q_star, float)), complex)
    cm = V220.commutant(U0, rng_seed=279)
    Cnorm, Coper = cm["norm"], cm["operator"]
    thetas = list(V220.THETAS)
    gaps = []
    for th in thetas:
        P = expm(-1j * th * Coper)
        gth = float(np.min(np.abs(np.angle(np.linalg.eigvals(P @ U0)))))
        gaps.append(gth)
    opening = gaps[-1] - gap0
    slope = (gaps[-1] - gaps[0]) / (thetas[-1] - thetas[0])
    data["commutant"] = {"obstruction": cm["obstruction"],
                         "norm": Cnorm,
                         "thetas": thetas, "gaps": gaps,
                         "opening": opening, "slope": slope}

    # ---- holonomy of the node (booked data) ------------------------------
    hol = None
    try:
        if m >= 2:
            plane = ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
            hol = V220.circle_holonomy(step, np.asarray(q_star, float),
                                       plane, 1e-2, m)
    except Exception as exc:                        # noqa: BLE001
        hol = {"error": str(exc)}
    data["holonomy"] = hol
    if mode.startswith("control"):
        data["transverse_exponent"] = V227.exponent_away(
            np.asarray(q_star, float), step)

    protected = bool(opening < 1e-6)
    status = ("CANDIDATE-PROTECTED" if protected else "CANDIDATE-MASSIVE")
    if mode.startswith("control"):
        status = "control:" + status
    data["status"] = status
    print("V279 {}: mode {}, m = {}, gap(0) {:.2e}, gap(1e-2) {:.2e}, "
          "opening {:.2e}, slope {:.3e}, commutant obstruction {:.1e}, "
          "norm {:.3f}".format(status, mode, m, gap0, gaps[-1], opening,
                               slope, cm["obstruction"], Cnorm), flush=True)
    checks.append({"name": "T1 protection test: opening at theta = 1e-2 "
                          "is {:.2e} (threshold 1e-6), slope {:.3e}, "
                          "invariant projector norm {:.3f}, obstruction "
                          "{:.1e}".format(opening, slope, Cnorm,
                                          cm["obstruction"]),
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    print("[V279 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V279_candidate_protection", "tier": 2,
           "property": "is the node produced by V277/V278 (or, as a "
                       "control, the star-module two-fold line) protected "
                       "against invariant on-site perturbations "
                       "(tier 2, O10)",
           "residual": True,
           "programme": "V279",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v279_candidate_protection.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V279 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
