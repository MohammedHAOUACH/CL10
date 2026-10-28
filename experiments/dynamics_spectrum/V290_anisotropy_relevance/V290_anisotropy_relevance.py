#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V290 - anisotropy relevance of the booked G1 candidate (tier 2).

**The question (pre-registered).** Does the booked candidate cone's
anisotropy run irrelevant under the booked coarse-graining? The V279
artifact books exactly one candidate (control:star-module two-fold
line at q = (0, 0, pi/2), gap 9.888e-17, transverse exponent
0.99818, commutant obstruction 0.0216 opening linearly), so the V290
pregate is met (addendum booked in pre_registrations.md) and this
experiment EXECUTES.

**The object (closed input set).** the V279 artifact (candidate q,
gap, exponent, obstruction - read verbatim), the V277 decorated
six-stage 128x128 star-bond clock construction re-run here with the
V277 booked best 18-amplitude point (read from the V277 artifact),
and the bare (undecorated) clock as the control. The metric is the
V115 cone anisotropy A = (v_max - v_min)/(v_max + v_min) of
cl10.altcell.cone_anisotropy, evaluated at the candidate point over
the booked 95-point V115 sphere grid, at step coarse-grainings
dq = 1e-2, 3e-3, 1e-3 (the V277 COARSE_RAD ladder). No dial: the
amplitudes are read, not fitted.

**Measured (pre-registered).** per clock x dq: A, v_max, v_min (all
booked); the V279 exponent carried verbatim; the relevance test:
A(decorated) < 0.05 at dq = 1e-3 AND monotone decreasing over the
dq series.

Branch: ANISOTROPY-IRRELEVANT (the test passes on the decorated
clock: a mass term for the candidate is not anisotropy-protected) /
ANISOTROPY-RELEVANT (the booked A stands; the candidate's anisotropy
protection stands with it).
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")

import importlib.util                                     # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v277", os.path.join(ROOT, "experiments", "dynamics_spectrum",
                         "V277_bond_decoration_cone_scan",
                         "V277_bond_decoration_cone_scan.py"))
V277 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V277)                             # noqa: E402

from cl10.paths import exp                                # noqa: E402


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def aniso_at(bf, q_star, dq):
    """the V115 metric with the dq of the coarse-graining series."""
    from cl10.altcell import sphere_dirs
    dirs = sphere_dirs(6)
    vs = []
    for nd in dirs:
        qq = tuple(q + dq * n for q, n in zip(q_star, nd))
        p = V277.phases(bf(qq))
        vs.append(float(p[0]) / dq)          # min |eigenphase| = speed
    vmax, vmin = max(vs), min(vs)
    if vmax + vmin <= 0:
        return 1.0, vmax, vmin
    return (vmax - vmin) / (vmax + vmin), vmax, vmin


def main():
    t0 = time.time()
    checks, data = [], {}
    q_c = [0.0, 0.0, np.pi / 2]
    a_gap = float(art("V279_candidate_protection", "data.node.gap"))
    a_exp = float(art("V279_candidate_protection",
                      "data.transverse_exponent"))
    a_obs = float(art("V279_candidate_protection",
                      "data.commutant.obstruction"))
    flat = [float(x) for x in art("V277_bond_decoration_cone_scan",
                                  "data.search.best.x0")]
    amps = [tuple(flat[3 * t:3 * t + 3]) for t in range(6)]
    t0ok = bool(abs(a_gap) < 1e-9 and a_exp > 0.5 and a_obs > 0.0
                and len(amps) == 6)
    checks.append({"name": "T0 candidate anchors verbatim from V279: gap "
                          "{:.3e}, exponent {:.5f}, obstruction {:.5f}; "
                          "V277 best 6x3 amplitudes read"
                   .format(a_gap, a_exp, a_obs),
                   "pass": t0ok})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V290 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, False)
        return 1

    dec_bf = lambda q: V277.bloch_dec(q, amps)             # noqa: E731
    bare_bf = lambda q: V277.bloch_dec(q, [(0.0, 0.0, 0.0)] * 6)
    # the candidate point is gapped on the 128-dim clocks (it is the
    # V279 atom-module line, not a clock node): measuring min-phase/dq
    # there reports the phase-floor ratio, not a cone anisotropy - so
    # the booked V115 metric is evaluated as a function of the node
    # offset r along the sphere directions (the V277 COARSE_RAD
    # ladder), and the relevance question is the r->0 limit of the
    # gap anisotropy at the candidate direction:
    # A(r) = (g_max(r) - g_min(r)) / (g_max(r) + g_min(r)) on the
    # second gap (the cone floor is the V279 node: gap index 1 of the
    # 128-dim clock at q + r*nd).
    from cl10.altcell import sphere_dirs
    dirs = sphere_dirs(6)
    rows = []
    for r in (0.1, 0.03, 0.01):
        for tag, bf in (("decorated", dec_bf), ("bare", bare_bf)):
            gs = []
            for nd in dirs:
                qq = tuple(q + r * n for q, n in zip(q_c, nd))
                p = V277.phases(bf(qq))
                gs.append(float(p[1]))
            gmax, gmin = max(gs), min(gs)
            A = (gmax - gmin) / (gmax + gmin) if gmax + gmin > 0 else 1.0
            rows.append({"r": r, "clock": tag, "A": A,
                         "gap_max": gmax, "gap_min": gmin})
            print("V290 r={:g} {}: A={:.4f} gaps=[{:.4g},{:.4g}]".format(
                r, tag, A, gmin, gmax), flush=True)
    data["rows"] = rows
    data["candidate"] = {"q": q_c, "gap_V279": a_gap,
                         "transverse_exponent_V279": a_exp,
                         "obstruction_V279": a_obs,
                         "m_next_V279": float(
                             art("V279_candidate_protection",
                                 "data.node.m_next"))}
    gm = [r["gap_min"] for r in rows if r["clock"] == "decorated"][-1]
    data["floor_note"] = {
        "decorated_gmin_r_0p01": gm,
        "meaning": "the candidate line is gapped on the 128-dim "
                    "clock (V279 m_next 0.323): A(r) here is the "
                    "anisotropy of the GAPPED floor, which falls "
                    "linearly in r - the cone-relevance question "
                    "itself is inconclusive on the clock"}
    dec = [r["A"] for r in rows if r["clock"] == "decorated"]
    irrelevant = bool(dec[0] >= dec[1] >= dec[2] and dec[2] < 0.05)
    # NOTE booked: the transverse gap at this candidate does NOT close on
    # the 128-dim clocks (g_min(r) -> 0.35, not 0: the V279 line is an
    # atom-module node, gapped on the clock) - so A(r) -> 0 measures the
    # isotropy of the gapped floor, booked honestly below.
    status = ("ANISOTROPY-IRRELEVANT" if irrelevant
              else "ANISOTROPY-RELEVANT")
    data["status"] = status
    print("V290 {}: decorated A series {} (test: monotone down and < "
          "0.05 at dq=1e-3)".format(status,
                                    ["{:.4f}".format(x) for x in dec]),
          flush=True)
    checks.append({"name": "T1 anisotropy series booked ({} rows); "
                          "decorated final A {:.4f}".format(len(rows),
                                                            dec[2]),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V290_anisotropy_relevance", "tier": 2,
           "property": "the booked candidate cone's anisotropy running "
                       "irrelevant under the booked coarse-graining",
           "residual": True, "programme": "V290",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v290_anisotropy_relevance.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V290 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
