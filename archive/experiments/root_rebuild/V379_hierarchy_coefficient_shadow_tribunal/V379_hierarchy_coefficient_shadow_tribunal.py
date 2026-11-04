#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V379 - the hierarchy-coefficient shadow tribunal (booked BEFORE the
run in review/pre_registrations_root.md, "V379_hierarchy_coefficient_
shadow_tribunal").

The sharpest hand-found contact after the module tribunal: the
observed hierarchy coefficient 2.076564111173832 sits 3.8e-05 from
3*exp(-1/e) - inside the booked 1e-04 band, zero fitted parameters.
Tribunal: (T1) the trap reading, recorded whatever follows;
(T2) the required-exponent census against every web integer;
(T3) the web-membership tribunal on e-type transcendentals (the V340
pslq pattern); (T4) the R near-miss map on the four 2I gapless
points - the shape of the missing third leg.

Verdicts: E-SHADOW-COEFFICIENT-NO-WEB-CONTACT / E-WEB-CONTACT-FOUND.
"""
import json
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)

import math  # noqa: E402
import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402

from cl10 import engine as E  # noqa: E402
from cl10 import supercell as SC  # noqa: E402
from cl10 import supercell_2i as S2  # noqa: E402
from cl10 import u10pred  # noqa: E402
from cl10.paths import exp  # noqa: E402

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, ok))


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def phases_np(M):
    return np.sort(np.angle(np.linalg.eigvals(M)))


def cluster_sizes(ph, tol=1e-9):
    out = []
    for p in sorted(ph):
        if out and abs(p - out[-1][0]) <= tol:
            out[-1][1] += 1
            out[-1][0] = p
        else:
            out.append([p, 1])
    return [c[1] for c in out]


def main():
    mp.mp.dps = 40
    data = {}

    # ---- T0 web anchors --------------------------------------------------
    v233 = json_art("V233_clock_arbitration")["data"]["readings"]
    obs = mp.mpf(repr(v233["observed"]))
    unc = mp.mpf(repr(v233["obs_uncertainty"]))
    v312d = json_art("V312_quaternionic_loop_ladder")["data"]
    runs = v312d["ladder_runs"]["0"]
    F_reads = {k: mp.mpf(repr(v["F_cross"])) for k, v in runs.items()}
    rho = mp.mpf(repr(json_art("V204_closed_ladder_limits_audit")
                             ["data"]["machinery"]["rho"]))
    phi = (1 + mp.sqrt(5)) / 2
    a_t0 = (abs(obs - mp.mpf("2.076564111173832")) == 0
            and abs(rho - mp.mpf("0.9780302704295731")) == 0)
    chk("anchors from artifacts: observed %r +- %s (V233 field), the "
        "three V312 F_cross reads, rho = 0.9780302704295731 (V312 "
        "field): %s" % (v233["observed"], v233["obs_uncertainty"], a_t0),
        a_t0)

    # ---- T1 the trap reading ---------------------------------------------
    A = 3 * mp.e ** (-1 / mp.e)
    dev = A - obs
    band = abs(dev) < unc
    chk("TRAP: 3*exp(-1/e) = %s ; |A - observed| = %.6e < the booked "
        "band 1.0e-04 - the hit is a REAL numerical fact (rel %.2e), "
        "recorded regardless of the tribunal: %s"
        % (mp.nstr(A, 17), abs(dev), abs(dev) / obs, band), band)
    data["three_exp_minus_one_over_e"] = mp.nstr(A, 25)
    data["deviation_to_observed"] = mp.nstr(dev, 8)

    # ---- T2 required-exponent census --------------------------------------
    ints = {30: "E8 Coxeter", 32: "spinor levels", 36: "module DIM",
            45: "adjoint", 64: "Casimir string", 208: "crossing n*",
            10: "DM rung", 17: "nu rung", 1: "I_1", 6: "I_3",
            23: "I_4?", 141: "I_6"}
    nreqs = {}
    ok_n = True
    for k, F in F_reads.items():
        if F == 0:
            continue
        nreq = mp.log(obs / 3) / mp.log(16 * mp.pi * F)
        dmin = min(abs(nreq - i) for i in ints)
        nreqs[k] = {"n_req": mp.nstr(nreq, 8),
                    "dist_to_nearest_web_int": mp.nstr(dmin, 4)}
        ok_n = ok_n and dmin > mp.mpf("0.05")
    chk("required-exponent census: n_req(observed, F_cross reads) = %s "
        "- no web integer (%s) sits within 0.05 of any required "
        "exponent: %s" % (json.dumps(nreqs), ",".join(map(str, sorted(
            ints))), ok_n), ok_n)
    data["n_req_by_F_read"] = nreqs

    # ---- T3 the web-membership tribunal on e ------------------------------
    X = mp.e ** (-1 / mp.e)
    Xinv = mp.e ** (1 / mp.e)
    loop = mp.mpf(repr(float(u10pred.LOOP)))
    basis = {"16pi": 16 * mp.pi, "1/pi": 1 / mp.pi, "phi": phi,
             "loop": loop, "rho": rho, "e": mp.e, "1/e": 1 / mp.e}
    for k, F in F_reads.items():
        basis["F_" + k] = F
    for w in ints:
        basis["int%d" % w] = mp.mpf(w)
    hits = []

    def relation(vec):
        """pslq plus the 40-dp verification: only EXACT relations fire."""
        c = mp.pslq(vec, maxcoeff=200)
        if c is None:
            return None
        # the target itself must participate (a relation among the
        # basis elements alone - e.g. -w*1 + 1*w for integer w - is
        # not a membership relation)
        if c[1] == 0 and (len(c) < 4 or c[3] == 0):
            return None
        if len(c) == 4:
            # for an INTEGER basis w the columns X and X*w are
            # collinear: the degenerate relation -w.X + 1.(Xw) has
            # both grouped coefficients (a + c w) and (b + d w)
            # vanish identically - membership requires a non-zero
            # grouped constant term (a + c w) != 0
            if c[0] + c[2] * w == 0:
                return None
        if abs(sum(ci * xi for ci, xi in zip(c, vec))) < mp.mpf("1e-25"):
            return [int(v) for v in c]
        return None

    for name, w in basis.items():
        for target, tname in ((X, "exp(-1/e)"), (Xinv, "exp(1/e)")):
            r1 = relation([mp.mpf(1), target, w])
            r2 = relation([mp.mpf(1), target, w, target * w])
            if r1 is not None or r2 is not None:
                hits.append((tname, name, r1, r2))
    a_t3 = len(hits) == 0
    chk("web-membership tribunal: pslq([1,X,w]) and pslq([1,X,w,Xw]) "
        "at maxcoeff 200 over BOTH e-targets and the full certified "
        "web basis (%d numbers): %d relations fired - e-type "
        "transcendentals are OFF the web: %s"
        % (len(basis), len(hits), a_t3), a_t3)
    data["web_basis_size"] = len(basis)
    data["web_hits"] = [[h[0], h[1]] for h in hits]

    # ---- T4 the R near-miss map -------------------------------------------
    pts = [(0.4, 2 * math.pi / 5), (0.8, math.pi / 5),
           (1.2, math.pi / 5), (1.6, 2 * math.pi / 5)]
    Upub = np.asarray(E.asnumpy(SC.internal(SC.G10)), dtype=np.complex128)
    deg_pub = []
    for q in SC.R_POINTS[:3]:
        deg_pub += cluster_sizes(phases_np(S2.step(q, Upub, SC.G10)))
    R_pub = sorted(set(deg_pub))
    sigs = []
    for om_pi, be in pts:
        U = S2.u_2i(om_pi * math.pi)
        deg = []
        for q in SC.R_POINTS[:3]:
            deg += cluster_sizes(phases_np(S2.step(q, U, be)))
        sigs.append((om_pi, sorted(set(deg))))
    a_t4 = R_pub == [32, 48] and all(s == [16, 48] for _, s in sigs)
    chk("R near-miss map: published signature %s; at all four 2I "
        "gapless points the signature is [16, 48] - the 48-cluster "
        "SURVIVES, the 32-cluster splits exactly dyadically: %s"
        % (R_pub, a_t4), a_t4)
    data["R_published"] = R_pub
    data["R_at_gapless_2I_points"] = [[p, s] for p, s in sigs]

    ok = all(c[1] for c in checks)
    if not hits:
        status = ("E-SHADOW-COEFFICIENT-NO-WEB-CONTACT (3*exp(-1/e) "
                  "reproduces the observed hierarchy coefficient to "
                  "%.6e - inside the observation band, zero fitted "
                  "parameters - but neither exp(1/e) nor exp(-1/e) "
                  "sits on ANY bounded integer relation to the "
                  "certified web (pslq, maxcoeff 200, %d basis "
                  "numbers), and no web integer reads the required "
                  "exponent: the cleanest number hit in the theory's "
                  "history is a SHADOW with no mechanism; F2 stands. "
                  "The shape constraint on any module-side solution "
                  "from the near-miss map: the clock-free core must "
                  "preserve the 32-cluster at R - the 48 leg already "
                  "survives, the M leg already survives, only the "
                  "dyadic 16+16 split of the 32 blocks the "
                  "three-leg carrier)"
                  % (abs(dev), len(basis)))
    else:
        status = ("E-WEB-CONTACT-FOUND (%s)"
                  % json.dumps([[h[0], h[1]] for h in hits]))
    out = {"id": "V379_hierarchy_coefficient_shadow_tribunal",
           "tier": 2,
           "residual": float(abs(dev)),
           "property": ("the 3*exp(-1/e) hit on the hierarchy "
                        "coefficient put to the web-membership "
                        "tribunal, required-exponent census, and the "
                        "R near-miss shape map, tier 2, O8"),
           "status": status,
           "pass": bool(ok),
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir,
                           "v379_hierarchy_coefficient_shadow_tribunal"
                           ".json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V379 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
