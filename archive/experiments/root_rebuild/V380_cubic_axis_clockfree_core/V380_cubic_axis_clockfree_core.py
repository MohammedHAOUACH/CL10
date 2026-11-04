#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V380 - the cubic-axis clock-free core (booked BEFORE the run in
review/pre_registrations_root.md, "V380_cubic_axis_clockfree_core").

The one carrier class the tribunal missed: a finite-order twist of
the SPATIAL-plane Spin(3) about a CUBIC axis - clock-free closure
(fixed order), cubic symmetry preserved at R by construction, and
the crystallographic transport untouched. If any point of the
booked 10 carries gapless M + the published R signature + closure
+ unitarity, the crystalline dilemma is BROKEN and the clock
arbitration (F2) must be re-run on this carrier, because closure
no longer consumes the clock at any g.

Verdicts: CUBIC-CORE-STILL-BROKEN / DILEMMA-BROKEN.
"""
import json
import math
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from cl10 import engine as E  # noqa: E402
from cl10 import supercell as SC  # noqa: E402
from cl10 import supercell_2i as S2  # noqa: E402

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, ok))


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


def unitarity(M):
    return float(np.max(np.abs(M @ M.conj().T - np.eye(M.shape[0]))))


def closure_k(M, cap=240, tol=1e-8):
    Mk = np.eye(M.shape[0], dtype=np.complex128)
    for k in range(1, cap + 1):
        Mk = Mk @ M
        e = min(float(np.max(np.abs(Mk - np.eye(M.shape[0])))),
                float(np.max(np.abs(Mk + np.eye(M.shape[0])))))
        if e <= tol:
            return k, e
    return None, None


def main():
    data = {}
    r3 = math.sqrt(3.0)

    # ---- T0 self-tests ----------------------------------------------------
    Sx, Sy, Sz = S2.spin_triple()
    I32 = np.eye(32, dtype=np.complex128)
    sq = max(float(np.max(np.abs(X @ X - I32 / 4))) for X in (Sx, Sy, Sz))
    anti = max(float(np.max(np.abs(a @ b + b @ a)))
               for a, b in ((Sx, Sy), (Sy, Sz), (Sx, Sz)))
    a_t0b = sq < 1e-12 and anti < 1e-12
    chk("Spin(3) spatial triple re-certified: S_a^2 = I/4 (max %.1e), "
        "anticommutators (%.1e): %s" % (sq, anti, a_t0b), a_t0b)
    Upub = np.asarray(E.asnumpy(SC.internal(SC.G10)), dtype=np.complex128)
    gaps_pub = min(float(np.min(np.abs(phases_np(
        S2.step(q, Upub, SC.G10))))) for q in SC.M_POINTS)
    deg_pub = []
    for q in SC.R_POINTS[:3]:
        deg_pub += cluster_sizes(phases_np(S2.step(q, Upub, SC.G10)))
    R_ref = sorted(set(deg_pub))
    a_t0a = gaps_pub < 1e-9 and R_ref == [32, 48]
    chk("published anchor through this machinery: M-gap %.1e < 1e-9, "
        "R signature %s == [32, 48]: %s" % (gaps_pub, R_ref, a_t0a),
        a_t0a)

    # ---- T1/T2 the ten-point cubic class ----------------------------------
    class_pts = ([((1.0, 0.0, 0.0), math.pi), ((0.0, 1.0, 0.0), math.pi),
                  ((0.0, 0.0, 1.0), math.pi),
                  ((1.0, 0.0, 0.0), math.pi / 2),
                  ((0.0, 1.0, 0.0), math.pi / 2),
                  ((0.0, 0.0, 1.0), math.pi / 2),
                  ((1 / r3, 1 / r3, 1 / r3), 2 * math.pi / 3),
                  ((1 / r3, -1 / r3, 1 / r3), 2 * math.pi / 3)]
                 + [((1.0, 0.0, 0.0), 2 * math.pi / 3),
                    ((0.0, 0.0, 1.0), 2 * math.pi / 3)])
    rows = []
    n_broken = 0
    for n, om in class_pts:
        U = S2.u_2i(om, n)
        ki, kerr = closure_k(U, cap=24, tol=1e-8)
        unit = max(unitarity(S2.step(q, U, SC.G10))
                   for q in (SC.GAMMA, SC.M_POINTS[0], SC.R_POINTS[0]))
        gapM = min(float(np.min(np.abs(phases_np(
            S2.step(q, U, SC.G10))))) for q in SC.M_POINTS)
        deg = []
        for q in SC.R_POINTS[:3]:
            deg += cluster_sizes(phases_np(S2.step(q, U, SC.G10)))
        r_sig = sorted(set(deg))
        kf, ferr = closure_k(S2.step(SC.GAMMA, U, SC.G10))
        ph = np.sort(np.angle(np.linalg.eigvals(U)))
        mir = float(np.max(np.abs(ph + ph[::-1])))
        t1 = (unit <= 1e-11 and ki is not None and ki <= 12
              and kf is not None and kf <= 240 and mir <= 1e-10)
        t2 = t1 and gapM <= 1e-9 and r_sig == R_ref
        if t2:
            n_broken += 1
        rows.append({"axis": list(n), "omega": om, "closure_internal": ki,
                     "closure_step": kf, "unitarity": unit, "gapM": gapM,
                     "R": r_sig, "mirror": mir, "T1": t1, "T2": t2})
    data["class_rows"] = rows
    all_t1 = sum(1 for r in rows if r["T1"])
    chk("T1 gates over the class: %d/10 points carry unitarity + "
        "internal closure <= 12 + step closure <= 240 + mirror <= 1e-10 "
        "(per-point table in data): recorded" % all_t1, True)
    per_gap = ["%s w=%.4f gapM=%.4e R=%s"
               % ("/".join("%g" % v for v in r["axis"]), r["omega"],
                  r["gapM"], ",".join(map(str, r["R"]))) for r in rows]
    data["per_point"] = per_gap

    ok = all(c[1] for c in checks)
    if n_broken == 0:
        status = ("CUBIC-CORE-STILL-BROKEN (the cubic-axis clock-free "
                  "class - the one carrier the tribunal missed - carries "
                  "no point with gapless M AND the published R "
                  "signature AND closure together; %d/10 pass T1; the "
                  "geometric obstruction survives on cubic axes too; the "
                  "crystalline dilemma stands as the named surviving "
                  "structure)" % all_t1)
    else:
        br = [r for r in rows if r["T2"]]
        status = ("DILEMMA-BROKEN (%d carrier point(s): %s) - follow-up "
                  "obligation booked: F2 clock re-arbitration on this "
                  "carrier, closure consumes no clock at any g"
                  % (n_broken, json.dumps(br)))
    out = {"id": "V380_cubic_axis_clockfree_core",
           "tier": 2,
           "residual": min(r["gapM"] for r in rows
                           if r["R"] == R_ref and r["gapM"] > 1e-9)
           if any(r["R"] == R_ref and r["gapM"] > 1e-9 for r in rows)
           else min(r["gapM"] for r in rows),
           "property": ("the three-leg carrier search: finite-order "
                        "spatial-plane twists about cubic axes - "
                        "clock-free closure x gapless M x published R "
                        "x crystallographic ladder, tier 2, O8/O9"),
           "status": status,
           "pass": bool(ok),
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v380_cubic_axis_clockfree_core.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    for p in per_gap:
        log("  point " + p)
    log("[V380 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
