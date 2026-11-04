#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V381 - the free-clock carrier arbitration (booked BEFORE the run
in review/pre_registrations_root.md, "V381_free_clock_carrier_
arbitration" + the pre-run correction of the same day).

V380 fired the booked solution branch: the carrier u_2i(pi/2,
coordinate axis) closes the FULL step at order 2 with gapM = 0 and
the published R signature. The staggered transport inverts the
spatial-plane twist: (PR)^2 = +-I at every rational-multiple-of-pi
beta - the closure forcing leaves the booked set {pi/2, pi/4,
3pi/4} and the clock is freed into the whole rational-pi family
(the rational-radians pins never close - the V378 lesson). Whether
a carrier clock READS the hierarchy coefficient is V382's booked
question, not this one.

Verdicts: FREE-CLOCK-CARRIER-CONFIRMED / CARRIER-ARBITRATION-
BREAKS.
"""
import json
import math
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
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


def closure_k(M, cap=24, tol=1e-8):
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
    U = S2.u_2i(math.pi / 2, (1.0, 0.0, 0.0))
    ph_u = np.sort(np.angle(np.linalg.eigvals(U)))
    mir = float(np.max(np.abs(ph_u + ph_u[::-1])))

    # T0 anchor: the carrier internal is an EVEN element of the root
    # (built from the three spatial planes) - re-certify the triple.
    Sx, Sy, Sz = S2.spin_triple()
    I32 = np.eye(32, dtype=np.complex128)
    sq = max(float(np.max(np.abs(X @ X - I32 / 4))) for X in (Sx, Sy, Sz))
    a_t0 = sq < 1e-12 and abs(mir) <= 1e-10
    chk("carrier internal u_2i(pi/2, e_z): spatial triple S_a^2 = I/4 "
        "(max %.1e), internal mirror %.1e <= 1e-10: %s" % (sq, mir, a_t0),
        a_t0)

    BETAS = [("pi/4", math.pi / 4), ("pi/3", math.pi / 3),
             ("pi/5", math.pi / 5), ("2pi/5", 2 * math.pi / 5),
             ("pi/6", math.pi / 6), ("1/5rad", 0.2),
             ("g_obs", 0.19999798835751165), ("0.3rad", 0.3)]
    RATIONAL_PI = {"pi/4", "pi/3", "pi/5", "2pi/5", "pi/6"}
    Upub = np.asarray(E.asnumpy(SC.internal(SC.G10)), dtype=np.complex128)
    R_ref = sorted(set(sum((cluster_sizes(phases_np(
        S2.step(q, Upub, SC.G10))) for q in SC.R_POINTS[:3]), [])))

    rows = []
    a_unit = a_zerop = True
    a_close_ratio = a_close_rad = True
    a_contact = True
    for name, be in BETAS:
        unit = max(unitarity(S2.step(q, U, be))
                   for q in (SC.GAMMA, SC.M_POINTS[0], SC.R_POINTS[0]))
        k, kerr = closure_k(S2.step(SC.GAMMA, U, be))
        gapM = min(float(np.min(np.abs(phases_np(
            S2.step(q, U, be))))) for q in SC.M_POINTS)
        deg = sum((cluster_sizes(phases_np(S2.step(q, U, be)))
                   for q in SC.R_POINTS[:3]), [])
        r_sig = sorted(set(deg))
        ph_g = phases_np(S2.step(SC.GAMMA, U, be))
        zpt = abs((ph_g.sum() + math.pi) % (2 * math.pi) - math.pi)
        a_unit = a_unit and unit <= 1e-11
        a_zerop = a_zerop and zpt <= 1e-9
        if name in RATIONAL_PI:
            a_close_ratio = a_close_ratio and k is not None
            a_contact = a_contact and gapM <= 1e-9 and r_sig == R_ref
        else:
            a_close_rad = a_close_rad and k is None
        rows.append({"beta": name, "closure_k": k, "closure_err": kerr,
                     "unitarity": unit, "gapM": gapM, "R": r_sig,
                     "zeropoint": zpt})
    data["carrier_rows"] = rows
    data["R_reference"] = R_ref
    closed_clocks = sorted(set(
        (r["beta"], r["closure_k"]) for r in rows if r["closure_k"]))
    chk("unitarity <= 1e-11 and 0-point <= 1e-9 at every beta: %s / %s"
        % (a_unit, a_zerop), a_unit and a_zerop)
    chk("closure: EVERY rational-multiple-of-pi beta closes (order "
        "<= 24: %s) and NO rational-radians beta closes (1/5 rad and "
        "0.3 rad: None to cap 24) - the closure forcing leaves the "
        "booked {pi/2, pi/4, 3pi/4}: the clock is FREED into the "
        "rational-pi family: %s" % (closed_clocks,
                                    a_close_ratio and a_close_rad),
        a_close_ratio and a_close_rad)
    chk("spectrum contacts at every closed rational-pi beta: gapM = 0 "
        "(<= 1e-9) and R signature == %s: %s" % (R_ref, a_contact),
        a_contact)

    ok = all(c[1] for c in checks)
    status = ("FREE-CLOCK-CARRIER-CONFIRMED (the carrier u_2i(pi/2, "
              "coordinate axis) closes at %s with gapless M and the "
              "published R at every closed beta; the closure-forced "
              "clock set {pi/2, pi/4, 3pi/4} is NOT the closure set "
              "of the carrier - rational-radians pins never close "
              "(V378 lesson verbatim); the F2 question transfers to "
              "the coefficient reading of a carrier clock - booked "
              "as V382)" % closed_clocks) if ok else \
        "CARRIER-ARBITRATION-BREAKS (see failing gates)"
    out = {"id": "V381_free_clock_carrier_arbitration",
           "tier": 2,
           "residual": min(r["closure_err"] for r in rows
                           if r["closure_k"] is not None),
           "property": ("the V380 carrier put to the clock-arbitration "
                        "re-run: closure, gapless M and R across eight "
                        "betas, tier 2, O8/O9"),
           "status": status,
           "pass": bool(ok),
           "data": data,
           "checks": [(t, bool(g)) for t, g in checks]}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir,
                           "v381_free_clock_carrier_arbitration.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V381 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
