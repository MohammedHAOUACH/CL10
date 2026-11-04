#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V378 - the crystalline dilemma, witnessed (booked BEFORE the run in
review/pre_registrations_root.md, "V378_crystalline_dilemma_witness").

V377 named the surviving structure of the module-replacement search:
the ladder objects exist only on crystallographic shift modules, and
clock-free closure only off them. This experiment puts that to a
direct witness test on the repository's only clock-free closure
engine - the 2I engine (cl10/supercell_2i): over the full 5-fold
class omega = pi*j/5 (j = 1..9) x beta in {G10, pi/5, 2pi/5, 1/5},
does ANY point carry the full conjunction closure(GAMMA, order <=
240) + gapless M + published R signature + internal mirror +
unitarity? DILEMMA-BROKEN (booked as the only module-side solution
carrier) / CRYSTALLINE-DILEMMA-WITNESSED.
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
from cl10.paths import exp  # noqa: E402

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


def internal_mirror(U):
    ph = np.sort(np.angle(np.linalg.eigvals(U)))
    return float(np.max(np.abs(ph + ph[::-1])))


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
    data = {"backend": "numpy (host; V337 ran cupy GPU1)"}

    # ---- T0b internal: the Spin(3) triple and the order-10 element ------
    Sx, Sy, Sz = S2.spin_triple()
    I32 = np.eye(32, dtype=np.complex128)
    sq = max(float(np.max(np.abs(X @ X - I32 / 4))) for X in (Sx, Sy, Sz))
    anti = max(float(np.max(np.abs(a @ b + b @ a)))
               for a, b in ((Sx, Sy), (Sy, Sz), (Sx, Sz)))
    U10 = S2.u_2i(2 * math.pi / 5)
    fifth = float(np.max(np.abs(np.linalg.matrix_power(U10, 5) + I32)))
    a_t0b = sq < 1e-12 and anti < 1e-12 and fifth < 1e-10
    chk("Spin(3) triple from the three spatial planes: S_a^2 = I/4 "
        "(max %.1e), anticommutators (%.1e), U_2I(2pi/5)^5 + I = %.1e "
        "(order 10) - the V337 T0b structure re-certified host-side: %s"
        % (sq, anti, fifth, a_t0b), a_t0b)

    # ---- T0a published anchor through THIS machinery ---------------------
    Upub = np.asarray(E.asnumpy(SC.internal(SC.G10)), dtype=np.complex128)
    gaps_pub = min(float(np.min(np.abs(phases_np(
        S2.step(q, Upub, SC.G10))))) for q in SC.M_POINTS)
    deg_pub = []
    for q in SC.R_POINTS[:3]:
        deg_pub += cluster_sizes(phases_np(S2.step(q, Upub, SC.G10)))
    R_ref = sorted(set(deg_pub))
    mir_pub = internal_mirror(Upub)
    v337 = json.load(open(exp("V337_2i_engine", "results",
                              "v337_2i_engine.json")))["data"]
    R_art = sorted(v337["T0a"]["pub_R_cluster_sizes"])
    a_t0a = (gaps_pub < 1e-9 and mir_pub < 1e-12 and R_ref == R_art)
    chk("published step re-measured through the V378 machinery: M-gap "
        "%.1e < 1e-9, mirror %.1e < 1e-12, R cluster signature %s == the "
        "V337 T0a field %s (backend reference match): %s"
        % (gaps_pub, mir_pub, R_ref, R_art, a_t0a), a_t0a)

    # ---- T1/T2 the 36-point class ----------------------------------------
    OMEGAS = [math.pi * j / 5 for j in range(1, 10)]
    BETAS = [SC.G10, math.pi / 5, 2 * math.pi / 5, 0.2]
    rows = []
    for om in OMEGAS:
        U = S2.u_2i(om)
        mir = internal_mirror(U)
        for be in BETAS:
            unit = max(unitarity(S2.step(q, U, be))
                       for q in (SC.GAMMA, SC.M_POINTS[0], SC.R_POINTS[0]))
            gapM = min(float(np.min(np.abs(phases_np(
                S2.step(q, U, be))))) for q in SC.M_POINTS)
            deg = []
            for q in SC.R_POINTS[:3]:
                deg += cluster_sizes(phases_np(S2.step(q, U, be)))
            r_sig = sorted(set(deg))
            k, cerr = closure_k(S2.step(SC.GAMMA, U, be))
            ki, _ = closure_k(U)
            holds = (unit <= 1e-11 and gapM <= 1e-9 and mir <= 1e-10
                     and r_sig == R_ref)
            rows.append({"omega_over_pi": om / math.pi, "beta": be,
                         "unitarity": unit, "gapM": gapM, "mirror": mir,
                         "R_holds": r_sig == R_ref, "closure_k": k,
                         "closure_internal_k": ki,
                         "closure_err": cerr, "holds": holds})
    n_hold = sum(1 for r in rows if r["holds"])
    n_closed = sum(1 for r in rows if r["closure_k"] is not None)
    n_int = sum(1 for r in rows if r["closure_internal_k"] is not None)
    a_closed = n_int == len(rows) and n_closed == 27
    chk("clock side: the 2I INTERNAL layer closes at every class point "
        "(%d/%d, order <= %d) - closure is clock-free in the internal "
        "layer; the FULL step closes at %d/36 = every rational-multiple-"
        "of-pi beta and never at beta = 1/5 rad (transport phase 0.2 "
        "rad/stage is incommensurate with 2pi): %s"
        % (n_int, len(rows),
           max(r["closure_internal_k"] for r in rows), n_closed,
           a_closed), a_closed)

    # overlap reproduction against the V337 cupy rows (backend cross-check)
    # (the V337 artifact stores the omega/pi and beta/pi RATIOS)
    worst = 0.0
    art = v337["grid"]
    n_ov = 0
    for r in rows:
        for p in art:
            if (abs(r["omega_over_pi"] - p["omega"]) < 1e-9
                    and abs(r["beta"] / math.pi - p["beta"]) < 1e-9):
                n_ov += 1
                worst = max(worst, abs(r["gapM"] - p["gapM"]))
    a_repro = worst <= 1e-8
    chk("the %d V337 grid points inside the booked class reproduce "
        "host-side vs the V337 cupy artifact at max |d gapM| %.1e <= "
        "1e-8: %s" % (n_ov, worst, a_repro), a_repro)

    rhold_gaps = [r["gapM"] for r in rows if r["R_holds"]]
    min_rhold = min(rhold_gaps)
    gapzeros = [(round(r["omega_over_pi"], 4), r["beta"]) for r in rows
                if r["gapM"] <= 1e-9]
    a_floor = min_rhold >= 0.1047 - 1e-6
    chk("the R-holding points of the class (beta = G10 = the published "
        "class) all carry gapM >= 0.1047 - measured floor %.6f (booked "
        "from the V337 row); gapM = 0 fires ONLY where R fails, at %s: %s"
        % (min_rhold, gapzeros, a_floor), a_floor)

    data["class_rows"] = rows
    data["n_holds"] = n_hold
    data["n_closed"] = n_closed
    data["R_reference"] = R_ref
    data["gapM_floor_R_holding"] = min_rhold
    data["gapM_zero_points"] = gapzeros

    ok = all(c[1] for c in checks)
    if n_hold == 0:
        status = ("CRYSTALLINE-DILEMMA-WITNESSED (0 of %d class points "
                  "carry the full conjunction: the internal 2I layer "
                  "closes at every class point at order <= %d - the "
                  "clock is born FREE there - yet the gapless M and the "
                  "published R multiplets NEVER co-occur: gapM = 0 only "
                  "where R fails, and R holds only at gapM >= %.6f; on "
                  "this engine the clock can live only by killing the "
                  "gapless spectrum or its R structure: no built module "
                  "carries both sides of the dilemma)"
                  % (len(rows),
                     max(r["closure_internal_k"] for r in rows),
                     min_rhold))
    else:
        br = [r for r in rows if r["holds"]]
        status = ("DILEMMA-BROKEN (%d point(s) carry the full "
                  "conjunction: %s)"
                  % (n_hold, json.dumps(br[:4])))
    out = {"id": "V378_crystalline_dilemma_witness",
           "tier": 2,
           "residual": min_rhold,
           "property": ("conjunction census of the 2I engine class - "
                        "clock-free closure x gapless M x published R x "
                        "mirror, the V377 named surviving structure put "
                        "to a witness test, tier 2, O8/O9"),
           "status": status,
           "pass": bool(ok),
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir,
                           "v378_crystalline_dilemma_witness.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V378 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
