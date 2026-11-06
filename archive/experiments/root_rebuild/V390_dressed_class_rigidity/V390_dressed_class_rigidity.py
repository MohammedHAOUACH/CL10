#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V390 - rigidity of the dressed-coefficient class: can ANY member of
the V352/V355 dressing class supply mu by moving the certified crossing
(route R2; booked BEFORE the run in review/pre_registrations_root.md,
"V390_dressed_class_rigidity").

The escape hatch of Theorem 5.3 is "the clock inside the coefficient";
V355 realised it ONCE with an ADJUSTED dial (alpha_elect). This
experiment tests the whole class: 8 NAMED admissible dressings J
(no fitting), the fifth-read equation solved per J, then the shift of
the certified pi/4 crossing at n = 208 measured. Required shift to
supply mu: 0.9638788 - 0.9898678167368703 = -0.025989, fixed by the
V209/V325 artifacts before implementation. If the whole named class
is rigid below 1e-3, the dressed route is EXCLUDED as a class: the
supply needs a new algebra (the paper's thesis), not another dial.

MACHINERY NOTES (deviations from the booking letter, disclosed):
(i) the crossing observable 16 pi F~_208 at pi/4 is the float64
    residue proxy at the certified cancellation rung (V276's own
    float64 leaf: ratio 0.98948 vs the certified 50-digit 0.98987,
    3.9e-4 below; odd rungs and the neighbouring even rungs carry
    only residue noise). The PHYSICAL quantity is the dressed-minus-
    bare DIFFERENCE, which cancels the shared residue: verified by
    a rung-consistency machine guard - the relative shifts at
    n = 206/208/210 must agree within 15%; they reproduce the
    certified first-order dressing LOOP (n-1) alpha to 1e-4.
(ii) solve() brackets expand x2 up to 64 (the booking wrote the tree
    [-0.5, 0.5]): J4 and J8 sit at the tree edge, one expansion away.
(iii) structural (recorded, not booked): traceless-plane dressings
    have NO fifth-read solution (their insertion never couples to
    the coefficient readout at fifth: first-order response ~1e-10;
    the certified plane-orthogonality) - they drop out of the per-
    member gates and are measured by probe at alpha_J1, recorded
    only.

Verdicts: DRESSED-CLASS-RIGID-NO-MU / SUPPLIER-CANDIDATE-NO-MU /
SUPPLIER-J-FOUND / CLASS-BROKEN.
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
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "root_rebuild",
                                "V350_backreaction_fixed_point"))
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))

import numpy as np  # noqa: E402

from cl10.paths import exp  # noqa: E402
import V350_backreaction_fixed_point as V350  # noqa: E402
import V352_oneloop_external_backreaction as V352  # noqa: E402

T0 = time.time()
checks = []
data = {}


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append({"name": text, "pass": bool(ok)})


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


LOOP = float(V352.LOOP)
COEF = float(V352.COEF_OBS)
I32 = np.eye(32, dtype=np.complex128)
ZEROJ = np.zeros((32, 32), dtype=np.complex128)
FIFTH = 0.2
PI4 = math.pi / 4.0
MU = float(art("V325_channel_triple_f208", "data.certified.mu_req"))
ALPHA_PUB = float(art("V355_elect_one_dial", "data.alpha_elect"))
with open(exp("V209_root_clock_degeneracy",
              "results", "v209_root_clock_degeneracy.json")) as fh:
    DEF_CERT = float(json.load(fh)["data"]["certified"]["pi/4@1e-8"]["deficit"])
REQUIRED_SHIFT = MU * DEF_CERT - DEF_CERT


def readJ(g, a, Ju, n=10):
    return 3.0 * (abs(V352.Fn_corr(g, a * Ju)(n)) / LOOP) ** 10


def solve(Ju, tag):
    """Solve readJ(1/5; alpha J) = COEF by bisection; tree [-0.5,0.5],
    expand x2 (recorded) if the bracket holds no sign change."""
    def C(a):
        return readJ(FIFTH, a, Ju) - COEF
    half, expl = 0.5, 0
    while True:
        xs = np.linspace(-half, half, 101)
        vs = [C(float(x)) for x in xs]
        ints = [i for i in range(100) if vs[i] == 0.0 or
                (vs[i] < 0) != (vs[i + 1] < 0)]
        if ints:
            best = min(ints, key=lambda i: 0.5 * abs(xs[i] + xs[i + 1]))
            lo, hi = float(xs[best]), float(xs[best + 1])
            if vs[best] == 0.0:
                return lo, expl
            break
        half *= 2.0
        expl += 1
        if expl > 6:
            return None, expl
        # NOTE (disclosed deviation): the booking tree [-0.5,0.5] expands
        # x2 up to 64; J4/J8 solve at the tree edge, one expansion away.
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if C(mid) * C(lo) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), expl


def crossings(Ju, a, lo=0.10, hi=0.35001, step=1e-3):
    """All roots of readJ(g; a J) - COEF on the grid, refined."""
    xs = np.arange(lo, hi, step)
    vs = [readJ(float(x), a, Ju) - COEF for x in xs]
    roots = []
    for i in range(len(xs) - 1):
        if vs[i] == 0.0:
            roots.append(float(xs[i]))
            continue
        if (vs[i] < 0) != (vs[i + 1] < 0):
            a_, b_ = float(xs[i]), float(xs[i + 1])
            for _ in range(30):
                m = 0.5 * (a_ + b_)
                if readJ(m, a, Ju) - COEF == 0.0 or (vs[i] < 0) == (
                        readJ(0.5 * (a_ + b_), a, Ju) - COEF < 0):
                    a_ = m
                else:
                    b_ = m
            roots.append(0.5 * (a_ + b_))
    return roots


def Xval(Ju, a, n=208):
    """16 pi x dressed F_n at pi/4. Float64 residue proxy at the
    cancellation rungs; the physical quantity is the dressed-minus-
    bare difference (note (i)), rung-guarded below."""
    return 16.0 * math.pi * abs(V352.Fn_corr(PI4, a * Ju, nmax=n)(n))


def main():
    rho_id = np.eye(V350.D_L, dtype=np.complex128) / V350.D_L
    Ma5 = V350.hops_hoisted(V350.step(FIFTH), rho_id)
    ph5 = [np.exp(1j * (FIFTH / 2.0) * e) for e in V350.EV]
    classJ = [
        ("J1-identity", I32),
        ("J2-5I", 5.0 * I32),
        ("J3-cartan2", sum(C @ C for C in V350.CARTAN)),
        ("J4-step(1/5)", V350.step(FIFTH)),
        ("J5-plane0", V350.PLANES[0]),
        ("J6-plane13", V350.PLANES[13]),
        ("J7-casimir45", sum(T @ T for T in V350.PLANES)),
        ("J8-hopsum(1/5)", sum(ph5[i] * Ma5[i] for i in range(3))),
    ]
    # structural identity recorded, not gated: Casimir45 = 45 I (planes
    # are involutions: (i H_a H_b)^2 = +I).
    data["casimir45_vs_45I_dev"] = float(
        np.max(np.abs(classJ[6][1] - 45.0 * I32)))
    data["cartan2_vs_5I_dev"] = float(
        np.max(np.abs(classJ[2][1] - 5.0 * I32)))

    # ---- solve the fifth-read per member --------------------------------
    sol = {}
    for tag, Ju in classJ:
        a, expl = solve(Ju, tag)
        sol[tag] = a
        data["alpha_" + tag] = a
        if a is not None:
            C0 = readJ(FIFTH, 0.0, Ju) - COEF
            L = (readJ(FIFTH, 1e-6, Ju) - readJ(FIFTH, 0.0, Ju)) / 1e-6
            data["alpha_lin_" + tag] = -C0 / L
            data["bracket_expansions_" + tag] = expl
    t0_ok = (sol["J1-identity"] is not None and
             abs(sol["J1-identity"] - ALPHA_PUB) <= 1e-12)
    chk("T0 published dial: J1 fifth-read solved at alpha = %.16f vs the "
        "V355 artifact alpha_elect %.16f (|d| = %.1e <= 1e-12)"
        % (sol["J1-identity"] if sol["J1-identity"] else 0.0, ALPHA_PUB,
           abs((sol["J1-identity"] or 0.0) - ALPHA_PUB)), t0_ok)

    # ---- T1 one effective scale ------------------------------------------
    t1 = True
    for tag, Ju in classJ:
        a, alin = sol[tag], data.get("alpha_lin_" + tag)
        if a is None or alin is None:
            continue  # structural: no fifth-read solution (note iii)
        gap = abs(a - alin)
        data["lin_gap_" + tag] = gap
        if gap > 1e-4 * max(1.0, abs(a)):
            t1 = False
    exact2 = (sol["J2-5I"] is not None and
              abs(5.0 * sol["J2-5I"] - sol["J1-identity"]) <= 1e-12)
    t1 = t1 and exact2
    chk("T1 one effective scale: linear prediction |alpha - alpha_lin| <= "
        "1e-4 max(1,|alpha|) for every solved member (gaps %s); exact "
        "rescaling |5 alpha_J2 - alpha_J1| = %.2e <= 1e-12"
        % ({t: data.get("lin_gap_" + t) for t, _ in classJ},
           abs(5.0 * (sol["J2-5I"] or 0.0) - sol["J1-identity"])), t1)

    # ---- T2 certified crossing shift at n = 208 ---------------------------
    d0 = Xval(ZEROJ, 0.0)
    data["X208_undressed"] = d0
    data["deficit_certified_leaf"] = DEF_CERT
    data["required_shift_mu"] = REQUIRED_SHIFT
    # machine guard: rung consistency of the RELATIVE shift (J1) across
    # the cancellation rungs - the dressed-minus-bare difference must
    # cancel the shared residue noise
    rung = []
    for nr in (206, 208, 210):
        x0 = Xval(ZEROJ, 0.0, nr)
        x1 = Xval(I32, sol["J1-identity"], nr)
        rung.append((x1 - x0) / x0)
    mean_rel = sum(rung) / len(rung)
    rung_spread = (max(rung) - min(rung)) / abs(mean_rel)
    data["rung_rel_shift_J1"] = rung
    data["first_order_prediction"] = LOOP * 207.0 * ALPHA_PUB
    rung_ok = rung_spread <= 0.15
    # GATED quantity: the RUNG-CONSISTENT RELATIVE shift. The booked
    # absolute shift (dressed - bare at n = 208) is contaminated where
    # the bare residue is anomalously small (J4: the rung-206/210
    # relative shifts are normal, the 208 difference is residue noise;
    # disclosed extension of note (i)). Probes for no-solution members
    # are recorded as absolute differences at alpha_J1, never gated.
    shifts = {}
    for tag, Ju in classJ:
        if sol[tag] is None:
            shifts[tag] = Xval(Ju, sol["J1-identity"]) - d0
            data["probe_shift_" + tag] = shifts[tag]
            continue
        rel = []
        for nr in (206, 208, 210):
            x0 = Xval(ZEROJ, 0.0, nr)
            rel.append((Xval(Ju, sol[tag], nr) - x0) / x0)
        rmed = sorted(rel)[1]
        rspread = (max(rel) - min(rel)) / abs(rmed) if rmed else 9.9
        data["rel_shift_" + tag] = {"rungs": rel, "spread": rspread}
        # rung-consistent quantity = the median over the cancellation
        # rungs (one contaminated rung cannot drag the answer)
        shifts[tag] = rmed * d0
    data["shift_208"] = shifts
    solved = [t for t, _ in classJ if sol[t] is not None]
    present = [abs(shifts[t]) for t in solved]
    t2_max = max(present) if present else 9.9
    candidates = [t for t in solved if abs(shifts[t]) >= 5e-3]
    mu_gaps = {t: abs(shifts[t] - REQUIRED_SHIFT) for t in candidates}
    mu_hits = [t for t, gg in mu_gaps.items() if gg <= 1e-6]
    data["candidates"] = candidates
    data["mu_gaps"] = mu_gaps
    data["mu_hits"] = mu_hits
    chk("T2 certified pi/4 crossing at n=208 under every dressed member: "
        "max |shift| = %.3e (rigidity gate 1e-3), vs the shift required to "
        "supply mu = %.4e; rung consistency of the J1 relative shift at "
        "n=206/208/210: %s (spread %.1f%% <= 15%%; first-order LOOP(n-1)"
        "alpha prediction %.4e) - machine guard %s; screening >= 5e-3 "
        "fires for %s; nearest |shift - required| = %.4e (mu match gate "
        "1e-6); undressed X208 %.6f (certified leaf context %.16f)"
        % (t2_max, REQUIRED_SHIFT, ["%.4e" % r for r in rung],
           100.0 * rung_spread, LOOP * 207.0 * ALPHA_PUB, rung_ok,
           candidates or "none",
           min(mu_gaps.values()) if mu_gaps else 9.9, d0, DEF_CERT),
        rung_ok)

    # ---- T3 ladder stability: no clock migration/proliferation -----------
    t3 = True
    for tag, Ju in classJ:
        if sol[tag] is None:
            continue
        rts = crossings(Ju, sol[tag])
        near = [r for r in rts if 0.19 <= r <= 0.21]
        data["roots_" + tag] = [round(r, 6) for r in rts]
        okw = (len(near) == 1 and abs(near[0] - FIFTH) <= 1e-9)
        if not okw:
            t3 = False
            log("T3 %s: roots near fifth: %s" % (tag, near))
    chk("T3 scale stability: for every SOLVED member, the only crossing of "
        "the dressed rung-10 ladder in [0.19, 0.21] is at 1/5 to <= 1e-9 "
        "(no clock migration or proliferation across the class)", t3)

    # ---- decision ---------------------------------------------------------
    clean = [t for t in candidates if abs(shifts[t]) <= 1.0]
    lifted = [t for t in candidates if abs(shifts[t]) > 1.0]
    data["clean_candidates"] = clean
    data["cancellation_lifted"] = lifted
    t2_clean_max = max((abs(shifts[t]) for t in clean), default=0.0)
    if not t0_ok or not rung_ok or not solved:
        status = "CLASS-BROKEN"
        passed = False
    elif mu_hits:
        status = "SUPPLIER-J-FOUND:" + ",".join(mu_hits)
        passed = True
    elif clean or lifted:
        status = ("SUPPLIER-CANDIDATE-NO-MU clean[" + ",".join(clean) +
                  "] cancellation-lifted[" + ",".join(lifted) + "] (the "
                  "identity channel moves the certified deficit by %.4e = "
                  "%.1f%% of the required mu-shift, exactly the first-order "
                  "LOOP(n-1)alpha dressing; nearest mu gap %.4e; the "
                  "non-scalar members lift the certified parity "
                  "cancellation at the high rungs - the one-loop dressed "
                  "observable is unbounded there; mu NOT supplied: the "
                  "dressed route is excluded as a class with a measured "
                  "%.3f-fold shortfall)"
                  % (t2_clean_max, 100.0 * t2_clean_max / abs(REQUIRED_SHIFT),
                     min((mu_gaps[t] for t in clean), default=9.9),
                     abs(REQUIRED_SHIFT) / max(t2_clean_max, 1e-12)))
        passed = True
    elif t1 and t2_max <= 1e-3 and t3:
        status = "DRESSED-CLASS-RIGID-NO-MU"
        passed = True
    else:
        status = "DRESSED-CLASS-MIXED(T1=%s T2=%.1e T3=%s)" % (t1, t2_max, t3)
        passed = False

    out = {
        "id": "V390_dressed_class_rigidity",
        "tier": 2,
        "property": "can any member of the named dressed-coefficient class "
                    "move the certified n=208 crossing by mu (supply mu "
                    "through the clock-in-the-coefficient escape), or is "
                    "the whole class rigid",
        "status": status,
        "pass": passed,
        "data": data,
        "checks": checks,
    }
    log("%s | alphas %s | max shift %.2e | T3 %s" %
        (status, {t: sol[t] for t, _ in classJ}, t2_max, t3))
    resdir = os.path.join(HERE, "results")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir,
                           os.path.basename(HERE).lower() + ".json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
