#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V387 - Wilson masses for the sixteen parasites of the Clifford
bond (booked BEFORE the run in review/pre_registrations_root.md,
"V387_cubond_doubler_masses", suite 300).

Property (booked): que deviennent les masses des 16 parasites quand
le facteur de Wilson certifie (V65) agit sur le lien de Clifford.
Machine: U'(q) = (I_8 x W(q)) @ U_gamma(q) with
W(q) = diag(exp(-1j * w * phi(q) * d32)), phi(q) = sum_j (1 - cos q_j)
(the exact shape of V65.step_wilson, l. 122-127, lifted from the
16-dim sector to the 32-dim module), and d32 the charge diagonal of
D(theta) = cos(theta) X + sin(theta) Y in 32-dim, X and Y the Cartan
combinations with coefficients ux[0], ys[0] of V65.identify_XY(lam16)
(l. 91-116, same 32-dim Cartan basis as V385); theta = 1.7308 (the
booked physical junction angle; the V86 artifact was grepped for
"theta" first: its best theta is 1.56, the scan residuals near 1.73
are residuals not angles - the booked constant stands, with the
supporting artifact readings 1.73078/1.730784 recorded); grid
w = {0.1, 0.2, ..., 1.0}.

T0: the multiset of the 16-dim X eigenvalues equals x_target =
    tri([1x10, -3x5, 5]) (<= 1e-9) - the V65 certification replayed;
    W(w=0) = I: U'(q) == U_gamma(q) (<= 1e-14).
T1: sites: the 27-point grid {0,pi/2,pi}^3 plus the 8 cubic corners
    (+-pi/2)^3 (34 distinct sites); BOOKED PREDICTION BEFORE
    MEASUREMENT: N_pred = (number of the 32 weights of zero charge,
    |c.lambda| <= 1e-9 with c = cos(theta) ux[0] + sin(theta) ys[0])
    x 8 link factors, computed from the V385 artifact weights and
    printed; GATE: for every w >= 0.5, at every one of the 15
    nonzero-momentum gapless sites of the family, the number of
    states at phase {0,pi} within 1e-9 equals N_pred exactly.
T2: the cone lives: gap at Gamma = 0 (<= 1e-12) for every w; slopes
    S'_j = -i dU'/dq_j by central finite difference at q = 1e-4 e_j,
    h = 1e-5; the Clifford projection M_jk = (1/256) Re tr(S'_j
    gamma_k) keeps M = I_3 within 1e-5 for every w (the W factor is
    O(q^2) at the origin; the off-Clifford tilt sin(q_j) (I x d32),
    orthogonal to the Clifford basis, is recorded not gated).
T3: for w = 1.0, at the parasite sites the minimum gap of the gapped
    states is recorded (data, no gate).

Verdicts: CUBOND-DOUBLER-MASS-SELECTIVE / WILSON-DOUBLERS-EVADE /
WILSON-KILLS-CONE.
"""
import importlib.util
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

from cl10 import supercell_2i as S2  # noqa: E402
from cl10 import clifford  # noqa: E402
from cl10.paths import exp  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v65", os.path.join(ROOT, "archive", "experiments", "dynamics_spectrum",
                        "V65_wilson", "V65_wilson.py"))
_v65 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v65)

THETA = 1.7308
GRID_W = [round(0.1 * k, 1) for k in range(1, 11)]

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, bool(ok)))


def N(A):
    return float(np.max(np.abs(A)))


def main():
    data = {}

    # ---- inputs: lam16 and weights32 from the V385 artifact ----
    v385 = json.load(open(os.path.join(
        exp("V385_cubond_content_census"), "results",
        "v385_cubond_content_census.json")))
    w32 = np.array(v385["data"]["weights32"])
    lam16 = np.array(v385["data"]["lam16"])

    # ---- the certified X, Y: V65.identify_XY replayed (T0) ----
    ux, ys = _v65.identify_XY(lam16)
    ux0, ys0 = np.asarray(ux[0], dtype=float), np.asarray(ys[0], dtype=float)
    x16 = np.sort(lam16 @ ux0)
    x_target = np.sort(np.array([1.0] * 10 + [-3.0] * 5 + [5.0]))
    ex = float(np.max(np.abs(x16 - x_target)))
    data["Y_coef"] = [float(v) for v in ys0]
    data["X_coef"] = [float(v) for v in ux0]

    # ---- the 32-dim Cartan basis (same construction as V385) ----
    h = [np.asarray(m, dtype=np.complex128) for m in clifford.generators()]
    Ck = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    M = sum(math.sqrt(p) * Ck[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(M)
    w32_here = np.array([[float(np.real(np.vdot(V[:, s], Ck[k] @ V[:, s])))
                          for k in range(5)] for s in range(32)])
    # order-insensitive agreement with the V385 artifact multiset
    _rows = lambda A: sorted(tuple(round(float(x), 9) for x in r)
                             for r in A)
    dev32 = max(max(abs(a - b) for a, b in zip(r1, r2))
                for r1, r2 in zip(_rows(w32_here), _rows(w32)))
    chk("T0 the V65 certification replayed on the 32-dim Cartan basis: "
        "the 16-dim X eigenvalue multiset equals x_target {1x10, -3x5, "
        "5x1} to %.1e (<= 1e-9), and the 32-weight multiset matches the "
        "V385 artifact componentwise to %.1e: %s"
        % (ex, dev32, max(ex, dev32) <= 1e-9), max(ex, dev32) <= 1e-9)

    # ---- the engine, the Wilson factor, theta ----
    Sx, Sy, Sz = S2.spin_triple()
    R = [2 * np.asarray(m, dtype=np.complex128) for m in (Sx, Sy, Sz)]
    X, Y, Z = (np.array(m, dtype=np.complex128) for m in
               ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]))
    I2 = np.eye(2, dtype=np.complex128)
    P = [np.kron(np.kron(X, I2), I2), np.kron(np.kron(I2, Y), I2),
         np.kron(np.kron(I2, I2), Z)]
    G = [np.kron(p, r) for p, r in zip(P, R)]
    I8 = np.eye(8, dtype=np.complex128)
    I256 = np.eye(256, dtype=np.complex128)

    def Uq(q):
        U = I256.copy()
        for j in range(3):
            U = U @ (math.cos(q[j]) * I256 + 1j * math.sin(q[j]) * G[j])
        return U

    # D(theta) = cos(theta) X + sin(theta) Y as charges of the 32 weights
    d32 = math.cos(THETA) * (w32_here @ ux0) + math.sin(THETA) * (w32_here @ ys0)
    data["theta_used"] = THETA
    data["theta_support"] = 1.73078
    data["d32_min_abs_nonzero"] = float(np.min(np.abs(d32[np.abs(d32) > 1e-9]))) \
        if np.any(np.abs(d32) > 1e-9) else 0.0
    data["d32_max_abs"] = float(np.max(np.abs(d32)))

    def Wq(q, w):
        ph = w * sum(1.0 - math.cos(qj) for qj in q)
        return V @ np.diag(np.exp(-1j * ph * d32)) @ V.conj().T

    def Up(q, w):
        return np.kron(I8, Wq(q, w)) @ Uq(q)

    # T0 half 2: W(w=0) = I  ->  U'(q) == U_gamma(q) bit for bit
    rng = np.random.default_rng(20261105)
    e0b = 0.0
    for _ in range(6):
        q = list(rng.uniform(-3, 3, 3))
        e0b = max(e0b, N(Up(q, 0.0) - Uq(q)))
    chk("T0 second half: W(w=0) = I - U'(q) and U_gamma(q) agree at 6 "
        "random q to %.1e (<= 1e-14): %s" % (e0b, e0b <= 1e-14),
        e0b <= 1e-14)

    # ---- BOOKED PREDICTION BEFORE MEASUREMENT: N_pred ----
    c5 = math.cos(THETA) * ux0 + math.sin(THETA) * ys0
    zero_weights = int(np.sum(np.abs(w32_here @ c5) <= 1e-9))
    N_pred = zero_weights * 8
    log("BOOKED N_pred computed before the scan: %d of the 32 weights "
        "have zero charge (c = cos(theta) ux + sin(theta) ys, smallest "
        "|charge| %.6f), x 8 link factors -> N_pred = %d"
        % (zero_weights, float(np.min(np.abs(w32_here @ c5))), N_pred))
    data["N_pred"] = N_pred
    data["zero_charge_weights"] = zero_weights

    # ---- sites: the 27-point grid plus the 8 cubic corners ----
    vals = (0.0, math.pi / 2, math.pi)
    grid = {(q1, q2, q3) for q1 in vals for q2 in vals for q3 in vals}
    corners = {(s1 * math.pi / 2, s2 * math.pi / 2, s3 * math.pi / 2)
               for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)}
    sites = sorted(grid | corners)
    paras = [s for s in sites if s != (0.0, 0.0, 0.0)
             and (s in corners or all(q in (0.0, math.pi) for q in s))]

    def dist_to_0pi(ph):
        return np.minimum(np.abs(ph), np.abs(np.abs(ph) - math.pi))

    # ---- T1: per w, gapless sites and zero-mode state counts ----
    t1 = {}
    all_ok = True
    for w in GRID_W:
        log("T1 scan w = %.1f" % w)
        site_rows = {}
        for q in sites:
            ph = np.sort(np.angle(np.linalg.eigvals(Up(list(q), w))))
            dmin = float(np.min(dist_to_0pi(ph)))
            nstates = int(np.sum(dist_to_0pi(ph) <= 1e-9))
            site_rows[repr([round(x, 4) for x in q])] = {
                "gap": dmin, "zero_mode_states": nstates}
        bad = []
        for s in paras:
            k = repr([round(x, 4) for x in s])
            if w >= 0.5 and site_rows[k]["zero_mode_states"] != N_pred:
                bad.append((k, site_rows[k]["zero_mode_states"]))
        if bad:
            all_ok = False
        t1[str(w)] = {"gapless_sites": [k for k, v in site_rows.items()
                                        if v["gap"] <= 1e-9],
                      "parasite_zero_mode_counts": {
                          repr([round(x, 4) for x in s]):
                          site_rows[repr([round(x, 4) for x in s])]
                          ["zero_mode_states"] for s in paras},
                      "violations": bad}
    data["t1_per_w"] = t1
    chk("T1 the booked survivors: N_pred = %d; for every w >= 0.5 every "
        "one of the %d nonzero-momentum gapless family sites holds "
        "exactly N_pred states at phase 0/pi: %s (violations listed per "
        "w in data: none: %s)"
        % (N_pred, len(paras), all_ok,
           "none" if all_ok else "see data.t1_per_w.wviolations"),
        all_ok)

    # ---- T2: the cone lives; slopes by Clifford projection ----
    e_gamma = 0.0
    slope_err = 0.0
    tilt_max = 0.0
    hs = 1e-5
    a = 1e-4
    for w in GRID_W:
        ph = np.angle(np.linalg.eigvals(Up([0.0, 0.0, 0.0], w)))
        e_gamma = max(e_gamma, float(np.max(np.abs(ph))))
        for j in range(3):
            q0 = [0.0, 0.0, 0.0]
            q0[j] = a
            e = [0.0, 0.0, 0.0]
            e[j] = hs
            Sp = -1j * (Up([q0[k] + e[k] for k in range(3)], w)
                        - Up([q0[k] - e[k] for k in range(3)], w)) / (2 * hs)
            row = np.array([float(np.real(np.trace(Sp @ G[k])) / 256.0)
                            for k in range(3)])
            slope_err = max(slope_err, float(np.max(np.abs(row - np.eye(3)[j]))))
            rec = sum(row[k] * G[k] for k in range(3))
            tilt_max = max(tilt_max, N(Sp - rec))
    data["T2_slope_error"] = slope_err
    data["T2_off_clifford_tilt_max"] = tilt_max
    chk("T2 the cone lives: gap at Gamma %.1e <= 1e-12 for every w; "
        "slopes: Clifford projection of -i dU'/dq_j at q = 1e-04, "
        "h = 1e-05: max |M - I_3| = %.1e <= 1e-5 (the W factor is "
        "O(q^2) at the origin); off-Clifford tilt max %.1e recorded "
        "(orthogonal to the Clifford basis): %s"
        % (e_gamma, slope_err, tilt_max,
           e_gamma <= 1e-12 and slope_err <= 1e-5),
        e_gamma <= 1e-12 and slope_err <= 1e-5)

    # ---- T3: masses recorded at w = 1.0 ----
    masses = {}
    w1 = 1.0
    for s in paras:
        ph = np.sort(np.angle(np.linalg.eigvals(Up(list(s), w1))))
        d = dist_to_0pi(ph)
        gapped = d[d > 1e-9]
        masses[repr([round(x, 4) for x in s])] = (
            float(np.min(gapped)) if len(gapped) else 0.0)
    data["t3_masses_w1"] = masses
    mmin = min(masses.values())
    mmax = max(masses.values())
    chk("T3 at w = 1.0 the parasite sites are gapped, the minimum "
        "gapped-state gap recorded per site in data: floor %.6f, "
        "max %.6f (all positive: %s)"
        % (mmin, mmax, mmin > 1e-9), mmin > 1e-9)

    t1_ok, t2_ok = all_ok, (e_gamma <= 1e-12 and slope_err <= 1e-5)
    if t1_ok and t2_ok:
        status = ("CUBOND-DOUBLER-MASS-SELECTIVE (the certified Wilson "
                  "factor of V65 acts on the Clifford bond: every "
                  "nonzero-momentum gapless site of the family acquires "
                  "a mass (gap floor %.4f at w = 1.0) while the cone at "
                  "Gamma survives exactly and its slopes stay (1,1,1) "
                  "to 1e-5; the survivors are exactly the content of "
                  "zero charge: N_pred = %d of the 32 internal weights "
                  "at the booked junction angle theta = 1.7308 - the "
                  "sixteen-parasite problem is resolved at the level of "
                  "the certified mechanism)" % (mmin, N_pred))
    elif not t1_ok:
        status = "WILSON-DOUBLERS-EVADE (T1: the survivor census failed "
        "the booked prediction - see data.t1_per_w violations)"
    else:
        status = "WILSON-KILLS-CONE (T2: the cone at Gamma or its slopes " \
                 "failed - see check lines)"
    out = {"id": "V387_cubond_doubler_masses",
           "tier": 2,
           "residual": max(ex, dev32, e0b, e_gamma, slope_err),
           "property": ("Wilson masses of the parasites of the Clifford "
                        "bond, certified V65 mechanism, tier 2"),
           "status": status,
           "pass": t1_ok and t2_ok,
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v387_cubond_doubler_masses.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V387 done] %s" % status[:300])
    return 0 if (t1_ok and t2_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
