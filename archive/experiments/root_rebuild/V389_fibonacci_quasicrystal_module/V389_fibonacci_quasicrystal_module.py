#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V389 - the quasicrystal Z[phi] module: can the Fibonacci module
supply mu as a NEW CERTIFIED INVARIANT (route R1; booked BEFORE the
run in review/pre_registrations_root.md, "V389_fibonacci_quasicrystal_module").

The certified ladder object is replaced by its golden-inflation
analogue: the composite stage S_phi(g) is the ordered product of the
certified 45-plane step over the Fibonacci word, clocks m_j = 1 (L)
and phi (S). The certified functional is untouched:
F_n^phi = cartan_norm(Bn^phi)/cartan_norm(S^n). Two legs are tested
against targets fixed by artifacts BEFORE implementation:
  T1 closure leg     - clock-free closure over five clocks (the V356
                       ico-scale leg);
  T2 invariant leg   - rho_phi = rho x mu^(1/208) (the new invariant);
  T3 crossing leg    - first rung crossing deficit = mu;
  T4 content leg     - the 32 internal content is untouched.

Verdicts: FIBONACCI-SUPPLIES-MU / FIBONACCI-PARTIAL / FIBONACCI-NULL.
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

import numpy as np  # noqa: E402

from cl10.paths import exp  # noqa: E402
import V350_backreaction_fixed_point as V350  # noqa: E402

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


def N(A):
    return float(np.max(np.abs(A)))


I32 = np.eye(32, dtype=np.complex128)
PHI = (1.0 + math.sqrt(5.0)) / 2.0
CLOCK = {"pi/4": math.pi / 4.0, "1/5": 0.2, "1/6": 1.0 / 6.0,
         "1/7": 1.0 / 7.0, "pi/6": math.pi / 6.0}
MU = float(art("V325_channel_triple_f208", "data.certified.mu_req"))
RHO = float(art("V276_certified_ladder_208", "data.machinery.rho"))
TARGET = RHO * MU ** (1.0 / 208.0)
COEF_AT_FIFTH = float(art("V206_clock_tension",
                          "data.data_clock.coef_at_fifth"))
LOOP = float(V350.LOOP)


def fib_words():
    """Fibonacci words w_{k+1} = w_k w_{k-1}; return all up to |w|>=377."""
    v, w = "L", "LS"
    words = [v, w]
    while len(w) < 377:
        w, v = w + v, w
        words.append(w)
    return words


def stage(word, g):
    S = I32.copy()
    for c in word:
        S = S @ V350.step(g * (1.0 if c == "L" else PHI))
    return S


def ladder(S, g, nmax=400):
    """The certified chain on the composite stage (Fn_corr structure,
    J = 0): B_{n} = B_{n-1} (M_phi S), Gn = S^n, F_n = norm ratio."""
    Ma = V350.hops_hoisted(
        S, np.eye(V350.D_L, dtype=np.complex128) / V350.D_L)
    ph = [np.exp(1j * (g / 2.0) * e) for e in V350.EV]
    M = sum(ph[i] * Ma[i] for i in range(3))
    MG = M @ S
    B, Gn = [I32], [I32]
    for _ in range(nmax):
        B.append(B[-1] @ MG)
        Gn.append(Gn[-1] @ S)
    return [1.0] + [float(np.real(V350.cartan_norm(B[n], V350.CARTAN)
                                  / V350.cartan_norm(Gn[n], V350.CARTAN)))
                    for n in range(1, nmax + 1)]


def main():
    words = fib_words()
    wA, wB, wC = words[-3], words[-2], words[-1]
    lens = [len(w) for w in words]
    data["word_lengths"] = lens[-4:]

    # ---- T0 construction guards ----------------------------------------
    ok_len = (len(wA), len(wB), len(wC)) == (144, 233, 377)
    chk("T0a exact Fibonacci lengths: last three words have lengths "
        "(144, 233, 377) = Fibonacci integers (exact string recurrence)",
        ok_len)

    S12 = stage(wC, CLOCK["pi/4"])
    dev_u = N(S12.conj().T @ S12 - I32)
    data["unitarity_dev"] = dev_u
    chk("T0b unitarity: max|S_phi^(144->377)(pi/4)^dag S - I| = %.1e "
        "(<= 1e-10)" % dev_u, dev_u <= 1e-10)

    Su = stage("L" * 377, CLOCK["pi/4"])
    P = I32.copy()
    Gp = V350.step(CLOCK["pi/4"])
    for _ in range(377):
        P = P @ Gp
    dev_id = N(Su - P)
    data["uniform_identity_dev"] = dev_id
    chk("T0c algebra identity: uniform modulation (m = 1) gives "
        "S = step(g)^377 to %.1e (<= 1e-12)" % dev_id, dev_id <= 1e-12)

    Fn_host = V350.Fn_of(0.2, nmax=10)
    c15 = 3.0 * (abs(Fn_host[10]) / LOOP) ** 10
    dev_read = abs(c15 - COEF_AT_FIFTH)
    data["undressed_fifth_read"] = c15
    chk("T0d host machinery: undressed fifth read %.6f vs V206 "
        "coef_at_fifth artifact %.6f (|d| %.1e <= 1e-9)"
        % (c15, COEF_AT_FIFTH, dev_read), dev_read <= 1e-9)

    t0_ok = ok_len and dev_u <= 1e-10 and dev_id <= 1e-12 and dev_read <= 1e-9

    # ---- T1 closure leg -------------------------------------------------
    clos = {}
    for tag, g in CLOCK.items():
        S = S12 if tag == "pi/4" else stage(wC, g)
        qv, Q = None, I32.copy()
        for q in range(1, 13):
            Q = Q @ S
            d = min(N(Q - I32), N(Q + I32))
            if d <= 1e-10:
                qv = q
                break
        clos[tag] = qv
    data["closure_q"] = {k: (v if v else "none<=12") for k, v in clos.items()}
    clockfree = all(v is not None for v in clos.values())
    t1_word = ("CLOSURE-CLOCK-FREE q=%d on all five clocks"
               % sorted(clos.values())[0]) if clockfree else \
        ("CLOCK-CONSUMING: no q <= 12 closes simultaneously (per-clock: "
         + ", ".join("%s:%s" % (k, v if v else "none")
                     for k, v in clos.items()) + ")")
    chk("T1 closure leg: " + t1_word, clockfree)

    # ---- T2 invariant leg -----------------------------------------------
    rs, Fsets = {}, {}
    for tag, w in (("F12=144", wA), ("F13=233", wB), ("F14=377", wC)):
        S = S12 if tag == "F14=377" else stage(w, CLOCK["pi/4"])
        F = ladder(S, CLOCK["pi/4"], 400)
        Fsets[tag] = F
        rs[tag] = math.sqrt(F[24] / F[22])
    spread = max(rs.values()) - min(rs.values())
    data["ratio_k"] = rs
    data["ratio_spread"] = spread
    rho_phi = sum(rs.values()) / 3.0
    data["rho_phi"] = rho_phi
    data["gap_to_target"] = abs(rho_phi - TARGET)
    data["gap_to_rho"] = abs(rho_phi - RHO)
    data["gap_to_rho_mu"] = abs(rho_phi - RHO * MU)
    data["target"] = TARGET
    t2_exists = spread <= 2e-3
    t2_hit = t2_exists and abs(rho_phi - TARGET) <= 1e-6
    chk("T2 invariant leg: convergence spread %.1e (<= 2e-3: %s); "
        "rho_phi = %.10f vs target rho x mu^(1/208) = %.10f "
        "(|d| = %.2e <= 1e-6: %s); gaps to rho %.2e, rho*mu %.2e"
        % (spread, t2_exists, rho_phi, TARGET, abs(rho_phi - TARGET),
           t2_hit, abs(rho_phi - RHO), abs(rho_phi - RHO * MU)), t2_hit)

    # ---- T3 crossing leg ------------------------------------------------
    cross = {}
    two_pi = 16.0 * math.pi
    for tag, g in CLOCK.items():
        F = ladder(stage(wC, g), g, 400)
        ns = [n for n in range(1, 401) if two_pi * F[n] <= 1.0]
        if ns:
            nst = ns[0]
            deficit = two_pi * F[nst]
            cross[tag] = (nst, deficit)
        else:
            cross[tag] = None
    data["crossings"] = {k: (v if v else "no-crossing<=400")
                         for k, v in cross.items()}
    t3_hits = [(k, v) for k, v in cross.items()
               if v and abs(v[1] - MU) <= 1e-6]
    best = min((abs(v[1] - MU), k) for k, v in cross.items() if v) \
        if any(cross.values()) else (9.9, "none")
    data["nearest_mu_gap"] = best[0]
    chk("T3 crossing leg: first 16pi F_n^phi crossing n*/deficit per clock "
        "%s; nearest |deficit - mu| = %.3e (gate 1e-6: %s)"
        % ({k: (v if v else "none") for k, v in cross.items()},
           best[0], bool(t3_hits)), bool(t3_hits))

    # ---- T4 content leg --------------------------------------------------
    ev = np.linalg.eigvals(S12)
    dev_mod = float(np.max(np.abs(np.abs(ev) - 1.0)))
    data["eig_modulus_dev"] = dev_mod
    th = np.sort(np.angle(ev))
    Pw = I32.copy()
    for _ in range(377):
        Pw = Pw @ Gp
    th0 = np.sort(np.angle(np.linalg.eigvals(Pw)))
    disp = float(np.mean(np.abs((th - th0 + np.pi) % (2 * np.pi) - np.pi)))
    data["spectral_shift_vs_certified"] = disp
    chk("T4 content leg: spectrum of S_phi^(377)(pi/4) on the unit circle "
        "(max ||lambda|-1| = %.1e <= 1e-10); sorted eigenphase shift vs "
        "certified step^377 = %.4f rad/level (recorded only)"
        % (dev_mod, disp), dev_mod <= 1e-10)

    # ---- decision --------------------------------------------------------
    data["ratio_nonconv"] = not t2_exists
    legs = int(t2_hit) + int(bool(t3_hits))
    if legs == 2:
        status = "FIBONACCI-SUPPLIES-MU"
    elif legs == 1:
        status = "FIBONACCI-PARTIAL"
    else:
        status = "FIBONACCI-NULL"
    passed = t0_ok

    out = {
        "id": "V389_fibonacci_quasicrystal_module",
        "tier": 1,
        "property": "can the quasicrystal Z[phi] (Fibonacci) module carry "
                    "both legs - clock-free closure and rho x mu^(1/208) "
                    "with a mu crossing - i.e. supply mu as a new certified "
                    "invariant",
        "status": status,
        "pass": passed,
        "data": data,
        "checks": checks,
    }
    log("%s | T0 ok=%s closure=%s rho-hit=%s crossing-hits=%d" %
        (status, t0_ok, clockfree, t2_hit, len(t3_hits)))
    resdir = os.path.join(HERE, "results")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir,
                           os.path.basename(HERE).lower() + ".json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
