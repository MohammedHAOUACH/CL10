#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V197 - the rung arithmetic: which exponent does every scale require (O2/O8)

The ladder masses sit on integer rungs exactly: `ladder_scale(n) = M_P * LOOP**n`
with n = RUNG_A the singlet A, n = RUNG_B the singlet B (the dark matter), and
n = 2*RUNG_B - RUNG_A the bare neutrino mass m_D^2 / M_N. The rung is a protected
pole: the one-step cell matrix factorises as M(q) = kron(P(q), U_cell^6) (V196),
so no bond dressing can move a rung mass; every scale that is NOT on an integer
rung can only be paid by a loop weight (O8), never by a mass shift.

This programme prints the rung exponent x(m) = log(m / M_P) / log(LOOP) of every
scale the paper quotes, and states the joint falsifier: rung-integer structure,
the dark matter at rung RUNG_B, and a measured neutrino mass near the 1/3
candidate are mutually exclusive - the pair fixes the seesaw partner off any
integer rung. No number is fitted: every entry is computed from cl10/config.yaml
inputs and the ladder definition.

Two conventions for a dimensionless phase are printed side by side: the rung
exponent against M_P (same ruler as the masses) and the bare loop exponent of
the phase itself.
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

from cl10 import cell, config, u10pred                        # noqa: E402

LOOP = u10pred.LOOP
MP = u10pred.M_PLANCK
LOG_RATIO = math.log(LOOP) / math.log(MP)   # log(LOOP)/log(M_P)


def rung_exp(m_gev):
    """x with m = M_P * LOOP**x, so ladder_scale(x) = m exactly."""
    return math.log(m_gev / MP) / math.log(LOOP)


def bare_exp(x_gev):
    """Exponent of a dimensionless phase measured against the loop alone."""
    return math.log(float(x_gev)) / math.log(LOOP)


# The two desert thresholds of the coupling repair (paper Sec. 13, the fitted
# pair; parameter_audit lists them as coincidence-class inputs).
DESERT_PAIR = [("M_10 (10-pair threshold)", 1.2e7),
               ("M_(1,2) (1-2 pair threshold)", 1.0e12)]

# The trace-phase suppression of the cosmological constant (paper Sec. 12.4 /
# O1): a dimensionless phase, both conventions printed.
LAMBDA_PHASE = 1.6e-120


def main():
    t0 = time.time()
    out = {"experiment": "V197_rung_arithmetic", "data": {}}
    checks = []

    masses = cell.masses()
    nl = u10pred.neutrino_ladder()

    # ---- T1 the ladder entries land on integer rungs exactly --------------
    rows_ladder = []
    exact_max = 0.0
    for name, m, rung in (
            ("singlet A (seesaw)", masses["m_A_GeV"], u10pred.RUNG_A),
            ("singlet B (dark matter)", masses["m_B_GeV"], u10pred.RUNG_B),
            ("bare neutrino m_D^2/M_N", masses["m_nu_eV_identity"] * 1e-9,
             2 * u10pred.RUNG_B - u10pred.RUNG_A),
            ("Planck (rung 0)", MP, 0)):
        x = rung_exp(m)
        exact_max = max(exact_max, abs(x - rung))
        rows_ladder.append({"name": name, "value_gev": m,
                            "rung_exp": round(x, 6), "rung": rung,
                            "off": x - rung})
    checks.append(("T1 the ladder entries are integer-exact on the rung ruler: "
                   "A %.6f, B %.6f, bare neutrino %.6f, Planck 0 (max off "
                   "%.1e) - and the rung is a protected pole (M(q) = kron "
                   "(P(q), U_cell^6), V196), so no bond dressing can move a "
                   "rung mass; off-rung scales can only be paid by loop "
                   "weights"
                   % (rows_ladder[0]["rung_exp"], rows_ladder[1]["rung_exp"],
                      rows_ladder[2]["rung_exp"], exact_max),
                   exact_max < 1e-9))

    # ---- T2 the off-ladder measured scales --------------------------------
    rows_off = []
    for name, m in (("electroweak v", u10pred.V_EW),
                    ("desert unification M_X", float(config.CFG.desert.m_x_gev))
                    ) + tuple(DESERT_PAIR):
        x = rung_exp(m)
        nearest = int(round(x))
        rows_off.append({"name": name, "value_gev": m,
                         "rung_exp": round(x, 4), "nearest_rung": nearest,
                         "off": round(x - nearest, 4)})
    offs = [abs(r["off"]) for r in rows_off]
    checks.append(("T2 off-ladder scales on the same ruler: %s; every one is "
                   "off by at least %.2f rungs, i.e. no integer rung supplies "
                   "it - the residue is loop-weight territory (O8), not a "
                   "shifted mass"
                   % ("; ".join("%s %.2f (nearest rung %d, off %+.2f)" %
                                (r["name"], r["rung_exp"], r["nearest_rung"],
                                 r["off"]) for r in rows_off),
                      min(offs)),
                   all(o > 1e-9 for o in offs)))

    # ---- T3 the neutrino pair and the trace phase --------------------------
    m_atm = float(config.CFG.neutrino.m_atm_ev) * 1e-9      # GeV
    m_cand = nl["m_nu_eV_third"] * 1e-9                     # 1/3 candidate
    x_atm, x_cand = rung_exp(m_atm), rung_exp(m_cand)
    x_bare = rung_exp(masses["m_nu_eV_identity"] * 1e-9)    # exactly 17
    x_lam_rung = rung_exp(LAMBDA_PHASE)
    x_lam_bare = bare_exp(LAMBDA_PHASE)
    checks.append(("T3 neutrino and phase exponents: bare ladder m_nu = %.6f "
                   "(rung %d exact); measured m_atm = %.4f eV needs rung "
                   "%.4f (off %+.2f from rung %d); the 1/3 candidate %.4f eV "
                   "needs %.4f (off %+.2f); the trace-phase suppression "
                   "1.6e-120 sits at rung %.4f against M_P, bare loop "
                   "exponent %.4f - nearest integer rung %d, %.2f rungs away: "
                   "no integer rung supplies either (O1/O2)"
                   % (x_bare, 2 * u10pred.RUNG_B - u10pred.RUNG_A,
                      config.CFG.neutrino.m_atm_ev, x_atm,
                      x_atm - round(x_atm), round(x_atm),
                      nl["m_nu_eV_third"], x_cand, x_cand - round(x_cand),
                      x_lam_rung, x_lam_bare, int(round(x_lam_rung)),
                      abs(x_lam_rung - round(x_lam_rung))),
                   abs(x_bare - (2 * u10pred.RUNG_B - u10pred.RUNG_A)) < 1e-9
                   and abs(x_atm - round(x_atm)) > 0.05
                   and abs(x_lam_rung - round(x_lam_rung)) > 0.05))

    # ---- T4 the joint falsifier -------------------------------------------
    # rung-integer structure + m_B on rung B + m_nu measured at the candidate:
    # the seesaw partner m_A required is m_B^2 / m_nu, whose rung is forced.
    m_pair = masses["m_B_GeV"] ** 2 / m_cand
    x_A_req = rung_exp(m_pair)
    A_int = int(round(x_A_req))
    checks.append(("T4 joint falsifier: with the dark matter at rung %d and a "
                   "measured sum m_nu at the 1/3 candidate %.4f eV, the "
                   "seesaw partner m_B^2/m_nu = %.3e GeV is forced to rung "
                   "%.4f - not an integer (%.2f rungs from rung %d, where the "
                   "built A sits): the three facts [integer rungs, m_B at rung "
                   "%d, m_nu near 0.05 eV] are mutually exclusive and one must "
                   "give; a confirmation of both measurements refutes the rung "
                   "arithmetic and demands a named mechanism"
                   % (u10pred.RUNG_B, nl["m_nu_eV_third"],
                      m_pair, x_A_req, abs(x_A_req - A_int), A_int,
                      u10pred.RUNG_B),
                   abs(x_A_req - A_int) > 0.05))

    ok = all(c[1] for c in checks)
    out["data"]["ladder_rows"] = rows_ladder
    out["data"]["off_ladder_rows"] = rows_off
    out["data"]["neutrino"] = {"x_bare": x_bare, "x_m_atm": x_atm,
                               "x_third_candidate": x_cand,
                               "x_A_required_by_pair": x_A_req}
    out["data"]["lambda_phase"] = {"rung_convention": x_lam_rung,
                                   "bare_loop_convention": x_lam_bare}
    out["pass"] = bool(ok)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "the rung arithmetic of every quoted scale (computed, nothing fitted): "
        "integer rung %d for the singlet A, %d for the singlet B (the dark "
        "matter) and %d for the bare neutrino mass, exact to %.0e; the "
        "electroweak v sits %.2f rungs off its nearest rung, M_X %.2f off, "
        "the two desert thresholds %.2f and %.2f off, the measured m_atm "
        "%.2f off rung %d, the 1/3 neutrino candidate %.2f off, and the trace "
        "phase %.2f off its nearest rung (%.2f bare loop exponent). The rung "
        "is a protected pole of M(q) = kron(P(q), U_cell^6), so none of these "
        "residues can be a shifted mass: each is a loop weight (O8). Joint "
        "falsifier: integer rungs + m_B at rung %d + sum m_nu at 0.05 eV force "
        "m_A to rung %.4f, not an integer - one of the three must give, and a "
        "confirmation of the two measurements refutes the rung arithmetic "
        "(F3 strengthened, V197)"
        % (u10pred.RUNG_A, u10pred.RUNG_B, 2 * u10pred.RUNG_B - u10pred.RUNG_A,
           exact_max,
           abs(rows_off[0]["off"]), abs(rows_off[1]["off"]),
           abs(rows_off[2]["off"]), abs(rows_off[3]["off"]),
           abs(x_atm - round(x_atm)), round(x_atm),
           abs(x_cand - round(x_cand)),
           abs(x_lam_rung - round(x_lam_rung)), x_lam_bare,
           u10pred.RUNG_B, x_A_req))
    out["residual"] = True
    out["tier"] = 5
    out["property"] = ("every off-rung scale is a loop weight, never a shifted "
                       "mass; the rung-integer structure and a measured "
                       "0.05 eV neutrino mass cannot both stand (V197)")
    out["programme"] = "V197"
    out["wall_time"] = round(time.time() - t0, 3)

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v197_rung_arithmetic.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V197 - the rung arithmetic")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    for r in rows_ladder:
        print("  ladder %-28s %11.4e GeV  rung %8.6f (exact %d)"
              % (r["name"], r["value_gev"], r["rung_exp"], r["rung"]))
    for r in rows_off:
        print("  off    %-28s %11.4e GeV  rung %8.4f -> nearest %d (off %+.4f)"
              % (r["name"], r["value_gev"], r["rung_exp"],
                 r["nearest_rung"], r["off"]))
    print("  m_atm rung %.4f | 1/3 candidate %.4f | required m_A rung %.4f"
          % (x_atm, x_cand, x_A_req))
    print("  lambda phase: rung %.4f (M_P ruler), bare loop exponent %.4f"
          % (x_lam_rung, x_lam_bare))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
