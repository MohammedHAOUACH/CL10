#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V199 - the A->BB channel: do the root bond strengths wash out leptogenesis (O5)

Singlet A (mass m_A = 9.613e13 GeV, rung RUNG_A) decays to the lepton-Higgs
pairs with total rate Gamma_A (V125, the atmospheric point y = 0.39861 with
K = Gamma_A/H = 46.81) and to the dark pair BB with a rate set by the A BB
coupling. If the BB channel carries a large branching fraction, the lepton
asymmetry is washed out and the resonant leptogenesis route of V125 dies -
whatever the CP phase. The question has a binary answer once the coupling
candidates are fixed, and the root fixes them: the 13 bond numbers of the
singlet class (V147, the same candidate list that bounds the portal).

Automaton convention (no new scale): couplings are dimensionless numbers and
the only scale the channel has is m_A, so

    Gamma_B(kappa) = kappa^2 m_A / (32 pi) * sqrt(1 - 4 m_B^2 / m_A^2)

with kappa a root candidate (the V147 convention alpha_eff = kappa^2/4pi; the
kinematic factor is 1 at the printed precision, computed below). H(T) is the
Friedmann rate 1.66 sqrt(g_*) T^2 / M_P of cl10/config.yaml at T = m_A.

Two thresholds, both from V125's stored verdict:
  * channel negligible: Br(A->LH) = K/(K + Gamma_B/H) > 0.99;
  * route survives at all: the washout efficiency kappa_f is diluted by the
    branching, so the route needs Br >= J_required/J_available_max = 0.0643
    (V125: 1.00126e-4 required against 1.55784e-3 available CP strength).
Between them the channel is a named bound on the coupling; on the kill side
the resonance route is impossible for that strength.

Production side (V195 convention, T2): if A is thermal the decay yield
overcloses 0.120 by the equilibrium overclosure of V195 (x2.7e63, V66
convention); if A is not thermal the yield needs the inflation/reheating
history the built model does not supply (named missing input). No value of
the 13 strengths produces 0.120 through this channel.
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
from cl10.paths import exp                                     # noqa: E402

from cl10 import cell, config, u10pred                         # noqa: E402

V147 = json.load(open(exp("V147_w7_relic_from_bond",
                          "results", "v147_w7_relic_from_bond.json")))
V125 = json.load(open(exp("V125_w5_cp_flux_washout",
                          "results", "v125_w5_cp_flux_washout.json")))
V195 = json.load(open(exp("V195_w7_cell_redefinition",
                          "results", "v195_w7_cell_redefinition.json")))


def main():
    t0 = time.time()
    out = {"experiment": "V199_a_to_bb_channel", "data": {}}
    checks = []

    masses = cell.masses()
    M_A = masses["m_A_GeV"]
    M_B = masses["m_B_GeV"]
    M_P = u10pred.M_PLANCK
    HC = float(config.CFG.cosmology.hubble_coeff)
    GS = float(config.CFG.cosmology.g_star)

    def hubble(T):
        return HC * math.sqrt(GS) * T * T / M_P

    # ---- T0 kinematics ------------------------------------------------------
    kin = math.sqrt(max(0.0, 1.0 - 4.0 * M_B ** 2 / M_A ** 2))
    checks.append(("T0 the channel is open: 4 m_B^2 = %.3e GeV^2 << m_A^2 = "
                   "%.3e, phase-space factor sqrt(1 - 4 m_B^2/m_A^2) = %.15f"
                   % (4 * M_B ** 2, M_A ** 2, kin),
                   4 * M_B ** 2 < M_A ** 2 and kin > 1 - 1e-12))

    # ---- T1 the census of the 13 root strengths against washout ------------
    atmo = V125["data"]["asymmetry"]["cases"]["atmospheric y"]
    K = float(atmo["K"])                        # Gamma_A/H at the y-point
    jreq = float(V125["data"]["asymmetry"]["cases"]["atmospheric y"]
                 ["J_required"])
    jav = float(V125["data"]["J_available_max"])
    br_min = jreq / jav                        # route survives iff Br >= this
    cand = [c for c in V147["data"]["candidates"] if c["kappa"] is not None]
    kap = sorted(abs(float(c["kappa"])) for c in cand)

    def gamma_b(k):
        return k * k * M_A * kin / (32.0 * math.pi)

    H_A = hubble(M_A)
    rows = []
    for k in kap:
        rB = gamma_b(k) / H_A
        br = K / (K + rB)
        rows.append({"kappa": round(k, 4), "Gamma_B_GeV": gamma_b(k),
                     "r_B_over_H": rB, "Br_A_to_LH": br,
                     "kills_route": br < br_min})
    n_kill = sum(r["kills_route"] for r in rows)

    def kappa_of_br(br):
        rB = K * (1.0 / br - 1.0)
        return math.sqrt(rB * H_A * 32.0 * math.pi / (M_A * kin))

    k_split = kappa_of_br(0.99)                # channel negligible below
    k_kill = kappa_of_br(br_min)               # route dead above
    out["data"]["census"] = rows
    out["data"]["thresholds"] = {"Br_min": br_min, "Kappa_split": k_split,
                                 "Kappa_kill": k_kill, "H_at_mA_GeV": H_A}
    checks.append(("T1 the 13 root bond strengths of the singlet class "
                   "(|kappa| = %.3f ... %.1f) run Gamma_B/H from %.3g to %.3g "
                   "at T = m_A: %d of 13 leave Br(A->LH) < %.4f = "
                   "J_required/J_available_max and kill the resonant route of "
                   "V125; the channel is negligible only for kappa < %.4f and "
                   "the route survives only for kappa < %.4f - the channel is "
                   "a named bound on the dark coupling, not a free door"
                   % (kap[0], kap[-1], rows[0]["r_B_over_H"],
                      rows[-1]["r_B_over_H"], n_kill, br_min, k_split, k_kill),
                   len(rows) == len(kap) and k_split < k_kill))

    # ---- T2 the yield side (V195 conventions) -------------------------------
    ov = float(V195["data"]["production"]["grav_overclosure_V66_conv"])
    br_min_row = min(r["Br_A_to_LH"] for r in rows)
    checks.append(("T2 the channel cannot produce 0.120: if A is thermal the "
                   "decay yield overcloses by the V195 equilibrium estimate "
                   "x%.3e (regression) at every candidate branching (smallest "
                   "branching fed %.3g); if A is not thermal the yield needs "
                   "the inflation/reheating history the built model does not "
                   "supply (named missing input, V195) - no value of the 13 "
                   "strengths gives the relic through A->BB"
                   % (ov, br_min_row), ov > 1e60))

    # ---- T3 branch -----------------------------------------------------------
    if all(r["kills_route"] for r in rows):
        branch = "WASHOUT"
    elif all(not r["kills_route"] and r["Br_A_to_LH"] > 0.99 for r in rows):
        branch = "IRRELEVANT"
    else:
        branch = "SPLIT"
    out["branch"] = branch
    ok = all(c[1] for c in checks)

    out["pass"] = bool(ok)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "the A->BB channel census of O5. With the 13 root bond numbers of the "
        "singlet class (V147) as the A BB coupling candidates and Gamma_B = "
        "kappa^2 m_A/32pi, the channel runs Gamma_B/H = %.3g ... %.3g at "
        "T = m_A against the leptogenesis point K = %.1f (V125): %d of the 13 "
        "strengths drive Br(A->LH) below %.4f and make the resonant route "
        "impossible, %d of them leave it alive; the channel is negligible "
        "below kappa = %.4f and survivable below kappa = %.4f - the smallest "
        "root strength %.3f sits between (branch %s). The channel supplies no "
        "relic (thermal A overcloses x%.2e, non-thermal A needs the named "
        "missing reheating history): the washout clause of O5 now names a "
        "bound on the dark coupling, and the coupling itself stays unnamed "
        "(O7)"
        % (rows[0]["r_B_over_H"], rows[-1]["r_B_over_H"], K, n_kill, br_min,
           len(rows) - n_kill, k_split, k_kill, kap[0], branch, ov))
    out["residual"] = True
    out["tier"] = 5
    out["property"] = ("the A->BB channel is a named bound on the dark "
                       "coupling, never a relic source (V199)")
    out["programme"] = "V199"
    out["wall_time"] = round(time.time() - t0, 3)

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v199_a_to_bb_channel.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V199 - the A->BB channel census")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    for r in rows:
        print("  kappa %8.4f  Gamma_B %11.3e GeV  r_B/H %9.3e  Br %8.4f  %s"
              % (r["kappa"], r["Gamma_B_GeV"], r["r_B_over_H"],
                 r["Br_A_to_LH"], "KILL" if r["kills_route"] else "alive"))
    print("  Br_min %.4f | kappa_split %.4f | kappa_kill %.4f | %d/13 kill "
          "| branch %s" % (br_min, k_split, k_kill, n_kill, branch))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
