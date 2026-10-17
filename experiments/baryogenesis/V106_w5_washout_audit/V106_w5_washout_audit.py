#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V106 - O5 audit: the 11% Delta B agreement under the SAME Yukawa.

**Why (Section 12.2).** The product (28/79)(g0^4/16pi)(1/16pi)^2 =
6.8e-10 sits 11% above the measured 6.12e-10. The paper already shows the
out-of-equilibrium condition fails at the model's own couplings
(Gamma/H = 4.6 at y = 1/8; 23 at the atmospheric y; 68 at the ladder
m_D = 118 GeV). This experiment QUANTIFIES the consequence the product
formula ignores: the strong-washout efficiency kappa_f of the SAME decay.

Efficiency (Buchmuller-Di Bari-Plumacher, Ann. Phys. 315 (2005) 305,
hep-ph/0401240, strong-washout + branching fit):
    kappa_f(K) = 2/(z_B(K) K) * (1 - exp(-z_B(K) K / 2)),
    z_B(K) = -(5/2) W_{-1}(-(4/5) pi^{-1/5} K^{-2/5}).
Zero fitted parameters: every input is read from V57's JSON or is a
standard constant (g_* = 106.75).

PASS (as an audit) = the self-consistent y (atmospheric m_nu) gives
K >> 1 (in equilibrium, Sakharov 3 fails) and kappa_f far from the
guessed dilution D = (1/16pi)^2: the 11% product agreement is an
accident of two wrong factors, not a near-miss.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
from cl10 import eta_cp, u10pred
from cl10 import neutrino as nv
from cl10 import baryogenesis as BARY
from cl10 import cosmology as COS

M_P = u10pred.M_PLANCK     # GeV
GSTAR = COS.G_STAR
SPHALERON = eta_cp.SPHALERON
DB_OBS = eta_cp.DELTA_B_OBS


def main():
    v57 = json.load(open(exp("V57_neutrino",
                                      "results", "v57_neutrino.json")))["data"]
    M_N = v57["M_N_GeV"]
    m_D = v57["m_D_GeV"]
    H = COS.HUBBLE_COEFF * math.sqrt(GSTAR) * M_N ** 2 / M_P

    def K_of_y(y):
        return (y * y * M_N / (8.0 * math.pi)) / H

    # the physical BDP washout efficiency of the library (V106's own audit):
    # kappa_f(K) with z_B the real Lambert-W solution, in domain for K >= 3.93
    f_eff = BARY.bdp_washout_efficiency

    cases = {}
    # (a) the model's own reading g0 = 1/8
    y = 0.125
    K = K_of_y(y)
    f, zB, _ = f_eff(K)
    eps = y ** 4 / (16 * math.pi)
    cases["y=g0=1/8"] = {"y": y, "K": round(K, 2), "f": round(f, 4),
                         "z_B": round(zB, 2),
                         "DB_product": SPHALERON * eps
                         * u10pred.LOOP ** 2,
                         "DB_washed": SPHALERON * eps * f}
    # (b) the atmospheric Yukawa: m_nu = y^2 v^2 / M_N
    for v in (u10pred.V_EW / math.sqrt(2), u10pred.V_EW):
        y_atm = math.sqrt(nv.M_ATM_EV * 1e-9 * M_N) / v
        K = K_of_y(y_atm)
        f, zB, _ = f_eff(K)
        cases[f"y=atmospheric(v={v:.0f})"] = {
            "y": round(y_atm, 4), "K": round(K, 1), "f": round(f, 4),
            "z_B": round(zB, 2)}
    # (c) the ladder Yukawa m_D = 118 GeV
    for v in (u10pred.V_EW / math.sqrt(2), u10pred.V_EW):
        y_lad = m_D / v
        K = K_of_y(y_lad)
        f, zB, _ = f_eff(K)
        cases[f"y=ladder(v={v:.0f})"] = {"y": round(y_lad, 4),
                                         "K": round(K, 1),
                                         "f": round(f, 4),
                                         "z_B": round(zB, 2)}

    # y allowed by out-of-equilibrium K < 1, and the m_nu it caps
    y_ooeq = math.sqrt(H * 8.0 * math.pi / M_N)
    m_nu_max = y_ooeq ** 2 * (u10pred.V_EW / math.sqrt(2)) ** 2 / M_N * 1e9   # eV

    atm = cases["y=atmospheric(v=174)"]
    D_guessed = u10pred.LOOP ** 2
    ok = atm["K"] > 10 and atm["f"] < 0.05
    status = (f"AUDIT PASSES BY REPORTING ITS OWN FAILURE: at the SAME "
              f"Yukawa that gives m_nu = {nv.M_ATM_EV} eV (y = {atm['y']}), the "
              f"decay is IN equilibrium (K = {atm['K']} >> 1, Sakharov 3 "
              f"fails) and the physical washout efficiency is "
              f"kappa_f = {atm['f']} - {atm['f']/D_guessed:.0f}x the "
              f"guessed dilution D = (1/16pi)^2 = {D_guessed:.1e} the "
              f"product formula uses. Out-of-equilibrium decay needs "
              f"y < {y_ooeq:.3f}, capping m_nu at {m_nu_max:.3e} eV, "
              f"{nv.M_ATM_EV / m_nu_max:.0f}x below the atmospheric splitting. "
              f"O5's 11% agreement is an accident of two wrong factors, "
              f"not a near-miss: at M_N = 9.6e13 GeV the mechanism cannot "
              f"supply both the neutrino mass and the asymmetry.")
    out = {"experiment": "V106_w5_washout_audit", "pass": bool(ok),
           "status": status,
           "data": {"M_N_GeV": M_N, "H_GeV": H, "cases": cases,
                    "y_out_of_equilibrium_bound": y_ooeq,
                    "m_nu_max_eV": m_nu_max,
                    "fit": "Buchmuller-Di Bari-Plumacher, Ann. Phys. 315, "
                           "305 (2005), hep-ph/0401240, strong-washout "
                           "kappa_f = 2/(z_B K)(1-exp(-z_B K/2))"}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v106_w5_washout_audit.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    print(json.dumps(cases, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
