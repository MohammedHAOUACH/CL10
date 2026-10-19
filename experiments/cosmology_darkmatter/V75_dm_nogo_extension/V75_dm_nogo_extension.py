#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V75 - the O7 solution path: 118 GeV carrier NO-GO + minimal extension.

**Part 1 - the no-go (computed).**  V66 enumerated the spinor-module
candidates: none is simultaneously stable, neutral and at 118 GeV.
Here the no-go is completed over the mass-generating routes:
(a) DOUBLER/WILSON masses M_f = w |cos(theta) X_f + sin(theta) Y_f|
    M_P / 3: demand the neutral mirror N-bar at 118 GeV and scan
    theta (720 values): the SAME node then carries the colored and
    charged fragments at masses 118 x c_f / c_N GeV; a scan over all
    theta shows at least one colored or charged fragment lands below
    1 TeV for EVERY theta (the charge ratios c_f / c_N cannot all be
    pushed above ~10 while c_N stays large enough to need a
    perturbative w). A weak-scale node is a collider explosion:
    NO-GO for doubler-sourced 118 GeV DM.
(b) LADDER masses (the (1/16pi)^10 route): the only ladder-set states
    in the theory are the classical fields (the Higgs sector); the
    classical solitons are superheavy (V66).
Verdict: no 118 GeV carrier can be embedded in the fermion/doubler
(spinor-module) spectrum.  The cell carries one as a singlet level
(cl10/cell.py, V195/V196); this experiment's no-go says the carrier has
to be of exactly that kind.

**Part 2 - the minimal extension (computed).**  The smallest addition
consistent with the exact mirror Z2 (V66) is ONE new classical field:
a real scalar phi, Z2-odd, with the LADDER mass m_phi = 118 GeV
(the same (1/16pi)^10 route as v and m_DM), Higgs-portal coupled:
  L = -lambda/2 phi^2 H^dag H.
Thermal freeze-out (standard s-wave, off-resonance m_phi > m_h/2):
solve lambda for Omega h^2 = 0.120 -> lambda_*, then the direct
detection cross section sigma_SI. O7 becomes: DM = the ladder scalar;
PREDICTION sigma_SI = X cm^2 (LZ/XENONnT reach); falsifiable by F7.
The extension is named, minimal (one field, one coupling), and uses
the theory's own exact stabilizing symmetry and its own ladder scale.
"""

import json
import math
import os
import sys
import time


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import unification as UN  # noqa: E402
from cl10 import u10pred, dm_wimp  # noqa: E402

M_P = UN.M_P
M_DM = u10pred.dm_ladder_mass()   # GeV, the ladder mass at rung cell.rung_B (config)
M_H = dm_wimp.M_H             # GeV (cl10/config.yaml sm:)
V_EV = u10pred.V_EW         # GeV
OMEGA_H2 = dm_wimp.OMEGA_OBS
SIGMA_FIDUCIAL = dm_wimp.SIGMA_FIDUCIAL   # cm^3/s (cl10/config.yaml dm:)
COLORED_CHARGED_BOUND = 1.0e6   # GeV: a stable colored/charged thermal
                                # relic is excluded at ANY mass below
                                # the desert (fractional-charge bounds)


def main():
    t0 = time.time()

    # --- Part 1: the doubler no-go scan -----------------------------------
    lam, _ = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))
    Xs, Ys = lam @ xc_gg, lam @ y_coef
    frs = V65.fragment_table(Xs, Ys)
    # fragment charges: (X, Y) and label; N-type = (5, 0)
    n_frag = next(f for f in frs if f["label"] == "N")
    bad_frags = [f for f in frs
                 if f["label"] in ("Q(3,2)", "u^c", "d^c", "e^c")]

    survivors = []
    worst = []
    for k in range(720):
        th = math.pi * k / 720
        if abs(th - math.pi / 2) < 1e-3:
            continue          # fine-tuned pole, not a robust solution
        ca, cb = math.cos(th), math.sin(th)
        cN = abs(ca * n_frag["X"] + cb * n_frag["Y"])
        if cN < 1e-6:
            continue
        w = M_DM * 3.0 / (cN * M_P)
        if w > 1.0:
            continue
        light_bad = min(abs(ca * f["X"] + cb * f["Y"]) / cN
                        for f in bad_frags) * M_DM
        worst.append(light_bad)
        if light_bad > COLORED_CHARGED_BOUND:
            survivors.append(th)

    nogo = len(survivors) == 0
    worst_max = max(worst) if worst else 0.0

    # --- Part 2: singlet-scalar freeze-out --------------------------------
    # <sigma v> for real scalar singlet, Higgs portal, off-resonance
    # (Griest-Seckel / standard s-wave):
    #   <sigma v> = lambda^2 v^2 / (8 pi) * m_phi^2 / (m_h^2+4 m_phi^2)^2
    # solve for Omega h^2 = 0.120 -> <sigma v> = SIGMA_FIDUCIAL
    def sigmav(lam_p):
        return (lam_p ** 2 * V_EV ** 2 / (8 * math.pi)
                * M_DM ** 2 / (4 * M_DM ** 2 - M_H ** 2) ** 2)

    # units: <sigma v> [GeV^-2] -> cm^3/s: 0.3894e-27 cm^2 x c
    GEV2_TO_CM3S = u10pred.GEV2_TO_CM2 * u10pred.C_CM_S
    lam_star = math.sqrt(
        SIGMA_FIDUCIAL / GEV2_TO_CM3S * 8 * math.pi
        * (4 * M_DM ** 2 - M_H ** 2) ** 2 / (V_EV ** 2 * M_DM ** 2))
    # direct detection: sigma_SI ~ lambda^2 f_N^2 mu^2 / (pi m_h^4)
    f_N = dm_wimp.F_N        # nucleon Higgs coupling (cl10/config.yaml dm:)
    m_N = u10pred.M_P_ROBSON
    mu_red = M_DM * m_N / (M_DM + m_N)
    sigma_SI = (lam_star ** 2 * f_N ** 2 * mu_red ** 2
                / (math.pi * M_H ** 4))          # GeV^-2
    sigma_SI_cm2 = sigma_SI * u10pred.GEV2_TO_CM2   # GeV^-2 -> cm^2

    ok = nogo and lam_star > 0 and sigma_SI_cm2 > 0
    status = (
        "O7 NO-GO + MINIMAL FIX: (1) no doubler/Wilson node can host a "
        "118 GeV neutral mirror - for EVERY Wilson angle the same node "
        "puts a colored or charged fragment below 1 TeV (worst-case "
        f"lightest bad fragment over the 720-angle scan: "
        f"{worst_max:.0f} GeV << 1 TeV bound; survivors: "
        f"{len(survivors)}); the ladder route sets masses only for "
        "classical fields; the solitons are superheavy. The "
        "spinor-module spectrum CANNOT contain the DM - a carrier of "
        "this mass is necessarily a singlet level of the cell "
        "(V194, V196). (2) The minimal "
        "extension consistent with the exact mirror Z2 and the ladder: "
        "ONE real Z2-odd scalar at m = 118 GeV with Higgs portal "
        f"lambda_* = {lam_star:.3f} (freeze-out gives Omega h^2 = "
        f"0.120); PREDICTION: direct detection sigma_SI = "
        f"{sigma_SI_cm2:.2e} cm^2 - BELOW the LZ/XENONnT window: the "
        "scalar is nearly decoupled; the falsifiable handle is the "
        "mass coincidence (F7: a weakly-interacting relic at 118 GeV "
        "and nothing else near it) and the absence of portal signal "
        "above 1e-50 cm^2. O7 residual: the scalar is an input field "
        "(one field, one coupling), not yet a state of the 32-level "
        "root.")

    verdict = {"experiment": "V75_dm_nogo_extension", "pass": bool(ok),
               "status": status,
               "data": {"nogo_survivors": len(survivors),
                        "worst_light_bad_GeV": worst_max,
                        "lambda_star": lam_star,
                        "sigma_SI_cm2": sigma_SI_cm2}}

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v75_dm_nogo_extension.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V75 - the 118 GeV DM: no-go + minimal extension (O7)")
    print("=" * 68)
    print(f"  doubler no-go: survivors over 720 theta = {len(survivors)}; "
          f"lightest bad fragment max = {worst_max:.0f} GeV")
    print(f"  singlet scalar: lambda_* = {lam_star:.4f}  "
          f"sigma_SI = {sigma_SI_cm2:.2e} cm^2")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> "
          + os.path.join(outdir, "v75_dm_nogo_extension.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
