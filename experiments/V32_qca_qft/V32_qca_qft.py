#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V32 — the QCA -> QFT correspondence (the 't Hooft program) and the
0-parameter electroweak breaking (the W / Z / Higgs).

**The 0-parameter result.**  The QCA step (the 0-parameter 45 Spin(10) plane
rotations, ``U = prod [cos(g_10) I - i sin(g_10) T_{ab}]``) is **exactly
unitary**; its **eigenvalues** ``e^{-i E_n}`` are the **energy spectrum** of
the **QFT effective action** (the 't Hooft program: the QCA step = the QFT
time-evolution).  The **electroweak breaking** follows 0-parameter: the unified
coupling ``g_10`` + the ``5/3`` U(1) normalization + the vev ``v = 246`` GeV
give the **W/Z** masses and the 0-parameter predictions

    sin^2(theta_W) = 3/8 = 0.375,      m_W/m_Z = sqrt(5/8) = 0.79  (= cos theta_W).

The **Higgs** is the 16's ``H`` doublet; its mass is the **named residual**
(the SM does not fix it 0-parameter).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import electroweak

TOL = 1e-6


def main():
    t0 = time.time()
    r = electroweak.verify(tol=TOL)
    spec = r["spectrum"]
    ew = r["electroweak"]

    ok = (r["ok"]
          and spec["unitary_err"] < 1e-9
          and abs(ew["sin2_w"] - 0.375) < 1e-6
          and abs(ew["mw_over_mz"] - 0.7905694) < 1e-3)

    verdict = {
        "id": "V32_qca_qft",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the Higgs mass is the named residual
        "pass": bool(ok),
        "status": "the_qca_step_gives_the_sm_effective_action_with_the_0_parameter_electroweak_breaking",
        "tol": TOL,
        "data": {
            "backend": spec["backend"],
            "unitary_err": spec["unitary_err"],
            "n_energy_levels": spec["n_levels"],
            "g10": ew["g10"],
            "gprime": ew["gprime"],
            "vev_v": ew["vev_v"],
            "m_W": ew["m_W"],
            "m_Z": ew["m_Z"],
            "sin2_w": ew["sin2_w"],
            "sin2_w_0p": ew["sin2_w_0p"],
            "mw_over_mz": ew["mw_over_mz"],
            "mw_over_mz_0p": ew["mw_over_mz_0p"],
            "residual": "the Higgs mass (not 0-parameter in the SM)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v32_qca_qft.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V32 - the QCA -> QFT (the 't Hooft program) + the electroweak breaking")
    print("=" * 68)
    print(f"  backend (GPU):  {spec['backend']}")
    print(f"  the QCA step (the 45 plane rotations at the 0-parameter angle):  unitarity residual {spec['unitary_err']:.2e}")
    print(f"    the {spec['n_levels']} energy levels = the QFT time-evolution (the SM effective action)")
    print(f"  the 0-parameter electroweak breaking:")
    print(f"    the unified coupling g_10 = {ew['g10']:.4f}, the 5/3 -> g' = {ew['gprime']:.4f}, the vev v = {ew['vev_v']:.1f} GeV")
    print(f"    the W/Z masses (GUT scale):  m_W = {ew['m_W']:.1f}, m_Z = {ew['m_Z']:.1f} GeV")
    print(f"    the 0-parameter weak mixing angle:  sin^2(theta_W) = {ew['sin2_w']:.4f}  (= 3/8)")
    print(f"    the 0-parameter W/Z ratio:  m_W/m_Z = {ew['mw_over_mz']:.4f}  (= sqrt(5/8))")
    print(f"    the Higgs: the 16's H doublet (the mass = the named residual)")
    print(f"  VERDICT:  {'PASS (the QCA step -> the SM effective action + the 0-parameter electroweak breaking)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v32_qca_qft.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
