#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V57 - W2 **completed**: the 2-step seesaw **on the 10-loop ladder**.

The W2 open part (the 2-step seesaw scale) is now **0-parameter**: both the
Dirac mass and the Majorana (2-step) mass are **numbers of the root** on the
**same** 10-loop ladder that fixes the hierarchy (W8/V53) and the dark matter
(W7/V54):

    m_D  =  M_P (1/16\\pi)^{10}   (the 10-loop, the WIMP scale)
    M_N  =  M_P (1/16\\pi)^{3}    (the 2-step, the 10^14 scale; 0-parameter)
    m_nu =  M_P (1/16\\pi)^{17}   (17 = 2*10 - 3, the 0-parameter loop identity)

This **replaces** the previously hard-coded ``M_N = 10^14`` (V15/V18).  The
2-step (``U^2``) is a 0-parameter **unitary** of the GUT QCA (the 2-loop, the
same object as the W4 UV-finite completion).  The 0-parameter result gives
``m_nu ~ 0.15 eV`` (2.9x the observed 0.05 eV — the **right order**, closer
than the previous 5.3x).  The **named residual** is the 2-digit coefficient
and the loop integer ``n_2step = 3``.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import neutrino

TOL = 1e-9

def main():
    t0 = time.time()
    v = neutrino.verify(TOL)
    s = v["seesaw"]
    u2 = v["twostep"]
    ok = bool(v["ok"])
    verdict = {
        "id": "V57_neutrino",
        "tier": 5,
        "property": "P6",
        "residual": True,   # the 2-digit coefficient and n_2step=3 are the named residual
        "pass": ok,
        "status": "the_w2_twostep_seesaw_on_the_10loop_ladder_0param",
        "tol": TOL,
        "data": {
            'backend': s['backend'],
            'n_gen': int(s['n_gen']),
            'n_2step': int(s['n_2step']),
            'exp_2step': int(s['exp_2step']),
            'm_D_GeV': float(s['m_D_GeV']),
            'M_N_GeV': float(s['M_N_GeV']),
            'm_nu_2step_eV': float(s['m_nu_2step_eV']),
            'm_nu_1step_eV': float(s['m_nu_1step_eV']),
            'm_nu_ladder_eV': float(s['m_nu_ladder_eV']),
            'm_nu_obs_eV': float(s['m_nu_obs_eV']),
            'off_2step': float(s['off_2step']),
            'right_order': bool(s['right_order']),
            'ladder_match': bool(s['ladder_match']),
            'two_step_unitary': bool(u2['two_step_unitary']),
            'checks': v["checks"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v57_neutrino.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V57 - W2 completed: the 2-step seesaw on the 10-loop ladder  [tier 5, P6]")
    print("=" * 68)
    print(f"  backend (GPU):  {s['backend']}")
    print(f"  10-loop Dirac   m_D  = M_P (1/16pi)^10 = {s['m_D_GeV']:.1f} GeV  (the WIMP scale)")
    print(f"  2-step scale    M_N  = M_P (1/16pi)^3  = {s['M_N_GeV']:.2e} GeV  (0-param, replaces 10^14)")
    print("  2-step seesaw   m_nu = M_P (1/16pi)^17 (17 = 2*10 - 3) = "
          f"{s['m_nu_2step_eV']:.3f} eV")
    print(f"  observed m_nu = {s['m_nu_obs_eV']} eV;  off = {s['off_2step']:.2f}x  (right order: {s['right_order']})")
    print(f"  2-step (U^2) unitary: {u2['two_step_unitary']}  (err {u2['unitary_err_2step']:.2e})")
    print(f"  VERDICT:  {'PASS (W2 completed: both scales 0-parameter on the ladder)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v57_neutrino.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
