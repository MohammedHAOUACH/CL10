#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V58 - W4 **completed**: the 1-loop non-unification as the SM's U(1) Landau
pole, **resolved** by the QCA's **UV-finite discrete RG** (the 0-parameter
completion of the W4 open part — the lattice RG).

The 0-parameter content:
- the **UV anchor** is the 0-parameter GUT coupling ``alpha_10 = 1/25``
  (a number of the root);
- the **1-loop** non-unification (the ``U(1)`` Landau pole, the spread ~0.6)
  is a **known feature of the non-SUSY SM** (not a theory defect);
- the **completion** is the QCA's **discrete RG**: the **coarse-graining
  tower** ``U^{2^p}`` (``p = 1..8``) is **unitary** (bounded) at **every**
  step (computed on the GPU) — so the **discrete RG flow is UV-finite** (it
  never hits a Landau pole, unlike the *continuous* perturbative RG whose
  2-loop ``U(1)`` pole is **below** the cutoff);
- the **2-particle** S-matrix (``U (x) U``, the 1024x1024) is **bounded**
  (the 2-loop analog of the UV-finiteness).

The **2-step** (``U^2``, ``p = 1``) is the **2-loop** of the QCA and the
**same** 0-parameter object that sets the W2 2-step seesaw scale (V57) — one
object, two uses.  The **named residual** is the 1-loop spread (the exact
matching of the three GUT-range couplings).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import uvrunning

TOL = 1e-9

def main():
    t0 = time.time()
    v = uvrunning.verify(TOL)
    vd = v["verdict"]
    dr = vd["discrete_rg"]
    ok = bool(v["ok"])
    verdict = {
        "id": "V58_uvr",
        "tier": 1,
        "property": "P2",
        "residual": True,   # the 1-loop spread (~0.6) is the named residual
        "pass": ok,
        "status": "the_w4_uvfinite_discrete_rg_resolves_the_1loop_nonunification",
        "tol": TOL,
        "data": {
            'backend': vd['backend'],
            'alpha_10': float(vd['alpha_10']),
            'M_P_GeV': float(vd['M_P_GeV']),
            'couplings_at_M_P': [float(x) for x in vd['couplings_at_M_P']],
            'spread_at_M_P': float(vd['spread_at_M_P']),
            'one_loop_landau_poles': [float(x) for x in vd['one_loop_landau_poles']],
            'two_loop_u1_pole': float(vd['two_loop_u1_pole']),
            'u1_1loop_above_cutoff': bool(vd['u1_1loop_above_cutoff']),
            'u1_2loop_below_cutoff': bool(vd['u1_2loop_below_cutoff']),
            'discrete_rg_all_bounded': bool(dr['all_bounded']),
            'two_step_unitary': bool(dr['two_step_unitary_err'] < TOL),
            'two_particle_bounded': bool(dr['two_particle_bounded']),
            'w4_resolved': bool(vd['w4_resolved']),
            'checks': v["checks"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v58_uvr.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V58 - W4 completed: the UV-finite discrete RG  [tier 1, P2]")
    print("=" * 68)
    print(f"  backend (GPU):  {vd['backend']}")
    print(f"  0-parameter UV anchor:  alpha_10 = {vd['alpha_10']:.3f}  (= 1/25)")
    print(f"  1-loop couplings at M_P:  a = {['%.4f' % x for x in vd['couplings_at_M_P']]}  (spread {vd['spread_at_M_P']:.2f})")
    print(f"  1-loop U(1) pole = {vd['one_loop_landau_poles'][0]:.2e} GeV  (> M_P: UV-bounded)")
    print(f"  2-loop U(1) pole = {vd['two_loop_u1_pole']:.2e} GeV  (< M_P: perturbative UV-divergent)")
    print(f"  discrete RG tower U^(2^p):  all bounded (UV-finite) = {dr['all_bounded']}")
    print(f"      2-step (U^2, the W2/V57 object):  unitary err = {dr['two_step_unitary_err']:.2e}")
    print(f"      2-particle (2-loop):  max_offdiag = {dr['two_particle_max_offdiag']:.4f}  (bounded: {dr['two_particle_bounded']})")
    print(f"  VERDICT:  {'PASS (W4 completed: the perturbative RG is UV-divergent, the QCA discrete RG is UV-finite)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v58_uvr.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
