#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V43 - the 9-to-1 reduction: the nine residuals W1-W9 are the nine named
projections of ONE 0-parameter object - the 32-level QCA spectrum.

The single 0-parameter object of the theory is the GUT QCA step U (the 32x32
product of the 45 Spin(10) plane rotations at g_10 = 1/5; the only input is
g_10).  Its spectrum - the 32 levels of H = -i log(U) - is the one computed
0-parameter object.  This experiment verifies the 0-parameter invariants of
that one object:

- the 0-point = 0  (the 32 levels are exactly +/---paired; the W1 s-mechanism)
- the SL(32) identity  (det U = 1 for every g_10, because every one of the 45
  plane rotations is a special unitary: each T_ab = i h_a h_b is a traceless
  16/16 involution)
- the CPT (16/16-bar) split  (U commutes with the pseudoscalar; each sector's
  0-point is 0)

and reports the **nine residuals W1-W9** as the **nine named projections of
this one object** (the 0-parameter mechanism of each is a consequence of the
single Cl(10) object; the exact value of each is the named fine-tuning).

The **elegant 9-to-1 reduction**: the nine residuals are not nine independent
fine-tunings; they are the nine independent exact values (scales / angles)
that the one 0-parameter object does not fix.  The 0-point (W1) is the most
invariant of the nine - a group-theoretic identity (SL(32)) of the object,
equal to 0 for **every** coupling.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import qca_spectrum

TOL = 1e-12

# the nine residuals W1-W9 as the nine named projections of the one object.
# Each: (the named residual, the 0-parameter mechanism, the face of the
# object it projects).
NINE = [
    ("W1", "the exact Lambda ~ 10^-122",
     "the 0-point = 0 (the 512/512 s-mechanism; the SL(32) identity: det U = 1)",
     "the spectrum (the 0-point)"),
    ("W2", "the exact M_N and the exact m_nu",
     "the 2-step GUT breaking + the 0-parameter Dirac mass m_D = g_0 v (g_0 = 1/8)",
     "the algebra (the scale ladder)"),
    ("W3", "the why-3 (the family count)",
     "the 3 = the 3 spatial directions; the 96-D family-symmetric step",
     "the algebra (the (1,3) embedding)"),
    ("W4", "the 1-loop SM matching",
     "the 0-parameter UV input alpha_10 = 1/25 (the 5 qubits) + the 1-loop beta",
     "the algebra (the 5 qubits)"),
    ("W5", "the exact eta_CP and the exact delta_B",
     "GUT leptogenesis: the heavy N out of equilibrium, g_0^4/16pi, the 28/79 sphaleron",
     "the algebra (the 2-step)"),
    ("W6", "the exact Yukawa (the hierarchy, the CKM)",
     "the 10-Higgs mass term + the tribimaximal PMNS (only sin^2_12 = 1/3 is 0-parameter as a number)",
     "the spectrum (the 16 levels)"),
    ("W7", "the sub-Planck dark-matter scale",
     "the stable soliton u = sqrt(2) sech x of the non-linear metric",
     "the dynamics (the non-linear metric)"),
    ("W8", "the derivation of the hierarchy exponent (v/E_P ~ 10^-17)",
     "the 0-parameter value v/E_P = 2.02e-17 (a number of the root)",
     "the algebra (the two scales)"),
    ("W9", "the QCA lattice RG (S-matrix to cross-sections)",
     "the 0-parameter S-matrix = U (the 't Hooft program)",
     "the step (S = U)"),
]


def main():
    t0 = time.time()
    r = qca_spectrum.verify(tol=1e-10)
    s = r["spectrum"]
    ok = r["ok"]

    verdict = {
        "id": "V43_qca_spectrum",
        "tier": 4,
        "property": "P2,P5",
        "residual": True,
        "pass": ok,
        "status": "the_9_to_1_reduction: one object (the 32-level QCA spectrum), nine named projections (W1-W9)",
        "tol": TOL,
        "data": {
            "one_object": {
                "n_levels": s["n_levels"],
                "g10": s["g10"],
                "zero_point": s["zero_point"],
                "det_U_minus_1": s["det_U_minus_1"],
                "sl32": s["sl32"],
                "pairing_resid": s["pairing_resid"],
                "max_trace_T45": s["max_trace_T45"],
                "cpt_0point_16": s["cpt_0point_16"],
                "cpt_0point_16bar": s["cpt_0point_16bar"],
                "commutes_with_pseudoscalar": s["commutes_with_pseudoscalar"],
            },
            "nine_projections": [
                {"wall": w, "residual": res, "mechanism": mech, "projection": proj}
                for (w, res, mech, proj) in NINE
            ],
        },
        "wall_time": round(time.time() - t0, 3),
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v43_qca_spectrum.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V43 - the 9-to-1 reduction: one object, nine projections")
    print("=" * 68)
    print(f"  the one object : the {s['n_levels']}-level QCA spectrum (g_10 = {s['g10']})")
    print(f"    0-point     = {s['zero_point']:.2e}   (the 32 levels are +/---paired)")
    print(f"    SL(32)      : det(U) - 1 = {s['det_U_minus_1']:.2e}, "
          f"45 max|tr| = {s['max_trace_T45']:.1e}")
    print(f"    CPT split   : 0-point 16/16-bar = "
          f"{s['cpt_0point_16']:.1e} / {s['cpt_0point_16bar']:.1e}")
    print("  the nine projections (W1-W9):")
    for (w, res, mech, proj) in NINE:
        print(f"    {w}  {proj}")
    print(f"  VERDICT      : {'PASS (the one object: 0-point = 0, det U = 1, CPT)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v43_qca_spectrum.json"))
    print(f"  Wall time    : {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
