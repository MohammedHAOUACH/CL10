#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V38 — the numbers of the root **read as 0-parameter couplings** (C3):
the GUT coupling, the Dirac coupling, and the 1-loop running DOWN to M_Z (quantified).

**The gap (before).**  The GUT coupling ``alpha_10 = 1/25`` and the Dirac
coupling ``g_0 = 1/8`` were *chosen* rationals labelled "0-parameter".

**The 0-parameter selection rule.**  The 32-D module of the Cl(10) root is a
**2^5** space (the symplectic Pauli pattern: **5** qubits, 10 = 2x5
generators), so a **reading** of the root data gives the coupling:
$g_{10} = 1/n_{qubits} = 1/5$, $\alpha_{10} = g_{10}^2 = 1/25$ — the
*reading* is the audit's mapping convention, not a derivation: the clock
itself is **closed** by the closure axiom at $g_{10} = \pi/4$ (V131,
Theorem 6), and no 0-parameter **analytic** condition selects an angle —
the audit's Point 1).  The **Dirac** coupling (O2) is
$g_0 = 1/2^{n_{space}} = 1/2^3 = 1/8$, the **3** being the spatial
dimension of the (1,3) (the same 0-parameter number as the family count —
one number, two uses).

**The 1-loop running DOWN (quantified, the named residual).**  The QCA
lattice UV cutoff is the **Planck** scale with the 0-parameter UV input
``alpha_10 = 1/25``.  The standard 1-loop QFT extrapolation with the SM beta
coefficients **down** to $M_Z$ is reported with the deviation from the
measured ``alpha_i(M_Z)``: the weak coupling $a_2$ lands within a factor
$\\sim 1.6$, the strong coupling $a_3$ within an order of magnitude, and the
U(1) 1-loop extrapolation hits its pole **inside** the lattice (the
naive continuation below the pole is unphysical) — the **1-loop matching**
(O4) is the **named residual**, now **quantified**; the QCA **lattice RG**
is the open part.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import unification, uvrunning

TOL = 1e-9


def main():
    t0 = time.time()
    root = unification.g10_from_root()
    down = unification.running_down()
    meet = unification.pairwise_meeting_scales()

    a1, a2, a3 = down["couplings_at_M_Z"]
    m1, m2, m3 = down["measured_at_M_Z"]
    ok = (root["n_qubits"] == 5
          and root["g10_matches_sm_qca"]
          and root["alpha_10_matches_u10pred"]
          and root["n_space"] == 3
          and root["g0_matches_coupling"]
          and root["g10"] == 1.0 / 5.0
          and abs(root["alpha_10"] - 1.0 / 25.0) < 1e-12
          and abs(root["g0"] - 1.0 / 8.0) < 1e-12
          and down["a1_within_factor_2"]
          and down["su3_pole_inside_lattice"]
          # the paper's Section 13.1 statement: under pure SM running the three
          # pairs meet at three DIFFERENT scales, all above the input scale
          and all(v > unification.M_Z for v in meet.values())
          and len({round(v / max(meet.values()), 3) for v in meet.values()}) == 3)

    verdict = {
        "id": "V38_couplings_root",
        "tier": 3,
        "property": "P3",
        "residual": True,      # the 1-loop SM matching (O4) is the named residual, now quantified
        "pass": bool(ok),
        "status": "the_g10_and_g0_readings_are_root_numbers_the_1loop_down_running_is_the_named_residual_O4_quantified",
        "tol": TOL,
        "data": {
            "backend": root["backend"],
            "n_module": root["n_module"],
            "n_qubits": root["n_qubits"],
            "g10": root["g10"],
            "alpha_10": root["alpha_10"],
            "g10_matches_sm_qca": root["g10_matches_sm_qca"],
            "alpha_10_matches_u10pred": root["alpha_10_matches_u10pred"],
            "n_space": root["n_space"],
            "g0": root["g0"],
            "g0_matches_coupling": root["g0_matches_coupling"],
            "alpha_10_at_M_P": down["alpha_10_at_M_P"],
            "couplings_at_M_Z": down["couplings_at_M_Z"],
            "measured_at_M_Z": down["measured_at_M_Z"],
            "ratio_to_measured": down["ratio_to_measured"],
            "a1_within_factor_2": down["a1_within_factor_2"],
            "u1_up_pole_scale": down["u1_up_pole_scale"],
            "u1_pole_above_cutoff": down["u1_pole_above_cutoff"],
            "su3_down_pole_scale": down["su3_down_pole_scale"],
            "su3_pole_inside_lattice": down["su3_pole_inside_lattice"],
            "u1_pole_from_MZ_input_1loop_GeV": unification.landau_pole_u1(),
            "u1_pole_from_MZ_input_2loop_GeV": uvrunning.two_loop_pole_u1(),
            "pairwise_meeting_scales_sm_GeV": meet,
            "m_z_ref_GeV": unification.M_Z,
            "residual": down["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v38_couplings_root.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V38 - the numbers of the root read as 0-parameter couplings (C3)")
    print("=" * 68)
    print(f"  backend (GPU):  {root['backend']}")
    print(f"  (1) the GUT coupling from the root:  the module is 2^{root['n_qubits']} = {root['n_module']} (the 5 qubits)")
    print(f"      g_10 = 1/{root['n_qubits']} = {root['g10']:.4f}  ->  alpha_10 = {root['alpha_10']:.4f}  (= 1/25, matches V21/V23/V32)")
    print(f"  (2) the Dirac coupling from the root:  g_0 = 1/2^{root['n_space']} = 1/8 = {root['g0']:.4f}  (the 3 = the 3-space, matches O2)")
    print("  (3) the 1-loop running DOWN from alpha_10 at M_P to M_Z (the 0-parameter UV input, the named residual quantified):")
    print(f"      a1(M_Z) = {a1:.4f}  (measured {m1}, ratio {down['ratio_to_measured'][0]:.2f})")
    print(f"      a2(M_Z) = {a2:.4f}  (measured {m2}, ratio {down['ratio_to_measured'][1]:.2f}),  a3(M_Z) = {a3:.4f}  (measured {m3}; below its pole -> non-perturbative)")
    print(f"      the U(1) 1-loop pole of the UP-running:  {down['u1_up_pole_scale']:.2e} GeV (above the cutoff: {down['u1_pole_above_cutoff']})")
    print(f"      the SU(3) 1-loop pole of the DOWN-running:  {down['su3_down_pole_scale']:.2e} GeV (inside the lattice: {down['su3_pole_inside_lattice']})")
    print(f"      the U(1) poles of the perturbative QFT run UP from the measured M_Z input (Fig. 4): "
          f"1-loop {unification.landau_pole_u1():.2e} GeV, 2-loop {uvrunning.two_loop_pole_u1():.2e} GeV")
    print(f"      pure-SM pairwise meeting scales (1-loop from the measured M_Z inputs): "
          f"(1,2) {meet['12']:.3e}, (1,3) {meet['13']:.3e}, (2,3) {meet['23']:.3e} GeV")
    print(f"  named residual:  {down['residual']}")
    print(f"  VERDICT:  {'PASS (the readings g_10 = 1/5 and g_0 = 1/8 are root-data conventions; the clock is closed at pi/4; the 1-loop down-running is the quantified named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v38_couplings_root.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
