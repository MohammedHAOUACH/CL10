#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V33 — baryogenesis: the honest 0-parameter limit (the 3 Sakharov conditions).

**The 0-parameter result (the honest limit).**  The 0-parameter QCA step (the
45 Spin(10) plane rotations) is a **unitary** dynamics, and the 3 **Sakharov**
conditions are **not** met by it:
(1) **No baryon-number violation** in the dynamics -- the 45 **commute** with
the pseudoscalar ``P`` (they **preserve** the ``16/16-bar`` split, ``B-L`` is
**conserved** by the QCA step);
(2) **CPT-invariant** -- the ``16/16-bar`` are **CPT conjugates**;
(3) **in equilibrium** -- the QCA step is **unitary** (reversible).
So the **0-parameter** QCA **does not generate** the observed matter excess (a
**named residual**, not a wall).

**But the GUT structure contains the ingredients.**  The 45 **does** contain
the ``X/Y`` that mediate ``p -> e^+ pi^0`` (the 0-parameter falsifiable proton
decay, a **B-violating process**); the **heavy ``N``** (``M_R = M_GUT``) decay
is **out of equilibrium** (``Gamma/H ~ 6e-4`` -- a **kinematic** fact); and the
**PMNS CP phase** (V29) is the 0-parameter C/CP structure.  A **baryogenesis**
from these ingredients requires a **non-0-parameter** ingredient (the
out-of-equilibrium decay + the CP-violating rate) -- the **baryon asymmetry**
is the **named residual**.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import baryogenesis

TOL = 1e-9


def main():
    t0 = time.time()
    r = baryogenesis.verify(tol=TOL)
    sak = r["sakharov"]
    ing = r["ingredients"]

    ok = (r["ok"]
          and not sak["b_violation_met"]
          and sak["cpt_invariant"]
          and sak["in_equilibrium"]
          and ing["heavy_N_out_of_equilibrium"]
          and ing["proton_lifetime_yr"] > 1e34)

    verdict = {
        "id": "V33_baryogenesis",
        "tier": 2,
        "property": "P2",
        "residual": True,           # the baryon asymmetry is the named residual
        "pass": bool(ok),
        "status": "the_0_parameter_QCA_does_not_generate_the_baryon_asymmetry_the_GUT_ingredients_are_the_named_residual",
        "tol": TOL,
        "data": {
            "backend": sak["backend"],
            "b_violation_in_dynamics": sak["b_violation_met"],
            "max_45_pseudoscalar_commutator": sak["max_45_pseudoscalar_commutator"],
            "cpt_invariant": sak["cpt_invariant"],
            "in_equilibrium": sak["in_equilibrium"],
            "unitary_err": sak["unitary_err"],
            "xy_proton_decay_process": ing["xy_process"],
            "proton_lifetime_yr": ing["proton_lifetime_yr"],
            "heavy_N_out_of_equilibrium": ing["heavy_N_out_of_equilibrium"],
            "gamma_over_h": ing["gamma_over_h"],
            "pmns_cp_phase": ing["pmns_cp_phase"],
            "residual": "the baryon asymmetry (the 0-parameter QCA does not generate it; the GUT ingredients require a non-0-parameter out-of-equilibrium + CP-violating rate)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v33_baryogenesis.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V33 - baryogenesis: the honest 0-parameter limit (the 3 Sakharov conditions)")
    print("=" * 68)
    print(f"  backend (GPU):  {sak['backend']}")
    print(f"  (1) B-violation in the dynamics:  {'MET' if sak['b_violation_met'] else 'NOT MET'}  (the 45 preserve the 16/16-bar;  [T,P] max = {sak['max_45_pseudoscalar_commutator']:.1e})")
    print(f"  (2) CPT:  {'invariant' if sak['cpt_invariant'] else 'violated'}  (the 16/16-bar are CPT conjugates)")
    print(f"  (3) equilibrium:  {'IN equilibrium' if sak['in_equilibrium'] else 'out'}  (the QCA step is unitary/reversible)")
    print("  -> the 0-parameter QCA does NOT generate the baryon asymmetry (the NAMED RESIDUAL)")
    print(f"  the GUT structure contains the ingredients:  the X/Y (p -> e+ pi0;  tau ~ {ing['proton_lifetime_yr']:.1e} yr), the heavy-N out of equilibrium (Gamma/H = {ing['gamma_over_h']:.2e}), the PMNS CP phase")
    print(f"  VERDICT:  {'PASS (the 0-parameter QCA does not generate the baryon asymmetry; the GUT ingredients are the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v33_baryogenesis.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
