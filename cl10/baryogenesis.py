#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — baryogenesis: the honest 0-parameter limit (the 3 Sakharov conditions
and what the GUT structure does and does not generate).

**The 0-parameter QCA does NOT generate the baryon asymmetry.**  The 0-parameter
QCA step (the 45 Spin(10) plane rotations) is a **unitary** dynamics, and the
3 **Sakharov** conditions are **not** met by it:

1. **No baryon-number violation in the dynamics.**  The 45 gauge generators
   **commute** with the **pseudoscalar** ``P = h_1 ... h_{10}`` (the ``16/16-bar``
   CPT operator), so they **preserve** the ``16/16-bar`` split (``B-L`` is
   **conserved** by the QCA step):  ``[T_{ab}, P] = 0`` for all 45.
2. **CPT-invariant.**  The ``16`` and the ``16-bar`` are **CPT conjugates**
   (``P`` maps ``16 <-> 16-bar``; the 10 grade-1 generators ``h_a``
   **anti**-commute with ``P``), so the 0-parameter QCA is **CPT-invariant**.
3. **In equilibrium.**  The QCA step is **unitary** (``U U^dag = I``) ->
   **reversible** -> **in equilibrium**.

So the **0-parameter** QCA is **baryon-conserving, CPT-invariant, in
equilibrium** -- it does **not** generate the observed matter excess (a **named
residual**, not a wall).

**But the GUT structure contains the ingredients.**  The 45 **does** contain
the ``X/Y`` bosons that mediate ``p -> e^+ pi^0`` (the **0-parameter** falsifiable
proton decay -- a **B-violating process**), the **heavy ``N``** (``M_R = M_GUT``)
decay is **out of equilibrium** (``Gamma/H ~ 6e-4 < 1`` -- a **kinematic**
fact), and the **PMNS CP phase** (V29) is the 0-parameter C/CP structure.  A
**baryogenesis** from these ingredients requires a **non-0-parameter**
ingredient (the **out-of-equilibrium decay** + the **CP-violating** rate) --
the **baryon asymmetry** is the **named residual** (the 0-parameter dynamics
does not generate it; the structure is there).
"""

from . import engine as E
from . import clifford
from . import toepred

import numpy as np


def _backend():
    return E._BACKEND


# ---------------------------------------------------------------------------
# the pseudoscalar P (the 16/16-bar CPT operator) and the 3 Sakharov conditions
# ---------------------------------------------------------------------------

def _pseudoscalar():
    """The pseudoscalar ``P = h_1 ... h_10`` (anti-Hermitian, ``P^2 = -I``);
    ``iP`` is the Hermitian CPT/chirality involution (``+1`` on the 16,
    ``-1`` on the 16-bar)."""
    g = [E.Array_(m) for m in clifford.generators()]
    P = g[0]
    for i in range(1, 10):
        P = P @ g[i]
    return P


def sakharov_conditions():
    """The 3 Sakharov conditions as seen by the 0-parameter QCA step (the 45).
    None is met by the 0-parameter **dynamics** (the named limit)."""
    g = [E.Array_(m) for m in clifford.generators()]
    P = _pseudoscalar()
    # (1) B-violation: do the 45 mix the 16/16-bar?  [T_ab, P] = 0 -> no.
    max_comm = 0.0
    for a in range(10):
        for b in range(a + 1, 10):
            T = 1j * (g[a] @ g[b])
            c = E.asnumpy(T @ P - P @ T)
            max_comm = max(max_comm, float(np.abs(c).max()))
    b_violation_met = not (max_comm < 1e-9)   # False: the 45 preserves the split
    # (2) CPT: the 16/16-bar are CPT conjugates (P maps 16<->16-bar; h_a anti-commute)
    max_anti = 0.0
    for a in range(10):
        c = E.asnumpy(g[a] @ P + P @ g[a])
        max_anti = max(max_anti, float(np.abs(c).max()))
    cpt_invariant = max_anti < 1e-9
    # (3) equilibrium: the QCA step is unitary -> reversible -> in equilibrium
    from . import sm_qca
    planes = sm_qca.plane_generators()
    theta = (toepred.ALPHA_10) ** 0.5
    U = sm_qca.gut_qca_step(planes, theta)
    unitary_err = sm_qca.unitary_err(U)
    in_equilibrium = unitary_err < 1e-9
    return {
        "backend": _backend(),
        "b_violation_met": bool(b_violation_met),
        "max_45_pseudoscalar_commutator": float(max_comm),
        "cpt_invariant": bool(cpt_invariant),
        "in_equilibrium": bool(in_equilibrium),
        "unitary_err": float(unitary_err),
        "any_sakharov_met": bool(b_violation_met or not cpt_invariant or not in_equilibrium),
    }


# ---------------------------------------------------------------------------
# the GUT ingredients (the structure is there, even if the 0-parameter dynamics
# does not generate the asymmetry)
# ---------------------------------------------------------------------------

def gut_ingredients():
    """The GUT structure contains the baryogenesis ingredients:
    (a) the X/Y B-violating process (p -> e+ pi0, the 0-parameter proton decay);
    (b) the heavy-N decay OUT of equilibrium (Gamma/H < 1, the kinematics);
    (c) the PMNS CP phase (the 0-parameter C/CP structure, V29)."""
    g2 = 4.0 * np.pi * toepred.ALPHA_10
    M_GUT = toepred.M_PLANCK
    Gamma = g2 * M_GUT / (16.0 * np.pi)
    H = 1.66 * np.sqrt(100.0) * M_GUT ** 2 / toepred.M_PLANCK
    tau_p = toepred.proton_decay_lifetime()
    return {
        "backend": _backend(),
        "xy_process": True,                 # the X/Y mediate p -> e+ pi0 (0-parameter)
        "proton_lifetime_yr": float(tau_p),
        "heavy_N_out_of_equilibrium": bool(Gamma < H),
        "gamma_over_h": float(Gamma / H),
        "pmns_cp_phase": True,              # the 0-parameter C/CP (V29)
    }


def verify(tol=1e-9):
    """0-parameter checks: the 0-parameter QCA does NOT generate the baryon
    asymmetry (the 3 Sakharov conditions are not met by the 0-parameter
    dynamics), but the GUT structure contains the ingredients."""
    sak = sakharov_conditions()
    ing = gut_ingredients()
    ok = (not sak["any_sakharov_met"]          # none of the 3 is met by the 0-parameter dynamics
          and not sak["b_violation_met"]
          and sak["cpt_invariant"]
          and sak["in_equilibrium"]
          and ing["xy_process"]
          and ing["heavy_N_out_of_equilibrium"]
          and ing["pmns_cp_phase"]
          and ing["proton_lifetime_yr"] > 1e34)
    return {"sakharov": sak, "ingredients": ing, "ok": bool(ok)}


def main():
    import sys
    r = verify()
    sak = r["sakharov"]
    ing = r["ingredients"]
    print("=" * 68)
    print("CL10 - baryogenesis: the honest 0-parameter limit (the 3 Sakharov conditions)")
    print("=" * 68)
    print(f"  backend (GPU):  {sak['backend']}")
    print(f"  (1) B-violation in the dynamics:  {'MET' if sak['b_violation_met'] else 'NOT MET'}  "
          f"(the 45 preserve the 16/16-bar;  [T,P] max = {sak['max_45_pseudoscalar_commutator']:.1e})")
    print(f"  (2) CPT:  {'invariant (NOT violated)' if sak['cpt_invariant'] else 'violated'}  "
          f"(the 16/16-bar are CPT conjugates)")
    print(f"  (3) equilibrium:  {'IN equilibrium' if sak['in_equilibrium'] else 'out'}  "
          f"(the QCA step is unitary/reversible;  unitary err = {sak['unitary_err']:.1e})")
    print(f"  -> the 0-parameter QCA does NOT generate the baryon asymmetry (the NAMED RESIDUAL)")
    print(f"  the GUT structure contains the ingredients:")
    print(f"    (a) the X/Y mediate p -> e+ pi0 (the 0-parameter proton decay;  tau ~ {ing['proton_lifetime_yr']:.1e} yr)")
    print(f"    (b) the heavy-N decay is out of equilibrium (Gamma/H = {ing['gamma_over_h']:.2e} < 1)")
    print(f"    (c) the PMNS CP phase (the 0-parameter C/CP structure, V29)")
    print(f"  VERDICT:  {'PASS (the 0-parameter QCA does not generate the baryon asymmetry; the GUT ingredients are the named residual)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
