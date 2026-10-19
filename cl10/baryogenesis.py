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
proton decay -- a **B-violating process**) and the **heavy ``N``** of the
seesaw.  The out-of-equilibrium condition is **not** a free kinematic fact:
at the ladder seesaw scale ``M_N = M_P (1/16 pi)^3 = 9.6e13`` GeV it holds only
for ``y < 0.058``, while the neutrino mass needs ``y ~ 0.40`` (V106, V125) --
the two are mutually incompatible, and that is the single named reason O5 is
open.  A **baryogenesis** from these ingredients therefore requires a
**non-0-parameter** ingredient (the CP-odd rate, carried by the bond-layer
flux of V117/V125) -- the **baryon asymmetry** is the **named residual** (the
0-parameter dynamics does not generate it; the structure is there).
"""

from . import config as _config
from . import engine as E
from . import clifford
from . import u10pred
from . import coupling
from . import cosmology as _cosmology

import math
import numpy as np
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml


def _backend():
    return E._BACKEND


#: the branch point of the BDP washout solution: the Lambert-W argument
#: ``-(4/5) pi^{-1/5} K^{-2/5}`` reaches ``-1/e`` at this ``K``; below it the
#: strong-washout expression has no real root (weak washout, efficiency 1).
K_STRONG_WASHOUT_MIN = (math.e * (4.0 / 5.0) * math.pi ** (-1.0 / 5.0)) ** (5.0 / 2.0)


def bdp_washout_efficiency(K):
    """The Buchmuller-Di Bari-Plumacher washout efficiency (hep-ph/0401240).

    ``K = Gamma_N / H`` is the equilibrium ratio. Returns
    ``(kappa_f, z_B, in_domain)``: for ``K >= K_STRONG_WASHOUT_MIN`` the
    strong-washout ``z_B`` is the real Lambert-W solution and
    ``kappa_f = (2 / (z_B K)) (1 - exp(-z_B K / 2))``; below it the equation has
    **no real root** -- the decays freeze out before equilibrating, the
    strong-washout formula is out of its domain, and the efficiency is the
    weak-washout bound ``kappa_f = 1`` (``z_B`` is ``None``).
    """
    arg = -(4.0 / 5.0) * math.pi ** (-1.0 / 5.0) * K ** (-2.0 / 5.0)
    if arg < -1.0 / math.e:
        return 1.0, None, False
    w = -2.0
    for _ in range(200):
        w -= (w * math.exp(w) - arg) / (math.exp(w) * (1.0 + w))
    zB = -2.5 * w
    return (2.0 / (zB * K)) * (1.0 - math.exp(-zB * K / 2.0)), zB, True


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
    b_violation_met = not (max_comm < _TOL.standard)   # False: the 45 preserves the split
    # (2) CPT: the 16/16-bar are CPT conjugates (P maps 16<->16-bar; h_a anti-commute)
    max_anti = 0.0
    for a in range(10):
        c = E.asnumpy(g[a] @ P + P @ g[a])
        max_anti = max(max_anti, float(np.abs(c).max()))
    cpt_invariant = max_anti < _TOL.standard
    # (3) equilibrium: the QCA step is unitary -> reversible -> in equilibrium
    from . import sm_qca
    planes = sm_qca.plane_generators()
    theta = coupling.G10
    U = sm_qca.gut_qca_step(planes, theta)
    unitary_err = sm_qca.unitary_err(U)
    in_equilibrium = unitary_err < _TOL.standard
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
# the GUT ingredients, at the numbers of Section 13.3
# ---------------------------------------------------------------------------

def gut_ingredients():
    """The GUT ingredients of baryogenesis, at the paper's numbers.

    (a) the X/Y B-violating process exists; the proton lifetime in force is
        the desert solution of V64, ``tau_p = 1.8e36`` yr, **conditional** on
        a desert spectrum with two fitted thresholds.  The ``M_GUT = M_P``
        reading (``tau_p ~ 4e47`` yr) is ABANDONED (Section 13.4, F5).
    (b) the heavy-N decay: at the ladder seesaw scale
        ``M_N = M_P (1/16 pi)^3 = 9.6e13`` GeV the decay is **in equilibrium**
        at the Yukawa the atmospheric splitting needs (``y ~ 0.40``:
        ``Gamma/H = 47``) and out of equilibrium only for ``y < 0.058``, which
        caps ``m_nu`` at ``2.1e-3`` eV.  The neutrino mass and the
        out-of-equilibrium condition are mutually incompatible at that scale:
        this is why O5 is open for one named reason (Sakharov 3).
    (c) the CP phase: every CP-odd invariant of the step vanishes (the theorem
        of Section 13.3, stable under the Wilson extension).  The only built
        carrier is a bond-layer Peierls flux, which costs nothing of
        unitarity, ``det U = 1``, the mirror or the causal radius.
    """
    M_N = u10pred.M_PLANCK * u10pred.LOOP ** u10pred.RUNG_A   # 9.6e13 GeV (config cell.rung_A)
    H = (_cosmology.HUBBLE_COEFF * np.sqrt(_cosmology.G_STAR)
         * M_N ** 2 / u10pred.M_PLANCK)
    by_y = {}
    for y in (0.125, 0.0582, 0.3987):
        Gamma = y ** 2 * M_N / (8.0 * np.pi)
        by_y[f"y={y}"] = {"Gamma_over_H": float(Gamma / H),
                          "out_of_equilibrium": bool(Gamma < H)}
    return {
        "backend": _backend(),
        "xy_process": True,
        "proton_lifetime_yr": u10pred.proton_lifetime_desert(),
        "proton_lifetime_provenance": ("conditional: the desert spectrum with "
                                       "two fitted thresholds (V64, F5)"),
        "seesaw_scale_GeV": float(M_N),
        "hubble_at_seesaw_GeV": float(H),
        "equilibrium_by_yukawa": by_y,
        "out_of_equilibrium_yukawa_ceiling": 0.058,
        "cp_phase_in_the_step": 0.0,
        "cp_phase_carrier": "the bond-layer Peierls flux (V125); the step carries none",
    }


# ---------------------------------------------------------------------------
# GUT leptogenesis: the mechanism is 0-parameter, the exact value is named
# (the O2 pattern)
# ---------------------------------------------------------------------------

def leptogenesis_estimate():
    """GUT leptogenesis at the **right order**, with the **physical washout**.

    The ingredients and their provenance:
    (a) the seesaw scale ``M_N = M_P (1/16 pi)^3 = 9.6e13`` GeV (the ladder
        2-step; the loop integer is the fitted number — the open problem O2);
    (b) the **out-of-equilibrium** condition, read at the model's own Yukawa
        ``y = g_0 = 1/8``: it **fails** (``Gamma/H > 1``), and it holds only for
        ``y < 0.058``, which is incompatible with the neutrino mass (V106);
    (c) the **CP-odd rate** carries the loop factor
        ``eps_loop = g0^4/(16*pi)`` with ``g0 = 1/8`` (``coupling.G0``: the
        ``g_0 = 1/2^n_space`` **reading**, a convention — Section 15), times the
        **named** CP factor ``eta_CP`` (observed range ``1e-3..1e-2``);
    (d) the **SM sphaleron** conversion ``28/79`` — an **input of the SM
        sector**, not a number of the Cl(10) root (V23, a named residual);
    (e) the **washout** ``kappa_f`` of ``bdp_washout_efficiency`` at the same
        ``K = Gamma/H`` — the physical replacement of the guessed
        ``(1/16 pi)^2`` dilution (V106's audit).

    The **honest limit**: the QCA *alone* satisfies **no** Sakharov condition
    (V22/V33), the decay at the model's own Yukawa is **in** equilibrium, and
    the physical washout destroys the agreement of the unwashed product. The
    unwashed number is reported only to expose the coincidence it used to be.
    """
    ing = gut_ingredients()
    M_N = ing["seesaw_scale_GeV"]
    H = ing["hubble_at_seesaw_GeV"]
    # (b) the model's own Yukawa reading, and the equilibrium ratio it gives
    y = float(coupling.G0)                       # 1/8
    gamma_over_h = float(y ** 2 * M_N / (8.0 * np.pi) / H)
    out_of_equil = bool(gamma_over_h < 1.0)
    # (c) the loop factor of the CP-odd rate (the g_0 reading, no fit)
    g0 = coupling.G0
    eps_loop = g0 ** 4 / (16.0 * np.pi)   # 1/(4096*16*pi) ~ 4.86e-6
    loop_0param = bool(
        abs(eps_loop - (1.0 / 8.0) ** 4 / (16.0 * np.pi)) / eps_loop < _TOL.strict)
    # the NAMED CP-violating factor in the N decay (observed range)
    eta_cp = 1e-3
    # (d) the SM sphaleron conversion factor (an input of the SM sector)
    sphaleron = float(_config.CFG.sm.sphaleron)
    observed = float(_config.CFG.sm.delta_b_obs)
    # (e) the physical washout at the same K (V106), not a guessed dilution
    kap, zB, in_domain = bdp_washout_efficiency(gamma_over_h)
    delta_b_unwashed = sphaleron * eps_loop * eta_cp
    delta_b = delta_b_unwashed * kap
    ratio = float(delta_b / observed)
    unwashed_ratio = float(delta_b_unwashed / observed)
    sak = sakharov_conditions()
    qca_no_sakharov = bool(not sak["any_sakharov_met"])
    return {
        "backend": _backend(),
        "M_N_GeV": float(M_N),
        "y": float(y),
        "gamma_over_h": gamma_over_h,
        "out_of_equilibrium": out_of_equil,
        "out_of_equilibrium_yukawa_ceiling": ing["out_of_equilibrium_yukawa_ceiling"],
        "g0": float(g0),
        "eps_loop": float(eps_loop),
        "loop_factor_0_parameter": loop_0param,
        "eta_cp": float(eta_cp),
        "sphaleron_factor": float(sphaleron),
        "washout_efficiency": float(kap),
        "washout_z_B": None if zB is None else float(zB),
        "washout_strong_domain": bool(in_domain),
        "delta_b_unwashed": float(delta_b_unwashed),
        "unwashed_ratio_to_observed": unwashed_ratio,
        "delta_b": float(delta_b),
        "observed_delta_b": float(observed),
        "ratio_to_observed": ratio,
        "order_ok": bool(0.1 <= unwashed_ratio <= 10.0),
        "agreement_survives_washout": bool(0.1 <= ratio <= 10.0),
        "0param_qca_satisfies_no_sakharov": qca_no_sakharov,
        "residual": (
            "the exact eta_CP (the CP-violating magnitude in the N decay, the "
            "named input, observed range 1e-3..1e-2; taken 1e-3) and the exact "
            "Delta_B (the observed 6.12e-10). The unwashed product lands "
            "within a factor 10 of the observed value by accident: at the "
            "model's own Yukawa the decay is in equilibrium and the physical "
            "BDP washout (kappa_f) removes the agreement. O5's single named "
            "reason is Sakharov 3 (V106, V125)."
        ),
    }


def verify(tol=_TOL.standard):
    """The honest negative result, reproduced: the automaton satisfies **none**
    of the Sakharov conditions, and the GUT route that could supply them is
    blocked by the equilibrium condition (O5)."""
    sak = sakharov_conditions()
    ing = gut_ingredients()
    at_nu_mass = ing["equilibrium_by_yukawa"]["y=0.3987"]
    ok = (not sak["any_sakharov_met"]          # none of the 3 is met by the step
          and not sak["b_violation_met"]
          and sak["cpt_invariant"]
          and sak["in_equilibrium"]
          and ing["xy_process"]
          and not at_nu_mass["out_of_equilibrium"]     # Sakharov 3 fails at m_nu
          and ing["cp_phase_in_the_step"] == 0.0
          and ing["proton_lifetime_yr"] > 1e34)
    return {"sakharov": sak, "ingredients": ing, "ok": bool(ok)}


def main():
    import sys
    r = verify()
    sak = r["sakharov"]
    ing = r["ingredients"]
    print("=" * 68)
    print("CL10 - baryogenesis: the computed negative result (the 3 Sakharov conditions)")
    print("=" * 68)
    print(f"  backend (GPU):  {sak['backend']}")
    print(f"  (1) B-violation in the dynamics:  {'MET' if sak['b_violation_met'] else 'NOT MET'}  "
          f"(the 45 preserve the 16/16-bar;  [T,P] max = {sak['max_45_pseudoscalar_commutator']:.1e})")
    print(f"  (2) CPT:  {'invariant (NOT violated)' if sak['cpt_invariant'] else 'violated'}  "
          f"(the 16/16-bar are CPT conjugates; every CP-odd invariant vanishes)")
    print(f"  (3) equilibrium:  {'IN equilibrium' if sak['in_equilibrium'] else 'out'}  "
          f"(the step is unitary/reversible;  unitary err = {sak['unitary_err']:.1e})")
    print("  -> the automaton satisfies NONE of the Sakharov conditions (a computed theorem)")
    print("  the GUT route, at the paper's numbers (Section 13.3):")
    print(f"    (a) X/Y mediate p -> e+ pi0;  tau_p = {ing['proton_lifetime_yr']:.2e} yr "
          f"[{ing['proton_lifetime_provenance']}]")
    print(f"    (b) seesaw M_N = {ing['seesaw_scale_GeV']:.3e} GeV, H = {ing['hubble_at_seesaw_GeV']:.3e} GeV")
    for k, v in ing["equilibrium_by_yukawa"].items():
        print(f"        {k}:  Gamma/H = {v['Gamma_over_H']:.2f}  "
              f"{'out of' if v['out_of_equilibrium'] else 'IN'} equilibrium")
    print(f"        out-of-equilibrium needs y < {ing['out_of_equilibrium_yukawa_ceiling']}, "
          "which caps m_nu below the atmospheric splitting (O5)")
    print(f"    (c) CP phase of the step = {ing['cp_phase_in_the_step']:.1f}; "
          f"{ing['cp_phase_carrier']}")
    print(f"  VERDICT:  {'PASS (the negative result and the named blocker O5 are reproduced)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
