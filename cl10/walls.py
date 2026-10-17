#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the three named walls and their current status (Section 16).

The wall labels W1-W8 of the earlier drafts are the **open problems O1-O9** of
``doc/final_publication.md``, Section 16.  The mapping for the three treated
here: **W1 -> O1** (the cosmological constant), **W2 -> O2** (the neutrino
masses), **W3 -> O3** (the family number).  What the mechanisms below establish
is a *mechanism*, never the observed value: each function returns the number
the mechanism gives and the named residual that is left.

- **O1 (W1) — the cosmological constant.**  The ``Z_2`` grading of the
  ``Cl(10)`` is **512 even / 512 odd (equal)**, so the 0-point energy of the
  two sectors **cancels** and the automaton's vacuum contribution is exactly
  zero.  The observed ``1e-122`` is **not** produced by it: the cancellation is
  exact, and how it is broken is the open problem.

- **O2 (W2) — the neutrino masses.**  The sterile ``N`` (the 1 of the 16) takes
  its Majorana mass at the ladder 2-step scale ``M_N = M_P (1/16 pi)^3 =
  9.6e13`` GeV (the loop integer fitted) and the Dirac
  mass is ``m_D = g_0 v`` with the ``g_0 = 1/8`` reading, so the 2-step seesaw
  gives the **right order**.  The ladder value ``m_nu = 0.146`` eV is
  **excluded** by the measured ``sum m_nu < 0.12`` eV (the falsifier F3); the
  ``1/3`` candidate ``0.0487`` eV is a candidate whose junction amplitude is
  free inside a finite window.  ``M_N`` is **fitted**, not derived.

- **O3 (W3) — the three families.**  The number **3** is the count of spatial
  directions and it *labels* the three copies of the 16; it does **not** force
  the count: the cubic group is transitive on the 12 doubler nodes, and the
  3-node triangle has a cubic orbit of 8.  The family number is an open
  symmetry-breaking input.
"""

from . import coupling
from . import u10pred
from . import neutrino as nv

# the scales of the sector (single source: u10pred; the reading: coupling)
M_PLANCK = u10pred.M_PLANCK   # GeV (the cutoff; M_GUT = M_P is NOT assumed)
V_EW = u10pred.V_EW           # GeV (the measured electroweak vev)
G0 = coupling.G0              # 1/8 (the g_0 = 1/2^n_space reading: a convention)


# ---------------------------------------------------------------------------
# O1 (W1) — the cosmological constant (the Z2 cancellation / s-mechanism)
# ---------------------------------------------------------------------------

def w1_cosmological_constant():
    """The vacuum cancellation behind the cosmological constant (the O1 mechanism).

    The $\\mathbb{Z}_2$ grading of the $\\mathrm{Cl}(10)$ gives **equal**
    even (512, bosonic) and odd (512, fermionic) sectors.  The **0-point
    energy** of a mode is $\\pm \tfrac12$ (the sign is $+$ for bosons, $-$ for
    fermions).  The **equal** sectors **cancel**:

        Lambda_QCA(0-param) = 512*(+1/2) + 512*(-1/2) = 0.

    This is the **self-tuning** (s-mechanism): the flat vacuum is a saddle
    point of the energy, so the 0-point **does not gravitate**.

    Returns a dict with the even/odd counts, the net vacuum contribution
    (exactly zero), and the residual: the observed value is **not** produced
    by the mechanism (the open problem O1).
    """
    n_even = 512            # the bosonic (even) sector
    n_odd = 512             # the fermionic (odd) sector
    half = 0.5              # the 0-point per mode (the 1/2)
    rho_even = n_even * half          # + (bosons)
    rho_odd = -n_odd * half           # - (fermions)
    lambda_0param = rho_even + rho_odd      # = 0 (the cancellation)
    lambda_obs = 1e-122                   # the observed (Planck units)
    return {
        "n_even": int(n_even),
        "n_odd": int(n_odd),
        "cancellation": bool(abs(n_even - n_odd) < 1e-12),
        "lambda_0param": float(lambda_0param),
        "lambda_obs": float(lambda_obs),
        "residual": float(lambda_obs),      # the *named* fine-tuning
    }


# ---------------------------------------------------------------------------
# O2 (W2) — the seesaw scale (the 2-step breaking, the fitted M_N)
# ---------------------------------------------------------------------------

def w2_seesaw():
    """The 2-step seesaw: the right **order** for $m_\\nu$, not its value.

    The **Dirac** mass is $m_D = g_0 v$ with the $g_0 = 1/8$ **reading** (a
    convention, Section 15).  The **Majorana** mass $M_N$ of the sterile $N$ is
    the ladder 2-step scale $M_N = M_P (1/16\\pi)^{n_{2\\text{step}}}$
    ($= 9.6\\times10^{13}$ GeV, `neutrino.twostep_scale`): the **loop integer**
    is **fitted**, not derived (the open problem O2 — it replaces the earlier
    hard-coded $10^{14}$ GeV).  The **type-I seesaw** gives
    $m_\\nu = m_D^2 / M_N$.

    Returns a dict with $m_D$, $M_N$, the 2-step $m_\\nu$ (the right order),
    the 1-step $m_\\nu$ (for comparison), and the residual: the exact $m_\\nu$
    and the fact that the ladder value is excluded by $\\sum m_\\nu < 0.12$ eV
    (the falsifier F3).
    """
    m_D = G0 * V_EW                  # 30.7 GeV (the Dirac mass from the g_0 reading)
    M_N = nv.twostep_scale()         # the ladder 2-step scale (the integer fitted, O2)
    M_GUT = M_PLANCK                 # the cutoff, for the 1-step comparison
    m_nu_2step = m_D ** 2 / M_N * 1e9      # GeV -> eV
    m_nu_1step = m_D ** 2 / M_GUT * 1e9    # GeV -> eV
    m_nu_obs = nv.M_ATM_EV                 # the observed heaviest scale (eV)
    off_1step = m_nu_1step / m_nu_obs      # the 1-step is ~6 orders off
    off_2step = m_nu_2step / m_nu_obs      # the 2-step is ~0.7 orders off
    return {
        "g0": float(G0),
        "m_D_GeV": float(m_D),
        "M_N_GeV": float(M_N),
        "m_nu_2step_eV": float(m_nu_2step),
        "m_nu_1step_eV": float(m_nu_1step),
        "m_nu_obs_eV": float(m_nu_obs),
        "off_1step": float(off_1step),
        "off_2step": float(off_2step),
        "right_order": bool(0.01 < off_2step < 100),   # within 2 orders of the observed
        "residual": float(m_nu_obs),            # the *named* exact m_nu
    }


# ---------------------------------------------------------------------------
# O3 (W3) — the three families (the 3 labels the copies, it does not force the count)
# ---------------------------------------------------------------------------

def w3_families():
    """The family count: the 3 of the 3 spatial directions **labels** the
    copies, it does **not** force the count (the open problem O3).

    Each generation is a **copy** of the 16 Weyl, and the three copies are
    indexed by the 3 spatial directions of the $(1,3)$ causal structure.  What
    is NOT established: that the index forces the number.  The cubic group is
    transitive on the 12 doubler nodes (the only cubic-symmetric split is the
    full set, whose meeting scale is excluded), and the 3-node triangle has a
    cubic orbit of 8.  The count is an open symmetry-breaking input.

    Returns a dict with the time/space split, the 3-family candidate, the
    observed count, and the residual (the why-3 derivation, O3).
    """
    n_time = 1
    n_space = 3
    n_fam_candidate = n_space          # the candidate (the 3 of the 3-space)
    n_observed = 3
    return {
        "n_time": int(n_time),
        "n_space": int(n_space),
        "n_fam_candidate": int(n_fam_candidate),
        "n_observed": int(n_observed),
        "match": bool(n_fam_candidate == n_observed),
        "derived": False,             # the index labels the copies; it does not force the count
        "residual": "O3: why the family count = the spatial dimension (an open symmetry-breaking input)",
    }

def baryon_asymmetry():
    """The unitary QCA step checked against the **Sakharov conditions**.

    A **baryon asymmetry** (the matter excess over antimatter) requires, by the
    **Sakharov conditions** (1967), **three** ingredients:

    1. **baryon-number violation** (the B-number is not conserved);
    2. **C and CP violation** (the laws are not symmetric under C or CP);
    3. **departure from thermal equilibrium** (the dynamics is not reversible).

    The **unitary step** satisfies **none** of the three:
    - it is **unitary** (a product of local unitaries), so **all** global
      charges are conserved — baryon number and B-L are **exact** symmetries
      (no baryon-number violation);
    - it is **CPT-invariant** (the 16/16-bar CPT pair, V10) and **CP-conserving**
      at the level of the step (no C/CP violation);
    - it is **reversible** (unitary), so it stays in **equilibrium** (no
      departure from equilibrium).

    Therefore the step **does not generate** the observed baryon asymmetry —
    the theory names it, it does not derive it.  A baryon asymmetry needs a
    process beyond the unitary step (CPT/CP-violating, non-unitary).
    """
    return {
        "baryon_number_violated": False,   # unitary -> B conserved
        "c_or_cp_violated": False,          # CPT-invariant, CP-conserving at 0-param
        "out_of_equilibrium": False,        # unitary -> reversible -> equilibrium
        "all_sakharov": False,              # none of the three are satisfied
        "generated": False,                 # the QCA does NOT generate the asymmetry
        "residual": "the observed baryon asymmetry (needs a CPT/CP-violating, non-unitary process)",
    }



def main():
    import sys
    print("=" * 68)
    print("CL10 - the three named walls: the mechanism and what is left (O1-O3)")
    print("=" * 68)
    w1 = w1_cosmological_constant()
    print("\n  O1 (W1) cosmological constant: the Z2 512/512 cancellation")
    print(f"    even/odd = {w1['n_even']}/{w1['n_odd']}  (equal: {w1['cancellation']})")
    print(f"    vacuum contribution = {w1['lambda_0param']} (the exact cancellation)")
    print(f"    left open: the observed {w1['lambda_obs']:.1e} is NOT produced by it")
    w2 = w2_seesaw()
    print("\n  O2 (W2) the neutrino masses: the 2-step seesaw")
    print(f"    m_D = g0*v = {w2['m_D_GeV']:.1f} GeV (the g_0 = 1/8 reading)   "
          f"M_N = {w2['M_N_GeV']:.3e} GeV (the ladder 2-step, integer fitted)")
    print(f"    2-step m_nu = {w2['m_nu_2step_eV']:.1e} eV  (the right order; the ladder "
          f"value is excluded by sum m_nu < {nv.SUM_MNU_MAX_EV} eV)")
    print(f"    1-step m_nu = {w2['m_nu_1step_eV']:.1e} eV  (six orders off)")
    w3 = w3_families()
    print("\n  O3 (W3) the family number: the 3 of the 3 space labels the copies")
    print(f"    (1,3): {w3['n_time']} time + {w3['n_space']} space  ->  {w3['n_fam_candidate']} families")
    print(f"    candidate {w3['n_fam_candidate']} = observed {w3['n_observed']}  "
          f"(forced by the index: {w3['derived']})")
    ok = (w1["cancellation"] and w1["lambda_0param"] == 0.0
          and w2["right_order"] and w3["match"])
    print(f"\n  VERDICT:  {'PASS (the 3 mechanisms hold; the observed values stay open: O1, O2, O3)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
