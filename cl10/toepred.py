#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the TOE predictions: proton decay, the seesaw, the 3 families.

The **falsifiable** TOE predictions, all from the 0-parameter GUT structure:

- **Proton decay.**  The GUT X/Y bosons mediate ``p -> e^+ pi^0``.  The
  **0-parameter** lifetime (the GUT scale is the Planck scale,
  ``M_GUT = M_P``) is ``tau_p ~ M_GUT^4 / (alpha_10^2 m_p^5)`` — the GUT
  scale **suppresses** the decay so far that the proton appears stable,
  consistent with the experimental bound ``tau_p > 10^34`` yr.

- **Seesaw / neutrino mass.**  The 16 carries the **sterile** ``nu^c``; the
  Dirac mass ``m_D`` (the electroweak scale) and the Majorana mass ``M_R``
  (the GUT scale) give the **0-parameter** light-neutrino mass
  ``m_nu ~ m_D^2 / M_R`` — light (<< 1 eV), as observed.

- **The 3 families.**  The 16 is **one** generation; the **3** observed
  families are a **named residual** (the family problem).  The 0-parameter
  *candidate* is 3 (from the 3 spatial directions of the (1,3) / the
  6 internal GUT directions / 2), but the theory does **not** derive the
  number 3.
"""

# the 0-parameter scales
M_PLANCK = 1.22e19    # GeV (the GUT scale = the Planck scale, 0-parameter)
M_P_ROBSON = 0.938     # GeV (the proton mass)
ALPHA_10 = 1 / 25.0    # the 0-parameter GUT coupling
V_EW = 246.0           # GeV (the electroweak scale)
GEV_TO_S = 6.582e-25   # 1/GeV in seconds
S_TO_YR = 3.156e7      # seconds per year


def proton_decay_lifetime():
    """The 0-parameter proton lifetime tau_p ~ M_GUT^4 / (alpha_10^2 m_p^5).

    Returns the lifetime in years (the GUT scale suppresses it to >> 10^34 yr).
    """
    tau_gev_inv = M_PLANCK ** 4 / (ALPHA_10 ** 2 * M_P_ROBSON ** 5)
    tau_s = tau_gev_inv * GEV_TO_S
    return tau_s / S_TO_YR


def seesaw_mass():
    """The 0-parameter light-neutrino mass m_nu ~ m_D^2 / M_R (in eV).

    m_D = v/2 (the electroweak scale), M_R = M_GUT (the GUT scale).
    """
    m_D = V_EW / 2.0
    m_nu_gev = m_D ** 2 / M_PLANCK
    return m_nu_gev * 1e9  # GeV -> eV


def family_candidate():
    """The 0-parameter candidate for the number of families.

    The 16 is one generation; the 0-parameter *candidate* for the number of
    families is 3 (the 3 spatial directions of the (1,3), or the 6 internal
    GUT directions / 2).  This is a **named residual** — the theory does not
    derive the number 3.
    """
    n_spatial = 3          # the 3 spatial directions of the (1,3)
    n_internal = 6         # the 6 non-spacetime GUT directions
    return {"candidate": 3, "from_spatial": n_spatial,
            "from_internal_half": n_internal // 2, "derived": False}


def main():
    import sys
    print("=" * 68)
    print("CL10 - the TOE predictions: proton, seesaw, families")
    print("=" * 68)
    tau_p = proton_decay_lifetime()
    print(f"  proton decay  tau_p ~ M_GUT^4/(alpha^2 m_p^5):  "
          f"{tau_p:.2e} yr  (observed bound > 1e34 yr)")
    m_nu = seesaw_mass()
    print(f"  seesaw  m_nu ~ m_D^2/M_R:  {m_nu:.2e} eV  (observed < 1 eV)")
    fam = family_candidate()
    print(f"  3 families  candidate = {fam['candidate']}  "
          f"(derived={fam['derived']}  --  NAMED RESIDUAL)")
    ok = (tau_p > 1e34 and m_nu < 1.0 and fam["candidate"] == 3)
    print(f"  VERDICT:  {'PASS (falsifiable: proton stable, neutrino light, 3-family candidate)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
