#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 - the headline numbers, each with its provenance class.

(The module name is the pre-CL10 label of the theory; the theory is CL10.  The
name is kept because the constants below are imported across the repository.)

Every number this module produces carries the provenance class of Section 15 of
the publication: **algebraic** (fixed by the root), **convention** (a mapping or
scale identification), **fitted** (read from data), **conditional** (built from
measured numbers on a named input), **open** (not fixed).

Two readings that the publication has ABANDONED are deliberately not computed
here, so that no PASS can ever be printed for them:

* the proton lifetime at ``M_GUT = M_P`` (``tau_p ~ 4e47`` yr).  Section 13.4:
  that margin "rest[s] entirely on that abandoned convention"; the lifetime in
  force is the desert solution, ``tau_p = 1.8e36`` yr (V64, falsifier F5);
* the seesaw with ``m_D = v/2`` and ``M_R = M_P``.  The ladder of Section 12.3
  puts ``m_D = M_P (1/16pi)^10 = 118`` GeV and
  ``M_N = M_P (1/16pi)^3 = 9.6e13`` GeV, and the bare ladder value 0.146 eV is
  already excluded by ``sum m_nu < 0.12`` eV (falsifier F3).
"""

from . import config as _config
from . import coupling

# ---------------------------------------------------------------------------
# the scales and units (the same constants the experiments import).  Their
# values come from `cl10/config.yaml`: one edit there moves the whole ladder.
# ---------------------------------------------------------------------------
M_PLANCK = float(_config.CFG.physics.m_planck_gev)   # GeV (the Planck scale; M_GUT = M_P is NOT assumed - see proton_lifetime_desert)
M_P_ROBSON = float(_config.CFG.physics.proton_mass_gev)   # GeV (the proton mass)
ALPHA_10 = coupling.ALPHA_10      # 1/25 = g^2 at g = 1/5: the alpha_10 = g^2 READING (a convention, Section 1.4)
V_EW = float(_config.CFG.physics.v_electroweak_gev)  # GeV (the measured electroweak vev)
GEV_TO_S = float(_config.CFG.units.gev_to_s)   # 1/GeV in seconds
S_TO_YR = float(_config.CFG.units.s_to_yr)     # seconds per year
GEV2_TO_CM2 = float(_config.CFG.units.gev2_to_cm2)   # (hbar c)^2: 1 GeV^-2 in cm^2
C_CM_S = float(_config.CFG.units.c_cm_s)             # the speed of light, cm/s

LOOP = float(_config.CFG.physics.loop)     # the loop factor of the ladder (Section 12)
# the rung indices of the cell singlets (config block `cell`; cl10/cell.py):
# rung_B is the dark-matter/Dirac rung, rung_A the seesaw rung - one source
RUNG_A = int(_config.CFG.cell.rung_A)
RUNG_B = int(_config.CFG.cell.rung_B)
ALPHA_EM_INV = float(_config.CFG.sm.alpha_em_inv)   # the electromagnetic coupling at M_Z
EV_PER_GEV = float(_config.CFG.units.ev_per_gev)  # GeV -> eV (config units)

# the desert solution of V64 (two FITTED thresholds; the coupling is the root value)
M_X_DESERT = float(_config.CFG.desert.m_x_gev)   # GeV
ALPHA_X_DESERT = float(_config.CFG.desert.alpha_x)
# the same content at the closed angle (V133): alpha_X = pi/64, threshold ordering destroyed
M_X_CLOSED = float(_config.CFG.desert.m_x_closed_gev)   # GeV
ALPHA_X_CLOSED = coupling.alpha_std(coupling.G_CLOSED)


def proton_lifetime(M_X, alpha_X):
    """``tau_p ~ M_X^4 / (alpha_X^2 m_p^5)`` in years (dimensional estimate;
    the renormalisation factor, the hadronic matrix element and the X/Y mass
    splitting are NOT included - Section 13.4)."""
    tau_gev_inv = M_X ** 4 / (alpha_X ** 2 * M_P_ROBSON ** 5)
    return tau_gev_inv * GEV_TO_S / S_TO_YR


def proton_lifetime_desert():
    """The lifetime in force: the desert solution, ``M_X = 1.8e16`` GeV at
    ``alpha_X = 1/25`` (V64).  Provenance **conditional**: the spectrum carries
    two fitted thresholds, and the identification ``alpha_X = 1/25`` is the
    ``alpha_10 = g^2`` convention."""
    return proton_lifetime(M_X_DESERT, ALPHA_X_DESERT)


def proton_lifetime_closed():
    """The same content at the closed angle ``alpha_X = pi/64 = 0.0491`` (V133):
    the meeting moves to ``M_X = 1.235e17`` GeV and the lifetime to
    ``2.8e39`` yr, with the ``(1,2)`` threshold above ``M_X`` - the desert
    threshold ordering is destroyed."""
    return proton_lifetime(M_X_CLOSED, ALPHA_X_CLOSED)


def ladder_scale(n):
    """One rung of the ladder: ``M_P (1/16pi)^n`` in GeV (Section 12.4).
    Provenance of the integer ``n``: **fitted** (read from the measured scale)."""
    return M_PLANCK * LOOP ** n


def dm_ladder_mass():
    """The dark-matter mass SCALE: ``M_P (1/16pi)^10 = 118`` GeV (V54, F7).
    Provenance **fitted** (the coefficient 1.18 against the canonical 100 GeV).
    The carrier: the cell's singlet level B occupies this rung n = 10
    (Section 12.2, cl10/cell.py, V195/V196); its portal strength to the Higgs
    is an input (Section 15, O7)."""
    return ladder_scale(RUNG_B)


def neutrino_ladder():
    """The neutrino sector of the ladder (Section 12.3, V57/V103/V114).

    ``m_D = M_P (1/16pi)^10``, ``M_N = M_P (1/16pi)^3``,
    ``m_nu = m_D^2 / M_N = M_P (1/16pi)^17``; the exponent identity
    ``2*10 - 3 = 17`` is exact, the choice ``n_2step = 3`` is **fitted** (and is
    not a link-word length: the ``Delta L = 2`` operator lives in Clifford grade
    5 while every adjoint link word is even-grade, V136)."""
    m_D = ladder_scale(RUNG_B)
    M_N = ladder_scale(RUNG_A)
    m_nu_gev = m_D ** 2 / M_N
    m_nu_eV = m_nu_gev * EV_PER_GEV
    return {
        "m_D_GeV": m_D,
        "M_N_GeV": M_N,
        "m_nu_eV_ladder": m_nu_eV,
        "m_nu_eV_third": m_nu_eV / 3.0,          # the 1/3 candidate (sin^2 theta_12 = 1/3)
        "sum_m_nu_candidate_eV": m_nu_eV / 3.0,  # normal ordering, two massless states
        "status_ladder": "EXCLUDED by sum m_nu < 0.12 eV (F3)",
        "status_third": "candidate, not a result: the junction amplitude is free inside a finite window (O2)",
        "provenance_exponent": "fitted (n_2step = 3); the structure 2*10-3 = 17 is exact",
    }


def hierarchy(F10=None):
    """The electroweak/Planck hierarchy of Section 12.1 (V85/V89/V108/V113, V206).

    The closed traced loops at the rung-10 clock ``g = 1/5`` measure
    ``F_10 = 0.019186`` (against ``1/16pi = 0.0199``), so

        v / M_P = 3 F_10^10 = 2.0881 x (1/16pi)^10 ,

    99.44% of the observed coefficient 2.0765.  At the closed clock the
    crossing moves to ``n* = 208`` and the coefficient there reads 2.702
    (+30%).  V206 measures the clock tension: the closure class reads
    {2.2995, 2.7018, 2.7048} while the observed value is read by the
    non-closed clock ``1/5`` to relative 1e-5 - one tension with the O4
    meeting coupling, not two.  The factor 3 is the one-split-node-per-
    axis-family choice, so the provenance is **conditional** on that node
    set (O8).  ``F10`` defaults to ``physics.f_10_measured`` of
    ``config.yaml``."""
    F10 = float(_config.CFG.physics.f_10_measured) if F10 is None \
        else float(F10)
    ratio = 3.0 * F10 ** RUNG_B            # v / M_P
    coefficient = ratio / LOOP ** RUNG_B   # the quoted two-digit number
    return {"F10": F10,
            "v_over_MPlanck": ratio,
            "coefficient": coefficient,
            "observed_v_over_MPlanck":
                float(_config.CFG.physics.v_over_mplanck_observed),
            "observed_coefficient":
                float(_config.CFG.physics.hierarchy_coefficient_observed)}


def family_number():
    """The family number: ``N_gen = 3 = n_space`` (Section 4, V26/V91/V116/V127).
    Provenance **open**: the index labels the three copies but does not force
    the count - the cubic group is transitive on the 12 doubler nodes, the only
    cubic-symmetric split is the full set of 12 (whose meeting scale is
    excluded), and the 3-node triangle has a cubic orbit of 8 (O3)."""
    return {"value": 3, "identification": "3 = n_space = C(3,2) = rank of the spatial metric",
            "provenance": "open (a symmetry-breaking input, O3)"}


# ---------------------------------------------------------------------------
# the coherence check: the numbers this module prints are the numbers the
# publication quotes.  This is NOT a test against nature - the provenance
# class says what each number is.
# ---------------------------------------------------------------------------
def verify(tol=0.02):
    """Recompute every headline number and compare it with the value quoted in
    the publication.  ``ok`` means code-and-paper agree, nothing more.

    The two proton rows carry a wider tolerance: the publication quotes V64,
    which fits two thresholds, while this module evaluates the same formula at
    the convention inputs (``alpha_GUT = 1/40``, ``M_X = 1.8e16`` GeV)."""
    tau_desert = proton_lifetime_desert()
    tau_closed = proton_lifetime_closed()
    nu = neutrino_ladder()
    hier = hierarchy()
    rows = {
        "tau_p_desert_yr": (tau_desert, 1.64e36, 0.2),
        "tau_p_closed_yr": (tau_closed, 2.68e39, 0.2),
        "m_DM_GeV": (dm_ladder_mass(), 118.5, tol),
        "m_nu_ladder_eV": (nu["m_nu_eV_ladder"], 0.146, tol),
        "m_nu_third_eV": (nu["m_nu_eV_third"], 0.0487, tol),
        "M_N_GeV": (nu["M_N_GeV"], 9.6e13, tol),
        "hierarchy_coefficient": (hier["coefficient"], 2.0881, tol),
        "v_over_MPlanck": (hier["v_over_MPlanck"], 2.0167e-17, tol),
    }
    dev, ok = {}, True
    for k, (computed, quoted, t) in rows.items():
        dev[k] = abs(computed - quoted) / quoted
        ok = ok and dev[k] < t
    return {"values": {k: v[0] for k, v in rows.items()},
            "quoted": {k: v[1] for k, v in rows.items()},
            "relative_deviation": dev,
            "ok": ok}


def main():
    import sys
    r = verify()
    print("=" * 68)
    print("CL10 - the headline numbers, with their provenance (Section 15)")
    print("=" * 68)
    print(f"  angle branch: {coupling.branch()}  g = {coupling.G10:.6g}   "
          f"[{coupling.provenance()}]")
    print(f"  proton  tau_p (desert, V64):        {r['values']['tau_p_desert_yr']:.2e} yr   "
          f"[conditional: 2 fitted thresholds]")
    print(f"  proton  tau_p (closed angle, V133): {r['values']['tau_p_closed_yr']:.2e} yr   "
          f"[conditional: threshold ordering destroyed]")
    print(f"  dark matter  M_P(1/16pi)^10:        {r['values']['m_DM_GeV']:.1f} GeV   "
          f"[fitted scale; the carrier is OPEN (O7)]")
    print(f"  seesaw    M_N:                      {r['values']['M_N_GeV']:.2e} GeV   [fitted n = 3]")
    print(f"  neutrino  ladder:                   {r['values']['m_nu_ladder_eV']:.3f} eV   "
          f"[EXCLUDED by sum m_nu < 0.12 eV (F3)]")
    print(f"  neutrino  1/3 candidate:            {r['values']['m_nu_third_eV']:.4f} eV   "
          f"[candidate, structure protected (O2)]")
    print(f"  hierarchy 3 F_10^10:                {r['values']['hierarchy_coefficient']:.4f} x (1/16pi)^10 "
          f"= {r['values']['v_over_MPlanck']:.3e}   "
          f"[conditional on the one-node-per-family set (O8)]")
    print(f"  families  N_gen = 3:                {family_number()['provenance']}")
    print(f"  code-and-paper agreement:           "
          f"max relative deviation {max(r['relative_deviation'].values()):.2e}")
    print(f"  VERDICT:  {'PASS (the module reproduces the published numbers)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
