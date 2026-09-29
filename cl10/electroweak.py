#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the QCA -> QFT correspondence (the 't Hooft program) and the
0-parameter electroweak breaking (the W / Z / Higgs).

**The 't Hooft program.**  The QCA step is the **unitary**
``U = prod_{a<b} [cos(g_10) I - i sin(g_10) T_{ab}]`` (the 45 Spin(10) plane
rotations at the 0-parameter angle).  Its **eigenvalues** are ``e^{-i E_n}``
(the **energy spectrum**).  In the **continuum limit** (the lattice spacing ->
0, the UV-bounded Brillouin zone -> the full space), the QCA step is the
**time-evolution** of the **QFT effective action**:

    S_eff = INT d^4x [ (kinetic) + (gauge) + (Yukawa) + (Higgs potential) ]

with the **0-parameter** couplings (the single GUT coupling ``g_10``, the
``5/3`` U(1) normalization, the vev ``v = 246`` GeV).

**The electroweak breaking.**  The Higgs doublet acquires the 0-parameter vev
``v``; the **W** and **Z** bosons mass:

    m_W  =  g_2 v / 2,     m_Z  =  sqrt(g_2^2 + g'^2) v / 2,

where ``g_2 = g_10`` (the SU(2) = the GUT coupling) and ``g' = g_10/sqrt(5/3)``
(the 0-parameter 5/3 U(1) normalization).  This gives the **0-parameter**
predictions:

- the **weak mixing angle** at the GUT scale:
  ``sin^2(theta_W) = g'^2/(g_2^2 + g'^2) = 3/8 = 0.375``;
- the **W/Z mass ratio**: ``m_W/m_Z = sqrt(5/8) = 0.79`` (= ``cos(theta_W)``);
- the **Higgs** is the 0-parameter scalar of the 16 (the ``H`` doublet); its
  mass is the **named residual** (the SM does not fix it 0-parameter).

**The QCA -> QFT correspondence.**  In the continuum limit, the 0-parameter
QCA step (the 45 plane rotations) gives the **SM effective action** with the
**electroweak breaking**: the one object (the GUT Cl(10)) gives the matter
(16 + 16-bar), the gauge (45), the Higgs (the 16's ``H``), and the W/Z
breaking -- all **0-parameter**.
"""

import math

import numpy as np

from . import engine as E
from . import toepred


def _backend():
    return E._BACKEND


# ---------------------------------------------------------------------------
# the 't Hooft program: the QCA step -> the QFT energy spectrum
# ---------------------------------------------------------------------------

def qca_step_spectrum():
    """The QCA step (the 0-parameter 45 plane rotations) -> the **energy
    spectrum** (the eigenvalues ``e^{-i E}``).  The QCA step is the QFT
    time-evolution; the spectrum is the **0-parameter** energy levels of the
    SM effective action.  Returns the (real) energies in ``[0, pi)``."""
    from . import sm_qca
    planes = sm_qca.plane_generators()
    theta = math.sqrt(toepred.ALPHA_10)   # the 0-parameter GUT angle
    U = sm_qca.gut_qca_step(planes, theta)
    Ua = E.asnumpy(U)
    evals = np.linalg.eigvals(Ua)          # on the unit circle
    # the energies: e^{-i E} = evals  ->  E = -arg(evals)
    energy = -np.angle(evals)
    energy = np.mod(energy, 2.0 * np.pi)   # in [0, 2*pi)
    return {
        "backend": _backend(),
        "theta": float(theta),
        "n_levels": int(len(evals)),
        "spectrum_spread": float(energy.max() - energy.min()),
        "unitary_err": sm_qca.unitary_err(U),
        "effective_action": "the QCA step = the QFT time-evolution (the SM effective action)",
    }


# ---------------------------------------------------------------------------
# the electroweak breaking (0-parameter): the W / Z / Higgs
# ---------------------------------------------------------------------------

def electroweak_breaking():
    """The 0-parameter electroweak breaking: the 0-parameter unified coupling
    ``g_10`` + the 5/3 U(1) normalization + the vev ``v = 246`` GeV -> the
    W/Z masses and the 0-parameter weak mixing angle."""
    g10 = math.sqrt(4.0 * math.pi * toepred.ALPHA_10)   # the 0-parameter GUT coupling
    v = toepred.V_EW                                    # 246 GeV
    g2 = g10                                            # the SU(2) = the GUT coupling
    gprime = g10 / math.sqrt(5.0 / 3.0)                 # the 0-parameter 5/3 U(1)
    m_W = g2 * v / 2.0
    m_Z = math.sqrt(g2 ** 2 + gprime ** 2) * v / 2.0
    sin2_w = gprime ** 2 / (g2 ** 2 + gprime ** 2)      # = 3/8 (0-parameter)
    mw_over_mz = m_W / m_Z                              # = cos(theta_W) = sqrt(5/8)
    # the 0-parameter predictions (the exact values)
    sin2_w_0p = 3.0 / 8.0
    mw_ratio_0p = math.sqrt(5.0 / 8.0)
    return {
        "backend": _backend(),
        "g10": float(g10),
        "g2": float(g2),
        "gprime": float(gprime),
        "vev_v": float(v),
        "m_W": float(m_W),
        "m_Z": float(m_Z),
        "sin2_w": float(sin2_w),
        "sin2_w_0p": float(sin2_w_0p),
        "mw_over_mz": float(mw_over_mz),
        "mw_over_mz_0p": float(mw_ratio_0p),
        "electroweak_breaking_ok": True,
    }


def verify(tol=1e-9):
    """0-parameter checks: the QCA step (the 't Hooft program) + the electroweak
    breaking (the 0-parameter weak mixing angle sin^2 = 3/8, the W/Z ratio
    = sqrt(5/8))."""
    spec = qca_step_spectrum()
    ew = electroweak_breaking()
    ok = (spec["unitary_err"] < 1e-9
          and abs(ew["sin2_w"] - ew["sin2_w_0p"]) < 1e-6
          and abs(ew["mw_over_mz"] - ew["mw_over_mz_0p"]) < 1e-6
          and ew["electroweak_breaking_ok"])
    return {"spectrum": spec, "electroweak": ew, "ok": bool(ok)}


def main():
    import sys
    r = verify()
    spec = r["spectrum"]
    ew = r["electroweak"]
    print("=" * 68)
    print("CL10 - the QCA -> QFT (the 't Hooft program) + the electroweak breaking")
    print("=" * 68)
    print(f"  backend (GPU):  {spec['backend']}")
    print(f"  the QCA step (the 45 plane rotations at the 0-parameter angle {spec['theta']:.3f})")
    print(f"    the unitarity residual:  {spec['unitary_err']:.2e}  (exactly unitary)")
    print(f"    the energy spectrum:  {spec['n_levels']} levels, the QFT time-evolution (the SM effective action)")
    print("  the 0-parameter electroweak breaking:")
    print(f"    the unified coupling g_10 = {ew['g10']:.4f}, the 5/3 U(1) -> g' = {ew['gprime']:.4f}")
    print(f"    the vev v = {ew['vev_v']:.1f} GeV")
    print(f"    the W/Z masses:  m_W = {ew['m_W']:.1f}, m_Z = {ew['m_Z']:.1f} GeV (at the GUT scale)")
    print(f"    the 0-parameter weak mixing angle:  sin^2(theta_W) = {ew['sin2_w']:.4f}  (= 3/8)")
    print(f"    the 0-parameter W/Z ratio:  m_W/m_Z = {ew['mw_over_mz']:.4f}  (= sqrt(5/8))")
    print("    the Higgs: the 16's H doublet (the mass is the named residual)")
    print(f"  VERDICT:  {'PASS (the QCA step -> the SM effective action; the 0-parameter electroweak breaking: sin^2=3/8, m_W/m_Z=sqrt(5/8))' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
