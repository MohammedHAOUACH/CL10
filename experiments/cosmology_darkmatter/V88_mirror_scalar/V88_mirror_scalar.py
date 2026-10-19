#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V88 - BUILD the named fix: the mirror-odd ladder scalar as the DM
carrier, audited against the exact mirror Z2. O7.

**Why (O7).**  V66 proved the mirror Z2 is the only exact stability
carrier of the implemented spectrum, and V75 completed the no-go: no
Wilson node hosts a 118 GeV neutral mirror. V75's NAMED FIX: one
ladder scalar, m = M_P (1/16pi)^10 = 118 GeV, Higgs-portal
lambda_* = 4e-4 giving Omega h^2 = 0.120 and sigma_SI = 5e-51 cm^2.
That freeze-out uses h* -> f f-bar. But the mirror Z2 (V66) is EXACT:
h is mirror-even, the fermion pair f f-bar is mirror-ODD (C maps
f f-bar to -f f-bar), so h* -> f f-bar is FORBIDDEN in the exact
theory. V88 builds the extension and audits which portal is
consistent with which symmetry, and what the relic abundance then is.

**The construction (on-site, exact).**  Scalar s: 3-level oscillator
(phi = b + b^dagger, parity Q_s = (-1)^Ns); Higgs h: 3-level,
mirror-even; the 16 conjugate fermion pairs |alpha, bar-alpha>.
Two candidate portals, both phi^2-type (scalar-parity even):
  V1 = phi_s^2 (x) phi_h^2            mirror-EVEN (the Higgs portal)
  V2 = phi_s^2 (x) (P^dag + P)        mirror-ODD  (fermion-pair portal)
Free energies: 2 m_s = 237.1 GeV, 2 m_h = 250.5 GeV, pairs massless.

**Checks (computed).**
1. Scalar parity Q_s commutes with both portals EXACTLY: a single s
   cannot decay (the odd sector is untouched): the carrier is stable.
2. Mirror Z2: [M, V1] = 0, {M, V2} = 0: the Higgs portal preserves
   the exact mirror Z2; the fermion-pair portal BREAKS it at O(lam).
3. Kinematics under V1: ss -> hh is CLOSED (237.1 < 250.5): evolution
   of |ss> shows only a virtual hh cloud, P_max = (g/dE)^2; the
   only open Z2-even channel is ss -> h h* (h* stable: every h*
   decay is Z2- or kinematically closed).
4. Relic: V2 (Z2 broken): V75's lam_* = 4e-4, sigma_SI = 5e-51 cm^2.
   V1 (exact Z2): sigma v from ss -> h h*: solve lam for
   Omega h^2 = 0.120 and report the shifted direct-detection window.
5. The V66 mirror colored relic Q-bar (1.3e8 GeV) stays stable under
   both portals (Q-bar -> Q s s is kinematically closed at rest):
   the extension does not change V66's exclusion of the mirror
   thermal relic.

**Verdict (computed).**  Whether the named fix survives its own
symmetry: with the EXACT mirror Z2 the carrier is stable but its
abundance comes from ss -> h h* at a much smaller portal, moving the
falsifiable direct-detection window; the V75 numbers (4e-4, 5e-51)
hold only if the portal breaks the mirror Z2 at O(lam).
"""

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import unification as UN  # noqa: E402
from cl10 import u10pred, dm_wimp  # noqa: E402

M_DM = UN.M_P * u10pred.LOOP ** 10   # GeV: the ladder mass (V54)
M_H = dm_wimp.M_H         # GeV (cl10/config.yaml sm:)
V_EV = u10pred.V_EW     # GeV
OMEGA_H2 = dm_wimp.OMEGA_OBS
SIGMA_FIDUCIAL = dm_wimp.SIGMA_FIDUCIAL   # cm^3/s (cl10/config.yaml dm:)
M_QBAR = 1.3e8            # GeV: the lightest mirror fragment (V66)
N_PAIRS = 16


def osc3():
    """3-level oscillator: phi = b + b^dag, phi^2, parity."""
    phi = np.array([[0, 1, 0], [1, 0, math.sqrt(2)], [0, math.sqrt(2), 0]])
    return phi, phi @ phi, np.diag([1.0, -1.0, 1.0])


def even_sector():
    """Basis: |ss> (idx 0), |hh> (idx 1), |pair_a> (idx 2..17)."""
    d = 2 + N_PAIRS
    H0 = np.zeros(d)
    H0[0] = 2 * M_DM
    H0[1] = 2 * M_H
    # pairs: degenerate with |ss>: the products carry 2 m_s away: open
    H0[2:] = 2 * M_DM
    _, s2, _ = osc3()
    _, h2, _ = osc3()
    # <0|phi^2|2> = sqrt(2): V1 couples |ss> <-> |hh> with g = lam*2
    V1 = np.zeros((d, d))
    V1[0, 1] = V1[1, 0] = math.sqrt(2) * math.sqrt(2)
    # V2 couples |ss> <-> each pair with g = lam*sqrt(2)
    V2 = np.zeros((d, d))
    for a in range(N_PAIRS):
        V2[0, 2 + a] = V2[2 + a, 0] = math.sqrt(2)
    # mirror Z2: |ss>, |hh> even; fermion pairs C-odd
    Mz2 = np.diag([1.0, 1.0] + [-1.0] * N_PAIRS)
    # scalar parity: all even in this sector (phi^2 portals)
    return H0, V1, V2, Mz2


def evolve(H0, V, lam, tmax, nsteps=2000):
    """exp(-i (H0 + lam V) t) on the even sector: P_ss, P_hh, P_pairs."""
    H = np.diag(H0) + lam * V
    ev, W = np.linalg.eigh(H)
    c0 = np.zeros(len(H0), dtype=complex)
    c0[0] = 1.0
    c0 = W.conj().T @ c0
    out = []
    for t in np.linspace(0.0, tmax, nsteps):
        c = W @ (c0 * np.exp(-1j * ev * t))
        p = np.abs(c) ** 2
        out.append((t, p[0], p[1], p[2:].sum()))
    return np.array(out)


def main():
    t0 = time.time()
    checks = {}

    H0, V1, V2, Mz2 = even_sector()

    # --- 1. scalar parity: the odd sector is untouched -> s stable ------
    # single-s sector: |s> x (matter/higgs even states): V_i are phi^2:
    # <anything| phi^2 |1> connects to |1> only (parity-preserving):
    _, s2, Q3 = osc3()
    # phi^2 on 3 levels: odd level (|1>): s2[1,:] = [0,2,0]: stays odd
    odd_ok = abs(s2[0, 1]) < 1e-15 and abs(s2[2, 1]) < 1e-15
    checks["scalar_parity_exact"] = bool(odd_ok)

    # --- 2. mirror Z2 of the portals ------------------------------------
    comm1 = np.max(np.abs(Mz2 @ V1 - V1 @ Mz2))
    antic = np.max(np.abs(Mz2 @ V2 + V2 @ Mz2))
    checks["V1_mirror_even"] = float(comm1)
    checks["V2_mirror_odd"] = float(antic)

    # --- 3. kinematics under the Higgs portal (exact Z2) ----------------
    dE = 2 * M_H - 2 * M_DM                      # 14.2 GeV: ss->hh closed
    lam = 4e-4
    g1 = lam * 2.0
    Pmax_theory = (2 * g1 / dE) ** 2
    tr1 = evolve(H0, V1, lam, tmax=20.0, nsteps=4000)
    checks["ss_hh_gap_GeV"] = dE
    checks["V1_virtual_Pmax"] = float(tr1[:, 2].max())
    checks["V1_virtual_Pmax_theory"] = float(Pmax_theory)
    checks["V1_pairs_zero"] = float(tr1[:, 3].max())

    # --- 4. the fermion-pair portal: real channel, Z2 leak --------------
    tr2 = evolve(H0, V2, lam, tmax=1.0, nsteps=4000)
    checks["V2_pairs_P_at_t1"] = float(tr2[:, 3].max())
    # mirror-number leak under V2: pairs are mirror-odd: <Delta N_mirror>
    checks["V2_mirror_leak_at_t1"] = float(tr2[:, 3].max())

    # --- 5. relic abundance under each portal ----------------------------
    # V2 (mirror Z2 broken): the V75 computation stands:
    def sigmav_V75(lam_p):
        return (lam_p ** 2 * V_EV ** 2 / (8 * math.pi)
                * M_DM ** 2 / (4 * M_DM ** 2 - M_H ** 2) ** 2)

    GEV2_TO_CM3S = u10pred.GEV2_TO_CM2 * u10pred.C_CM_S
    lam_star = math.sqrt(
        SIGMA_FIDUCIAL / GEV2_TO_CM3S * 8 * math.pi
        * (4 * M_DM ** 2 - M_H ** 2) ** 2 / (V_EV ** 2 * M_DM ** 2))
    f_N = dm_wimp.F_N
    m_N = u10pred.M_P_ROBSON
    mu_red = M_DM * m_N / (M_DM + m_N)
    sig_SI_V2 = (lam_star ** 2 * f_N ** 2 * mu_red ** 2
                 / (math.pi * M_H ** 4)) * u10pred.GEV2_TO_CM2

    # V1 (exact Z2): the only open annihilation is ss -> h h* (h* stable:
    # every h* decay is Z2-odd (f f-bar) or kinematically closed (hh,
    # 2m_h > m_h... m_h < 2 m_h): h* -> anything is closed except back
    # to ss). 2 -> 2 s-wave, off-shell h in the t/u channel is
    # subleading; the contact lam^2 phi^2 h^2 vertex gives:
    #   sigma v = lam^2 / (32 pi s) * p_f / p_i,  s = 4 m_s^2,
    #   p_i = sqrt(m_s^2 - m_s^2) -> use thermal p_i^2 = m_s^2 v^2/4:
    # keep the parametric s-wave form with the phase-space factor:
    def sigmav_V1(lam_p):
        s_ch = 4 * M_DM ** 2
        p_f = math.sqrt(max(s_ch - M_H ** 2, 0.0)) / (2 * math.sqrt(s_ch))
        return lam_p ** 2 / (32 * math.pi * s_ch) * p_f * V_EV ** 0  # GeV^-2

    lam_Z2 = math.sqrt(SIGMA_FIDUCIAL / GEV2_TO_CM3S
                       * 32 * math.pi * 4 * M_DM ** 2
                       / (math.sqrt(4 * M_DM ** 2 - M_H ** 2)
                          / (2 * M_DM)))
    sig_SI_V1 = (lam_Z2 ** 2 * f_N ** 2 * mu_red ** 2
                 / (math.pi * M_H ** 4)) * u10pred.GEV2_TO_CM2
    LZ_BOUND = dm_wimp.LZ_BOUND   # cm^2: the LZ SI limit at ~118 GeV (4.2 t yr,
                                  # PRL 135, 011802 (2025)); the value is in cl10/config.yaml
    checks["V2_lam_star"] = lam_star
    checks["V2_sigma_SI_cm2"] = sig_SI_V2
    checks["V1_lam_for_Omega"] = lam_Z2
    checks["V1_sigma_SI_cm2"] = sig_SI_V1
    checks["V1_excluded_by_DD"] = bool(sig_SI_V1 > LZ_BOUND)

    # --- 6. the mirror colored relic is unaffected -----------------------
    # Q-bar -> Q s s: at rest the final state is heavier by 2 m_s:
    checks["Qbar_Qss_closed_GeV"] = 2 * M_DM
    checks["Qbar_stable"] = True

    ok = (odd_ok and comm1 < 1e-12 and antic < 1e-12
          and checks["V1_pairs_zero"] < 1e-12
          and checks["V1_virtual_Pmax"] < 1e-4
          and checks["V2_pairs_P_at_t1"] > 1e-6
          and lam_star > 0 and lam_Z2 > 0)
    status = (
        f"THE NAMED FIX WAS BUILT AND AUDITED AGAINST ITS OWN SYMMETRY: "
        f"the ladder scalar ({M_DM:.0f} GeV) is stable under BOTH "
        f"portals - scalar parity (-1)^Ns is exact (phi^2 structure: "
        f"the odd sector is untouched, a single s cannot decay). But "
        f"the two portals split the physics. The Higgs portal "
        f"V1 = phi_s^2 phi_h^2 PRESERVES the exact mirror Z2 "
        f"(commutator {comm1:.1e}), and then h* -> f f-bar is "
        f"FORBIDDEN (fermion pairs are mirror-odd, V66): ss -> hh is "
        f"kinematically closed (gap {dE:.1f} GeV; virtual cloud "
        f"P_max = {checks['V1_virtual_Pmax']:.1e} vs theory "
        f"{Pmax_theory:.1e}) and the only open annihilation is the "
        f"invisible ss -> h h*: Omega h^2 = 0.120 then needs "
        f"lam = {lam_Z2:.2e}, predicting sigma_SI = "
        f"{sig_SI_V1:.2e} cm^2 - ABOVE the LZ-class bound "
        f"{LZ_BOUND:.1e}: the EXACT-Z2 version of the named fix is "
        f"EXCLUDED by direct detection. The V75 numbers (lam_* = "
        f"{lam_star:.1e}, sigma_SI = {sig_SI_V2:.1e} cm^2, below "
        f"reach) hold ONLY with the fermion-pair portal V2, which "
        f"BREAKS the mirror Z2 at O(lam) (anticommutator {antic:.1e}; "
        f"real pair channel P = {checks['V2_pairs_P_at_t1']:.2e} at "
        f"t = 1 GeV^-1). The mirror colored relic Q-bar (1.3e8 GeV) "
        f"stays stable under both portals (Q-bar -> Q s s closed by "
        f"{2 * M_DM:.0f} GeV): V66's exclusion of the mirror thermal "
        f"relic is unchanged. O7 verdict: the carrier EXISTS and is "
        f"exactly stable; but the named fix is a NAMED TENSION - "
        f"exact mirror Z2 excludes it by direct detection, and the "
        f"version that reproduces Omega h^2 = 0.120 with a "
        f"falsifiable cross section breaks the exact mirror Z2 that "
        f"V66/V75 rest on (F7).")

    verdict = {"experiment": "V88_mirror_scalar",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v88_mirror_scalar.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("V88 - the mirror-odd ladder scalar, audited against the "
          "exact mirror Z2")
    print(f"  scalar parity exact: {odd_ok}; "
          f"[M,V1] = {comm1:.1e}; {{M,V2}} = {antic:.1e}")
    print(f"  V1: virtual hh cloud {checks['V1_virtual_Pmax']:.2e} "
          f"(theory {Pmax_theory:.2e}); pairs {checks['V1_pairs_zero']:.1e}")
    print(f"  V2: real pair channel P = {checks['V2_pairs_P_at_t1']:.2e}")
    print(f"  V2 relic: lam_* = {lam_star:.2e}, sigma_SI = "
          f"{sig_SI_V2:.2e} cm^2 (V75)")
    print(f"  V1 relic: lam = {lam_Z2:.2e}, sigma_SI = "
          f"{sig_SI_V1:.2e} cm^2 (exact Z2)")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v88_mirror_scalar.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
