#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the sub-Planck dark-matter scale: the 0-parameter 10-loop WIMP scale
(V54, O7).

The O7 residual (Section 13): the **sub-Planck** dark-matter scale that gives
the **observed** abundance.  V42 showed the 0-parameter **Planck-mass soliton**
DM ($m \sim M_P$, the stable soliton of the non-linear metric) is
**under**-abundant ($\Omega_{\rm 0param}/\Omega_{\rm obs} \sim 10^{-52}$) — the
observed DM is *not* the Planck-mass soliton.

**The sub-Planck scale this module computes (O7).**  The scale is the **10-loop
suppression of the Planck scale**:

    m_DM  =  M_P \times (1/16\pi)^{10}

where $(1/16\pi)^{10}$ is the **same 0-parameter 10-loop factor** as V53
(the O8 hierarchy): the **10** = the 0-parameter number of
$\mathrm{Cl}(10)$ generators ($n_{\rm gen}$), the **1-loop factor**
$1/(16\pi) \approx 0.0199$.  So

    m_DM  =  1.22e19 \times 9.7e-18  \approx  118 GeV

— the **WIMP / weak scale**.

**What the WIMP miracle adds (an input, not a derivation).**  A thermally
decoupled relic of mass $\sim 100$ GeV with a weak-scale annihilation
cross-section gives $\Omega h^2 \approx 0.12$: that statement is **standard
cosmology**, not a result of the automaton.  The ladder scale sits within a
factor $1.18$ of that reference scale.  This module does **not** conclude that
the observed DM mass is the weak scale: **no state of the implemented spectrum
carries the scale** (the paper's O7/F6), so the coincidence is a scale
statement, not the prediction of a particle.

**The unification with V53.**  The **same** 10-loop factor $(1/16\pi)^{10}$
$= (n_{\rm gen}$ loop suppressions) gives **both** the hierarchy exponent -17
(V53: $v/E_P \sim (1/16\pi)^{10}$) **and** the sub-Planck DM mass
(V54: $m_{\rm DM} = M_P \times (1/16\pi)^{10}$).  One 0-parameter number
(unified): the 10-loop suppression of the single Planck scale.

**What is derived, what is input, what is open.**  Derived: the ladder
arithmetic $M_P (1/16\pi)^{10} = 118$ GeV and its contrast with the V42
Planck-mass soliton.  Input: the observed $\Omega h^2 = 0.12$, the $100$ GeV
reference, and the relic convention ``relic_abundance`` (normalised to the
observed value at the reference scale **by construction**, so it is a contrast
device, not an abundance prediction).  Open (**O7**): which state carries the
scale, and the cross-section (the RG running, O9).  The **V42 gap is not
resolved** here: it is relocated, from "the soliton is under-abundant" to "the
right scale has no carrier".
"""

from . import config as _config
from . import engine as E
from . import u10pred
from . import cosmology3

# the 0-parameter constants (the values live in `cl10/config.yaml`)
M_P = u10pred.M_PLANCK          # GeV
N_GEN = int(_config.CFG.physics.n_gen)   # the 0-parameter number of Cl(10) generators
LOOP = float(_config.CFG.physics.loop)   # the 0-parameter 1-loop QFT factor
WEAK_SCALE = float(_config.CFG.physics.weak_scale_gev)   # GeV (the 0-parameter WIMP-miracle reference)
OMEGA_OBS = float(_config.CFG.physics.omega_dm_h2)       # the observed dark-matter density (Omega h^2)


def _backend():
    return E._BACKEND


def F10():
    r"""The 0-parameter 10-loop factor $(1/16\pi)^{10}$ ($10 = n_{\rm gen}$)."""
    return float(LOOP ** N_GEN)


def relic_abundance(m_GeV):
    r"""The 0-parameter relic abundance for a mass $m$ (the V42 convention):
    $\Omega h^2 = \Omega_{\rm obs} \times (100/m)^3$, which gives
    $\Omega_{\rm obs}$ at the WIMP-miracle scale ($100$ GeV) by construction
    and the V42 Planck-soliton gap at $M_P$."""
    return float(OMEGA_OBS * (WEAK_SCALE / m_GeV) ** 3)


def dm_wimp():
    r"""The O7 scale: $M_P (1/16\pi)^{10} = 118$ GeV, the comparison with the
    $100$ GeV WIMP-miracle reference, the convention-normalised relic contrast,
    and the V42 Planck-soliton gap (for contrast).  **No state of the
    implemented spectrum carries this scale** - the carrier is the open problem
    O7 - and ``relic_abundance`` is normalised to the observed density at the
    reference scale by construction, so the abundance numbers are a contrast,
    not a prediction."""
    f10 = F10()
    m_dm = float(M_P * f10)                 # ~ 118 GeV
    m_over_weak = float(m_dm / WEAK_SCALE)  # ~ 1.18
    omega_dm = relic_abundance(m_dm)        # ~ 0.07
    omega_ratio = float(omega_dm / OMEGA_OBS)

    # the V42 Planck-soliton (for contrast): the 10^-52 gap
    v42 = cosmology3.dm_mass_scale()
    planck_omega = relic_abundance(M_P)
    planck_gap = float(planck_omega / OMEGA_OBS)

    one_line = (
        "the sub-Planck ladder scale M_P*(1/16pi)^10 = 118 GeV lands within a "
        "factor 1.18 of the WIMP/weak scale (the same 10-loop factor as the O8 "
        "hierarchy); no state of the implemented spectrum carries it - the "
        "carrier is the open problem O7 - and the V42 Planck-soliton gap "
        "(1e-52) is relocated, not resolved"
    )
    return {
        "backend": _backend(),
        "F10": f10,
        "n_gen": int(N_GEN),
        "M_P_GeV": float(M_P),
        "m_DM_GeV": m_dm,
        "weak_scale_GeV": float(WEAK_SCALE),
        "m_over_weak": m_over_weak,
        "omega_DM": omega_dm,
        "omega_obs": float(OMEGA_OBS),
        "omega_ratio_to_observed": omega_ratio,
        "v42_planck_omega": float(planck_omega),
        "v42_planck_gap": planck_gap,
        "v42_quantified_gap": float(v42["abundance_gap"]),
        "carrier_status": "open (O7): no state of the implemented spectrum carries the scale",
        "relic_convention": "Omega(m) = Omega_obs*(100 GeV/m)^3: normalised to the observed value at the reference scale by construction (a contrast device, not an abundance prediction)",
        "one_line": one_line,
    }


def verify(tol=_config.CFG.tolerances.strict):
    r"""The arithmetic of the O7 scale: (i) the ladder identity
    $m = M_P(1/16\pi)^{n_{\rm gen}}$; (ii) the scale lands in the WIMP window
    and within a factor 3 of the $100$ GeV reference (the reported coincidence);
    (iii) the relic contrast is the convention it claims to be (the identity
    $\\Omega(m)=\\Omega_{\rm obs}(100/m)^3$, so no independent abundance is
    asserted); (iv) the V42 gap agrees with ``cosmology3.dm_mass_scale`` across
    the modules.  Nothing here tests the DM **carrier**: that is O7."""
    r = dm_wimp()
    ladder = M_P * LOOP ** N_GEN
    omega_identity = abs(r["omega_DM"] - OMEGA_OBS * (WEAK_SCALE / r["m_DM_GeV"]) ** 3)
    cross_module = abs(r["v42_planck_gap"] - r["v42_quantified_gap"])
    ok = (
        abs(r["m_DM_GeV"] - ladder) <= tol * ladder
        and 10.0 <= r["m_DM_GeV"] <= 1000.0
        and (1.0 / 3.0) <= r["m_over_weak"] <= 3.0
        and omega_identity <= tol
        and cross_module <= tol * r["v42_planck_gap"]
    )
    return {"r": r, "ok": bool(ok)}


def main():
    import sys
    r = verify()
    d = r["r"]
    print("=" * 72)
    print("CL10 - the sub-Planck ladder scale (the 10-loop suppression of M_P); the carrier is open (O7)")
    print("=" * 72)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  10-loop factor (1/16pi)^10 = {d['F10']:.3e}  (10 = n_gen)")
    print(f"  ladder scale = M_P * (1/16pi)^10 = {d['m_DM_GeV']:.1f} GeV   (lands at the WIMP/weak scale)")
    print(f"  WIMP-miracle reference = {d['weak_scale_GeV']:.0f} GeV   ladder/reference = {d['m_over_weak']:.3f}")
    print(f"  relic contrast Omega h^2 = {d['omega_DM']:.4f}   (observed {d['omega_obs']}; normalised at the reference scale by construction, not a prediction)")
    print(f"  contrast ratio:  {d['omega_DM']:.3f} vs {d['omega_obs']}  (ratio {d['omega_ratio_to_observed']:.2f})")
    print(f"  V42 contrast:  the Planck-mass soliton gives Omega h^2 = {d['v42_planck_omega']:.2e}  "
          f"(gap {d['v42_planck_gap']:.1e}; V42 quantified {d['v42_quantified_gap']:.1e})")
    print(f"  carrier:  {d['carrier_status']}")
    print(f"  VERDICT:  {'PASS (the ladder scale M_P(1/16pi)^10 = 118 GeV; no state of the implemented spectrum carries it - the open problem O7)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
