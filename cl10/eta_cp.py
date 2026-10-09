#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the baryon asymmetry as a product of 0-parameter **factors**
(V55, W5).

The W5 residual (Section 13): the **CP-violating factor** $\eta_{CP}$ in
heavy-$N$ decay (V40 takes $\eta_{CP} = 10^{-3}$, the observed range
$10^{-3}\text{–}10^{-2}$).  This module presents the baryon asymmetry as
a **product of 0-parameter factors** (the physical GUT-leptogenesis
structure):

    \Delta B \;=\; \frac{28}{79}\;\varepsilon_{CP}\;D

with $\frac{28}{79}$ = the **SM sphaleron** conversion — an **input of the SM
sector** (it counts the SM left-handed degrees of freedom, $n_{\rm gen}
$\times$ the $\mathrm{SU}(2)$ doublets), **not** a number of the Cl(10) root:
it is 0-parameter **only conditional** on the full derivation of the SM
content (V23, a named residual),
$\varepsilon_{CP} = g_0^4/16\pi$ = the **1-loop CP-odd rate**
($g_0 = 1/8 = 1/2^3$, the 0-parameter GUT coupling -- its *reading* is the
audit's mapping convention), and $D = (1/16\pi)^2$ = the **2-loop
dilution** (the washout in the 2-step GUT breaking).  The **factors** are
0-parameter (root / QFT / SM), but the **result is not a 0-parameter
prediction**:

- In the **0-parameter QCA alone** the CP-odd theorem gives
  $\varepsilon_{CP} = 0$ (every CP-odd invariant vanishes: the QCA
  satisfies none of the 3 Sakharov conditions) -- the asymmetry is
  nonzero **only** through the **named irreducible CP phase** $\delta$
  (an $\mathcal{O}(1)$ input, the minimal irreducible input the theorem
  demands).
- The 2-digit **ratio** $1.12$ against the observed $\Delta B =
  6.1e-10$ is a **named residual** whose value is **read from the
  data**: the "12% match" depends entirely on the non-derived choice
  $g_0 = 1/8$.

Result:

    \Delta B \;=\; \frac{28}{79}\,\frac{g_0^4}{16\pi}\,(1/16\pi)^2
    \;\approx\; 6.8e-10

**Named residuals**: the CP phase $\delta$ ($\mathcal{O}(1)$) and the
2-digit ratio 1.12.  The **factors** of W5 are resolved (loop factors,
not a bare named input); their **values** remain named.
"""

import numpy as np

from . import engine as E
from . import walls

# the 0-parameter constants
G0 = walls.G0                  # 1/8 (the 0-parameter GUT coupling)
LOOP = 1.0 / (16.0 * np.pi)   # the 0-parameter 1-loop QFT factor
SPHALERON = 28.0 / 79.0        # the SM sphaleron conversion: an input of the SM
#                                sector (the LH d.o.f. count), NOT a root number;
#                                0-parameter only conditional on V23 (named residual)
DELTA_B_OBS = 6.12e-10         # observed baryon asymmetry (Planck 2018; PDG 2024)
ETA_CP_V40 = 1.0e-3            # the V40 named value (the observed range 1e-3..1e-2)


def _backend():
    return E._BACKEND


def cp_odd_loop_factor():
    r"""The 0-parameter CP-odd **rate** loop factor $g_0^4/(16\pi)$
    ($g_0 = 1/8$)."""
    return float(G0 ** 4 / (16.0 * np.pi))


def eta_cp():
    r"""The W5 result: the baryon asymmetry as a **product of factors of three
    provenances** — root / QFT / SM (the physical GUT-leptogenesis structure).

    $$\Delta B \;=\; \frac{28}{79}\;\varepsilon_{CP}\;D$$

    with the **physical factors** (three provenance classes, NOT all of the
    root):

    - $\frac{28}{79}$: the **SM sphaleron** conversion — provenance **SM**
      (an input of the SM sector: the left-handed d.o.f. count; it is a root
      number only if the SM content is derived, V23 = named residual).  From the SM gauge
      group).
    - $\varepsilon_{CP} = g_0^4/(16\pi)$: the **1-loop CP-odd decay
      rate**.  $g_0 = 1/8 = 1/2^3$ (the 0-parameter GUT coupling, Theorem
      6; the **3** = the spatial dimension; the *reading* is the audit's
      mapping convention) is the CP-odd **4-coupling** product
      ($g_0^4$); $16\pi$ = the 4D 1-loop factor.
    - $D = (1/16\pi)^2$: the **2-loop dilution** (the washout of the
      asymmetry in the 2-step GUT breaking).

    The **0-parameter kinematics**: the out-of-equilibrium decay
    $\Gamma/H \sim 6\times10^{-4} < 1$ and the 2-step GUT breaking
    ($M_{N1}, M_{N2} \sim 10^{14}$ GeV).

    The 0-parameter QCA **alone** gives $\varepsilon_{CP} = 0$ (the
    CP-odd theorem): the asymmetry is nonzero **only** through the
    **named** irreducible CP phase $\delta$ ($\mathcal{O}(1)$).

    $$\Delta B \;=\; \frac{28}{79}\,\frac{g_0^4}{16\pi}\,(1/16\pi)^2
    \;\approx\; 6.8e-10,$$

    a **ratio 1.12** to the observed $6.1e-10$ (a **12%** match) -- the
    2-digit ratio is a **named residual** (read from the data), as is
    the CP phase.  The factors are 0-parameter in three **different** senses:
    numbers of the root ($16\pi$, the loop counts), QFT loop factors
    ($\varepsilon_{CP}$, $D$), or an **input of the SM sector** ($28/79$,
    conditional on V23).

    Returns the physical structure, the $\Delta B$ match, and the V40
    comparison.
    """
    eps_cp = cp_odd_loop_factor()                 # g0^4/(16pi) ~ 4.86e-6 (1-loop)
    d_dilution = float(LOOP ** 2)                 # (1/16pi)^2 ~ 3.96e-4 (2-loop)
    delta_b_0param = float(SPHALERON * eps_cp * d_dilution)
    ratio = float(delta_b_0param / DELTA_B_OBS)
    g0_is_1over2cubed = bool(abs(G0 - (1.0 / 2.0) ** 3) < 1e-12)  # 3 = spatial

    # the V40 comparison (the named eta_CP = 1e-3)
    delta_b_v40 = float(SPHALERON * eps_cp * ETA_CP_V40)
    ratio_v40 = float(delta_b_v40 / DELTA_B_OBS)

    one_line = (
        "the baryon asymmetry Delta_B = (28/79) x eps_CP x D is a product of "
        "FACTORS of three provenances (28/79 the SM sphaleron = an input of the "
        "SM sector, a root number only conditional on V23; eps_CP = g_0^4/16pi the "
        "1-loop CP-odd rate, g_0 = 1/2^3 a mapping convention; D = (1/16pi)^2 "
        "the 2-loop dilution): 6.8e-10, ratio 1.12 to the observed 6.1e-10.  "
        "The 0-parameter QCA alone gives eps_CP = 0 (the CP-odd theorem); the "
        "NAMED residuals are the CP phase (O(1)) and the 2-digit ratio 1.12"
    )
    return {
        "backend": _backend(),
        "loop_factor": float(LOOP),
        "structure": {
            "sphaleron": float(SPHALERON),                 # 28/79 (the SM)
            "eps_CP": eps_cp,                              # g0^4/(16pi) (1-loop)
            "D_dilution": d_dilution,                      # (1/16pi)^2 (2-loop)
            "g0": float(G0),
            "g0_is_1_over_2_cubed": g0_is_1over2cubed,      # 3 = the spatial dim
            "loop_count": 3,                                # 1 (CP) + 2 (dilution)
        },
        "eta_CP_0param": d_dilution,
        "eta_CP_observed": float(ETA_CP_V40),
        "eta_CP_over_observed": float(d_dilution / ETA_CP_V40),
        "g0": float(G0),
        "cp_odd_loop_factor": eps_cp,
        "sphaleron": float(SPHALERON),
        "delta_B_0param": delta_b_0param,
        "delta_B_observed": float(DELTA_B_OBS),
        "delta_B_ratio": ratio,
        "delta_B_v40": delta_b_v40,
        "delta_B_ratio_v40": ratio_v40,
        "one_line": one_line,
    }


def verify(tol=1e-12):
    """0-parameter checks on the FACTORS: (i) the 0-parameter eta_CP loop
    factor is within a factor 3 of the V40 value (1e-3); (ii) the
    product gives the right ORDER (the 2-digit ratio 1.12 -- the named
    residual); (iii) it **improves** on V40 (the ratio is smaller)."""
    r = eta_cp()
    ok = (
        (1.0 / 3.0) <= r["eta_CP_over_observed"] <= 3.0
        and (1.0 / 3.0) <= r["delta_B_ratio"] <= 3.0
        and r["delta_B_ratio"] < r["delta_B_ratio_v40"]
    )
    return {"r": r, "ok": bool(ok)}


def main():
    import sys
    r = verify()
    d = r["r"]
    print("=" * 72)
    print("CL10 - the baryon asymmetry as a product of 0-parameter factors (W5)")
    print("=" * 72)
    st = d["structure"]
    print(f"  backend (GPU):  {d['backend']}")
    print("  Delta_B = (28/79) x eps_CP x D, every factor 0-parameter:")
    print(f"    - sphaleron   = 28/79 = {st['sphaleron']:.4f}   (the SM)")
    print(f"    - eps_CP      = g_0^4/16pi = {st['eps_CP']:.3e}   (1-loop; g_0 = 1/2^3)")
    print(f"    - D (dilution)= (1/16pi)^2 = {st['D_dilution']:.3e}   (2-loop)")
    print(f"    - loop counts = 1 (CP) + 2 (dilution) = {st['loop_count']}")
    print(f"  product Delta_B = {d['delta_B_0param']:.3e}   (observed {d['delta_B_observed']:.1e};  2-digit ratio {d['delta_B_ratio']:.3f} = named residual)")
    print(f"  V40 comparison:  Delta_B = {d['delta_B_v40']:.3e}  (ratio {d['delta_B_ratio_v40']:.3f})")
    print(f"  VERDICT:  {'PASS (factors 0-parameter; right order 1.12; QCA alone gives eps_CP = 0 by the CP-odd theorem; the CP phase + 2-digit ratio are the named residuals)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
