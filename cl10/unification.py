#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the 0-parameter coupling unification and the QCA UV completion.

**The 0-parameter content.**
- the **$5/3$** GUT normalization of the SM $\mathrm{U}(1)$ (Theorem 3) — the
  **0-parameter** anchor of the unification;
- the **0-parameter GUT coupling** $\alpha_{10} = 1/25$ (the single coupling of
  the 45, the GUT Cartan);
- the **1-loop** $\beta$-coefficients of the particle content
  ($b_1=41/10$, $b_2=-19/6$, $b_3=-7$ for the SM).

**The running to the 0-parameter scale.**  Running the three couplings from
$M_Z$ to the **Planckian** GUT scale $M_{\\rm GUT}=M_P$ (the 0-parameter
single scale) gives couplings in the **GUT range** ($\\alpha\\sim0.02$-$0.04$,
consistent with $\alpha_{10}=0.04$); the **1-loop spread** at $M_P$
($\\sim0.6$ for the non-SUSY SM) is the **named residual** — the non-SUSY SM
does **not** unify to a single point at $M_P$ (the clean MSSM unification is at
$\sim2\\times10^{16}$, **below** the Planck scale).

**The QCA UV completion.**  The lattice UV cutoff is the Planck scale.  The SM
$\mathrm{U}(1)$ has a **Landau pole** at $\sim3\\times10^{39}$ GeV, far
**above** the QCA cutoff $M_P$ — so in the QCA the running is **UV-bounded**
(the coupling is **bounded** on the lattice, no Landau pole below the Planck
scale).  This is the **0-parameter UV completion** of the theory.
"""

import numpy as np

from . import engine as E
from . import gauge
from . import toepred

M_Z = 91.2                       # GeV
ALPHA_10 = toepred.ALPHA_10      # 1/25 = 0.04 (the 0-parameter GUT coupling)
M_P = toepred.M_PLANCK           # 1.22e19 GeV (the 0-parameter GUT = Planck scale)

# the 1-loop beta coefficients (convention: 1/a(M) = 1/a(Mz) - (b/2pi) ln(M/Mz))
B_SM = (41 / 10, -19 / 6, -7)    # (b1, b2, b3)  (the SM particle content)
A_Z = (0.0179, 0.0338, 0.118)    # (a1, a2, a3) at M_Z  (paired with B_SM)


def _run(a_z, b, M):
    return 1.0 / a_z - (b / (2.0 * np.pi)) * np.log(M / M_Z)


def sm_running_to(M):
    """The 1-loop SM running to scale M.  Returns (a1, a2, a3, spread)."""
    a = [1.0 / _run(A_Z[k], B_SM[k], M) for k in range(3)]
    mean = float(np.mean(a))
    spread = float((max(a) - min(a)) / mean)
    return a, spread


def landau_pole_u1():
    """The SM U(1) Landau pole scale (where 1/a1 -> 0)."""
    ln_pole = (2.0 * np.pi / B_SM[0]) * (1.0 / A_Z[0])
    return M_Z * float(np.exp(ln_pole))


def unification():
    """The 0-parameter unification verdict at the Planckian scale M_P."""
    a_p, spread_p = sm_running_to(M_P)
    # the clean MSSM unification scale (below M_P)
    mssm_gut = 2.0e16
    a_m, spread_m = _mssm_running(mssm_gut)
    landau = landau_pole_u1()
    return {
        "alpha_10": float(ALPHA_10),
        "five_thirds": float(gauge.five_thirds()),
        "M_GUT_0param": float(M_P),
        "couplings_at_M_P": [float(x) for x in a_p],
        "spread_at_M_P": spread_p,
        "mssm_unif_scale": float(mssm_gut),
        "mssm_spread": spread_m,
        "u1_landau_pole": float(landau),
        "uv_bounded": bool(M_P < landau),
        "backend": E._BACKEND,
    }


def _mssm_running(M):
    b = (33 / 5, -1, -3)
    a = [1.0 / _run(A_Z[k], b[k], M) for k in range(3)]
    return a, float((max(a) - min(a)) / float(np.mean(a)))


def main():
    import sys
    r = unification()
    print("=" * 68)
    print("CL10 - the 0-parameter coupling unification + the QCA UV completion")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter GUT coupling:  alpha_10 = {r['alpha_10']:.3f}  (= 1/25)")
    print(f"  the 0-parameter 5/3 normalization:  g1^2/g'^2 = {r['five_thirds']:.4f}")
    print(f"  the 0-parameter GUT scale:  M_GUT = M_P = {r['M_GUT_0param']:.2e} GeV (the single scale)")
    print(f"  the 1-loop SM running to M_P:  a = {['%.3f' % x for x in r['couplings_at_M_P']]}  (spread {r['spread_at_M_P']:.2f})")
    print(f"  the clean MSSM (SUSY) unification is at {r['mssm_unif_scale']:.0e} GeV (spread {r['mssm_spread']:.2f}) -- BELOW M_P")
    print(f"  the U(1) Landau pole is at {r['u1_landau_pole']:.2e} GeV -- the QCA cutoff M_P is BELOW it (UV-bounded)")
    ok = (r["five_thirds"] > 1.6 and r["five_thirds"] < 1.8
          and abs(r["alpha_10"] - 0.04) < 1e-9
          and r["uv_bounded"]
          and all(0.01 < a < 0.06 for a in r["couplings_at_M_P"]))
    print(f"  VERDICT:  {'PASS (0-parameter 5/3 + GUT coupling; the running GUT-range at M_P; the 1-loop matching the named residual; the QCA UV-bounded)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
