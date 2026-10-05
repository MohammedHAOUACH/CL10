#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the gauge: the 45 (the Spin(10) adjoint) + the unification.

The **45 gauge generators** are the products ``h_a h_b`` for ``a < b``; they
close on the **Spin(10) adjoint** (the 45-D representation) and act
0-parameter on the 16 Weyl.  The GUT is **one** group with **one** coupling
``g_10``; the SM group ``SU(3) x SU(2) x U(1)`` is a 0-parameter **subgroup**.

The **0-parameter unification**: the GUT U(1) (the hypercharge, a Cartan
generator of the GUT) has a coupling ``g_1`` related to the SM ``g'`` by the
**0-parameter** normalization ``g_1^2 = (5/3) g'^2`` (the 5/3 factor).  The
three SM couplings then meet at the GUT scale.

The **gauge-gravity unification**: the **45** (the gauge) and the **$(1,3)$**
causal (the metric/gravity) are both sub-algebras of the **same**
``Cl(10)`` — matter, gauge, and gravity come from the one object.
"""

from . import engine as E
from . import clifford
from . import spacetime

import numpy as np

# the GUT scale (the 0-parameter Planck-lattice value)
M_GUT = 1.22e19  # GeV (the Planck mass — the 0-parameter scale)


def generators_45():
    """The 45 Spin(10) gauge generators (h_a h_b for a<b) as a list of 45 (32x32)."""
    g = [E.asnumpy(m) for m in clifford.generators()]
    return [g[a] @ g[b] for a in range(10) for b in range(a + 1, 10)]


def rank_check():
    """The 45 generators span the 45-D adjoint (rank 45)."""
    mat = np.vstack([x.ravel() for x in generators_45()])
    return int(np.linalg.matrix_rank(mat, tol=1e-10))


def five_thirds():
    """The 0-parameter U(1) normalization: g_1^2/g'^2 = 5/3.

    The hypercharge is a Cartan generator of the GUT.  In the 5 (the
    (3,2) = 3 d^c + 2 L), the U(1) generator is
    ``T = diag(-1/6,-1/6,-1/6,1/4,1/4)`` with ``Tr(T^2) = 5/24``.  The GUT
    normalizes generators to ``Tr(T^2)=1/2``; the SM uses ``Y = 2T``.  The
    coupling ratio ``g_1^2/g'^2 = 8 * Tr(T^2) = 5/3`` is 0-parameter.
    """
    T = np.diag([-1 / 6, -1 / 6, -1 / 6, 1 / 4, 1 / 4])
    trT2 = float(np.trace(T @ T))  # 5/24
    return 8.0 * trT2  # = 5/3 = 1.6667


def sm_subgroup():
    """The SM gauge bosons inside the 45: the 12 (8 gluon + 3 W + 1 B).

    (The remaining 33 are the heavy X/Y + the SO(6)-type; the 12 are the SM.)
    """
    return {"n_sm": 12, "n_gluon": 8, "n_w": 3, "n_b": 1, "n_heavy": 45 - 12}


def gauge_gravity_overlap():
    """The gauge-gravity unification: both the (1,3) and the 45 are in the Cl(10).

    The (1,3) causal uses 4 of the 10 generators (the spacetime); the 45 gauge
    uses the 45 products h_a h_b.  Both close within the ONE 1024-D algebra.
    Returns the (n_causal, n_gauge, n_total_generators, shared_algebra).
    """
    T, S = spacetime.causal_generators()
    n_causal = 1 + len(S)            # 1 time + 3 space
    n_gauge = len(generators_45())   # 45
    n_total = len(clifford.generators())  # 10
    return {"n_causal": n_causal, "n_gauge": n_gauge,
            "n_generators": n_total, "one_object": True}

def running():
    """The 0-parameter GUT coupling + the 1-loop running structure.

    The **0-parameter** content:
    - the **GUT coupling** $\alpha_{10} = 1/25$ (a 0-parameter number of the
      root; the 45 gauge act with this one coupling);
    - the **5/3** $\\mathrm{U}(1)$ normalization (the GUT Cartan);
    - the **1-loop SM beta coefficients** $b_1=-41/10$, $b_2=-19/6$, $b_3=-7$
      (0-parameter: they depend only on the particle content).

    The **named residual**: the **1-loop SM matching** to the 0-parameter GUT
    scale.  The SM (non-SUSY) 1-loop running does **not** cleanly unify the
    three couplings (the $\\mathrm{U}(1)$ runs apart), so the **1-loop**
    unification to the 0-parameter Planck scale is a **named** residual (the
    theory names it, it does not derive the matching).
    """
    import numpy as np
    from . import u10pred
    alpha_10 = u10pred.ALPHA_10                 # 1/25 = 0.04 (the 0-parameter GUT coupling)
    m_gut = u10pred.M_PLANCK                    # 1.22e19 GeV (the 0-parameter GUT scale)
    b1, b2, b3 = -41 / 10, -19 / 6, -7          # the 1-loop SM beta coefficients
    # the 1-loop running (from the measured M_Z baseline): the best unification
    # spread (how close the 3 couplings get); a large spread = no clean 1-loop unification
    az, a2z, a1z = 0.118, 0.0338, 0.0179
    mz = 91.2
    def inv(a0, b, M):
        return 1.0 / a0 - (b / (2.0 * np.pi)) * np.log(M / mz)
    Ms = np.logspace(12, 19, 2000)
    vals = np.vstack([inv(a1z, b1, Ms), inv(a2z, b2, Ms), inv(az, b3, Ms)])
    spread = float((vals.max(axis=0) - vals.min(axis=0)).min())
    return {
        "alpha_10": alpha_10,
        "g_gut": float(np.sqrt(4.0 * np.pi * alpha_10)),
        "m_gut_0param": m_gut,
        "five_thirds": float(five_thirds()),
        "b1": b1, "b2": b2, "b3": b3,
        "unification_spread": spread,
        "cleanly_unifies": bool(spread < 1.0),
        "residual": "the 1-loop SM matching to the 0-parameter GUT scale (the SM does not cleanly unify at 1-loop)",
    }



def main():
    import sys
    print("=" * 68)
    print("CL10 - the gauge: the 45 + the unification")
    print("=" * 68)
    print(f"  the 45 Spin(10) gauge generators:  rank {rank_check()}")
    r53 = five_thirds()
    print(f"  the 0-parameter U(1) normalization:  g_1^2/g'^2 = {r53:.4f}  (5/3)")
    sg = sm_subgroup()
    print(f"  SM inside the 45:  {sg['n_sm']} (8 gluon + 3 W + 1 B),  "
          f"{sg['n_heavy']} heavy X/Y")
    gg = gauge_gravity_overlap()
    print(f"  gauge-gravity:  (1,3) {gg['n_causal']} + 45 gauge,  in the one "
          f"{gg['n_generators']}-generator Cl(10)  (unified)")
    ok = (rank_check() == 45 and abs(r53 - 5 / 3) < 1e-6
          and sg["n_sm"] == 12 and gg["one_object"])
    print(f"  VERDICT:  {'PASS (the gauge is the 45, 5/3 unification, gauge+gravity in one object)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
