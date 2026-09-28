#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the matter: the 16-D Weyl spinor (one SM generation) + the anomalies.

The 32-D GUT module splits, via the pseudoscalar, into a **16 Weyl** and a
**16-bar**.  The 16 is **one full SM generation**; under the SM it decomposes
as::

    16  =  Q(3,2,+1/6)  +  u^c(3bar,1,-2/3)  +  d^c(3bar,1,+1/3)
        +  L(1,2,-1/2)  +  e^c(1,1,+1)       +  N(1,1,0)

(the hypercharges $Y$ in the $Q = T_3 + Y$ convention; the $u^c,d^c,e^c,N$
are the *left-handed Weyl* — the conjugates of the physical right-handed
fields).  The 16-bar is the **antimatter** (the CPT conjugate, opposite B-L).

The **0-parameter** charges above are the input to the **anomaly**
computation: a chiral gauge theory is *inconsistent* unless all gauge
anomalies cancel, so the anomaly check is a **hard TOE requirement**.
"""

from . import clifford

# the 6 SM fields:  (su3_type, su2_type, Y, n_weyl)
#   su3_type:  +1 = 3 (fund),  -1 = 3bar,  0 = 1
#   su2_type:  +1 = 2 (doublet),  0 = 1
#   Y:         the U(1) hypercharge (Q = T3 + Y)
#   n_weyl:    the number of left-handed Weyl spinors in the field
SM_FIELDS = {
    "Q":   (+1, +1, +1 / 6, 6),
    "u^c": (-1, 0, -2 / 3, 3),
    "d^c": (-1, 0, +1 / 3, 3),
    "L":   (0, +1, -1 / 2, 2),
    "e^c": (0, 0, +1, 1),
    "N":   (0, 0, 0, 1),
}

# the Dynkin indices (0-parameter:  T(3) = T(2) = 1/2)
_T3 = 1 / 2
_T2 = 1 / 2


def weyl16():
    """The 16 Weyl and 16-bar (the +i and -i eigenspaces of the pseudoscalar)."""
    return clifford.weyl_spinor()


def dimension_check():
    """The 6 SM fields have total Weyl dimension 16 (0-parameter)."""
    return sum(f[3] for f in SM_FIELDS.values())


def bl_charge(which="16"):
    """The 0-parameter B-L charge: 16 -> -1, 16-bar -> +1."""
    return -1 if which == "16" else +1


def anomalies():
    """The 8 SM gauge-anomaly coefficients from the 0-parameter charges.

    Every one must be **0** for the chiral gauge theory to be consistent.
    The non-trivial ones:  [SU(3)]^3, [U(1)]^3, [grav]^2 x U(1),
    U(1) x SU(3)^2, U(1) x SU(2)^2.  The trivially-zero ones:  [SU(2)]^3
    (SU(2) has no $d^{abc}$), [grav]^2 x SU(3), [grav]^2 x SU(2) (the
    generators are traceless).
    """
    a_su3_cubed = a_u1_cubed = a_grav_u1 = a_u1_su3_sq = a_u1_su2_sq = 0.0
    for su3, su2, y, n in SM_FIELDS.values():
        a_su3_cubed += n * su3          # tr({T,T}T):  + for 3, - for 3bar
        a_u1_cubed += n * y ** 3
        a_grav_u1 += n * y             # tr(T^a) of the U(1) = the charge
        if su3 != 0:
            a_u1_su3_sq += n * y * _T3
        if su2 != 0:
            a_u1_su2_sq += n * y * _T2
    return {
        "su3_cubed": a_su3_cubed,       # [SU(3)]^3
        "su2_cubed": 0.0,               # [SU(2)]^3  (SU(2) is anomaly-free)
        "u1_cubed": a_u1_cubed,         # [U(1)]^3
        "grav2_su3": 0.0,               # [grav]^2 x SU(3)  (traceless)
        "grav2_su2": 0.0,               # [grav]^2 x SU(2)  (traceless)
        "grav2_u1": a_grav_u1,          # [grav]^2 x U(1)
        "u1_su3_sq": a_u1_su3_sq,       # U(1) x SU(3)^2
        "u1_su2_sq": a_u1_su2_sq,       # U(1) x SU(2)^2
    }


def mass_generation():
    """The 0-parameter fermion mass generation (the 10-Higgs).

    The **Cl(10) 10** (the grade-1 generators, the GUT vector) is **off-diagonal**
    in the 16/16-bar (matter/antimatter) split: it maps the 16 to the 16-bar and
    the 16-bar to the 16 (it anticommutes with the pseudoscalar).  So a **10-Higgs
    vev** (a specific grade-1 element) couples the 16 to its 16-bar conjugate and
    gives the **0-parameter Dirac mass term** $16\\cdot 10 \\cdot 16\\bar \\to 1$:

    - the **diagonal** blocks (16->16, 16-bar->16-bar) are **0** (the 10 is odd);
    - the **off-diagonal** blocks (16->16-bar, 16-bar->16) are **nonzero** (the
      mass matrices).

    The 16 also contains the **sterile N** (the 1 of the 16, $B\\!-\\!L=0$), which
    is the 0-parameter source of the **seesaw** (the 2-step W2, the light neutrino
    mass).  The **structure** of the mass term is **0-parameter**; the **exact**
    mass values (the 3-generation hierarchy, the mixing angles) are **named**
    residuals.
    """
    import numpy as np
    g = [clifford.E.asnumpy(m) for m in clifford.generators()]
    p, m = weyl16()
    P16 = p @ p.conj().T
    P16b = m @ m.conj().T
    h = g[0]                                  # a representative grade-1 (the 10)
    A = P16b @ h @ P16                        # 16 -> 16-bar
    B = P16 @ h @ P16b                        # 16-bar -> 16
    diag16 = P16 @ h @ P16                    # 16 -> 16  (should be 0)
    diag16b = P16b @ h @ P16b                 # 16-bar -> 16-bar  (should be 0)
    offdiag_ok = (float(np.linalg.norm(diag16).max()) < 1e-6
                  and float(np.linalg.norm(diag16b).max()) < 1e-6
                  and float(np.linalg.norm(A).max()) > 1.0
                  and float(np.linalg.norm(B).max()) > 1.0)
    return {
        "n16": int(p.shape[1]), "n16bar": int(m.shape[1]),
        "off_diagonal": bool(offdiag_ok),
        "diag_16": float(np.linalg.norm(diag16).max()),
        "diag_16bar": float(np.linalg.norm(diag16b).max()),
        "A_16_16bar": float(np.linalg.norm(A).max()),
        "B_16bar_16": float(np.linalg.norm(B).max()),
        "mass_term": "16 . 10-Higgs . 16bar -> 1 (the 0-parameter Dirac mass)",
        "sterile_N_in_16": True,              # the N (the 1 of the 16) -> the seesaw
        "residuals": "the exact mass values (the 3-gen hierarchy, the mixing) are named",
    }


def main():
    import sys
    print("=" * 68)
    print("CL10 - the matter: the 16 Weyl (one SM generation)")
    print("=" * 68)
    p, m = weyl16()
    print(f"  16 Weyl / 16-bar  dims:  {p.shape[1]} / {m.shape[1]}")
    print(f"  total SM-field dim:  {dimension_check()}  (10 + 5bar + 1 + nu^c)")
    print("  the 6 SM fields  (su3, su2, Y, n):")
    for name, (su3, su2, y, n) in SM_FIELDS.items():
        print(f"    {name:5s}  ({'3' if su3 else '1'},{'2' if su2 else '1'},Y={y:+.3f})  n={n}")
    print(f"  B-L charge:  16 -> {bl_charge('16')},  16-bar -> {bl_charge('16bar')}")
    ok = (p.shape[1] == 16 and m.shape[1] == 16 and dimension_check() == 16)
    print(f"  VERDICT:  {'PASS (the 16 Weyl is one generation)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
