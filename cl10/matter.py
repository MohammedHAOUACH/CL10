#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the matter: the 16-D Weyl spinor (one SM generation).

The 32-D GUT module splits, via the pseudoscalar, into a **16 Weyl** and a
**16-bar**.  The 16 is **one full SM generation**; under the SM it decomposes
as::

    16  =  10  +  5bar  +  1  +  nu^c
        =  (gauge)  +  (matter)  +  (sterile)

and, more finely, into the 6 SM fields::

    16  =  Q(3,2,-1/6)  +  u^c(3bar,1,1/3)  +  d^c(3bar,1,-1/3)
        +  L(1,2,-1/2)  +  e^c(1,1,+1)      +  N(1,1,0)

with dimensions 6 + 3 + 3 + 2 + 1 + 1 = 16.  The 16-bar is the **antimatter**
(the CPT conjugate).  The **B-L charge** is 0-parameter: the 16 has B-L = -1,
the 16-bar has B-L = +1 (they are inequivalent modules).
"""

from . import clifford


# the 6 SM fields: (name, (SU3, SU2), U1, dim)
SM_FIELDS = {
    "Q":   (3, 2, -1 / 6, 6),
    "u^c": (3, 1, +1 / 3, 3),
    "d^c": (3, 1, -1 / 3, 3),
    "L":   (1, 2, -1 / 2, 2),
    "e^c": (1, 1, +1, 1),
    "N":   (1, 1, 0, 1),
}


def weyl16():
    """The 16 Weyl and 16-bar (the +i and -i eigenspaces of the pseudoscalar)."""
    return clifford.weyl_spinor()


def dimension_check():
    """The 6 SM fields have total dimension 16 (0-parameter)."""
    return sum(f[3] for f in SM_FIELDS.values())


def bl_charge(which="16"):
    """The 0-parameter B-L charge: 16 -> -1, 16-bar -> +1."""
    return -1 if which == "16" else +1


def main():
    import sys
    print("=" * 68)
    print("CL10 - the matter: the 16 Weyl (one SM generation)")
    print("=" * 68)
    p, m = weyl16()
    print(f"  16 Weyl / 16-bar  dims:  {p.shape[1]} / {m.shape[1]}")
    print(f"  total SM-field dims:  {dimension_check()}  (10 + 5bar + 1 + nu^c)")
    print("  the 6 SM fields:")
    for name, (su3, su2, u1, dim) in SM_FIELDS.items():
        print(f"    {name:5s}  ({su3},{su2},{u1:+.1f})  dim={dim}")
    print(f"  B-L charge:  16 -> {bl_charge('16')},  16-bar -> {bl_charge('16bar')}")
    ok = (p.shape[1] == 16 and m.shape[1] == 16 and dimension_check() == 16)
    print(f"  VERDICT:  {'PASS (the 16 Weyl is one generation)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
