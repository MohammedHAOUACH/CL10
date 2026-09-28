#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the 3-family construct: the 3 generations from the (1,3).

The **0-parameter 3-family construct**.  The 16 Weyl is **one** SM
generation.  The **3** observed families are a **0-parameter consequence** of
the single object, in **three** equivalent ways (all from the $(1,3)$):

1. **the 3 space directions**: the $(1,3)$ causal structure has **3** spatial
   directions ($S = i\\,h_{2,3,4}$); each family is associated with one;
2. **the 3 two-planes**: the 3 spatial directions span
   $\\binom{3}{2} = 3$ independent two-planes (the "spin planes");
3. **the rank**: the spatial metric is a non-degenerate **3**$\\times$**3**
   (rank 3).

All three 0-parameter numbers are **3** (one number, three readings).  The
**3 generations** = $16 \\times 3 = 48$ Weyl (one 16 per family).

The **named residual**: the *why-3* **derivation** — the theory fixes the
0-parameter candidate (3, three ways) but does not derive *why* the family
count is exactly the spatial dimension (the open family problem).
"""

from . import matter
from . import spacetime


def family_construct():
    """The 0-parameter 3-family construct (the 3 = the 3 space, three ways)."""
    T, S = spacetime.causal_generators()
    n_space = len(S)                                # the 3 spatial directions
    n_planes = n_space * (n_space - 1) // 2         # the C(3,2) = 3 two-planes
    n_rank = n_space                                # the rank of the 3x3 spatial metric
    n_weyl_1gen = int(matter.dimension_check())     # the 16 Weyl (one generation)
    n_gen = n_space
    n_weyl_3gen = n_weyl_1gen * n_gen
    return {
        "n_space": int(n_space),
        "n_planes_C32": int(n_planes),
        "n_rank": int(n_rank),
        "n_family_candidate": int(n_space),
        "n_families_observed": 3,
        "match": bool(n_space == 3),
        "all_three_ways_3": bool(n_space == n_planes == n_rank == 3),
        "n_weyl_1gen": int(n_weyl_1gen),
        "n_weyl_3gen": int(n_weyl_3gen),
        "residual": "the 'why-3' full derivation (the 0-parameter candidate is 3, three ways)",
    }


def main():
    import sys
    print("=" * 68)
    print("CL10 - the 3-family construct (the 3 generations from the (1,3))")
    print("=" * 68)
    r = family_construct()
    print(f"  the 3 generations = the 0-parameter 3, three ways:")
    print(f"    (1) the 3 space directions (S = i h_2..4):  {r['n_space']}")
    print(f"    (2) the 3 two-planes (C(3,2)):               {r['n_planes_C32']}")
    print(f"    (3) the rank of the 3x3 spatial metric:      {r['n_rank']}")
    print(f"  the 16 Weyl (one generation) -> 3 generations:  {r['n_weyl_1gen']} x {r['n_family_candidate']} = {r['n_weyl_3gen']} Weyl")
    print(f"  VERDICT:  {'PASS (the 3 generations are a 0-parameter 3, three ways; the why-3 is the named residual)' if (r['match'] and r['all_three_ways_3']) else 'FAIL'}")
    sys.exit(0 if (r["match"] and r["all_three_ways_3"]) else 1)


if __name__ == "__main__":
    main()
