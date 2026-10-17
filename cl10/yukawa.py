#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the 0-parameter fermion-mass structure.

Two layers:

- **The structure** (V28, `yukawa_matrix`): the **10-Higgs 16->16-bar** (the
  grade-1 GUT vector is **off-diagonal** in the 16/16-bar split) gives the
  0-parameter **Dirac** mass term $16\,10\,16\bar$ — the **where** of the mass
  (all 16 Weyl get a Dirac mass).  With the **3 families**, the mass term gives
  a **3 x 3** matrix per fermion type.  At the **0-parameter** (family-symmetric)
  level the 3 families are **equal**, so the 3x3 is the **democratic** base
  $|1\rangle\langle1|$ : **rank 1**, spectrum $\{3,0,0\}$ — **one heavy** level
  + a **2-fold light** level (the 0-parameter **2-level** structure).

  The 3x3 is therefore **NOT** non-singular: a rank-1 3x3 has two zero
  eigenvalues.  The **three distinct** masses are **not** a 0-parameter
  statement — they require the **10-Higgs vev** to break the 2-fold light
  degeneracy (the O6 named residual).  This is the same statement as
  `cl10.yukawa_3x3` (V45/V49) and as Section 10 / Theorem 13 of the
  publication; the two modules share **one** implementation (this module
  delegates to `yukawa_3x3.mass_matrix_3x3`), so the code cannot drift from
  the paper again.

- **The hierarchy + mixing** (V29, `democratic_matrix`,
  `tribimaximal_mixing`, `hierarchy_residual`): the three families are the
  three two-planes of the $(1,3)$ 3-space.  Diagonalizing the **democratic**
  base gives the **tribimaximal** mixing, which **diagonalizes** it.  Of its
  angles only $\sin^2(12) = 1/3$ is **derived 0-parameter** (the unique heavy
  eigenvector); $\sin^2(23)$ and $\sin^2(13)$ sweep $[0, 2/3]$ under the $O(2)$
  freedom of the 2-fold light subspace.  The **mass hierarchy** (the three
  observed levels) is the **named residual** (the 10-Higgs vev breaking).
"""

import numpy as np

from . import engine as E
from . import family
from . import matter
from . import yukawa_3x3


def yukawa_matrix():
    """The 0-parameter 3x3 Yukawa/mass-matrix structure (the democratic base).

    Derivation: the 3-fold family symmetry (the 3 spatial directions $S_a$ are
    equivalent, see `yukawa_3x3.family_symmetry`) forces the 10-Higgs to couple
    the 3 families **equally**: the 0-parameter base is the **democratic** matrix
    $|1\\rangle\\langle1|$, of **rank 1** and spectrum $\\{3,0,0\\}$.  A rank-1
    3x3 matrix is **singular** (two zero eigenvalues): the three distinct masses
    are NOT 0-parameter.
    """
    mg = matter.mass_generation()
    fam = family.family_construct()
    n_gen = int(fam["n_family_candidate"])          # the 3 generations
    mm = yukawa_3x3.mass_matrix_3x3()               # the democratic base |1><1|
    offdiag = bool(mg["off_diagonal"])
    diag_zero = (mg["diag_16"] < 1e-6 and mg["diag_16bar"] < 1e-6)
    return {
        "n_gen": n_gen,
        "n_dim": int(mm["n"]),
        "rank": int(mm["rank"]),                    # 1 (the democratic base)
        "spectrum": list(mm["eigenvalues"]),        # {3, 0, 0}
        "n_distinct_levels": int(mm["n_distinct_levels"]),   # 2
        "n_heavy": int(mm["n_heavy"]),              # 1
        "n_light_2fold": int(mm["n_light_2fold"]),  # 2
        "democratic": bool(mm["democratic"]),
        "singular": bool(mm["rank"] < n_gen),       # True: rank 1 < 3
        "off_diagonal_Dirac": offdiag,
        "diag_zero": bool(diag_zero),
        "mass_term": mg["mass_term"],
        "residual": ("the 3-level mass hierarchy (the 10-Higgs vev breaking the "
                     "2-fold light degeneracy; the observed ratios) + the exact "
                     "mixing angles (CKM/PMNS) — the 0-parameter object fixes the "
                     "2-level structure, not the 3 levels"),
    }


def democratic_matrix():
    """The 0-parameter democratic 3x3 mass matrix (the 3 equal families).

    The three families (the three two-planes of the $(1,3)$ 3-space) are
    **equal** at the 0-parameter level: the GUT structure gives the
    **democratic** (all-ones) base, $|1\\rangle\\langle1|$, **rank 1**.
    """
    n = 3
    return E.Array_(np.ones((n, n), dtype=complex))


def tribimaximal_mixing():
    """The tribimaximal PMNS of the democratic base, with the honest status.

    The TBM diagonalises the democratic base **exactly** (correction: the
    previous implementation of `family.tribimaximal_3x3` had a wrong sign and did
    not diagonalise).  The honest status of each angle is **measured** by
    `yukawa_3x3.mixing_from_democratic`:

    - sin^2(12) = 1/3 is **genuinely derived** 0-parameter (the heavy eigenvector
      (1,1,1)/sqrt(3) is **unique**; the angle does not move under the O(2)
      freedom of the light subspace);
    - sin^2(23) and sin^2(13) **sweep [0, 2/3]** under the O(2) freedom of the
      2-fold degenerate subspace: the TBM values (1/2, 0) are the choice
      **adapted to the transposition (23)** of S3, not a consequence.
    """
    m = yukawa_3x3.mixing_3x3()
    obs = (0.307, 0.579, 0.022)     # the observed PMNS (world average)
    return {
        "sin2_12": m["sin2_12"], "sin2_23": m["sin2_23"], "sin2_13": m["sin2_13"],
        "obs": obs, "backend": E._BACKEND,
        "diagonalizes_democratic": m["diagonalizes_democratic"],
        "sin2_12_is_0param_derived": m["sin2_12_is_0param_derived"],
        "sin2_23_is_0param_derived": m["sin2_23_is_0param_derived"],
        "sin2_13_is_0param_derived": m["sin2_13_is_0param_derived"],
        "sin2_23_O2_range": m["sin2_23_O2_range"],
        "sin2_13_O2_range": m["sin2_13_O2_range"],
        "mixing_status": m["mixing_status"],
    }


def hierarchy_residual():
    """The 0-parameter base spectrum (2-level) -> the 3-level hierarchy is the residual.

    The democratic base is **rank 1**: spectrum {3, 0, 0} — one heavy level and
    a **2-fold degenerate light** level.  The **three-level** mass hierarchy is
    **not** 0-parameter: it requires the **10-Higgs vev** breaking.  The exact
    ratios are the **named residual**.
    """
    D = E.asnumpy(democratic_matrix())
    base = np.sort(np.abs(np.linalg.eigvals(D)))
    return {
        "base_spectrum": [float(x) for x in base],
        "degenerate": bool(base[1] < 1e-9 or abs(base[0] - base[1]) < 1e-9),
        "residual": "the 3-level mass hierarchy (the 10-Higgs vev breaking; the exact ratios)",
    }


def main():
    import sys
    y = yukawa_matrix()
    mix = tribimaximal_mixing()
    h = hierarchy_residual()
    obs = mix["obs"]
    err = (abs(mix["sin2_12"] - obs[0]) + abs(mix["sin2_23"] - obs[1])
           + abs(mix["sin2_13"] - obs[2])) / 3
    print("=" * 68)
    print("CL10 - the 0-parameter mass structure (V28 structure + V29 mixing/hierarchy)")
    print("=" * 68)
    print(f"  backend (GPU):  {mix['backend']}")
    print(f"  [V28] the {y['n_gen']}x{y['n_gen']} democratic mass matrix:  rank {y['rank']} "
          f"(spectrum {['%.3g' % x for x in y['spectrum']]}: {y['n_heavy']} heavy + "
          f"{y['n_light_2fold']}-fold light), Dirac off-diagonal={y['off_diagonal_Dirac']}")
    print("  [V29] the 0-parameter democratic base -> the TBM that diagonalizes it")
    print("        (only sin^2(12)=1/3 is derived 0-parameter; sin^2(23), sin^2(13) are an O(2) choice):")
    print(f"    sin^2(12)={mix['sin2_12']:.3f} (obs {obs[0]:.3f}), "
          f"sin^2(23)={mix['sin2_23']:.3f} (obs {obs[1]:.3f}), "
          f"sin^2(13)={mix['sin2_13']:.3f} (obs {obs[2]:.3f})")
    print(f"  the mean 0-parameter vs observed mixing deviation:  {err:.3f}")
    print(f"  the 0-parameter base spectrum (2-level):  {['%.3g' % x for x in h['base_spectrum']]}")
    print(f"  named residual:  {h['residual']}")
    ok = (y["rank"] == 1 and y["n_distinct_levels"] == 2 and y["off_diagonal_Dirac"]
          and err < 0.1 and h["degenerate"])
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
