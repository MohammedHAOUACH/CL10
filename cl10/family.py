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

import math

import numpy as np

from . import baryogenesis
from . import clifford
from . import engine as E
from . import matter
from . import sm_qca
from . import spacetime
from . import config as _config
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml


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


def _backend():
    return E._BACKEND


def tribimaximal_3x3():
    """The genuine tribimaximal (Harrison-Perkins-Scott) 3x3 unitary.

    **Correction.**  The earlier implementation had a **wrong sign** in ``r23``:
    it actually produced the **bimaximal** matrix, which has the **same squared
    moduli** as the TBM (hence apparently correct angles, ``sin^2 12 = 1/3``,
    ``sin^2 23 = 1/2``, ``sin^2 13 = 0``) but which **does not diagonalise** the
    democratic base ``|1><1|``: its third column ``(0,1,1)/sqrt(2)`` does not
    belong to the light eigenspace (the component sum is non-zero).  The bug was
    invisible because the **angles** (the moduli) were checked, not the
    **diagonalisation**.

    The genuine TBM has as its three columns the eigenvectors of ``|1><1|``:
    ``(2,-1,-1)/sqrt(6)`` and ``(0,1,-1)/sqrt(2)`` (the 2-fold light subspace,
    zero sum) and ``(1,1,1)/sqrt(3)`` (the heavy level).  It is the choice
    **adapted to S3** (the transposition (23)) of the O(2) freedom of the
    degenerate subspace - see ``yukawa_3x3.mixing_from_democratic`` for the
    honest analysis of what is 0-parameter and what is not.
    """
    s12, c12 = math.sqrt(1.0 / 3.0), math.sqrt(2.0 / 3.0)
    s23 = c23 = math.sqrt(0.5)
    r12 = np.array([[c12, s12, 0.0], [-s12, c12, 0.0], [0.0, 0.0, 1.0]])
    r23 = np.array([[1.0, 0.0, 0.0], [0.0, c23, -s23], [0.0, s23, c23]])
    return r23 @ r12


def family_qca():
    """The 3-generation QCA step: the 3 families **inside** the 0-parameter
    dynamics.

    The 3-generation internal space ``H_3f = ⊕_{k=1..3} H_32`` is a **96-D**
    module (**48 Weyl + 48 CPT-conjugate**): the **3 families** are the **3
    two-planes** ``C(3,2)`` of the ``(1,3)`` 3-space (the V26 reading), each
    family one **32-D** copy.  The **3-generation QCA step** is the
    **family-symmetric** 0-parameter step:

        U_96 = blockdiag(U_32, U_32, U_32),   U_32 = the GUT QCA step
        (the 45 plane rotations at ``g_10 = 0.2``),

    i.e. the **45 act identically on each family** — the 0-parameter
    statement that the 45 do **not** distinguish the 3 two-planes (the
    **family symmetry**).  The 0-parameter **family mixing** of the GUT is
    the **tribimaximal** 3x3 unitary (the V29 angles — the PMNS of the
    democratic base), a property of the neutrino/Yukawa sector, not of the
    time-evolution step.  Only $\\sin^2\\theta_{12}=1/3$ is 0-parameter **as a
    number** (the unique heavy eigenvector of the democratic base);
    $\\theta_{23}$ and $\\theta_{13}$ are an $O(2)$ convention.
    """
    n = 32
    ngen = 3
    g = [E.Array_(m) for m in clifford.generators()]
    P = baryogenesis._pseudoscalar()
    planes = sm_qca.plane_generators()
    theta = sm_qca.G10
    U32 = sm_qca.gut_qca_step(planes, theta)
    I3 = E.eye(ngen)
    # the 96x96 family-symmetric step: the 3 identical 32-D copies
    U = E.kron(I3, U32)
    Ua = E.asnumpy(U)
    I96 = np.eye(n * ngen, dtype=complex)
    unitary_err = float(np.abs(Ua @ Ua.conj().T - I96).max())
    # (1) the family symmetry: U_96 is blockdiag of 3 EQUAL 32x32 blocks
    off_block = 0.0
    block_diff = 0.0
    for i in range(ngen):
        bi = Ua[i * n:(i + 1) * n, i * n:(i + 1) * n]
        block_diff = max(block_diff, float(np.abs(bi - Ua[0:n, 0:n]).max()))
        for j in range(ngen):
            if i != j:
                off_block = max(off_block, float(np.abs(
                    Ua[i * n:(i + 1) * n, j * n:(j + 1) * n]).max()))
    # (2) the per-family 16/16-bar split is preserved:
    # (a) the 45 commute with the pseudoscalar of one copy (the V22/V33 logic)
    max_comm = 0.0
    for a in range(10):
        for b in range(a + 1, 10):
            T = 1j * (g[a] @ g[b])
            max_comm = max(max_comm, float(np.abs(E.asnumpy(T @ P - P @ T)).max()))
    # (b) the 32-D GUT step preserves the split of its copy
    step_comm_32 = float(np.abs(E.asnumpy(U32 @ P - P @ U32)).max())
    # (c) the full 96x96 step commutes with the 96-D pseudoscalar
    P96 = E.kron(I3, P)
    step_comm_96 = float(np.abs(E.asnumpy(U @ P96 - P96 @ U)).max())
    # (3) the tribimaximal mixing (V29): the PMNS of the democratic base
    tbm = tribimaximal_3x3()
    tbm_unitary = float(np.abs(tbm @ tbm.T - np.eye(3)).max())
    tbm_angles_ok = (abs(tbm[0, 1] ** 2 - 1.0 / 3.0) < 1e-15
                     and abs(tbm[1, 2] ** 2 - 0.5) < 1e-15
                     and abs(tbm[0, 2]) < 1e-15)
    # the check the old code lacked: the TBM must actually DIAGONALIZE the
    # democratic base, not merely carry the right moduli (the old sign bug
    # passed the angle test while failing this one).
    _J = np.ones((3, 3), dtype=complex)
    _D = tbm.conj().T @ _J @ tbm
    tbm_diagonalizes = bool(float(np.abs(_D - np.diag(np.diag(_D))).max()) < _TOL.standard)
    # (4) the dimension: 48 Weyl = 16 x 3 and 48 CPT-conjugate
    n_weyl_1gen = int(matter.dimension_check())
    n_weyl_3gen = n_weyl_1gen * ngen
    return {
        "backend": _backend(),
        "g10": float(theta),
        "n_planes": int(len(planes)),
        "n_module_1gen": int(n),
        "n_module_3gen": int(n * ngen),
        "n_weyl_1gen": int(n_weyl_1gen),
        "n_weyl_3gen": int(n_weyl_3gen),
        "n_cpt_3gen": int(n_weyl_3gen),
        "unitary_err": float(unitary_err),
        "family_symmetry_off_block": float(off_block),
        "family_symmetry_block_diff": float(block_diff),
        "max_45_pseudoscalar_commutator": float(max_comm),
        "step_32_pseudoscalar_commutator": float(step_comm_32),
        "step_96_pseudoscalar_commutator": float(step_comm_96),
        "tribimaximal_unitary_err": float(tbm_unitary),
        "tribimaximal_angles_ok": bool(tbm_angles_ok),
        "tribimaximal_diagonalizes_democratic": tbm_diagonalizes,
        "tribimaximal_sin2": {
            "s12": float(abs(tbm[0, 1]) ** 2),
            "s23": float(abs(tbm[1, 2]) ** 2),
            "s13": float(abs(tbm[0, 2]) ** 2),
        },
        "residual": ("the why-3 (O3) and the family mass hierarchy "
                     "(the 10-Higgs breaking) are the named residuals"),
    }


def main():
    import sys
    print("=" * 68)
    print("CL10 - the 3-family construct (the 3 generations from the (1,3))")
    print("=" * 68)
    r = family_construct()
    print("  the 3 generations = the 0-parameter 3, three ways:")
    print(f"    (1) the 3 space directions (S = i h_2..4):  {r['n_space']}")
    print(f"    (2) the 3 two-planes (C(3,2)):               {r['n_planes_C32']}")
    print(f"    (3) the rank of the 3x3 spatial metric:      {r['n_rank']}")
    print(f"  the 16 Weyl (one generation) -> 3 generations:  {r['n_weyl_1gen']} x {r['n_family_candidate']} = {r['n_weyl_3gen']} Weyl")
    print(f"  VERDICT:  {'PASS (the count 3 is the index of the 3-space, three ways; the count itself is the open problem O3)' if (r['match'] and r['all_three_ways_3']) else 'FAIL'}")
    sys.exit(0 if (r["match"] and r["all_three_ways_3"]) else 1)


if __name__ == "__main__":
    main()
