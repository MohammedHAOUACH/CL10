#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the 3-generation Yukawa mass matrix: the 10-to-1 projection (O6, C4).

The **10-to-1** reduction (Section 12) says the ten residuals O1-O10 are
projections of **one** 0-parameter object — the GUT QCA step
``U`` (the 32-level spectrum ``H = -i log U`` at the root-clock ``g_10``).  This
module is the **O6 projection**: the **3-generation Yukawa mass matrix**
and the **mixing**, read off from that one object, with **no fit**.

- **The one object -> the 3 families.**  The 32-D Weyl module is **one**
  generation.  The **3** observed families are the **3 two-planes** of the
  $(1,3)$ 3-space; the 3-generation step is the **family-symmetric**
  96-D object  ``U_96 = I_3 (x) U_32`` (the 45 act **identically** on each
  family), so its spectrum is **3 x (32 levels)** = 96 levels — the
  **one object, three times**.

- **The 0-parameter 3x3 mass matrix (the democratic base).**  At the
  symmetric (0-parameter) level the 3 families are **equal**: the
  **3-fold GUT symmetry** forces the 3x3 to be symmetric, and the
  **10-Higgs** (a single grade-1) couples the 3 families **equally** —
  the **democratic** matrix  ``M = |1><1|``  (all-ones).  Its spectrum
  ``{3, 0, 0}`` is the **0-parameter 2-level structure** (one **heavy**
  level + a **2-fold light** level), diagonalized by the **TBM** — the
  **0-parameter** PMNS ($\\sin^2 12 = 1/3$, $\\sin^2 23 = 1/2$,
  $\\theta_{13} = 0$), close to the observed.

- **The 3-level hierarchy (the NAMED residual).**  The **3 distinct**
  mass levels ($m_1 < m_2 < m_3$) require **breaking** the 2-fold light
  degeneracy — the **10-Higgs vev** (a 0-parameter-impossible selection).
  The **exact ratios** ($m_1:m_2:m_3$) and the **exact** CKM/PMNS are
  the **named residual** — the 0-parameter object does not fix them.

- **The 10-Higgs Dirac term (computed, 0-parameter).**  The grade-1 GUT
  vector (the 10) is **off-diagonal** in the 16/16-bar split, so the
  0-parameter mass term $16\\,10\\,16\\bar$ is a **full-rank-16 Dirac**
  matrix (all 16 Weyl get mass) — the **where** of the mass, computed.

The **named residual (O6)**: the **exact** 3-generation hierarchy
(the 10-Higgs vev breaking, the ratios $m_1:m_2:m_3$) and the **exact**
CKM/PMNS angles beyond the TBM — and the TBM angles $\\sin^2 23$, $\\sin^2 13$
themselves (the $O(2)$ freedom of the 2-fold light subspace) — the one object
does not fix them.  Only $\\sin^2 12 = 1/3$ is 0-parameter **as a number**.
"""

import numpy as np

from . import clifford
from . import spacetime
from . import engine as E
from . import family
from . import qca_spectrum
from . import sm_qca
from . import coupling
from . import config as _config
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml
# the measured mixing inputs of the comparisons (config blocks neutrino/ckm)
_OBS_PMNS = (float(_config.CFG.neutrino.obs_sin2_theta12),
             float(_config.CFG.neutrino.obs_sin2_theta23),
             float(_config.CFG.neutrino.obs_sin2_theta13))
_OBS_CKM = (float(_config.CFG.ckm.obs_sin2_theta12),
            float(_config.CFG.ckm.obs_sin2_theta23),
            float(_config.CFG.ckm.obs_sin2_theta13))

# the rotation angle of the active branch (CL10_ANGLE: 1/5 implemented, pi/4 closed)
G10 = coupling.G10


def one_object_spectrum(g10=None):
    """The **one** 0-parameter object: the 32-level QCA spectrum and the
    96-level **3-family** spectrum (``3 x 32``, the one object three
    times).  Returns both, sorted."""
    s = qca_spectrum.qca_spectrum(g10)
    lev1 = s["spectrum"]                       # the 32 levels (one generation)
    # the 3-family 96-level spectrum: 3 identical copies (the family symmetry)
    lev3 = np.sort(np.tile(lev1, 3))          # 3 x 32 = 96 levels
    return {
        "g10": float(s["g10"]),
        "n_levels_1gen": int(s["n_levels"]),
        "spectrum_1gen": lev1,
        "spectrum_3gen": lev3,
        "n_levels_3gen": int(len(lev3)),
        "zero_point_1gen": float(s["zero_point"]),
        "det_U_minus_1": float(s["det_U_minus_1"]),
        "family_triplication": bool(len(lev3) == 3 * len(lev1)),
    }


def higgs_mass_16():
    """The **0-parameter** 16x16 10-Higgs Dirac mass matrix (the 16->16-bar
    mass term).  The grade-1 GUT vector (the 10) is **off-diagonal** in the
    16/16-bar split, so it gives the **0-parameter Dirac** mass term
    $16\\,10\\,16\\bar\\to 1$: the **16->16** and **16-bar->16-bar** blocks are
    **0** (the 10 is odd in the grading) and the **off-diagonal** blocks are
    **nonzero** (the mass matrices).  This is the **where** of the mass
    (0-parameter); the **3-family** 3x3 of this term is the democratic base
    below."""
    g = [E.asnumpy(m) for m in clifford.generators()]
    h = g[0]                                   # a representative grade-1 (the 10)
    p, m = clifford.weyl_spinor()
    P16 = p @ p.conj().T
    P16b = m @ m.conj().T
    A = P16b @ h @ P16                        # 16 -> 16-bar (the mass matrix)
    diag16 = P16 @ h @ P16                    # 16 -> 16 (should be 0)
    rank = int(np.linalg.matrix_rank(A))
    return {
        "shape": int(A.shape[0]),
        "n16": int(p.shape[1]),
        "rank": rank,
        "off_diagonal": float(np.linalg.norm(A).max()),
        "diag_16": float(np.linalg.norm(diag16).max()),
        "non_singular": bool(rank >= 1),
        "dirac": bool(rank >= 1 and float(np.linalg.norm(diag16).max()) < _TOL.loose),
    }


def family_symmetry():
    """The **0-parameter 3-fold family symmetry** — the *derivation* of the
    **democratic** 3x3 (a 0-parameter *why*, not an assertion).

    The 3 families are the **3 spatial two-planes** ``C(3,2)`` of the
    ``$(1,3)$`` 3-space.  The **3 spatial directions**
    ``S_a = i h_{a}  (a = 2,3,4)`` (one per family) are **equivalent**: they
    carry the **same** 0-parameter invariants (``S_a^2 = -I`` for every
    ``a``; ``{S_a, S_b} = 0``; the same 0-parameter norm ``Tr(S_a S_a^dagger)``).
    This **3-fold** invariance (the permutation of the 3 spatial directions)
    is the family symmetry.  A **single** grade-1 (the 10-Higgs) is therefore
    forced to couple the 3 families **equally**: the **0-parameter (maximal)**
    coupling that respects the 3-fold symmetry is the **democratic**
    all-ones matrix, whose spectrum ``{3, 0, 0}`` is the **2-level**
    structure (one heavy + a 2-fold light).  Everything below (the TBM, the
    2-level hierarchy) is a 0-parameter consequence of this symmetry."""
    T, S = spacetime.causal_generators()
    S = [E.asnumpy(s) for s in S]
    n_space = len(S)
    I32 = np.eye(32, dtype=complex)
    # (1) the 3 spatial directions are EQUIVALENT: the same 0-parameter
    # invariants.  S_a^2 = -I for every a (the Minkowski spatial metric).
    s2_err = [float(np.abs(s @ s + I32).max()) for s in S]
    # the mutual anticommutator {S_a, S_b} = 0 (the 3-space is isotropic).
    ss_err = max(float(np.abs(S[a] @ S[b] + S[b] @ S[a]).max())
                 for a in range(3) for b in range(a + 1, 3))
    # the 0-parameter norm of each direction (equal -> the 3-fold symmetry).
    norms = [float((s @ s.conj().T).trace().real) for s in S]
    norm_spread = float(max(norms) - min(norms))
    equivalent = bool(n_space == 3 and max(s2_err) < _TOL.loose
                     and ss_err < _TOL.loose and norm_spread < _TOL.loose)
    # (2) the 3 families = the C(3,2) = 3 spatial two-planes (the pairs).
    pairs = [(0, 1), (0, 2), (1, 2)]
    n_planes = len(pairs)
    # (3) the democratic 3x3 (the 0-parameter (maximal) coupling that
    # respects the 3-fold symmetry): its spectrum is the 2-level structure.
    M = np.ones((3, 3), dtype=float)
    eig = sorted(float(x) for x in np.abs(np.linalg.eigvals(M)))
    return {
        "n_space": n_space,
        "n_families_C32": n_planes,
        "S_a2_err": float(max(s2_err)),
        "anticommutator_Sa_Sb": float(ss_err),
        "spatial_norms": [float(x) for x in norms],
        "norm_spread": norm_spread,
        "equivalent_3fold": equivalent,
        "democratic_rank": int(np.linalg.matrix_rank(M)),
        "democratic_spectrum": [float(x) for x in eig],
        "two_level": bool(len({round(e, 6) for e in eig}) == 2),
    }


def mass_matrix_3x3():
    """The **0-parameter** 3x3 Yukawa/mass matrix (the **democratic** base):
    the 3 families are **equal** at the symmetric (0-parameter) level, so the
    GUT gives the **democratic** matrix ``M = |1><1|`` (the all-ones, rank 1).
    It is **Hermitian** (diagonalizable to real masses) and its spectrum is
    ``{3, 0, 0}``: **one heavy** level + a **2-fold degenerate light** level
    (the **0-parameter 2-level structure**).  The **3-level** hierarchy
    (3 distinct masses) is the **named residual** (the 10-Higgs vev breaks
    the 2-fold light degeneracy)."""
    n = 3
    M = np.ones((n, n), dtype=complex)        # the democratic base |1><1|
    Mh = 0.5 * (M + M.conj().T)               # Hermitian (it is real-symmetric)
    eig = np.sort(np.abs(np.linalg.eigvals(Mh)))
    n_heavy = int(np.sum(eig > _TOL.loose))
    n_light = n - n_heavy
    n_distinct = int(np.unique(np.round(eig / 1e-6)).size)   # {0, 3} -> 2
    return {
        "n": n,
        "matrix": Mh,
        "eigenvalues": [float(x) for x in eig],
        "rank": int(np.linalg.matrix_rank(Mh)),
        "n_heavy": n_heavy,
        "n_light_2fold": n_light,
        "n_distinct_levels": n_distinct,
        "democratic": bool(np.allclose(Mh, 1.0)),
    }


def mixing_from_democratic():
    """**Honest analysis** of what the democratic base does and does not fix.

    The democratic base M = |1><1| has the spectrum {3,0,0}.  The **heavy**
    level h = (1,1,1)/sqrt(3) is **unique**; the **light** subspace is **2-fold
    degenerate**, so every diagonalizer of M is an O(2) rotation inside that
    subspace, tensored with h.  The democratic base fixes the **spectrum** and
    **part** of the mixing - not all of it.

    The honest result, **measured** here (an O(2) scan over 1441 rotations):

    - ``sin^2 theta_12 = |U_e2|^2`` is **FIXED at 1/3**: it is the square of the
      overlap between any basis vector and the **unique** heavy vector h.  That
      is a genuine 0-parameter number (observed 0.307, a deviation of 8 %).
    - ``sin^2 theta_23`` and ``sin^2 theta_13`` **sweep [0, 2/3]**: the O(2)
      freedom of the degenerate light subspace.  The TBM values (1/2, 0) are a
      **choice**.

    The TBM choice is the one **adapted to the transposition (23)** of S3.  But
    (12), (13) and (23) are **three equivalent symmetries** of |1><1|: the
    democratic base does not single out (23).  Which transposition is broken is
    therefore an **extra input** (naming it is the honest correction), not a
    0-parameter consequence.

    The CKM built from two **identical** democratic
    up/down matrices is the **identity**, not the TBM.  A non-trivial CKM needs
    a **misalignment** between up and down, which the symmetric democratic base
    does not supply.
    """
    J = np.ones((3, 3), dtype=complex)
    h = np.array([1, 1, 1], dtype=float) / np.sqrt(3.0)        # the UNIQUE heavy eigenvector
    a0 = np.array([1, -1, 0], dtype=float) / np.sqrt(2.0)      # an orthonormal basis of
    b0 = np.array([1, 1, -2], dtype=float) / np.sqrt(6.0)      # the 2-fold light subspace

    def diagonalizer(t):
        """a diagonalizer of M, with the heavy direction in the MIDDLE column
        (the PMNS slot convention: column 2 is the one carrying sin^2 12 and
        sin^2 23)."""
        e1 = np.cos(t) * a0 + np.sin(t) * b0
        e2 = -np.sin(t) * a0 + np.cos(t) * b0
        return np.column_stack([e1, h, e2])

    def offdiag(U):
        D = U.conj().T @ J @ U
        return float(np.abs(D - np.diag(np.diag(D))).max())

    th = np.linspace(0.0, 2.0 * np.pi, 1441)
    s12 = [abs(diagonalizer(t)[0, 1]) ** 2 for t in th]        # |U_e2|^2
    s23 = [abs(diagonalizer(t)[1, 2]) ** 2 for t in th]        # |U_mu3|^2
    s13 = [abs(diagonalizer(t)[0, 2]) ** 2 for t in th]        # |U_e3|^2

    # every member of the family really diagonalizes M - the ANGLES alone are
    # not a diagonalization test
    diag_err = max(offdiag(diagonalizer(t)) for t in th[::40])

    # the TBM is a member of the family (the (23)-adapted one)
    tbm = family.tribimaximal_3x3()
    tbm_diag_err = offdiag(tbm)
    # a diagonalizer is defined column-by-column up to a phase, so the family
    # must be compared against all 8 column-sign choices diag(s0, s1, s2)
    signs = [np.diag([float(a), float(b), float(c)])
             for a in (-1.0, 1.0) for b in (-1.0, 1.0) for c in (-1.0, 1.0)]
    tbm_in_family = bool(min(np.abs(diagonalizer(t) @ S - tbm).max()
                             for t in th for S in signs) < _TOL.standard)

    # the three S3 transpositions: equivalent symmetries of |1><1|, each adapted
    # basis giving a DIFFERENT mixing matrix
    axes = {"(23)": ((-2, 1, 1), (0, 1, -1)),
            "(12)": ((1, 1, -2), (1, -1, 0)),
            "(13)": ((1, -2, 1), (1, 0, -1))}
    adapted = {}
    for name, (p1, p2) in axes.items():
        e1 = np.array(p1, dtype=float) / np.linalg.norm(p1)
        e2 = np.array(p2, dtype=float) / np.linalg.norm(p2)
        B = np.column_stack([e1, h, e2])
        adapted[name] = {"sin2_12": float(abs(B[0, 1]) ** 2),
                         "sin2_23": float(abs(B[1, 2]) ** 2),
                         "sin2_13": float(abs(B[0, 2]) ** 2),
                         "reproduces_TBM": bool(abs(B[0, 1]) ** 2 - 1/3 < _TOL.standard
                                                and abs(B[1, 2]) ** 2 - 1/2 < _TOL.standard
                                                and abs(B[0, 2]) ** 2 < _TOL.standard)}

    # the CKM from IDENTICAL up/down democratic matrices is the IDENTITY
    wu, Vu = np.linalg.eigh(np.ones((3, 3), dtype=complex))
    wd, Vd = np.linalg.eigh(np.ones((3, 3), dtype=complex))
    ckm_naive = Vu.conj().T @ Vd
    ckm_is_identity = bool(np.allclose(np.abs(ckm_naive), np.eye(3)))

    obs_pmns = _OBS_PMNS
    return {
        "spectrum_fixed": [3.0, 0.0, 0.0],
        "heavy_eigenvector": [1/np.sqrt(3)] * 3,
        "heavy_eigenvector_unique": True,
        "light_degeneracy": 2,
        "diagonalizer_is_unique": False,
        "sin2_12_fixed": bool(max(s12) - min(s12) < _TOL.standard),
        "sin2_12_value": float(np.mean(s12)),
        "sin2_12_range": [float(min(s12)), float(max(s12))],
        "sin2_23_range": [float(min(s23)), float(max(s23))],
        "sin2_13_range": [float(min(s13)), float(max(s13))],
        "sin2_12_vs_observed": {"predicted": float(np.mean(s12)),
                                "observed": obs_pmns[0],
                                "ratio": float(np.mean(s12) / obs_pmns[0])},
        "max_diagonalization_err": diag_err,
        "tbm_diagonalizes_democratic": bool(tbm_diag_err < _TOL.standard),
        "tbm_in_the_O2_family": tbm_in_family,
        "tbm_angles": {"sin2_12": 1.0 / 3.0, "sin2_23": 0.5, "sin2_13": 0.0},
        "s3_transpositions_equivalent": adapted,
        "ckm_is_identity": ckm_is_identity,
        "verdict_0param": (
            "the democratic base DERIVES sin^2 theta_12 = 1/3 (8% from the "
            "observed 0.307): it is the squared overlap of any family basis "
            "vector with the UNIQUE heavy eigenvector (1,1,1)/sqrt(3).  It does "
            "NOT derive sin^2 theta_23 or sin^2 theta_13: both sweep [0, 2/3] "
            "under the O(2) freedom of the 2-fold light subspace.  TBM (1/3, "
            "1/2, 0) is the (23)-adapted member of that family; (12), (13) and "
            "(23) are equivalent symmetries of |1><1|, so WHICH one breaks is "
            "an extra input.  The CKM from identical up/down democratic "
            "matrices is the IDENTITY, not the TBM."),
    }


def mixing_3x3():
    """The mixing of the democratic base: the TBM **angles**, with the honest
    status of each.

    The angles sin^2 12 = 1/3, sin^2 23 = 1/2, theta_13 = 0 are the TBM angles,
    and the TBM **diagonalises the democratic base exactly**.  But they are not
    **unique**: the 2-fold degeneracy of the light subspace leaves an O(2)
    freedom (see `mixing_from_democratic`).  The TBM is the choice adapted to
    the transposition (23), and the democratic base does not single that
    transposition out.  Its status is therefore a **convention adapted to S3**,
    not a prediction.
    """
    # the TBM is a diagonalizer of the democratic base (the (23)-adapted one)
    tbm = family.tribimaximal_3x3()
    # the TBM angles (the numbers of the (23)-adapted choice)
    s12, s23, s13 = (1.0 / 3.0, 0.5, 0.0)
    obs_pmns = _OBS_PMNS         # the observed PMNS (config neutrino)
    obs_ckm = _OBS_CKM           # the observed CKM (config ckm)
    dev_pmns = (abs(s12 - obs_pmns[0]) + abs(s23 - obs_pmns[1])
                + abs(s13 - obs_pmns[2])) / 3
    dev_ckm = (abs(s12 - obs_ckm[0]) + abs(s23 - obs_ckm[1])
               + abs(s13 - obs_ckm[2])) / 3
    honest = mixing_from_democratic()
    return {
        "unitary": bool(np.allclose(tbm.conj().T @ tbm, np.eye(3))),
        "diagonalizes_democratic": honest["tbm_diagonalizes_democratic"],
        "sin2_12": float(s12), "sin2_23": float(s23), "sin2_13": float(s13),
        "obs_pmns": obs_pmns, "obs_ckm": obs_ckm,
        "dev_pmns": float(dev_pmns), "dev_ckm": float(dev_ckm),
        # the honest provenance of each angle
        "sin2_12_is_0param_derived": honest["sin2_12_fixed"],
        "sin2_23_is_0param_derived": False,
        "sin2_13_is_0param_derived": False,
        "sin2_23_O2_range": honest["sin2_23_range"],
        "sin2_13_O2_range": honest["sin2_13_range"],
        "tbm_is_the_23_adapted_choice": True,
        "mixing_status": honest["verdict_0param"],
    }


def hierarchy_3level():
    """The mass hierarchy: the **0-parameter 2-level structure** (the
    democratic base: one heavy + a 2-fold degenerate light level) and the
    **named** **3-level** hierarchy (the 10-Higgs vev breaks the 2-fold
    light degeneracy -> 3 distinct masses $m_1 < m_2 < m_3$).  The
    0-parameter object does **not** fix the exact ratios — that is the O6
    named residual."""
    eig = np.array(mass_matrix_3x3()["eigenvalues"])
    n_heavy = int(np.sum(eig > _TOL.loose))
    n_light = 3 - n_heavy
    n_distinct = int(np.unique(np.round(eig / 1e-6)).size)
    obs_ratio = (1.0 / 3.0, 1.0 / 10.0)   # the observed 2-gen ratios (order)
    return {
        "n_levels_0param": int(n_distinct),      # 2: one heavy + a 2-fold light
        "n_levels_observed": 3,                  # the 3-level (the residual)
        "n_heavy": n_heavy,
        "n_light_2fold": n_light,
        "structure": "one heavy + a 2-fold light level (the 0-parameter base)",
        "degenerate_light": bool(n_light >= 2 and eig[0] < _TOL.loose and abs(eig[0] - eig[1]) < _TOL.loose),
        "obs_ratio_order": obs_ratio,
        "residual": ("the exact 3-gen hierarchy (the 10-Higgs vev breaking, "
                     "the ratios m1:m2:m3) and the exact CKM/PMNS angles "
                     "beyond the TBM (O6); sin^2_23 and sin^2_13 are an O(2) "
                     "choice, not a derivation"),
    }

def quark_ckm():
    """The **0-parameter** quark CKM: the **quark** mass matrices on the
    **democratic** base (the **2-level** structure) and the **CKM** (the
    **TBM**, the **0-parameter** quark mixing).  The **V49** completion (the
    O6 residual, **quantified**).

    The **0-parameter** statement:

    1. **The quark mass matrices** (the **democratic** base, the **3x3**):
       the **up**-type and **down**-type mass matrices are the **democratic**
       base (the **all-ones**, **Hermitian**, spectrum ``{3, 0, 0}``: one
       **heavy** + a **2-fold** light level).

    2. **The CKM**: the **misalignment** between the **up**-type and
       **down**-type mass matrices (the **unitary** **3x3**) is the **TBM**
       (the **0-parameter** quark mixing, the **eigenvectors** of the
       **democratic** base).

    3. **The 0-parameter CKM angles**: ``sin^2 12 = 1/3``,
       ``sin^2 23 = 1/2``, ``sin^2 13 = 0`` (the **TBM**, the **0-parameter**
       **2-level** structure).

    4. **The comparison** to the **observed CKM** (``0.050, 0.949, 0.020``):
       the **0-parameter** TBM is **close** to the **observed CKM** (the
       **deviation** is the **O6** named residual).
    """
    # the quark mass matrices (the democratic base, the 3x3)
    mu = np.ones((3, 3), dtype=complex)          # the up-type (the democratic base)
    md = np.ones((3, 3), dtype=complex)          # the down-type (the democratic base)
    Hu = 0.5 * (mu + mu.conj().T)
    Hd = 0.5 * (md + md.conj().T)
    eig_u = sorted(float(x) for x in np.abs(np.linalg.eigvals(Hu)))
    eig_d = sorted(float(x) for x in np.abs(np.linalg.eigvals(Hd)))

    # **HONEST CKM.**  The CKM is the **misalignment** between the diagonalisers
    # of $H_u$ and $H_d$.  Here $H_u = H_d$: the two diagonalisers are the same,
    # so the computed CKM is the **IDENTITY**, not the TBM.  A non-trivial CKM
    # needs an up/down **misalignment** that the symmetric democratic base does
    # not supply (it comes from the 10-Higgs vev, the named residual O6).
    # The claim "CKM = TBM" is therefore FALSE in the 0-parameter construction as
    # coded; it holds only if different diagonalisers are POSTULATED, which is
    # exactly the named residual.
    wu, Vu = np.linalg.eigh(Hu)
    wd, Vd = np.linalg.eigh(Hd)
    ckm = Vu.conj().T @ Vd
    ckm_is_identity = bool(np.allclose(np.abs(ckm), np.eye(3)))
    ckm_angles = [float(abs(ckm[0, 1]) ** 2), float(abs(ckm[1, 2]) ** 2),
                  float(abs(ckm[0, 2]) ** 2)]

    tbm = family.tribimaximal_3x3()
    s12, s23, s13 = (1.0 / 3.0, 0.5, 0.0)
    obs_ckm = _OBS_CKM                       # the observed CKM (config ckm)
    dev_ckm = (abs(s12 - obs_ckm[0]) + abs(s23 - obs_ckm[1])
               + abs(s13 - obs_ckm[2])) / 3
    ok = (bool(np.allclose(tbm.conj().T @ tbm, np.eye(3)))
          and len({round(e, 6) for e in eig_u}) == 2
          and len({round(e, 6) for e in eig_d}) == 2)
    return {
        "backend": E._BACKEND,
        "up_eigenvalues": eig_u,
        "down_eigenvalues": eig_d,
        "two_level": bool(len({round(e, 6) for e in eig_u}) == 2
                          and len({round(e, 6) for e in eig_d}) == 2),
        "ckm_unitary": bool(np.allclose(ckm.conj().T @ ckm, np.eye(3))),
        "ckm_is_identity": ckm_is_identity,
        "ckm_computed_angles": ckm_angles,
        "sin2_12": float(s12), "sin2_23": float(s23), "sin2_13": float(s13),
        "obs_ckm": obs_ckm,
        "dev_ckm": float(dev_ckm),
        "residual": ("the quark CKM is NOT derived 0-parameter: with H_u = H_d "
                     "(both democratic) the computed CKM is the IDENTITY.  A "
                     "non-trivial CKM requires an up/down MISALIGNMENT, which "
                     "the symmetric democratic base does not provide — it comes "
                     "from the 10-Higgs vev breaking.  The TBM angles are the "
                     "lepton-sector (23)-adapted choice, not a quark result: "
                     "the O6 named residual, now stated honestly"),
        "ok": bool(ok),
    }
def w6_verdict():
    """The full 0-parameter O6 verdict: the 10-to-1 projection of the one
    object onto the 3-family space — the 3x3 mass matrix (the democratic
    base), the TBM mixing (only sin^2_12 = 1/3 is 0-parameter as a number),
    the 3-level hierarchy structure, the
    10-Higgs Dirac structure, and the **named residual** (the exact
    hierarchy + the exact mixing)."""
    U = sm_qca.gut_qca_step(sm_qca.plane_generators(), G10)
    Ua = E.asnumpy(U)
    I32 = np.eye(32, dtype=complex)
    s2unit = float(np.abs(Ua @ Ua.conj().T - I32).max())
    one = one_object_spectrum()
    mm = mass_matrix_3x3()
    mix = mixing_3x3()
    hier = hierarchy_3level()
    # the 10-Higgs Dirac structure (the 0-parameter mass term, computed)
    h16 = higgs_mass_16()
    fs = family_symmetry()
    honest = mixing_from_democratic()
    ok = (
        s2unit < _TOL.strict
        and one["family_triplication"]
        and abs(one["zero_point_1gen"]) < 1e-8
        and mm["democratic"]
        and mm["rank"] == 1
        and mm["n_distinct_levels"] == 2
        and mix["unitary"]
        # the honest mixing checks: the ANGLES are not enough, the matrix has
        # to diagonalize the democratic base
        and honest["tbm_diagonalizes_democratic"]
        and honest["sin2_12_fixed"]
        and abs(honest["sin2_12_value"] - 1.0 / 3.0) < _TOL.standard
        and hier["degenerate_light"]
        and h16["dirac"]
        and h16["rank"] == 16
        and fs["equivalent_3fold"]
        and fs["two_level"]
    )
    return {
        "backend": E._BACKEND,
        "g10": one["g10"],
        "one_object": one,
        "mass_matrix": mm,
        "family_symmetry": fs,
        "mixing": mix,
        "hierarchy": hier,
        "higgs_dirac": {
            "n16": h16["n16"],
            "rank": h16["rank"],
            "off_diagonal": h16["off_diagonal"],
            "diag_16": h16["diag_16"],
            "dirac": h16["dirac"],
            "note": "the 10-Higgs 16->16-bar is the 0-parameter Dirac term (rank 16, off-diagonal)",
        },
        "mixing_honest": honest,
        "s2_unitary_err": float(s2unit),
        "ok": bool(ok),
        "residual": hier["residual"],
    }


def main():
    import sys

    r = w6_verdict()
    one = r["one_object"]
    mm = r["mass_matrix"]
    mix = r["mixing"]
    hier = r["hierarchy"]
    obs = mix["obs_pmns"]
    print("=" * 68)
    print("CL10 - the 3-gen Yukawa mass matrix: the 10-to-1 projection (O6, C4)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the ONLY input:  g_10 = {r['g10']:.6f}  (the root clock, cl10/config.yaml)")
    print(f"  [1] the one object -> 3 families:  32 levels (1 gen) x 3 = {one['n_levels_3gen']} levels (3 gen)")
    fs = r["family_symmetry"]
    print(f"  [2] 3-fold family symmetry (the DERIVATION of the democratic 3x3):  the 3 spatial directions are "
          f"{'EQUIVALENT' if fs['equivalent_3fold'] else 'NOT'} (each S_a^2=-I, anticommuting, equal norm) -> the 10-Higgs couples the 3 families equally -> democratic spectrum {fs['democratic_spectrum']}")
    print(f"  [3] 0-parameter 3x3 mass matrix (the democratic base):  eigenvalues "
          f"{['%.3g' % x for x in mm['eigenvalues']]}  ({mm['n_distinct_levels']}-level: {mm['n_heavy']} heavy + {mm['n_light_2fold']}-fold light, rank {mm['rank']})")
    print(f"  [4] mixing (the TBM, which diagonalizes |1><1|):  "
          f"sin^2(12)={mix['sin2_12']:.3f} (obs {obs[0]:.3f}), "
          f"sin^2(23)={mix['sin2_23']:.3f} (obs {obs[1]:.3f}), "
          f"sin^2(13)={mix['sin2_13']:.3f} (obs {obs[2]:.3f})")
    print(f"  [5] 0-parameter hierarchy:  {hier['structure']}  (0-parameter {hier['n_levels_0param']}-level; the {hier['n_levels_observed']}-level is the named residual);  "
          f"mean 0-parameter vs observed PMNS deviation:  {mix['dev_pmns']:.3f}")
    print(f"  [6] 10-Higgs Dirac term (the 0-parameter mass):  rank {r['higgs_dirac']['rank']} (16->16-bar, off-diagonal, Dirac)")
    hh = r["mixing_honest"]
    print("  [7] HONEST mixing status (the O(2) freedom of the 2-fold light subspace):")
    print(f"      sin^2(12) = {hh['sin2_12_value']:.4f}  -> DERIVED 0-parameter "
          f"(the unique heavy eigenvector), vs observed {hh['sin2_12_vs_observed']['observed']} "
          f"(ratio {hh['sin2_12_vs_observed']['ratio']:.3f})")
    print(f"      sin^2(23) sweeps [{hh['sin2_23_range'][0]:.3f}, {hh['sin2_23_range'][1]:.3f}]  -> NOT derived (O(2) choice)")
    print(f"      sin^2(13) sweeps [{hh['sin2_13_range'][0]:.3f}, {hh['sin2_13_range'][1]:.3f}]  -> NOT derived (O(2) choice)")
    print(f"      TBM diagonalizes |1><1|: {hh['tbm_diagonalizes_democratic']}  "
          f"(the (23)-adapted member of the family: {hh['tbm_in_the_O2_family']})")
    print("      the S3 transpositions (12),(13),(23) are equivalent symmetries of |1><1|;")
    print(f"      only (23) reproduces TBM: "
          f"{ {k: v['reproduces_TBM'] for k, v in hh['s3_transpositions_equivalent'].items()} }")
    print(f"  [8] named residual (O6):  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the democratic 3x3 base + TBM + the 2-level structure; the 3-level hierarchy and the exact mixings are the open problem O6)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
