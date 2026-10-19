"""The engine's invariants as an executable guard.

The publication states these as facts about the implemented step (Sections 5,
8.3, 8.4, 9). Until now each was measured only inside one scattered
experiment, so a change to the stage order, the bond convention or the sector
basis could break one and nothing would notice. These are the invariants the
headline claims rest on. They run on the numpy core so they never contend with
a GPU run, and they assert the numbers the paper prints, not their signs.
"""
import os

os.environ["CL10_BACKEND"] = "numpy"

import numpy as np  # noqa: E402

from cl10 import supercell as SC  # noqa: E402

QS = [SC.GAMMA, SC.X_POINTS[0], SC.M_POINTS[0], SC.R_POINTS[0],
      (0.31, -0.77, 0.44), (0.1, 0.2, 0.3)]


def _two_angle(a, b, q):
    """The six-stage bond product with angle `a` on the intra-block stages and
    `b` on the inter-block ones (the axis-symmetric family of V192)."""
    M = np.eye(SC.NSUB, dtype=complex)
    for even, j in SC.STAGES:
        M = SC.stage_matrix(even, j, q[j], a if even else b) @ M
    return M


def test_step_exactly_unitary_at_every_momentum():
    _, _, U16 = SC.sectors()
    for q in QS:
        M = SC.bloch(q, U16)
        assert np.abs(M @ M.conj().T - np.eye(len(M))).max() < 1e-12
        assert abs(abs(np.linalg.det(M)) - 1.0) < 1e-12


def test_matter_antimatter_spectral_mirror():
    """The mirror is between the two sectors, not inside one: the `16`-sector
    one-step spectrum is the exact negative of the `16bar` one at every
    momentum (Theorem 4)."""
    from cl10 import clifford, engine as E
    B, U32, _ = SC.sectors()
    w = E.asnumpy(clifford.pseudoscalar()).astype(np.complex128)
    ev, evec = np.linalg.eig(w)
    Bb = np.linalg.qr(evec[:, np.argsort(np.imag(ev))[:16]])[0]
    V16 = np.kron(np.eye(SC.NSUB), B)
    V16b = np.kron(np.eye(SC.NSUB), Bb)
    for q in QS:
        M = SC.bloch(q, U32)
        p16 = np.sort(np.angle(np.linalg.eigvals(V16.conj().T @ M @ V16)))
        p16b = np.sort(-np.angle(
            np.linalg.eigvals(V16b.conj().T @ M @ V16b)))
        assert np.abs(p16 - p16b).max() < 1e-9


def test_16_sector_closed_by_the_full_step():
    """The bond stages act on the sub-site index and the internal layer is
    even, so the full step leaves the 16 sector invariant (Proposition 2)."""
    B, U32, _ = SC.sectors()
    V16 = np.kron(np.eye(SC.NSUB), B)
    P = V16 @ V16.conj().T
    M = SC.bloch((0.31, -0.77, 0.44), U32)
    assert np.abs((np.eye(P.shape[0]) - P) @ M @ P).max() < 1e-12


def test_front_speed_exactly_one_planck_cell_per_stage():
    """Strict causality, and the front reaching the bound: three blocks per
    step, one Planck cell per elementary stage, nothing outside."""
    r = SC.run(L=8, n_steps=3, x0=(4, 4, 4), sector=True)
    assert abs(r["final_norm"] - 1.0) < 1e-12
    for s in r["steps"]:
        assert s["n_outside"] == 0
        assert s["max_radius"] == 3 * s["k"]


def test_bond_product_covariant_at_gamma_broken_elsewhere():
    """Section 8.4 rests on both halves of this: the six-stage product is
    exactly O_h-equivariant at the zone centre and genuinely broken at every
    other momentum - the breaking is the non-commutativity of the stages, not
    a numerical artifact."""
    eq = SC.cubic_equivariance()
    assert eq[SC.GAMMA] < 1e-14
    for q in [SC.X_POINTS[0], SC.M_POINTS[0], SC.R_POINTS[0],
              (0.31, -0.77, 0.44)]:
        assert eq[tuple(round(v, 4) for v in q)] > 1e-2


def test_two_angle_family_is_covariant_at_gamma():
    """The premise of V192: every axis-symmetric two-angle product keeps the
    zone centre exactly covariant, whatever the two angles are - so the node
    analysis at Gamma applies to the whole family, and what separates its
    members is the spectrum, not the symmetry."""
    W = [W for _, W in SC.oh_subsite_matrices()]
    rng = np.random.default_rng(0)
    for _ in range(6):
        a, b = rng.uniform(0.05, np.pi - 0.05, 2)
        P = _two_angle(a, b, SC.GAMMA)
        assert max(np.abs(Wk @ P @ Wk.T - P).max() for Wk in W) < 1e-14


def test_stage_order_is_the_causal_one():
    """The stage schedule is not decoration: it is three intra-block and three
    inter-block stages, one per axis in each parity, which is what makes one
    step a radius-three operator on the block lattice."""
    assert len(SC.STAGES) == 6
    assert sorted(j for even, j in SC.STAGES if even) == [0, 1, 2]
    assert sorted(j for even, j in SC.STAGES if not even) == [0, 1, 2]
