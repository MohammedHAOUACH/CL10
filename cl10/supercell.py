#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 - **the automaton**: a quantum cellular automaton on a cubic lattice of
**2x2x2 blocks** (the *supercell*), one **time-independent** step, on-site
space $\mathbb C^{8}\otimes\mathbb C^{32}$.

**Why the block is the cell.**  The internal layer is the product of the 45
even Spin(10) plane rotations, so it commutes with the pseudoscalar $w$ and
the 32 splits as $16\oplus\overline{16}$ *dynamically*, not just algebraically.
The transport layer is a **staggered** swap-rotation: at each stage it acts on
the bonds of one parity along one axis, and the parity alternates.  A signal
therefore sees a spatial pattern that repeats after **two** cells in each
direction, and the one-step operator has period **six** in time (two parities
$\times$ three axes).  The translation group of the automaton is thus not the
one-cell group: it is the group of the **2x2x2 block**.  The automaton is
written directly in that language - the cell is the block, the step is the
six-stage cycle folded into **one** time-independent operator, and the
staggering is not a time dependence but the **sub-lattice structure of the
cell**.

**The cell.**  Eight sub-sites $e = (e_1,e_2,e_3) \in \{0,1\}^3$ (the position
inside the block, in units of the Planck spacing; the sub-site index is the
integer $e = e_1 + 2 e_2 + 4 e_3$, so the bit of axis $j$ is $2^j$) carrying
the 36-D cell module $M = S \oplus R_A \oplus R_B$ of Axiom 1 (cl10/cell.py):
the 32-D $\mathrm{Cl}(10,0)$ Dirac module $S$ plus the two neutral singlet
sectors, each a conjugate pair carrying $\pm$ its rung phase (the ladder);
on-site space $\mathbb C^{288}$.  The 256- and 128-dimensional matrices below
are the spinor
block, which carries all of the geometry.  The
internal layer $U_{32} = \prod_{a<b}[\cos g_{10} - i \sin g_{10}\, T_{ab}]$
acts on every sub-site.  A **stage** of the step is a swap-rotation
$\cos g_{10} - i\sin g_{10} P$ on a set of bonds:

* an **intra-block** stage (parity 0 along axis $j$): the four bonds
  $(e, e \oplus 2^j)$ with $e_j = 0$, all inside the block - a $256\times256$
  **one-cell** operator;
* an **inter-block** stage (parity 1 along axis $j$): the bonds joining
  sub-site $e_j = 1$ of block $R$ to sub-site $e_j = 0$ of block $R + \hat
  e_j$ - still a two-cell operator, now between neighbouring blocks.

The step is the six stages in the cycle ``STAGES``
(intra-$x$, inter-$y$, intra-$z$, inter-$x$, intra-$y$, inter-$z$) with the
internal layer applied once per stage.  It is a product of exactly unitary
two-cell operators, so the step is **exactly unitary** and **strictly causal**:
one block of budget per stage, hence a light cone of Manhattan radius $3k$
blocks after $k$ steps - which is exactly **one Planck cell per elementary
step**, the same front speed as the unfolded description.

**Momentum.**  Because the step is one operator, the Bloch decomposition is
exact at **every** step: $M(q)$ is the $256\times256$ (or $128\times128$ on
the 16 sector) one-step matrix at block momentum $q$, and the eight sub-sites
are its internal structure.  A hop across a block face carries the phase
$\exp(\pm iq)$, so $q$ is $2\pi$-periodic and the zone is
$q\in[-\pi,\pi]^3$.  The cubic high-symmetry momenta are the carriers of the
spectrum:

    GAMMA  (0,0,0)                    gapped
    X      (+-pi/2, 0, 0)             star 6,  little group order 8
    M      (+-pi/2, +-pi/2, 0)        star 12, little group C2v (4) - the gapless nodes
    R      (+-pi/2, +-pi/2, +-pi/2)   star 8,  little group C3v (6) - exactly
                                      3-fold degenerate multiplets at every corner

What the unfolded one-cell description calls the Nielsen-Ninomiya doublers is,
in the cell of the automaton, the band structure at the 12 **M** momenta, with
paired chiralities (their sum is 0, as the theorem requires).  The 8 **R**
corners - one $O_h$ star - carry exactly 3-fold degenerate multiplets: the
only exactly massless multiplets of the step, and the only exact degeneracy of
the zone.

**The one number.**  The only parameter of the whole step is the rotation angle
$g_{10}$, and the closure axiom fixes it (Sections 1.4 and 9): the physical clock
is the **closed** root $g_{10} = \pi/4$, the order-8 step, where the zone
centre is pinned to the Floquet pair and $M(0)^2 = I$.  The reference angle
$g_{10} = 1/5$ (the convention $1/n_{\rm qubits}$, the dispersive Wilson-family
member of the Section 8 searches) and the degenerate second root $\pi/2$ stay
selectable with ``CL10_ANGLE=implemented|closed|half``; the default is ``closed``.
The cell, the cycle and the staggering are fixed by the algebra.
"""

import math

import numpy as np

from . import clifford
from . import config as _config
from . import engine as E
from . import sm_qca
from . import coupling
from . import cell
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml

# ---------------------------------------------------------------------------
# the cell
# ---------------------------------------------------------------------------
NSUB = int(_config.CFG.physics.n_subsites)      # sub-sites per block (2x2x2)
NINT = int(_config.CFG.physics.n_internal)      # the Cl(10,0) module per sub-site
NSECTOR = int(_config.CFG.physics.n_sector)     # the Weyl sector (the +i eigenspace of w)
NONSITE = NSUB * NINT          # 256
NSECTOR_ONSITE = NSUB * NSECTOR  # 128

_config.enforce_angle()   # a CL10_ANGLE pin set after config was imported is
                         # invisible here: fail instead of binding the wrong clock
G10 = coupling.G10        # the active clock: angle.default_branch in
                          # cl10/config.yaml, overridden by CL10_ANGLE only when
                          # it is set BEFORE cl10 is imported for the first time

# the six stages of the one step: (intra-block?, axis).
# The order is not free: the product of the six stage matrices is exactly
# O_h-equivariant at GAMMA and broken at every other momentum (see
# cubic_equivariance below, and Section 8.4); the stage SUM is invariant at
# all momenta, so the breaking is the non-commutativity of the stages.
STAGES = [(True, 0), (False, 1), (True, 2),
          (False, 0), (True, 1), (False, 2)]

# ---------------------------------------------------------------------------
# the Brillouin zone of the block lattice: q in [-pi, pi]^3, the crystal
# momentum conjugate to a translation by one cell (the half-zone points
# pi/2 are the X, M, R high-symmetry momenta)
# ---------------------------------------------------------------------------
PI2 = math.pi / 2
GAMMA = (0.0, 0.0, 0.0)
X_POINTS = [q for p in ((PI2, 0.0, 0.0), (0.0, PI2, 0.0), (0.0, 0.0, PI2))
            for q in (p, tuple(-v for v in p))]
M_POINTS = [q for sx in (1, -1) for sy in (1, -1)
            for q in ((sx * PI2, sy * PI2, 0.0),
                      (sx * PI2, 0.0, sy * PI2),
                      (0.0, sx * PI2, sy * PI2))]          # 12, star 12, C2v
R_POINTS = [(sx * PI2, sy * PI2, sz * PI2)
            for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]  # 8, star 8
HIGH_SYM = {"GAMMA": [GAMMA], "X": X_POINTS, "M": M_POINTS, "R": R_POINTS}


def _signed_permutations():
    """The 48 elements of O_h as index/sign maps on the momentum components."""
    import itertools
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            out.append((perm, signs))
    return out


_OH = _signed_permutations()


def star(q):
    r"""The O_h orbit (the star) of the block momentum ``q``, with the
    components folded back into $[-\pi,\pi)$ (the block momentum is
    $2\pi$-periodic: the Bloch phase of a hop across a block face is
    $\exp(\pm i q)$)."""
    def norm(qq):
        return tuple(round(((v + math.pi) % (2 * math.pi)) - math.pi, 9)
                     for v in qq)
    pts = set()
    for perm, signs in _OH:
        pts.add(norm(tuple(signs[i] * q[perm[i]] for i in range(3))))
    return sorted(pts)


def stabilizer_order(q):
    """The order of the little group of ``q`` in O_h (48 / |star|)."""
    return 48 // len(star(q))


# the stars of the zone: the node momenta form a 12-point star (little group
# C2v, order 4); the 8 signed corners are one star (little group C3v, order 6).
GAMMA_STAR = [GAMMA]
X_STAR = star(X_POINTS[0])
M_STAR = star(M_POINTS[0])
R_STAR = star(R_POINTS[0])


# ---------------------------------------------------------------------------
# the internal layer and the sector basis
# ---------------------------------------------------------------------------

def internal(g10=G10):
    """The 32x32 internal unitary on the active backend."""
    return sm_qca.gut_qca_step(sm_qca.plane_generators(), g10)


def sectors(g10=G10):
    """The 16 / 32 internal unitaries and the 16 basis (the +i eigenspace of w).

    Returns ``(B, U32, U16)`` with ``U16 = B^dag U32 B`` (the 16 sector is
    exactly invariant: the 45 generators are even, so they commute with w).
    """
    w = E.asnumpy(clifford.pseudoscalar()).astype(np.complex128)
    ev, evec = np.linalg.eig(w)
    B = np.linalg.qr(evec[:, np.argsort(-np.imag(ev))[:16]])[0]
    U32 = E.asnumpy(internal(g10)).astype(np.complex128)
    return B, U32, B.conj().T @ U32 @ B


# ---------------------------------------------------------------------------
# the bond stages
# ---------------------------------------------------------------------------

def stage_matrix(even, j, q, g10=G10):
    """The 8x8 sub-site matrix of one bond stage at block momentum ``q``.

    An intra-block stage (``even=True``) carries no phase; an inter-block
    stage carries $\\exp(\\pm i q_j)$ for a hop in $\\pm\\hat e_j$."""
    c, s = math.cos(g10), math.sin(g10)
    S = np.eye(NSUB, dtype=complex)
    bit = 1 << j
    for e in range(NSUB):
        if (e & bit) == (0 if even else bit):
            o = e ^ bit
            ph = 1.0 if even else (np.exp(1j * q) if (e & bit) else np.exp(-1j * q))
            S[e, e] = c
            S[o, o] = c
            S[e, o] = -1j * s * ph
            S[o, e] = -1j * s * np.conj(ph)
    return S


def bond_product(q, g10=G10, stages=STAGES):
    """The 8x8 bond part P(q) of the one-step Bloch matrix: the six stage
    matrices multiplied in stage order. The full step is ``bloch(q, U)``;
    when the internal angle equals the bond angle it factorizes as
    P(q) (x) U(g10)**ns."""
    M = np.eye(NSUB, dtype=complex)
    for even, j in stages:
        M = stage_matrix(even, j, q[j], g10) @ M
    return M


def oh_subsite_matrices():
    """The 48 elements of O_h as permutation matrices on the eight sub-sites:
    sub-site index e has bit j = its coordinate along axis j, and a sign flip
    is the reflection about the block centre, x -> 1 - x."""
    out = []
    for perm, signs in _OH:
        W = np.zeros((NSUB, NSUB))
        for e in range(NSUB):
            x = [((e >> j) & 1) * 2 - 1 for j in range(3)]
            y = [signs[i] * x[perm[i]] for i in range(3)]
            e2 = sum(int((v + 1) // 2) << j for j, v in enumerate(y))
            W[e2, e] = 1.0
        out.append(((perm, signs), W))
    return out


def cubic_equivariance(qs=(GAMMA, X_POINTS[0], M_POINTS[0], R_POINTS[0],
                           (0.31, -0.77, 0.44)), g10=G10, stages=STAGES):
    """The defect of the O_h equivariance of the bond product,
    ``W P(q) W^T = P(R q)``, maximized over the 48 group elements, per
    momentum.

    Measured (V191, V192): the engine's six-stage product is exact at
    GAMMA (1.1e-16) and broken at every other momentum (1.07e-1 generic,
    1.40e-1 M, 1.46e-1 R), while its stage SUM is invariant to 1.1e-15 at all
    momenta. Section 8.4 uses exactly this: the node analysis at GAMMA is a
    statement about an O_h-invariant matrix, and it is valid there."""
    out = {}
    for q in qs:
        P = bond_product(q, g10, stages)
        d = 0.0
        for (perm, signs), W in oh_subsite_matrices():
            rq = tuple(signs[i] * q[perm[i]] for i in range(3))
            d = max(d, float(np.abs(W @ P @ W.T
                                    - bond_product(rq, g10, stages)).max()))
        out[tuple(round(v, 4) for v in q)] = d
    return out


def bloch(q, U, g10=G10, stages=STAGES):
    """The exact one-step Bloch matrix of the automaton at block momentum
    ``q``: ``(8*n, 8*n)`` for an ``(n,n)`` internal unitary (128 on the 16
    sector, 256 on the full module).

    The build goes through ``E.kron``/``E.matmul``, so it runs where ``U``
    lives: host for the 128/256 matrices of the momentum sweeps (measured
    faster there, and the sweeps are parallelised by ``cl10/parpool.py``),
    device if a caller hands in a device ``U``."""
    M = np.eye(NSUB * U.shape[0], dtype=complex)
    for even, j in stages:
        S = stage_matrix(even, j, q[j], g10)
        M = E.matmul(E.kron(S, U), M)
    return M


def phases(M):
    """The sorted eigenphases.  The eigendecomposition follows the engine
    routing (``routing.eig_gpu_min_dim``): at 128/256 the host is measured
    faster, above the threshold the device is used."""
    return np.sort(np.angle(E.asnumpy(E.eigvals(M))))


def gap(q, U, g10=G10, stages=STAGES):
    return float(np.min(np.abs(phases(bloch(q, U, g10, stages)))))


def degeneracy(q, U, g10=G10, stages=STAGES, tol=_TOL.standard):
    """The eigenphase multiplet structure at ``q``: [(gap, multiplicity)]."""
    a = np.sort(np.abs(phases(bloch(q, U, g10, stages))))
    out = []
    for v in a:
        if out and abs(v - out[-1][0]) < tol:
            out[-1][1] += 1
        else:
            out.append([float(v), 1])
    return [(v, m) for v, m in out]


# ---------------------------------------------------------------------------
# the real-space step (runs on the active backend: cupy / jax / numpy)
# ---------------------------------------------------------------------------

def internal_layer(fields, U):
    """The internal layer: ``U`` on every sub-site of every block.
    ``fields`` has shape ``(L,L,L,8,32)`` (or ``(L,L,L,8,16)``)."""
    return E.einsum("...ai,ij->...aj", fields, U)


def bond_layer(fields, even, j, theta=G10):
    """One bond stage: the swap-rotation on the bonds of parity ``even``
    along axis ``j``.  Intra-block stages act inside the cell; inter-block
    stages join sub-site $e_j=1$ of a block to sub-site $e_j=0$ of the next
    block (a two-cell operator, one block of budget)."""
    c, s = math.cos(theta), math.sin(theta)
    bit = 1 << j
    if even:
        # the stage matrix stays on the host: einsum then runs where the FIELD
        # lives (a device field pulls it over), so a host field never becomes a
        # device field through this constant.
        S = stage_matrix(True, j, 0.0, theta)
        return E.einsum("ab,...bi->...ai", S, fields)
    hi = [e for e in range(NSUB) if e & bit]
    lo = [e for e in range(NSUB) if not (e & bit)]
    f_plus = E.roll(fields, -1, j)          # the field of the next block
    f_minus = E.roll(fields, 1, j)          # the field of the previous block
    out = c * fields
    out[..., hi, :] = out[..., hi, :] - 1j * s * f_plus[..., [e ^ bit for e in hi], :]
    out[..., lo, :] = out[..., lo, :] - 1j * s * f_minus[..., [e ^ bit for e in lo], :]
    return out


def stage(fields, U, t, theta=G10):
    """One elementary **stage** of the cycle: the internal layer plus the bond
    stage ``STAGES[t % 6]``.  Six stages make one step; a stage moves a signal
    at most one block along its own axis (one Planck cell of budget)."""
    even, j = STAGES[t % len(STAGES)]
    return bond_layer(internal_layer(fields, U), even, j, theta)


def step(fields, U, theta=G10, stages=STAGES):
    """One automaton **step** (six stages, each = internal layer + bond stage).
    ``fields`` has shape ``(L,L,L,8,n)`` on a torus."""
    f = fields
    for even, j in stages:
        f = internal_layer(f, U)
        f = bond_layer(f, even, j, theta)
    return f


def run(L=8, n_steps=4, x0=(4, 4, 4), e0=0, theta=G10, seed=0,
        tol=_TOL.strict,
        sector=True):
    """Run the automaton from a state localised on **one sub-site** of one
    block; return the per-step support and its light-cone radius (in blocks).

    ``sector=True`` evolves the 16 sector (128 on-site), ``sector=False`` the
    full module (256)."""
    n = NSECTOR if sector else NINT
    B, U32, U16 = sectors(theta)
    U = E.Array_(U16 if sector else U32)
    rng = np.random.default_rng(seed)
    f = E.zeros((L, L, L, NSUB, n))
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    f[x0 + (e0,)] = E.Array_(v / np.linalg.norm(v))
    steps = []
    for k in range(1, n_steps + 1):
        f = step(f, U, theta)
        a = E.asnumpy(f)
        nrm = np.sqrt(np.sum(np.abs(a) ** 2, axis=(-2, -1)))
        support = np.argwhere(nrm > tol)
        d = []
        for (x, y, z) in support:
            dd = 0
            for a0, a1 in ((x, x0[0]), (y, x0[1]), (z, x0[2])):
                delta = abs(a0 - a1)
                dd += min(delta, L - delta)
            d.append(dd)
        steps.append({"k": k, "n_support": int(len(support)),
                      "max_radius": int(max(d)) if d else 0,
                      "n_outside": int(sum(1 for r in d if r > 3 * k))})
    return {"backend": E._BACKEND, "L": int(L), "n_steps": int(n_steps),
            "theta": float(theta), "g10": float(G10), "sector": bool(sector),
            "steps": steps,
            "final_norm": E.field_norm(f)}


# ---------------------------------------------------------------------------
# the invariants of the cell
# ---------------------------------------------------------------------------

def verify(tol=_TOL.standard):
    """The invariants of the automaton as written on the block.

    1. the one-step Bloch matrix is **exactly unitary** at generic momenta
       and at every high-symmetry point;
    2. the step is **strictly causal**: the support stays inside the light
       cone of 3 blocks per step (= one Planck cell per elementary stage) and
       the front reaches it exactly;
    3. the **16 sector is closed by the step** (the bond stages act on the
       sub-site index, the internal layer is even): the 16 -> 16bar leakage
       is 0, and the 16bar spectrum is the exact conjugate of the 16;
    4. the **zone structure of the active branch** (Sections 8.3 and 9):
       closed (g = pi/4, the root clock) - GAMMA pinned to the Floquet pair
       {0, pi} with a gapless band, the 12 **M** channels exactly gapless
       (64+64 of 128 quasi-energies, V203), and the zone-centre step an
       involution, M(0)^2 = I;
       implemented (g = 1/5, the reference Wilson-family member) - GAMMA gapped
       (0.0365), the 12 **M** momenta near-gapless (2.0e-4), every **R** corner
       an exactly 3-fold multiplet;
       half (g = pi/2) - M(0) = I on the whole 128-dimensional space: the bond
       layer is a pure translation (V134);
    5. the stars and little groups: the M momenta form a 12-point star (little
       group C2v, order 4), the 8 corners one star (little group C3v, order 6).
    6. the **cubic equivariance of the bond product**: exact at GAMMA, broken
       at every other momentum (Section 8.4).
    7. the **cell module of Axiom 1** (cl10/cell.py): the 36-level step - the
       spinor block plus the singlet conjugate pairs at their rung phases -
       is exactly unitary with +- paired spectrum and det = 1, and the
       B-parity commutes with it and with all 45 plane generators extended
       by zero (the dark matter is stable by this invariant, Section 12.2).
    """
    B, U32, U16 = sectors()
    r = {"backend": E._BACKEND, "g10": float(G10), "branch": coupling.branch()}

    # 1. unitarity of the one-step matrix
    uerr = 0.0
    for q in [GAMMA, (0.31, -0.77, 0.44), M_POINTS[0], R_POINTS[0],
              X_POINTS[0], (0.1, 0.2, 0.3)]:
        M = bloch(q, U16)
        uerr = max(uerr, float(np.abs(M @ M.conj().T - np.eye(len(M))).max()))
    r["step_unitary_err"] = uerr

    # 2. causality on the block lattice
    rr = run(L=8, n_steps=3, x0=(4, 4, 4), sector=True)
    r["causal_err"] = abs(rr["final_norm"] - 1.0)
    r["causal_ok"] = all(s["n_outside"] == 0 and s["max_radius"] == 3 * s["k"]
                         for s in rr["steps"])
    r["front_radius"] = [s["max_radius"] for s in rr["steps"]]
    r["support_growth"] = [s["n_support"] for s in rr["steps"]]

    # 3. the 16 sector is closed; the mirror is exact
    qg = (0.31, -0.77, 0.44)
    M256 = bloch(qg, U32)
    # the 16bar block: Bbar = the -i eigenspace of w
    w = E.asnumpy(clifford.pseudoscalar()).astype(np.complex128)
    ev, evec = np.linalg.eig(w)
    Bb = np.linalg.qr(evec[:, np.argsort(np.imag(ev))[:16]])[0]
    V16 = np.kron(np.eye(NSUB), B)        # 256 x 128: the 16 sector of the cell
    V16b = np.kron(np.eye(NSUB), Bb)      # the 16bar sector
    # the leakage: the 16 <-> 16bar off-diagonal block of the one-step matrix
    leak = float(np.abs(V16.conj().T @ M256 @ V16b).max())
    # the 16-sector step IS the 128-level Bloch matrix
    block_err = float(np.abs(V16.conj().T @ M256 @ V16
                             - bloch(qg, U16)).max())
    U16b = Bb.conj().T @ U32 @ Bb
    sp16 = np.sort(phases(bloch(qg, U16)))
    sp16b = np.sort(phases(bloch(qg, U16b)))
    r["sector_leakage"] = leak
    r["sector_block_err"] = block_err
    r["spectral_mirror_err"] = float(np.abs(sp16 + sp16b[::-1]).max())

    # 4. the zone structure
    r["gap_gamma"] = gap(GAMMA, U16)
    mg = [gap(q, U16) for q in M_POINTS]
    r["n_M_gapless"] = int(sum(1 for g in mg if g < 1e-3))
    r["max_M_gap"] = float(max(mg))
    r["M_star_gaps"] = [gap(q, U16) for q in M_STAR]
    r["R_degeneracy"] = [[v, m] for v, m in degeneracy(R_POINTS[0], U16)[:4]]
    r["R_first_mult"] = int(degeneracy(R_POINTS[0], U16)[0][1])
    r["R_mult_all"] = sorted({degeneracy(q, U16)[0][1] for q in R_POINTS})
    # the closure of the one-step matrix at the zone centre (Section 9):
    # M(0)^2 = I at the closed angle, M(0) = I at the half angle
    M0 = bloch(GAMMA, U16)
    r["M0_sq_minus_I"] = float(np.abs(M0 @ M0 - np.eye(len(M0))).max())
    r["M0_minus_I"] = float(np.abs(M0 - np.eye(len(M0))).max())

    # 5. the little groups of the block zone
    r["little_M"] = stabilizer_order(M_POINTS[0])
    r["little_R"] = stabilizer_order(R_POINTS[0])
    r["star_M"] = len(star(M_POINTS[0]))
    r["star_R"] = len(star(R_POINTS[0]))

    # 6. the cubic equivariance of the bond product (Section 8.4)
    eq = cubic_equivariance()
    r["cubic_equivariance"] = {str(k): v for k, v in eq.items()}
    r["equivariant_at_gamma"] = eq[GAMMA] < tol

    # 7. the cell module of Axiom 1 (cl10/cell.py): the 36-level step and its
    # exact B-parity - unitary, and odd on the dark-matter sector alone
    Ucell = cell.cell_step_np(G10)
    Gcell = cell.parity_np()
    r["cell_unitary_err"] = float(np.abs(Ucell @ Ucell.conj().T - np.eye(cell.DIM)).max())
    r["cell_parity_step_comm"] = float(np.abs(Gcell @ Ucell - Ucell @ Gcell).max())
    r["cell_parity_gauge_max"] = float(max(
        np.abs(Gcell @ X - X @ Gcell).max() for X in cell.gauge_ext_np()))

    zone = {
        "implemented": (r["gap_gamma"] > 1e-3 and r["n_M_gapless"] == 12
                        and r["max_M_gap"] < 1e-3 and r["R_first_mult"] == 3),
        "closed": (r["gap_gamma"] < _TOL.strict and r["n_M_gapless"] == 12
                   and r["max_M_gap"] < _TOL.strict and r["M0_sq_minus_I"] < _TOL.strict),
        "half": (r["M0_minus_I"] < _TOL.strict and r["n_M_gapless"] == 12),
    }[coupling.branch()]
    ok = (uerr < tol and r["causal_ok"] and r["causal_err"] < tol
          and leak < tol and block_err < tol
          and r["spectral_mirror_err"] < _TOL.standard
          and zone and r["little_M"] == 4
          and r["little_R"] == 6 and r["star_M"] == 12 and r["star_R"] == 8
          and r["equivariant_at_gamma"]
          and r["cell_unitary_err"] < tol and r["cell_parity_step_comm"] < tol
          and r["cell_parity_gauge_max"] < tol)
    r["ok"] = bool(ok)
    return r


def main():
    import sys
    r = verify()
    print("=" * 68)
    print("CL10 - the automaton on the 2x2x2 block: one step, one cell")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  on-site space:  {NSUB} sub-sites x {NINT} internal = {NONSITE}")
    print(f"  the angle:  branch {r['branch']}   g_10 = {r['g10']:.6g}")
    print(f"  one-step unitarity:  {r['step_unitary_err']:.2e}")
    print(f"  light cone (blocks after k steps):  {r['front_radius']}  "
          f"(3 blocks = 6 Planck cells per step)")
    print(f"  support growth:  {r['support_growth']}")
    print(f"  16 sector closed by the step:  leakage {r['sector_leakage']:.2e}")
    print(f"  16bar spectrum = conjugate of 16:  {r['spectral_mirror_err']:.2e}")
    print(f"  gap at GAMMA:  {r['gap_gamma']:.4f}   gapless M points: "
          f"{r['n_M_gapless']}/12 (max gap {r['max_M_gap']:.2e})")
    print(f"  zone-centre closure:  |M(0)^2 - I| = {r['M0_sq_minus_I']:.2e}   "
          f"|M(0) - I| = {r['M0_minus_I']:.2e}")
    print(f"  R-point multiplet:  {r['R_first_mult']}-fold (all R: {r['R_mult_all']})")
    print(f"  block-zone stars:  |star(M)| = {r['star_M']} (little group {r['little_M']}), "
          f"|star(R)| = {r['star_R']} (little group {r['little_R']})")
    print(f"  cell (36-level, cl10/cell.py):  unitarity {r['cell_unitary_err']:.2e}, "
          f"[G,U] {r['cell_parity_step_comm']:.2e}, [G,45] {r['cell_parity_gauge_max']:.2e}")
    print(f"  VERDICT:  {'PASS' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
