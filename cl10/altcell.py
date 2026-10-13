#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 alt-cell exploration: the generalizations the alt_cell census (V149-V160)
needs, WITHOUT touching ``clifford``/``sm_qca``/``supercell``.

Two independent axes are generalized here:

1. **The signature.**  ``signature_generators(q)`` is the Cl(10-q, q) set:
   the 10 qubit generators of ``clifford.generators()`` with the first ``q``
   multiplied by ``i`` (anti-Hermitian, square -I).  Everything downstream
   (plane generators, the 45-plane product, the pseudoscalar, the 16/16-bar
   split) is recomputed per signature.

2. **The cell.**  ``cell(shape)`` is the mixed-radix block (the bit trick of
   ``supercell`` is the shape (2,2,2) special case - the regression anchor);
   ``lattice_cell(kind)`` is an FCC/diamond basis cell with nearest-neighbour
   bonds and image labels; ``one_factorizations`` enumerates the ordered
   partitions of the bond set into perfect matchings (the stage cycles).

A cell is fully described by its **bond stages**: a stage is a list of bonds
``(i, j, v)`` - the amplitude of sub-site ``j`` of the block displaced by the
integer cell vector ``v`` hops to sub-site ``i`` of the block at 0, carrying
the Bloch phase ``exp(+i q . v)``.  This one representation covers the block
stagger (intra: v = 0; inter: v = +-e_j) and the lattice matchings (v = the
image label of the NN bond).

The module is exploration-only (tier 4): it never enters the hard-tier
verdict and no existing module imports it.
"""

import itertools
import math

import numpy as np

from . import clifford
from . import supercell as SC

# ---------------------------------------------------------------------------
# 1. the signature axis: Cl(10-q, q)
# ---------------------------------------------------------------------------


def signature_generators(q=0):
    """The 10 Cl(10-q, q) generators: ``h_a^2 = eps_a I`` with
    ``eps_a = -1`` for ``a < q``, ``+1`` otherwise (host complex128)."""
    h = [np.asarray(m, dtype=np.complex128) for m in clifford.generators()]
    return [1j * m if a < q else m for a, m in enumerate(h)]


def signature_squares(q=0):
    """The 10 signs eps_a of the signature Cl(10-q, q)."""
    return [-1.0] * q + [1.0] * (10 - q)


def plane_generators_sig(q=0):
    """The 45 ``T_ab = i h_a h_b`` of the signature: Hermitian iff the two
    generators share their square-class, anti-Hermitian otherwise;
    ``T_ab^2 = eps_a eps_b I``."""
    h = signature_generators(q)
    return [1j * (h[a] @ h[b]) for a in range(10) for b in range(a + 1, 10)]


def plane_product(planes, g):
    """The implemented 45-plane product ``Pi [cos g - i sin g T_ab]``
    (the ``sm_qca.gut_qca_step`` order, host matrices)."""
    n = planes[0].shape[0]
    I = np.eye(n, dtype=np.complex128)
    U = I.copy()
    c, s = math.cos(g), math.sin(g)
    for T in planes:
        U = U @ (c * I - 1j * s * T)
    return U


def unitary_err(U):
    """``||U U^dag - I||_max``."""
    return float(np.abs(U @ U.conj().T - np.eye(U.shape[0])).max())


def exact_order(U, Nmax=64):
    """The smallest exact period N <= Nmax of ``U`` (V131 method):
    ``(N, err)``; ``(0, best_err)`` if none reaches 1e-9."""
    I = np.eye(U.shape[0], dtype=np.complex128)
    best = (0, 9.9)
    P = I.copy()
    for N in range(1, Nmax + 1):
        P = P @ U
        e = float(np.abs(P - I).max())
        if e < best[1]:
            best = (N, e)
        if e < 1e-9:
            break
    return best


def pseudoscalar_sig(q=0):
    """``w = h_1 ... h_10`` of the signature; ``w^2 = -prod(eps_a)``."""
    h = signature_generators(q)
    w = h[0].copy()
    for m in h[1:]:
        w = w @ m
    return w


def weyl_split_sig(q=0):
    """The 16/16-bar of the signature: ``(B_plus, B_minus, w2)``.

    ``w^2 = -(-1)^q``: for even q the eigenvalues are +-i (the q=0 case,
    the +i eigenspace is the 16); for odd q they are +-1 (the +1
    eigenspace is the 16).  Both are 16/16 (w is traceless)."""
    w = pseudoscalar_sig(q)
    w2 = float(np.real(np.trace(w @ w)) / w.shape[0])
    ev, evec = np.linalg.eig(w)
    key = np.imag(ev) if w2 < 0 else np.real(ev)
    order = np.argsort(-key)
    B = np.linalg.qr(evec[:, order[:16]])[0]
    Bb = np.linalg.qr(evec[:, order[16:]])[0]
    return B, Bb, w2


def sectors_sig(q=0, g=None):
    """The internal step of the signature on the 16 sector:
    ``(B, U32, U16)`` with ``U32`` the 45-plane product at ``g``
    (default ``SC.G10``) and ``U16 = B^dag U32 B``.

    The planes are even (grade 2), so they commute with ``w`` for every
    signature: the 16 split is invariant whenever it exists."""
    g = SC.G10 if g is None else g
    B, _Bb, _w2 = weyl_split_sig(q)
    U32 = plane_product(plane_generators_sig(q), g)
    return B, U32, B.conj().T @ U32 @ B


# ---------------------------------------------------------------------------
# 2. the cell axis: mixed-radix blocks
# ---------------------------------------------------------------------------


def cell(shape):
    """The mixed-radix block cell ``shape = (n1, n2, n3)``:
    ``nsub`` sub-sites, ``digit(e, j)`` the j-th mixed-radix digit of ``e``,
    ``site_index(digits)`` its inverse (the bit trick at shape (2,2,2))."""
    shape = tuple(int(n) for n in shape)
    radix = [int(np.prod(shape[:j])) for j in range(3)]
    nsub = int(np.prod(shape))

    def digit(e, j):
        return (e // radix[j]) % shape[j]

    def site_index(digits):
        return sum(int(d) * r for d, r in zip(digits, radix))

    return {"kind": "block", "shape": shape, "nsub": nsub, "radix": radix,
            "digit": digit, "site_index": site_index,
            "sites": [tuple(digit(e, j) for j in range(3))
                      for e in range(nsub)]}


def block_bonds(c, j, phase):
    """The bond list of one block stage: every sub-site with digit
    ``d_j = phase`` hops to digit ``d_j + 1 mod n_j``; a wrap (``d + 1 =
    n_j``) is an inter-block bond with image label ``+e_j`` (the Bloch
    phase ``exp(+i q_j)`` - the ``supercell.stage_matrix`` convention)."""
    n_j = c["shape"][j]
    r = c["radix"][j]
    bonds = []
    for e in range(c["nsub"]):
        if c["digit"](e, j) == phase:
            if phase + 1 < n_j:
                bonds.append((e, e + r, (0, 0, 0)))
            else:
                bonds.append((e, e - (n_j - 1) * r, (0, 0, 0)))
                bonds[-1] = (e, e - (n_j - 1) * r,
                             tuple(1 if k == j else 0 for k in range(3)))
    return bonds


def stages_shape(shape):
    """The canonical stagger of a block: phase-major, ``sum_j n_j`` stages
    ``[(j, phase), ...]`` (for (2,2,2) this is the phase-major cycle, NOT
    the implemented one - see ``stages_baseline``)."""
    shape = tuple(shape)
    return [(j, r) for r in range(max(shape)) for j in range(3)
            if r < shape[j]]


def stages_baseline():
    """The implemented ``supercell.STAGES`` as ``(j, phase)`` pairs
    (phase 0 = intra-block, phase 1 = inter-block)."""
    return [(j, 0 if even else 1) for even, j in SC.STAGES]


def cell_stage_bonds(c, stage):
    """The bonds of one stage of one cell (block: ``(j, phase)``)."""
    j, phase = stage
    return block_bonds(c, j, phase)


# ---------------------------------------------------------------------------
# 3. the cell axis: FCC / diamond basis cells with NN bonds
# ---------------------------------------------------------------------------

FCC_PRIMITIVE = ((0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0))


def lattice_cell(kind):
    """A basis cell of a lattice: ``sites`` = Cartesian positions (units of
    the cubic constant a), ``lattice`` = the primitive cell vectors,
    ``bonds`` = ALL nearest-neighbour bonds ``(i, j, v)`` with the integer
    image label ``v`` in lattice coordinates (``v`` in units of the lattice
    vectors; the real-space displacement is ``r_j - r_i + v . lattice``).

    - ``"fcc"``: the 4-site conventional cubic cell, 24 NN bonds.
    - ``"diamond"``: the 2-site FCC primitive cell (basis +- the (1/4,1/4,1/4)
      shift), 4 NN bonds.
    """
    if kind == "fcc":
        sites = [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5),
                 (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
        lattice = np.eye(3)                      # conventional cubic cell
        nn = math.sqrt(0.5)
    elif kind == "diamond":
        lattice = np.array(FCC_PRIMITIVE)
        sites = [(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)]
        nn = math.sqrt(3.0) / 4.0
    else:
        raise ValueError(f"unknown lattice {kind!r}")
    sites = [tuple(s) for s in sites]
    bonds = []
    for i in range(len(sites)):
        for j in range(len(sites)):
            for n in itertools.product((-1, 0, 1), repeat=3):
                ri = np.array(sites[i])
                rj = np.array(sites[j]) + np.asarray(n, float) @ lattice
                d = float(np.linalg.norm(rj - ri))
                if abs(d - nn) < 1e-9 and i < j:
                    bonds.append((i, j, n))
    return {"kind": "lattice", "lattice_kind": kind, "sites": sites,
            "lattice": lattice, "nsub": len(sites), "bonds": bonds,
            "nn": nn}


def one_factorizations(lc, cap=None):
    """All ordered partitions of the cell bond set into perfect matchings
    (backtracking, lexicographic bond order; each matching is built
    greedily from the lowest uncovered sub-site and covers every sub-site
    exactly once, which fixes the order of the matchings inside a
    factorization).

    Returns ``(factorizations, skipped)``: each factorization is a list of
    stages (each stage a list of bonds); ``skipped`` counts the branches
    cut by ``cap``.  A bond set with no perfect matching yields none."""
    bonds = sorted(lc["bonds"])
    all_sites = {i for (i, _j, _v) in bonds} | \
        {j for (_i, j, _v) in bonds}
    out = []
    skipped = [0]

    def rec(remaining, stages):
        if len(out) == (cap if cap is not None else float("inf")):
            skipped[0] += 1
            return
        if not remaining:
            out.append([list(s) for s in stages])
            return
        uncovered = {i for (i, _j, _v) in remaining} | \
                    {j for (_i, j, _v) in remaining}
        if len(uncovered) % 2 or not uncovered:
            skipped[0] += 1
            return
        v0 = min(uncovered)
        for k, b in enumerate(remaining):
            if v0 not in (b[0], b[1]):
                continue
            i, j, _v = b
            if i not in uncovered or j not in uncovered:
                continue
            stage = [b]
            rest = remaining[:k] + remaining[k + 1:]
            # greedily complete the perfect matching of this stage:
            # every sub-site appears in exactly one bond of the stage
            stack = [(stage, {b[0], b[1]}, rest)]
            while stack:
                if cap is not None and len(out) >= cap:
                    skipped[0] += 1
                    stack.clear()
                    return
                st, cov, rs = stack.pop()
                need = all_sites - cov
                if not need:
                    rec(rs, stages + [st])
                    continue
                w0 = min(need)
                for m, b2 in reversed(list(enumerate(rs))):
                    if w0 not in (b2[0], b2[1]):
                        continue
                    if b2[0] in cov or b2[1] in cov:
                        continue
                    stack.append((st + [b2], cov | {b2[0], b2[1]},
                                  rs[:m] + rs[m + 1:]))

    rec(bonds, [])
    return out, skipped[0]


def lattice_stage_bonds(lc, stage):
    """The bonds of one lattice stage (a perfect matching): the stage is
    already a bond list."""
    return list(stage)


# ---------------------------------------------------------------------------
# 4. one bond representation for every cell: stage matrices, Bloch, real space
# ---------------------------------------------------------------------------


def stage_matrix_bonds(bonds, q, g10=SC.G10, nsub=8):
    """The ``nsub x nsub`` matching matrix of one stage: identity on the
    unpaired sub-sites; on a bond ``(i, j, v)`` the 2x2 rotation
    ``[[c, -i s e^{+i q.v}], [-i s e^{-i q.v}, c]]`` (the
    ``supercell.stage_matrix`` convention)."""
    c, s = math.cos(g10), math.sin(g10)
    S = np.eye(nsub, dtype=complex)
    for (i, j, v) in bonds:
        # q is the full momentum vector, or the scalar momentum along the
        # bond's own axis (the supercell.stage_matrix convention: the image
        # label v is a single +-1 on that axis, so q * sum(v) = q).
        qv = (q * sum(v)) if np.isscalar(q) else np.dot(q, v)
        ph = np.exp(1j * qv)
        S[i, i] = c
        S[j, j] = c
        S[i, j] = -1j * s * ph
        S[j, i] = -1j * s * np.conj(ph)
    return S


def stage_matrix_cell(c, j, phase, q, g10=SC.G10):
    """The generalized ``supercell.stage_matrix`` (regression anchor:
    exact equality at shape (2,2,2), phase 0 = intra, phase 1 = inter)."""
    return stage_matrix_bonds(block_bonds(c, j, phase), q, g10, c["nsub"])


def bloch_cell(c, U, g10, stages, q):
    """The exact one-step Bloch matrix of a generalized cell:
    ``Pi kron(S_stage, U)`` (the ``supercell.bloch`` product structure)."""
    M = np.eye(c["nsub"] * U.shape[0], dtype=complex)
    for st in stages:
        S = stage_matrix_bonds(cell_stage_bonds(c, st), q, g10, c["nsub"])
        M = np.kron(S, U) @ M
    return M


def bloch_lattice(lc, U, g10, factorization, q):
    """The one-step Bloch matrix of a lattice cell: one stage per matching
    of the factorization."""
    M = np.eye(lc["nsub"] * U.shape[0], dtype=complex)
    for stage in factorization:
        S = stage_matrix_bonds(stage, q, g10, lc["nsub"])
        M = np.kron(S, U) @ M
    return M


def bond_layer_bonds(fields, bonds, g10=SC.G10):
    """One real-space bond stage on a torus: ``fields`` has shape
    ``(L, L, L, nsub, n)``; a bond ``(i, j, v)`` (v in **cell** units,
    integer) sends the field of cell ``x + v`` to sub-site ``i`` of cell
    ``x`` (the ``supercell.bond_layer`` convention).  The Hermitian
    partner ``(j, i, -v)`` is applied too: one listed bond is one edge
    of the matching, seen once."""
    c, s = math.cos(g10), math.sin(g10)
    out = fields.copy()
    for (i, j, v) in bonds:
        # each sub-site appears in at most one bond of a stage: the
        # paired endpoints get the 2x2 rotation, the rest stay identity
        fi = fields[..., i, :].copy()
        fj = fields[..., j, :]
        src = np.roll(fj, shift=tuple(-int(x) for x in v),
                      axis=(0, 1, 2))
        src2 = np.roll(fi, shift=tuple(int(x) for x in v),
                       axis=(0, 1, 2))
        out[..., i, :] = c * fi - 1j * s * src
        out[..., j, :] = c * fj - 1j * s * src2
    return out


def front_speed(c, U, stages, g10=SC.G10, Lb=6, n_stages=12, seed=0,
                tol=1e-12):
    """The causal-radius test (V36 method) on a generalized cell: run the
    automaton (internal layer + bond stages, cycling ``stages``) from a
    state localized on one sub-site; the support radius is measured in
    **Planck cells** of the underlying lattice (cell coordinate
    ``n_j R_j + d_j`` for a block; ``R_j + s_j`` for a lattice cell with
    integer site offsets).  Returns the per-stage radii and the front
    speed (max radius / stages)."""
    rng = np.random.default_rng(seed)
    n = U.shape[0]
    shape = c.get("shape", (1, 1, 1))
    f = np.zeros((Lb, Lb, Lb, c["nsub"], n), dtype=complex)
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x0 = (Lb // 2, Lb // 2, Lb // 2)
    e0 = 0
    f[x0 + (e0,)] = v / np.linalg.norm(v)
    R = np.arange(Lb)
    axes = []
    for j in range(3):
        if c["kind"] == "block":
            d = np.array([c["digit"](e, j) for e in range(c["nsub"])])
            axes.append(shape[j] * R.reshape(
                [Lb if k == j else 1 for k in range(3)] + [1]) + d.reshape(
                    [1, 1, 1, c["nsub"]]))
        else:
            s = np.array([int(round(c["sites"][e][j]))
                          for e in range(c["nsub"])])
            axes.append(R.reshape(
                [Lb if k == j else 1 for k in range(3)] + [1]) + s.reshape(
                    [1, 1, 1, c["nsub"]]))
    cx, cy, cz = np.broadcast_arrays(*axes)
    if c["kind"] == "block":
        origin = tuple(shape[k] * x0[k] + c["digit"](e0, k) for k in range(3))
    else:
        origin = tuple(x0[k] + int(round(c["sites"][e0][k]))
                       for k in range(3))
    L = [shape[k] * Lb if c["kind"] == "block" else Lb for k in range(3)]

    def _md(d, Lk):
        return np.minimum(np.abs(d), Lk - np.abs(d))

    dist = (_md(cx - origin[0], L[0]) + _md(cy - origin[1], L[1])
            + _md(cz - origin[2], L[2]))
    radii = []
    outside = 0
    for k in range(1, n_stages + 1):
        f = np.einsum("...ai,ij->...aj", f, U)
        st = stages[(k - 1) % len(stages)]
        bonds = cell_stage_bonds(c, st) if c["kind"] == "block" \
            else lattice_stage_bonds(c, st)
        f = bond_layer_bonds(f, bonds, g10)
        nrm = np.sqrt(np.sum(np.abs(f) ** 2, axis=-1))
        sup = np.argwhere(nrm > tol)
        d = [int(dist[tuple(p)]) for p in sup]
        outside += int(sum(1 for r in d if r > k))
        radii.append(max(d) if d else 0)
    return {"radii": radii, "n_outside_lightcone": outside,
            "front_speed": radii[-1] / n_stages if radii else 0.0,
            "final_norm": float(np.linalg.norm(f.ravel()))}


# ---------------------------------------------------------------------------
# 4b. the motion axis: per-stage internal factors and bond angles
# ---------------------------------------------------------------------------


def stage_bonds(c, st):
    """The bonds of one stage in either representation: a block stage
    ``(j, phase)`` (resolved through ``cell_stage_bonds``) or an explicit
    bond list (a perfect matching, as produced by ``one_factorizations``)."""
    if isinstance(st, tuple):
        return cell_stage_bonds(c, st)
    return list(st)


def _broadcast(Us, thetas, n):
    """Normalize (Us, thetas) to per-stage lists of length ``n``: a single
    matrix / scalar is repeated, a list must have length ``n``."""
    Us_l = [Us] * n if isinstance(Us, np.ndarray) else list(Us)
    th_l = [thetas] * n if np.isscalar(thetas) else list(thetas)
    if len(Us_l) != n or len(th_l) != n:
        raise ValueError(f"per-stage lists must have length {n}: got "
                         f"{len(Us_l)}, {len(th_l)}")
    return Us_l, th_l


def bloch_stages(c, Us, thetas, stages, q):
    """The exact one-step Bloch matrix with a PER-STAGE internal matrix
    ``Us[t]`` and a PER-STAGE bond angle ``thetas[t]``:
    ``M = prod_t kron(S_t(thetas[t], q), Us[t])``.  Regression anchor:
    exact equality with ``bloch_cell`` when every entry is equal."""
    Us_l, th_l = _broadcast(Us, thetas, len(stages))
    M = np.eye(c["nsub"] * Us_l[0].shape[0], dtype=complex)
    for t, st in enumerate(stages):
        S = stage_matrix_bonds(stage_bonds(c, st), q, th_l[t], c["nsub"])
        M = np.kron(S, Us_l[t]) @ M
    return M


def step_stages(fields, Us, thetas, stages):
    """One real-space automaton step with per-stage internal factors and
    bond angles (the ``supercell.step`` semantics: internal layer then bond
    stage, per stage).  Anchor: equals ``supercell.step`` on the (2,2,2)
    cell with uniform entries."""
    Us_l, th_l = _broadcast(Us, thetas, len(stages))
    f = fields
    for t, st in enumerate(stages):
        f = SC.internal_layer(f, Us_l[t])
        if isinstance(st, tuple):
            j, phase = st
            f = SC.bond_layer(f, phase == 0, j, th_l[t])
        else:
            f = bond_layer_bonds(f, st, th_l[t])
    return f


def front_speed_stages(c, Us, thetas, stages, Lb=6, n_stages=12, seed=0,
                       tol=1e-12):
    """``front_speed`` generalized to per-stage internal factors and bond
    angles (same support-growth test, same Planck-cell radius convention);
    stages may be block ``(j, phase)`` pairs or explicit bond lists."""
    Us_l, th_l = _broadcast(Us, thetas, len(stages))
    rng = np.random.default_rng(seed)
    n = Us_l[0].shape[0]
    shape = c.get("shape", (1, 1, 1))
    f = np.zeros((Lb, Lb, Lb, c["nsub"], n), dtype=complex)
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x0 = (Lb // 2, Lb // 2, Lb // 2)
    e0 = 0
    f[x0 + (e0,)] = v / np.linalg.norm(v)
    R = np.arange(Lb)
    axes = []
    for j in range(3):
        if c["kind"] == "block":
            d = np.array([c["digit"](e, j) for e in range(c["nsub"])])
            axes.append(shape[j] * R.reshape(
                [Lb if k == j else 1 for k in range(3)] + [1]) + d.reshape(
                    [1, 1, 1, c["nsub"]]))
        else:
            s = np.array([int(round(c["sites"][e][j]))
                          for e in range(c["nsub"])])
            axes.append(R.reshape(
                [Lb if k == j else 1 for k in range(3)] + [1]) + s.reshape(
                    [1, 1, 1, c["nsub"]]))
    cx, cy, cz = np.broadcast_arrays(*axes)
    if c["kind"] == "block":
        origin = tuple(shape[k] * x0[k] + c["digit"](e0, k) for k in range(3))
    else:
        origin = tuple(x0[k] + int(round(c["sites"][e0][k]))
                       for k in range(3))
    L = [shape[k] * Lb if c["kind"] == "block" else Lb for k in range(3)]

    def _md(d, Lk):
        return np.minimum(np.abs(d), Lk - np.abs(d))

    dist = (_md(cx - origin[0], L[0]) + _md(cy - origin[1], L[1])
            + _md(cz - origin[2], L[2]))
    radii = []
    outside = 0
    for k in range(1, n_stages + 1):
        f = SC.internal_layer(f, Us_l[(k - 1) % len(stages)])
        st = stages[(k - 1) % len(stages)]
        f = bond_layer_bonds(f, stage_bonds(c, st),
                             th_l[(k - 1) % len(stages)])
        nrm = np.sqrt(np.sum(np.abs(f) ** 2, axis=-1))
        sup = np.argwhere(nrm > tol)
        d = [int(dist[tuple(p)]) for p in sup]
        outside += int(sum(1 for r in d if r > k))
        radii.append(max(d) if d else 0)
    return {"radii": radii, "n_outside_lightcone": outside,
            "front_speed": radii[-1] / n_stages if radii else 0.0,
            "final_norm": float(np.linalg.norm(f.ravel()))}


# ---------------------------------------------------------------------------
# 5. stars and the node census
# ---------------------------------------------------------------------------

_OH = [(perm, signs) for perm in itertools.permutations(range(3))
       for signs in itertools.product((1, -1), repeat=3)]


def oh_star(q):
    """The O_h orbit of a Cartesian momentum, folded into [-pi, pi)^3
    (the ``supercell.star`` convention: 2pi-periodic cubic units)."""
    def norm(qq):
        return tuple(round(((v + math.pi) % (2 * math.pi)) - math.pi, 9)
                     for v in qq)
    return sorted({norm(tuple(signs[i] * q[perm[i]] for i in range(3)))
                   for perm, signs in _OH})


def fcc_star(q):
    """The O_h orbit of a Cartesian momentum (units 2pi/a) modulo the FCC
    reciprocal lattice (the BCC lattice: integer vectors of even component
    sum), canonicalized in the reciprocal basis."""
    B = np.array([(-1, 1, 1), (1, -1, 1), (1, 1, -1)], float)  # 2pi/a basis
    Binv = np.linalg.inv(B)

    def canon(qq):
        c = Binv @ np.array(qq, float)
        c = (c + 1e-9) % 1.0
        return tuple(np.round(B @ c, 9))

    return sorted({canon(tuple(signs[i] * q[perm[i]] for i in range(3)))
                   for perm, signs in _OH})


def sphere_dirs(n=6):
    """The cone directions of the V115 scan (the same grid)."""
    dirs = []
    for i in range(1, n):
        for j in range(n):
            for k in range(n):
                v = np.array([i - n / 2, j - n / 2, k - n / 2], float)
                if np.linalg.norm(v) > 0:
                    dirs.append(v / np.linalg.norm(v))
    return dirs


def gap_at(M):
    """min |eigenphase| of the one-step matrix (exact, small matrices)."""
    return float(np.min(np.abs(np.sort(np.angle(np.linalg.eigvals(M))))))


def cone_anisotropy(bloch_fn, q_star, dirs=None, dq=1e-3):
    """The V115 anisotropy metric at a gapless point: the cone speed
    ``min|omega| / dq`` over the sphere directions,
    ``A = (v_max - v_min) / (v_max + v_min)``."""
    dirs = sphere_dirs(6) if dirs is None else dirs
    vs = []
    for nd in dirs:
        qq = tuple(q + dq * n for q, n in zip(q_star, nd))
        vs.append(gap_at(bloch_fn(qq)) / dq)
    vmax, vmin = max(vs), min(vs)
    if vmax + vmin <= 0:
        return 1.0, vmax, vmin
    return (vmax - vmin) / (vmax + vmin), vmax, vmin


def chirality_at(bloch_fn, q_star, dq=1e-4):
    """The V63 node chirality: the sign of the product of the axis slopes
    of the lowest branch at the node (nonzero axes only)."""
    sl = []
    for ax in range(3):
        if abs(q_star[ax]) < 1e-12:
            continue
        qa = list(q_star); qa[ax] -= dq
        qb = list(q_star); qb[ax] += dq
        pa = np.angle(np.linalg.eigvals(bloch_fn(tuple(qa))))
        pb = np.angle(np.linalg.eigvals(bloch_fn(tuple(qb))))
        sl.append(float((pb[np.argmin(np.abs(pb))]
                         - pa[np.argmin(np.abs(pa))]) / (2 * dq)))
    if not sl:
        return 0.0, sl
    return float(np.sign(np.prod(sl))), sl


# ---------------------------------------------------------------------------
# 6. the (1+3) frame of a signature (V150)
# ---------------------------------------------------------------------------


def frame_residual(q, time_idx, space_idx):
    """The Minkowski-subalgebra residual of one frame choice, in the
    ``spacetime.causal_generators`` convention: time ``T = h_t``, space
    ``S_a = i h_{s_a}`` (so the baseline frame squares to ``T^2 = +I``,
    ``S_a^2 = -I``).  The rotations ``J_a = S_b S_c`` (cyclic) and boosts
    ``K_a = T S_a`` must satisfy the so(1,3) form

        [J_a, J_b] = 2 eps_abc J_c,   [J_a, K_b] = 2 eps_abc K_c,
        [K_a, K_b] = -2 eps_abc J_c.

    The actual structure constants carry the signature factors
    (``[J_a,J_b] = 2 eps_{s_c} J_c`` etc.); the residual is the max
    deviation from the target (0 iff the frame is exactly (1+3))."""
    h = signature_generators(q)
    T = h[time_idx]
    S = [1j * h[k] for k in space_idx]
    J = [S[1] @ S[2], S[2] @ S[0], S[0] @ S[1]]
    K = [T @ S[0], T @ S[1], T @ S[2]]
    basis = J + K
    actual = np.zeros((6, 6, 6))
    for a in range(6):
        for b in range(6):
            C = basis[a] @ basis[b] - basis[b] @ basis[a]
            for c in range(6):
                actual[a, b, c] = float(np.real(
                    np.trace(C @ basis[c].conj().T)
                    / np.trace(basis[c] @ basis[c].conj().T)))
    target = np.zeros((6, 6, 6))
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        target[a, b, c] = 2.0
        target[b, a, c] = -2.0
        target[a + 3, b + 3, c] = -2.0
        target[b + 3, a + 3, c] = 2.0
        target[a, b + 3, c + 3] = 2.0
        target[b + 3, a, c + 3] = -2.0
        target[b, a + 3, c + 3] = -2.0
        target[a + 3, b, c + 3] = 2.0
    return float(np.abs(actual - target).max())


def main():
    import sys
    # the regression anchor: the generalized stage matrix equals
    # supercell.stage_matrix on the (2,2,2) cell, every axis, both phases.
    c = cell((2, 2, 2))
    err = 0.0
    for j in range(3):
        for p in (0, 1):
            err = max(err, float(np.abs(
                stage_matrix_cell(c, j, p, 0.3)
                - SC.stage_matrix(p == 0, j, 0.3)).max()))
    print(f"altcell anchor: stage_matrix_cell vs supercell.stage_matrix "
          f"max err = {err:.2e}")
    lc = lattice_cell("diamond")
    fac, sk = one_factorizations(lc)
    print(f"diamond: {len(lc['bonds'])} bonds, {len(fac)} factorizations, "
          f"{sk} skipped")
    lc = lattice_cell("fcc")
    fac, sk = one_factorizations(lc, cap=50)
    print(f"fcc: {len(lc['bonds'])} bonds, {len(fac)} factorizations "
          f"(cap 50), {sk} skipped")
    # the motion anchors: bloch_stages / step_stages / front_speed_stages
    # with uniform entries equal the implemented motion exactly.
    _B, _U32, U16 = sectors_sig(0)
    st = stages_baseline()
    err_bloch = 0.0
    for q in ((0.3, -1.1, 2.0), (math.pi / 2, math.pi / 2, 0.0)):
        err_bloch = max(err_bloch, float(np.abs(
            bloch_stages(c, [U16] * 6, [SC.G10] * 6, st, q)
            - bloch_cell(c, U16, SC.G10, st, q)).max()))
    print(f"altcell anchor: bloch_stages vs bloch_cell "
          f"max err = {err_bloch:.2e}")
    rng = np.random.default_rng(1)
    f0 = (rng.standard_normal((4, 4, 4, 8, 16))
          + 1j * rng.standard_normal((4, 4, 4, 8, 16)))
    err_step = float(np.abs(
        step_stages(f0, [U16] * 6, [SC.G10] * 6, st)
        - SC.step(f0, U16, SC.G10, SC.STAGES)).max())
    print(f"altcell anchor: step_stages vs supercell.step "
          f"max err = {err_step:.2e}")
    fs1 = front_speed(c, U16, st, SC.G10, Lb=6, n_stages=12)
    fs2 = front_speed_stages(c, [U16] * 6, [SC.G10] * 6, st,
                             Lb=6, n_stages=12)
    err_fs = abs(fs1["front_speed"] - fs2["front_speed"]) + \
        abs(fs1["n_outside_lightcone"] - fs2["n_outside_lightcone"])
    print(f"altcell anchor: front_speed_stages vs front_speed "
          f"diff = {err_fs:.2e}")
    sys.exit(0 if err < 1e-14 and err_bloch == 0.0
             and err_step < 1e-14 and err_fs == 0.0 else 1)


if __name__ == "__main__":
    main()
