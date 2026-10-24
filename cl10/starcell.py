"""The star module: the committed automaton step with the icosahedral star
substituted for the cubic shift module.

The implemented construction (cl10/altcell.py, cl10/supercell.py) is a
staggered matching step on a mixed-radix block cell: each stage is a perfect
matching of the sub-sites, a bond (i, j, v) applies the 2x2 rotation

    [[c, -i s e^{+i q.v}], [-i s e^{-i q.v}, c]],      c = cos g, s = sin g

and the one-step matrix is the product of `kron(S_stage, U)` over the stages,
with U the committed internal map on the 16 sector.  The only thing that makes
that construction cubic is the list of bond displacements: the block cell hops
along the three axes e_1, e_2, e_3, and the Bloch phase is exp(i q . v) with v
the image label.

V212 and V213 close that construction off: the one locus of its Brillouin zone
where symmetry forces isotropy has a quadratic touching, every locus that
touches linearly has its anisotropy forced by symmetry, and the tangent
anisotropy is marginal (exponent -0.0001 on a radius ladder, 0.0000 on a
properly rescaled block ladder), so no coarse-graining inside it produces a
rotational infrared.  The alternative the root's own Clifford algebra points at
is the icosahedral star (cl10/continuum.py: `icosahedral_star`, `star_module`):
twelve signed directions, six axes, six equiangular lines at 1/sqrt 5, an
isotropic fourth moment where the cubic axes have an anisotropic one, and a
Z-module of rank 6 -- which is why the module is handled on its own footing
here rather than as a 3D lattice cell.

The 6D footing is exact, not an approximation: a rank-6 module is the image of
Z^6 under the physical projection pi carrying the i-th basis vector to the i-th
module generator, and a step that hops along the module generators is a
staggered matching step on the periodic Z^6 cell with Bloch phases
exp(i q . pi(e_i)).  That is the same code shape as `altcell.bloch_cell` with
the integer image label replaced by the physical image of the module axis.

Everything internal is reused unchanged: the internal map is the committed
`supercell.sectors()[2]`, the clock comes from cl10/config.yaml, and the only
new input is the module itself.  What this module does not claim is that the
star construction is the physical one: V214, V215 and V216 measure what it
gives -- its invariants, its front speed, its gapless locus, and whether its
tangent algebra is Clifford, the property V187 proved mutually exclusive with
gaplessness inside the cubic module.

The step matrices are built and diagonalised through `cl10.engine`, so a
1024-dimensional module step runs on the GPU when the engine commits to a GPU
core (the module cell is the smallest one that carries six periodic axes) and
on the host otherwise; the host-only helpers here are bookkeeping.
"""
import math

import numpy as np

from cl10 import continuum as K
from cl10 import engine as E

NDIM = 6


def module_cell(shape=(2,) * NDIM):
    """The mixed-radix periodic cell of the module: `len(shape)` module axes,
    `nsub = prod(shape)` sub-sites, `digit(e, i)` the i-th mixed-radix digit of
    e and `site_index(digits)` its inverse -- the same bookkeeping as
    `altcell.cell`, in six directions instead of three."""
    shape = tuple(int(n) for n in shape)
    radix = [int(np.prod(shape[:i])) for i in range(len(shape))]

    def digit(e, i):
        return (e // radix[i]) % shape[i]

    def site_index(digits):
        return sum(int(d) * int(r) for d, r in zip(digits, radix))

    return {"kind": "module", "shape": shape, "ndim": len(shape),
            "nsub": int(np.prod(shape)), "radix": radix, "digit": digit,
            "site_index": site_index,
            "sites": [tuple((e // r) % n for r, n in zip(radix, shape))
                      for e in range(int(np.prod(shape)))]}


def module_bonds(c, i, phase):
    """One stage of the module cell: every sub-site whose i-th digit equals
    `phase` hops to digit `phase + 1 mod n_i`.  A wrap is an inter-cell hop,
    and it is the only one that carries a Bloch phase -- the physical image of
    the i-th module axis (see `stage_matrix`)."""
    shape, nsub = c["shape"], c["nsub"]
    n = shape[i]
    bonds = []
    for e in range(nsub):
        if c["digit"](e, i) != phase:
            continue
        digits = list(c["sites"][e])
        nxt = list(digits)
        nxt[i] = (phase + 1) % n
        wraps = bool(phase + 1 == n)
        bonds.append((e, c["site_index"](nxt), i, wraps))
    return bonds


def module_stages(shape=(2,) * NDIM, order="phase-major"):
    """The stagger of the module cell.  `phase-major` alternates intra- and
    inter-cell stages per axis the way `supercell.STAGES` does for the cubic
    block; `block-major` runs all intra-cell stages first.  The choice is part
    of the construction, not a consequence of it, so V214 measures both."""
    shape = tuple(int(n) for n in shape)
    if order == "phase-major":
        return [(i, phase) for phase in range(max(shape))
                for i, n in enumerate(shape) if phase < n]
    if order == "block-major":
        return [(i, phase) for i, n in enumerate(shape) for phase in range(n)]
    raise ValueError("unknown stage order " + str(order))


def stage_matrix(bonds, q, axes, g10, nsub, bond_g10=None):
    """The matching matrix of one module stage: identity on the unpaired
    sub-sites, and on a bond (i, j, axis, wraps) the 2x2 rotation

        [[c, -i s e^{+i q.v}], [-i s e^{-i q.v}, c]]

    with v the physical image of that module axis -- the `altcell` convention
    with the integer image label replaced by the module generator.
    ``bond_g10`` is the decoupled bond angle beta (V239/V244); ``None`` is
    the published step, bit-identically."""
    b = g10 if bond_g10 is None else bond_g10
    c, s = math.cos(b), math.sin(b)
    S = np.eye(nsub, dtype=complex)
    for (a, b_i, axis, wraps) in bonds:
        ph = np.exp(1j * float(np.dot(q, axes[axis]))) if wraps else 1.0 + 0j
        S[a, a] = c
        S[b_i, b_i] = c
        S[a, b_i] = -1j * s * ph
        S[b_i, a] = -1j * s * np.conj(ph)
    return S


def bloch_star(cell, U, g10, stages, q, axes, bond_g10=None):
    """The one-step Bloch matrix of the module: `Pi kron(S_stage, U)` over the
    stagger, with U the committed internal map.  ``bond_g10`` decouples the
    bond angle beta from the internal angle g10 (V244); ``None`` is the
    published step.  The product runs on the device when the engine is on a
    GPU core, which is where the 1024-dimensional module step belongs; `U`
    is moved once and the stage matrices ride along."""
    U = E.to_device(U)
    M = E.eye(cell["nsub"] * U.shape[0], complex)
    for (i, phase) in stages:
        S = stage_matrix(module_bonds(cell, i, phase), q, axes, g10,
                         cell["nsub"], bond_g10)
        M = E.matmul(E.kron(E.to_device(S), U), M)
    return M


def eigenphases(M):
    """The eigenphases of a step matrix, from the routed eigensolver."""
    return np.angle(E.asnumpy(E.eigvals(M)))


def gap_star(cell, U, g10, stages, q, axes, bond_g10=None):
    """min |eigenphase| of the module's one-step matrix."""
    return float(np.min(np.abs(eigenphases(
        bloch_star(cell, U, g10, stages, q, axes, bond_g10)))))


def parity(cell):
    """The reflection of the module cell, digit by digit (d -> -d mod n): the
    operator the mirror Z2 has to be compared against."""
    P = np.zeros((cell["nsub"], cell["nsub"]))
    for e, digits in enumerate(cell["sites"]):
        refl = [(-d) % n for d, n in zip(digits, cell["shape"])]
        P[cell["site_index"](refl), e] = 1.0
    return P


def lr_radius(axes):
    """The Lieb-Robinson radius per step: the longest module generator, since a
    step hops along one of them (in units where the icosahedron's vertex
    directions are the generators)."""
    return float(max(np.linalg.norm(np.asarray(a, float)) for a in axes))


def word_displacement(coefs, axes):
    """The physical displacement of a word in module steps: the number of steps
    is sum |coefs|, the reach is the physical projection."""
    coefs = np.asarray(coefs, float)
    A = np.array([np.asarray(a, float) for a in axes])
    d = coefs @ A
    return {"steps": int(np.abs(coefs).sum()),
            "displacement": [float(x) for x in d],
            "norm": float(np.linalg.norm(d))}


def star_axes(shape=None):
    """The module generators: the independent axes of the icosahedral star, as
    physical displacements.  This is the shift module."""
    return K.module_generators()["vectors"]
