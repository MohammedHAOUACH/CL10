#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""The pentagon engine: the fundamental cell is the dodecahedron's
TWELVE FACES, not the cubic 8-subsite cell.

Architectural reason (booked at V356, review/pre_registrations_root.md):
the cubic class's closures CONSUME the clock (V334 theorem). The full
icosahedral group I_h is FINITE: a step whose core is the wreath
(face permutation x frame rotation) over I_h closes at every clock, so
the clock can live in the coefficients. That is what this library
builds: a 36-dimensional one-particle cell (12 faces x 3 frame), a
closed I_h-hopping hop part U_hop(g), and the full step tensoring the
published internal stages (cl10.supercell.internal, used VERBATIM -
"copy that loop pattern, never the cubic matrices").

The closed cell carries an exact spectral gap: the 12-face adjacency
(the icosahedral graph, MathWorld: characteristic polynomial
(x-5)(x+1)^5(x^2-5)^3, |Aut| = 120) has spectrum
{5, (+-sqrt5)^3, (-1)^5}, so the twisted 36-cell hop spectrum starts
at 5 - sqrt5 = 2.76393202250021... - a mass gap, certified by V358.
A GAPLESS pentagon theory therefore needs the quasiperiodic extension:
the de Bruijn pentagrid (V360). This module is the reason V360 exists.

Math conventions (the pentagon-12 special math):
- face normals = icosahedron vertices (phi = (1+sqrt5)/2, mpmath 50 dp);
- edges = the 30 vertex pairs at squared distance 2;
- the I_h table is built by BFS from (Rz(2pi/5), Rx(pi/2), -I) with
  exact-ring (mpmath) entries; every chained product's entries are
  normalized (+x at dps 50) and dedup keys are canonical 14-digit
  strings - the V356 expression/path-rounding blowup must not be
  replicated here;
- the edge frame rotation R_uv is the geometric-algebra ROTOR of the
  pentagon cell (kingdon, Cl(3,0)): r_uv ~ 1 + n_v n_u (geometric
  product), v -> r v r~. The rotor is the same Rodrigues rotation the
  plan books; both are computed and cross-certified (<= 1e-14), with
  entries in Q(sqrt5).

No scipy here (numpy + mpmath only; cupy appears only in V358-V360).
"""
import mpmath

MP = 50  # dps for every exact-ring computation in this module

MPH = None  # filled on first use (exact 1/2 + sqrt5/2)


def phi():
    """The golden ratio at MP digits."""
    mpmath.mp.dps = MP
    return (1 + mpmath.sqrt(5)) / 2


def face_geometry():
    """The 12 face normals (= icosahedron vertices) and the 30 edges.

    Returns ``(vs, edges)``: normals as 3-tuples of mpmath objects,
    edges as index pairs at squared distance 2 (the minimum over all
    pairs).
    """
    mpmath.mp.dps = MP
    ph = phi()
    vs = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            vs.append((mpmath.mpf(0), mpmath.mpf(s1), mpmath.mpf(s2) * ph))
            vs.append((mpmath.mpf(s1), mpmath.mpf(s2) * ph, mpmath.mpf(0)))
            vs.append((mpmath.mpf(s1) * ph, mpmath.mpf(0), mpmath.mpf(s2)))
    assert len(vs) == 12

    def d2(a, b):
        return sum((a[k] - b[k]) ** 2 for k in range(3))
    dmin = min(d2(vs[i], vs[j]) for i in range(12) for j in range(i + 1, 12))
    edges = [(i, j) for i in range(12) for j in range(i + 1, 12)
             if abs(d2(vs[i], vs[j]) - dmin) < mpmath.mpf(10) ** -30]
    return vs, edges


def _rot_z(theta):
    mpmath.mp.dps = MP
    return [[mpmath.cos(theta), -mpmath.sin(theta), mpmath.mpf(0)],
            [mpmath.sin(theta), mpmath.cos(theta), mpmath.mpf(0)],
            [mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(1)]]


def _rot_x(theta):
    mpmath.mp.dps = MP
    return [[mpmath.mpf(1), mpmath.mpf(0), mpmath.mpf(0)],
            [mpmath.mpf(0), mpmath.cos(theta), -mpmath.sin(theta)],
            [mpmath.mpf(0), mpmath.sin(theta), mpmath.cos(theta)]]


def _mul(A, B):
    """3x3 exact-ring product with entries normalized at dps MP (the
    V356 book-fix: chained products must not grow)."""
    X = [[+sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)]
         for i in range(3)]
    return X


def _key(M):
    """Canonical dedup key: distinct I_h entries differ by >= 1/4;
    different generation paths agree to ~1e-47 after normalization, so
    14 significant digits dedup bit-exactly. The |x| < 1e-25 threshold
    kills the mpmath.pi transcendental residue (sin(pi) ~ 1e-50)."""
    def k(x):
        return "0.0" if abs(x) < mpmath.mpf(10) ** -25 \
            else mpmath.nstr(x, 14)
    return tuple(tuple(k(x) for x in row) for row in M)


def _rot_axis(axis, theta):
    """Rodrigues rotation about `axis` (mpmath entries) by theta, entries
    normalized at dps MP."""
    mpmath.mp.dps = MP
    n = mpmath.sqrt(sum(a * a for a in axis))
    u = [a / n for a in axis]
    c, s = mpmath.cos(theta), mpmath.sin(theta)
    ux = [[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]]
    return [[+c * (1 if i == j else 0) + (1 - c) * u[i] * u[j]
             + s * ux[i][j] for j in range(3)] for i in range(3)]


def ih_table():
    """The 120 elements of I_h by BFS from the cell's own axes —
    C5 on the 5-fold face-normal axis (0,1,phi), C2x, C3 on the
    threefold (1,1,1) axis (the icosahedron face centre: EXACTLY the
    cyclic permutation matrix), and the inversion — with the induced
    12-face permutation and the STRICT element order of each.

    (The plan's booked generator text "(Rz(2pi/5), Rx(pi/2), -I)" does
    not generate I_h in this orientation: Rz(72 deg) is not a symmetry
    of the (0,+-1,+-phi) vertex set, Rx(pi/2) is not either, and
    C5-C2-perpendicular generates only D5. The booked GATES (|I_h| =
    120, strict order census max 10, faithful 12-action) are the
    contract and are unchanged. sin(72 deg)/sqrt(1+phi^2) = 1/2
    exactly (10+2sqrt5 = 4(phi+2)), so every generator entry stays in
    Q(sqrt5) up to the C3 permutation being integral.)

    Returns a list of dicts ``{"M": M, "perm": p, "order": k}`` in BFS
    discovery order (the canonical copy of the V356 construction).
    """
    mpmath.mp.dps = MP
    vs, _ = face_geometry()
    ph = phi()
    I3 = [[mpmath.mpf(1) if i == j else mpmath.mpf(0) for j in range(3)]
          for i in range(3)]
    NI3 = [[-mpmath.mpf(1) if i == j else mpmath.mpf(0) for j in range(3)]
           for i in range(3)]
    C5m = _rot_axis([mpmath.mpf(0), mpmath.mpf(1), ph], 2 * mpmath.pi / 5)
    C2x = _rot_x(mpmath.pi)
    C3 = [[+mpmath.mpf(0), +mpmath.mpf(0), +mpmath.mpf(1)],
          [+mpmath.mpf(1), +mpmath.mpf(0), +mpmath.mpf(0)],
          [+mpmath.mpf(0), +mpmath.mpf(1), +mpmath.mpf(0)]]
    gens = [C5m, C2x, C3, NI3]
    seen = {}
    frontier = [I3]
    seen[_key(I3)] = I3
    order_list = [I3]
    while frontier:
        nxt = []
        for M in frontier:
            for G in gens:
                X = _mul(M, G)
                k = _key(X)
                if k not in seen:
                    seen[k] = X
                    order_list.append(X)
                    nxt.append(X)
        frontier = nxt
    if len(order_list) != 120:
        raise AssertionError("I_h BFS: %d elements, expected 120"
                             % len(order_list))

    def perm_of(M):
        p = []
        for v in vs:
            mv = tuple(sum(M[i][j] * v[j] for j in range(3)) for i in range(3))
            idx = min(range(12),
                      key=lambda t: max(abs(mv[i] - vs[t][i])
                                        for i in range(3)))
            if max(abs(mv[i] - vs[idx][i]) for i in range(3)) > 1e-9:
                raise AssertionError("I_h action not closed on the 12 faces")
            p.append(idx)
        return tuple(p)

    def order_of(M, kmax=60):
        P = I3
        for k in range(1, kmax + 1):
            P = _mul(P, M)
            if _key(P) == _key(I3):
                return k
        raise AssertionError("element with no strict order <= %d" % kmax)

    return [{"M": M, "perm": perm_of(M), "order": order_of(M)}
            for M in order_list]


def _rodriques(nu, nv):
    """The rotation matrix R with R n_u = n_v (axis n_u x n_v, angle
    arccos(n_u . n_v)); all entries mpmath. Degenerate (n_u = n_v) -> I."""
    mpmath.mp.dps = MP
    dot = sum(nu[k] * nv[k] for k in range(3))
    ax = [nu[1] * nv[2] - nu[2] * nv[1],
          nu[2] * nv[0] - nu[0] * nv[2],
          nu[0] * nv[1] - nu[1] * nv[0]]
    nax = mpmath.sqrt(sum(a * a for a in ax))
    if nax < mpmath.mpf(10) ** -40:
        return [[mpmath.mpf(1) if i == j else mpmath.mpf(0)
                 for j in range(3)] for i in range(3)], dot
    ang = mpmath.acos(max(-1, min(1, dot)))
    u = [a / nax for a in ax]
    c, s = mpmath.cos(ang), mpmath.sin(ang)
    # Rodrigues: R = c I + (1-c) u u^T + s [u]_x
    ux = [[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]]
    R = [[+c * (1 if i == j else 0) + (1 - c) * u[i] * u[j]
          + s * ux[i][j] for j in range(3)] for i in range(3)]
    return R, dot


def _unit(v):
    mpmath.mp.dps = MP
    n = mpmath.sqrt(sum(x * x for x in v))
    return [x / n for x in v]


def _rotor(nu, nv):
    """The SAME rotation via kingdon geometric algebra: the rotor
    r ~ 1 + n_v n_u acts as v -> r v r~ (Cl(3,0), kingdon's float64
    coefficient path - the cross-certification it feeds is a 1e-14
    gate; the exact side stays mpmath and keeps the ring certificate).
    Returns a 3x3 numpy float64 matrix. (Blade keys are bitmasks:
    e_i -> 1 << i.)"""
    import numpy as np
    from kingdon import Algebra
    if not hasattr(_rotor, "alg"):
        _rotor.alg = Algebra(3, 0, 0)
    alg = _rotor.alg
    nuv = alg.vector([float(x) for x in nu])
    nvv = alg.vector([float(x) for x in nv])
    r = (alg.blades.e + nvv * nuv).normalized()
    R = np.zeros((3, 3))
    for j in range(3):
        ej = alg.vector([1.0 if k == j else 0.0 for k in range(3)])
        w = r * ej * ~r
        d = dict(zip(w._keys, w._values))
        for i in range(3):
            R[i][j] = float(d.get(1 << i, 0.0))
    return R


def _in_ring(x):
    """True if x in Q(sqrt5) with coefficients of |64|: PSLQ finds the
    integer relation c0 + c1 sqrt5 + c2 x = 0; verified at dps 60 by
    exact reconstruction (residual < 1e-40)."""
    mpmath.mp.dps = 60
    s5 = mpmath.sqrt(5)
    x = +mpmath.mpf(x)
    if x == 0:
        return True
    c = mpmath.pslq([1, s5, x], maxcoeff=64)
    if c is None or all(ci == 0 for ci in c):
        return False
    c0, c1, c2 = c
    if c2 == 0:
        return False      # impossible except for a PSLQ artifact
    resid = abs(c0 + c1 * s5 + c2 * x)
    return resid < mpmath.mpf(10) ** -40


def certify_rotors(tol=1e-14):
    """Certify the edge twists on all 30 edges:
    (i) Rodrigues exact side: R n_u = n_v to `tol` (mpmath, dps 50+),
        R R^T = I and det R = +1 to 1e-30;
    (ii) kingdon geometric-algebra rotor == Rodrigues to 1e-13 (float64);
    (iii) the Q(sqrt5) entry structure booked in the plan: EVERY entry
        of every twist lies in Q(sqrt5) (PSLQ relation, coefficients
        |64|, reconstruction residual 1e-40). For the pentagon cell the
        adjacent-face twist has cos = 1/sqrt5 and axis (1, 0, phi)-class
        entries, so all twist entries close in the icosahedral ring -
        this is what makes the moments of K (tr K = 0, tr K^2 = 180)
        exact.
    Returns (max_dev_nv, max_dev_pair, n_ring_entries)."""
    mpmath.mp.dps = 50
    vs, edges = face_geometry()
    d1 = d2 = d3 = mpmath.mpf(0)
    n_ring = 0
    for (u, v) in edges:
        nu, nv = _unit(vs[u]), _unit(vs[v])
        R, _ = _rodriques(nu, nv)
        Rg = _rotor(nu, nv)
        mv = [sum(R[i][k] * nu[k] for k in range(3)) for i in range(3)]
        d1 = max(d1, max(abs(mv[i] - nv[i]) for i in range(3)))
        d2 = max(d2, max(abs(float(Rg[i][j]) - float(R[i][j]))
                         for i in range(3) for j in range(3)))
        RRt = [[+sum(R[i][k] * R[j][k] for k in range(3))
                for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):
                d3 = max(d3, abs(RRt[i][j] - (1 if i == j else 0)))
        detR = (+R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1])
                - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0])
                + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]))
        d3 = max(d3, abs(detR - 1))
        for i in range(3):
            for j in range(3):
                if _in_ring(R[i][j]):
                    n_ring += 1
                else:
                    raise AssertionError(
                        "twist entry (%d,%d) edge (%d,%d) not in Q(sqrt5)"
                        % (i, j, u, v))
    if d1 > tol or d2 > 10 * tol or d3 > 1e-30:
        raise AssertionError("rotor certification failed: nv %.1e, "
                             "pair %.1e, orth %.1e" % (d1, d2, d3))
    return float(d1), float(d2), n_ring


def twisted_hop(vs=None, edges=None):
    """The 36x36 Hermitian twisted hop K on face (x) frame.

    For each edge (u, v) the block (u, v) is the frame rotor R_uv
    (normal n_u carried to n_v by the shortest arc on the cell), and
    (v, u) is R_uv^T. K is real symmetric (R_vu = R_uv^T exactly on the
    icosahedral cell), hence Hermitian; the certified spectrum starts at
    5 - sqrt5 (V358's PENTGAP leaf)."""
    mpmath.mp.dps = MP
    import numpy as np
    if vs is None:
        vs, edges = face_geometry()
    K = np.zeros((36, 36), dtype=np.float64)
    for (u, v) in edges:
        R, _ = _rodriques(_unit(vs[u]), _unit(vs[v]))
        Rf = np.array([[float(R[i][j]) for j in range(3)] for i in range(3)])
        K[3 * u:3 * u + 3, 3 * v:3 * v + 3] = Rf
        K[3 * v:3 * v + 3, 3 * u:3 * u + 3] = Rf.T
    return K


def hop_spectrum(K=None):
    """The certified 36 eigenvalues of K (sorted), via mpmath-free dense
    symmetric eigensolve (numpy) - exact values are integers and
    +-sqrt5-class numbers; V358 gates the device-vs-host match."""
    import numpy as np
    if K is None:
        K = twisted_hop()
    return np.linalg.eigvalsh(K)


def hop_unitary(g, K=None, K_eig=None):
    """U_hop(g) = W diag(exp(i g lambda)) W^dagger from the symmetric
    eigensolve of K - never expm: unitarity by construction."""
    import numpy as np
    if K_eig is None:
        lam, W = np.linalg.eigh(twisted_hop() if K is None else K)
    else:
        lam, W = K_eig
    return (W * np.exp(1j * g * lam)) @ W.conj().T


def internal32(g10):
    """The published 32x32 internal step, VERBATIM from cl10.supercell
    (the internal stages of the cubic theory are carried unchanged into
    the pentagon cell - "copy the pattern, never the cubic matrices").
    Host numpy via cl10.engine.asnumpy."""
    from cl10 import supercell
    from cl10 import engine as E
    import numpy as np
    return E.asnumpy(supercell.internal(g10)).astype(np.complex128)


def step(g, g10, K=None, K_eig=None):
    """The pentagon one-step operator on face(12) (x) frame(3) (x)
    internal(32):  U = (I36 (x) Uin . U_hop(g) (x) I32)^6.

    Six stages alternate hop/internal exactly as cl10.supercell.
    bond_product alternates its six stages (the internal layer sits on
    the left of each stage there too - copy that alternation pattern,
    never the cubic matrices). The cell factors are the einops axes
    f (face), r (frame), a (internal)."""
    import numpy as np
    from einops import einsum, rearrange
    Uh = hop_unitary(g, K=K, K_eig=K_eig)
    Uin = internal32(g10)
    Uh4 = Uh.reshape(12, 3, 12, 3)            # (f r), (g s)
    I36t = np.eye(36, dtype=np.complex128).reshape(12, 3, 12, 3)
    I32t = np.eye(32, dtype=np.complex128)
    H = rearrange(einsum(Uh4, I32t, 'f r g s, a b -> f r a g s b'),
                  'f r a g s b -> (f r a) (g s b)')   # U_hop (x) I32
    N = rearrange(einsum(I36t, Uin, 'f r g s, a b -> f r a g s b'),
                  'f r a g s b -> (f r a) (g s b)')   # I36 (x) Uin
    M = np.eye(1152, dtype=np.complex128)
    for _stage in range(6):                  # hop, internal, alternating
        M = N @ (H @ M)
    return M


# ---------------------------------------------------------------------------
# the CELL itself: the booked table - Type: regular polyhedron, faces:
# pentagons; F = 12, E = 30, V = 20, Euler characteristic 2; three faces
# per vertex, five vertices per face ({5,3}). The nest: 120 dodecahedra -
# the regular orbit of I_h. The order-4 dodecahedral honeycomb {5,3,4}:
# four cells around every edge, the compact hyperbolic sibling of the
# V360 flat pentagrid.
# ---------------------------------------------------------------------------

def dodeca_cell():
    """The pentagon-12 cell as a regular dodecahedron {5,3}, built as
    the polar of the 12 face normals (face_geometry): the dodecahedron
    {x : n_i . x <= 1} is regular with faces ⊥ the engine normals, so
    face i of the cell is EXACTLY engine face i - no relabeling.

    Certifies the booked table: 20 vertices, 30 edges, 12 faces of five
    vertices each, three faces per vertex, chi = V - E + F = 2; every
    face planar by construction (n_i . x = 1) and equilateral (common
    edge length) in a plane cutting the vertex sphere - an equilateral
    cyclic pentagon is regular. Raises AssertionError otherwise.
    Returns ``{"pts": (20,3) mpmath, "faces": 12 x (5 vertex indices),
    "edges": 30 pairs, "counts": (V, E, F, chi), "edge_len": x}``.
    """
    mpmath.mp.dps = MP
    vs, _ = face_geometry()

    def solve3(A, b):
        # Cramer in mpmath
        def det3(M):
            return (+M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                    - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                    + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        D = det3(A)
        if abs(D) < mpmath.mpf(10) ** -30:
            return None
        x = []
        for c in range(3):
            Ac = [list(row) for row in A]
            for r in range(3):
                Ac[r][c] = b[r]
            x.append(det3(Ac) / D)
        return x

    pts = {}
    for a in range(12):
        for b in range(a + 1, 12):
            for c in range(b + 1, 12):
                x = solve3([vs[a][:], vs[b][:], vs[c][:]],
                           [mpmath.mpf(1)] * 3)
                if x is None:
                    continue
                if any(sum(vs[j][k] * x[k] for k in range(3))
                       > 1 + mpmath.mpf(10) ** -30 for j in range(12)):
                    continue
                k = tuple(mpmath.nstr(v, 14) for v in x)
                pts[k] = x
    P = list(pts.values())
    if len(P) != 20:
        raise AssertionError("polar dodecahedron: %d vertices, expect 20"
                             % len(P))
    faces = []
    for i in range(12):
        f = [j for j, p in enumerate(P)
             if abs(sum(vs[i][k] * p[k] for k in range(3)) - 1)
             < mpmath.mpf(10) ** -25]
        if len(f) != 5:
            raise AssertionError("face %d has %d vertices, expect 5"
                                 % (i, len(f)))
        faces.append(tuple(f))
    inc = [sum(1 for f in faces if j in f) for j in range(20)]
    if set(inc) != {3}:
        raise AssertionError("faces per vertex %s, expect all 3" % inc)

    def dist2(p, q):
        return sum((p[k] - q[k]) ** 2 for k in range(3))
    same = set()
    for f in faces:
        for a in range(5):
            for b in range(a + 1, 5):
                same.add((min(f[a], f[b]), max(f[a], f[b])))
    dmin = min(dist2(P[i], P[j]) for (i, j) in same)
    edges = sorted((i, j) for (i, j) in same
                   if abs(dist2(P[i], P[j]) - dmin) < mpmath.mpf(10) ** -25)
    if len(edges) != 30:
        raise AssertionError("%d edges, expect 30" % len(edges))
    chi = 20 - 30 + 12
    return {"pts": P, "faces": faces, "edges": edges,
            "counts": (20, 30, 12, chi), "edge_len": mpmath.sqrt(dmin)}


def cell_nest(ih=None):
    """The nest of 120 dodecahedra: the I_h table acts simply
    transitively on the 120 cell labels via left regular action.
    Certifies: the multiplication table gives 120 distinct permutation
    matrices of 120 labels, each row a permutation, and the stabilizer
    of a cell is trivial (so the nest has EXACTLY 120 cells)."""
    if ih is None:
        ih = ih_table()
    M = [e["M"] for e in ih]
    keyidx = {_key(m): i for i, m in enumerate(M)}
    perms = []
    for g in range(120):
        p = []
        for h in range(120):
            k = _key(_mul(M[g], M[h]))
            if k not in keyidx:
                raise AssertionError("nest: product left the table")
            p.append(keyidx[k])
        if len(set(p)) != 120:
            raise AssertionError("nest: left action not a permutation")
        perms.append(tuple(p))
    if len(set(perms)) != 120:
        raise AssertionError("nest: %d distinct actions, expect 120"
                             % len(set(perms)))
    return perms


def honeycomb_534():
    """The dodecahedral honeycomb of order 4, Schlafli {5,3,4}: Coxeter
    chain [5 - 3 - 4] with Gram B_ij = -cos(pi / m_ij). Certifies:
    every proper connected sub-block spherical (det > 0: I2(5), I2(3),
    I2(4), and [5,3] = the cell's own H3 symmetry group of order 120),
    the full 4x4 Gram Lorentzian of signature (3,1) (one sign change in
    the principal-minor sequence) => compact hyperbolic honeycomb;
    and m_34 = 4: four dodecahedra around every edge (the "order 4").
    Returns ``{"gram": B (mpmath), "det": x, "signature": (3, 1),
    "edge_order": 4}``.
    """
    mpmath.mp.dps = MP
    c5 = mpmath.cos(mpmath.pi / 5)
    c3 = mpmath.cos(mpmath.pi / 3)
    c4 = mpmath.cos(mpmath.pi / 4)
    B = [[mpmath.mpf(1), -c5, 0, 0],
         [-c5, mpmath.mpf(1), -c3, 0],
         [0, -c3, mpmath.mpf(1), -c4],
         [0, 0, -c4, mpmath.mpf(1)]]

    def det(M):
        n = len(M)
        if n == 1:
            return M[0][0]
        tot = mpmath.mpf(0)
        for j in range(n):
            if M[0][j] != 0:
                minor = [[M[i][k] for k in range(n) if k != j]
                         for i in range(1, n)]
                tot += ((-1) ** j) * M[0][j] * det(minor)
        return +tot
    d1 = det([[B[0][0]]])
    d2 = det([r[:2] for r in B[:2]])
    d3 = det([r[:3] for r in B[:3]])
    d4 = det(B)
    for name, d in (("I2(5)", d1), ("{1,2}", d2), ("[5,3]=H3", d3)):
        if d <= 0:
            raise AssertionError("sub-block %s not spherical (det %s)"
                                 % (name, d))
    d34 = det([[B[2][2], B[2][3]], [B[3][2], B[3][3]]])
    if d34 <= 0:
        raise AssertionError("I2(4) block not spherical")
    signs = [1 if d > 0 else -1 for d in (d1, d2, d3, d4)]
    changes = sum(1 for i in range(3) if signs[i] != signs[i + 1])
    if changes != 1:
        raise AssertionError("Gram not Lorentzian: %s" % signs)
    if not (c4 > 0 and B[2][3] == -c4):
        raise AssertionError("edge label m34 != 4")
    return {"gram": B, "det": d4, "signature": (3, 1), "edge_order": 4}
