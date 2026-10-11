#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V140 - the dressed gauge gate on the directed-link network, and the
dispersion of the quintet.

**Where this stands.**  V138 closed the qubit-link route to a lattice spin 2:
the 8-dimensional bond factor is A1g + T2g + A2u + T1u - it carries T2g but
NEVER Eg, so the symmetric-traceless quintet Eg + T2g is unavailable in the
qubit link at any coupling, and the gapless GAMMA multiplets are flat bands.
The channel that DOES carry the quintet is the symmetric rank-2 space of the
LINK DIRECTIONS (V138 T5) - the directional bond layer of Section 10. V79
built that layer (one Spin(10) adjoint link per bond, exact Gauss algebra)
and named its open step: the naive bond gate K = exp(i g10 J.L) preserves
the bond diagonal charge but LEAKS the single-site Gauss charge - "the
gauge-invariant gate needs Kogut-Susskind dressing (named open step)".
This experiment BUILDS the dressed gate on the directed-link network and
then runs the dispersion of the quintet.

**Part A - the dressed gate (the directed-link network).**  A bond is a
DIRECTED link x -> x': the tail carries the link generator +L^A, the head
the conjugate-representation generator -Lbar^A. Computed:
- T1 the naive sign convention C_x = T_x + L_left - L_right (the V79 chain
  note) does NOT close the Gauss algebra at a vertex with two links (the
  same-link commutator carries s^2 = +1); the closing convention is tail
  +L, head -Lbar, exact on the two-link star.
- T2 the bare V79 gate leaks the vertex charge (re-measured on the directed
  link, leak O(1)).
- T3 the Kogut-Susskind ELECTRIC gate K_E = exp(i g10 C_x^2), C_x the
  directed vertex charge, is EXACTLY gauge-invariant (commutes with all 45
  C_x^A), exactly unitary, and its spectrum is the V79 content theorem:
  phases g10 * {45, 85, 117} with multiplicities 32, 288, 1120.
- T4 the two-site DRESSED HOPPING is obstructed on the fixed-representation
  link: the flux-dressed current-current C_x . C_x' and its
  transporter-completed version both leak the vertex charges (measured);
  an exactly invariant hopping needs the link in the REGULAR representation
  (a group-algebra link with a true parallel transporter). The named open
  step of V79 is therefore split by computation: electric part CLOSED,
  hopping part obstructed with a named fix.

**Part B - the dispersion of the quintet (the directional network).**  On
the directional network the link flux is the link ORIENTATION (a c-number
sign), so the dressed hop exists exactly: the hop along axis j carries the
vertex dressing exp(-i g10 sin(q_j) J_j), J_k the spin-2 generators on the
symmetric rank-2 space sym^2 (the lattice home of the metric, V138 T5).
The one-step Bloch matrix M(q) = exp(-i g10 sum_j sin(q_j) J_j) is exactly
unitary and O_h-covariant (the 24 proper rotations exactly; the 24
improper ones map the step to its inverse - the directed network is
chiral). Computed:
- T6 at GAMMA the step is the identity: SIX gapless bands, content
  A1g + Eg + T2g - the quintet IS there, where the qubit bond factor
  never had Eg (V138 T1).
- T7 the quintet dispersion is MASSLESS and LINEAR with the helicity
  ladder speeds {0, g10, 2 g10} (the m = 0 branch flat - the longitudinal
  mode does not propagate; the A1g trace band flat too), EXACTLY
  isotropic to leading order (the spin-2 spectrum theorem: n.J has
  eigenvalues {-2,-1,0,1,2} for every direction n - direction spread
  O(eps^2) at |q| = eps), and UV-bounded (max 2 g10 sqrt(3)).
This is the kinematics the qubit link never had (V138: flat at GAMMA;
V139: anisotropy 0.857 at the nodal nodes): the escape route of V138 T5
works kinematically on the directional layer.

**What is NOT derived (stated).**  The gauge group of the directional
network is the lattice direction group O_h (exact covariance), not
Spin(10); the Spin(10)-invariant dressed hopping is obstructed on the
fixed adjoint link (T4); the ladder speeds g10 * {1, 2} = {0.2, 0.4} are
not 1 - the dressing angle per hop is an input of the construction; the
coupling of the quintet to the matter layer and the non-linear completion
are not constructed.

Machinery: the root's own 45 generators -> f^{ABC} -> adjoint link (V79);
the ten O_h irreps built from the defining representation (V138); the
spin-2 generators on sym^2 constructed here and validated against su(2)
and the spectrum theorem. Zero new parameters: g10 = 1/5 everywhere.
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10 import engine as E                                          # noqa: E402
from cl10 import sm_qca                                               # noqa: E402
from cl10 import supercell as SC                                      # noqa: E402

ETA = SC.G10                       # g10 = 1/5, the root angle
TOL_MACH = 1e-10
TOL_GAP = 1e-9


# ---------------------------------------------------------------------------
# the root: 45 generators -> f^{ABC} -> adjoint link operators (V79 route)
# ---------------------------------------------------------------------------
def root_generators():
    return [E.asnumpy(T).astype(np.complex128)
            for T in sm_qca.plane_generators()]


def structure_constants(Ts):
    n = len(Ts)
    k = float(np.real(np.trace(Ts[0] @ Ts[0])))
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(a + 1, n):
            comm = Ts[a] @ Ts[b] - Ts[b] @ Ts[a]
            for c in range(n):
                v = float(np.real(-1j * np.trace(comm @ Ts[c]) / k))
                f[a, b, c] = v
                f[a, c, b] = -v
                f[b, a, c] = -v
    return f, k


def adjoint_link_operators(f):
    return [np.array(-1j * f[A], dtype=np.complex128) for A in range(45)]


# ---------------------------------------------------------------------------
# O_h: signed permutations, classes, the ten irreps built from the defining
# representation, and the symmetric rank-2 (metric) representation (V138)
# ---------------------------------------------------------------------------
def oh_elements():
    return [(p, s) for p in itertools.permutations(range(3))
            for s in itertools.product((1, -1), repeat=3)]


OH = oh_elements()


def m3(perm, signs):
    A = np.zeros((3, 3))
    for i in range(3):
        A[i, perm[i]] = signs[i]
    return A


def perm_sign(perm):
    P = np.array([[1.0 if j == perm[i] else 0.0 for j in range(3)]
                  for i in range(3)])
    return round(float(np.linalg.det(P)))


def classes():
    seen, out = set(), []
    for g in OH:
        if g in seen:
            continue
        cls = set()
        for b in OH:
            B = m3(*b)
            M = B @ m3(*g) @ np.linalg.inv(B)
            p = tuple(int(np.argmax(np.abs(M[i]))) for i in range(3))
            s = tuple(int(M[i, p[i]]) for i in range(3))
            cls.add((p, s))
        seen |= cls
        out.append(sorted(cls))
    return out


CLS = classes()

IDX = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def sym2_basis():
    """The six orthonormal symmetric 3x3 matrices of the IDX basis."""
    B = []
    for (i, j) in IDX:
        S = np.zeros((3, 3))
        if i == j:
            S[i, i] = 1.0
        else:
            S[i, j] = S[j, i] = 1.0 / math.sqrt(2.0)
        B.append(S)
    return B


SB = sym2_basis()


def sym2_rep(perm, signs):
    """The 6-dimensional O_h action on symmetric rank-2 tensors."""
    A = m3(perm, signs)
    D = np.zeros((6, 6))
    for b in range(6):
        Sp = A @ SB[b] @ A.T
        for a in range(6):
            D[a, b] = float(np.trace(SB[a].T @ Sp))
    return D


def build_irreps():
    idv = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]) / math.sqrt(3.0)
    K = np.eye(6) - np.outer(idv, idv)
    C2 = [((0, 1, 2), (-1, -1, 1)), ((0, 1, 2), (-1, 1, -1)),
          ((0, 1, 2), (1, -1, -1))]
    p = sum(K @ sym2_rep(*v) @ K
            for v in [((0, 1, 2), (1, 1, 1))] + C2) / 4.0
    w, V = np.linalg.eigh((p + p.conj().T) / 2)
    Eg = V[:, w > 0.5]
    Q, _ = np.linalg.qr(np.column_stack([idv, Eg]))
    rest = K - Q @ Q.conj().T
    w2, V2 = np.linalg.eigh((rest + rest.conj().T) / 2)
    T2 = V2[:, w2 > 0.5]

    def one(x):
        return np.array([[x]])

    reps = {
        "A1g": lambda g: one(1.0),
        "A2g": lambda g: one(perm_sign(g[0])),
        "Eg": lambda g: Eg.conj().T @ sym2_rep(*g) @ Eg,
        "T2g": lambda g: T2.conj().T @ sym2_rep(*g) @ T2,
        "T1g": lambda g: np.linalg.det(m3(*g)) * m3(*g),
        "A1u": lambda g: one(np.linalg.det(m3(*g))),
        "A2u": lambda g: one(perm_sign(g[0]) * np.linalg.det(m3(*g))),
        "Eu": lambda g: np.linalg.det(m3(*g)) * (Eg.conj().T @ sym2_rep(*g) @ Eg),
        "T1u": lambda g: m3(*g),
        "T2u": lambda g: np.linalg.det(m3(*g)) * (T2.conj().T @ sym2_rep(*g) @ T2),
    }
    chi = {r: [float(np.real(np.trace(f(c[0])))) for c in CLS]
           for r, f in reps.items()}
    return reps, chi, {"Eg": Eg, "T2g": T2, "idv": idv}


REPS, CHI, _IR = build_irreps()


def inner(a, b):
    return sum(len(c) * x * y for c, x, y in zip(CLS, a, b)) / 48.0


def decompose(chi_v):
    return {r: round(inner(chi_v, CHI[r]), 4) for r in CHI}


def chi_of(repfun):
    return [float(np.real(np.trace(np.atleast_2d(repfun(*c[0]))))) for c in CLS]


# ---------------------------------------------------------------------------
# the spin-2 generators on sym^2: J_k(S) = i (A_k S - S A_k), A_k the
# antisymmetric rotation generators of the vector representation
# ---------------------------------------------------------------------------
def rotation_generators():
    """A_k real antisymmetric: (A_k)_{ij} = -eps_{kij}."""
    lev = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
           (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
    A = []
    for k in range(3):
        M = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                M[i, j] = -lev.get((k, i, j), 0)
        A.append(M)
    return A


def spin2_generators():
    """The 6x6 Hermitian generators J_k on the orthonormal sym^2 basis."""
    A = rotation_generators()
    J = []
    for k in range(3):
        M = np.zeros((6, 6), dtype=complex)
        for b in range(6):
            act = 1j * (A[k] @ SB[b] - SB[b] @ A[k])
            for a in range(6):
                M[a, b] = np.trace(SB[a].T @ act)
        J.append(M)
    return J


J6 = spin2_generators()
TRACE_V = np.zeros(6, dtype=complex)
for a in range(6):
    TRACE_V[a] = np.trace(SB[a]) / math.sqrt(3.0)


def M_dir(q):
    """The dressed directional one-step Bloch matrix on sym^2."""
    Hq = sum(math.sin(qj) * J for qj, J in zip(q, J6))
    return expm(-1j * ETA * Hq)


def bands(q):
    return np.sort(np.abs(np.angle(np.linalg.eigvals(M_dir(q)))))


def fib_dirs(n=100):
    ga = (1 + math.sqrt(5)) / 2
    out = []
    for i in range(n):
        z = 1 - 2 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1 - z * z))
        th = 2 * math.pi * i / ga
        out.append((r * math.cos(th), r * math.sin(th), z))
    return out


# ---------------------------------------------------------------------------
# Part A helpers: vertex charges as matvecs
# ---------------------------------------------------------------------------
def apply_T(v3, TA):
    """(T^A on the matter factor) v, v3 shape (32, 32, 45)."""
    return (TA @ v3.reshape(32, 32 * 45)).reshape(32, 32, 45)


def apply_L(v3, LA, link):
    """(L^A on link factor `link`) v."""
    if link == 1:
        return np.einsum('bce,de->bcd', v3, LA)
    return np.einsum('abe,de->abd', v3, LA)


def main():
    t0 = time.time()
    out = {"experiment": "V140_dressed_gate_quintet", "data": {}}
    checks = []

    # ---- T0 the machinery is self-validating -------------------------------
    Ts = root_generators()
    f, k32 = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    rng = np.random.default_rng(7)
    alg_err = 0.0
    for _ in range(60):
        a, b = (int(x) for x in rng.integers(0, 45, 2))
        if a == b:
            continue
        comm = Ls[a] @ Ls[b] - Ls[b] @ Ls[a]
        rhs = 1j * sum(f[a, b, c] * Ls[c] for c in range(45)
                       if abs(f[a, b, c]) > 1e-12)
        alg_err = max(alg_err, float(np.abs(comm - rhs).max()))
    cas64 = float(np.abs(sum(L @ L for L in Ls) - 64.0 * np.eye(45)).max())
    cas45 = float(np.abs(sum(T @ T for T in Ts) - 45.0 * np.eye(32)).max())
    chi_sym2 = chi_of(sym2_rep)
    sym2_dec = decompose(chi_sym2)
    sum_d2 = sum(int(round(CHI[r][0])) ** 2 for r in CHI)
    su2_err = max(float(np.abs(J6[i] @ J6[j] - J6[j] @ J6[i]
                               - 1j * J6[l]).max())
                  for (i, j, l) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
    herm_err = max(float(np.abs(J - J.conj().T).max()) for J in J6)
    tr_ann = max(float(np.abs(J @ TRACE_V).max()) for J in J6)
    spec_err = 0.0
    for nhat in fib_dirs(20):
        nn = math.sqrt(sum(x * x for x in nhat))
        nJ = sum(c * J for c, J in zip(nhat, J6)) / nn
        ev = np.sort(np.linalg.eigvalsh(nJ))
        spec_err = max(spec_err, float(np.abs(ev - np.array(
            [-2, -1, 0, 0, 1, 2])).max()))
    t0ok = (alg_err < TOL_MACH and cas64 < TOL_MACH and cas45 < TOL_MACH
            and sym2_dec["A1g"] > 0.99 and sym2_dec["Eg"] > 0.99
            and sym2_dec["T2g"] > 0.99 and sum_d2 == 48
            and su2_err < 1e-12 and herm_err < 1e-12 and tr_ann < 1e-12
            and spec_err < 1e-10)
    out["data"]["machinery"] = {
        "link_algebra_err": alg_err, "casimir64_err": cas64,
        "casimir45_err": cas45, "sym2_decomposition": sym2_dec,
        "sum_d_squared": sum_d2, "su2_err": su2_err,
        "trace_annihilated_err": tr_ann,
        "spin2_spectrum_err": spec_err}
    checks.append(("T0 machinery validates: adjoint link algebra to %.1e, "
                   "Casimirs 64/45 to %.1e/%.1e, sym^2 = %s, sum d^2 = %d, "
                   "spin-2 generators: su(2) to %.1e, trace annihilated, "
                   "n.J spectrum = {-2,-1,0,0,1,2} over 20 directions to "
                   "%.1e" % (alg_err, cas64, cas45, sym2_dec, sum_d2,
                             su2_err, spec_err), t0ok))

    # ---- T1 the directed-link Gauss convention -----------------------------
    # two-link star at x: link 1 arrives (x is its head), link 2 leaves (x
    # is its tail). Naive convention: C = T - L1 + L2. Closing convention:
    # the head carries the CONJUGATE generator: C = T - Lbar1 + L2.
    Lbar = [L.conj() for L in Ls]

    def C_star_naive(v, A):      # v shape (32, 45, 45): matter, link1, link2
        return (np.einsum('ab,bxy->axy', Ts[A], v)
                - np.einsum('ax,mxy->may', Ls[A], v)
                + np.einsum('mxy,ay->mxa', v, Ls[A]))

    def C_star_close(v, A):
        return (np.einsum('ab,bxy->axy', Ts[A], v)
                - np.einsum('ax,mxy->may', Lbar[A], v)
                + np.einsum('mxy,ay->mxa', v, Ls[A]))

    rng2 = np.random.default_rng(11)
    leak_naive, leak_close = 0.0, 0.0
    for _ in range(30):
        A, B = (int(x) for x in rng2.integers(0, 45, 2))
        if A == B:
            continue
        v = (rng2.standard_normal((32, 45, 45))
             + 1j * rng2.standard_normal((32, 45, 45)))
        for Cf, acc in ((C_star_naive, "naive"), (C_star_close, "close")):
            lhs = Cf(Cf(v, B), A) - Cf(Cf(v, A), B)
            rhs = 1j * sum(f[A, B, c] * Cf(v, c) for c in range(45)
                           if abs(f[A, B, c]) > 1e-12)
            e = float(np.linalg.norm(lhs - rhs)) / float(np.linalg.norm(v))
            if acc == "naive":
                leak_naive = max(leak_naive, e)
            else:
                leak_close = max(leak_close, e)
    t1ok = leak_close < 1e-10 and leak_naive > 1e-2
    out["data"]["gauss_convention"] = {
        "naive_sign_leak": leak_naive, "conjugate_leak": leak_close,
        "note": "same-link commutator carries s^2 = +1: the head must "
                "carry the conjugate generator -Lbar, not -L"}
    checks.append(("T1 the directed-link Gauss charge closes ONLY with the "
                   "conjugate head convention (tail +L, head -Lbar): "
                   "[C^A,C^B] = i f C on the two-link star to %.1e, while "
                   "the naive V79 chain sign convention leaks %.1e"
                   % (leak_close, leak_naive), t1ok))

    # ---- T2 the bare V79 gate leaks the vertex charge ----------------------
    # directed link x -> x': tail x (C_x = T_x + L), head x' (C_x' = T_x' -
    # Lbar), on (32, 32, 45).
    def apply_bondH(v3):         # the full bare gate H = g10 sum_A J^A L^A
        out_ = np.zeros_like(v3)
        for A in range(45):
            w = np.einsum('bce,de->bcd', v3, Ls[A])
            out_ = out_ + (Ts[A] @ w.reshape(32, 32 * 45)
                           ).reshape(32, 32, 45)
            out_ = out_ + (Ts[A] @ w.transpose(1, 0, 2).reshape(32, 32 * 45)
                           ).reshape(32, 32, 45).transpose(1, 0, 2)
        return ETA * out_

    def C_tail(v3, A):
        return apply_T(v3, Ts[A]) + apply_L(v3, Ls[A], 1)

    leak_bare = 0.0
    for A in range(0, 45, 5):
        v3 = (rng2.standard_normal((32, 32, 45))
              + 1j * rng2.standard_normal((32, 32, 45)))
        lhs = apply_bondH(C_tail(v3, A))
        rhs = C_tail(apply_bondH(v3), A)
        leak_bare = max(leak_bare, float(np.linalg.norm(lhs - rhs))
                        / float(np.linalg.norm(v3)))
    t2ok = leak_bare > 1e-2
    out["data"]["bare_gate"] = {"vertex_leak": leak_bare}
    checks.append(("T2 the bare V79 bond gate K = exp(i g10 J.L) leaks the "
                   "directed vertex charge (relative leak %.2e): the V79 "
                   "named open step re-measured on the directed link"
                   % leak_bare, t2ok))

    # ---- T3 the dressed KS electric gate: exactly invariant ----------------
    # C_x^A = T^A (x) I + I (x) L^A on 32 x 45 (the tail vertex of one
    # directed link); K_E = exp(i g10 sum_A C^A C^A). The generator
    # expands through the two Casimirs of T0:
    #   sum_A C^A C^A = (45 + 64) I + 2 sum_A T^A (x) L^A,
    # a single 1440 x 1440 matrix; the invariance test uses the charge as
    # a matvec, never as a stored 45-fold stack.
    C2op = 109.0 * np.eye(1440, dtype=np.complex128)
    for T, L in zip(Ts, Ls):
        C2op += 2.0 * np.kron(T, L)
    herm_C2 = float(np.abs(C2op - C2op.conj().T).max())
    KE = expm(1j * ETA * C2op)
    unit_E = float(np.abs(KE @ KE.conj().T - np.eye(1440)).max())

    def C_matvec(v, A):          # v shape (32, 45)
        return (np.einsum('ab,be->ae', Ts[A], v)
                + np.einsum('be,ae->ab', Ls[A], v))
    inv_E = 0.0
    for A in range(45):
        for _ in range(2):
            v = (rng2.standard_normal((32, 45))
                 + 1j * rng2.standard_normal((32, 45)))
            vf = v.reshape(1440)
            lhs = C_matvec((KE @ vf).reshape(32, 45), A).reshape(1440)
            rhs = KE @ C_matvec(v, A).reshape(1440)
            inv_E = max(inv_E, float(np.linalg.norm(lhs - rhs))
                        / float(np.linalg.norm(vf)))
    ev = np.linalg.eigvalsh(C2op)
    clusters = []
    for e in ev:
        if clusters and abs(e - clusters[-1][0]) < 1e-6 * max(1.0, abs(e)):
            clusters[-1][1] += 1
        else:
            clusters.append([float(e), 1])
    content_ok = ([c[1] for c in clusters] == [32, 288, 1120]
                  and all(abs(c[0] - t) < 1e-6
                          for c, t in zip(clusters, (45.0, 85.0, 117.0))))
    phases = sorted({round((ETA * c[0]) % (2 * math.pi), 6)
                     for c in clusters})
    t3ok = (herm_C2 < 1e-12 and unit_E < 1e-12 and inv_E < 1e-10
            and content_ok)
    out["data"]["dressed_electric_gate"] = {
        "generator_hermiticity_err": herm_C2, "unitarity_err": unit_E,
        "gauss_invariance_err_max45": inv_E,
        "C2_levels": [round(c[0], 4) for c in clusters],
        "multiplicities": [c[1] for c in clusters],
        "gate_phases_mod2pi": phases}
    checks.append(("T3 the Kogut-Susskind ELECTRIC gate K_E = exp(i g10 "
                   "C_x^2) is EXACTLY gauge-invariant: [C^A, K_E] = 0 for "
                   "all 45 A to %.1e, unitary to %.1e, and its spectrum is "
                   "the V79 content theorem: C2 levels %s with "
                   "multiplicities %s (phases g10 x {45,85,117})"
                   % (inv_E, unit_E, [c[0] for c in clusters],
                      [c[1] for c in clusters]), t3ok))

    # ---- T4 the two-site dressed hopping is obstructed ---------------------
    # candidates on (32, 32, 45): (a) flux-dressed current-current
    # G1 = sum_A C_tail^A C_head^A; (b) transporter-completed
    # G2 = sum_A C_tail^A (I x I x V) C_head^A (I x I x V^dag), V = e^{i g10 L}.
    def C_head(v3, A):
        return apply_T(v3.transpose(1, 0, 2), Ts[A]).transpose(1, 0, 2) \
            - apply_L(v3, Lbar[A], 1)

    def G1(v3):
        """sum_A (C_tail^A C_head^A + C_head^A C_tail^A)/2 v - Hermitian."""
        s = np.zeros_like(v3)
        for A in range(45):
            s = s + 0.5 * (C_tail(C_head(v3, A), A) + C_head(C_tail(v3, A), A))
        return s

    Vs = [expm(1j * ETA * L) for L in Ls]

    def G2(v3):
        """sum_A C_tail^A (I x I x V) C_head^A (I x I x V^dag) v."""
        s = np.zeros_like(v3)
        for A in range(45):
            w = np.einsum('bcd,de->bce', v3, Vs[A].conj().T)
            w = C_head(w, A)
            w = np.einsum('bcd,de->bce', w, Vs[A])
            s = s + C_tail(w, A)
        return s

    leak1_t, leak1_h, leak2_t, leak2_h = 0.0, 0.0, 0.0, 0.0
    for A in range(0, 45, 5):
        v3 = (rng2.standard_normal((32, 32, 45))
              + 1j * rng2.standard_normal((32, 32, 45)))
        g1 = G1(v3)
        leak1_t = max(leak1_t, float(np.linalg.norm(
            C_tail(g1, A) - G1(C_tail(v3, A)))) / np.linalg.norm(v3))
        leak1_h = max(leak1_h, float(np.linalg.norm(
            C_head(g1, A) - G1(C_head(v3, A)))) / np.linalg.norm(v3))
        g2 = G2(v3)
        leak2_t = max(leak2_t, float(np.linalg.norm(
            C_tail(g2, A) - G2(C_tail(v3, A)))) / np.linalg.norm(v3))
        leak2_h = max(leak2_h, float(np.linalg.norm(
            C_head(g2, A) - G2(C_head(v3, A)))) / np.linalg.norm(v3))
    t4ok = min(leak1_t, leak1_h, leak2_t, leak2_h) > 1e-3
    out["data"]["hopping_obstruction"] = {
        "current_current_leak_tail": leak1_t,
        "current_current_leak_head": leak1_h,
        "transporter_leak_tail": leak2_t,
        "transporter_leak_head": leak2_h,
        "note": "no exactly invariant two-site gate on the fixed adjoint "
                "link; the invariant hopping needs a regular-representation "
                "(group-algebra) link with a true parallel transporter"}
    checks.append(("T4 the two-site DRESSED HOPPING is obstructed on the "
                   "fixed adjoint link: the flux-dressed current-current "
                   "leaks the tail/head vertex charges (%.2e / %.2e) and "
                   "the transporter-completed version leaks (%.2e / %.2e) "
                   "- V79's open step split by computation: electric part "
                   "CLOSED (T3), hopping part obstructed, named fix = "
                   "regular-representation link"
                   % (leak1_t, leak1_h, leak2_t, leak2_h), t4ok))

    # ---- T5 the dressed directional gate: unitary + O_h-covariant ----------
    rng3 = np.random.default_rng(3)
    unit_dir, cov_proper, cov_improper_inv = 0.0, 0.0, 0.0
    qs = [(0.0, 0.0, 0.0)] + [tuple(rng3.uniform(-math.pi, math.pi, 3))
                              for _ in range(8)]
    for q in qs:
        M = M_dir(q)
        unit_dir = max(unit_dir, float(np.abs(M @ M.conj().T
                                             - np.eye(6)).max()))
        for g in OH:
            R = m3(*g)
            D = sym2_rep(*g)
            gq = tuple(float(np.dot(R[i], q)) for i in range(3))
            Mgq = M_dir(gq)
            lhs = D @ M @ D.conj().T
            det = round(float(np.linalg.det(R)))
            if det > 0:
                cov_proper = max(cov_proper, float(np.abs(lhs - Mgq).max()))
            else:
                cov_improper_inv = max(
                    cov_improper_inv,
                    float(np.abs(lhs - Mgq.conj().T).max()))
    t5ok = unit_dir < 1e-12 and cov_proper < 1e-12 and cov_improper_inv < 1e-12
    out["data"]["directional_gate"] = {
        "unitarity_err": unit_dir, "covariance_proper_err": cov_proper,
        "covariance_improper_as_inverse_err": cov_improper_inv,
        "note": "the 24 improper elements map the step to its inverse: "
                "the directed network is chiral"}
    checks.append(("T5 the dressed directional step M(q) = exp(-i g10 "
                   "sin(q).J) is exactly unitary (%.1e) and O_h-covariant: "
                   "D(g) M(q) D(g)^dag = M(gq) for the 24 proper rotations "
                   "to %.1e; the 24 improper ones give M(gq)^dag instead "
                   "(to %.1e) - the directed network is chiral"
                   % (unit_dir, cov_proper, cov_improper_inv), t5ok))

    # ---- T6 the GAMMA content: the quintet IS there ------------------------
    M0 = M_dir((0.0, 0.0, 0.0))
    gap0 = float(np.abs(M0 - np.eye(6)).max())
    chi_dir = chi_of(sym2_rep)
    dir_dec = decompose(chi_dir)
    t6ok = (gap0 < 1e-14 and dir_dec["Eg"] > 0.99 and dir_dec["T2g"] > 0.99
            and dir_dec["A1g"] > 0.99)
    out["data"]["gamma_content"] = {
        "M0_minus_I": gap0, "n_gapless": 6,
        "decomposition": dir_dec}
    checks.append(("T6 at GAMMA the directional step is the identity "
                   "(dev %.1e): SIX gapless bands with content %s - the "
                   "quintet Eg+T2g IS present in the directional layer, "
                   "where the qubit bond factor A1g+T2g+A2u+T1u never had "
                   "Eg (V138 T1)" % (gap0, dir_dec), t6ok))

    # ---- T7 the dispersion of the quintet ----------------------------------
    eps_list = (1e-2, 1e-3, 1e-4)
    axes = {"G-X": (1, 0, 0), "G-M": (1, 1, 0), "G-R": (1, 1, 1)}
    disp = {}
    for name, d in axes.items():
        nn = math.sqrt(sum(x * x for x in d))
        rows = {}
        for eps in eps_list:
            q = tuple(eps * v / nn for v in d)
            b = bands(q)
            rows["%.0e" % eps] = [round(float(x), 8) for x in b]
        # local exponent of the top band
        b2 = rows["1e-04"][-1] / 1e-4
        expo = math.log(rows["1e-04"][-1] / rows["1e-02"][-1]) / math.log(1e-2)
        disp[name] = {"bands_by_eps": rows, "top_speed": b2,
                      "local_exponent": expo}
    # the helicity ladder: speeds of the 6 bands at eps = 1e-4 along (1,0,0)
    ladder = [round(float(x) / 1e-4, 4) for x in bands((1e-4, 0.0, 0.0))]
    ladder_ok = (all(abs(a) < 1e-3 for a in ladder[:2])
                 and all(abs(a - ETA) < 1e-3 for a in ladder[2:4])
                 and all(abs(a - 2 * ETA) < 1e-3 for a in ladder[4:]))
    expo_ok = all(abs(v["local_exponent"] - 1.0) < 0.01
                  for v in disp.values())
    # isotropy: top-band speed over 100 directions at fixed |q| = eps
    iso = {}
    for eps in (1e-2, 1e-3):
        sp = [bands(tuple(eps * v for v in n))[-1] / eps for n in fib_dirs(100)]
        iso["%.0e" % eps] = {"min": round(min(sp), 6),
                             "max": round(max(sp), 6),
                             "spread": round(max(sp) - min(sp), 8)}
    iso_ok = iso["1e-03"]["spread"] < 1e-4 and iso["1e-02"]["spread"] < 1e-2
    # UV bound: max over a 13^3 zone grid
    grid = np.linspace(-math.pi, math.pi, 13)
    top_max = 0.0
    for qx in grid:
        for qy in grid:
            for qz in grid:
                top_max = max(top_max, float(bands((qx, qy, qz))[-1]))
    uv_pred = 2 * ETA * math.sqrt(3.0)
    uv_ok = abs(top_max - uv_pred) < 0.05
    t7ok = ladder_ok and expo_ok and iso_ok and uv_ok
    out["data"]["dispersion"] = {
        "bands_by_direction": disp, "helicity_ladder_speeds": ladder,
        "predicted_speeds": [0.0, 0.0, ETA, ETA, 2 * ETA, 2 * ETA],
        "isotropy": iso, "zone_top_max": top_max,
        "uv_bound_2g10sqrt3": uv_pred}
    checks.append(("T7 the quintet dispersion is MASSLESS LINEAR with the "
                   "helicity ladder speeds %s (predicted {0,0,g10,g10,"
                   "2g10,2g10} = %s): local exponent %s, direction spread "
                   "of the top speed %s at |q| = eps (isotropy exact to "
                   "leading order, the spin-2 spectrum theorem), UV-bounded "
                   "by 2 g10 sqrt(3) = %.4f (zone max %.4f)"
                   % (ladder, [round(ETA * m, 3) for m in (0, 0, 1, 1, 2, 2)],
                      {k: round(v["local_exponent"], 3)
                       for k, v in disp.items()}, iso, uv_pred, top_max),
                  t7ok))

    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["pass"] = ok
    out["status"] = (
        "THE DRESSED GAUGE GATE ON THE DIRECTED-LINK NETWORK, AND THE "
        "QUINTET DISPERSION. Part A (the V79 open step, split by "
        "computation): the directed-link Gauss charge closes ONLY with the "
        "conjugate head convention (tail +L, head -Lbar) - the naive chain "
        f"sign convention leaks {leak_naive:.1e} on the two-link star "
        f"(T1); the bare bond gate leaks the vertex charge {leak_bare:.2f} "
        "(T2); the Kogut-Susskind ELECTRIC gate exp(i g10 C_x^2) is "
        f"exactly gauge-invariant (max over all 45 charges {inv_E:.1e}) "
        "and its spectrum IS the V79 content theorem, phases g10 x "
        "{45, 85, 117} with multiplicities 32, 288, 1120 (T3); every "
        "two-site dressed-hopping candidate on the fixed adjoint link "
        f"leaks (min {min(leak1_t, leak1_h, leak2_t, leak2_h):.2e}) - the "
        "invariant hopping needs a regular-representation link (T4). "
        "Part B (the V138 T5 escape route, run): on the directional "
        "network the flux is the link orientation (a c-number), the "
        "dressed hop exists exactly, and the step M(q) = exp(-i g10 "
        "sin(q).J) is unitary and O_h-covariant (improper elements invert "
        "the step: the network is chiral, T5). At GAMMA it is the "
        "identity: six gapless bands with content A1g + Eg + T2g - the "
        "quintet is present where the qubit link never had Eg (T6). Its "
        "dispersion is massless, linear (exponent 1.00), EXACTLY isotropic "
        "to leading order (direction spread O(eps^2)), with the helicity "
        "ladder speeds {0, g10, 2 g10} - the m = 0 branch and the A1g "
        "trace do not propagate - and UV-bounded by 2 g10 sqrt(3) (T7). "
        "VERDICT: the directional bond layer supplies the kinematics the "
        "qubit link was proved not to have (V138: no Eg, flat bands; "
        "V139: anisotropy 0.857): a massless, linear, isotropic quintet "
        "at the zone centre. What is NOT derived: the directional "
        "network's gauge group is the lattice direction group O_h, not "
        "Spin(10); the Spin(10)-invariant dressed hopping is obstructed "
        "on the fixed adjoint link (T4); the ladder speeds g10 x {1,2} = "
        "{0.2, 0.4} are not 1 (the dressing angle per hop is an input); "
        "the coupling of the quintet to the matter layer and the "
        "non-linear completion are not constructed.")

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v140_dressed_gate_quintet.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V140 - the dressed gate on the directed-link network + the "
          "quintet dispersion")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    print("  wall time %.1f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
