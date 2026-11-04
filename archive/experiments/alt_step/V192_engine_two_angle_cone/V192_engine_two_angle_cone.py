#!/usr/bin/env python
"""V192 - can the engine's OWN six-stage product reach a conical node at Gamma?

V191 settled the two halves of the cone question: no single-parity multiplet at
an inversion-invariant momentum can be conical (the parity no-go), and a
parity-degenerate A1g/T1u pair gives an exactly isotropic non-abelian cone - but
built on exp(-i g A), the exponential of a hopping, which is exponentially local
and so not the engine's finite-depth circuit.

V191 also measured where the symmetry lives: the six-stage product is exactly
O_h-equivariant AT GAMMA (1.11e-16) and broken at every other momentum
(1.07e-1 generic, 1.40e-1 M, 1.46e-1 R). So the obstacle at the node is not the
symmetry, it is the spectrum: at the single angle the root uses the four bond
levels are +6, -2, -6, +2 and the A1g and T1u levels do not pair.

The largest family of six-stage products that keeps Gamma exactly covariant is
the axis-symmetric two-angle one: one angle a on the three intra-block stages,
one angle b on the three inter-block ones. This programme searches that family,
together with the internal angle, for a four-state gapless node and measures
what sits there.

Three points about this search, each easy to get wrong.

*(i) The pairing condition.* The step is the Kronecker product of the bond
factor and the internal factor, so its 128 phases are the pairwise sums of the 8
bond phases and the 16 internal ones, and a state is gapless exactly when one
such sum vanishes. A four-state node needs TWO such cancellations at the same
point. The internal spectrum is +-symmetric, which supplies the partner of a
cancellation but at the opposite internal phase, not at one half-turn away: two
bond levels at opposite phases therefore do NOT pair for free. A search for a
degeneracy of the two bond phases through a root-finder on the wrapped phase
difference converges on the wrap discontinuity, and it measures the cone on a
two-level subspace it has assumed gapless.

*(ii) The subspace.* Everything here is measured on the gapless eigenspace of
the full one-step matrix, orthonormalised (eig returns an arbitrary, non-
orthogonal basis for a degenerate eigenvalue even of a unitary matrix), and its
O_h content is decomposed by character inner product rather than assumed.

*(iii) The internal phases.* They are computed from the sixth power of the
internal step itself, never from a generator: C.U is a PRODUCT of exponentials
of the 45 plane rotations, not the exponential of their sum, so the eigenphases
of U^6 are not linear in the angle: they are never read off a small-angle
expansion: that would be wrong by 4e-3 at theta = 2.36 - the whole size of the
gap being asked about.
"""
import itertools
import json
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402
from scipy.optimize import brentq, minimize  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10 import altcell as A                       # noqa: E402
from cl10 import supercell as SC                    # noqa: E402
from cl10.altcell import FactorizedCore             # noqa: E402

C = FactorizedCore((2, 2, 2), "baseline")
NS = C.nsub
GAM = np.zeros(3)
DIRS = C.dirs


def coord(e):
    return np.array([2 * C.c["digit"](e, j) - 1.0 for j in range(3)])


G48 = []
for _perm in itertools.permutations(range(3)):
    for _sg in itertools.product((1, -1), repeat=3):
        R3 = np.zeros((3, 3))
        for i, (p, s) in enumerate(zip(_perm, _sg)):
            R3[i, p] = s
        G48.append(R3)


def _perm_of(R3):
    out = []
    for e in range(NS):
        y = R3 @ coord(e)
        out.append(C.c["site_index"](tuple(int(round((v + 1) / 2)) for v in y)))
    return np.array(out)


W48 = []
for _P in [_perm_of(R3) for R3 in G48]:
    W = np.zeros((NS, NS))
    for e in range(NS):
        W[_P[e], e] = 1.0
    W48.append(W)


def class_of(R3):
    tr, det = round(np.trace(R3)), int(round(np.linalg.det(R3)))
    if det > 0:
        if tr == 3:
            return "E"
        if tr == 0:
            return "C3"
        if tr == 1:
            return "C4"
        w, vecs = np.linalg.eigh(R3)
        ax = vecs[:, int(np.argmax(w))]
        return "C2p" if sum(abs(x) > 1e-9 for x in ax) == 1 else "C2"
    if tr == -3:
        return "i"
    if tr == 0:
        return "S6"
    if tr == -1:
        return "S4"
    w, vecs = np.linalg.eigh(R3)
    nm = vecs[:, int(np.argmin(w))]
    return "sh" if sum(abs(x) > 1e-9 for x in nm) == 1 else "sd"


CLS = [class_of(R3) for R3 in G48]
CHAR = {
    # O_h, classes E, 8C3 (C3), 6C2' edge axes (C2), 3C2 face axes (C2p = C4^2),
    # 6C4, i, 8S6, 6sigma_d (sd), 3sigma_h (sh), 6S4. The improper entries
    # follow from i C2p = sigma_h, i C2 = sigma_d, i C4 = S4, i C3 = S6, with
    # one extra minus for the u irreps. Get those entries wrong and the
    # 8-dimensional bond module decomposes into seven irreps of total dimension
    # 13 instead of the correct A1g + T2g + A2u + T1u of dimension 8.
    "A1g": dict(E=1, C3=1, C2=1, C2p=1, C4=1, i=1, S6=1, sd=1, sh=1, S4=1),
    "A2g": dict(E=1, C3=1, C2=-1, C2p=1, C4=-1, i=1, S6=1, sd=-1, sh=1, S4=-1),
    "Eg": dict(E=2, C3=-1, C2=0, C2p=2, C4=0, i=2, S6=-1, sd=0, sh=2, S4=0),
    "T1g": dict(E=3, C3=0, C2=-1, C2p=-1, C4=1, i=3, S6=0, sd=-1, sh=-1, S4=1),
    "T2g": dict(E=3, C3=0, C2=1, C2p=-1, C4=-1, i=3, S6=0, sd=1, sh=-1, S4=-1),
    "A1u": dict(E=1, C3=1, C2=1, C2p=1, C4=1, i=-1, S6=-1, sd=-1, sh=-1, S4=-1),
    "A2u": dict(E=1, C3=1, C2=-1, C2p=1, C4=-1, i=-1, S6=-1, sd=1, sh=-1, S4=1),
    "Eu": dict(E=2, C3=-1, C2=0, C2p=2, C4=0, i=-2, S6=1, sd=0, sh=-2, S4=0),
    "T1u": dict(E=3, C3=0, C2=-1, C2p=-1, C4=1, i=-3, S6=0, sd=1, sh=1, S4=-1),
    "T2u": dict(E=3, C3=0, C2=1, C2p=-1, C4=-1, i=-3, S6=0, sd=-1, sh=1, S4=1),
}


def projector(name):
    d = CHAR[name]["E"]
    P = np.zeros((NS, NS), dtype=complex)
    for k in range(48):
        P += CHAR[name][CLS[k]] * W48[k]
    return d * P / 48


BASIS = {}
for _n in CHAR:
    w, v = np.linalg.eigh(projector(_n).real)
    BASIS[_n] = v[:, w > 0.5]

NI = C.U(np.zeros(45)).shape[0]
W48_ON = [np.kron(W, np.eye(NI, dtype=complex)) for W in W48]


def internal_phases(th):
    """the 16 phases of the internal step. Computed exactly, not from a
    generator: C.U is a PRODUCT of exponentials of the 45 plane rotations, not
    the exponential of their sum, so the eigenphases of its sixth power are not
    linear in the angle. Reading them off a small-angle expansion and
    extrapolating is wrong by 4e-3 at theta = 2.36, which is the whole size of
    the gap being asked about."""
    return np.angle(np.linalg.eigvals(
        np.linalg.matrix_power(C.U(np.full(45, th)), C.ns)))


def wrap(x):
    return np.angle(np.exp(1j * x))


def P6(q, a, b):
    """the engine's six-stage bond product: angle a on intra-block stages,
    angle b on inter-block ones. At q = 0 this family is exactly
    O_h-covariant."""
    M = np.eye(NS, dtype=complex)
    for even, j in SC.STAGES:
        M = SC.stage_matrix(even, j, q[j], a if even else b) @ M
    return M


def step(a, b, th, q=GAM):
    """the full one-step matrix: bond factor times internal factor."""
    Vm = np.linalg.matrix_power(C.U(np.full(45, th)), C.ns)
    return np.kron(P6(q, a, b), Vm)


def bond_levels(a, b):
    """the phase P6(0) carries on A1g and on T1u; both isotypic blocks are
    scalar at Gamma, so the normalised trace is exact."""
    P = P6(GAM, a, b)
    return tuple(np.trace(BASIS[n].conj().T @ P @ BASIS[n])
                 / BASIS[n].shape[1] for n in ("A1g", "T1u"))


def pairing(a, b):
    """Im(lambda_A1g conj(lambda_T1u)): continuous, and zero exactly when the
    two bond phases coincide or are opposite; the sign of the real part tells
    the two apart. Searching the wrapped phase DIFFERENCE instead puts a
    root-finder on the discontinuity."""
    z1, z3 = bond_levels(a, b)
    return float(np.imag(z1 * np.conj(z3)))


def is_opposite(a, b):
    z1, z3 = bond_levels(a, b)
    return float(np.real(z1 * np.conj(z3))) < 0.0


def all_phases(a, b, th):
    """the 128 phases of the step at Gamma: the eigenvalues of a Kronecker
    product are the pairwise products, so the 8 bond phases and the 16 internal
    ones give them all without touching the 128 by 128 matrix."""
    bp = np.angle(np.linalg.eigvals(P6(GAM, a, b)))
    return wrap(bp[:, None] + internal_phases(th)[None, :]).ravel()


def four_smallest(a, b, th):
    """the objective of the search: the sum of the four smallest phase
    magnitudes. Zero exactly when four on-site states sit at phase zero, which
    is what a four-dimensional node at Gamma needs."""
    return float(np.sort(np.abs(all_phases(a, b, th)))[:4].sum())


def oh_content(Q):
    """the O_h content of a subspace with orthonormal basis Q, by the character
    inner product over the 48 elements."""
    chi = [np.trace(Q.conj().T @ W @ Q).real for W in W48_ON]
    out = {}
    for n, ch in CHAR.items():
        m = sum(ch[CLS[k]] * chi[k] for k in range(48)) / 48.0
        if m > 0.4:
            out[n] = int(round(m))
    return out


def cone_at(a, b, th, causal=True):
    """the node the six-stage product carries at Gamma, measured on the
    gapless eigenspace of the full one-step matrix."""
    M0 = step(a, b, th)
    w, V = np.linalg.eig(M0)
    ph = np.angle(w)
    gap = float(np.abs(wrap(ph)).min())
    # a window of 1e-6 catches near-degenerate neighbours and turns the
    # eigenspace into a non-invariant subspace; 1e-12 keeps it an eigenspace
    Q = V[:, np.abs(wrap(ph)) < 1e-12]
    # eig gives an arbitrary non-orthogonal basis for a degenerate eigenvalue
    # even of a unitary matrix; characters and tangent operators need an
    # orthonormal one
    Q, R = np.linalg.qr(Q)
    ng = int(np.sum(np.abs(np.diag(R)) > 1e-8))
    Q = Q[:, :ng]
    if ng == 0:
        return dict(gap=gap, gapless_onsite=0, onsite_dim=int(ph.size),
                    A=1.0, rho=0.0, v=0.0, chirality=0.0, content={},
                    front_speed=float("nan"), n_outside=-1)
    DQ = 1e-5
    Vm = np.linalg.matrix_power(C.U(np.full(45, th)), C.ns)
    S = []
    for k in range(3):
        e = np.zeros(3)
        e[k] = DQ
        dM = (np.kron(P6(e, a, b), Vm) - np.kron(P6(-e, a, b), Vm)) / (2 * DQ)
        S.append(1j * Q.conj().T @ M0.conj().T @ dM @ Q)
    nrm = [float(np.linalg.norm(s)) for s in S]
    com = max(float(np.linalg.norm(S[i] @ S[j] - S[j] @ S[i]))
              for i in range(3) for j in range(i + 1, 3))
    rho = com / (sum(n * n for n in nrm) / 3.0) if max(nrm) > 0 else 0.0
    sp = np.array([np.sort(np.abs(np.linalg.eigvalsh(
        d[0] * S[0] + d[1] * S[1] + d[2] * S[2]))) for d in DIRS])
    vmax = sp.max(axis=1)
    An = float((vmax.max() - vmax.min()) / vmax.max()) if vmax.max() > 0 else 1.0
    chi = (float(np.imag(np.trace(S[0] @ S[1] @ S[2])) /
                 (nrm[0] * nrm[1] * nrm[2])) if min(nrm) > 0 else 0.0)
    if not causal:
        return dict(gap=gap, gapless_onsite=ng, onsite_dim=int(ph.size),
                    A=An, rho=rho, v=float(vmax.max()), chirality=chi,
                    content=oh_content(Q), front_speed=float("nan"),
                    n_outside=-1)
    fr = A.front_speed_stages(C.c, [C.U(np.full(45, th))] * len(SC.STAGES),
                              [th] * len(SC.STAGES), SC.STAGES, Lb=6,
                              n_stages=12)
    return dict(gap=gap, gapless_onsite=ng, onsite_dim=int(ph.size),
                A=An, rho=rho, v=float(vmax.max()), chirality=chi,
                content=oh_content(Q),
                front_speed=float(fr["front_speed"]),
                n_outside=int(fr["n_outside_lightcone"]))


def main():
    t0 = time.time()
    out = {}

    # ---- (1) the family keeps Gamma covariant, and nothing else -----------
    def equi(qs, a, b):
        d = 0.0
        for R3, W in zip(G48, W48):
            for q in qs:
                qq = np.array(q)
                d = max(d, float(np.abs(W @ P6(qq, a, b) @ W.T
                                        - P6(R3 @ qq, a, b)).max()))
        return d

    eq = {"Gamma": equi([GAM], 0.9, 1.4),
          "M": equi([(np.pi / 2, np.pi / 2, 0.0)], 0.9, 1.4),
          "generic": equi([(0.3, -0.2, 0.7), (0.1, 0.4, -0.6)], 0.9, 1.4)}
    print("(1) two-angle family, O_h equivariance: " +
          "  ".join(f"{k}={v:.2e}" for k, v in eq.items()))
    out["equivariance"] = eq

    def P6a(angs, q):
        M = np.eye(NS, dtype=complex)
        for (even, j), g in zip(SC.STAGES, angs):
            M = SC.stage_matrix(even, j, q[j], g) @ M
        return M

    def defG(angs):
        P = P6a(angs, GAM)
        return max(float(np.abs(W @ P @ W.T - P).max()) for W in W48)

    rng = np.random.default_rng(1)
    gen = min(defG(rng.uniform(0.1, 3.0, 6)) for _ in range(200))
    peraxis = defG([0.3, 0.3, 0.9, 0.9, 1.7, 1.7])
    print(f"(1) not enlargeable: a generic six-angle product misses Gamma "
          f"covariance by {gen:.2f} (best of 200 draws), a per-axis "
          f"generalisation by {peraxis:.2f}")
    out["maximality"] = dict(generic_six_angle_best=gen, per_axis=peraxis)

    # ---- (2) the bond levels at Gamma: what has to pair -------------------
    print("\n(2) the bond levels at Gamma, at the angle the root uses "
          "(g = 0.2):")
    for (a, b) in ((0.2, 0.2), (0.9, 1.4)):
        bp = np.sort(np.angle(np.linalg.eigvals(P6(GAM, a, b))))
        z1, z3 = bond_levels(a, b)
        print(f"    a={a:.3f} b={b:.3f}  phases "
              f"{np.round(bp, 4).tolist()}  A1g {np.angle(z1):+.4f} "
              f"T1u {np.angle(z3):+.4f}  product "
              f"{z1 * np.conj(z3):+.4f}")

    # ---- (3) three angles, four cancellations ----------------------------
    # A four-state node at Gamma needs two bond levels and two internal levels
    # whose phases cancel in pairs: two conditions on three angles, so the
    # solution set, if it exists, is a curve. Grid, then polish.
    print("\n(3) searching (a, b, theta) for four simultaneous cancellations:")
    cand = []
    for a in np.linspace(0.05, np.pi - 0.05, 26):
        for b in np.linspace(0.05, np.pi - 0.05, 26):
            for th in np.linspace(0.02, np.pi - 0.02, 50):
                cand.append((four_smallest(a, b, th), a, b, th))
    cand.sort()
    print(f"    {len(cand)} grid points; the four best residuals "
          f"{[f'{c[0]:.3e}' for c in cand[:4]]}")
    pol = []
    for f0, a, b, th in cand[:60]:
        r = minimize(lambda v: four_smallest(v[0], v[1], v[2]), [a, b, th],
                     method="Nelder-Mead",
                     options=dict(xatol=1e-13, fatol=1e-16, maxiter=6000))
        if r.fun < 1e-10 and 0.01 < r.x[0] < np.pi - 0.01 \
                and 0.01 < r.x[1] < np.pi - 0.01 and 0.01 < r.x[2] < np.pi - 0.01:
            pol.append((float(r.fun), float(r.x[0]), float(r.x[1]),
                        float(r.x[2])))
    pol.sort()
    uniq = []
    for f, a, b, th in pol:
        if not any(abs(a - u[1]) < 1e-4 and abs(b - u[2]) < 1e-4
                   and abs(th - u[3]) < 1e-4 for u in uniq):
            uniq.append((f, a, b, th))
    print(f"    {len(pol)} polished points reach four exact cancellations; "
          f"{len(uniq)} distinct")
    out["polished"] = [dict(residual=f, a=a, b=b, theta=th)
                       for f, a, b, th in uniq]

    # ---- (4) what sits at those points -----------------------------------
    rows = []
    print(f"\n(4) {'a':>7} {'b':>7} {'theta':>8} {'gap':>9} {'gapl':>7} "
          f"{'A':>9} {'rho':>6} {'slope':>9} {'front':>6}  content")
    for i, (f, a, b, th) in enumerate(uniq):
        r = cone_at(a, b, th, causal=(i < 3))
        r.update(a=a, b=b, theta=th)
        r["v_planck_per_stage"] = 2.0 * r["v"] / len(SC.STAGES)
        rows.append(r)
        if i < 3 or i % 10 == 0:
            print(f"{a:7.4f} {b:7.4f} {th:8.5f} {r['gap']:9.2e} "
                  f"{r['gapless_onsite']:3d}/{r['onsite_dim']:<3d} "
                  f"{r['A']:9.2e} {r['rho']:6.3f} {r['v']:9.2e} "
                  f"{r['front_speed']:6.2f}  {r['content']}")
    out["rows"] = rows

    good = [r for r in rows if r["gap"] < 1e-8 and r["A"] < 1e-6
            and r["rho"] > 0.5 and r["gapless_onsite"] >= 4
            and r["n_outside"] == 0 and r["v"] > 1e-6
            and r["front_speed"] <= 1.0 + 1e-12]
    def single_parity(r):
        return len({c[-1] for c in r["content"]}) == 1

    spp = [r for r in rows if single_parity(r)]
    mxp = [r for r in rows if not single_parity(r)]
    flat = [r for r in rows if r["gap"] < 1e-8 and r["v"] < 1e-9]
    print(f"\n(5) conical nodes of the six-stage product: {len(good)}")
    print(f"    single-parity gapless spaces: {len(spp)}, largest slope "
          f"{max([r['v'] for r in spp], default=0.0):.3e} - flat, which is "
          f"what the parity no-go forces")
    print(f"    parity-mixed gapless spaces: {len(mxp)}, slope "
          f"{min([r['v'] for r in mxp], default=0.0):.3e}-"
          f"{max([r['v'] for r in mxp], default=0.0):.3e}, rho_comm at most "
          f"{max([r['rho'] for r in mxp], default=0.0):.3f}, anisotropy "
          f"{min([r['A'] for r in mxp], default=1.0):.3f}-"
          f"{max([r['A'] for r in mxp], default=1.0):.3f}")
    print(f"    gapless and flat: {len(flat)} of {len(rows)}")
    out["conical_nodes"] = len(good)
    out["flat_nodes"] = len(flat)
    out["single_parity_points"] = len(spp)
    out["parity_mixed_points"] = len(mxp)
    out["single_parity_max_slope"] = max([r["v"] for r in spp], default=0.0)
    out["parity_mixed_slope"] = [min([r["v"] for r in mxp], default=0.0),
                                 max([r["v"] for r in mxp], default=0.0)]
    out["parity_mixed_max_rho"] = max([r["rho"] for r in mxp], default=0.0)

    # ---- (6) the root's rule: ONE angle for bond and internal step -------
    diag = []
    dg = np.linspace(0.02, np.pi - 0.02, 600)

    def pairing_diag(x):
        return pairing(x, x)

    sp = np.array([pairing_diag(x) for x in dg])
    for i in range(len(dg) - 1):
        if sp[i] == 0 or sp[i] * sp[i + 1] < 0:
            xr = brentq(pairing_diag, dg[i], dg[i + 1], xtol=1e-14)
            if not is_opposite(xr, xr):
                continue
            d = cone_at(xr, xr, xr)
            d["angle"] = float(xr)
            diag.append(d)
    print(f"(6) single-angle rule: {len(diag)} angles where the A1g and T1u "
          f"bond phases are opposite on the diagonal a = b, measured with the "
          f"same angle on the internal step")
    for d in diag:
        print(f"    angle {d['angle']:.5f}  gap {d['gap']:9.2e}  "
              f"gapless {d['gapless_onsite']:4d}/{d['onsite_dim']}  "
              f"A {d['A']:9.2e}  rho {d['rho']:6.3f}  slope {d['v']:9.2e}  "
              f"{d['content']}")
    conical = [d for d in diag if d["gap"] < 1e-8 and d["A"] < 1e-6
               and d["v"] > 1e-6]
    print(f"    of those, conical: {len(conical)}")
    out["diagonal"] = diag

    # what this experiment claims is the no-go, not a cone: the engine's own
    # circuit carries no conical node, and its exactly gapless points split
    # exactly as the parity theorem forces - the single-parity ones flat, the
    # parity-mixed ones Abelian
    ok = (len(rows) > 0 and len(good) == 0 and len(conical) == 0
          and out["single_parity_max_slope"] < 1e-9
          and out["parity_mixed_max_rho"] < 1e-6)
    verdict = {
        "id": "V192_engine_two_angle_cone",
        "tier": 4,
        "property": "can the engine's own six-stage product carry a conical "
                    "node at Gamma: the axis-symmetric two-angle stage family "
                    "(a on the three intra-block stages, b on the three "
                    "inter-block ones, the largest six-stage family keeping "
                    "Gamma exactly O_h-covariant) searched jointly with the "
                    "internal angle for four simultaneous phase cancellations",
        "residual": None,
        "programme": "V192",
        "equivariance": eq,
        "grid_points": len(cand),
        "exact_four_cancellation_points": len(uniq),
        "conical_nodes": len(good),
        "flat_nodes": len(flat),
        "single_parity_points": len(spp),
        "parity_mixed_points": len(mxp),
        "best": (min(rows, key=lambda r: r["A"]) if rows else None),
        "largest_slope": (max(r["v"] for r in rows) if rows else 0.0),
        "diagonal_opposite_phases": len(diag),
        "single_angle_cone": len(conical),
        "single_parity_max_slope": out["single_parity_max_slope"],
        "parity_mixed_slope": out["parity_mixed_slope"],
        "parity_mixed_max_rho": out["parity_mixed_max_rho"],
        "max_step_gap": max((r["gap"] for r in rows), default=0.0),
        "max_polish_residual": max((f for f, _, _, _ in uniq), default=0.0),
        "n_directions": len(DIRS),
        "status": (f"{len(cand)} points of the (a, b, theta) grid searched, "
                   f"{len(uniq)} carry exactly gapless on-site states: "
                   f"{len(spp)} single-parity (flat to "
                   f"{out['single_parity_max_slope']:.1e}, as the parity "
                   f"no-go forces) and {len(mxp)} parity-mixed (slope "
                   f"{out['parity_mixed_slope'][0]:.2f}-"
                   f"{out['parity_mixed_slope'][1]:.2f} but Abelian, rho "
                   f"<= {out['parity_mixed_max_rho']:.1e}); conical nodes "
                   f"{len(good)}"),
        "pass": bool(ok),
        "data": out,
        "wall_time": time.time() - t0,
    }
    print(f"\nverdict: pass={ok}  conical nodes of the causal six-stage "
          f"product: {len(good)}  flat gapless points: {len(flat)}  "
          f"single-angle cone: {len(conical)}")
    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v192_engine_two_angle_cone.json"), "w") as fh:
        json.dump(verdict, fh, indent=2)
    print(f"wall {verdict['wall_time']:.2f} s")


if __name__ == "__main__":
    main()
