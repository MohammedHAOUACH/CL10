#!/usr/bin/env python
"""V191 - the parity no-go, and the Dirac cone it leaves open.

Every earlier programme in this line (V184, V188, V189, and the twelve W9
rounds that preceded this one) measured the same thing at every gapless point
it could find: the touching is FLAT.  The exponent of omega(q) about the node
is 2, the cone anisotropy A is O(1), and the tangent operators are numerically
zero.  Those programmes reported the measurement and treated the flatness as a
property of the search.  It is not.  It is a theorem about the cell, and the
theorem says exactly what has to change.

    (1) The 8 sub-sites of the 2x2x2 cell form the O_h module

            A1g (+) T2g (+) A2u (+) T1u        (multiplicity free)

        so every O_h-invariant bond Hamiltonian A(0) is a scalar on each of the
        four components, and every gapless multiplet at Gamma is a sum of them.

    (2) NO irrep of O_h admits a linear term.  The linear term sum_k q_k S_k in
        the effective Hamiltonian of a multiplet D exists iff
        Hom_O_h(T1u, D (x) D*) is non-zero - the three S_k must form a polar
        vector.  But D (x) D* is always EVEN under inversion, while q is ODD
        (T1u).  The multiplicity is therefore 0 for all ten irreps, measured
        below by the character inner product and again by direct linear algebra
        on the T1u block of the actual step (dim of the solution space: 0).
        This is why every scan found exponent 2: at a momentum whose little
        group contains inversion, a single-parity multiplet CANNOT be conical.

    (3) The escape is the one real Dirac semimetals use: the gapless space must
        be REDUCIBLE and mix parities.  The O_h products give the admissible
        pairs in this cell: A1g/T1u (A1g (x) T1u = T1u) and T1u/T2g
        (T1u (x) T2g = Eg (+) T1u (+) T2u).  The intertwiner is unique by
        Schur, hence the cone is isotropic by construction: for the strongest
        invariant hopping the coupling matrix obeys T T^dag = 16 I exactly
        (isotropy defect 4e-16).

    (4) Constructed and measured: the O_h-invariant bond Hamiltonian built from
        two hopping orbits (weights -24 and +8) puts A1g and T1u at the SAME
        level (-192) with T2g and A2u 384 away.  P(q) = exp(-i g A(q)) is
        exactly equivariant, P(Wq) = W P(q) W^T to 1.2e-14 over all 48
        elements.  With the internal angle tuned so one internal eigenphase
        cancels the common phase, the automaton is gapless at Gamma to 1e-11
        with a 4-dimensional multiplet A1g (+) T1u, cone anisotropy
        A = 0.00000 over the 215 directions, tangent algebra non-abelian
        (rho_comm = 1/sqrt 2), velocity v = g |T|.

    (5) What it costs.  The bond Hamiltonian of (4) is NOT the engine's step.
        The engine multiplies six stage matrices; the cubic group permutes the
        stages, and the product is not invariant under a permutation of
        non-commuting factors: the measured O_h defect of the implemented step
        is 1.07e-1, against 1.1e-15 for the stage SUM.  And with equal stage
        weights the stage-sum spectrum is pinned at four distinct levels
        (A1g +6, T1u +2, T2g -2, A2u -6), so the A1g/T1u degeneracy of (4) is
        unreachable inside the present stage structure.  A conical node
        therefore needs a bond layer redesigned around an O_h-invariant
        hopping whose stage factorisation is itself group-symmetric.

        One caveat carried by (4): P(q) = exp(-i g A(q)) is the exponential of
        a nearest-neighbour hopping, so it is exponentially local, not a
        finite-depth circuit, and its cone velocity 96 g is unbounded.  A
        range-one circuit realisation would additionally need 96 g <= 1 (the
        Lieb-Robinson bound on the front speed); the family allows it, the
        factorisation is the open piece.

The verdict is not "W9 is solved".  It is sharper: W9 is impossible for any
single-parity node, and it is solved - exactly, isotropically, with a
non-abelian tangent algebra - by a parity-degenerate node; the remaining work
is a stage factorisation of that bond layer, not another search.
"""
import importlib.util
import itertools
import json
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

V189 = f"{ROOT}/experiments/alt_step/V189_cell_step_joint_search/" \
       "V189_cell_step_joint_search.py"
_spec = importlib.util.spec_from_file_location("V189", V189)
V = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V)

C = V.Core((2, 2, 2), "baseline")
NS = C.nsub
GAM = np.zeros(3)
DQ = 1e-3
DIRS = C.dirs

# ---------------------------------------------------------------- the group --
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
        x = np.array([2 * C.c["digit"](e, j) - 1.0 for j in range(3)])
        y = R3 @ x
        out.append(C.c["site_index"]((int(round((y[0] + 1) / 2)),
                                      int(round((y[1] + 1) / 2)),
                                      int(round((y[2] + 1) / 2)))))
    return np.array(out)


PERM = [_perm_of(R3) for R3 in G48]
W48 = []
for P in PERM:
    W = np.zeros((NS, NS))
    for e in range(NS):
        W[P[e], e] = 1.0
    W48.append(W)


def class_of(R3):
    """the O_h class of a signed permutation, from trace, det and axis."""
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


def wrap(x):
    return np.angle(np.exp(1j * x))


def main():
    t0 = time.time()
    out = {}

    # ---- (1) the site module ------------------------------------------------
    chi_site = [int(np.trace(W)) for W in W48]
    mult = {n: sum(CHAR[n][CLS[k]] * chi_site[k] for k in range(48)) / 48
            for n in CHAR}
    present = {n: m for n, m in mult.items() if abs(m) > 1e-9}
    proj = {}
    idem = 0.0
    for n in present:
        proj[n] = projector(n)
        idem = max(idem, float(np.abs(proj[n] @ proj[n] - proj[n]).max()))
    print(f"(1) site module of the 2x2x2 cell: "
          f"{' + '.join(f'{n} x{int(round(m))}' for n, m in present.items())}")
    print(f"    projector idempotence {idem:.2e}; "
          f"ranks {[int(np.linalg.matrix_rank(proj[n].real)) for n in present]}")
    out["site_module"] = {n: int(round(m)) for n, m in present.items()}
    out["projector_idempotence"] = idem

    # ---- (2) the parity no-go ----------------------------------------------
    # Hom_O_h(T1u, D (x) D*) = <chi_{D (x) D*}, chi_T1u> = sum chi_D^2 chi_T1u / 48
    hom = {n: sum(CHAR[n][CLS[k]] ** 2 * CHAR["T1u"][CLS[k]] for k in range(48))
           / 48 for n in CHAR}
    print("(2) dim Hom_O_h(T1u, D x D*) for every irrep D:")
    print("    " + "  ".join(f"{n}:{int(round(v))}" for n, v in hom.items()))
    out["linear_term_multiplicity"] = {n: int(round(v)) for n, v in hom.items()}

    # the same statement by direct linear algebra on the T1u block of the step
    def A_stage_sum(q):
        tot = np.zeros((NS, NS), dtype=np.complex128)
        for bonds in C.bonds:
            for (i, j, v) in bonds:
                ph = np.exp(1j * (q[0] * v[0] + q[1] * v[1] + q[2] * v[2]))
                tot[i, j] += ph
                tot[j, i] += np.conj(ph)
        return tot

    ev0, evec0 = np.linalg.eigh(A_stage_sum(GAM))
    idx = [k for k in range(NS) if abs(ev0[k] - 2.0) < 1e-9]
    Qt, _ = np.linalg.qr(evec0[:, idx])
    Dl = [Qt.conj().T @ W @ Qt for W in W48]
    d = 3
    rows = []
    for R3, D in zip(G48, Dl):
        kr = np.kron(D, D.conj())
        for j in range(d):
            M2 = np.zeros((d * d, 3 * d * d), dtype=complex)
            for m in range(d):
                for n in range(d):
                    for mp in range(d):
                        for npp in range(d):
                            M2[m * d + n, j * d * d + mp * d + npp] += \
                                kr[m * d + n, mp * d + npp]
            M = np.zeros((d * d, 3 * d * d), dtype=complex)
            for i in range(d):
                if R3[i, j] != 0:
                    for m in range(d):
                        for n in range(d):
                            M[m * d + n, i * d * d + m * d + n] -= R3[i, j]
            rows.append(M2 - M)
    sv = np.linalg.svd(np.vstack(rows), compute_uv=False)
    dim_lin = int(sum(s < 1e-10 for s in sv))
    print(f"    direct solve on the T1u block of the stage sum: "
          f"dim {dim_lin} (smallest sigma {sv[-1]:.3f})")
    out["linear_term_dim_direct"] = dim_lin

    # ---- (3) the invariant hopping family and the admissible pairs ----------
    terms = [(i, j, v) for i in range(NS) for j in range(NS)
             for v in itertools.product((-1, 0, 1), repeat=3)]

    def key(t):
        i, j, v = t
        best = None
        for R3, P in zip(G48, PERM):
            vv = tuple(int(round(x)) for x in R3 @ np.array(v))
            k = (int(P[i]), int(P[j]), vv)
            if best is None or k < best:
                best = k
        return best

    orbs = {}
    for t in terms:
        orbs.setdefault(key(t), []).append(t)
    orb_list = sorted(orbs.values(), key=lambda m: key(m[0]))

    def bloch(mem, q):
        M = np.zeros((NS, NS), dtype=complex)
        for (i, j, v) in mem:
            M[i, j] += np.exp(1j * np.dot(q, v))
        return M + M.conj().T

    comp = list(present)
    basis = {}
    for n in comp:
        w, v = np.linalg.eigh(proj[n].real)
        basis[n] = v[:, w > 0.5]
    lam = {n: np.array([np.trace(basis[n].conj().T @ bloch(m, GAM) @ basis[n]).real
                        / basis[n].shape[1] for m in orb_list]) for n in comp}
    print(f"(3) {len(orb_list)} O_h-invariant hopping orbits "
          f"(displacements |v| <= 1); level carried per component")

    pairs = [("A1g", "T1u"), ("T1u", "T2g")]
    ncou = {}
    for (a, b) in pairs:
        cnt = 0
        best = (0.0, None, None)
        for oi, mem in enumerate(orb_list):
            dA = [np.zeros((NS, NS), dtype=complex) for _ in range(3)]
            for (i, j, v) in mem:
                for k in range(3):
                    dA[k][i, j] += 1j * v[k]
            T = np.stack([basis[a].conj().T @ dd @ basis[b] for dd in dA])
            nrm = float(np.sqrt((np.abs(T) ** 2).sum()))
            if nrm > 1e-12:
                cnt += 1
            if nrm > best[0]:
                Tm = T.reshape(3, -1)          # rows = direction k
                best = (nrm, oi, np.linalg.eigvalsh(Tm @ Tm.conj().T))
        iso = float(np.abs(best[2] - best[2].mean()).max() / best[2].mean())
        print(f"    {a}/{b}: {cnt} orbits couple linearly; strongest orbit "
              f"{best[1]} TT^dag eigenvalues {np.round(best[2], 6)}, "
              f"isotropy defect {iso:.2e}")
        ncou[f"{a}/{b}"] = {"orbits": cnt, "isotropy_defect": iso,
                            "TTdag": [float(x) for x in best[2]]}
    out["linear_coupling"] = ncou

    # ---- (4) the construction -----------------------------------------------
    best = None
    for a in range(len(orb_list)):
        for b in range(a + 1, len(orb_list)):
            dd = lam["T1u"][a] - lam["A1g"][a]
            ee = lam["T1u"][b] - lam["A1g"][b]
            if abs(dd - ee) < 1e-9:
                continue
            w = np.zeros(len(orb_list))
            w[a], w[b] = ee, -dd
            lv = {n: float(w @ lam[n]) for n in comp}
            sep = min(abs(lv[n] - lv["A1g"]) for n in comp
                      if n not in ("A1g", "T1u"))
            if sep < 0.5:
                continue
            dA = [np.zeros((NS, NS), dtype=complex) for _ in range(3)]
            for k in range(3):
                for x in (a, b):
                    for (i, j, vv) in orb_list[x]:
                        dA[k][i, j] += w[x] * 1j * vv[k]
            T = np.stack([basis["A1g"].conj().T @ dd2 @ basis["T1u"]
                          for dd2 in dA]).reshape(3, 3)
            ev = np.linalg.eigvalsh(T @ T.conj().T)
            if ev.mean() <= 1e-6:
                continue
            iso = float(np.abs(ev - ev.mean()).max() / ev.mean())
            if ev.mean() > 1e-6 and (best is None
                                     or ev.mean() / (1 + 10 * iso) > best[0]):
                best = (float(ev.mean() / (1 + 10 * iso)), w, lv, iso,
                        float(np.sqrt(ev.mean())), (a, b))
    score, w, lv, iso, tscale, pair = best

    def Aq(q):
        o = np.zeros((NS, NS), dtype=complex)
        for x in np.nonzero(np.abs(w) > 0)[0]:
            o += w[x] * bloch(orb_list[x], q)
        return o

    def P_sym(q, g):
        return expm(-1j * g * Aq(q))

    eqd = 0.0
    for R3, W in zip(G48, W48):
        for q in ([0.3, -0.2, 0.7], [0.9, 0.9, 0.1], [0.1, 0.4, -0.6]):
            qq = np.array(q)
            eqd = max(eqd, float(np.abs(W @ P_sym(qq, 0.7) @ W.T
                                        - P_sym(R3 @ qq, 0.7)).max()))
    print(f"(4) construction: orbits {pair} weights "
          f"{w[pair[0]]:+.1f}, {w[pair[1]]:+.1f}; levels "
          + " ".join(f"{n}={lv[n]:+.0f}" for n in comp))
    print(f"    equivariance P(Wq) = W P(q) W^T defect {eqd:.2e}; "
          f"coupling |T| {tscale:.1f}, isotropy defect {iso:.2e}")

    rows = []
    for g in (0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0):
        target = g * lv["A1g"]

        def mismatch(th):
            phv = np.angle(np.linalg.eigvals(np.linalg.matrix_power(
                C.U(np.full(45, th)), C.ns)))
            return min(abs(wrap(target - p)) for p in phv)
        grid = np.linspace(0.01, np.pi - 0.01, 1500)
        vals = np.array([mismatch(t) for t in grid])
        k = int(np.argmin(vals))
        r = minimize_scalar(mismatch, bracket=(grid[max(0, k - 1)], grid[k],
                                               grid[min(len(grid) - 1, k + 1)]),
                            method="brent", options={"xtol": 1e-15})
        th = float(r.x)
        Vm = np.linalg.matrix_power(C.U(np.full(45, th)), C.ns)
        phv = np.angle(np.linalg.eigvals(Vm))
        php = np.angle(np.linalg.eigvals(P_sym(GAM, g)))
        gap = float(np.abs(wrap(php[:, None] + phv[None, :])).min())

        Qg = np.hstack([basis["A1g"], basis["T1u"]])
        P0 = P_sym(GAM, g)
        S = []
        for k2 in range(3):
            e = np.zeros(3)
            e[k2] = DQ
            dP = (P_sym(GAM + e, g) - P_sym(GAM - e, g)) / (2 * DQ)
            S.append(1j * Qg.conj().T @ P0.conj().T @ dP @ Qg)
        herm = float(max(np.abs(s - s.conj().T).max() for s in S))
        nrm = [float(np.linalg.norm(s)) for s in S]
        com = max(float(np.linalg.norm(S[i] @ S[j] - S[j] @ S[i]))
                  for i in range(3) for j in range(i + 1, 3))
        rho = com / (sum(n * n for n in nrm) / 3.0)
        sp = np.array([np.sort(np.abs(np.linalg.eigvalsh(
            dd[0] * S[0] + dd[1] * S[1] + dd[2] * S[2]))) for dd in DIRS])
        vmax = sp.max(axis=1)
        A = float((vmax.max() - vmax.min()) / vmax.max())
        # the in-band spectrum: {0, 0, +-v|q|} for a spin-1 Dirac pair
        band = sp.min(axis=1)
        rows.append(dict(g=g, theta=th, gap=gap, mult=int(Qg.shape[1]),
                         module="A1g+T1u", A=A, rho_comm=rho,
                         v=float(vmax.max()), flat_band=float(np.abs(band).max()),
                         hermiticity=herm))
        print(f"    g={g:5.2f} theta={th:8.5f} gap={gap:8.1e} mult={Qg.shape[1]} "
              f"A={A:8.5f} rho={rho:6.3f} v={vmax.max():8.4f} "
              f"flat={np.abs(band).max():.1e}")
    out["construction"] = dict(orbits=list(pair), weights=[float(w[pair[0]]),
                                                           float(w[pair[1]])],
                               levels={n: lv[n] for n in comp},
                               equivariance_defect=eqd, T=tscale,
                               isotropy_defect=iso, rows=rows)

    # ---- (5) the cost: the engine's own step -------------------------------
    def P_engine(q):
        return V.SC.bond_product(q)   # the engine's own bond product

    MOM = {"Gamma": [(0.0, 0.0, 0.0)],
           "M": [(np.pi / 2, np.pi / 2, 0.0)],
           "R": [(np.pi / 2, np.pi / 2, np.pi / 2)],
           "generic": [(0.3, -0.2, 0.7), (0.9, 0.9, 0.1), (0.1, 0.4, -0.6)]}

    def defect(fn, qs):
        d = 0.0
        for R3, W in zip(G48, W48):
            for q in qs:
                qq = np.array(q)
                d = max(d, float(np.abs(W @ fn(qq) @ W.T - fn(R3 @ qq)).max()))
        return d

    ed = {k: defect(P_engine, v) for k, v in MOM.items()}
    sd = {k: defect(lambda q: expm(-1j * 0.7 * A_stage_sum(q)), v)
          for k, v in MOM.items()}
    lv_stage = {}
    for n in comp:
        B = basis[n]
        lv_stage[n] = float(np.trace(B.conj().T @ A_stage_sum(GAM) @ B).real
                            / B.shape[1])
    print("(5) implemented 6-stage product, O_h defect: " +
          "  ".join(f"{k}={ed[k]:.2e}" for k in MOM))
    print("    stage SUM, same momenta: " +
          "  ".join(f"{k}={sd[k]:.2e}" for k in MOM))
    print("    stage-sum levels: " + " ".join(f"{n}={v:+.3f}"
                                              for n, v in lv_stage.items())
          + f"  -> A1g/T1u splitting {abs(lv_stage['A1g'] - lv_stage['T1u']):.3f}")
    out["engine"] = dict(product_o_h_defect=ed, sum_o_h_defect=sd,
                         stage_levels=lv_stage)

    ok = (all(v == 0 for v in out["linear_term_multiplicity"].values())
          and out["linear_term_dim_direct"] == 0
          and all(r["A"] < 1e-6 and r["gap"] < 1e-8 and r["mult"] == 4
                  for r in rows))
    verdict = {
        "id": "V191_parity_dirac_cone",
        "tier": 4,
        "property": "why every gapless node of the 2x2x2 cell found so far "
                    "is flat, and the node that is not: the linear-term "
                    "multiplicity Hom_O_h(T1u, D (x) D*) for all ten "
                    "irreps at an inversion-invariant momentum, and the "
                    "cone on a parity-degenerate A1g+T1u pair",
        "residual": None,
        "programme": "V191",
        "question": "can a gapless node of the 2x2x2 cell be conical?",
        "no_go": "no linear term exists for any single irrep of O_h at a "
                 "momentum whose little group contains inversion: "
                 "Hom(T1u, D x D*) = 0 for all ten irreps",
        "construction": "a parity-degenerate A1g+T1u pair at Gamma gives an "
                        "exactly gapless, exactly isotropic, non-abelian cone",
        "cone_anisotropy": min(r["A"] for r in rows),
        "gap": min(r["gap"] for r in rows),
        "rho_comm": rows[0]["rho_comm"],
        "velocity": [r["v"] for r in rows],
        "engine_o_h_defect": ed,
        "status": (f"parity no-go: linear-term multiplicity 0 for all ten "
                   f"irreps (direct solve on the T1u block: dim "
                   f"{out['linear_term_dim_direct']}); the parity-degenerate "
                   f"A1g+T1u pair gives an isolated cone - min gap "
                   f"{min(r['gap'] for r in rows):.1e}, min anisotropy "
                   f"{min(r['A'] for r in rows):.1e}, rho_comm "
                   f"{rows[0]['rho_comm']:.3f}; the implemented 6-stage "
                   f"product is O_h-equivariant at Gamma to "
                   f"{ed['Gamma']:.1e} and broken elsewhere "
                   f"(M {ed['M']:.1e}, R {ed['R']:.1e}, generic "
                   f"{ed['generic']:.1e})"),
        "pass": bool(ok),
        "data": out,
        "wall_time": time.time() - t0,
    }
    print(f"\nverdict: pass={ok}  min A = {verdict['cone_anisotropy']:.6f}  "
          f"min gap = {verdict['gap']:.2e}  "
          f"engine O_h defect: Gamma {ed['Gamma']:.1e}, generic {ed['generic']:.1e}")
    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v191_parity_dirac_cone.json"), "w") as fh:
        json.dump(verdict, fh, indent=2)
    print(f"wall {verdict['wall_time']:.2f} s")


if __name__ == "__main__":
    main()
