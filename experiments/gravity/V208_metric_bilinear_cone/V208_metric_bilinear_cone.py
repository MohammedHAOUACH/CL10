#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V208 - the metric is a bilinear of the bond layer, and its cone is a module fact.

The paper's Section 11 takes a symmetric rank-2 field and G = 1 as inputs. This
experiment settles, on the source's own representation content, how much of that
input the bond layer actually forces, and how much it does not. Three numbers
define the gravitational sector - the field content, the light cone, Newton's
constant - and they are treated separately here.

**T1 the content is forced, and it is a bilinear.** The polarisation space of a
symmetric tensor on the cubic lattice is Sym^2(R^3), whose O_h content is
A1g + Eg + T2g (the trace plus the five-component lattice spin-2 quintet). The
carrier the automaton actually has - one amplitude on each of the six oriented
bonds - is a different 6-dimensional module: A1g + Eg + T1u. The quintet's T2g
is absent and an odd T1u sits in its place: **a single bond amplitude cannot be
the metric, by characters, exactly.** The metric content first appears in the
bilinear module bond (x) Sym^2 (36 states), whose O_h content and invariant
algebra are computed here: the invariant algebra has dimension 34 = sum of the
squared multiplicities, i.e. 34 independent gravitational self-couplings are
allowed and none is chosen.

**T2 the hop cannot move the metric, by parity.** The velocity of the hop (the
derivative of the bond phase at q = 0) is a POLAR vector operator, the spin-2
dressing an AXIAL one; both checks are exact. The metric sector is the even half
of the bilinear module (18 states), and the hop's velocity block inside it is
zero to 2.7e-16: the hop does not move the metric at first order. What is left
in the even sector is the dressing, whose spectrum is g x {0, +-1, +-2} in every
direction (the spin-2 spectrum theorem) - isotropic exactly, with cone radius
2g. This is where the implemented lattice graviton's speed comes from.

**T3 an exactly lightlike cone does exist in the module, with the quintet on
it.** The space of O-covariant Hermitian velocity triples (H_1,H_2,H_3) on the
bilinear module is 160-dimensional. Inside it sits a triple, committed in
cone_generator.json and verified here rather than searched at run time, whose
acoustic symbol is exactly Clifford-regular: {H_i,H_j} = 2 delta_ij S with
S = H_1^2 = H_2^2 = H_3^2, defect 1.4e-13. Its 36-branch spectrum is identical
in all 44 sampled directions to 5e-15; the cone carries 8 branches (content
E + T1 + T2: the metric quintet appaired with the odd spin-connection vector -
the first-order, Palatini pairing) at slope exactly 1, and 4 branches at slope
exactly 0 (the constrained pieces, including the trace A1g, which does not
propagate). The step M(q) = exp(i sum sin(q_k) H_k) is unitary to 1e-14 and
M(0) = I exactly: within this class no mass term for the metric can be written,
which is the emergent gauge invariance, as a number. The step is quasi-local:
the real-space kernel falls from 9.1e-1 on site to 1.1e-1 at one bond, 3.4e-2,
5.4e-3, 7.9e-4, 7.2e-5, and below 7e-6 beyond five bonds.

**T4 the implemented step is not that triple.** The committed dressing of the
directional layer (V140/V141) explains 0.9% of the cone triple, and the hop
0.3%; the exact-cone family is 53-dimensional (Jacobian rank 107 of 160) and
none of its members is the implemented generator. So the light cone of the
lattice graviton is not forced by the module - it is one of a 53-parameter
family that the dynamics does not pick. The implemented generator gives
c_gw = 2 g = 1.5708 in units of the matter cone (V141's committed
dressing_only_speed), against |c_gw/c - 1| < 1e-15 from GW170817; the clock that
would give exactly 1 (g = 1/2) is inside V206's dead zone (n_star = 3, deficit
1.44e-3, no crossing of the closed class).

**T5 Newton's constant, as a number of the root.** The bond entropy is
committed: S = ln 45 = 3.8066624897703196 nats per Planck bond (V129, V101,
V181). Identifying that bond entropy with the gravitational entropy of one
Planck area, S = A/4G, gives G = 1/(4 ln 45) = 6.5667e-2 in cell units,
16 pi G = 3.3012, and a cell size a / l_P = 2 sqrt(ln 45) = 3.9022. The
Bekenstein gap the paper already names, 4 ln 45 = 15.2266, is exactly the ratio
by which the naive bond count overshoots A/4 at G = 1. The identification is the
input; the number is not.
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


def art(area, exp, key):
    pth = os.path.join(ROOT, "experiments", area, exp, "results",
                       "%s.json" % exp.lower())
    with open(pth) as fh:
        j = json.load(fh)
    for k in key.split("."):
        j = j[k]
    return j


# ---------------------------------------------------------------- group theory
def build_Oh():
    OHS = [(p, s) for p in itertools.permutations(range(3))
           for s in itertools.product((1, -1), repeat=3)]
    Rv = []
    for p, s in OHS:
        R = np.zeros((3, 3))
        for i in range(3):
            R[i, p[i]] = s[i]
        Rv.append(R)
    return OHS, Rv


CLASSES = {  # class key -> (name, size, characters of the ten O_h irreps)
    ((3., 3., 1., True)): ("E", 1, [1, 1, 2, 3, 3, 1, 1, 2, 3, 3]),
    ((0., 0., 1., False)): ("8C3", 8, [1, 1, -1, 0, 0, 1, 1, -1, 0, 0]),
    ((-1., 3., 1., True)): ("3C2", 3, [1, 1, 2, -1, -1, 1, 1, 2, -1, -1]),
    ((1., -1., 1., False)): ("6C4", 6, [1, -1, 0, 1, -1, 1, -1, 0, 1, -1]),
    ((-1., 3., 1., False)): ("6C2p", 6, [1, -1, 0, -1, 1, 1, -1, 0, -1, 1]),
    ((-3., 3., -1., True)): ("i", 1, [1, 1, 2, 3, 3, -1, -1, -2, -3, -3]),
    ((0., 0., -1., False)): ("8S6", 8, [1, 1, -1, 0, 0, -1, -1, 2, 0, 0]),
    ((1., 3., -1., True)): ("3sh", 3, [1, 1, 2, -1, -1, -1, -1, -2, 1, 1]),
    ((-1., -1., -1., False)): ("6S4", 6, [1, -1, 0, 1, -1, -1, 1, 0, -1, 1]),
    ((1., 3., -1., False)): ("6sd", 6, [1, -1, 0, -1, 1, -1, 1, 0, 1, -1])}
IRREPS = ["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"]


def ckey(R):
    k = (round(np.trace(R).real, 9), round(np.trace(R @ R).real, 9),
         round(np.linalg.det(R).real, 9))
    diag = bool(np.max(np.abs(R - np.diag(np.diag(R)))) < 1e-9)
    return k + (diag,)


def decompose(chi, keys):
    out = {}
    for i, r in enumerate(IRREPS):
        m = sum(CLASSES[k][2][i] * chi[j] for j, k in enumerate(keys)) / 48.
        if abs(m) > 0.4:
            out[r] = int(round(m))
    return out


def main():
    t0 = time.time()
    out = {"experiment": "V208_metric_bilinear_cone", "data": {}}
    checks = []

    # ---- the module: 6 oriented bonds x Sym^2(R^3) -------------------------
    OHS, Rv = build_Oh()
    keys = [ckey(R) for R in Rv]
    IDX = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    SB = []
    for k, l in IDX:
        S = np.zeros((3, 3))
        S[k, l] += 1.0
        if l != k:
            S[l, k] += 1.0
        SB.append(S / math.sqrt(2.0 if l != k else 1.0))
    DIRS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    Dp, Db = [], []
    for R in Rv:
        Mm = np.zeros((6, 6))
        for a in range(6):
            T = R @ SB[a] @ R.T
            for b in range(6):
                Mm[b, a] = np.trace(SB[b].T @ T)
        Dp.append(Mm)
        Pm = np.zeros((6, 6))
        for j, d in enumerate(DIRS):
            Rd = tuple(int(round(x)) for x in R @ np.array(d))
            Pm[DIRS.index(Rd), j] = 1.0
        Db.append(Pm)
    Dm = [np.kron(Pb, Pp) for Pb, Pp in zip(Db, Dp)]

    chi_bond = [np.trace(P) for P in Db]
    chi_pol = [np.trace(P) for P in Dp]
    chi_mod = [a * b for a, b in zip(chi_bond, chi_pol)]
    d_bond = decompose(chi_bond, keys)
    d_pol = decompose(chi_pol, keys)
    d_mod = decompose(chi_mod, keys)
    comm_dim = sum(c * c for c in d_mod.values())
    sizes_ok = all(CLASSES[k][1] == keys.count(k) for k in CLASSES)
    t1ok = (sizes_ok and d_bond == {"A1g": 1, "Eg": 1, "T1u": 1}
            and d_pol == {"A1g": 1, "Eg": 1, "T2g": 1}
            and d_mod.get("Eg") == 3 and d_mod.get("T2g") == 2
            and d_bond.get("T2g", 0) == 0 and d_pol.get("T1u", 0) == 0
            and comm_dim == 34)
    out["data"]["content"] = {
        "bond_module_6": d_bond, "polarisation_Sym2": d_pol,
        "bilinear_module_36": d_mod, "invariant_algebra_dim": comm_dim,
        "single_bond_has_T2g": int(d_bond.get("T2g", 0)),
        "polarisation_has_T1u": int(d_pol.get("T1u", 0))}
    checks.append(("T1 content: bond(6 oriented) = %s while Sym^2 = %s - the "
                   "quintet's T2g is absent from a single bond (%d copies) and "
                   "T1u is absent from Sym^2 (%d), so the metric is a bilinear: "
                   "bond (x) Sym^2 = %s, whose invariant algebra has dimension "
                   "%d (= sum m^2: %d independent self-couplings, none chosen)"
                   % (d_bond, d_pol, d_bond.get("T2g", 0), d_pol.get("T1u", 0),
                      d_mod, comm_dim, comm_dim), t1ok))

    # ---- the two velocity triples and the parity wall ----------------------
    lev = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
           (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
    Bb = [np.zeros((6, 6)) for _ in range(3)]
    for k in range(3):
        for j, d in enumerate(DIRS):
            Bb[k][j, j] = d[k]
    Kp = []
    for k in range(3):
        Aa = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                Aa[i, j] = -lev.get((k, i, j), 0)
        Mm = np.zeros((6, 6), dtype=complex)
        for b in range(6):
            act = 1j * (Aa @ SB[b] - SB[b] @ Aa)
            for a in range(6):
                Mm[a, b] = np.trace(SB[a].T @ act)
        Kp.append(Mm)
    Bop = [np.kron(Bb[k], np.eye(6)) for k in range(3)]
    Kop = [np.kron(np.eye(6), Kp[k]) for k in range(3)]

    viol_p = 0.0
    viol_a = 0.0
    for R, D in zip(Rv, Dm):
        det = round(np.linalg.det(R).real)
        for a in range(3):
            lhs = D @ Bop[a] @ D.conj().T
            rhs = sum(R[b, a] * Bop[b] for b in range(3))
            viol_p = max(viol_p, np.abs(lhs - rhs).max())
            lhs = D @ Kop[a] @ D.conj().T
            rhs = sum(det * R[b, a] * Kop[b] for b in range(3))
            viol_a = max(viol_a, np.abs(lhs - rhs).max())
    iota = [i for i, R in enumerate(Rv) if np.allclose(R, -np.eye(3))][0]
    Pe = (np.eye(36) + Dm[iota]) / 2.0
    hop_even = [float(np.linalg.norm(Pe @ Bop[k] @ Pe)) for k in range(3)]
    drs_even = [float(np.linalg.norm(Pe @ Kop[k] @ Pe)) for k in range(3)]

    dirs44 = [np.array(v, dtype=float) for v in
              [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1),
               (1, 1, 1), (1, 1, -1), (1, 2, 3), (1, 2, -2), (0.31, -0.22, 0.77),
               (0.9, 0.1, -0.3), (-0.4, 0.75, 0.5), (0.05, 0.99, 0.12),
               (1, 3, 1), (2, 2, 1), (5, 1, 2), (1, 1, 4), (0.6, -0.8, 0.1),
               (0.2, 0.3, 0.93)]]
    rng = np.random.default_rng(2024)
    for _ in range(24):
        dirs44.append(rng.standard_normal(3))
    dirs44 = [n / np.linalg.norm(n) for n in dirs44]

    g_root = math.pi / 4.0
    drs_spread = 0.0
    for n in dirs44:
        ev = np.sort(np.linalg.eigvalsh(sum(n[k] * Kop[k] for k in range(3))))
        drs_spread = max(drs_spread, float(np.abs(ev - np.array(
            [-2.] * 6 + [-1.] * 6 + [0.] * 12 + [1.] * 6 + [2.] * 6)).max()))
    t2ok = (viol_p < 1e-13 and viol_a < 1e-13 and max(hop_even) < 1e-13
            and min(drs_even) > 1.0 and drs_spread < 1e-12
            and round(np.trace(Pe).real) == 18)
    out["data"]["parity"] = {
        "hop_polar_vector_violation": viol_p,
        "dressing_axial_vector_violation": viol_a,
        "even_sector_rank": int(round(np.trace(Pe).real)),
        "hop_velocity_block_in_even_sector": hop_even,
        "dressing_block_in_even_sector": drs_even,
        "dressing_spectrum_spread_over_directions": drs_spread,
        "implemented_cone_radius_2g": 2.0 * g_root}
    checks.append(("T2 parity: the hop velocity is a polar vector (covariance "
                   "%.1e) and the dressing an axial one (%.1e); the even "
                   "(metric) half of the bilinear module has rank %d and the "
                   "hop's velocity block inside it vanishes (%s), while the "
                   "dressing survives (%s); the dressing's spectrum is exactly "
                   "g x {0,+-1,+-2} in all %d directions (max deviation %.1e), "
                   "so the implemented cone radius is 2g = %.6f"
                   % (viol_p, viol_a, int(round(np.trace(Pe).real)),
                      ["%.1e" % v for v in hop_even],
                      ["%.3f" % v for v in drs_even], len(dirs44),
                      drs_spread, 2 * g_root), t2ok))

    # ---- the committed exact-cone generator --------------------------------
    gen = json.load(open(os.path.join(HERE, "cone_generator.json")))
    H = [np.array(gen["A"][k]) + 1j * np.array(gen["S"][k]) for k in range(3)]
    herm = max(np.abs(h - h.conj().T).max() for h in H)
    PROP = [i for i, R in enumerate(Rv) if np.linalg.det(R) > 0]
    cov = 0.0
    for i in PROP:
        R, D = Rv[i], Dm[i]
        for a in range(3):
            lhs = D @ H[a] @ D.conj().T
            rhs = sum(R[b, a] * H[b] for b in range(3))
            cov = max(cov, np.abs(lhs - rhs).max())
    S3 = (sum(h @ h for h in H) / 3.0)
    S3 = (S3 + S3.conj().T) / 2.0
    def1 = [H[i] @ H[j] + H[j] @ H[i] - (2 * S3 if i == j else 0.0)
            for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]]
    cliff = float(np.linalg.norm(np.concatenate(
        [np.concatenate([z.real.reshape(-1), z.imag.reshape(-1)])
         for z in def1])))
    bands = []
    for n in dirs44:
        bands.append(np.sort(np.linalg.eigvalsh(sum(n[k] * H[k]
                                                   for k in range(3)))))
    bands = np.array(bands)
    bs = float(np.abs(bands - bands[0]).max())
    w0, U0 = np.linalg.eigh(S3)
    top = float(w0[-1])
    cone_sel = [i for i in range(36) if abs(w0[i] - top) < 1e-8]
    zero_sel = [i for i in range(36) if abs(w0[i]) < 1e-9]
    # content is measured in the chiral group O (the triple is covariant under
    # the 24 proper rotations; the directed network is chiral by construction)
    CH24 = {"A1": [1, 1, 1, 1, 1], "A2": [1, 1, 1, -1, -1],
            "E": [2, -1, 2, 0, 0], "T1": [3, 0, -1, 1, -1],
            "T2": [3, 0, -1, -1, 1]}
    def cls24(R):
        t = round(np.trace(R).real, 9)
        t2 = round(np.trace(R @ R).real, 9)
        di = bool(np.max(np.abs(R - np.diag(np.diag(R)))) < 1e-9)
        if abs(t - 3) < 1e-9:
            return 0
        if abs(t) < 1e-9:
            return 1
        if abs(t + 1) < 1e-9 and abs(t2 - 3) < 1e-9 and di:
            return 2
        if abs(t - 1) < 1e-9 and abs(t2 + 1) < 1e-9:
            return 3
        return 4
    lab24 = [cls24(Rv[i]) for i in PROP]
    D24 = [Dm[i] for i in PROP]
    def mults(vals):
        Q = U0[:, vals]
        Pi = Q @ Q.conj().T
        cont = {}
        for r, ch in CH24.items():
            m = sum(ch[lab24[j]] * np.trace(Pi @ D24[j])
                    for j in range(24)) / 24.
            if abs(m.real) > 0.4:
                cont[r] = int(round(m.real))
        return cont, Pi
    cone_content, Pi_top = mults(cone_sel)
    zero_content, Pi_0 = mults(zero_sel)
    PE = 2.0 * sum(CH24["E"][lab24[j]] * D24[j] for j in range(24)) / 24.
    PT2 = 3.0 * sum(CH24["T2"][lab24[j]] * D24[j] for j in range(24)) / 24.
    quintet = float(np.trace(Pi_top @ (PE + PT2)).real)
    n_pos = int(np.sum(np.abs(bands[0] - 1.0) < 1e-9))
    n_neg = int(np.sum(np.abs(bands[0] + 1.0) < 1e-9))

    def Mstep(q):
        return expm(1j * sum(np.sin(q[k]) * H[k] for k in range(3)))
    unit = max(np.abs(Mstep(q) @ Mstep(q).conj().T - np.eye(36)).max()
               for q in [(0.3, -0.7, 1.1), (2.0, 2.0, 0.4),
                         (math.pi, math.pi, math.pi), (1.3, -0.4, 2.1)])
    at0 = float(np.abs(Mstep((0.0, 0.0, 0.0)) - np.eye(36)).max())
    N = 15
    Qs = 2.0 * np.pi * np.arange(N) / N
    Kk = np.zeros((N, N, N, 36, 36), dtype=complex)
    for a in range(N):
        for b in range(N):
            for c in range(N):
                Kk[a, b, c] = Mstep((Qs[a], Qs[b], Qs[c]))
    ker = np.fft.ifftn(Kk, axes=(0, 1, 2))
    prof = {}
    for i in range(N):
        for j in range(N):
            for l in range(N):
                r = max(min(i, N - i), min(j, N - j), min(l, N - l))
                prof[r] = max(prof.get(r, 0.0), float(np.abs(ker[i, j, l]).max()))
    t3ok = (herm < 1e-13 and cov < 1e-12 and cliff < 1e-11 and bs < 1e-12
            and abs(top - 1.0) < 1e-9 and len(cone_sel) == 8
            and quintet > 4.99 and len(zero_sel) == 4
            and unit < 1e-13 and at0 < 1e-15
            and prof[1] < 0.2 and prof[5] < 1e-4
            and max(v for r, v in prof.items() if r >= 6) < 1e-5)
    out["data"]["cone"] = {
        "clifford_defect": cliff, "O_covariance_violation": cov,
        "band_spectrum_spread_over_%d_directions" % len(dirs44): bs,
        "cone_subspace_dim": len(cone_sel),
        "cone_slope": math.sqrt(top),
        "cone_branches_at_plus_1": n_pos, "cone_branches_at_minus_1": n_neg,
        "cone_content": cone_content, "quintet_weight_on_cone": quintet,
        "non_propagating_subspace_dim": len(zero_sel),
        "non_propagating_content": zero_content,
        "step_unitarity": float(unit), "step_at_q0_minus_identity": at0,
        "kernel_profile": {str(r): prof[r] for r in sorted(prof)},
        "generator_file": "cone_generator.json"}
    checks.append(("T3 the committed triple (cone_generator.json) is Hermitian "
                   "(%.1e) and O-covariant (%.1e), is exactly Clifford-regular "
                   "({H_i,H_j} = 2 delta_ij S, S = H_1^2 = H_2^2 = H_3^2, defect "
                   "%.1e), its 36-branch spectrum is direction-independent over "
                   "%d directions to %.1e, the cone subspace has dim %d at slope "
                   "%.12f (%d branches at +1, %d at -1) with content %s - the "
                   "quintet weight on it is %.4f, i.e. E + T2 of the metric "
                   "appaired with the odd T1 connection vector - the non-"
                   "propagating subspace has dim %d with content %s (the trace "
                   "A1 does not propagate), the step exp(i sum sin(q) H) is "
                   "unitary to %.1e with M(0) = I to %.1e (no mass term exists "
                   "in this class), and the hopping kernel falls %.1e -> %.1e "
                   "-> %.1e -> %.1e -> %.1e -> %.1e -> below %.1e"
                   % (herm, cov, cliff, len(dirs44), bs, len(cone_sel),
                      math.sqrt(top), n_pos, n_neg, cone_content, quintet,
                      len(zero_sel), zero_content, unit, at0, prof[0], prof[1],
                      prof[2], prof[3], prof[4], prof[5],
                      max(v for r, v in prof.items() if r >= 6)), t3ok))

    # ---- the implemented step is not a cone generator ----------------------
    def i3(A, B):
        return sum(np.trace(A[a].conj().T @ B[a]).real for a in range(3))
    fB = (i3(H, Bop) / i3(Bop, Bop)) ** 2 * i3(Bop, Bop) / i3(H, H)
    fK = (i3(H, Kop) / i3(Kop, Kop)) ** 2 * i3(Kop, Kop) / i3(H, H)
    # the size of the search space, and of the exact-cone family through the
    # committed generator: a covariant triple is an O-intertwiner V (x) M -> M,
    # so its dimension is sum_r m_r(V (x) M) m_r(M) by characters; the family
    # dimension is the kernel of the linearised Clifford-regularity map.
    def mult24(chi):
        return {r: sum(CH24[r][lab24[j]] * chi[j] for j in range(24)) / 24.
                for r in CH24}
    mM = mult24([np.trace(Dm[i]) for i in PROP])
    mVM = mult24([np.trace(Dm[i]) * np.trace(Rv[i]) for i in PROP])
    cov_dim = int(round(sum(mVM[r] * mM[r] for r in CH24)))
    srng = np.random.default_rng(20261021)
    def rherm():
        A = (srng.standard_normal((36, 36))
             + 1j * srng.standard_normal((36, 36)))
        return (A + A.conj().T) / 2.0
    def cov_project(Ts):
        outp = [np.zeros((36, 36), dtype=complex) for _ in range(3)]
        for i in PROP:
            D, R = Dm[i], Rv[i]
            DT = [D @ Tb @ D.conj().T for Tb in Ts]
            for a in range(3):
                for b in range(3):
                    outp[a] += R[a, b] * DT[b]
        return [o / 24. for o in outp]
    def vec3(Ts):
        return np.concatenate([t.real.reshape(-1) for t in Ts]
                              + [t.imag.reshape(-1) for t in Ts])
    rows, rank = [], 0
    for _ in range(4 * cov_dim + 60):
        v = vec3(cov_project([rherm(), rherm(), rherm()]))
        r = int(np.linalg.matrix_rank(np.array(rows + [v]), tol=1e-9))
        if r > rank:
            rows.append(v)
            rank = r
        if rank == cov_dim:
            break
    _, sv, Vt = np.linalg.svd(np.array(rows), full_matrices=False)
    basis = Vt[sv > 1e-9]
    def unvec(vv):
        b = vv.reshape(6, 36, 36)
        return [b[k] + 1j * b[3 + k] for k in range(3)]
    def sym(A, B):
        return A @ B + B @ A
    def linearised(dH):
        d1, d2, d3 = dH
        return [sym(H[0], d2) + sym(d1, H[1]), sym(H[0], d3) + sym(d1, H[2]),
                sym(H[1], d3) + sym(d2, H[2]),
                sym(H[0], d1) - sym(H[1], d2),
                sym(H[0], d1) - sym(H[2], d3)]
    cols = []
    for v in basis:
        ms = linearised(unvec(v))
        cols.append(np.concatenate([m.real.reshape(-1) for m in ms]
                                   + [m.imag.reshape(-1) for m in ms]))
    Jm = np.array(cols).T
    jac_rank = int(np.linalg.matrix_rank(Jm, tol=1e-9))
    family_dim = cov_dim - jac_rank
    drs_speed = art("gauged_layer", "V141_regular_link_quintet_coupling",
                    "data.trilemma.dressing_only_speed")
    dz = art("scales_loopfactor", "V206_clock_tension", "data.dead_zone")
    half = {k: v for k, v in dz.items() if k == "1/2"}
    gw_over_c = drs_speed
    t4ok = (fB < 0.02 and fK < 0.02 and abs(gw_over_c - 2 * g_root) < 1e-7
            and "1/2" in half and half["1/2"]["n_star"] <= 3
            and cov_dim == 160 and rank == cov_dim
            and jac_rank == 107 and family_dim == 53)
    out["data"]["implemented"] = {
        "hop_fraction_of_cone_generator": float(fB),
        "dressing_fraction_of_cone_generator": float(fK),
        "implemented_cone_speed_from_V141": drs_speed,
        "GW170817_bound_abs_cgw_over_c_minus_1": 1e-15,
        "clock_giving_cone_1": 0.5,
        "V206_dead_zone_entry_for_1_2": half,
        "covariant_velocity_space_dim": cov_dim,
        "clifford_constraint_jacobian_rank": jac_rank,
        "exact_cone_family_dimension": family_dim}
    checks.append(("T4 the implemented layer is not a cone generator: the "
                   "covariant velocity space is %d-dimensional by characters, "
                   "the Clifford-regularity Jacobian at the committed generator "
                   "has rank %d in it, so the exact-cone family through it is "
                   "%d-dimensional, and neither implemented triple lies in it - "
                   "the hop accounts for %.4f of the committed generator and "
                   "the dressing for %.4f; the implemented cone speed is %.6f "
                   "(V141 dressing_only_speed = 2 g_root, g_root = pi/4) "
                   "against the matter cone 1 and |c_gw/c - 1| < 1e-15 "
                   "(GW170817); the clock that would give exactly 1, g = 1/2, "
                   "sits in V206's dead zone %s"
                   % (cov_dim, jac_rank, family_dim, fB, fK, drs_speed, half),
                   t4ok))

    # ---- Newton's constant from the bond entropy ---------------------------
    ln45 = art("gauged_layer", "V129_gauged_layer_gap", "data.entanglement.ln45")
    G_lat = 1.0 / (4.0 * ln45)
    out["data"]["newton"] = {
        "bond_entropy_ln45_from_V129": ln45,
        "G_in_cell_units_1_over_4ln45": G_lat,
        "sixteen_pi_G": 16.0 * math.pi * G_lat,
        "cell_over_planck_2_sqrt_ln45": 2.0 * math.sqrt(ln45),
        "bekenstein_deficit_factor_4ln45": 4.0 * ln45}
    t5ok = abs(ln45 - math.log(45.0)) < 1e-12 and 0.065 < G_lat < 0.067
    checks.append(("T5 the bond entropy is committed as ln 45 = %.16f (V129); "
                   "identifying one Planck bond with one Planck area gives "
                   "G = 1/(4 ln 45) = %.6f in cell units, 16 pi G = %.5f, and a "
                   "cell of 2 sqrt(ln 45) = %.5f Planck lengths; the Bekenstein "
                   "gap the paper names is the factor 4 ln 45 = %.5f"
                   % (ln45, G_lat, 16 * math.pi * G_lat, 2 * math.sqrt(ln45),
                      4 * ln45), t5ok))

    ok = all(p for _, p in checks)
    for name, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", name))
    out["checks"] = [n for n, _ in checks]
    out["pass"] = ok
    out["status"] = ("THE METRIC IS A BILINEAR OF THE BOND LAYER AND ITS CONE IS "
                     "A MODULE FACT, NOT A MODULE VERDICT. (i) Content: Sym^2 = "
                     "A1g+Eg+T2g while the six oriented bonds carry A1g+Eg+T1u - "
                     "one bond amplitude is not the metric, by exact characters; "
                     "the metric content first appears in the 36-state bilinear, "
                     "whose invariant algebra has dimension 34. (ii) Parity: the "
                     "hop's velocity is polar and vanishes inside the even "
                     "(metric) sector to 2.7e-16, so the metric's first-order "
                     "motion can only come from the axial dressing (isotropic "
                     "exactly, radius 2g) or from pairing the even and odd "
                     "halves. (iii) Such a pairing exists and is exactly "
                     "lightlike: the committed triple has Clifford defect 1.4e-13, "
                     "a direction-independent 36-branch spectrum (5e-15 over 44 "
                     "directions), the quintet E+T2 appaired with the odd "
                     "connection vector T1 on the cone at slope 1, four branches "
                     "at slope exactly 0 (the trace among them), M(0) = I exactly "
                     "(no mass term in the class), and a quasi-local kernel. "
                     "(iv) It is not the implemented generator (hop 0.3%% of it, "
                     "dressing 0.9%%; the exact-cone family is "
                     "53-dimensional), and the implemented one gives c_gw = "
                     "%.6f against the matter cone 1 - excluded by GW170817 at "
                     "1e-15 - while the clock that would fix it, g = 1/2, lies in "
                     "the closed-class dead zone. So the content of the "
                     "gravitational field is supplied by the root; its light cone "
                     "is permitted in a 53-parameter family and not selected by "
                     "anything in the current construction." % drs_speed
                     if ok else "see failed checks")
    out["wall_time"] = round(time.time() - t0, 3)
    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v208_metric_bilinear_cone.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("  VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
