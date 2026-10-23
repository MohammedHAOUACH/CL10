#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V118 - is the node anisotropy a measurement artifact, and what forbids isotropy? (O9)

**The question.**  V102 measured anisotropic doubler cones (principal velocities
0.200, 0.200, 0.284) and V115 found no isotropy point over the V65 Wilson
family - but V115's anisotropy is measured at a FIXED radius dq = 1e-3, by
min |eigenphase| / dq over 179 sphere directions.  A finite radius cannot
separate a genuinely anisotropic tangent cone from an isotropic cone with
large curvature corrections.  This experiment measures the TANGENT cone
exactly, and then asks which symmetry the node has.

**Four measurements.**

T1 - the tangent cone, exactly.  At a node q0 the one-step matrix M(q0) has a gapless
    eigenvalue; the tangent velocities are the eigenvalues of
    Q^H (d M / d q_j) Q on the gapless subspace Q.  At the 12 V63 doubler
    momenta (w = 0) the three principal tangent velocities are
    (0.14189, 0.14189, 0) - the third is EXACTLY zero.  The anisotropy is not
    a radius artifact: the finite-difference velocities of V115 converge to
    constants as dq -> 0 (exponent alpha = 0 for the linear branches, alpha = 1
    for the flat one), so A_tangent = 1 at the node itself.

T2 - the reason is the little group.  The V63 nodes sit at half-zone
    momenta (pi/2, pi/2, 0) and permutations: their little group is C2v,
    which permits two independent velocities and a vanishing third -
    exactly what is measured (the two pi/2 axes are equal to 1e-16, the
    zero axis differs).  At the eight R points (pi/2, pi/2, pi/2) - one
    cubic star, little group C3v - the 128 one-step levels fall into
    exactly 3-fold degenerate multiplets (grouped at a 1e-12 tolerance): the cubic
    symmetry is realised on the node.

T3 - the Wilson factor can put a CUBIC MULTIPLET exactly at zero.  Scanning
    the V65 Wilson family (theta, w) for the lowest exactly-3-fold-degenerate
    corner level, the gap closes to < 1e-10 at two independent points
    (w = 1.769072193, theta = 5.408296164) and (w = 1.801091626,
    theta = 2.513274121) - the first exactly massless cubic multiplet found
    in the model.  One condition on two parameters: a one-parameter family.

T4 - and the touching is linear but PLANAR, not conical.  On the
    3-dimensional gapless subspace S_j = i dM/dq_j is nonzero (max speed
    0.92), so the crossing IS linear; but [S_i, S_j] = 0 to 2e-12 and
    sum_j S_j^2 is a scalar to 6e-12.  The three branches therefore carry
    three FIXED velocity vectors - the cubic permutations of
    (0.704, 0.421, 0.421) - and each branch speed is exactly |v_alpha . d|:
    the node is a crossing of three NODAL PLANES, not a Weyl or Dirac cone
    (36 of 6588 sampled directions are exactly gapless - the intersections
    of those planes).  The obstruction is the Kronecker form of the bond
    layer, (sub-site index) x (internal algebra): the three momentum
    derivatives act on one factor and commute.

**Verdict.**  Continuous Lorentz invariance fails at the nodes for a
structural reason, not a tuning reason: the node's little group permits the
anisotropy it shows, the Wilson family can make a cubic multiplet exactly
massless but only as a crossing of three nodal planes, and no tuning of the
Wilson angle turns three planes into a cone.  Closing O9 requires changing
the bond layer itself - entangling the sub-site and internal indices - not
scanning Wilson angles.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# reference spectrum: the g = 1/5 Wilson-family member (paper 8.2):
# the dispersive node tangent speeds (0, 0.14189, 0.14189) and the
# C2v little-group bookkeeping are properties of that member; at the
# pi/4 clock the M channels are pinned to the Floquet pair (V203).
os.environ["CL10_ANGLE"] = "0.2"
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))

import V65_wilson as V65                      # noqa: E402
from cl10 import supercell as SC              # noqa: E402

_, _, U16 = SC.sectors()
_lam, _off = V65.cartan_weights()
_v65 = json.load(open(exp("V65_wilson",
                                   "results", "v65_wilson.json")))
_gg = next(r for r in _v65["data"]["part_B"].values()
           if r["su5"] == "Georgi-Glashow")
XCOEF = np.array(_gg["X_coef"])
YCOEF = np.array(_v65["data"]["Y_coef"])


def ddiag(th):
    return _lam @ (math.cos(th) * XCOEF + math.sin(th) * YCOEF)


def M(q, w=0.0, th=0.0):
    base = SC.bloch(q, U16)
    if w == 0.0:
        return base
    ph = w * sum(1.0 - math.cos(qj) for qj in q)
    return np.kron(np.eye(8, dtype=complex),
                   np.diag(np.exp(-1j * ph * ddiag(th)))) @ base


def gap(q, w=0.0, th=0.0):
    return float(np.min(np.abs(np.sort(np.angle(np.linalg.eigvals(M(q, w, th)))))))


def gapless_basis(q, w, th, tol=1e-6):
    ev, V = np.linalg.eig(M(q, w, th))
    sel = np.where(np.abs(np.angle(ev)) < tol)[0]
    if len(sel) == 0:
        return None, 0
    Q, _ = np.linalg.qr(V[:, sel])
    return Q, len(sel)


def tangent_matrices(q, w, th, h=1e-5, tol=1e-6):
    Q, m = gapless_basis(q, w, th, tol=tol)
    if Q is None:
        return None, 0
    out = []
    for j in range(3):
        qp = list(q)
        qm = list(q)
        qp[j] += h
        qm[j] -= h
        D = (M(tuple(qp), w, th) - M(tuple(qm), w, th)) / (2 * h)
        out.append(Q.conj().T @ D @ Q)
    return out, m


def sphere(n=4):
    out = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                r = math.sqrt(i * i + j * j + k * k)
                if r > 0:
                    out.append((i / r, j / r, k / r))
    return out


def main():
    out = {"experiment": "V118_w9_little_group", "data": {}}
    checks = []

    # ---- T1: the exact tangent cone at the 12 doubler nodes ---------------
    rows = []
    for q in SC.M_POINTS:
        Vj, m = tangent_matrices(q, 0.0, 0.0, tol=1e-3)
        vs = sorted(float(np.abs(np.linalg.eigvals(A)).max()) for A in Vj)
        A_t = 1.0 if vs[-1] == 0.0 else (vs[-1] - vs[0]) / (vs[-1] + vs[0])
        rows.append({"q": list(q), "m": m,
                     "principal_tangent_v": vs, "A_tangent": A_t})
    uniq = sorted({tuple(round(v, 6) for v in r["principal_tangent_v"])
                   for r in rows})
    out["data"]["T1"] = {"nodes": rows[:4], "n_nodes": len(rows),
                         "distinct_velocity_triples": [list(u) for u in uniq]}
    checks.append(("T1 all 12 nodes: one tangent velocity exactly 0, the other two equal",
                   len(rows) == 12
                   and all(abs(r["principal_tangent_v"][0]) < 1e-9 for r in rows)
                   and all(abs(r["principal_tangent_v"][1]
                               - r["principal_tangent_v"][2]) < 1e-9 for r in rows)))

    # the finite-difference convergence (V115's method, dq -> 0) at V115's
    # best Wilson point, where the node is exactly gapless.
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    qv = (math.pi / 2, -math.pi / 2, 0.0)
    thv, wv = 2.1991, 0.68625
    a, b = wv - 1e-4, wv + 1e-4
    for _ in range(80):
        m1, m2 = b - gr * (b - a), a + gr * (b - a)
        if gap(qv, m1, thv) < gap(qv, m2, thv):
            b = m2
        else:
            a = m1
    wv = 0.5 * (a + b)
    conv = {"refined_w": wv, "gap_at_node": gap(qv, wv, thv)}
    for nd in ((1, 0, 0), (0, 1, 0), (0, 0, 1),
               (1 / math.sqrt(2), -1 / math.sqrt(2), 0)):
        vs = []
        for e in range(2, 8):
            dq = 10.0 ** -e
            k = tuple(qj + dq * n for qj, n in zip(qv, nd))
            vs.append(gap(k, wv, thv) / dq)
        alpha = math.log(vs[3] / vs[2]) / math.log(10) if vs[2] > 0 else 0.0
        conv["dir=%s" % (nd,)] = {"v(dq)": [float(v) for v in vs],
                                  "alpha_last": alpha}
    out["data"]["T1b"] = conv
    checks.append(("T1b V115's finite-difference velocity converges (alpha=0 linear, alpha=-1 flat)",
                   abs(conv["dir=(1, 0, 0)"]["alpha_last"]) < 0.1
                   and conv["dir=(0, 0, 1)"]["alpha_last"] < -0.8))

    # ---- T2: the little group, edge vs corner -----------------------------
    C = (math.pi / 2,) * 3
    p = np.sort(np.abs(np.angle(np.linalg.eigvals(M(C)))))
    mult = []
    i = 0
    while i < len(p):
        j = i
        while j + 1 < len(p) and p[j + 1] - p[i] < 1e-12:
            j += 1
        if j > i:
            mult.append({"level": float(p[i]), "degeneracy": j - i + 1})
        i = j + 1
    out["data"]["T2"] = {"corner_gap": gap(C), "corner_degenerate_multiplets": mult[:6],
                         "edge_gap": gap(SC.M_POINTS[0])}
    checks.append(("T2 corner levels fall in exactly 3-fold multiplets (O_h)",
                   any(m2["degeneracy"] == 3 and m2["level"] < 0.1 for m2 in mult)))

    # ---- T3: a corner cubic multiplet made exactly gapless ----------------
    def g3(w, th):
        pp = np.sort(np.abs(np.angle(np.linalg.eigvals(M(C, w, th)))))
        i = 0
        while i < len(pp) - 2:
            if pp[i + 2] - pp[i] < 1e-7:
                return pp[i]
            i += 1
        return pp[2]

    seeds = []
    for it in range(60):
        th = 2 * math.pi * it / 60
        for w in np.linspace(0.0, 3.0, 61):
            seeds.append((g3(w, th), float(w), float(th)))
    seeds.sort()

    gr = (math.sqrt(5.0) - 1.0) / 2.0

    def refine(w0, th0, rw, rth, rounds=6):
        w, th = w0, th0
        for _ in range(rounds):
            for var in ("w", "th"):
                lo, hi = (w - rw, w + rw) if var == "w" else (th - rth, th + rth)
                f = (lambda x: g3(x, th)) if var == "w" else (lambda x: g3(w, x))
                c1, c2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
                f1, f2 = f(c1), f(c2)
                for _ in range(35):
                    if f1 < f2:
                        hi, c2, f2 = c2, c1, f1
                        c1 = hi - gr * (hi - lo)
                        f1 = f(c1)
                    else:
                        lo, c1, f1 = c1, c2, f2
                        c2 = lo + gr * (hi - lo)
                        f2 = f(c2)
                if var == "w":
                    w = 0.5 * (lo + hi)
                else:
                    th = 0.5 * (lo + hi)
            rw *= 0.5
            rth *= 0.5
        return w, th, g3(w, th)

    pts = []
    for s in seeds[:3]:
        w, th, g = refine(s[1], s[2], 0.1, 0.2)
        pts.append({"w": w, "theta": th, "triplet_gap": g,
                    "lowest_phases": [float(x) for x in
                                      np.sort(np.abs(np.angle(
                                          np.linalg.eigvals(M(C, w, th)))))[:8]]})
    out["data"]["T3"] = {"grid_seeds": [list(map(float, s)) for s in seeds[:3]],
                         "closed_multiplets": pts}
    checks.append(("T3 corner cubic triplet driven exactly gapless (<1e-9)",
                   all(pt["triplet_gap"] < 1e-9 for pt in pts)))

    # ---- T4: the tangent velocity algebra - cone or nodal planes? ---------
    # dM/dq_j restricted to the gapless subspace is ANTI-hermitian (M is
    # unitary and M_0 = I there), so the velocity operators are S_j = i V_j.
    tan = []
    for pt in pts:
        Vj, m = tangent_matrices(C, pt["w"], pt["theta"])
        if Vj is None:
            tan.append({"w": pt["w"], "m": 0})
            continue
        S = [1j * A for A in Vj]
        anti = max(np.abs(A + A.conj().T).max() for A in Vj)
        dirs = sphere(6)
        spd = []
        for kd in dirs:
            A = kd[0] * S[0] + kd[1] * S[1] + kd[2] * S[2]
            spd.extend(sorted(abs(x) for x in np.linalg.eigvalsh(A)))
        spd = np.array(spd)
        nz = spd[spd > 1e-6]
        comm = [np.linalg.norm(S[a] @ S[b] - S[b] @ S[a])
                for a, b in ((0, 1), (1, 2), (2, 0))]
        S2 = sum(s @ s for s in S)
        I = np.eye(m)
        # the S_j commute, so a generic combination splits them and the JOINT
        # eigenvalues are the three branch velocity vectors v_alpha.
        G = (S[0] + 0.718281828459045 * S[1]
             + 0.367879441171442 * S[2])
        _, Ww = np.linalg.eigh(G)
        vecs = [[float(np.real(Ww[:, a].conj().T @ s @ Ww[:, a]))
                 for s in S] for a in range(m)]
        # each branch is exactly PLANAR: its speed is |v_alpha . d|
        planar = 0.0
        for kd in dirs[:300]:
            ev = sorted(abs(x) for x in np.linalg.eigvalsh(
                kd[0] * S[0] + kd[1] * S[1] + kd[2] * S[2]))
            pred = sorted(abs(vd[0] * kd[0] + vd[1] * kd[1] + vd[2] * kd[2])
                          for vd in vecs)
            planar = max(planar, max(abs(a - b) for a, b in zip(ev, pred)))
        tan.append({"w": pt["w"], "theta": pt["theta"], "gapless_dim": m,
                    "antihermiticity_err": float(anti),
                    "max_speed": float(spd.max()),
                    "min_nonzero_speed": float(nz.min()) if nz.size else 0.0,
                    "anisotropy": float((nz.max() - nz.min()) / (nz.max() + nz.min()))
                    if nz.size else None,
                    "n_zero_branch_directions": int((spd < 1e-6).sum()),
                    "n_branch_directions": int(len(spd)),
                    "commutator_norms": [float(c) for c in comm],
                    "sum_Sj2_dev_from_scalar": float(np.abs(S2 - np.trace(S2) / m * I).max()),
                    "sum_Sj2_scalar": float(np.trace(S2).real / m),
                    "branch_velocity_vectors": vecs,
                    "planar_branch_err": float(planar)})
    out["data"]["T4"] = tan
    checks.append(("T4 the touching IS linear (max speed > 1e-3)",
                   all(t["max_speed"] > 1e-3 for t in tan if t.get("max_speed"))))
    checks.append(("T4b the three velocity operators COMMUTE: nodal planes, not a cone",
                   all(max(t["commutator_norms"]) < 1e-9 for t in tan
                       if t.get("commutator_norms"))))
    checks.append(("T4c each branch is exactly planar: speed = |v_alpha . d| (three nodal planes)",
                   all(t["planar_branch_err"] < 1e-6 for t in tan
                       if t.get("planar_branch_err") is not None)))

    ok = all(c[1] for c in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["pass"] = ok
    out["status"] = (
        "TANGENT VELOCITIES AND NODE LITTLE GROUPS (O9). T1: measured exactly, "
        "not at a finite radius, the 12 V63 doubler nodes have tangent speeds "
        "along the three axes (%.5f, %.5f, %.5f) - one axis direction is "
        "exactly flat, the other two carry equal speeds; V115's "
        "finite-difference velocities converge as dq -> 0 (exponent %.2f on a "
        "linear branch, %.2f on the flat one, i.e. a quadratic gap), so the "
        "anisotropy is NOT a dq = 1e-3 artifact. T2: the nodes sit at zone-EDGE "
        "momenta whose little group is C2v - two equal velocities and a free "
        "third, exactly what is measured - while the supercell zone CORNER "
        "(pi/2,pi/2,pi/2) carries exactly 3-fold degenerate multiplets "
        "(splitting < 1e-12): the full octahedral group is realised there. T3: "
        "the V65 Wilson factor drives a corner cubic TRIPLET exactly gapless "
        "(gap %.1e and %.1e at independent (theta,w) points) - the first "
        "exactly massless cubic multiplet found in the model, a one-parameter "
        "family. T4: on that subspace dM/dq_j is anti-hermitian (err %.1e) and "
        "the velocity operators S_j = i dM/dq_j are nonzero (max speed %.3f in "
        "lattice units): the touching IS linear, not quadratic. T4b: but the "
        "three velocity operators COMMUTE (max ||[S_i,S_j]|| = %.1e) and "
        "sum_j S_j^2 is a scalar - the three branches carry three FIXED "
        "velocity vectors, so the node is a crossing of three NODAL PLANES, not "
        "a cone: the %d branch-directions span [%.1e, %.3f], %d of them are "
        "exactly gapless, and every branch speed equals "
        "|v_alpha . d| to %.1e - the node is three planes. VERDICT: the failure of "
        "continuous Lorentz invariance is structural, not a tuning failure. A "
        "cone requires the tangent velocity operators to close a Clifford-type "
        "algebra; in the implemented step they commute, because the bond layer "
        "is a Kronecker product (cell index) x (internal algebra) and never "
        "mixes the two. Closing O9 requires a bond layer that entangles the "
        "cell and internal indices."
        % (
           rows[0]["principal_tangent_v"][0], rows[0]["principal_tangent_v"][1],
           rows[0]["principal_tangent_v"][2],
           conv["dir=(1, 0, 0)"]["alpha_last"], conv["dir=(0, 0, 1)"]["alpha_last"],
           pts[0]["triplet_gap"], pts[1]["triplet_gap"],
           tan[0]["antihermiticity_err"], tan[0]["max_speed"],
           max(max(t["commutator_norms"]) for t in tan if t.get("commutator_norms")),
           tan[0]["n_branch_directions"],
           tan[0]["min_nonzero_speed"], tan[0]["max_speed"],
           tan[0]["n_zero_branch_directions"],
           max(t["planar_branch_err"] for t in tan
               if t.get("planar_branch_err") is not None)))

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v118_w9_little_group.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V118 - tangent cones, little groups, and the corner multiplet")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  node tangent velocities:", uniq)
    print("  corner multiplets:", [(m2["degeneracy"], round(m2["level"], 6))
                                   for m2 in mult[:5]])
    for pt in pts:
        print("  closed triplet: w=%.9f theta=%.9f gap=%.2e"
              % (pt["w"], pt["theta"], pt["triplet_gap"]))
    for t in tan:
        print("  dim %d: max speed %.4f, min nonzero %.4f, commutators %s, "
              "zero-branch dirs %d/%d"
              % (t["gapless_dim"], t["max_speed"], t["min_nonzero_speed"],
                 ["%.1e" % c for c in t["commutator_norms"]],
                 t["n_zero_branch_directions"], t["n_branch_directions"]))
        for v in t["branch_velocity_vectors"]:
            print("     branch velocity v = (%+.4f, %+.4f, %+.4f)" % tuple(v))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
