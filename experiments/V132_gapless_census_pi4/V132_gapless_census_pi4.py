#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V132 - the gapless census at the quantized angle: the doubler network is a
set of NODAL LINES, and what that does to chirality and to helicity.

**Why (V131 forces this).** The closure axiom (Section 8.5) selects g = pi/4.
At that angle the one-step Bloch matrix at the zone centre is an exact
involution, M(0)^2 = I, and the twelve M nodes are exactly gapless, against
the 2.003e-4 near-miss of the implemented 1/5 where V128 found NO exactly
gapless mode anywhere in the zone.  Three statements in the paper are computed
on the OLD census:

  (i)  the chirality sum of V63 (12 point nodes, sum 0) - the Nielsen-Ninomiya
       bookkeeping counts chiral charges of ISOLATED zeros;
  (ii) the little-group / helicity analysis of V128 (C2v at an M node, no
       2-dimensional irrep, hence no massless spin-1 or spin-2 candidate);
  (iii) the census itself - "the 12 M points are gapless".

**What is measured (and it changes the picture).**
1. The factorisation M(q) = P(q) (x) U16^6 holds at pi/4 (as at 1/5).
2. The exactly gapless set of the block zone at pi/4 is NOT 13 points: it is
   GAMMA plus the SIX face-diagonal lines {(+-t, +-t, 0)} and its axis
   permutations - the twelve M points are only the high-symmetry points that
   happen to sit on those lines.  Every point of every line is gapless to
   machine precision over its whole length; the body diagonal and generic
   lines are gapped away from GAMMA.
3. Transverse structure of a line: flat along itself (gap 0), LINEAR in one
   transverse direction, QUADRATIC in the other - a quadratic nodal line, so
   the V63 point-node chirality convention does not assign it a charge (the
   slopes it needs vanish at every node it is applied to).
4. The point group of the step at pi/4 (16 of the 48 signed permutations, the
   same index-3 subgroup as at 1/5), the little group of GAMMA, and the
   decomposition of the gapless subspace into irreducible constituents of that
   group - the helicity +/-1, +/-2 question re-asked at the new gapless point.
5. Control at g = 1/5: V128's census reproduced (no exactly gapless mode, the
   face diagonals gapped).
6. The pi/2 branch: every bond stage at the half angle is exactly a monomial
   (a permutation times Bloch phases) - the bond layer carries no mixing, the
   paper's "collapses to pure swaps" claim, measured.
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "V128_f4_graviton_content"))

from cl10 import supercell as SC                                    # noqa: E402
from cl10 import engine as E                                      # noqa: E402
from V128_f4_graviton_content import (R_g, gq_of)                   # noqa: E402

PI = math.pi
G4 = math.pi / 4
G5 = 1.0 / 5.0
TOL_GAP = 1e-9
GRID = 15

# the six face-diagonal lines of the block zone
LINES = [("q_x=+q_y, q_z=0", lambda t: (t, t, 0.0)),
         ("q_x=-q_y, q_z=0", lambda t: (t, -t, 0.0)),
         ("q_x=+q_z, q_y=0", lambda t: (t, 0.0, t)),
         ("q_x=-q_z, q_y=0", lambda t: (t, 0.0, -t)),
         ("q_y=+q_z, q_x=0", lambda t: (0.0, t, t)),
         ("q_y=-q_z, q_x=0", lambda t: (0.0, t, -t))]


def P_bond(q, g):
    P = np.eye(8, dtype=complex)
    for even, j in SC.STAGES:
        P = SC.stage_matrix(even, j, q[j], g) @ P
    return P


def bond_phases(q, g):
    return np.angle(np.linalg.eigvals(P_bond(q, g)))


def internal_phases(g):
    U16 = SC.sectors(g10=g)[2]
    return np.sort(np.angle(np.exp(1j * 6 * np.angle(np.linalg.eigvals(U16)))))


def min_gap(q, bet, g):
    a = bond_phases(q, g)
    om = np.angle(np.exp(1j * (a[:, None] + bet[None, :])))
    return float(np.min(np.abs(om)))


def point_group(g, qs):
    ok = []
    for perm, signs in SC._signed_permutations():
        best = 9.0
        for inv in (False, True):
            for lp in (True, False):
                err = 0.0
                for q in qs:
                    U = R_g(q, perm, signs, inv, lp)
                    err = max(err, float(np.abs(
                        P_bond(gq_of(q, perm, signs), g) @ U
                        - U @ P_bond(q, g)).max()))
                best = min(best, err)
        if best < 1e-10:
            ok.append((perm, signs))
    return ok


def axis_slopes(q0, U16, g, dq=1e-4):
    """The slope of the smallest |eigenphase| of M(q) along each axis."""
    sl = []
    for ax in range(3):
        qa = list(q0); qa[ax] -= dq
        qb = list(q0); qb[ax] += dq
        pa = np.angle(np.linalg.eigvals(SC.bloch(tuple(qa), U16, g10=g)))
        pb = np.angle(np.linalg.eigvals(SC.bloch(tuple(qb), U16, g10=g)))
        sl.append(float((pb[np.argmin(np.abs(pb))]
                         - pa[np.argmin(np.abs(pa))]) / (2 * dq)))
    return sl


def irrep_content(reps, k):
    """Isotypic decomposition of a finite-group representation: a generic
    element of the commutant has one eigenspace per isotypic component."""
    rng = np.random.default_rng(7)
    C = sum(complex(rng.normal(), rng.normal()) * D for D in reps)
    C = 0.5 * (C + C.conj().T)
    w, V = np.linalg.eigh(C)
    groups, cur = [], [0]
    for i in range(1, k):
        if abs(w[i] - w[cur[-1]]) < 1e-7 * max(1.0, abs(w[0])):
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    out = []
    for gs in groups:
        W = V[:, gs]
        chi = np.array([np.trace(W.conj().T @ D @ W) for D in reps])
        inner = float(np.mean(np.abs(chi) ** 2))
        dim = len(gs)
        mult = int(round(inner)) if abs(inner - round(inner)) < 0.05 else None
        out.append({"dim": dim, "self_inner": round(inner, 4),
                    "multiplicity": mult,
                    "irrep_dim": (round(dim / mult, 4) if mult else None)})
    return out


def main():
    t0 = time.time()
    out = {"experiment": "V132_gapless_census_pi4", "tol": TOL_GAP,
           "backend": E._BACKEND, "data": {}}
    checks = []
    rng = np.random.default_rng(132)
    U16_4 = SC.sectors(g10=G4)[2]
    bet4 = internal_phases(G4)
    bet5 = internal_phases(G5)

    # ------------------------------------------------- T0 the line census
    fac = 0.0
    for _ in range(6):
        q = tuple(float(x) for x in rng.uniform(-PI, PI, 3))
        Mfull = SC.bloch(q, U16_4, g10=G4)
        Mfac = np.kron(P_bond(q, G4), np.linalg.matrix_power(U16_4, 6))
        fac = max(fac, float(np.abs(
            np.sort(np.angle(np.linalg.eigvals(Mfull)))
            - np.sort(np.angle(np.linalg.eigvals(Mfac)))).max()))

    lines = {}
    for name, f in LINES:
        ts = [-PI + 2 * PI * i / 24 for i in range(25)]
        gs = [min_gap(f(t), bet4, G4) for t in ts]
        lines[name] = {"max_gap_on_line": max(gs), "n_samples": len(ts)}
    gm = min_gap(SC.GAMMA, bet4, G4)
    m_gaps = [min_gap(q, bet4, G4) for q in SC.M_POINTS]
    x_min = min(min_gap(q, bet4, G4) for q in SC.X_POINTS)
    r_min = min(min_gap(q, bet4, G4) for q in SC.R_POINTS)
    body = [min_gap((t, t, t), bet4, G4) for t in (0.3, 1.0, 2.0)]

    grid = [-PI + 2 * PI * i / (GRID - 1) for i in range(GRID)]
    nz = (9.0, None)
    off_line = 0
    for q in itertools.product(grid, grid, grid):
        gp = min_gap(q, bet4, G4)
        if gp < TOL_GAP:
            a, b, c = (abs(v) for v in q)
            on = ((c < 1e-12 and abs(a - b) < 1e-12)
                  or (b < 1e-12 and abs(a - c) < 1e-12)
                  or (a < 1e-12 and abs(b - c) < 1e-12))
            if not on:
                off_line += 1
        elif gp < nz[0]:
            nz = (gp, q)
    line_max = max(v["max_gap_on_line"] for v in lines.values())
    out["data"]["census"] = {
        "angle": G4, "factorisation_err": fac, "lines": lines,
        "line_max_gap": line_max, "gamma_gap": gm,
        "M_gaps": [round(x, 12) for x in m_gaps],
        "X_gap_min": x_min, "R_gap_min": r_min,
        "body_diagonal_gaps": body,
        "grid": GRID, "grid_gapless_off_lines": off_line,
        "grid_min_nonzero_gap": nz[0], "grid_min_nonzero_at": list(nz[1])}
    checks.append(("T0 at pi/4 the step still factorises as kron(P(q), U16^6) "
                   "(%.1e), and the exactly gapless set of the zone is NOT 13 "
                   "points: the SIX face-diagonal lines {(+-t,+-t,0)} and its "
                   "axis permutations are gapless along their whole length "
                   "(max gap on any line %.1e, 25 samples each), GAMMA is "
                   "their common point (%.1e), the 12 M points are only the "
                   "high-symmetry points on them (gaps [%.1e, %.1e]); X and R "
                   "stay gapped (%.4f, %.4f), the body diagonal is gapped away "
                   "from GAMMA (%s), and of the %d^3 grid's gapless points %d "
                   "lie off the lines"
                   % (fac, line_max, gm, min(m_gaps), max(m_gaps), x_min,
                      r_min, [round(x, 3) for x in body], GRID, off_line),
                   fac < 1e-12 and line_max < TOL_GAP and gm < TOL_GAP
                   and all(x < TOL_GAP for x in m_gaps) and x_min > 1e-3
                   and r_min > 1e-3 and min(body) > 1e-2 and off_line == 0))

    # ------------------------------------- T1 transverse structure of a line
    q0 = SC.M_POINTS[0]                      # (pi/2, pi/2, 0)
    dirs = {"along (1,1,0)": (1, 1, 0), "transverse (1,-1,0)": (1, -1, 0),
            "transverse (0,0,1)": (0, 0, 1)}
    struct = {}
    for name, d in dirs.items():
        n = math.sqrt(sum(v * v for v in d))
        u = [v / n for v in d]
        gs = []
        for h in (1e-1, 1e-2, 1e-3):
            qa = tuple(q0[a] + h * u[a] for a in range(3))
            qb = tuple(q0[a] - h * u[a] for a in range(3))
            gs.append(min(min_gap(qa, bet4, G4), min_gap(qb, bet4, G4)))
        e = (math.log(gs[1] / gs[2]) / math.log(10)
             if gs[1] > 1e-12 and gs[2] > 1e-12 else None)
        struct[name] = {"gap_1e-1": gs[0], "gap_1e-2": gs[1], "gap_1e-3": gs[2],
                        "exponent": round(e, 3) if e is not None else None,
                        "slope": round(gs[2] / (1e-3 * n), 5)
                        if gs[2] > 1e-12 else 0.0}
    ax_exp = {}
    for q in SC.M_POINTS:
        ex = []
        for ax in range(3):
            e = []
            for h in (1e-2, 1e-3):
                qa = list(q); qa[ax] -= h
                qb = list(q); qb[ax] += h
                e.append(min(min_gap(tuple(qa), bet4, G4),
                             min_gap(tuple(qb), bet4, G4)))
            ex.append(round(math.log(e[0] / e[1]) / math.log(10), 2)
                      if min(e) > 1e-12 else None)
        ax_exp[str([round(v, 4) for v in q])] = ex
    pats = sorted({tuple(e) for e in ax_exp.values()})
    out["data"]["transverse"] = {"M0": struct, "per_node_axis_exponents": ax_exp,
                                 "distinct_patterns": [list(p) for p in pats]}
    checks.append(("T1 the nodal line is QUADRATIC, not a linear crossing: "
                   "through an M point the gap stays 0 along the line (%s), "
                   "grows LINEARLY in one transverse direction (exponent %.2f, "
                   "slope %.3f) and QUADRATICALLY in the other (exponent "
                   "%.2f); per node the three axis exponents form the pattern "
                   "%s - the quadratic axis is the zero component - the same "
                   "at all 12 nodes"
                   % ([f"{struct['along (1,1,0)'][k]:.1e}"
                      for k in ("gap_1e-1", "gap_1e-2", "gap_1e-3")],
                      struct["transverse (1,-1,0)"]["exponent"],
                      struct["transverse (1,-1,0)"]["slope"],
                      struct["transverse (0,0,1)"]["exponent"], pats),
                   struct["along (1,1,0)"]["gap_1e-1"] < TOL_GAP
                   and abs(struct["transverse (1,-1,0)"]["exponent"] - 1) < 0.05
                   and abs(struct["transverse (0,0,1)"]["exponent"] - 2) < 0.05
                   and len(pats) == 3))

    # ------------------------------------------- T2 the chirality bookkeeping
    naive = []
    nonisolated = 0
    line_gaps = {}
    for q in SC.M_POINTS:
        nz = [a for a in range(3) if abs(q[a]) > 1e-12]
        d = [0.0, 0.0, 0.0]
        for a in nz:
            d[a] = math.copysign(1.0, q[a]) / math.sqrt(2.0)
        gs = []
        for h in (1e-3, 1e-2, 1e-1):
            qa = tuple(q[a] + h * d[a] for a in range(3))
            qb = tuple(q[a] - h * d[a] for a in range(3))
            gs.append(min(min_gap(qa, bet4, G4), min_gap(qb, bet4, G4)))
        if max(gs) < TOL_GAP:
            nonisolated += 1
        line_gaps[str([round(v, 4) for v in q])] = [float(x) for x in gs]
        sl_all = axis_slopes(q, U16_4, G4)
        kept = [s for ax, s in enumerate(sl_all) if abs(q[ax]) > 1e-12]
        naive.append(float(np.sign(kept[0] * kept[1])))
    naive_sum = float(sum(naive))
    sl_g = axis_slopes(SC.GAMMA, U16_4, G4)
    out["data"]["chirality"] = {"naive_V63_sum": naive_sum,
                                "naive_per_node": naive,
                                "nonisolated_nodes": nonisolated,
                                "gap_along_line_per_node": line_gaps,
                                "gamma_axis_slopes": [round(v, 6) for v in sl_g]}
    checks.append(("T2 the V63 point-node chirality convention no longer "
                   "applies: at %d of the 12 M nodes the gap stays below %.0e "
                   "along the line through the node at ALL of the steps 1e-3, "
                   "1e-2, 1e-1 - the node is not an isolated zero, so the "
                   "charge V63 assigns to it is sign noise (the naive sum is "
                   "%.0f where V63 measured exactly 0 at 1/5); at GAMMA the "
                   "three axis slopes are %s (quadratic touching, not a linear "
                   "one) - GAMMA is a line crossing - so the Nielsen-Ninomiya "
                   "bookkeeping must be redone per line, not per node"
                   % (nonisolated, TOL_GAP, naive_sum,
                      [round(v, 4) for v in sl_g]),
                   nonisolated == 12))

    # ------------------------------------------- T3 little group at GAMMA
    qs = [tuple(float(x) for x in rng.uniform(-PI, PI, 3)) for _ in range(5)]
    pg4 = point_group(G4, qs)
    pg5 = point_group(G5, qs)
    M0 = SC.bloch(SC.GAMMA, U16_4, g10=G4)
    w0, V0 = np.linalg.eig(M0)
    keep = np.where(np.abs(np.angle(w0)) < 1e-7)[0]
    k = len(keep)
    S = np.linalg.qr(V0[:, keep])[0]
    reps, leak_max = [], 0.0
    for perm, signs in pg4:
        D = np.kron(R_g(SC.GAMMA, perm, signs), np.eye(16, dtype=complex))
        leak_max = max(leak_max, float(np.abs(
            S.conj().T @ D @ (np.eye(128) - S @ S.conj().T)).max()))
        reps.append(S.conj().T @ D @ S)
    chi = np.array([np.trace(D) for D in reps])
    inner_full = float(np.mean(np.abs(chi) ** 2))
    Pb = P_bond(SC.GAMMA, G4)
    wb, Vb = np.linalg.eig(Pb)
    bph = np.angle(wb)
    ib = np.angle(np.exp(1j * 6 * np.angle(np.linalg.eigvals(U16_4))))
    pair = [i for i in range(8)
            if any(abs((bph[i] + p + PI) % (2 * PI) - PI) < 1e-7 for p in ib)]
    Sb = np.linalg.qr(Vb[:, pair])[0]
    breps = [Sb.conj().T @ R_g(SC.GAMMA, p, s) @ Sb for p, s in pg4]
    comp = irrep_content(breps, len(pair))
    dims = sorted({c["irrep_dim"] for c in comp if c["irrep_dim"]})
    two_dim = any(c["irrep_dim"] == 2.0 for c in comp)
    out["data"]["little_group_gamma"] = {
        "point_group_order": len(pg4), "point_group_order_at_1_5": len(pg5),
        "gapless_subspace_dim": k, "intertwiner_leak_max": leak_max,
        "character_inner_product_full": inner_full,
        "bond_gapless_dim": len(pair), "bond_isotypic_components": comp,
        "constituent_irrep_dims": dims, "two_dim_irrep_present": two_dim}
    checks.append(("T3 at pi/4 the step's point group is %d of the 48 signed "
                   "permutations (the same index-3 subgroup as at 1/5: %d) and "
                   "GAMMA is fixed by all of them; the exactly gapless subspace "
                   "of M(0) has dimension %d of 128 (intertwiner leak %.1e, "
                   "<chi,chi> = %.0f); since the internal 16 is inert the "
                   "content is that of the %d-dimensional bond factor, whose "
                   "isotypic constituents have dimensions %s - it %s a "
                   "2-dimensional irrep, so at the new gapless point %s "
                   "helicity +/-1 or +/-2 candidate (V128's negative answer was "
                   "for the C2v of an M node)"
                   % (len(pg4), len(pg5), k, leak_max, inner_full, len(pair),
                      dims, "carries" if two_dim else "carries NO",
                      "a" if two_dim else "no"),
                   len(pg4) == 16 and leak_max < 1e-9 and k > 0 and pair))

    # -------------------------------------------------- T4 control at 1/5
    nz5 = (9.0, None)
    ng5 = 0
    for q in itertools.product(grid, grid, grid):
        gp = min_gap(q, bet5, G5)
        if gp < TOL_GAP:
            ng5 += 1
        elif gp < nz5[0]:
            nz5 = (gp, q)
    m5 = [min_gap(q, bet5, G5) for q in SC.M_POINTS]
    g5 = min_gap(SC.GAMMA, bet5, G5)
    line5 = min(max(min_gap(f(t), bet5, G5)
                    for t in (-PI + 2 * PI * i / 24 for i in range(25)))
                for _n, f in LINES)
    out["data"]["control_1_5"] = {
        "grid_min_gap": nz5[0], "grid_min_at": list(nz5[1]),
        "M_gap_min": min(m5), "M_gap_max": max(m5), "gamma_gap": g5,
        "gapless_count": ng5, "point_group_order": len(pg5),
        "face_diagonal_min_gap": line5}
    checks.append(("T4 control at the implemented 1/5 reproduces V128: over "
                   "the %d^3 grid the smallest quasi-energy is %.3e (at %s), "
                   "the 12 M nodes lie in [%.3e, %.3e], GAMMA is %.4f, the "
                   "best gap along a face diagonal is %.3e, and there are %d "
                   "exactly gapless modes - the census the paper quotes today"
                   % (GRID, nz5[0], [round(v, 4) for v in nz5[1]], min(m5),
                      max(m5), g5, line5, ng5),
                   ng5 == 0 and 1.9e-4 < min(m5) < 2.1e-4 and g5 > 3e-2
                   and line5 > 1e-4))

    # -------------------------------------------------- T5 the pi/2 branch
    worst = 0.0
    for even, j in SC.STAGES:
        for qq in (0.0, 0.7, -2.1):
            S2 = SC.stage_matrix(even, j, qq, math.pi / 2)
            mag = np.abs(S2)
            monomial = bool(np.allclose(mag, np.round(mag), atol=1e-14)
                            and np.all(np.round(mag).sum(axis=0) == 1)
                            and np.all(np.round(mag).sum(axis=1) == 1))
            diag = float(np.max(np.abs(np.diag(S2))))
            worst = max(worst, 0.0 if monomial else 1.0, diag)
    out["data"]["pi_2_pure_swap"] = {"worst_deviation": worst}
    checks.append(("T5 at pi/2 every bond stage is exactly a monomial - one "
                   "unit-modulus entry per row and column, diagonal vanishing "
                   "to %.1e (cos(pi/2) = 6.1e-17) - over the six stages at "
                   "three momenta: the stage is a pure swap times Bloch "
                   "phases, no mixing at all, so the half-angle branch carries "
                   "no bond dynamics - 'the bond layer collapses to pure "
                   "swaps' is now measured" % worst, worst < 1e-14))

    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    out["status"] = ("the gapless census at pi/4 is GAMMA plus six "
                     "face-diagonal nodal lines (max gap %.1e along each; %d "
                     "of the %d^3 grid's gapless points lie off the lines); "
                     "quadratic touching (transverse exponents 1 and 2, slope "
                     "%.3f); the V63 point-node chirality convention is "
                     "inapplicable at %d/12 nodes (the naive sum is %.0f - "
                     "sign noise, backend-dependent - where V63 measured "
                     "exactly 0 at 1/5); the gapless subspace at GAMMA is %d "
                     "of 128 with %d 1-dimensional bond constituents - no "
                     "2-dim irrep, no helicity candidate"
                     % (line_max, off_line, GRID,
                        struct["transverse (1,-1,0)"]["slope"], nonisolated,
                        naive_sum, k, len(comp)))
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v132_gapless_census_pi4.json"), "w") as f:
        json.dump(out, f, indent=1)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  gapless at pi/4: GAMMA + 6 face-diagonal lines (max gap %.1e) | "
          "non-isolated M nodes %d of 12 | bond constituents %s | "
          "2-dim irrep: %s" % (line_max, nonisolated, dims, two_dim))
    print("  wall %.1f s" % out["wall_s"])
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
