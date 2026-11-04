#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V330 - is the family number 3 forced by a size-3 orbit anywhere?
(tier 3, O3, theory-hunt closure programme.)

**The question (pre-registered).** V91 decided the count 3 is a
symmetry-breaking input on the 12 nodes. Does ANY booked group-action
object carry an orbit of size exactly 3, in every row of the class -
is 3 forced anywhere in the booked machinery?

**Closed class (booked pre-run, verbatim).** Groups: the 48 signed
permutations V114 group_elements() and the 24-element chiral subgroup
V114 chiral_elements() (imported from V114_junction_r by path); the
cubic node action is V91's construction verbatim (its cub signed-
permutation matrices, its fold(), its image(); V91 nests them inside
main(), so the short bodies are copied verbatim - hunt precedent).
Objects and actions:
 (a) the 12 doubler nodes (M_POINTS): orbits under the 48 group and
     under the chiral 24, by V91's fold/image node action;
 (b) the 6 chiral stabilizer elements of the V76 modulation on (Z/4)^3
     (V114's profile-stabilizer body verbatim): orbits under conjugation
     by the 48 group (elements as 3x3 signed permutation matrices,
     conjugation P_h -> P_g P_h P_g^T);
 (c) the 16-level single-generation basis: the 16x16 block B(g) of
     cell_step_np(g10 = g) on the Weyl-16 (the 16-subspace of the 32-D
     Weyl module from clifford.weyl_spinor; the step preserves the
     16/16-bar split - V66), at each of the five scan_clocks; the
     monomiality of B(g) tested and booked at every clock; a permutation
     row exists only where B(g) is monomial, with the orbit sizes of
     the induced permutation. Closed-clock anchor booked alongside:
     angle.branches.closed = pi/4 and the V239 grid_labels_by_order
     {'2': 1, '4': 0, '8': 0, '16': 0}.
"Forced 3": an orbit of size exactly 3 in EVERY census row (nodes under
both groups, the stabilizer-6 under conjugation, and every existing
16-level permutation row). The full orbit-size histogram booked always.

**T0 (pre-registered).** V91 data.orbit_sizes == [12] and
data.triangle_orbit_size == 8 (artifact leaves; the construction
re-derived verbatim and compared); V114 data.stabilizer_size == 12 and
data.stabilizer_chiral == 6 (re-derived). Failure: MACHINERY-DRIFT.

**Branch (pre-registered).** GEN3-ORBIT-FORCED (the orbit booked with
its object, group and clock) / GEN3-NOT-FORCED-DECLARED (the histogram
booked; no size-3 orbit in every row).
"""
import itertools
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
sys.path.insert(0, exp("V114_junction_r"))

import V114_junction_r as V114                            # noqa: E402
from cl10 import cell as CE                               # noqa: E402
from cl10 import clifford as CL                           # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import supercell as SC                          # noqa: E402

CFG = _config.CFG
TWO_PI = 2.0 * math.pi
SCAN = [float(g) for g in CFG.pillar_surgery.scan_clocks]

# ---- V91's construction, bodies verbatim (nested in V91 main) -----------
SEL = [(math.pi / 2, math.pi / 2, 0.0),
       (math.pi / 2, 0.0, math.pi / 2),
       (0.0, math.pi / 2, math.pi / 2)]


def art(name, key):
    """A sibling verdict JSON, read by experiment name (house art())."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def fold(q):
    return tuple(round(v % TWO_PI, 9) for v in q)


def main():
    t0 = time.time()
    checks, data = [], {}

    G48 = V114.group_elements()
    GH24 = V114.chiral_elements()
    assert len(G48) == 48 and len(GH24) == 24

    def matrix(gp):
        """V114's action (x -> (s_k x_{p_k})) as a 3x3 signed permutation
        matrix; identical to V91's cub construction R[i, p] = s."""
        p, s = gp
        R = np.zeros((3, 3))
        for i in range(3):
            R[i, p[i]] = float(s[i])
        return R

    # ---- (a) the 12 doubler nodes under both groups ---------------------
    nodes = [fold(q) for q in SC.M_POINTS]

    def image(q, R):
        Rq = fold(tuple(np.dot(R, q)))
        j = min(range(len(nodes)),
                key=lambda k: sum(abs(a - b) for a, b in zip(nodes[k], Rq)))
        d = sum(abs(a - b) for a, b in zip(nodes[j], Rq))
        return j if d < 1e-9 else None

    closed_all = all(image(q, matrix(g)) is not None for q in nodes
                     for g in G48)

    def node_orbits(gens):
        mats = [matrix(g) for g in gens]
        seen, orbits = set(), []
        for i0 in range(len(nodes)):
            if i0 in seen:
                continue
            orb, stack = {i0}, [i0]
            while stack:
                i = stack.pop()
                for R in mats:
                    j = image(nodes[i], R)
                    if j is not None and j not in orb:
                        orb.add(j)
                        stack.append(j)
            seen |= orb
            orbits.append(sorted(orb))
        return [len(o) for o in orbits]

    orb_nodes_48 = node_orbits(G48)
    orb_nodes_24 = node_orbits(GH24)

    # V91's triangle-orbit re-derivation (its verbatim body) --------------
    cub = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, (p, s) in enumerate(zip(perm, signs)):
                R[i, p] = s
            cub.append(R)
    sel_idx = []
    tri_ok = True
    for q in SEL:
        j = image(fold(q), np.eye(3))
        if j is None:
            tri_ok = False
        sel_idx.append(j)
    tris = set()
    for R in cub:
        img = tuple(sorted(image(q, R) for q in [nodes[i]
                                                for i in sel_idx]))
        if all(v is not None for v in img):
            tris.add(img)

    a_orb = art("V91_count3_forcing", "data.orbit_sizes")
    a_tri = int(art("V91_count3_forcing", "data.triangle_orbit_size"))
    t0_v91 = bool(a_orb == [12] and a_tri == 8
                  and orb_nodes_48 == [12] and len(tris) == 8
                  and closed_all and tri_ok)

    # ---- (b) the 6 chiral stabilizer elements (V114 body verbatim) -----
    pts = list(itertools.product(range(4), repeat=3))
    idx = {p: i for i, p in enumerate(pts)}

    def act(x, g):
        p, s = g
        return tuple((s[k] * x[p[k]]) % 4 for k in range(3))

    v76 = np.array([math.cos(math.pi * (x[1] - x[2]) / 2)
                    + math.cos(math.pi * (x[0] - x[2]) / 2)
                    + math.cos(math.pi * (x[0] - x[1]) / 2)
                    for x in pts], float)
    stab = []
    for g in G48:
        wp = np.array([v76[idx[act(x, g)]] for x in pts])
        if np.abs(wp - v76).max() < 1e-12:
            stab.append(g)
    stab_chiral = [g for g in stab if g in GH24]

    a_stab = int(art("V114_junction_r", "data.stabilizer_size"))
    a_stabc = int(art("V114_junction_r", "data.stabilizer_chiral"))
    t0_v114 = bool(a_stab == 12 and a_stabc == 6
                   and len(stab) == 12 and len(stab_chiral) == 6)

    # conjugation orbits of the 6, by the 48 group
    def canonical(R):
        return tuple(tuple(round(float(v), 9) for v in row)
                     for row in R)

    mats48 = [matrix(g) for g in G48]
    key_to_mat = {canonical(R): R for R in mats48}
    orbits_stab = []
    seen = set()
    for h0 in stab_chiral:
        key = canonical(matrix(h0))
        if key in seen:
            continue
        orb, stack = {key}, [key]
        while stack:
            h = key_to_mat[stack.pop()]
            for g in mats48:
                c = canonical(g @ h @ g.T)
                if c not in orb:
                    orb.add(c)
                    stack.append(c)
        seen |= orb
        orbits_stab.append(len(orb))

    # ---- (c) the 16-level step permutation at the five clocks ----------
    p16, _m16 = CL.weyl_spinor()
    p16 = np.asarray(p16).astype(complex)

    a_grid = art("V239_two_dial_step_family",
                 "data.closure.grid_labels_by_order")
    t0_v239 = bool(a_grid == {"2": 1, "4": 0, "8": 0, "16": 0})

    rows_c = []
    for g in SCAN:
        U = CE.cell_step_np(g10=g)
        U32 = np.asarray(U)[:32, :32]
        B = p16.conj().T @ U32 @ p16
        mmax = np.abs(B).max(axis=0)
        second = np.sort(np.abs(B), axis=0)[-2]
        monomial = bool(mmax.min() > 1.0 - 1e-9
                        and second.max() < 1e-9)
        row = {"clock": g, "monomial": monomial,
               "min_col_max": float(mmax.min()),
               "max_col_second": float(second.max())}
        if monomial:
            pi = np.argmax(np.abs(B), axis=0)
            seen_l, ors = set(), []
            for i0 in range(16):
                if i0 in seen_l:
                    continue
                orb, x = set(), i0
                while x not in orb:
                    orb.add(x)
                    seen_l.add(x)
                    x = int(pi[x])
                ors.append(len(orb))
            row["permutation_orbit_sizes"] = sorted(ors)
        rows_c.append(row)

    # ---- the census and forced-3 decision -------------------------------
    rows = {
        "nodes12_x_G48": {"orbit_sizes": orb_nodes_48},
        "nodes12_x_chiral24": {"orbit_sizes": orb_nodes_24},
        "stabilizer6_x_conj_G48": {"orbit_sizes": sorted(orbits_stab)},
    }
    has3 = [any(s == 3 for s in r["orbit_sizes"])
            for r in rows.values()]
    for row in rows_c:
        if row["monomial"]:
            rows["levels16_x_step_k=%r" % row["clock"]] = {
                "orbit_sizes": row["permutation_orbit_sizes"]}
            has3.append(any(s == 3
                            for s in row["permutation_orbit_sizes"]))
    forced3 = bool(has3) and all(has3)
    data["census"] = {"rows": rows, "levels16_rows": rows_c,
                      "row_has_size3_orbit": has3,
                      "closed_clock_pi_over_4":
                          float(CFG.angle.branches.closed),
                      "v239_grid_labels_by_order": a_grid}
    hist = {}
    for r in rows.values():
        for s in r["orbit_sizes"]:
            hist[str(s)] = hist.get(str(s), 0) + 1
    data["orbit_size_histogram"] = hist

    data["t0"] = {"v91_orbit_sizes_artifact": a_orb,
                  "nodes_orbit_sizes_rederived": orb_nodes_48,
                  "v91_triangle_orbit_artifact": a_tri,
                  "triangle_orbit_rederived": len(tris),
                  "cubic_action_closed": closed_all,
                  "triangle_is_node_set": tri_ok,
                  "v114_stabilizer_artifact": a_stab,
                  "v114_stabilizer_chiral_artifact": a_stabc,
                  "stabilizer_rederived": len(stab),
                  "stabilizer_chiral_rederived": len(stab_chiral)}
    checks.append({"name": "T0 V91 leaves {} (triangle {}) reproduced by "
                           "the verbatim construction ({}, {}) and V114 "
                           "stabilizer {}/{} re-derived as {}/{}; V239 "
                           "closure labels {}".format(
                               a_orb, a_tri, orb_nodes_48, len(tris),
                               a_stab, a_stabc, len(stab),
                               len(stab_chiral), a_grid),
                   "pass": bool(t0_v91 and t0_v114 and t0_v239)})
    if not (t0_v91 and t0_v114 and t0_v239):
        data["status"] = "MACHINERY-DRIFT"
        print("V330 MACHINERY-DRIFT: the V91/V114/V239 anchors do not "
              "reproduce; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    branch = ("GEN3-ORBIT-FORCED" if forced3
              else "GEN3-NOT-FORCED-DECLARED")
    data["status"] = branch
    data["wall_time"] = time.time() - t0
    checks.append({"name": "branch {}: orbit-size histogram {} over rows "
                           "{}; a size-3 orbit present per row: {}"
                   .format(branch, hist, list(rows), has3),
                   "pass": True})
    print("V330 {}: histogram {}; nodes rows {} / {}; stab orbits {}; "
          "monomial clocks {}; wall {:.0f} s.".format(
              branch, hist, orb_nodes_48, orb_nodes_24, orbits_stab,
              [r["clock"] for r in rows_c if r["monomial"]],
              data["wall_time"]), flush=True)
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V330_generation_forcing_census", "tier": 3,
           "property": "does any booked group-action object carry an "
                       "orbit of size exactly 3 in every row - is the "
                       "family number 3 forced anywhere (tier 3, O3)",
           "residual": True,
           "programme": "theory-hunt V330",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v330_generation_forcing_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V330 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
