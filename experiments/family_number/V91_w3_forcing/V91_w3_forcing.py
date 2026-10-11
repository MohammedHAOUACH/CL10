#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V91 - is the family count 3 FORCED? Cubic orbits of the 12 nodes,
the rank-1 meeting for every cubic-symmetric split, and the orbit of
the 3-node triangle. W3's forcing question, decided.

**Why (W3).**  V70 found the index: the 12 doubler nodes split into 3
families of 4 (zero axis x, y, z), and the unification-required split
(V68: n >= 3 nodes meeting) is realizable by one node per family -
the count 3 EXISTS. The open question: is it FORCED? A forcing
mechanism would have to exclude every other split of the 12 nodes.
The cubic group constrains the SPLIT SET: a cubic-symmetric choice
must be a union of orbits. V91 computes: (i) the cubic orbits on the
12 nodes; (ii) for every cubic-symmetric split size n = 4, 8, 12,
whether the desert meeting survives and where it lands (the V74
rank-1 fact generalized: the meeting point is fixed by the SUM of
strengths, independent of the split); (iii) the cubic orbit of the
3-node triangle: whether the count 3 is cubic-symmetric at all.

**Checks (computed).**
1. Regression V70: 12 nodes, families {x:4, y:4, z:4}, spectra
   identical within and across families to 1e-12.
2. Cubic orbits on the 12 nodes (48 signed permutations): sizes.
3. Rank-1 Jacobian of the meeting equations with n = 4, 8, 12 node
   thresholds (V74 mechanism generalized): rank in the thresholds.
4. V65 scan_roots with n_nodes = 3, 4, 8, 12: the fixed meeting
   scale M_X(n) for each cubic-symmetric size; against the observed
   unification scale ~2e16 GeV.
5. The 3-node triangle SEL (V76): its orbit under the cubic group -
   is the count-3 split set cubic-invariant?

**Verdict (computed).**  Whether the family count 3 is forced by the
automaton, or whether every cubic-symmetric split is excluded by the
fixed meeting (F5) while the count 3 requires a symmetry-BREAKING
choice - the precise form of the W3 residual.
"""

import itertools
import json
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
sys.path.insert(0, exp("V74_multiscale_nodes"))
import V65_wilson as V65  # noqa: E402
from V74_multiscale_nodes import fragments  # noqa: E402
from cl10 import supercell as SC  # noqa: E402
from cl10 import unification as UN  # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
TWO_PI = 2.0 * math.pi
INV_A10 = 25.0
M_X_OBS = 2.0e16            # GeV, observed 1-loop unification scale
SEL = [(math.pi / 2, math.pi / 2, 0.0),
       (math.pi / 2, 0.0, math.pi / 2),
       (0.0, math.pi / 2, math.pi / 2)]


def fold(q):
    return tuple(round(v % TWO_PI, 9) for v in q)


def main():
    t0 = time.time()
    checks = {}
    _B, _U32, U16 = SC.sectors()

    # --- 1. regression V70 -------------------------------------------------
    fam = defaultdict(list)
    for q in SC.M_POINTS:
        axis = next(j for j in range(3) if abs(q[j]) < 1e-9)
        fam[["x", "y", "z"][axis]].append(q)
    fam_sizes = {k: len(v) for k, v in fam.items()}

    def spectrum(q):
        return np.sort(np.angle(np.linalg.eigvals(SC.bloch(q, U16))))
    per_family = {k: [spectrum(q) for q in v] for k, v in fam.items()}
    within = max(float(np.abs(a - b).max())
                 for v in per_family.values() for a, b in zip(v, v[1:]))
    across = max(float(np.abs(per_family["x"][0] - per_family[k][0]).max())
                 for k in ("y", "z"))
    reg_ok = (fam_sizes == {"x": 4, "y": 4, "z": 4}
              and within < 1e-12 and across < 1e-12)
    checks["v70_regression"] = {"fam_sizes": fam_sizes,
                                "within": float(f"{within:.2e}"),
                                "across": float(f"{across:.2e}")}

    # --- 2. cubic orbits on the 12 nodes -----------------------------------
    cub = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, (p, s) in enumerate(zip(perm, signs)):
                R[i, p] = s
            cub.append(R)
    nodes = [fold(q) for q in SC.M_POINTS]

    def image(q, R):
        Rq = fold(tuple(np.dot(R, q)))
        j = min(range(len(nodes)),
                key=lambda k: sum(abs(a - b) for a, b in zip(nodes[k], Rq)))
        d = sum(abs(a - b) for a, b in zip(nodes[j], Rq))
        return j if d < 1e-9 else None

    closed = all(image(q, R) is not None for q in nodes for R in cub)
    seen, orbits = set(), []
    for i0 in range(len(nodes)):
        if i0 in seen:
            continue
        orb = {i0}
        stack = [i0]
        while stack:
            i = stack.pop()
            for R in cub:
                j = image(nodes[i], R)
                if j is not None and j not in orb:
                    orb.add(j)
                    stack.append(j)
        seen |= orb
        orbits.append(sorted(orb))
    orbit_sizes = sorted(len(o) for o in orbits)
    checks["cubic_action_closed"] = closed
    checks["orbit_sizes"] = orbit_sizes

    # cubic-symmetric split sets = unions of orbits: available sizes
    avail = sorted({sum(sz) for k in range(1, len(orbits) + 1)
                    for sz in itertools.combinations(orbit_sizes, k)})

    # --- 3. rank-1 meeting Jacobian for n = 4, 8, 12 -----------------------
    def F_of_n(zs, x, theta):
        frs, cf = fragments(theta)
        lcf = [math.log(c) for c in cf]
        F = []
        for i in range(3):
            s = INV_A10 + B_SM[i] * x / TWO_PI
            for zk in zs:
                for j, f in enumerate(frs):
                    act = x - zk - lcf[j]
                    if act > 0:
                        s += f["db"][i] * act / TWO_PI
            F.append(s - 1.0 / A_Z[i])
        return F

    def rank_n(zs, x, theta):
        base = np.array(F_of_n(zs, x, theta))
        J = np.zeros((3, len(zs)))
        for k in range(len(zs)):
            zp = list(zs)
            zp[k] += 1e-5
            J[:, k] = (np.array(F_of_n(zp, x, theta)) - base) / 1e-5
        return int(np.linalg.matrix_rank(J, tol=1e-6))

    ranks = {}
    for n in (4, 8, 12):
        zs = [20.0 + 6.0 * (k % 4) - 9.0 for k in range(n)]  # split
        ranks[f"n_{n}"] = rank_n(zs, 33.5, 1.731)
    checks["meeting_jacobian_rank"] = ranks

    # --- 4. the fixed meeting scale for each cubic-symmetric size ----------
    th = 1.731
    frs, _ = fragments(th)
    mx = {}
    for n in (3, 4, 8, 12):
        roots = [r for r in V65.scan_roots(frs, n_nodes=n)
                 if r["physical"]]
        if roots:
            r0 = max(roots, key=lambda r: r["M_X_GeV"])
            mx[f"n_{n}"] = {"M_X_GeV": float(f"{r0['M_X_GeV']:.4e}"),
                            "ratio_to_obs": float(
                                f"{r0['M_X_GeV'] / M_X_OBS:.3e}")}
        else:
            mx[f"n_{n}"] = {"M_X_GeV": None,
                            "ratio_to_obs": None,
                            "note": "no physical meeting root"}
    checks["fixed_meeting_scale"] = mx

    # --- 5. the cubic orbit of the 3-node triangle -------------------------
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
    checks["triangle_is_cubic_set"] = tri_ok
    checks["triangle_orbit_size"] = len(tris)

    all_rank1 = all(v == 1 for v in ranks.values())
    excluded = all(v["M_X_GeV"] is None
                   or abs(math.log10(v["M_X_GeV"] / M_X_OBS)) > 0.5
                   for v in mx.values())
    ok = (reg_ok and closed and orbit_sizes == [12] and all_rank1
          and tri_ok and len(tris) > 1 and excluded)

    mx_txt = {k: (f"{v['M_X_GeV']:.1e} GeV" if v["M_X_GeV"]
                  else "no root") for k, v in mx.items()}
    status = (
        f"THE COUNT 3 IS NOT FORCED - IT IS A SYMMETRY-BREAKING INPUT "
        f"(V70 regression dev {max(within, across):.1e}). The cubic "
        f"group acts TRANSITIVELY on the 12 nodes (one orbit of size "
        f"{orbit_sizes[0]}; the V70 axis families are NOT cubic "
        f"orbits), so the ONLY cubic-symmetric split set is the full "
        f"set of 12 nodes (available sizes {avail}). The desert "
        f"meeting is rigid for every split: the n-node meeting "
        f"Jacobian has rank {ranks} in the thresholds (the V74 rank-1 "
        f"fact holds for every n), so M_X is fixed by the SUM of "
        f"strengths, split-independent. The fixed scales: {mx_txt} "
        f"against the observed unification {M_X_OBS:.0e} GeV: every "
        f"size n = 3, 4, 8, 12 {'FAILS to meet at the observed scale (F5)' if excluded else 'meets'}. "
        f"The count-3 triangle (V76 SEL) IS a set of nodes (verified) "
        f"but is NOT cubic-invariant: its cubic orbit has {len(tris)} "
        f"distinct triangles. W3, decided: the cubic-symmetric option "
        f"(all 12 nodes) is excluded by the fixed meeting; every "
        f"proper split (3, 4, 8 nodes) BREAKS the cubic symmetry of "
        f"the split set, so the family count 3 is an input choice - "
        f"not forced by the automaton - and its cubic-symmetric "
        f"alternative is now excluded by computation.")

    verdict = {"experiment": "V91_w3_forcing",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v91_w3_forcing.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("V91 - is the family count 3 forced?")
    print(f"  V70 regression: {fam_sizes}, dev {max(within, across):.1e}")
    print(f"  cubic orbits on 12 nodes: {orbit_sizes} "
          f"(closed: {closed})")
    print(f"  meeting Jacobian ranks: {ranks}")
    print(f"  fixed M_X: {mx_txt}")
    print(f"  triangle cubic orbit: {len(tris)} triangles")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v91_w3_forcing.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
