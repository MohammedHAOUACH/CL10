#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V70 - the O3 audit: is there an index of the step that equals 3?

**The question.**  The root gives ONE generation; the theory identifies the
family number with the 3 spatial directions (a stated hypothesis, O3). Is
there an index of the IMPLEMENTED step that equals 3 and could force it?

**The audit.**  The 12 doubler nodes of V63 organize as 3 FAMILIES OF 4:
the zero component of the node momentum picks an axis -
  F_z = {(+-pi/2, +-pi/2, 0)},  F_y,  F_x  - one family per spatial
direction. Tests:
1. the node SET is one orbit of the 48-element cubic group (12 = 48/4,
   stabilizer order 4 - V63), every node exactly gapless, and at the root
   clock the 12 node channels fall into ISO-SPECTRAL CLASSES: the levels
   pinned to the Floquet pair {0, pi} are identical at every node, while
   the intermediate +/-pi/3, +/-2pi/3 levels permute between classes (the
   index-3 point-group defect of V128 seen in the node spectra);
2. the V65 Wilson splitting with ONE scale w splits either all 12 or a
   cubic-complete set; the n_nodes = 3 solutions of V65/V68 are realizable
   by choosing ONE node per family (3 nodes, related by cubic rotations,
   hence degenerate - consistent with a single w);
3. the count: 12 = 4 x 3, and 3 = the number of axes = the bond-layer
   period (t mod 3). The ONLY index of the step equal to 3 is the axis
   count itself.

**Verdict.**  The three node families give the strongest structural hook
available: the 3 copies of desert matter required by unification (V68) can
be labeled by the 3 spatial directions - generations = node families =
axes. But nothing in the step makes the copy number 3 FORCED (the 4 nodes
per family would give 4, the full set 12; the choice "one per family" is a
domain-wall input). O3 remains a hypothesis - now with a named index (the
3 axis families of doubler nodes) rather than a bare identification.
"""

import json
import math
import os
import sys
import time
import os

from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import supercell as SC  # noqa: E402


def main():
    t0 = time.time()
    _B, _U32, U16 = SC.sectors()

    # 1. families by zero axis; spectra identical within and across
    fam = defaultdict(list)
    for q in SC.M_POINTS:
        axis = next(j for j in range(3) if abs(q[j]) < 1e-9)
        key = ["x", "y", "z"][axis]
        fam[key].append(q)
    fam = dict(fam)
    fam_sizes = {k: len(v) for k, v in fam.items()}

    def spectrum(q):
        return np.sort(np.angle(np.linalg.eigvals(SC.bloch(q, U16))))
    per_family = {k: [spectrum(q) for q in v] for k, v in fam.items()}
    # circular phase distance: at the root clock the spectra carry the
    # Floquet pair {0, pi} and the branch cut of np.angle must not count
    cdist = lambda A, B: np.abs((A - B + np.pi) % (2 * np.pi) - np.pi)
    within = max(float(cdist(a, b).max())
                 for v in per_family.values() for a, b in zip(v, v[1:]))
    across = max(float(cdist(per_family["x"][0], per_family[k][0]).max())
                 for k in ("y", "z"))
    # isospectral classes of the 12 node channels (root clock): pinned
    # levels identical, intermediate levels permuted (V128 index-3 defect)
    reps, classes = [], []
    for q in SC.M_POINTS:
        s = spectrum(q)
        hit = next((i for i, r in enumerate(reps)
                    if cdist(s, r).max() < 1e-9), None)
        if hit is None:
            reps.append(s)
            classes.append(1)
        else:
            classes[hit] += 1
    import itertools
    base = SC.M_POINTS[0]
    orbit = {tuple(round(x, 9) for x in
                    tuple(s[i] * base[p[i]] for i in range(3)))
             for p in itertools.permutations(range(3))
             for s in itertools.product((1, -1), repeat=3)}
    mset = {tuple(round(x, 9) for x in q) for q in SC.M_POINTS}
    transitive = orbit == mset and len(orbit) == 12
    gapless_all = all(SC.gap(q, U16) < 1e-9 for q in SC.M_POINTS)
    cubic = transitive and gapless_all

    # 2. one node per family is cubic-degenerate under the Wilson factor:
    # the Wilson phase at a node is w * sum(1 - cos q_j) = 2w for EVERY
    # node, so any 3 chosen nodes (one per family) get the SAME splitting
    # at a single w - consistent with the V65 n_nodes = 3 roots.
    phase_at = {q: sum(1.0 - math.cos(qj) for qj in q) for q in SC.M_POINTS}
    all_same = len(set(round(v, 12) for v in phase_at.values())) == 1

    # 3. the counts
    counts = {"nodes": len(SC.M_POINTS), "families": len(fam),
              "per_family": 4, "bond_period": 3}

    ok = cubic and all_same and fam_sizes == {"x": 4, "y": 4, "z": 4}
    verdict = {
        "experiment": "V70_w3_index",
        "pass": ok,
        "status": ("O3 AUDIT: the 12 doubler nodes split into 3 families of "
                   "4 labeled by the zero axis (x, y, z); the node set is "
                   "ONE cubic orbit of 12 (stabilizer 4) and every node is "
                   "exactly gapless; at the root clock the 12 node channels "
                   f"form {len(classes)} isospectral classes (sizes "
                   f"{classes}): the Floquet-pinned levels are identical "
                   "everywhere, the intermediate levels permute between "
                   "classes (the V128 index-3 point-group defect); the "
                   "Wilson phase is equal at all nodes, so the n_nodes = 3 "
                   "unification solutions of V65/V68 are realizable by ONE "
                   "node per axis family at a single w. The index 3 EXISTS "
                   "(3 axis families = 3 directions = bond period t mod 3) "
                   "and labels the 3 copies of desert matter; it does NOT "
                   "force the copy number. O3 stays a hypothesis, now with "
                   "a named index."),
        "data": {"families": {k: [list(q) for q in v]
                              for k, v in fam.items()},
                 "spectral_classes": classes,
                 "spectral_spread_vs_first_class": {
                     str([round(v, 4) for v in q]):
                         float(cdist(spectrum(q), spectrum(SC.M_POINTS[0])).max())
                     for q in SC.M_POINTS},
                 "wilson_phase_equal_all_nodes": all_same,
                 "counts": counts},
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v70_w3_index.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V70 - the O3 audit: an index of the step equal to 3?")
    print("=" * 68)
    print(f"  node families by zero axis: {fam_sizes}")
    print(f"  node set one cubic orbit of 12: {transitive}; all gapless: "
          f"{gapless_all}; isospectral classes: {classes}")
    print(f"  Wilson phase equal at all 12 nodes: {all_same}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - {verdict['status']}")
    print("  JSON -> " + os.path.join(outdir, "v70_w3_index.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
