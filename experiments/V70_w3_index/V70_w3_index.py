#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V70 - the W3 audit: is there an index of the step that equals 3?

**The question.**  The root gives ONE generation; the theory identifies the
family number with the 3 spatial directions (a stated hypothesis, W3). Is
there an index of the IMPLEMENTED step that equals 3 and could force it?

**The audit.**  The 12 doubler nodes of V63 organize as 3 FAMILIES OF 4:
the zero component of the node momentum picks an axis -
  F_z = {(+-pi/2, +-pi/2, 0)},  F_y,  F_x  - one family per spatial
direction. Tests:
1. the cubic group permutes the families transitively (the Floquet spectra
   of the 4 nodes of one family are identical, and the three families carry
   the same spectrum);
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
domain-wall input). W3 remains a hypothesis - now with a named index (the
3 axis families of doubler nodes) rather than a bare identification.
"""

import json
import math
import os
import sys
import time

from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import brillouin as BR  # noqa: E402


def main():
    t0 = time.time()
    _B, _U32, U16 = BR.blocks()

    # 1. families by zero axis; spectra identical within and across
    fam = defaultdict(list)
    for q in BR.DOUBLERS:
        axis = next(j for j in range(3) if abs(q[j]) < 1e-9)
        key = ["x", "y", "z"][axis]
        fam[key].append(q)
    fam = dict(fam)
    fam_sizes = {k: len(v) for k, v in fam.items()}

    def spectrum(q):
        return np.sort(np.angle(np.linalg.eigvals(BR.floquet(q, U16))))
    per_family = {k: [spectrum(q) for q in v] for k, v in fam.items()}
    within = max(float(np.abs(a - b).max())
                 for v in per_family.values() for a, b in zip(v, v[1:]))
    across = max(float(np.abs(per_family["x"][0] - per_family[k][0]).max())
                 for k in ("y", "z"))
    cubic = within < 1e-12 and across < 1e-12

    # 2. one node per family is cubic-degenerate under the Wilson factor:
    # the Wilson phase at a node is w * sum(1 - cos q_j) = 2w for EVERY
    # node, so any 3 chosen nodes (one per family) get the SAME splitting
    # at a single w - consistent with the V65 n_nodes = 3 roots.
    phase_at = {q: sum(1.0 - math.cos(qj) for qj in q) for q in BR.DOUBLERS}
    all_same = len(set(round(v, 12) for v in phase_at.values())) == 1

    # 3. the counts
    counts = {"nodes": len(BR.DOUBLERS), "families": len(fam),
              "per_family": 4, "bond_period": 3}

    ok = cubic and all_same and fam_sizes == {"x": 4, "y": 4, "z": 4}
    verdict = {
        "experiment": "V70_w3_index",
        "pass": ok,
        "status": ("W3 AUDIT: the 12 doubler nodes split into 3 families of "
                   "4 labeled by the zero axis (x, y, z); the cubic group "
                   "acts transitively (identical Floquet spectra, dev "
                   f"{max(within, across):.1e}); the Wilson phase is equal "
                   "at all nodes, so the n_nodes = 3 unification solutions "
                   "of V65/V68 are realizable by ONE node per axis family "
                   "at a single w. The index 3 EXISTS (3 axis families = 3 "
                   "directions = bond period t mod 3) and labels the 3 "
                   "copies of desert matter; it does NOT force the copy "
                   "number (4 per family, 12 total are equally available). "
                   "W3 stays a hypothesis, now with a named index."),
        "data": {"families": {k: [list(q) for q in v]
                              for k, v in fam.items()},
                 "spectra_identical_within": within,
                 "spectra_identical_across": across,
                 "wilson_phase_equal_all_nodes": all_same,
                 "counts": counts},
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v70_w3_index.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V70 - the W3 audit: an index of the step equal to 3?")
    print("=" * 68)
    print(f"  node families by zero axis: {fam_sizes}")
    print(f"  spectra identical within {within:.1e}, across {across:.1e} "
          f"-> cubic transitive: {cubic}")
    print(f"  Wilson phase equal at all 12 nodes: {all_same}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - {verdict['status']}")
    print("  JSON -> " + os.path.join(outdir, "v70_w3_index.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
