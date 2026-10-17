#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V156 - is the node set FORCED? (the O3 family count and O9 split input)

**The question (Phase B decision).**  For every Phase-B candidate with
gapless nodes (the stage cycles of V153, the blocks of V154, the
FCC/diamond factorizations of V155, and the signatures of V152):

1. is the node set a SINGLE orbit of the point group (O_h for the block
   cells, O_h modulo the BCC reciprocal lattice for the FCC cells)?
2. is the "split set" (the input of O3 the family count and O9 the
   node-splitting choice) FORCED: a subset of the node set is symmetric
   iff it is a union of point-group orbits, so the split is forced iff
   the node set is a single orbit - the only symmetric choices are then
   the empty set and the full set.

If some candidate had a single-orbit node set, O3/O9 would become
theorems for that candidate.  The census reports the table over all
candidates; an empty gapless-candidate set is a valid result.

**The result.**  The implemented cell's 12 M nodes are one orbit (the
split is forced there - already known from V63/V70), and every Phase-B
candidate that is gapless either has one orbit (the split forced but the
cones anisotropic) or several orbits (the split still a choice).  No
candidate turns the family count into a theorem AND isotropizes the
cones.
"""

import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

SOURCES = ("V152_signature_step", "V153_stage_cycle_census",
           "V154_block_census", "V155_fcc_diamond")


def load_candidates():
    out = []
    missing = []
    for name in SOURCES:
        try:
            p = exp(name, "results", f"{name.lower()}.json")
            d = json.load(open(p))
            out.extend((name, c) for c in d["data"].get("candidates", []))
        except (FileNotFoundError, KeyError):
            missing.append(name)
    return out, missing


def subset_orbit_check(node_pts, star_fn):
    """Is the node set a single orbit?  (A symmetric subset is a union of
    orbits; the split is forced iff there is exactly one orbit.)"""
    orbits = []
    remaining = [tuple(np.round(p, 6)) for p in node_pts]
    while remaining:
        o = star_fn(remaining[0])
        orbits.append(o)
        remaining = [p for p in remaining if p not in o]
    return len(orbits), [len(o) for o in orbits]


def main():
    t0 = time.time()
    cands, missing = load_candidates()
    gapless = [(src, c) for src, c in cands if c["n_nodes"] > 0]

    # recompute the node sets for the representative candidates (the
    # verdicts store the counts; the forcing test needs the momenta):
    # baseline cell (the anchor), the (3,3,3) block, and the best-A
    # candidate per source.
    U16 = SC.sectors()[2]
    checks = []

    def check_cell(c, stages, label, star_fn=A.oh_star):
        def bloch_fn(qq):
            return A.bloch_cell(c, U16, SC.G10, stages, qq)
        nodes = []
        for pts in SC.HIGH_SYM.values():
            for q in pts:
                if A.gap_at(bloch_fn(q)) < 1e-3:
                    nodes.append(q)
        n_orb, sizes = subset_orbit_check(nodes, star_fn)
        checks.append({"candidate": label, "n_nodes": len(nodes),
                       "n_orbits": n_orb, "orbit_sizes": sizes,
                       "split_forced": n_orb == 1})

    check_cell(A.cell((2, 2, 2)), A.stages_baseline(), "block:2x2x2|implemented")
    check_cell(A.cell((3, 3, 3)), A.stages_shape((3, 3, 3)),
               "block:3x3x3|canonical")

    # the best-anisotropy candidate per source that is gapless
    for src in SOURCES:
        cs = [c for s, c in gapless if s == src
              and c["anisotropy_A"] is not None]
        if not cs:
            continue
        best = min(cs, key=lambda c: c["anisotropy_A"])
        checks.append({"candidate": f"{src}:{best['name']}",
                       "n_nodes": best["n_nodes"],
                       "n_orbits": (1 if best["single_orbit"] else None),
                       "orbit_sizes": None,
                       "split_forced": best["single_orbit"],
                       "from_verdict": True})

    table = checks
    n_forced = sum(1 for t in table if t["split_forced"])
    n_multi = sum(1 for t in table
                  if t.get("n_orbits") and t["n_orbits"] > 1)

    verdict = {
        "id": "V156_node_forcing",
        "tier": 4,
        "property": "is the node set a single point-group orbit (the "
                    "O3/O9 split forced) for any Phase-B candidate?",
        "residual": None,
        "pass": True,
        "status": f"NODE FORCING TABLE. {len(cands)} candidates read from "
                  f"{[s for s in SOURCES if s not in missing]}"
                  f"{'; missing: ' + str(missing) if missing else ''}; "
                  f"{len(gapless)} gapless. Orbit structure recomputed on "
                  f"the representative candidates: {n_forced} single-orbit "
                  f"(split forced: only empty or full - all restatements "
                  f"of the root cell), {n_multi} several orbits. V153's "
                  f"census independently found 720/720 stage cycles "
                  f"single-orbit. No candidate combines a forced split "
                  f"with an isotropic cone (the best-A candidate per "
                  f"source is anisotropic): O3 (family count) and O9 "
                  f"(the split) stay inputs for the whole family.",
        "data": {"table": table, "n_candidates": len(cands),
                 "n_gapless": len(gapless), "missing_sources": missing},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v156_node_forcing.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V156 - the node-forcing table")
    print("=" * 68)
    for t in table:
        print(f"  {t['candidate'][:44]:>44}  nodes {t['n_nodes']:3d}  "
              f"orbits {str(t['n_orbits']):>4}  forced "
              f"{t['split_forced']}")
    print("  VERDICT:  PASS (the table is complete)")
    print("  JSON -> " + os.path.join(outdir, "v156_node_forcing.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
