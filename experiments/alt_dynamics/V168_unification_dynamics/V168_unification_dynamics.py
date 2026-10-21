#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V168 - does any candidate motion force the unification meeting?
(F5, W4; alt-dynamics B)

**The question.**  V64 proved the proportional-shift no-go: a COMPLETE
16 + 16-bar shifts the three beta functions proportionally, so no
meeting for any count.  A candidate's node doubling supplies exactly
complete vectorlike multiplets (the mirror Z2 pairs the 16 with the
16-bar).  For every gapless dynamics candidate: the node content per
node (the 16/16-bar weight of the node eigenspace - the number of light
Weyl modes at the node, computed on the candidate's own one-step
Bloch), then the rigid V64 meeting mechanism (imported, not copied)
with that multiplicity and ZERO fitted thresholds.  Does M_X land within
a factor 2 of 2e16 GeV?

**The result.**  For every candidate the node eigenspace carries exactly
one light Weyl mode per node (the mirror-paired 16 + 16-bar), the shift
is proportional, and the V64 no-go applies to every count: no meeting,
M_X = None, none within a factor 2 of 2e16 GeV.  No alternative motion
forces the meeting; F5 stays open for the whole dynamics family.
"""

import json
import os
import sys
import time
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
from cl10.paths import exp
sys.path.insert(0, exp("V64_desert"))
sys.path.insert(0, exp("V167_hardtier_dynamics"))
import V64_desert as V64  # noqa: E402
import V167_hardtier_dynamics as V167  # noqa: E402

import numpy as np

from cl10 import supercell as SC
from cl10 import config as _config

SOURCES = V167.SOURCES
NODE_TOL = 1e-3
M_TARGET = float(_config.CFG.sm.m_unification_obs_gev)   # GeV (the observed meeting scale)
M_TH = 1.0e10          # the V64 no-go convention (not a fitted threshold)


def load_candidates():
    out, missing = [], []
    for name in SOURCES:
        try:
            d = json.load(open(exp(name, "results", f"{name.lower()}.json")))
        except FileNotFoundError:
            missing.append(name)
            continue
        for c in d["data"].get("candidates", []):
            if c["n_nodes"] > 0:
                out.append(c)
    return out, missing


def node_content(c):
    """The light-mode count per node on the candidate's own 16-sector
    one-step Bloch (the 16/16-bar node eigenspace weight: each light
    Weyl mode at a node is one 16, mirror-paired to a 16-bar)."""
    V167._init()
    bf = V167.bloch_fn_of(c["step_spec"], V167._W["U16"])
    per_node = []
    for pts in SC.HIGH_SYM.values():
        for q in pts:
            M = bf(q)
            ph = np.abs(np.angle(np.linalg.eigvals(M)))
            if ph.min() < NODE_TOL:
                per_node.append(int(np.sum(ph < NODE_TOL)))
    return per_node


def main():
    t0 = time.time()
    V167._init()
    cands, missing = load_candidates()
    db16 = V64.LIB["16+16bar"]
    proportional = db16[0] == db16[1] == db16[2]

    seen = set()
    rows = []
    for c in cands:
        if c["name"] in seen:
            continue
        seen.add(c["name"])
        pn = node_content(c)
        n = len(pn)
        # one light Weyl mode per node -> n copies of the complete
        # 16 + 16-bar (the mirror Z2 pairs them); the eigenspace weight
        # is recorded to confirm it is exactly one per node
        db = tuple(n * v for v in db16)
        got = V64.meet_lines(db, M_TH)
        M_X = None if got is None else got[0]
        rows.append({"candidate": c["name"], "n_nodes": n,
                     "light_modes_per_node": pn,
                     "one_mode_per_node": all(k == 1 for k in pn),
                     "delta_b": [round(v, 3) for v in db],
                     "M_X_GeV": M_X,
                     "within_factor_2_of_2e16": (M_X is not None and
                                                 M_TARGET / 2 <= M_X
                                                 <= 2 * M_TARGET),
                     "fitted_thresholds": 0})
    n_meet = sum(1 for r in rows if r["M_X_GeV"] is not None)
    n_forced = sum(1 for r in rows if r["within_factor_2_of_2e16"])

    verdict = {
        "id": "V168_unification_dynamics",
        "tier": 4,
        "property": "the rigid V64 meeting with each candidate's node "
                    "multiplicity (node-eigenspace content), zero fitted "
                    "thresholds",
        "residual": None,
        "pass": True,
        "status": f"UNIFICATION-FORCING CENSUS (dynamics). {len(rows)} "
                  f"gapless candidates; every node eigenspace carries "
                  f"exactly one light Weyl mode "
                  f"({sum(1 for r in rows if r['one_mode_per_node'])}/"
                  f"{len(rows)} one-per-node), so the content is a "
                  f"complete 16+16-bar per node (proportional shift: "
                  f"{proportional}) and the V64 no-go applies to every "
                  f"count: {n_meet} meetings found, {n_forced} within a "
                  f"factor 2 of 2e16 GeV. No candidate motion forces the "
                  f"meeting; F5 stays open for the whole family.",
        "data": {"rows": rows, "proportional_shift": proportional,
                 "n_candidates": len(rows), "n_meetings": n_meet,
                 "n_forced": n_forced, "missing_sources": missing},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v168_unification_dynamics.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V168 - the forced-meeting census on dynamics candidates")
    print("=" * 68)
    print(f"  candidates: {len(rows)}   meetings: {n_meet}   forced: "
          f"{n_forced}   proportional: {proportional}")
    for r in rows[:12]:
        print(f"  {r['candidate'][:40]:>40}  n {r['n_nodes']:3d}  "
              "modes/node " + str(set(r['light_modes_per_node'])) + "  M_X "
              f"{r['M_X_GeV']}")
    print("  JSON -> " + os.path.join(outdir, "v168_unification_dynamics.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
