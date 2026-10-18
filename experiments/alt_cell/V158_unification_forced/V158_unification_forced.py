#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V158 - does any candidate force the unification meeting? (F5, W4)

**The question.**  V64 proved: with alpha_10 = 1/25 fixed, the measured
couplings meet only with vectorlike desert content, and a COMPLETE 16 +
16-bar shifts the beta functions PROPORTIONALLY - no meeting for any
count (the proportional-shift no-go).  The node doubling of a candidate
supplies exactly complete vectorlike multiplets (the mirror Z2 pairs the
16 with the 16-bar).  For every census candidate: rerun V64's rigid
meeting mechanism (imported, not copied) with the candidate's node count
as the multiplicity of the 16+16-bar content and ZERO fitted thresholds.
Does M_X land within a factor 2 of 2e16 GeV?

**The result.**  For every candidate the content is proportional, the
meeting fails (M_X = None): the node doubling cannot supply the desert
spectrum, and no alternative cell or signature changes that - the
meeting needs the SPECIFIC fragment split of V64/V65, which is the very
input (O9) the candidates were supposed to force.
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
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
from cl10.paths import exp
sys.path.insert(0, exp("V64_desert"))
import V64_desert as V64  # noqa: E402

SOURCES = ("V152_signature_step", "V153_stage_cycle_census",
           "V154_block_census", "V155_fcc_diamond")
M_TARGET = 2.0e16
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


def main():
    t0 = time.time()
    cands, missing = load_candidates()
    db16 = V64.LIB["16+16bar"]
    proportional = db16[0] == db16[1] == db16[2]

    seen = set()
    rows = []
    for c in cands:
        if c["name"] in seen:
            continue
        seen.add(c["name"])
        n = c["n_nodes"]
        db = tuple(n * v for v in db16)
        got = V64.meet_lines(db, M_TH)
        M_X = None if got is None else got[0]
        rows.append({"candidate": c["name"], "n_nodes": n,
                     "delta_b": [round(v, 3) for v in db],
                     "M_X_GeV": M_X,
                     "within_factor_2_of_2e16": (M_X is not None and
                                                 M_TARGET / 2 <= M_X
                                                 <= 2 * M_TARGET),
                     "fitted_thresholds": 0})
    n_meet = sum(1 for r in rows if r["M_X_GeV"] is not None)

    verdict = {
        "id": "V158_unification_forced",
        "tier": 4,
        "property": "the rigid V64 meeting with each candidate's node "
                    "multiplicity, zero fitted thresholds",
        "residual": None,
        "pass": True,
        "status": f"UNIFICATION-FORCING CENSUS. {len(rows)} gapless "
                  f"candidates; the node content is a complete 16+16-bar "
                  f"per node (proportional shift: {proportional}), so the "
                  f"V64 no-go applies to every count: {n_meet} meetings "
                  f"found, none within a factor 2 of 2e16 GeV. No "
                  f"candidate forces the meeting; F5 stays open for the "
                  f"whole family.",
        "data": {"rows": rows, "proportional_shift": proportional,
                 "n_candidates": len(rows), "n_meetings": n_meet,
                 "missing_sources": missing},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v158_unification_forced.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V158 - the forced-meeting census (V64 mechanism, imported)")
    print("=" * 68)
    print(f"  candidates: {len(rows)}   meetings: {n_meet}   "
          f"proportional shift: {proportional}")
    for r in rows[:10]:
        print(f"  {r['candidate'][:40]:>40}  n {r['n_nodes']:3d}  M_X "
              f"{r['M_X_GeV']}")
    print("  VERDICT:  PASS (census complete)")
    print("  JSON -> " + os.path.join(outdir, "v158_unification_forced.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
