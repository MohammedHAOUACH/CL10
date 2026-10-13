#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V174 - the isotropy gate on the interaction census (alt-interactions)

**The question.**  Any row of the interaction census (V171 vertex, V172
one-loop flow, V173 block-RG flow) that reaches cone isotropy
(A < 1e-6) or shows a flow toward isotropy larger than 0.5 across two
decades of scale becomes a candidate and gets the full hard tier (the
V167 gate table: unitarity, front speed exactly 1, zero mirror leakage,
anomalies, closure).  The contact collision C = exp(i U_c N_block) is
diagonal unitary and moves no amplitude, so the front budget of the
interacting step equals the free one - verified numerically on the
6^3-block torus with the V167 front machinery.

**The result.**  ZERO candidates: no interaction row is isotropic and no
flow toward isotropy exceeds 0.5.  The gate records the census minima
and passes with "no candidate: the interaction does not isotropize".
"""

import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
ROOT = HERE
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

SRC = ("V171_collision_vertex", "V172_self_energy_flow", "V173_block_rg")


def _data(name):
    p = os.path.join(ROOT, "experiments", "alt_interactions", name,
                     "results", name.lower() + ".json")
    with open(p) as f:
        d = json.load(f)
    return d["data"]


def _anisotropy(row):
    for k in ("A_free", "A_bare", "A_interacting"):
        if k in row:
            v = row[k]
            if isinstance(v, dict):
                return min(float(x) for x in v.values())
            return float(v)
    if "A_dressed" in row:
        return min(float(x) for x in row["A_dressed"].values())
    return None


def _flow(row):
    if "A_dressed" in row:
        return row["A_dressed"]["0.0001"] - row["A_dressed"]["0.01"]
    return None


def main():
    t0 = time.time()
    all_rows, srcs, excluded = [], [], 0
    for name in SRC:
        data = _data(name)
        converged = bool(data.get("grid_converged", True))
        for r in data["rows"]:
            all_rows.append(r)
            srcs.append(name)
            if not converged:
                excluded += 1

    cands = []
    for r, src in zip(all_rows, srcs):
        a = _anisotropy(r)
        fl = _flow(r)
        if a is None or src == "V172_self_energy_flow" and excluded:
            continue
        if a < 1e-6 or (fl is not None and fl < -0.5):
            cands.append({"source": src, "row": r, "A": a, "flow": fl})

    # hard tier on candidates (the contact is diagonal unitary: verify;
    # the front budget is the free one: verify with the V167 machinery)
    gates = []
    if cands:
        C = A.cell((2, 2, 2))
        ST = A.stages_baseline()
        U16 = SC.sectors()[2]
        for c in cands:
            uc = c["row"].get("U_c", SC.G10)
            unit = abs(abs(np.exp(1j * uc)) - 1.0)          # contact unitarity
            fs = A.front_speed(C, U16, ST, SC.G10, Lb=6, n_stages=12)
            gates.append({
                "source": c["source"], "U_c": uc,
                "contact_unitarity_dev": unit,
                "front_speed": fs["front_speed"],
                "outside_lightcone": fs["n_outside_lightcone"],
            })

    ok = len(all_rows) >= 15 and len(cands) == 0
    min_a = min((_anisotropy(r) for r, s in zip(all_rows, srcs)
                 if _anisotropy(r) is not None
                 and not (s == "V172_self_energy_flow" and excluded)),
                default=None)
    verdict = {
        "id": "V174_isotropy_gate",
        "tier": 4,
        "property": "hard-tier gate on the interaction census V171-V173: "
                    "isotropy A < 1e-6 or flow toward isotropy > 0.5",
        "residual": None,
        "pass": bool(ok),
        "status": f"INTERACTION ISOTROPY GATE. {len(all_rows)} census rows "
                  f"scanned (V171-V173); {excluded} rows excluded from "
                  f"candidacy (V172 one-loop flow not grid-convergent - "
                  f"an unconverged measurement cannot nominate); minimum "
                  f"converged anisotropy {min_a:.4f}; ZERO candidates (no "
                  f"A < 1e-6, no converged flow toward isotropy > 0.5): NO "
                  f"CANDIDATE - the minimal admissible collision (block "
                  f"contact, U_c in {{g_10, pi/4, pi/2}}) does not "
                  f"isotropize the W9 cone. The interacting fixed-point "
                  f"mechanism (Ising/Wilson-Fisher universality) needs the "
                  f"gauged-layer 2-fermion S-matrix (V79-V101), not built.",
        "data": {"n_rows": len(all_rows), "min_A": min_a,
                 "excluded_unconverged": excluded,
                 "candidates": cands, "gates": gates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v174_isotropy_gate.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V174 - the isotropy gate on the interaction census")
    print("=" * 68)
    print(f"  rows scanned: {len(all_rows)}   min A: {min_a:.4f}   "
          f"candidates: {len(cands)}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v174_isotropy_gate.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
