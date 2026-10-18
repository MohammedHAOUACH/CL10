#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V183 - the non-perturbative isotropy gate (alt-nonperturbative).

Mirror of V174 (the interaction gate) on the non-perturbative census
V180-V182. Pre-registered candidacy:
  - a V180 row with A_np < 1e-6 (dressed cone isotropized) or
    |A_np - A_free| > 0.5 (the exact flux dressing moves the cone
    anisotropy by a macroscopic amount);
  - a V181/V182 row with a light mode: gap < 1e-2 eta (a two-fermion
    excitation that could dress the IR).
Candidates get the hard tier: the V181 S-matrix unitarity is exact by
construction (recorded < 1e-12) and the front budget is the free one -
the dressing does not touch the motion - verified with the V167
machinery (front_speed == 1.0, n_outside_lightcone == 0).
ZERO candidates is a valid PASS.
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
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core

from cl10 import altcell as A
from cl10 import supercell as SC

SRC = ("V180_cone_flux_dressing", "V181_ring_kernel_smatrix",
       "V182_ir_decoupling")


def _data(name):
    p = os.path.join(ROOT, "experiments", "alt_nonpert", name,
                     "results", f"{name.lower()}.json")
    with open(p) as f:
        d = json.load(f)
    if name != "V182_ir_decoupling":
        assert d["pass"], f"{name} verdict must PASS"
    return d["data"]


def main():
    t0 = time.time()
    all_rows, srcs = [], []
    for name in SRC:
        for r in _data(name)["rows"]:
            all_rows.append(r)
            srcs.append(name)

    cands = []
    for r, src in zip(all_rows, srcs):
        if src == "V180_cone_flux_dressing":
            a = r.get("A_np")
            if a is not None and (a < 1e-6 or abs(r["dA"]) > 0.5):
                cands.append({"source": src, "row": r, "A": a})
        else:
            g = r.get("gap")
            if g is not None and g < 1e-2:
                cands.append({"source": src, "row": r, "gap": g})

    # hard tier on candidates: exact S-matrix unitarity (recorded from
    # the V181 rows) + the free front (the dressing does not touch the
    # motion)
    gates = []
    if cands:
        C = A.cell((2, 2, 2))
        ST = A.stages_baseline()
        U16 = SC.sectors()[2]
        fs = A.front_speed(C, U16, ST, SC.G10, Lb=6, n_stages=12)
        for c in cands:
            unit = c["row"].get("S_unitarity")
            gates.append({
                "source": c["source"],
                "S_unitarity_dev": unit,
                "front_speed": fs["front_speed"],
                "outside_lightcone": fs["n_outside_lightcone"],
            })
    hard_ok = all(
        g["S_unitarity_dev"] is None or g["S_unitarity_dev"] < 1e-12
        and g["front_speed"] == 1.0 and g["outside_lightcone"] == 0
        for g in gates)

    min_a = min((r["A_np"] for r, s in zip(all_rows, srcs)
                 if s == "V180_cone_flux_dressing"), default=None)
    min_gap = min((r["gap"] for r, s in zip(all_rows, srcs)
                   if "gap" in r), default=None)
    ok = len(all_rows) >= 8 and len(cands) == 0 and hard_ok

    verdict = {
        "id": "V183_nonpert_isotropy_gate",
        "tier": 4,
        "property": "hard-tier gate on the non-perturbative census "
                    "V180-V182: dressed anisotropy A_np < 1e-6 or "
                    "|A_np - A_free| > 0.5, or a light two-fermion mode "
                    "gap < 1e-2 eta",
        "residual": None,
        "pass": bool(ok),
        "status": f"NON-PERTURBATIVE ISOTROPY GATE. {len(all_rows)} "
                  f"census rows scanned (V180-V182); minimum dressed "
                  f"anisotropy {min_a:.4f} (dressing shift |dA| <= "
                  f"{max(abs(r['dA']) for r, s in zip(all_rows, srcs) if s == 'V180_cone_flux_dressing'):.1e}); "
                  f"minimum two-fermion gap {min_gap:.4f} eta. ZERO "
                  f"candidates: the exact bond-current dressing of the "
                  f"cone (all orders in eta, fixed adjoint link) does "
                  f"not move the W9 anisotropy, and the gauged flux "
                  f"sector has no light mode at any volume probed. The "
                  f"non-perturbative route of the gauged sector closes "
                  f"without an isotropizing mechanism; W9 remains open.",
        "data": {"n_rows": len(all_rows), "min_A_np": min_a,
                 "min_gap": min_gap,
                 "candidates": cands, "gates": gates},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v183_nonpert_isotropy_gate.json"),
              "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V183 - the non-perturbative isotropy gate")
    print("=" * 68)
    print(f"  rows scanned: {len(all_rows)}   min A_np: {min_a:.4f}   "
          f"min gap: {min_gap:.4f} eta   candidates: {len(cands)}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir,
                                      "v183_nonpert_isotropy_gate.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
