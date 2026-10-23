#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V36 - the automaton: the causal radius **implemented** (the 't Hooft
program, C1).

**The 0-parameter result.**  The step of CL10 is a **product of local
(2-cell) unitaries** on the 2x2x2 block lattice:

- the **internal layer** (the 0-parameter 45 Spin(10) plane rotations at
  ``g_10``) on each of the 8 sub-sites of the cell;
- the **bond stage** (the x, y, z directions **cycle** through the six stages
  of the step - one direction per stage, so a signal gets **one cell of budget
  per stage**; the 't Hooft staggered even / odd bonds, now read as the
  sub-lattice structure of the cell: the even bonds are inside the block, the
  odd bonds join neighbouring blocks): the **0-parameter swap-rotation**
  ``U_bond = cos(g_10) I - i sin(g_10) P`` (``P`` = the 0-parameter swap of the
  two cell states, a permutation - no angle).

The **only** number of the step is ``g_10`` (the 0-parameter GUT coupling).
Therefore:

- the step is **unitary** (the norm is preserved);
- a **localized** excitation stays **inside** the **light cone** (the
  Manhattan ball of radius k in Planck cells after k stages) - **causality is
  a property of the dynamics, not an axiom**;
- the **front** of the support **reaches** the light cone **exactly** (the
  **front speed is c = 1 cell per stage** - the light cone is **realized**);
- the light cone is a **structural** property (the **same** cone is realized
  for different 0-parameter angles - it is carried by the **2-cell operator
  structure**, not by the value of the angle).
"""

import sys
import os
import json
import time
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import engine as E
from cl10 import supercell as SC

TOL = 1e-9


def cone(Lb=8, n_stages=12, x0=(4, 4, 4), e0=0, theta=SC.G10, bonds=True,
         seed=0, tol=1e-12):
    """Run the automaton from a state localised on one sub-site, stage by
    stage; the support radius is measured in **Planck cells** (the cell of the
    automaton is the block: cell coordinate = 2R + e)."""
    rng = np.random.default_rng(seed)
    U = E.Array_(SC.sectors(theta)[1])
    f = E.zeros((Lb, Lb, Lb, SC.NSUB, SC.NINT))
    v = rng.standard_normal(SC.NINT) + 1j * rng.standard_normal(SC.NINT)
    f[x0 + (e0,)] = E.Array_(v / np.linalg.norm(v))

    # the cell coordinate of every (block, sub-site)
    R = np.arange(Lb)
    sub = np.arange(SC.NSUB)
    ex, ey, ez = sub & 1, (sub >> 1) & 1, (sub >> 2) & 1
    cx = 2 * R[:, None, None, None] + ex
    cy = 2 * R[None, :, None, None] + ey
    cz = 2 * R[None, None, :, None] + ez
    L = 2 * Lb

    def _md(d):
        return np.minimum(np.abs(d), L - np.abs(d))

    sx = 2 * x0[0] + (e0 & 1)
    sy = 2 * x0[1] + ((e0 >> 1) & 1)
    sz = 2 * x0[2] + ((e0 >> 2) & 1)
    dist = _md(cx - sx) + _md(cy - sy) + _md(cz - sz)

    steps = []
    for k in range(1, n_stages + 1):
        f = SC.stage(f, U, k - 1, theta) if bonds else SC.internal_layer(f, U)
        nrm = np.sqrt(np.sum(np.abs(E.asnumpy(f)) ** 2, axis=-1))
        support = np.argwhere(nrm > tol)
        d = [int(dist[tuple(p)]) for p in support]
        steps.append({"k": k, "n_support": int(len(support)),
                      "max_radius": max(d) if d else 0,
                      "n_outside": int(sum(1 for r in d if r > k)),
                      "outside_amp": float(nrm[dist > k].max())})
    return {"steps": steps,
            "final_norm": float(np.linalg.norm(E.asnumpy(f).ravel()))}


def main():
    t0 = time.time()
    k_struct = 9
    r_full = cone()
    r_int = cone(bonds=False)
    r_long = cone(Lb=32, x0=(16, 16, 16), n_stages=30)
    r_alt = cone(theta=0.3, n_stages=k_struct)

    unitary_err = abs(r_full["final_norm"] - 1.0)
    causal_ok = all(s["n_outside"] == 0 for s in r_full["steps"] + r_long["steps"])
    saturates_ok = all(s["max_radius"] == s["k"] for s in r_full["steps"])
    front_speed = (r_full["steps"][-1]["max_radius"]
                   / r_full["steps"][-1]["k"])
    structural_ok = all(
        r_alt["steps"][i]["n_outside"] == 0
        and r_alt["steps"][i]["max_radius"] == r_full["steps"][i]["max_radius"]
        == r_alt["steps"][i]["k"]
        for i in range(k_struct))
    nontrivial_ok = (r_full["steps"][-1]["n_support"] > 1
                     and r_int["steps"][-1]["n_support"] == 1)

    ok = (unitary_err < TOL and causal_ok and saturates_ok
          and abs(front_speed - 1.0) < 1e-9
          and structural_ok and nontrivial_ok)

    verdict = {
        "id": "V36_local_qca",
        "tier": 2,
        "property": "P2,P4",
        "residual": False,   # the automaton step is a 0-parameter PASS (no new residual)
        "pass": bool(ok),
        "status": "the_causal_radius_1_is_implemented_the_light_cone_is_realized",
        "tol": TOL,
        "data": {
            "backend": E._BACKEND,
            "g10": float(SC.G10),
            "cell": f"{SC.NSUB} sub-sites x {SC.NINT} internal",
            "stages_per_step": len(SC.STAGES),
            "blocks_per_step": 3,
            "unitary_err": float(unitary_err),
            "causal_ok": bool(causal_ok),
            "saturates_ok": bool(saturates_ok),
            "max_outside_amplitude": max(
                s["outside_amp"]
                for s in r_full["steps"] + r_long["steps"] + r_alt["steps"]),
            "front_speed": float(front_speed),
            "structural_ok": bool(structural_ok),
            "nontrivial_ok": bool(nontrivial_ok),
            "support_growth": [s["n_support"] for s in r_full["steps"]],
            "radius_growth": [s["max_radius"] for s in r_full["steps"]],
            "long_run_Lb": 32,
            "long_run_stages": 30,
            "long_run_radius": [s["max_radius"] for s in r_long["steps"]],
            "long_run_outside": sum(s["n_outside"] for s in r_long["steps"]),
            "long_run_unitary_err": abs(r_long["final_norm"] - 1.0),
            "internal_only_growth": [s["n_support"] for s in r_int["steps"]],
            "theta_tested": 0.3,
            "residual": "none new (the QCA lattice RG remains the named open "
                        "problem, as before)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v36_local_qca.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V36 - the automaton: the causal radius IMPLEMENTED (C1)")
    print("=" * 68)
    print(f"  backend (GPU):  {E._BACKEND}")
    print(f"  the cell:  {SC.NSUB} sub-sites x {SC.NINT} internal = {SC.NONSITE};  "
          f"g_10 = {SC.G10}  (the 45 rotations + the swap-rotations)")
    print(f"  the unitary (the norm is preserved):  residual = {unitary_err:.2e}")
    print(f"  the LIGHT CONE is REALIZED:  the front speed c = {front_speed:.3f}  "
          f"(1 Planck cell / stage, {len(SC.STAGES)} stages = 1 step)")
    print(f"    the support after k stages:  {verdict['data']['support_growth']}")
    print(f"    the radius after k stages (cells):  {verdict['data']['radius_growth']}")
    print(f"    long run (64-cell torus, 30 stages): radius "
          f"{verdict['data']['long_run_radius'][-1]} cells, outside "
          f"{verdict['data']['long_run_outside']}, largest amplitude outside "
          f"the cone {verdict['data']['max_outside_amplitude']:.2e}")
    print(f"    the internal layer alone (no 2-cell coupling):  "
          f"{verdict['data']['internal_only_growth']}")
    print(f"  the light cone is STRUCTURAL (independent of g_10, tested at 0.3):  "
          f"{structural_ok}")
    print(f"  VERDICT:  {'PASS (the causal radius 1 is a property of the dynamics, not an axiom)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v36_local_qca.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
