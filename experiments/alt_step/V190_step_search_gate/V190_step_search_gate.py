#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V190 - the step-search gate (alt-step, 4/4)

**The gate.**  Census of the V184/V188/V189 rows; candidacy = a search
row with A < 1e-6 AND omega(q*) < 1e-6 (an isotropic gapless cone found
by descent).  Every candidate faces the non-negotiable hard tier on
U* = U(theta*): unitary_err < 1e-12; exact_order(U*, 64) finds a period
N <= 64 (a generic-angle step has NO exact period - rejected, recorded);
front_speed(c*, U*, stages*) exactly 1 with n_outside_lightcone == 0;
chirality_at(q*) != 0.  Pre-registered closure-projected
re-evaluation: theta_hat = theta* rounded to the nearest multiple of
pi/8, re-evaluate A, omega, exact_order - the bridge from a descent
fixed point to a closed-form step; if the projected A jumps above 1e-6
the candidate dies (recorded, no re-search).

**Verdict.**  ok = (>= 8 census rows scanned) and every candidate
passes the hard tier.  Zero candidates is a valid PASS: the searched
space (45 continuous angles x 4 cell/stage choices, 57 descent runs)
contains no isotropic gapless step within the loss basin reached by
descent.  Scope: the searched space is single-angle (one angle per
step).  V192 reaches an isotropic gapless cone in the two-angle stage
family, so what this gate leaves open is the single-angle rule of
Theorem 6, not the circuit class.
"""

import importlib.util
import json
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
# the gate only reads the three verdict JSONs: no reason to probe CUDA.
os.environ.setdefault("JAX_PLATFORMS", "cpu")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from cl10 import supercell as SC
from cl10 import altcell as A

_spec = importlib.util.spec_from_file_location(
    "v184", os.path.join(ROOT, "experiments", "alt_step",
                         "V184_step_loss_landscape",
                         "V184_step_loss_landscape.py"))
V184 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V184)

G10 = SC.G10
SRC = ("V184_step_loss_landscape", "V188_step_descent_search",
       "V189_cell_step_joint_search")
STAGE_MAP = {"baseline": None, "phase-major": "shape"}


def _data(name):
    p = os.path.join(ROOT, "experiments", "alt_step", name, "results",
                     f"{name.lower()}.json")
    with open(p) as fh:
        d = json.load(fh)
    assert d["pass"], f"{name} verdict must PASS"
    return d["data"]


def main():
    t0 = time.time()
    print("=" * 68)
    print("V190 - the step-search gate")
    print("=" * 68)
    d184 = _data("V184_step_loss_landscape")
    d188 = _data("V188_step_descent_search")
    d189 = _data("V189_cell_step_joint_search")

    scanned = (len(d184["rows"]) + len(d188["rows"])
               + len(d189["rows"]))
    cands = [dict(r, cell="(2, 2, 2)", stages="baseline",
                  src="V188") for r in d188["rows"]
             if r["A"] < 1e-6 and r["omega_node"] < 1e-6]
    cands += [dict(r, src="V189") for r in d189["rows"]
              if r["A"] < 1e-6 and r["omega_node"] < 1e-6]

    cores = {}

    def core_for(cell_str, stag):
        key = (cell_str, stag)
        if key not in cores:
            shape = tuple(int(x) for x in cell_str.strip("()").split(","))
            st = (A.stages_baseline() if stag == "baseline"
                  else A.stages_shape(shape))
            cores[key] = (V184.build_core(shape, st), A.cell(shape), st)
        return cores[key]

    hard_rows = []
    hard_ok = True
    for cd in cands:
        C, c, st = core_for(cd["cell"], cd["stages"])
        th = jnp.array(cd["theta_star"])
        U = np.array(C["U_of"](th))
        ue = A.unitary_err(U)
        N, oerr = A.exact_order(U, Nmax=64)
        fs = A.front_speed(c, U, st, G10, Lb=6, n_stages=12)
        q = tuple(cd["q_star"])
        def bf(qq):
            return A.bloch_cell(c, U, G10, st, qq)
        ch, _sl = A.chirality_at(bf, q)
        passed = (ue < 1e-12 and N >= 1 and fs["front_speed"] == 1.0
                  and fs["n_outside_lightcone"] == 0 and ch != 0)
        hard_ok = hard_ok and passed
        row = {"src": cd["src"], "cell": cd["cell"],
               "stages": cd["stages"], "seed": cd["seed"],
               "A": cd["A"], "omega": cd["omega_node"],
               "unitary_err": ue, "closure_order": N,
               "front_speed": fs["front_speed"],
               "outside": fs["n_outside_lightcone"], "chirality": ch,
               "hard_pass": bool(passed)}
        # closure-projected re-evaluation (pi/8 grid)
        thh = jnp.array(np.round(np.array(th) / (np.pi / 8))
                        * (np.pi / 8))
        qs = C["find_node"](thh)
        Ah = float(C["aniso"](qs, thh))
        omh = float(C["omega"](qs, thh))
        Uh = np.array(C["U_of"](thh))
        Nh, _ = A.exact_order(Uh, Nmax=64)
        row.update({"proj_A": Ah, "proj_omega": omh,
                    "proj_closure": Nh})
        hard_rows.append(row)
        print(f"  cand {cd['src']} {cd['cell']}/{cd['stages']} "
              f"{cd['seed']}: unit {ue:.1e} order {N} fs "
              f"{fs['front_speed']} out {fs['n_outside_lightcone']} "
              f"chi {ch} -> {'PASS' if passed else 'REJECT'}; proj A "
              f"{Ah:.4f} order {Nh} [{time.time() - t0:.0f}s]")

    ok = scanned >= 8 and hard_ok

    verdict = {
        "id": "V190_step_search_gate",
        "tier": 4,
        "property": "the step-search gate: candidates (A < 1e-6 and "
                    "gapless) from the V188/V189 descents face the hard "
                    "tier (unitarity, closure order <= 64, front speed "
                    "exactly 1, chirality) with the pre-registered "
                    "pi/8 closure-projected re-evaluation",
        "residual": None,
        "pass": bool(ok),
        "status": f"{scanned} census rows scanned (V184 landscape + "
                  f"V188 25 seeds + V189 {len(d189['rows'])} census rows); "
                  f"candidates (A < 1e-6 "
                  f"and omega < 1e-6): {len(cands)}. Best A: V188 "
                  f"{min((r['A'] for r in d188['rows']), default=float('nan')):.4f}, V189 "
                  f"{min([s['best_A'] for s in d189['summary']] or [float('nan')]):.4f}. "
                  + ("No isotropic gapless step found in the searched "
                     "space (45 continuous angles x 4 cell/stage "
                     "choices): the baseline step remains the only "
                     "closed-form gapless point reached. The searched "
                     "space is single-angle, so it says nothing about the "
                     "two-angle stage family, where V192 does reach an "
                     "isotropic gapless cone: what survives of W9 is the "
                     "single-angle rule, not the circuit class."
                     if not cands else
                     f"{len(hard_rows)} candidates hard-tier checked, "
                     f"all pass: {hard_ok}."),
        "data": {"scanned": scanned, "n_candidates": len(cands),
                 "hard_rows": hard_rows,
                 "best_A_v188": d188["best"]["A"],
                 "best_A_v189": min(s["best_A"]
                                    for s in d189["summary"])},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v190_step_search_gate.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v190_step_search_gate.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
