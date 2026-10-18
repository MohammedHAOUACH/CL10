#!/usr/bin/env python
"""V189 - does the merged library reproduce the search's own objective?

V189's verdict is the 2026-10-18 00:40 GPU run of the joint cell x step
search; the search itself costs four hours of GPU and is not re-executed.
What can be re-executed cheaply is its OBJECTIVE: every row of the verdict
stores the (theta*, q*) the descent reported, so the anisotropy and the node
frequency can be re-evaluated at those exact points under the current library.
A match means the carried numbers are the numbers the current library
produces; the descent is deterministic given the objective.

This is the check, not the search.  It takes seconds per cell on the card and
minutes on the CPU core (the (3,3,3) cell is M = 432).

    ./run_reproduce_check.sh                  # the (3,3,3) cell, on the card
    ./run_reproduce_check.sh --all            # the four cells
    CL10_JAX_CPU=1 ./run_reproduce_check.sh   # the CPU core instead

Verdict: ok = every cell's aniso deviation < 1e-3 (the anchor gate V189
itself uses) and omega(q*) still exactly 0.
"""
import importlib.util
import json
import os
import sys

os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
if os.environ.get("CL10_JAX_CPU", "0") == "1":
    # JAX_PLATFORMS, not the retired JAX_PLATFORM_NAME: under the latter JAX
    # still initializes the CUDA plugin, and the engine's pinned card index
    # then reaches for a device id the CPU-only platform does not have.
    os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)
_DEV = jax.devices()[0]
if _DEV.platform != "gpu" and os.environ.get("CL10_JAX_CPU", "0") != "1":
    sys.exit(f"V189-check: JAX resolved to {_DEV} - the card is unavailable. "
             f"Free it, or set CL10_JAX_CPU=1 to check on the CPU core.")

from cl10 import altcell as A  # noqa: E402
from cl10 import supercell as SC  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "V184", os.path.join(ROOT, "experiments", "alt_step",
                         "V184_step_loss_landscape",
                         "V184_step_loss_landscape.py"))
V184 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V184)

VERDICT = os.path.join(HERE, "results", "v189_cell_step_joint_search.json")
TOL = 1e-3  # the anchor gate V189 uses against V188's baseline

CELLS = [((2, 2, 2), "baseline", A.stages_baseline()),
         ((2, 2, 2), "phase-major", A.stages_shape((2, 2, 2))),
         ((2, 2, 3), "phase-major", A.stages_shape((2, 2, 3))),
         ((3, 3, 3), "phase-major", A.stages_shape((3, 3, 3)))]


def main():
    argv = sys.argv[1:]
    cells = CELLS if "--all" in argv else [c for c in CELLS if c[0] == (3, 3, 3)]
    d = json.load(open(VERDICT))
    rows = d["data"]["rows"]

    print("=" * 68)
    print(f"V189 objective reproduction - merged library, "
          f"JAX {jax.__version__} on {_DEV}, G10 = {SC.G10}")
    print("=" * 68)
    worst, ok = 0.0, True
    for shape, stagname, stages in cells:
        sel = [r for r in rows
               if r["cell"] == str(shape) and r["stages"] == stagname]
        best = min(sel, key=lambda r: r["A"])
        C = V184.build_core(shape, stages)
        th, q = np.array(best["theta_star"]), np.array(best["q_star"])
        a, om = float(C["aniso"](q, th)), float(C["omega"](q, th))
        dev = abs(a - best["A"])
        worst = max(worst, dev)
        good = dev < TOL and abs(om) < 1e-12
        ok = ok and good
        print(f"{str(shape):>9} {stagname:>12} {best['seed']:>11}  M={C['c']['nsub']*16:<5}"
              f" A {best['A']:.12f} -> {a:.12f}  dev {dev:.2e}"
              f"  omega {best['omega_node']:.1e} -> {om:.1e}"
              f"  gap2 {float(C['gap2'](q, th)):.1e}  [{'ok' if good else 'FAIL'}]")
    print("-" * 68)
    print(f"REPRODUCE: {'ok' if ok else 'FAIL'} - worst aniso deviation "
          f"{worst:.2e} against the {TOL:.0e} anchor gate; the candidacy "
          f"threshold is 1e-6.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
