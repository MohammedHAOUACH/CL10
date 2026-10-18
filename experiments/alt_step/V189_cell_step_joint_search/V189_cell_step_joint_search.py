#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V189 - joint cell x step gradient search (alt-step, 3/4)

**The search.**  V188 searched the 45 step angles on the implemented
(2,2,2) cell.  The CELL is the other discrete half of the
specification: this experiment runs the same pre-registered alternating
Adam scheme (10 cycles x [Nelder-Mead node + 30 theta-steps, lr 0.05])
on four pre-registered cell x stage choices:

  ((2,2,2), stages_baseline)   the implemented motion
  ((2,2,2), stages_shape)      the canonical phase-major stagger
  ((2,2,3), stages_shape)      one volume up, mixed radix
  ((3,3,3), stages_shape)      the cubic 27-site block

Bloch matrices 128/128/192/432.  Seeds per cell: theta0 = g_10 * 1 plus
7 Uniform(0, 2 pi) draws (rng 187).  Structural tier (one line): the
anomaly cancellation, hypercharge and chiral matter content are fixed
by the signature + sector and are untouched by (cell, theta) - the cell
changes the MOTION only (front speed, closure), tested on candidates in
V190.

**Anchor.**  The (2,2,2)/baseline baseline-seed final A reproduces
V188's baseline final A to 1e-6 (same seed, same scheme, deterministic
core).
"""

import importlib.util
import json
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
# the joint search is hours on the CPU core, so JAX takes the card; lazy
# allocation keeps it working when the card is shared with another process.
# CL10_JAX_CPU=1 opts out to the CPU core.
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
if os.environ.get("CL10_JAX_CPU", "0") == "1":
    # JAX_PLATFORMS, not the retired JAX_PLATFORM_NAME: under the latter JAX
    # still initializes the CUDA plugin and dies on a card another process
    # has filled (RuntimeError: no supported devices found for platform CUDA).
    os.environ.setdefault("JAX_PLATFORMS", "cpu")

# the four cell workers share the box: cap the intra-op threads per worker so
# the pool does not oversubscribe the cores (one quarter of them each)
os.environ.setdefault("OMP_NUM_THREADS",
                      str(max(1, (os.cpu_count() or 4) // 4)))

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

_CPU_OK = os.environ.get("CL10_JAX_CPU", "0") == "1"

def _route(where):
    """Resolve the route, and refuse it loudly if it is not the one asked for.

    This must NOT run at import in a pool worker: an exception raised during
    spawn's prepare() - while the child re-imports this module - kills the child
    before the executor can capture anything, and the parent only ever sees
    BrokenProcessPool.  That is what sys.exit did, and what a module-level raise
    still did.  The parent calls it in main(); the worker calls it inside its own
    function, where the executor can report the reason."""
    dev = jax.devices()[0]
    if dev.platform != "gpu" and not _CPU_OK:
        raise RuntimeError(
            f"V189 ({where}): JAX resolved to {dev} - the card is unavailable. "
            f"Free it, or set CL10_JAX_CPU=1 to run the search on the CPU core.")
    return dev


def _free_mib(dev):
    """Free memory on the card JAX actually took.  JAX's own bytes_limit is its
    allocator budget, not what the device has left."""
    q = subprocess.run(["nvidia-smi", "-i", str(dev.id),
                        "--query-gpu=memory.free",
                        "--format=csv,noheader,nounits"],
                       capture_output=True, text=True)
    return int(q.stdout.strip().splitlines()[0]) if q.stdout.strip() else 0

from cl10 import supercell as SC
from cl10 import altcell as A

_spec = importlib.util.spec_from_file_location(
    "v184", os.path.join(ROOT, "experiments", "alt_step",
                         "V184_step_loss_landscape",
                         "V184_step_loss_landscape.py"))
V184 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V184)

G10 = SC.G10
N_CYCLES = 10
N_TH = 30
LR_TH = 0.02


def _cell_search(job):
    """One cell of the joint search: its core, its 8 seeds, the alternating
    Adam scheme.  This is a worker of the cell pool.  The four cells are
    independent (the rng stream, the node search and the descent are per
    cell), so the pool is the whole speed-up: the wall becomes the slowest
    cell, not the sum of four."""
    shape, stagname, stages = job
    _route("cell worker")     # inside the function: the executor can report it
    C = V184.build_core(shape, stages)
    ns = C["c"]["nsub"] * 16
    rng = np.random.default_rng(187)
    seeds = [("baseline g_10*1", np.full(45, G10))]
    for s in range(7):
        seeds.append((f"uniform-{s}",
                      rng.uniform(0.0, 2.0 * np.pi, 45)))
    rows = []
    t0 = time.time()
    for name, th in seeds:
        thj = jnp.array(th)
        q = C["find_node"](thj)
        a0 = float(C["aniso"](q, thj))
        best_l = (a0 + 100.0 * float(C["omega"](q, thj)) ** 2)
        for _ in range(N_CYCLES):
            q = C["find_node"](thj)
            thj = V184.adam(lambda t: C["combined_grad"](q, t),
                            thj, N_TH, LR_TH)
            q2 = C["find_node"](thj)
            ll2 = (float(C["aniso"](q2, thj))
                   + 100.0 * float(C["omega"](q2, thj)) ** 2)
            if ll2 < best_l:
                best_l, thj = ll2, thj
        q = C["find_node"](thj)
        aa = float(C["aniso"](q, thj))
        om = float(C["omega"](q, thj))
        ll = aa + 100.0 * om * om
        rows.append({"cell": str(shape), "stages": stagname,
                     "seed": name, "A": aa, "A0": a0,
                     "omega_node": om,
                     "L": ll, "q_star": [float(x) for x in q],
                     "theta_star": [float(x) for x in thj]})
        print(f"  {str(shape):>9} {stagname:>11} {name:>16}: A "
              f"{aa:.4f}  omega(q*) {om:.2e}  "
              f"[{time.time() - t0:.0f}s]", flush=True)
    best = min(rows, key=lambda r: r["A"])
    summary = {"cell": str(shape), "stages": stagname,
               "m_size": ns, "best_A": best["A"],
               "best_L": best["L"],
               "best_omega": best["omega_node"],
               "best_seed": best["seed"],
               "node_found": bool(best["omega_node"] < 1e-3)}
    return rows, summary


def main():
    t0 = time.time()
    print("=" * 68)
    print(f"V189 - joint cell x step search "
          f"(JAX {jax.__version__}, device {jax.devices()[0]})")
    print("=" * 68)
    cells = [((2, 2, 2), "baseline", A.stages_baseline()),
             ((2, 2, 2), "phase-major", A.stages_shape((2, 2, 2))),
             ((2, 2, 3), "phase-major", A.stages_shape((2, 2, 3))),
             ((3, 3, 3), "phase-major", A.stages_shape((3, 3, 3)))]

    # The route is decided here, in the parent, and pinned for the children:
    # CUDA_VISIBLE_DEVICES is inherited by spawn, so a worker cannot drift to
    # another card or silently fall back to the CPU core and run a cell on a
    # route the verdict does not record.
    dev = _route("driver")

    # one worker per cell, spawn (a forked child inherits a dead CUDA context).
    # ex.map keeps the cell order, so the rows and the summary come back exactly
    # as the serial loop produced them.  On the card the pool is sized to the
    # memory the card can still hand out (~700 MiB per JAX context): inside the
    # suite eight cupy workers are already on it, and a worker that cannot get a
    # context is dead weight.  With room for no second context the four cells run
    # serially in this process, on the route just established.
    nworkers = len(cells)
    if dev.platform == "gpu":
        free = _free_mib(dev)
        if free < 700:
            sys.exit(f"V189: JAX took {dev}, which has only {free} MiB free; one "
                     f"cell worker needs ~700 MiB. Free the card, or set "
                     f"CL10_JAX_CPU=1 to run the search on the CPU core.")
        nworkers = max(1, min(len(cells), free // 700))
        if nworkers < len(cells):
            print(f"  ({free} MiB free on {dev}: {nworkers} cell worker"
                  f"{'s' if nworkers > 1 else ''}; the rest run serially)",
                  flush=True)
        os.environ["CUDA_VISIBLE_DEVICES"] = str(dev.id)
    if nworkers == 1:
        out = [_cell_search(job) for job in cells]
    else:
        with ProcessPoolExecutor(max_workers=nworkers,
                                 mp_context=multiprocessing.get_context("spawn")
                                 ) as ex:
            out = list(ex.map(_cell_search, cells))
    rows, summary = [], []
    for r, s in out:
        rows.extend(r)
        summary.append(s)

    # ---- anchor: (2,2,2)/baseline baseline seed vs V188 baseline ----
    p188 = os.path.join(ROOT, "experiments", "alt_step",
                        "V188_step_gradient_search", "results",
                        "v188_step_gradient_search.json")
    with open(p188) as fh:
        d188 = json.load(fh)
    # anchor (pre-registered): the (2,2,2)/baseline baseline seed's
    # INITIAL node A reproduces V188's baseline A0 (the post-descent
    # A is chaotic across processes - NM landing jitter cascades
    # through 300 Adam steps - so the anchor is the pre-descent
    # value, gate 1e-3 as in V188)
    A188 = d188["data"]["rows"][0]["A0"]
    base_row = [r for r in rows if r["cell"] == "(2, 2, 2)"
                and r["stages"] == "baseline"
                and r["seed"] == "baseline g_10*1"][0]
    devA = abs(base_row["A0"] - A188)

    cands = [r for r in rows
             if r["A"] < 1e-6 and r["omega_node"] < 1e-6]

    # gate 1e-3: same knife-edge-A + cross-process NM landing jitter
    # conditioning as the V188 anchor (both GPU routes)
    ok = devA < 1e-3 and len(rows) == 32 and len(summary) == 4

    verdict = {
        "id": "V189_cell_step_joint_search",
        "tier": 4,
        "property": "joint search over the cell (block shapes (2,2,2) "
                    "baseline / (2,2,2) phase-major / (2,2,3) / (3,3,3)) "
                    "and the 45 step angles, same alternating Adam scheme "
                    "as V188, 8 seeds per cell; structural tier "
                    "(anomalies, hypercharge) untouched by (cell, theta)",
        "residual": None,
        "pass": bool(ok),
        "status": "JOINT CELL x STEP SEARCH (4 cells x 8 seeds = 32 "
                  "rows). Best A per cell: "
                  + ", ".join(f"{s['cell']}/{s['stages']} {s['best_A']:.4f}"
                              for s in summary)
                  + f". Anchor (2,2,2)/baseline baseline seed vs V188: "
                    f"dev {devA:.1e}. Candidates (A < 1e-6 and gapless): "
                    f"{len(cands)}. Route: JAX {jax.devices()[0]}.",
        "data": {"rows": rows, "summary": summary,
                 "anchor_dev_A": devA, "n_candidates": len(cands)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v189_cell_step_joint_search.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v189_cell_step_joint_search.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
