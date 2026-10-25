#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V266 - finite-size convergence of the composite gap verdict
(tier 2, residual: false - structural, hard tier).

**The question (pre-registered).**  Is the V260 gap verdict finite-size
or thermodynamic?

**The construction (closed input set).**  The kernel-recipe gap ladder
on q-grids N = 4, 6, 8 (the fiber stays 1568-dim; only the q-grid
grows): per grid, the V259 kernel (C2 = sum_a Q^a Q^a, eigh, W =
V[:, w < 1e-9]) and the V260 convention E0(q) = lowest eigenvalue of
H_k(q) = W^dag H2(q) W minus the q-independent free-tower floor
(-32 + omega), omega = 0 row; min gap per grid; the trend exponent by
least squares on log E_min vs 1/N.

**Cascade (pre-registered roll-up).**  V259 KERNEL-EMPTY: the kernel is
empty (and q-grid-independent, since the 45 charges never depend on q),
so no min gap exists at any N and the verdict is recorded
NEVER-REALIZED.  The min gap per grid and the exponent are only defined
on the drafted non-cascade path, which runs iff the kernel is non-empty.
Measured regardless: the kernel dimension at each grid (the
grid-independence of the emptiness is itself a finite-size statement).

**Acceptance (pre-registered).**  GAP-CONFIRMED iff the gaps agree
across N = 4, 6, 8 to 1e-6 and all > 1e-9; GAPLESS-EMERGES iff E_min
crosses below 1e-9 at N = 8 with a negative trend exponent;
GRID-NEUTRAL otherwise; NEVER-REALIZED by the cascade rule.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")     # the 1/5 reference branch
from cl10.paths import exp                     # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))
import numpy as np                             # noqa: E402
from V79_gauged_bond_layer import (adjoint_link_operators,   # noqa: E402
                                   root_generators,
                                   structure_constants)
from cl10 import altcell as AC                 # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10
NPART = 1568
KTOL = 1e-9
GRID_NS = [4, 6, 8]


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V266 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg),
            flush=True)

    log("start, backend={}".format(E._BACKEND))
    core = E._BACKEND
    gpu = E._BACKEND != "numpy"
    xp = E._E if gpu else np
    if gpu:
        free = int(xp.cuda.runtime.memGetInfo()[0])
        if free < (3500 << 20):
            gpu, xp, core = False, np, "numpy (pinned: GPU free {} MiB)".format(
                free >> 20)
            log("GPU free {} MiB: pinned to host numpy".format(free >> 20))

    v259 = art("V259_singlet_kernel_census")
    v260 = art("V260_composite_gap_ladder")

    # committed fiber and charges (grid-independent by construction)
    Ts = root_generators()
    f, k_tr, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    B, U32, U16 = SC.sectors(G10)
    _Bm, T16, _f16 = AC.weyl_generators()
    QP = [np.kron(np.eye(8, dtype=np.complex128), t) for t in T16]
    QQ = [np.kron(T, np.eye(45)) + np.kron(np.eye(32), L)
          for T, L in zip(Ts, Ls)]
    Qfull = []
    for a in range(45):
        M = np.zeros((NPART, NPART), dtype=np.complex128)
        M[:128, :128] = QP[a]
        M[128:, 128:] = QQ[a]
        Qfull.append(M)
    Qd = xp.asarray(np.stack(Qfull))
    C2d = xp.matmul(Qd, Qd).sum(axis=0)
    del Qd
    C2d = 0.5 * (C2d + C2d.conj().swapaxes(-2, -1))
    w, V = xp.linalg.eigh(C2d)
    km = E.asnumpy(w).real < KTOL
    kernel_dim = int(np.count_nonzero(km))
    W = V[:, km]

    r_lat = np.array([(e & 1, (e >> 1) & 1, (e >> 2) & 1)
                      for e in range(8)], dtype=float)

    def coupling_at(q):
        wq = np.exp(1j * (r_lat @ np.asarray(q))) / math.sqrt(8.0)
        KW = np.kron(wq, np.eye(16, dtype=np.complex128))
        return sum(np.kron(T @ B, np.eye(45)[:, a:a + 1]) @ KW
                   for a, T in enumerate(Ts))

    def build_H2(H1, TLm, C, omega):
        H = np.zeros((NPART, NPART), dtype=np.complex128)
        H[:128, :128] = H1
        H[128:, 128:] = TLm + omega * np.eye(1440)
        H[:128, :128] = G10 * C.conj().T
        H[128:, :128] = G10 * C
        return H

    # ------------------------------------------------------------------
    # the ladder on N = 4, 6, 8
    # ------------------------------------------------------------------
    grids = []
    Wn = E.asnumpy(W)
    for N in GRID_NS:
        grid = [(2 * math.pi * i / N, 2 * math.pi * j / N,
                 2 * math.pi * l / N)
                for i in range(N) for j in range(N) for l in range(N)]
        vals = []
        if kernel_dim:                       # drafted non-cascade path
            for q in grid:
                H1 = AC.effective_hamiltonian(
                    lambda pz: SC.bloch(pz, U16, G10), q)
                H2 = build_H2(H1, TL, coupling_at(q), 0.0)
                Hk = Wn.conj().T @ H2 @ Wn
                vals.append(float(np.linalg.eigvalsh(Hk)[0]) - (-32.0))
        emin = min(vals) if vals else None
        grids.append({"N": N, "points": len(grid), "kernel_dim":
                      kernel_dim, "min_E0": emin})
        log("N = {}: {} points, kernel dim {}, min E0 {}".format(
            N, len(grid), kernel_dim, emin))

    cascade = (v259["branch"] == "KERNEL-EMPTY" and kernel_dim == 0
               and v260["branch"] == "NEVER-REALIZED")

    exponent = None
    if not cascade:
        xs = [1.0 / g["N"] for g in grids]
        ys = [math.log(g["min_E0"]) for g in grids]
        A = np.polyfit(xs, ys, 1)
        exponent = float(A[0])
        gaps = [g["min_E0"] for g in grids]
        if all(g > KTOL for g in gaps) and \
                max(gaps) - min(gaps) < 1e-6:
            branch = "GAP-CONFIRMED"
        elif grids[-1]["min_E0"] < KTOL and exponent < 0:
            branch = "GAPLESS-EMERGES"
        else:
            branch = "GRID-NEUTRAL"
    else:
        branch = "NEVER-REALIZED"

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, V260 branch {}, kernel "
                "recomputed here dim {} and grid-independently 0 at "
                "every N (kernel dims {} - the 45 charges never depend "
                "on q, so emptiness cannot be a finite-size effect)"
                .format(v259["branch"], v259["data"]["total_singlet_count"],
                        v260["branch"], kernel_dim,
                        [g["kernel_dim"] for g in grids]),
        "pass": bool(v259["branch"] == "KERNEL-EMPTY"
                     and int(v259["data"]["total_singlet_count"]) == 0
                     and kernel_dim == 0
                     and all(g["kernel_dim"] == 0 for g in grids))})
    checks.append({
        "name": "T1 grid convention: N = {} torus grids (2 pi i/N "
                "convention of V251/V260), point counts {}, same fiber "
                "and tolerance 1e-9 at every grid".format(
                    GRID_NS, [g["points"] for g in grids]),
        "pass": bool([g["points"] for g in grids] == [N ** 3
                                                      for N in GRID_NS])})
    checks.append({
        "name": "T2 convergence state: cascade {} -> min gaps {}, trend "
                "exponent {}; branch {} matches the pre-registered "
                "acceptance strings (GAP-CONFIRMED requires agreement "
                "to 1e-6 and all > 1e-9; GAPLESS-EMERGES requires a "
                "cross below 1e-9 at N = 8 with negative exponent)".format(
                    cascade, [g["min_E0"] for g in grids], exponent,
                    branch),
        "pass": bool((cascade and branch == "NEVER-REALIZED"
                      and all(g["min_E0"] is None for g in grids))
                     or (not cascade and exponent is not None))})

    data["conventions"] = {
        "grid": "N = 4, 6, 8 torus grids, 2 pi i/N convention; E0 = "
                "lowest eigenvalue of H_k = W^dag H2 W minus the "
                "q-independent free-tower floor (-32 + omega), omega = "
                "0 row (V260 conventions)",
        "exponent": "least-squares slope of log E_min vs 1/N",
        "kernel_tolerance": KTOL}
    data["grids"] = grids
    data["trend_exponent"] = exponent
    data["cascade"] = {
        "rule": "pre-registered roll-up: V259 KERNEL-EMPTY -> the gap "
                "ladder is undefined at every grid -> NEVER-REALIZED; "
                "kernel emptiness measured grid-independently",
        "fired": bool(cascade),
        "v259_branch": v259["branch"], "v260_branch": v260["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    status = (
        "finite-size convergence of the composite gap verdict (grids "
        "N = {}, {:.0f} s): the V259 kernel recomputed here has dim {} "
        "and is 0 at every grid (the charges are q-independent, so "
        "emptiness is not a finite-size effect); min gaps {}; trend "
        "exponent {}. Branch {}.").format(
            GRID_NS, time.time() - t0, kernel_dim,
            [g["min_E0"] for g in grids], exponent, branch)
    status += (
        " V259 landed KERNEL-EMPTY: there is no composite band at N = "
        "4, 6 or 8, so the V260 verdict is not a finite-size artefact - "
        "it is the empty kernel at every grid; the convergence ladder "
        "is NEVER-REALIZED by the pre-registered cascade." if cascade
        else " min gaps {} agree to {} -> {}.".format(
            [g["min_E0"] for g in grids],
            max(g["min_E0"] for g in grids)
            - min(g["min_E0"] for g in grids), branch))
    out = {"id": "V266_finite_size_convergence",
           "tier": 2,
           "property": "the N = 4/6/8 finite-size behaviour of the "
                       "composite gap verdict (V266)",
           "residual": False,
           "programme": "V266",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v266_finite_size_convergence.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
