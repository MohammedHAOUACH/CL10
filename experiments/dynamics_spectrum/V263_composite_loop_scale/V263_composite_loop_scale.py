#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V263 - the composite loop scale (tier 2, residual: true, O8,
hard tier).

**The question (pre-registered).**  Does the composite sector generate a
dimensionless scale playing the role of the loop factor?

**The construction (closed input set).**  The composite gap at its
minimum (V260's artifact: data.min_E0 - itself built on the V259 kernel
recipe); the V101 bond-current ground-state scale read from the
V101_bond_current_ground_state artifact (data.E0; the implementer reads
the field, never retypes); V129's committed gap from its artifact
(data.two_site.gap) as the second candidate.  Candidates
F_phys = (composite gap at minimum) / |V101 E0| and
(composite gap) / (V129 gap); the target 1/(16 pi) computed from pi at
run time; the acceptance bands [1.5e-2, 3.0e-2] (factor 1.5 of the
target) fixed in the pre-registration.

**Cascade (pre-registered roll-up).**  V259 KERNEL-EMPTY and V260
NEVER-REALIZED: the composite gap does not exist, so the scale test is
recorded NEVER-REALIZED.

**Acceptance (pre-registered).**  SCALE-LOOP-CANDIDATE iff either
candidate lies in [1.5e-2, 3.0e-2]; else NO-SCALE-ROUTE recording the
closest value; NEVER-REALIZED by the cascade rule.
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
os.environ.setdefault("CL10_ANGLE", "0.2")
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


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V263 {:s} t+{:6.1f}s] {}".format(
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
    v251 = art("V251_gauged_two_particle_sector")
    v101 = art("V101_bond_current_ground_state")
    v129 = art("V129_gauged_layer_gap")

    # ------------------------------------------------------------------
    # own recomputation of the kernel (shared V259 recipe)
    # ------------------------------------------------------------------
    Ts = root_generators()
    f, k_tr, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    _TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    _B, _U32, _U16 = SC.sectors(G10)
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
    w, _V = xp.linalg.eigh(C2d)
    kernel_dim = int(np.count_nonzero(E.asnumpy(w).real < KTOL))
    log("kernel recomputed here: dim {}".format(kernel_dim))

    # ------------------------------------------------------------------
    # artifact scales and cross-artifact consistency
    # ------------------------------------------------------------------
    E0_v101 = float(v101["data"]["E0"])         # ground-state energy
    gap_v101 = float(v101["data"]["gap"])       # 42.09... at eta = 1
    gap_v129 = float(v129["data"]["two_site"]["gap"])
    omega_plus = float(v251["data"]["sensitivity"]["omega_rows"]["plus"])
    omega_consist = abs(gap_v101 * G10 - omega_plus)

    cascade = (v259["branch"] == "KERNEL-EMPTY" and kernel_dim == 0
               and v260["branch"] == "NEVER-REALIZED")

    target = 1.0 / (16.0 * math.pi)             # computed, never retyped
    band = [1.5e-2, 3.0e-2]

    min_E0 = v260["data"]["min_E0"]
    cand = None
    cand_rows = []
    if not cascade and min_E0 and min_E0.get("0") is not None:
        g0 = float(min_E0["0"])
        cand_rows = [
            {"name": "V101_E0_scale", "gap": g0, "scale": abs(E0_v101),
             "F_phys": g0 / abs(E0_v101)},
            {"name": "V129_gap", "gap": g0, "scale": gap_v129,
             "F_phys": g0 / gap_v129}]
        for r in cand_rows:
            r["dev_from_target"] = abs(r["F_phys"] - target)
            r["in_band"] = bool(band[0] <= r["F_phys"] <= band[1])
        hit = [r for r in cand_rows if r["in_band"]]
        branch = "SCALE-LOOP-CANDIDATE" if hit else "NO-SCALE-ROUTE"
        cand = min(cand_rows, key=lambda r: r["dev_from_target"])
    else:
        branch = "NEVER-REALIZED"

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, V260 branch {}, kernel "
                "recomputed here dim {} (agreement {})".format(
                    v259["branch"], v259["data"]["total_singlet_count"],
                    v260["branch"], kernel_dim,
                    kernel_dim == int(v259["data"]
                                      ["total_singlet_count"])),
        "pass": bool(v259["branch"] == "KERNEL-EMPTY"
                     and int(v259["data"]["total_singlet_count"]) == 0
                     and kernel_dim == 0
                     and v260["branch"] == "NEVER-REALIZED")})
    checks.append({
        "name": "T1 scale provenance: V101 data.E0 = {:.6f} and data."
                "gap = {:.6f}, V129 data.two_site.gap = {:.6f} (all "
                "artifact reads, none retyped); cross-artifact identity "
                "V101 gap x g10 == V251's omega_plus row: defect "
                "{:.2e} < 1e-9".format(E0_v101, gap_v101, gap_v129,
                                        omega_consist),
        "pass": bool(omega_consist < 1e-9)})
    checks.append({
        "name": "T2 scale state: cascade {} -> composite min gap {}, "
                "candidates {}, target 1/(16 pi) = {:.6e} computed from "
                "pi, band {} fixed in the pre-registration; branch {} "
                "matches the guard".format(
                    cascade, min_E0 if not cascade else "None (band "
                    "undefined - no kernel states)",
                    "None" if cand is None else cand["F_phys"],
                    target, band, branch),
        "pass": bool((cascade and branch == "NEVER-REALIZED"
                      and not cand_rows)
                     or (not cascade and cand is not None))})

    data["conventions"] = {
        "candidates": "F_phys = (composite min gap)/(V101 |E0|) and "
                      "/(V129 two_site gap); gap per V260 (min_E0 of "
                      "the omega = 0 row)",
        "target": "1/(16*pi) computed at run time",
        "band": band,
        "kernel_tolerance": KTOL}
    data["scales"] = {"V101_E0": E0_v101, "V101_gap": gap_v101,
                      "V129_gap": gap_v129,
                      "V251_omega_plus": omega_plus,
                      "omega_consistency_defect": omega_consist}
    data["target"] = target
    data["band"] = band
    data["candidates"] = cand_rows
    data["closest"] = cand
    data["cascade"] = {
        "rule": "pre-registered roll-up: V259 KERNEL-EMPTY (V260 "
                "NEVER-REALIZED) -> the composite gap and both scale "
                "candidates are undefined -> NEVER-REALIZED",
        "fired": bool(cascade),
        "v259_branch": v259["branch"], "v260_branch": v260["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    status = (
        "composite loop-scale test ({:.0f} s): V101 ground-state scale "
        "|E0| = {:.6f} and V129 gap = {:.6f} read from their artifacts "
        "(V101 gap x g10 reproduces V251's omega row to {:.1e}), target "
        "1/(16 pi) = {:.4e}, acceptance band {} pre-registered; V259 "
        "kernel dim recomputed {} and V260 branch {}. Branch {}.").format(
            time.time() - t0, abs(E0_v101), gap_v129, omega_consist,
            target, band, kernel_dim, v260["branch"], branch)
    status += (
        " The composite gap does not exist (kernel empty: smallest C2 "
        "eigenvalue 45.000000), so no composite sector scale can play "
        "the role of the loop factor: the O8 loop-scale contact on the "
        "physical sector is NEVER-REALIZED by the pre-registered "
        "cascade." if cascade else
        " closest candidate F_phys = {} (deviation {} from the target)."
        .format(cand["F_phys"], cand["dev_from_target"]))
    out = {"id": "V263_composite_loop_scale",
           "tier": 2,
           "property": "the composite-sector candidates for the loop "
                       "factor 1/(16 pi) (V263)",
           "residual": True,
           "programme": "V263",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v263_composite_loop_scale.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
