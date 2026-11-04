#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V261 - the composite F9 contact (tier 5, residual: true, F9).

**The question (pre-registered).**  At composite gapless loci, does the
velocity anisotropy meet the SME ceiling read from
_config.CFG.bounds.sme_electron_velocity_aniso?

**The construction (closed input set).**  V260's gapless loci (its
artifact data.gapless_loci); the same 12-ray finite-difference gradient
V254 uses (AC.sphere_dirs(6), dq = 1e-3) with V235's normalization
dvv = (v_max - v_min)/v_max, applied to the composite band of the
V259 kernel recipe (H_k(q) = W^dag H2(q) W, drafted below exactly as in
V260); the bound read from config at run time, never retyped.

**Cascade (pre-registered roll-up + V261 entry).**  V259 landed
KERNEL-EMPTY and V260 NEVER-REALIZED: there is no composite gapless
locus, so the branch is F9-COMPOSITE-NO-LOCUS (recorded, exit 0).  The
dvv path below is the pre-registered non-cascade branch and runs iff
V260 landed COMPOSITE-GAPLESS.

**Acceptance (pre-registered).**  F9-COMPOSITE-CLOSED iff best dvv < the
config bound; F9-COMPOSITE-FIRED with the ratio; F9-COMPOSITE-NO-LOCUS
if V260 landed GAPPED or V259 KERNEL-EMPTY.
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
from cl10 import config as _config             # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10
NPART = 1568
KTOL = 1e-9
DQ = 1e-3                     # V254's finite-difference step


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V261 {:s} t+{:6.1f}s] {}".format(
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

    v251 = art("V251_gauged_two_particle_sector")
    v259 = art("V259_singlet_kernel_census")
    v260 = art("V260_composite_gap_ladder")
    bound = float(_config.CFG.bounds.sme_electron_velocity_aniso)
    GRID_N = int(v251["data"]["grid"]["axes"])
    omegas = {str(r): float(v) for r, v in
              v251["data"]["sensitivity"]["omega_rows"].items()}
    loci = v260["data"]["gapless_loci"]

    # ------------------------------------------------------------------
    # own recomputation of the kernel (shared V259 recipe)
    # ------------------------------------------------------------------
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
    kernel_dim = int(np.count_nonzero(E.asnumpy(w).real < KTOL))
    W = V[:, E.asnumpy(w).real < KTOL]
    log("kernel recomputed here: dim {}".format(kernel_dim))

    cascade = (v259["branch"] == "KERNEL-EMPTY"
               or v260["branch"] != "COMPOSITE-GAPLESS"
               or not loci)

    dvv_rows = []
    best = None
    if not cascade:
        # ------- drafted non-cascade path (composite gapless loci) ----
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
            H[:128, 128:] = G10 * C.conj().T
            H[128:, :128] = G10 * C
            return H

        Wn = E.asnumpy(W)

        def E0_at(qq, om):
            H1 = AC.effective_hamiltonian(
                lambda pz: SC.bloch(pz, U16, G10), qq)
            H2 = build_H2(H1, TL, coupling_at(qq), om)
            Hk = Wn.conj().T @ H2 @ Wn
            return (float(np.linalg.eigvalsh(Hk)[0]) - (-32.0 + om))

        dirs = AC.sphere_dirs(6)
        for row in loci:
            q0 = np.asarray(row["q"], float)
            om = omegas[row["row"]]
            vs = [abs(E0_at(tuple(a + DQ * x for a, x in zip(q0, d)), om)
                      - row["E0"]) / DQ for d in dirs]
            vmax, vmin = max(vs), min(vs)
            dvv_rows.append({**row, "v_max": vmax, "v_min": vmin,
                             "dvv": (vmax - vmin) / vmax if vmax > 0
                                    else 1.0})
        best = min(dvv_rows, key=lambda r: r["dvv"])
        fired = best["dvv"] >= bound
        branch = "F9-COMPOSITE-FIRED" if fired else "F9-COMPOSITE-CLOSED"
    else:
        branch = "F9-COMPOSITE-NO-LOCUS"

    ratio = None if best is None else best["dvv"] / bound

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, V260 branch {} with {} gapless "
                "loci, kernel recomputed here dim {} (agreement {})"
                .format(v259["branch"],
                        v259["data"]["total_singlet_count"],
                        v260["branch"], len(loci), kernel_dim,
                        kernel_dim == int(v259["data"]
                                          ["total_singlet_count"])),
        "pass": bool(v259["branch"] == "KERNEL-EMPTY"
                     and int(v259["data"]["total_singlet_count"]) == 0
                     and kernel_dim == 0)})
    checks.append({
        "name": "T1 bound provenance: sme_electron_velocity_aniso read "
                "from _config.CFG.bounds at run time = {:.3e} (never "
                "retyped); 12-ray loop AC.sphere_dirs(6), dq = 1e-3, "
                "V235 normalization (V254's convention)".format(bound),
        "pass": bool(bound > 0)})
    checks.append({
        "name": "T2 contact state: cascade {} -> best dvv {}, ratio {}, "
                "dvv rows {} (== V260 loci count {}); branch {} matches "
                "the guard".format(
                    cascade, "None" if best is None else best["dvv"],
                    "None" if ratio is None else "{:.3e}".format(ratio),
                    len(dvv_rows), len(loci), branch),
        "pass": bool((cascade and not dvv_rows and best is None
                      and branch == "F9-COMPOSITE-NO-LOCUS")
                     or (not cascade and best is not None))})

    data["conventions"] = {
        "bound_source": "config bounds.sme_electron_velocity_aniso",
        "dvv": "V235 normalization (v_max - v_min)/v_max, 12 rays "
               "AC.sphere_dirs(6), dq = 1e-3 (V254's loop), applied to "
               "the composite band one-sided gradients at the V260 loci",
        "grid_axes": GRID_N, "omega_rows": omegas,
        "kernel_tolerance": KTOL}
    data["loci_from_V260"] = loci
    data["dvv"] = dvv_rows
    data["best"] = best
    data["bound"] = bound
    data["ratio_best"] = ratio
    data["cascade"] = {
        "rule": "V260 loci empty or V259 KERNEL-EMPTY -> "
                "F9-COMPOSITE-NO-LOCUS (pre-registered)",
        "fired": bool(cascade),
        "v259_branch": v259["branch"], "v260_branch": v260["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    status = (
        "composite F9 contact (grid {}^3, {:.0f} s): the V259 recipe "
        "recomputed here gives kernel dim {}, V259 branch {}, V260 "
        "branch {} with {} gapless loci, the SME bound reads {:.3e} "
        "from config. Branch {}.").format(
            GRID_N, time.time() - t0, kernel_dim, v259["branch"],
            v260["branch"], len(loci), bound, branch)
    status += (
        " There is no composite gapless locus - the 1568-dim fiber "
        "kernel is empty (smallest C2 eigenvalue 45.000000 per V259) - "
        "so the F9 velocity-anisotropy contact has nothing to measure "
        "on the physical sector: NO-LOCUS by the pre-registered rule; "
        "the 12-ray dvv path ran not a line." if cascade else
        " best dvv {} vs bound {}, ratio {}.".format(
            best["dvv"], bound, ratio))
    out = {"id": "V261_composite_f9_contact",
           "tier": 5,
           "property": "the SME velocity-anisotropy contact measured on "
                       "the composite gapless loci (V261)",
           "residual": True,
           "programme": "V261",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v261_composite_f9_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
