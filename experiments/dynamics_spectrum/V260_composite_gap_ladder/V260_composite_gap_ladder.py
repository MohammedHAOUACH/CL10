#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V260 - the composite gap ladder on the singlet kernel (tier 2,
residual: false - structural, hard tier).

**The question (pre-registered).**  Is the composite projected band
gapless anywhere in the zone, at the published N = 6 grid?

**The construction (closed input set; the V259 kernel recipe).**  V251's
fiber, charges and build_H2 verbatim (grid axes and omega rows read from
V251's artifact); C2 = sum_a Q_G^a Q_G^a, kernel basis W = V[:, w < 1e-9];
composite band E0(q) = lowest eigenvalue of H_k(q) = W^dag H2(q) W above
the q-independent ground offset: the free-tower floor (-32 + omega)
(the only exactly q-independent spectral floor of H2 - recorded
explicitly here and in data.conventions; E0 < 1e-9 therefore means
"the composite band touches or undercuts the free flux-tower floor").
Gapless loci = q with E0 < 1e-9 (repo-wide tolerance); the velocity
tensor at those loci by V235's normalization dvv = (v_max - v_min)/v_max
along the 12 star rays (AC.sphere_dirs(6)) at dq = 1e-3, V254's
measurement loop; V129's committed gap read from its artifact.

**Measured (pre-registered).**  min E0 over the grid per omega row; the
gapless loci; dvv at those loci; the ratio of min gap to V129's gap.
Cross-check (new-behavior proof): at q = 0 the composite kernel dimension
and the C2 kernel eigenvalue count must equal the number of singlet
copies computed tower-sector by tower-sector (restricted C2 spectra on
the T.L eigenspaces) - asserted equal below 1e-9; a wrong projection
shows up here.

**Cascade (pre-registered roll-up).**  V259 landed KERNEL-EMPTY: the
kernel is empty, so the composite band does not exist and the gap
ladder is recorded NEVER-REALIZED.  The full ladder path below is the
drafted non-cascade branch of the pre-registration: it runs iff the
kernel is non-empty.
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
NPART = 1568                  # 128 + 1440
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
        print("[V260 {:s} t+{:6.1f}s] {}".format(
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
    v129 = art("V129_gauged_layer_gap")
    GRID_N = int(v251["data"]["grid"]["axes"])
    omegas = {str(r): float(v) for r, v in
              v251["data"]["sensitivity"]["omega_rows"].items()}
    v129_gap = float(v129["data"]["two_site"]["gap"])

    # ------------------------------------------------------------------
    # the V259 kernel recipe, recomputed here (shared construction)
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
    w_np = E.asnumpy(w).real
    ker_mask = w_np < KTOL
    kernel_dim = int(np.count_nonzero(ker_mask))
    W = V[:, ker_mask]
    log("kernel recomputed here: dim {}".format(kernel_dim))

    # ------------------------------------------------------------------
    # q = 0 cross-check (new-behavior proof): kernel count == C2 kernel
    # eigenvalue count == singlet copies resolved tower-sector by
    # tower-sector (restricted C2_Q spectra on the T.L eigenspaces)
    # ------------------------------------------------------------------
    c2q = np.matmul(np.stack(QQ), np.stack(QQ)).sum(axis=0)
    c2q = 0.5 * (c2q + c2q.conj().T)
    tw_ev, tw_V = np.linalg.eigh(TL)
    sector_rows = []
    n_sector = 0
    leak_max = 0.0
    for s in (-32.0, -12.0, 4.0):
        idx = np.abs(tw_ev - s) < 1e-9
        Vs = tw_V[:, idx]                        # 1440 x m, orthonormal
        # the tower sector must be C2_Q-invariant for the restriction to
        # be a genuine block: leakage of (I - Vs Vs^dag) C2_Q Vs
        cVs = c2q @ Vs
        leak = float(np.linalg.norm(cVs - Vs @ (Vs.conj().T @ cVs)))
        leak_max = max(leak_max, leak)
        R = Vs.conj().T @ cVs
        R = 0.5 * (R + R.conj().T)
        rw = np.linalg.eigvalsh(R)
        cnt = int(np.count_nonzero(rw < KTOL))
        sector_rows.append({"sector": s, "dim": int(idx.sum()),
                            "singlet_copies": cnt,
                            "min_restricted_C2": float(rw[0])})
        n_sector += cnt
    xcheck_ok = (kernel_dim == n_sector and leak_max < KTOL)
    log("q=0 cross-check: kernel {} vs tower-sector singlet sum {} "
        "(max sector leakage {:.1e})".format(kernel_dim, n_sector,
                                             leak_max))

    # ------------------------------------------------------------------
    # cascade (pre-registered roll-up): V259 KERNEL-EMPTY
    # ------------------------------------------------------------------
    cascade = (v259["branch"] == "KERNEL-EMPTY" and kernel_dim == 0)
    consistency = (v259["branch"] == "KERNEL-EMPTY") == (
        int(v259["data"]["total_singlet_count"]) == 0)
    v259_dims_all_zero = all(r["kernel_dim"] == 0
                             for r in v259["data"]["kernel_dims"])

    grid = [(2 * math.pi * i / GRID_N, 2 * math.pi * j / GRID_N,
             2 * math.pi * l / GRID_N)
            for i in range(GRID_N) for j in range(GRID_N)
            for l in range(GRID_N)]

    min_E0 = {}
    loci_rows = []
    dvv_rows = []
    if not cascade:
        # ------- drafted non-cascade path (kernel non-empty) ----------
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
        for rname, om in omegas.items():
            floor = -32.0 + om                  # q-independent offset
            vals = []
            for q in grid:
                H1 = AC.effective_hamiltonian(
                    lambda qq: SC.bloch(qq, U16, G10), q)
                H2 = build_H2(H1, TL, coupling_at(q), om)
                Hk = Wn.conj().T @ H2 @ Wn
                vals.append(float(np.linalg.eigvalsh(Hk)[0]) - floor)
            min_E0[rname] = min(vals)
            for q, e in zip(grid, vals):
                if e < KTOL:
                    loci_rows.append({"q": [float(x) for x in q],
                                      "E0": e, "row": rname})
        dirs = AC.sphere_dirs(6)
        for row in loci_rows:
            q0 = np.asarray(row["q"], float)

            def E0_at(qq):
                H1 = AC.effective_hamiltonian(
                    lambda pz: SC.bloch(pz, U16, G10), qq)
                H2 = build_H2(H1, TL, coupling_at(qq),
                              omegas[row["row"]])
                Hk = Wn.conj().T @ H2 @ Wn
                return (float(np.linalg.eigvalsh(Hk)[0])
                        - (-32.0 + omegas[row["row"]]))

            vs = [abs(E0_at(tuple(a + DQ * x for a, x in zip(q0, d)))
                      - row["E0"]) / DQ for d in dirs]
            vmax, vmin = max(vs), min(vs)
            dvv_rows.append({**row,
                             "v_max": vmax, "v_min": vmin,
                             "dvv": (vmax - vmin) / vmax if vmax > 0
                                    else 1.0})
        branch = "COMPOSITE-GAPLESS" if loci_rows else "COMPOSITE-GAPPED"
    else:
        branch = "NEVER-REALIZED"

    ratio = None
    if not cascade and min_E0.get("0") is not None:
        ratio = min_E0["0"] / v129_gap

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, per-q kernel dims all zero {}, "
                "kernel recomputed here dim {} (agreement {})".format(
                    v259["branch"], v259["data"]["total_singlet_count"],
                    v259_dims_all_zero, kernel_dim,
                    kernel_dim == int(v259["data"]["total_singlet_count"])),
        "pass": bool(consistency and v259_dims_all_zero
                     and kernel_dim == int(v259["data"]
                                           ["total_singlet_count"]))})
    checks.append({
        "name": "T1 q = 0 projection cross-check: composite kernel dim "
                "{} == singlet copies resolved tower-sector by "
                "tower-sector {} (sectors {}, max T.L-eigenspace "
                "C2_Q-leakage {:.1e} < 1e-9: each sector is an exact "
                "block, so the restricted count is a genuine "
                "independent census)".format(
                    kernel_dim, n_sector, sector_rows, leak_max),
        "pass": bool(xcheck_ok)})
    checks.append({
        "name": "T2 ladder state: cascade {} -> composite band undefined "
                "(no kernel state to disperse); min E0 {}, gapless loci "
                "{}, V129 committed gap {:.6f} booked for the ledger, "
                "ratio {}".format(
                    cascade, "None" if cascade else
                    {k: round(v, 12) for k, v in min_E0.items()},
                    len(loci_rows), v129_gap,
                    "None" if ratio is None else "{:.6g}".format(ratio)),
        "pass": bool(cascade or min_E0)})

    data["conventions"] = {
        "grid_axes_from_V251_artifact": GRID_N,
        "grid_points": len(grid),
        "omega_rows": omegas,
        "kernel_tolerance": KTOL,
        "offset": "E0(q) = lowest eigenvalue of H_k(q) = W^dag H2(q) W "
                  "minus the q-independent free-tower floor (-32 + "
                  "omega) - the only exactly q-independent spectral "
                  "floor of H2; gapless = E0 < 1e-9 (touch or undercut)",
        "dvv": "V235 normalization (v_max - v_min)/v_max along the 12 "
               "rays AC.sphere_dirs(6), dq = 1e-3 (V254's loop)"}
    data["kernel"] = {"dim_recomputed_here": kernel_dim,
                      "dim_V259_artifact":
                          int(v259["data"]["total_singlet_count"]),
                      "tower_sector_singlets": sector_rows,
                      "xcheck_max_sector_leakage": leak_max}
    data["min_E0"] = min_E0 if not cascade else None
    data["gapless_loci"] = loci_rows
    data["dvv"] = dvv_rows
    data["v129_gap"] = v129_gap
    data["gap_to_v129_ratio"] = ratio
    data["cascade"] = {
        "rule": "pre-registered roll-up: V259 KERNEL-EMPTY (no q-point "
                "carries a singlet) -> V260-V266 recorded NEVER-REALIZED",
        "fired": bool(cascade),
        "v259_branch": v259["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    c2_first = v259["data"]["c2_spectrum"]["first_nonzero"]
    status = (
        "composite gap ladder on the V259 singlet kernel (grid {}^3, "
        "omega rows from V251, tolerance 1e-9, {:.0f} s): the kernel "
        "recomputed here has dim {}, the q = 0 cross-check kernel == "
        "tower-sector singlet sum ({}) holds (max sector leakage "
        "{:.1e}), the V129 committed gap is {:.6f}. Branch {}.").format(
            GRID_N, time.time() - t0, kernel_dim, n_sector, leak_max,
            v129_gap, branch)
    status += (
        " V259 landed KERNEL-EMPTY: the 1568-dim coupled fiber carries "
        "no gauge singlet (smallest C2 eigenvalue {:.6f} from V259), so "
        "no composite (matter + flux) band exists to disperse - the gap "
        "ladder, its gapless loci and the dvv velocity tensor are "
        "NEVER-REALIZED on this construction; the drafted ladder path "
        "above is the pre-registered non-cascade branch and ran not a "
        "line.".format(c2_first) if cascade else
        " min E0 per omega row {}, {} gapless loci, ratio to the V129 "
        "gap {}.".format(min_E0, len(loci_rows), ratio))
    out = {"id": "V260_composite_gap_ladder",
           "tier": 2,
           "property": "the composite band on the singlet kernel: min "
                       "gap, gapless loci and the dvv velocity tensor "
                       "at them (V260)",
           "residual": False,
           "programme": "V260",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v260_composite_gap_ladder.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
