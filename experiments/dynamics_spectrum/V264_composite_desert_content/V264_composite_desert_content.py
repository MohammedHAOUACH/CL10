#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V264 - the composite desert content (tier 3, residual: true, F2,
hard tier).

**The question (pre-registered).**  Does the composite singlet
multiplicity table supply desert content that V237's threshold rule can
realize?

**The construction (closed input set).**  The composite multiplicity
table: the joint-Cartan decomposition of the V259 singlet kernel (the
V262 Cartan machinery: P-block V67-convention planes + the flux-block
Cartan as the centralizer of a generic element in the committed V79
structure constants, seed 0x62).  Kernel states are annihilated by all
45 total charges (C2 = sum Q^2 is PSD and w^dag C2 w = 0), so every
kernel state is an internal Spin(10) singlet: it contributes
db = (0, 0, 0) to the desert beta functions and carries no SM (1,2)
branching - both measured, not assumed (the Cartan weights of each
kernel basis vector are recorded).  The rule: V237's threshold matching
copied VERBATIM from experiments/couplings_unification/
V237_desert_threshold_matching/V237_desert_threshold_matching.py
l. 142-225 (the closed V64 content library, UN.M_Z/A_Z/B_SM,
alpha_X = pi/64, the pinned M_X from config desert.m_x_closed_gev, the
fixed-M_X two-row solve, the (1,2)-pair demand dz12 against the
log-spacing spac12).  Acceptance (V237's own inequality, verbatim): the
demand is realizable iff |dz12| <= spac12.

**Cascade (pre-registered roll-up + V264 entry).**  V259 KERNEL-EMPTY:
the multiplicity table has no rows at all (no kernel states), so the
desert content is empty and the branch is NO-CONTENT.

**Acceptance (pre-registered).**  F2-COMPOSITE-CLOSED iff the (1,2)-
shift demand <= the log-spacing to M_X; else F2-COMPOSITE-FIRED with
the shift ratio; else the NO-CONTENT cascade.
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
sys.path.insert(0, exp("V64_desert"))
import V64_desert as V64                       # noqa: E402
from cl10 import altcell as AC                 # noqa: E402
from cl10 import config as CFG                 # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import unification as UN             # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10
NPART = 1568
KTOL = 1e-9
TWO_PI = 2.0 * math.pi


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V264 {:s} t+{:6.1f}s] {}".format(
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
    v237 = art("V237_desert_threshold_matching")

    # ------------------------------------------------------------------
    # the V259 recipe (shared construction) + Cartan controls
    # ------------------------------------------------------------------
    Ts = root_generators()
    f, k_tr, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    _TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
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
    log("kernel recomputed here: dim {}".format(kernel_dim))

    # flux-block Cartan (V262 convention): centralizer of a generic
    # element in the committed structure constants
    rng = np.random.default_rng(0x62)
    c0 = rng.standard_normal(45)
    adz = np.einsum("a,abc->bc", c0, f)
    _u, s_ad, vt_ad = np.linalg.svd(adz)
    cart_vecs = vt_ad[np.abs(s_ad) < 1e-9]
    rank = len(cart_vecs)
    CartanQ = [sum(v[a] * QQ[a] for a in range(45))
               for v in cart_vecs[:5]]

    # ------------------------------------------------------------------
    # the composite multiplicity table (drafted; runs iff non-empty)
    # ------------------------------------------------------------------
    multiplicity_table = []
    if kernel_dim:
        Wn = E.asnumpy(W)
        weights = []
        for i in range(Wn.shape[1]):
            v = Wn[:, i]
            weights.append(tuple(round(float(np.real(
                np.vdot(v, Cq @ v))), 9) for Cq in CartanQ))
        counts = {}
        for wt in weights:
            counts[wt] = counts.get(wt, 0) + 1
        multiplicity_table = [{"cartan_weight": list(wt), "mult": m,
                               "db": [0.0, 0.0, 0.0],
                               "note": "kernel states are annihilated "
                                       "by all 45 total charges (ker "
                                       "C2): Spin(10) singlets, db "
                                       "identically zero, no SM (1,2) "
                                       "branching"}
                              for wt, m in sorted(counts.items())]
    n_12 = sum(r["mult"] for r in multiplicity_table
               if any(abs(x) > KTOL for x in r["cartan_weight"]))
    nB = n_12   # (1,2)-carrying multiplets: none can come from singlets

    # ------------------------------------------------------------------
    # V237's rule, copied verbatim (l. 142-225), run on the composite
    # content AND (as the fidelity control) on V237's own content
    # ------------------------------------------------------------------
    def v237_rule(dA, dB, nA, nBc):
        MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
        AX = math.pi / 64.0
        MXF = CFG.CFG["desert"]["m_x_closed_gev"]
        x_pref = [math.log(1.2e7 / MZ), math.log(1.0e12 / MZ)]
        xfix = math.log(MXF / MZ)
        C = [(1.0 / A_Z[i] - 1.0 / AX) * TWO_PI
             - xfix * (B_SM[i] + nA * dA[i] + nBc * dB[i])
             for i in range(3)]
        R = [[-nA * dA[i], -nBc * dB[i]] for i in range(3)]
        Mm = [[R[0][k] - R[1][k] for k in (0, 1)],
              [R[1][k] - R[2][k] for k in (0, 1)]]
        v = [C[0] - C[1], C[1] - C[2]]
        zA, zB = np.linalg.solve(np.array(Mm), np.array(v))
        dz10 = zA - x_pref[0]
        dz12 = zB - x_pref[1]
        spac12 = xfix - x_pref[1]
        return {"zA_GeV": MZ * math.exp(zA), "zB_GeV": MZ * math.exp(zB),
                "dz10": dz10, "dz12": dz12, "spac12": spac12,
                "shift_over_spacing_12": abs(dz12) / spac12,
                "MXF": MXF, "alpha_X": AX}

    # fidelity control (runs always): the copied rule on V237's OWN
    # closed content must reproduce V237's committed artifact
    dA_lib, dB_lib = V64.LIB["10+10bar"], V64.LIB["(1,2)"]
    fid = v237_rule(dA_lib, dB_lib, 3, 1)
    ref = v237["data"]["matching"]
    fid_dev = max(
        abs(fid["dz12"] - ref["dz12_efolds"]),
        abs(fid["spac12"] - ref["spacing12"]),
        abs(fid["dz10"] - ref["dz10_efolds"]),
        abs(fid["shift_over_spacing_12"] -
            ref["shift_over_spacing_12"]),
        abs(fid["zB_GeV"] - ref["zB_solved_GeV"]) /
        ref["zB_solved_GeV"])
    log("V237 rule fidelity: max deviation vs artifact {:.2e}".format(
        fid_dev))

    cascade = (v259["branch"] == "KERNEL-EMPTY" and kernel_dim == 0)
    shift_demand = None
    ratio = None
    if not cascade:
        if nB == 0:
            branch = "NO-CONTENT"             # no (1,2) carrier exists
        else:
            comp = v237_rule(dA_lib, dB_lib, 0, nB)
            shift_demand = comp["dz12"]
            ratio = comp["shift_over_spacing_12"]
            branch = ("F2-COMPOSITE-CLOSED"
                      if abs(shift_demand) <= comp["spac12"]
                      else "F2-COMPOSITE-FIRED")
    else:
        branch = "NO-CONTENT"

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, V260 branch {}, kernel "
                "recomputed here dim {} -> composite multiplicity table "
                "has {} rows".format(
                    v259["branch"], v259["data"]["total_singlet_count"],
                    v260["branch"], kernel_dim, len(multiplicity_table)),
        "pass": bool(v259["branch"] == "KERNEL-EMPTY"
                     and int(v259["data"]["total_singlet_count"]) == 0
                     and kernel_dim == 0 and not multiplicity_table)})
    checks.append({
        "name": "T1 V237 rule fidelity: the copied threshold-matching "
                "code (l. 142-225 conventions: V64 closed content, "
                "UN.M_Z/A_Z/B_SM, alpha_X = pi/64, M_X pinned at the "
                "committed desert.m_x_closed_gev = {:.4e}) reproduces "
                "V237's own artifact (dz12, dz10, spacing12, "
                "shift/spacing, zB) to max deviation {:.2e} < 1e-9; the "
                "flux Cartan rank is {}".format(
                    fid["MXF"], fid_dev, rank),
        "pass": bool(fid_dev < 1e-9 and rank == 5)})
    checks.append({
        "name": "T2 content state: cascade {} -> multiplets {}, (1,2)-"
                "carrying multiplets {}, shift demand {}, ratio {}, "
                "branch {} matches the guard (F2 inequality: |dz12| <= "
                "spacing12 = {:.3f} e-folds, V237's own)".format(
                    cascade, len(multiplicity_table), nB, shift_demand,
                    ratio, branch, ref["spacing12"]),
        "pass": bool((cascade and branch == "NO-CONTENT"
                      and shift_demand is None)
                     or (not cascade and nB > 0
                         and ratio is not None))})

    data["conventions"] = {
        "rule_source": "experiments/couplings_unification/"
                       "V237_desert_threshold_matching/V237_desert_"
                       "threshold_matching.py l. 142-225 copied verbatim",
        "content": "composite multiplicity table: joint Cartan weights "
                   "of the singlet kernel (V262 Cartan convention); "
                   "kernel states are Spin(10) singlets (annihilated "
                   "by all 45 charges) -> db = (0, 0, 0), no (1,2) "
                   "branching",
        "inequality": "realizable iff |dz12| <= spacing12 (V237's own)",
        "kernel_tolerance": KTOL}
    data["composite_multiplicity_table"] = multiplicity_table
    data["n_12_carrying"] = nB
    data["shift_demand_dz12"] = shift_demand
    data["ratio"] = ratio
    data["fidelity_control"] = {"deviation": fid_dev,
                                "recomputed": {k: fid[k] for k in
                                               ("dz12", "dz10",
                                                "spac12",
                                                "shift_over_spacing_12")}}
    data["cascade"] = {
        "rule": "V259 KERNEL-EMPTY (empty multiplicity table) -> "
                "NO-CONTENT (pre-registered cascade branch)",
        "fired": bool(cascade),
        "v259_branch": v259["branch"], "v260_branch": v260["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    status = (
        "composite desert content vs V237's threshold rule ({:.0f} s): "
        "the copied rule reproduces V237's own artifact to {:.1e}; the "
        "V259 kernel recomputed here has dim {}, so the composite "
        "multiplicity table has {} rows and the (1,2)-carrying count is "
        "{}. Branch {}.").format(
            time.time() - t0, fid_dev, kernel_dim,
            len(multiplicity_table), nB, branch)
    status += (
        " The singlet kernel is empty: the composite sector supplies no "
        "desert content at all (kernel states, had they existed, are "
        "internal Spin(10) singlets with db = (0, 0, 0) and no SM (1,2) "
        "branching), so there is no threshold the (1,2)-shift demand "
        "could act on: NO-CONTENT by the pre-registered cascade; V237's "
        "original verdict (NOT REALIZABLE BY ITS OWN CONTENT) stands "
        "unchanged and the copied rule was proven faithful to it."
        if cascade else
        " composite (1,2) content {} -> demand {} e-folds vs spacing "
        "{}.".format(nB, shift_demand, ref["spacing12"]))
    out = {"id": "V264_composite_desert_content",
           "tier": 3,
           "property": "the composite singlet content fed to V237's "
                       "desert threshold rule (V264)",
           "residual": True,
           "programme": "V264",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v264_composite_desert_content.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
