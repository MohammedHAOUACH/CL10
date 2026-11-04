#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V265 - the composite F8 contact (tier 5, residual: true, F8).

**The question (pre-registered).**  Does the composite sector carry a
spin-2 channel propagating at c = 1?

**The construction (closed input set).**  V229's T2g stress vertex:
Emat (polarization tensors, l. 128) and spin2_generators (l. 154)
imported from experiments/gauged_layer/V229_diagonal_bond_quintet/
V229_diagonal_bond_quintet.py (module-level, __main__-guarded; the
quintet basis q5 and the T2g off_rows come from the same source),
applied along the star axes: the committed set is SIX star axes
(cl10.starcell.star_axes = continuum.module_generators vectors, "six
diagonal star bonds", V229 T2 / V214; V229's rays add the +- sign
reversals = 12 directed rays, l. 234).  The plan's phrase "15 star
axes" has no committed referent in the repo (grep: no other file
mentions 15 axes); following the plan's never-guess rule the V229
source is used and the discrepancy recorded here.
Drafted path (runs iff the kernel is non-empty): composite
ground/first-excited states = eigenvectors of H_k(q) = W^dag H2(q) W
(V259 recipe, V260 conventions); the stress vertex along axis a in
polarization channel g is V229's scalar d^T E^gamma d (T1 convention,
l. 195) times the matter-block axis quadrupole Phi_a =
diag_e((r_e . u_a)^2) kron I_16 on the P block, zero on the flux block
(zero-parameter lift, recorded here); the channel is T2g-selected by
the matrix elements <psi_1| W^dag Phi_a W |psi_0> times the Emat
coefficients; the composite c_gw is the slope of the linear part of
E_1(q) - E_0(q) along the T2g-selected star axis (delta = 1e-3..5e-3
least squares).  The acceptance band |c_gw - 1| < 1e-15 is V258's
(gw_band = 1e-15, "V208's GW170817 contact",
V258_falsifier_conjunction.py l. 96-97).

**Cascade (pre-registered roll-up).**  V259 KERNEL-EMPTY: the kernel
carries no ground and no first-excited composite state, so the spin-2
channel and c_gw are recorded NEVER-REALIZED.

**Acceptance (pre-registered).**  F8-COMPOSITE-OK iff |c_gw - 1| <
1e-15; else F8-COMPOSITE-FIRED with the ratio (channel-absent recorded
as deviation undefined).
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
sys.path.insert(0, exp("V229_diagonal_bond_quintet"))
from V229_diagonal_bond_quintet import (Emat, spin2_generators,   # noqa
                                        q5, diag_rows, off_rows)
from cl10 import altcell as AC                 # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import starcell as X                 # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10
NPART = 1568
KTOL = 1e-9
GW_BAND = 1e-15               # V258 l. 96-97: V208's GW170817 band


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V265 {:s} t+{:6.1f}s] {}".format(
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

    # ------------------------------------------------------------------
    # the V259 recipe (shared construction)
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
    km = E.asnumpy(w).real < KTOL
    kernel_dim = int(np.count_nonzero(km))
    W = V[:, km]
    log("kernel recomputed here: dim {}".format(kernel_dim))

    # ------------------------------------------------------------------
    # V229 quintet machinery controls (run regardless of the cascade)
    # ------------------------------------------------------------------
    gram = np.array([[float(np.real(np.trace(Emat(g).conj().T
                                             @ Emat(hh))))
                      for hh in range(5)] for g in range(5)])
    gram_dev = float(np.abs(gram - np.eye(5)).max())
    q5_def = float(np.abs(q5 @ q5.conj().T - np.eye(5)).max())
    J6 = spin2_generators()
    eps3 = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
            (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
    j_algebra = 0.0
    for k in range(3):
        for l in range(3):
            cm = J6[k] @ J6[l] - J6[l] @ J6[k]
            tgt = sum(1j * eps3.get((k, l, m), 0) * J6[m]
                      for m in range(3))
            j_algebra = max(j_algebra, float(np.abs(cm - tgt).max()))
    star = [np.asarray(a, float) for a in X.star_axes()]
    axes = [a / np.linalg.norm(a) for a in star]
    n_axes = len(axes)
    log("quintet machinery: Emat gram {:.1e}, q5 {:.1e}, J so(3) "
        "algebra {:.1e}, star axes {}".format(gram_dev, q5_def,
                                              j_algebra, n_axes))

    # ------------------------------------------------------------------
    # drafted composite contact (runs iff the kernel is non-empty)
    # ------------------------------------------------------------------
    cascade = (v259["branch"] == "KERNEL-EMPTY" and kernel_dim == 0)
    c_gw = None
    deviation = None
    channel = []
    r_lat = np.array([(e & 1, (e >> 1) & 1, (e >> 2) & 1)
                      for e in range(8)], dtype=float)
    if not cascade:
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

        def bands_at(q):
            H1 = AC.effective_hamiltonian(
                lambda pz: SC.bloch(pz, U16, G10), q)
            H2 = build_H2(H1, TL, coupling_at(q), 0.0)
            Hk = Wn.conj().T @ H2 @ Wn
            ww, VV = np.linalg.eigh(Hk)
            return ww, VV

        ww0, V0 = bands_at((0.0, 0.0, 0.0))
        psi0, psi1 = V0[:, 0], V0[:, 1] if kernel_dim >= 2 else None
        for g in off_rows:                    # T2g quintet components
            for ia, u in enumerate(axes):
                coeff = float(u @ Emat(g) @ u)
                if abs(coeff) < 1e-12:
                    continue
                Phi = np.zeros((NPART, NPART), dtype=np.complex128)
                rphi = np.diag([(r_lat[e] @ u) ** 2 for e in range(8)])
                Phi[:128, :128] = np.kron(rphi, np.eye(16))
                mel = (abs(np.vdot(psi1, Wn.conj().T @ Phi @ Wn @ psi0))
                       if psi1 is not None else 0.0)
                channel.append({"g": int(g), "axis": ia,
                                "coeff": coeff, "matrix_element": mel})
        sel = [c for c in channel if c["matrix_element"] > 1e-12]
        if not sel:
            branch = "F8-COMPOSITE-FIRED"     # channel-absent
            deviation = None
        else:
            top = max(sel, key=lambda c: c["matrix_element"])
            u = axes[top["axis"]]
            ds = np.arange(1, 6) * 1e-3
            dE = [float(bands_at(tuple(d * x for x in u))[0][1]
                        - bands_at(tuple(d * x for x in u))[0][0])
                  for d in ds]
            c_gw = float(np.polyfit(ds, dE, 1)[0])
            deviation = abs(c_gw - 1.0)
            branch = ("F8-COMPOSITE-OK" if deviation < GW_BAND
                      else "F8-COMPOSITE-FIRED")
    else:
        branch = "NEVER-REALIZED"

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, V260 branch {}, kernel "
                "recomputed here dim {} -> {} composite states to "
                "sandwich (ground/first-excited undefined)".format(
                    v259["branch"], v259["data"]["total_singlet_count"],
                    v260["branch"], kernel_dim, kernel_dim),
        "pass": bool(v259["branch"] == "KERNEL-EMPTY"
                     and int(v259["data"]["total_singlet_count"]) == 0
                     and kernel_dim == 0)})
    checks.append({
        "name": "T1 V229 quintet machinery (imported and live): Emat "
                "polarization gram = I5 to {:.1e}, q5 orthonormality "
                "{:.1e}, spin2_generators J_k close the so(3) algebra "
                "to {:.1e}, T2g off_rows {}, diag_rows {}, {} "
                "star axes = {} directed rays (V229's X.star_axes / "
                "l. 234 convention)".format(
                    gram_dev, q5_def, j_algebra, off_rows, diag_rows,
                    n_axes, 2 * n_axes),
        "pass": bool(gram_dev < 1e-9 and q5_def < 1e-9
                     and j_algebra < 1e-9 and len(off_rows) == 3
                     and len(diag_rows) == 2 and n_axes == 6)})
    checks.append({
        "name": "T2 contact state: cascade {} -> kernel state count {}, "
                "c_gw {}, band {:.0e} (V258 l. 96-97 provenance: V208's "
                "GW170817 contact), T2g channel entries {}; branch {} "
                "matches the guard".format(
                    cascade, kernel_dim, c_gw, GW_BAND, len(channel),
                    branch),
        "pass": bool((cascade and branch == "NEVER-REALIZED"
                      and c_gw is None and not channel)
                     or (not cascade and branch.startswith(
                         "F8-COMPOSITE")))})

    data["conventions"] = {
        "vertex": "V229 Emat polarization tensors (l. 128) and scalar "
                  "d^T E^gamma d convention (V229 T1); T2g = off_rows; "
                  "spin2_generators J_k (l. 154)",
        "axes": "six star axes cl10.starcell.star_axes (V229 l. 175, "
                "V214's six generators; +- reversals = 12 directed "
                "rays, V229 l. 234); plan wording '15 star axes' has "
                "no committed referent - recorded",
        "lift": "matter-block axis quadrupole diag_e((r_e.u)^2) kron "
                "I16 on the P block, zero on the flux block "
                "(zero-parameter, recorded)",
        "band": GW_BAND,
        "band_provenance": "V258_falsifier_conjunction.py l. 96-97 "
                           "(gw_band = 1e-15, V208's GW170817 contact)",
        "kernel_tolerance": KTOL}
    data["machinery_controls"] = {"emat_gram_dev": gram_dev,
                                 "q5_dev": q5_def,
                                 "j_algebra_dev": j_algebra,
                                 "off_rows": off_rows,
                                 "diag_rows": diag_rows,
                                 "n_axes": n_axes}
    data["kernel_state_count"] = kernel_dim
    data["channel"] = channel
    data["c_gw"] = c_gw
    data["deviation"] = deviation
    data["cascade"] = {
        "rule": "pre-registered roll-up: V259 KERNEL-EMPTY -> no "
                "composite ground/first-excited kernel states -> "
                "NEVER-REALIZED",
        "fired": bool(cascade),
        "v259_branch": v259["branch"], "v260_branch": v260["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    status = (
        "composite F8 contact ({:.0f} s): V229's quintet machinery "
        "verified live on import (Emat gram {:.1e}, q5 {:.1e}, J so(3) "
        "algebra {:.1e}, T2g rows {}, {} star axes = {} directed "
        "rays), V259 kernel dim recomputed {}, V260 branch {}. "
        "Branch {}.").format(
            time.time() - t0, gram_dev, q5_def, j_algebra, off_rows,
            n_axes, 2 * n_axes, kernel_dim, v260["branch"], branch)
    status += (
        " The singlet kernel is empty: there is no composite ground or "
        "first-excited state to carry a spin-2 channel, so the F8 "
        "c = 1 contact on the physical sector is NEVER-REALIZED (the "
        "band is |c_gw - 1| < 1e-15, V258's); the drafted vertex-"
        "sandwich path ran not a line." if cascade else
        " c_gw = {} (deviation {}) vs band {}.".format(
            c_gw, deviation, GW_BAND))
    out = {"id": "V265_composite_f8_contact",
           "tier": 5,
           "property": "the GW170817 c = 1 contact measured on the "
                       "composite T2g spin-2 channel (V265)",
           "residual": True,
           "programme": "V265",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v265_composite_f8_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
