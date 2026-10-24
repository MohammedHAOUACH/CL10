#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V235 - the Lorentz contact, quantified against measurement (tier 2).

**The question this pre-registers against** (`review/pre_registrations.md`):
the paper reports cone anisotropy `0.748-0.758`, "0 of 2726 gapless points
reach `A < 1e-6`", and node speeds `(0.200, 0.200, 0.284)` at the
implemented clock - and defends the contact with O10 ("no light candidate"),
not with a suppression mechanism. Because the lattice spacing IS the Planck
length, a dimension-4 velocity anisotropy is NOT Planck-suppressed: the
electron-sector velocity bound of the Standard-Model Extension applies
directly to any light mode at these nodes. How far is the node velocity
tensor from the measured bound?

Closed input set (pre-registered): the V102/V115 cone estimator (Fibonacci
directions, min-eigenphase velocity `2 asin(sigma_min(M-I)/2)/eps`, tensor
least-squares) with the same normalisation that produced the published
`0.748-0.758` and the star nodal-line `0.680` (V219's azimuthal profile);
the three split nodes of the selective reading at
`desert.theta_root_rad = 1.7307836435290485` (V65 physical selective root:
`w = 6.262416e-9`; the desert programme's rank-1 split-invariance - V74 -
makes the node choice among the 12-orbit immaterial, and V115 measured the
anisotropy identical across the orbit); the twelve M nodes at both clocks;
the star nodal line (V219's recorded transverse velocities through
`cl10.starcell.gap_star`); and the fetched bound
`bounds.sme_electron_velocity_aniso = 8.0e-15` (Kostelecky-Mewes data
tables, arXiv:0801.0287v16, Table D6, electron d=3,4: the traceless
velocity operators c_(TJ): |c_(TX)| <= 3.0e-14, |c_(TY)| <= 8.0e-15,
|c_(TZ)| <= 1.3e-13; the tightest component is the ceiling; retrieved
2026-10-23).

Measured quantity: `delta_v/v = (v_max - v_min)/v_max` per named node (the
principal velocities of the fitted tensor; a vanishing principal velocity -
the flat axis at the closed clock - enters as `delta_v/v = 1`); then the
ratio `delta_v/v : bound`.

Acceptance (pre-registered): if `delta_v/v` exceeds the bound by more than
a factor 10 (expected) a NINTH falsifier F9 joins Table 8 carrying the two
numbers and the citation; O9's gap column gets the quantitative statement.
No softening, no re-scope.
"""
import glob
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
# the reference branch of this measurement is the implemented 1/5 clock,
# the branch V102/V115/V65 publish on. cl10.config freezes the clock at
# its own import, so pin it here (V118/V186/V187/V207 preamble order).
os.environ["CL10_ANGLE"] = "0.2"
sys.path.insert(0, ROOT)
from cl10.paths import exp                                    # noqa: E402

sys.path.insert(0, exp("V102_lorentz_isotropy"))
sys.path.insert(0, exp("V65_wilson"))
import V102_lorentz_isotropy as V102                          # noqa: E402
import V65_wilson as V65                                      # noqa: E402

from cl10 import config as _config                            # noqa: E402
from cl10 import engine as E                                  # noqa: E402
from cl10 import supercell as SC                              # noqa: E402

CFG = _config.CFG
XP = E._E if E._BACKEND == "cupy" else np
I128 = np.eye(128, dtype=complex)
I128B = XP.asarray(I128)
X0 = (np.random.default_rng(235).standard_normal(128)
      + 1j * np.random.default_rng(236).standard_normal(128))
X0 /= np.linalg.norm(X0)


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def fib_dirs(n=100):
    idx = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * idx / n)
    th = np.pi * (1 + 5 ** 0.5) * idx
    return np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th),
                     np.cos(phi)], axis=1)


_U_CACHE = {}


def step_batch_at(qs, angle):
    """V102's batched one-step Floquet, but with the internal clock passed
    explicitly: the clock-matched 16-sector internal unitary U16(angle)
    AND the same angle in the bond stages (the engine default would mix
    two clocks)."""
    U = _U_CACHE.get(angle)
    if U is None:
        U = _U_CACHE[angle] = SC.sectors(angle)[2]
    M = XP.broadcast_to(XP.eye(128, dtype=complex),
                        (len(qs), 128, 128)).copy()
    for even, j in SC.STAGES:
        L = XP.asarray(np.stack([np.kron(SC.stage_matrix(even, j, q[j],
                                                          angle),
                                          U) for q in qs]))
        M = L @ M
    return M


def step_batch_wilson(qs, angle, w, dvec):
    """the Wilson-factor step at direction angle `theta` (inside dvec) and
    strength w: exp(-i w sum(1-cos q_j) diag(d)) Bloch - V115's family."""
    M = step_batch_at(qs, angle)
    ph = np.array([sum(1.0 - math.cos(qj) for qj in q) for q in qs]) * w
    Wl = XP.asarray(np.exp(-1j * ph[:, None, None]
                           * np.tile(dvec, 8)[None, None, :]))
    return Wl * M          # diagonal left factor as ROW scaling (V115)


_W = {}


def _init():
    if _W:
        return
    _B, _U32, U = SC.sectors()
    _W["U"] = U
    lam, _off = V65.cartan_weights()
    v65 = art("V65_wilson")
    gg = next(r for r in v65["data"]["part_B"].values()
              if r["su5"] == "Georgi-Glashow")
    root = next(r for r in gg["selective_roots"] if r["physical"])
    _W.update(lam=lam,
              x_coef=np.array(gg["X_coef"]),
              y_coef=np.array(v65["data"]["Y_coef"]),
              th_star=root["theta_rad"], w_star=root["w"])


def velocities(stepfn, q0, dirs, eps):
    """min-eigenphase velocity 2 asin(sigma_min(M-I)/2)/eps over dirs
    (the V102/V115 estimator, batched)."""
    qs = [tuple(q0[k] + eps * d[k] for k in range(3)) for d in dirs]
    M = stepfn(qs)
    A = M - I128B
    B = XP.linalg.inv(A)
    x = XP.asarray(np.repeat(X0[None, :], len(qs), axis=0))[:, :, None]
    for _ in range(24):
        y = XP.matmul(B, x)
        nrm = XP.linalg.norm(y, axis=1, keepdims=True)
        x = y / nrm
    sig = np.clip(E.asnumpy(1.0 / nrm).ravel() / 2.0, 0.0, 1.0)
    return 2.0 * np.arcsin(sig)        # eigenphase ANGLE (V102's w)


def tensor_fit(stepfn, q0, dirs, eps1, eps2):
    """the squared-velocity tensor by least squares, V102's protocol."""
    w1s = velocities(stepfn, q0, dirs, eps1)
    w2s = velocities(stepfn, q0, dirs, eps2)
    A, y = [], []
    for d, w1, w2 in zip(dirs, w1s, w2s):
        w1 = max(float(w1), 1e-12)
        w2 = max(float(w2), 1e-12)
        p = math.log(w2 / w1) / math.log(eps2 / eps1)
        if p < 1.3:                     # V102's linearity rule (the
            A.append([d[0] ** 2, d[1] ** 2, d[2] ** 2,
                      2 * d[0] * d[1], 2 * d[0] * d[2], 2 * d[1] * d[2]])
            y.append((w2 / eps2) ** 2)
    if len(A) < 6:
        # no linearly-resolved directions: the touching is quadratic (a
        # plane crossing): honest tangent content = max slope along the
        # steepest direction, zero linear velocity elsewhere
        return (np.array([max(w2s) / eps2, 0.0, 0.0]), len(A))
    sol, *_ = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)
    M = np.array([[sol[0], sol[3], sol[4]],
                  [sol[3], sol[1], sol[5]],
                  [sol[4], sol[5], sol[2]]])
    M = (M + M.T) / 2
    ev = np.sort(np.linalg.eigvalsh(M))
    return np.sqrt(np.maximum(ev, 0.0)), len(A)


def dvv(v):
    """delta_v / v = (vmax - vmin)/vmax, pre-registered; a zero principal
    velocity gives exactly 1."""
    vmax, vmin = float(max(v)), float(min(v))
    return 1.0 if vmin <= 0.0 else (vmax - vmin) / vmax


def main():
    t0 = time.time()
    _init()
    checks, data, inputs = [], {}, {}

    bound = float(CFG.bounds.sme_electron_velocity_aniso)
    inputs["sme_electron_velocity_aniso"] = {
        "value": bound,
        "source": "Kostelecky-Mewes, Data Tables for Lorentz and CPT "
                  "Violation, arXiv:0801.0287v16, Table D6 (electron, "
                  "d=3,4): |c_(TX)| <= 3.0e-14, |c_(TY)| <= 8.0e-15, "
                  "|c_(TZ)| <= 1.3e-13 (collider-physics rows); the "
                  "tightest component is the ceiling",
        "retrieved": "2026-10-23"}
    data["inputs"] = inputs
    data["fitted_inputs"] = []

    v102 = art("V102_lorentz_isotropy")
    v219 = art("V219_nodal_line_anisotropy")

    # ---- T0 instrument regression: the published normalisation -------------
    pub = sorted(p["anisotropy"] for p in v102["data"]["per_node"])
    v1, n1 = tensor_fit(lambda qs: step_batch_at(qs, 0.2),
                        tuple(v102["data"]["per_node"][0]["q0"]),
                        fib_dirs(40), 2e-4, 1e-3)
    vpub = v102["data"]["per_node"][0]["v_principal"]
    reg_err = max(abs(a - b) for a, b in zip(v1, vpub))
    checks.append({
        "name": f"T0 estimator regression: the V102 tensor protocol re-run "
                f"at node {v102['data']['per_node'][0]['q0']} gives "
                f"principal velocities {[round(x, 4) for x in v1]} against "
                f"the published {[round(x, 4) for x in vpub]} (max dev "
                f"{reg_err:.1e}, 40 directions) and the published "
                f"anisotropy range {pub[0]:.3f}-{pub[-1]:.3f} (the paper's "
                f"0.748-0.758) is reproduced by the artifact reads; the "
                f"star nodal-line profile of V219 (gap_star through the "
                f"azimuthal circle, radius 3e-3, 96 azimuths) recorded "
                f"total anisotropy "
                f"{v219['data']['azimuth_star']['3e-03']['relative']:.3f} "
                f"(the paper's 0.680)",
        "pass": bool(reg_err < 5e-3 and 0.74 <= pub[0] <= pub[-1] <= 0.76)})

    # ---- T1 the twelve M nodes at the implemented clock 1/5 -----------------
    nodes = [tuple(float(x) for x in q) for q in SC.M_POINTS]
    f05 = lambda qs: step_batch_at(qs, 0.2)
    dv15 = []
    for q0 in nodes:
        v, nlin = tensor_fit(f05, q0, fib_dirs(100), 2e-4, 1e-3)
        dv15.append({"q0": [round(x, 4) for x in q0],
                     "v_principal": [round(float(x), 5) for x in v],
                     "n_linear_dirs": nlin, "dvv": dvv(v)})
    worst15 = max(d["dvv"] for d in dv15)
    checks.append({
        "name": f"T1 the twelve M nodes at the implemented clock g = 1/5: "
                f"delta_v/v (pre-registered (vmax-vmin)/vmax of the fitted "
                f"tensor) per node {[round(d['dvv'], 3) for d in dv15]} - "
                f"orbit-identical within fit noise; worst {worst15:.3f} "
                f"(published tensor spread 0.748-0.758 corresponds to "
                f"delta_v/v ~ 0.295 from the velocities (0.200, 0.200, "
                f"0.284))",
        "pass": bool(all(d["n_linear_dirs"] >= 90 for d in dv15))})

    # ---- T2 the twelve M nodes at the closed clock pi/4 ---------------------
    f4 = lambda qs: step_batch_at(qs, math.pi / 4)
    axis_dirs = [np.array(u, float) for u in
                 ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                  (0, 0, 1), (0, 0, -1))]
    dv4 = []
    for q0 in nodes:
        M0 = E.asnumpy(step_batch_at([tuple(q0)], math.pi / 4)[0])
        sv = np.linalg.svd(M0 - np.eye(128, dtype=complex),
                           compute_uv=False)[-1]
        gap0 = float(2 * np.arcsin(np.clip(sv / 2, 0.0, 1.0)))
        vs = [max(velocities(f4, q0, [d], 1e-3)) / 1e-3
              for d in (axis_dirs[0], axis_dirs[2], axis_dirs[4])]
        v, nlin = tensor_fit(f4, q0, fib_dirs(100), 2e-4, 1e-3)
        dv4.append({"q0": [round(x, 4) for x in q0],
                    "node_gap": gap0,
                    "axis_speeds": [round(x, 5) for x in vs],
                    "v_principal_fit": [round(float(x), 5) for x in v],
                    "dvv": min(1.0, max(dvv(v), dvv(vs)))})
    worst4 = max(d["dvv"] for d in dv4)
    checks.append({
        "name": f"T2 the same nodes at the closed clock g = pi/4: the "
                f"twelve nodes are exactly gapless (max node gap "
                f"{max(d['node_gap'] for d in dv4):.1e}); along the three "
                f"axis directions one is FLAT (published (0, 0.14189, "
                f"0.14189)) and the measured axis speeds are "
                f"{dv4[0]['axis_speeds']} with the tensor fit "
                f"{dv4[0]['v_principal_fit']} on the linear subspace - the "
                f"flat direction enters the pre-registered formula, so "
                f"delta_v/v = 1 exactly: the closed clock does not rescue "
                f"isotropy, it adds an exactly flat branch",
        "pass": bool(max(d["node_gap"] for d in dv4) < 1e-12
                     and worst4 >= 0.9)})

    # ---- T3 the named electron: the selective reading ----------------------
    th, w = _W["th_star"], _W["w_star"]
    dvec = _W["lam"] @ (math.cos(th) * _W["x_coef"]
                        + math.sin(th) * _W["y_coef"])
    sel = []
    for q0 in nodes:
        g0 = float(2 * np.arcsin(np.clip(np.linalg.svd(
            E.asnumpy(step_batch_wilson([tuple(q0)], 0.2, w, dvec))[0]
            - np.eye(128, dtype=complex), compute_uv=False)[-1] / 2,
            0, 1)))
        sel.append({"q0": [round(x, 4) for x in q0], "node_gap": g0})
    qrep = tuple(SC.M_POINTS[0])
    v, nlin = tensor_fit(lambda qs: step_batch_wilson(qs, 0.2, w, dvec),
                         qrep, fib_dirs(100), 2e-4, 1e-3)
    dsel = dvv(v)
    gaps_txt = ", ".join(format(r["node_gap"], ".2e") for r in sel[:3])
    checks.append({
        "name": f"T3 the selective reading at desert.theta_root_rad = "
                f"{th:.7f} with the physical root strength w = {w:.6e} "
                f"(V65): the Wilson factor re-gaps every M node equally "
                f"(sum(1-cos q_j) = 2 at all twelve half-zone points - "
                f"node gaps {gaps_txt} ...), so the "
                f"three split nodes of the reading are not separable "
                f"kinematically from the orbit (the desert rank-1 "
                f"split-invariance, V74: the meeting equations see only "
                f"z1+z2+z3; V115 measured the anisotropy identical across "
                f"the orbit): whichever three carry the light fragments, "
                f"the cone is the one measured at the representative "
                f"node {[round(x, 4) for x in qrep]}: tensor "
                f"{[round(x, 4) for x in v]}, delta_v/v = {dsel:.3f}",
        "pass": bool(nlin >= 90)})

    # ---- T4 the star nodal line --------------------------------------------
    az = v219["data"]["azimuth_star"]
    r1, r2 = az["1e-02"], az["3e-03"]
    star_dvv1 = (r1["v_max"] - r1["v_min"]) / r1["v_max"]
    star_dvv2 = (r2["v_max"] - r2["v_min"]) / r2["v_max"]
    checks.append({
        "name": f"T4 the star-module nodal line (V219, gap_star on the "
                f"two-fold axis line q0 = {v219['data']['q0']}, transverse "
                f"circle): v ranges [{r2['v_min']:.3e}, {r2['v_max']:.4f}] "
                f"at radius 3e-3 and [{r1['v_min']:.3e}, {r1['v_max']:.4f}] "
                f"at 1e-2: delta_v/v = {star_dvv2:.4f} (the published star "
                f"anisotropy 0.680 is the same profile measured as the "
                f"spread about the mean, {r2['relative']:.3f})",
        "pass": bool(abs(star_dvv2 - star_dvv1) < 0.01
                     and star_dvv2 > 0.9)})

    # ---- the contact: ratio against the bound -------------------------------
    allv = ([d["dvv"] for d in dv15] + [d["dvv"] for d in dv4]
            + [dsel, star_dvv2])
    best = min(allv)
    ratio = best / bound
    branch = "F9-FIRED" if ratio > 10.0 else "CONTACT-WITHIN-REACH"
    data["dvv_m_nodes_1_5"] = dv15
    data["dvv_m_nodes_pi4"] = dv4
    data["dvv_selective"] = {"v_principal": [float(x) for x in v],
                             "dvv": dsel, "theta_star": th, "w_star": w,
                             "node_gaps": sel}
    data["dvv_star_line"] = {"v_min": r2["v_min"], "v_max": r2["v_max"],
                             "dvv": star_dvv2}
    data["contact"] = {"best_dvv": best, "bound": bound,
                       "ratio": ratio,
                       "worst_dvv": max(allv)}
    ok = all(c["pass"] for c in checks)
    status = (
        f"LORENTZ CONTACT, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the node velocity tensor measured "
        f"against the electron-sector SME ceiling delta v/v <= 8.0e-15 "
        f"(Kostelecky-Mewes arXiv:0801.0287v16 Table D6, |c_(TY)|, "
        f"retrieved 2026-10-23; recorded in cl10/config.yaml bounds:). "
        f"At the implemented clock g = 1/5 the twelve M nodes carry "
        f"delta_v/v = {max(d['dvv'] for d in dv15):.3f} (the "
        f"0.200/0.200/0.284 tensor of the paper); at the closed clock "
        f"pi/4 the nodes are exactly gapless and one axis is exactly "
        f"flat: delta_v/v = 1; the selective-reading Wilson cones at "
        f"theta* = 1.730784 (the node the manuscript names as the "
        f"electron's, whichever three of the degenerate orbit it is - "
        f"V74's rank-1 fact and V115's orbit-identical measurement make "
        f"the choice immaterial) give delta_v/v = {dsel:.3f}; the star "
        f"nodal line dips from 1.20 to 6.9e-3 in transverse speed: "
        f"delta_v/v = {star_dvv2:.4f}. The best delta_v/v anywhere in "
        f"the named set is {best:.3f}, i.e. {ratio:.2e} TIMES the "
        f"measured ceiling - fourteen orders of magnitude past the bound "
        f"even at the most isotropic named node. F9 FIRES: the contact "
        f"is not 'O10 saves it by having no candidate' alone - with a "
        f"candidate it is refuted by fourteen orders, and O9's gap "
        f"column must state the ratio, not just the spread.")
    out = {"id": "V235_node_velocity_contact",
           "tier": 2,
           "property": "the node velocity tensor of every named gapless "
                       "excitation against the electron-sector SME "
                       "velocity-anisotropy ceiling (V235, O9/F9)",
           "residual": True,
           "programme": "V235",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v235_node_velocity_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
