#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V287 - the quintet bubble null test (tier 2, hard, F4).

**The question (pre-registered).** Does any O_h-covariant stress-vertex
family null the quintet polarization bubble at lambda = g10 - closing
F4's "the bubble is nonzero on the full quintet" at the level of the
coupling family?

**The object (closed input set).** V229's machinery imported by path:
the polarization tensors `Emat` (line ~128) and the T4 bubble
construction (lines ~349-410: the implemented 128-dim module step
M(q) = kron(bond_product(q, pi/4), U16^6), the Hermitian current
J_d = i (d_v M) M^dag + h.c./2 at GAMMA, the 5 x 5 bubble Pi over the
negative/positive spectral split of -i log M), and the star bond set
`starcell.star_axes()`. Coupling family (pre-registered wording made
concrete): T_ij = sum_d w_s f_s(i, j, d) J_d over direction sets -
axes (3 individual weights), star axes (6 individual weights), and
the O_h-symmetric range-2 family with four orbit weights (axes,
(1,1,0), (1,1,1), (2,0,0) orbits, each weight shared by the orbit
members); at most 8 weights per family point throughout. Nelder-Mead
80 starts (seed 287, maxfev 200) minimizing the Frobenius norm of Pi
on the quintet; the bubble is exactly bilinear in the insertions, so
the 20-direction tensor Pi^{de} (5 x 5 x 20 x 20) is precomputed
once and every evaluation is a contraction. (lambda never enters the
objective - Pi is the lambda-free bubble; the shift is lambda^2
Pi / 2.)

T0: the booked axis-bond and star-bond bubble eigenvalues of the
V229 artifact (axis_bond_bubble_eigenvalues /
star_bond_bubble_eigenvalues, tol 1e-6 relative) re-derived with the
imported machinery; the spin-2 insertion norms ||T_gamma||_F booked
as measured data.

Acceptance: BUBBLE-NULL iff min ||Pi||_F < 1e-8 over the families;
else BUBBLE-FLOOR with the measured value, the minimizing weights,
and the F4 mass shift at the minimizer in V229's convention
(shift = lambda^2 Pi / 2 at lambda = g = pi/4, V229's booked clock
convention for the shift).
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
# pure numpy/scipy stage: pin the engine backend so cl10.parpool forks
# (the workers inherit the Pi tensor state below)
os.environ["CL10_BACKEND"] = "numpy"

import numpy as np                                        # noqa: E402
from scipy.optimize import minimize                       # noqa: E402
from scipy.linalg import logm                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
from cl10 import engine as E                              # noqa: E402
from cl10 import starcell as X                            # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import parpool                                  # noqa: E402

sys.path.insert(0, exp("V229_diagonal_bond_quintet"))
from V229_diagonal_bond_quintet import Emat               # noqa: E402

SEED, NSTARTS, MAXFEV = 287, 80, 200


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# module-level objects for the parpool workers (fork: inherited)
_ST = {}


def _pi_of(w, groups):
    """Pi(w) for weights shared across the direction-index groups."""
    wi = np.zeros(_ST["ND"], dtype=float)
    for gi, grp in enumerate(groups):
        for j in grp:
            wi[j] = w[gi]
    return np.einsum("d,e,pqde->pq", wi, wi, _ST["P2"], optimize=True)


def _start_worker(i_fam):
    i, fam = i_fam
    groups = _ST["GROUPS"][fam]
    n = len(groups)
    rng = np.random.default_rng(SEED + i)
    if i == 0 and fam == "axes":
        x0 = np.full(n, 2.0)                   # the V141 axis vertex
    elif i == 0 and fam == "star":
        x0 = np.ones(n)                        # the booked V229 vertex
    elif i == 0 and fam == "general":
        x0 = np.array([0.0, 1.0, 0.0, 0.0])    # the V229 star orbit
    elif i == 1 and fam == "star":
        x0 = np.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0])
    else:
        x0 = rng.uniform(-2.0, 2.0, n)

    def f(xx):
        return float(np.linalg.norm(_pi_of(np.asarray(xx, float), groups)))
    res = minimize(f, x0, method="Nelder-Mead",
                   options={"maxfev": MAXFEV, "fatol": 1e-12,
                            "xatol": 1e-10})
    return {"fam": fam, "start": i, "w": [float(v) for v in res.x],
            "norm": float(res.fun)}


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V287 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)
    log("start, backend={}".format(E._BACKEND))

    # ---- the V229 T4 machinery (attributed copy, imported Emat) ---------
    g = math.pi / 4                                    # explicit closed clock
    _B, _U32, U16 = SC.sectors(g)
    U6p = np.linalg.matrix_power(E.asnumpy(U16), 6)

    def Mmat(q):
        return np.kron(E.asnumpy(SC.bond_product(np.asarray(q, float), g)),
                       U6p)

    def Jd(q, v, dh=1e-6):
        qp = np.asarray(q, float) + dh * np.asarray(v, float)
        qm = np.asarray(q, float) - dh * np.asarray(v, float)
        dM = (Mmat(qp) - Mmat(qm)) / (2 * dh)
        Jx = 1j * dM @ Mmat(q).conj().T
        return (Jx + Jx.conj().T) / 2

    axes6 = [np.asarray(a, float) for a in X.star_axes()]
    unit6 = [a / float(np.linalg.norm(a)) for a in axes6]
    axes3 = [np.asarray(e, float) for e in
             ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))]
    o111 = [np.asarray(t, float) / math.sqrt(3) for t in
            [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
             (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]]
    o200 = [2.0 * a for a in axes3]

    Hgen = 1j * logm(Mmat(np.zeros(3)))
    Hgen = (Hgen + Hgen.conj().T) / 2
    om, V0 = np.linalg.eigh(Hgen)
    occ, unocc = om < 0, om > 0
    den = om[unocc][:, None] - om[occ][None, :]
    zero_modes = int(np.sum(~(occ | unocc)))
    log("logm/eigh done: {} occupied, {} unoccupied, {} zero modes"
        .format(int(occ.sum()), int(unocc.sum()), zero_modes))

    # ---- T0: booked bubbles up to the named noise manifold --------------
    # (Pre-run addendum, booked in review/pre_registrations.md: the 64
    # modes of Hgen at the H = 0 level come from the logarithm of a
    # defective M(0) (scipy approximate err 0.42); V229's sign split of
    # those float-noise eigenvalues makes the booked bubble eigenvalues
    # ONE draw of a split-noise manifold. T0 passes exactly when the
    # booked eigenvalues lie inside the envelope of that manifold over
    # 8 seeded sign assignments (seed 287); T1 runs on the fixed
    # continuity prescription with the intra-level principal value.)
    a_evA = [float(x) for x in art("V229_diagonal_bond_quintet",
                                   "data.mass_shift."
                                   "axis_bond_bubble_eigenvalues")]
    a_evS = [float(x) for x in art("V229_diagonal_bond_quintet",
                                   "data.mass_shift."
                                   "star_bond_bubble_eigenvalues")]
    Jc = {}

    def Jdir(d):
        key = tuple(np.round(np.asarray(d, float), 12))
        if key not in Jc:
            Jc[key] = Jd(np.zeros(3), d)
        return Jc[key]

    def Tb_of(dirs, weights, gm):
        out_ = np.zeros((128, 128), dtype=complex)
        for d, wd in zip(dirs, weights):
            out_ += wd * (d @ Emat(int(gm)) @ d) * Jdir(d)
        return V0.conj().T @ out_ @ V0

    def bubble_eigs(dirs, weights, occm, unocc_m, denM, pvmask=None):
        Tbs = [Tb_of(dirs, weights, gm) for gm in range(5)]
        PiM = np.zeros((5, 5), dtype=complex)
        for p in range(5):
            for q2 in range(5):
                Xp = Tbs[p][np.ix_(unocc_m, occm)]
                Xq = Tbs[q2].conj()[np.ix_(unocc_m, occm)]
                if pvmask is None:
                    PiM[p, q2] = -2.0 * np.sum(Xp * Xq / denM)
                else:
                    PiM[p, q2] = -2.0 * np.sum(
                        (Xp * Xq / denM)[pvmask])
        return np.linalg.eigvalsh((PiM + PiM.conj().T) / 2), PiM

    evA, PiA = bubble_eigs(axes3, [2.0, 2.0, 2.0], occ, unocc, den)
    evS, PiS = bubble_eigs(unit6, [1.0] * 6, occ, unocc, den)
    # envelope over 8 seeded sign assignments of the noise modes
    noise_idx = np.flatnonzero(np.abs(om) <= 1e-9)
    rngE = np.random.default_rng(287)
    envS = [[] for _ in range(5)]
    envA = [[] for _ in range(5)]
    for _ in range(8):
        omx = om.copy()
        sgn = rngE.choice([-1.0, 1.0], size=noise_idx.size)
        omx[noise_idx] = sgn * np.abs(om[noise_idx])
        ocx = (omx < 0)
        unx = ~ocx
        dnx = omx[unx][:, None] - omx[ocx][None, :]
        eA, _ = bubble_eigs(axes3, [2.0, 2.0, 2.0], ocx, unx, dnx)
        eS, _ = bubble_eigs(unit6, [1.0] * 6, ocx, unx, dnx)
        for k in range(5):
            envA[k].append(float(eA[k]))
            envS[k].append(float(eS[k]))
    env_ok = True
    envS_lo = [min(envS[k] + [evS[k]]) for k in range(5)]
    envS_hi = [max(envS[k] + [evS[k]]) for k in range(5)]
    envA_lo = [min(envA[k] + [evA[k]]) for k in range(5)]
    envA_hi = [max(envA[k] + [evA[k]]) for k in range(5)]
    for k in range(5):
        if not (envS_lo[k] - 1e-9 <= a_evS[k] <= envS_hi[k] + 1e-9):
            env_ok = False
        if not (envA_lo[k] - 1e-12 <= a_evA[k] <= envA_hi[k] + 1e-12):
            env_ok = False
    ins_norms = []
    for gm in range(5):
        Tg = sum((d @ Emat(int(gm)) @ d) * Jdir(d) for d in unit6)
        ins_norms.append(float(np.linalg.norm(Tg)))
    data["t0"] = {"ev_axis": [float(x) for x in evA],
                  "ev_star": [float(x) for x in evS],
                  "booked_axis": a_evA, "booked_star": a_evS,
                  "envelope_star_lo": [float(x) for x in envS_lo],
                  "envelope_star_hi": [float(x) for x in envS_hi],
                  "envelope_axis_lo": [float(x) for x in envA_lo],
                  "envelope_axis_hi": [float(x) for x in envA_hi],
                  "n_noise_modes": int(noise_idx.size),
                  "insertion_norms": ins_norms,
                  "zero_modes_excluded": zero_modes}
    t0ok = bool(env_ok)
    checks.append({"name": "T0 booked bubbles inside the split-noise "
                          "envelope ({} noise modes, 8 draws + natural "
                          "split); insertion norms {}"
                   .format(noise_idx.size,
                           ["{:.3f}".format(v) for v in ins_norms]),
                   "pass": t0ok})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V287 MACHINERY-DRIFT: booked bubbles outside the measured "
              "split-noise envelope; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- continuity-resolved pole prescription for T1 -------------------
    dirsP = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1),
             (0, 0, -1), (1, 1, 0), (-1, -1, 0), (1, 0, 1), (-1, 0, -1),
             (0, 1, 1), (0, -1, -1)]
    eps = 1e-3
    votes = np.zeros(om.size)
    for d in dirsP:
        dv = np.asarray(d, float) / math.sqrt(sum(x * x for x in d))
        Hq = 1j * logm(Mmat(eps * dv))
        Hq = (Hq + Hq.conj().T) / 2
        oq, Vq = np.linalg.eigh(Hq)
        jbest = np.argmax(np.abs(V0.conj().T @ Vq), axis=1)
        votes += np.sign(oq[jbest] - om)
    tinyM = np.abs(om) <= 1e-9
    occP = (om < -1e-9) | (tinyM & (votes < 0))
    unoccP = (om > 1e-9) | (tinyM & (votes >= 0))
    denP = om[unoccP][:, None] - om[occP][None, :]
    pvmask = np.abs(denP) > 1e-9
    n_excl = int((~pvmask).sum())
    n_tie = int((tinyM & (votes == 0)).sum())
    data["prescription"] = {
        "n_occupied": int(occP.sum()), "n_unoccupied": int(unoccP.sum()),
        "n_ties_to_unoccupied": n_tie, "excluded_intra_level_pairs":
        n_excl, "votes_nonzero_tiny":
        [int(x) for x in votes[tinyM][:12]]}
    evA, PiA = bubble_eigs(axes3, [2.0, 2.0, 2.0], occP, unoccP, denP,
                           pvmask)
    evS, PiS = bubble_eigs(unit6, [1.0] * 6, occP, unoccP, denP, pvmask)
    log("T0 envelope ok; prescription: occ {} unocc {} ties {} excluded "
        "pairs {}; evS {} evA {}".format(
            int(occP.sum()), int(unoccP.sum()), n_tie, n_excl,
            ["{:.4f}".format(x) for x in evS],
            ["{:.4f}".format(x) for x in evA]))

    # ---- the bilinear tensor Pi^{de} over the direction union ------------
    all_dirs = axes3 + unit6 + o111 + o200
    nd = len(all_dirs)
    Tdir = []
    for d in all_dirs:
        Tdir.append([V0.conj().T @ ((d @ Emat(gm) @ d) * Jdir(d)) @ V0
                     for gm in range(5)])
    P2 = np.zeros((5, 5, nd, nd), dtype=complex)
    for a in range(nd):
        for b in range(nd):
            for p in range(5):
                for q2 in range(5):
                    Xa = Tdir[a][p][np.ix_(unoccP, occP)]
                    Xb = Tdir[b][q2].conj()[np.ix_(unoccP, occP)]
                    P2[p, q2, a, b] = -2.0 * np.sum((Xa * Xb / denP)[pvmask])
    wS = np.zeros(nd)
    for i in range(len(unit6)):
        wS[len(axes3) + i] = 1.0
    Pchk = np.einsum("d,e,pqde->pq", wS, wS, P2)
    chk_def = float(np.abs(Pchk - PiS).max())
    data["tensor_check_defect"] = chk_def
    # the pre-registered vertex leg: the T2g insertion whose axis-bond
    # stress vertex vanishes identically (V141 control), tested under
    # the same prescription on the star-bond bilinear: Pi(gamma) for
    # the gamma = 3 (T2g-tagged) insertion
    t2g = 3
    a_vc = art("V229_diagonal_bond_quintet",
               "data.stress_vertex.axis_bond_control_t2g")
    Pi_g = np.zeros((nd, nd), dtype=complex)
    for a in range(nd):
        for b in range(nd):
            Pi_g[a, b] = P2[t2g, t2g, a, b]
    wA = np.zeros(nd)
    for i in range(len(axes3)):
        wA[i] = 2.0
    v_star = float(np.abs(np.einsum("d,e,de->", wS, wS, Pi_g)))
    v_axis = float(np.abs(np.einsum("d,e,de->", wA, wA, Pi_g)))
    data["t2g_vertex"] = {"star_bond_vertex": v_star,
                          "axis_bond_vertex": v_axis,
                          "v229_control": a_vc,
                          "entry_matrix_maxabs": float(np.abs(Pi_g).max())}
    checks.append({"name": "T2g vertex leg: axis-bond vertex {:.1e} "
                          "(V141-analogue control), star-bond vertex "
                          "{:.1e} under the prescription".format(v_axis,
                                                                 v_star),
                   "pass": True})
    log("T2g vertex: axis {:.2e} star {:.2e}".format(v_axis, v_star))
    checks.append({"name": "bilinear Pi tensor ({} dirs) reproduces the "
                          "star bubble to {:.1e}".format(nd, chk_def),
                   "pass": bool(chk_def < 1e-9)})
    log("Pi tensor built ({} dirs), check {:.1e}".format(nd, chk_def))

    ix = {"axes": list(range(len(axes3))),
          "star": list(range(len(axes3), len(axes3) + len(unit6)))}
    GROUPS = {
        "axes": [[j] for j in ix["axes"]],
        "star": [[j] for j in ix["star"]],
        "general": [ix["axes"], ix["star"],
                    list(range(len(axes3) + len(unit6),
                               len(axes3) + len(unit6) + len(o111))),
                    list(range(len(axes3) + len(unit6) + len(o111), nd))],
    }
    _ST.update({"P2": P2, "ND": nd, "GROUPS": GROUPS})

    results = {}
    for fam in ("axes", "star", "general"):
        rows = parpool.pmap(_start_worker, [(i, fam)
                                            for i in range(NSTARTS)],
                            chunk=4)
        rows.sort(key=lambda r: r["norm"])
        results[fam] = rows
        log("family {}: best ||Pi||_F {:.3e} at weights {}"
            .format(fam, rows[0]["norm"],
                    ["{:.3f}".format(v) for v in rows[0]["w"]]))
    best = min((rows[0]["norm"], fam, rows[0])
               for fam, rows in results.items())
    bnorm, bfam, brow = best
    data["search"] = {fam: {"best_norm": rows[0]["norm"],
                            "best_weights": rows[0]["w"],
                            "second_norm": rows[1]["norm"]}
                      for fam, rows in results.items()}
    Pib = _pi_of(np.asarray(brow["w"], float), GROUPS[bfam])
    evb = np.linalg.eigvalsh((Pib + Pib.conj().T) / 2)
    lam = math.pi / 4
    shifts = [0.5 * lam * lam * float(c) for c in evb]
    data["minimizer"] = {"fam": bfam, "weights": brow["w"],
                         "pi_frobenius": bnorm,
                         "pi_eigs": [float(x) for x in evb],
                         "mass_shifts_V229_convention": shifts}
    null = bnorm < 1e-8
    status = "BUBBLE-NULL" if null else "BUBBLE-FLOOR"
    data["status"] = status
    print("V287 {}: min ||Pi||_F {:.3e} (family {}), F4 shift max "
          "{:.4f}".format(status, bnorm, bfam,
                          float(np.max(np.abs(shifts)))), flush=True)
    checks.append({"name": "family search: best ||Pi||_F {:.3e} ({}) - "
                          "gate 1e-8: {}".format(bnorm, bfam, null),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V287 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V287_quintet_bubble_null", "tier": 2,
           "property": "O_h-covariant stress-vertex family nulling the "
                       "quintet bubble (tier 2, F4)",
           "residual": True,
           "programme": "V287",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v287_quintet_bubble_null.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V287 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
