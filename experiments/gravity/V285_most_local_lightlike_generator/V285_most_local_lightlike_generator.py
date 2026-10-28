#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V285 - the most local lightlike generator (tier 2, hard, F8/O9).

**The question (pre-registered).** Within the 53-dimensional family
of exact-lightcone generators (V208: the Clifford-regularity
Jacobian at the committed generator has rank 107 in the 160-dim
covariant space, family dimension 53), is the committed generator
selected by locality - or is it a boundary accident of the family?

Object (V208's construction region re-derived in-script by path,
attributed copies): the bilinear module 6 (x) Sym^2 (36-dim), the
O(48) representation Dm / Rv, the committed triple from
V208's cone_generator.json, the covariant 160-space basis
(cov_project sampling, as V208 builds it), the linearised
Clifford-regularity Jacobian (rank 107, nullity 53 - re-measured),
and V208's hopping-kernel profile machinery (N=15 FFT grid).
Objective: the tail functional T(delta) = sum_{r>=2} kernel_r^2 of
the kernel profile at the family point.

Family parameterization (pre-reg wording made concrete, booked):
tangent coordinates delta in R^53 along the Jacobian nullspace
directions v_a; the point H + sum delta_a v_a is retracted onto the
exact Clifford-regularity manifold by fixed-Jacobian Gauss-Newton in
the covariant 160 coordinates (residual: the six V208 def1
anticommutator identities AND the exact-cone scalar - the top
eigenvalue of the common S = H_a^2 equals 1, which is what makes the
band slopes exactly +-1/0, i.e. exact-cone membership maintained to
machine precision after each retract); the retraction is the log/
exponential correction of the tangent step at quadratic order.
Minimize T over the 53 coordinates by Nelder-Mead, 50 starts
(seed 285, maxfev 200 each), box |delta_a| <= sqrt(i3(H,H)/3) per
coordinate (booked choice), workers via cl10.parpool.

Measured: T0 anchors (V208 artefact): kernel profile entries
(0.91086..., 0.10705..., 0.034139, 0.0054357, 7.8786e-4, 7.1579e-5)
tol 1e-9 relative; hop projection 0.0028977782769928693 and
dressing projection 0.0092074659831694 (tol 1e-9 relative); Jacobian
rank 107 and family dim 53 re-derived. Then: T_committed, the
minimizer T*, its hop/dressing projections and i3-distance from the
committed generator (relative), the 53-gradient norm at T*, the
number of starts whose solution touches the box, the spread of the
best values.

Acceptance: T0 fail -> MACHINERY-DRIFT. CONE-LOCALLY-SELECTED iff
an interior strict minimum with T* < T_committed and hop projection
> 0.0092; else CONE-SELECTION-ABSENT with the boundary data booked.
"""
import itertools
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
# the whole stage is numpy/scipy; pin the engine backend to numpy so
# cl10.parpool forks (the workers inherit the retraction state below)
os.environ["CL10_BACKEND"] = "numpy"

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
from cl10 import parpool                                  # noqa: E402

from scipy.optimize import minimize                       # noqa: E402

SEED, NSTARTS, MAXFEV = 285, 50, 200


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# ---------------------------------------------------------------------------
# V208 machinery (attributed copies of the construction region:
# build_Oh/module/bilinear reps lines 92-166, velocity triples lines
# 191-213, kernel-profile lines 351-375, covariant basis + Jacobian
# lines 405-455 - same conventions, no numbers retyped)
# ---------------------------------------------------------------------------
def build_Oh():
    OHS = [(p, s) for p in itertools.permutations(range(3))
           for s in itertools.product((1, -1), repeat=3)]
    Rv = []
    for p, s in OHS:
        R = np.zeros((3, 3))
        for i in range(3):
            R[i, p[i]] = s[i]
        Rv.append(R)
    return OHS, Rv


IDX = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
DIRS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def module_reps():
    OHS, Rv = build_Oh()
    SB = []
    for k, l in IDX:
        S = np.zeros((3, 3))
        S[k, l] += 1.0
        if l != k:
            S[l, k] += 1.0
        SB.append(S / math.sqrt(2.0 if l != k else 1.0))
    Dp, Db = [], []
    for R in Rv:
        Mm = np.zeros((6, 6))
        for a in range(6):
            T = R @ SB[a] @ R.T
            for b in range(6):
                Mm[b, a] = np.trace(SB[b].T @ T)
        Dp.append(Mm)
        Pm = np.zeros((6, 6))
        for j, d in enumerate(DIRS):
            Rd = tuple(int(round(x)) for x in R @ np.array(d))
            Pm[DIRS.index(Rd), j] = 1.0
        Db.append(Pm)
    Dm = [np.kron(Pb, Pp) for Pb, Pp in zip(Db, Dp)]
    return Rv, Dm


def velocity_ops():
    lev = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
           (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
    Bb = [np.zeros((6, 6)) for _ in range(3)]
    for k in range(3):
        for j, d in enumerate(DIRS):
            Bb[k][j, j] = d[k]
    SB = []
    for k, l in IDX:
        S = np.zeros((3, 3))
        S[k, l] += 1.0
        if l != k:
            S[l, k] += 1.0
        SB.append(S / math.sqrt(2.0 if l != k else 1.0))
    Kp = []
    for k in range(3):
        Aa = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                Aa[i, j] = -lev.get((k, i, j), 0)
        Mm = np.zeros((6, 6), dtype=complex)
        for b in range(6):
            act = 1j * (Aa @ SB[b] - SB[b] @ Aa)
            for a in range(6):
                Mm[a, b] = np.trace(SB[a].T @ act)
        Kp.append(Mm)
    Bop = [np.kron(Bb[k], np.eye(6)) for k in range(3)]
    Kop = [np.kron(np.eye(6), Kp[k]) for k in range(3)]
    return Bop, Kop


def vec3(Ts):
    return np.concatenate([t.real.reshape(-1) for t in Ts]
                          + [t.imag.reshape(-1) for t in Ts])


def unvec(vv):
    b = vv.reshape(6, 36, 36)
    return [b[k] + 1j * b[3 + k] for k in range(3)]


def sym(A, B):
    return A @ B + B @ A


def i3(A, B):
    return sum(np.trace(A[a].conj().T @ B[a]).real for a in range(3))


def triple_ops(tr):
    """Residual pieces of one triple: def1 six identities and the
    common S."""
    S3 = (sum(h @ h for h in tr) / 3.0)
    S3 = (S3 + S3.conj().T) / 2.0
    def1 = [sym(tr[0], tr[0]) - 2 * S3, sym(tr[1], tr[1]) - 2 * S3,
            sym(tr[2], tr[2]) - 2 * S3,
            sym(tr[0], tr[1]), sym(tr[0], tr[2]), sym(tr[1], tr[2])]
    return def1, S3


CONE_DIRS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1),
             (0, 1, 1), (1, 1, 1), (1, 2, 3)]


def profile(tr, N=15):
    """V208's kernel profile: M(q) = exp(i sum sin(q) H) on the N^3
    grid, kernel = ifftn, profile = max abs over Chebyshev rings."""
    Qs = 2.0 * np.pi * np.arange(N) / N
    Kk = np.zeros((N, N, N, 36, 36), dtype=complex)
    for a in range(N):
        for b in range(N):
            for c in range(N):
                nrm = math.sin(Qs[a]) * tr[0] + math.sin(Qs[b]) * tr[1] \
                    + math.sin(Qs[c]) * tr[2]
                w, U = np.linalg.eigh(nrm)
                Kk[a, b, c] = (U * np.exp(1j * w)) @ U.conj().T
    ker = np.fft.ifftn(Kk, axes=(0, 1, 2))
    prof = {}
    for i in range(N):
        for j in range(N):
            for l in range(N):
                r = max(min(i, N - i), min(j, N - j), min(l, N - l))
                prof[r] = max(prof.get(r, 0.0),
                              float(np.abs(ker[i, j, l]).max()))
    return prof


# ---------------------------------------------------------------------------
# module-level state for the parpool workers (fork inherits these)
# ---------------------------------------------------------------------------
_ST = {}


def _resid(c):
    tr = unvec(c @ _ST["Bmat"])
    def1, S3 = triple_ops(tr)
    out = []
    for z in def1:
        out.append(z.real.reshape(-1))
        out.append(z.imag.reshape(-1))
    w = np.linalg.eigvalsh(S3)
    # pin the booked slope-1 cone block: the top eight eigenvalues of
    # the common S sit exactly at 1 (the 8+4/4-branch structure that
    # T0 re-derives live), the rest of the spectrum is free
    out.append(np.array([w[-1] - 1.0, w[-8] - 1.0]))
    return np.concatenate(out)


def _retract(c0, maxit=3):
    """fixed-J Gauss-Newton onto the exact manifold."""
    c = c0.copy()
    Q, R = _ST["QR"]
    r = _resid(c)
    for _ in range(maxit):
        b = -np.linalg.lstsq(R, Q.T @ r, rcond=None)[0]
        c = c + b
        r = _resid(c)
        if np.linalg.norm(r) < 1e-10:
            break
    return c, r


def _tail_of(tr):
    prof = profile(tr)
    return sum(prof[r] ** 2 for r in range(2, 8))


def _obj(delta):
    d = np.asarray(delta, float)
    box = _ST["BOX"]
    d = np.clip(d, -box, box)
    c0 = _ST["cH"] + d @ _ST["NULL"]
    c, r = _retract(c0)
    nr = float(np.linalg.norm(r))
    if not np.isfinite(nr) or nr > 1e-8:
        return 1e3
    tr = unvec(c @ _ST["Bmat"])
    if not np.all(np.isfinite(tr)):
        return 1e3
    # cone condition after retraction: the unit-direction spectra must
    # keep exactly the booked branch structure at slope 1 (counts
    # re-derived at T0 from the live generator and the V208 artefact,
    # tol 1e-6); the rest of the 36 is non-propagating content
    for nd in CONE_DIRS:
        nn = np.sqrt(sum(x * x for x in nd))
        M = sum(x / nn * h for x, h in zip(nd, tr))
        M = (M + M.conj().T) / 2
        w = np.linalg.eigvalsh(M)
        n1 = int(np.sum(np.abs(w - 1.0) < 1e-6))
        nm1 = int(np.sum(np.abs(w + 1.0) < 1e-6))
        if n1 != _ST["BR"][0] or nm1 != _ST["BR"][1]:
            return 1e3
    try:
        return _tail_of(tr)
    except np.linalg.LinAlgError:
        return 1e3


def _start_worker(i):
    rng = np.random.default_rng(SEED + i)
    box = _ST["BOX"]
    if i == 0:
        x0 = np.zeros(53)                      # the committed generator
    else:
        x0 = rng.uniform(-box, box, 53)

    def f(xx):
        try:
            return _obj(xx)
        except (np.linalg.LinAlgError, FloatingPointError):
            return 1e3
    res = minimize(f, x0, method="Nelder-Mead",
                   options={"maxfev": MAXFEV, "fatol": 1e-10,
                            "xatol": 1e-8})
    x = np.clip(res.x, -box, box)
    return {"start": i, "x": [float(v) for v in x], "T": float(res.fun),
            "fev": int(res.nfev)}


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V285 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    Rv, Dm = module_reps()
    Bop, Kop = velocity_ops()
    gen = json.load(open(os.path.join(
        exp("V208_metric_bilinear_cone"), "cone_generator.json")))
    H = [np.array(gen["A"][k]) + 1j * np.array(gen["S"][k])
         for k in range(3)]
    PROP = [i for i, R in enumerate(Rv) if np.linalg.det(R) > 0]

    def cov_project(Ts):
        outp = [np.zeros((36, 36), dtype=complex) for _ in range(3)]
        for i in PROP:
            D, R = Dm[i], Rv[i]
            DT = [D @ Tb @ D.conj().T for Tb in Ts]
            for a in range(3):
                for b in range(3):
                    outp[a] += R[a, b] * DT[b]
        return [o / 24. for o in outp]

    # covariant basis (V208 sampling convention)
    srng = np.random.default_rng(20261021)

    def rherm():
        A = (srng.standard_normal((36, 36))
             + 1j * srng.standard_normal((36, 36)))
        return (A + A.conj().T) / 2.0
    rows, rank = [], 0
    for _ in range(4 * 160 + 60):
        v = vec3(cov_project([rherm(), rherm(), rherm()]))
        r = int(np.linalg.matrix_rank(np.array(rows + [v]), tol=1e-9))
        if r > rank:
            rows.append(v)
            rank = r
        if rank == 160:
            break
    _, sv, Vt = np.linalg.svd(np.array(rows), full_matrices=False)
    Bmat = Vt[sv > 1e-9]                        # (160, 1152)
    log("covariant basis {} x {} (cov_dim {})".format(*Bmat.shape,
                                                      Bmat.shape[0]))

    def linearised(dH):
        d1, d2, d3 = dH
        return [sym(H[0], d2) + sym(d1, H[1]),
                sym(H[0], d3) + sym(d1, H[2]),
                sym(H[1], d3) + sym(d2, H[2]),
                sym(H[0], d1) - sym(H[1], d2),
                sym(H[0], d1) - sym(H[2], d3)]
    cols = []
    for v in Bmat:
        ms = linearised(unvec(v))
        cols.append(np.concatenate([m.real.reshape(-1) for m in ms]
                                   + [m.imag.reshape(-1) for m in ms]))
    Jm = np.array(cols).T                       # (1152 x 160)
    jac_rank = int(np.linalg.matrix_rank(Jm, tol=1e-9))
    family_dim = 160 - jac_rank

    # ---- T0 anchors against the V208 artefact ----------------------------
    profH = profile(H)
    a_prof = art("V208_metric_bilinear_cone", "data.cone.kernel_profile")
    pdef = max(abs(profH[r] - float(a_prof[str(r)])) / float(a_prof[str(r)])
               for r in range(5))
    fB = (i3(H, Bop) / i3(Bop, Bop)) ** 2 * i3(Bop, Bop) / i3(H, H)
    fK = (i3(H, Kop) / i3(Kop, Kop)) ** 2 * i3(Kop, Kop) / i3(H, H)
    a_fB = float(art("V208_metric_bilinear_cone",
                     "data.implemented."
                     "hop_fraction_of_cone_generator"))
    a_fK = float(art("V208_metric_bilinear_cone",
                     "data.implemented."
                     "dressing_fraction_of_cone_generator"))
    a_b1 = int(art("V208_metric_bilinear_cone",
                   "data.cone.cone_branches_at_plus_1"))
    a_bm1 = int(art("V208_metric_bilinear_cone",
                    "data.cone.cone_branches_at_minus_1"))
    wH = np.linalg.eigvalsh((H[0] + H[0].conj().T) / 2)
    b1 = int(np.sum(np.abs(wH - 1.0) < 1e-6))
    bm1 = int(np.sum(np.abs(wH + 1.0) < 1e-6))
    t0ok = (pdef < 1e-9 and abs(fB - a_fB) < 1e-9 * abs(a_fB)
            and abs(fK - a_fK) < 1e-9 * abs(a_fK)
            and b1 == a_b1 and bm1 == a_bm1
            and jac_rank == 107 and family_dim == 53)
    T_committed = sum(profH[r] ** 2 for r in range(2, 8))
    data["t0"] = {"profile_defect": pdef, "fB": fB, "fK": fK,
                  "fB_art": a_fB, "fK_art": a_fK, "rank": jac_rank,
                  "family_dim": family_dim, "T_committed": T_committed,
                  "branches_at_1": b1, "branches_at_minus_1": bm1,
                  "branches_art": [a_b1, a_bm1],
                  "profile": {str(r): profH[r] for r in range(8)}}
    checks.append({"name": "T0 V208 anchors live: profile {:.1e}, "
                          "hop {:.10f} vs {:.10f}, dressing {:.10f} vs "
                          "{:.10f}, rank {} family {}; T_committed "
                          "{:.6e}".format(pdef, fB, a_fB, fK, a_fK,
                                          jac_rank, family_dim,
                                          T_committed),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V285 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1
    _, svJ, Vh = np.linalg.svd(Jm, full_matrices=False)
    NULLc = Vh[jac_rank:]                      # (53, 160) coefficients
    log("Jacobian rank {} family dim {}".format(jac_rank, family_dim))

    # ---- retract machinery (module state shared with workers) -----------
    cH = np.zeros(160)
    vH = vec3(H)
    for j in range(160):
        cH[j] = float(np.dot(vH, Bmat[j]))
    trH = unvec(cH @ Bmat)
    log("committed triple recovered from cov coords: defect {:.1e}"
        .format(np.abs(np.array(trH) - np.array(H)).max()))

    # fixed Jacobian of the residual at the committed point
    BOX = math.sqrt(i3(H, H) / 3.0)
    _ST.update({"Bmat": Bmat, "cH": cH, "NULL": NULLc, "BOX": BOX,
                "BR": (b1, bm1)})
    c_ref = cH.copy()
    J0 = np.zeros((len(_resid(c_ref)), 160))
    eps = 1e-6
    for j in range(160):
        cp = c_ref.copy()
        cp[j] += eps
        J0[:, j] = (_resid(cp) - _resid(c_ref)) / eps
    Q0, R0 = np.linalg.qr(J0)
    _ST.update({"J0": J0, "QR": (Q0, R0)})
    log("retract Jacobian built ({}x{}); box {:.4f}; start search "
        "{} x maxfev {}".format(J0.shape[0], J0.shape[1], BOX, NSTARTS,
                                MAXFEV))

    rows_s = []
    for b0 in range(0, NSTARTS, 5):
        rows_s += parpool.pmap(_start_worker,
                               list(range(b0, min(b0 + 5, NSTARTS))),
                               chunk=1)
        log("search {:2d}/{} starts, best T {:.6e}".format(
            len(rows_s), NSTARTS, min(r["T"] for r in rows_s)))
    rows_s.sort(key=lambda r: r["T"])
    best = rows_s[0]
    x = np.array(best["x"])
    c_star, rres = _retract(_ST["cH"] + x @ _ST["NULL"])
    tr_star = unvec(c_star @ _ST["Bmat"])
    fB_s = (i3(tr_star, Bop) / i3(Bop, Bop)) ** 2 * i3(Bop, Bop) \
        / i3(tr_star, tr_star)
    fK_s = (i3(tr_star, Kop) / i3(Kop, Kop)) ** 2 * i3(Kop, Kop) \
        / i3(tr_star, tr_star)
    dist = math.sqrt(i3(np.array(tr_star) - np.array(H),
                        np.array(tr_star) - np.array(H)) / i3(H, H))
    # 53-dim gradient at the minimizer (finite differences, retracted)
    g = np.zeros(53)
    hstep = 1e-5
    for a in range(53):
        xp = x.copy()
        xm = x.copy()
        xp[a] += hstep
        xm[a] -= hstep
        g[a] = (_obj(xp) - _obj(xm)) / (2 * hstep)
    n_box = int(np.sum(np.abs(x) > BOX - 1e-6))
    top5 = [r["T"] for r in rows_s[:5]]
    data["search"] = {"best_T": best["T"], "T_committed": T_committed,
                      "hop_star": fB_s, "dressing_star": fK_s,
                      "distance_rel": dist, "grad_norm": float(np.linalg.norm(g)),
                      "box_face_coords": n_box, "top5": top5,
                      "retract_defect": float(np.linalg.norm(rres)),
                      "x_best": [float(v) for v in x],
                      "all_T": [r["T"] for r in rows_s]}
    interior = (float(np.linalg.norm(g)) < 1e-6 * max(1.0, T_committed)
                and n_box == 0)
    found = (best["T"] < T_committed and fB_s > 0.0092 and interior)
    status = "CONE-LOCALLY-SELECTED" if found else "CONE-SELECTION-ABSENT"
    data["status"] = status
    print("V285 {}: T* {:.6e} vs committed {:.6e}, hop* {:.6f} (gate "
          "0.0092), dist {:.4f}, grad {:.2e}, box coords {}".format(
              status, best["T"], T_committed, fB_s, dist,
              np.linalg.norm(g), n_box), flush=True)
    checks.append({"name": "T1 minimization: {} starts, best T* {:.6e} "
                          "(committed {:.6e}), hop projection at the "
                          "minimizer {:.6f}, distance {:.4f}, interior "
                          "diagnostics grad {:.2e} box {}".format(
                              NSTARTS, best["T"], T_committed, fB_s,
                              dist, np.linalg.norm(g), n_box),
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V285 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V285_most_local_lightlike_generator", "tier": 2,
           "property": "locality (hopping-kernel tail) selection inside "
                       "the 53-dim exact-lightcone family through the "
                       "committed generator (tier 2, F8/O9)",
           "residual": True,
           "programme": "V285",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v285_most_local_lightlike_generator.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V285 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
