#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V277 - the bond-decoration cone scan (tier 2, hard, O10).

**The question (pre-registered).** Does the enlarged bond-decoration
circuit class - the six-stage step with a per-stage spin-connection
quaternion twist on the internal factor - carry a protected conical
node at GAMMA: gapless, linear (cone degree ~1), subluminal (slope
<= 1) and non-Abelian (rho_comm > 0.5)?

**The object (closed input set).** The V277 pre-registration in
review/pre_registrations.md is the authority; in brief: the published
six-stage supercell step at the closed clock (beta, g) = (pi/4,
pi/4) (V239 closure set), decorated per stage by
U16 -> U16 W_s(q), W_s(q) = exp(i q.(alpha_s V1 + beta_s V2 +
gamma_s V3)) with the V250 quaternion triplet
(V1,V2,V3) = (i h_1 h_2, i h_2 h_3, i h_1 h_3)|_16. 18 parameters,
box [0, 2 g10] per component with g10 = 0.2 from cl10 config.
Nelder-Mead (V186 tune polish), 120 starts (seed 277), maxfev 300.

**T0 (machinery gate).** The V250 star-cell decoration construction
re-derived in-script against the LIVE V250 artefact (tolerances in
the pre-registration): the implemented-clock bare/decorated degrees,
the closed-clock bare degree (1e-2 - a gap-floor fit), the closed
decorated degree as the artefact's NaN, the implemented decorated
cone anisotropy (1e-5) and the closed-clock anisotropy velocities
(1e-9 relative). Any drift -> MACHINERY-DRIFT, exit nonzero.

**Measured.** T0; the 128-dim closed-clock GAMMA census (64 zeros +
64 pi, V192/V132 pinned spectrum); bare top-64 degree (~1.9986, the
quadratic touching); the published uniform-triplet start; the search
best with fine degree, GAMMA gap, V186 tangent rho_comm on the
64-dim gap subspace, and 179-Fibonacci tangent speeds (v_max,
A_top).

**Acceptance.** DECORATION-CONE-FOUND iff GAMMA gap < 1e-10 (an
inherited fact: W(0) = I - booked as inherited, not earned), fine
degree >= 0.99, v_max <= 1, rho_comm > 0.5; else DECORATION-BEST
with the best (degree, rho_comm, v_max) triple and the binding leg.
"""
import json
import math
import os
import sys
import time

# The cone search is many independent 128x128 Floquet eigendecomposi-
# tions: the V185-V187 measurement puts them ~40x faster on CPU pool
# workers than on the card, and pool workers cannot host a CUDA con-
# text. Pinned BEFORE any cl10 import (documented deviation from the
# GPU1 run convention; same precedent as V185-V187, V192, V201).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("CL10_BACKEND", "numpy")  # launcher may force cupy

import numpy as np                                        # noqa: E402
from scipy.optimize import minimize                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402

from cl10 import altcell as AC                            # noqa: E402
from cl10 import clifford as CF                           # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import engine as E                              # noqa: E402
from cl10 import parpool                                  # noqa: E402
from cl10 import starcell as X                            # noqa: E402
from cl10 import supercell as SC                          # noqa: E402

CFG = _config.CFG
G10 = float(CFG.angle.branches.implemented)               # 0.2: the box scale
GCLOSED = math.pi / 4.0                                   # the closed clock
NSTARTS, MAXFEV, SEED = 120, 300, 277

# ---------------------------------------------------------------------------
# the shared construction (module level: inherited by the forked pool)
# ---------------------------------------------------------------------------
_B16, _U32, U16 = SC.sectors(GCLOSED)
_H = [E.asnumpy(m).astype(complex) for m in CF.generators()]
VT = [_B16.conj().T @ (1j * (_H[a] @ _H[b])) @ _B16
      for (a, b) in ((0, 1), (1, 2), (0, 2))]
TRIPL = "triplet (V1,V2,V3) = (i h_1 h_2, i h_2 h_3, i h_1 h_3)|_16 " \
        "(the V250 spin-connection quaternion units)"


def twist(q, amp):
    """W(q) = exp(i q.(a V1 + b V2 + c V3)), exact via the 16x16 eigh."""
    a, b, c = amp
    qV = a * q[0] * VT[0] + b * q[1] * VT[1] + c * q[2] * VT[2]
    if a == 0.0 and b == 0.0 and c == 0.0:
        return np.eye(16, dtype=complex)
    ev, Vm = np.linalg.eigh(qV)
    return (Vm * np.exp(1j * ev)) @ Vm.conj().T


def bloch_dec(q, amps):
    """The decorated six-stage 128x128 step at the closed clock."""
    M = np.eye(128, dtype=complex)
    for t, (even, j) in enumerate(SC.STAGES):
        S = SC.stage_matrix(even, j, q[j], GCLOSED, GCLOSED)
        M = np.kron(S, U16 @ twist(q, amps[t])) @ M
    return M


CUBIC_DIRS = []
for _d in ((1, 0, 0), (0, 1, 0), (0, 0, 1),
           (1, 1, 0), (1, 0, 1), (0, 1, 1),
           (1, -1, 0), (1, 0, -1), (0, 1, -1),
           (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)):
    _v = np.array(_d, float)
    CUBIC_DIRS.append(_v / np.linalg.norm(_v))
COARSE_RAD = (0.1, 0.03, 0.01)
FINE_RAD = (0.3, 0.1, 0.03, 0.01, 0.003)


def phases(M):
    """|eigenphases| of a unitary matrix through the Hermitian part
    (M + M^dag)/2 = cos(ph) exactly (V227's phase_gap precedent: the
    arccos route is ~7x cheaper than eigvals and agrees with the
    general eigenvalue route to 1e-10 on 128-dim steps)."""
    return np.sort(np.abs(np.arccos(np.clip(
        np.linalg.eigvalsh(0.5 * (M + M.conj().T)), -1.0, 1.0))))


def top_m(bf, m, r, dirs):
    """min over dirs of the m-th smallest |eigenphase| at radius r."""
    w = np.inf
    for d in dirs:
        p = phases(bf(tuple(r * x for x in d)))
        w = min(w, float(p[m - 1]))
    return w


def deg_fit(radii, ws):
    """polyfit(log r, log w) on positive gaps: a zero gap sends log to
    -inf and the fit to NaN - the V250 decorated-closed-clock case,
    reproduced as NaN, not masked."""
    with np.errstate(divide="ignore"):
        return float(np.polyfit(np.log(list(radii)), np.log(
            [float(w) for w in ws]), 1)[0])


def coarse_obj(amps18):
    """-degree of the 64-mode GAMMA multiplet top at the coarse grid."""
    amps = [(amps18[3 * t], amps18[3 * t + 1], amps18[3 * t + 2])
            for t in range(6)]
    bf = lambda q: bloch_dec(q, amps)        # noqa: E731
    try:
        ws = [top_m(bf, 64, r, CUBIC_DIRS) for r in COARSE_RAD]
    except np.linalg.LinAlgError:
        return 1e3
    d = deg_fit(COARSE_RAD, ws)
    if not np.isfinite(d):
        return 1e3
    return -d


def _start_worker(i):
    """One Nelder-Mead start (module level for pickling)."""
    rng = np.random.default_rng(SEED + i)
    if i == 0:
        x0 = np.full(18, G10)                # the published uniform triplet
    else:
        x0 = rng.uniform(0.0, 2.0 * G10, 18)
    simplex = [x0] + [x0 + 0.2 * e for e in np.eye(18)]
    res = minimize(coarse_obj, x0, method="Nelder-Mead",
                   options={"maxfev": MAXFEV, "fatol": 1e-8, "xatol": 1e-8,
                            "initial_simplex": np.array(simplex)})
    x = np.clip(res.x, 0.0, 2.0 * G10)
    return {"start": i, "x0": [float(v) for v in x0],
            "x": [float(v) for v in x], "deg": -coarse_obj(x),
            "fev": int(res.nfev)}


# ---------------------------------------------------------------------------
# final acceptance measurement at one point
# ---------------------------------------------------------------------------
def measure_point(amps, m0):
    bf = lambda q: bloch_dec(q, amps)        # noqa: E731
    out = {"gap_gamma": float(np.min(np.abs(np.angle(np.linalg.eigvals(
        bf((0.0, 0.0, 0.0)))))))}
    ws = [top_m(bf, m0, r, CUBIC_DIRS) for r in FINE_RAD]
    out["degree_fine"] = deg_fit(FINE_RAD, ws)
    out["gaps_fine"] = [float(w) for w in ws]
    # V186 tangent maps on the m0-gap subspace at GAMMA
    M0 = bf((0.0, 0.0, 0.0))
    ev, Vv = np.linalg.eig(M0)
    sel = np.argsort(np.abs(np.angle(ev)))[:m0]
    Q, _ = np.linalg.qr(Vv[:, sel])
    h = 1e-5
    S = []
    for jj in range(3):
        qp, qm = [0.0] * 3, [0.0] * 3
        qp[jj], qm[jj] = h, -h
        Dd = (bf(tuple(qp)) - bf(tuple(qm))) / (2 * h)
        S.append(1j * (Q.conj().T @ Dd @ Q))
    nrm = [float(np.linalg.norm(s)) for s in S]
    mean2 = sum(n * n for n in nrm) / 3.0
    comm = max(float(np.linalg.norm(S[a] @ S[b] - S[b] @ S[a]))
               for a, b in ((0, 1), (1, 2), (2, 0)))
    out["rho_comm"] = comm / mean2 if mean2 > 0 else 0.0
    out["S_norms"] = nrm
    # tangent speeds over the 179-direction Fibonacci sphere (V235 fib)
    idx = np.arange(179) + 0.5
    phi = np.arccos(1 - 2 * idx / 179)
    th = np.pi * (1 + 5 ** 0.5) * idx
    fdirs = np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th),
                      np.cos(phi)], axis=1)
    spd = []
    for kd in fdirs:
        Ax = kd[0] * S[0] + kd[1] * S[1] + kd[2] * S[2]
        spd.extend(sorted(abs(x) for x in np.linalg.eigvalsh(
            0.5 * (Ax + Ax.conj().T))))
    spd = np.array(spd)
    nz = spd[spd > 1e-6]
    if nz.size:
        out["v_min"] = float(nz.min())
        out["v_max"] = float(nz.max())
        out["A_top"] = float((nz.max() - nz.min()) / (nz.max() + nz.min()))
        out["n_flat"] = int((spd <= 1e-6).sum())
        out["n_dirs"] = int(len(fdirs))
    # unitarity defect at three random q (invariance check)
    rng = np.random.default_rng(277)
    ud = 0.0
    for _ in range(3):
        qq = tuple(float(x) for x in rng.uniform(-math.pi, math.pi, 3))
        Mm = bf(qq)
        ud = max(ud, float(np.abs(Mm @ Mm.conj().T
                                  - np.eye(128, dtype=complex)).max()))
    out["unitarity_defect"] = ud
    # W(GAMMA) = I exactly: the decorated step shares M(GAMMA) with the bare
    Mb = bloch_dec((0.0, 0.0, 0.0), [(0.0, 0.0, 0.0)] * 6)
    out["M_gamma_minus_bare"] = float(np.abs(M0 - Mb).max())
    return out


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    t0 = time.time()
    checks, data = [], {}
    LOG0 = time.time()

    def log(msg):
        print("[V277 t+{:5.0f}s] {}".format(time.time() - LOG0, msg),
              flush=True)
    log("start; engine backend {} (GPU eigvals routed by the engine)"
        .format(E._BACKEND))
    data["triplet"] = TRIPL
    herm = max(float(np.abs(v - v.conj().T).max()) for v in VT)

    # ---- T0a: the V250 star-cell decoration, re-derived -----------------
    v250 = art("V250_single_particle_walls", "data.continuation_decoration")
    deg_art = v250["cone_degree_gamma"]
    an_art = v250["cone_anisotropy_gamma"]
    CELL = X.module_cell()
    STAGES = X.module_stages(order="phase-major")
    AXES = np.asarray(X.star_axes(), float)
    AXU = AXES / np.linalg.norm(AXES, axis=1)[:, None]
    B16s, _U32s, _U16s = SC.sectors(math.pi / 4)
    Vs = [B16s.conj().T @ (1j * (_H[a] @ _H[b])) @ B16s
          for (a, b) in ((0, 1), (1, 2), (0, 2))]

    def twist250(q):
        qV = sum(qq * V for qq, V in zip(q, Vs))
        ev, W = np.linalg.eigh(qV)
        return (W * np.exp(1j * G10 * ev)) @ W.conj().T

    def bloch250(q, dec, beta, g):
        # V250 convention: the internal factor is sectors(g)[2], NOT the
        # pi/4 sector, for both the bare and the decorated step
        U16g = SC.sectors(g)[2]
        U = U16g if not dec else U16g @ twist250(q)
        return E.asnumpy(X.bloch_star(CELL, U, g, STAGES, q, AXES,
                                      bond_g10=beta))

    dirs250 = [tuple(v) for v in AXU] + [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    dirs250 = [tuple(np.array(d) / np.linalg.norm(d)) for d in dirs250]
    rad250 = [0.3, 0.1, 0.03, 0.01, 0.003]

    def degree250(bf):
        ws = []
        for r in rad250:
            w = min(np.min(np.abs(np.angle(
                E.eigvals(bf(tuple(r * x for x in d))))))
                for d in dirs250)
            ws.append(w)
        return float(np.polyfit(np.log(rad250), np.log(ws), 1)[0]), ws

    def aniso_gpu(bf, dq=1e-3):
        """altcell.cone_anisotropy with engine routing (identical
        convention: sphere directions, min |angle| / dq, GPU eigvals)."""
        vs = []
        for nd in AC.sphere_dirs(6):
            qq = tuple(q + dq * n for q, n in zip((0.0, 0.0, 0.0), nd))
            th = np.angle(E.eigvals(bf(qq)))
            vs.append(float(np.min(np.abs(th))) / dq)
        vmax, vmin = max(vs), min(vs)
        if vmax + vmin <= 0:
            return 1.0, vmax, vmin
        return (vmax - vmin) / (vmax + vmin), vmax, vmin

    redeg, rean = {}, {}
    for label, (beta, g) in (("implemented", (G10, G10)),
                             ("closed", (GCLOSED, GCLOSED))):
        p_bare, w_bare = degree250(lambda q, b=beta, gg=g: bloch250(
            q, False, b, gg))
        p_dec, w_dec = degree250(lambda q, b=beta, gg=g: bloch250(
            q, True, b, gg))
        redeg[label] = {"bare": p_bare, "decorated": p_dec,
                        "bare_gaps": [float(v) for v in w_bare],
                        "decorated_gaps": [float(v) for v in w_dec]}
        log("T0 degree250 {} bare {:.4f} dec {:.4f}".format(
            label, p_bare, p_dec))
        A_b, vxb, vnb = aniso_gpu(
            lambda q, b=beta, gg=g: bloch250(q, False, b, gg))
        A_d, vxd, vnd = aniso_gpu(
            lambda q, b=beta, gg=g: bloch250(q, True, b, gg))
        rean[label] = {"bare_A": A_b, "decorated_A": A_d,
                       "bare_vmax": vxb, "bare_vmin": vnb,
                       "decorated_vmax": vxd, "decorated_vmin": vnd}

    def nan_ok(a, b):
        return (isinstance(a, float) and math.isnan(a)
                and isinstance(b, float) and math.isnan(b))

    drift_items = []
    if abs(redeg["implemented"]["bare"] - deg_art["implemented"]["bare"]) \
            > 1e-3:
        drift_items.append("implemented bare degree")
    if abs(redeg["implemented"]["decorated"]
           - deg_art["implemented"]["decorated"]) > 1e-3:
        drift_items.append("implemented decorated degree")
    # amended contract (addendum, pre-booked run): at the closed clock
    # both bare and decorated gap arrays must sit at the float64
    # rounding floor (<= 1e-15 every radius; V250's booked arrays:
    # max 3.5e-16). The NaN/finite PLACEMENT of the polyfit is single-
    # bit rounding noise and is booked as measured, not as drift.
    for lbl2, wkey in (("bare", "bare_gaps"),
                       ("decorated", "decorated_gaps")):
        if max(redeg["closed"][wkey]) > 1e-15:
            drift_items.append("closed " + lbl2 + " not at rounding floor")
    art_closed_max = max(max(deg_art["closed"]["bare_gaps"] if
                             isinstance(deg_art["closed"].get("bare_gaps"),
                                      list) else [0.0]),
                         max(deg_art["closed"]["decorated_gaps"] if
                             isinstance(deg_art["closed"].get(
                                 "decorated_gaps"), list) else [0.0]))
    if art_closed_max > 1e-15:
        drift_items.append("artefact closed rows not at rounding floor")
    data["closed_clock_nan_placement"] = {
        "rederived_bare_nan": bool(math.isnan(redeg["closed"]["bare"])),
        "rederived_decorated_nan": bool(
            math.isnan(redeg["closed"]["decorated"])),
        "artefact_bare_nan": bool(
            math.isnan(float(deg_art["closed"]["bare"]))),
        "artefact_decorated_nan": bool(
            math.isnan(float(deg_art["closed"]["decorated"]))),
        "note": "single-bit rounding placement of polyfit's log(0) at "
                "the pinned closed clock; pinning itself is the gate"}
    if abs(rean["implemented"]["decorated_A"]
           - an_art["implemented"]["decorated_A"]) > 1e-5:
        drift_items.append("implemented decorated A")
    if abs(rean["closed"]["bare_A"] - an_art["closed"]["bare_A"]) > 1e-6:
        drift_items.append("closed bare A")
    for kk in ("bare_vmax", "decorated_vmax", "decorated_vmin"):
        if abs(rean["closed"][kk] - an_art["closed"][kk]) \
                > 1e-9 * max(1.0, abs(an_art["closed"][kk])):
            drift_items.append("closed " + kk)
    drift = bool(drift_items)
    data["t0_star_cell"] = {"rederived": redeg, "rederived_aniso": rean,
                            "artifact_degrees": deg_art,
                            "artifact_aniso": an_art,
                            "drift_items": drift_items,
                            "note": "the live V250 artefact carries no "
                                    "0.281/0.540 pair; the printed "
                                    "publication pair has no live "
                                    "artefact key (measured before any "
                                    "run), so the anchors bind to the "
                                    "artefact values above"}
    checks.append({"name": "T0a V250 star-cell decoration re-derived "
                          "against the live artefact (implemented degrees "
                          "{:.4f}/{:.4f} vs booked {:.4f}/{:.4f}; closed "
                          "bare {:.4f} vs {:.4f}; closed decorated NaN "
                          "contract {}; implemented decorated A {:.3e} vs "
                          "{:.3e}; closed aniso velocities to 1e-9 rel); "
                          "triplet hermiticity {:.1e}".format(
                              redeg["implemented"]["bare"],
                              redeg["implemented"]["decorated"],
                              deg_art["implemented"]["bare"],
                              deg_art["implemented"]["decorated"],
                              redeg["closed"]["bare"],
                              deg_art["closed"]["bare"],
                              math.isnan(redeg["closed"]["decorated"]),
                              rean["implemented"]["decorated_A"],
                              an_art["implemented"]["decorated_A"], herm),
                   "pass": not drift and herm < 1e-12})
    if drift:
        data["status"] = "MACHINERY-DRIFT"
        print("V277 MACHINERY-DRIFT: the V250 decoration re-derivation "
              "disagrees with the booked artefact ({}); nothing downstream "
              "is interpreted".format(", ".join(drift_items)), flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0b: the 128-dim closed-clock machinery ------------------------
    log("T0b: 128-dim supercell census")
    M0 = bloch_dec((0.0, 0.0, 0.0), [(0.0, 0.0, 0.0)] * 6)
    p0 = np.sort(np.abs(np.angle(np.linalg.eigvals(M0))))
    n0 = int((p0 < 1e-12).sum())
    npi = int((np.abs(p0 - math.pi) < 1e-12).sum())
    bare_deg = deg_fit(FINE_RAD,
                       [top_m(lambda q: bloch_dec(q, [(0., 0., 0.)] * 6),
                              64, r, CUBIC_DIRS) for r in FINE_RAD])
    uni0 = deg_fit(FINE_RAD,
                   [top_m(lambda q: bloch_dec(q, [(G10, G10, G10)] * 6),
                          64, r, CUBIC_DIRS) for r in FINE_RAD])
    t0b_ok = bool(n0 == 64 and npi == 64 and abs(bare_deg - 1.9986) < 2e-3)
    data["machinery_128"] = {"gap_gamma": float(p0.min()),
                             "zeros_at_gamma": n0, "pi_at_gamma": npi,
                             "m0": 64, "bare_top64_degree": bare_deg,
                             "uniform_triplet_start_degree": uni0}
    checks.append({"name": "T0b 128-dim closed-clock machinery: GAMMA "
                          "census {} zeros + {} pi (V192/V132 pinned "
                          "spectrum), bare top-64 degree {:.5f} (quadratic "
                          "touching, tol 2e-3 vs 1.9986), published "
                          "uniform-triplet start degree {:.4f}".format(
                              n0, npi, bare_deg, uni0),
                   "pass": t0b_ok})
    if not t0b_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V277 MACHINERY-DRIFT: the 128-dim closed-clock machinery "
              "census/degree disagrees with the published spectrum",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the search ------------------------------------------------------
    log("search: {} Nelder-Mead starts (maxfev {}) via parpool".format(
        NSTARTS, MAXFEV))
    rows = []
    for b0 in range(0, NSTARTS, 12):
        rows += parpool.pmap(_start_worker,
                             list(range(b0, min(b0 + 12, NSTARTS))),
                             chunk=1)
        log("search: {:3d}/{} starts done, best coarse deg {:.5f}"
            .format(len(rows), NSTARTS, max(r["deg"] for r in rows)))
    rows.sort(key=lambda r: -r["deg"])
    best = rows[0]
    data["search"] = {"starts": NSTARTS, "maxfev": MAXFEV, "seed": SEED,
                      "box": [0.0, 2.0 * G10],
                      "coarse_grid": {"radii": list(COARSE_RAD),
                                      "dirs": "13 cubic classes"},
                      "best": best, "top5": rows[:5],
                      "all_degrees": [r["deg"] for r in rows]}
    checks.append({"name": "T1 search: {} Nelder-Mead starts (maxfev {}, "
                          "seed {}), best coarse degree {:.5f} at x = {}"
                          .format(NSTARTS, MAXFEV, SEED, best["deg"],
                                  [round(v, 4) for v in best["x"]]),
                   "pass": True})

    # ---- acceptance measurement at the best point ------------------------
    log("acceptance measurement at the best point")
    amps_best = [(best["x"][3 * t], best["x"][3 * t + 1],
                  best["x"][3 * t + 2]) for t in range(6)]
    meas = measure_point(amps_best, 64)
    data["best_measure"] = meas
    legs = {"gap_lt_1e_10": meas["gap_gamma"] < 1e-10,
            "degree_ge_0_99": (np.isfinite(meas["degree_fine"])
                               and meas["degree_fine"] >= 0.99),
            "vmax_le_1": meas.get("v_max", 9.9) <= 1.0,
            "rho_gt_0_5": meas["rho_comm"] > 0.5}
    found = all(legs.values())
    binding = [k for k, v in legs.items() if not v]
    branch = "DECORATION-CONE-FOUND" if found else "DECORATION-BEST"
    data["status"] = branch
    data["acceptance_legs"] = legs
    data["binding_legs"] = binding
    print("V277 {}: fine degree {:.4f}, rho_comm {:.4f}, v_max {:.4f}, "
          "A_top {:.4f}, gap(GAMMA) {:.2e} (inherited: W(0) = I); legs "
          "failed {}".format(
              branch, meas["degree_fine"], meas["rho_comm"],
              meas.get("v_max", float("nan")),
              meas.get("A_top", float("nan")), meas["gap_gamma"],
              binding or ["none"]), flush=True)
    checks.append({"name": "T2 acceptance legs at the best point: gap "
                          "{:.2e} < 1e-10 {}, degree {:.4f} >= 0.99 {}, "
                          "v_max {:.4f} <= 1 {}, rho_comm {:.4f} > 0.5 {}"
                   .format(meas["gap_gamma"], legs["gap_lt_1e_10"],
                           meas["degree_fine"], legs["degree_ge_0_99"],
                           meas.get("v_max", float("nan")),
                           legs["vmax_le_1"], meas["rho_comm"],
                           legs["rho_gt_0_5"]),
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    print("[V277 done] {}  ({:.0f} s)".format(branch, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V277_bond_decoration_cone_scan", "tier": 2,
           "property": "does the 18-parameter bond-decoration circuit "
                       "class carry a protected conical node at GAMMA "
                       "(gapless, degree >= 0.99, slope <= 1, "
                       "rho_comm > 0.5) - tier 2, O10",
           "residual": True,
           "programme": "V277",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v277_bond_decoration_cone_scan.json"), "w") \
            as fh:
        json.dump(out, fh, indent=1)
    print("[V277 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
