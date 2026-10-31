#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V329 - does any selection functional have an interior minimum in r?
(tier 3, O2, theory-hunt closure programme.)

**The question (pre-registered).** V303/V307 booked the junction
amplitude r = kappa/w0 as energy-driven to the boundary r* = 2 (tree-
level r-flatness by the trace theorem, one-loop monotone to the
singularity). Is that booked edge the whole stationary structure: over
the r-analogue of V317's selection-functional census, does ANY booked
functional have an interior minimum of the energy object on
r in [1.0, 1.999], or is the only stationary structure the boundary
r = 2 singularity (with the r = 1 endpoint degeneracy)?

**Closed class (booked pre-run, verbatim).** The object is V303/V307's
canonical junction data, their bodies verbatim: NODES = M_POINTS (12),
the pi/2-support edge list, the adjacency A_S of every non-empty subset
of the 12 nodes (4095 split sets), canonical levels |1 + (r/2) mu_k|
and one-loop dE = sum log lv^2 with V307's massless convention
(any lv < 1e-14 gives -inf). V317's four functionals take no r argument
(booked honestly; V317 is imported by path and its weight
G = 1/(4 ln 45) used verbatim), so the census runs the four functional
SHAPES on the r-dependent object:
  F1(r) = dE(A12, r)                       (T*/membership analogue)
  F2(r) = G * (E2(r) - E1(r)), E1 <= E2 the two lowest of the 4095
          subset energies (weighted-gap analogue; grid points where
          E1 or E2 is -inf excluded, count booked)
  F3(r) = sum of the three smallest of the twelve levels(A12, r)
          (fold/channel-energy analogue)
  F4(r) = min over nodes of 1/|(M(r)^-1)_ii|, M(r) = 2 I + r A12
          (node p_A analogue; M singular on the grid booked as gap 0)
Grid: r = 1.0, 1.0001, ..., 1.9990 (step 1e-4, 9991 points). Interior
minimum at index j iff the forward difference changes sign < 0 -> > 0
at j AND |F(j+1) - 2F(j) + F(j-1)| >= 1e-12. The grid is the booked
class: no widening.

**T0 (pre-registered).** V307 data.t3.r_singularity =
1.9999999999999791 re-derived by V307's own bisection (|d| <= 1e-3);
V307 data.t2.r_atm = 1.9191154379808832 re-derived as
loop_ladder(17)*EV_PER_GEV/M_ATM_EV - 1 (|d| <= 1e-9); V305
data.spectra.selected.r = 2.0 exactly. Failure: MACHINERY-DRIFT.

**Branch (pre-registered).** R-INTERIOR-MIN-FOUND (some functional has
an interior minimum: r*, the functional, distances to r_atm and to 2
booked) / R-INTERIOR-EMPTY-DECLARED (all four counts zero: the only
stationary structure is the booked boundary - r stays declared).
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

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
sys.path.insert(0, exp("V317_selection_conjunction"))

import V317_selection_conjunction as V317                 # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402

PI2 = math.pi / 2
SUPPORT = [(0, 1, -1), (1, 0, -1), (1, -1, 0)]
NODES = SC.M_POINTS
G_ENTROPY = V317.G_ENTROPY          # 1/(4 ln 45), the V317 weight verbatim
SUPPORT_SET = set()
for _s in SUPPORT:
    SUPPORT_SET.add(tuple(_s))
    SUPPORT_SET.add(tuple(-v for v in _s))


def art(name, key):
    """A sibling verdict JSON, read by experiment name (house art())."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# ---- V303/V307 machinery, bodies verbatim --------------------------------
def node_index(q):
    for i, p in enumerate(NODES):
        if all(abs(a - b) < 1e-12 for a, b in zip(q, p)):
            return i
    raise KeyError(q)


def edge_list():
    edges = []
    for i, j in itertools.combinations(range(12), 2):
        d = tuple(round((NODES[i][a] - NODES[j][a]) / PI2) for a in range(3))
        if d in SUPPORT_SET:
            edges.append((i, j))
    return edges


def adjacency(S, edges):
    idx = {v: k for k, v in enumerate(S)}
    A = np.zeros((len(S), len(S)))
    for i, j in edges:
        if i in idx and j in idx:
            A[idx[i], idx[j]] += 1.0
            A[idx[j], idx[i]] += 1.0
    return A


def levels(A, r):
    return np.abs(1.0 + 0.5 * r * np.linalg.eigvalsh(A))


def dE(A, r):
    lv = levels(A, r)
    if np.any(lv < 1e-14):
        return -math.inf
    return float(np.sum(np.log(lv ** 2)))


def main():
    t0 = time.time()
    checks, data = [], {}
    edges = edge_list()
    A12 = adjacency(list(range(12)), edges)
    mu12 = np.linalg.eigvalsh(A12)
    w12, V12 = np.linalg.eigh(A12)
    P2 = np.abs(V12) ** 2                     # |eigenvector components|^2

    # ---- T0 gates -------------------------------------------------------
    lo, hi = 1.5, 2.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.min(levels(A12, mid)) > 1e-14:
            lo = mid
        else:
            hi = mid
    r_sing = 0.5 * (lo + hi)
    a_sing = float(art("V307_onset_runaway_structure",
                       "data.t3.r_singularity"))
    t0a = abs(r_sing - a_sing)
    S17 = nv.loop_ladder(17) * nv.EV_PER_GEV
    r_atm = S17 / nv.M_ATM_EV - 1.0
    a_ratm = float(art("V307_onset_runaway_structure", "data.t2.r_atm"))
    t0b = abs(r_atm - a_ratm)
    a_rsel = float(art("V305_onset_neutrino_contacts",
                       "data.spectra.selected.r"))
    data["t0"] = {"r_singularity_rederived": r_sing,
                  "r_singularity_artifact": a_sing, "abs_diff": t0a,
                  "r_atm_rederived": r_atm, "r_atm_artifact": a_ratm,
                  "r_atm_abs_diff": t0b,
                  "selected_r_artifact": a_rsel}
    checks.append({"name": "T0 the booked onset anchors: singularity "
                          "re-derived {:.16g} vs artifact {:.16g} (|d| "
                          "{:.1e} <= 1e-3); r_atm re-derived {:.16g} vs "
                          "{:.16g} (|d| {:.1e} <= 1e-9); V305 selected "
                          "r = {:.4f} == 2.0".format(
                              r_sing, a_sing, t0a, r_atm, a_ratm, t0b,
                              a_rsel),
                   "pass": bool(t0a <= 1e-3 and t0b <= 1e-9
                                and a_rsel == 2.0)})
    if not (t0a <= 1e-3 and t0b <= 1e-9 and a_rsel == 2.0):
        data["status"] = "MACHINERY-DRIFT"
        print("V329 MACHINERY-DRIFT: the V307/V305 onset anchors do not "
              "re-derive; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the r grid; the 4095 subset spectra are r-independent (levels
    # are |1 + (r/2) mu|), so eigenvalues are booked once ----------------
    grid = np.arange(1.0, 1.9990 + 0.5e-4, 1e-4)
    n_grid = grid.size
    MU = np.zeros((4095, 12))                 # 0-padding: lv = 1, log = 0
    m = 0
    for mask in range(1, 1 << 12):
        S = [v for v in range(12) if (mask >> v) & 1]
        MU[m, :len(S)] = np.linalg.eigvalsh(adjacency(S, edges))
        m += 1

    F1 = np.empty(n_grid)
    F3 = np.empty(n_grid)
    F4 = np.empty(n_grid)
    E1 = np.empty(n_grid)
    E2 = np.empty(n_grid)
    for j, r in enumerate(grid):
        lv12 = np.abs(1.0 + 0.5 * r * mu12)
        F1[j] = (-math.inf if np.any(lv12 < 1e-14)
                 else float(np.sum(np.log(lv12 ** 2))))
        F3[j] = float(np.sort(lv12)[:3].sum())
        den = 2.0 + r * w12
        if np.any(np.abs(den) < 1e-12):
            F4[j] = 0.0
        else:
            F4[j] = float(np.min(1.0 / np.abs((P2 / den).sum(axis=1))))
        lv = np.abs(1.0 + 0.5 * r * MU)
        with np.errstate(divide="ignore"):
            Eall = 2.0 * np.log(lv).sum(axis=1)
        bad = (lv < 1e-14).any(axis=1)
        Eall[bad] = -math.inf
        part = np.partition(Eall, 1)
        E1[j], E2[j] = part[0], part[1]
    with np.errstate(invalid="ignore"):
        gap = E2 - E1
    F2 = G_ENTROPY * gap                # -inf rows propagate as excluded

    # ---- the interior-minimum test (booked rule) ------------------------
    def minima(F, name):
        count_nonfinite = int(np.count_nonzero(~np.isfinite(F)))
        events = []
        mins = []
        for j in range(1, n_grid - 1):
            if not (np.isfinite(F[j - 1]) and np.isfinite(F[j])
                    and np.isfinite(F[j + 1])):
                continue
            d_prev = F[j] - F[j - 1]
            d_next = F[j + 1] - F[j]
            if d_prev < 0.0 < d_next:
                curv = abs(F[j + 1] - 2.0 * F[j] + F[j - 1])
                if curv >= 1e-12:
                    mins.append({"r_star": float(grid[j]),
                                 "curvature": float(curv),
                                 "F_at_star": float(F[j]),
                                 "dist_to_r_atm": float(abs(grid[j]
                                                            - r_atm)),
                                 "dist_to_2": float(abs(grid[j] - 2.0))})
                else:
                    events.append({"r": float(grid[j]),
                                   "type": "sign change below curvature "
                                           "floor",
                                   "curvature": float(curv)})
            elif d_prev > 0.0 > d_next:
                events.append({"r": float(grid[j]), "type": "interior "
                               "maximum", "curvature":
                               float(abs(F[j + 1] - 2.0 * F[j]
                                         + F[j - 1]))})
        decim = [{"r": float(grid[j]), "F": (float(F[j])
                                             if np.isfinite(F[j]) else
                                             "-inf")}
                 for j in range(0, n_grid, 200)]
        return name, {"grid": "1.0..1.999 step 1e-4 (n = %d)" % n_grid,
                      "n_excluded_nonfinite": count_nonfinite,
                      "n_interior_minima": len(mins), "minima": mins,
                      "other_stationarity_events": events,
                      "decimated_grid_every_200": decim}

    rows = [minima(F1, "F1_full_set_one_loop_energy"),
            minima(F2, "F2_weighted_best_second_gap"),
            minima(F3, "F3_three_channel_energy_fold"),
            minima(F4, "F4_node_local_junction_gap")]
    data["functionals"] = {name: row for name, row in rows}
    data["sites_selected_per_functional"] = {
        name: int(row["n_interior_minima"]) for name, row in rows}
    total_mins = sum(row["n_interior_minima"] for _, row in rows)

    branch = ("R-INTERIOR-MIN-FOUND" if total_mins > 0
              else "R-INTERIOR-EMPTY-DECLARED")
    data["status"] = branch
    data["wall_time"] = time.time() - t0
    checks.append({"name": "branch {}: interior minima per functional "
                           "{} (grid {}); the boundary singularity at "
                           "r = 2 and the r = 1 endpoint degeneracy are "
                           "the booked edge; total minima {}".format(
                               branch,
                               data["sites_selected_per_functional"],
                               rows[0][1]["grid"], total_mins),
                   "pass": True})
    print("V329 {}: per-functional interior-minimum counts {}; wall "
          "{:.0f} s.".format(branch,
                             data["sites_selected_per_functional"],
                             data["wall_time"]), flush=True)
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V329_junction_interior_minimality", "tier": 3,
           "property": "does any booked selection functional have an "
                       "interior minimum of the junction energy object on "
                       "r in [1.0, 1.999] - or is the boundary r = 2 the "
                       "only stationary structure (tier 3, O2)",
           "residual": True,
           "programme": "theory-hunt V329",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v329_junction_interior_minimality.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V329 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
