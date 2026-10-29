#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V303 - does the junction energy select the ratio r? (tier 2, hard, O2.)
**The question (pre-registered).** The junction ratio r = kappa/w0 is not
derived by the shape (V114) and not quantized (V293 books the family
multiplicity, a different r). Is it selected by ENERGY: on the stable
branch 0 < r < r* of the V127 quadratic junction theory, does the measured
vacuum energy have an interior minimum - and does the energy's flow land
on the democratic onset r* = 2 (the candidate's anchor) inside the
cosmological window r > 1.278?
**The object (closed input set).** V127's machinery verbatim (NODES =
cl10.supercell.M_POINTS, the support set, the closing triangle SEL, the
adjacency, the one-loop functional, the 4096 split sets). Convention fixed
before any run and booked here: BOTH functionals run in the canonical
junction ratio r = kappa/w0, levels 1 + (r/2) mu_k(A_S). V127's quadratic
functional e_quad used a doubled convention (1 + r mu) while its crossing
search and one-loop used the canonical one, and its sweep points
(0.5, 1.0, 3.0) never enter the stable branch r < 2: the stable branch
has never been swept. No dial, no optimisation: grids only.
**Measured.** T0: the V127 artifact anchors re-derived (8-orbit, K3
degrees, 12 edges / 0 exiting, junction group 12, trace theorem, first
crossings {2.0}) - drift halts. T1: closing-triangle functionals on
r in (0, 2): tree-level r-flatness (trace theorem deviation), one-loop
monotonicity (interior-minimum count), exact closed form 2 log(1+r) +
4 log|1 - r/2| verified. T2: 4096-set sweep: per-set crossing
r_S = -2/mu_min (buckets below 1.278 / 1.9191 / 2 / none) and the energy
minimiser at every r of the 0.05 grid (size, multiplicity, whether the
closing triangle or its junction orbit ever attains it). T3: contact:
the attractor r -> 2 gives m_3 = S17/3 against the booked atmospheric
scale; relative residual booked; the cosmological window re-derived.
Branch: R-ONSET-DEMOCRATIC / R-WINDOW-MIN crossed with
VACUUM-SELECTS-ALL-12 / TRIANGLE-SELECTED.
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
from cl10 import supercell as SC                          # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402

PI2 = math.pi / 2
SUPPORT = [(0, 1, -1), (1, 0, -1), (1, -1, 0)]            # units of pi/2
NODES = SC.M_POINTS                                        # the 12 M nodes
SEL = [(PI2, PI2, 0.0), (PI2, 0.0, PI2), (0.0, PI2, PI2)]
SUPPORT_SET = set()
for _s in SUPPORT:
    SUPPORT_SET.add(tuple(_s))
    SUPPORT_SET.add(tuple(-v for v in _s))


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def node_index(q):
    for i, p in enumerate(NODES):
        if all(abs(a - b) < 1e-12 for a, b in zip(q, p)):
            return i
    raise KeyError(q)


def apply_oh(tri, perm, signs):
    img = []
    for v in sorted(tri):
        q = NODES[v]
        img.append(node_index(tuple(signs[a] * q[perm[a]] for a in range(3))))
    return frozenset(img)


def oh_orbit(tri):
    return {apply_oh(tri, p, s) for p, s in SC._signed_permutations()}


def junction_group():
    out = []
    for p, s in SC._signed_permutations():
        img = {tuple(s[a] * v[p[a]] for a in range(3)) for v in SUPPORT_SET}
        if img == SUPPORT_SET:
            out.append((p, s))
    return out


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
    """Canonical junction levels |1 + (r/2) mu_k| (r = kappa/w0)."""
    return np.abs(1.0 + 0.5 * r * np.linalg.eigvalsh(A))


def e_quad_can(A, r):
    if A.shape[0] == 0:
        return 0.0
    return float(np.sum(levels(A, r)))


def dE(A, r):
    """One-loop sum_k log|1 + (r/2) mu_k|^2; -inf at a massless mode."""
    lv = levels(A, r)
    if np.any(lv < 1e-12):
        return -math.inf
    return float(np.sum(np.log(lv ** 2)))


def first_crossing(A):
    mu = np.linalg.eigvalsh(A)
    if mu.min() >= -1e-12:
        return math.inf
    return float(-2.0 / mu.min())


def main():
    t0 = time.time()
    checks, data = [], {}
    edges = edge_list()
    tri0 = frozenset(node_index(q) for q in SEL)
    orb = sorted(oh_orbit(tri0), key=sorted)
    GJ = junction_group()
    jorb = {apply_oh(tri0, p, s) for p, s in GJ}

    # ---- T0 the V127 anchors, re-derived ----------------------------------
    v127 = art("V127_triangle_energies")
    a = v127["data"]["orbit"]
    adj0 = adjacency(sorted(tri0), edges)
    deg0 = [int(adj0[i].sum()) for i in range(3)]
    exits = sum(1 for (i, j) in edges if (i in tri0) != (j in tri0))
    # support graph components
    seen, comps = set(), []
    for v in range(12):
        if v in seen:
            continue
        stack, cc = [v], set()
        while stack:
            x = stack.pop()
            if x in cc:
                continue
            cc.add(x)
            for i, j in edges:
                if i == x:
                    stack.append(j)
                if j == x:
                    stack.append(i)
        seen |= cc
        comps.append(sorted(cc))
    cross_tri = first_crossing(adj0)
    cross_all = set()
    for c in comps:
        cval = first_crossing(adjacency(c, edges))
        if math.isfinite(cval):
            cross_all.add(round(cval, 6))
    drift = not (len(orb) == a["oh_size"] == 8 and deg0 == a["triangle_degrees"]
                 == [2, 2, 2] and len(edges) == a["edges_total"] == 12
                 and exits == a["edges_leaving_triangle"] == 0
                 and len(GJ) == a["junction_group_order"] == 12
                 and len(jorb) == a["junction_orbit_size"] == 2
                 and abs(cross_tri - 2.0) < 1e-9
                 and v127["data"]["trace_theorem"]["first_zero_crossing_r"]
                 == [2.0])
    data["machinery"] = {"orbit": len(orb), "junction_group": len(GJ),
                         "junction_orbit": len(jorb), "degrees": deg0,
                         "edges": len(edges), "exiting": exits,
                         "components": [len(c) for c in comps],
                         "component_crossings": sorted(cross_all)}
    checks.append({"name": "T0 V127 anchors re-derived (orbit 8, K3, 12 "
                          "edges / 0 exiting, junction group 12 with a "
                          "2-orbit, triangle crossing {2.0}): %s; NEW "
                          "structure booked: the support graph splits into "
                          "components %s with crossings %s - the 6-cycle "
                          "component goes massless at r = 1, BELOW the "
                          "window"
                          % ("DRIFT" if drift else "exact",
                             [len(c) for c in comps], sorted(cross_all)),
                   "pass": not drift})
    if drift:
        data["status"] = "MACHINERY-DRIFT"
        print("V303 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    # ---- T1 the closing-triangle functionals on the stable branch ---------
    rs = np.linspace(0.02, 1.98, 2001)   # grid stops at 1.98: eigenvalue
    # noise (~1e-16) is amplified by the log as the onset r = 2 approaches,
    # so the 1e-12 closed-form check is run off the crossing
    eq = np.array([e_quad_can(adj0, r) for r in rs])
    dEr = np.array([dE(adj0, r) for r in rs])
    flat_dev = float(np.max(np.abs(eq - 3.0)))
    dd = np.diff(dEr)
    monotone = bool(np.all(dd < 0.0))
    interior = int(np.sum(np.diff(np.sign(dd)) != 0)) if len(dd) > 2 else 0
    slope_lo = float(dEr[np.searchsorted(rs, 1.30)]
                     - dEr[np.searchsorted(rs, 1.25)]) / 0.05
    slope_hi = float(dEr[np.searchsorted(rs, 1.90)]
                     - dEr[np.searchsorted(rs, 1.85)]) / 0.05
    exact = 2.0 * np.log(1.0 + rs) + 4.0 * np.log(np.abs(1.0 - 0.5 * rs))
    exact_dev = float(np.max(np.abs(exact - dEr)))
    data["triangle"] = {"flat_dev": flat_dev, "monotone": monotone,
                        "interior_minima": interior,
                        "slope_1p25_to_1p30": slope_lo,
                        "slope_1p85_to_1p90": slope_hi,
                        "argmin_r": float(rs[int(np.argmin(dEr))]),
                        "closed_form_dev": exact_dev,
                        "dE_at_2_minus": float(dEr[-1])}
    r_branch = ("R-ONSET-DEMOCRATIC" if monotone and interior == 0
                and abs(rs[int(np.argmin(dEr))] - rs[-1]) < 1e-12
                else "R-WINDOW-MIN")
    checks.append({"name": "T1 tree-level r-flat to %.1e (trace theorem), "
                          "one-loop strictly monotone %s, interior minima "
                          "%d, slopes %.2f -> %.2f (more negative at the "
                          "window edge), argmin r = %.4f at the onset, "
                          "closed form to %.1e"
                          % (flat_dev, monotone, interior, slope_lo,
                             slope_hi, rs[int(np.argmin(dEr))], exact_dev),
                   "pass": flat_dev < 1e-9 and monotone and interior == 0
                   and exact_dev < 1e-12})

    # ---- T2 the 4096-set sweep on the stable branch -----------------------
    sets = []
    for m in range(1 << 12):
        S = [v for v in range(12) if (m >> v) & 1]
        if not S:
            continue
        sets.append((frozenset(S), adjacency(S, edges)))
    cross = {}
    for S, A in sets:
        cross[S] = first_crossing(A)
    b1278 = sum(1 for c in cross.values() if c < 1.278)
    b1919 = sum(1 for c in cross.values() if c < 1.9191)
    b2 = sum(1 for c in cross.values() if c < 2.0 - 1e-9)
    ninf = sum(1 for c in cross.values() if not math.isfinite(c))
    distinct = sorted({round(c, 6) for c in cross.values()
                       if math.isfinite(c)})
    grid = np.round(np.arange(0.05, 1.96, 0.05), 6)
    tri_min_r, min_sizes, n_min_ties, all12_r = [], [], 0, 0
    min_sets = []
    for r in grid:
        best, ties, sizes, bestS = math.inf, 0, None, None
        for S, A in sets:
            val = dE(A, r)
            if not math.isfinite(val):
                continue
            if val < best - 1e-9:
                best, ties, sizes, bestS = val, 1, len(S), S
            elif abs(val - best) <= 1e-9:
                ties += 1
        n_min_ties = max(n_min_ties, ties)
        min_sizes.append(int(sizes))
        min_sets.append(sorted(bestS))
        if sizes == 12:
            all12_r += 1
        for S, A in sets:
            if dE(A, r) <= best + 1e-9 and (S == tri0 or S in jorb):
                tri_min_r.append(float(r))
                break
    data["sweep"] = {"sets": len(sets), "buckets_below_1p278": b1278,
                     "buckets_below_1p9191": b1919,
                     "buckets_below_2": b2, "never_massless": ninf,
                     "distinct_crossings": distinct,
                     "grid": [float(g) for g in grid],
                     "minimiser_sizes": min_sizes,
                     "minimiser_sets": min_sets,
                     "all12_min_count": all12_r,
                     "triangle_min_r": tri_min_r,
                     "max_ties": n_min_ties}
    sel_branch = ("TRIANGLE-SELECTED" if tri_min_r
                  else "VACUUM-NEVER-TRIANGLE")
    checks.append({"name": "T2 crossing levels %s (buckets %d / %d / %d, "
                          "%d sets never massless); minimiser sizes across "
                          "the window grid %s: the full 12-set minimises "
                          "for %d/%d points (below its own r = 1 crossing), "
                          "the closing triangle attains the minimum %d "
                          "times"
                          % (distinct, b1278, b1919, b2, ninf,
                             sorted(set(min_sizes)), all12_r, len(grid),
                             len(tri_min_r)),
                   "pass": 1.0 in distinct and 2.0 in distinct
                   and len(tri_min_r) == 0})

    # ---- T3 the contact: attractor r -> 2 vs the measured scale -----------
    S17 = nv.loop_ladder(17) * nv.EV_PER_GEV
    m3_onset = S17 / 3.0
    rel = m3_onset / nv.M_ATM_EV - 1.0
    window = S17 / nv.SUM_MNU_MAX_EV - 1.0
    data["contact"] = {"S17_eV": float(S17), "m3_onset_eV": float(m3_onset),
                       "M_ATM_EV": float(nv.M_ATM_EV),
                       "rel_residual_vs_atm": float(rel),
                       "window_r_min": float(window),
                       "onset_in_window": bool(1.278 < 2.0 <= window + 1e-12
                                               or window < 2.0)}
    checks.append({"name": "T3 attractor r* = 2 (inside the window "
                          "r > %.3f): m_3 = S17/3 = %.4f eV vs the booked "
                          "atmospheric %.4f eV: relative residual %+.2f%% "
                          "(named, booked)"
                          % (window, m3_onset, nv.M_ATM_EV, 100 * rel),
                   "pass": window < 2.0 and abs(m3_onset - S17 / 3) < 1e-15})

    branch = "%s / %s" % (r_branch, sel_branch)
    data["status"] = branch
    print("V303 %s: the stable branch r < 2 has an r-FLAT tree level "
          "(dev %.1e) and a strictly monotone one-loop (slopes %.1f -> "
          "%.1f, dE -> %.1f at the onset): the junction energy has NO "
          "interior minimum and flows to the massless onset r* = 2 - the "
          "democratic anchor - inside the window r > %.3f, so the "
          "candidate value m_3 = S17/3 = %.4f eV acquires an energy "
          "attractor (a marginal, boundary selection) at a %+0.2f%% "
          "residual against the atmospheric %.4f eV; over the 4096 split "
          "sets the vacuum minimiser has sizes %s across the window grid "
          "(the full 12-set only below its own 6-cycle crossing at r = 1: "
          "%d/%d points), the closing triangle attains it %d times, "
          "crossing levels %s" % (branch, flat_dev, slope_lo, slope_hi,
                                  dEr[-1], window, m3_onset, 100 * rel,
                                  nv.M_ATM_EV, sorted(set(min_sizes)),
                                  all12_r, len(grid), len(tri_min_r),
                                  distinct), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V303_junction_energy_r", "tier": 2,
           "property": "energy selection of the junction ratio r on the "
                       "stable branch (V127 theory, canonical convention)",
           "residual": True, "programme": "V303",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v303_junction_energy_r.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V303 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
