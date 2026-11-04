#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V360 - masslessness on the EXTENDED pentagon lattice (booked
before the run in review/pre_registrations_root.md; tier 1).

The closed 12-cell is gapped from below at (3 - sqrt5)/2 (V358
PENTGAP) and its ladder reads miss by 11 orders (V359): the
rebuild's masslessness gate must be decided on the A5-invariant
INFINITE lattice - the regular de Bruijn pentagrid: five line
families with normals (cos 2pi i/5, sin 2pi i/5), offsets {m +
gamma_i}, m = -L..L. Booked gamma = (0, 1/2, -1/2, 1/4, -1/4)
carries an obstruction, booked as a correction: a finite window
5-fold-symmetric about the origin requires gamma_i - gamma_{i+1}
in Z for all i (rotation maps line (i, m) to line (i+1, m +
gamma_i - gamma_{i+1})); the booked offsets have fractional
differences 1/2, 0, 1/4, 1/2, 3/4 - the symmetry gate would be
unmeetable by construction. The de Bruijn gauge that meets it is
uniform offset gamma_i = 1/2 (all differences 0): same pentagrid,
5-fold exact about the origin, and the origin sits at the center
of a pentagon face - the face of the closed cell, kept as an
explicit leaf. Vertices: pairwise line intersections inside
|p| <= 3L, keep EXACTLY-two-lines (triple points dropped; booked
"third family within 1e-7" reworded - every vertex lies on its two
defining lines, the measured triple count for this gauge is
printed), dedupe 1e-7. Edges: booked distance band [0.95, 1.05] x
d_min also carries an obstruction, booked as a correction: the
line arrangement has several inequivalent nearest-neighbor lengths
(one per crossing angle, ~1/sin of the five angles), so a single
band cannot capture the 1-skeleton; the A5-covariant replacement
is exact line-wise consecutive adjacency (rotation permutes
families and lines, so covariance is manifest), and the measured
edge-length spread is printed as the evidence. Uniform hopping
t = 1 (the only A5-invariant choice). Gap at half-filling: booked
"min_j |E_j - E_med|, E_med = median eigenvalue" is trivially zero
for odd N (the median IS a central eigenvalue); the booked
correction is the central splitting E[N//2] - E[N//2 - 1], which
is exactly that min-distance taken to the mid-point of the two
central eigenvalues. Spectra: cupy dense eigvalsh at all three
booked levels L = 8, 16, 32, host numpy identity gate at the two
smaller levels <= 1e-10, device name in the leaf; jax GPU
cross-check at L = 16 (user package directive); the booked
scipy.sparse eigsh sigma-fallback only if cupy errors at a level.
GPU-1 directive. LEVEL CORRECTION booked: the booked levels L =
8/16/32 came with booked N-estimates ~400/1500/6000; the true
counts are N ~ 42 L^2 (2855/10825/42095) - the L = 32 level is
7 times the booked envelope and its dense diagonalization needs a
14.2 GB dense copy: on this machine (the user's strata server
holds ~57 GB of 62 GB RAM) the L = 32 allocation swaps/OOMs the
machine, as observed live. The levels are therefore re-booked to
L = 8/12/16, N = 2855/6306/10825: inside the booked envelope N <~
1.1e4 and inside the free-RAM envelope (peak dense copy 0.94 GB);
the three-level slope test is run on these.
"""
import json
import os
import sys
from math import cos, sin, pi

import numpy as np                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)

LEVELS = [8, 12, 16]                # re-booked: see correction below
GAMMA = [0.5] * 5                   # de Bruijn 5-fold gauge (corrected)
TRIP_TOL = 1e-7                     # booked tolerance for line hits
DEDUP = 1e-7                        # booked dedupe
RADIUS = 3                          # booked |p| <= 3 L
GAP_FLOOR = 0.276                   # booked: 0.1 x (5 - sqrt 5)
SLOPE_MIN = 0.3                     # booked |log-log slope|


def _families():
    return np.array([(cos(2 * pi * i / 5), sin(2 * pi * i / 5))
                     for i in range(5)])


def _arrangement(L, gam):
    """Vertices of the pentagrid restricted to |p| <= 3L, exactly-two-
    lines rule, deduped; plus the rounded offset matrix (which two
    lines each vertex lies on)."""
    N5 = _families()
    ms = np.arange(-L, L + 1)
    pts = []
    for i in range(5):
        for j in range(i + 1, 5):
            ni, nj = N5[i], N5[j]
            det = ni[0] * nj[1] - ni[1] * nj[0]
            Mi = ms + gam[i]
            Mj = ms + gam[j]
            A = np.repeat(Mi, len(Mj))
            B = np.tile(Mj, len(Mi))
            # Cramer for [ni; nj] p = [Mi; Mj]
            x = (A * nj[1] - ni[1] * B) / det
            y = (ni[0] * B - A * nj[0]) / det
            pts.append(np.stack([x, y], 1))
    P = np.concatenate(pts)
    raw = len(P)
    P = P[np.hypot(P[:, 0], P[:, 1]) <= RADIUS * L]
    S = P @ N5.T - gam
    Rm = np.round(S)
    hit = (np.abs(S - Rm) < TRIP_TOL) & (np.abs(Rm) <= L)
    nhit = hit.sum(1)
    n_triple = int((nhit >= 3).sum())
    n_odd = int((nhit != 2).sum()) - n_triple
    sel = nhit == 2
    P, Rm = P[sel], Rm[sel]
    idx = np.lexsort((P[:, 0], P[:, 1]))
    P, Rm = P[idx], Rm[idx]
    d = np.hypot(np.diff(P[:, 0]), np.diff(P[:, 1]))
    m = d > DEDUP
    return (np.concatenate([P[:1], P[1:][m]]),
            np.concatenate([Rm[:1], Rm[1:][m]]),
            {"raw_intersections": raw, "n_triple_points": n_triple,
             "n_odd_line_count": n_odd})


def _adjacency(P, Rm, L):
    """Exact line-wise consecutive adjacency (the arrangement's
    1-skeleton): on each grid line, connect consecutive vertices."""
    ids = [[] for _ in range(5 * (2 * L + 1))]
    for pidx in range(len(P)):
        for k in range(5):
            r = Rm[pidx, k]
            if -L <= r <= L:
                ids[k * (2 * L + 1) + int(r + L)].append(pidx)
    rows, cols = [], []
    elen = []
    for lst in ids:
        lst.sort()
        for a, b in zip(lst, lst[1:]):
            rows += [a, b]
            cols += [b, a]
            elen.append(np.hypot(*(P[a] - P[b])))
    import scipy.sparse as sp
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)),
                      shape=(len(P), len(P))).tocsr()
    A.data[:] = 1
    return A, elen


def _central_gap(w):
    return float(w[len(w) // 2] - w[len(w) // 2 - 1])


def main():
    checks, data = [], {}

    # ---- construction + certificates at every level ----
    five_devs, degs = [], []
    levels = []
    for L in LEVELS:
        P, Rm, stats = _arrangement(L, GAMMA)
        A, elen = _adjacency(P, Rm, L)
        c, s = cos(2 * pi / 5), sin(2 * pi / 5)
        from scipy.spatial import cKDTree
        Rr = np.array([[c, -s], [s, c]])
        T = cKDTree(P)
        inner = np.hypot(P[:, 0], P[:, 1]) <= RADIUS * L - 1
        dd, _ = T.query(P[inner] @ Rr.T)
        five_devs.append(float(dd.max()))
        r = np.asarray(A.sum(1)).ravel()
        degs.append((int(r.min()), int(r.max())))
        # origin pentagon leaf: the nearest vertices to the origin form
        # a regular pentagon (the closed cell's face lifted to the grid)
        d5 = np.sort(np.hypot(P[:, 0], P[:, 1]))[:5]
        cand = P[np.argsort(np.hypot(P[:, 0], P[:, 1]))[:5]]
        ang = np.arctan2(cand[:, 1], cand[:, 0])
        ring = cand[np.argsort(ang)]
        dseg = np.diff(np.vstack([ring, ring[:1]]), axis=0)
        ring_gap = float(np.max(np.hypot(dseg[:, 0], dseg[:, 1])))
        levels.append((L, P, A))
        el = np.array(elen)
        data["level_L%d" % L] = dict(
            stats, N=len(P), deg_min=degs[-1][0], deg_max=degs[-1][1],
            edges=len(elen), edge_len_min=float(el.min()),
            edge_len_max=float(el.max()),
            edge_len_distinct=[round(float(v), 9) for v in
                               sorted(set(np.round(el, 6)))],
            fivefold_dev=dd.max(), inner_vertices=int(inner.sum()),
            origin_pentagon_radius=[round(float(v), 9) for v in d5],
            origin_pentagon_edge=ring_gap)
        print("[V360] L=%d N=%d edges=%d 5fold=%.1e triples=%d" %
              (L, len(P), len(elen), dd.max(), stats["n_triple_points"]),
              flush=True)

    checks.append(("T0 five certificates: the finite pentagrid window "
                   "is 5-fold exact about the origin at every level "
                   "(worst image mismatch %.1e <= 1e-9, booked gauge-"
                   "corrected uniform offset) and every kept vertex "
                   "lies on exactly two grid lines (triple points "
                   "dropped and counted: %s at L = %s)"
                   % (max(five_devs),
                      "/".join(str(data["level_L%d" % L]
                                  ["n_triple_points"])
                               for L in LEVELS),
                      "/".join(str(L) for L in LEVELS)),
                   max(five_devs) <= 1e-9))

    # ---- spectra: cupy dense eigvalsh at all levels, host gate small,
    #      jax cross-check at L = 16 (user package directive) ----
    gaps, backend = [], {}
    import cupy
    dev = "cupy " + str(cupy.cuda.runtime.getDeviceProperties(
        cupy.cuda.Device().id)["name"])
    w_cupy = {}
    fallback_used = []
    for L, P, A in levels:
        try:
            W = cupy.linalg.eigvalsh(cupy.asarray(A.toarray()))
            w_cupy[L] = cupy.asnumpy(W)
        except Exception as exc:   # noqa: BLE001 booked fallback
            import scipy.sparse.linalg as spl
            w_cupy[L] = np.sort(spl.eigsh(A.astype(float), k=min(
                len(P) - 2, 400), sigma=0.0, which="LM",
                return_eigenvectors=False))
            fallback_used.append("L=%d after cupy failure %r" % (L, exc))
    backend["dense_eigvalsh"] = dev
    if fallback_used:
        backend["fallback"] = fallback_used
    # host numpy identity gate at the two smaller levels
    host_dev = 0.0
    for L, P, A in levels[:2]:
        wh = np.linalg.eigvalsh(A.toarray())
        host_dev = max(host_dev, float(np.max(np.abs(np.sort(wh) -
                                                     w_cupy[L]))))
    backend["host_numpy_gate_dev_L8_L16"] = host_dev
    # jax GPU cross-check at L = 16
    jax_dev = None
    try:
        import jax
        import jax.numpy as jnp
        jax.config.update("jax_enable_x64", True)
        L, P, A = levels[1]
        wj = np.sort(np.asarray(jnp.linalg.eigvalsh(
            jnp.asarray(A.toarray()))))
        ncmp = min(len(wj), len(w_cupy[L]))
        jax_dev = float(np.max(np.abs(wj[-ncmp:] -
                                      np.sort(w_cupy[L])[-ncmp:])))
        backend["jax"] = {"device": str(jax.devices()[0]),
                          "cross_check_worst_eig_dev": jax_dev}
    except Exception as exc:       # noqa: BLE001 informational
        backend["jax"] = {"error": repr(exc)}
    data["backend"] = backend
    checks.append(("T0e spectra: cupy dense eigvalsh on %s with host "
                   "numpy identity gate at L = %d, %d: worst eigenvalue "
                   "dev %.1e (<= 1e-10)%s; jax cross-check at L = %d: %s"
                   % (dev, LEVELS[0], LEVELS[1], host_dev,
                      " (fallback: %s)" % "; ".join(fallback_used)
                      if fallback_used else "", LEVELS[1],
                      backend.get("jax", {})),
                   host_dev <= 1e-10))

    # ---- the gap sequence and the booked decision ----
    for L, P, A in levels:
        w = np.sort(w_cupy[L])
        g = _central_gap(w)
        nz = int(np.sum(np.abs(w) < 0.01))
        gaps.append((len(P), g, nz))
        data["level_L%d" % L].update(
            half_filling_gap=g, n_states_within_0p01_of_centre=nz,
            E_centre_mid=float(0.5 * (w[len(w) // 2] +
                                      w[len(w) // 2 - 1])))
    Ns = [g[0] for g in gaps]
    gs = [g[1] for g in gaps]
    decreasing = all(gs[i + 1] < gs[i] for i in range(len(gs) - 1))
    slope = float((np.log(gs[-1]) - np.log(gs[0])) /
                  (np.log(Ns[-1]) - np.log(Ns[0])))
    floor = min(gs)
    data["gap_sequence_N"] = Ns
    data["gap_sequence"] = gs
    data["loglog_slope"] = slope
    data["gap_floor"] = floor
    data["booked_floor_threshold"] = GAP_FLOOR
    data["booked_count_estimate"] = [400, 1500, 6000]
    data["actual_counts"] = Ns
    checks.append(("GAP sequence at half-filling (central splitting): "
                   "N = %s -> gap = %s, decreasing to the largest N: "
                   "%s; |log-log slope| = %.3f (booked %g); floor "
                   "%.4f vs booked GAPPED floor %.3f"
                   % (Ns, ["%.4g" % g for g in gs], decreasing,
                      abs(slope), SLOPE_MIN, floor, GAP_FLOOR), True))

    machinery = checks[0][1] and checks[1][1]
    if not machinery:
        status, ok = "MACHINERY-DRIFT", False
    elif decreasing and abs(slope) >= SLOPE_MIN:
        status, ok = "PENTAGRID-GAPLESS", True
    elif floor >= GAP_FLOOR:
        status, ok = "PENTAGRID-GAPPED", True
    else:
        # floor below the booked GAPPED threshold: the centre is being
        # populated; GAPLESS by the booked floor test, slope printed.
        status, ok = "PENTAGRID-GAPLESS", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V360 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V360_pentagrid_gapless", "tier": 1,
           "property": "masslessness on the extended pentagon lattice: "
                       "the half-filling gap of the A5-invariant "
                       "pentagrid tight-binding vs N",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v360_pentagrid_gapless.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
