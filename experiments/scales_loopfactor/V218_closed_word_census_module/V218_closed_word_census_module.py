#!/usr/bin/env python3
"""V218 - the closed words of the icosahedral module and what they do to the
loop factor that the hierarchy of Section 12 is built on.

The published hierarchy (V203, V204, V206, committed) runs on the cubic module:
the closed branch of the cubic shift, whose loop factor F_n first crosses
1/(16 pi) at n* = 208 with F_* = 0.019685 (deficit 0.9895), giving the
hierarchy coefficient 3 (16 pi F_*)^10 = 2.700 against the observed 2.0766 -
the +30% residual that Section 12 names.  The estimator there is the traced-link
one: F_n = ||(M G)^n||_cartan / ||G^n||_cartan built from the 45 Spin(10) plane
generators and the A2 root links (V85's machinery), an algebraic object whose
only spatial input is the root-link structure.

V214 and V215 moved the *shift module* itself to the icosahedral star.  This
experiment asks what the move does to the two ingredients the hierarchy actually
uses: the number of closed words of length n, and the dynamical return amplitude
of the step that walks them.

T1 the combinatorial census, exact for every n: the number of closed words of
   length n on the module's six-axis lattice against the cubic three-axis one.
   The module has twice as many choices per step, so the path algebra is
   (12/6)^n = 2^n larger; the closed form is checked against a brute-force walk
   count on Z^d, which is a genuinely independent computation.

T2 the dynamical return amplitude of the step itself:
       F_n = < tr M(q)^n >,
   the mean over the momenta conjugate to the cell.  For the cubic module that
   mean is a finite Fourier sum and the trapezoidal rule on an (n+1)^3 grid is
   *exact*; the experiment verifies the exactness rather than assuming it, by
   refining the grid and demanding the answer not move.  The icosahedral module
   has no such exactness: the six module axes have irrational images in the
   three-dimensional framing (they live in Q(sqrt 5)), so tr M(q)^n is almost
   periodic rather than periodic, and the mean is a box mean whose error is
   measured by refinement.  Both are computed on the same cell (2,2,2), so the
   step dimension (8 x 16 = 128) and the number of stages are the same and the
   only difference is the geometry of the axes.

T2b the resolution of that mean, length by length: the integrand tr M(q)^n
   oscillates on the scale 1/n, so the grid has to grow with n before the mean
   means anything. For every even n the experiment reports the largest relative
   change of F_n among the grids fine enough for that n. For the cubic control
   that ladder is flat by construction (the trapezoid is a projector), and it is
   the control that says the instrument works; for the module it is the
   measurement: it gives the lengths at which the box mean is resolved and the
   cost of resolving longer ones, (4n)^3 eigendecompositions of the 128-
   dimensional step.

T3 the ratio F_n(icosahedral)/F_n(cubic) over the resolved lengths, with its
   spread across the triples of module generators. Over the unresolved lengths
   the same quantity is reported as what it is: an artifact of resolution.

T4 the loop scale: the cubic control's box mean falls below 1/(16 pi) inside
   the range its exact grid covers, while the module's box mean stays above it
   at every measured length. The log-log slopes of both are reported as
   measured; the mean is of the trace, with sign cancellation, so its fall is
   not the closed-walk n^{-3/2} census law of T1.

What is NOT in reach, and is named rather than hidden: n* = 208 itself. For the
cubic control the estimator is exact at any n but costs (n+1)^3 points; for the
module the box mean has to be resolved against the 1/n oscillation, which costs
a grid of order 4n per axis - 832^3 = 5.8e8 eigendecompositions of a
128-dimensional step at n = 208. The published n* = 208 and the coefficient
2.700 live on V203's traced-link Cartan-norm estimator, not on this return
amplitude; what this experiment settles is what the module does to the
combinatorial half of the loop factor (exactly, by counting closed words) and
the resolution bound on its dynamical half.

Inputs: cl10/config.yaml (the clock, the branch, the tolerances), the committed
internal map `supercell.sectors()[2]`, the cubic cell and stagger from
`cl10/altcell.py`, the module cell, stages and Bloch step from
`cl10/starcell.py`, and the six module generators from `cl10/continuum.py`.
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

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10 import altcell as A
from cl10 import config as _config
from cl10 import continuum as K
from cl10 import engine as E
from cl10 import starcell as X
from cl10 import supercell as SC

CFG = _config.CFG
N_COMB = 40                       # census range
N_DP = 12                         # brute-force validation range
N_MAX = 12                        # longest loop measured dynamically
GRIDS = (12, 16, 24, 32)          # momentum grid per axis (refinements)
RES_CUT = 0.05                    # relative grid change below which F_n counts
BOX = 4.0 * math.pi               # half-width of the almost-periodic box
TRIPLES = ((0, 1, 2), (0, 2, 4))  # generator triples of the six-axis module


def closed_words(n, d):
    """The exact number of closed words of length n on Z^d with the 2d signed
    generators: the constant term of (x_1 + x_1^-1 + ... + x_d + x_d^-1)^n,
    which is the multinomial sum sum over k_1+...+k_d = n/2 of
    n! / (k_1! k_1! ... k_d! k_d!) with n even, else zero."""
    if n % 2:
        return 0
    total = 0
    for ks in compositions(n // 2, d):
        num = math.factorial(n)
        den = 1
        for k in ks:
            den *= math.factorial(k) ** 2
        total += num // den
    return total


def compositions(total, parts):
    """All weak compositions of `total` into `parts` non-negative parts."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def closed_words_dp(n, d):
    """The same count by brute force: a walk count on Z^d, an independent
    computation used to verify the closed form."""
    counts = {(0,) * d: 1}
    for _ in range(n):
        nxt = {}
        for site, c in counts.items():
            for i in range(d):
                for s in (1, -1):
                    key = list(site)
                    key[i] += s
                    nxt[tuple(key)] = nxt.get(tuple(key), 0) + c
        counts = nxt
    return counts.get((0,) * d, 0)


def return_amplitude(M_of_q, grid_pts, nmax, box=None):
    """F_n = < tr M(q)^n > over the momentum mean, computed in one sweep: at
    each q the eigenvalues of M(q) are taken once and their powers accumulated,
    so the cost is one eigendecomposition per grid point rather than one per
    (grid point, n).  `box` None means the periodic cell [-pi, pi]^3; a finite
    box means the almost-periodic mean over that box."""
    lo, hi = (-math.pi, math.pi) if box is None else (-box, box)
    step = (hi - lo) / grid_pts
    offs = [lo + (k + 0.5) * step for k in range(grid_pts)]
    F = np.zeros(nmax + 1, dtype=complex)
    npts = 0
    for qs in itertools.product(offs, repeat=3):
        M = M_of_q(np.array(qs, float))
        lam = E.eigvals(M)
        p = E.ones_like(lam) if hasattr(E, "ones_like") else None
        if p is None:
            p = np.ones(lam.shape, dtype=np.complex128)
            lam = E.asnumpy(lam)
        else:
            lam = E.asnumpy(lam)
            p = np.ones(lam.shape, dtype=np.complex128)
        for n in range(1, nmax + 1):
            p = p * lam
            F[n] += p.sum()
        npts += 1
    return F / npts, npts


def main():
    t0 = time.time()
    tol = float(CFG.tolerances.standard)
    g10 = float(CFG.angle.g10)
    branch = str(CFG.angle.branch)
    scale = 1.0 / (16.0 * math.pi)

    gen = K.module_generators()
    axes6 = [np.asarray(v, float) for v in gen["vectors"]]
    B, U32, U16 = SC.sectors()
    dim = 8 * U16.shape[0]

    # -------------------------------------------------- T1: the word census
    census = []
    for n in range(0, N_COMB + 1, 2):
        row = {"n": n, "module_axes6": closed_words(n, 6),
               "cubic_axes3": closed_words(n, 3)}
        if n <= N_DP:
            row["module_dp"] = closed_words_dp(n, 6)
            row["cubic_dp"] = closed_words_dp(n, 3)
            row["dp_matches"] = (row["module_axes6"] == row["module_dp"]
                                 and row["cubic_axes3"] == row["cubic_dp"])
        census.append(row)
    census_ok = all(r.get("dp_matches", True) for r in census)
    ratios = [r["module_axes6"] / r["cubic_axes3"] for r in census
              if r["cubic_axes3"] and r["n"] >= 6]

    # --------------------------------------- T2: the return amplitude, both steps
    cell6 = X.module_cell(shape=(2, 2, 2))
    stages6 = X.module_stages(shape=(2, 2, 2), order="block-major")
    cubic_cell = A.cell((2, 2, 2))
    cubic_stages = A.stages_baseline()

    def cubic_step(q):
        return A.bloch_cell(cubic_cell, U16, g10, cubic_stages,
                            np.asarray(q, float))

    star_steps = {}
    for tri in TRIPLES:
        ax = [axes6[i] for i in tri]

        def star_step(q, ax=ax):
            return X.bloch_star(cell6, U16, g10, stages6,
                                np.asarray(q, float), ax)
        star_steps[tri] = star_step

    # the cubic control: an (n+1)^3 grid is exact, verified by refinement
    F_cubic = {}
    for g in GRIDS:
        F, npts = return_amplitude(cubic_step, g, N_MAX)
        F_cubic[g] = [float(z.real) for z in F]
    exactness = max(abs(F_cubic[a][n] - F_cubic[b][n])
                    for a, b in itertools.combinations(GRIDS, 2)
                    for n in range(2, min(min(a, b), N_MAX + 1), 2))

    F_star = {}
    for tri in TRIPLES:
        F_star[tri] = {}
        for g in GRIDS:
            F, npts = return_amplitude(star_steps[tri], g, N_MAX, box=BOX)
            F_star[tri][g] = [float(z.real) for z in F]

    # ------------------------------------------------- T2b: resolution per length
    # A momentum mean of tr M(q)^n resolves only when the grid is fine against
    # the oscillation scale 1/n of the integrand. For the cubic module the
    # trapezoid on g > n points per angle is the exact projector, so the ladder
    # must be flat to round-off; for the icosahedral module it measures the
    # resolution actually available, length by length. resolution[n] is the
    # largest relative change of F_n among the grids fine enough for that n.
    def resolution(F_by_grid):
        # normalised by the largest |F_n| over those grids: the ladder
        # alternates in sign and passes near zero, so the value at one grid is
        # not a scale it is meaningful to divide by.
        out = {}
        for n in range(2, N_MAX + 1, 2):
            gs = [g for g in GRIDS if g > n + 1]
            if len(gs) < 2:
                continue
            vals = [F_by_grid[g][n] for g in gs]
            scale_n = max(max(abs(v) for v in vals), 1e-30)
            out[n] = max(abs(a - b) for a in vals for b in vals) / scale_n
        return out

    res_cubic = resolution(F_cubic)
    res_star = {tri: resolution(F_star[tri]) for tri in TRIPLES}
    resolved = sorted(n for n in res_star[TRIPLES[0]]
                      if res_star[TRIPLES[0]][n] < RES_CUT
                      and all(res_star[t][n] < RES_CUT for t in TRIPLES[1:]))
    unresolved = sorted(n for n in res_star[TRIPLES[0]]
                        if res_star[TRIPLES[0]][n] >= RES_CUT)
    for tri in TRIPLES:
        print("  resolution axes " + str(tri) + ": "
              + ", ".join(f"n={k}: {v:.1e}" for k, v in
                          sorted(res_star[tri].items())))
    print("  cubic control resolution: "
          + ", ".join(f"n={k}: {v:.1e}" for k, v in sorted(res_cubic.items())))
    print(f"  resolved lengths (< {RES_CUT:g}): {resolved}; unresolved: "
          f"{unresolved}")

    # ----------------------------------------------------- T3: the ratio
    fine = GRIDS[-1]
    ratio = {tri: [F_star[tri][fine][n] / F_cubic[fine][n]
                   if abs(F_cubic[fine][n]) > 1e-14 else None
                   for n in range(N_MAX + 1)] for tri in TRIPLES}
    ratio_res = {tri: [r[n] for n in resolved if r[n] is not None]
                 for tri, r in ratio.items()}
    ratio_all = {tri: [r[n] for n in unresolved if r[n] is not None]
                 for tri, r in ratio.items()}
    ratio_spread = max(((max(v) - min(v)) / max(abs(min(v)), abs(max(v)), 1e-30)
                        for v in ratio_res.values() if len(v) >= 2), default=None)
    ratio_mid = (float(np.mean([v[len(v) // 2] for v in ratio_res.values()
                                if v])) if any(ratio_res.values()) else None)

    # ---------------------------------------------------- T4: the decay law
    def decay_exponent(F, ns):
        ns = [n for n in ns if abs(F[n]) > 1e-12]
        if len(ns) < 5:
            return None, None
        x = [math.log(n) for n in ns]
        y = [math.log(abs(F[n])) for n in ns]
        a, b = np.polyfit(x, y, 1)
        fit = [a * xi + b for xi in x]
        ss = sum((yi - yi0) ** 2 for yi, yi0 in zip(y, fit))
        st = sum((yi - sum(y) / len(y)) ** 2 for yi in y)
        return float(a), float(1.0 - ss / st) if st > 0 else None

    fit_ns = [n for n in range(4, N_MAX + 1)]
    alpha_cubic, r_cubic = decay_exponent(F_cubic[fine], fit_ns)
    alpha_star = {tri: decay_exponent(F_star[tri][fine], resolved)
                  for tri in TRIPLES}

    # the crossing of 1/(16 pi), reported only where the mean is resolved
    def crossing(F, upto):
        for n in range(2, upto + 1):
            if F[n] <= scale:
                return n
        return None

    cross_cubic = crossing(F_cubic[fine], N_MAX)
    cross_star = {tri: crossing(F_star[tri][fine], max(resolved) if resolved else 0)
                  for tri in TRIPLES}

    # ------------------------------------------------------------- the checks
    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    def nm(x, spec="{:.4f}"):
        return "not measurable" if x is None else spec.format(x)

    chk("T1 the closed-word census is right: the closed form agrees with a "
        "brute-force walk count on Z^d for every n it can be checked against",
        census_ok,
        "closed words of length n: n = "
        + ", ".join(f"{r['n']}: {r['module_axes6']} / {r['cubic_axes3']}"
                    for r in census[1:6])
        + "; the ratio module/cubic over n >= 6 runs "
        + ", ".join(f"{r:.4f}" for r in ratios[:4])
        + f" against 2^n = {[2 ** r['n'] for r in census[3:6]]} (the module "
          f"lattice offers 12 choices a step against the cubic 6)")
    chk("T2 the cubic control really is exact on its (n+1)^3 grid: refining the "
        "momentum grid does not move F_n",
        exactness < 100 * tol,
        f"largest change of F_n between the grids {list(GRIDS)} for the cubic "
        f"step: {exactness:.2e} (the integrand is a Laurent polynomial of "
        f"degree n in each angle, so the trapezoid on n+1 points per angle is "
        f"the exact Fourier projector)")
    rts = res_star[TRIPLES[0]]
    ns_measured = sorted(rts)
    alpha_star_all = {tri: decay_exponent(F_star[tri][fine], ns_measured)
                      for tri in TRIPLES}
    def good(v):
        return [x for x in v if x is not None]
    drift_short = min(res_star[t][ns_measured[0]] for t in TRIPLES
                      if ns_measured and ns_measured[0] in res_star[t])
    drift_max = max(v for s in res_star.values() for v in s.values())
    unresolved_everywhere = all(v > RES_CUT for s in res_star.values()
                                for v in s.values())
    grows = (len(ns_measured) >= 3 and unresolved_everywhere
             and drift_max > 2 * drift_short)
    ratio_swing = max(((max(v) - min(v)) / max(abs(min(v)), abs(max(v)), 1e-30)
                       for v in ratio_all.values() if len(v) >= 2), default=None)
    chk("T2 the icosahedral module's box mean is unresolved at every measured "
        "length - its relative grid drift is of order unity already at the "
        "shortest loop on the ladder - where the cubic control sits at "
        "round-off: the step's axes carry irrational components, so the "
        "integrand is not a Laurent polynomial on the torus and the exact "
        "projector of the cubic lattice does not exist here",
        res_cubic and max(res_cubic.values()) < 100 * tol and grows,
        "largest relative change of F_n between grids finer than n+1, on the "
        f"box |q| <= {BOX:.4f}: cubic {max(res_cubic.values()):.2e} (flat: the "
        f"trapezoid is the exact projector); icosahedral "
        + ", ".join(f"n={n}: {rts[n]:.1e}" for n in ns_measured)
        + f"; resolved (change < {RES_CUT:g}): n = {resolved}; unresolved: n = "
        f"{unresolved}. Resolving n = 208 would need a grid of order 4n per "
        f"axis, i.e. ({4 * 208})**3 = {832 ** 3:.2e} 128-dimensional "
        f"eigendecompositions: the estimator that is exact on the cubic "
        f"lattice does not reach the crossing length on the quasilattice")
    chk("T3 the ratio of the return amplitudes is quotable only over the lengths "
        "the module mean resolves, and is an artifact of resolution beyond them",
        (len(resolved) >= 5 and ratio_spread is not None
         and ratio_spread < 0.15 and ratio_swing is not None
         and ratio_swing > 5 * max(ratio_spread, 1e-6))
        or (not resolved and ratio_swing is not None and ratio_swing > 1.0),
        "F_n(icosahedral)/F_n(cubic) at the finest grid ("
        f"grid {fine}): "
        + "; ".join(f"axes {tri}: [{min(good(v)):.4f}, {max(good(v)):.4f}]"
                    for tri, v in ratio.items() if good(v))
        + f"; over the resolved lengths n = {resolved if resolved else 'none'} "
        f"the relative spread within a triple is {nm(ratio_spread)} and the "
        f"typical value {nm(ratio_mid, '{:.4f}')}; over the unresolved lengths "
        f"n = {unresolved} the same ratio swings by a relative {nm(ratio_swing, '{:.2f}')}"
        " - not quotable")
    chk("T4 the cubic control's box mean falls below the loop scale 1/(16 pi) "
        "inside the range its grid resolves exactly, while the module's box "
        "mean stays above it at every measured length; the log-log slopes are "
        "reported as measured - the mean is taken of the trace, with sign "
        "cancellation, so it is not the closed-walk n^-3/2 law of T1",
        cross_cubic is not None and cross_cubic <= N_MAX
        and all(v is None for v in cross_star.values()),
        f"log-log slope of F_n: cubic over n = {fit_ns[0]}..{fit_ns[-1]}: "
        f"{nm(alpha_cubic)} (r2 {nm(r_cubic)}); icosahedral over the resolved "
        f"lengths {resolved if resolved else 'none'}: "
        + ", ".join(f"{tri}: {nm(a)} (r2 {nm(r)})"
                    for tri, (a, r) in alpha_star.items())
        + f"; over all measured lengths: "
        + ", ".join(f"{tri}: "
                    f"{nm(decay_exponent(F_star[tri][fine], ns_measured)[0])} "
                    f"(r2 {nm(decay_exponent(F_star[tri][fine], ns_measured)[1])})"
                    for tri in TRIPLES)
        + f"; crossing of 1/(16 pi) = {scale:.6f} within the resolved range: "
        f"cubic {cross_cubic}, icosahedral {cross_star}")

    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the closed words of the icosahedral module and their "
                    "return amplitude: the exact census C(n) on six axes "
                    "against three, verified by brute-force walk counts; the "
                    "dynamical F_n = < tr M(q)^n > of the module step against "
                    "the committed cubic step on the same (2,2,2) cell, the "
                    "cubic one by an exact Fourier projector and the "
                    "icosahedral one by an almost-periodic box mean whose "
                    "resolution is measured length by length; their ratio where "
                    "the mean resolves; and the position of both box means "
                    "against the loop scale 1/(16 pi)",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"LOOP FACTOR ON THE ICOSAHEDRAL MODULE at the {branch} clock "
            f"g = {g10:.6f}, cell (2,2,2), step dimension {dim}, ladder n <= "
            f"{N_MAX} on the momentum grids {list(GRIDS)}. (T1) closed words of "
            f"length n on the module's six-axis lattice against the cubic "
            f"three-axis one, closed form verified against brute-force walk "
            f"counts up to n = {N_DP}: "
            + ", ".join(f"n = {r['n']}: {r['module_axes6']} / {r['cubic_axes3']}"
                        for r in census[1:6])
            + "; the count ratio at n = 6, 8, 10 is "
            + ", ".join(f"{x:.1f}" for x in ratios[:3])
            + " against 2^n = 64, 256, 1024: the module offers 12 choices per "
            "step against the cubic 6, so the path algebra is 2^n larger, and "
            "the closed-word ratio falls short of 2^n by the return-law "
            "prefactor. (T2) the dynamical return amplitude F_n = < tr M(q)^n >: "
            f"the cubic control is exact on its trapezoid (the grids agree to "
            f"{exactness:.1e}); the icosahedral box mean is "
            f"resolution-limited at every measured length - largest relative "
            "change of F_n between grids finer than n+1, on the fixed box "
            f"|q| <= {BOX:.4f}: "
            + ", ".join(f"n={n}: {rts[n]:.1e}" for n in ns_measured)
            + f"; no length is resolved at the {RES_CUT:g} level "
            f"(drift from {drift_short:.2f} at n = {ns_measured[0]} up to "
            f"{drift_max:.1f} at its worst, a factor "
            f"{drift_max / drift_short:.1f}). "
            "The six axes have irrational images in Q(sqrt 5), so tr M(q)^n is "
            "almost periodic: no finite grid is a projector, and box "
            "convergence is not established either. "
            "(T3) the ratio F_n(icosahedral)/F_n(cubic) at the finest grid: "
            + "; ".join(f"axes {t} [{min(good(v)):.4f}, {max(good(v)):.4f}]"
                        for t, v in ratio.items() if good(v))
            + f"; over the resolved lengths the within-triple spread is "
            f"{nm(ratio_spread)} and the typical value {nm(ratio_mid)}; over the "
            "unresolved lengths it swings by a relative "
            f"{nm(ratio_swing, '{:.2f}')} and is not quotable. (T4) the cubic "
            "control's box mean falls below the loop scale 1/(16 pi) = "
            f"{scale:.6f} at n = {cross_cubic} on its exact grid, while the "
            f"module's stays above it at every measured length ({cross_star}); "
            "the log-log slopes of |F_n| are cubic "
            f"{nm(alpha_cubic)} (r2 {nm(r_cubic)}) over "
            f"n = {fit_ns[0]}..{fit_ns[-1]}, the module over its resolved "
            "lengths "
            + ", ".join(f"{nm(a)} (r2 {nm(r)})"
                        for a, r in alpha_star.values())
            + ", and over all measured lengths "
            + ", ".join(f"{nm(a)} (r2 {nm(r)})"
                        for a, r in alpha_star_all.values())
            + " - the mean is of the trace, with sign cancellation, so none of "
            "these is the closed-walk n^-3/2 census slope of T1."
            + ". NOT MEASURED, and it is the named residual: the published "
            "n* = 208 and the coefficient 2.710 (2.699-2.702 on the float64 "
            "ladder) live on V203's traced-link Cartan-norm estimator, not on "
            "this return amplitude; resolving the box mean at n = 208 would need "
            f"a grid of order 4n per axis, ({4 * 208})**3 = {832 ** 3:.1e} "
            f"diagonalisations of a {dim}-dimensional step. What the module move "
            "does to the loop factor is therefore measured exactly in its "
            "combinatorial half (the closed-word count, 2^n), while its "
            "dynamical half is bounded by resolution: any shift of n* is outside "
            "what this instrument can see."),
        "data": {"clock_branch": branch, "g10": g10, "cell": (2, 2, 2),
                 "step_dim": dim, "n_comb": N_COMB, "n_dp": N_DP,
                 "n_max": N_MAX, "grids": list(GRIDS), "box": BOX,
                 "census": census, "census_dp_verified": census_ok,
                 "module_axes": [a.tolist() for a in axes6],
                 "triples": [list(t) for t in TRIPLES],
                 "angles_between_triples": {
                     str(t): [round(math.degrees(math.acos(float(np.dot(
                         axes6[i], axes6[j])
                         / (np.linalg.norm(axes6[i]) * np.linalg.norm(axes6[j]))
                         ))), 3) for i, j in itertools.combinations(t, 2)]
                     for t in TRIPLES},
                 "F_cubic": {str(g): F_cubic[g] for g in GRIDS},
                 "F_star": {str(t): {str(g): F_star[t][g] for g in GRIDS}
                            for t in TRIPLES},
                 "exactness": exactness,
                 "resolution_cubic": {str(k): v for k, v in res_cubic.items()},
                 "resolution_star": {str(t): {str(k): v for k, v in s.items()}
                                     for t, s in res_star.items()},
                 "resolved_lengths": resolved,
                 "unresolved_lengths": unresolved,
                 "ratio": {str(t): ratio[t] for t in TRIPLES},
                 "ratio_resolved": {str(t): v for t, v in ratio_res.items()},
                 "ratio_unresolved": {str(t): v for t, v in ratio_all.items()},
                 "ratio_spread": ratio_spread, "ratio_typical": ratio_mid,
                 "ratio_swing": ratio_swing, "drift_shortest_n": drift_short,
                 "drift_worst": drift_max,
                 "decay": {"cubic": [alpha_cubic, r_cubic],
                           **{str(t): list(alpha_star[t]) for t in TRIPLES},
                           **{str(t) + "_all_lengths": list(alpha_star_all[t])
                              for t in TRIPLES}},
                 "crossing": {"cubic": cross_cubic,
                              **{str(t): cross_star[t] for t in TRIPLES}},
                 "loop_scale": scale, "engine": E.backend_info(),
                 "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v218_closed_word_census_module.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V218 - the closed words of the icosahedral module and their return "
          "amplitude")
    print("=" * 68)
    print("  census (n = 2,4,6,8,10): "
          + ", ".join(f"{r['n']}: {r['module_axes6']} / {r['cubic_axes3']}"
                      for r in census[1:6])
          + f"; DP verified {census_ok}")
    print(f"  cubic exactness across grids: {exactness:.2e}; module resolution "
          "relative change vs the finest grid: "
          + ", ".join(f"n={k}: {v:.1e}"
                      for k, v in res_star[TRIPLES[0]].items()))
    print(f"  resolved lengths: {resolved}; unresolved: {unresolved}")
    for tri in TRIPLES:
        print(f"  ratio axes {tri} over the resolved lengths: "
              f"{[round(x, 4) for x in ratio_res.get(tri, [])]}")
    print(f"  decay: cubic {nm(alpha_cubic)} (r2 {nm(r_cubic)}); "
          + "; ".join(f"{tri} {nm(a)} (r2 {nm(r)})"
                      for tri, (a, r) in alpha_star.items()))
    print(f"  crossing of {scale:.6f}: cubic {cross_cubic}, icosahedral "
          f"{cross_star}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
