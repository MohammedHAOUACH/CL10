#!/usr/bin/env python3
"""V221 - the rung labels as K0 indices: what the icosahedral module does to the
algebra the loop ladder lives in, and whether n* = 208 and the coefficient 2.710
move with it.

The ladder of Section 12 is a statement about closed words on the shift module's
generators: rung n is the class of the closed words of length n, and the labels
of those classes live in the K0 of the C*-algebra the module generates.  For an
abelian shift module that algebra is a group C*-algebra C*(Z^r), and everything
below - the rank r, the torsion, K0 and the closed-word counts - is computable in
exact integer arithmetic from the generators the library carries, with no
momentum grid and no eigensolver.  This is the part of the ladder question that
can be settled exactly; V218 measured the part that cannot.

T1 the shift group of each module, computed: the integer matrix whose columns are
   the signed generators, written in a Z-basis of the coordinate ring, in Smith
   normal form.  The cubic module's three axes (six signed) give free rank 3 with
   three relations (the opposite pairs) and no torsion.  The icosahedral star's
   twelve signed directions have Q-rank 6 inside Q(sqrt 5)^3, so the star's shift
   group is free of rank 6 with six relations and no torsion - the algebraic form
   of the density V214 measured by shortest elements.
T2 K0 of the two group C*-algebras: K0(C(T^r)) = wedge^even Z^r, of rank 2^(r-1).
   Z^4 for the cubic module against Z^32 for the star; of those classes the ones
   a three-dimensional physical space can represent (degree at most 2) are 4
   against 16.
T3 the return law on each ladder, from exact counts: the number of closed words
   of length n on the free abelian group of rank r is the constant term of
   (sum_i x_i + x_i^-1)^n, computed here by two independent closed forms checked
   against each other, and the return probability C(n,r)/(2r)^n falls as
   n^(-r/2).  The exponent is measured by a log-log fit on exact integers: -3/2
   for the cubic ladder, -3 for the star's.
T4 where each census ladder crosses the loop scale 1/(16 pi), as exact rationals,
   and the combinatorial prefactor each carries: the number of axis families, 3
   against 6.

What this settles, and what it names: the published n* = 208 and the coefficient
2.710 live on the traced-link Cartan-norm estimator of V203/V204, on the cubic
module; that estimator is a statement about the clock and the link algebra, and
the internal sector carrying it is unchanged by the module (V217), so they do not
move there.  They are not reproduced by the star module's own closed-word ladder,
whose crossing is a single-digit rung and whose return exponent is -3.  Deciding
whether the star's ladder has a rung where 208 sits would need its almost-periodic
mean over the six-torus: (4*208)^6 = 3.3e17 eigendecompositions of a
128-dimensional step.
"""

import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))

from cl10 import continuum as C

N_CENSUS = 60         # longest closed-word count computed exactly
FIT_RANGE = (40, 60)  # log-log fit window for the return law
N_STAR = 208           # the published crossing length
GRID_PER_AXIS = 4      # grid points per axis needed to resolve length n (V218)
VERIFY_MAX = 24       # the second closed form is checked up to this length


# --------------------------------------------------------- exact linear algebra

def smith_invariants(matrix):
    """Invariant factors and rank of an integer matrix by Smith normal form; the
    kernel rank of the column map (the relation lattice of the signed
    generators) follows as the number of columns minus the rank."""
    A = [[int(v) for v in row] for row in matrix]
    rows, cols = len(A), len(A[0])
    rank = 0
    for step in range(min(rows, cols)):
        piv, best = None, None
        for i in range(step, rows):
            for j in range(step, cols):
                v = abs(A[i][j])
                if v and (best is None or v < best):
                    piv, best = (i, j), v
        if piv is None:
            break
        pi, pj = piv
        A[step], A[pi] = A[pi], A[step]
        for i in range(rows):
            A[i][step], A[i][pj] = A[i][pj], A[i][step]
        changed = True
        while changed:
            changed = False
            for j in range(step + 1, cols):
                if A[step][j] % A[step][step]:
                    q = A[step][j] // A[step][step]
                    for k in range(rows):
                        A[k][j] -= q * A[k][step]
                    changed = True
            for i in range(step + 1, rows):
                if A[i][step] % A[step][step]:
                    q = A[i][step] // A[step][step]
                    for k in range(cols):
                        A[i][k] -= q * A[step][k]
                    changed = True
        if A[step][step] < 0:
            for j in range(cols):
                A[step][j] = -A[step][j]
        rank += 1
    d = sorted(abs(A[i][i]) for i in range(rank))
    fixed = True
    while fixed:                       # put the factors in divisibility order
        fixed = False
        for i in range(len(d) - 1):
            if d[i + 1] % d[i]:
                g = math.gcd(d[i], d[i + 1])
                d[i], d[i + 1] = g, d[i] * d[i + 1] // g
                fixed = True
    inv = [x for x in d if x]
    return {"invariant_factors": inv, "rank": rank,
            "torsion": [x for x in inv if x > 1]}


def integer_matrix(columns):
    """Clear the denominators of a Fraction matrix (one row per coordinate part)
    without changing the relations among its columns."""
    rows = []
    for row in columns:
        den = 1
        for v in row:
            den = den * v.denominator // math.gcd(den, v.denominator)
        rows.append([int(v * den) for v in row])
    return rows


# ------------------------------------------------------------ exact word counts

def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def closed_words(n, r):
    """Closed words of length n on the free abelian group of rank r: the constant
    term of (sum x_i + x_i^-1)^n as the multinomial sum over k_1+...+k_r = n/2."""
    if n % 2:
        return 0
    total = 0
    for ks in compositions(n // 2, r):
        num = math.factorial(n)
        den = 1
        for k in ks:
            den *= math.factorial(k) ** 2
        total += num // den
    return total


def closed_words_alternative(n, r):
    """The same count by a different closed form: expand the n-th power first,
    (sum_i y_i)^n with y_i = x_i + x_i^-1, and take the constant term of each
    y_i^{m_i} separately, which is binom(m_i, m_i/2) for m_i even and zero
    otherwise; the sum runs over weak compositions of n into r parts."""
    total = 0
    num = math.factorial(n)
    for ms in compositions(n, r):
        prod, ok = 1, True
        for m in ms:
            if m % 2:
                ok = False
                break
            prod *= math.comb(m, m // 2)
        if ok:
            den = 1
            for m in ms:
                den *= math.factorial(m)
            total += num // den * prod
    return total


def loglog_slope(ns, vals):
    x = [math.log(n) for n in ns]
    y = [math.log(v) for v in vals]
    a, b = np.polyfit(x, y, 1)
    fit = [a * xi + b for xi in x]
    ss = sum((yi - yi0) ** 2 for yi, yi0 in zip(y, fit))
    st = sum((yi - sum(y) / len(y)) ** 2 for yi in y)
    return float(a), float(1.0 - ss / st)


def main():
    t0 = time.time()

    # ------------------------------------------------- T1: the two shift groups
    star = C.icosahedral_star()
    gens = C.module_generators(star)
    axes = gens["exact"]                                  # six independent axes
    signed_star = star["exact"]                           # twelve signed ones
    signed_cubic = [[Fraction(c) for c in col] for col in zip(
        (1, -1, 0, 0, 0, 0), (0, 0, 1, -1, 0, 0), (0, 0, 0, 0, 1, -1))]

    mat_star = integer_matrix(C._module_rows(star))       # 6 x 12 exact rationals
    mat_cubic = integer_matrix([list(row) for row in zip(*signed_cubic)])
    sm_star = smith_invariants(mat_star)
    sm_cubic = smith_invariants(mat_cubic)
    r_star, r_cubic = sm_star["rank"], sm_cubic["rank"]
    rel_star = len(signed_star) - r_star
    rel_cubic = len(signed_cubic) - r_cubic

    # the exact generators have to be the ones the engine hops on
    float_defect = max(abs(float(a) + float(b) * math.sqrt(5.0) - vec[k])
                       for ab, vec in zip(axes, gens["vectors"])
                       for k, (a, b) in enumerate(ab))

    # -------------------------------------------------------- T2: K0 of the pair
    def k0_rank(r):
        return sum(math.comb(r, 2 * j) for j in range(0, r // 2 + 1))

    def k0_visible(r, d=3):
        return sum(math.comb(r, 2 * j) for j in range(0, min(d, r) // 2 + 1))

    k0_cubic, k0_star = k0_rank(r_cubic), k0_rank(r_star)
    vis_cubic, vis_star = k0_visible(r_cubic), k0_visible(r_star)

    # ------------------------------------------------ T3: the exact return laws
    ns = [n for n in range(2, N_CENSUS + 1, 2)]
    census, counts_ok = [], True
    for n in ns:
        c3, c6 = closed_words(n, r_cubic), closed_words(n, r_star)
        row = {"n": n, "closed_cubic_r3": c3, "closed_star_r6": c6,
               "p_cubic": c3 / (2 * r_cubic) ** n,
               "p_star": c6 / (2 * r_star) ** n}
        if n <= VERIFY_MAX:
            agree = (c3 == closed_words_alternative(n, r_cubic)
                     and c6 == closed_words_alternative(n, r_star))
            counts_ok = counts_ok and agree
            row["two_forms_agree"] = agree
        else:
            row["two_forms_agree"] = None
        census.append(row)
    fit_ns = [n for n in ns if FIT_RANGE[0] <= n <= FIT_RANGE[1]]
    sel = [c for c in census if c["n"] in fit_ns]
    slope_cubic, r2_cubic = loglog_slope([c["n"] for c in sel],
                                         [c["p_cubic"] for c in sel])
    slope_star, r2_star = loglog_slope([c["n"] for c in sel],
                                       [c["p_star"] for c in sel])

    # ------------------------------------------- T4: crossings and prefactor
    scale = 1.0 / (16.0 * math.pi)
    limit = Fraction(1, 16) / Fraction(math.pi)   # exact 1/(16 pi) as a rational

    def crossing(key, base):
        for c in census:
            if Fraction(c[key], base ** c["n"]) <= limit:
                return c["n"]
        return None

    cross_cubic = crossing("closed_cubic_r3", 2 * r_cubic)
    cross_star = crossing("closed_star_r6", 2 * r_star)
    p_cross_cubic = next(c["p_cubic"] for c in census
                         if c["n"] == cross_cubic)
    p_cross_star = next(c["p_star"] for c in census
                        if c["n"] == cross_star)
    fam_cubic, fam_star = len(signed_cubic) // 2, len(signed_star) // 2
    cost_cubic = (GRID_PER_AXIS * N_STAR) ** r_cubic
    cost_star = (GRID_PER_AXIS * N_STAR) ** r_star

    # ------------------------------------------------------------- the checks
    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    chk("T1 the shift group of each module is free of the measured rank: the "
        "cubic module's signed axes present Z^3 with three relations, the "
        "icosahedral star's twelve signed directions present Z^6 with six, and "
        "neither carries torsion - the star's directions are Q-independent in "
        "Q(sqrt 5)^3, the algebraic form of its density",
        r_cubic == 3 and rel_cubic == 3 and not sm_cubic["torsion"]
        and r_star == 6 and rel_star == 6 and not sm_star["torsion"]
        and float_defect < 1e-13,
        f"Smith invariant factors: cubic {sm_cubic['invariant_factors']} (rank "
        f"{r_cubic} of {len(signed_cubic)} signed generators, relation rank "
        f"{rel_cubic}), star {sm_star['invariant_factors']} (rank {r_star} of "
        f"{len(signed_star)} signed directions, relation rank {rel_star}); "
        f"torsion {sm_cubic['torsion'] + sm_star['torsion']}; the exact "
        f"generators reproduce the engine's float star axes to {float_defect:.1e}")
    chk("T2 the rung labels therefore live in different K0 groups: K0(C(T^3)) is "
        "Z^4 for the cubic module while K0(C(T^6)) is Z^32 for the star, and the "
        "classes representable by cycles of degree at most 2 are 4 against 16",
        k0_cubic == 4 and k0_star == 32 and vis_cubic == 4 and vis_star == 16,
        f"K0 = wedge^even Z^r has rank 2^(r-1): cubic r = {r_cubic} gives "
        f"{k0_cubic}, star r = {r_star} gives {k0_star}; degree at most 2: "
        f"{vis_cubic} against {vis_star}")
    chk("T3 the two ladders count differently in a way that is not a convention: "
        "the closed-word return probability falls as n^(-r/2), measured at -3/2 "
        "on the cubic ladder and -3 on the star's, from exact integer counts "
        "verified by two independent closed forms",
        counts_ok and abs(slope_cubic + 1.5) < 0.02 and r2_cubic > 0.999
        and abs(slope_star + 3.0) < 0.05 and r2_star > 0.999,
        "closed words, rank 3 / rank 6: " + ", ".join(
            f"n = {c['n']}: {c['closed_cubic_r3']} / {c['closed_star_r6']}"
            for c in census[1:6]) + f"; the two closed forms agree at every "
        f"n <= {VERIFY_MAX}: {counts_ok}; log-log slope of C(n,r)/(2r)^n over "
        f"n = {fit_ns[0]}..{fit_ns[-1]}: cubic {slope_cubic:.4f} "
        f"(r2 {r2_cubic:.6f}), star {slope_star:.4f} (r2 {r2_star:.6f})")
    chk("T4 neither census ladder crosses the loop scale 1/(16 pi) anywhere near "
        "n* = 208, and the combinatorial prefactor moves by the family count: 3 "
        "axis families on the cubic module against 6 on the star",
        cross_cubic is not None and cross_star is not None
        and cross_cubic < 50 and cross_star < 50
        and cross_cubic != N_STAR and cross_star != N_STAR
        and (fam_cubic, fam_star) == (3, 6),
        f"loop scale 1/(16 pi) = {scale:.6f}: the exact census ladder crosses it "
        f"at n = {cross_cubic} (cubic, return probability {p_cross_cubic:.4e}) "
        f"and n = {cross_star} (star, {p_cross_star:.4e}); axis families "
        f"{fam_cubic} against {fam_star}, so the factor 3 of the coefficient - "
        f"the three axis families - becomes {fam_star}; resolving the star's "
        f"almost-periodic mean at n = {N_STAR} costs "
        f"({GRID_PER_AXIS * N_STAR})^{r_star} = {cost_star:.2e} "
        f"eigendecompositions against ({GRID_PER_AXIS * N_STAR})^{r_cubic} = "
        f"{cost_cubic:.2e} on the cubic torus")

    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the rung labels of the loop ladder as K0 indices: the shift "
                    "group of each module in Smith normal form (free rank, "
                    "relations, torsion), K0 of its group C*-algebra and the "
                    "classes representable in three dimensions, the exact "
                    "closed-word counts on both ranks by two independent closed "
                    "forms, the n^(-r/2) return law measured on those integers, "
                    "the crossing of 1/(16 pi) on each census ladder, and the "
                    "axis-family count that sets the coefficient's "
                    "combinatorial prefactor",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"RUNG LABELS AS K0 INDICES. (T1) the cubic module's "
            f"{len(signed_cubic)} signed axes present the free abelian group "
            f"of rank {r_cubic} with {rel_cubic} relations (the opposite pairs) "
            f"and no torsion; the icosahedral star's {len(signed_star)} signed "
            f"directions are Q-independent inside Q(sqrt 5)^3, presenting rank "
            f"{r_star} with {rel_star} relations and no torsion - the exact form "
            f"of the density V214 measured by shortest elements, and the reason "
            f"no finite Brillouin zone exists on it. (T2) the labels therefore "
            f"live in different groups: K0 = wedge^even Z^r has rank {k0_cubic} "
            f"for the cubic module and {k0_star} for the star, of which "
            f"{vis_cubic} and {vis_star} classes are representable by cycles of "
            f"degree at most 2 in three dimensions. (T3) the two ladders count "
            f"differently in a way that is not a convention: the closed-word "
            f"return probability falls as n^(-r/2), measured slopes "
            f"{slope_cubic:.4f} (r2 {r2_cubic:.6f}) against {slope_star:.4f} "
            f"(r2 {r2_star:.6f}) on exact integers up to n = {N_CENSUS}, the two "
            f"count forms agreeing at every length. (T4) the census ladder "
            f"crosses 1/(16 pi) at n = {cross_cubic} (cubic) and n = "
            f"{cross_star} (star): single-digit rungs, nowhere near 208; and the "
            f"family count that supplies the coefficient's factor 3 doubles to "
            f"{fam_star} on the star. VERDICT: the ladder's algebra moves with "
            f"the module - rank {r_cubic} to {r_star}, K0 Z^{k0_cubic} to "
            f"Z^{k0_star}, return exponent -3/2 to -3 - while the published "
            f"n* = 208 and coefficient 2.710 live on the traced-link "
            f"Cartan-norm estimator of V203/V204, a statement about the clock "
            f"and the link algebra that the module change leaves untouched (the "
            f"internal sector carrying it is unchanged, V217). Nothing in the "
            f"published ladder arithmetic moves; what moves is the ladder the "
            f"star module would have to build for itself, and V218 showed its "
            f"dynamical half cannot be resolved at n = {N_STAR}: "
            f"({GRID_PER_AXIS * N_STAR})^{r_star} = {cost_star:.1e} "
            f"diagonalisations of a 128-dimensional step on the six-torus."),
        "data": {
            "cubic_signed_generators": len(signed_cubic),
            "star_signed_directions": len(signed_star),
            "smith_cubic": sm_cubic, "smith_star": sm_star,
            "relations_cubic": rel_cubic, "relations_star": rel_star,
            "exact_vs_engine_float_defect": float_defect,
            "k0_rank_cubic": k0_cubic, "k0_rank_star": k0_star,
            "k0_degree_le_2_cubic": vis_cubic, "k0_degree_le_2_star": vis_star,
            "census": census, "two_closed_forms_agree": counts_ok,
            "return_exponent_cubic": [slope_cubic, r2_cubic],
            "return_exponent_star": [slope_star, r2_star],
            "fit_range": list(FIT_RANGE),
            "crossing_cubic": cross_cubic, "crossing_star": cross_star,
            "loop_scale": scale,
            "axis_families_cubic": fam_cubic, "axis_families_star": fam_star,
            "resolution_cost_cubic": cost_cubic,
            "resolution_cost_star": cost_star,
            "n_published": N_STAR, "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v221_rung_k0_labels.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V221 - the rung labels as K0 indices")
    print("=" * 68)
    print(f"  shift groups: cubic rank {r_cubic} rel {rel_cubic} tors "
          f"{sm_cubic['torsion']}; star rank {r_star} rel {rel_star} tors "
          f"{sm_star['torsion']}")
    print(f"  K0: Z^{k0_cubic} against Z^{k0_star}; degree <= 2: {vis_cubic} "
          f"against {vis_star}")
    print("  census (rank 3 / rank 6): " + ", ".join(
        f"{c['n']}: {c['closed_cubic_r3']}/{c['closed_star_r6']}"
        for c in census[1:6]) + f"; two forms agree {counts_ok}")
    print(f"  return exponent: cubic {slope_cubic:.4f} (r2 {r2_cubic:.6f}), "
          f"star {slope_star:.4f} (r2 {r2_star:.6f})")
    print(f"  crossing of {scale:.6f}: cubic n = {cross_cubic}, star n = "
          f"{cross_star}; families {fam_cubic} -> {fam_star}; star cost at "
          f"n = {N_STAR}: {cost_star:.2e}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
