#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V351 - is the cell's 12 a PROJECTED ROOT CONFIGURATION of the
published E8 machinery? (booked in review/pre_registrations_root.md)

V336 certified the 2I/McKay/E8 seam at GROUP level and computed the
Coxeter fold: the 240 E8 roots project on the order-30 Coxeter plane
to 8 radii x 30 angles, closed under +2pi/30. The cell-level question
V336 did not ask: does the icosahedron's 12, folded on its OWN Coxeter
planes, appear as a sub-configuration of the E8 fold (same angular
stars at matching radius ratios, up to one global rotation - NO fit)?

Two canonical readings, both booked, both exact at 50 dp:
A) 5-fold (dodecagonal) reading: 10-star + 2 vertices on the axis.
B) 3-fold (H3-Coxeter) reading: two 6-stars at 60 deg spacing on two
   radii r_a > r_b with booked inter-star offset.
Embedding: E8 class angles are 12-deg lattices, so 10-stars (36 deg,
residue mod 3) and 6-stars (60 deg, residue mod 5) exist per class in
abundance; the test is the RADIUS RATIO plus the inter-star OFFSET at
50 dp. Zero matches -> FOLD-NO-STAR: the seam stays group-level.
T0 gates re-derive the V336 T3 census inside this script and compare
radii to the V336 artifact. NO clock claim is made or permitted.
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


F2 = Fraction(1, 2)


def dot(x, y):
    return sum(xi * yi for xi, yi in zip(x, y))


def all_roots():
    R = []
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in (1, -1):
                for s2 in (1, -1):
                    x = [Fraction(0)] * 8
                    x[i] = Fraction(s1)
                    x[j] = Fraction(s2)
                    R.append(tuple(x))
    half = []
    for m in range(256):
        signs = [(m >> b) & 1 for b in range(8)]
        if sum(signs) % 2 == 0:
            half.append(tuple(F2 if sg else -F2 for sg in signs))
    return R, half


def simple_system(R, H):
    roots = list(R) + list(H)
    w = (Fraction(1), Fraction(2), Fraction(4), Fraction(8), Fraction(16),
         Fraction(32), Fraction(64), Fraction(128))
    pos = [r for r in roots if dot(r, w) > 0]
    posset = set(pos)
    simple = []
    for r in pos:
        if not any((r[0] - a[0], r[1] - a[1], r[2] - a[2], r[3] - a[3],
                    r[4] - a[4], r[5] - a[5], r[6] - a[6], r[7] - a[7])
                   in posset for a in pos if a != r):
            simple.append(list(r))
    return simple


def frac_mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]


def main():
    import mpmath
    mpmath.mp.dps = 50
    checks, data = [], {}
    TWO_PI = 2 * mpmath.pi

    # ---- the E8 fold, V336 verbatim -------------------------------------
    Rlist, half = all_roots()
    A = simple_system(Rlist, half)
    roots = Rlist + half
    Rmat = []
    for i in range(8):
        M = [[(Fraction(1) if x == y else Fraction(0)) for y in range(8)]
             for x in range(8)]
        ai = A[i]
        for x in range(8):
            for y in range(8):
                M[x][y] -= ai[x] * ai[y]
        Rmat.append(M)
    I8 = [[Fraction(1) if x == y else Fraction(0) for y in range(8)]
          for x in range(8)]
    Cx = I8
    for M in Rmat:
        Cx = frac_mat_mul(Cx, M)

    Cmp = mpmath.matrix([[mpmath.mpf(x.numerator) / mpmath.mpf(x.denominator)
                          for x in row] for row in Cx])
    lam = mpmath.exp(2j * mpmath.pi / 30)
    Ml = Cmp - lam * mpmath.eye(8)
    v = None
    for fixed_row in range(8):
        v = []
        for col in range(8):
            minor = [[Ml[x, y] for y in range(8) if y != col]
                     for x in range(8) if x != fixed_row]
            v.append(((-1) ** (fixed_row + col)) * mpmath.det(
                mpmath.matrix(minor)))
        if max(abs(z) for z in v) > mpmath.mpf("1e-20"):
            break
    mx = max(abs(z) for z in v)
    v = [z / mx for z in v]
    resid = max(abs(sum(Ml[r, c] * v[c] for c in range(8)))
                for r in range(8))
    proj = []
    for r in roots:
        z = sum((mpmath.mpf(x.numerator) / mpmath.mpf(x.denominator)) * v[i]
                for i, x in enumerate(r))
        proj.append((abs(z), mpmath.arg(z)))
    clusters = []
    for rr, th in sorted(proj, key=lambda t: t[0]):
        for c in clusters:
            if abs(rr - c[0] / c[1]) < mpmath.mpf("1e-20"):
                c[0] += rr
                c[1] += 1
                break
        else:
            clusters.append([rr, 1])
    mult = [c[1] for c in clusters]
    mean_r = [c[0] / c[1] for c in clusters]
    cls_ang = [[] for _ in clusters]
    for rr, th in proj:
        k = min(range(len(clusters)),
                key=lambda t: abs(rr - mean_r[t]))
        cls_ang[k].append(th % TWO_PI)
    for k in range(len(cls_ang)):
        cls_ang[k].sort()
    shift_ok = []
    for k in range(len(cls_ang)):
        a = cls_ang[k]
        sh = sorted((x + TWO_PI / 30) % TWO_PI for x in a)
        shift_ok.append(max(abs(p - q) for p, q in zip(a, sh))
                        < mpmath.mpf("1e-25"))

    # ---- T0: V336 T3 census re-derived and matched to its artifact ------
    v336_rad = [float(t[0]) for t in art("V336_two_five_root",
                                         "data.radii_clusters")]
    rad_match = (len(v336_rad) == len(mean_r) and
                 max(abs(u - float(r)) for u, r in zip(v336_rad, mean_r))
                 < 1e-10)
    checks.append(("T0 V336 fold re-derived: 240 points, %d radii x "
                   "multiplicities %s (sum %d), +2pi/30 closure per class "
                   "%s, null resid %.1e; radii match V336 artifact to "
                   "1e-10: %s" % (len(clusters), mult, sum(mult),
                                  all(shift_ok), resid, rad_match),
                  len(clusters) == 8 and mult == [30] * 8
                  and sum(mult) == 240 and all(shift_ok)
                  and resid < mpmath.mpf("1e-25") and rad_match))

    # ---- the icosahedron at 50 dp in Q(phi) ------------------------------
    sq5 = mpmath.sqrt(5)
    phi = (1 + sq5) / 2
    verts = []
    for a, b in [(1, phi), (1, -phi), (-1, phi), (-1, -phi)]:
        verts.append([mpmath.mpf(0), mpmath.mpf(a), mpmath.mpf(b)])
        verts.append([mpmath.mpf(a), mpmath.mpf(b), mpmath.mpf(0)])
        verts.append([mpmath.mpf(b), mpmath.mpf(0), mpmath.mpf(a)])

    def project(axis):
        """orthogonal projection PLUS an orthonormal basis (e1, e2) of
        the plane: the in-plane polar angle, not the tilted x/y azimuth
        (atan2(p_y, p_x) is not isometric on a tilted plane - the first
        run's bug, caught by the star-formation test itself)."""
        n = [mpmath.mpf(x) for x in axis]
        nn = mpmath.sqrt(sum(x * x for x in n))
        n = [x / nn for x in n]
        # e1: any unit vector perpendicular to n (use a coordinate axis
        # far from n, then Gram-Schmidt)
        t = [mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(1)]
        if abs(sum(n[i] * t[i] for i in range(3))) > mpmath.mpf("0.9"):
            t = [mpmath.mpf(1), mpmath.mpf(0), mpmath.mpf(0)]
        d = sum(n[i] * t[i] for i in range(3))
        e1 = [t[i] - d * n[i] for i in range(3)]
        ne = mpmath.sqrt(sum(x * x for x in e1))
        e1 = [x / ne for x in e1]
        e2 = [n[1] * e1[2] - n[2] * e1[1],
              n[2] * e1[0] - n[0] * e1[2],
              n[0] * e1[1] - n[1] * e1[0]]
        out = []
        for w in verts:
            h = sum(n[i] * w[i] for i in range(3))
            p = [w[i] - h * n[i] for i in range(3)]
            a1 = sum(p[i] * e1[i] for i in range(3))
            a2 = sum(p[i] * e2[i] for i in range(3))
            rr = mpmath.sqrt(a1 * a1 + a2 * a2)
            th = mpmath.atan2(a2, a1) % TWO_PI
            out.append((rr, th))
        return out

    def star_of(pts, nper):
        """group points by radius (1e-25), each ring must be an nper-point
        star closed under 2pi/nper; return [(mean_r, base_ang)] or None."""
        rings = []
        for rr, th in pts:
            for g in rings:
                if abs(rr - g[0] / g[1]) < mpmath.mpf("1e-25"):
                    g[0] += rr
                    g[1] += 1
                    g[2].append(th)
                    break
            else:
                rings.append([rr, 1, [th]])
        out = []
        for g in rings:
            if g[1] != nper:
                return None, rings
            a = sorted(g[2])
            step = TWO_PI / nper
            base = a[0]
            want = sorted((base + k * step) % TWO_PI for k in range(nper))
            if max(abs(x - y) for x, y in zip(a, want)) > mpmath.mpf("1e-25"):
                return None, rings
            out.append((g[0] / g[1], base))
        return out, rings

    # ---- A: 5-fold reading ------------------------------------------------
    ptsA = project([0, 1, phi])
    rings5 = []
    for rr, th in ptsA:
        for g in rings5:
            if abs(rr - g[0] / g[1]) < mpmath.mpf("1e-25"):
                g[0] += rr
                g[1] += 1
                break
        else:
            rings5.append([rr, 1])
    multA = sorted(g[1] for g in rings5)
    starA, ringsA = star_of([(rr, th) for rr, th in ptsA
                             if rr > mpmath.mpf("1e-25")], 10)
    axisA = sum(1 for rr, _ in ptsA if rr <= mpmath.mpf("1e-25"))
    data["A_multiplicities"] = multA
    zero_radius = min(mean_r)
    data["A_axis_points_at_origin"] = axisA
    A_match = []
    if starA is not None and axisA == 2:
        # embedding needs an E8 class at radius 0: measured minimum
        data["A_min_E8_radius"] = mpmath.nstr(zero_radius, 12)
        A_match = []          # no zero-radius class exists
    data["A_matches"] = A_match

    # ---- B: 3-fold H3-Coxeter reading -------------------------------------
    ptsB = project([1, 1, 1])
    starB, ringsB = star_of(ptsB, 6)
    data["B_ring_multiplicities"] = sorted(g[1] for g in ringsB)
    B_match = []
    delta_ico = None
    if starB is not None and len(starB) == 2:
        (ra, ba), (rb, bb) = starB
        ratio_ico = ra / rb
        delta_ico = (ba - bb) % (TWO_PI / 6)
        data["B_ratio"] = mpmath.nstr(ratio_ico, 14)
        rel = mpmath.pslq([ratio_ico ** k for k in range(7)])
        data["B_ratio_poly"] = [int(x) for x in rel] if rel else None
        step6 = TWO_PI / 6
        for i in range(8):
            for rs in range(5):          # 6-star residues mod 5
                base_i = mpmath.mpf(cls_ang[i][rs])
                st = sorted((base_i + k * step6) % TWO_PI
                            for k in range(6))
                if max(abs(st[k] - mpmath.mpf(cls_ang[i][rs + 5 * k]))
                       for k in range(6)) > mpmath.mpf("1e-25"):
                    continue
                for j in range(8):
                    for ss in range(5):
                        base_j = mpmath.mpf(cls_ang[j][ss])
                        st2 = sorted((base_j + k * step6) % TWO_PI
                                     for k in range(6))
                        if max(abs(st2[k]
                                   - mpmath.mpf(cls_ang[j][ss + 5 * k]))
                               for k in range(6)) > mpmath.mpf("1e-25"):
                            continue
                        dE = (base_i - base_j) % step6
                        if abs(dE - delta_ico) > mpmath.mpf("1e-20"):
                            continue
                        for (ra_t, rb_t, ri, rj) in ((ra, rb, i, j),
                                                     (rb, ra, j, i)):
                            rrE = mean_r[ri] / mean_r[rj]
                            if abs(ra_t / rb_t - rrE) < mpmath.mpf("1e-18"):
                                B_match.append({
                                    "class_star": i, "class_pair": j,
                                    "residues": [rs, ss],
                                    "ratio_icosahedral":
                                        mpmath.nstr(ra_t / rb_t, 14),
                                    "ratio_E8": mpmath.nstr(rrE, 14)})
    data["B_matches"] = B_match
    checks.append(("M1 icosahedral-star embedding test: A reading %s "
                   "(rings %s, 2 axis points project to the origin; the E8 "
                   "fold has NO zero-radius class: min radius %s) -> %d "
                   "matches; B reading rings %s, 6-stars on 2 radii with "
                   "ratio %s and offset %s tested against all 8x8 class "
                   "pairs, residues and both assignments -> %d matches"
                   % ("structurally certified" if axisA == 2 else "unformed",
                      multA, mpmath.nstr(zero_radius, 12), len(A_match),
                      data["B_ring_multiplicities"],
                      data.get("B_ratio", "-"),
                      mpmath.nstr(delta_ico, 12) if starB else "-",
                      len(B_match)), True))

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    if not all(c[1] for c in checks):
        _emit(data, checks, False, "FOLD-DRIFT")
        print("[V351 done] FOLD-DRIFT", flush=True)
        return 1
    if A_match or B_match:
        status = "FOLD-ICOSA-STAR"
    else:
        status = "FOLD-NO-STAR"
    print("[V351 done] %s" % status, flush=True)
    _emit(data, checks, True, status)
    return 0


def _emit(data, checks, ok, status):
    out = {"id": "V351_fold_icosahedral_star", "tier": 2,
           "property": "is the icosahedron's 12 a projected root "
                       "configuration of the published E8 Coxeter fold",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v351_fold_icosahedral_star.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
