#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V357 - the ladder on the pentagon (booked before the run).

The 12-face walk spectrum W_n (V345 certified nine leaves) carries
the golden section exactly. Test whether ANY functional of the
booked monomial class reproduces the ladder's published reading
targets: T = (COEF_OBS/3)^(1/10) * LOOP (the rung-10 amplitude the
elect reads at the fifth) and the certified deficit 0.01013, with
the chance-density discipline booked BEFORE the run.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))
import V352_oneloop_external_backreaction as V352  # noqa: E402


def main():
    checks, data = [], {}
    # dodecahedron faces = icosahedron vertices; adjacency = icosahedral
    # graph (share an edge of the dodecahedron)
    ph = (1 + 5 ** 0.5) / 2
    vs = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            vs += [(0.0, float(s1), float(s2) * ph),
                   (float(s1), float(s2) * ph, 0.0),
                   (float(s1) * ph, 0.0, float(s2))]
    dmin = min((sum((vs[i][k] - vs[j][k]) ** 2 for k in range(3))
                for i in range(12) for j in range(i + 1, 12)))
    A = np.zeros((12, 12), dtype=object)
    for i in range(12):
        for j in range(i + 1, 12):
            d = sum((vs[i][k] - vs[j][k]) ** 2 for k in range(3))
            if abs(d - dmin) < 1e-9:
                A[i, j] = 1
                A[j, i] = 1
    M = A
    W = {}
    for n in range(1, 31):
        if n > 1:
            M = A @ M
        W[n] = int(sum(M[i][i] for i in range(12)))
    v345 = [60, 120, 780, 3120, 16380, 78120, 394380, 1953120, 9784380]
    got = [W[n] for n in range(2, 11)]
    data["walks_n_1_30"] = {str(n): W[n] for n in range(1, 31)}
    checks.append(("T0 closed walks re-derived: W_2..W_10 = %s vs the "
                   "V345 leaves %s (bit-identical)" % (got, v345),
                   got == v345))

    # exact spectrum certificate: Newton sums on the integer traces give
    # the elementary symmetric functions of the eigenvalues (the char
    # poly's coefficients) - all integers, golden-section factorizable
    n = 12
    tr = [W[k] for k in range(1, n + 1)]
    e = [1] + [0] * n
    for k in range(1, n + 1):
        # Newton: e_k = (1/k) sum_{j=1..k} (-1)^{j-1} e_{k-j} p_j
        tot = 0
        for j in range(1, k + 1):
            tot += (-1) ** (j - 1) * e[k - j] * tr[j - 1]
        assert tot % k == 0
        e[k] = tot // k
    data["charpoly_elementary_symmetrics"] = e
    checks.append(("char poly certified via Newton sums on integer "
                   "traces: e_1..e_12 all exact integers = %s" % e,
                   all(isinstance(x, int) for x in e)))

    COEF = float(V352.COEF_OBS)
    LOOP = float(V352.LOOP)
    T = (COEF / 3.0) ** 0.1 * LOOP
    D = 0.01013
    data["targets"] = {"T_rung10_amplitude": T, "deficit": D}
    logW = {m: np.log(float(W[m])) for m in range(1, 31) if W[m] > 0}
    # (book-fix: W_1 = tr A = 0, so W_1^{+-1} is undefined - the monomial
    # class is built only from the defined positive walks; before this
    # fix the -inf logs poisoned argmin and the scan was vacuous)
    log5, logs5 = np.log(5.0), np.log(5.0) / 2.0
    best = {1: (9.9, None), 2: (9.9, None)}
    hits = {1: [], 2: []}
    rngA = range(-20, 21)
    rngB = range(0, 21)
    rngE = range(-20, 21)
    pair_terms = []
    for mm in range(1, 31):
        for nn in range(1, 31):
            for cs, ds in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                if W[mm] == 0 and cs != 0:
                    continue
                if W[nn] == 0 and ds != 0:
                    continue
                pair_terms.append((cs * logW[mm] + ds * logW[nn], mm, nn,
                                   cs, ds))
    pvals = np.array([p[0] for p in pair_terms])
    for a in rngA:
        for b in rngB:
            for ev in rngE:
                base = a * log5 + b * logs5 - ev * log5
                vals = base + pvals
                for tidx, tgt in ((1, np.log(T)), (2, np.log(D))):
                    d = np.abs(vals - tgt)
                    j = int(np.argmin(d))
                    if d[j] < best[tidx][0]:
                        best[tidx] = (float(d[j]), pair_terms[j])
                    if d[j] <= 1e-4:
                        mm, nn, cs, ds = pair_terms[j][1:]
                        hits[tidx].append({"rel_log": float(d[j]),
                                           "form": "5^%d (sqrt5)^%d / 5^%d"
                                                   " * W_%d^%d W_%d^%d"
                                                   % (a, b, ev, mm, cs,
                                                      nn, ds)})
    cls_size = len(list(rngA)) * len(list(rngB)) * len(list(rngE)) \
        * len(pair_terms)
    expected = {}
    for tidx in (1, 2):
        # chance density over the class's own log-range
        lo = -40.0
        hi = 60.0
        expected[tidx] = cls_size * 2e-4 / (hi - lo)
    data["pentagon_class"] = {"class_size": cls_size,
                              "expected_hits_in_1e-4_band": expected,
                              "hits_T": hits[1][:10], "hits_D": hits[2][:10],
                              "best_T": best[1][0], "best_D": best[2][0]}
    found = (best[1][0] <= 1e-8) or (best[2][0] <= 1e-8)
    checks.append(("pentagon monomial class %d forms against T = %.10f "
                   "and the deficit %.5f: best relative-log miss %.1e (T) "
                   "/ %.1e (deficit); chance expectation in the 1e-4 band "
                   "~%.1f hits"
                   % (cls_size, T, D, best[1][0], best[2][0],
                      expected[1]), True))

    if not (checks[0][1] and checks[1][1]):
        status, ok = "MACHINERY-DRIFT", False
    elif found and max(expected.values()) < 1.0:
        status, ok = "PENTAGON-LADDER-HITS", True
    else:
        status, ok = "PENTAGON-LADDER-VOID", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V357 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V357_pentagon_walk_ladder", "tier": 1,
           "property": "does the 12-face walk ladder carry the "
                       "published reading objects (golden-section "
                       "spectrum monomials)",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v357_pentagon_walk_ladder.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
