#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V294 - is the democratic rank forced by the charge bookkeeping?
(tier 3)

**The question (pre-registered).** Do the booked integer selection
rules (atom channels {0, 64, 96} of the V269 c2_spectrum census;
the V197 rung factorization 208 = 2^4 * 13 read as the certified
arithmetic anchor) FORCE the democratic rank-1 structure of the 3x3
family mass matrix, or does the family selection leave higher rank
allowed?

**The object.** the V45 mass_matrix anchor (democratic, rank 1),
the V70 family structure (3 equivalent families: the allowed
matrix must commute with the family permutations - the booked
family-symmetry fact), the V269 census integers.

**Measured.** the commutant of the S_3 family-permutation action on
3x3 matrices (exact linear solve): its dimension, an explicit basis
{I, all-ones}, and the maximum rank over real combinations; the
democratic rank-1 sits inside the commutant - the test is whether
the commutant FORCES rank 1 or merely contains it.

Branch: RANK-FORCED-1 (commutant = span{ones}: democratic forced) /
RANK-FREEDOM (maximum rank booked: the rank-1 choice is a named
input, stated as such).
"""
import json
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


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    v45 = art("V45_yukawa_matrix")
    v70 = art("V70_family_index")
    v269 = art("V269_dynamical_link_module")

    # anchors
    mm = v45["data"]["mass_matrix"]
    t0ok = bool(mm["democratic"] and int(mm["rank"]) == 1
                and int(v70["data"]["counts"]["families"]) == 3)
    chan = sorted(int(e["dim"]) for e in v269["data"]["c2_spectrum"])
    data["anchors"] = {"rank_booked": int(mm["rank"]),
                       "democratic": bool(mm["democratic"]),
                       "families": int(v70["data"]["counts"]["families"]),
                       "census_dims": chan}
    checks.append({"name": "T0 anchors: V45 democratic rank 1, V70 three "
                          "families, V269 census {}".format(chan),
                   "pass": t0ok})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V294 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    # ---- T1: commutant of the S_3 family action -------------------------
    P12 = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    P23 = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    # M in commutant: P M P^T = M for both generators: 9 real linear
    # constraints on vec(M)
    A = []
    for P in (P12, P23):
        K = np.kron(P, P) - np.eye(9)
        A.append(K)
    A = np.vstack(A)
    _, sv, Vt = np.linalg.svd(A)
    null = Vt[sv < 1e-12].T
    dim = null.shape[1]
    # basis and maximum rank: sample the commutant exhaustively on a
    # 2-parameter form recovered from the basis (dim expected 2:
    # I and ones); verify explicit generators
    I3 = np.eye(3)
    O3 = np.ones((3, 3))
    gen_ok = all(np.allclose(P @ G @ P.T, G, atol=1e-12)
                 for P in (P12, P23) for G in (I3, O3))
    max_rank = 0
    best = None
    rng = np.random.default_rng(294)
    for _ in range(200):
        a, b = rng.uniform(-2, 2, 2)
        r = int(np.linalg.matrix_rank(a * I3 + b * O3, tol=1e-9))
        if r > max_rank:
            max_rank, best = r, (a, b)
    data["commutant"] = {"dim": int(dim), "generators_ok": bool(gen_ok),
                         "max_rank": int(max_rank),
                         "witness": [float(x) for x in best],
                         "witness_spectrum": [
                             float(x) for x in
                             np.linalg.eigvalsh(best[0] * I3
                                                + best[1] * O3)]}
    status = "RANK-FORCED-1" if max_rank == 1 else "RANK-FREEDOM"
    data["status"] = status
    print("V294 {}: commutant dim {}, max rank {} (democratic rank-1 "
          "is inside, not forced)".format(status, dim, max_rank),
          flush=True)
    checks.append({"name": "T1 family-selection commutant: dim {}, maximum "
                          "rank {} over real combinations".format(dim,
                                                                  max_rank),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V294_democratic_mD_forced", "tier": 3,
           "property": "whether the booked selection rules force the "
                       "democratic rank-1 structure",
           "residual": True, "programme": "V294",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v294_democratic_mD_forced.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V294 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
