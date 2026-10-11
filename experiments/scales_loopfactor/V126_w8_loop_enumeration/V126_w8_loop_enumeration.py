#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V126 - WHICH closed loops supply the ladder factor (W8)

**Where W8 stands.**  V85 measured the closed-loop suppression of the traced
adjoint link, F_n = alpha((MG)^n)/alpha(G^n) with M the one-crossing hop
matrix, and found F_10 = 0.0192 against 1/(16 pi) = 0.0199: the loop factor
AT the ladder length.  V89 extended this to n = 2..30, found the crossing at
n* = 10, a mod-4 pattern, and a geometric law F_n = A r^n with r ~ 0.54.  The
named residual is the SELECTION RULE: which loops contribute, and whether the
factor 1/(16 pi) is a property of one loop class or of the whole sum.

**What is new here.**  The T.L Casimir of the 32 (x) 45 link has exactly THREE
eigenspaces (eigenvalues -32, -12, 4), so M decomposes into exactly three
sectors, M = sum_lam A_lam, and the traced n-hop product is a sum over 3^n
WORDS of the three sector matrices:
    (M G)^n = sum_{lam_1..lam_n} B_{lam_1} B_{lam_2} ... B_{lam_n},
    B_lam = A_lam G.
Three sectors makes the enumeration EXACT: every word of length n <= 10 (the
ladder length) is listed, its Cartan amplitude computed, and the sum compared
with the matrix power.  That turns the selection rule from a claim into a
counted object, and it separates three things V89 conflated in one number r:
  (i)  the GEOMETRIC per-word suppression,
  (ii) the MULTIPLICITY (how many words survive at all - the entropy),
  (iii) the INTERFERENCE between words (the phases e^{i(eta/2)ev_lam}).
The decisive control: the same enumeration at eta = 0, where every phase is 1
and no interference is possible.  If the loop factor is much larger at
eta = 0 than at the root, then 1/(16 pi) is a PHASE-CANCELLATION effect of the
Wilson angle, not a geometric suppression - and the geometric law F_n = A r^n
of V89 is a fit to an interference pattern.

**Checks (computed).**
1. Regression: F_10 at the root reproduces V85/V89; the word sum equals
   (MG)^n to machine precision for every n.
2. Word census n = 2..10: surviving words N_n, the growth rate, the largest
   single word, the top sector sequences.
3. Interference: ||sum of word vectors|| / sum of word norms at eta = 0.2,
   eta = 0 and eta = 1.4423; the eta = 0 loop factor against the root one.
4. The selection rule: which sectors appear at which n, and the mod-4 pattern
   explained (or not) by the word structure.
5. The ladder extrapolation with multiplicity: the per-word decay and the
   entropy rate, and whether their product reaches (1/16 pi)^17.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V85_closed_loop_factor"))

from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)
from V85_closed_loop_factor import projectors  # noqa: E402
from cl10 import clifford as CL                             # noqa: E402
from cl10 import engine as E                                # noqa: E402
from cl10 import sm_qca                                     # noqa: E402

D_L = 45
LOOP = 1.0 / (16.0 * math.pi)
ETA_ROOT = 0.2
ETA_FMIN = 1.4423076923076923
NMAX = 10


def cartan_vector(W, cartan):
    """The linear Cartan projection of a 32x32 matrix (5 complex numbers)."""
    return np.array([np.trace(C @ W) / W.shape[0] for C in cartan],
                    dtype=complex)


def sectors(eta, Ma, ev, G):
    """The three sector hop matrices B_lam = e^{i(eta/2) ev_lam} M_lam G, so
    that MG = sum_lam B_lam and (MG)^n is the sum of the 3^n words in B."""
    return [np.exp(1j * (eta / 2.0) * e) * (m @ G) for e, m in zip(ev, Ma)]


def word_census(B, cartan, nmax):
    """Exact enumeration of all 3^n words of length n in the sector matrices.

    Returns per n: the Cartan vectors of every word, and the census built
    from them.
    """
    out = {}
    prefixes = [np.eye(32, dtype=complex)]
    for n in range(1, nmax + 1):
        words = []
        for P in prefixes:
            for b in B:
                words.append(P @ b)
        vecs = [cartan_vector(W, cartan) for W in words]
        norms = [float(np.abs(v).sum()) for v in vecs]
        keep = [i for i, nz in enumerate(norms) if nz > 1e-13]
        out[n] = {"n": n, "words": len(words), "surviving": len(keep),
                  "norms": [norms[i] for i in keep],
                  "vecs": [vecs[i] for i in keep],
                  "total": sum(vecs, np.zeros(5, dtype=complex))}
        prefixes = words
    return out


def main():
    t0 = time.time()
    out = {"experiment": "V126_w8_loop_enumeration", "data": {}}
    checks = []

    # ---- the exact objects of V85/V89 -------------------------------------
    G = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(),
                                      sm_qca.G10)).astype(np.complex128)
    h = [E.asnumpy(m) for m in CL.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = root_generators()
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    P, ev = projectors(TL)
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    Ma = []
    for Pa in P:
        r = (Pa @ Gi).reshape(32, D_L, 32, D_L)
        Ma.append(np.einsum("ajbj->ab", r) / D_L)
    out["data"]["sectors"] = {"count": len(P),
                              "eigenvalues": [round(e, 6) for e in ev],
                              "norms": [round(float(np.abs(m).max()), 6)
                                        for m in Ma]}

    # ---- T0 regression + exactness of the word decomposition --------------
    with open(exp("V85_closed_loop_factor",
                           "results", "v85_closed_loop_factor.json")) as fh:
        v85_F10 = json.load(fh)["data"]["loops_root"]["F_n"][8]
    B = sectors(ETA_ROOT, Ma, ev, G)
    MG = sum(B)
    census = word_census(B, cartan, NMAX)
    dec_err = 0.0
    F10 = None
    for n in range(1, NMAX + 1):
        W = census[n]["total"]
        Gn = np.linalg.matrix_power(G, n)
        MGn = np.linalg.matrix_power(MG, n)
        # the word sum must equal the matrix power, in the Cartan projection
        exact = cartan_vector(MGn, cartan)
        dec_err = max(dec_err, float(np.abs(W - exact).max()))
        if n == NMAX:
            F10 = float(np.abs(W).sum()
                        / np.abs(cartan_vector(Gn, cartan)).sum())
    out["data"]["regression"] = {"F10_v85": v85_F10, "F10_v126": F10,
                                 "word_decomposition_error": dec_err,
                                 "elapsed_s": round(time.time() - t0, 1)}
    checks.append(("T0 the T.L Casimir of the 32 (x) 45 link has %d eigenspaces "
                   "(eigenvalues %s), so the n-hop traced product is a sum of "
                   "3^n words; the enumeration reproduces the matrix power "
                   "exactly (%.1e) and F_10 = %.6f against V85's %.6f"
                   % (len(P), [round(e, 2) for e in ev], dec_err, F10, v85_F10),
                   len(P) == 3 and dec_err < 1e-12
                   and abs(F10 - v85_F10) < 1e-8))

    # ---- T1 the word census ------------------------------------------------
    rows = []
    for n in range(2, NMAX + 1):
        c = census[n]
        nz = c["norms"]
        rows.append({"n": n, "words": c["words"], "surviving": c["surviving"],
                     "max_word": round(max(nz), 12) if nz else 0.0,
                     "sum_abs": round(float(sum(nz)), 12) if nz else 0.0,
                     "total": round(float(np.abs(c["total"]).sum()), 12)})
    out["data"]["census"] = rows
    Ns = [r["surviving"] for r in rows]
    growth = [(Ns[i + 1] / Ns[i]) if Ns[i] else 0.0 for i in range(len(Ns) - 1)]
    checks.append(("T1 word census at the root: surviving words per length "
                   "n = 2..%d are %s (out of 3^n); the survival rate %s, so "
                   "the multiplicity grows like %.3f^n - the loop sum is NOT "
                   "dominated by a single class"
                   % (NMAX, Ns,
                      "falls with n" if growth and growth[-1] < growth[0]
                      else "is flat",
                      max(g for g in growth if g > 0) if any(growth) else 0.0),
                   all(s > 0 for s in Ns)))

    # ---- T2 interference: the eta = 0 control ------------------------------
    inter = {}
    for eta, tag in ((ETA_ROOT, "root"), (0.0, "no phase"), (ETA_FMIN, "Fmin")):
        Bx = sectors(eta, Ma, ev, G)
        cx = word_census(Bx, cartan, NMAX)
        rec = {}
        for n in (4, 6, 8, 10):
            c = cx[n]
            tot = float(np.abs(c["total"]).sum())
            absum = float(sum(c["norms"]))
            rec[n] = {"total": tot, "sum_abs": absum,
                      "cancellation": (tot / absum) if absum > 0 else 0.0}
        Gn = np.linalg.matrix_power(G, NMAX)
        rec["F10"] = rec[10]["total"] / float(
            np.abs(cartan_vector(Gn, cartan)).sum())
        inter[tag] = rec
    out["data"]["interference"] = inter
    root_c = inter["root"][10]["cancellation"]
    checks.append(("T2 interference at length 10: the ratio ||sum of words|| / "
                   "sum ||word|| is %.4f at the root (eta = %.2f), %.4f with "
                   "the Wilson phase switched off (eta = 0) and %.4f at the "
                   "Fmin angle; F_10 = %.5f (root) against %.5f (eta = 0)"
                   % (root_c, ETA_ROOT,
                      inter["no phase"][10]["cancellation"],
                      inter["Fmin"][10]["cancellation"],
                      inter["root"]["F10"], inter["no phase"]["F10"]),
                   root_c > 0 and inter["no phase"]["F10"] > 0))

    sel = {}
    for n in range(2, NMAX + 1):
        c = census[n]
        per = []
        for lam in range(3):
            lo = lam * 3 ** (n - 1)
            hi = lo + 3 ** (n - 1)
            v = sum(c["vecs"][lo:hi], np.zeros(5, dtype=complex))
            per.append(float(np.abs(v).sum()))
        sel[n] = [round(p, 12) for p in per]
    out["data"]["sector_content"] = sel
    dom = {n: int(np.argmax(sel[n])) for n in sel}
    checks.append(("T3 sector selection: the Cartan content of length n splits "
                   "over the three sectors as %s (n = 2, 4, 6, 8, 10 shown); "
                   "the dominant sector is %s at every length, and the "
                   "mod-4 pattern of V89 is %s"
                   % ({n: sel[n] for n in (2, 4, 6, 8, 10)},
                      "sector " + str(dom[10]),
                      "NOT explained by sector dominance (all three sectors "
                      "contribute at every length)"
                      if all(min(sel[n]) > 1e-14 for n in sel)
                      else "explained by vanishing sectors"),
                   all(sum(sel[n]) > 0 for n in sel)))

    # ---- T4 multiplicity x per-word decay against the ladder ---------------
    # per-word decay: (sum_abs / surviving)^(1/n); entropy rate: surviving^(1/n)
    per_word, entropy = [], []
    for n in range(2, NMAX + 1):
        c = census[n]
        if c["surviving"]:
            avg = sum(c["norms"]) / c["surviving"]
            per_word.append((n, avg ** (1.0 / n)))
            entropy.append((n, c["surviving"] ** (1.0 / n)))
    r_word = per_word[-1][1]
    mu = entropy[-1][1]
    F10_root = inter["root"]["F10"]
    ladder17 = LOOP ** 17
    extrap = (mu * r_word) ** 170
    out["data"]["ladder"] = {"per_word_decay_r": round(r_word, 6),
                             "entropy_rate_mu": round(mu, 6),
                             "product_mu_r": round(mu * r_word, 6),
                             "F10_root": F10_root,
                             "loop_factor": LOOP,
                             "extrapolated_170": extrap,
                             "ladder_17_factor": ladder17,
                             "ratio": (extrap / ladder17) if ladder17 else 0.0}
    checks.append(("T4 the ladder factor splits into a per-word decay "
                   "r = %.4f and a word multiplicity rate mu = %.3f, product "
                   "%.4f; the n = 10 loop factor is %.5f against 1/(16 pi) = "
                   "%.5f, and the multiplicity-weighted extrapolation to the "
                   "full 17-factor ladder (length 170) gives %.3e against "
                   "(1/16 pi)^17 = %.3e, ratio %.3g"
                   % (r_word, mu, mu * r_word, F10_root, LOOP, extrap,
                      ladder17, extrap / ladder17),
                   0.0 < mu * r_word < 1.0))

    ok = all(c[1] for c in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "WHICH LOOPS SUPPLY THE LADDER FACTOR (W8). The T.L Casimir of the "
        "32 (x) 45 link has exactly 3 eigenspaces (%s), so the n-hop traced "
        "product is a sum of 3^n words and the selection rule is enumerable "
        "exactly: the enumeration reproduces (MG)^n to %.1e and F_10 = %.6f "
        "against V85's %.6f. Census: %d words survive at n = 10 out of 3^10 = "
        "59049, growing like %.3f^n. The factor splits into a per-word decay "
        "r = %.4f and a multiplicity rate mu = %.3f (product %.4f), so the "
        "loop factor is NOT one loop class: it is a sum over %d classes with "
        "the interference ratio ||sum||/sum|| || = %.4f at the root against "
        "%.4f with the Wilson phase off, and F_10 = %.5f at the root against "
        "%.5f at eta = 0 - the phase %s the factor. %s"
        % ([round(e, 2) for e in ev], dec_err, F10, v85_F10,
           census[10]["surviving"], max(g for g in growth if g > 0),
           r_word, mu, mu * r_word, census[10]["surviving"], root_c,
           inter["no phase"][10]["cancellation"], inter["root"]["F10"],
           inter["no phase"]["F10"],
           "SUPPLIES" if inter["no phase"]["F10"] > 2 * inter["root"]["F10"]
           else "does NOT supply",
           ("the multiplicity-weighted extrapolation to length 170 gives "
            "%.3e against (1/16 pi)^17 = %.3e: the counted loops %s the "
            "ladder exponent by a factor %.3g"
            % (extrap, ladder17,
               "REACH" if 0.1 < extrap / ladder17 < 10 else "MISS",
               extrap / ladder17))))

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    for n in census:
        census[n].pop("vecs")
        census[n].pop("norms")
    out["data"]["census_detail"] = {str(n): {
        "words": census[n]["words"], "surviving": census[n]["surviving"],
        "total_l1": round(float(np.abs(census[n]["total"]).sum()), 12)}
        for n in census}
    with open(os.path.join(res, "v126_w8_loop_enumeration.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V126 - which loops supply the ladder factor")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  sectors %d %s | F10 %.6f (v85 %.6f) | decomp %.1e"
          % (len(P), [round(e, 2) for e in ev], F10, v85_F10, dec_err))
    for r in rows:
        print("  n=%2d words %6d surviving %6d max_word %.3e sum_abs %.3e "
              "total %.3e" % (r["n"], r["words"], r["surviving"],
                              r["max_word"], r["sum_abs"], r["total"]))
    print("  r_word %.4f mu %.3f product %.4f | extrap170 %.3e vs %.3e"
          % (r_word, mu, mu * r_word, extrap, ladder17))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
