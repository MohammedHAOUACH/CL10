#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V136 - the MINIMAL adjoint word that carries the Majorana number (O2).

**Where O2 stands.**  The seesaw of Section 10.3 puts the heavy Majorana scale
at M_N = M_P (1/16 pi)^{n_2step} with n_2step = 3, and Table 5 classes that 3
as a FITTED coefficient: nothing in the automaton says the Majorana mass is a
3-loop object rather than a 1-, 2- or 4-loop one.  V114 closed the junction
route to the coefficient r = kappa/w_0 (no cubic-symmetric closed junction on
the period-4 lattice); the EXPONENT was never attacked.

**The question, made finite.**  A Majorana mass for the neutrino is an operator
that carries the weight 2 lambda_nu (it maps the antineutrino weight to the
neutrino weight: Delta L = 2).  The automaton's own algebra is the full
M_32 = the 1024-D Clifford algebra of the root, graded by Clifford degree, and
the Cartan of Spin(10) grades it by weight.  So:

  (i)  does an operator of weight 2 lambda_nu exist in the algebra at all,
       and in which Clifford grade (which Spin(10) irrep) does it sit?
  (ii) what is the MINIMAL number k* of adjoint (45) factors whose product
       spans it?  k* is the number of link crossings of the object, i.e. the
       loop exponent the ladder must carry.

If k* = 3, the exponent n_2step = 3 is a number of the root, not a fit.  The
answer is a structural NO-GO, and it is the finding: Clifford multiplication
changes the grade by +-2, so every product of grade-2 (adjoint) factors is
EVEN-grade, while the Delta L = 2 operator sits in ODD grade 5 (the 126-type
channel) and the 10-Higgs weight in odd grade 1.  Neither is a gauge-link
composite at any length: the seesaw exponent is not a link-word length.

**Checks (computed).**
1. The weight system of the root: the 16 weights (norm^2 5), the 45 adjoint
   weights (40 roots at norm^2 8 + 5 Cartan at 0), the 243 weights of the
   whole algebra (norm^2 <= 20), the Cartan-commutant dimension 32.
2. The Majorana operator: the weight space of 2 lambda_nu in the algebra has
   dimension exactly 1, and it is a grade-5 element - the 126-type channel of
   Spin(10) exists inside the root's own algebra, no new field.
3. The SM bookkeeping: 2 lambda_nu is NOT an SM singlet (X-charge -6,
   Y-charge -1), and NO vector (10-Higgs) weight can neutralise it: the norm
   mismatch |2 lambda_nu| = sqrt(20) vs 2 x 2 is exact.
4. The parity theorem: the multiplication table maps grade g to g +- 2 only
   (exact on all 45 x 1024 table entries), so the cumulative span S_k of
   products of the 45 has zero odd-grade content and never contains the
   Majorana direction: k* = None.
5. The same obstruction for the vector-10 (Higgs) weight, grade 1.
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
sys.path.insert(0, exp("V65_wilson"))

from cl10 import clifford as CL                              # noqa: E402
from cl10 import engine as E                                 # noqa: E402
from cl10 import sm_qca                                      # noqa: E402
import V65_wilson as V65                                     # noqa: E402

TOL = 1e-9


def clifford_basis():
    """The 1024 ordered products, with their grade."""
    g = [E.asnumpy(m) for m in CL.generators()]
    I = np.eye(32, dtype=np.complex128)
    out = [I]
    grades = [0]
    import itertools
    for r in range(1, 11):
        for combo in itertools.combinations(range(10), r):
            prod = I
            for a in combo:
                prod = prod @ g[a]
            out.append(prod)
            grades.append(r)
    return out, np.array(grades), g


def main():
    t0 = time.time()
    out = {"experiment": "V136_w2_majorana_word", "backend": E._BACKEND,
           "tol": TOL, "data": {}}
    checks = []

    # ------------------------------------------------ 1. the weight system
    B, grades, h = clifford_basis()
    B = [np.asarray(m, dtype=np.complex128) for m in B]
    Bd = np.array([m.conj().T / 32 for m in B])

    def coeffs(M):
        return np.einsum('ijk,jk->i', Bd, M)

    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Mk = np.zeros((5, 1024, 1024), dtype=np.complex128)
    for k, C in enumerate(cartan):
        for j, Ej in enumerate(B):
            Mk[k, :, j] = coeffs(C @ Ej - Ej @ C)
    herm = max(float(np.abs(Mk[k] - Mk[k].conj().T).max()) for k in range(5))
    comm = max(float(np.abs(Mk[a] @ Mk[b] - Mk[b] @ Mk[a]).max())
               for a in range(5) for b in range(5))
    A = sum(math.sqrt(p) * Mk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, U = np.linalg.eigh(A)
    rj = np.array([np.real(np.diag(U.conj().T @ (Mk[k]) @ U))
                   for k in range(5)]).T
    offdiag = max(float(np.abs(U.conj().T @ Mk[k] @ U
                              - np.diagflat(rj[:, k])).max())
                  for k in range(5))
    rr = np.round(rj, 5)
    weights = {}
    for i in range(1024):
        weights[tuple(rr[i])] = weights.get(tuple(rr[i]), 0) + 1
    n2set = sorted({round(float(np.dot(v, v)), 4) for v in weights})
    zero_dim = weights.get(tuple([0.0] * 5), 0)

    lam, lam_off = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))
    Xc = lam @ xc_gg
    Yc = lam @ y_coef
    lam_n2 = np.round(np.einsum('ij,ij->i', lam, lam), 6)
    Ts = [E.asnumpy(T).astype(np.complex128)
          for T in sm_qca.plane_generators()]
    C45 = np.array([coeffs(T) for T in Ts])
    # the 45 as weight vectors of the Cartan
    g45 = []
    for a in range(45):
        v = U.conj().T @ C45[a]
        g45.append(tuple(rr[int(np.argmax(np.abs(v)))]))
    roots45 = [g for g in g45 if np.linalg.norm(g) > 1e-8]
    n_roots = len({tuple(np.round(g, 5)) for g in roots45})
    n_cartan = 45 - len(roots45)
    root_n2 = sorted({round(float(np.dot(g, g)), 4) for g in roots45})

    checks.append(("T1 the weight system of the root: the 16 weights all have "
                   "norm^2 %s (spread %.1e), the 45 adjoint factors split into "
                   "%d root-weight vectors at norm^2 %s plus %d Cartan ones, "
                   "and the full 1024-D algebra carries %d weights with "
                   "norm^2 in %s and a %d-dimensional Cartan commutant (ad "
                   "Hermitian to %.1e, commuting to %.1e, joint diagonalisation "
                   "off-diagonal %.1e)"
                   % (sorted(set(lam_n2)), float(lam_n2.max() - lam_n2.min()),
                      len(roots45), root_n2, n_cartan, len(weights), n2set,
                      zero_dim, herm, comm, offdiag),
                   bool(lam_n2.max() - lam_n2.min() < 1e-8
                        and len(roots45) == 40
                        and n_cartan == 5 and root_n2 == [8.0]
                        and n2set[-1] == 20.0 and zero_dim == 32
                        and herm < 1e-10 and comm < 1e-10 and offdiag < 1e-9)))

    # ------------------------------------------- 2. the Majorana operator
    dbl = [s for s in range(16)
           if abs(Yc[s] + 0.5) < 1e-6 and abs(Xc[s] + 3) < 1e-6]
    nu = dbl[0]
    w_maj = 2.0 * lam[nu]
    cols = [i for i in range(1024)
            if np.linalg.norm(rr[i] - np.round(w_maj, 5)) < 1e-4]
    u = U[:, cols[0]] if cols else np.zeros(1024, dtype=np.complex128)
    grade_content = np.bincount(grades[np.abs(u) > 1e-9], minlength=11)
    g_maj = int(np.argmax(grade_content))
    checks.append(("T2 the Delta L = 2 Majorana operator EXISTS in the root's "
                   "own algebra: the weight 2 lambda_nu = %s (norm^2 %.1f) has "
                   "an algebra weight space of dimension exactly %d, and it "
                   "sits in Clifford grade %d - the %d-form channel of "
                   "Spin(10) (the 126-type object of the type-II seesaw), with "
                   "no field outside the root. The lepton doublet is the two "
                   "weights with X = -3, Y = -1/2: %s"
                   % (np.round(w_maj, 3).tolist(), float(w_maj @ w_maj),
                      len(cols), g_maj, g_maj, dbl),
                   bool(len(dbl) == 2 and len(cols) == 1
                        and grade_content[g_maj] == 32 and g_maj == 5)))

    # ------------------------------------- 3. the SM bookkeeping of the weight
    chX = float(w_maj @ xc_gg)
    chY = float(w_maj @ y_coef)
    # the 10-Higgs: a vector weight (norm^2 4) of the algebra
    # can a vector (10-Higgs) weight neutralise the Majorana weight?
    vec_all = [rr[i] for i in range(1024)
               if abs(float(np.dot(rr[i], rr[i])) - 4.0) < 1e-6]
    best = min((float(np.linalg.norm(w_maj + 2.0 * v)) for v in vec_all))
    checks.append(("T3 the Majorana weight is NOT an SM singlet and the 10-Higgs "
                   "cannot make it one: its X-charge is %.3f and its Y-charge "
                   "%.3f, and over all %d vector weights (norm^2 4) the closest "
                   "two-insertion neutralisation leaves a residual weight of "
                   "norm %.3f - the norm mismatch is exact (|2 lambda_nu| = "
                   "%.3f against 2 x 2 = 4), so the Delta L = 2 neutralisation "
                   "must come from the grade-5 (126-type) channel itself, not "
                   "from the Higgs of Section 11.5"
                   % (chX, chY, len(vec_all), best, float(np.linalg.norm(w_maj))),
                   bool(abs(chX) > 1e-6 and abs(chY) > 1e-6 and len(vec_all) > 0
                        and best > 1e-3)))

    # ------------------------------------------ 4. the minimal adjoint word
    import itertools
    combos = []
    for r in range(0, 11):
        for c in itertools.combinations(range(10), r):
            combos.append(c)
    index = {c: i for i, c in enumerate(combos)}

    def mul_right(S, a, b):
        seq = list(S) + [a, b]
        sign = 1.0
        for i in range(1, len(seq)):
            j = i
            while j > 0 and seq[j - 1] > seq[j]:
                seq[j - 1], seq[j] = seq[j], seq[j - 1]
                sign = -sign
                j -= 1
        red = []
        for x in seq:
            if red and red[-1] == x:
                red.pop()
            else:
                red.append(x)
        return 1j * sign, index[tuple(red)]

    pairs = [(a, b) for a in range(10) for b in range(a + 1, 10)]
    order_ok = all(np.allclose(Ts[k], 1j * (h[a] @ h[b]))
                   for k, (a, b) in enumerate(pairs))
    MAP = []
    for (a, b) in pairs:
        dst = np.zeros(1024, dtype=np.int64)
        cf = np.zeros(1024, dtype=np.complex128)
        for i, S in enumerate(combos):
            c, j = mul_right(S, a, b)
            dst[i] = j
            cf[i] = c
        MAP.append((dst, cf))

    def right_apply(vec, Aidx):
        dst, cf = MAP[Aidx]
        o = np.zeros(1024, dtype=np.complex128)
        np.add.at(o, dst, vec * cf)
        return o

    def step(Q, S):
        """Orthonormal basis of the new products {x T_A} modulo span S."""
        cand = np.array([right_apply(row, A)
                         for row in Q for A in range(45)])
        add = []
        for i in range(0, cand.shape[0], 4000):
            blk = cand[i:i + 4000]
            blk = blk - blk @ S.conj().T @ S
            qq, rr_ = np.linalg.qr(blk.T)
            d = np.abs(np.diag(rr_))
            a = qq[:, d > 1e-8].conj().T
            if a.shape[0]:
                add.append(a)
                S = np.vstack([S, a])
            if S.shape[0] >= 1024:
                break
        return (np.vstack(add) if add else np.zeros((0, 1024),
                                                    dtype=np.complex128)), S

    q, _ = np.linalg.qr(C45.T)
    P = q.conj().T                                   # exactly-1 span (45)
    S = P.copy()                                     # cumulative span
    dims = {}
    content = {}
    odd = {}
    kstar = None
    for k in range(1, 7):
        dims[k] = int(S.shape[0])
        content[k] = float(np.abs(S @ u.conj()).max())
        odd[k] = float(np.abs(S[:, grades % 2 == 1]).max())
        if content[k] > 1e-8:
            kstar = k
            break
        P, S = step(P, S)

    # the exact parity statement on the multiplication table itself
    parity_ok = True
    for Aidx in range(45):
        dst, cf = MAP[Aidx]
        keep = np.abs(cf) > 0
        if not np.all((grades[dst[keep]] - grades) % 2 == 0):
            parity_ok = False
            break

    # control 2: the vector-10 (Higgs) weight, norm^2 4 - is it an adjoint
    # composite at all?
    w_vec = next((rr[i] for i in range(1024)
                  if abs(float(np.dot(rr[i], rr[i])) - 4.0) < 1e-6), None)
    k_vec = None
    if w_vec is not None:
        cd = [i for i in range(1024)
              if np.linalg.norm(rr[i] - w_vec) < 1e-4]
        Ws = U[:, cd]                                  # the whole weight space
        gc_v = np.zeros(11)
        for j in range(Ws.shape[1]):
            gc_v += np.bincount(grades[np.abs(Ws[:, j]) > 1e-9], minlength=11)
        g_vec = sorted(int(g) for g in range(11) if gc_v[g] > 1e-9)
        q, _ = np.linalg.qr(C45.T)
        Pv = q.conj().T
        Sv = Pv.copy()
        for k in range(1, 7):
            if float(np.abs(Sv @ Ws).max()) > 1e-8:
                k_vec = k
                break
            Pv, Sv = step(Pv, Sv)

    checks.append(("T4 the parity theorem: the Clifford multiplication table "
                   "maps grade g to grade g +- 2 only (exact on all 45 x 1024 "
                   "entries: %s), so every product of adjoint (grade-2) factors "
                   "is EVEN-grade; the cumulative span S_k has dimensions %s "
                   "with odd-grade content %s and content along the Majorana "
                   "direction %s, hence k* = %s - the Delta L = 2 operator "
                   "(grade %d, odd) is NEVER an adjoint composite"
                   % (parity_ok, [dims[k] for k in sorted(dims)],
                      ["%.1e" % odd[k] for k in sorted(odd)],
                      ["%.1e" % content[k] for k in sorted(content)], kstar,
                      g_maj),
                   bool(order_ok and parity_ok and kstar is None
                        and max(odd.values()) < 1e-12)))

    checks.append(("T5 the same obstruction for the Higgs: the vector-10 weight "
                   "(norm^2 4) has a weight space of dimension %d carried "
                   "entirely by ODD Clifford grades %s, and it is reached at "
                   "k = %s, while the adjoint channel (grade 2, even) is at "
                   "k = 1 by construction - the two odd channels of the ladder, "
                   "the 10-Higgs and the 126-type Majorana operator, are "
                   "matter-layer objects, not gauge-link composites; the fitted "
                   "n_2step = 3 is therefore NOT a link-word length"
                   % (len(cd), g_vec, k_vec),
                   bool(kstar is None and k_vec is None and len(g_vec) >= 1
                        and all(g % 2 == 1 for g in g_vec))))

    n2step_fitted = 3
    derived = (kstar == n2step_fitted)
    out["data"] = {
        "weights": {"n_weights": len(weights), "norm2_set": n2set,
                    "cartan_commutant_dim": zero_dim,
                    "roots45": n_roots, "cartan45": n_cartan,
                    "root_norm2": root_n2},
        "majorana": {"weight": np.round(w_maj, 6).tolist(),
                     "norm2": float(w_maj @ w_maj),
                     "weight_space_dim": len(cols),
                     "grade": g_maj, "grade_components": int(grade_content[g_maj]),
                     "X_charge": chX, "Y_charge": chY,
                     "doublet_indices": dbl},
        "word": {"dims_Sk": dims, "content_Sk": content, "odd_Sk": odd,
                 "parity_table_ok": bool(parity_ok), "k_star": kstar,
                 "k_vector10": k_vec, "grades_vector10": g_vec,
                 "weight_space_dim_vector10": len(cd),
                 "generator_order_ok": bool(order_ok)},
        "n_2step_fitted": n2step_fitted,
        "n_2step_derived": bool(derived),
    }
    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    out["status"] = ("THE DELTA L = 2 OPERATOR IS NOT A GAUGE-LINK COMPOSITE "
                     "(O2). The operator exists inside the root: the weight "
                     "2 lambda_nu (norm^2 %.1f) has a 1-dimensional algebra "
                     "weight space sitting in Clifford grade %d - the 126-type "
                     "channel - with no field outside the root. But Clifford "
                     "multiplication changes grade by +-2, so products of the "
                     "45 adjoint link factors are always EVEN-grade (exact on "
                     "the whole table) and the span S_k (dims %s) never contains "
                     "it: k* = %s for every k, and the vector-10 Higgs weight "
                     "(grades %s, all odd) is likewise unreachable (k = %s). "
                     "The ladder exponent n_2step = 3 is therefore NOT the "
                     "length of an adjoint word: the seesaw operator belongs "
                     "to the odd (matter) half of the Clifford module, and O2 "
                     "stays open"
                     % (float(w_maj @ w_maj), g_maj,
                        [dims[k] for k in sorted(dims)], kstar, g_vec, k_vec))
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v136_w2_majorana_word.json"), "w") as f:
        json.dump(out, f, indent=1)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  k* = %s | dims %s | content %s | wall %.1f s"
          % (kstar, dims, {k: round(content[k], 6) for k in content},
             out["wall_s"]))
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
