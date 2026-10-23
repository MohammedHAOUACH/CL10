#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V84 - the EQUAL-SHIFT theorem: why no automaton content can supply
the desert spectrum, and exactly what kind of content the data demand.

**Why (O4, the structural closure).**  V74 excluded the node mechanism,
V81+V82 excluded the composite content of the gauged bond layer under
every natural mass ansatz, and found the numerical identity
sum db_1 = sum db_2 = 144.056 over 16 (x) 45. This experiment finds the
THEOREM behind that number and its exact domain of validity:

    over the matter 16 itself,  (3/5) sum Y^2 = sum T_2 = 10/3 * 3/5 = 2,
    and the cross terms of 16 (x) H vanish because the 16 weights sum to
    zero in every Cartan direction; so the equality over 16 (x) H holds
    if and only if it holds over the flux set H.

The flux set of the gauged layer (the 45 weights in the V81/V82
normalization) DOES satisfy it - measured here to 1e-12 - so every
fermionic content of the automaton (full 16's tensored with any link
flux of the breaking sector) shifts U(1) and SU(2) EQUALLY.

**What the data demand.**  The identity is UNIVERSAL over complete
multiplets - verified here on the 16, the flux 45, the SO(10) vector 10
and the SU(5) 10 inside the 16 (all ratio db1/db2 = 1 to 1e-6): it is
the index-proportionality of Tr Y^2 and Tr T3^2, both proportional to
the representation index with the ratio fixed by the fundamental,
(3/5) Tr_5 Y^2 = Tr_5 T3^2. Therefore the mismatch the meeting
equations require (x12 != x23) can ONLY come from SU(5)-FRAGMENT
content or from split thresholds inside multiplets. Computed from the
same charge operators: the V64 repair's 10 + 10bar pairs shift equally
(mismatch contribution 0); its WHOLE mismatch is carried by the (1,2)
pair (db1/db2 = 3/5), which is not an SU(5) multiplet.

**Verdict.**  O4 acquires a structural form: every complete-multiplet
content at common thresholds - in particular every fermionic content
the automaton supplies (full 16's x any breaking-sector flux) - shifts
U(1) and SU(2) equally; split thresholds fail under every natural
ansatz (V82); the automaton has no fragment Weyl states. The desert
spectrum is provably not automaton content - a named input, and a
fragment-producing mechanism the named open problem.
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
os.environ["CL10_ANGLE"] = "0.2"  # built on the V65/V81/V82 desert machinery: its anchors are 1/5-branch, so pin the clock BEFORE cl10 latches it (V65 pins the same)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
sys.path.insert(0, exp("V81_link_desert_content"))
sys.path.insert(0, exp("V82_composite_desert_exclusion"))

from V82_composite_desert_exclusion import content  # noqa: E402


def main():
    t0 = time.time()
    checks = {}

    lam, gg, yc, al, comp, db_all, sm_idx = content()
    from V81_link_desert_content import sm_roots
    bL, b1, b2 = sm_roots(al, lam, gg, yc)
    bL2 = float(bL @ bL)

    def db_of(w):
        from V81_link_desert_content import c2_labels
        Ys = w @ yc
        m2 = (w @ bL) / bL2
        c2 = c2_labels(w, b1, b2)
        return np.stack([(2 / 3) * (3 / 5) * Ys ** 2,
                         (2 / 3) * m2 ** 2,
                         (2 / 3) * c2 / 8.0], axis=1)

    # ---- 1. the identity over the 16 itself -------------------------------
    db16 = db_of(lam)
    s16 = db16.sum(axis=0)
    checks["over_16"] = {"sum_db1": round(float(s16[0]), 10),
                         "sum_db2": round(float(s16[1]), 10),
                         "sum_db3": round(float(s16[2]), 10),
                         "diff_12": float(abs(s16[0] - s16[1]))}

    # ---- 2. zero-sum of the 16 weights (cross-term killer) ---------------
    checks["zero_sum_16"] = {
        "sum_Y": float(abs((lam @ yc).sum())),
        "sum_m": float(abs((lam @ bL / bL2).sum()))}

    # ---- 3. the identity over the flux set (the 45 weights) --------------
    dbal = db_of(al)
    sa = dbal.sum(axis=0)
    checks["over_flux45"] = {"sum_db1": round(float(sa[0]), 10),
                             "sum_db2": round(float(sa[1]), 10),
                             "diff_12": float(abs(sa[0] - sa[1]))}

    # ---- 4. the composite 16 (x) 45: expansion check ---------------------
    # sum over 16 x H = |H| * sum_16 + 16 * sum_H + cross (cross = 0)
    dbt = db_all.sum(axis=0)
    pred1 = 45 * s16[0] + 16 * sa[0]
    pred2 = 45 * s16[1] + 16 * sa[1]
    checks["composite_expansion"] = {
        "sum_db1_measured": round(float(dbt[0]), 8),
        "sum_db1_predicted": round(float(pred1), 8),
        "sum_db2_measured": round(float(dbt[1]), 8),
        "sum_db2_predicted": round(float(pred2), 8),
        "cross_term_db1": float(dbt[0] - pred1),
        "cross_term_db2": float(dbt[1] - pred2)}

    # ---- 5. universality: the vector 10 of SO(10) (5 + 5bar of SU(5)) ---
    # weights of the SO(10) vector: +/- e_i in the Cartan basis: the 45
    # roots are 2(e_i +- e_j): e_i = (r_a + r_b)/4.
    roots = al[np.abs(al).sum(axis=1) > 1e-12]
    e = []
    used = np.zeros(5, dtype=bool)
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            s = (roots[i] + roots[j]) / 4
            if np.count_nonzero(np.abs(s) > 1e-9) == 1 and \
               np.abs(np.abs(s).max() - 1) < 1e-9:
                k = int(np.argmax(np.abs(s)))
                if not used[k]:
                    e.append(np.sign(s[k]) * s)
                    used[k] = True
            if used.all():
                break
        if used.all():
            break
    vec = np.array(e + [-x for x in e])
    db10 = db_of(vec)
    s10 = db10.sum(axis=0)
    checks["vector10"] = {
        "sum_db1": round(float(s10[0]), 10),
        "sum_db2": round(float(s10[1]), 10),
        "ratio_db1_over_db2": round(float(s10[0] / s10[1]), 6),
        "note": "the identity holds on the vector 10 too: it is the "
                "index-proportionality of Tr Y^2 and Tr T3^2 over any "
                "complete SU(5)/Spin(10) multiplet"}

    # ---- 6. the V64 repair content: where the mismatch actually lives ---
    # identify the 5bar + 1 inside the 16 by labels: d^c: c2 != 0, m = 0,
    # Y = 1/3; L: c2 = 0, m = +-1/2; nu^c: c2 = 0, m = 0, Y = 0
    Ys16 = lam @ yc
    m16 = (lam @ bL) / bL2
    from V81_link_desert_content import c2_labels
    c216 = c2_labels(lam, b1, b2)
    is_d = (np.abs(m16) < 1e-9) & (c216 > 1e-9) & (Ys16 > 0.2)
    is_L = (c216 < 1e-9) & (np.abs(np.abs(m16) - 0.5) < 1e-9)
    is_n = (c216 < 1e-9) & (np.abs(m16) < 1e-9) & (np.abs(Ys16) < 1e-9)
    keep10 = ~(is_d | is_L | is_n)
    w10 = lam[keep10]
    db_10 = db_of(w10)
    s_10 = db_10.sum(axis=0)
    # (1,2) pair: the L weights and their conjugates - NOT an SU(5) multiplet
    w12 = lam[is_L]
    db_12 = db_of(w12)
    s_12 = db_12.sum(axis=0)
    tot64 = 6 * s_10 + 2 * s_12
    checks["V64_content"] = {
        "n_10": int(keep10.sum()),
        "db_10": [round(float(x), 6) for x in s_10],
        "ratio_db1_over_db2_for_10": round(float(s_10[0] / s_10[1]), 6),
        "db_(1,2)pair": [round(float(x), 6) for x in 2 * s_12],
        "ratio_db1_over_db2_for_(1,2)": round(float(s_12[0] / s_12[1]), 6),
        "total_3x10+3x10bar+1x(1,2)": [round(float(x), 6) for x in tot64],
        "mismatch_from_10s": float(6 * (s_10[0] - s_10[1])),
        "mismatch_from_(1,2)": float(2 * (s_12[0] - s_12[1])),
        "note": "the 10 + 10bar pairs shift U(1) and SU(2) EQUALLY "
                "(complete SU(5) multiplets): the whole db1 - db2 "
                "mismatch of the V64 repair is carried by the (1,2) "
                "pair, which is not an SU(5) multiplet"}

    # ---- 7. the mismatch the content must supply -------------------------
    from cl10 import unification as UN
    TWO_PI = 2 * math.pi
    x12 = ((1 / UN.A_Z[0] - 1 / UN.A_Z[1])
           / ((UN.B_SM[0] - UN.B_SM[1]) / TWO_PI))
    x23 = ((1 / UN.A_Z[1] - 1 / UN.A_Z[2])
           / ((UN.B_SM[1] - UN.B_SM[2]) / TWO_PI))
    checks["required_mismatch"] = {
        "x12_no_content": round(float(x12), 3),
        "x23_no_content": round(float(x23), 3),
        "note": "a content with sum db1 = sum db2 at common thresholds "
                "contributes equally to the two meeting equations and "
                "cannot reconcile x12 != x23"}

    ok = (checks["over_16"]["diff_12"] < 1e-12
          and checks["over_flux45"]["diff_12"] < 1e-12
          and abs(checks["composite_expansion"]["cross_term_db1"]) < 1e-10
          and abs(checks["composite_expansion"]["cross_term_db2"]) < 1e-10
          and checks["zero_sum_16"]["sum_Y"] < 1e-10
          and checks["zero_sum_16"]["sum_m"] < 1e-10
          and abs(checks["vector10"]["ratio_db1_over_db2"] - 1.0) < 1e-6
          and abs(checks["V64_content"]["ratio_db1_over_db2_for_10"]
                  - 1.0) < 1e-6
          and abs(checks["V64_content"]["ratio_db1_over_db2_for_(1,2)"]
                  - 0.6) < 1e-6
          and checks["V64_content"]["n_10"] == 10
          and abs(checks["V64_content"]["mismatch_from_10s"]) < 1e-12
          and abs(checks["V64_content"]["mismatch_from_(1,2)"] + 0.266666)
          < 1e-4)

    status = (
        "THE EQUAL-SHIFT THEOREM IS VERIFIED AND UNIVERSAL: over the "
        f"matter 16 sum db1 = sum db2 = {checks['over_16']['sum_db1']:.6f}; "
        "the identity also holds over the flux 45, over the SO(10) "
        "vector 10, and over the SU(5) 10 inside the 16 - it is the "
        "index-proportionality of Tr Y^2 and Tr T3^2: both are "
        "proportional to the representation index and their ratio is "
        "fixed by the fundamental, (3/5) Tr_5 Y^2 = Tr_5 T3^2. Hence "
        "EVERY complete-multiplet content at common thresholds - in "
        "particular every fermionic content of the automaton (full 16's "
        "x any breaking-sector flux) - shifts U(1) and SU(2) EQUALLY "
        "and cannot reconcile x12 = "
        f"{checks['required_mismatch']['x12_no_content']} with x23 = "
        f"{checks['required_mismatch']['x23_no_content']}. The V64 "
        "repair's mismatch is carried ENTIRELY by its (1,2) pair "
        "(db1/db2 = 3/5, mismatch_from_10s = 0, mismatch_from_(1,2) = "
        f"{checks['V64_content']['mismatch_from_(1,2)']:.4f}): the desert "
        "content MUST contain SU(5)-fragment Weyl states (or split "
        "thresholds inside multiplets, excluded for every natural "
        "ansatz by V82). The automaton supplies only complete Spin(10) "
        "multiplets: O4 is a structural no-go - the desert spectrum "
        "is provably not automaton content, and a fragment-producing "
        "mechanism is the named open problem.")

    verdict = {"experiment": "V84_equal_shift_theorem",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v84_equal_shift_theorem.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V84 - the equal-shift theorem")
    print("=" * 68)
    print(f"  over 16:      db1 = {checks['over_16']['sum_db1']}, "
          f"db2 = {checks['over_16']['sum_db2']}  "
          f"(diff {checks['over_16']['diff_12']:.1e})")
    print(f"  over flux45:  db1 = {checks['over_flux45']['sum_db1']}, "
          f"db2 = {checks['over_flux45']['sum_db2']}  "
          f"(diff {checks['over_flux45']['diff_12']:.1e})")
    print(f"  cross terms 16x45: {checks['composite_expansion']['cross_term_db1']:.1e}, "
          f"{checks['composite_expansion']['cross_term_db2']:.1e}")
    print(f"  vector 10:    db1/db2 = "
          f"{checks['vector10']['ratio_db1_over_db2']}")
    print(f"  V64 10:       db1/db2 = "
          f"{checks['V64_content']['ratio_db1_over_db2_for_10']}  "
          f"(n_10 = {checks['V64_content']['n_10']})")
    print(f"  V64 (1,2):    db1/db2 = "
          f"{checks['V64_content']['ratio_db1_over_db2_for_(1,2)']}  "
          f"mismatch = "
          f"{checks['V64_content']['mismatch_from_(1,2)']:.4f}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v84_equal_shift_theorem.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
