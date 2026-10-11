#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V135 - the coupling-selection tension at the closed angle.

**Why (V131 forces this).** Section 14 states the one free parameter's
status: the selection rule g_10 = 1/n_qubits = 1/5 is a convention the data
reject - scanning alpha = 1/n^2 against the one-loop GUT coupling 0.0231 at
M_P, the best fit is n = 7 (gap 0.115), not n = 5 (gap 0.734).  The closure
axiom (V131) replaces the free angle by g = pi/4, and the same coupling
convention alpha_10 = g^2 then carries pi^2/16 in rotation units, pi/64 =
0.0491 in the SM normalisation.  Nobody has asked what the closed angle does
to that scan.  This experiment asks it.

**What is measured.**
1. V56's scan reproduced: the one-loop mean at M_P, the gaps of 1/n^2 for
   n = 2..8, best n = 7.
2. The closed value pi/64 against the same target: its gap, and the factor by
   which the closure WORSENS the coupling tension.
3. Whether the closed value is of the 1/n^2 form at all (it is not: n = 4.51).
4. The same question in angle space: the closed g = pi/4 against the rule
   g = 1/n, n = 2..8.
5. The three-way verdict: closure (algebraic), the GUT scan (data), and the
   kappa_s desert fit (V133, which absorbs any coupling) - with the numbers
   that say the closure buys the algebra and pays in coupling.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cl10 import unification as UN                                    # noqa: E402

TWO_PI = 2.0 * math.pi
MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
M_P = UN.M_P
G4 = math.pi / 4
A_QUANT = math.pi / 64.0                 # the closed angle in SM units
A_ROOT = 1.0 / 25.0                      # g = 1/5 in the same units


def main():
    t0 = time.time()
    out = {"experiment": "V135_coupling_selection_closed", "data": {}}
    checks = []

    # --------------------------------------- T1 the one-loop target and scan
    xs = math.log(M_P / MZ)
    a_mp = [1.0 / (1.0 / A_Z[i] - B_SM[i] * xs / TWO_PI) for i in range(3)]
    mean_gut = float(np.mean(a_mp))
    gaps = {n: abs(1.0 / n ** 2 / mean_gut - 1.0) for n in range(2, 9)}
    best = min(gaps, key=gaps.get)
    out["data"]["scan"] = {"a_at_M_P": a_mp, "mean_GUT": mean_gut,
                           "gaps_1_over_n2": {str(k): v for k, v in gaps.items()},
                           "best_n": best, "gap_best": gaps[best],
                           "gap_root_n5": gaps[5]}
    checks.append(("T1 V56's scan reproduced: the one-loop couplings at M_P "
                   "are %s, mean %.6f; the candidates alpha = 1/n^2 have gaps "
                   "%s, best n = %d (gap %.3f) against the root n = 5 (gap "
                   "%.3f)" % ([round(v, 6) for v in a_mp], mean_gut,
                              {k: round(v, 3) for k, v in gaps.items()},
                              best, gaps[best], gaps[5]),
                   abs(mean_gut - 0.023068) < 1e-4 and best == 7
                   and abs(gaps[5] - 0.734) < 1e-2))

    # ------------------------------------------- T2 the closed value's gap
    gap_q = abs(A_QUANT / mean_gut - 1.0)
    gap_r = abs(A_ROOT / mean_gut - 1.0)
    out["data"]["closed_value"] = {
        "alpha_10_rotation_units": G4 ** 2, "alpha_SM_units": A_QUANT,
        "gap_to_mean_GUT": gap_q, "gap_of_1_over_25": gap_r,
        "worsening_factor": gap_q / gap_r,
        "ratio_to_data_optimal": A_QUANT / (1.0 / best ** 2)}
    checks.append(("T2 at the closed angle the same convention gives alpha_10 "
                   "= (pi/4)^2 = %.6f in rotation units, pi/64 = %.6f in the "
                   "SM normalisation; its gap to the one-loop mean %.6f is "
                   "%.3f against %.3f for the old root value 1/25 - the "
                   "closure WORSENS the coupling tension by a factor %.2f and "
                   "sits %.2fx above the data-optimal 1/%d = %.6f"
                   % (G4 ** 2, A_QUANT, mean_gut, gap_q, gap_r,
                      gap_q / gap_r, A_QUANT * best ** 2, best,
                      1.0 / best ** 2),
                   abs(A_QUANT - 0.049087) < 1e-6 and gap_q > gap_r))

    # ------------------------------------------ T3 is pi/64 of the 1/n^2 form
    n_real = math.sqrt(1.0 / A_QUANT)
    near = {n: abs(1.0 / n ** 2 / A_QUANT - 1.0) for n in (4, 5, 6, 7)}
    out["data"]["form_check"] = {"n_real": n_real,
                                 "gaps_relative_to_closed": {str(k): v
                                                             for k, v in near.items()}}
    checks.append(("T3 the closed value is not of the 1/n^2 form: pi/64 = "
                   "1/n^2 needs n = %.4f, and the nearest integers give "
                   "relative gaps %s - the selection rule and the closure "
                   "axiom cannot both be true"
                   % (n_real, {k: round(v, 3) for k, v in near.items()}),
                   abs(n_real - 4.5135) < 1e-3
                   and abs(n_real - round(n_real)) > 0.4
                   and min(near.values()) > 0.15))

    # ------------------------------------------------ T4 the same in g-space
    gr = {n: abs(1.0 / n - G4) for n in range(2, 9)}
    bestg = min(gr, key=gr.get)
    out["data"]["angle_space"] = {"g_closed": G4, "gaps_1_over_n": {
        str(k): v for k, v in gr.items()}, "best_n": bestg,
        "ratio_closed_to_nearest": G4 / (1.0 / bestg),
        "ratio_closed_to_data_preferred_1_7": G4 * 7.0}
    checks.append(("T4 in angle space the closure selects g = pi/4 = %.6f; the "
                   "rule g = 1/n for n = 2..8 comes nearest at n = %d "
                   "(distance %.4f, the closed angle is %.2fx larger) and is "
                   "%.2fx the data-preferred 1/7 = %.4f - the two selection "
                   "principles disagree by a factor of order 5"
                   % (G4, bestg, gr[bestg], G4 * bestg, G4 * 7, 1.0 / 7),
                   bestg == 2 and abs(G4 * 7 - 5.4978) < 1e-3))

    # ------------------------------------------------------ T5 the verdict
    out["data"]["three_way"] = {
        "closure_says_g": G4, "closure_alpha_SM": A_QUANT,
        "gut_scan_says_alpha": 1.0 / best ** 2, "gut_scan_n": best,
        "kappa_s_absorbs_any_coupling": True,
        "gap_closed_vs_scan": gap_q, "gap_root_vs_scan": gap_r}
    checks.append(("T5 three statements about g_10 now coexist and pull apart: "
                   "the closure axiom fixes g = pi/4 (alpha = %.4f, gap %.3f "
                   "to the one-loop mean); the coupling scan prefers alpha = "
                   "1/%d^2 = %.4f (gap %.3f); and the desert meeting survives "
                   "at the closed value only through the six free kappa_s "
                   "(V133), which absorb any coupling. The closure buys the "
                   "algebra - M(0)^2 = I, no continuous phase in the dynamics "
                   "- and pays %.2fx in coupling tension: g_10 remains the "
                   "one free parameter, now with a named algebraic candidate "
                   "the data disfavours"
                   % (A_QUANT, gap_q, best, 1.0 / best ** 2, gaps[best],
                      gap_q / gap_r),
                   gap_q / gap_r > 1.5))

    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    out["status"] = ("the closed coupling pi/64 = 0.0491 has relative gap "
                     "1.128 to the one-loop mean 0.023068 where 1/25 had "
                     "0.734 (1.54x worse), is not of the 1/n^2 form "
                     "(n = 4.51), and pi/4 is 1.57x the nearest 1/n and "
                     "5.50x the data-preferred 1/7")
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v135_coupling_selection_closed.json"),
              "w") as f:
        json.dump(out, f, indent=1)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  closed pi/64 gap %.3f vs root 1/25 gap %.3f (worsening %.2fx) | "
          "scan best n=%d" % (gap_q, gap_r, gap_q / gap_r, best))
    print("  wall %.1f s" % out["wall_s"])
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
