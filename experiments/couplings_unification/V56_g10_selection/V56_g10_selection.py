#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V56 — the 0-parameter GUT coupling g_10: the 1-loop data check (W6).

The W6 residual (the "0-parameter is partially conventional" point): the choice
to read $1/n_{\rm qubits}$ as the GUT coupling is a **mapping convention**.  This
experiment makes the **1-loop data check** of that convention explicit and
honest.

The 0-parameter $g_{10} = 1/5$ (a **number of the root**: the 32-D module,
5 qubits) gives $\alpha_{10} = 1/25 = 0.04$.  The 1-loop running up to $M_P$
gives a GUT coupling (mean of the three) of $\sim 0.023$ (the 1-loop does not
fully unify — the spread is the W4 named residual).  The 0-parameter
$\alpha_{10} = 1/25$ is **within a factor 1.73** of this 1-loop GUT coupling —
the **right order**.

The 0-parameter content: the $g_{10} = 1/5$ (a root number).  The named
residual: the 1-loop matching (W4) and the exact GUT coupling.

**The honest scan (this experiment states it explicitly).**  Scanning the
0-parameter candidates $\\alpha = 1/n^2$ ($n = 2\\ldots8$) against the 1-loop
GUT coupling, the **data-optimal** value is $n = 7$ ($1/49 = 0.0204$), **not**
the root $n = 5$ ($1/25 = 0.04$): the data-optimal $n$ is $6.4\times$ closer to
the measured coupling than the root value.  The data check therefore confirms
the **order only** — it does **not** select $n$.  This is precisely what makes
the reading $1/n_{\rm qubits}$ a **mapping convention** rather than a fit:
$n = 5$ is chosen because it is the number of qubits of the 32-D module
($10 = 2 \times 5$), not because the data prefer it.
"""
import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import g10_selection

TOL = 1e-12


def main():
    t0 = time.time()
    r = g10_selection.verify()
    d = r["r"]
    ok = r["ok"]
    g = d["gut_1loop"]

    verdict = {
        "id": "V56_g10_selection",
        "tier": 6,
        "property": "P5,P6",
        "residual": True,   # the 1-loop matching (W4) and the exact GUT coupling are the named residual
        "pass": bool(ok),
        "status": (f"the_0_parameter_GUT_coupling_g_10=1/5_is_a_number_of_the_root_the_32D_module_"
                   f"5_qubits_giving_alpha_10=1/25=0.04_within_a_factor_{d['ratio_to_1loop']:.2f}_of_the_1_loop_GUT_"
                   f"coupling_{d['gut_1loop']['mean_GUT']:.3f}_at_M_P_the_1_loop_matching_is_the_W4_named_residual"),
        "tol": TOL,
        "data": {
            "backend": d["backend"],
            "n_qubits": d["n_qubits"],
            "g10": d["g10"],
            "alpha_10_0param": d["alpha_10_0param"],
            "g10_matches_sm_qca": d["g10_matches_sm_qca"],
            "a_at_M_P": g["a_at_M_P"],
            "mean_GUT": g["mean_GUT"],
            "spread": g["spread"],
            "ratio_to_1loop": d["ratio_to_1loop"],
            "n_closest_to_1loop": d["n_closest_to_1loop"],
            "best_n_1loop": d["best_n_1loop"],
            "candidate_scan_gap": d["candidate_scan_gap"],
            "gap_root_n5": d["gap_root_n5"],
            "gap_best_n": d["gap_best_n"],
            "data_check_selects_n": d["data_check_selects_n"],
            "residual": ("the 1-loop matching (W4) and the exact GUT coupling — "
                         "the 0-parameter g_10 = 1/5 is a root number in the "
                         "right order (within a factor 2 of the 1-loop GUT "
                         "coupling)"),
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v56_g10_selection.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V56 - the 0-parameter GUT coupling g_10 (the 1-loop data check, W6)")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  0-parameter g_10 = 1/{d['n_qubits']} = {d['g10']:.3f}  (a number of the root: the 32-D module, 5 qubits)")
    print(f"  0-parameter alpha_10 = 1/25 = {d['alpha_10_0param']:.4f}")
    print(f"  1-loop GUT coupling (mean at M_P) = {g['mean_GUT']:.4f}   (a1={g['a_at_M_P'][0]:.4f} a2={g['a_at_M_P'][1]:.4f} a3={g['a_at_M_P'][2]:.4f}; spread {g['spread']:.2f})")
    print(f"  ratio 0-parameter / 1-loop = {d['ratio_to_1loop']:.3f}   (within a factor ~2 = the right order)")
    print(f"  0-parameter candidate scan:  n closest to the 1-loop GUT coupling = {d['n_closest_to_1loop']}")
    print(f"    -> the DATA-optimal n is {d['best_n_1loop']} (gap {d['gap_best_n']:.3f}), "
          f"NOT the root n = {d['n_qubits']} (gap {d['gap_root_n5']:.3f})")
    print(f"    -> the data check confirms the ORDER only; it does NOT select n "
          f"(selects n: {d['data_check_selects_n']})")
    print("    -> reading 1/n_qubits as the coupling is therefore a MAPPING CONVENTION "
          "(n = 5 is the number of qubits of the 32-D module, 10 = 2 x 5)")
    print(f"  named residual (W6):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the 0-parameter g_10 = 1/5 is a root number in the right order; the 1-loop matching is the W4 named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v56_g10_selection.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
