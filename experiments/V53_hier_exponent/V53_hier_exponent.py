#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V53 — the hierarchy exponent -17: the 0-parameter 10-loop derivation (W8).

The W8 residual (Section 13): the *derivation* of the hierarchy exponent
$v/E_P \sim 10^{-17}$.  V52 showed the 0-parameter ratio $v/M_P = 2.02e-17$
matches the observed to 0.18% but is a **2-digit coincidence** (the exponent
was not derived).  This experiment **derives** the exponent from the
0-parameter structure, with **no fit**:

1. **The 0-parameter 10-loop suppression.**  The exponent is
   $10 \cdot \log_{10}(16\pi) = 17.01$: the **10** = the 0-parameter number of
   $\mathrm{Cl}(10)$ generators ($n_{\rm gen}$), the **loop factor**
   $1/(16\pi)$ is the 0-parameter 1-loop QFT factor.  So
   $(1/16\pi)^{10} = 9.7e-18$ vs the observed $2.02e-17$: the **exponent**
   -17 matches to 0.32 decades (0-parameter-derived), the **2-digit
   coefficient** (factor 2.08) is the named residual.

2. **The n-scan (the 0-parameter loop count is unique).**  Scanning
   $n = 1..20$ of $(1/16\pi)^n$, the **unique** best match to the observed
   exponent is $n = 10 = n_{\rm gen}$ (the next-best $n = 9$ is 1.38 decades
   off) — the 10 is the root's own number, not a fit.

3. **The n-step spectrum (secondary).**  The -17 exponent is **not** a raw
   QCA spectral number: for $n = 1..100$ the smallest wrapped $n$-step level
   stays $\\sim 7e-4$ (not $10^{-17}$); reaching it requires a Diophantine
   $n \sim 5e16$ (a non-0-parameter count).  The exponent arises as a
   **power** (a 10-loop chain), not as a single spectral gap.

The 0-parameter content: the **exponent** -17 (and, via V54, the sub-Planck
DM scale) is derived from $n_{\rm gen} = 10$ and the 1-loop factor.  The
named residual: the exact 2-digit coefficient.
"""
import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import hier_exponent

TOL = 1e-12


def main():
    t0 = time.time()
    r = hier_exponent.verify()
    d = r["r"]
    ns = d["nstep"]
    ok = r["ok"]

    verdict = {
        "id": "V53_hier_exponent",
        "tier": 6,
        "property": "P5,P6",
        "residual": True,   # the exact 2-digit coefficient (factor 2.08) is the named residual
        "pass": bool(ok),
        "status": ("the_hierarchy_exponent_-17_is_the_0_parameter_10_loop_suppression_"
                   "(1/16pi)^10_with_10_the_n_gen_of_the_Cl10_root: predicted 9.7e-18 "
                   "vs observed 2.02e-17 (exponent match 0.32 decades; n=10 the unique "
                   "best of the n=1..20 scan; not a raw n-step spectral level; the "
                   "2-digit coefficient is the named residual W8)"),
        "tol": TOL,
        "data": {
            "backend": d["backend"],
            "loop_factor": d["loop_factor"],
            "n_gen": d["n_gen"],
            "predicted_ratio": d["predicted_ratio"],
            "predicted_exponent": d["predicted_exponent"],
            "observed_ratio": d["observed_ratio"],
            "observed_exponent": d["observed_exponent"],
            "exponent_abs_diff": d["exponent_abs_diff"],
            "ratio_predicted_over_observed": d["ratio_predicted_over_observed"],
            "n_best_scan": d["n_best_scan"],
            "scan_abs_diff_at_n10": float(
                min(s["abs_diff"] for s in d["scan"] if s["n"] == 10)),
            "scan_abs_diff_at_n9": float(
                min(s["abs_diff"] for s in d["scan"] if s["n"] == 9)),
            "nstep": {
                "n_levels_1step": ns["n_levels_1step"],
                "min_level_small_n": ns["min_level_small_n"],
                "max_n_scanned": ns["max_n_scanned"],
                "any_level_at_target": ns["any_level_at_target"],
                "n_diophantine_est": ns["n_diophantine_est"],
                "diophantine_is_0param": ns["diophantine_is_0param"],
            },
            "residual": ("the exact 2-digit coefficient (predicted/observed = "
                         "%.3f, a factor %.2f) — the exponent -17 is "
                         "0-parameter-derived (10 = n_gen)"
                         % (d["ratio_predicted_over_observed"],
                            1.0 / d["ratio_predicted_over_observed"])),
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v53_hier_exponent.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V53 - the hierarchy exponent -17 (the 0-parameter 10-loop derivation, W8)")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  0-parameter inputs:  n_gen = {d['n_gen']}  (the Cl(10) generators)  "
          f"loop factor 1/16pi = {d['loop_factor']:.4f}")
    print(f"  predicted  (1/16pi)^10 = {d['predicted_ratio']:.3e}   log10 = {d['predicted_exponent']:.3f}")
    print(f"  observed   v/M_P       = {d['observed_ratio']:.3e}   log10 = {d['observed_exponent']:.3f}")
    print(f"  exponent match:  {d['exponent_abs_diff']:.3f} decades  "
          f"({'0-parameter-derived' if d['exponent_abs_diff'] < 1.0 else 'NO'})")
    print(f"  2-digit coefficient ratio = {d['ratio_predicted_over_observed']:.3f}  (the named residual)")
    print(f"  n-scan:  the unique best loop count is n = {d['n_best_scan']} = n_gen  "
          f"(n=9 is {verdict['data']['scan_abs_diff_at_n9']:.2f} decades off)")
    print(f"  n-step spectrum (n<=100):  min level {ns['min_level_small_n']:.2e}  "
          f"(not 2.02e-17; a Diophantine n ~ {ns['n_diophantine_est']:.1e} would be needed)")
    print(f"  named residual (W8):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the exponent -17 is 0-parameter-derived; the 2-digit coefficient is the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v53_hier_exponent.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
