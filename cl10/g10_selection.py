#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the rotation angle and the coupling readings: the 1-loop data check (V56).

The step has **one** parameter, the angle $g_{10}$.  The implemented branch
reads it as $1/n_{\rm qubits} = 1/5$ (the 32-D module of the $\mathrm{Cl}(10)$
root is a $2^5$ space: **5** qubits, $10 = 2\times 5$ generators) and the
coupling as $\alpha_{10} = g_{10}^2 = 1/25$.  Both readings are **mapping
conventions** (Section 15): what the algebra fixes is the *form* of the step —
45 plane rotations sharing one angle — not its value.  This module makes the
**data check** of that convention explicit and honest.

**The 1-loop data check.**  Running the three measured $M_Z$ couplings up to
the cutoff $M_P$ at 1-loop gives a mean GUT-range coupling $\\sim 0.023$ with a
spread of $\\sim 0.5$: the three couplings do **not** meet (the open problem
O4).  The convention value $\\alpha_{10} = 0.04$ sits within a factor **1.73**
of that mean — the **right order**, nothing more.

**What the scan shows.**  Scanning the reading $\\alpha = 1/n^2$ over $n$, the
value closest to the 1-loop mean is $n = 7$, **not** the root $n = 5$: the data
confirm the **order** only, they do not select the reading.  The exact coupling
values stay with the open problem O4 — the meeting needs the six threshold
coefficients as inputs (Section 10.4).
"""

import numpy as np

from . import config as _config
from . import engine as E
from . import unification

# the Planck scale: the UV cutoff (NOT identified with M_GUT)
M_P = unification.M_P
M_Z = unification.M_Z


def _backend():
    return E._BACKEND


def _a10_candidates(n_min=2, n_max=8):
    r"""The 0-parameter candidate $\alpha_{10} = 1/n^2$ for $n = 2\ldots 8$
    (the possible 0-parameter qubit counts of a $2^n$ module)."""
    return {n: 1.0 / (n * n) for n in range(n_min, n_max + 1)}


def gut_1loop():
    """The 1-loop running **up** to $M_P$: the mean GUT-range coupling and the
    (non-unification) spread (the open problem O4)."""
    a, spread = unification.sm_running_to(M_P)
    mean = float(np.mean(a))
    return {
        "a_at_M_P": [float(x) for x in a],
        "mean_GUT": mean,
        "spread": float(spread),
    }


def g10_selection():
    """The V56 result: the implemented angle $g_{10} = 1/5$ (the root reading)
    against the 1-loop GUT coupling (the data check).  Returns the angle and
    the $\\alpha_{10} = g^2$ reading, the 1-loop GUT coupling, the ratio, and
    the candidate scan over the reading $\\alpha = 1/n^2$."""
    root = unification.g10_from_root()
    gut = gut_1loop()
    a10_0param = float(root["alpha_10"])            # 1/25
    mean_gut = gut["mean_GUT"]
    ratio = float(a10_0param / mean_gut)            # ~ 1.73

    # the candidate scan over the reading alpha = 1/n^2: which n is closest to the 1-loop mean
    cands = _a10_candidates()
    scan = {n: float(abs(a10 / mean_gut - 1.0)) for n, a10 in cands.items()}
    best_n = min(scan, key=scan.get)
    n_root = int(root["n_qubits"])
    # the honesty: the scan points at n = 7 (1/49 = 0.0204), NOT n = 5
    # (1/25 = 0.04).  n = 5 is NOT the best fit: it is the number of the root
    # (the 32-D module, 10 = 2 x 5).  The data-check therefore does not fix n:
    # it confirms only the ORDER.  That is exactly what makes the 1/n_qubits
    # reading a MAPPING CONVENTION and not a prediction.
    gap_root = scan[n_root]
    gap_best = scan[best_n]
    fit_prefers_root = bool(best_n == n_root)

    one_line = (
        "the implemented angle g_10 = 1/5 (the 1/n_qubits reading of the root "
        "module: the 32-D space, 5 qubits) gives alpha_10 = g^2 = 1/25 = 0.04, "
        "within a factor %.2f of the 1-loop GUT-range mean (%.3f at M_P; the "
        "three couplings do not meet, the open problem O4).  HONEST: the "
        "candidate scan picks n = %d (alpha = 1/%d = %.4f, gap %.3f), NOT the "
        "root n = %d (gap %.3f) — the data check confirms the ORDER only, it "
        "does not select n; reading 1/n_qubits as the coupling is a mapping "
        "convention, not a fit" % (ratio, mean_gut, best_n, best_n, cands[best_n],
                                    gap_best, n_root, gap_root)
    )
    return {
        "backend": _backend(),
        "n_qubits": n_root,
        "g10": float(root["g10"]),
        "alpha_10_0param": a10_0param,
        "g10_matches_sm_qca": bool(root["g10_matches_sm_qca"]),
        "gut_1loop": gut,
        "ratio_to_1loop": ratio,
        "best_n_1loop": int(best_n),
        "n_closest_to_1loop": int(best_n),
        "candidate_scan_gap": {int(k): float(v) for k, v in scan.items()},
        "gap_root_n5": float(gap_root),
        "gap_best_n": float(gap_best),
        "fit_prefers_root": fit_prefers_root,
        "data_check_selects_n": fit_prefers_root,
        "one_line": one_line,
    }


def verify(tol=_config.CFG.tolerances.strict):
    """Checks: (i) the implemented angle is the root reading $1/n_{\rm qubits}$
    (the 32-D module, 5 qubits) and matches the SM QCA step; (ii) the
    $\alpha_{10} = g^2$ reading is within a factor 3 of the 1-loop GUT-range
    mean (the order check; the values are the open problem O4)."""
    r = g10_selection()
    ok = (
        r["g10_matches_sm_qca"]
        and (1.0 / 3.0) <= r["ratio_to_1loop"] <= 3.0
        and r["n_qubits"] == 5
    )
    return {"r": r, "ok": bool(ok)}


def main():
    import sys
    r = verify()
    d = r["r"]
    print("=" * 72)
    print("CL10 - the rotation angle and the coupling readings (the 1-loop data check, O4)")
    print("=" * 72)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  the implemented angle g_10 = 1/{d['n_qubits']} = {d['g10']:.3f}   "
          f"(the 1/n_qubits reading: a CONVENTION, Section 15)")
    print(f"  the rotation-squared reading alpha_10 = g^2 = 1/25 = {d['alpha_10_0param']:.4f}   (a convention)")
    g = d["gut_1loop"]
    print(f"  1-loop GUT coupling (mean at M_P) = {g['mean_GUT']:.4f}   (spread {g['spread']:.2f}, the open problem O4)")
    print(f"  ratio implemented / 1-loop = {d['ratio_to_1loop']:.3f}   (within a factor ~2 = the right order)")
    print(f"  candidate scan alpha = 1/n^2:  n closest to the 1-loop GUT coupling = {d['n_closest_to_1loop']}")
    print(f"    -> the DATA-optimal n is {d['best_n_1loop']} (gap {d['gap_best_n']:.3f}), "
          f"NOT the root n = {d['n_qubits']} (gap {d['gap_root_n5']:.3f})")
    print(f"    -> the data check confirms the ORDER only; it does NOT select n "
          f"(selects n: {d['data_check_selects_n']}) — reading 1/n_qubits is a MAPPING CONVENTION")
    print(f"  VERDICT:  {'PASS (the implemented angle has the right ORDER; the data prefer n = 7, so 1/n_qubits is a convention and the coupling VALUES are the open problem O4)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
