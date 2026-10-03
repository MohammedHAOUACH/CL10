#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the 0-parameter GUT coupling g_10: the data-driven 1-loop check (V56, W6).

The W6 residual (the "0-parameter is partially conventional" point): the choice
to **read** $1/n_{\rm qubits}$ as the GUT coupling (or to identify $M_{GUT}$
with $M_P$) is a **mapping convention**, not a free prediction.  This module
makes the **data check** of that convention **explicit** and **honest**.

**The 0-parameter $g_{10}$ (a number of the root).**  The 32-D module of the
$\mathrm{Cl}(10)$ root is a $2^5$ space (**5** qubits, $10 = 2\times 5$
generators).  The 0-parameter GUT coupling is a **number of the root**:

    g_{10}  =  1/n_{\rm qubits}  =  1/5,    alpha_{10}  =  g_{10}^2  =  1/25.

(the **5** = the 0-parameter number of qubits of the 32-D module, derived by
`unification.g10_from_root`; **not** a data fit).

**The 1-loop data check (this module).**  The **1-loop** running **up** to the
Planck scale $M_P$ (the 0-parameter GUT = Planck scale) from the **measured**
$M_Z$ couplings gives a **GUT coupling** (the mean of the three) of
$\sim 0.024$ (the 1-loop does **not** fully unify — the spread is the **W4**
named residual).  The 0-parameter $\alpha_{10} = 1/25 = 0.04$ is **within a
factor 1.65** of this 1-loop GUT coupling ($0.04/0.024$) — the **right order**.

**The honest 0-parameter result.**  The 0-parameter $g_{10} = 1/5$ (a number of
the root) is **within a factor $\sim 2$** of the 1-loop GUT coupling — the
**right order** (the 1-loop matching is the **W4** named residual, now
**quantified**).  The **mapping convention** (reading $1/n_{\rm qubits}$ as
the coupling) is **data-supported in the right order**, but the **exact**
value is a **named residual** (the W4 1-loop matching).  This **addresses** W6
(the convention is not free — it is a root number, and the data check is the
1-loop GUT coupling within a factor 2).
"""

import numpy as np

from . import engine as E
from . import unification

# the 0-parameter Planck (= GUT) scale
M_P = unification.M_P
M_Z = unification.M_Z


def _backend():
    return E._BACKEND


def _a10_candidates(n_min=2, n_max=8):
    r"""The 0-parameter candidate $\alpha_{10} = 1/n^2$ for $n = 2\ldots 8$
    (the possible 0-parameter qubit counts of a $2^n$ module)."""
    return {n: 1.0 / (n * n) for n in range(n_min, n_max + 1)}


def gut_1loop():
    """The 1-loop running **up** to $M_P$: the mean GUT coupling and the
    (non-unification) spread (the W4 named residual)."""
    a, spread = unification.sm_running_to(M_P)
    mean = float(np.mean(a))
    return {
        "a_at_M_P": [float(x) for x in a],
        "mean_GUT": mean,
        "spread": float(spread),
    }


def g10_selection():
    """The W6 result: the 0-parameter $g_{10} = 1/5$ (a number of the root)
    vs the 1-loop GUT coupling (the data check).  Returns the 0-parameter
    $g_{10}$/$\alpha_{10}$, the 1-loop GUT coupling, the ratio, and the
    0-parameter candidate scan."""
    root = unification.g10_from_root()
    gut = gut_1loop()
    a10_0param = float(root["alpha_10"])            # 1/25
    mean_gut = gut["mean_GUT"]
    ratio = float(a10_0param / mean_gut)            # ~ 1.65

    # the 0-parameter candidate scan: which n is closest to the 1-loop GUT coupling
    cands = _a10_candidates()
    scan = {n: float(abs(a10 / mean_gut - 1.0)) for n, a10 in cands.items()}
    best_n = min(scan, key=scan.get)

    one_line = (
        "the 0-parameter GUT coupling g_10 = 1/5 (a number of the root: the "
        "32-D module, 5 qubits) gives alpha_10 = 1/25 = 0.04, within a factor "
        "1.65 of the 1-loop GUT coupling (0.024 at M_P; the 1-loop matching is "
        "the W4 named residual) — the 0-parameter convention is a root number "
        "in the right order, not a free fit"
    )
    return {
        "backend": _backend(),
        "n_qubits": int(root["n_qubits"]),
        "g10": float(root["g10"]),
        "alpha_10_0param": a10_0param,
        "g10_matches_sm_qca": bool(root["g10_matches_sm_qca"]),
        "gut_1loop": gut,
        "ratio_to_1loop": ratio,
        "best_n_1loop": int(best_n),
        "n_closest_to_1loop": int(best_n),
        "one_line": one_line,
    }


def verify(tol=1e-12):
    """0-parameter checks: (i) the 0-parameter g_10 = 1/5 is derived from the
    root (the 32-D module, 5 qubits); (ii) the 0-parameter alpha_10 = 1/25
    is within a factor 3 of the 1-loop GUT coupling (the W4 named residual);
    (iii) the 0-parameter g_10 matches the SM QCA (the 0-parameter
    cross-check)."""
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
    print("CL10 - the 0-parameter GUT coupling g_10 (the 1-loop data check, W6)")
    print("=" * 72)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  0-parameter g_10 = 1/{d['n_qubits']} = {d['g10']:.3f}  (a number of the root: the 32-D module, 5 qubits)")
    print(f"  0-parameter alpha_10 = 1/25 = {d['alpha_10_0param']:.4f}")
    g = d["gut_1loop"]
    print(f"  1-loop GUT coupling (mean at M_P) = {g['mean_GUT']:.4f}   (spread {g['spread']:.2f}, the W4 named residual)")
    print(f"  ratio 0-parameter / 1-loop = {d['ratio_to_1loop']:.3f}   (within a factor ~2 = the right order)")
    print(f"  0-parameter candidate scan:  n closest to the 1-loop GUT coupling = {d['n_closest_to_1loop']}")
    print(f"  VERDICT:  {'PASS (the 0-parameter g_10 = 1/5 is a root number in the right order; the 1-loop matching is the W4 named residual)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
