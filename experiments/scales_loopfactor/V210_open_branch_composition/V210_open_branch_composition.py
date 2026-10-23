#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V210 - the composition law, measured where the functional is non-degenerate.

Section 12.1 separates two objects that a geometric extrapolation conflates:
$(MG)^{170}$, ONE chain of 170 hops on ONE particle line, and the 17 radiative
insertions of the neutrino ladder, which sit on DIFFERENT lines and therefore
multiply.  The published numbers for both were read at the root clock
$g = \pi/4$: the chain $F_{170} = 1.57$ against $(F_{10})^{17} = 3.9\times10^{29}$,
a ratio $R_{17} = 4.0\times10^{-30}$.

V209 killed that arithmetic where it stands: at a closure-class angle every
$n \equiv 2 \pmod 4$ ratio is $0/0$ in exact integer arithmetic, and approaching
the angle the ratio is unbounded with poles accumulating there ($F_{210} =$
$3.7\times10^{7}$ at $g = \pi/4 + 10^{-8}$ against the float64 floor ratio
$0.71$).  $F_2$, $F_6$, $F_{10}$, $F_{170}$, $(F_{10})^{17}$, $R_{17}$ and the
"$33\times$ branch separation" are all floor artefacts.

This experiment re-reads the same two objects at the implemented OPEN clock
$g = 1/5$, where the traced-link phases of the $\mathbf{288}$ and $\mathbf{1120}$
Casimir sectors are non-degenerate (V209 T0) and the functional is a plain
measurement rather than a limit:

1. **One measurement, two precisions.** The whole ladder $F_n$, $n \le 172$, is
   evaluated at 30 and at 50 decimal digits; the relative drift is the error
   bar, and it is $10^{-13}$-level: this ladder is not a floor.
2. **The two compositions are different numbers.** The chain value is
   $F_{170} = 1.31\times10^{-21}$; the 17-fold product of the measured rung-10
   factor is $(F_{10})^{17} = (0.0191863910262)^{17} = 6.47\times10^{-30}$.
   One chain of 170 hops suppresses $2.0\times10^{8}$ times LESS than the
   product of its 17 ten-hop pieces: consecutive hops on one line interact
   through the propagation, so a chain is not a stack of factors.
3. **What the ladder actually carries at rung 17.** The multiplicative object
   is $6.4694\times10^{-30} = \mathrm{deficit}^{17}\times(1/16\pi)^{17}$ with
   the measured clock deficit $0.9644132$ - the 17th power amplifies the 3.6%
   clock deficit into a factor $0.537$ on the rung-17 suppression.
4. **The odd-$n$ selection rule is a closure-class property.** At $g = 1/5$ the
   odd chains do NOT cancel ($F_1 = 0.154$, $F_3 = 0.242987$, $F_5 = 0.104$,
   $F_7 = 0.368$), reproducing the $1/5$ entry of V206's odd-selection record;
   at the closure class they cancel exactly (V209 T1, integer arithmetic).

**Verdict.** The composition law survives, with new numbers and in the only
regime where it is measurable: chain $1.31\times10^{-21}$, product
$6.47\times10^{-30}$, ratio $2.02\times10^{8}$.  Every version of it quoted at
the root clock is a ratio of round-off floors with no limit there (V209).
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                   # noqa: E402
from cl10 import config as _config                           # noqa: E402
sys.path.insert(0, exp("V209_root_clock_degeneracy"))
from V209_root_clock_degeneracy import (LOOP, mp_ladder,     # noqa: E402
                                        skeleton)

REF_G = str(float(_config.CFG.angle.branches.implemented))   # the open clock, g = 1/5
N_CHAIN = 172                       # one chain longer than the 170-hop ladder
RUNGS = 17                          # insertions in the neutrino ladder


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published number has a
    source, it is not a decimal restated here."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


ODD_REF_V206 = art("V206_clock_tension", "data.odd_selection.1/5")
F10_REF_V85 = art("V209_root_clock_degeneracy", "data.open_branch.F_10_30dp")


def main():
    t0 = time.time()
    checks, data = [], {}
    _N_int, blocks, dens, _D, cart, planes, _TL = skeleton()

    lad = {dps: mp_ladder(REF_G, N_CHAIN, dps, planes, blocks, dens, cart)
           for dps in (30, 50)}
    assert all(lad[dps][n]["F"] is not None for dps in lad
               for n in range(1, N_CHAIN + 1))

    # ---- T1 the open-clock ladder is a measurement, not a limit ----------
    drift = max(abs(float(lad[30][n]["F"] / lad[50][n]["F"]) - 1.0)
                for n in range(1, N_CHAIN + 1))
    F10 = float(lad[50][10]["F"])
    checks.append({
        "name": f"T1 the ladder at the open clock g = 1/5 is one measurement: 30"
                f" and 50 decimal digits agree over every length n <= {N_CHAIN}"
                f" (max relative drift {drift:.1e}), and F_10 ="
                f" {F10:.13f} reproduces the V85/V203/V209 value"
                f" {F10_REF_V85:.13f} to {abs(F10 / F10_REF_V85 - 1):.1e}"
                f" - unlike the closure class (V209) nothing here is a ratio of"
                f" round-off floors",
        "pass": bool(drift < 1e-11 and abs(F10 / F10_REF_V85 - 1) < 1e-11)})
    data["precision"] = {
        "dps": [30, 50], "max_relative_drift_n_le_172": drift,
        "F_n_dps50": {str(n): float(lad[50][n]["F"])
                      for n in range(1, N_CHAIN + 1)},
        "F_n_dps30": {str(n): float(lad[30][n]["F"])
                      for n in range(1, N_CHAIN + 1)}}

    # ---- T2 chain versus product -----------------------------------------
    F170 = float(lad[50][170]["F"])
    F170_30 = float(lad[30][170]["F"])
    prod = F10 ** RUNGS
    ratio = F170 / prod
    checks.append({
        "name": f"T2 the two compositions are different measured numbers: one"
                f" chain of 170 hops on one line gives F_170 = {F170:.5e}"
                f" (30-digit value {F170_30:.5e}, drift"
                f" {abs(F170 / F170_30 - 1):.1e}), while 17 insertions on"
                f" different lines multiply to (F_10)^17 = {prod:.5e}; the"
                f" chain suppresses {ratio:.4g} times LESS than the product of"
                f" its 17 ten-hop pieces, so consecutive hops interact through"
                f" the propagation and a chain is not a stack of factors"
                f" (the root-clock entries F_170 = 1.57, (F_10)^17 = 3.9e29,"
                f" R_17 = 4.0e-30 are floor ratios, not measurements: n = 170"
                f" and n = 10 are both"
                f" 4k+2 lengths, and those ratios have no limit at the closure"
                f" angle, V209)",
        "pass": bool(abs(F170 / F170_30 - 1) < 1e-11 and ratio > 1e7
                     and F170 < 1.0 and prod < LOOP ** RUNGS * 0.7)})
    data["composition"] = {"F_170_chain": F170, "F_10": F10,
                           "product_17_x_F10": prod,
                           "chain_over_product": ratio,
                           "odd_chain_values_dps50": {
                               str(n): float(lad[50][n]["F"])
                               for n in range(1, 12, 2)}}

    # ---- T3 what the ladder carries at rung 17 ----------------------------
    deficit = F10 / LOOP
    loop17 = LOOP ** RUNGS
    ladder_factor = prod / loop17
    checks.append({
        "name": f"T3 the rung-17 suppression of the ladder, measured:"
                f" (F_10(1/5))^17 = {prod:.6e} = {ladder_factor:.5f} x"
                f" (1/16 pi)^17 ({loop17:.6e}); the identity"
                f" prod/(1/16pi)^17 = deficit^17 holds to"
                f" {abs(ladder_factor / deficit ** RUNGS - 1):.1e} with the"
                f" measured clock deficit {deficit:.7f} - the 17th power turns"
                f" the {(1 - deficit) * 100:.2f}% clock deficit into a"
                f" {(1 - ladder_factor) * 100:.1f}% correction of the rung-17"
                f" exponent ruler",
        "pass": bool(abs(ladder_factor / deficit ** RUNGS - 1) < 1e-12
                     and 0.4 < ladder_factor < 0.7)})
    data["rung17"] = {"deficit_clock": deficit, "loop_factor_pow_17": loop17,
                      "measured_product_pow_17": prod,
                      "correction_of_the_ruler": ladder_factor}

    # ---- T4 the odd-n rule is a closure-class property --------------------
    odd = {n: float(lad[50][n]["F"]) for n in range(1, 12, 2)}
    checks.append({
        "name": f"T4 the odd-n cancellation is a closure-class property, not a"
                f" law of the ladder: at g = 1/5 the odd chains are finite and"
                f" non-zero ({', '.join(f'F_{n} = {odd[n]:.4f}' for n in (1, 3, 5, 7))}),"
                f" and F_3 = {odd[3]:.16f} reproduces the 1/5 entry of V206's"
                f" odd-selection record ({ODD_REF_V206}) to"
                f" {abs(odd[3] / ODD_REF_V206 - 1):.1e}; at the closure class"
                f" the same traces vanish exactly in integer arithmetic (V209"
                f" T1)",
        "pass": bool(abs(odd[3] / ODD_REF_V206 - 1) < 1e-9
                     and all(v > 1e-6 for v in odd.values()))})
    data["odd_open_branch"] = {str(n): v for n, v in odd.items()}

    status = (f"COMPOSITION LAW, MEASURED ON THE OPEN CLOCK (zero free"
              f" parameters, {time.time() - t0:.0f} s). At g = 1/5 the traced-"
              f"loop Cartan ladder is precision-stable end to end: 30- and 50-"
              f"digit evaluations of F_n agree to {drift:.0e} for every n <="
              f" {N_CHAIN}. The one-line chain of 170 hops gives F_170 ="
              f" {F170:.4e}; the 17 separate insertions multiply to (F_10)^17 ="
              f" {prod:.4e}; a chain therefore suppresses {ratio:.3g} times less"
              f" than the product of its pieces, which is why the ladder's 17"
              f" factors must be insertions on different lines and why a"
              f" geometric fit of (MG)^170 cannot supply them. The rung-17"
              f" suppression carried by the ladder is {prod:.4e} ="
              f" {ladder_factor:.4f} x (1/16 pi)^17: the clock deficit"
              f" {deficit:.6f} raised to the 17th power. The odd-length"
              f" cancellation used elsewhere is specific to the closure class"
              f" (odd chains at 1/5: F_3 = {odd[3]:.4f}). VERDICT: the"
              f" composition law stands on these numbers; the root-clock"
              f" entries (F_170 = 1.57, (F_10)^17 = 3.9e29, R_17 = 4.0e-30)"
              f" are floor ratios with no limit at the closure angle (V209).")
    out = {"id": "V210_open_branch_composition",
           "tier": 4,
           "property": "one-chain and multi-leg loop compositions are distinct "
                       "measured numbers on the open clock (V210)",
           "residual": True,
           "programme": "V210",
           "branch": "W-NAMED",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": all(c["pass"] for c in checks)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v210_open_branch_composition.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
