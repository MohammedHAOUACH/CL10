#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V204 - the closed branch is the branch: ladder concepts and limits, audited.

The closure axiom (V131) makes the clock angle a root element: $g = \pi/4$,
order 8, zero continuous unfixed numbers.  The loop machinery at that angle
(V200, V202, V203) carries an exact odd-$n$ selection rule and crosses
$1/(16\pi)$ at $n^* = 208$.  This experiment audits every ladder concept and
every stated limit of the theory under those numbers, in one place:

1. The machinery: odd-$n$ cancellation to machine zero, angle-independent word
   norms, the first crossing at $n^* = 208$, the interference deficit there.
2. The hierarchy coefficient at the crossing: the per-channel factor
   $(16\pi F_{208})^{10}$, the natural node-set readings, and the named
   residual against the observed $v / M_P$.
3. The rung census under the selection rule: the odd rungs ($3$, $17$) carry
   EXACTLY zero dynamical loop weight - the seesaw and neutrino scales are
   cell phases (V196), never loop-dynamical products; the even rungs supply a
   factor only at $n^*$; the rung masses are the cell's root fact (V196),
   reproduced from the module.
4. The seven falsifiers recomputed at the closed angle: F2 forced to
   $\alpha_X = \pi/64$, F5 moved to $\tau_p = 2.8\times10^{39}$ yr, F3/F6/F7
   verdicts unchanged (their inputs carry no angle argument), F4's channel
   defect persists at $g\{1,2\} = \pi/2$.

**Verdict.** The loop MECHANISM and the $n^* = 208$ ladder length stand with
zero free parameters; the hierarchy coefficient reads $+30\%$ over the
observed value as a named node-set residual (O8); the desert threshold
ordering is destroyed ($\tau_p = 2.8\times10^{39}$ yr, F5 out of any planned
reach); the rung masses, the neutrino tension, the baryon verdict and the
dark-matter mass are branch-independent facts of the cell and the running.
"""
import json
import math
import os
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from cl10 import baryogenesis  # noqa: E402
from cl10 import config as _config  # noqa: E402
from cl10 import coupling  # noqa: E402
from cl10 import u10pred  # noqa: E402

G10_REF = coupling.G_CLOSED          # pi/4 - the physical branch (V131)
F208_REF = 0.0196864                 # V203 regression: F at first crossing
RHO_REF = 0.978030                   # V200/V203: angle-independent norm
G_OPEN = float(_config.CFG.angle.branches.implemented)   # 1/5: the open clock


def main():
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    checks, data = [], {}

    # ---- 1. the loop machinery at the root angle --------------------------
    G, cartan, Ma, ev = hop_machinery(G10_REF)
    F, rho = loop_ratios(G, cartan, Ma, ev, G10_REF, 260)
    odd_max = max(abs(F[n]) for n in range(1, 12, 2))
    crossing = next(n for n in range(2, 261, 2) if F[n] < LOOP)
    deficit = F[crossing] / LOOP
    checks.append({
        "name": f"T1 the loop machinery at the root angle g = pi/4: odd-n"
                f" ratios cancel to {odd_max:.1e} (n <= 11), the word norm"
                f" is the angle-independent {rho:.6f}, and the FIRST"
                f" crossing of 1/(16 pi) is n* = {crossing} with"
                f" F_* = {F[crossing]:.6f} (deficit {deficit:.4f}) - the"
                f" crossing is on the n = 4k class, and its certified value by"
                f" continuity at 50 digits is F_* = 0.0196928,"
                f" deficit 0.989868 (V209 T5): the float64 reading is low by"
                f" 4e-4 relative",
        "pass": bool(odd_max < 1e-13 and crossing == 208
                     and abs(F[crossing] - F208_REF) < 1e-3 * F208_REF
                     and abs(rho - RHO_REF) < 1e-5)})
    data["machinery"] = {"odd_max_n11": float(odd_max), "rho": rho,
                         "crossing_n": int(crossing),
                         "F_cross": float(F[crossing]),
                         "deficit_vs_loop": float(deficit)}

    # ---- 2. the hierarchy coefficient at the crossing ----------------------
    base = (16 * math.pi * F[crossing]) ** 10
    readings = {"single_channel": base, "one_per_family": 3.0 * base,
                "all_twelve": 12.0 * base}
    required = COEF_OBS / base
    over = 3.0 * base / COEF_OBS - 1.0
    checks.append({
        "name": f"T2 the hierarchy coefficient of the closed ladder: per-"
                f"channel factor (16 pi F_*)^10 = {base:.5f}; node-set"
                f" readings single {base:.3f} / one-per-family"
                f" {3 * base:.3f} / all-twelve {12 * base:.1f} against"
                f" the observed {COEF_OBS:.4f}: the observed value requires"
                f" sum factor {required:.3f}, no natural set realizes it,"
                f" and the selective reading sits {100 * over:+.1f}% over -"
                f" the named O8 residual",
        "pass": bool(base < COEF_OBS < 3.0 * base
                     and 2.2 < required < 2.4 and 0.28 < over < 0.32)})
    data["coefficient"] = {"per_channel": base, "readings": readings,
                           "required_sum_factor": required,
                           "selective_reading_over": over,
                           "coef_observed": COEF_OBS}

    # ---- 3. the rung census under the selection rule ----------------------
    w3, w17 = abs(F[3]), abs(F[17])
    f10 = F[10]
    G_o, cart_o, Ma_o, ev_o = hop_machinery(G_OPEN)
    F_o, _rho_o = loop_ratios(G_o, cart_o, Ma_o, ev_o, G_OPEN, 12)
    f10_open = F_o[10]
    dyn_open = math.log10((f10_open / LOOP) ** 10)
    dyn_ratio_10 = math.log10((f10 / LOOP) ** 10)
    rung3 = u10pred.ladder_scale(3)
    rung10 = u10pred.dm_ladder_mass()
    nu = u10pred.neutrino_ladder()
    m_nu_ladder = float(nu["m_nu_eV_ladder"])
    checks.append({
        "name": f"T3 rung census under the odd-n rule: the odd rungs"
                f" carry EXACTLY zero dynamical loop weight"
                f" (|F_3| = {w3:.1e}, |F_17| = {w17:.1e}) - the seesaw"
                f" rung 3 ({rung3:.3e} GeV) and the neutrino rung 17 are"
                f" cell phases (V196), never loop-dynamical products -"
                f" while the rung-10 ratio at the CLOSED angle is not a number:"
                f" n = 10 is a 4k+2 length, both Cartan norms vanish exactly"
                f" there, and the float64 F_10 = {f10:.2f} (10^{dyn_ratio_10:.1f}"
                f" decades) is a floor ratio with no limit (V209). On the OPEN"
                f" clock the rung-10 factor IS measured,"
                f" F_10(1/5) = {f10_open:.7f} = {f10_open / LOOP:.5f} x"
                f" 1/(16 pi), i.e. {dyn_open:+.3f} decades from the ladder mass"
                f" {rung10:.2f} GeV: the dynamical loop factor lives on the open"
                f" branch (and, at the closed angle, only at n* = {crossing} on"
                f" the 4k class); the rung masses themselves are reproduced from"
                f" the module (m_B = {rung10:.2f} GeV, m_nu ladder"
                f" {m_nu_ladder:.4f} eV)",
        "pass": bool(w3 < 1e-12 and w17 < 1e-12
                     and abs(rung3 - 9.6132e13) < 1e-3 * 9.6132e13
                     and abs(rung10 - 118.57) < 5e-2
                     and abs(m_nu_ladder - 0.1462) < 1e-3
                     and abs(f10_open - 0.01918639) < 1e-6
                     and abs(dyn_open) < 0.2)})
    data["rungs"] = {"weight_F3": float(F[3]), "weight_F17": float(F[17]),
                     "F_10_closed_float64_floor_ratio": float(f10),
                     "F_10_open_branch": float(f10_open),
                     "dyn_decades_off_rung10_open": dyn_open,
                     "rung3_gev": rung3, "rung10_gev": rung10,
                     "m_nu_ladder_eV": m_nu_ladder}

    # ---- 4. the falsifiers at the closed angle ----------------------------
    a_root = coupling.alpha_std(coupling.G_CLOSED)      # pi/64
    a_obs = float(_config.CFG.desert.alpha_x)           # 1/25 (the meeting)
    tau_closed = u10pred.proton_lifetime_closed()
    lep = baryogenesis.leptogenesis_estimate()
    win = [float(v) for v in _config.CFG.physics.wimp_window_gev]
    f6_unwashed = float(lep["unwashed_ratio_to_observed"])
    f6_washout = bool(lep["agreement_survives_washout"])
    f6_ooeq = bool(lep["out_of_equilibrium"])
    checks.append({
        "name": f"T4 the seven limits at the closed branch: F1 anomaly-"
                f"free content is charge algebra, no angle argument (V26/"
                f"V91); F2 the root coupling is alpha_X = pi/64 ="
                f" {a_root:.5f}, a factor {a_root / a_obs:.3f} above the"
                f" observed 1/25 meeting - O4 tension unchanged; F3 the"
                f" ladder m_nu = {m_nu_ladder:.3f} eV stays excluded by"
                f" sum m_nu < {_config.CFG.neutrino.sum_mnu_max_ev} eV, the"
                f" 1/3 candidate unchanged (V103);"
                f" F4 the quintet channel defect persists at"
                f" g{{1,2}} = {2 * coupling.G_CLOSED:.3f} != 1 (V138);"
                f" F5 tau_p = {tau_closed:.2e} yr - the desert ordering is"
                f" destroyed and the Hyper-Kamiokande channel closes (V133);"
                f" F6 the decay at the model's own y is in equilibrium"
                f" (out_of_equilibrium = {f6_ooeq}): the unwashed value is"
                f" {f6_unwashed:.2f} x observed and the washed agreement"
                f" (ratio in [0.1, 10]: {f6_washout}) survives only as the"
                f" product of two wrong factors (V106) - no angle input;"
                f" F7 m_B = {rung10:.2f} GeV inside the"
                f" WIMP window {win[0]:.0f}-{win[1]:.0f} GeV, branch-"
                f"independent",
        "pass": bool(abs(a_root - math.pi / 64) < 1e-12
                     and a_root / a_obs > 1.2
                     and m_nu_ladder > _config.CFG.neutrino.sum_mnu_max_ev
                     and abs(2 * coupling.G_CLOSED - 1.0) > 0.1
                     and 1e39 < tau_closed < 1e40
                     and 2.5 < f6_unwashed < 3.1 and not f6_ooeq
                     and f6_washout
                     and win[0] <= rung10 <= win[1])})
    data["limits"] = {"alpha_X_root": a_root,
                      "alpha_X_root_over_observed": a_root / a_obs,
                      "tau_p_closed_yr": tau_closed,
                      "f6_unwashed_ratio_to_observed": f6_unwashed,
                      "f6_out_of_equilibrium": f6_ooeq,
                      "f6_washout_survives": f6_washout,
                      "m_B_gev": rung10}

    status = (
        f"THE CLOSED BRANCH AUDITED END TO END (the clock angle is the root"
        f" element g = pi/4, V131). (i) Machinery: odd-n cancellation"
        f" {odd_max:.1e}, norm {rho:.6f} angle-independent, first crossing"
        f" n* = {crossing} (F_* = {F[crossing]:.6f}, deficit {deficit:.4f})."
        f" (ii) Hierarchy: per-channel factor {base:.5f}, readings"
        f" {base:.3f} / {3 * base:.3f} / {12 * base:.1f} against the"
        f" observed {COEF_OBS:.4f} - the required sum factor {required:.3f}"
        f" is no natural node set, the selective reading sits"
        f" {100 * over:+.1f}% over (the named O8 residual). (iii) Rungs:"
        f" |F_3| = {w3:.1e} and |F_17| = {w17:.1e} - the odd rungs carry"
        f" zero dynamical weight, the seesaw and neutrino scales are cell"
        f" phases (V196); the rung-10 ratio at the closed angle is a floor"
        f" ratio of two exactly vanishing norms and has NO limit (4k+2 class,"
        f" V209), while on the open clock the rung-10 factor sits"
        f" {dyn_open:+.3f} decades from the ladder mass - the dynamical factor"
        f" lives on the open branch, and at the closed angle only at n* = 208"
        f" on the 4k class. (iv) Limits: F2 forced"
        f" to alpha_X = pi/64 (x{a_root / a_obs:.3f} the meeting),"
        f" F5 tau_p = {tau_closed:.2e} yr (Hyper-K channel closed), F3/F6/"
        f" F7 verdicts branch-independent, F1 algebraic, F4 defect persists."
        f" VERDICT: the ladder MECHANISM is a closed-root theorem with zero"
        f" free parameters; what nature adds at the rung-10 scale is the"
        f" node-set sum factor, and it is named, not fitted.")
    out = {"id": "V204_closed_ladder_limits_audit",
           "tier": 4,
           "property": "the closed-branch ladder concepts and the seven"
                       " limits audited at g = pi/4 (V204)",
           "residual": True,
           "programme": "V204",
           "branch": "W-NAMED",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": all(c["pass"] for c in checks)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v204_closed_ladder_limits_audit.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
