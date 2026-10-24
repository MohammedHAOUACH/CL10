#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V234 - the rung exponents: derived or declared (tier 5).

**The question this pre-registers against** (`review/pre_registrations.md`):
the rung labels are config integers (`cl10/cell.py`, `cell: {rung_A: 3,
rung_B: 10}`; `rung_phase` applies `int()`, so a non-integer rung is
unrepresentable); V196 shows the six-stage cell product reproduces the stage
phases `loop^3` and `loop^10`; F3's seesaw partner of a rung-10 dark matter
with a `0.05 eV` neutrino sits at rung `2.7265` - not an integer. Are
`(3, 10)` derived or declared?

Measured quantities (pre-registered):
(a) which `n` (integer grid `0..30` + the real values obtained by inverting
    every cell-step level phase) give a level phase `loop^n` matching an
    actual 36-level cell level;
(b) whether the ladder identity `m_nu = m_B^2 / m_A` holds for any pair
    `(n_A, n_B)` OTHER than `(3, 10)` within the measured scale contacts
    (`neutrino.m_atm_ev = 0.0501`, `m_sol_ev = 0.0086`, `sum_mnu_max_ev =
    0.0642`, dark matter `M_P*loop^10 = 118.571 GeV` inside
    `physics.wimp_window_gev`);
(c) whether `n = 2.72` is realized as an exact phase.

Acceptance (pre-registered): `branch: RUNGS-DERIVED` only if `(3, 10)` is
the UNIQUE pair satisfying (a)+(b); otherwise `branch: RUNGS-DECLARED` and
the rung assignment stays declared. If (c) is empty, F3 is recorded as fired
for the integer-rung ladder.
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
from cl10.paths import exp                                    # noqa: E402

from cl10 import config as _config                            # noqa: E402
from cl10 import cell as _cell                                # noqa: E402
from cl10 import u10pred as _u10                              # noqa: E402

CFG = _config.CFG
LOOP = _u10.LOOP
MP = float(CFG.physics.m_planck_gev)


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def rung_n(value):
    """the rung a scale sits at: n = log_loop(value / M_P), real-valued."""
    return math.log(value / MP, LOOP)


def main():
    t0 = time.time()
    checks, data = [], {}
    data["fitted_inputs"] = []

    v196 = art("V196_cell_definition")
    v197 = art("V197_rung_arithmetic")

    # ---- T0 instruments: rung machinery exact ------------------------------
    rp_err = max(abs(_cell.rung_phase(n) - LOOP ** n) for n in range(31))
    ms = _cell.masses()
    assert abs(ms["m_B_GeV"] - _u10.dm_ladder_mass()) == 0.0
    m_nu_id = ms["m_nu_eV_identity"]
    six_d = abs(float(v196["data"]["six_stage_phase_A"])
                - 6.0 * _cell.phase_A())
    six_ok = six_d <= 5e-20          # 2 float64 ulp at 4.7e-05
    checks.append({
        "name": f"T0 rung machinery: rung_phase(n) == loop^n exactly for "
                f"all integers 0..30 (max err {rp_err:.1e}); the two singlet "
                f"masses are ladder facts m_A = {ms['m_A_GeV']:.6e} GeV "
                f"(rung {_cell.RUNG_A}), m_B = {ms['m_B_GeV']:.6f} GeV "
                f"(rung {_cell.RUNG_B}); the identity m_B^2/m_A = "
                f"{m_nu_id:.6f} eV; the six-stage product phase is 6x the "
                f"stage phase (V196, {'OK' if six_ok else 'ERR'}) and the "
                f"stage phases are loop^3 = {_cell.phase_A():.6e}, loop^10 "
                f"= {_cell.phase_B():.6e}",
        "pass": bool(rp_err == 0.0 and six_ok
                     and _cell.RUNG_A == 3 and _cell.RUNG_B == 10)})
    data["instruments"] = {"rung_phase_max_err": rp_err,
                           "m_A_GeV": ms["m_A_GeV"], "m_B_GeV": ms["m_B_GeV"],
                           "m_nu_identity_eV": m_nu_id}

    # ---- (a) which n are realized as cell-step level phases ----------------
    ev = np.linalg.eigvals(_cell.cell_step_np())
    raw = [abs(float(a)) for a in np.angle(ev) if abs(float(a)) > 1e-9]
    angles = sorted({round(a, 9) for a in raw})
    # match every integer rung n = 0..30 against the level phases
    realized, inversions = {}, []
    for n in range(31):
        p = LOOP ** n
        d = min(abs(a - p) for a in raw)
        if d <= 1e-15:
            realized[n] = d
    for a in angles:
        inversions.append({"level_phase": a, "rung": rung_n(a * MP)})
    # resolution floor: loop^n below float64 phase resolution (2.2e-16) is
    # invisible on the step: rung 10's phase 9.71e-18 is certified only by
    # the exact construction (u10 ladder + V197's exact float rows)
    v197_rows = {r["name"]: r for r in v197["data"]["ladder_rows"]}
    rung10_off = abs(float(v197_rows["singlet B (dark matter)"]["rung_exp"])
                     - 10.0)
    n_invisible = min(n for n in range(31) if LOOP ** n < 2.3e-16)
    data["realized_rungs_float64"] = realized
    data["level_phase_inversions"] = inversions
    checks.append({
        "name": f"T1 (a): among the 36-level cell-step phases {angles}, "
                f"loop^3 (7.87e-06) IS realized, matched to "
                f"{realized.get(3, float('nan')):.1e} absolute; loop^10 "
                f"(9.71e-18) is BELOW the float64 phase resolution "
                f"(2.2e-16): every loop^n for n >= {n_invisible} folds onto "
                f"angle 0, so rung 10 is certified only by the exact "
                f"construction (u10 ladder_scale, V197 rung_exp = 10.0 to "
                f"{rung10_off:.0e}); NO other integer n in 0..30 matches "
                f"any level phase (the remaining level phases are the gauge "
                f"clock phases pi/4 and 3pi/4, inverting to rungs "
                f"{inversions[0]['rung']:.3f} and {inversions[-1]['rung']:.3f} "
                f"- not rungs of the ladder)",
        "pass": bool(3 in realized and all(
            n not in realized for n in range(31) if n not in (3,))
            and rung10_off == 0.0)})

    # ---- (b) the pair search against the measured contacts ------------------
    lo, hi = CFG.neutrino.m_sol_ev, CFG.neutrino.sum_mnu_max_ev
    wlo, whi = CFG.physics.wimp_window_gev
    pairs = []
    for nb in range(31):
        mB = MP * LOOP ** nb
        if not (wlo <= mB <= whi):
            continue
        for na in range(31):
            mnu = (mB * mB / (MP * LOOP ** na)) * 1e9   # eV
            pairs.append({
                "pair": [na, nb], "m_B_GeV": mB, "m_nu_eV": mnu,
                "in_contacts": bool(lo <= mnu <= hi),
                "vs_atm_factor": mnu / float(CFG.neutrino.m_atm_ev),
                "excluded_by_sum": bool(mnu > hi)})
    hits = [p for p in pairs if p["in_contacts"]]
    nl = _u10.neutrino_ladder()
    data["pair_search"] = {"dm_window_runks": sorted(
        {p["pair"][1] for p in pairs}),
        "candidate_pairs": pairs, "pairs_in_contacts": hits,
        "third_generation_candidate_eV": nl["m_nu_eV_third"],
        "third_status": nl["status_third"],
        "ladder_status": nl["status_ladder"]}
    checks.append({
        "name": f"T2 (b): the dark-matter contact (wimp window "
                f"[{wlo:.0f}, {whi:.0f}] GeV) pins rung_B = 10 uniquely "
                f"(rung 9 = {MP * LOOP ** 9 / 1e3:.2f} TeV, rung 11 = "
                f"{MP * LOOP ** 11:.2f} GeV); for rung_B = 10 the identity "
                f"gives m_nu = loop^(20 - n_A) M_P: n_A = 3 -> 0.146248 eV "
                f"(2.92x the atmospheric scale, ABOVE the cosmological sum "
                f"bound 0.0642 eV - the published EXCLUDED ladder, V197/F3), "
                f"n_A = 2 -> 0.002910 eV (below the solar scale 0.0086); NO "
                f"integer pair (n_A, n_B) in 0..30 lands inside the "
                f"measured contacts - not even (3, 10) itself; the only "
                f"ladder-adjacent value in contact is the third-generation "
                f"candidate 0.048749 eV, which is a free junction amplitude "
                f"(O2), not a pair identity",
        "pass": bool(len(hits) == 0 and len(pairs) <= 31)})

    # ---- (c) is n = 2.72 realizable? ----------------------------------------
    need = rung_n((MP * LOOP ** 10) ** 2 / (float(CFG.neutrino.m_atm_ev)
                                           * 1e-9) )
    trunc = _cell.rung_phase(2.72)          # int() truncation
    exact_272 = LOOP ** 2.72
    data["n_272"] = {"required_rung_A_for_atm_at_rungB10": need,
                     "rung_phase_2_72_returns": trunc,
                     "loop_pow_2_72": exact_272,
                     "truncation_to_loop_2": abs(trunc - LOOP ** 2) == 0.0,
                     "config_runks_are_ints": [CFG.cell.rung_A,
                                               CFG.cell.rung_B]}
    checks.append({
        "name": f"T3 (c): the seesaw partner of a rung-10 dark matter that "
                f"puts the identity at the atmospheric 0.0501 eV sits at "
                f"rung {need:.4f} (exactly log_loop(m_B^2 / m_atm / M_P)); "
                f"rung_phase(2.72) returns loop^{int(2.72)} (int() "
                f"truncation), NOT loop^2.72 (a {exact_272 / LOOP ** 2:.3f}x"
                f" factor away), and cell.rung_A/rung_B are config ints "
                f"(3, 10) - a non-integer rung is UNREPRESENTABLE: (c) is "
                f"EMPTY and F3 is recorded as fired for the integer-rung "
                f"ladder",
        "pass": bool(abs(trunc - LOOP ** 2) == 0.0
                     and abs(need - 2.726527059612546) < 1e-9)})

    # ---- decision -----------------------------------------------------------
    unique = (len(hits) == 1 and hits[0]["pair"] == [3, 10])
    branch = "RUNGS-DERIVED" if unique else "RUNGS-DECLARED"
    ok = all(c["pass"] for c in checks)
    status = (
        f"RUNGS {branch} (zero fitted parameters, {time.time() - t0:.0f} s): "
        f"(a) on the 36-level cell step only loop^3 is numerically realized "
        f"(matched to 1e-16 absolute); loop^10 = 9.71e-18 is below the "
        f"float64 phase resolution and certified only by the exact ladder "
        f"construction; no other integer 0..30 matches any level phase. "
        f"(b) the dark-matter window pins rung_B = 10 uniquely, and rung_B = "
        f"10 admits NO rung_A whose identity m_B^2/m_A lands inside the "
        f"measured neutrino contacts: rung_A = 3 gives 0.146 eV (excluded "
        f"by the cosmological sum bound 0.0642 eV, the published F3 "
        f"EXCLUDED ladder), rung_A = 2 gives 0.0029 eV below the solar "
        f"scale; the set of pairs satisfying (a)+(b) is EMPTY, so (3, 10) "
        f"is not the unique solution of anything - the rungs are config "
        f"declarations (cell: rung_A 3, rung_B 10; u10 already records the "
        f"exponent provenance as 'fitted (n_2step = 3)'). (c) the partner "
        f"rung 2.7265 demanded by the atmospheric scale is not realizable: "
        f"rung_phase truncates its argument to int, so loop^2.72 cannot be "
        f"produced; F3 is recorded as FIRED for the integer-rung ladder. "
        f"VERDICT: branch RUNGS-DECLARED; Table 6's fitted count stands "
        f"(the rungs were already fitted/declared); F3 fires for the "
        f"integer ladder and the neutrino mass lives on the free third-"
        f"generation junction amplitude (O2), not on a rung identity.")
    out = {"id": "V234_rung_phase_spectrum",
           "tier": 5,
           "property": "which rung exponents are realized cell-level phases "
                       "and whether the ladder identity has any integer "
                       "pair solution inside the measured scale contacts "
                       "(V234, F3)",
           "residual": True,
           "programme": "V234",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v234_rung_phase_spectrum.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
