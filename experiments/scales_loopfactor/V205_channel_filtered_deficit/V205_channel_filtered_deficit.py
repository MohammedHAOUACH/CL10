#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V205 - is the O8 deficit reachable by the pinned channel? Measured no.

**Why (O8, Section 12.1).** The closed-branch hierarchy coefficient is
C = w (16 pi F_208)^10 with measured per-channel factor (16 pi F_208)^10
= 0.89968 (V203, V204): the selective reading w = 3 sits +30.0% over the
observed 2.0765, and the observed value needs per-channel 0.69219 - an
insertion amplitude 0.963879, a 3.61% deficit, against the measured
crossing deficit 1.05%. The residual lives in ONE measured number. The
one remaining structural route (never computed): the pinned M-channel
phases at pi/4, {0^32, +-pi/3^16, +-2pi/3^16, pi^32} (V203) - a Z_3 of
Floquet phases that could in principle re-phase the loop interference
and deepen the deficit with zero free parameters. This experiment
measures whether any channel-side filter can move the deficit.

**What is measured.**
1. Phase quantization at the closed clocks: the sector phases of the
   traced word sum, e^{i(g/2) ev} with ev in {-32, -12, +4}, are roots
   of unity of order dividing 4 at both closed clocks ({1, -1, -1} at
   pi/2, {1, i, i} at pi/4) while at the 1/5 reference they are
   irrational rotations (no order <= 64): at the closed clock the
   interference is a finite fourth-root filter and nothing else - the
   generic-angle interference that carried the old 0.56% agreement is
   structurally absent at every closed clock (V202).
2. Channel-side Z_3 filters on the pinned spectrum: the propagator
   trace and every nontrivial Z_3-filtered trace of the M-node channel
   are EXACTLY zero (32 + 32 e^{ipi} + 16 x 2 cos(pi/3) + 16 x 2
   cos(2pi/3): the classes cancel by the measured multiplicities):
   every filter supported by the pinned channel has empty kernel.
3. Cartan-direction census: the insertion amplitude at the crossing
   length is measured per Cartan direction; it is bimodal ({0.90,
   0.99}) - the direction choice moves the deficit to 1% or to 10%,
   never to the required 3.6%.
4. The arithmetic: the closed deficit, the required one, the 0.055%
   fingerprint of the 1/5-branch deficit, and the adjacent-length
   bracket showing no length reading realizes the observed coefficient.

**Verdict.** The O8 deficit is channel-blind: no structural filter
inside the built automaton moves it. The closed branch predicts the
hierarchy coefficient 3 x 0.989484^10 = 2.699 against the observed
2.0765; the observed value is 3 x (the generic-clock deficit)^10 to
0.055%. The named residual is the interference deficit of the crossing
length at the closed clock - nothing else.
"""
import json
import math
import os
import sys

os.environ["CL10_BACKEND"] = "numpy"  # the V203 census amplitudes are read at host-numpy precision: pin the core

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V85_closed_loop_factor"))
sys.path.insert(0, exp("V203_closed_branch_coefficient"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios, REF_G, F10_REF)
from cl10 import supercell as SC
from cl10 import u10pred

PI3 = math.pi / 3.0


def art(name, key):
    """A sibling verdict JSON, read by experiment name (cl10.paths resolves the
    category): the number another experiment measured is an input with a source."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


F208_REF = art("V203_closed_branch_coefficient", "data.closed.F_even_2_260.208")
RHO_REF = art("V131_algebraic_closure", "data.loop_machinery_pi_4.even_decay_per_step")
COEF_ONE_PER_FAMILY = art("V204_closed_ladder_limits_audit",
                          "data.coefficient.readings.one_per_family")


def main():
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    checks, data = [], {}

    # ---- 1. sector-phase quantization at the closed clocks -------------
    G4, cartan, Ma4, ev = hop_machinery(math.pi / 4)
    ev = [float(e) for e in ev]
    phase_info = {}
    for name, g in (("pi/2", math.pi / 2), ("pi/4", math.pi / 4),
                    ("1/5", REF_G)):
        ph = np.exp(1j * (g / 2.0) * np.array(ev))
        r4 = float(np.max(np.abs(ph ** 4 - 1.0)))
        minq = min(float(np.max(np.abs(ph ** q - 1.0)))
                   for q in range(1, 65))
        phase_info[name] = {"max_dev_ph4_1": r4, "min_order64_dev": minq,
                            "phases": [[round(float(p.real), 4),
                                        round(float(p.imag), 4)]
                                       for p in np.round(ph, 4)]}
    checks.append({
        "name": f"T1 the traced-word sector phases e^(i(g/2)ev), ev = "
                f"{[round(e, 1) for e in ev]}, are exactly 4th roots of "
                f"unity at the closed clocks (max |ph^4-1|: pi/2 "
                f"{phase_info['pi/2']['max_dev_ph4_1']:.1e}, pi/4 "
                f"{phase_info['pi/4']['max_dev_ph4_1']:.1e}; the values "
                f"are {{1,-1,-1}} and {{1,i,i}}) and irrational at the "
                f"reference angle (no order <= 64 within "
                f"{phase_info['1/5']['min_order64_dev']:.2f}): at the "
                f"closed clock the sector interference is a finite "
                f"fourth-root filter - the generic-angle interference "
                f"that carried the old 0.56% agreement cannot be "
                f"re-created at any closed clock (V202)",
        "pass": bool(phase_info["pi/2"]["max_dev_ph4_1"] < 1e-13
                     and phase_info["pi/4"]["max_dev_ph4_1"] < 1e-13
                     and phase_info["1/5"]["min_order64_dev"] > 1e-2)})
    data["sector_phases"] = phase_info

    # ---- 2. closed-branch crossing regression (V203 numbers) ------------
    F4, rho4 = loop_ratios(G4, cartan, Ma4, ev, math.pi / 4, 260)
    crossings = [n for n in range(2, 261) if 1e-9 < F4[n] < LOOP]
    first = crossings[0]
    F_first = F4[first]
    deficit_closed = F_first / LOOP
    checks.append({
        "name": f"T2 the closed branch regresses V203: first crossing of"
                f" 1/(16 pi) at n* = {first}, F = {F_first:.6f}, measured"
                f" deficit {deficit_closed:.6f} (rho = {rho4:.6f})",
        "pass": bool(first == 208 and abs(F_first - F208_REF) < 5e-5
                     and abs(rho4 - RHO_REF) < 1e-5)})
    data["closed_crossing"] = {"n": first, "F": float(F_first),
                               "deficit_vs_loop": float(deficit_closed),
                               "rho": float(rho4)}

    # ---- 3. pinned-channel Z_3 filters: empty kernels -------------------
    q0 = tuple(float(v) for v in SC.M_POINTS[0])
    _B, _U32, U16c = SC.sectors(math.pi / 4)
    w = SC.phases(SC.bloch(np.array(q0), U16c, g10=math.pi / 4))
    resid = [int(round(float(np.real(x)) / PI3)) % 3 for x in w]
    mult = {}
    for r in resid:
        mult[r] = mult.get(r, 0) + 1
    filt = {}
    for k in range(3):
        s = sum(np.exp(1j * float(np.real(x)) + 2j * math.pi * k * r / 3.0)
                for x, r in zip(w, resid))
        filt[f"z3_filter_k{k}"] = float(abs(s))
    checks.append({
        "name": f"T3 the pinned M-channel spectrum at pi/4 (V203:"
                f" {{0^32, +-pi/3^16, +-2pi/3^16, pi^32}}) filters to zero:"
                f" residue multiplicities {sorted(mult.items())} and the"
                f" propagator trace with every nontrivial Z_3 filter is"
                f" {filt['z3_filter_k0']:.2e}, {filt['z3_filter_k1']:.2e},"
                f" {filt['z3_filter_k2']:.2e} - every channel-side filter"
                f" supported by the pinned Floquet phases has EMPTY"
                f" KERNEL: the channel cannot re-phase the loop"
                f" interference it does not trace over",
        "pass": bool(sum(mult.values()) == 128
                     and all(v < 1e-12 for v in filt.values()))})
    data["pinned_channel"] = {"residue_multiplicities":
                              {str(k): v for k, v in sorted(mult.items())},
                              "z3_filtered_traces": filt}

    # ---- 4. per-Cartan deficit at the crossing length -------------------
    Bn, Gn = np.eye(32, dtype=np.complex128), np.eye(32, dtype=np.complex128)
    M = sum(p * m for p, m in
            zip(np.exp(1j * (math.pi / 4 / 2) * np.array(ev)), Ma4))
    MG = M @ G4
    for _ in range(first):
        Bn = Bn @ MG
        Gn = Gn @ G4
    per_c = []
    for j, C in enumerate(cartan):
        num = abs(np.trace(C @ Bn)) / 32
        den = abs(np.trace(C @ Gn)) / 32
        per_c.append(float(num / den / LOOP))
    x_req_t4 = float((COEF_OBS / 3.0) ** 0.1)
    near = min(abs(p - x_req_t4) for p in per_c)
    checks.append({
        "name": f"T4 the per-Cartan insertion amplitudes at n* = {first}:"
                f" {[round(d, 4) for d in per_c]} - bimodal {{0.90, 0.99}}:"
                f" the choice of Cartan insertion direction puts the"
                f" deficit at 1% or at 10%, never at the required"
                f" {100 * (1 - x_req_t4):.2f}% (closest direction misses by"
                f" {near:.4f} in amplitude): the Cartan-direction freedom"
                f" does not realize the observed coefficient",
        "pass": bool(all(0.85 < p < 1.02 for p in per_c) and near > 0.02)})
    data["per_cartan_insertion_at_crossing"] = per_c

    # ---- 5. the arithmetic of the named residual -------------------------
    x_closed = float(deficit_closed)
    x_req = float((COEF_OBS / 3.0) ** 0.1)
    x_branch = float(F10_REF / LOOP)
    F204 = F_first * rho4 ** -4
    F212 = F_first * rho4 ** 4
    bracket = {str(n): float(3.0 * (16 * math.pi * F) ** 10)
               for n, F in ((first - 4, F204), (first, F_first),
                            (first + 4, F212))}
    checks.append({
        "name": f"T5 the named residual, localized: closed deficit"
                f" {x_closed:.6f} -> C(3) = {3.0 * x_closed ** 10:.4f}"
                f" vs observed {COEF_OBS:.4f}; the observed value needs"
                f" insertion amplitude {x_req:.6f} (deficit"
                f" {100 * (1 - x_req):.3f}%) against the measured"
                f" {100 * (1 - x_closed):.3f}% - a {100 * (x_closed / x_req - 1):.2f}%"
                f" gap; the 1/5-branch deficit {x_branch:.6f} matches the"
                f" required value to {100 * abs(x_req / x_branch - 1):.3f}%"
                f" (fingerprint: nature still reads the generic clock);"
                f" adjacent lengths bracket the observed coefficient"
                f" between {sorted(bracket.values())[1]:.3f} and"
                f" {max(bracket.values()):.3f} - no length reading realizes"
                f" it; and by T1-T4 no channel or Cartan filter moves it",
        "pass": bool(abs(3.0 * x_closed ** 10 - COEF_ONE_PER_FAMILY) < 0.05
                     and abs(x_req / x_branch - 1) < 1e-3
                     and min(bracket.values()) < COEF_OBS
                     < max(bracket.values()))})
    data["arithmetic"] = {
        "coef_observed": COEF_OBS,
        "closed_insertion_amplitude": x_closed,
        "required_insertion_amplitude": x_req,
        "branch_1_5_insertion_amplitude": x_branch,
        "fingerprint_deviation": float(x_req / x_branch - 1),
        "adjacent_length_coefficients": bracket}

    coef3 = 3.0 * x_closed ** 10
    closed_phase_dev = max(phase_info["pi/2"]["max_dev_ph4_1"],
                           phase_info["pi/4"]["max_dev_ph4_1"])
    status = (
        f"THE O8 DEFICIT IS CHANNEL-BLIND, MEASURED. (i) The closed-clock"
        f" sector phases are exactly 4th roots of unity at both closed"
        f" angles (max |ph^4-1| {closed_phase_dev:.1e}; {{1,-1,-1}} and"
        f" {{1,i,i}}) and irrational at 1/5: at the closed clock the"
        f" sector interference is a finite fourth-root filter the pinned"
        f" channel cannot re-weight. (ii) The pinned M-channel Z_3 of"
        f" Floquet phases {{0, +-pi/3, +-2pi/3, pi}} has every filtered trace zero"
        f" (max {max(filt.values()):.1e}): the only re-phasing the pinned"
        f" channel supports has empty kernel. (iii) The per-Cartan"
        f" insertion amplitudes at n* = {first} are bimodal (0.90, 0.99):"
        f" no Cartan direction supplies the required {x_req:.6f} (closest"
        f" miss {near:.4f}). (iv) Arithmetic: closed C(3) = {coef3:.4f} vs"
        f" observed {COEF_OBS:.4f} ({100 * (coef3 / COEF_OBS - 1):+.1f}%);"
        f" required insertion {x_req:.6f} vs measured {x_closed:.6f} - a"
        f" {100 * (x_closed / x_req - 1):.2f}% gap; the required value"
        f" equals the 1/5-branch deficit to {100 * abs(x_req / x_branch - 1):.3f}%"
        f"; adjacent lengths bracket, none matches. VERDICT: no measured"
        f" structure inside the built automaton closes O8 at the closed"
        f" clock; the named residual is the crossing-length deficit itself"
        f" ({100 * (1 - x_req):.2f}% required, {100 * (1 - x_closed):.2f}%"
        f" measured), and the only measured fact about it is the"
        f" fingerprint: nature's coefficient is the generic-clock deficit"
        f" to {100 * abs(x_req / x_branch - 1):.3f}%.")
    out = {"id": "V205_channel_filtered_deficit",
           "tier": 4,
           "property": "channel-side filters of the O8 hierarchy deficit "
                       "at the closed clock: measured blind (V205)",
           "residual": True,
           "programme": "V205",
           "branch": "W-NAMED",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": all(c["pass"] for c in checks)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v205_channel_filtered_deficit.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
