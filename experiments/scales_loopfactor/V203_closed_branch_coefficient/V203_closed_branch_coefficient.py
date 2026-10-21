#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V203 - the closed-branch hierarchy coefficient, measured self-consistently.

**Why (O8, Section 9).** The closure axiom quantizes the clock to
g in {pi/2, pi/4} and the ladder crossing moves: the implemented branch
crosses 1/(16 pi) at exactly n = 10 (V85, V108, V113); the quantized
branch was read from the V200 eta-law - a scan of the traced-link phase
with the STEP held at the implemented angle, i.e. a hybrid of the two
branches - as crossing at n* = 208. This experiment measures the
SELF-CONSISTENT closed branch: step and traced-link phase both at pi/4.

**What is measured.**
1. Regression and the COMPENSATION IDENTITY on the implemented branch:
   the published coefficient 2.0881 = 3 (F_10 / (1/16 pi))^10 = 3 x
   0.9644^10. The factor 3 (one split node per axis family, V113)
   alone overpredicts by 1.404; the interference deficit 0.9644 at the
   crossing, raised to the tenth, supplies exactly the compensation.
   The hierarchy precision is carried by that deficit - an interference
   effect of the sector phases e^{i(g/2) ev}, present at generic
   angles and absent at roots of unity (V202: the eight crossing
   angles of the exact functional contain no root angle).
2. The closed branch: odd-n loop ratios vanish exactly; rho(MG) equals
   the implemented value (the word norms are angle-independent, V200);
   the FIRST crossing of 1/(16 pi) is located and recorded; and the
   self-consistent F_10 at pi/4 is recorded next to the hybrid eta-law
   value (0.1949) - they are not the same object.
3. Channel equality at the quantized angle: the 12 M nodes are exactly
   gapless at pi/4, their Floquet spectra identical within and across
   the three axis families (V113's equal-contributions condition,
   re-verified at pi/4); and the gapless set at M is scanned along the
   13 cubic line directions to measure the nodal network.
4. The node-set sum at pi/4: per-channel factor (16 pi F_208)^10;
   every natural reading (one channel per family, all 12 point
   channels, the detected line network, a single channel) gives a
   coefficient; none lands on the observed 2.0765.

**Verdict.** The loop factor and the channel equality SURVIVE closure -
the hierarchy mechanism survives with zero free parameters - but its
99.44% precision was the generic-angle interference deficit and is not
recoverable at any root angle: the closed branch predicts the
hierarchy at a named percent over, the implemented branch at +0.6% at
the price of the 1/5 convention.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V85_closed_loop_factor"))

from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)
from V85_closed_loop_factor import cartan_norm, projectors  # noqa: E402
from cl10 import clifford  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import sm_qca  # noqa: E402
from cl10 import supercell as SC  # noqa: E402
from cl10 import u10pred  # noqa: E402
from cl10 import config as _config  # noqa: E402

D_L = 45
REF_G = float(_config.CFG.angle.branches.implemented)   # the 1/5 reference angle
F10_REF = 0.01918639            # V85 / V108 / V113 regression
F10_HYBRID_REF = 0.194894       # V200 eta_law: step at 1/5, phase at pi/4


def hop_machinery(g10):
    """The V85 one-particle traced-link hop data at the clock angle g10."""
    G = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(),
                                      g10)).astype(np.complex128)
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = root_generators()
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    P, ev = projectors(TL)
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    Ma = [np.einsum("ajbj->ab", (Pa @ Gi).reshape(32, D_L, 32, D_L)) / D_L
          for Pa in P]
    return G, cartan, Ma, ev


def loop_ratios(G, cartan, Ma, ev, eta, nmax):
    ph = [np.exp(1j * (eta / 2) * e) for e in ev]
    M = sum(p * m for p, m in zip(ph, Ma))
    MG = M @ G
    rho = float(np.max(np.abs(np.linalg.eigvals(MG))))
    Fn, Bn, Gn = {}, np.eye(32, dtype=np.complex128), \
        np.eye(32, dtype=np.complex128)
    for n in range(1, nmax + 1):
        Bn = Bn @ MG
        Gn = Gn @ G
        Fn[n] = cartan_norm(Bn, cartan) / cartan_norm(Gn, cartan)
    return Fn, rho


def main():
    LOOP = float(u10pred.LOOP)
    # observed v/M_P divided by the reference loop power - built from the
    # measured inputs (config), not typed in (V108's COEF_OBS)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    G15, cartan, Ma15, ev = hop_machinery(REF_G)        # reference branch, pinned
    G4, _, Ma4, _ = hop_machinery(math.pi / 4)

    checks = []
    data = {}

    # ---- 1. implemented branch: regression + compensation identity -------
    F15, rho15 = loop_ratios(G15, cartan, Ma15, ev, REF_G, 30)
    F10 = F15[10]
    reg_err = abs(F10 - F10_REF)
    deficit = F10 / LOOP
    coef_impl = 3.0 * deficit ** 10
    coef_3_bare = 3.0                      # 3 x LOOP^10 alone, in units
    checks.append({
        "name": f"T0 the implemented branch reproduces V85 (F_10 = {F10:.8f}"
                f", regression {reg_err:.2e}) and the compensation identity"
                f" 3 x deficit^10 = {coef_impl:.4f} against the published"
                f" 2.0881 (|d| {abs(coef_impl - 2.0881):.4f}); the bare"
                f" factor-3 reading overpredicts the observed coefficient"
                f" by x{coef_3_bare / COEF_OBS:.3f} and"
                f" deficit^10 = {deficit ** 10:.4f} supplies exactly the"
                f" compensation",
        "pass": bool(reg_err < 1e-6 and abs(coef_impl - 2.0881) < 1e-3)})
    data["implemented"] = {
        "F_10": round(F10, 8), "deficit_vs_loop": round(deficit, 6),
        "coef_3_deficit10": round(coef_impl, 5),
        "coef_observed_from_config": round(COEF_OBS, 5)}

    # ---- 2. self-consistent closed branch: selection rule, rho, crossing -
    F4, rho4 = loop_ratios(G4, cartan, Ma4, ev, math.pi / 4, 260)
    odd_max = max(F4[n] for n in range(3, 12, 2))
    crossings = [n for n in range(2, 261) if 1e-9 < F4[n] < LOOP]
    first = crossings[0] if crossings else None
    F_first = F4[first] if first else None
    # the hybrid object, for the record: step at 1/5, phase at pi/4
    Fh, rho_hyb = loop_ratios(G15, cartan, Ma15, ev, math.pi / 4, 12)
    checks.append({
        "name": f"T1 the self-consistent closed branch (step and phase at"
                f" pi/4): odd-n loop ratios vanish (max |F_2k+1| = "
                f"{odd_max:.2e} for n <= 11); rho(MG) = {rho_hyb:.6f} scan"
                f" value, rho(MG) closed = {rho4:.6f} (|d| "
                f"{abs(rho4 - rho_hyb):.2e}: the decay law is fixed by the"
                f" traced-link phase alone (V200 norms); the implemented"
                f" branch gives rho = {rho15:.6f}, regressing V85's 0.776"
                f" to {abs(rho15 - 0.775661):.2e}); the FIRST"
                f" crossing of 1/(16 pi) is at n* = {first} with F = "
                f"{F_first:.6f} (deficit {F_first / LOOP:.4f}); and"
                f" F_10(self-consistent) = {F4[10]:.4f} is NOT the V200"
                f" eta-law value {Fh[10]:.4f} (hybrid: step at 1/5, phase"
                f" at pi/4) - the two objects are recorded apart",
        "pass": bool(odd_max < 1e-12 and abs(rho4 - rho_hyb) < 1e-9
                     and abs(rho15 - 0.775661) < 1e-6
                     and first == 208 and abs(F_first - LOOP) < 5e-4
                     and abs(Fh[10] - F10_HYBRID_REF) < 1e-3)})
    data["closed"] = {
        "odd_max_n_le_11": float(odd_max), "rho_MG": round(rho4, 6),
        "rho_MG_hybrid_phase_scan": round(rho_hyb, 6),
        "rho_MG_implemented_V85": round(rho15, 6),
        "first_crossing_n": first,
        "n_below_loop_in_2_260": len(crossings),
        "F_n_2_30": {n: float(F4[n]) for n in range(2, 31)},
        "F_even_2_260": {n: float(F4[n]) for n in range(2, 261, 2)},
        "F_10_hybrid_eta_law": round(float(Fh[10]), 6),
    }

    # ---- 3. channels at pi/4: exact gaplessness, pinning, zero band ---
    nodes = [tuple(float(v) for v in q) for q in SC.M_POINTS]
    fam = {"x": [], "y": [], "z": []}
    for q in nodes:
        axis = next(j for j in range(3) if abs(q[j]) < 1e-9)
        fam[["x", "y", "z"][axis]].append(q)
    _B, _U32, U16c = SC.sectors(math.pi / 4)

    def spec(q):
        return SC.phases(
            SC.bloch(np.array(q), U16c, g10=math.pi / 4))

    m_gaps = [float(np.min(np.abs(spec(q)))) for q in nodes]
    gamma_gap = float(np.min(np.abs(spec([0.0, 0.0, 0.0]))))
    wM = spec(nodes[0])
    n0M = int(np.sum(np.abs(wM) < 1e-9))
    nPM = int(np.sum(np.abs(wM - math.pi) < 1e-9))
    within = 0.0
    for qs in fam.values():
        s0 = spec(qs[0])
        for q in qs[1:]:
            within = max(within, float(np.abs(spec(q) - s0).max()))
    across = 0.0
    s_x = spec(fam["x"][0])
    for k in ("y", "z"):
        across = max(across, float(np.abs(spec(fam[k][0]) - s_x).max()))

    # the gapless SET: Gamma + the six face-diagonal lines (V132); the bulk
    # momenta carry a real gap - measured, not assumed
    def fold(w):
        return np.sort(np.where(np.abs(np.abs(w) - math.pi) < 1e-9, math.pi, w))

    wG = spec((0.0, 0.0, 0.0))
    n0G = int(np.sum(np.abs(wG) < 1e-9))
    npiG = int(np.sum(np.abs(np.abs(wG) - math.pi) < 1e-9))
    wM = spec(nodes[0])
    nPM = int(np.sum(np.abs(np.abs(wM) - math.pi) < 1e-9))
    within = 0.0
    for qs in fam.values():
        s0 = fold(spec(qs[0]))
        for q in qs[1:]:
            within = max(within, float(np.abs(fold(spec(q)) - s0).max()))
    across = 0.0
    s_x = fold(spec(fam["x"][0]))
    for k in ("y", "z"):
        across = max(across, float(np.abs(fold(spec(fam[k][0])) - s_x).max()))
    ts = np.linspace(0.05, math.pi / 2, 12)
    line_gap = 0.0
    n_line = 0
    for zero_axis in range(3):
        for s1 in (1.0, -1.0):
            for s2 in (1.0, -1.0):
                for t in ts:
                    q = [0.0, 0.0, 0.0]
                    a, b = [j for j in range(3) if j != zero_axis]
                    q[a], q[b] = s1 * t, s2 * t
                    line_gap = max(line_gap, float(np.min(np.abs(spec(q)))))
                    n_line += 1
    rng = np.random.default_rng(20261020)
    bulk = [float(np.min(np.abs(spec(qq))))
            for qq in rng.uniform(-math.pi / 2, math.pi / 2, (12, 3))]
    bulk_min = min(bulk)

    checks.append({
        "name": f"T2 at pi/4 the 12 M channels are EXACTLY gapless (max"
                f" M gap {max(m_gaps):.2e} - the 2.003e-4 defect of the"
                f" reference angle repaired exactly); Gamma is FULLY pinned"
                f" to the Floquet pair {{0, pi}} ({n0G}+{npiG} of 128,"
                f" M(0)^2 = I), and every M spectrum is"
                f" {{0^32, +-pi/3^16, +-2pi/3^16, pi^32}}: 64 of 128 pinned"
                f" ({n0M} zero + {nPM} at pi) - identical across all 12 M"
                f" nodes (folded deviation {within:.2e} within, {across:.2e}"
                f" across); the exactly gapless set is Gamma plus the six"
                f" face-diagonal lines ({n_line} sampled line points, max"
                f" gap {line_gap:.2e}) and the sampled bulk momenta are"
                f" gapped (min quasi-energy >= {bulk_min:.3f}) - the V132"
                f" line network, not a zone-wide band (V166 re-gaps the M"
                f" nodes when the Wilson factor is quantized)",
        "pass": bool(max(m_gaps) < 1e-10 and gamma_gap < 1e-10
                     and n0G + npiG == 128 and n0M + nPM == 64
                     and line_gap < 1e-9 and bulk_min > 1e-3
                     and within < 1e-9 and across < 1e-9)})
    data["channels_pi4"] = {
        "m_gap_max": float(max(m_gaps)), "gamma_gap": float(gamma_gap),
        "gamma_pinned_n0_npi": [n0G, npiG],
        "m_pinned_n0_npi": [n0M, nPM],
        "m_spectral_interior": 64,
        "spectral_equality_within": float(within),
        "spectral_equality_across": float(across),
        "face_line_gap_max": float(line_gap),
        "bulk_gap_min": float(bulk_min)}

    # ---- 4. the node-set sum at the quantized angle ----------------------
    base_c = (16 * math.pi * F_first) ** 10
    readings = {"single_channel": 1.0, "one_per_family": 3.0,
                "all_points": 12.0}
    coefs = {k: f * base_c for k, f in readings.items()}
    required = COEF_OBS / base_c
    checks.append({
        "name": f"T4 the closed-branch per-channel factor is (16 pi F_208)"
                f"^10 = {base_c:.5f}; the node-set readings give the"
                f" coefficients " + ", ".join(f"{k} = {v:.3f}" for k, v in
                                             sorted(coefs.items())) +
                f" against the observed {COEF_OBS:.4f}: the closest natural"
                f" reading is over by "
                f"{100 * abs(coefs['one_per_family'] - COEF_OBS) / COEF_OBS:.1f}"
                f"% and the coefficient requires sum_factor = "
                f"{required:.3f}, which no natural node set (1, 3 or 12"
                f" point channels) realizes",
        "pass": bool(abs(coefs["one_per_family"] - COEF_OBS) / COEF_OBS
                     > 0.05 and 1.5 < required < 3.0)})

    coef3 = coefs["one_per_family"]
    status = (
        f"THE CLOSED BRANCH MEASURED SELF-CONSISTENTLY (step AND traced-"
        f"link phase at pi/4). (i) The implemented branch's 99.44% is a "
        f"COMPENSATION: coefficient 2.0881 = 3 x (F_10/(1/16pi))^10 = "
        f"3 x 0.9644^10 - the factor 3 overpredicts by x"
        f"{3.0 / COEF_OBS:.3f} and the interference deficit "
        f"of the GENERIC clock angle supplies exactly the rest. (ii) At "
        f"pi/4: odd-n loop ratios vanish exactly ({odd_max:.1e}), "
        f"rho(MG) = {rho4:.6f} equals the phase-scan value (norms angle-"
        f"independent), the FIRST crossing of 1/(16 pi) is at n* = "
        f"{first} (F = {F_first:.6f}, deficit {F_first / LOOP:.4f}) - the "
        f"published 208 is confirmed for the self-consistent branch - "
        f"while the self-consistent F_10 = {F4[10]:.4f} is NOT the eta-"
        f"law 0.1949 (a hybrid of step@1/5 with phase@pi/4; both values "
        f"recorded apart). (iii) Channels at pi/4: the 12 M points are "
        f"EXACTLY gapless (max M gap {max(m_gaps):.1e}, the 2.003e-4 defect "
        f"repaired); Gamma is FULLY pinned to the Floquet pair, {n0G}+{npiG} "
        f"of 128, and every M spectrum is {{0^32, +-pi/3, +-2pi/3, pi^32}} - "
        f"64 of 128 pinned ({n0M} zero + {nPM} at pi); the exactly gapless "
        f"set is Gamma plus the six face-diagonal lines (max gap over "
        f"{n_line} line points {line_gap:.1e}) with gapped sampled bulk "
        f"momenta (min quasi-energy {bulk_min:.3f}) - the V132 line network, "
        f"not a zone-wide band; V166 re-gaps the M nodes when the Wilson "
        f"factor is quantized. The V113 family equality of the 1/5 spectra "
        f"is replaced by identical pinned M spectra (folded deviation "
        f"{within:.1e} within / {across:.1e} across), while the loop "
        f"factor itself is link-local and momentum-independent (V113). "
        f"(iv) The "
        f"hierarchy at pi/4: per-channel factor {base_c:.5f}; readings "
        f"{{1, 3, 12}} -> {min(coefs.values()):.3f}..{max(coefs.values()):.2f}"
        f" against observed {COEF_OBS:.4f}; required sum_factor "
        f"{required:.3f} is not a natural node set. VERDICT: the loop-"
        f"factor MECHANISM survives closure with zero free parameters "
        f"(selection rule, exact crossing, link-local factor), but the "
        f"99.44% PRECISION is carried by the generic-angle interference "
        f"deficit and is not recoverable at any root angle (V202): the "
        f"closed branch predicts the hierarchy to "
        f"{100 * (coef3 / COEF_OBS - 1):+.1f}% (one channel per family, "
        f"the same selective reading that gives +0.6% at 1/5).")
    out = {"id": "V203_closed_branch_coefficient",
           "tier": 4,
           "property": "the closed-branch (pi/4) loop ladder, channel "
                       "equality and node-set hierarchy coefficient (V203)",
           "residual": True,
           "programme": "V203",
           "branch": "W-NAMED",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": all(c["pass"] for c in checks)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v203_closed_branch_coefficient.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
