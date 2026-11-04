#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V315 - the closure-quantizer census (tier 2, hard, O4/O8, theory-hunt).

**The question (pre-registered).** Is there a booked-class deformation of
the six-stage step whose closure class contains a RATIONAL clock pi*a/b
whose order-k closure coexists with the exact odd-n selection rule and
whose rung-10 read lands on the observed hierarchy coefficient inside
pillar_surgery.coef_band_abs?  1/5 = pi*(1/(5 pi)) is not rational in pi:
the search is for the closest rational closure AND the coefficient contact
at it, not for 1/5; |b - 10| is booked as the clock-order distance.

**The class (closed, fixed pre-run).** D-family: 2^6 = 64 per-stage twist
decorations on/off at the booked amplitude (alpha_10, alpha_10, alpha_10),
alpha_10 = angle.alpha_10, on V277's bloch_dec with the fixed GCLOSED
replaced by the scan clock g in both stage_matrix slots.  P-family: all
720 permutations of SC.STAGES (SC.bond_product's stages= parameter is the
engine's blessing).  Total 784.  BOOKED PRE-RUN: V277 books W(0) = I
(inherited), so every D-variant has the identity step at q = 0; and the
V308 traced_ladder body has NO stage-wise U16 link (its link layer is the
32-dim GUT hop and adjoint link operators), so the rung-10 ladder leg is
the plain traced ladder at the variant clock for EVERY variant - the table
is normalized: per-clock ladder reads + per-variant closure rows, the
join is exact by the class semantics, not a post-run reclassification.

**Measured (pre-registered).** Closure census in V311's pattern: every
g = pi*a/b, 0 < a < b <= closure_den_max (79 distinct), smallest
k <= closure_k_max with ||U(0)^k -/+ I||_max < CLOSE_EPS (V311's own 1e-12
definition line copied verbatim).  Per closure clock: ladder leg at that
clock (traced_ladder nmax=12), odd-rule check odd_max = max|F_n| over odd
n <= 11 vs tolerances.strict (V312's test), rung-10 coefficient
read_coef(F_10, LOOP) against COEF_OBS (V308's loader formula).
**T0 (pre-registered).** identity variant closes the clock set {pi/2,
pi/4, 3pi/4} (anchor: V311 p=0 artifact row keys verbatim; the 128-dim
fact M(0)^2 = I at pi/4 from V203 T2 sets the expected k = 2 - the k
difference vs V311's 32-dim k = 4 is booked, not hidden); the plain
rung-10 read at 1/5 reproduces V206's coef_at_fifth to < 1e-9.
Failure: MACHINERY-DRIFT, nothing interpreted.

**Acceptance (pre-registered).** QUANTIZER-CLOSURE-FOUND iff some variant
has a closure clock with odd-n exact AND |coef - COEF_OBS| <=
coef_band_abs; QUANTIZER-NEGATIVE otherwise, with the full census booked:
per-variant closure rows, per-clock ladder contacts, closest row named
with its clock-order distance |b - 10|, and the no-closure counts per
family.
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from cl10 import clifford as CF                            # noqa: E402
from cl10 import config as _config                         # noqa: E402
from cl10 import engine as E                               # noqa: E402
from cl10 import supercell as SC                           # noqa: E402
from cl10 import u10pred                                   # noqa: E402

CFG = _config.CFG
DEN_MAX = int(CFG.pillar_surgery.closure_den_max)
K_MAX = int(CFG.pillar_surgery.closure_k_max)
CLOSE_EPS = 1e-12                # V311's own definition line, verbatim
BAND = float(CFG.pillar_surgery.coef_band_abs)
TOL_STRICT = float(CFG.tolerances.strict)
G10 = float(CFG.angle.branches.implemented)   # 1/5, the rung-10 clock
ALPHA10 = float(CFG.angle.alpha_10)           # the booked decoration amp
GCLOSED = math.pi / 4.0                      # V277's sector-build clock
NMAX_LEG = 12                                 # the pre-registered leg

# ---------------------------------------------------------------------
# the shared construction (V277's twist/VT block, copied verbatim and
# re-parameterized: the stage_matrix clocks become the scan clock g)
# ---------------------------------------------------------------------
_B16, _U32, U16C = SC.sectors(GCLOSED)
_H = [E.asnumpy(m).astype(complex) for m in CF.generators()]
VT = [_B16.conj().T @ (1j * (_H[a] @ _H[b])) @ _B16
      for (a, b) in ((0, 1), (1, 2), (0, 2))]


def twist(q, amp):
    """W(q) = exp(i q.(a V1 + b V2 + c V3)), exact via the 16x16 eigh
    (V277 lines 94-98 copied verbatim)."""
    a, b, c = amp
    qV = a * q[0] * VT[0] + b * q[1] * VT[1] + c * q[2] * VT[2]
    if a == 0.0 and b == 0.0 and c == 0.0:
        return np.eye(16, dtype=complex)
    ev, Vm = np.linalg.eigh(qV)
    return (Vm * np.exp(1j * ev)) @ Vm.conj().T


def dec_step(q, g, amps, stages):
    """The decorated six-stage 128x128 step at the SCAN clock g: V277's
    bloch_dec with GCLOSED -> g in both stage_matrix slots (the class)."""
    M = np.eye(128, dtype=complex)
    for t, (even, j) in enumerate(stages):
        S = SC.stage_matrix(even, j, q[j], g, g)
        M = np.kron(S, U16C @ twist(q, amps[t])) @ M
    return M


def traced_ladder(g, nmax=NMAX_LEG):
    """The V203 traced-link ladder at the clock g: F_n, rho (V308's
    traced_ladder body, the plain leg the class pre-books for every
    variant)."""
    G, cartan, Ma, ev = hop_machinery(g)
    F, rho = loop_ratios(G, cartan, Ma, ev, g, nmax)
    return F, rho


def read_coef(Fv, LOOP):
    """The published channel reading 3 (16 pi F)^10 = 3 (F/LOOP)^10
    (V308's read_coef, copied)."""
    return 3.0 * (Fv / LOOP) ** 10


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def closure_row(U):
    """Smallest k <= K_MAX with ||U^k -/+ I||_max < CLOSE_EPS, or None."""
    Uk = np.eye(128, dtype=complex)
    I128 = np.eye(128, dtype=complex)
    for k in range(1, K_MAX + 1):
        Uk = Uk @ U
        dev = min(float(np.abs(Uk - I128).max()),
                  float(np.abs(Uk + I128).max()))
        if dev < CLOSE_EPS:
            return {"k": k, "dev": dev}
    return None


def main():
    t0 = time.time()
    checks, data = [], {}
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)

    # ---- clock scan (V311's pattern: 79 distinct pi*a/b clocks) ------
    clocks = []
    seen = set()
    for b in range(2, DEN_MAX + 1):
        for a in range(1, b):
            g = np.pi * a / b
            key = round(g, 15)
            if key in seen:
                continue
            seen.add(key)
            clocks.append({"key": "pi*%d/%d" % (a, b), "g": float(g),
                           "a": a, "b": b})
    data["clocks_tested"] = len(clocks)

    # ---- T0a: identity variant closes the anchored clock set ----------
    anchor = art("V311_quaternionic_root_step",
                 "data.signatures.0.closure_clocks")
    anchor_keys = sorted(anchor.keys())
    t0a_rows = {}
    t0a_ok = True
    for ck in anchor_keys:
        g = math.pi * int(ck.split("*")[1].split("/")[0]) \
            / int(ck.split("/")[1])
        row = closure_row(dec_step((0.0, 0.0, 0.0), g,
                                   [(0.0, 0.0, 0.0)] * 6, SC.STAGES))
        t0a_rows[ck] = {"rederived_k": (None if row is None else
                                        row["k"]),
                        "rederived_dev": (None if row is None
                                          else row["dev"]),
                        "v311_32dim_k": anchor[ck]["k"],
                        "v311_32dim_dev": anchor[ck]["dev"]}
        if row is None:
            t0a_ok = False
    # the 128-dim booked fact at pi/4 is M(0)^2 = I (V203 T2): k = 2
    if t0a_rows.get("pi*1/4", {}).get("rederived_k") != 2:
        t0a_ok = False
    data["t0a_identity_closure"] = t0a_rows
    checks.append({"name": "T0a identity decoration + default order "
                          "closes the anchored clock set {} (V311 p=0 "
                          "artifact keys verbatim; 128-dim k at pi/4 = 2 "
                          "per V203's booked M(0)^2 = I - the k "
                          "difference vs V311's 32-dim signature row "
                          "(pi*1/4: k=4) is booked, all devs < {:.0e})"
                   .format(anchor_keys, CLOSE_EPS),
                   "pass": bool(t0a_ok)})
    if not t0a_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V315 MACHINERY-DRIFT: the identity six-stage step no "
              "longer closes the anchored clocks; nothing interpreted",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0b: the plain rung-10 read reproduces V206 coef_at_fifth ---
    coef_fifth_art = float(art("V206_clock_tension",
                               "data.data_clock.coef_at_fifth"))
    F5, _rho5 = traced_ladder(G10)
    coef_fifth = read_coef(F5[10], LOOP)
    t0b = abs(coef_fifth - coef_fifth_art)
    data["t0b_coef_fifth"] = {"read": coef_fifth,
                              "artifact": coef_fifth_art,
                              "abs_diff": t0b}
    checks.append({"name": "T0b plain rung-10 read at g = 1/5 reproduces "
                          "V206's booked coef_at_fifth {:.15f} (|d| "
                          "{:.1e} < 1e-9)".format(coef_fifth_art, t0b),
                   "pass": bool(t0b < 1e-9)})
    if t0b >= 1e-9:
        data["status"] = "MACHINERY-DRIFT"
        print("V315 MACHINERY-DRIFT: the rung-10 machinery no longer "
              "reproduces V206's booked clock read; nothing interpreted",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the per-clock ladder leg (the class books it variant-plain:
    # W(0) = I inherited on D; the traced-link ladder carries no
    # stage-wise U16 link - booked in the pre-registration) ------------
    print("[V315] per-clock ladder legs ({} clocks, nmax = {})"
          .format(len(clocks), NMAX_LEG), flush=True)
    clock_leg = {}
    for c in clocks:
        g = c["g"]
        F, rho = traced_ladder(g)
        odd_max = max(abs(F[n]) for n in range(1, NMAX_LEG, 2))
        coef = read_coef(F[10], LOOP)
        clock_leg[c["key"]] = {"g": g, "a": c["a"], "b": c["b"],
                               "F_10": float(F[10]),
                               "odd_max": float(odd_max),
                               "odd_rule_exact": bool(odd_max
                                                      < TOL_STRICT),
                               "coef": coef,
                               "residual": abs(coef - COEF_OBS),
                               "rho": float(rho)}
    data["clock_leg"] = clock_leg
    data["coefficient_observed"] = COEF_OBS
    data["band_coef"] = BAND

    # ---- the 784-variant census --------------------------------------
    print("[V315] variant census: 64 D-decorations + 720 P-orderings",
          flush=True)
    amp_on = (ALPHA10, ALPHA10, ALPHA10)
    amp_off = (0.0, 0.0, 0.0)
    variants = []
    n_found_rows = 0
    best_rows = []
    d_no_closure = 0
    p_no_closure = 0
    for fam, amps, stages in (
            [("D", list(itertools.product((amp_off, amp_on), repeat=6)),
              [SC.STAGES] * 64)]
            + [("P", [(amp_off,) * 6] * 720,
                list(itertools.permutations(SC.STAGES)))]):
        for amps6, st in zip(amps, stages):
            closures = {}
            for c in clocks:
                U = dec_step((0.0, 0.0, 0.0), c["g"], amps6, st)
                row = closure_row(U)
                if row is not None:
                    closures[c["key"]] = row
            if not closures:
                if fam == "D":
                    d_no_closure += 1
                else:
                    p_no_closure += 1
            variants.append({
                "family": fam,
                "amps": ["on" if a is amp_on else "off" for a in amps6],
                "stage_order": ["T{}".format(j) if e else "F{}"
                                .format(j) for (e, j) in st],
                "n_closures": len(closures),
                "closures": {k: {"k": v["k"], "dev": v["dev"]}
                             for k, v in closures.items()}})
            for k, v in closures.items():
                leg = clock_leg[k]
                if leg["odd_rule_exact"] and leg["residual"] <= BAND:
                    n_found_rows += 1
                if leg["odd_rule_exact"]:
                    best_rows.append({"family": fam, "clock": k,
                                      "k": v["k"],
                                      "order_distance_b_minus_10":
                                          abs(leg["b"] - 10),
                                      "coef": leg["coef"],
                                      "residual": leg["residual"]})
    best_rows.sort(key=lambda r: (r["residual"],
                                  r["order_distance_b_minus_10"]))
    data["variants"] = variants
    data["no_closure_counts"] = {"D": d_no_closure, "P": p_no_closure}
    data["odd_exact_closure_rows"] = best_rows
    data["closest_odd_closure_row"] = (best_rows[0] if best_rows
                                       else None)

    # ---- branch -------------------------------------------------------
    branch = ("QUANTIZER-CLOSURE-FOUND" if n_found_rows
              else "QUANTIZER-NEGATIVE")
    data["status"] = branch
    checks.append({"name": "branch {}: {} variants censused (64 D + 720 "
                          "P), clocks with >= 1 closure variant: {}, "
                          "odd-exact closure rows: {}, closest odd-exact "
                          "closure {} (|coef - COEF_OBS| {:.4e}, band "
                          "{:.0e}, |b-10| {}), no-closure counts D {} P "
                          "{}; table normalized per the class (per-clock "
                          "plain leg + per-variant closure rows - the "
                          "D-family is identity-degenerate at q=0 by "
                          "V277's booked W(0)=I)".format(
                              branch, len(variants),
                              sum(1 for k in clock_leg
                                  if any(k in v["closures"]
                                         for v in variants)),
                              len(best_rows),
                              "none" if not best_rows else
                              best_rows[0]["clock"],
                              best_rows[0]["residual"] if best_rows
                              else float("nan"), BAND,
                              best_rows[0]["order_distance_b_minus_10"]
                              if best_rows else "n/a",
                              d_no_closure, p_no_closure),
                   "pass": True})
    print("V315 {}: {} variants; odd-exact rational closures {} at the "
          "coefficient band (closest odd-exact closure {}, |d| {:.3e}, "
          "|b-10| {}); no-closure D {} P {}.".format(
              branch, len(variants), "FOUND" if n_found_rows else "miss "
              "the band", "none" if not best_rows else best_rows[0]
              ["clock"], best_rows[0]["residual"] if best_rows else -1.0,
              best_rows[0]["order_distance_b_minus_10"] if best_rows
              else "n/a", d_no_closure, p_no_closure), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V315_closure_quantizer_census", "tier": 2,
           "property": "closure-quantizer census of the six-stage step: "
                       "64 twist decorations x 720 stage orderings, "
                       "rational closure clocks, odd-rule coexistence and "
                       "the rung-10 coefficient contact (tier 2, hard, "
                       "O4/O8)",
           "residual": True, "programme": "theory-hunt V315",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v315_closure_quantizer_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V315 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
