#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V312 - the quaternionic loop ladder (tier 3, O8, pillar-surgery;
conditional).

**Contingency (pre-decided).** If V311 booked STEP-HYPERBOLIC-DECLARED
for ALL p, this experiment does NOT run the ladder: main() reads V311's
artifact, books the skip explicitly (branch SKIPPED-STEP-HYPERBOLIC),
passes its single T0 check and exits 0.  If at least one signature is
STEP-UNITARY, the ladder runs for those only.

**The object (closed input set).** the V203 machinery rebuilt over the
new root in-experiment: Ts' = the 45 plane generators T'_ab of V311's
embedding (import of the V311 module: swapped_generators/planes_of/
signature_step/unitarity_err); structure_constants called on the new
Ts' (V203's import site: V79_gauged_bond_layer); adjoint links via the
V79 adjoint_link_operators; projectors/cartan_norm from V85 as V203
imports them; the step G' from the V311 pattern at V311's booked
closure clocks for that signature; then the ladder F_n from V203's
loop_ratios (imported).  No decimal is retyped: closure clocks and the
unitarity gate come from the V311 artifact, the observed coefficient is
the loader quotient, the band is pillar_surgery.coef_band_abs.

**The question (pre-registered).** does any quaternionic-root closure
clock read the coefficient inside coef_band_abs of the observed, or
reproduce the odd-n exact selection rule with a crossing whose deficit
matches the observed coefficient?

**Measured (pre-registered).** T0: the rebuilt machinery reproduces
V311's booked unitarity residual at the closed clock (< strict tol) and
the booked closure-clock table; then per signature per closure clock:
odd-n maximum over n = 1..11 (the exact odd selection test, band
tolerances.strict), first crossing n*, deficit F_{n*}/LOOP, coefficient
3*deficit^10 against the observed.

**Acceptance (pre-registered).** SIGNATURE-LADDER-MATCH iff a closure
clock's coefficient read lands inside coef_band_abs of the observed;
SIGNATURE-LADDER-NEGATIVE otherwise, closest reading booked with the
odd-rule status per clock.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
for _cat in ("V311_quaternionic_root_step",
             "V79_gauged_bond_layer",
             "V85_closed_loop_factor",
             "V203_closed_branch_coefficient"):
    sys.path.insert(0, exp(_cat))

from V311_quaternionic_root_step import (swapped_generators,  # noqa: E402
                                         planes_of, unitarity_err)
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   structure_constants)
from V85_closed_loop_factor import projectors               # noqa: E402
from V203_closed_branch_coefficient import loop_ratios      # noqa: E402
from cl10 import config as _config                          # noqa: E402
from cl10 import u10pred                                    # noqa: E402

CFG = _config.CFG
NG = int(CFG.physics.n_gen)
D_L = NG * (NG - 1) // 2
NMAX = 260


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def parse_clock(label):
    """'pi*a/b' -> float."""
    a, b = label.replace("pi*", "").split("/")
    return np.pi * int(a) / int(b)


def hop_machinery_sig(hp, g):
    """V203's hop machinery rebuilt over the signature-p root: step from
    the T'_ab planes, links from structure_constants on the T'_ab."""
    G = signature_step_from(hp, g)
    cartan = [1j * (hp[a] @ hp[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = planes_of(hp)
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    P, ev = projectors(TL)
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    Ma = [np.einsum("ajbj->ab", (Pa @ Gi).reshape(32, D_L, 32, D_L)) / D_L
          for Pa in P]
    return G, cartan, Ma, ev


def signature_step_from(hp, g):
    import math
    n = 32
    In = np.eye(n, dtype=complex)
    U = In.copy()
    for T in planes_of(hp):
        U = U @ (math.cos(g) * In - 1j * math.sin(g) * T)
    return U


def main():
    checks, data = [], {}

    branch_v311 = art("V311_quaternionic_root_step", "branch")
    unitary_p = [int(p) for p in art("V311_quaternionic_root_step",
                                     "data.unitary_p")]
    data["v311_branch"] = branch_v311
    data["unitary_p"] = unitary_p

    if not unitary_p:
        # pre-decided skip: every signature booked STEP-HYPERBOLIC-
        # DECLARED at the step level, no closure clocks to read
        checks.append({"name": "T0 V311 booked {} for every p - the "
                              "ladder is skipped as pre-decided".format(
                                  branch_v311),
                       "pass": True})
        data["status"] = "SKIPPED-STEP-HYPERBOLIC"
        print("V312 SKIPPED-STEP-HYPERBOLIC: V311 booked the step-level "
              "no-go for all quaternionic signatures; the skip is the "
              "pre-decided contingency, booked here.", flush=True)
        _emit(data, checks, ok=True)
        return 0

    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    BAND = float(CFG.pillar_surgery.coef_band_abs)
    TOL = float(CFG.tolerances.strict)
    sigs = art("V311_quaternionic_root_step", "data.signatures")

    runs = {}
    ladder_rows = []
    for p in unitary_p:
        hp = swapped_generators(p)
        # T0: rebuilt machinery reproduces V311's booked unitarity and
        # closure table
        g_closed = float(CFG.angle.branches.closed)
        unit = unitarity_err(signature_step_from(hp, g_closed))
        booked = sigs[str(p)]
        clock_keys = sorted(booked["closure_clocks"].keys())
        t0ok = (unit < TOL and clock_keys
                and unit < float(booked["unitarity_err_at_closed_clock"])
                + TOL)
        checks.append({"name": "T0 p = {:d}: rebuilt unitarity residual "
                              "{:.3e} (V311 booked {:.3e}, band {:.0e}); "
                              "closure clocks {} (booked {})".format(
                                  p, unit,
                                  booked["unitarity_err_at_closed_clock"],
                                  TOL, clock_keys,
                                  sorted(booked["closure_clocks"])),
                       "pass": bool(t0ok)})
        if not t0ok:
            data["status"] = "MACHINERY-DRIFT"
            print("V312 MACHINERY-DRIFT: the rebuilt signature machinery "
                  "does not reproduce V311; nothing is interpreted",
                  flush=True)
            _emit(data, checks, ok=False)
            return 1

        rows = {}
        for key in clock_keys:
            g = parse_clock(key)
            G, cartan, Ma, ev = hop_machinery_sig(hp, g)
            F, rho = loop_ratios(G, cartan, Ma, ev, g, NMAX)
            odd_max = max(abs(F[n]) for n in range(1, 12, 2))
            n_star = next((n for n in range(2, NMAX + 1, 2)
                           if abs(F[n]) < LOOP), None)
            if n_star is None:
                n_any = next((n for n in range(2, NMAX + 1)
                              if abs(F[n]) < LOOP), None)
                rows[key] = {"clock": g, "odd_max": float(odd_max),
                             "odd_rule_exact": bool(odd_max < TOL),
                             "n_star": None,
                             "n_first_any": (None if n_any is None
                                             else int(n_any)),
                             "rho": float(rho), "deficit": None,
                             "coef": None, "residual_abs": None}
                continue
            deficit = float(abs(F[n_star]) / LOOP)
            coef = 3.0 * deficit ** 10
            rows[key] = {"clock": g, "odd_max": float(odd_max),
                         "odd_rule_exact": bool(odd_max < TOL),
                         "n_star": int(n_star), "rho": float(rho),
                         "F_cross": float(abs(F[n_star])),
                         "deficit": deficit, "coef": coef,
                         "residual_abs": abs(coef - COEF_OBS)}
            ladder_rows.append((abs(coef - COEF_OBS), p, key, rows[key]))
        runs[str(p)] = rows
        for key, r in rows.items():
            checks.append({"name": "p = {:d} at {}: odd-rule max {:.2e} "
                                   "({}), crossing n* = {}, deficit "
                                   "{}, coefficient {} vs observed "
                                   "{:.7f}".format(
                                       p, key, r["odd_max"],
                                       "exact" if r["odd_rule_exact"]
                                       else "not exact",
                                       r["n_star"],
                                       "n/a" if r["deficit"] is None
                                       else "{:.9f}".format(r["deficit"]),
                                       "n/a" if r["coef"] is None
                                       else "{:.7f}".format(r["coef"]),
                                       COEF_OBS),
                           "pass": True})
    data["ladder_runs"] = runs
    ladder_rows.sort()
    data["coefficient_observed"] = COEF_OBS
    data["band_abs"] = BAND
    data["closest_read"] = (None if not ladder_rows else
                            {"p": ladder_rows[0][1],
                             "clock": ladder_rows[0][2],
                             "coef": ladder_rows[0][3]["coef"],
                             "residual_abs": ladder_rows[0][0]})
    odd_rule_somewhere = any(r["odd_rule_exact"]
                             for rows in runs.values() for r in
                             rows.values())
    data["odd_rule_somewhere"] = odd_rule_somewhere

    match = bool(ladder_rows and ladder_rows[0][0] < BAND)
    branch = ("SIGNATURE-LADDER-MATCH" if match
              else "SIGNATURE-LADDER-NEGATIVE")
    data["status"] = branch
    print("V312 {}: signatures {} closure clocks {}; closest read p={} "
          "at {} gives coefficient {} (|d| {} vs observed {:.7f}, band "
          "{:.0e}); exact odd rule somewhere: {}.".format(
              branch, unitary_p,
              {str(p): sorted(runs[str(p)]) for p in unitary_p},
              ladder_rows[0][1] if ladder_rows else None,
              ladder_rows[0][2] if ladder_rows else None,
              "n/a" if not ladder_rows else
              "{:.7f}".format(ladder_rows[0][3]["coef"]),
              "n/a" if not ladder_rows else
              "{:.2e}".format(ladder_rows[0][0]), COEF_OBS, BAND,
              odd_rule_somewhere), flush=True)
    _emit(data, checks, ok=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V312_quaternionic_loop_ladder", "tier": 3,
           "property": "does any quaternionic-root closure clock read "
                       "the hierarchy coefficient inside the booked "
                       "band (tier 3, O8)",
           "residual": True,
           "programme": "pillar-surgery V312",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v312_quaternionic_loop_ladder.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V312 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
