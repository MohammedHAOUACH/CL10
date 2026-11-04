#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V308 - coefficient variants at the fifth clock (tier 3, O8, pillar-surgery).

**The question (pre-registered).** the rung-10 clock 1/5 reads the
coefficient at 2.088106865906787 (V206 `coef_at_fifth`) vs the observed
2.0765641 - a gap of -0.55%.  Does any *variant reading of the loop
machinery at that clock* close it?

**The object (closed input set).** V203's hop_machinery(g10) and
loop_ratios imported verbatim (the same import path V204/V276 use);
V206's booked data_clock (coef_at_fifth); V218's return_amplitude ladder
function reused on the icosahedral module step.  The observed
coefficient is the loader quotient (V_EW/M_P)/LOOP^10, never a typed
decimal.  Band: pillar_surgery.coef_band_abs from config.

Variants, all at the rung-10 clock g = angle.branches.implemented:
  (a) plain rung-10 read: 3 (F_10/LOOP)^10 from the traced-link ladder;
  (b) certified-continuation read: the V276 T1 geometric continuation
      recipe - (rho, F_cross) re-derived at g inside this experiment,
      F_10_geom = F_cross * rho**(-(n_star - 10)), read 3 (F_10_geom/
      LOOP)^10;
  (c) icosahedral-module closed-word ladder: F_10 = < tr M(q)^10 > box
      mean of V218 on both booked generator triples, read
      3 (F_10/LOOP)^10, with V218's grid-resolution test at n = 10
      recorded next to it.

**Measured (pre-registered).** T0 gate: the plain rung-10 read must
reproduce V206's booked coef_at_fifth to 1e-9.  Each variant's
coefficient, its residual against the observed, and (for c) whether the
F_10 is grid-resolved.  If a variant lands inside coef_band_abs, it is
re-evaluated by a second independent call and required to differ from
the plain rung-10 read.

**Acceptance (pre-registered).** COEF-VARIANT-CLOSED iff some variant
is within coef_band_abs of the observed (supplied variant booked, twice
reproduced); COEF-GAP-PERSISTS otherwise, with the five smallest variant
residuals and their provenance booked.
"""
import json
import math
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))
sys.path.insert(0, exp("V218_closed_word_census_module"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from V218_closed_word_census_module import return_amplitude  # noqa: E402
from cl10 import config as _config  # noqa: E402
from cl10 import continuum as K  # noqa: E402
from cl10 import starcell as X  # noqa: E402
from cl10 import supercell as SC  # noqa: E402
from cl10 import u10pred  # noqa: E402

CFG = _config.CFG
G5 = float(CFG.angle.branches.implemented)      # the rung-10 clock 1/5
NMAX = 260                                      # V276's continuation range
GRIDS = (12, 16, 24, 32)                        # V218's grids, unchanged
BOX = 4.0 * math.pi                             # V218's box half-width
TRIPLES = ((0, 1, 2), (0, 2, 4))                # V218's booked triples
RES_CUT = 0.05                                  # V218's resolution cut


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def traced_ladder(g, nmax=NMAX):
    """The V203 traced-link ladder at the clock g: F_n, rho."""
    G, cartan, Ma, ev = hop_machinery(g)
    F, rho = loop_ratios(G, cartan, Ma, ev, g, nmax)
    return F, rho


def star_ladders(g):
    """V218's icosahedral return-amplitude ladder at the clock g, both
    booked triples, all booked grids: {triple: {grid: [F_0..F_nmax]}}."""
    axes6 = [np.asarray(v, float) for v in K.module_generators()["vectors"]]
    _B, _U32, U16 = SC.sectors()
    cell6 = X.module_cell(shape=(2, 2, 2))
    stages6 = X.module_stages(shape=(2, 2, 2), order="block-major")
    out = {}
    for tri in TRIPLES:
        ax = [axes6[i] for i in tri]

        def star_step(q, ax=ax):
            return X.bloch_star(cell6, U16, g, stages6,
                                np.asarray(q, float), ax)
        out[tri] = {}
        for grid in GRIDS:
            F, _npts = return_amplitude(star_step, grid, 10, box=BOX)
            out[tri][grid] = [float(z.real) for z in F]
    return out


def star_resolution(F_by_grid, n=10):
    """V218's T2b resolution at length n over the grids fine enough."""
    gs = [g for g in GRIDS if g > n + 1]
    vals = [F_by_grid[g][n] for g in gs]
    scale_n = max(max(abs(v) for v in vals), 1e-30)
    return (max(abs(a - b) for a in vals for b in vals) / scale_n, gs)


def read_coef(Fv, LOOP):
    """The published channel reading 3 (16 pi F)^10 = 3 (F/LOOP)^10."""
    return 3.0 * (Fv / LOOP) ** 10


def main():
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    BAND = float(CFG.pillar_surgery.coef_band_abs)
    checks, data = [], {}

    coef_fifth_art = float(art("V206_clock_tension",
                               "data.data_clock.coef_at_fifth"))

    # ---- T0 gate: the plain rung-10 read reproduces V206 --------------
    F, rho = traced_ladder(G5)
    coef_a = read_coef(F[10], LOOP)
    t0_def = abs(coef_a - coef_fifth_art)
    data["t0"] = {"g": G5, "F_10": float(F[10]), "rho": float(rho),
                  "coef_plain_read": coef_a,
                  "coef_at_fifth_artifact": coef_fifth_art,
                  "abs_diff": t0_def}
    checks.append({"name": "T0 plain rung-10 read at g = 1/5 reproduces "
                          "V206's booked coef_at_fifth {:.15f} (|d| "
                          "{:.1e} < 1e-9)".format(coef_fifth_art, t0_def),
                   "pass": bool(t0_def < 1e-9)})
    if t0_def >= 1e-9:
        data["status"] = "MACHINERY-DRIFT"
        print("V308 MACHINERY-DRIFT: the rung-10 machinery no longer "
              "reproduces V206's booked clock read; nothing downstream "
              "is interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- variant (a): plain rung-10 read ------------------------------
    # ---- variant (b): certified continuation at g = 1/5 (V276 T1) -----
    crossing = next((n for n in range(2, NMAX + 1, 2) if F[n] < LOOP),
                    None)
    degenerate = (crossing is None or crossing == 10)
    F_cross = float(F[crossing]) if crossing is not None else None
    F_geom = (F_cross * rho ** (-(crossing - 10.0))
              if crossing is not None else None)
    coef_b = read_coef(F_geom, LOOP) if F_geom is not None else None
    data["continuation"] = {
        "crossing_n": (None if crossing is None else int(crossing)),
        "F_cross": F_cross, "rho": float(rho),
        "F_10_geom": F_geom,
        "degenerate_at_this_clock": bool(degenerate),
        "direction": ("F_n decays with n (rho = {:.16f} < 1): the V276 "
                      "continuation from (rho, F_cross) at n* = {} back "
                      "to n = 10 multiplies F_cross by rho**(-({}-10)) = "
                      "{:.6f}".format(rho, crossing, crossing,
                                      rho ** (-(crossing - 10.0))
                                      if crossing is not None else
                                      float("nan")))}

    # ---- variant (c): icosahedral closed-word ladder ------------------
    F_star = star_ladders(G5)
    star_rows = {}
    for tri in TRIPLES:
        res, gs = star_resolution(F_star[tri])
        F10 = F_star[tri][GRIDS[-1]][10]
        star_rows[str(tri)] = {
            "F_10_finest_grid": F10,
            "coef": read_coef(F10, LOOP),
            "resolution_at_n10": res, "grids": list(gs),
            "resolved": bool(res < RES_CUT)}
        data["continuation"]["star_ladder_" + str(tri)] = {
            "F_by_grid": {str(g): F_star[tri][g] for g in GRIDS}}

    variants = {
        "a_plain_rung10": {"coef": coef_a,
                           "provenance": "V203 traced-link F_10 at g=1/5, "
                                         "read 3(F/LOOP)^10 (the booked "
                                         "V206 fifth-clock read)"},
        "b_certified_continuation": {
            "coef": coef_b,
            "provenance": "V276 T1 geometric continuation "
                          "F_10_geom = F_cross * rho**(-(n*-10)), "
                          "(rho, F_cross) re-derived at g=1/5 inside "
                          "this experiment{}".format(
                              " - DEGENERATE: n* = 10 at this clock, the "
                              "continuation is the identity"
                              if degenerate else "")},
    }
    for tri in TRIPLES:
        variants["c_icosahedral_" + str(tri)] = {
            "coef": star_rows[str(tri)]["coef"],
            "provenance": "V218 return-amplitude F_10 box mean on triple "
                          "{} at g=1/5, grid {}, resolution {:.2e} "
                          "({})".format(tri, GRIDS[-1],
                                        star_rows[str(tri)]
                                        ["resolution_at_n10"],
                                        "resolved"
                                        if star_rows[str(tri)]["resolved"]
                                        else "UNRESOLVED - see booking")}
    data["variants"] = {}
    scored = []
    for name, v in variants.items():
        res_abs = None if v["coef"] is None else abs(v["coef"] - COEF_OBS)
        data["variants"][name] = {
            "coef": v["coef"], "provenance": v["provenance"],
            "residual_abs": res_abs,
            "residual_rel": (None if v["coef"] is None else
                             v["coef"] / COEF_OBS - 1.0)}
        if res_abs is not None:
            scored.append((res_abs, name))
        checks.append({"name": "variant {}: coefficient {} vs observed "
                              "{:.7f} (|d| {})".format(
                                  name,
                                  "n/a" if v["coef"] is None
                                  else "{:.7f}".format(v["coef"]),
                                  COEF_OBS,
                                  "n/a" if res_abs is None
                                  else "{:.2e}".format(res_abs)),
                       "pass": True})
    scored.sort()
    data["five_smallest_residuals"] = [
        {"name": n, "residual_abs": r,
         "provenance": variants[n]["provenance"]}
        for r, n in scored[:5]]
    data["coefficient_observed"] = COEF_OBS
    data["band_abs"] = BAND

    # ---- branch ----
    hits = [(r, n) for r, n in scored if r < BAND]
    supplied = None
    if hits:
        # the supplied variant owes a second reproduction and must differ
        # from the plain rung-10 read
        r_hit, n_hit = hits[0]
        if n_hit == "a_plain_rung10":
            F2, _ = traced_ladder(G5)
            ok2 = abs(read_coef(F2[10], LOOP) - coef_a) < 1e-12
            checks.append({"name": "supplied variant (a) second call "
                                  "reproduces: {}".format(ok2),
                           "pass": bool(ok2)})
            supplied = n_hit if ok2 else None
        elif n_hit.startswith("b_"):
            F2, rho2 = traced_ladder(G5)
            cr2 = next((n for n in range(2, NMAX + 1, 2)
                        if F2[n] < LOOP), None)
            coef2 = (read_coef(F2[cr2] * rho2 ** (-(cr2 - 10.0)), LOOP)
                     if cr2 is not None else None)
            ok2 = coef2 is not None and abs(coef2 - coef_b) < 1e-12
            checks.append({"name": "supplied variant (b) second call "
                                  "reproduces: {}".format(ok2),
                           "pass": bool(ok2)})
            supplied = n_hit if ok2 else None
        else:
            tri = tuple(int(x) for x in n_hit.split("_")[-1]
                        .strip("()").split(","))
            axes6 = [np.asarray(v, float)
                     for v in K.module_generators()["vectors"]]
            _B, _U32, U16 = SC.sectors()
            cell6 = X.module_cell(shape=(2, 2, 2))
            stages6 = X.module_stages(shape=(2, 2, 2), order="block-major")
            ax = [axes6[i] for i in tri]

            def star_step2(q, ax=ax):
                return X.bloch_star(cell6, U16, G5, stages6,
                                    np.asarray(q, float), ax)
            F2, _ = return_amplitude(star_step2, GRIDS[-1], 10, box=BOX)
            coef2 = read_coef(float(F2[10].real), LOOP)
            ref = star_rows[str(tri)]["coef"]
            ok2 = abs(coef2 - ref) < 1e-9
            differs = abs(ref - coef_a) > 1e-6
            checks.append({"name": "supplied variant (c) second call "
                                  "reproduces: {}, differs from plain "
                                  "rung-10 read: {}".format(ok2, differs),
                           "pass": bool(ok2 and differs)})
            supplied = n_hit if (ok2 and differs) else None

    branch = ("COEF-VARIANT-CLOSED" if supplied else "COEF-GAP-PERSISTS")
    data["status"] = branch
    data["supplied_variant"] = supplied
    gap = coef_a / COEF_OBS - 1.0
    print("V308 {}: the fifth-clock plain read sits at {:.9f}, gap {:+.4%} "
          "against the observed {:.7f}; the certified continuation at this "
          "clock is {} and reads {}; the icosahedral ladder reads {}; five "
          "smallest residuals {}.".format(
              branch, coef_a, gap, COEF_OBS,
              "degenerate (n* = 10)" if degenerate
              else "non-degenerate (n* = {})".format(crossing),
              "n/a" if coef_b is None else "{:.7f}".format(coef_b),
              {t: "{:.7f}".format(star_rows[str(t)]["coef"])
               for t in TRIPLES},
              [(n, "{:.2e}".format(r)) for r, n in scored[:5]]),
          flush=True)
    _emit(data, checks, ok=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V308_coefficient_variants_fifth_clock", "tier": 3,
           "property": "does a variant reading of the loop machinery at "
                       "the rung-10 clock close the -0.55% coefficient "
                       "gap (tier 3, O8)",
           "residual": True,
           "programme": "pillar-surgery V308",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v308_coefficient_variants_fifth_clock.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V308 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
