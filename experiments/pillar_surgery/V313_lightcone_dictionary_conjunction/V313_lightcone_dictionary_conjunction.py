#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V313 - the lightcone dictionary conjunction (tier 3, F8/O10,
pillar-surgery).

**The question (pre-registered).** the exactly lightlike graviton cone
exists but only in a 53-parameter family nothing selects, whose known
repair clock g = 1/2 collapses the loop ladder (paper Section 11.6,
V208/V285/V229).  Re-parameterize the family and test the conjunction
member-wise.

**The object (closed input set).** V285's machinery reused wholesale by
import (build_Oh/module_reps/velocity_ops, the covariant 160-space
basis, the Clifford-regularity Jacobian with rank 107 / nullity 53,
the tangent coordinates delta in R^53, the Gauss-Newton retraction and
the Nelder-Mead locality functional T with pillar_surgery.nelder_starts
starts, seed 285 as V285 books).  The ladder leg is V203's
hop_machinery/loop_ratios run at every scan clock.  NOTE booked before
the run: the Jacobian, the family and T have NO clock parameter in
their construction - they are properties of the V208 committed
generator - so the family re-measured at every scan_clocks entry is
the same object; the scan's content is clock-blindness of the
dictionary side plus the ladder intactness of each clock.

**Measured (pre-registered).**
1. T0 gate (repo convention: it outranks interpretation): the
   Jacobian rank/nullity at pi/4 must reproduce V285's booked 107/53,
   with V208's profile/hop/dressing anchors (the V285 T0 pattern); on
   failure the verdict is MACHINERY-DRIFT and the scan does not run.
2. the family dimension re-measured at every clock in
   pillar_surgery.scan_clocks (0.125, 0.2, 0.5, pi/4, pi/2).
3. ladder intactness per clock: the V203 first crossing n* and deficit
   F_{n*}/LOOP at that clock, the band
   |deficit(g)/deficit(pi/4) - 1| <= ladder_intact_band against the
   V206-booked pi/4 deficit; the pi/4 deficit itself is re-derived and
   anchored to the booked value (a vanishing odd length is not a
   crossing where the exact odd rule holds - the V276/V204 scan
   convention).
4. the T minimization over the 53 coordinates (nelder_starts starts,
   seed 285, maxfev 200, box as V285): T_committed, T*, interior
   diagnostics (53-gradient norm, box-face coordinates), distinct
   minima count, hop/dressing projections and i3-distance of the
   minimizer.
5. the bond-entropy weight G = 1/(4 ln 45) reported as a weight only
   (never refit): G*T* and G*T_committed.

**Acceptance (pre-registered).** LIGHTCONE-SELECTED iff a member at a
ladder-intact clock uniquely minimizes a booked selector functional
(T* < T_committed strictly, interior, one distinct minimum);
LIGHTCONE-CLOCK-CONFLICT otherwise - families empty at every
ladder-intact clock, or the selector never unique: the bounded no-go
with the per-clock family dimensions booked.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V285_most_local_lightlike_generator"))
sys.path.insert(0, exp("V203_closed_branch_coefficient"))

import numpy as np                                    # noqa: E402

import V285_most_local_lightlike_generator as V285    # noqa: E402
from V203_closed_branch_coefficient import (hop_machinery,   # noqa: E402
                                            loop_ratios)
from cl10 import config as _config                     # noqa: E402
from cl10 import engine as E                           # noqa: E402
from cl10 import parpool                               # noqa: E402
from cl10 import u10pred                               # noqa: E402

CFG = _config.CFG
NSTARTS = int(CFG.pillar_surgery.nelder_starts)
MAXFEV = 200                                          # V285's booked cap
BAND_INTACT = float(CFG.pillar_surgery.ladder_intact_band)
SCAN = [float(g) for g in CFG.pillar_surgery.scan_clocks]


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def ladder_at(g, LOOP, tol_strict, nmax=260):
    """V203's first crossing of LOOP at clock g and its deficit.  The
    scan skips the identically-vanishing odd lengths on clocks where the
    exact odd selection rule holds (V276/V204 pattern) - on those
    clocks a zero is not a crossing."""
    G, cartan, Ma, ev = hop_machinery(g)
    F, rho = loop_ratios(G, cartan, Ma, ev, g, nmax)
    odd_max = max(abs(F[n]) for n in range(1, 12, 2))
    odd_rule = odd_max < tol_strict
    rng = range(2, nmax + 1, 2) if odd_rule else range(2, nmax + 1)
    n_star = next((n for n in rng if abs(F[n]) < LOOP), None)
    deficit = float(abs(F[n_star]) / LOOP) if n_star else None
    return {"clock": g, "n_star": (None if n_star is None else
                                   int(n_star)),
            "deficit": deficit, "rho": float(rho),
            "odd_rule_exact": bool(odd_rule),
            "odd_max": float(odd_max)}


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V313 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    # ---- V285 machinery rebuilt in-experiment ----------------------------
    Rv, Dm = V285.module_reps()
    Bop, Kop = V285.velocity_ops()
    gen = json.load(open(os.path.join(
        exp("V208_metric_bilinear_cone"), "cone_generator.json")))
    H = [np.array(gen["A"][k]) + 1j * np.array(gen["S"][k])
         for k in range(3)]
    PROP = [i for i, R in enumerate(Rv) if np.linalg.det(R) > 0]

    def cov_project(Ts):
        outp = [np.zeros((36, 36), dtype=complex) for _ in range(3)]
        for i in PROP:
            D, R = Dm[i], Rv[i]
            DT = [D @ Tb @ D.conj().T for Tb in Ts]
            for a in range(3):
                for b in range(3):
                    outp[a] += R[a, b] * DT[b]
        return [o / 24. for o in outp]

    srng = np.random.default_rng(20261021)      # V285's basis seed

    def rherm():
        A = (srng.standard_normal((36, 36))
             + 1j * srng.standard_normal((36, 36)))
        return (A + A.conj().T) / 2.0
    rows, rank = [], 0
    for _ in range(4 * 160 + 60):
        v = V285.vec3(cov_project([rherm(), rherm(), rherm()]))
        r = int(np.linalg.matrix_rank(np.array(rows + [v]), tol=1e-9))
        if r > rank:
            rows.append(v)
            rank = r
        if rank == 160:
            break
    _, sv, Vt = np.linalg.svd(np.array(rows), full_matrices=False)
    Bmat = Vt[sv > 1e-9]

    def linearised(dH):
        d1, d2, d3 = dH
        return [V285.sym(H[0], d2) + V285.sym(d1, H[1]),
                V285.sym(H[0], d3) + V285.sym(d1, H[2]),
                V285.sym(H[1], d3) + V285.sym(d2, H[2]),
                V285.sym(H[0], d1) - V285.sym(H[1], d2),
                V285.sym(H[0], d1) - V285.sym(H[2], d3)]
    cols = []
    for v in Bmat:
        ms = linearised(V285.unvec(v))
        cols.append(np.concatenate([m.real.reshape(-1) for m in ms]
                                   + [m.imag.reshape(-1) for m in ms]))
    Jm = np.array(cols).T
    jac_rank = int(np.linalg.matrix_rank(Jm, tol=1e-9))
    family_dim = 160 - jac_rank

    # ---- T0: V208/V285 anchors, pi/4 gate --------------------------------
    profH = V285.profile(H)
    a_prof = art("V208_metric_bilinear_cone", "data.cone.kernel_profile")
    pdef = max(abs(profH[r] - float(a_prof[str(r)])) / float(a_prof[str(r)])
               for r in range(5))
    fB = (V285.i3(H, Bop) / V285.i3(Bop, Bop)) ** 2 * V285.i3(Bop, Bop) \
        / V285.i3(H, H)
    fK = (V285.i3(H, Kop) / V285.i3(Kop, Kop)) ** 2 * V285.i3(Kop, Kop) \
        / V285.i3(H, H)
    a_fB = float(art("V208_metric_bilinear_cone",
                     "data.implemented."
                     "hop_fraction_of_cone_generator"))
    a_fK = float(art("V208_metric_bilinear_cone",
                     "data.implemented."
                     "dressing_fraction_of_cone_generator"))
    a_b1 = int(art("V208_metric_bilinear_cone",
                   "data.cone.cone_branches_at_plus_1"))
    a_bm1 = int(art("V208_metric_bilinear_cone",
                    "data.cone.cone_branches_at_minus_1"))
    wH = np.linalg.eigvalsh((H[0] + H[0].conj().T) / 2)
    b1 = int(np.sum(np.abs(wH - 1.0) < 1e-6))
    bm1 = int(np.sum(np.abs(wH + 1.0) < 1e-6))
    T_committed = sum(profH[r] ** 2 for r in range(2, 8))
    a_Tcommit = float(art("V285_most_local_lightlike_generator",
                          "data.t0.T_committed"))
    gate = (jac_rank == 107 and family_dim == 53 and pdef < 1e-9
            and abs(fB - a_fB) < 1e-9 * abs(a_fB)
            and abs(fK - a_fK) < 1e-9 * abs(a_fK)
            and b1 == a_b1 and bm1 == a_bm1
            and abs(T_committed - a_Tcommit) < 1e-12)
    data["t0"] = {"rank": jac_rank, "family_dim": family_dim,
                  "profile_defect": pdef, "fB": fB, "fK": fK,
                  "fB_art": a_fB, "fK_art": a_fK,
                  "branches": [b1, bm1], "branches_art": [a_b1, a_bm1],
                  "T_committed": T_committed,
                  "T_committed_artifact": a_Tcommit}
    checks.append({"name": "T0 pi/4 gate: Jacobian rank {} / nullity {} "
                          "must be 107/53; V208 anchors profile {:.1e}, "
                          "hop {:.10f}/{:.10f}, dressing {:.10f}/{:.10f}, "
                          "T_committed {:.6e} (artifact {:.6e})".format(
                              jac_rank, family_dim, pdef, fB, a_fB,
                              fK, a_fK, T_committed, a_Tcommit),
                   "pass": bool(gate)})
    if not gate:
        data["status"] = "MACHINERY-DRIFT"
        print("V313 MACHINERY-DRIFT: the pi/4 Jacobian/anchors do not "
              "reproduce V285's booked 107/53; the scan does not run",
              flush=True)
        _emit(data, checks, ok=False)
        return 1
    log("T0 gate passed")

    _, svJ, Vh = np.linalg.svd(Jm, full_matrices=False)
    NULLc = Vh[jac_rank:]

    # ---- family re-measured at every scan clock ---------------------------
    per_clock = {}
    for g in SCAN:
        # the Jacobian has no clock parameter (booked note): the same
        # re-measurement applies at every clock label
        rk = int(np.linalg.matrix_rank(Jm, tol=1e-9))
        per_clock[repr(g)] = {"rank": rk, "family_dim": 160 - rk}
    data["family_per_clock"] = per_clock
    checks.append({"name": "family dimension at the scan clocks {}: {} "
                           "(the construction is clock-blind - booked "
                           "in the pre-reg)".format(
                               SCAN, {k: v["family_dim"]
                                      for k, v in per_clock.items()}),
                   "pass": True})

    # ---- ladder intactness per scan clock ---------------------------------
    LOOP = float(u10pred.LOOP)
    def_pi4 = float(art("V206_clock_tension",
                        "data.closure_class.pi/4.deficit"))
    tol_strict = float(CFG.tolerances.strict)
    ladders = {}
    for g in SCAN:
        row = ladder_at(g, LOOP, tol_strict)
        row["intact"] = (row["deficit"] is not None
                         and abs(row["deficit"] / def_pi4 - 1.0)
                         <= BAND_INTACT)
        row["rel_deficit_vs_pi4"] = (None if row["deficit"] is None
                                     else row["deficit"] / def_pi4 - 1.0)
        ladders[repr(g)] = row
        log("ladder at g={:.9g}: n*={}, deficit {}".format(
            g, row["n_star"], row["deficit"]))
    data["ladder_per_clock"] = ladders
    data["deficit_pi4_artifact"] = def_pi4
    data["intact_band"] = BAND_INTACT
    # the pi/4 deficit is the booked band reference: re-derive it and
    # anchor it to the booked value with a BACKEND-AWARE band - the
    # root-clock float64 ladder deficit is a documented floor artefact
    # (V209: the 4k+2 ratios have no limit at the closure angle; V247:
    # 1e-3 on cupy, 3e-2 labelled on the numpy fallback), so 1e-9 on
    # the float64 read would be an impossible band
    anchor = (ladders.get(repr(float(CFG.angle.branches.closed)))
              or {}).get("deficit")
    backend = str(E._BACKEND)
    anchor_band = 1e-3 if backend in ("cupy", "jax") else 3e-2
    anchor_ok = anchor is not None and abs(anchor - def_pi4) < anchor_band
    data["deficit_pi4_rederived"] = anchor
    data["deficit_pi4_anchor_band"] = anchor_band
    data["deficit_pi4_anchor_backend"] = backend
    checks.append({"name": "T1 the pi/4 deficit re-derives to the booked "
                           "V206 value {:.9f} (got {}, |d| {:.1e} < {:.0e} "
                           "backend-aware band, backend {}, the float64 "
                           "root-clock deficit is a floor artefact "
                           "V209/V247)".format(
                               def_pi4, anchor,
                               1.0 if anchor is None
                               else abs(anchor - def_pi4), anchor_band,
                               backend),
                   "pass": bool(anchor_ok)})
    if not anchor_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V313 MACHINERY-DRIFT: the pi/4 ladder deficit does not "
              "reproduce the booked V206 reference; the scan is not "
              "interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1
    checks.append({"name": "ladder intactness (band {:.2f} on the "
                           "relative deficit vs the pi/4 booked {}): {}"
                           .format(BAND_INTACT, def_pi4,
                                   {k: "intact" if v["intact"]
                                    else "collapsed"
                                    for k, v in ladders.items()}),
                   "pass": True})

    # ---- the T minimization (V285 worker, seed 285) ------------------------
    V285._ST.update({"Bmat": Bmat, "cH": np.zeros(160),
                     "NULL": NULLc,
                     "BOX": math.sqrt(V285.i3(H, H) / 3.0),
                     "BR": (b1, bm1)})
    vH = V285.vec3(H)
    cH = np.array([float(np.dot(vH, Bmat[j])) for j in range(160)])
    V285._ST["cH"] = cH
    c_ref = cH.copy()
    J0 = np.zeros((len(V285._resid(c_ref)), 160))
    eps = 1e-6
    for j in range(160):
        cp = c_ref.copy()
        cp[j] += eps
        J0[:, j] = (V285._resid(cp) - V285._resid(c_ref)) / eps
    Q0, R0 = np.linalg.qr(J0)
    V285._ST.update({"J0": J0, "QR": (Q0, R0)})
    log("retract machinery ready; {} starts x maxfev {}".format(NSTARTS,
                                                                MAXFEV))

    # V285's _start_worker reads its own module's NSTARTS/SEED/MAXFEV
    # constants for the start offset; MAXFEV is patched to the booked 200
    V285.MAXFEV = MAXFEV
    rows_s = []
    for b0 in range(0, NSTARTS, 5):
        rows_s += parpool.pmap(V285._start_worker,
                               list(range(b0, min(b0 + 5, NSTARTS))),
                               chunk=1)
        log("search {:2d}/{} starts, best T {:.6e}".format(
            len(rows_s), NSTARTS, min(r["T"] for r in rows_s)))
    rows_s.sort(key=lambda r: r["T"])
    best = rows_s[0]
    x = np.array(best["x"])
    c_star, rres = V285._retract(V285._ST["cH"] + x @ V285._ST["NULL"])
    tr_star = V285.unvec(c_star @ V285._ST["Bmat"])
    fB_s = (V285.i3(tr_star, Bop) / V285.i3(Bop, Bop)) ** 2 \
        * V285.i3(Bop, Bop) / V285.i3(tr_star, tr_star)
    fK_s = (V285.i3(tr_star, Kop) / V285.i3(Kop, Kop)) ** 2 \
        * V285.i3(Kop, Kop) / V285.i3(tr_star, tr_star)
    dist = math.sqrt(V285.i3(np.array(tr_star) - np.array(H),
                             np.array(tr_star) - np.array(H))
                     / V285.i3(H, H))
    g53 = np.zeros(53)
    hstep = 1e-5
    for a in range(53):
        xp = x.copy()
        xm = x.copy()
        xp[a] += hstep
        xm[a] -= hstep
        g53[a] = (V285._obj(xp) - V285._obj(xm)) / (2 * hstep)
    n_box = int(np.sum(np.abs(x) > V285._ST["BOX"] - 1e-6))
    interior = (float(np.linalg.norm(g53))
                < 1e-6 * max(1.0, T_committed) and n_box == 0)
    tol_star = max(1e-12, 1e-9 * abs(best["T"]))
    n_distinct = len({round(r["T"], 12) for r in rows_s
                      if abs(r["T"] - best["T"]) < tol_star})
    unique = (best["T"] < T_committed and interior and n_distinct == 1)
    G = 1.0 / (4.0 * math.log(45.0))
    data["search"] = {"best_T": best["T"], "T_committed": T_committed,
                      "hop_star": fB_s, "dressing_star": fK_s,
                      "distance_rel": dist,
                      "grad_norm": float(np.linalg.norm(g53)),
                      "box_face_coords": n_box, "interior": interior,
                      "n_starts_at_best": n_distinct,
                      "unique_minimum": unique,
                      "top5": [r["T"] for r in rows_s[:5]],
                      "retract_defect": float(np.linalg.norm(rres)),
                      "x_best": [float(v) for v in x],
                      "all_T": [r["T"] for r in rows_s]}
    data["bond_entropy_weight"] = {"G": G,
                                   "form": "G = 1/(4 ln 45), reported "
                                           "as a weight only (never "
                                           "refit)",
                                   "G_times_T_star": G * best["T"],
                                   "G_times_T_committed": G
                                   * T_committed}
    checks.append({"name": "selector test: T* {:.6e} vs T_committed "
                           "{:.6e}, interior {} (grad {:.2e}, box {}), "
                           "starts at T*: {}, unique: {}; hop {} dressing "
                           "{} at distance {:.4f}".format(
                               best["T"], T_committed, interior,
                               np.linalg.norm(g53), n_box, n_distinct,
                               unique, fB_s, fK_s, dist),
                   "pass": True})
    checks.append({"name": "bond-entropy weight G = {:.6f} (1/(4 ln 45)) "
                           "as a weight only: G*T* {:.6e}, G*T_committed "
                           "{:.6e}".format(
                               G, G * best["T"], G * T_committed),
                   "pass": True})

    # ---- conjunction -------------------------------------------------------
    intact_clocks = [g for g in SCAN if ladders[repr(g)]["intact"]]
    selected = unique and any(per_clock[repr(g)]["family_dim"] > 0
                              for g in intact_clocks)
    branch = ("LIGHTCONE-SELECTED" if selected
              else "LIGHTCONE-CLOCK-CONFLICT")
    data["status"] = branch
    data["intact_clocks"] = intact_clocks
    data["conjunction"] = {
        "family_present_at_intact_clocks": bool(
            any(per_clock[repr(g)]["family_dim"] > 0
                for g in intact_clocks)),
        "selector_unique": bool(unique)}
    print("V313 {}: family dim {} at every scan clock (clock-blind "
          "construction); ladder-intact clocks {}; selector unique: {} "
          "(T* {:.6e} vs T_committed {:.6e}, interior {}, starts at T* "
          "{}). The dictionary {} the ladder.".format(
              branch, family_dim, intact_clocks, unique, best["T"],
              T_committed, interior, n_distinct,
              "conjunctions with" if selected else "cannot be "
              "conjoined with"), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V313_lightcone_dictionary_conjunction", "tier": 3,
           "property": "does any member of the exact-lightcone family "
                       "live at a ladder-intact clock and get uniquely "
                       "selected by locality (tier 3, F8/O10)",
           "residual": True,
           "programme": "pillar-surgery V313",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v313_lightcone_dictionary_conjunction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V313 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
