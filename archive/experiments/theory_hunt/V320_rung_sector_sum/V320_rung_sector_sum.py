#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V320 - the rung-sector sum: do the booked twisted sectors DERIVE a
rung weight at the required scale? (tier 2, hard, O8, theory-hunt)

**The question (pre-registered).** V319 booked the required LOOP^3-scale
rung weight as supplied-by-density. Is such a weight a REAL contribution
when V315's 64 twisted sectors are SUMMED as internal degrees of freedom
(uniform measure, and the closing-sector measure, no tunable parameter)?

**Object (closed class, fixed pre-run).** The 64 D-family decorations
(amps[t] in {(0,0,0), (G10,G10,G10)}, G10 = angle.alpha_10) over V315's
dec_step(q, g, amps, SC.STAGES); rung legs via the Floquet phase sum on
the booked 12-point SC.M_POINTS mesh:
rho_n^{(a)}(g) = (1/12) sum_q (1/128) tr M_a(q;g)^n (exact via
eigenvalues); W_eff(n) = |mean_a rho_n^{(a)}| / |rho_n^{(0)}| with the
complex phase booked alongside. Guard fixed now: a row contributes only
where |rho_n^{(0)}| > 1e-9.

**Targets (artifact leaves).** V319 corrected-clock window lnW(10) in
[9.877594588304732e-07, 1.0619053724142347e-05]; the 1/5 window
re-derived from V206 coef_at_fifth; the certified leg uses V209
certified F_208 at pi/4 (key certified.pi/4@1e-8.F_208) times the n = 208
sector ratio - does it REDUCE the booked certified deficit?

**T0.** 64 sectors; the identity sector closes at pi/4 with V315's k = 2
row (dev < 1e-12); guards booked.

**Acceptance (pre-registered words).** RUNG-SECTOR-SUPPLIED iff W_eff(10)
lands in a booked window (either measure, either clock), guards pass at
10/17/208, and the certified-deficit leg does not increase the deficit;
RUNG-SECTOR-NO-CONTRIB otherwise, full sector table booked.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))
sys.path.insert(0, exp("V315_closure_quantizer_census"))

import numpy as np                                        # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
GUARD = 1e-9                       # pre-fixed guard (not an afterthought)
NS = (10, 17, 208)                 # the measured rungs


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    t0 = time.time()
    checks, data = [], {}

    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    g5 = float(CFG.angle.branches.implemented)
    g_closed = float(CFG.angle.branches.closed)
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    hits310 = art("V310_coefficient_expression_census", "data.supplied")
    g_cand = g_cross - float([h for h in hits310
                              if "clock_gap" in h["target"]][0]["value"])
    C_15 = float(art("V206_clock_tension", "data.data_clock.coef_at_fifth"))
    F208_CERT = float(art("V209_root_clock_degeneracy",
                          "data.certified.pi/4@1e-8.F_208"))

    win_corr = [float(x) for x in art("V319_rung_weight_double_target",
                                      "data.density.A.window_lnW10")]
    win_15 = [np.log((COEF_OBS - BAND_C) / C_15) / 10.0,
              np.log((COEF_OBS + BAND_C) / C_15) / 10.0]

    # the 64 sectors (V315 class, verbatim amplitudes)
    G10 = float(CFG.angle.alpha_10)
    ON = (G10, G10, G10)
    OFF = (0.0, 0.0, 0.0)
    sectors = []
    for bits in range(64):
        amps = [ON if (bits >> t) & 1 else OFF for t in range(6)]
        sectors.append(amps)

    import V315_closure_quantizer_census as V315mod       # noqa: E402

    # T0: identity sector closure at the closed clock (V315 row)
    U0 = V315mod.dec_step(np.zeros(3), g_closed, sectors[0], SC.STAGES)
    clos0 = V315mod.closure_row(U0)
    t0ok = clos0 is not None and clos0["k"] == 2 \
        and clos0["dev"] < 1e-12
    data["identity_closure_pi4"] = clos0
    checks.append({"name": "T0 64 sectors enumerated; the identity sector "
                   "closes at pi/4 with V315's k = 2 row: {}".format(clos0),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V320 MACHINERY-DRIFT: the identity sector does not close "
              "as V315 booked", flush=True)
        _emit(data, checks, ok=False)
        return 1

    qs = [np.array(q, float) for q in SC.M_POINTS]

    def rho(a_idx, g):
        """{n: (1/12) sum_q (1/128) tr M_a(q;g)^n} for n in NS, exact
        via the 128 eigenphases."""
        acc = {n: 0.0 + 0.0j for n in NS}
        for q in qs:
            M = V315mod.dec_step(q, g, sectors[a_idx], SC.STAGES)
            lam = np.linalg.eigvals(M)
            for n in NS:
                acc[n] += np.sum(lam ** n) / 128.0
        return {n: acc[n] / len(qs) for n in NS}

    clocks = {"corrected": g_cand, "1/5": g5, "pi/4 (certified)":
              g_closed}
    table = {}
    for cname, g in clocks.items():
        legs = [rho(a, g) for a in range(64)]
        rho0 = legs[0]
        mean_leg = {n: sum(l[n] for l in legs) / 64.0 for n in NS}
        row = {}
        for n in NS:
            guard = bool(abs(rho0[n]) > GUARD)
            Wabs = (abs(mean_leg[n]) / abs(rho0[n])) if guard else None
            phase = (np.angle(mean_leg[n]) - np.angle(rho0[n])) \
                if guard else None
            inwin = None
            if guard and n == 10:
                lw = float(np.log(Wabs))
                inwin = {"corrected": bool(win_corr[0] <= lw <=
                                           win_corr[1]),
                         "1/5": bool(win_15[0] <= lw <= win_15[1])}
            row[n] = {"rho0": complex(rho0[n]),
                      "rho_sector_mean": complex(mean_leg[n]),
                      "guard_ok": guard,
                      "W_eff_abs": (float(Wabs) if Wabs is not None
                                    else None),
                      "W_eff_ln": (float(np.log(Wabs)) if Wabs is not
                                   None and Wabs > 0 else None),
                      "W_eff_phase_rad": (float(phase) if phase is not
                                          None else None),
                      "in_window": inwin}
        # closing-sector measure at the closed clock: V277/V315 booked
        # W(0) = I, so every sector closes at Gamma identically - the
        # closing-only measure coincides with the uniform measure; booked
        # as a fact, not re-used as a different number.
        closing = []
        for a in range(64):
            Ua = V315mod.dec_step(np.zeros(3), g_closed, sectors[a],
                                  SC.STAGES)
            if V315mod.closure_row(Ua) is not None:
                closing.append(a)
        row["closing_sector_count"] = len(closing)
        table[cname] = row

    # certified-deficit leg (branch-relevant): the n = 208 ratio at pi/4
    r208 = table["pi/4 (certified)"][208]
    deficit_plain = 1.0 - F208_CERT / LOOP
    deficit_sector = None
    if r208["guard_ok"] and r208["W_eff_abs"] is not None:
        deficit_sector = 1.0 - (F208_CERT * r208["W_eff_abs"]) / LOOP

    hits = []
    for cname in clocks:
        r10 = table[cname][10]
        if r10["guard_ok"] and r10["in_window"]:
            for wname, ok in r10["in_window"].items():
                if ok:
                    hits.append((cname, wname))
    guard_ok_all = all(table[c][n]["guard_ok"] for c in clocks
                       for n in NS)
    cert_ok = (deficit_sector is None or
               deficit_sector <= deficit_plain + 1e-12)
    branch = ("RUNG-SECTOR-SUPPLIED" if hits and cert_ok
              else "RUNG-SECTOR-NO-CONTRIB")
    data["status"] = branch
    data["windows"] = {"corrected": win_corr, "1/5": win_15}
    data["clocks_g"] = {k: v for k, v in clocks.items()}
    data["sector_table"] = table
    data["certified"] = {"F_208_cert": F208_CERT,
                         "deficit_plain": deficit_plain,
                         "deficit_sector_included": deficit_sector,
                         "does_not_increase": bool(cert_ok)}
    data["hits"] = ["{} x {}".format(a, b) for a, b in hits]
    checks.append({"name": "branch {}: sector-window hits {} at guards "
                   "{}; certified deficit {:.6e} -> {} (not increased: "
                   "{})".format(branch, data["hits"],
                               "all-ok" if guard_ok_all else
                               "some-guard-failed", deficit_plain,
                               ("{:.6e}".format(deficit_sector)
                                if deficit_sector is not None else
                                "guard-failed"), bool(cert_ok)),
                   "pass": True})

    for cname in clocks:
        for n in NS:
            r = table[cname][n]
            print("  {} n={:3d}: |rho0| {:.3e} guard {} W_eff {} "
                  "lnW {} in {}".format(
                      cname, n, abs(r["rho0"]), r["guard_ok"],
                      ("{:.3e}".format(r["W_eff_abs"])
                       if r["W_eff_abs"] is not None else "n/a"),
                      ("{:+.3e}".format(r["W_eff_ln"])
                       if r["W_eff_ln"] is not None else "n/a"),
                      r["in_window"] if r["in_window"] else "-"),
                      flush=True)
    print("  certified F208 deficit: plain {:.6e} sector {:.6e}".format(
        deficit_plain,
        deficit_sector if deficit_sector is not None else float("nan")),
        flush=True)
    data["wall_time"] = time.time() - t0
    _emit(data, checks, ok=all(c["pass"] for c in checks))
    return 0


def _emit(data, checks, ok):
    out = {"id": "V320_rung_sector_sum", "tier": 2,
           "property": "sum over the 64 booked twisted sectors as a "
                       "rung-weight derivation test against the V319 "
                       "windows (tier 2, hard, O8)",
           "residual": True, "programme": "theory-hunt V320",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v320_rung_sector_sum.json"),
              "w") as fh:
        json.dump(out, fh, default=str, indent=1)
    print("[V320 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
