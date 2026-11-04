#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V258 - the falsifier conjunction of the interacting programme
(tier 5, residual: true).

**The question (pre-registered).**  Does any single theory point
satisfy (F9 dressed AND F8 ratio AND F2 both bands AND O8 ladder)
simultaneously?  The programme's roll-up as an experiment.

**The object (closed input set).**  Only the V249-V257 verdict
artifacts (cl10.paths.exp), plus the F8-ratio input computed here:
the dressed fermion cone slope (the V254 b = 16 rows) against V208's
cone machinery re-evaluated on the dressed effective step -
cone_generator.json re-diagonalized (its Clifford triple, Hermiticity
and cone subspace reproduced here, the 8-dim cone-subspace check),
c_gw = 2 g (V141's committed dressing_only_speed reading, the V208
convention), against |c_gw/c - 1| < 1e-15 (GW170817, the V208 band).
The ratio, not the absolute, is the GW170817 contact.  No other
artifact's number is recomputed.

**Measured (pre-registered).**  The conjunction table (per experiment
leg x per theory point g in {pi/4, pi/2, 0.2} at beta = pi/4; the
C_star point is never-realized per V249) and the dressed F8 ratio.

**Acceptance (pre-registered).**  Branch CONJUNCTION-ALIVE (a named
point passes all four legs; the falsifiers become decidable by a next
programme) or CONJUNCTION-DEAD (the per-experiment table is the
artifact; F2/F8/F9/O8 become theorems of the whole construction, both
single-particle and interacting).
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")
from cl10.paths import exp                                     # noqa: E402

from cl10 import coupling                                      # noqa: E402

G10 = coupling.G10


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    arts = {}
    for i in ("V249_star_module_closure_census", "V250_single_particle_walls",
              "V251_gauged_two_particle_sector",
              "V252_dynamical_self_energy_loop", "V253_block_rg_window",
              "V254_dressed_f9_contact", "V255_dressed_f2_contact",
              "V256_ladder_rung_reconstruction",
              "V257_interacting_ledger"):
        arts[i.split("_")[0]] = art(i)

    # ---- F8 dressed ratio: V208's cone machinery re-evaluated ------------
    gen = json.load(open(os.path.join(exp("V208_metric_bilinear_cone"),
                                      "cone_generator.json")))
    H = [np.array(gen["A"][k]) + 1j * np.array(gen["S"][k])
         for k in range(3)]
    herm = max(np.abs(h - h.conj().T).max() for h in H)
    S3 = sum(h @ h for h in H) / 3.0
    S3 = (S3 + S3.conj().T) / 2.0
    def1 = [H[i] @ H[j] + H[j] @ H[i] - (2 * S3 if i == j else 0.0)
            for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]]
    cliff = float(np.linalg.norm(np.concatenate(
        [np.concatenate([z.real.reshape(-1), z.imag.reshape(-1)])
         for z in def1])))
    w0 = np.linalg.eigvalsh(S3)
    top = float(w0[-1])
    cone_dim_gen = int(sum(abs(w - top) < 1e-8 for w in w0))
    slope_gen = math.sqrt(top)

    v254 = arts["V254"]
    m_rows = [r for r in v254["data"]["rows"] if r["class"] == "M"]
    qstar_row = min(m_rows, key=lambda r: math.hypot(
        r["q"][0] - math.pi / 2, r["q"][1] - math.pi / 2, r["q"][2]))
    c_dressed = qstar_row["v_max"]
    c_gw = 2.0 * G10                       # V141 dressing_only_speed,
    gw_band = 1e-15                        # V208's GW170817 contact
    ratio = c_gw / c_dressed
    f8_closed = abs(ratio - 1.0) < gw_band

    # ---- legs from the landed artifacts -----------------------------------
    f9_closed = v254["branch"] == "F9-DRESSED-CLOSED"
    o8_closed = arts["V256"]["branch"] == "LADDER-CONSISTENT"

    gs = [math.pi / 4, math.pi / 2, 0.2]
    table = []
    for g in gs:
        f2_rows = [r for r in arts["V255"]["data"]["rows"]
                   if abs(r["g"] - g) < 1e-12]
        f2_closed = any(r["in_ladder_band"]
                        and (r["in_alpha_band_std"] or r["in_alpha_band_10"])
                        for r in f2_rows)
        table.append({"beta": math.pi / 4, "g": g,
                      "c_star_point": "never-realized (V249)",
                      "F9_dressed": bool(f9_closed),
                      "F8_ratio": bool(f8_closed),
                      "F2_both_bands": bool(f2_closed),
                      "O8_ladder": bool(o8_closed),
                      "all_four": bool(f9_closed and f8_closed
                                       and f2_closed and o8_closed)})
    alive = any(t["all_four"] for t in table)

    # ---------------- checks ---------------------------------------------
    checks.append({
        "name": "T0 V208's committed cone generator reproduced here: "
                "Hermitian to {:.1e}, Clifford-regular to {:.1e}, cone "
                "subspace dim {} at slope {:.12f} (the 8-dim check and "
                "slope 1 are V208's landed numbers)".format(
                    herm, cliff, cone_dim_gen, slope_gen),
        "pass": bool(herm < 1e-13 and cliff < 1e-11
                     and cone_dim_gen == 8 and abs(slope_gen - 1) < 1e-9)})
    checks.append({
        "name": "T1 dressed F8 ratio (ratio, not absolute): c_gw = 2g = "
                "{:.6f} over the dressed matter cone c = {:.6f} (V254 "
                "b = 16 row at QSTAR) gives |c_gw/c - 1| = {:.4f} "
                "against the GW170817 band {:.0e}".format(
                    c_gw, c_dressed, abs(ratio - 1.0), gw_band),
        "pass": True})
    checks.append({
        "name": "T2 every leg is a landed verdict, not a new measurement: "
                "F9 {}, F2 per-point from V255 rows, O8 {} - F8 is the "
                "only input computed here".format(
                    arts["V254"]["branch"], arts["V256"]["branch"]),
        "pass": bool(all(a["pass"] for a in arts.values()))})
    checks.append({
        "name": "T3 the table covers exactly the pre-registered points: "
                "{} theory points at beta = pi/4, C_star recorded "
                "never-realized".format(len(table)),
        "pass": bool(len(table) == 3)})

    branch = "CONJUNCTION-ALIVE" if alive else "CONJUNCTION-DEAD"
    data["table"] = table
    data["f8"] = {"c_gw": c_gw, "c_dressed": c_dressed, "ratio": ratio,
                  "band": gw_band, "generator_herm": herm,
                  "generator_cliff": cliff, "cone_dim_gen": cone_dim_gen,
                  "slope_gen": slope_gen}
    data["leg_sources"] = {"F9": "V254", "F2": "V255", "O8": "V256",
                           "F8": "computed here from V208's generator "
                                 "+ V254's dressed slope"}
    data["ledger_branch"] = arts["V257"]["branch"]
    data["fitted_inputs"] = []
    ok = all(c["pass"] for c in checks)
    status = (
        "falsifier conjunction of the interacting programme ({:.0f} s): "
        "{} theory points x 4 legs. Dressed F8: |c_gw/c - 1| = {:.4f} "
        "(c_gw = {:.4f}, dressed cone c = {:.4f}) vs {:.0e}; F9 {}; O8 "
        "{}. Branch {}: {}.".format(
            time.time() - t0, len(table), abs(ratio - 1.0), c_gw,
            c_dressed, gw_band, arts["V254"]["branch"],
            arts["V256"]["branch"], branch,
            "a named point passes all four legs - the falsifiers become "
            "decidable by a next programme" if alive else
            "no point passes all four - F2/F8/F9/O8 stand as theorems "
            "of the whole construction, single-particle and interacting"))
    out = {"id": "V258_falsifier_conjunction",
           "tier": 5,
           "property": "the interacting programme's falsifier "
                       "conjunction (F9 dressed and F8 ratio and F2 "
                       "bands and O8 ladder) as one experiment (V258)",
           "residual": True,
           "programme": "V258",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v258_falsifier_conjunction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
