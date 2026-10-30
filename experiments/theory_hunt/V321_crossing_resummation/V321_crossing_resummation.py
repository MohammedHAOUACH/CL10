#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V321 - the crossing resummation: canonical finite part of the 4k+2
power divergence at the root clock (tier 2, hard, O8, theory-hunt).

Fact (V209, booked): the certified 4k+2 class diverges at pi/4 (F_210
certified = 3.885e48) while the detuned class blows up as a POWER in the
detuning (slope ~11/decade). A power divergence has a canonical finite
part: subtract c*delta^-m with (c, m) MEASURED at the detuned points;
targets enter only the decision.

Object (closed): F_n(delta) = the exact-arithmetic mp_ladder of V209
(90 dps, nmax 260) at g = pi/4 + delta over the fixed grid {1e-1, 1e-2,
1e-3, 1e-4, 1e-5}; n in {2,6,10,14,18,22} plus tail probes {250,254,
258}; guards: slope agreement < 5% (and m > 0), c-spread < 1%,
fp-spread < 1%, with the subtraction constant c_n FIXED at the finest
grid point (else booked NON-POLYDIVERGENT / NON-CANONICAL, no fp).
Anchor: the convergent 4k class untouched (values booked at pi/4+1e-8,
never subtracted); the certified F_208 leaf (V209, pi/4+1e-8)
reproduced to 1e-12.

T0 (amended at run 2, booked in the note): machinery IDENTITY -
V209's exact detuned calls (nmax 12, 55 dps) replayed exactly (rel
< 1e-12); the certified F_208 anchor at pi/4+1e-8 to 1e-12; the
certified F_210 leaf quoted as fact. The 90-dps production grid is
booked as MORE precise than the 55-dps leaf (the 1e-5 leaf carries
~8 digits at 1e52; a 1e-9 leaf-agreement gate is unsatisfiable).

Acceptance (pre-registered words). RESUM-CROSSING-TIED iff n = 10 is
canonical AND |3(fp_10/LOOP)^10 - COEF_OBS| <= coef_band_abs AND the 4k
anchor holds - the first zero-fitted algebraic supply of the hierarchy
coefficient, the canonically tamed divergent class; RESUM-NO-TIE
otherwise with the full per-n table booked. One tied rung-contact does
not close the ledger (rung-17 floor, V276 residue) - the booking says so.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V209_root_clock_degeneracy"))

from cl10 import config as _config                        # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
DELAS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
NS = [2, 6, 10, 14, 18, 22, 250, 254, 258]


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    import mpmath as mp
    from V209_root_clock_degeneracy import mp_ladder, skeleton
    t0 = time.time()
    checks, data = [], {}
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)

    N, blocks, dens, D, cart, planes, TL = skeleton()
    DPS = 90                       # booked precision: the 4k+2 tail
    #                probes reach ~1e65 at delta = 1e-5; V209 ran 55 dps
    #                for n <= 12 and 50 dps at 212 - the raise is booked,
    #                machinery-only, no class change
    NMAX = 260
    exprs = {d: "mp.pi/4 + mp.mpf('{}')".format(repr(d)) for d in DELAS}

    Fs = {}
    for d in DELAS:
        r = mp_ladder(exprs[d], NMAX, DPS, planes, blocks, dens, cart)
        Fs[d] = {n: r[n]["F"] for n in NS + [4, 8, 208]}
    r_anchor = mp_ladder("mp.pi/4 + mp.mpf('1e-8')", NMAX, DPS,
                         planes, blocks, dens, cart)

    # ---- T0 (gate-hygiene amended, booked at the second run) ----------
    # the booked V209 detuned leaf at delta = 1e-5 is itself a 55-dps
    # value (true 90-dps value differs by 5.1e-8 relative - the leaf
    # carries ~8 significant digits at 1e52), so a 1e-9 agreement gate
    # against the leaf is unsatisfiable in principle. T0 is therefore
    # machinery IDENTITY: replay V209's exact call (nmax 12, 55 dps) and
    # require exact agreement; the production grid runs 90 dps because
    # the finite-part subtraction needs absolute error far below the
    # fp scale at e50-class terms (55 dps would give ~e44 absolute
    # error - garbage). Decision gate untouched.
    v209_det = art("V209_root_clock_degeneracy", "data.branch_structure"
                   ".detuned")
    t0a = {}
    for key, d in (("1e-1", 1e-1), ("1e-3", 1e-3), ("1e-5", 1e-5)):
        booked = float(v209_det[key]["2"])
        r_id = mp_ladder("mp.pi/4 + mp.mpf('{}')".format(repr(d)), 12, 55,
                         planes, blocks, dens, cart)
        got = float(abs(r_id[2]["F"]))
        t0a[key] = abs(got - booked) / booked
    cert208 = float(art("V209_root_clock_degeneracy",
                        "data.certified.pi/4@1e-8.F_208"))
    cert210 = art("V209_root_clock_degeneracy",
                  "data.certified.pi/4@1e-8.F_210")
    d208 = abs(float(abs(r_anchor[208]["F"])) - cert208)
    t0ok = all(v < 1e-12 for v in t0a.values()) and d208 < 1e-12
    data["t0_detuned_F2_rel"] = t0a
    data["t0_F208_anchor_at_1e-8_absdiff"] = d208
    data["certified_leaves"] = {"F_208": cert208, "F_210": cert210}
    data["anchor_4k_values_at_pi4_plus_1e-8"] = {
        "F_4": float(abs(r_anchor[4]["F"])),
        "F_8": float(abs(r_anchor[8]["F"])),
        "F_208": float(abs(r_anchor[208]["F"]))}
    checks.append({"name": "T0 machinery identity: V209's exact detuned "
                   "calls (nmax 12, 55 dps) replayed exactly (rel {}) and "
                   "the certified F_208 anchor at pi/4+1e-8 reproduced to "
                   "{:.1e} < 1e-12; certified F_210 quoted {:.4e}".format(
                       {k: "{:.1e}".format(v) for k, v in t0a.items()},
                       d208, cert210), "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V321 MACHINERY-DRIFT: the detuned grid does not reproduce "
              "V209", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- canonical finite parts ---------------------------------------
    # guards in mp; the subtraction constant c_n is FIXED at the finest
    # grid point, never re-fitted per point (pre-registration)
    table = {}
    for n in NS:
        pts = [(d, abs(Fs[d][n])) for d in DELAS]
        fin = [(d, v) for d, v in pts if mp.isfinite(v)]
        row = {"values": {repr(d): float(v) for d, v in fin}}
        if len(fin) < 3:
            row["status"] = "NON-CANONICAL-OVERFLOW"
            table[n] = row
            continue
        (d1, f1), (d2, f2), (d3, f3) = (fin[0], fin[2], fin[4]) \
            if len(fin) >= 5 else tuple(fin[::2])
        m_a = (mp.log(f2) - mp.log(f1)) / mp.log(d1 / d2)
        m_b = (mp.log(f3) - mp.log(f2)) / mp.log(d2 / d3)
        m = (m_a + m_b) / 2
        slope_ok = bool(abs(m_a - m_b) / max(abs(m), mp.mpf("1e-12"))
                        < 0.05 and m > 0)
        dfin, vfin = fin[-1]
        c = vfin * dfin ** m
        cs = [v * d ** m for d, v in fin]
        c_spread = float((max(cs) - min(cs)) / abs(c)) if c != 0 else 1.0
        fps = [v - c * d ** (-m) for d, v in fin]
        canonical = bool(slope_ok and c_spread < 0.01)
        fp_spread = float("nan")
        fp = None
        if canonical:
            fp = sum(fps) / len(fps)
            fp_spread = float((max(fps) - min(fps))
                              / max(abs(fp), mp.mpf("1e-30")))
            canonical = bool(fp_spread < 0.01)
        row.update({"m_3pt": [float(m_a), float(m_b)], "m": float(m),
                    "slope_ok": slope_ok, "c_spread": c_spread,
                    "fp": (float(fp) if fp is not None else None),
                    "fp_spread": fp_spread,
                    "status": ("CANONICAL" if canonical else
                               ("NON-POLYDIVERGENT" if not slope_ok
                                else "NON-CANONICAL"))})
        if fp is not None:
            row["coef_read"] = 3.0 * (abs(fp) / LOOP) ** 10
            row["coef_norm"] = abs(row["coef_read"] - COEF_OBS) / BAND_C
        table[n] = row

    r10 = table[10]
    anchor = bool(d208 < 1e-12)
    tie = (r10.get("status") == "CANONICAL"
           and r10.get("coef_norm", 1e30) <= 1.0 and anchor)
    branch = "RESUM-CROSSING-TIED" if tie else "RESUM-NO-TIE"
    data["status"] = branch
    data["table"] = table
    data["COEF_OBS"] = COEF_OBS
    data["F_req_rung10"] = LOOP * (COEF_OBS / 3.0) ** 0.1
    data["anchor_4k_untouched"] = anchor
    checks.append({"name": "branch {}: n=10 status {} coef-read {} "
                   "(norm {}), anchor {}; per-n canonical: {}".format(
                       branch, r10.get("status"),
                       ("{:.7f}".format(r10["coef_read"])
                        if "coef_read" in r10 else "n/a"),
                       ("{:.3f}".format(r10["coef_norm"])
                        if "coef_norm" in r10 else "n/a"), anchor,
                       {n: table[n]["status"] for n in NS}),
                   "pass": True})

    for n in NS:
        r = table[n]
        print("  n={:3d} {}: m {} fp {} coef {} norm {}".format(
            n, r["status"],
            ("{:.3f}".format(r["m"]) if "m" in r else "-"),
            ("{:.6e}".format(r["fp"]) if r.get("fp") is not None
             else "-"),
            ("{:.7f}".format(r["coef_read"]) if "coef_read" in r
             else "-"),
            ("{:.3f}".format(r["coef_norm"]) if "coef_norm" in r
             else "-")), flush=True)
    print("  anchor F_208@1e-8 mp {:.15g} | cert {:.15g} | F_req(10) "
          "{:.9f}".format(float(abs(r_anchor[208]["F"])), cert208,
                          LOOP * (COEF_OBS / 3.0) ** 0.1), flush=True)
    data["wall_time"] = time.time() - t0
    _emit(data, checks, ok=all(c["pass"] for c in checks))
    return 0


def _emit(data, checks, ok):
    out = {"id": "V321_crossing_resummation", "tier": 2,
           "property": "canonical finite part of the 4k+2 power "
                       "divergence at the root clock, tie test against "
                       "the hierarchy coefficient (tier 2, hard, O8)",
           "residual": True, "programme": "theory-hunt V321",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v321_crossing_resummation.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V321 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
