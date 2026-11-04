#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V331 - does any (kappa, portal channel) row deliver the asymmetry?
(tier 3, O5, theory-hunt closure programme.)

**The question (pre-registered).** V199 bounded the A->BB channel
against washout; V298 booked the rung-portal channel coefficients. Over
the closed (kappa, channel) product, does ANY row put the baryon yield
inside 5% of the observed delta_b_obs = CFG.sm.delta_b_obs = 6.12e-10?

**Closed class (booked pre-run, verbatim).** kappa = the 13 booked bond
strengths of V199 data.census UNION a 400-point log grid
geomspace(0.01, 100); rows with kappa < Kappa_split or kappa >
Kappa_kill rejected (thresholds re-derived). Channels = the three key
names of V298_rung_portal_ledger data.ledger with their artifact
values. Row formula (V199's nested bodies verbatim - they are defined
inside V199 main(), so the short blocks are copied verbatim; V199 is
imported by path for lineage and its numpy-backend pin):
  Gamma_B(k) = k^2 m_A sqrt(1 - 4 m_B^2/m_A^2)/(32 pi), m_A, m_B from
  cl10.cell.masses(); r_B/H = Gamma_B/H_A with H_A = hubble(m_A) =
  hubble_coeff sqrt(g_*) T^2/M_P at T = m_A (V199's hubble);
  feed fraction Br_A_to_LH(k) = K/(K + r_B/H), K the V125
  atmospheric-point equilibrium ratio (artifact leaf).
  Y_B(k, channel) = Br_A_to_LH(k) * kappa_f(k) * coef_channel,
  kappa_f = cl10.baryogenesis.bdp_washout_efficiency(k)[0] (the BDP
  efficiency at the row's own equilibrium ratio; weak-washout domain
  efficiency 1 per the booked library convention).
Contact test: |Y_B/delta_b_obs - 1| <= 0.05 (ledger_scale_rel_tol).

**T0 (pre-registered).** Reproduce all 13 booked rows of V199
data.census: Gamma_B_GeV and r_B_over_H each to 1e-6 relative;
re-derive the thresholds against V199 data.thresholds {Br_min,
Kappa_split, Kappa_kill} to 1e-6 relative, same keys. Failure:
MACHINERY-DRIFT.

**Branch (pre-registered).** ASYM-CHANNEL-CONTACT (the (kappa*,
channel) booked with its Y_B and residual) / ASYM-NO-CHANNEL-DECLARED
(the closest (kappa, channel), its Y_B/6.12e-10 ratio and the band
multiple booked).
"""
import json
import math
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"   # V199's own convention, booked

import numpy as np                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                              # noqa: E402
sys.path.insert(0, exp("V199_a_to_bb_channel"))

import V199_a_to_bb_channel as V199_mod                 # noqa: E402,F401

from cl10 import baryogenesis as BG                     # noqa: E402
from cl10 import cell as cell                           # noqa: E402
from cl10 import config as _config                      # noqa: E402
from cl10 import u10pred                                # noqa: E402

CFG = _config.CFG


def art(name, key):
    """A sibling verdict JSON, read by experiment name (house art())."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# ---- V199's formula blocks, bodies verbatim (nested in V199 main) ------
def main():
    t0 = time.time()
    checks, data = [], {}

    masses = cell.masses()
    M_A = masses["m_A_GeV"]
    M_B = masses["m_B_GeV"]
    M_P = u10pred.M_PLANCK
    HC = float(CFG.cosmology.hubble_coeff)
    GS = float(CFG.cosmology.g_star)

    def hubble(T):
        return HC * math.sqrt(GS) * T * T / M_P

    kin = math.sqrt(max(0.0, 1.0 - 4.0 * M_B ** 2 / M_A ** 2))

    def gamma_b(k):
        return k * k * M_A * kin / (32.0 * math.pi)

    H_A = hubble(M_A)

    K = float(art("V125_cp_flux_washout",
                  "data.asymmetry.cases.atmospheric y")["K"])
    jreq = float(art("V125_cp_flux_washout",
                     "data.asymmetry.cases.atmospheric y")["J_required"])
    jav = float(art("V125_cp_flux_washout", "data.J_available_max"))
    br_min = jreq / jav

    def kappa_of_br(br):
        rB = K * (1.0 / br - 1.0)
        return math.sqrt(rB * H_A * 32.0 * math.pi / (M_A * kin))

    # ---- T0 gates: the 13 booked rows and the booked thresholds --------
    # V199 stores its kappa leaves ROUNDED (round(k, 4)) but computed the
    # booked Gamma_B/r_B_over_H at the ORIGINAL V147 candidate strengths
    # (its verbatim convention) - the re-run uses the originals and pairs
    # each with its booked rounded leaf (the mapping booked in data); the
    # census below then uses the rounded booked leaves, as pre-registered.
    rows199 = art("V199_a_to_bb_channel", "data.census")
    th199 = art("V199_a_to_bb_channel", "data.thresholds")
    cand = [c for c in art("V147_relic_from_bond", "data.candidates")
            if c["kappa"] is not None]
    kap_orig = sorted(abs(float(c["kappa"])) for c in cand)
    worst = 0.0
    kappa_map = []
    for r, k in zip(rows199, kap_orig):
        g = gamma_b(k)
        rB = g / H_A
        dg = abs(g / float(r["Gamma_B_GeV"]) - 1.0)
        dr = abs(rB / float(r["r_B_over_H"]) - 1.0)
        worst = max(worst, dg, dr)
        assert round(k, 4) == float(r["kappa"]), (k, r["kappa"])
        kappa_map.append({"booked_kappa_leaf": float(r["kappa"]),
                          "original_V147_strength": k})
    data["kappa_original_vs_leaf"] = kappa_map
    k_split = kappa_of_br(0.99)
    k_kill = kappa_of_br(br_min)
    dth = max(abs(k_split / float(th199["Kappa_split"]) - 1.0),
              abs(k_kill / float(th199["Kappa_kill"]) - 1.0),
              abs(br_min / float(th199["Br_min"]) - 1.0))
    data["t0"] = {"rows199": len(rows199),
                  "worst_rel_row_reproduction": worst,
                  "thresholds_booked": th199,
                  "thresholds_rederived": {"Br_min": br_min,
                                           "Kappa_split": k_split,
                                           "Kappa_kill": k_kill,
                                           "H_at_mA_GeV": H_A},
                  "worst_rel_threshold_reproduction": dth,
                  "K_atmo": K, "m_A_GeV": M_A, "m_B_GeV": M_B,
                  "kinematic": kin}
    checks.append({"name": "T0 all {:d} booked V199 rows re-run to {:.1e} "
                           "relative (gate 1e-6) and the booked thresholds "
                           "Br_min {:.6g} / Kappa_split {:.6g} / Kappa_kill "
                           "{:.6g} re-derived to {:.1e} relative (same "
                           "keys)".format(len(rows199), worst,
                                          float(th199["Br_min"]),
                                          float(th199["Kappa_split"]),
                                          float(th199["Kappa_kill"]),
                                          dth),
                   "pass": bool(worst <= 1e-6 and dth <= 1e-6
                                and k_split < k_kill)})
    if not (worst <= 1e-6 and dth <= 1e-6 and k_split < k_kill):
        data["status"] = "MACHINERY-DRIFT"
        print("V331 MACHINERY-DRIFT: the V199 washout machinery does not "
              "re-run; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the closed (kappa, channel) census -----------------------------
    delta = float(CFG.sm.delta_b_obs)
    band = float(CFG.pillar_surgery.ledger_scale_rel_tol)   # 5%
    channels = art("V298_rung_portal_ledger", "data.ledger")
    chans = {name: float(entry["value"])
             for name, entry in channels.items()}
    kap_set = sorted({float(r["kappa"]) for r in rows199}
                     | set(float(v)
                           for v in np.geomspace(0.01, 100.0, 400)))
    kept, rejected, best = 0, 0, None
    rows_out = []
    for k in kap_set:
        if k < k_split or k > k_kill:
            rejected += 1
            continue
        kept += 1
        rB = gamma_b(k) / H_A
        br = K / (K + rB)
        kap_f, _zB, _dom = BG.bdp_washout_efficiency(k)
        for name, coef in chans.items():
            yb = br * kap_f * coef
            resid = abs(yb / delta - 1.0)
            row = {"kappa": k, "channel": name, "Br_A_to_LH": br,
                   "kappa_f": kap_f, "coef": coef, "Y_B": yb,
                   "Y_over_delta": yb / delta,
                   "residual_over_band": resid / band,
                   "contact": bool(resid <= band)}
            rows_out.append(row)
            if best is None or resid < best["residual_over_band"] * band:
                best = row
    contacts = [r for r in rows_out if r["contact"]]
    top10 = sorted(rows_out, key=lambda r: abs(r["Y_over_delta"] - 1.0))[:10]
    data["census"] = {"n_kappa_rows": len(kap_set), "kept": kept,
                      "rejected_outside_thresholds": rejected,
                      "channels": chans, "band_rel": band,
                      "n_rows": len(rows_out), "n_contacts": len(contacts),
                      "contacts": contacts,
                      "closest": best,
                      "closest_band_multiples": best["residual_over_band"],
                      "closest_10_rows": top10}

    branch = ("ASYM-CHANNEL-CONTACT" if contacts
              else "ASYM-NO-CHANNEL-DECLARED")
    data["status"] = branch
    data["wall_time"] = time.time() - t0
    checks.append({"name": "branch {}: {} rows over {} kept kappa x {} "
                           "channels; {} contacts at the {:.0%} band; "
                           "closest Y_B/delta_b_obs = {:.3e} (kappa {}, "
                           "channel {}) = {:.0f} band-multiples".format(
                               branch, len(rows_out), kept, len(chans),
                               len(contacts), band, best["Y_over_delta"],
                               best["kappa"], best["channel"],
                               best["residual_over_band"]),
                   "pass": True})
    print("V331 {}: closest Y_B = {:.3e} at kappa {:.4g} channel {} "
          "({:.3e} x delta, {:.0f} band-multiples); kept {} rejected {}; "
          "wall {:.0f} s.".format(branch, best["Y_B"], best["kappa"],
                                  best["channel"], best["Y_over_delta"],
                                  best["residual_over_band"], kept,
                                  rejected, data["wall_time"]),
          flush=True)
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V331_channel_asymmetry_closure", "tier": 3,
           "property": "does any booked (kappa, portal-channel) row land "
                       "the baryon yield on delta_b_obs = 6.12e-10 within "
                       "the 5% band (tier 3, O5)",
           "residual": True,
           "programme": "theory-hunt V331",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v331_channel_asymmetry_closure.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V331 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
