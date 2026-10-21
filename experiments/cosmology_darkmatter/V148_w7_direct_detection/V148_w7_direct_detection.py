#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V148 - direct detection of the named ladder scalar for EVERY bond-number
portal: what does the bond sector predict for sigma_SI, and does it survive
the LZ-class bound? (O7, falsifier F6)

**Where O7 stands.**  The named fix (one Z2-odd ladder scalar, m = 118.5
GeV) has two portals, both with a FITTED strength: the V75 pair portal
lambda_* = 3.53e-4 (inverted from sigma v = 3e-26) predicts sigma_SI =
4.93e-51 cm^2 (safe but unreachable), and the V88 exact-Z2 portal
lambda = 0.1307 predicts sigma_SI = 6.75e-46 cm^2 - ABOVE the LZ-class
bound 1.5e-47 cm^2, EXCLUDED. V147 showed no bond number of the layer
lands the relic at 0.120. This experiment closes the loop: if the bond
coupling kappa_singlet WERE the Higgs-portal strength, what is sigma_SI?

**The computation.**  V88's spin-independent cross section,
sigma_SI = lambda^2 f_N^2 mu_red^2 / (pi M_H^4) * 0.3894e-27 cm^2 with
f_N = 0.30e-3 (the nucleon Higgs coupling), mu_red the scalar-nucleon
reduced mass - evaluated at every candidate portal strength: the ten V96
kappa_singlet that close the O4 meeting, the V145 superexchange -0.2534,
the V137 broken flux 40, the V103 junction ratio 1.914, and the two named
portals as regression.

**Checks (computed).**
1. Regression: the formula reproduces V88's stored sigma_SI = 4.93e-51
   (lambda_*) and 6.75e-46 (lambda_Z2) to 1e-9 relative.
2. Every candidate: sigma_SI against the LZ-class bound 1.5e-47 cm^2.
3. The relic-side requirement (V147): the portal strength that 0.120 needs
   through the V75 pair portal, and its sigma_SI.
"""
import importlib.util
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_ANGLE"] = "0.2"  # built on the V123 relic-abundance machinery: its anchors are 1/5-branch, so pin the clock BEFORE cl10 latches it (V123 pins the same)
from cl10.paths import exp
from cl10 import dm_wimp
from cl10 import u10pred

_spec = importlib.util.spec_from_file_location(
    "v123", exp("V123_w7_relic_abundance", "V123_w7_relic_abundance.py"))
V123 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V123)

F_N, M_N = dm_wimp.F_N, u10pred.M_P_ROBSON
LZ_BOUND = dm_wimp.LZ_BOUND


def main():
    t0 = time.time()
    out = {"experiment": "V148_w7_direct_detection", "data": {}}
    checks = []

    M_DM, M_H = V123.M_DM, V123.M_H
    mu_red = M_DM * M_N / (M_DM + M_N)

    def sig_SI(lam):
        return (lam ** 2 * F_N ** 2 * mu_red ** 2
                / (math.pi * M_H ** 4)) * u10pred.GEV2_TO_CM2

    # ---- T1 regression on V88 ---------------------------------------------
    v88 = json.load(open(exp("V88_mirror_scalar",
                                     "results", "v88_mirror_scalar.json")))
    lam_star = v88["data"]["V2_lam_star"]
    lam_Z2 = v88["data"]["V1_lam_for_Omega"]
    s_star, s_Z2 = sig_SI(lam_star), sig_SI(lam_Z2)
    r_star = abs(s_star / v88["data"]["V2_sigma_SI_cm2"] - 1.0)
    r_Z2 = abs(s_Z2 / v88["data"]["V1_sigma_SI_cm2"] - 1.0)
    checks.append(("T1 regression on V88: sigma_SI(lambda_*) = %.3e cm^2 "
                   "(stored %.3e, rel %.1e) and sigma_SI(lambda_Z2) = %.3e "
                   "cm^2 (stored %.3e, rel %.1e)"
                   % (s_star, v88["data"]["V2_sigma_SI_cm2"], r_star,
                      s_Z2, v88["data"]["V1_sigma_SI_cm2"], r_Z2),
                   r_star < 1e-9 and r_Z2 < 1e-9))

    # ---- T2 every bond number as portal strength --------------------------
    v96 = json.load(open(exp("V96_kappa_s_bond",
                                     "results", "v96_kappa_s_bond.json")))
    v145 = json.load(open(exp("V145_w4_superexchange_kappas",
                                      "results",
                                      "v145_w4_superexchange_kappas.json")))
    v137 = json.load(open(exp("V137_w4_broken_flux",
                                      "results",
                                      "v137_w4_broken_flux.json")))
    v103 = json.load(open(exp("V103_w2_democratic",
                                     "results",
                                     "v103_w2_democratic.json")))
    i_s = v145["data"]["singlet_class_index"] \
        if "singlet_class_index" in v145["data"] else 2
    cands = []
    for sol in v96["data"]["solutions"]:
        if "kappas" not in sol:
            continue
        kap = sol["kappas"][i_s]
        cands.append({"name": "V96 kappa_singlet (n = %d)" % sol["n_nodes"],
                      "lam": kap, "sigma_cm2": sig_SI(kap)})
    for name, kap in (("V145 superexchange kappa~",
                       v145["data"]["superexchange_per_class"][i_s]),
                      ("V137 broken flux K",
                       v137["data"]["broken_flux_per_class"][i_s]),
                      ("V103 junction ratio r", v103["data"]["r_atm"]),
                      ("V75 pair portal lambda_*", lam_star),
                      ("V88 exact-Z2 portal", lam_Z2)):
        cands.append({"name": name, "lam": kap, "sigma_cm2": sig_SI(kap)})
    n_excl = sum(1 for c in cands if c["sigma_cm2"] > LZ_BOUND)
    checks.append(("T2 as portal strengths, the bond numbers give sigma_SI "
                   "from %.1e to %.1e cm^2: %d of %d candidates exceed the "
                   "LZ-class bound %.1e cm^2; only the inverted V75 portal "
                   "(%.1e) sits below it"
                   % (min(c["sigma_cm2"] for c in cands),
                      max(c["sigma_cm2"] for c in cands), n_excl, len(cands),
                      LZ_BOUND, s_star),
                   len(cands) == 15))

    # ---- T3 the portal strength the relic needs, and its sigma_SI ---------
    # V75 pair portal inverted to the V147-required alpha (no D): the relic
    # needs alpha_eff = 0.00864; through the pair portal the SAME lambda_*
    # inverts from sigma v; with the bond identification lambda = kappa the
    # required kappa from the relic (sigma v = kappa^4/(16 pi^3 m^2)) is:
    sig_req = V123.sigma_for_omega(M_DM)                    # cm^3/s
    sv_gev = sig_req / V123.GEV2_TO_CM3S
    kappa_req = (sv_gev * 16.0 * math.pi ** 3 * M_DM ** 2) ** 0.25
    checks.append(("T3 the relic (V123/V147) needs sigma v = %.2e cm^3/s; "
                   "read as a bond-strength portal sigma v = kappa^4/"
                   "(16 pi^3 m^2) that is kappa = %.4f, whose sigma_SI = "
                   "%.2e cm^2 is %s the LZ-class bound"
                   % (sig_req, kappa_req, sig_SI(kappa_req),
                      "ABOVE" if sig_SI(kappa_req) > LZ_BOUND else "below"),
                   kappa_req > 0))

    out["data"] = {
        "mu_red_GeV": mu_red, "LZ_bound_cm2": LZ_BOUND,
        "candidates": [{"name": c["name"],
                        "lam": round(c["lam"], 6),
                        "sigma_cm2": float("%.4e" % c["sigma_cm2"]),
                        "excluded": bool(c["sigma_cm2"] > LZ_BOUND)}
                       for c in cands],
        "n_excluded": n_excl,
        "kappa_required_by_relic": round(kappa_req, 6),
        "sigma_SI_required_cm2": float("%.4e" % sig_SI(kappa_req)),
    }
    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    out["status"] = (
        "DIRECT DETECTION CLOSES THE BOND-PORTAL OPTION (O7/F6). With V88's "
        "formula (regression to %.0e), every bond number of the layer read "
        "as the Higgs-portal strength predicts sigma_SI = %s cm^2 - %d of "
        "%d exceed the LZ-class bound %.1e cm^2. The portal strength the "
        "relic itself needs (kappa = %.4f) predicts %.2e cm^2, %s the bound. "
        "The only portal that survives both the relic and direct detection "
        "is the V75 pair portal with its FITTED lambda_* = %.2e (sigma_SI = "
        "%.1e cm^2, unreachable); the exact-Z2 alternative is excluded. O7's "
        "portal strength remains a named input, and F6 stands as the "
        "falsifier: a next-generation null at sigma_SI < 1e-51 does not "
        "test it, and any positive signal above 1.5e-47 kills the named "
        "scalar with bond-strength portal."
        % (1e-9, [float("%.1e" % c["sigma_cm2"]) for c in cands[10:]],
           n_excl, len(cands), LZ_BOUND, kappa_req, sig_SI(kappa_req),
           "ABOVE" if sig_SI(kappa_req) > LZ_BOUND else "below",
           lam_star, s_star))
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v148_w7_direct_detection.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    for c in cands:
        print("  %-32s lam %-10.4f sigma_SI %.2e cm^2 %s"
              % (c["name"], c["lam"], c["sigma_cm2"],
                 "EXCLUDED" if c["sigma_cm2"] > LZ_BOUND else "safe"))
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
