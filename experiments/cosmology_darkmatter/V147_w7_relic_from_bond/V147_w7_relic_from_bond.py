#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V147 - the relic abundance from the BOND coupling of the singlet class:
is the same bond sector that closes (or fails to close) the O4 meeting also
the source of the dark-matter annihilation rate? (O7)

**Where O7 stands.**  V123 computed the relic abundance at the ladder mass
m = 118.4 GeV: the measured Omega h^2 = 0.120 requires alpha_eff = 0.086 in
sigma v = alpha^2/(pi m^2), while the root supplies only the loop factor
1/16pi = 0.0199 (Omega h^2 = 0.0246) or g10^2/4pi = 0.0032 (0.798); with the
nodal-plane final-state enhancement D = 1485 they move to 0.081 and 0.156 -
never 0.120. The named portals (V75 pair portal, V88 exact-Z2) reproduce
0.120 only because their lambda was inverted from the fiducial cross
section.

**The candidate.**  The DM scalar is the singlet fragment of the 16: the
class (X, Y) = (1, 0) of multiplicity 1 - the same class whose bond number
kappa_s enters the V96 meeting (index 2 of the six classes, broken flux
K = 40 by V137, superexchange kappa~ = -0.2534 by V145). If the scalar
annihilates through its bond coupling into the automaton's own light
(nodal-plane) sector, the annihilation strength is

    alpha_eff = kappa_singlet^2 / (4 pi)

(the standing convention alpha = g^2/4pi), sigma v = alpha_eff^2/(pi m^2),
and the relic follows from V123's freeze-out integral with the measured
D. The test is DECISIVE because kappa_singlet is not free: it is V96's
stored value for each of the ten node counts that close the meeting, and
the two derived patterns kappa~ (V145) and K (V137). A hit - Omega h^2
within 0.01 of 0.120 with the doubler enhancement for the same n that
closes the meeting - converts the six named bond numbers into one sector
with two independent hits.

**Checks (computed).**
1. Regression: V123's freeze-out solver reproduces the WIMP miracle
   (Omega h^2 = 0.10 at 100 GeV, sigma v = 3e-26) and the required
   alpha_eff = 0.086 at the ladder mass.
2. The doubler enhancement D recomputed from the implemented step equals
   V123's stored 1485.
3. The singlet class is identified: (X, Y) = (1, 0), multiplicity 1,
   K_singlet = 40, kappa~_singlet = -0.2534.
4. The relic for every candidate coupling: the ten V96 kappa_singlet,
   kappa~, K, the V103 junction ratio r = 1.914, the loop factor, the root
   coupling - with and without D; and the alpha_eff that 0.120 requires.
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
from cl10.paths import exp

_spec = importlib.util.spec_from_file_location(
    "v123", exp("V123_w7_relic_abundance", "V123_w7_relic_abundance.py"))
V123 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V123)


def main():
    t0 = time.time()
    out = {"experiment": "V147_w7_relic_from_bond", "data": {}}
    checks = []

    M_DM, OMEGA_OBS = V123.M_DM, V123.OMEGA_OBS
    G2CM3 = V123.GEV2_TO_CM3S

    # ---- T0 regression: solver + required alpha ---------------------------
    om100 = V123.omega_of_sigma(100.0, V123.SIGMA_FIDUCIAL)
    sig_req = V123.sigma_for_omega(M_DM)
    alpha_req = math.sqrt(sig_req / G2CM3 * math.pi * M_DM ** 2)
    checks.append(("T0 regression on V123: Omega h^2 = %.3f at 100 GeV with "
                   "sigma v = 3e-26, and alpha_eff = %.4f required at the "
                   "ladder mass m = %.1f GeV" % (om100, alpha_req, M_DM),
                   abs(om100 - 0.10) < 0.030
                   and abs(alpha_req - 0.00864) < 0.001))

    # ---- T1 the doubler enhancement D -------------------------------------
    rows = V123.node_velocities()
    speeds = [b for r in rows for b in r["branch_speeds"] if b > 1e-6]
    D = sum(1.0 / s ** 3 for s in speeds)
    v123 = json.load(open(exp("V123_w7_relic_abundance",
                                      "results",
                                      "v123_w7_relic_abundance.json")))
    Dstored = v123["data"]["doublers"]["D_planar"]
    checks.append(("T1 the nodal-plane enhancement recomputed here: %d nodes, "
                   "D = %.1f, matching V123's stored %.1f to %.1e"
                   % (len(rows), D, Dstored, abs(D - Dstored)),
                   len(rows) == 12 and abs(D - Dstored) < 1.0))

    # ---- T2 the singlet class ---------------------------------------------
    v145 = json.load(open(exp("V145_w4_superexchange_kappas",
                                      "results",
                                      "v145_w4_superexchange_kappas.json")))
    classes = v145["data"]["classes"]
    kapc = v145["data"]["superexchange_per_class"]
    v137 = json.load(open(exp("V137_w4_broken_flux",
                                      "results",
                                      "v137_w4_broken_flux.json")))
    Kc = v137["data"]["broken_flux_per_class"]
    sing = [i for i, (cc, m) in enumerate(classes)
            if m == 1 and abs(cc[1]) < 1e-9 and abs(cc[0] - 1.0) < 1e-9]
    i_s = sing[0]
    checks.append(("T2 the singlet fragment class (X, Y) = (1, 0) is class "
                   "index %d of six: multiplicity 1, K_singlet = %.1f, "
                   "kappa~_singlet = %.4f"
                   % (i_s, Kc[i_s], kapc[i_s]),
                   len(sing) == 1 and abs(Kc[i_s] - 40.0) < 1e-6
                   and abs(kapc[i_s] + 0.253431) < 1e-4))

    # ---- T3 the relic for every candidate coupling ------------------------
    v96 = json.load(open(exp("V96_kappa_s_bond",
                                      "results", "v96_kappa_s_bond.json")))
    v103 = json.load(open(exp("V103_w2_democratic",
                                      "results",
                                      "v103_w2_democratic.json")))
    r_atm = v103["data"]["r_atm"]

    def relic(alpha, with_D=True):
        sv = alpha ** 2 / (math.pi * M_DM ** 2)
        if with_D:
            sv *= D
        return V123.freeze_out(M_DM, sv)[1]

    cands = []
    for sol in v96["data"]["solutions"]:
        if "kappas" not in sol:
            continue
        kap = sol["kappas"][i_s]
        a = kap ** 2 / (4 * math.pi)
        cands.append({"name": "V96 kappa_singlet (n = %d)" % sol["n_nodes"],
                      "kappa": kap, "alpha": a,
                      "Omega": relic(a, False), "Omega_D": relic(a, True)})
    for name, kap in (("V145 superexchange kappa~", kapc[i_s]),
                      ("V137 broken flux K", Kc[i_s]),
                      ("V103 junction ratio r", r_atm)):
        a = kap ** 2 / (4 * math.pi)
        cands.append({"name": name, "kappa": kap, "alpha": a,
                      "Omega": relic(a, False), "Omega_D": relic(a, True)})
    for name, a in (("loop factor 1/16pi", V123.LOOP),
                    ("root g10^2/4pi", V123.G10 ** 2 / (4 * math.pi))):
        cands.append({"name": name, "kappa": None, "alpha": a,
                      "Omega": relic(a, False), "Omega_D": relic(a, True)})
    hits = [c for c in cands if abs(c["Omega_D"] - OMEGA_OBS) < 0.01]
    best = min(cands, key=lambda c: abs(c["Omega_D"] - OMEGA_OBS))
    # the alpha_eff that 0.120 needs, with and without D
    alpha_req_D = alpha_req / math.sqrt(D)
    checks.append(("T3 the relic for every candidate: %d candidates evaluated; "
                   "closest to 0.120 with the doubler enhancement: %s "
                   "(alpha = %.4f, Omega h^2 = %.4f, |dOmega| = %.4f); %d "
                   "candidates hit within 0.01; the measured abundance needs "
                   "alpha_eff = %.4f without D, %.5f with D"
                   % (len(cands), best["name"], best["alpha"],
                      best["Omega_D"], abs(best["Omega_D"] - OMEGA_OBS),
                      len(hits), alpha_req, alpha_req_D),
                   len(cands) == 15))

    out["data"] = {
        "m_DM_GeV": M_DM, "D": D, "alpha_required": alpha_req,
        "alpha_required_with_D": alpha_req_D,
        "singlet_class_index": i_s,
        "candidates": [{k: (round(v, 6) if isinstance(v, float) else v)
                        for k, v in c.items()} for c in cands],
        "n_hits": len(hits),
        "best": {k: (round(v, 6) if isinstance(v, float) else v)
                 for k, v in best.items()},
    }
    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    if hits:
        status = ("THE BOND SECTOR HITS TWICE: kappa_singlet = %s gives "
                  "Omega h^2 = %s with the nodal-plane enhancement D = %.0f, "
                  "within 0.01 of the measured 0.120, for the same node "
                  "counts whose meeting V96 closes. The annihilation rate is "
                  "a number of the bond layer."
                  % ([round(c["kappa"], 4) for c in hits],
                     [round(c["Omega_D"], 4) for c in hits], D))
    else:
        status = ("THE BOND COUPLING DOES NOT SET THE RELIC ABUNDANCE: every "
                  "candidate kappa_singlet - the ten V96 values that close "
                  "the meeting (%s), the superexchange -0.2534, the broken "
                  "flux 40, the junction ratio %.3f - gives, through alpha = "
                  "kappa^2/4pi and sigma v = alpha^2/(pi m^2) with D = %.0f, "
                  "an Omega h^2 between %.3g and %.3g; the closest is %s at "
                  "%.4f (|dOmega| = %.3f). The measured 0.120 needs alpha_eff "
                  "= %.4f without D or %.5f with D - no bond number of the "
                  "layer lands there. O7's portal remains named, not derived."
                  % ([round(c["kappa"], 3) for c in cands[:10]], r_atm, D,
                     min(c["Omega_D"] for c in cands),
                     max(c["Omega_D"] for c in cands),
                     best["name"], best["Omega_D"],
                     abs(best["Omega_D"] - OMEGA_OBS), alpha_req, alpha_req_D))
    out["status"] = status
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v147_w7_relic_from_bond.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    for c in cands:
        print("  %-32s kappa %-9s alpha %.4f Omega %.3g (D: %.3g)"
              % (c["name"], ("%.4f" % c["kappa"]) if c["kappa"] is not None
                 else "-", c["alpha"], c["Omega"], c["Omega_D"]))
    print("  closest: %s | hits %d/15 | wall %.1f s"
          % (best["name"], len(hits), out["wall_s"]))
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
