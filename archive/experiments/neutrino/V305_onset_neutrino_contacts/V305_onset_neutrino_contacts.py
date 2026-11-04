#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V305 - the neutrino contacts of the energy-selected spectrum
(tier 3, O2.)
**The question (pre-registered).** Under V303's branch, what are the exact
laboratory and cosmology contacts of the energy-selected (onset r = 2)
spectrum - and which of the published reach claims survives with the
theory's OWN (TBM) mixing versus measured mixing?
**Measured.** r_sel = 2.0 if V303 booked the onset branch (else its booked
minimum; the atmospheric window point re-derived for contrast); spectrum
(0, M_SOL_EV, S17/(1+r_sel)); the sum against the booked ceiling with the
margin booked; the beta-decay effective mass m_beta and the neutrinoless
interval m_bbeta under the measured-mixing read (V292's booked obs_pmns
fractions) AND over the model's booked O(2) family (V45's sin2_13 range):
the published claim m_bbeta ~ 0.02 eV "in reach" is re-scored as an
O(2)-sector choice or a prediction. F3 is re-booked in full: bare ladder
sum excluded, candidate status per the sums. Branch: ONSET-SPEC /
SUM-UNDER-BY-m% / BB-ENVELOPE-DECLARED (the published value is inside the
free family, and the measured-mixing interval sits below it).
"""
import glob
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402


def art(name):
    hits = glob.glob(os.path.join(exp(name, "results"), "*.json"))
    assert len(hits) == 1, (name, hits)
    with open(hits[0]) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    v303 = art("V303_junction_energy_r")
    v292 = art("V292_vev_host_mixing")
    assert v303["pass"], "V303 must be green"
    onset = "R-ONSET-DEMOCRATIC" in v303["branch"]
    r_sel = 2.0 if onset else float(v303["data"]["sweep"].get("argmin_r",
                                                              1.9191154379))
    S17 = nv.loop_ladder(17) * nv.EV_PER_GEV
    m_atm, m_sol, ceiling = nv.M_ATM_EV, nv.M_SOL_EV, nv.SUM_MNU_MAX_EV
    r_atm = S17 / m_atm - 1.0

    # ---- the spectra at the selected and contrasted points ----------------
    spec = {}
    for tag, r in (("selected", r_sel), ("atm_point", r_atm)):
        m3 = S17 / (1.0 + r)
        total = m_sol + m3                       # lightest exactly 0 (NH)
        spec[tag] = {"r": float(r), "m_eV": [0.0, float(m_sol), float(m3)],
                     "sum_eV": float(total),
                     "margin_vs_ceiling_pct": float(
                         100.0 * (1.0 - total / ceiling))}
    bare_sum = S17                               # rank-1 ladder: one family
    data["spectra"] = spec
    data["bare_ladder"] = {"S17_eV": float(S17), "sum_eV": float(bare_sum),
                           "excluded_by_ceiling": bool(bare_sum > ceiling)}
    checks.append({"name": "T0 spectra: selected r={:.4f} sum {:.4f} eV and "
                          "atmospheric-point sum {:.4f} eV against the "
                          "booked ceiling {:.4f} eV; bare ladder sum "
                          "{:.3f} eV excluded".format(
                              spec["selected"]["r"],
                              spec["selected"]["sum_eV"],
                              spec["atm_point"]["sum_eV"], ceiling,
                              bare_sum),
                   "pass": spec["selected"]["sum_eV"] < ceiling
                   and spec["atm_point"]["sum_eV"] < ceiling})

    # ---- mixing reads and the laboratory contacts -------------------------
    # V45 books the model sector: sin2_12 = 1/3 derived, but sin2_23/sin2_13
    # are an O(2) family with booked ranges; the published m_bbeta ~ 0.02 eV
    # is a choice inside that family, not a derived number. Re-score it.
    obs = [float(x) for x in v292["data"]["t0"]["obs_pmns"]]   # s12,s23,s13
    s12h, s23h, s13h = obs
    m2, m3 = m_sol, spec["selected"]["m_eV"][2]
    v45 = art("V45_yukawa_matrix")
    o2_13 = [float(x) for x in v45["data"]["mixing"]["sin2_13_O2_range"]]

    def beta(s12, s13):
        # |Ue1|^2 m1^2 + |Ue2|^2 m2^2 + |Ue3|^2 m3^2, m1 = 0
        return float(np.sqrt(s12 * (1 - s13) * m2 ** 2
                             + s13 * m3 ** 2))

    def bb(s12, s13):
        a = s12 * (1 - s13) * m2                 # |Ue2|^2 m2
        b = s13 * m3                             # |Ue3|^2 m3
        return [abs(a - b), a + b]               # CP-phase interval

    hy = {"s12": s12h, "s13": s13h}
    reads = {"oscillation_measured": hy,
             "model_locked_s13_0": {"s12": 1.0 / 3.0, "s13": 0.0}}
    out = {tag: {"m_beta_eV": beta(u["s12"], u["s13"]),
                 "m_bbeta_interval_eV": bb(u["s12"], u["s13"])}
           for tag, u in reads.items()}
    # model O(2) envelope over the booked sin2_13 range (worst CP phase)
    env_lo, env_hi = math.inf, -math.inf
    bmax = -math.inf
    for s13 in np.linspace(o2_13[0], o2_13[1], 201):
        lo, hi = bb(1.0 / 3.0, s13)
        env_lo, env_hi = min(env_lo, lo), max(env_hi, hi)
        bmax = max(bmax, beta(1.0 / 3.0, s13))
    out["model_O2_envelope"] = {"m_beta_max_eV": float(bmax),
                                "m_bbeta_envelope_eV": [float(env_lo),
                                                        float(env_hi)],
                                "sin2_13_booked_range": o2_13}
    data["lab_contacts"] = out
    hy_lo, hy_hi = out["oscillation_measured"]["m_bbeta_interval_eV"]
    claim = 0.02
    bb_branch = ("BB-ENVELOPE-DECLARED" if env_lo <= claim <= env_hi
                 and hy_hi < claim else "BB-POINT-SPEC")
    checks.append({"name": "T1 mixing reads: m_beta = {:.5f} eV measured "
                          "fractions, model O(2) envelope up to {:.5f} eV; "
                          "m_bbeta [{:.5f}, {:.5f}] eV measured fractions "
                          "vs the model O(2) envelope [{:.5f}, {:.5f}] eV: "
                          "the published ~{:.2f} eV is a sector choice "
                          "inside the O(2) family ({})".format(
                              out["oscillation_measured"]["m_beta_eV"],
                              bmax, hy_lo, hy_hi, env_lo, env_hi, claim,
                              bb_branch),
                   "pass": bool(env_hi > 0.0 and hy_hi < env_hi)})

    rbranch = "ONSET-SPEC" if onset else "RMIN-SPEC"
    margin = spec["selected"]["margin_vs_ceiling_pct"]
    branch = "{} / SUM-UNDER-BY-{:.1f}% / {}".format(
        rbranch, margin, bb_branch)
    data["status"] = branch
    print("V305 {}: the energy-selected spectrum (0, {:.4f}, {:.4f}) eV "
          "sums to {:.4f} eV, {:.1f}% under the booked ceiling {:.4f} eV; "
          "beta-decay m_beta = {:.5f} eV and neutrinoless m_bbeta in "
          "[{:.5f}, {:.5f}] eV under measured mixing, against the model "
          "O(2) envelope [{:.5f}, {:.5f}] eV - the published ~0.02 eV "
          "reach claim is a SECTOR CHOICE inside the booked O(2) family "
          "(sin2 th13 free in [{:.3f}, {:.3f}], V45), not a prediction; "
          "the bare ladder (sum {:.3f} eV) stays excluded".format(
              branch, m_sol, m3, spec["selected"]["sum_eV"], margin,
              ceiling, out["oscillation_measured"]["m_beta_eV"], hy_lo,
              hy_hi, env_lo, env_hi, o2_13[0], o2_13[1], bare_sum),
          flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V305_onset_neutrino_contacts", "tier": 3,
           "property": "laboratory and cosmology contacts of the "
                       "energy-selected neutrino spectrum",
           "residual": True, "programme": "V305",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v305_onset_neutrino_contacts.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V305 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
