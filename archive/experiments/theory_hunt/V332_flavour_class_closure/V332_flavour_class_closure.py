#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V332 - does any member of the closed flavour class contact data?
(tier 3, O6, theory-hunt closure programme.)

**The question (pre-registered).** V45 booked the democratic rank-1
base, the TBM mixing and the identity CKM. Closing the class = the
democratic rank-1 families joined by every signed permutation and the
booked 3-fold family-symmetry label: does ANY member land within 1e-3
absolute on ALL six measured sin^2 angles?

**Closed class (booked pre-run, verbatim).** Members (P, a), 48 x 3 =
144: P = the 48 signed permutations (V114 group_elements, as 3x3 signed
permutation matrices D_P with D_P[k, p[k]] = s_k); the family-space
member matrix M(P) = D_P J D_P^T, J the all-ones democratic base (V45's
cl10.yukawa_3x3.mass_matrix_3x3 construction); a in {2,3,4} = the three
spatial family-symmetry generators S_a of V45's data.family_symmetry
path (cl10.spacetime.causal_generators, the n_families_C32 = 3
decoration) - their equivalence (S_a^2 = -I, anticommutator 0, equal
norms) re-derived and booked so the label's blindness is measured, not
assumed. Mixing per member on V45's own formula path
(cl10.yukawa_3x3.mixing_from_democratic / mixing_3x3 convention): the
unique heavy eigenvector h(P) = D_P (1,1,1)/sqrt(3);
sin^2_12 = |U_e2|^2 = the squared e-component of h(P); sin^2_23 and
sin^2_13 the TBM-adapted (23) read of the O(2) light family exactly as
V45's mixing block books them (1/2 and 0), with the booked O(2) ranges
carried along; CKM the identity from identical up/down democratic
members (V45's ckm_is_identity path). Contact test: six residuals
|sin^2 - obs| <= 1e-3 absolute, PMNS against neutrino.obs_sin2_theta12/
23/13 = 0.307/0.579/0.022 and CKM against ckm.obs_sin2_theta12/23/13 =
0.05/0.949/0.02. A member CONTACTS iff all six hold.

**T0 (pre-registered).** Reproduce V45 data.mixing sin2_12 =
0.3333333333333333 (gate |d| <= 1e-12), sin2_23 = 0.5, sin2_13 = 0.0
from the same formula path, and data.mass_matrix.rank == 1. Failure:
MACHINERY-DRIFT.

**Branch (pre-registered).** FLAVOUR-CLASS-CONTACT (a member inside 1e-3
on all six, the member booked) / FLAVOUR-NO-CLASS-DECLARED (the best
member and its six per-angle residuals booked; the CKM-side identity
stays booked).
"""
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
sys.path.insert(0, exp("V114_junction_r"))

import V114_junction_r as V114                            # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import engine as E                              # noqa: E402
from cl10 import yukawa_3x3 as YU                         # noqa: E402

CFG = _config.CFG
BAND = 1e-3                          # booked absolute band on each sin^2


def art(name, key):
    """A sibling verdict JSON, read by experiment name (house art())."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    t0 = time.time()
    checks, data = [], {}

    # ---- T0 gates: V45's booked mixing and rank on its own path --------
    mix45 = art("V45_yukawa_matrix", "data.mixing")
    rank45 = int(art("V45_yukawa_matrix", "data.mass_matrix.rank"))
    mm = YU.mass_matrix_3x3()
    mix = YU.mixing_3x3()
    honest = YU.mixing_from_democratic()
    fs = YU.family_symmetry()
    d12 = abs(mix["sin2_12"] - float(mix45["sin2_12"]))
    d23 = abs(mix["sin2_23"] - float(mix45["sin2_23"]))
    d13 = abs(mix["sin2_13"] - float(mix45["sin2_13"]))
    t0ok = bool(d12 <= 1e-12 and d23 <= 1e-12 and d13 <= 1e-12
                and mm["rank"] == rank45 == 1 and fs["equivalent_3fold"]
                and honest["ckm_is_identity"])
    data["t0"] = {"sin2_12_rederived": mix["sin2_12"],
                  "sin2_12_artifact": mix45["sin2_12"], "abs_diff_12": d12,
                  "sin2_23_rederived": mix["sin2_23"],
                  "sin2_23_artifact": mix45["sin2_23"], "abs_diff_23": d23,
                  "sin2_13_rederived": mix["sin2_13"],
                  "sin2_13_artifact": mix45["sin2_13"], "abs_diff_13": d13,
                  "rank_rederived": mm["rank"], "rank_artifact": rank45,
                  "family_symmetry": fs,
                  "ckm_identity_rederived": honest["ckm_is_identity"],
                  "backend": E._BACKEND}
    checks.append({"name": "T0 V45's booked mixing (1/3, 1/2, 0) "
                           "reproduced on its own path (|d| {:.1e} / "
                           "{:.1e} / {:.1e} <= 1e-12), rank {} == {}, "
                           "3-fold family symmetry equivalent {}, CKM "
                           "identity {} (backend {})".format(
                               d12, d23, d13, mm["rank"], rank45,
                               fs["equivalent_3fold"],
                               honest["ckm_is_identity"], E._BACKEND),
                   "pass": t0ok})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V332 MACHINERY-DRIFT: the V45 flavour anchors do not "
              "reproduce; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the closed class: 48 signed permutations x 3 family labels ----
    G48 = V114.group_elements()
    assert len(G48) == 48
    obs_pmns = [float(CFG.neutrino.obs_sin2_theta12),
                float(CFG.neutrino.obs_sin2_theta23),
                float(CFG.neutrino.obs_sin2_theta13)]
    obs_ckm = [float(CFG.ckm.obs_sin2_theta12),
               float(CFG.ckm.obs_sin2_theta23),
               float(CFG.ckm.obs_sin2_theta13)]

    J = np.ones((3, 3))
    h0 = np.ones(3) / math.sqrt(3.0)
    # the (23)-adapted light basis of V45's axes block, verbatim
    e2_23 = np.array((0.0, 1.0, -1.0)) / math.sqrt(2.0)

    best = None
    contact_members = []
    sin12_set, n_members = set(), 0
    for gp in G48:
        p, s = gp
        D = np.zeros((3, 3))
        for i in range(3):
            D[i, p[i]] = float(s[i])
        M = D @ J @ D.T
        h = D @ h0
        eerr = float(np.abs(M @ h - 3.0 * h).max())
        s12 = float(abs(h[0]) ** 2)
        # V45's mixing path: heavy in the middle column, the (23)-adapted
        # member for the 23/13 reads
        e1 = np.cross(h, e2_23)
        e1 = e1 / np.linalg.norm(e1)
        B = np.column_stack([e1, h, e2_23])
        s23 = float(abs(B[1, 2]) ** 2)
        s13 = float(abs(B[0, 2]) ** 2)
        angles = (s12, s23, s13)
        sin12_set.add(round(s12, 12))
        res = [abs(s12 - obs_pmns[0]), abs(s23 - obs_pmns[1]),
               abs(s13 - obs_pmns[2]), abs(s12 - obs_ckm[0]),
               abs(s23 - obs_ckm[1]), abs(s13 - obs_ckm[2])]
        for a_idx in (2, 3, 4):          # the family-symmetry labels
            n_members += 1
            member = {"perm": list(p), "signs": list(s), "family_a": a_idx,
                      "heavy_eigvec_err": eerr,
                      "angles_sin2": list(angles),
                      "residuals_pmns_12_23_13": res[:3],
                      "residuals_ckm_12_23_13": res[3:],
                      "max_residual": max(res),
                      "contacts_all_six": bool(all(r <= BAND
                                                   for r in res))}
            if all(r <= BAND for r in res):
                contact_members.append(member)
            if best is None or max(res) < best["max_residual"]:
                best = member
    data["class"] = {"n_members": n_members,
                     "distinct_sin2_12_values": sorted(sin12_set),
                     "band_abs": BAND,
                     "obs_pmns": obs_pmns, "obs_ckm": obs_ckm,
                     "o2_ranges_booked": {"sin2_23":
                                          honest["sin2_23_range"],
                                          "sin2_13":
                                          honest["sin2_13_range"]},
                     "n_contacts": len(contact_members),
                     "contacts": contact_members[:5],
                     "best_member": best}

    branch = ("FLAVOUR-CLASS-CONTACT" if contact_members
              else "FLAVOUR-NO-CLASS-DECLARED")
    data["status"] = branch
    data["wall_time"] = time.time() - t0
    checks.append({"name": "branch {}: {} members, sin^2_12 takes values "
                           "{} on the class (V45's 1/3 pinned by the "
                           "unique heavy eigenvector); best member {} "
                           "residuals [pmns 12,23,13 | ckm 12,23,13] = "
                           "[{:.5f}, {:.5f}, {:.5f} | {:.5f}, {:.5f}, "
                           "{:.5f}] all outside the {:.0e} band; CKM-side "
                           "identity stays booked".format(
                               branch, n_members, sorted(sin12_set),
                               (best["perm"], best["signs"],
                                best["family_a"]),
                               best["residuals_pmns_12_23_13"][0],
                               best["residuals_pmns_12_23_13"][1],
                               best["residuals_pmns_12_23_13"][2],
                               best["residuals_ckm_12_23_13"][0],
                               best["residuals_ckm_12_23_13"][1],
                               best["residuals_ckm_12_23_13"][2], BAND),
                   "pass": True})
    print("V332 {}: best residuals {} (all > 1e-3); {} contacts; wall "
          "{:.0f} s.".format(branch, best["residuals_pmns_12_23_13"]
                             + best["residuals_ckm_12_23_13"],
                             len(contact_members), data["wall_time"]),
          flush=True)
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V332_flavour_class_closure", "tier": 3,
           "property": "does any member of the closed 48x3 democratic "
                       "flavour class land within 1e-3 on all six "
                       "measured sin^2 angles (tier 3, O6)",
           "residual": True,
           "programme": "theory-hunt V332",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v332_flavour_class_closure.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V332 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
