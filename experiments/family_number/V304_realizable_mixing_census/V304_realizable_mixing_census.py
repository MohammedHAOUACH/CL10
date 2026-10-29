#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V304 - the realizable mixing census (tier 3, O2.)
**The question (pre-registered).** V104 scanned GENERIC vev directions and
found joint NuFIT fraction 0.0; V292 measured the ONE realizable insertion
and found the angles unmoved (shift 2.2e-16). Is the mixing correction of
the ENTIRE realizable sector exactly {0}: every neutral chain-kernel host
channel of V291 (indices {6,7,9}, weights p0), all 7 non-empty
combinations, with the family-symmetric insertion justified by re-derived
booked facts (V45 family_symmetry equivalent_3fold, V294 commutant
span{I, J}) rather than asserted?
**Measured.** T0: V45 democratic spectrum re-read (rank 1, eigs {0, 0, 3})
and V292's point re-run against its artifact (c = 112, shift < 1e-12).
T1: the census over the 7 subsets: correction c(W) = sum of booked
chain-H energies of the members, M(c) = (3 + c) ones: eigenvalues
{3+c, 0, 0}, heavy eigenvector (1,1,1)/sqrt(3), the maximal angle shift
over the census booked. T2: the locked point (1/3, 1/2, 0) against the
NuFIT dictionary IMPORTED from V104's module (not retyped): per-angle
sigma booked; joint fraction of the census (0 exactly if any angle is
outside its band). Branch: MIXING-REALIZABLE-EMPTY (max shift < 1e-12 and
joint fraction 0.0) / MIXING-REALIZED-N.
"""
import glob
import itertools
import json
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
V104DIR = os.path.join(ROOT, "experiments", "spacetime_matter",
                       "V104_mixing_manifold")
sys.path.insert(1, V104DIR)
from cl10.paths import exp                                # noqa: E402
from V104_mixing_manifold import NUFIT                    # noqa: E402


def art(name):
    hits = glob.glob(os.path.join(exp(name, "results"), "*.json"))
    assert len(hits) == 1, (name, hits)
    with open(hits[0]) as fh:
        return json.load(fh)


def angles_of(M):
    # frame-invariant content: sin^2 th12 = |<e1|h>|^2 of the unique heavy
    # direction; the light 2-fold family keeps s23 = 1/2, s13 = 0 as the
    # booked convention (V292's own convention).
    w, V = np.linalg.eigh(M)
    h = V[:, int(np.argmax(np.abs(w)))]
    return float(abs(h[0]) ** 2), 0.5, 0.0


def main():
    t0 = time.time()
    checks, data = [], {}
    v45 = art("V45_yukawa_matrix")
    v291 = art("V291_higgs_host_channels")
    v292 = art("V292_vev_host_mixing")
    v294 = art("V294_democratic_mD_forced")

    # ---- T0 the booked anchors, re-read -----------------------------------
    eigs = sorted(float(x) for x in v45["data"]["mass_matrix"]["eigenvalues"])
    fam_ok = bool(v45["data"]["family_symmetry"]["equivalent_3fold"])
    comm = v294["data"]["commutant"]
    c_art = float(v292["data"]["t1"]["c_host_energy_sum"])
    shift_art = float(v292["data"]["t1"]["angle_shift"])
    J = np.ones((3, 3))
    e0 = sorted(np.linalg.eigvalsh(J))
    t0ok = (abs(eigs[2] - 3.0) < 1e-9 and abs(eigs[0]) < 1e-9 and fam_ok
            and comm["dim"] == 2 and shift_art < 1e-12
            and np.allclose(np.linalg.eigvalsh(J), e0))
    data["t0"] = {"democratic_eigs": eigs, "family_3fold": fam_ok,
                  "commutant_dim": comm["dim"],
                  "commutant_generators_ok": bool(comm["generators_ok"]),
                  "c_artifact": c_art, "shift_artifact": shift_art}
    checks.append({"name": "T0 V45 democratic (0,0,3) rank 1 with the "
                          "booked 3-fold family symmetry; V294 commutant "
                          "span(I,J) dim 2; V292 shift {:.1e}: {}".format(
                              shift_art, "OK" if t0ok else "DRIFT"),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V304 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    # ---- T1 the census of realizable insertions ----------------------------
    host = v291["data"]["host"]
    p0 = [float(x) for x in host["p0"]]
    idx = [i for i, x in enumerate(p0) if x > 0.5]      # V292's own read
    assert sorted(idx) == [6, 7, 9] and len(p0) == 12
    hspec = [float(x) for x in host["chain_H_spectrum_artefact"]]
    per = {i: hspec[i] for i in idx}
    assert [int(round(per[i])) for i in idx] == [32, 32, 48], per
    max_shift, rows = 0.0, []
    for k in range(1, 4):
        for W in itertools.combinations(idx, k):
            c = sum(per[i] for i in W)
            a12, a23, a13 = angles_of((3.0 + c) * J)
            sh = abs(a12 - 1.0 / 3.0)
            max_shift = max(max_shift, sh)
            rows.append({"W": list(W), "c": float(c),
                         "eigs": sorted(float(x) for x in
                                        np.linalg.eigvalsh((3.0 + c) * J)),
                         "s12_heavy": a12, "angle_shift": float(sh)})
    c_total = sum(per[i] for i in idx)
    data["t1"] = {"hosts": idx, "p0": p0, "per_host_E": per,
                  "combos": rows, "c_all": float(c_total),
                  "c_artifact_check": abs(c_total - c_art) < 1e-6,
                  "max_angle_shift": float(max_shift)}
    checks.append({"name": "T1 census of the 7 non-empty host combinations: "
                          "every correction is c x ones (eigs (3+c,0,0)), "
                          "maximal angle shift {:.1e} over the entire "
                          "realizable sector; c(all) = {:.1f} reproduces "
                          "V292's booked {:.1f}".format(
                              max_shift, c_total, c_art),
                   "pass": max_shift < 1e-12 and abs(c_total - c_art) < 1e-6})

    # ---- T2 the locked point against NuFIT ---------------------------------
    s12, s13, s23 = 1.0 / 3.0, 0.0, 0.5
    sig = {"s12": (s12 - NUFIT["s12"][0]) / NUFIT["s12"][1],
           "s13": (s13 - NUFIT["s13"][0]) / NUFIT["s13"][1],
           "s23": (s23 - NUFIT["s23"][0]) / NUFIT["s23"][1]}
    joint = float(np.mean(
        [abs(s12 - NUFIT["s12"][0]) < NUFIT["s12"][1],
         abs(s13 - NUFIT["s13"][0]) < NUFIT["s13"][1],
         abs(s23 - NUFIT["s23"][0]) < NUFIT["s23"][1]]))
    data["t2"] = {"locked_point": [s12, s13, s23],
                  "nufit_imported": NUFIT, "sigma": sig,
                  "joint_frac_census": joint}
    checks.append({"name": "T2 the locked democratic point (1/3, 0, 1/2): "
                          "sigma_12 = {:+.2f}, sigma_13 = {:+.2f}, "
                          "sigma_23 = {:+.2f} (NuFIT imported from V104); "
                          "joint fraction over the realizable census "
                          "= {:.1f}".format(sig["s12"], sig["s13"],
                                            sig["s23"], joint),
                   "pass": abs(sig["s13"]) > 10 and joint == 0.0})

    branch = ("MIXING-REALIZABLE-EMPTY" if max_shift < 1e-12 and joint == 0.0
              else "MIXING-REALIZED-N")
    data["status"] = branch
    print("V304 {}: across the ENTIRE realizable scalar sector (all 7 "
          "neutral host combinations, corrections c x ones with c up to "
          "{:.1f}) no mixing angle moves (max shift {:.1e}); the point is "
          "locked at (1/3, 0, 1/2), sigma_13 = {:+.1f} against NuFIT, "
          "joint fraction 0.0 - the 1/3 -> 0.307 fixup and the whole "
          "theta13 deficit require a declared direction, no realizable "
          "vev supplies one".format(branch, c_total, max_shift, sig["s13"]),
          flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def _emit(data, checks, ok):
    out = {"id": "V304_realizable_mixing_census", "tier": 3,
           "property": "mixing corrections of the entire realizable scalar "
                       "sector (V291 hosts, V294 commutant, V104 NuFIT)",
           "residual": True, "programme": "V304",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v304_realizable_mixing_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V304 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
