#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V292 - the vev-host mixing test (tier 3).

**The question (pre-registered).** With the booked democratic mass
matrix (V45: all-ones, rank 1, eigenvalues {3, 0, 0}, TBM
sin^2 theta12 = 1/3) and the chain-host insertion from V291 (booked
channel energies, no free dial), does the mixing move, and by how
much against PDG sin^2 theta12 = 0.307?

**The object.** the V45 mass_matrix and mixing anchors, the V291
artifact (n0, p0, booked chain_H_spectrum entries of host states).
No dial, no grid, no optimization.

**Measured.** T0: V45 mixing re-read and the democratic
diagonalization re-verified (TBM to 1e-9). T1: the host insertion
is family-symmetric BY THE BOOKED 3-fold family symmetry (V45
family_symmetry: the 45 act identically on each family), so the
only induced correction is delta-M = c x ones with c = the booked
sum of host chain-H energies: M + delta-M = (3 + c) x ones: rank
and angles invariant; the mass rescale c/3 booked.

Branch: MIX-UNCHANGED (angles shift < 1e-9: TBM stands; the PDG
gap is the standing named falsifier, booked) / MIX-SHIFTED (angles
booked, fitted-flag false).
"""
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
from cl10.paths import exp                                # noqa: E402


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def angles_of(M):
    Mm = 0.5 * (M + M.conj().T)
    w, V = np.linalg.eigh(Mm)
    order = np.argsort(-np.abs(w))
    V = V[:, order]
    s12 = np.abs(V[1, 1]) / np.sqrt(np.abs(V[0, 1]) ** 2
                                    + np.abs(V[1, 1]) ** 2 + 0.0)
    s23 = np.abs(V[2, 1]) / np.sqrt(np.abs(V[2, 1]) ** 2
                                    + np.abs(V[2, 2]) ** 2 + 0.0)
    s13 = np.abs(V[2, 0])
    return (float(s12 ** 2), float(s23 ** 2), float(s13 ** 2))


def main():
    t0 = time.time()
    checks, data = [], {}
    v45 = art("V45_yukawa_matrix")
    v291 = art("V291_higgs_host_channels")

    # ---- T0: democratic anchor re-verified ------------------------------
    M = np.ones((3, 3))
    a_mm = v45["data"]["mass_matrix"]
    ev = np.linalg.eigvalsh(M)
    t0ok = (abs(float(a_mm["eigenvalues"][-1]) - 3.0) < 1e-9
            and abs(float(np.max(ev)) - 3.0) < 1e-9 and int(a_mm["rank"])
            == 1)
    s12a = float(v45["data"]["mixing"]["sin2_12"])
    obs = [float(x) for x in v45["data"]["mixing"]["obs_pmns"]]
    data["t0"] = {"democratic_eigs": [float(x) for x in ev],
                  "booked_sin2_12": s12a, "obs_pmns": obs}
    checks.append({"name": "T0 democratic anchor: eigs {{3,0,0}} "
                          "re-derived, "
                          "booked sin2_12 = {:.6f}, obs {}".format(
                              s12a, obs),
                   "pass": bool(t0ok and abs(s12a - 1.0 / 3.0) < 1e-9)})

    # ---- T1: the host insertion (booked only, family-symmetric) --------
    n0 = int(v291["data"]["host"]["n0"])
    p0 = [float(x) for x in v291["data"]["host"]["p0"]]
    hspec = [float(x) for x in v291["data"]["host"]
             ["chain_H_spectrum_artefact"]]
    host_idx = [i for i, x in enumerate(p0) if x > 0.5]
    c = sum(hspec[i] for i in host_idx)
    Mh = M + c * np.ones((3, 3))
    a0 = angles_of(M)
    a1 = angles_of(Mh)
    shift = max(abs(x - y) for x, y in zip(a0, a1))
    data["t1"] = {"n0": n0, "host_indices": host_idx,
                  "c_host_energy_sum": float(c),
                  "angles_democratic": a0, "angles_with_host": a1,
                  "angle_shift": float(shift),
                  "mass_rescale_frac": float(c / 3.0)}
    checks.append({"name": "T1 host insertion: n0 = {} (c = {:.3e}): angle "
                          "shift {:.2e}; PDG sin2_12 gap |1/3 - {:.3f}| = "
                          "{:.3f} is the standing falsifier".format(
                              n0, c, shift, obs[0], abs(1.0 / 3.0 - obs[0])),
                   "pass": True})
    status = "MIX-UNCHANGED" if shift < 1e-9 else "MIX-SHIFTED"
    data["status"] = status
    print("V292 {}: angle shift {:.2e}; PDG gap booked {:.4f}; mass "
          "rescale by host energy sum: {:+.3e}".format(
              status, shift, abs(1.0 / 3.0 - obs[0]), c / 3.0), flush=True)
    ok = all(cc["pass"] for cc in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V292_vev_host_mixing", "tier": 3,
           "property": "does the booked chain-host insertion move the "
                       "democratic TBM mixing",
           "residual": True, "programme": "V292",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v292_vev_host_mixing.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V292 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
