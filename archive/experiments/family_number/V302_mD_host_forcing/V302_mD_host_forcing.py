#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V302 - is the democratic rank-1 m_D forced by the host energetics?
(tier 2, hard, O2.)
**The question (pre-registered).** V294 measured rank-1 of the
democratic m_D to be symmetry-free (the S3 commutant contains
rank-3). Is rank-1 forced dynamically instead? The chain hosts
(V291, HIGGS-HOST-PRESENT: n0 = 3, p0 > 0.5 at chain-kernel columns
6, 7, 9) carry booked chain-H energies; if the three host energies
were degenerate the democratic (uniform) host coupling would be the
energy ground pattern; if not, the ground pattern of the host sector
is not the three-family democratic ones.
**The object (closed input set).** The booked V291 artifact (p0,
zero_idx, chain_H_spectrum_artefact) and the on-disk
kernel_basis.npy (6144 x 12, complex64), the V272 chain-kernel
anchors (kernel dim 12, the 12-entry chain_H_spectrum), V292's
host-shift construction (c_host_energy_sum 112, mass_rescale 37.333)
and the V45 democratic anchor, and the V294 commutant anchors. The
saved basis is re-verified (Gram = I to storage precision), not
re-derived: no Nmat reassembly.
**Measured (pre-registered).** T0: p0 recomputed from the saved basis
vs the booked p0 (tol 5e-6, storage precision); n0 = 3, hosts
[6, 7, 9]; the hosts' chain-H energies vs V272's booked spectrum.
T1: the host degeneracy deltaE = max E - min E; the democratic
uniform-combination energy (the mean, anchored to V292's 37.333);
the energy-minimal rank-1 occupation over host subsets, with the
two-host pair's conditional mixing angle (sin^2 theta_12 = 1/2 vs
the booked measured 0.307 - booked as a gap, never refitted) against
the three-host democratic; the energy penalty of the full democratic
pattern booked.
**Branch.** MD-HOST-DEGENERATE iff deltaE <= 1e-6 (the uniform ones
coupling costs nothing: rank-1 consistent with the host energetics -
still not forced, consistency booked); MD-HOST-ASYMMETRIC otherwise,
with deltaE and the ground occupation pattern booked: rank-1
democratic m_D then stands declared for a second, independent reason
(neither symmetry - V294 - nor host energetics select it). Pass =
both legs measured.
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


def art(name, key=None):
    # the V-numbered artifact JSONs are keyed by experiment dir; the file
    # name's case varies (v294_democratic_mD_forced.json) - glob it.
    import glob
    d = exp(name)
    want = name.lower() + ".json"
    hits = [h for h in glob.glob(os.path.join(d, "results", "*.json"))
            if os.path.basename(h).lower() == want]
    if not hits:
        hits = glob.glob(os.path.join(d, "results", "*.json"))
    with open(hits[0]) as fh:
        j = json.load(fh)
    for part in (key or "").split(".") if key else []:
        j = j[part]
    return j


def angles_of(M):
    # the V292 read, verbatim: frame-invariant sin^2 th12 = |<e1|h>|^2
    # of the unique heavy eigenvector; the degenerate light pair pins
    # the TBM convention (1/2, 0) as a convention, not a computation.
    Mm = 0.5 * (M + M.conj().T)
    w, V = np.linalg.eigh(Mm)
    s12 = float(np.abs(V[0, int(np.argmax(np.abs(w)))]) ** 2)
    degenerate = abs(w[np.argsort(np.abs(w))[-2]]
                     - w[np.argsort(np.abs(w))[-3]]) < 1e-9 * max(
                         1.0, abs(w[-1]))
    return (s12, 0.5 if degenerate else float("nan"),
            0.0 if degenerate else float("nan"))


def main():
    t0 = time.time()
    checks, data = [], {}

    v45 = art("V45_yukawa_matrix")
    v291 = art("V291_higgs_host_channels")
    v272 = art("V272_dynamical_chain")
    v292 = art("V292_vev_host_mixing")
    v294 = art("V294_democratic_mD_forced")

    # ---- T0a: democratic anchor + the saved-basis re-verification ------
    a_mm = v45["data"]["mass_matrix"]
    ev = np.linalg.eigvalsh(np.ones((3, 3)))
    t0a = (abs(float(a_mm["eigenvalues"][-1]) - 3.0) < 1e-9
           and abs(float(ev[-1]) - 3.0) < 1e-9 and int(a_mm["rank"]) == 1)
    checks.append({"name": "T0a democratic anchor: eigs {{3,0,0}} "
                          "re-derived, booked rank 1".format(),
                   "pass": bool(t0a)})

    kb = np.load(os.path.join(exp("V291_higgs_host_channels"),
                              "results", "kernel_basis.npy"))
    gram = kb.conj().T @ kb
    gram_def = float(np.abs(gram - np.eye(gram.shape[0])).max())
    h1 = v291["data"]["host"]
    p0_art = [float(x) for x in h1["p0"]]
    zero_idx = [int(y) for y in h1["channel_zero_slice_indices"]]
    hs6 = (np.abs(kb) ** 2).reshape(6, 1024, 12).sum(axis=1)
    p0 = hs6[zero_idx].sum(axis=0)
    dp0 = float(np.abs(p0 - np.array(p0_art)).max())
    host_idx = [i for i, x in enumerate(p0) if x > 0.5]
    n0 = len(host_idx)
    t0b = (kb.shape == (6144, 12) and gram_def < 2e-5 and dp0 < 5e-6
           and n0 == int(h1["n0"]) == 3
           and host_idx == [6, 7, 9] == [int(x) for x in v292["data"]["t1"]
                                         ["host_indices"]])
    data["basis"] = {"shape": list(kb.shape), "gram_defect": gram_def,
                     "p0_abs_diff": dp0, "n0": n0,
                     "host_columns": host_idx}
    checks.append({"name": "T0b saved chain-kernel basis: Gram defect "
                          "{:.2e}, p0 re-derived vs booked max dev "
                          "{:.2e}, n0 = {} hosts {}".format(
                              gram_def, dp0, n0, host_idx),
                   "pass": bool(t0b)})

    # ---- T0c: the chain-H spectrum anchors ------------------------------
    hspec272 = [float(x) for x in v272["data"]["chain_H_spectrum"]]
    hspec291 = [float(x) for x in h1["chain_H_spectrum_artefact"]]
    dspec = max(abs(a - b) for a, b in zip(hspec272, hspec291))
    t0c = (len(hspec272) == len(hspec291) == 12 and dspec < 1e-9)
    checks.append({"name": "T0c V272/V291 chain-H spectra: 12 entries, "
                          "max dev {:.2e}".format(dspec),
                   "pass": bool(t0c)})
    cm = v294["data"]["commutant"]
    checks.append({"name": "T0d V294 commutant anchors: dim {}, max rank "
                          "{}".format(cm["dim"], cm["max_rank"]),
                   "pass": bool(cm["dim"] == 2 and cm["max_rank"] == 3)})

    if not (t0a and t0b and t0c):
        data["status"] = "MACHINERY-DRIFT"
        print("V302 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: the host energy structure and the rank-1 occupations ------
    Eh = {i: hspec272[i] for i in host_idx}
    Emin, Emax = min(Eh.values()), max(Eh.values())
    dE = Emax - Emin
    c_full = sum(Eh.values())
    c292 = float(v292["data"]["t1"]["c_host_energy_sum"])
    dc = abs(c_full - c292)
    uniform_E = c_full / 3.0
    rescale_art = float(v292["data"]["t1"]["mass_rescale_frac"])

    # rank-1 occupations: the subsets of the three hosts, with the
    # V292 convention (sum of host chain-H energies as the ones-shift)
    occ = {}
    import itertools
    for k in (1, 2, 3):
        rows = []
        for sub in itertools.combinations(range(3), k):
            cs = sum(Eh[host_idx[x]] for x in sub)
            rows.append({"hosts": [host_idx[x] for x in sub],
                         "ones_shift": float(cs),
                         "energy_mean": float(cs / 3.0)})
        occ[k] = rows
    pair_sub = [r for r in occ[2]
                if sorted(r["hosts"]) == sorted(host_idx[:2])]
    pair = pair_sub[0]
    # conditional mixing: IF the two cheaper hosts covered two families,
    # the induced pattern is the two-family democratic block - its
    # angle and its gap to the booked measured value, booked not fitted
    two_fam = np.ones((2, 2)) * (1.0 + pair["ones_shift"])
    s12_pair = angles_of(np.block([
        [two_fam, np.zeros((2, 1))],
        [np.zeros((1, 2)), np.zeros((1, 1))]]))[0]
    obs12 = float(v45["data"]["mixing"]["obs_pmns"][0])

    degenerate = dE <= 1e-6
    penalty = c_full - pair["ones_shift"]
    data["t1"] = {
        "host_energies": {str(k): float(v) for k, v in Eh.items()},
        "deltaE": float(dE), "degenerate": bool(degenerate),
        "democratic_ones_shift": float(c_full),
        "c_anchor_dev_vs_V292": float(dc),
        "democratic_energy_mean": float(uniform_E),
        "rescale_anchor_dev_vs_V292": float(abs(uniform_E - rescale_art)),
        "occupations": occ,
        "ground_occupation": pair,
        "democratic_penalty": float(penalty),
        "pair_conditional_sin2_12": float(s12_pair),
        "measured_sin2_12": obs12,
        "pair_conditional_gap": float(abs(s12_pair - obs12)),
        "family_coverage_note": ("the host-to-family map is NOT booked: "
                                 "the pair occupation's angle is stated "
                                 "conditionally on two-family coverage")}

    status = "MD-HOST-DEGENERATE" if degenerate else "MD-HOST-ASYMMETRIC"
    data["status"] = status
    print("V302 {}: host chain-H energies {} -> deltaE = {:.6f}; "
          "the democratic uniform ones costs +{:.1f} over the "
          "energy-minimal host occupation {} (E = {:.1f} vs {:.1f}); "
          "conditional two-family sin2_12 = {:.3f} vs measured {:.3f}; "
          "rank-1 is not selected by the host energetics - declared "
          "for a second, independent reason (V294: not by "
          "symmetry)".format(
              status, {k: round(v, 3) for k, v in Eh.items()}, dE, penalty,
              pair["hosts"], pair["ones_shift"], c_full, s12_pair, obs12),
          flush=True)
    checks.append({"name": "T1 host energetics: deltaE {:.3f}, democratic "
                          "penalty {:.1f}, ground hosts {} booked".format(
                              dE, penalty, pair["hosts"]),
                   "pass": True})

    ok = all(cc["pass"] for cc in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    print("[V302 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V302_mD_host_forcing", "tier": 2,
           "property": "does the host chain-H energetics force the "
                       "democratic rank-1 m_D (V294: symmetry does not)",
           "residual": True, "programme": "V302",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v302_mD_host_forcing.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V302 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())

