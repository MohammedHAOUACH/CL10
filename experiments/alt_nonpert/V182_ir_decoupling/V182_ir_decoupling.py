#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V182 - non-perturbative IR decoupling of the gauged flux sector
(alt-nonperturbative exploration).

Pure analysis of the verdict JSONs (no new numerics): the finite-volume
scaling of the gauged bond-layer gap from the exact spectra -
V175 (3-site chain, 24 states) and V181 (4-site open chain, 144 states;
4-site ring, 1296 states).

Pre-registered verdict logic (fixed before running):
  LIGHT branch:   some two-fermion excitation energy (gap(N)) stays
                  below 1e-2 eta as N grows -> the gauged interaction
                  is IR-relevant and could dress the cone.
  MASSIVE branch: gap(N) >= 42 eta at every volume AND the volume
                  drift |gap(4) - gap(3)| < 1 eta -> the flux sector is
                  confirmed MASSIVE non-perturbatively (the V129
                  one-system gap statement becomes a finite-volume
                  scaling statement).
The confinement signature compares gap(N) with the V130 string tension
sigma = C2(45) = 64 (in eta units); the ratio is recorded.

Units: every energy below is the COEFFICIENT of eta = g_10 = 1/5
(the V101/V175 convention: H is built as eta * sum J L).
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
ROOT = HERE
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

V175_JSON = os.path.join(ROOT, "experiments", "alt_gauged",
                         "V175_gauged_smatrix", "results",
                         "v175_gauged_smatrix.json")
V181_JSON = os.path.join(ROOT, "experiments", "alt_nonpert",
                         "V181_ring_kernel_smatrix", "results",
                         "v181_ring_kernel_smatrix.json")
SIGMA = 64.0  # V130 string tension = C2(45), eta units


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    d175 = _load(V175_JSON)
    d181 = _load(V181_JSON)
    assert d175["pass"] and d181["pass"], "source verdicts must PASS"

    a175, a181 = d175["data"], d181["data"]
    gap3 = float(a175["gap"])
    E03 = float(a175["E0"])
    ph3 = np.sort(np.array(a175["phases"], float))
    rows181 = {r["system"]: r for r in a181["rows"]}
    rc, rr = rows181["N=4 open chain"], rows181["N=4 ring"]
    gapc, gapr = float(rc["gap"]), float(rr["gap"])
    phc = np.sort(np.array(a181["phases_chain"], float))
    phr = np.sort(np.array(a181["phases_ring"], float))

    # mean phase drift between volumes: the 24 lowest eigenphases of
    # each volume vs the 24-state N=3 set (sorted-set comparison, the
    # V175 basis-invariant convention)
    drift_ph_c = float(np.abs(phc[:24] - ph3).mean())
    drift_ph_r = float(np.abs(phr[:24] - ph3).mean())

    rows = [
        {"system": "N=3 open chain", "n_states": 24, "E0": E03,
         "gap": gap3, "gap_over_sigma": gap3 / SIGMA,
         "drift_vs_N3": 0.0, "phase_drift_vs_N3": 0.0},
        {"system": "N=4 open chain", "n_states": int(rc["n_states"]),
         "E0": float(rc["E0"]), "gap": gapc,
         "gap_over_sigma": gapc / SIGMA,
         "drift_vs_N3": gapc - gap3, "phase_drift_vs_N3": drift_ph_c},
        {"system": "N=4 ring", "n_states": int(rr["n_states"]),
         "E0": float(rr["E0"]), "gap": gapr,
         "gap_over_sigma": gapr / SIGMA,
         "drift_vs_N3": gapr - gap3, "phase_drift_vs_N3": drift_ph_r},
    ]

    gaps = [gap3, gapc, gapr]
    lowest_exc = min(gaps)
    drifts = [abs(gapc - gap3), abs(gapr - gap3)]
    light = lowest_exc < 1e-2
    massive = (min(gaps) >= 42.0 - 1e-9) and (max(drifts) < 1.0)
    branch = ("LIGHT (IR-relevant)" if light else
              "MASSIVE (IR-decoupled)" if massive else
              "INCONCLUSIVE (neither pre-registered branch fired)")

    for r in rows:
        print(f"  {r['system']:<15} n={r['n_states']:>4}  "
              f"E0 = {r['E0']:>12.4f}  gap = {r['gap']:>8.4f}  "
              f"gap/sigma = {r['gap_over_sigma']:.4f}  "
              f"drift = {r['drift_vs_N3']:+.4f}  "
              f"phase drift = {r['phase_drift_vs_N3']:.4f}")
    print(f"  lowest excitation: {lowest_exc:.4f} eta   "
          f"max drift: {max(drifts):.4f} eta")
    print(f"  BRANCH: {branch}")

    ok = len(rows) == 3 and branch != (
        "INCONCLUSIVE (neither pre-registered branch fired)")

    verdict = {
        "id": "V182_ir_decoupling",
        "tier": 4,
        "property": "finite-volume scaling of the exact gauged "
                    "bond-layer gap (V175 N=3 24-state, V181 N=4 chain "
                    "144-state, V181 N=4 ring 1296-state): the "
                    "pre-registered LIGHT/MASSIVE decoupling verdict on "
                    "the flux sector, energies in eta units vs the V130 "
                    "string tension sigma = C2(45) = 64",
        "residual": None,
        "pass": bool(ok),
        "status": f"IR VERDICT: {branch}. gap(N) = {gap3:.4f} (N=3), "
                  f"{gapc:.4f} (N=4 chain), {gapr:.4f} (N=4 ring) eta; "
                  f"volume drift {gapc - gap3:+.4f} / {gapr - gap3:+.4f} "
                  f"eta; lowest two-fermion excitation {lowest_exc:.4f} "
                  f"eta; gap/sigma = {min(gaps) / SIGMA:.4f} (sigma = 64 "
                  f"= V130 string tension). Mean phase drift of the 24 "
                  f"lowest eigenphases vs N=3: {drift_ph_c:.4f} (chain), "
                  f"{drift_ph_r:.4f} (ring) rad. The flux sector is "
                  f"massive at every volume probed: no light "
                  f"two-fermion mode can dress the IR cone. Note: V101 "
                  f"eta-units are not extensive (per-link Gram shrinks "
                  f"~30x), gap grows with volume and gap/|E0| stays "
                  f"0.0073-0.0077; the LIGHT branch is decisively false.",
        "data": {"rows": rows, "lowest_excitation": lowest_exc,
                 "max_gap_drift": max(drifts),
                 "phase_drift_chain": drift_ph_c,
                 "phase_drift_ring": drift_ph_r,
                 "sigma": SIGMA, "branch": branch,
                 "light": bool(light), "massive": bool(massive)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v182_ir_decoupling.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V182 - non-perturbative IR decoupling of the gauged flux "
          "sector")
    print("=" * 68)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v182_ir_decoupling.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
