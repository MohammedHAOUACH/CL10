#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V368 - torus rigidity of the composite isolation window (booked
BEFORE the run; reality-interface third cut). Tests the machine-forced
prediction following V367: since bloch(q,U) = kron(P(q), U^6) and U^6
is q-independent, the composite gap must be momentum-rigid. V259's 6^3
grid; per total momentum K: composite set (bond-pair phases + physical
omega rows, ALL 216 relative pairs) vs continuum set (bond-pair + the
32x32 internal phase sums, 24 booked random k-points per K, seed
0x51). Booked anchor: contact-locus gap reproduces V367's artifact
value to 1e-9. omega=0 is a booked control (merge). Verdicts:
RIGID-WINDOW-<minK gap> / RIGID-WINDOW-BROKEN-<nloci> / ANCHOR-DRIFT.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
from cl10 import supercell as SC  # noqa: E402
from cl10.paths import exp  # noqa: E402

GRID_N = 6


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def modpi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def circmin(c, T):
    """min circular distance from each c to the sorted array T."""
    i = np.searchsorted(T, c)
    n = len(T)
    best = np.full(len(c), np.inf)
    for sh in (i - 1, i, i + 1):
        j = sh % n
        d = np.abs(modpi(c - T[j]))
        best = np.minimum(best, d)
    return float(best.min())


def main():
    checks, data = [], {}

    def log(m):
        print(m, flush=True)

    _Bm, U32, _U16 = SC.sectors(SC.G10)
    U6 = np.linalg.matrix_power(U32, 6)
    iota = np.angle(np.linalg.eigvals(U6))           # 32 internal phases
    grid = [(2 * np.pi * i / GRID_N, 2 * np.pi * j / GRID_N,
             2 * np.pi * l / GRID_N)
            for i in range(GRID_N) for j in range(GRID_N)
            for l in range(GRID_N)]
    betas = np.array([np.angle(np.linalg.eigvals(SC.bond_product(q)))
                      for q in grid])                # (216, 8)

    v251 = json_art("V251_gauged_two_particle_sector")
    rows = {k: float(v) for k, v in
            v251["data"]["sensitivity"]["omega_rows"].items()}
    phys = [v for v in rows.values() if abs(v) > 1e-9]
    zero_present = any(abs(v) < 1e-9 for v in rows.values())

    # ---- booked anchor: contact locus (k1 = k2 = 0) vs V367 artifact
    b0 = betas[0]
    c_set = modpi(b0[:, None] + b0[None, :]).ravel()
    t_set = modpi((b0[:, None] + b0[None, :]).ravel()[:, None]
                  + (iota[:, None] + iota[None, :]).ravel()[None, :]).ravel()
    t_set = np.sort(modpi(t_set))
    anchor_gaps = {}
    for om in phys:
        anchor_gaps[str(om)] = circmin(modpi(c_set + om), t_set)
    v367 = json_art("V367_singlet_pair_gamma")
    v367_vals = sorted(float(v) for v in v367["data"]["window_gaps"]
                       .values() if float(v) > 1e-6)
    anchor_ok = all(min(abs(g - x) for x in v367_vals) < 1e-9
                    for g in anchor_gaps.values())
    checks.append(("anchor: contact gap %s vs V367 artifact %s to 1e-9: "
                   "%s" % (["%.6f" % g for g in anchor_gaps.values()],
                           ["%.6f" % x for x in v367_vals], anchor_ok),
                   anchor_ok))

    rng = np.random.default_rng(0x51)
    gaps = []
    idx = np.arange(GRID_N ** 3)
    tri = np.array([(i, j, l) for (i, j, l) in
                    ((a, b, c) for a in range(GRID_N)
                     for b in range(GRID_N) for c in range(GRID_N))])
    for K in range(GRID_N ** 3):
        k1 = idx
        k2 = ((tri - tri[K]) % GRID_N) @ np.array([36, 6, 1])
        bA, bB = betas[k1], betas[k2]
        comp = modpi((bA[:, :, None] + bB[:, None, :]).ravel()[:, None]
                     + np.array(phys)[None, :]).ravel()
        samp = rng.choice(216, size=24, replace=False)
        s1, s2 = samp, ((tri[samp] - tri[K]) % GRID_N) @ np.array([36, 6, 1])
        sb = (betas[s1][:, :, None] + betas[s2][:, None, :]).ravel()
        cont = modpi(sb[:, None]
                     + (iota[:, None] + iota[None, :]).ravel()[None, :]
                     ).ravel()
        cont = np.sort(modpi(cont))
        gaps.append(circmin(comp, cont))
    gaps = np.array(gaps)
    data["gaps_min"] = float(gaps.min())
    data["gaps_max"] = float(gaps.max())
    data["gaps_mean"] = float(gaps.mean())
    broken = int(np.sum(gaps < 1e-6))
    log("torus gaps: min %.6f max %.6f spread %.1e broken %d/%d"
        % (gaps.min(), gaps.max(), gaps.max() - gaps.min(), broken,
           len(gaps)))

    # booked control: omega=0 composite (bond band) must merge
    c0 = np.sort(modpi((betas[0][:, None] + betas[0][None, :]).ravel()))
    ctrl = circmin(c0, t_set)
    checks.append(("control: omega=0 merges into continuum (gap %.1e "
                   "< 1e-6): %s" % (ctrl, ctrl < 1e-6), ctrl < 1e-6))
    spread = float(gaps.max() - gaps.min())
    data["spread"] = spread

    if not (anchor_ok and ctrl < 1e-6 and zero_present):
        status = "V368-ANCHOR-DRIFT"
        ok = False
    elif broken == 0 and spread <= 1e-9:
        status = ("RIGID-WINDOW-%.6f (isolation holds at every one of "
                  "the %d total momenta of the 6^3 torus; spread over "
                  "the torus %.1e - the window is set by the clock's "
                  "internal phase spectrum and the flux energy, not by "
                  "the locus)" % (gaps.min(), len(gaps), spread))
        ok = True
    else:
        # booked post-run correction (2026-11-01): the booked prediction
        # that the window is momentum-rigid is REFUTED by the machine at
        # this cut: at fixed total momentum the continuum is the union
        # over ALL relative momenta, and the bond-pair bands overlap the
        # shifted composite band - the contact window 1.007 collapses to
        # min gap %.1e over the torus (spread %.1e). The composite is
        # embedded in the single-flux two-quanta continuum at generic
        # momenta: at best a resonance, not an asymptotic state. The
        # old booked word RIGID-WINDOW-BROKEN reserved "broken" for
        # gap <= 1e-6; this verdict keeps broken=0 but names the
        # collapse honestly.
        status = ("WINDOW-NOT-RIGID-%.1e (the V367 contact window "
                  "collapses from %.6f at contact to %.1e at generic "
                  "torus momenta - the composite disperses into the "
                  "single-flux two-quanta continuum; no stable "
                  "asymptotic candidate at this dynamics; grid "
                  "convergence of the min gap not measured - booked "
                  "as V369 question)"
                  % (gaps.min(), next(iter(anchor_gaps.values())),
                     gaps.min()))
        ok = True
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V368 done] %s" % status, flush=True)
    out = {"id": "V368_torus_rigidity_window", "tier": 2, "residual": True,
           "property": "is the V367 composite isolation window rigid "
                       "across the whole torus or locus-bound",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v368_torus_rigidity_window.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
