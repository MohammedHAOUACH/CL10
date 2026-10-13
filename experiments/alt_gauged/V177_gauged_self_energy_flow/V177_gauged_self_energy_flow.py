#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V177 - one-loop flow of the cone anisotropy from the GAUGED vertex
(alt-gauged, step 3)

**The question.**  V172 ran the one-loop Fock flow with the contact
vertex: MARGINAL (alpha = 1), not grid-convergent (the flat axis, V118).
V176 showed the gauged bond-current vertex is NOT a scalar: it carries
Eg + T2g (spin-2) weight on the cone.  Does the one-loop self-energy of
the GAUGED vertex (V79/V101 bond current, coupling eta = g_10 = 1/5 - the
only admissible strength) drive the W9 anisotropy toward isotropy?  Same
protocol as V172 (static Fock, normal-ordered tadpole subtraction, all 12
M-node singularities excluded and ray-subtracted), vertex replaced by
``current_vertex=True`` in ``altcell.fock_self_energy``.

**The result.**  Recorded: A_dressed(dq) at three scales for two grids
(16, 24 - the grid drift is a RESULT, V172 precedent), the cubic-harmonic
content and |k| exponent of the gauged self-energy, and the flow dA.
"""

import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
UC = SC.G10                      # eta = g_10: the only admissible strength
DQS = (1e-2, 1e-3, 1e-4)
QSTAR = (math.pi / 2, math.pi / 2, 0.0)
GRIDS = (16, 24)


def bf(q):
    return A.bloch_cell(C, SC.sectors()[2], SC.G10, ST, q)


def _fit(ks, vals):
    lg = [(math.log(k), math.log(max(v, 1e-16)))
          for k, v in zip(ks, vals) if v > 1e-13]
    if len(lg) < 2:
        return None
    return float(np.polyfit([a for a, _ in lg], [b for _, b in lg], 1)[0])


def dressed_A(grids, gn, dirs, D26, ncone, ncub, blk):
    """A_dressed(dq) and the Sigma harmonic content on grid gn."""
    G = grids[gn]
    Ad, wts, aniso = {}, {}, {}
    for di, dq in enumerate(DQS):
        vs = []
        for i, d in enumerate(dirs):
            q = tuple(a + dq * b for a, b in zip(QSTAR, d))
            Heff, _ = A.lowdin_downfold(bf, q, QSTAR, n_low=1)
            S = UC * (G[1 + di * blk + i] - G[0])
            ev = np.linalg.eigvalsh(Heff + S)
            vs.append(float(np.min(np.abs(ev))) / dq)
        Ad[dq] = (max(vs) - min(vs)) / (max(vs) + min(vs))
        tr = np.array([float(np.real(
            (UC * (G[1 + di * blk + ncone + j] - G[0]))[0, 0]))
                       for j in range(ncub)])
        wts[dq] = A.cubic_harmonic_weights(tr, D26)
        aniso[dq] = float(np.abs(tr - tr.mean()).max())
    return Ad, wts, aniso


def main():
    t0 = time.time()
    print("=" * 68)
    print("V177 - one-loop flow with the gauged bond-current vertex")
    print("=" * 68)
    dirs = A.sphere_dirs(6)
    D26 = A.oh26_dirs()
    ncone, ncub = len(dirs), len(D26)
    blk = ncone + ncub
    kappas = [(0.0, 0.0, 0.0)]
    for dq in DQS:
        for d in dirs:
            kappas.append(tuple(dq * x for x in d))
        for d in D26:
            kappas.append(tuple(dq * x for x in d))

    grids = {gn: A.fock_self_energy(bf, QSTAR, 1.0, gn, kappas,
                                    sub=(0.15, 0.6), nodes=SC.M_POINTS,
                                    current_vertex=True)
             for gn in GRIDS}
    print(f"  self-energy grids {GRIDS} done [{time.time() - t0:.0f}s]")

    # discriminating: Sigma not identically zero (vertex reaches the cone)
    nz = max(float(np.abs(grids[gn][1] - grids[gn][0]).max())
             for gn in GRIDS)
    print(f"  discriminating: max |Sigma(kappa) - tadpole| = {nz:.3e}")

    # bare downfolded cone speeds
    bare = {}
    for dq in DQS:
        for i, d in enumerate(dirs):
            q = tuple(a + dq * b for a, b in zip(QSTAR, d))
            _, ev = A.lowdin_downfold(bf, q, QSTAR, n_low=1)
            bare[(dq, i)] = float(np.min(np.abs(ev))) / dq
    A_bare = {}
    for dq in DQS:
        vs = [bare[(dq, i)] for i in range(ncone)]
        A_bare[dq] = (max(vs) - min(vs)) / (max(vs) + min(vs))

    rows = []
    for gn in GRIDS:
        Ad, wts, aniso = dressed_A(grids, gn, dirs, D26, ncone, ncub, blk)
        rows.append({
            "grid": gn, "U_c": UC,
            "A_dressed": {str(dq): Ad[dq] for dq in DQS},
            "A_bare": {str(dq): A_bare[dq] for dq in DQS},
            "dA_vs_free": Ad[1e-4] - A_bare[1e-4],
            "w_A1g_Sigma": {str(dq): wts[dq]["A1g"] for dq in DQS},
            "w_aniso_Sigma": {str(dq): wts[dq]["Eg"] + wts[dq]["T2g"]
                              for dq in DQS},
            "alpha_aniso_Sigma": _fit(list(DQS), [aniso[dq]
                                                  for dq in DQS]),
            "isotropic": bool(min(Ad.values()) < 1e-6),
        })
        print(f"  grid {gn}: A_dressed = "
              + "  ".join(f"{Ad[dq]:.4f}" for dq in DQS)
              + f"   alpha_aniso(Sigma) = {rows[-1]['alpha_aniso_Sigma']:.2f}")

    # 2-route: grid drift on A_dressed(1e-3)
    devs = [abs(rows[0]["A_dressed"]["0.001"] - rows[1]["A_dressed"]["0.001"])]

    anchor_ok = abs(A_bare[1e-3] - 0.998) < 1e-3
    ok = anchor_ok and len(rows) == 2 and nz > 1e-12

    r24 = rows[1]
    verdict = {
        "id": "V177_gauged_self_energy_flow",
        "tier": 4,
        "property": "static Fock self-energy of the GAUGED bond-current "
                    "vertex (Gamma^A = I_8 (x) T^A_16, U_c = g_10) on the "
                    "M-node cone: dressed anisotropy A(dq), cubic-harmonic "
                    "content and |k| scaling of Sigma, grids 16 vs 24",
        "residual": True,
        "pass": bool(ok),
        "status": f"ONE-LOOP GAUGED-VERTEX FLOW (2 rows, grids 16/24, "
                  f"normal-ordered, 12 M-node singularities ray-subtracted "
                  f"- V172 protocol, vertex replaced by the Spin(10) "
                  f"current). Sigma is NOT zero (max |Sigma - tadpole| = "
                  f"{nz:.2e}): the gauged vertex reaches the cone. Its "
                  f"anisotropic part scales as |k|^"
                  f"{r24['alpha_aniso_Sigma']:.2f} (V172 contact: 1.00 - "
                  f"MARGINAL) with aniso weight "
                  f"{r24['w_aniso_Sigma']['0.0001']:.4f} at dq=1e-4. "
                  f"A_dressed(1e-4) = {r24['A_dressed']['0.0001']:.4f} vs "
                  f"bare {r24['A_bare']['0.0001']:.4f}: flow dA = "
                  f"{r24['dA_vs_free']:+.4f}. Grid drift 16 vs 24 on "
                  f"A_dressed(1e-3) = {devs[0]:.1e} (IR-sensitive, V172 "
                  f"flat-axis precedent - drift recorded, not a failure).",
        "data": {"rows": rows, "node": list(QSTAR), "dqs": list(DQS),
                 "route_dev": devs, "sigma_nonzero": nz,
                 "candidates": [r for r in rows if r["isotropic"]]},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v177_gauged_self_energy_flow.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  bare A(1e-3) = {A_bare[1e-3]:.4f} (V63/V115 0.998)")
    print(f"  2-route grid drift: {devs[0]:.1e}   isotropic rows: "
          f"{sum(r['isotropic'] for r in rows)}   "
          f"VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir,
                                      "v177_gauged_self_energy_flow.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
