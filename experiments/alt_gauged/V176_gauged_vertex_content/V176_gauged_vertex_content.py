#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V176 - the O_h content of the gauged vertex on the cone (alt-gauged,
step 2: the discriminating measurement)

**The question.**  V171 proved the minimal contact collision vertex is an
O_h SCALAR (w_A1g = 1.0000): it cannot rotate the W9 cone.  The gauged
bond-current interaction (V79/V101, S-matrix built in V175) has a
DIFFERENT vertex: the Spin(10) current ``Gamma^A = I_8 (x) T^A_16`` on the
128-D Bloch space (T^A_16 the Weyl generators of ``altcell.weyl_generators``
- anchored to exponentiate to the implemented step).  The adjoint 45
restricted to the cubic lattice carries the spin-2 pieces (V140 measured
sym^2 content A1g + Eg + T2g), so this vertex CAN carry quadrupole (l = 2)
structure.  Measure the node-mode current form factor

``F(kappa n) = sum_A |<u(q*)| Gamma^A |u(q* + kappa n)>|^2``

along the 26 O_h directions for |kappa| in {0.05, 0.1, 0.2, 0.4}, its
cubic-harmonic weights, and the low-kappa exponents.  This is the
discriminating comparison against V171's scalar contact.

**The result.**  Recorded: the weights (w_Eg + w_T2g nonzero or exactly
zero - a structural fact either way) and the exponents.  Anchors: the
Casimir sum ``sum_A Gamma^A Gamma^A = 45 I_128`` and the 32-D-module
recomputation of the form factor (route 2).
"""

import json
import math
import os
import sys
import time
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
os.environ["CL10_ANGLE"] = "0.2"  # the W9 cone and its vertex anchors are 1/5-branch (V115/V123) numbers: pin the clock, the closed default_branch angle must not leak in

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
QSTAR = (math.pi / 2, math.pi / 2, 0.0)
KS = (0.05, 0.1, 0.2, 0.4)


def bf(q):
    return A.bloch_cell(C, SC.sectors()[2], SC.G10, ST, q)


def _fit(ks, vals):
    lg = [(math.log(k), math.log(max(v, 1e-14)))
          for k, v in zip(ks, vals) if v > 1e-13]
    if len(lg) < 2:
        return None
    return float(np.polyfit([a for a, _ in lg], [b for _, b in lg], 1)[0])


def main():
    t0 = time.time()
    print("=" * 68)
    print("V176 - the gauged current vertex on the W9 cone")
    print("=" * 68)
    B, T16, f = A.weyl_generators()
    Gammas = [np.kron(np.eye(8), T) for T in T16]
    T32 = A.plane_generators_sig(0)
    Gammas32 = [np.kron(np.eye(8), T) for T in T32]

    # anchor: Casimir on the Bloch space
    cas = sum(G @ G for G in Gammas)
    err_cas = float(np.abs(cas - 45.0 * np.eye(128)).max())
    print(f"  anchor sum_A Gamma^A Gamma^A = 45 I_128: err {err_cas:.2e}")

    D26 = A.oh26_dirs()
    u0 = A.node_modes(bf, QSTAR, n=1)[0][:, 0]

    def F_at(q):
        uq = A.node_modes(bf, q, n=1)[0][:, 0]
        return float(sum(abs(v.conj() @ uq) ** 2 for v in
                         (G @ u0 for G in Gammas)))

    F0 = F_at(QSTAR)
    # route 2: 32-D module recomputation at 3 sampled (kappa, dir)
    v0 = np.kron(np.eye(8), B) @ u0
    devs = []
    for kap, di in ((0.0, 0), (0.2, 5), (0.4, 17)):
        q = tuple(a + kap * b for a, b in zip(QSTAR, D26[di]))
        uq = A.node_modes(bf, q, n=1)[0][:, 0]
        vq = np.kron(np.eye(8), B) @ uq
        F16 = sum(abs(v.conj() @ uq) ** 2 for v in (G @ u0 for G in Gammas))
        F32 = sum(abs(v.conj() @ vq) ** 2 for v in (G @ v0 for G in Gammas32))
        devs.append(abs(F16 - F32))
    print(f"  route 2 (32-D module) devs: {max(devs):.2e}")

    rows = []
    Fmean0 = None
    aniso_vals, iso_vals, kappas_used = [], [], []
    for kap in KS:
        F = np.array([F_at(tuple(a + kap * b
                                 for a, b in zip(QSTAR, d))) for d in D26])
        w = A.cubic_harmonic_weights(F, D26)
        aniso = float(np.abs(F - F.mean()).max())
        if Fmean0 is None:
            Fmean0 = F.mean()
        iso = float(abs(F.mean() - F0))
        aniso_vals.append(aniso)
        iso_vals.append(iso)
        kappas_used.append(kap)
        rows.append({
            "kappa": kap, "F_mean": float(F.mean()),
            "w_A1g": w["A1g"], "w_Eg": w["Eg"], "w_T2g": w["T2g"],
            "aniso_max": aniso, "iso_shift": iso,
        })
        print(f"  |k| = {kap:4.2f}: F_mean = {F.mean():.6f}  A1g = "
              f"{w['A1g']:.6f}  Eg = {w['Eg']:.6f}  T2g = {w['T2g']:.6f}  "
              f"aniso = {aniso:.3e}")
    a_an = _fit(kappas_used, aniso_vals)
    a_iso = _fit(kappas_used, iso_vals)
    w_last = rows[-1]
    spin2 = w_last["w_Eg"] + w_last["w_T2g"]
    print(f"  exponents: alpha_aniso = {a_an}, alpha_iso = {a_iso}")
    print(f"  discriminating: w_Eg + w_T2g at |k| = 0.4: {spin2:.4f}")

    ok = (err_cas < 1e-9 and max(devs) < 1e-12 and len(rows) == 4
          and F0 > 0.0 and all(math.isfinite(r["F_mean"]) for r in rows))

    verdict = {
        "id": "V176_gauged_vertex_content",
        "tier": 4,
        "property": "O_h (cubic-harmonic) content of the gauged "
                    "bond-current vertex Gamma^A = I_8 (x) T^A_16 on the "
                    "W9 node cone: form factor sum_A |<u(q*)|Gamma^A|"
                    "u(q*+k n)>|^2 over 26 directions, 4 |k| values",
        "residual": None,
        "pass": bool(ok),
        "status": f"GAUGED VERTEX CONTENT ON THE CONE (4 rows). The "
                  f"on-node current expectation F(0) = {F0:.6f} (the cone "
                  f"mode carries Spin(10) charge). The form factor "
                  f"decomposes as A1g = {w_last['w_A1g']:.4f}, Eg = "
                  f"{w_last['w_Eg']:.4f}, T2g = {w_last['w_T2g']:.4f} at "
                  f"|k| = 0.4: the gauged vertex is "
                  f"{'NOT a scalar - it carries the spin-2 (Eg+T2g) pieces the contact vertex lacked' if spin2 > 1e-6 else 'EXACTLY an O_h scalar on the cone (structural no-go: the cone mode is invisible to the spin-2 pieces)'}. "
                  f"Anisotropic part scales as k^{a_an:.2f}, the A1g shift "
                  f"as k^{a_iso:.2f} (V171 contact: alpha_aniso = 2, "
                  f"A1g = 1.0000). Casimir anchor 45 I_128 err "
                  f"{err_cas:.1e}; 2-route (16-D vs 32-D module) dev "
                  f"{max(devs):.1e}.",
        "data": {"rows": rows, "F0": F0, "alpha_aniso": a_an,
                 "alpha_iso": a_iso, "casimir_err": err_cas,
                 "route_dev": max(devs), "spin2_weight_k04": spin2},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v176_gauged_vertex_content.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v176_gauged_vertex_content.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
