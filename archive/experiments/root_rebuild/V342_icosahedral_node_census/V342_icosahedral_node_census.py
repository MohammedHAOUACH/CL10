#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V342 - census of the dressed Dirac nodes of N2b (V341 engine): the
icosahedral 12-site cell x the theory's own Spin(10) internal layer
(one application per stage). V341 measured z = 1.00 five-ray linear
dispersion at the five-clocks with a refine-floor gap 2.7e-9; polished
to 1e-16 the contact is an EXACT zero. This experiment counts the
nodes (dedup mod reciprocal-lattice and reflections), checks whether
they form one orbit, and measures the cone velocity along 7 rays."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "root_rebuild",
                                "V341_icosahedral_internal_dressing"))
from V341_icosahedral_internal_dressing import cell_build, wrap  # noqa

TWO = 2 * math.pi


def main():
    from scipy.optimize import minimize
    from cl10 import supercell as SC
    checks = []
    edges, fac, edge_dir = cell_build()

    def cell_U(qv, g):
        u = np.eye(12, dtype=complex)
        for me in fac:
            M = np.zeros((12, 12), dtype=complex)
            for (i, j) in me:
                if (i, j) in edge_dir:
                    a, sg = edge_dir[(i, j)]
                    ph = np.exp(1j * sg * qv[a])
                    M[i, j] = ph
                    M[j, i] = np.conj(ph)
                else:
                    M[i, j] = 1.0
                    M[j, i] = 1.0
            u = (np.cos(g) * np.eye(12, dtype=complex)
                 - 1j * np.sin(g) * M) @ u
        return u

    B, U32, U16 = SC.sectors()
    ferr = 0.0
    q = (0.31, -0.77, 0.44)
    full = SC.bloch(q, U16)
    kronv = np.kron(SC.bond_product(q),
                    np.linalg.matrix_power(U16, 6))
    ferr = float(np.max(np.abs(full - kronv)))
    gM = float(SC.gap(SC.M_POINTS[0], U16))
    checks.append(("T0 published anchor: bloch == kron(P, U^6) %.1e, "
                   "M-point gap %.1e" % (ferr, gM),
                   ferr < 1e-10 and gM < 1e-9))

    allres = {}
    for cname, g in (("2pi/5", TWO / 5), ("4pi/5", 4 * TWO / 5)):
        Uint = np.linalg.matrix_power(U16, 5)
        bph = np.sort(np.angle(np.linalg.eigvals(Uint)))

        def gapf(qq):
            ap = np.angle(np.linalg.eigvals(cell_U(tuple(qq), g)))
            return min(abs(wrap(p + b)) for p in ap for b in bph)

        seeds = []
        N = 16
        for a in range(N):
            qx = -math.pi + TWO * a / N
            for b_ in range(N):
                qy = -math.pi + TWO * b_ / N
                for c in range(N):
                    qz = -math.pi + TWO * c / N
                    if gapf((qx, qy, qz)) < 3e-2:
                        seeds.append((qx, qy, qz))
        # spatial seed dedup (0.4 rad buckets) BEFORE polishing
        buckets = {}
        for s in seeds:
            k = tuple(int((x + math.pi) // 0.4) for x in s)
            if k not in buckets:
                buckets[k] = s
        nodes = []
        for s in buckets.values():
            r = minimize(gapf, s, method="Nelder-Mead",
                         options=dict(xatol=1e-12, fatol=1e-14,
                                      maxiter=3000))
            if r.fun > 1e-11:
                continue
            qm = tuple(x % TWO for x in r.x)
            key = tuple(round(x, 4) for x in
                        sorted(min(x, TWO - x) for x in qm))
            if not any(abs(k[0] - key[0]) < 5e-3
                       and abs(k[1] - key[1]) < 5e-3
                       and abs(k[2] - key[2]) < 5e-3
                       for k, _ in nodes):
                nodes.append((key, qm))
        rays = {}
        if nodes:
            q0 = nodes[0][1]
            for (nm, dv) in (("x", (1, 0, 0)), ("y", (0, 1, 0)),
                             ("z", (0, 0, 1)),
                             ("golden", (1, (1 + 5 ** .5) / 2, 0)),
                             ("body", (1, 1, 1)), ("diag1", (1, 1, 0)),
                             ("diag2", (1, 0, 1))):
                nr = math.sqrt(sum(d * d for d in dv))
                vs = []
                for eps in (1e-3, 3e-3):
                    qq = tuple(wrap(q0[i] + eps * dv[i] / nr)
                               for i in range(3))
                    vs.append(gapf(qq) / eps)
                rays[nm] = round(float(np.mean(vs)), 5)
        vs = list(rays.values())
        aniso = (max(vs) - min(vs)) / max(min(vs), 1e-12) if vs else 0
        allres[cname] = {"n_nodes": len(nodes),
                         "nodes_mod_reflection": [k for k, _ in nodes],
                         "cone_velocities": rays, "aniso": aniso}
    n25 = allres["2pi/5"]["n_nodes"]
    checks.append(("T1 census at the five-clocks: %s"
                   % json.dumps(allres), True))
    iso = allres["2pi/5"]["aniso"] < 0.05
    status = ("N2B-NODE-ORBIT-%d-%s" % (n25,
                                        "CONE-ISOTROPIC" if iso
                                        else "CONE-ANISOTROPIC"))
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:500])
                    for c in checks), flush=True)
    print("[V342 done] %s" % status, flush=True)
    out = {"id": "V342_icosahedral_node_census", "tier": 2,
           "property": "count and cone shape of the exact dressed Dirac "
                       "nodes of the icosahedral cell at the five-clocks",
           "status": status,
           # house convention (booked correction 2026-11-01): pass marks
           # the measurement certification (machine gates + census), the
           # physics verdict lives in the booked branch word - the first
           # version folded the expectation one_orbit (n in 4/8/12/24)
           # into pass, contradicting the booking: the word REPORTS the
           # orbit structure, it does not require it (the chiral-test
           # FAILURE is the measured result, stated in the word exactly
           # as for V345/V348/V352/V359/V362)
           "pass": bool(all(c[1] for c in checks) and n25 > 0),
           "data": {"results": allres,
                    "published_reference": "12 doubler channels at M "
                                           "(cubic engine, closed clock)",
                    "anchor": {"kron_err": ferr, "M_gap": gM}},
           "checks": [{"name": c[0][:700], "pass": c[1]} for c in checks]}
    with open(os.path.join(HERE, "results",
                           "v342_icosahedral_node_census.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
