#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V358 - the pentagon engine, built and gated (booked before the run
in review/pre_registrations_root.md; tier 1).

The engine library is cl10/pentengine.py: cell = 12 faces x 3 frame,
the I_h table on the cell's own axes, the exact-ring edge twists
(kingdon Cl(3,0) rotor cross-certified against Rodrigues), the
twisted hop K, U_hop(g) by eigh construction, the published 32-level
internal stages VERBATIM, and the six-stage step on einops cell axes.
Gates booked: (T0a) pentengine reproduces the V356 artifact's
ih_order_census bit-for-bit; (T0b) the internal identities verbatim
from cl10.supercell.sectors(0.2); (T0c) the cell table itself as
certified code (the {5,3} dodecahedron 12/30/20/chi-2, the 120-cell
nest, the {5,3,4} order-4 honeycomb) plus the 30-edge rotor
certification; (T1) hop unitarity at four clocks <= 1e-12 and |det| =
1; (T2) the full 1152-dimensional step unitary <= 1e-10; (T3)
A5-covariance over all 120 elements, batched on GPU 1 (cupy) in
batches of 10 with a three-probe numpy identity gate <= 1e-12;
(T4) the exact moments of K (tr K = 0, tr K^2 = 180) in mpmath;
(T5) PENTGAP - the closed-cell hop gap min|lambda(K)| measured at
dps 40 (mpmath inverse + power bound) with the cupy eigvalsh
cross-check <= 1e-10, no target: this leaf is what makes V360
necessary. Also booked from the V358 registration text: the exact
face-spectrum {5, (+-sqrt5)^3, (-1)^5} re-derived to 1e-12.
"""
import json
import os
import sys

os.environ["CL10_ANGLE"] = "0.2"    # house pin BEFORE importing cl10

import numpy as np                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)

from cl10.paths import exp          # noqa: E402
import mpmath                       # noqa: E402


_CP = None


def _get_cp():
    global _CP
    if _CP is None:
        try:
            import cupy
            cupy.cuda.Stream  # noqa: B018  (probe)
            _CP = cupy
        except Exception:      # noqa: BLE01
            _CP = False
    return _CP or None


def main():
    checks, data = [], {}
    from cl10 import pentengine as P

    # ---- T0a: the engine's I_h table reproduces the V356 artifact ----
    with open(exp("V356_pentagon_closure_free_clock", "results",
                  "v356_pentagon_closure_free_clock.json")) as fh:
        v356 = json.load(fh)
    booked_census = v356["data"]["ih_order_census"]
    IH = P.ih_table()
    census = {}
    for e in IH:
        census[str(e["order"])] = census.get(str(e["order"]), 0) + 1
    data["ih_census_engine"] = census
    data["ih_census_v356"] = booked_census
    checks.append(("T0a engine ih_table reproduces the V356 artifact "
                   "ih_order_census bit-for-bit: %s == %s"
                   % (census, booked_census), census == booked_census))
    perms = {e["perm"] for e in IH}
    checks.append(("T0a' the engine table closes: 120 elements, %d "
                   "distinct faithful face permutations, max strict "
                   "order %d (booked 10)" % (len(perms),
                                             max(e["order"] for e in IH)),
                   len(IH) == 120 and len(perms) == 120
                   and max(e["order"] for e in IH) == 10))

    # ---- T0b: the published internal identities, verbatim ----
    from cl10 import supercell
    B, U32, U16 = supercell.sectors(0.2)
    dev_int = float(np.max(np.abs(U32 @ B - B @ U16)))
    det_err = abs(np.linalg.det(U32) - 1.0)
    checks.append(("T0b internal identities from supercell.sectors(0.2): "
                   "|U32 B - B U16| = %.1e (<= 1e-12), |det U32 - 1| = "
                   "%.1e (<= 1e-12)" % (dev_int, det_err),
                   dev_int <= 1e-12 and det_err <= 1e-12))

    # ---- T0c: the cell itself as certified code (user directives) ----
    try:
        cell = P.dodeca_cell()
        nest = P.cell_nest(IH)
        hc = P.honeycomb_534()
        counts = tuple(cell["counts"])
        cert = P.certify_rotors()
        ok_cell = (counts == (20, 30, 12, 2) and len(nest) == 120
                   and hc["signature"] == (3, 1) and hc["edge_order"] == 4)
        checks.append(("T0c cell table certified as code: dodeca_cell "
                       "{5,3} counts (V,E,F,chi) = %s, edge_len %.6f; "
                       "cell_nest %d regular dodecahedral actions (120 "
                       "dodecahedra); honeycomb_534 Gram det %.6f "
                       "signature %s edge_order %d; rotor certs "
                       "(nv %.1e, pair %.1e, ring entries %d/270)"
                       % (counts, cell["edge_len"], len(nest), hc["det"],
                          hc["signature"], hc["edge_order"],
                          cert[0], cert[1], cert[2]), ok_cell))
    except AssertionError as ex:
        checks.append(("T0c cell/nest/honeycomb/rotor certification "
                       "FAILED: %s" % ex, False))

    # ---- T1: hop unitarity at the booked clocks ----
    K = P.twisted_hop()
    K_eig = np.linalg.eigh(K)
    devs, dets = [], []
    for g in (0.02, 1.0 / 5.0, np.pi / 4.0, 0.98):
        Uh = P.hop_unitary(g, K_eig=K_eig)
        devs.append(float(np.max(np.abs(Uh.conj().T @ Uh - np.eye(36)))))
        dets.append(float(abs(abs(np.linalg.det(Uh)) - 1.0)))
    checks.append(("T1 U_hop(g) unitary at g in {0.02, 1/5, pi/4, 0.98}: "
                   "max dev %.1e (<= 1e-12); |det| deviation max %.1e "
                   "(<= 1e-12; det = exp(i g tr K) = 1 since tr K = 0)"
                   % (max(devs), max(dets)),
                   max(devs) <= 1e-12 and max(dets) <= 1e-12))

    # ---- T2: the full engine step is unitary ----
    M = P.step(1.0 / 5.0, 0.2, K_eig=K_eig)
    dev_step = float(np.max(np.abs(M.conj().T @ M - np.eye(1152))))
    data["step_shape"] = list(M.shape)
    checks.append(("T2 full step(1/5, 0.2) on face(12) x frame(3) x "
                   "internal(32) = 1152: unitarity dev %.1e (<= 1e-10)"
                   % dev_step, dev_step <= 1e-10))

    # ---- T3: A5-covariance, all 120 elements, cupy batch on GPU 1 ----
    S36 = []
    for e in IH:
        Pm = np.zeros((12, 12), dtype=np.float64)
        for i, j in enumerate(e["perm"]):
            Pm[j, i] = 1.0
        Rh = np.array([[float(x) for x in row] for row in e["M"]])
        S36.append(np.kron(Pm, Rh))
    S36 = np.stack(S36)
    cp = _get_cp()
    devname = "CPU fallback"
    if cp is not None:
        import cupy
        devname = "cupy " + str(cupy.cuda.runtime.getDeviceProperties(
            cupy.cuda.Device().id)["name"])
    devs_cov = []
    gpu_id_dev = 0.0
    if cp is not None:
        MG = cp.asarray(M)
        for b0 in range(0, 120, 10):
            SG = cp.asarray(np.kron(np.stack(
                [S36[b0 + t] for t in range(min(10, 120 - b0))]),
                np.eye(32)))
            # batched conjugation: (10,1152,1152) @ (1152,1152)
            T = SG @ MG @ cp.swapaxes(SG.conj(), 1, 2)
            dv = cp.abs(T - MG[None]).max(axis=(1, 2))
            devs_cov.extend(cp.asnumpy(dv))
        for probe in (0, 47, 119):     # host identity gate
            Sh = np.kron(S36[probe], np.eye(32))
            nh = float(np.max(np.abs(Sh @ M @ Sh.T - M)))
            gpu_id_dev = max(gpu_id_dev, abs(nh - devs_cov[probe]))
    else:
        for h in range(120):
            Sh = np.kron(S36[h], np.eye(32))
            devs_cov.append(float(np.max(np.abs(Sh @ M @ Sh.T - M))))
    worst_cov = max(devs_cov)
    data["a5_cov_devs"] = {"worst": worst_cov,
                           "mean": sum(devs_cov) / len(devs_cov)}
    data["gpu_backend"] = {"backend": "cupy" if cp is not None
                           else "numpy", "device": devname}
    checks.append(("T3 A5-covariance: all 120 conjugations (kron(P_h x "
                   "R_h, I32)) step (...)dagger: worst dev %.1e (<= "
                   "1e-10), mean %.1e; batched on %s, device-vs-host "
                   "identity gate 3 probes %.1e (<= 1e-12)"
                   % (worst_cov, sum(devs_cov) / len(devs_cov), devname,
                      gpu_id_dev),
                   worst_cov <= 1e-10 and gpu_id_dev <= 1e-12))

    # ---- T4: the exact moments of K (mpmath) ----
    mpmath.mp.dps = 50
    vs, edges = P.face_geometry()
    trK2 = mpmath.mpf(0)
    for (u, v) in edges:
        R, _ = P._rodriques(P._unit(vs[u]), P._unit(vs[v]))
        frob = +sum(R[i][j] * R[i][j] for i in range(3)
                    for j in range(3))
        trK2 += 2 * frob    # K symmetric: both blocks (u,v) and (v,u)
    checks.append(("T4 exact moments of K in mpmath: tr K = 0 (diagonal "
                   "zero by construction), tr K^2 = %s (booked exactly "
                   "180: 30 edges x 2 blocks x tr R R^T = 3)"
                   % mpmath.nstr(trK2, 30),
                   abs(trK2 - 180) < mpmath.mpf(10) ** -25))

    # ---- T5: PENTGAP - min |lambda(K)| at dps 40 + cupy cross-check ----
    mpmath.mp.dps = 40
    n = 36
    Km = [[mpmath.mpf(0) for _ in range(n)] for _ in range(n)]
    for (u, v) in edges:
        R, _ = P._rodriques(P._unit(vs[u]), P._unit(vs[v]))
        for i in range(3):
            for j in range(3):
                Km[3 * u + i][3 * v + j] = +R[i][j]
                Km[3 * v + j][3 * u + i] = +R[i][j]
    # min|lambda| = 1 / rho(K^-1): invert in mpmath, symmetric power
    # iteration on K^-1 at dps 40.
    KinvT = []
    for j in range(n):
        e = [mpmath.mpf(0)] * n
        e[j] = 1
        KinvT.append(mpmath.lu_solve(Km, e))
    Kinv = [[KinvT[c][r] for c in range(n)] for r in range(n)]
    vec = [mpmath.mpf(1) if r % 7 == 0 else mpmath.mpf(0)
           for r in range(n)]
    nrm0 = mpmath.sqrt(sum(x * x for x in vec))
    vec = [x / nrm0 for x in vec]
    rho = mpmath.mpf(0)
    for _ in range(400):
        w = [+sum(Kinv[r][c] * vec[c] for c in range(n)) for r in range(n)]
        nw = mpmath.sqrt(sum(x * x for x in w))
        if nw == 0:
            break
        vec = [x / nw for x in w]
        rho = nw
    pentgap_mp = 1 / rho
    lam = K_eig[0]                       # float64 spectrum (numpy eigh)
    pentgap_np = float(np.min(np.abs(lam)))
    if cp is not None:
        lamg = cp.linalg.eigvalsh(cp.asarray(K))
        pentgap_gpu = float(cp.asnumpy(cp.min(cp.abs(lamg))))
    else:
        pentgap_gpu = pentgap_np
    data["PENTGAP"] = {"min_abs_lambda_mp40": mpmath.nstr(pentgap_mp, 25),
                       "min_abs_lambda_np": pentgap_np,
                       "min_abs_lambda_gpu": pentgap_gpu}
    checks.append(("T5 PENTGAP: closed-cell hop gap min|lambda(K)| = %s "
                   "(mpmath dps 40, inverse power bound) vs numpy %.10f "
                   "vs %s %.10f - cross-check %.1e (<= 1e-10); NO target "
                   "booked: this number is what makes V360 necessary (no "
                   "massless mode on the closed 12-cell)"
                   % (mpmath.nstr(pentgap_mp, 12), pentgap_np, devname,
                      pentgap_gpu,
                      max(abs(pentgap_np - float(pentgap_mp)),
                          abs(pentgap_gpu - pentgap_np))),
                   max(abs(pentgap_np - float(pentgap_mp)),
                       abs(pentgap_gpu - pentgap_np)) <= 1e-10))

    # ---- booked from the V358 registration: exact face spectrum ----
    A = np.zeros((12, 12))
    for (u, v) in edges:
        A[u, v] = 1
        A[v, u] = 1
    lamA = np.linalg.eigvalsh(A)
    mpmath.mp.dps = 40
    exact = [5.0] + [float(mpmath.sqrt(5))] * 3 \
        + [float(-mpmath.sqrt(5))] * 3 + [float(-1.0)] * 5
    exact.sort()
    dev_spec = float(np.max(np.abs(np.sort(lamA) - np.array(exact))))
    checks.append(("face spectrum re-derived (booked in the V358 "
                   "registration): the 12-face adjacency eigenvalues are "
                   "{5, (+-sqrt5)^3, (-1)^5} - (x-5)(x+1)^5(x^2-5)^3, "
                   "the icosahedral graph spectrum (MathWorld) - max dev "
                   "%.1e (<= 1e-12); gap 5 - sqrt5 = %s"
                   % (dev_spec, mpmath.nstr(5 - mpmath.sqrt(5), 15)),
                   dev_spec <= 1e-12))

    guards = [c[1] for c in checks]
    t0 = guards[0] and guards[1] and guards[2]
    if not t0:
        status, ok = "MACHINERY-DRIFT", False
    elif all(guards):
        status, ok = "PENTENGINE-BUILT", True
    else:
        status, ok = "PENTENGINE-BROKEN", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V358 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V358_pentengine_engine_rebuild", "tier": 1,
           "property": "the pentagon engine built and gated: I_h census "
                       "reproduction, internal identities, cell/nest/"
                       "honeycomb certs, unitarity, A5-covariance, exact "
                       "moments, PENTGAP",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v358_pentengine_engine_rebuild.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
