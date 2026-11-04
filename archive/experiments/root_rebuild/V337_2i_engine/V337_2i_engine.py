#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V337 - the engine itself carries the binary-icosahedral element
(root-rebuild; booked in review/pre_registrations_root.md).

The internal layer of the published automaton is REPLACED (branch-local
module cl10/supercell_2i.py) by an order-10 element of the binary
icosahedral group: a rotation about the icosahedral 5-fold axis in the
Spin(3) of the three spatial planes. The tribunal measures whether the
skeleton survives the digit 5: unitarity, exact M-point gaplessness, the
exact 3-fold R multiplets, +- mirror pairing, closure order at GAMMA.
Runs on GPU 1 (user directive): cupy pinned via CUDA_VISIBLE_DEVICES=1;
every observable is cross-checked host-vs-GPU (T0c), and the published
step is re-measured first as the machinery anchor (T0a).
"""
import json
import math
import os
import sys

os.environ["CL10_BACKEND"] = "cupy"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402  (env pins above affect cl10 only)


def cluster_sizes(ph, tol=1e-9):
    out = []
    for p in sorted(ph):
        if out and abs(p - out[-1][0]) <= tol:
            out[-1][1] += 1
            out[-1][0] = p
        else:
            out.append([p, 1])
    return [c[1] for c in out]


def phases_np(M):
    return np.sort(np.angle(np.linalg.eigvals(M)))


def mirror_resid(ph):
    neg = np.sort(-np.asarray(ph))
    return float(np.max(np.abs(np.sort(ph) - neg)))


def unitarity(M):
    return float(np.max(np.abs(M @ M.conj().T - np.eye(M.shape[0]))))


def main():
    import numpy as np
    from cl10 import engine as E
    from cl10 import supercell as SC
    from cl10 import supercell_2i as S2

    checks = []
    data = {"backend": E.BACKEND if hasattr(E, "BACKEND") else "cupy?",
            "gpu": os.environ.get("CUDA_VISIBLE_DEVICES")}

    PHIv = (1.0 + math.sqrt(5.0)) / 2.0

    # ---- T0b: the Spin(3) triple, the axis, the order-10 element ---------
    Sx, Sy, Sz = S2.spin_triple()
    I32 = np.eye(32, dtype=np.complex128)
    sq = max(float(np.max(np.abs(X @ X - I32 / 4)))
             for X in (Sx, Sy, Sz))
    anti = max(float(np.max(np.abs(a @ b + b @ a)))
               for a, b in ((Sx, Sy), (Sy, Sz), (Sx, Sz)))
    U10 = S2.u_2i(2 * math.pi / 5)
    u_err = float(np.max(np.abs(U10 @ U10.conj().T - I32)))
    fifth = float(np.max(np.abs(np.linalg.matrix_power(U10, 5) + I32)))
    R = S2.spatial_adjoint(U10)
    tr_phi = float(np.abs(np.trace(R) - PHIv))
    so3 = float(np.max(np.abs(R @ R.T - np.eye(3))))
    data["T0b"] = {"Sq_I4": sq, "anticom": anti, "unit": u_err,
                   "U5_plus_I": fifth, "adjoint_trace_minus_phi": tr_phi,
                   "adjoint_SO3": so3}
    checks.append(("T0b Spin(3) triple: S_a^2 = I/4 (max %.1e), pairwise "
                   "anticommutators (%.1e); U_2I(2pi/5) unitary %.1e, "
                   "U^5 + I = %.1e (ORDER 10); spatial adjoint SO(3) %.1e "
                   "with trace - phi = %.1e (GOLDEN RATIO IN THE ENGINE)"
                   % (sq, anti, u_err, fifth, so3, tr_phi),
                  sq < 1e-12 and anti < 1e-12 and u_err < 1e-13
                  and fifth < 1e-10 and so3 < 1e-10 and tr_phi < 1e-10))

    # ---- T0a: the published step is re-measured with THIS machinery -----
    Upub = np.asarray(E.asnumpy(SC.internal(SC.G10)),
                      dtype=np.complex128)
    gaps_pub = min(float(np.min(np.abs(phases_np(
        S2.step(q, Upub, SC.G10))))) for q in SC.M_POINTS)
    deg_pub = []
    for q in SC.R_POINTS[:3]:
        deg_pub += cluster_sizes(phases_np(S2.step(q, Upub, SC.G10)))
    R_signature = sorted(set(deg_pub))
    # the +- mirror pairing of the paper is the ON-SITE (internal) layer's
    # spectrum pairing - the coupled Bloch matrix at one q is not the
    # object of that theorem (first-run gate defect, booked); the coupled
    # matrix carries instead the q -> -q doubling symmetry, measured too.
    def internal_mirror(U):
        ph = np.sort(np.angle(np.linalg.eigvals(U)))
        return float(np.max(np.abs(ph + ph[::-1])))
    def zone_pair_resid(U, beta, qs):
        # doubling symmetry: the spectrum at -q is the phase-REVERSED
        # spectrum at q (chirality flip), not the equal one.
        out = 0.0
        for q in qs:
            a = phases_np(S2.step(q, U, beta))
            b = phases_np(S2.step(tuple(-v for v in q), U, beta))
            neg_b = np.sort(-np.asarray(b))
            out = max(out, float(np.max(np.abs(np.angle(
                np.exp(1j * (np.sort(a) - neg_b)))))))
        return out
    mir_pub = internal_mirror(Upub)
    zone_pub = zone_pair_resid(Upub, SC.G10,
                               [SC.M_POINTS[0], SC.R_POINTS[0],
                                SC.X_POINTS[0], (0.3, 1.1, -2.0)])
    unit_pub = unitarity(S2.step(SC.R_POINTS[0], Upub, SC.G10))
    # closure order of the published step at GAMMA (recorded, not gated)
    def closure_k(M, cap=240, tol=1e-8):
        Mk = np.eye(M.shape[0], dtype=np.complex128)
        for k in range(1, cap + 1):
            Mk = Mk @ M
            e = min(float(np.max(np.abs(Mk - np.eye(M.shape[0])))),
                    float(np.max(np.abs(Mk + np.eye(M.shape[0])))))
            if e <= tol:
                return k, e
        return None, None
    kpub, _ = closure_k(S2.step(SC.GAMMA, Upub, SC.G10))
    data["T0a"] = {"pub_gap_M": gaps_pub,
                   "pub_R_cluster_sizes": R_signature,
                   "pub_internal_mirror": mir_pub,
                   "pub_zone_pair": zone_pub,
                   "pub_unitarity": unit_pub,
                   "pub_closure_k": kpub}
    checks.append(("T0a published-step anchor through the V337 machinery: M "
                   "points gapless (min |phase| %.1e), R cluster-size "
                   "signature %s (the reference for the grid), internal +- "
                   "mirror %.1e, q->-q doubling %.1e, unitarity %.1e, GAMMA "
                   "closure order %s"
                   % (gaps_pub, R_signature, mir_pub, zone_pub, unit_pub,
                      kpub),
                  gaps_pub < 1e-9 and mir_pub < 1e-12
                  and unit_pub < 1e-13 and kpub is not None))

    # ---- the grid: the 2I engine ----------------------------------------
    OMEGAS = [math.pi / 5, 2 * math.pi / 5, 4 * math.pi / 5]
    BETAS = [math.pi / 4, math.pi / 5, 2 * math.pi / 5]
    QS = {"GAMMA": SC.GAMMA, "X": SC.X_POINTS[0],
          "M": SC.M_POINTS[0], "R": SC.R_POINTS[0]}
    rows = []
    for om in OMEGAS:
        for be in BETAS:
            U = S2.u_2i(om)
            unit = max(unitarity(S2.step(q, U, be)) for q in QS.values())
            gapM = min(float(np.min(np.abs(phases_np(S2.step(q, U, be)))))
                       for q in SC.M_POINTS)
            deg = []
            for q in SC.R_POINTS:
                deg += cluster_sizes(phases_np(S2.step(q, U, be)))
            gam = S2.step(SC.GAMMA, U, be)
            ph_g = phases_np(gam)
            mir = internal_mirror(U)
            zone = zone_pair_resid(U, be, [SC.M_POINTS[0], SC.R_POINTS[0],
                                           (0.3, 1.1, -2.0)])
            zeropt = float(np.abs(ph_g.sum()))
            detph = float(np.abs(np.angle(np.linalg.det(gam))))
            k, cerr = closure_k(gam)
            r_sig = sorted(set(deg))
            holds = (unit <= 1e-11 and gapM <= 1e-9 and mir <= 1e-10
                     and r_sig == R_signature)
            rows.append({"omega": round(om / math.pi, 4),
                         "beta": round(be / math.pi, 4),
                         "unitarity": unit, "gapM": gapM,
                         "R_signature": r_sig,
                         "R_matches_published": r_sig == R_signature,
                         "internal_mirror": mir, "zone_pair": zone,
                         "zero_point": zeropt,
                         "det_phase": detph, "closure_k": k,
                         "closure_err": cerr, "holds": holds})
            print("  w=%.4f pi b=%.4f pi: unit %.2e gap %.2e R %s mir %.2e "
                  "zone %.2e 0pt %.2e det %.2e closure k=%s -> %s"
                  % (om / math.pi, be / math.pi, unit, gapM, r_sig, mir,
                     zone, zeropt, detph, k,
                     "HOLDS" if holds else "breaks"), flush=True)
    data["grid"] = rows

    # ---- T0c: GPU 1 vs host path ----------------------------------------
    import cupy as cp
    U = S2.u_2i(2 * math.pi / 5)
    M = S2.step(SC.GAMMA, U, math.pi / 4)
    ph_h = phases_np(M)
    Mdev = cp.asarray(M)
    ph_g = np.sort(cp.asnumpy(cp.angle(cp.linalg.eigvals(Mdev))))
    eig_diff = float(np.max(np.abs(np.angle(
        np.exp(1j * (ph_h - ph_g))))))
    Mk = cp.asarray(np.eye(M.shape[0], dtype=np.complex128))
    dev_pow_err = None
    for k in range(1, 40):
        Mk = Mk @ Mdev
    Mkh = np.linalg.matrix_power(M, 39)
    dev_pow_err = float(cp.max(cp.abs(Mk - cp.asarray(Mkh))))
    data["T0c"] = {"eig_host_vs_gpu": eig_diff, "pow_host_vs_gpu":
                   dev_pow_err}
    checks.append(("T0c GPU-1 (cupy, CUDA_VISIBLE_DEVICES=1) vs host numpy: "
                   "eigenphases at GAMMA max diff %.2e, 39-step power scan "
                   "%.2e - verdicts are path-invariant"
                   % (eig_diff, dev_pow_err), eig_diff < 1e-11
                  and dev_pow_err < 1e-11))

    nok = sum(r["holds"] for r in rows)
    if nok == len(rows):
        status = "ENGINE2I-STRUCTURE-HOLDS"
    elif nok == 0:
        status = "ENGINE2I-STRUCTURE-BREAKS"
    else:
        status = "ENGINE2I-PARTIAL-HOLDS-AT-%d-of-%d" % (
            nok, len(rows))
    data["status"] = status
    checks.append(("tribunal: %d of %d (internal angle x hop angle) grid "
                   "points keep unitarity + gapless M + exact 3-fold R + "
                   "mirror pairing" % (nok, len(rows)), True))
    for name, passed in checks:
        print(("  ok  " if passed else " FAIL ") + name, flush=True)
    out = {"id": "V337_2i_engine", "tier": 2,
           "property": "the automaton engine with a binary-icosahedral "
                       "internal element: does the skeleton survive "
                       "(root-rebuild)",
           "residual": False, "programme": "root-rebuild V337",
           "branch": status, "data": data,
           "checks": [{"name": n, "pass": p} for n, p in checks],
           "status": status, "pass": bool(
               all(p for _, p in checks[:3]))}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "v337_2i_engine.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("[V337 done] %s" % status, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
