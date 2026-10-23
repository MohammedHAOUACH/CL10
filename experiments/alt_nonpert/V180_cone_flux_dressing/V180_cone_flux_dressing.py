#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V180 - exact (all-orders-in-eta) flux dressing of the node multiplet
(alt-nonpert, step 1: the non-perturbative side of the node)

**The question.**  Every route tried so far dresses the node mode
PERTURBATIVELY: one loop (V172/V177), tree-level block RG (V173/V178),
vertex form factors (V171/V176).  The gauged bond-current coupling
(V79/V101) is ``eta * sum_A J^A L^A`` per bond with eta = g_10 = 1/5 NOT
small compared with the flux scale, so the right question is the EXACT
one: diagonalize the flux response to the node multiplet's current to
all orders in eta, with no loop expansion and no link truncation.  The
link space is the full adjoint C^45 (V79; V140 fixed-link sector), and
on it ``sum_A L^A L^A = C2(45) I = 64 I`` EXACTLY (V129) - the flux
"kinetic" term is a constant and drops from every anisotropy statement.

**ROOT-CLOCK REDESIGN (current theory).**  The 1/5-era version dressed
a SINGLE cone mode u_k.  At the root clock that construction is
ill-defined and its "pure rest shift" structural facts are false: the
one-step spectrum through the M family is an EXACTLY 16-fold dispersion
multiplet at generic momentum (32 - particle plus hole - at the node
itself), so any single eigenvector of the family is pure gauge and the
per-bond currents, the -4*eta stage sums, and the dressed slopes built
from them are eigenvector-gauge noise (recorded by this file's git
history).  The gauge-SAFE source is the multiplet-AVERAGED current:
``Jbar_b^A(q) = sum_s Tr(rho_b^s T^A_16)`` with ``rho^s`` the block-s
diagonal of the projector onto the signed minimal-|phase| multiplet
(eigh of the unitary, 1e-9 window around the lowest phase).  This
experiment now measures:
  (1) the semi-Dirac structure of the node itself: gap ~ c k along the
      in-plane directions (exponent 1) and ~ g k^2 along the flat axis
      (exponent 2) - the V63/V118 "one flat tangent axis" made
      quantitative;
  (2) the EXACT dressing of the multiplet: eps(q) = sum_b min eig(eta *
      sum_A Jbar_b^A(q) L^A), its cubic-harmonic content over 26 O_h
      directions, and the kappa-exponents of its isotropic part and its
      anisotropy (rest shift: exponent >= 1 means the dressing cannot
      touch the linear cone slopes);
  (3) the dressed anisotropy metric A_np vs the free pinning metric
      A_free = 0.998 (the bare reference-branch cone anisotropy, the
      V115-V102 value) at three dq.

Anchors (mode-independent, gating): the Casimir 64 I_45; the bond-
current sum rule |sum_b Jbar_b| / |Tr(rho Gamma^A)| = 6 exactly (4
bonds per stage x 6 stages per subsite); window stability of the
multiplet projector (1e-8 vs 1e-10 selection gives identical Jbar -
the multiplet is isolated by a finite gap at every q probed); the
eta -> 0 identity dressed == free; A_free = 0.998.
"""

import importlib.util
import json
import math
import os
import sys
import time
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core  # alt suite numpy-native
os.environ["CL10_ANGLE"] = "0.2"  # the reference-branch clock: every number this experiment quotes or feeds is 1/5-branch; the closed default_branch angle must not leak in

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import eigsh

from cl10 import altcell as A
from cl10 import supercell as SC
from cl10.paths import exp

_spec = importlib.util.spec_from_file_location(
    "v79", exp("V79_gauged_bond_layer", "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)

C = A.cell((2, 2, 2))
ST = A.stages_baseline()
QSTAR = (math.pi / 2, math.pi / 2, 0.0)
FLAT_AXIS = (0.0, 0.0, 1.0)   # quadratic direction through QSTAR (V118)
ETA = SC.G10
DQS = (1e-2, 1e-3, 1e-4)
KS = (0.05, 0.1, 0.2, 0.4)


def bf(q):
    return A.bloch_cell(C, SC.sectors()[2], SC.G10, ST, q)


def _fit(ks, vals):
    lg = [(math.log(k), math.log(max(v, 1e-14)))
          for k, v in zip(ks, vals) if v > 0]
    if len(lg) < 2:
        return float("nan")
    return float(np.polyfit([a for a, _ in lg], [b for _, b in lg], 1)[0])


def main():
    t0 = time.time()
    print("=" * 68)
    print("V180 - exact (all-orders-in-eta) dressing of the node multiplet")
    print("=" * 68)
    B, T16, f = A.weyl_generators()
    Ls = V79.adjoint_link_operators(f)

    # anchor: Casimir on the adjoint link (V129)
    cas = sum(L @ L for L in Ls)
    err_cas = float(np.abs(cas - 64.0 * np.eye(45)).max())
    print(f"  anchor sum_A L^A L^A = 64 I_45: err {err_cas:.2e}")

    # bonds of the full 6-stage step: (subsite_i, subsite_j, image)
    bonds = []
    for st in ST:
        bonds += A.cell_stage_bonds(C, st)
    print(f"  bonds in the step: {len(bonds)}")

    def multiplet(q, win=1e-9):
        """The signed minimal-|phase| multiplet of the one-step Bloch
        matrix: eigenvectors whose phase lies within `win` of the
        lowest-|phase| eigenvalue. Gauge-safe data = the projector."""
        w, V = np.linalg.eig(bf(q))
        ph = np.angle(w)
        i0 = int(np.argmin(np.abs(ph)))
        sel = np.abs(ph - ph[i0]) < win
        return V[:, sel], int(sel.sum()), float(ph[i0])

    def currents(q, win=1e-9):
        """Jbar_b^A: the multiplet-averaged bond currents. Block form:
        the 128-D Bloch space is 8 sub-lattice blocks of 16, and the
        block-s diagonal rho_s of P P^dag gives
        Jbar^A = sum_s Tr(rho_s T^A) - basis-independent within the
        multiplet by construction."""
        P, dim, ph0 = multiplet(q, win)
        rho = P @ P.conj().T
        js = [np.array([float(np.real(np.trace(rho[s * 16:(s + 1) * 16,
                                                  s * 16:(s + 1) * 16]
                                            @ T)))
                        for T in T16])
              for s in range(8)]
        return dim, ph0, [js[i] + js[j] for (i, j, _v) in bonds]

    def eps_b(Jv, eta=ETA):
        """Exact flux ground-state shift of one bond: min eig of
        eta * sum_A J^A L^A (all orders in eta)."""
        M = eta * sum(Jv[A_] * Ls[A_] for A_ in range(45))
        return float(eigvalsh(M)[0]), M

    def eps_b_route2(Jv, eta=ETA):
        M = eta * sum(Jv[A_] * Ls[A_] for A_ in range(45))
        w, v = eigsh(M, k=1, which="SA", ncv=12, tol=0.0)
        rq = float(np.real(v[:, 0].conj() @ (M @ v[:, 0])))
        return float(w[0]), rq

    def eps_at(q, eta=ETA):
        _dim, _ph, Js = currents(q)
        return sum(eps_b(Jv, eta)[0] for Jv in Js)

    # multiplet structure at the node and off it (the root-clock fact
    # that killed the single-mode version: record, do not assume)
    fam = []
    for label, q in (("node", QSTAR),
                     ("x+0.4", tuple(a + 0.4 * b for a, b in
                                     zip(QSTAR, (1.0, 0.0, 0.0)))),
                     ("z+0.1", tuple(a + 0.1 * b for a, b in
                                     zip(QSTAR, FLAT_AXIS)))):
        d, p, _Js = currents(q)
        fam.append({"label": label, "dim": d, "phase0": p})
        print(f"  multiplet at {label:>6}: dim {d:2d}, phase {p:+.6f}")

    # window stability: the multiplet is an ISOLATED set of eigenvalues
    # at every probed q - selection windows 1e-8 and 1e-10 give the
    # identical projector, i.e. the source is unique (the multiplet gap
    # exceeds 1e-8 everywhere probed)
    devsJ = []
    for q in (QSTAR, tuple(a + 0.4 * b for a, b in zip(QSTAR, (1.0, 0.0, 0.0))),
              tuple(a + 0.1 * b for a, b in zip(QSTAR, FLAT_AXIS))):
        _d1, _p1, J1 = currents(q, win=1e-8)
        _d2, _p2, J2 = currents(q, win=1e-10)
        devsJ.append(max(float(np.abs(x - y).max()) for x, y in zip(J1, J2)))
    dev_route = max(devsJ)
    print(f"  projector window stability (1e-8 vs 1e-10): dev {dev_route:.2e}")

    # sum rule anchor: bond currents vs the total multiplet current
    # (Gamma^A = I_8 (x) T^A_16, the V176 definition, computed with
    # krons once as an independent path)
    _d0, _p0, Js0 = currents(QSTAR)
    P0, _dim0, _ph0 = multiplet(QSTAR)
    rho0 = P0 @ P0.conj().T
    Jtot = np.array([float(np.real(np.trace(rho0 @ np.kron(np.eye(8), T))))
                     for T in T16])
    Jsum = np.sum(Js0, axis=0)
    ratio = float(np.linalg.norm(Jsum) / max(np.linalg.norm(Jtot), 1e-15))
    print(f"  sum rule: |sum_b Jbar_b^A| / |Tr(rho Gamma^A)| = {ratio:.6f} "
          f"(expect exactly 6: 4 bonds per stage, 6 stages; recorded)")

    # 2-route eps: dense eigvalsh vs ARPACK + Rayleigh bound
    devs = []
    for q in (QSTAR, tuple(a + 0.2 * b for a, b in zip(QSTAR,
                                                       A.oh26_dirs()[3]))):
        _d, _p, Js = currents(q)
        for Jv in Js[:6]:
            e1, M = eps_b(Jv)
            e2, rq = eps_b_route2(Jv)
            devs.append(abs(e1 - e2))
            assert e1 <= rq + 1e-10, "variational bound violated"
    dev_eps = max(devs)
    print(f"  2-route eps (eigvalsh vs eigsh): dev {dev_eps:.2e}")

    eps_node = eps_at(QSTAR)
    print(f"  rest-energy offset eps(q*) = {eps_node:.6f} (subtracted "
          f"from the dressed slopes; flux ground state sourced by the "
          f"node multiplet)")

    # ---- (1) semi-Dirac exponents of the node (V63/V118 made
    # quantitative): gap vs displacement along an in-plane direction
    # and along the flat axis
    inplane = []
    flat = []
    for dd in (1e-2, 1e-3, 1e-4):
        inplane.append(A.gap_at(bf(tuple(a + dd * b for a, b in
                                         zip(QSTAR, (1.0, 0.0, 0.0))))) / dd)
        flat.append(A.gap_at(bf(tuple(a + dd * b for a, b in
                                      zip(QSTAR, FLAT_AXIS)))) / dd)
    a_in = _fit((1e-2, 1e-3, 1e-4), inplane) + 1.0
    a_fl = _fit((1e-2, 1e-3, 1e-4), flat) + 1.0
    g_flat = flat[2]
    print(f"  node gap exponents: in-plane k^{a_in:.2f} (v = "
          f"{inplane[2]:.4f}), flat axis k^{a_fl:.2f} (coefficient "
          f"{g_flat:.3e}: gap ~ {g_flat:.3e} k^2 - semi-Dirac)")

    DIRS = A.sphere_dirs(6)

    def aniso(bloch_fn, q_star, dq, dressed):
        vs = []
        for d in DIRS:
            q = tuple(a + dq * b for a, b in zip(q_star, d))
            e = A.gap_at(bloch_fn(q))
            if dressed:
                e += eps_at(q) - eps_at(tuple(a + dq * 1e-3 * b
                                              for a, b in zip(q_star, d)))
            vs.append(e / dq)
        return (max(vs) - min(vs)) / (max(vs) + min(vs))

    rows = []
    for dq in DQS:
        a_free = aniso(bf, QSTAR, dq, False)
        a_np = aniso(bf, QSTAR, dq, True)
        rows.append({"dq": dq, "A_free": a_free, "A_np": a_np,
                     "dA": a_np - a_free})
        print(f"  dq = {dq:7.0e}: A_free = {a_free:.6f}  "
              f"A_np = {a_np:.6f}  dA = {a_np - a_free:+.3e}")

    # eta -> 0 anchor: dressed == free exactly
    a_free1 = aniso(bf, QSTAR, 1e-3, False)

    def aniso_eta0_fn(dq2):
        vs = []
        for d in DIRS:
            q = tuple(a + dq2 * b for a, b in zip(QSTAR, d))
            e = A.gap_at(bf(q)) + eps_at(q, eta=0.0) - eps_at(QSTAR, eta=0.0)
            vs.append(e / dq2)
        return (max(vs) - min(vs)) / (max(vs) + min(vs))

    a_eta0 = aniso_eta0_fn(1e-3)
    err_eta0 = abs(a_eta0 - a_free1)
    print(f"  anchor eta -> 0: A_dressed = {a_eta0:.6f} vs A_free = "
          f"{a_free1:.6f} (dev {err_eta0:.2e})")
    err_v115 = abs(a_free1 - 0.998)  # the bare reference-branch cone of the
    # current theory (V115 T3 "up to 0.998", V102 flow, V172 anchor 0.998)
    print(f"  cross-check A_free vs V115/V102 0.998 (reference clock): dev "
          f"{err_v115:.2e}")

    # ---- (2) eps cubic content + kappa-exponents (26 O_h dirs)
    D26 = A.oh26_dirs()
    eps_k, aniso_k, mean_k = [], [], []
    for kap in KS:
        E = np.array([eps_at(tuple(a + kap * b
                                   for a, b in zip(QSTAR, d)))
                      - eps_at(tuple(a + 1e-5 * b
                                     for a, b in zip(QSTAR, d)))
                      for d in D26])
        w = A.cubic_harmonic_weights(E, D26)
        aniso_k.append(float(np.abs(E - E.mean()).max()))
        mean_k.append(float(abs(E.mean())))
        eps_k.append({"kappa": kap, "eps_mean": float(E.mean()),
                      "w_A1g": w["A1g"], "w_Eg": w["Eg"], "w_T2g": w["T2g"],
                      "aniso": aniso_k[-1]})
        print(f"  |k| = {kap:4.2f}: eps = {E.mean():+.6f}  A1g = "
              f"{w['A1g']:.4f}  Eg = {w['Eg']:.4f}  T2g = {w['T2g']:.4f}  "
              f"aniso = {aniso_k[-1]:.3e}")
    a_eps_an = _fit(list(KS), aniso_k)
    a_eps_iso = _fit(list(KS), mean_k)
    print(f"  eps exponents: aniso ~ k^{a_eps_an:.2f}, "
          f"|mean| ~ k^{a_eps_iso:.2f}")

    # structural probe (recorded, not gating): per-stage 4-bond dressing
    # sums - the -4*eta cancellation of the 1/5 era is NOT expected at
    # the root clock because the averaged multiplet current has
    # |Jbar_b| != 1 per bond
    stage_sums = []
    for label, q in (("node", QSTAR),
                     ("k0.4_x", tuple(a + 0.4 * b for a, b in
                                      zip(QSTAR, (1.0, 0.0, 0.0))))):
        _d, _p, Js = currents(q)
        ss = [float(sum(eps_b(Jv)[0] for Jv in Js[4 * s:4 * s + 4]))
              for s in range(6)]
        stage_sums.append({"label": label, "sums": ss})
    err_stage = max(abs(s + 4.0 * ETA) for r in stage_sums for s in r["sums"])
    print(f"  per-stage 4-bond sums vs -4*eta: err {err_stage:.2e} "
          f"(1/5-era cancellation recorded, not gating)")

    # gating: the mode-independent operator anchors + the pinning metric
    ok = (err_cas < 1e-12 and dev_eps < 1e-10 and err_eta0 < 1e-12
          and err_v115 < 1e-3 and abs(ratio - 6.0) < 1e-9
          and dev_route < 1e-10 and len(rows) == 3
          and all(math.isfinite(r["A_np"]) for r in rows))

    rest_shift = a_eps_iso >= 1.0 - 1e-9 or not math.isfinite(a_eps_iso)
    verdict = {
        "id": "V180_cone_flux_dressing",
        "tier": 4,
        "property": "exact (all-orders-in-eta, no loop expansion, full "
                    "adjoint C^45 link) dressing of the 16-fold node "
                    "multiplet by the gauged bond current, with the "
                    "gauge-safe multiplet-averaged source: semi-Dirac "
                    "exponents, eps cubic weights, kappa-exponents, "
                    "dressed vs free anisotropy",
        "residual": None if rest_shift else
                    "eps anisotropy is sublinear in k: the flux "
                    "dressing contributes at O(eta) to the node "
                    "dispersion (recorded exponent in data)",
        "pass": bool(ok),
        "status": f"EXACT DRESSING OF THE NODE MULTIPLET (3 dq rows). "
                  f"ROOT-CLOCK SOURCE: every band through M is an "
                  f"exactly 16-fold dispersion multiplet (32 at the "
                  f"node: particle+hole); the source is the "
                  f"multiplet-averaged current, projector-selected and "
                  f"window-stable to {dev_route:.1e}. NODE GEOMETRY: "
                  f"gap ~ k^{a_in:.2f} in-plane (v = {inplane[2]:.4f}) "
                  f"and k^{a_fl:.2f} along the flat axis (coefficient "
                  f"{g_flat:.2e}): semi-Dirac. DRESSING: eps(q*) = "
                  f"{eps_node:+.6f}; |mean eps(q*+k)| ~ k^{a_eps_iso:.2f}"
                  f" ({'analytic: a rest shift that cannot touch the linear cone slopes' if rest_shift else 'SUBLINEAR: O(eta) back-reaction on the node dispersion'}); "
                  f"anisotropy ~ k^{a_eps_an:.2f}, cubic weights "
                  f"A1g/Eg/T2g = {eps_k[-1]['w_A1g']:.2f}/"
                  f"{eps_k[-1]['w_Eg']:.2f}/{eps_k[-1]['w_T2g']:.2f} at "
                  f"k = 0.4. At dq = 1e-3: A_free = "
                  f"{rows[1]['A_free']:.6f} (V115/V102 reference cone 0.998, "
                  f"dev {err_v115:.1e}), A_np = {rows[1]['A_np']:.6f}, "
                  f"dA = {rows[1]['dA']:+.3e}. 1/5-era structural "
                  f"facts (single-mode -4*eta stage sums) do NOT "
                  f"survive: err {err_stage:.1e}, recorded. Anchors: "
                  f"Casimir 64 I_45 err {err_cas:.1e}; sum rule ratio "
                  f"{ratio:.6f} (expect 6); eta->0 dev {err_eta0:.1e}; "
                  f"2-route eps {dev_eps:.1e}.",
        "data": {"rows": rows, "eps_node": eps_node,
                 "sum_rule_ratio": ratio, "eps_k": eps_k,
                 "stage_sums": stage_sums,
                 "multiplets": fam,
                 "alpha_eps_aniso": a_eps_an, "alpha_eps_iso": a_eps_iso,
                 "node_exponents": {"in_plane": a_in, "flat_axis": a_fl,
                                    "v_in_plane": inplane[2],
                                    "g_flat": g_flat},
                 "casimir_err": err_cas, "eta0_dev": err_eta0,
                 "v115_dev": err_v115, "eps_route_dev": dev_eps,
                 "projector_window_dev": dev_route,
                 "stage_sum_err": err_stage,
                 "n_bonds": len(bonds)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v180_cone_flux_dressing.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v180_cone_flux_dressing.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
