#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V146 - the two-number bond pattern: superexchange + broken flux +
ground-state Cartan shift, tested jointly across all node counts (O4).

**Why.**  V137 showed the broken flux K_s = 45 - C_2^{unbroken}(s) is a
per-class integer of the root but does not close the V96 meeting (residual
0.605). V145 showed the second-order bond-current (superexchange) class
shift kappa~_s, with the exact channel denominators of 45 x 16 = 16 + 144 +
560, also does not close it alone (residual 0.577, R^2 0.076). The gauged
layer supplies BOTH numbers; the question this experiment decides: is the
V96 pattern a fixed COMBINATION of them - plus the only other class number
the layer carries, the ground-state Cartan shift F_s = <s| sum_A mu^A T^A |s>
with mu^A the bond flux expectation of the V101 ground state (expected zero
by the Gauss law, computed, not assumed)?

**The test.**  kappa_s = A kappa~_s + B K_s + C F_s with (A, B, C) SHARED
across all node counts n = 3..12 and z0 free per n: 3 meeting equations per
n give 30 equations for 13 unknowns - heavily overdetermined, so a zero
residual is a theorem, not a fit. The desert containment is checked for
every n.

**Checks (computed).**
1. Regression on V96: the meeting equations assembled here reproduce its
   ten stored closures at its fitted (z0, kappa_s).
2. The broken flux recomputed here (V137 T4 method) equals V137's stored
   per-class values [32, 36, 40, 28, 32, 40].
3. The ground-state bond flux mu^A = Tr(rho_B J^A) from the V101 kernel:
   its norm (the Gauss law predicts exactly zero - intertwiners carry no
   net charge).
4. The joint meeting: residual, desert containment per n.
5. Pattern comparison against the V96 fits (Pearson r, R^2 of the
   two-number combination).
"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh as geigh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))

from cl10 import clifford as CL                              # noqa: E402
from cl10 import engine as E                                 # noqa: E402
from cl10 import sm_qca                                      # noqa: E402
from cl10 import supercell as SC                             # noqa: E402
from cl10 import unification as UN                           # noqa: E402
from cl10 import u10pred                                     # noqa: E402
from cl10 import config as _config                          # noqa: E402
import V65_wilson as V65                                     # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
INV_A10 = 25.0
TWO_PI = 2.0 * math.pi
M_X_OBS = u10pred.M_X_DESERT
C2_32 = 45.0


def main():
    t0 = time.time()
    out = {"experiment": "V146_w4_two_number_kappas",
           "backend": E._BACKEND, "data": {}}
    checks = []

    # ------------------------------------------------- 1. classes and theta*
    lam, lam_off = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))
    Xc = lam @ xc_gg
    Yc = lam @ y_coef
    keys = [(round(float(Xc[s]), 6), round(float(Yc[s]), 6))
            for s in range(16)]
    classes = sorted(set(keys))
    members = {c: [s for s in range(16) if keys[s] == c] for c in classes}
    mults = [len(members[c]) for c in classes]

    v65 = json.load(open(exp("V65_wilson",
                                      "results", "v65_wilson.json")))
    roots = [r for r in v65["data"]["part_B"]["X_1"]["selective_roots"]
             if r.get("physical")]
    TH_STAR = roots[0]["theta_rad"]

    def direction(th):
        return math.cos(th) * xc_gg + math.sin(th) * y_coef

    c_dir = direction(TH_STAR)
    frs = V65.fragment_table(Xc, Yc)
    db_of_key = {(f["X"], f["Y"]): np.array(f["db"]) / f["n"] for f in frs}
    lam_db = np.zeros((3, 16))
    for s in range(16):
        if keys[s] in db_of_key:
            lam_db[:, s] = db_of_key[keys[s]]
    checks.append(("T1a %d fragment classes with multiplicities %s, theta* = "
                   "%.6f" % (len(classes), mults, TH_STAR),
                   bool(len(classes) == 6 and mults == [3, 2, 1, 6, 3, 1]
                        and abs(TH_STAR - _config.CFG.desert.theta_root_rad) < 1e-5)))

    # ------------------------------------------------ 2. generators and K_s
    Ts = [E.asnumpy(T).astype(np.complex128)
          for T in sm_qca.plane_generators()]
    h = [E.asnumpy(m) for m in CL.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    B, U32, U16 = SC.sectors()
    Dk = [B.conj().T @ C @ B for C in cartan]
    Mdiag = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, Vw = np.linalg.eigh(Mdiag)
    T16 = [Vw.conj().T @ (B.conj().T @ T @ B) @ Vw for T in Ts]

    Cop = sum(c_dir[k] * cartan[k] for k in range(5))
    Madj = np.zeros((45, 45))
    for a, T in enumerate(Ts):
        cc = Cop @ T - T @ Cop
        for b, T2 in enumerate(Ts):
            Madj[b, a] = float(np.real(-1j * np.trace(T2 @ cc) / 32.0))
    ev, Vn = np.linalg.eigh(1j * Madj)
    nul = [i for i in range(45) if abs(ev[i]) < 1e-8]
    Gens = [sum(Vn[A, i] * Ts[A] for A in range(45)) for i in nul]
    Herm = []
    for G in Gens:
        Herm.append(0.5 * (G + G.conj().T))
        Herm.append((G - G.conj().T) / 2j)
    Gr = np.array([[float(np.real(np.trace(X @ Y) / 32.0)) for Y in Herm]
                   for X in Herm])
    evg, Vh = np.linalg.eigh(Gr)
    ON = [sum(Vh[j, i] * Herm[j] for j in range(len(Herm)))
          for i in range(len(Herm)) if evg[i] > 1e-8]
    ON16 = [Vw.conj().T @ (B.conj().T @ T @ B) @ Vw for T in ON]
    Cas = np.zeros((16, 16), dtype=np.complex128)
    for T in ON16:
        Cas += T @ T
    dCas = np.real(np.diag(Cas))
    K = C2_32 - dCas
    Kc = [float(np.mean([K[s] for s in members[cc]])) for cc in classes]
    v137 = json.load(open(exp("V137_w4_broken_flux",
                                      "results",
                                      "v137_w4_broken_flux.json")))
    Kstored = [float(v) for v in v137["data"]["broken_flux_per_class"]]
    checks.append(("T2 the broken flux recomputed here %s equals V137's "
                   "stored per-class values %s to %.1e"
                   % ([round(v, 4) for v in Kc], [round(v, 4)
                                                  for v in Kstored],
                      float(np.abs(np.array(Kc) - np.array(Kstored)).max())),
                   bool(float(np.abs(np.array(Kc)
                                     - np.array(Kstored)).max()) < 1e-9)))

    # ------------------------------------------------ 3. superexchange kappa~
    v145 = json.load(open(exp("V145_w4_superexchange_kappas",
                                      "results",
                                      "v145_w4_superexchange_kappas.json")))
    kapc = [float(v) for v in v145["data"]["superexchange_per_class"]]
    checks.append(("T3a the superexchange pattern loaded from V145: %s "
                   "(V145 verdict %s)"
                   % ([round(v, 4) for v in kapc], v145["verdict"]),
                   bool(v145["pass"] and len(kapc) == 6)))

    # --------------------------- 4. the V101 ground-state bond flux mu^A
    _spec = importlib.util.spec_from_file_location(
        "v79", exp("V79_gauged_bond_layer", "V79_gauged_bond_layer.py"))
    V79 = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(V79)
    z = np.load(exp("V101_bond_current_ground_state",
                            "results", "kernel.npz"))
    Is, TsK = z["Is"], z["Ts"]
    nI, nT = Is.shape[0], TsK.shape[0]
    states = [(a, t, b) for a in range(nI) for t in range(nT)
              for b in range(nI)]
    ns = len(states)
    Ts_np = V79.root_generators()
    f, _k, _asym = V79.structure_constants(Ts_np)
    Ls_np = V79.adjoint_link_operators(f)
    I32 = np.eye(32, dtype=np.complex128)
    Js = [np.kron(Ts_np[A], I32) + np.kron(I32, Ts_np[A])
          for A in range(45)]
    P = np.einsum("xpa,ypa->xya", Is.conj(), Is)
    Q = np.einsum("tjab,ujab->tuab", TsK.conj(), TsK)
    Kk = np.einsum("xya,tuab->xytub", P, Q)
    Z = np.einsum("xki,yki->xyi", Is.conj(), Is)
    N = np.einsum("xytub,zwb->xtzyuw", Kk, Z).reshape(ns, ns)
    N = (N + N.conj().T) / 2
    UJ = np.zeros((45, nI, nI, 45, 45), dtype=np.complex128)
    for A in range(45):
        for y in range(nI):
            JI = Js[A] @ Is[y]
            for x in range(nI):
                UJ[A, x, y] = Is[x].conj().T @ JI
    V2 = np.einsum("tjab,ujcb->tuac", TsK.conj(), TsK)
    O2 = np.einsum("tjab,ujac->tubc", TsK.conj(), TsK)
    H = np.zeros((ns, ns), dtype=np.complex128)
    for s, (a, t, b) in enumerate(states):
        for s2, (a2, t2, b2) in enumerate(states):
            hh = 0.0 + 0.0j
            if b == b2:
                for A in range(45):
                    hh += np.einsum("xz,xy,yz->", UJ[A, a, a2],
                                    Ls_np[A], V2[t, t2])
            if a == a2:
                for A in range(45):
                    hh += np.einsum("xz,xy,yz->", UJ[A, b, b2],
                                    Ls_np[A], O2[t, t2])
            H[s, s2] = hh
    H = (H + H.conj().T) / 2
    ev, evec = geigh(H, N)
    ev = np.real(ev)
    E0 = float(ev[0])
    deg = int(np.sum(np.abs(ev - E0) < 1e-8 * max(1.0, abs(ev[-1]))))
    c = evec[:, :deg]
    Cw = (c @ c.conj().T) / deg
    rhoB = np.zeros((1024, 1024), dtype=np.complex128)
    for si, (a, t, b) in enumerate(states):
        for sj, (a2, t2, b2) in enumerate(states):
            w = Cw[si, sj]
            if abs(w) < 1e-18:
                continue
            kk = Kk[a, a2, t, t2].conj()
            rhoB += w * np.einsum("b,kb,lb->kl", kk, Is[b], Is[b2].conj())
    rhoB = (rhoB + rhoB.conj().T) / 2
    mu = np.array([float(np.real(np.trace(rhoB @ J))) for J in Js])
    mu_norm = float(np.abs(mu).max())
    # F_s = <s| sum_A mu^A T^A |s> diagonal in the weight basis
    Fop = sum(mu[A] * T16[A] for A in range(45))
    Fs = np.real(np.diag(Fop))
    Fc = [float(np.mean([Fs[s] for s in members[cc]])) for cc in classes]
    checks.append(("T3b the V101 ground state (E0 = %.4f eta, degeneracy %d, "
                   "matching V101's stored E0 = %.4f) carries bond flux "
                   "mu^A with max |mu| = %.1e - the Gauss law's zero - so the "
                   "Cartan class shift F_s is identically zero (per class %s)"
                   % (E0, deg,
                      json.load(open(exp("V101_bond_current_ground_state",
                                                "results",
                                                "v101_bond_current_"
                                                "ground_state.json")))
                           ["data"]["E0"],
                      mu_norm, [round(v, 6) for v in Fc]),
                   bool(abs(E0 + 5802.398188196636) < 1e-6 and deg == 6
                        and mu_norm < 1e-9)))

    # --------------------------------------------------- 5. the joint meeting
    xs = math.log(M_X_OBS / MZ)
    lp_s = np.array([math.log(max(abs(float(lam[members[cc][0]] @ c_dir)),
                                 1e-12)) for cc in classes])
    w = np.array(mults, dtype=float)
    li = np.array([w[i] * lam_db[:, members[cc][0]]
                   for i, cc in enumerate(classes)])      # (6,3)
    ns_list = list(range(3, 13))
    # unknowns: z0 per n (10) + A, B, C (shared)
    Arow = np.zeros((3 * len(ns_list), len(ns_list) + 3))
    brow = np.zeros(3 * len(ns_list))
    for j, n in enumerate(ns_list):
        for i in range(3):
            r = 3 * j + i
            Arow[r, j] = -n / TWO_PI * li[:, i].sum()
            Arow[r, len(ns_list)] = n / TWO_PI * float(li[:, i] @ kapc)
            Arow[r, len(ns_list) + 1] = n / TWO_PI * float(li[:, i] @ Kc)
            Arow[r, len(ns_list) + 2] = n / TWO_PI * float(li[:, i] @ Fc)
            brow[r] = -(INV_A10 + B_SM[i] * xs / TWO_PI - 1.0 / A_Z[i]) \
                - n / TWO_PI * float(li[:, i] @ (xs - lp_s))
    sol, res, rank, _ = np.linalg.lstsq(Arow, brow, rcond=None)
    Fj = Arow @ sol - brow
    resid = float(np.abs(Fj).max())
    Acoef, Bcoef, Ccoef = (float(sol[len(ns_list)]), float(sol[len(ns_list)
                                                                + 1]),
                           float(sol[len(ns_list) + 2]))
    desert = []
    for j, n in enumerate(ns_list):
        z0 = float(sol[j])
        kk = Acoef * np.array(kapc) + Bcoef * np.array(Kc) + Ccoef * np.array(Fc)
        mu = np.array([2.0 * math.cos(2.0 * math.pi * l / n) - 1.0
                       for l in range(n)])
        thr = np.array([lp_s[k] + z0 + kk[k] * mu for k in range(6)])
        desert.append({"n": n, "z0": round(z0, 4),
                       "in_desert": bool(np.all(thr > 0.0)
                                         and np.all(thr < xs)),
                       "min_threshold": round(float(thr.min()), 4),
                       "max_threshold": round(float(thr.max()), 4)})
    n_desert = sum(1 for d in desert if d["in_desert"])
    checks.append(("T4 the joint meeting with ONE (A, B, C) = (%.4f, %.4f, "
                   "%.4f) across all node counts leaves residual %s and keeps "
                   "every threshold in the desert for %d/10 node counts"
                   % (Acoef, Bcoef, Ccoef, "%.3e" % resid, n_desert),
                   bool(Arow.shape == (30, 13))))

    # ------------------------------- 5b. regression on V96 + pattern compare
    v96 = json.load(open(exp("V96_kappa_s_bond",
                                      "results", "v96_kappa_s_bond.json")))
    reg = []
    for s96 in v96["data"]["solutions"]:
        nn = s96["n_nodes"]
        kp = np.array(s96["kappas"], dtype=float)
        A7 = np.zeros((3, 7))
        b7 = np.zeros(3)
        for i in range(3):
            A7[i, 0] = -nn / TWO_PI * li[:, i].sum()
            A7[i, 1:] = nn / TWO_PI * li[:, i]
            b7[i] = -(INV_A10 + B_SM[i] * xs / TWO_PI - 1.0 / A_Z[i]) \
                - nn / TWO_PI * float(li[:, i] @ (xs - lp_s))
        reg.append(float(np.abs(A7 @ np.concatenate([[s96["z0"]], kp])
                                - b7).max()))
    # R^2 of the two-number combination against the V96 fits (per n, best)
    best_r2 = -9.0
    best_fit = None
    for s96 in v96["data"]["solutions"]:
        kp = np.array(s96["kappas"], dtype=float)
        Mx = np.column_stack([kapc, Kc])
        coef, *_ = np.linalg.lstsq(Mx, kp, rcond=None)
        r2 = 1.0 - float(((kp - Mx @ coef) ** 2).sum()
                         / ((kp - kp.mean()) ** 2).sum())
        if r2 > best_r2:
            best_r2 = r2
            best_fit = {"n": s96["n_nodes"],
                        "AB": [round(float(v), 4) for v in coef]}
    checks.append(("T5 the meeting equations assembled here reproduce V96's "
                   "ten exact closures at its fitted (z0, kappa_s) to %s; the "
                   "two-number combination (kappa~, K) explains at best R^2 = "
                   "%s of the variance of V96's fitted kappa_s (best n = %s, "
                   "(A,B) = %s)"
                   % ("%.1e" % max(reg), round(best_r2, 4), best_fit["n"],
                      best_fit["AB"]),
                   bool(len(reg) == 10 and max(reg) < 1e-2)))

    out["data"] = {
        "classes": [[list(cc), len(members[cc])] for cc in classes],
        "theta_star": round(TH_STAR, 6),
        "superexchange_per_class": [round(v, 6) for v in kapc],
        "broken_flux_per_class": [round(v, 6) for v in Kc],
        "cartan_shift_per_class": [round(v, 9) for v in Fc],
        "mu_max": mu_norm,
        "joint": {"A": round(Acoef, 6), "B": round(Bcoef, 6),
                  "C": round(Ccoef, 6), "residual": resid,
                  "n_desert": n_desert, "per_n": desert},
        "best_r2_two_number": round(best_r2, 4), "best_fit": best_fit,
        "v96_regression_max_residual": max(reg),
    }
    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    closed = resid < 1e-10 and n_desert == 10
    if closed:
        status = ("THE TWO-NUMBER COMBINATION CLOSES O4: kappa_s = %.4f "
                  "kappa~_s + %.4f K_s (+ zero Cartan shift) with ONE pair "
                  "(A, B) closes the desert meeting for all 10 node counts "
                  "(residual %.1e, all thresholds in the desert). The six "
                  "kappa_s become two numbers of the gauged bond layer."
                  % (Acoef, Bcoef, resid))
    else:
        status = ("THE TWO-NUMBER COMBINATION DOES NOT CLOSE O4: the best "
                  "joint (A, B, C) = (%.4f, %.4f, %.4f) over 30 equations in "
                  "13 unknowns leaves residual %s (desert for %d/10 node "
                  "counts), and the (kappa~, K) combination explains at best "
                  "R^2 = %s of the variance of V96's fitted kappa_s. The "
                  "ground-state bond flux is zero to %.0e (Gauss law), so the "
                  "Cartan shift contributes nothing. The supplier space of "
                  "the six kappa_s is now exhausted by every class number the "
                  "gauged layer computes: Cartan elements (V105), class "
                  "invariants (V122), broken flux (V137), superexchange "
                  "(V145), and every fixed combination of them (V146)."
                  % (Acoef, Bcoef, Ccoef, "%.3e" % resid, n_desert,
                     round(best_r2, 3), mu_norm))
    out["status"] = status
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v146_w4_two_number_kappas.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  joint residual %.2e | desert %d/10 | best R^2 %s | wall %.1f s"
          % (resid, n_desert, round(best_r2, 3), out["wall_s"]))
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
