#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V137 - the per-class bond strength FROM the gauged layer (O4).

**Where O4 stands.**  V96 closes the desert meeting exactly once the bond
layer is allowed a junction strength per fragment CLASS, kappa_s (s = 0..5):
    z_{s,l} = z_s(theta, z0) + kappa_s (2 cos(2 pi l / n) - 1).
Three meeting equations, six unknowns: the meeting is a SOLUTION FAMILY, and
Table 5 classes the six kappa_s as fitted inputs.  V105 proved they are not
obtainable from any on-site Cartan element; V110 proved the most general
Hermitian bond operator on the 80 root-vector pairs never closes the meeting
(floor 30.4).  What was never tried is the number the gauged layer itself
supplies: the flux of BROKEN generators through each fragment class.

**The candidate, with no free parameter.**  At the V65 selective root theta*
the Wilson direction c leaves a 13-dimensional centralizer in so(10) - the
SU(3) x SU(2) x U(1)^3 of V67.  For a fragment state s of the 16, the
quadratic Casimir of that unbroken subalgebra,
    C_s = sum_{A in unbroken, Tr T^A T^B / 32 = delta^{AB}} <s| (T^A)^2 |s>,
is a number of the root, and so is the broken flux
    K_s = C_2(32) - C_s = 45 - C_s .
K_s is the second-order response of the fragment to the adjoint link: the
part of the link fluctuation that the unbroken algebra does not screen.  The
three checks below decide whether it is the kappa_s pattern.

**Checks (computed).**
1. Regression: the six fragment classes of V96 (multiplicities 3,2,1,6,3,1)
   and the theta* of the V65 physical selective root.
2. The class-blindness theorem of the FULL adjoint coupling: sum over all 45
   generators of sum_t |<t|T^A|s>|^2 = C_2(32) = 45 for every weight s - the
   naive second-order link shift is exactly class-degenerate (this is the
   V93/V95 degeneracy, now an identity).
3. The unbroken subalgebra: the centralizer of the Wilson direction in the 45
   has dimension exactly 13, with a Hermitian orthonormal basis.
4. The broken flux per class: 45 - C_s is constant on each class and takes
   four distinct integer values over the six classes - the automaton does
   supply a per-class number with no free parameter.
5. The meeting test: kappa_s = A K_s (one overall link scale A, the link gap
   is not fixed by the root) fed into the three meeting equations with (z0, A)
   the only unknowns, for n = 3..12: the residual and whether every threshold
   stays in the desert.  Zero residual = O4 becomes a number of the root.
6. The regression and the pattern comparison: the meeting equations assembled
   here reproduce V96's ten exact closures at its fitted (z0, kappa_s), and the
   Pearson correlation / variance explained (R^2) of V96's fitted kappa_s
   against the broken flux.
"""
import json
import math
import os
import sys
import time
import os

import numpy as np

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
from cl10 import config as _config                         # noqa: E402
from cl10 import u10pred                                     # noqa: E402
import V65_wilson as V65                                     # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
INV_A10 = 25.0
TWO_PI = 2.0 * math.pi
M_X_OBS = u10pred.M_X_DESERT
C2_32 = 45.0


def main():
    t0 = time.time()
    out = {"experiment": "V137_broken_flux", "backend": E._BACKEND,
           "data": {}}
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

    checks.append(("T1 regression on V96: %d fragment classes with "
                   "multiplicities %s, theta* = %.6f rad from the V65 physical "
                   "selective root, weight norm^2 spread %.1e"
                   % (len(classes), mults, TH_STAR,
                      float(np.einsum('ij,ij->i', lam, lam).max()
                            - np.einsum('ij,ij->i', lam, lam).min())),
                   bool(len(classes) == 6 and mults == [3, 2, 1, 6, 3, 1]
                        and abs(TH_STAR - _config.CFG.desert.theta_root_rad) < 1e-5)))

    # ------------------------------------------ 2. the class-blindness identity
    Ts = [E.asnumpy(T).astype(np.complex128)
          for T in sm_qca.plane_generators()]
    h = [E.asnumpy(m) for m in CL.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    B, U32, U16 = SC.sectors()
    Dk = [B.conj().T @ C @ B for C in cartan]
    Mdiag = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, Vw = np.linalg.eigh(Mdiag)
    T16 = [Vw.conj().T @ (B.conj().T @ T @ B) @ Vw for T in Ts]
    tot = np.array([float(np.real(sum(np.sum(np.abs(T[:, s]) ** 2)
                                      for T in T16))) for s in range(16)])
    checks.append(("T2 the full adjoint coupling is exactly class-blind: "
                   "sum_A sum_t |<t|T^A|s>|^2 = %s for the 16 weights (against "
                   "C_2(32) = 45, spread %.1e) - the second-order link shift "
                   "from the COMPLETE adjoint is a Casimir, hence identical for "
                   "every fragment: this is the V93/V95 class degeneracy as an "
                   "identity, and it is why kappa_s cannot come from the "
                   "unbroken-plus-broken sum"
                   % (np.round(tot, 6).tolist(), float(tot.max() - tot.min())),
                   bool(float(np.abs(tot - C2_32).max()) < 1e-9)))

    # --------------------------------------------- 3. the unbroken subalgebra
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
    gram = max(abs(float(np.real(np.trace(ON[i] @ ON[j]) / 32.0))
                   - (1.0 if i == j else 0.0))
               for i in range(len(ON)) for j in range(len(ON)))
    checks.append(("T3 the Wilson direction at theta* leaves a centralizer of "
                   "dimension %d in so(10) - the SU(3) x SU(2) x U(1)^3 of V67 "
                   "- with a Hermitian basis orthonormal to %.1e in the "
                   "Tr/32 pairing" % (len(ON), gram),
                   bool(len(ON) == 13 and gram < 1e-9)))

    # ------------------------------------------------- 4. the broken flux
    ON16 = [Vw.conj().T @ (B.conj().T @ T @ B) @ Vw for T in ON]
    Cas = np.zeros((16, 16), dtype=np.complex128)
    for T in ON16:
        Cas += T @ T
    dCas = np.real(np.diag(Cas))
    offCas = float(np.abs(Cas - np.diagflat(dCas)).max())
    K = C2_32 - dCas
    Kc = [float(np.mean([K[s] for s in members[cc]])) for cc in classes]
    Cc = [float(np.mean([dCas[s] for s in members[cc]])) for cc in classes]
    within = max(float(max(K[members[cc]]) - min(K[members[cc]]))
                 for cc in classes)
    checks.append(("T4 the broken flux K_s = 45 - C_2^{unbroken}(s) is a "
                   "per-class integer number of the root: per class %s "
                   "(unbroken Casimir %s; spread within a class %.1e, "
                   "off-diagonal Casimir %.1e) - %d distinct values over the six "
                   "classes, so the gauged layer does distinguish the fragments "
                   "without any input"
                   % ([round(v, 4) for v in Kc],
                      [round(v, 4) for v in Cc], within, offCas,
                      len(set(round(v, 6) for v in Kc))),
                   bool(within < 1e-9 and offCas < 1e-9
                        and len(set(round(v, 6) for v in Kc)) >= 4)))

    # --------------------------------------------------- 5. the meeting test
    xs = math.log(M_X_OBS / MZ)
    lp_s = np.array([math.log(max(abs(float(lam[members[cc][0]] @ c_dir)),
                                 1e-12)) for cc in classes])
    w = np.array(mults, dtype=float)
    li = np.array([w[i] * lam_db[:, members[cc][0]]
                   for i, cc in enumerate(classes)])      # (6,3)

    Khat = np.array(Kc)
    Khat = Khat / np.abs(Khat).max()
    rows = []
    for n in range(3, 13):
        # unknowns (z0, A): kappa_s = A * Khat_s
        A_eq = np.zeros((3, 2))
        b_eq = np.zeros(3)
        for i in range(3):
            A_eq[i, 0] = -n / TWO_PI * li[:, i].sum()
            A_eq[i, 1] = n / TWO_PI * float(li[:, i] @ Khat)
            b_eq[i] = -(INV_A10 + B_SM[i] * xs / TWO_PI - 1.0 / A_Z[i]) \
                - n / TWO_PI * float(li[:, i] @ (xs - lp_s))
        sol, res, rank, _ = np.linalg.lstsq(A_eq, b_eq, rcond=None)
        F = A_eq @ sol - b_eq
        z0, A = float(sol[0]), float(sol[1])
        zs = lp_s + z0
        kap = A * Khat
        mu = np.array([2.0 * math.cos(2.0 * math.pi * l / n) - 1.0
                       for l in range(n)])
        thr = np.array([zs[k] + kap[k] * mu for k in range(6)])
        in_desert = bool(np.all(thr > 0.0) and np.all(thr < xs))
        rows.append({"n": n, "residual": float(np.linalg.norm(F)),
                     "z0": round(z0, 4), "A": round(A, 4),
                     "kappa_s": [round(float(v), 4) for v in kap],
                     "thresholds_in_desert": in_desert,
                     "min_threshold": round(float(thr.min()), 4),
                     "max_threshold": round(float(thr.max()), 4)})
    best = min(rows, key=lambda r: r["residual"])
    closed = [r for r in rows if r["residual"] < 1e-10
              and r["thresholds_in_desert"]]
    checks.append(("T5 the meeting test: kappa_s = A K_s with the broken-flux "
                   "pattern leaves a residual of %s (best over n = 3..12: n = "
                   "%d, A = %s, z0 = %s, thresholds in the desert: %s); %d of "
                   "10 node counts close exactly with every threshold in the "
                   "desert - the broken flux %s against the V96 free pattern"
                   % ("%.3e" % best["residual"], best["n"], best["A"],
                      best["z0"], best["thresholds_in_desert"], len(closed),
                      [round(v, 3) for v in Kc]),
                   bool(len(rows) == 10)))

    # ------------------------------- 6. the broken flux against V96's fit
    v96 = json.load(open(exp("V96_kappa_s_bond",
                                      "results", "v96_kappa_s_bond.json")))
    Kv = np.array(Kc, dtype=float)
    reg = []
    cors = []
    for sol in v96["data"]["solutions"]:
        nn = sol["n_nodes"]
        kp = np.array(sol["kappas"], dtype=float)
        A7 = np.zeros((3, 7))
        b7 = np.zeros(3)
        for i in range(3):
            A7[i, 0] = -nn / TWO_PI * li[:, i].sum()
            A7[i, 1:] = nn / TWO_PI * li[:, i]
            b7[i] = -(INV_A10 + B_SM[i] * xs / TWO_PI - 1.0 / A_Z[i]) \
                - nn / TWO_PI * float(li[:, i] @ (xs - lp_s))
        reg.append(float(np.abs(A7 @ np.concatenate([[sol["z0"]], kp])
                                - b7).max()))
        r = float(np.corrcoef(kp, Kv)[0, 1])
        Afit = float(kp @ Kv / (Kv @ Kv))
        r2 = 1.0 - float(((kp - Afit * Kv) ** 2).sum()
                         / ((kp - kp.mean()) ** 2).sum())
        cors.append({"n": sol["n_nodes"], "pearson_r": round(r, 4),
                     "R2_proportional": round(r2, 4),
                     "kappas": [round(float(v), 4) for v in kp]})
    best_r2 = max(c["R2_proportional"] for c in cors)
    checks.append(("T6 the regression on V96 and the pattern comparison: the "
                   "meeting equations assembled here reproduce V96's ten exact "
                   "closures at its fitted (z0, kappa_s) to %s - the accuracy "
                   "V96's stored values allow, they are rounded to 1e-4 - so "
                   "the test above measures the same three equations; and the "
                   "broken flux against the six kappa_s V96 FITS to close gives "
                   "Pearson r %s with variance explained by one scale times the "
                   "flux R^2 %s (best %s) - the fitted pattern %s (n = 3) "
                   "against the root's flux %s"
                   % ("%.1e" % max(reg), [c["pearson_r"] for c in cors],
                      [c["R2_proportional"] for c in cors], round(best_r2, 4),
                      cors[0]["kappas"], [round(v, 3) for v in Kc]),
                   bool(len(cors) == 10 and max(reg) < 1e-2)))

    out["data"] = {
        "classes": [[list(cc), len(members[cc])] for cc in classes],
        "theta_star": round(TH_STAR, 6),
        "class_blind_sum": np.round(tot, 8).tolist(),
        "unbroken_dim": len(ON), "gram_leak": gram,
        "casimir_unbroken_per_class": [round(v, 6) for v in Cc],
        "broken_flux_per_class": [round(v, 6) for v in Kc],
        "meeting": rows, "best": best, "n_closed": len(closed),
        "against_v96_fit": cors, "best_r2": round(best_r2, 4),
        "v96_regression_max_residual": max(reg),
    }
    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    out["status"] = ("THE GAUGED LAYER SUPPLIES A PER-CLASS BOND NUMBER, BUT IT "
                     "IS NOT THE ONE THAT CLOSES THE MEETING (O4). At theta* the "
                     "Wilson direction leaves a 13-dimensional centralizer, and "
                     "the flux of BROKEN generators through each fragment class, "
                     "K_s = 45 - C_2^{unbroken}(s) = %s, is an integer number of "
                     "the root with no free parameter (the FULL adjoint sum is "
                     "class-blind to %.0e, which is the V93/V95 degeneracy as an "
                     "identity). Fed into the V96 meeting as kappa_s = A K_s it "
                     "leaves a residual %s and closes in the desert for %d/10 "
                     "node counts; against the six kappa_s V96 fits it explains "
                     "at best R^2 = %s of their variance. The kappa_s pattern is "
                     "%s"
                     % ([round(v, 3) for v in Kc], float(tot.max() - tot.min()),
                        "%.2e" % best["residual"], len(closed),
                        round(best_r2, 3),
                        "SUPPLIED by the root" if closed else
                        "still NOT supplied by the broken flux"))
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v137_broken_flux.json"), "w") as f:
        json.dump(out, f, indent=1)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  broken flux %s | best residual %.2e (n=%d) | closed %d/10 | "
          "wall %.1f s" % ([round(v, 3) for v in Kc], best["residual"],
                           best["n"], len(closed), out["wall_s"]))
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
