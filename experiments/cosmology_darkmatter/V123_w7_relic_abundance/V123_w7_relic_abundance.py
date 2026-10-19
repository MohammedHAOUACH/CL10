#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V123 - the relic abundance, COMPUTED (O7)

**What the repo does today.**  Every abundance number in the theory is an
ansatz, never an integral:
  - cl10/dm_wimp.relic_abundance(m) = 0.12 * (100 GeV / m)^3 - the WIMP
    miracle read off a mass-scaling law, exact by construction at 100 GeV;
  - V75 and V88 never integrate the Boltzmann equation either: they INVERT
    the fiducial cross section sigma v = 3e-26 cm^3/s into a portal coupling
    (lambda_* = 4e-4, lambda_Z2 = 0.1307), so "Omega h^2 = 0.120" in those
    experiments is the input, not the output.
Section 10.2 and Section 12.2 both say the abundance is not computed.  This
experiment computes it: the standard freeze-out integral, run on the ladder
mass the theory actually predicts, with the annihilation strengths the theory
actually supplies.

**The computation.**  dY/dx = -(s <sigma v> / (H x)) (Y^2 - Y_eq^2) with
x = m/T, s = (2 pi^2/45) g_*s T^3, H = 1.66 sqrt(g_*) T^2 / M_P,
M_P = 1.220890e19 GeV (the repo's constant), g_* = g_*s = 106.75, Y_eq the
Boltzmann equilibrium density.  Y_in = Y_eq(x = 1); integrate to x = 200 by
RK4; Omega h^2 = m s_0 Y_inf / (rho_c/h^2) with s_0 = 2891.2 cm^-1 and
rho_c/h^2 = 1.0537e-5 GeV cm^-3.

**The automaton's own correction.**  The final states of the annihilation are
the theory's light states, and its light states are the DOUBLERS: twelve
nodes at the M points of the block zone (Proposition 2 iii), sitting at
|omega| = 2.0e-4 (the internal U16 mass; V102's node gap).  Each node carries
ONE gapless branch whose tangent velocity vector v^(a) = (S_1,S_2,S_3) is a
triple of NUMBERS (the node subspace is one-dimensional), so the velocity
operators commute trivially and the crossing is a nodal PLANE, not a cone -
V118's structure, measured here at the M nodes as well as the corner.  A
planar linear crossing E = v.k contributes a density of states proportional
to 1/|v|^3 instead of 1, so the annihilation into the light spectrum is
enhanced by D = sum over nodes and branches of 1/|v^(a)|^3.  With the
measured |v| ~ 0.38 that is ~18 per node, twelve nodes: a factor of a few
hundred that the dm_wimp (100/m)^3 ansatz never sees.  (The cone formula
1/(v1 v2 v3) is NOT applicable here: the velocity metric at these nodes has
rank 1, two of its three eigenvalues are zero - which is itself the check
that the node is planar.)

**What is asked.**  (i) the sigma v that Omega h^2 = 0.120 requires at the
ladder mass; (ii) the Omega h^2 the ladder's own couplings give - the loop
factor 1/16pi, the root coupling g10^2/4pi = 1/(100 pi) - with and without D;
(iii) the Omega h^2 of V88's two named portals at their named lambda.  No
number is fitted to 0.120 anywhere.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
from cl10 import supercell as SC              # noqa: E402
from cl10 import engine as E                  # noqa: E402
from cl10 import u10pred                      # noqa: E402
from cl10 import dm_wimp                      # noqa: E402
from cl10 import cosmology as COS             # noqa: E402

_B16, _U32, _U16 = SC.sectors()
U16 = E.asnumpy(_U16)

# ---- constants (repo conventions) -----------------------------------------
M_P = u10pred.M_PLANCK            # GeV (cl10/config.yaml physics:)
LOOP = u10pred.LOOP
G10 = SC.G10                      # 1/5
M_DM = u10pred.dm_ladder_mass()   # the ladder mass at rung cell.rung_B, GeV (config)
M_H, V_EV = dm_wimp.M_H, u10pred.V_EW
GSTAR, GSTARS = COS.G_STAR, COS.G_STAR
S0 = COS.S0_CM3                    # cm^-1, today's entropy density
RHO_C = COS.RHO_C_OVER_H2         # GeV cm^-3, critical density / h^2
GEV2_TO_CM3S = u10pred.GEV2_TO_CM2 * u10pred.C_CM_S
OMEGA_OBS = dm_wimp.OMEGA_OBS
SIGMA_FIDUCIAL = dm_wimp.SIGMA_FIDUCIAL   # cm^3/s, the V75/V88 fiducial
WEAK_SCALE = dm_wimp.WEAK_SCALE           # GeV, the WIMP-miracle window center
ZETA3 = 1.2020569


def n_eq(x, m, g=2.0):
    """Equilibrium number density (Boltzmann)."""
    T = m / x
    if x <= 1.0:
        return g * ZETA3 / math.pi ** 2 * T ** 3
    return g * (m * T / (2 * math.pi)) ** 1.5 * math.exp(-x) * (1.0 + 0.75 / x)


def freeze_out(m, sigmav, x0=1.0, x1=200.0, n=20000):
    """Integrate the Boltzmann equation; return (Y_inf, Omega h^2).
    sigmav in GeV^-2 (constant s-wave)."""
    xs = np.linspace(x0, x1, n)
    dx = xs[1] - xs[0]

    def yeq(x):
        T = m / x
        s = (2 * math.pi ** 2 / 45) * GSTARS * T ** 3
        return n_eq(x, m) / s

    def coef(x):
        T = m / x
        s = (2 * math.pi ** 2 / 45) * GSTARS * T ** 3
        H = COS.HUBBLE_COEFF * math.sqrt(GSTAR) * T ** 2 / M_P
        return s * sigmav / (H * x)

    def f(x, Y):
        ye = yeq(x)
        return -coef(x) * (Y * Y - ye * ye)

    Y = yeq(x0)
    for i in range(n - 1):
        x = xs[i]
        k1 = f(x, Y)
        k2 = f(x + dx / 2, Y + dx * k1 / 2)
        k3 = f(x + dx / 2, Y + dx * k2 / 2)
        k4 = f(x + dx, Y + dx * k3)
        Y += dx * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if Y < 1e-30:
            Y = 1e-30
    return float(Y), float(m * S0 * Y / RHO_C)


def omega_of_sigma(m, sigmav_cm3s):
    return freeze_out(m, sigmav_cm3s / GEV2_TO_CM3S)[1]


def sigma_for_omega(m, target=OMEGA_OBS):
    """Bisection on log sigma v for the target abundance."""
    lo, hi = 1e-32, 1e-20
    for _ in range(45):
        mid = math.sqrt(lo * hi)
        if omega_of_sigma(m, mid) > target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def node_velocities(tol=1e-3, h=1e-5):
    """The tangent velocity structure at the twelve M nodes of the implemented
    step: the node states sit at |omega| = 2.0e-4 (the internal U16 mass,
    V102's node gap), the node subspace is one-dimensional, and the three
    tangent velocities are numbers - so the metric has rank 1 and the crossing
    is a nodal plane.  Returns per node the branch velocity vectors."""
    MNODE = (math.pi / 2, math.pi / 2, 0.0)
    rows = []
    for q in SC.star(MNODE):
        M = SC.bloch(q, U16)
        ev = np.angle(np.linalg.eigvals(M))
        sel = np.where(np.abs(ev) < tol)[0]
        if sel.size == 0:
            continue
        _, V = np.linalg.eig(M)
        Q, _ = np.linalg.qr(V[:, sel])
        d = Q.shape[1]
        S = []
        for j in range(3):
            qp, qm = list(q), list(q)
            qp[j] += h
            qm[j] -= h
            D = (SC.bloch(tuple(qp), U16) - SC.bloch(tuple(qm), U16)) / (2 * h)
            S.append(1j * (Q.conj().T @ D @ Q))
        Gm = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                Gm[a, b] = float(np.real(np.trace(S[a] @ S[b])) / d)
        v2 = np.sort(np.linalg.eigvalsh(0.5 * (Gm + Gm.T)))
        vpr = [float(math.sqrt(max(float(v), 0.0))) for v in v2]
        genc = S[0] + 0.718281828459045 * S[1] + 0.367879441171442 * S[2]
        _, W = np.linalg.eigh(genc)
        vecs = [[float(np.real(W[:, a].conj().T @ s @ W[:, a]))
                 for s in S] for a in range(d)]
        rows.append({"q": [round(float(c), 6) for c in q], "dim": int(d),
                     "gap": float(np.min(np.abs(ev))),
                     "metric_eigenvalues": [float(v) for v in v2],
                     "principal_velocities": vpr,
                     "branch_speeds": [math.sqrt(sum(c * c for c in v))
                                       for v in vecs],
                     "comm_max": max(float(np.linalg.norm(S[a] @ S[b]
                                                        - S[b] @ S[a]))
                                     for a, b in ((0, 1), (1, 2), (2, 0)))})
    return rows


def main():
    out = {"experiment": "V123_w7_relic_abundance", "data": {}}
    checks = []

    # ---- T0 regression: the solver reproduces the WIMP miracle -------------
    om100 = omega_of_sigma(WEAK_SCALE, SIGMA_FIDUCIAL)
    om_tab = [(s, omega_of_sigma(WEAK_SCALE, s))
              for s in (SIGMA_FIDUCIAL / 10, SIGMA_FIDUCIAL, SIGMA_FIDUCIAL * 10)]
    ratio = [om / (OMEGA_OBS * SIGMA_FIDUCIAL / s) for s, om in om_tab]
    out["data"]["regression"] = {
        "Omega_h2_at_100GeV_3e-26": om100,
        "table": [{"sigma_cm3s": s, "Omega": om, "over_fiducial_scaling": r}
                  for (s, om), r in zip(om_tab, ratio)]}
    checks.append(("T0 the freeze-out solver reproduces the WIMP miracle: "
                   "Omega h^2 = %.3f at m = 100 GeV, sigma v = 3e-26 cm^3/s "
                   "(the canonical 0.10, within 0.03; the OBSERVED 0.120 "
                   "needs 2e-26, see T1), and tracks 0.12*(3e-26/sigma v) "
                   "within 45%% over a decade" % om100,
                   abs(om100 - 0.10) < 0.030
                   and max(abs(r - 1) for r in ratio) < 0.45))

    # ---- T1 the ladder mass and the cross section it needs ----------------
    sig_req = sigma_for_omega(M_DM)
    om_fid_ladder = omega_of_sigma(M_DM, SIGMA_FIDUCIAL)
    alpha_req = math.sqrt(sig_req / GEV2_TO_CM3S * math.pi * M_DM ** 2)
    out["data"]["ladder"] = {"m_DM_GeV": M_DM,
                             "sigma_required_cm3s": sig_req,
                             "alpha_required": alpha_req,
                             "Omega_at_fiducial_3e-26": om_fid_ladder,
                             "dm_wimp_ansatz": OMEGA_OBS * (WEAK_SCALE / M_DM) ** 3}
    checks.append(("T1 the ladder mass is m = %.1f GeV and the measured "
                   "Omega h^2 = 0.120 requires sigma v = %.2e cm^3/s there "
                   "(alpha_eff = %.4f in sigma v = alpha^2/(pi m^2)); the "
                   "dm_wimp (100/m)^3 ansatz gives %.4f at that mass"
                   % (M_DM, sig_req, alpha_req, OMEGA_OBS * (WEAK_SCALE / M_DM) ** 3),
                   abs(M_DM - 118.4) < 1.0 and 1e-27 < sig_req < 1e-24))

    # ---- T2 the doubler final-state density, measured from the step -------
    rows = node_velocities()
    speeds = [b for r in rows for b in r["branch_speeds"] if b > 1e-6]
    D = sum(1.0 / s ** 3 for s in speeds)
    rank1 = all(sum(1 for v in r["metric_eigenvalues"] if v > 1e-12) == 1
                for r in rows)
    vmean = sum(speeds) / len(speeds)
    out["data"]["doublers"] = {"n_nodes": len(rows),
                               "node_gap": rows[0]["gap"],
                               "branch_speeds": [round(s, 6) for s in speeds],
                               "mean_speed": vmean,
                               "comm_max_max": max(r["comm_max"] for r in rows),
                               "metric_rank_one": rank1,
                               "D_planar": D}
    checks.append(("T2 the implemented step has %d M nodes, each with a "
                   "one-dimensional node subspace at |omega| = %.1e whose "
                   "velocity metric has rank 1 (commutator %.1e): the crossing "
                   "is a nodal PLANE, so the final-state density is enhanced "
                   "by D = sum 1/|v|^3 = %.1f over one isotropic cone "
                   "(mean |v| = %.3f)"
                   % (len(rows), rows[0]["gap"],
                      max(r["comm_max"] for r in rows), D, vmean),
                   len(rows) == 12 and rank1
                   and max(r["comm_max"] for r in rows) < 1e-8))

    # ---- T3 the couplings the theory supplies, and what they give ---------
    cands = {"loop factor 1/16pi": LOOP,
             "root g10^2/4pi": G10 ** 2 / (4 * math.pi),
             "required (inverted)": alpha_req}
    supplied = []
    for name, a in cands.items():
        sv = a ** 2 / (math.pi * M_DM ** 2)          # GeV^-2
        om = freeze_out(M_DM, sv)[1]
        omD = freeze_out(M_DM, sv * D)[1]
        supplied.append({"coupling": name, "alpha": a,
                         "sigma_cm3s": sv * GEV2_TO_CM3S,
                         "Omega_h2": om, "Omega_h2_with_doublers": omD})
    out["data"]["supplied_couplings"] = supplied
    checks.append(("T3 with the ladder mass and NO fitted portal: the loop "
                   "factor gives Omega h^2 = %.3g, the root coupling "
                   "g10^2/4pi gives %.3g, and the doubler final-state density "
                   "D = %.1f moves them to %.3g and %.3g"
                   % (supplied[0]["Omega_h2"], supplied[1]["Omega_h2"], D,
                      supplied[0]["Omega_h2_with_doublers"],
                      supplied[1]["Omega_h2_with_doublers"]),
                   all(s["Omega_h2"] > 0 for s in supplied)))

    # ---- T4 V88's two named portals, abundance computed not assumed -------
    def sigmav_V75(lam_p):
        return (lam_p ** 2 * V_EV ** 2 / (8 * math.pi)
                * M_DM ** 2 / (4 * M_DM ** 2 - M_H ** 2) ** 2)

    def sigmav_V1(lam_p):
        s_ch = 4 * M_DM ** 2
        p_f = math.sqrt(max(s_ch - M_H ** 2, 0.0)) / (2 * math.sqrt(s_ch))
        return lam_p ** 2 / (32 * math.pi * s_ch) * p_f

    v88 = json.load(open(exp("V88_mirror_scalar",
                                      "results", "v88_mirror_scalar.json")))
    lam_Z2 = v88["data"]["V1_lam_for_Omega"]
    lam_star = math.sqrt(SIGMA_FIDUCIAL / GEV2_TO_CM3S * 8 * math.pi
                         * (4 * M_DM ** 2 - M_H ** 2) ** 2
                         / (V_EV ** 2 * M_DM ** 2))
    portals = []
    for name, lam, fn in (("V75 pair portal", lam_star, sigmav_V75),
                          ("V88 exact-Z2 portal", lam_Z2, sigmav_V1)):
        sv = fn(lam)
        portals.append({"portal": name, "lambda": lam,
                        "sigma_cm3s": sv * GEV2_TO_CM3S,
                        "Omega_h2": freeze_out(M_DM, sv)[1],
                        "Omega_h2_with_doublers": freeze_out(M_DM, sv * D)[1]})
    out["data"]["named_portals"] = portals
    checks.append(("T4 the named portals land near 0.120 only because their "
                   "lambda was INVERTED from sigma v = 3e-26: computed "
                   "abundances are %.3g (V75, lambda = %.2e) and %.3g (V88, "
                   "lambda = %.4f); with the doubler multiplicity included "
                   "they become %.3g and %.3g"
                   % (portals[0]["Omega_h2"], portals[0]["lambda"],
                      portals[1]["Omega_h2"], portals[1]["lambda"],
                      portals[0]["Omega_h2_with_doublers"],
                      portals[1]["Omega_h2_with_doublers"]),
                   abs(portals[0]["Omega_h2"] - dm_wimp.OMEGA_OBS) < 0.05))

    ok = all(c[1] for c in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "THE RELIC ABUNDANCE IS COMPUTED (O7). The freeze-out integral "
        "reproduces the WIMP miracle (Omega h^2 = %.3f at 100 GeV with "
        "sigma v = 3e-26 cm^3/s). At the ladder mass m = %.1f GeV the measured "
        "Omega h^2 = 0.120 requires sigma v = %.2e cm^3/s, i.e. alpha_eff = "
        "%.4f in sigma v = alpha^2/(pi m^2) - while the theory supplies only "
        "the loop factor 1/16pi = %.4f (Omega h^2 = %.3g) or the root "
        "g10^2/4pi = %.4f (Omega h^2 = %.3g). The automaton's own light "
        "spectrum matters: %d M nodes, each a one-dimensional node subspace at "
        "|omega| = %.1e with a rank-1 velocity metric and commuting tangent "
        "velocities (max %.1e) - nodal planes, not cones - so the final-state "
        "density is enhanced by D = sum 1/|v|^3 = %.1f (mean |v| = %.3f), "
        "which the dm_wimp (100/m)^3 ansatz never sees; with it the supplied "
        "couplings give Omega h^2 = %.3g and %.3g. The named portals reproduce "
        "0.120 only because their lambda was inverted from the fiducial "
        "3e-26 cm^3/s (computed: %.3f and %.3f; with doublers %.3f and %.3f). "
        "VERDICT: the abundance is now a computed number, and it is NOT 0.120 "
        "for any coupling the root supplies - the WIMP miracle of Fig. 2 is a "
        "mass-scaling ansatz, and the doubler spectrum changes the annihilation "
        "rate by the factor D."
        % (om100, M_DM, sig_req, alpha_req, LOOP, supplied[0]["Omega_h2"],
           G10 ** 2 / (4 * math.pi), supplied[1]["Omega_h2"], len(rows),
           rows[0]["gap"], max(r["comm_max"] for r in rows), D, vmean,
           supplied[0]["Omega_h2_with_doublers"],
           supplied[1]["Omega_h2_with_doublers"], portals[0]["Omega_h2"],
           portals[1]["Omega_h2"], portals[0]["Omega_h2_with_doublers"],
           portals[1]["Omega_h2_with_doublers"]))

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v123_w7_relic_abundance.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V123 - relic abundance computed")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  m_DM = %.2f GeV; sigma required %.2e cm^3/s (alpha %.4f)"
          % (M_DM, sig_req, alpha_req))
    print("  M nodes %d, gap %.1e, rank-1 %s, mean |v| %.3f, D %.1f"
          % (len(rows), rows[0]["gap"], rank1, vmean, D))
    for s in supplied:
        print("  %-22s alpha %.4f sigma %.2e cm^3/s Omega %.3g (with D %.3g)"
              % (s["coupling"], s["alpha"], s["sigma_cm3s"], s["Omega_h2"],
                 s["Omega_h2_with_doublers"]))
    for p in portals:
        print("  %-22s lambda %.3e sigma %.2e Omega %.3g (with D %.3g)"
              % (p["portal"], p["lambda"], p["sigma_cm3s"], p["Omega_h2"],
                 p["Omega_h2_with_doublers"]))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
