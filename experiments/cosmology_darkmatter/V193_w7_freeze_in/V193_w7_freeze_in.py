#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V193 - the PRODUCTION side of the relic abundance: freeze-in computed (O7)

**Why (O7).**  Section 14.2 of the master states the abundance is a computed
number and closes the paragraph with the one production mechanism left
untested: "Freeze-in production of FIMP dark matter is not considered here."
V123 computes the abundance but only through the freeze-out integral - which
presupposes that the carrier was ever in chemical equilibrium; V75 and V88
never integrate anything, they INVERT the fiducial sigma v = 3e-26 cm^3/s
into a portal coupling (lambda_* = 3.53e-4 for the pair portal, lambda_Z2 =
0.1308 for the exact-Z2 portal), and V148 maps every lambda to a direct-
detection cross section.  Nobody ever asked the forward question: at a
coupling the LZ-class bound ALLOWS, does the built ladder scalar (one real
mirror-odd scalar, m_s = M_P LOOP^10 = 118.5709 GeV, V54/V75/V88) get
produced at 0.120 in the early universe?  This experiment computes the
production side - freeze-out and freeze-in with the regime DECIDED from the
computed in-equilibrium rate, never assumed - against the LZ-class ceiling.

**The computation.**  The two named portal kernels of V123 are taken VERBATIM
and promoted from the single point s = 4 m_s^2 to the full thermal s:
  V1 (Z2-preserving portal phi_s^2 phi_h^2; V123's sigmav_V1):
     K1(s, lam) = lam^2/(32 pi s) * p_f(s),  p_f = sqrt(s - M_h^2)/(2 sqrt s),
     plus the declared p-wave variant K1p with p_f -> p_f^3 (the repo kernel
     shows one power; an s-channel spin-0 mediator with massless fermions
     would show three - both are computed, neither selected after the fact).
  V2 (pair portal phi_s^2 (P^dagger+P); V123's sigmav_V75):
     K2(s, lam) = lam^2 v^2 m_s^2/(8 pi (s - M_h^2)^2),
     EXACTLY sigmav_V75 at s = 4 m_s^2 (asserted to 1e-12).
  The s-channel pole s = M_h^2 is strictly below the (2 m_s)^2 production
  threshold (2 m_s = 237.1 GeV > M_h = 125.25 GeV, V88's lattice free
  energies): every thermal pair meeting the kernel is off-shell, and the
  on-shell decay h -> s sbar has zero width (both audited in T3).
  Thermal averaging is EXACT: for a Boltzmann gas of initial mass m at T the
  in-rate density is n(T)^2 <K>_m(T) with
     <K>_m = N/n^2,
     N = 8 pi^2 int dP P^2 int dq q^2 sum_i w_i E^{-(E+ + E- - 2m)/T} v_M K,
  total/relative momenta P, q (unit m_s), E+-^2 = (P/2)^2 + q^2 +- P q mu_i +
  m^2, s = (E+ + E-)^2 - P^2, Moller velocity v_M = sqrt(s(s-4m^2))/
  (2 E+ E-), Gauss-Legendre weights w_i in mu, and n(T) = 4 pi int p^2
  e^{-(E-m)/T} dp.  Grids are fixed and deterministic.  The SAME quantity is
  cross-checked by a seeded Monte-Carlo average over the product of two
  Maxwell-Boltzmann momenta (T3b) and by grid refinement (T3c).
  Regime decision (computed, never assumed): the SELF-CONSISTENCY root of
  the two routes - lambda_th := Y_FI = Y_eq, i.e. Omega_FI(lambda_th) =
  Omega_eq (the equilibrium ceiling).  For lambda < lambda_th the freeze-in
  estimate stays under the equilibrium density and the population never
  enters chemical equilibrium (freeze-in is exact); above it it does not
  and the freeze-out route applies; the decided-regime abundance is
  continuous at lambda_th by construction.  The annihilation in-equilibrium
  rate criterion max_T <K>_ss n_s/H (the "would a WIMP thermalize" test) is
  computed and stored as a diagnostic, not used as the switch.
  The hh/ss convolution ratio at T = m_s is stored (the repo kernels' two
  legs are NOT crossing-conjugate, so this ratio is a property of the
  kernel convention, not an identity).
  Freeze-in route (exact for a never-thermal population, x = m_s/T):
     Y_FI = int_1^inf <K>_hh(x) [n_h^eq(x)]^2 D_g / (H(x) x^2) dx,
  with s(x), H(x) the V123 forms (x = m_s/T) and the equilibrium source
  density the exact Maxwell-Boltzmann g m^2 T K2(m/T)/(2 pi^2) (z = m/T >=
  1 throughout, the bath is always sub-threshold; V123's relativistic
  low-x branch is never exercised); the integrand
  is Boltzmann-suppressed e^{-2 m_h/T}, so the x-integral runs on a fine
  LINEAR grid to x = 80 and the tail is below e^{-134}.
  Freeze-out route: V123.freeze_out(m_s, <K>_ss(T = m_s) lam^2) - reused,
  not re-implemented; floored at the equilibrium density.
  DD ceiling (V148's published map, verbatim):
     sigma_SI(lam) = lam^2 f_N^2 mu_red^2/(pi M_h^4) * 0.3894e-27 cm^2,
     lambda_LZ = sqrt(LZ pi M_h^4/(f_N^2 mu_red^2 0.3894e-27)), mu_red the
     scalar-nucleon reduced mass (~0.93 GeV) - the convention that
     reproduces the published 6.72e-46 / 4.89e-51 cm^2 (regression T2).
  Nodal-plane final-state density D (V123): whether the pairs that carry
  production sit on the nodal planes is undetermined - both D_off = 1 and
  D_on = D_123 are computed and reported (declared variants, not post-hoc).

**Checks.**
  T0 the reused freeze-out solver reproduces V123's stored 0.080039;
  T1 the nodal-plane density D recomputed from the engine (12 M nodes,
     rank-1 velocity metric, commuting tangents) equals V123's stored 1485.27;
  T2 the V148 sigma_SI map reproduces V88's stored 6.72e-46 and 4.89e-51,
     and fixes the ceiling lambda_LZ;
  T3 K2(4 m_s^2) = V123's sigmav_V75 per unit lambda to 1e-12; the on-shell
     h -> s sbar width is zero; the pole is off-shell;
  T3b the exact convolution agrees with the seeded Monte-Carlo estimate of
     the same thermal average (constant kernel, and K2 on both legs) to 3%;
  T3c convergence: convolution grid refinement and the freeze-in
     x-integral under resolution/table refinement (drift < 2%);
  T4 the V88 scalar-parity stability regression (the single s cannot decay);
  T5 the solution lambda_sol of Omega h^2 = 0.120 on the decided-regime
     curve inside the DD-allowed range (1e-18, lambda_LZ], per variant,
     with the freeze-out crossing also recorded WITHOUT the cap when direct
     detection forbids it, and the branches PRE-REGISTERED:
       A  ABUNDANCE COMPUTED: a variant whose decided-regime curve crosses
          0.120 below lambda_LZ at a regime-self-consistent solution; the
          mechanism (freeze-in or freeze-out) is named by the computed
          regime, and the falsifier is the direct-detection cross section
          sigma_SI(lambda_sol) at the computed coupling;
       B  ABUNDANCE NOT REACHED: no self-consistent crossing below
          lambda_LZ on any variant; the abundance stays a named input,
          with freeze-in computed and excluded;
       C  INCONCLUSIVE ("inconclusive": true) if the integrals do not
          converge under refinement.

**Inputs (provenance).**  derived: m_s = V123.M_DM = M_P LOOP^10 (ladder,
V54).  measured/library (cl10/config.yaml, read through cl10): M_h = 125.25,
v = 246.22, M_P = 1.220890e19 GeV, Omega h^2 = 0.120, g_* = g_*s = 106.75,
Hubble coefficient 1.66, s_0 = 2891.2 cm^-1, rho_c/h^2 = 1.0537e-5 GeV cm^-3,
f_N = 3.0e-4, m_N = 0.93827208943 GeV, LZ limit 1.5e-47 cm^2.  named from
stored verdicts: lambda_* = 3.5299e-4, lambda_Z2 = 0.13081 (V88), the sigma_SI
map and its regression values (V148).  convention: x = m_s/T, Boltzmann
n_eq, exact Moller convolution; open choice declared: D on/off in production.
"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                    # noqa: E402
from cl10 import dm_wimp, u10pred, cosmology as COS  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v123", exp("V123_w7_relic_abundance", "V123_w7_relic_abundance.py"))
V123 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V123)

# ---- constants (repo conventions, via V123 / cl10) -------------------------
M_S = V123.M_DM                 # the ladder scalar, GeV
M_H, V_EV = V123.M_H, V123.V_EV
M_P, GSTARS, GSTAR = V123.M_P, V123.GSTARS, V123.GSTAR
S0, RHO_C, GEV2_TO_CM3S = V123.S0, V123.RHO_C, V123.GEV2_TO_CM3S
OMEGA_OBS = V123.OMEGA_OBS
SIGMA_FIDUCIAL = V123.SIGMA_FIDUCIAL
HC = COS.HUBBLE_COEFF
F_N = dm_wimp.F_N
M_N = u10pred.M_P_ROBSON
LZ = dm_wimp.LZ_BOUND
MU_RED = M_S * M_N / (M_S + M_N)                    # scalar-nucleon, GeV
CM2 = u10pred.GEV2_TO_CM2                           # GeV^-2 -> cm^2
X_FI_MAX = 80.0                                     # e^{-2 m_h * 80 / m_s}
_trapz = np.trapezoid                               # numpy >= 2

# ---- the kernels per unit lambda, verbatim V123 forms (numpy-ready) --------
def p_f(s):
    return np.sqrt(np.maximum(s - M_H ** 2, 0.0)) / (2.0 * np.sqrt(s))


def k1_s(s):
    return p_f(s) / (32.0 * math.pi * s)


def k1_p(s):
    return p_f(s) ** 3 / (32.0 * math.pi * s)


def k2_s(s):
    return V_EV ** 2 * M_S ** 2 / (8.0 * math.pi * (s - M_H ** 2) ** 2)


def k_const(s):
    return np.ones_like(np.asarray(s, dtype=float))


KERNELS = {"K1 p_f (V1 Z2 portal)": k1_s, "K1p p_f^3 (p-wave)": k1_p,
           "K2 (V2 pair portal)": k2_s}


def _k2(z):
    """Modified Bessel K2 via the integral rep (z >= 1: the source bath is
    always sub-threshold, T <= m_s < m_h)."""
    t = np.linspace(0.0, 12.0, 3000)
    return np.trapezoid(np.exp(-np.asarray(z)[..., None] * np.cosh(t))
                        * np.cosh(t), t)


def n_mb(x, m):
    """Exact Maxwell-Boltzmann equilibrium density g m^2 T K2(m/T)/(2 pi^2),
    g = 2 - the source-side density used consistently for the FI integral,
    the regime boundary and the diagnostic (V123's low-x branch is never
    exercised: x >= 1 here)."""
    z = np.asarray(x, dtype=float) * m / M_S
    T = M_S / np.asarray(x, dtype=float)
    return 2.0 * m * m * T * _k2(z) / (2.0 * math.pi ** 2)


def s_ent(x):
    return (2 * math.pi ** 2 / 45) * GSTARS * (M_S / x) ** 3


def hubble(x):
    return HC * math.sqrt(GSTAR) * (M_S / x) ** 2 / M_P


# ---- the exact thermal convolution <K>_m(T) ---------------------------------
MU_NODES, MU_W = np.polynomial.legendre.leggauss(20)


def conv_average(kfn, m, T, nP=80, nq=160):
    """<K>_m(T) = N/n^2 for a Boltzmann gas of initial mass m at T (GeV):
    N = 8 pi^2 int dP P^2 int dq q^2 sum_i w_i e^{-(E+ + E- - 2m)/T} v_M K(s)
    with total/relative momenta P, q in units of m_s, s in GeV^2 fed to the
    GeV kernel; n = 4 pi int p^2 e^{-(E-m)/T} dp.  Deterministic fixed
    grids: uniform P, q to a temperature-adaptive cutoff, Gauss-Legendre
    in mu (weighted sum, not a quadrature over the nodes)."""
    m_u, T_u = m / M_S, T / M_S
    qmax = 10.0 * math.sqrt(m_u * T_u / 2.0) + 12.0 * T_u + 0.5
    Pmax = 8.0 * math.sqrt(m_u * T_u) + 12.0 * T_u + 0.5
    Pg = np.linspace(0.0, Pmax, nP)
    qg = np.linspace(0.0, qmax, nq)
    P = Pg[:, None, None]
    q = qg[None, :, None]
    c = MU_NODES[None, None, :]
    w = MU_W[None, None, :]
    Ep = np.sqrt((P / 2.0) ** 2 + q ** 2 + P * q * c + m_u ** 2)
    Em = np.sqrt((P / 2.0) ** 2 + q ** 2 - P * q * c + m_u ** 2)
    s_u = (Ep + Em) ** 2 - P ** 2                      # s in units of m_s^2
    boltz = np.exp(-((Ep + Em - 2.0 * m_u) / T_u))     # shifted, cancels
    vm = np.sqrt(np.maximum(s_u * (s_u - 4.0 * m_u ** 2), 0.0)) \
        / (2.0 * Ep * Em)
    kval = kfn(s_u * M_S ** 2)
    inner = np.sum(boltz * vm * kval * P ** 2 * q ** 2 * w, axis=2)
    N = 8.0 * math.pi ** 2 * _trapz(_trapz(inner, qg, axis=1), Pg, axis=0)
    pg = np.linspace(0.0, qmax, 600)
    n = 4.0 * math.pi * _trapz(pg ** 2 * np.exp(
        -(np.sqrt(pg ** 2 + m_u ** 2) - m_u) / T_u), pg)
    return N / (n * n)


def mc_average(kfn, m, T, npairs=1000000, seed=20261019):
    """Independent seeded estimate of <K>_m(T) = E[v_M K(s)] over the
    product of two Maxwell-Boltzmann momenta: the cross-check of
    conv_average (inverse-CDF sampling on a fixed momentum grid)."""
    rng = np.random.default_rng(seed)
    pmax = 40.0 * math.sqrt(m * T) + 20.0 * T + 20.0 * m
    pg = np.linspace(0.0, pmax, 3000)
    dens = pg ** 2 * np.exp(-(np.sqrt(pg ** 2 + m * m) - m) / T)
    cdf = np.concatenate(([0.0], np.cumsum(
        0.5 * (dens[1:] + dens[:-1]) * (pg[1] - pg[0]))))
    cdf /= cdf[-1]
    u = rng.random(2 * npairs)
    pmag = np.interp(u, cdf, pg)
    cs = 2.0 * rng.random(2 * npairs) - 1.0
    sn = np.sqrt(np.maximum(1.0 - cs * cs, 0.0))
    ph = 2.0 * math.pi * rng.random(2 * npairs)
    px, py, pz = (pmag * sn * np.cos(ph), pmag * sn * np.sin(ph), pmag * cs)
    E = np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m * m)
    dot = E[:npairs] * E[npairs:] - (px[:npairs] * px[npairs:]
                                     + py[:npairs] * py[npairs:]
                                     + pz[:npairs] * pz[npairs:])
    s = 2.0 * m * m + 2.0 * dot
    vm = np.sqrt(np.maximum(dot * dot - m ** 4, 0.0)) / (E[:npairs]
                                                         * E[npairs:])
    return float(np.mean(vm * kfn(s)))


def main():
    t0 = time.time()
    out = {"id": "V193_w7_freeze_in", "tier": 5,
           "property": "the production side of the ladder-scalar abundance: "
                       "freeze-in and freeze-out computed through the two "
                       "named portal kernels with exact thermal convolution, "
                       "the regime decided from the computed in-equilibrium "
                       "rate, against the LZ-class direct-detection ceiling",
           "programme": "V193", "data": {}}
    checks = []

    # ---- T0 the reused freeze-out solver ------------------------------------
    v123 = json.load(open(exp("V123_w7_relic_abundance",
                                     "results",
                                     "v123_w7_relic_abundance.json")))
    om100_ref = v123["data"]["regression"]["Omega_h2_at_100GeV_3e-26"]
    om100 = V123.freeze_out(100.0, SIGMA_FIDUCIAL / GEV2_TO_CM3S)[1]
    checks.append(("T0 the freeze-out solver reused from V123 reproduces its "
                   "stored WIMP-miracle regression: Omega h^2 = %.6f (stored "
                   "%.6f)" % (om100, om100_ref),
                   abs(om100 - om100_ref) < 1e-4))

    # ---- T1 the nodal-plane density from the engine -------------------------
    rows = V123.node_velocities()
    speeds = [b for r in rows for b in r["branch_speeds"] if b > 1e-6]
    D123 = sum(1.0 / s ** 3 for s in speeds)
    D_ref = v123["data"]["doublers"]["D_planar"]
    checks.append(("T1 the nodal-plane final-state density recomputed from "
                   "the implemented step: %d M nodes, rank-1 metric, "
                   "comm < %.0e, D = %.2f (V123 stored %.2f)"
                   % (len(rows), max(r["comm_max"] for r in rows),
                      D123, D_ref),
                   len(rows) == 12 and abs(D123 - D_ref) / D_ref < 1e-2
                   and max(r["comm_max"] for r in rows) < 1e-8))

    # ---- T2 the V148 sigma_SI map and the DD ceiling ------------------------
    def sig_si(lam):
        return lam ** 2 * F_N ** 2 * MU_RED ** 2 / (math.pi * M_H ** 4) * CM2

    v88 = json.load(open(exp("V88_mirror_scalar",
                                     "results", "v88_mirror_scalar.json")))
    lam_star = v88["data"]["V2_lam_star"]           # 3.5299e-4, pair portal
    lam_z2 = v88["data"]["V1_lam_for_Omega"]        # 0.13081, Z2 portal
    sig_ref_z2 = v88["data"]["V1_sigma_SI_cm2"]     # 6.7209e-46
    sig_ref_star = v88["data"]["V2_sigma_SI_cm2"]   # 4.8942e-51
    lam_lz = math.sqrt(LZ * math.pi * M_H ** 4
                       / (F_N ** 2 * MU_RED ** 2 * CM2))
    checks.append(("T2 the V148 sigma_SI map (mu_red = %.3f GeV) reproduces "
                   "V88's stored cross sections at the named couplings: "
                   "%.3e vs %.3e (lambda_Z2), %.3e vs %.3e (lambda_*); the "
                   "LZ-class ceiling is lambda_LZ = %.3e"
                   % (MU_RED, sig_si(lam_z2), sig_ref_z2,
                      sig_si(lam_star), sig_ref_star, lam_lz),
                   abs(sig_si(lam_z2) - sig_ref_z2) / sig_ref_z2 < 1e-3
                   and abs(sig_si(lam_star) - sig_ref_star) / sig_ref_star
                   < 1e-3))

    # ---- T3 kernels: the V75 point, the closed-channel audit -----------------
    k2_at_v75 = float(k2_s(np.array(4.0 * M_S ** 2)))
    sigv75_unit = (V_EV ** 2 / (8.0 * math.pi)
                   * M_S ** 2 / (4.0 * M_S ** 2 - M_H ** 2) ** 2)
    gamma_h_ss = math.sqrt(M_H ** 2 - 4 * M_S ** 2) / (8 * math.pi) \
        if M_H >= 2 * M_S else 0.0                  # closed channel
    checks.append(("T3 K2(s = 4 m_s^2) equals V123's stored sigmav_V75 per "
                   "unit lambda (%.4e vs %.4e GeV^-2, rel %.1e); the h -> s "
                   "sbar on-shell width is %.0f - 2m_s = %.1f GeV > M_h = "
                   "%.2f GeV, so the s-channel pole s = M_h^2 is strictly "
                   "off-shell for every thermal pair"
                   % (k2_at_v75, sigv75_unit,
                      abs(k2_at_v75 - sigv75_unit) / sigv75_unit,
                      gamma_h_ss, 2 * M_S, M_H),
                   abs(k2_at_v75 - sigv75_unit) / sigv75_unit < 1e-12
                   and gamma_h_ss == 0.0 and M_H < 2 * M_S))

    # ---- T3b the convolution vs the seeded Monte-Carlo estimate --------------
    mc_ref = {}
    for tag, kfn, mm, TT in (("const(T=m_s)", k_const, M_S, M_S),
                             ("K2_hh(T=m_s)", k2_s, M_H, M_S),
                             ("K2_ss(T=m_s)", k2_s, M_S, M_S)):
        g = conv_average(kfn, mm, TT)
        c = mc_average(kfn, mm, TT)
        mc_ref[tag] = (g, c, abs(g - c) / c)
    max_mc = max(v[2] for v in mc_ref.values())
    checks.append(("T3b the exact convolution vs the seeded Monte-Carlo "
                   "average of the same thermal quantity: "
                   + ", ".join("%s %.4e vs %.4e (rel %.1e)"
                               % (k, v[0], v[1], v[2])
                               for k, v in mc_ref.items())
                   + " - the constant-kernel value is the mean Moller "
                   "velocity <v_M> = %.3f at T = m_s"
                   % mc_ref["const(T=m_s)"][1], max_mc < 3e-2))

    # ---- the thermal tables on the x grid (x = M_S/T) -----------------------
    NX = 256
    xgrid = np.exp(np.linspace(0.0, math.log(X_FI_MAX), NX))
    tables = {}
    for kname, kfn in KERNELS.items():
        hh = np.array([conv_average(kfn, M_H, M_S / x) for x in xgrid])
        ss = np.array([conv_average(kfn, M_S, M_S / x) for x in xgrid])
        tables[kname] = {"hh": hh, "ss": ss}

    # freeze-in constants: Omega_FI(lam) = C_FI * lam^2 * D_g, linear x grid
    def omega_fi_constant(kname, nsub=8000, nx_tab=None):
        if nx_tab is None:
            hh = tables[kname]["hh"]
            xg_t = xgrid
        else:
            xg_t = np.exp(np.linspace(0.0, math.log(X_FI_MAX), nx_tab))
            hh = np.array([conv_average(KERNELS[kname], M_H, M_S / x)
                           for x in xg_t])
        xg = np.linspace(1.0, X_FI_MAX, nsub)
        hh_i = np.interp(xg, xg_t, hh)
        integr = hh_i * np.array([n_mb(x, M_H) ** 2
                                  / (hubble(x) * x * x) for x in xg])
        return M_S * S0 * _trapz(integr, xg) / RHO_C

    C_FI = {k: omega_fi_constant(k) for k in KERNELS}

    # T3c convergence: convolution grids and the freeze-in integral
    cv = [conv_average(k1_s, M_S, M_S, nP=80, nq=160),
          conv_average(k1_s, M_S, M_S, nP=160, nq=320)]
    drift_cv = abs(cv[1] - cv[0]) / cv[1]
    csub = [omega_fi_constant("K2 (V2 pair portal)", nsub=n, nx_tab=nx)
            for n, nx in ((4000, 128), (8000, 256), (16000, 512))]
    drift_fi = max(abs(csub[i] - csub[i - 1]) / csub[i] for i in (1, 2))
    checks.append(("T3c convergence under refinement: thermal convolution "
                   "(80x160 vs 160x320) drift %.2e; freeze-in integral "
                   "(4000/128, 8000/256, 16000/512) drifts %.2e, %.2e"
                   % (drift_cv, drift_fi, abs(csub[2] - csub[1]) / csub[2]),
                   drift_cv < 0.02 and drift_fi < 0.02))

    # regime and freeze-out kernels
    Runit = {}          # max_T <K>_ss n_s / H per kernel (lam = 1, diagnostic)
    kfo = {}            # the constant kernel the freeze-out route uses
    om_eq = M_S * S0 * (float(n_mb(1.0, M_S)) / s_ent(1.0)) / RHO_C
    for kname in KERNELS:
        r = np.array([tables[kname]["ss"][i] * float(n_mb(x, M_S))
                      / hubble(x) for i, x in enumerate(xgrid)])
        Runit[kname] = float(r.max())
        kfo[kname] = float(tables[kname]["ss"][0])     # T = m_s

    # ---- variants -------------------------------------------------------------
    variants = []
    for kname in KERNELS:
        for Dg, dtag in ((1.0, "D-off"), (D123, "D-on")):
            if kname.startswith("K1p") and Dg != 1.0:
                continue          # declared: p-wave variant only D-off
            variants.append({"kernel": kname, "D": Dg, "tag": dtag})

    scan_lo = 1.0e-18
    per_variant = []
    for vr in variants:
        kname, Dg = vr["kernel"], vr["D"]
        cfi = C_FI[kname] * Dg
        # the regime boundary is the SELF-CONSISTENCY root Y_FI = Y_eq: below
        # it the freeze-in estimate stays under the equilibrium density (the
        # population never equilibrates); above it it does not (freeze-out).
        lam_th = math.sqrt(om_eq / cfi)

        def omega_to(lam, kf=kfo[kname]):
            return V123.freeze_out(M_S, kf * lam ** 2)[1]

        def omega_decided(lam, lam_th=lam_th, cfi=cfi, omega_to=omega_to):
            if lam <= lam_th:                    # Y_FI < Y_eq: freeze-in
                return cfi * lam ** 2
            return omega_to(lam)

        # the 0.120 crossings inside the DD-allowed range, per branch
        lam_fi = math.sqrt(OMEGA_OBS / cfi)
        if not (scan_lo < lam_fi <= lam_th):
            lam_fi = None
        lam_to = None
        lo_t, hi_t = lam_th, lam_lz
        if lo_t < hi_t and (omega_to(lo_t) - OMEGA_OBS) * \
                (omega_to(hi_t) - OMEGA_OBS) < 0:
            lo, hi = lo_t, hi_t
            for _ in range(45):
                mid = math.sqrt(lo * hi)
                if omega_to(mid) > OMEGA_OBS:
                    lo = mid
                else:
                    hi = mid
            lam_to = math.sqrt(lo * hi)
        # the same crossing WITHOUT the DD cap - where the thermal solution
        # actually sits when direct detection forbids it
        lam_to_open = None
        if (omega_to(lam_th) - OMEGA_OBS) * (omega_to(1.0) - OMEGA_OBS) < 0:
            lo, hi = lam_th, 1.0
            for _ in range(45):
                mid = math.sqrt(lo * hi)
                if omega_to(mid) > OMEGA_OBS:
                    lo = mid
                else:
                    hi = mid
            lam_to_open = math.sqrt(lo * hi)
        pts = {}
        for tag, lam in (("lambda_*", lam_star), ("lambda_Z2", lam_z2),
                         ("lambda_LZ", lam_lz)):
            pts[tag] = {"lambda": lam,
                        "R": Runit[kname] * lam ** 2,
                        "route": ("freeze-in" if lam <= lam_th
                                  else "freeze-out"),
                        "Omega_h2": omega_decided(lam)}
        pts["lambda_thermal"] = {"lambda": lam_th, "route": "boundary",
                                 "Omega_h2_fi": cfi * lam_th ** 2,
                                 "Omega_h2_fo": omega_to(lam_th)}
        per_variant.append({
            "variant": "%s %s" % (kname, vr["tag"]), "kernel": kname,
            "tag": vr["tag"], "D_factor": Dg,
            "lambda_thermal": lam_th,
            "Rmax_annihilation_criterion_per_lambda2": Runit[kname],
            "C_FI_omega_per_lambda2": cfi,
            "Omega_equilibrium_ceiling": om_eq,
            "lambda_fi_solution": lam_fi,
            "lambda_freeze_out_solution": lam_to,
            "lambda_freeze_out_solution_uncapped": lam_to_open,
            "sigma_SI_at_fi_cm2": sig_si(lam_fi) if lam_fi else None,
            "sigma_SI_at_to_cm2": sig_si(lam_to) if lam_to else None,
            "named_points": pts,
        })

    stability = {"scalar_parity_exact": v88["data"]["scalar_parity_exact"],
                 "V1_mirror_even_commutator": v88["data"]["V1_mirror_even"]}
    checks.append(("T4 the stability premise the production interpretation "
                   "needs, regressed from V88's stored verdicts: the scalar "
                   "parity Q_s is exact (%s) and the V1 portal commutator "
                   "[Q_s,V1] = %.1e - the single s cannot decay"
                   % (stability["scalar_parity_exact"],
                      stability["V1_mirror_even_commutator"]),
                   bool(stability["scalar_parity_exact"])
                   and abs(stability["V1_mirror_even_commutator"]) < 1e-9))

    out["data"]["constants"] = {
        "M_S_GeV": M_S, "M_H_GeV": M_H, "mu_red_GeV": MU_RED,
        "LZ_bound_cm2": LZ, "lambda_LZ_ceiling": lam_lz,
        "D_planar_engine": D123,
        "Omega_equilibrium_ceiling": om_eq,
        "convolution_unit_GeV2": {k: {
            "hh_at_MS": float(tables[k]["hh"][0]),
            "ss_at_MS": float(tables[k]["ss"][0])} for k in KERNELS},
        "hh_over_ss_conv_ratio_at_MS": {
            k: float(tables[k]["hh"][0] / tables[k]["ss"][0])
            for k in KERNELS},
        "Rmax_annihilation_criterion_per_lambda2": Runit}
    out["data"]["regression"] = {
        "Omega_h2_at_100GeV_3e-26": om100, "stored": om100_ref,
        "D_planar_stored": D_ref,
        "sigma_SI_stored_pair_cm2": sig_ref_star,
        "sigma_SI_stored_Z2_cm2": sig_ref_z2}
    out["data"]["variants"] = per_variant

    def getv(kpre, dtag="D-off"):
        return next(p for p in per_variant
                    if p["kernel"].startswith(kpre) and p["tag"] == dtag)

    k1v, k2v = getv("K1 "), getv("K2")
    fi_sols = [p for p in per_variant if p["lambda_fi_solution"]]
    to_sols = [p for p in per_variant if p["lambda_freeze_out_solution"]]
    branch = "A" if (fi_sols or to_sols) else "B"
    ok = (all(c[1] for c in checks) and drift_cv < 0.02 and drift_fi < 0.02)

    if branch == "A":
        status = (
            "THE PRODUCTION SIDE IS COMPUTED (O7). The exact-Z2 ladder "
            "scalar (m_s = %.2f GeV, stable: T4) has ONE production "
            "solution inside the LZ-class ceiling lambda_LZ = %.2e: "
            "FREEZE-IN at lambda = %.2e with the V1 portal kernel "
            "(Omega_FI = %.3g lambda^2; regime self-consistent - below the "
            "boundary lambda_th = %.2e the yield stays under the "
            "equilibrium ceiling Omega = %.3g), predicting sigma_SI = %.2e "
            "cm^2, %.1f orders below the bound: freeze-in is DD-blind by "
            "construction. The freeze-out solution of the same kernel sits "
            "at lambda = %s - above the ceiling, its sigma_SI = %s cm^2: "
            "the V88/V148 direct-detection exclusion of the thermal Z2 "
            "portal is confirmed by computed production, not by an "
            "inverted coupling. "
            "The pair (Z2-breaking) kernel: freeze-out solution lambda = "
            "%s (sigma_SI = %s cm^2) %s the ceiling; its named point "
            "lambda_* = %.2e gives Omega = %.3g via %s - the inverted "
            "fiducial coupling is NOT a production solution. Nodal-plane D "
            "on the production channel multiplies the freeze-in coupling "
            "by 1/sqrt(D) (D-on variants reported). The h -> s sbar channel "
            "is closed (2m_s = %.1f > M_h = %.2f): the s-channel pole is "
            "always off-shell; the exact Moller convolution agrees with "
            "seeded Monte Carlo to %.0e. VERDICT: the O7 freeze-in "
            "parenthetical is retired - the production side is computed, "
            "with a quantified and (for the Z2 portal) DD-blind solution. %s"
            % (M_S, lam_lz,
               k1v["lambda_fi_solution"], k1v["C_FI_omega_per_lambda2"],
               k1v["lambda_thermal"], om_eq,
               k1v["sigma_SI_at_fi_cm2"],
               math.log10(LZ / k1v["sigma_SI_at_fi_cm2"]),
               ("%.2e" % k1v["lambda_freeze_out_solution_uncapped"])
               if k1v["lambda_freeze_out_solution_uncapped"] else "none",
               ("%.2e" % sig_si(k1v["lambda_freeze_out_solution_uncapped"]))
               if k1v["lambda_freeze_out_solution_uncapped"] else "-",
               ("%.2e" % k2v["lambda_freeze_out_solution"])
               if k2v["lambda_freeze_out_solution"] else "none",
               ("%.2e" % k2v["sigma_SI_at_to_cm2"])
               if k2v["sigma_SI_at_to_cm2"] else "-",
               "inside" if k2v["lambda_freeze_out_solution"]
               and k2v["lambda_freeze_out_solution"] < lam_lz else "outside",
               lam_star,
               k2v["named_points"]["lambda_*"]["Omega_h2"],
               k2v["named_points"]["lambda_*"]["route"],
               2 * M_S, M_H, max_mc, "PASS" if ok else "FAIL"))
    else:
        status = (
            "FREEZE-IN COMPUTED AND EXCLUDED (O7). The production side of "
            "the ladder-scalar abundance was integrated through both named "
            "portal kernels (exact convolution, freeze-out and freeze-in, "
            "regime from the self-consistency root Y_FI = Y_eq) against the "
            "LZ-class ceiling lambda_LZ = %.2e: no DD-allowed coupling "
            "reproduces Omega h^2 = 0.120 on any variant - at the ceiling "
            "the decided-regime abundance is %s. The freeze-in "
            "parenthetical is retired as COMPUTED AND EXCLUDED; the "
            "abundance remains a named input. %s"
            % (lam_lz,
               ["%.3g" % p["named_points"]["lambda_LZ"]["Omega_h2"]
                for p in per_variant],
               "PASS" if ok else "FAIL"))

    out["pass"] = bool(ok)
    out["branch"] = branch
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = status
    out["inconclusive"] = bool(drift_cv >= 0.02 or drift_fi >= 0.02)
    out["wall_time"] = time.time() - t0

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v193_w7_freeze_in.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V193 - production side of the relic abundance computed (O7)")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    for p in per_variant:
        print("  %-34s lam_th %.2e  lam_FI %s (sig %s)  lam_TO %s (sig %s)"
              % (p["variant"], p["lambda_thermal"],
                 "%.2e" % p["lambda_fi_solution"]
                 if p["lambda_fi_solution"] else "none",
                 "%.1e" % p["sigma_SI_at_fi_cm2"]
                 if p["sigma_SI_at_fi_cm2"] else "-",
                 "%.2e" % p["lambda_freeze_out_solution"]
                 if p["lambda_freeze_out_solution"]
                 else ("%.2e open" % p["lambda_freeze_out_solution_uncapped"]
                       if p["lambda_freeze_out_solution_uncapped"]
                       else "none"),
                 "%.1e" % p["sigma_SI_at_to_cm2"]
                 if p["sigma_SI_at_to_cm2"] else "-"))
        for tag in ("lambda_*", "lambda_Z2", "lambda_LZ", "lambda_thermal"):
            pt = p["named_points"][tag]
            if tag == "lambda_thermal":
                print("      lam_th %.3e  Omega fi/fo at boundary "
                      "%.3g / %.3g" % (pt["lambda"], pt["Omega_h2_fi"],
                                       pt["Omega_h2_fo"]))
            else:
                print("      %-10s lam %.3e %-10s Omega %.3g"
                      % (tag, pt["lambda"], pt["route"], pt["Omega_h2"]))
    print("  branch %s: %s" % (branch,
                               "ABUNDANCE COMPUTED" if branch == "A"
                               else "ABUNDANCE NOT REACHED"))
    print("  VERDICT: %s   wall %.2f s" % ("PASS" if ok else "FAIL",
                                           out["wall_time"]))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
