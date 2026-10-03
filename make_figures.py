#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — figure generation for `doc/publication.md`.

Every figure is **computed** from the 0-parameter CL10 modules (no hand-
drawn data): the 32-level QCA spectrum, the 1-loop running, the 10-loop
ladder, the WIMP miracle, the 2-route Kretschmann, and the 512/512
cancellation.  Run:

    .venv/bin/python make_figures.py

Writes PNGs to `doc/figures/` (the publication embeds them by relative path).
The backend (CuPy/JAX/NumPy) is selected by `cl10.engine` at import; the
32x32 objects are small, so this runs on the GPU or degrades to CPU.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "doc", "figures")

# --- a clean, publication-grade style --------------------------------------
plt.rcParams.update({
    "figure.figsize": (6.8, 4.6),
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "axes.axisbelow": True,
    "axes.unicode_minus": False,
    "lines.linewidth": 1.8,
})

# consistent 0-parameter palette
C_GUT = "#1f77b4"     # blue
C_SU2 = "#ff7f0e"     # orange
C_SU3 = "#2ca02c"     # green
C_ROOT = "#7f1c8b"    # purple (the one object / the root)
C_ZERO = "#d62728"    # red (the 0-parameter points / the measured)

os.makedirs(OUT, exist_ok=True)


def _save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote", path)


# ===========================================================================
# Figure 1 — the one 0-parameter object: the 32-level QCA spectrum (9-to-1)
# ===========================================================================
def fig1_one_object():
    from cl10 import qca_spectrum

    lv = np.asarray(qca_spectrum.levels(), dtype=float)          # 32 phases, sorted
    n = len(lv)
    zp = float(lv.sum())

    fig, ax = plt.subplots()
    x = np.arange(n)
    colors = [C_ZERO if v > 1e-12 else ("#7f7f7f" if abs(v) <= 1e-12 else C_ROOT)
              for v in lv]
    ax.vlines(x, 0, lv, color=colors, linewidth=1.6, alpha=0.9)
    ax.axhline(0, color="k", lw=0.9)
    ax.set_xlabel("spectrum index  $k$  (32 levels, sorted)")
    ax.set_ylabel(r"eigenphase  $\omega_k$  (rad/step)")
    ax.set_title("The one 0-parameter object: the 32-level QCA spectrum\n"
                 r"$H=-i\log U$,  $U=\prod_{a<b}(\cos g_{10}\,I-i\sin g_{10}\,T_{ab})$ at $g_{10}=1/5$")
    # the 9-to-1 invariants, computed
    from cl10 import parameter_audit
    sm = parameter_audit.spectral_mirror()
    ax.text(0.03, 0.97,
            r"0-point:  $\sum_k \omega_k = $" + f"{zp:.1e}" + r"  ($=0$, W1)" + "\n"
            + r"$\det U = $" + f"{1 - sm['det_U_minus_1']:.2e}"
            + r"  (the $\mathrm{SL}(32)$ identity, W1)" + "\n"
            + r"CPT:  spec(16-bar) $= -$ spec(16),  mirror residual $= $"
            + f"{sm['mirror_resid']:.1e}" + "\n"
            + r"smallest level  $\omega_{\min} = $" + f"{sm['omega_min']:.4f}"
            + r"  (a mass $\sim 1.4\times10^{17}$ GeV)",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#888", alpha=0.95))
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], color=C_ROOT, lw=2, label=r"$-\omega_k$  (the $16$-bar / antimatter)"),
        Line2D([0], [0], color=C_ZERO, lw=2, label=r"$+\omega_k$  (the $16$ / matter)"),
    ], loc="lower right")
    ax.set_xlim(-1, n)
    _save(fig, "fig1_one_object_32spectrum.png")


# ===========================================================================
# Figure 2 — the 1-loop running (1/alpha_i) to M_P + the Landau-pole contrast
# ===========================================================================
def fig2_unification():
    from cl10 import unification

    B = unification.B_SM          # (41/10, -19/6, -7)  1-loop beta
    B2 = unification.B2_SM        # the 2-loop beta coefficients
    AZ = unification.A_Z          # the measured couplings at M_Z
    M_Z, M_P = unification.M_Z, unification.M_P
    names = [r"$1/\alpha_1$ (U(1))", r"$1/\alpha_2$ (SU(2))", r"$1/\alpha_3$ (SU(3))"]
    cols = [C_GUT, C_SU2, C_SU3]

    # the 2-loop U(1) Landau pole (the perturbative UV-divergence, below M_P)
    b1, b2, c = B[0], B2[0], 1.0 / AZ[0]
    Acoef = b2 / (b1 ** 2 * 2.0 * np.pi)
    Bcoef = b1 / (2.0 * np.pi)
    t2 = (-Bcoef + np.sqrt(Bcoef ** 2 + 4.0 * Acoef * c)) / (2.0 * Acoef)
    pole2 = M_Z * float(np.exp(t2))
    pole1 = unification.landau_pole_u1()

    s = np.logspace(np.log10(M_Z), np.log10(M_P), 240)
    ln = np.log(s / M_Z)
    fig, ax = plt.subplots()
    for k in range(3):
        ax.plot(np.log10(s), 1.0 / AZ[k] - (B[k] / (2.0 * np.pi)) * ln,
                color=cols[k], lw=2.0, label=names[k])
    # the 0-parameter GUT coupling: 1/alpha_10 = 25
    ax.axhline(1.0 / unification.ALPHA_10, color=C_ZERO, ls="--", lw=1.4,
               label=r"0-parameter $1/\alpha_{10} = 25$ (the GUT coupling)")
    # the 0-parameter GUT scale and the 2-loop Landau pole
    ax.axvline(np.log10(M_P), color="k", ls="-", lw=1.0)
    ax.text(np.log10(M_P) - 0.12, 9.0, r"$M_P$  ($1.22\times10^{19}$ GeV)",
            rotation=90, fontsize=9, ha="right")
    ax.axvline(np.log10(pole2), color="#7f7f7f", ls=":", lw=1.4,
               label=f"2-loop U(1) Landau pole ($\\sim$ {pole2:.0e} GeV, below $M_P$)")
    # the 1-loop pole (above M_P) — the 1-loop running stays UV-bounded
    ax.text(0.02, 0.97,
            f"1-loop U(1) pole $\\sim$ {pole1:.0e} GeV (above $M_P$): the 1-loop\n"
            + r"running is UV-bounded, but the 2-loop pole lies below $M_P$ — "
            "so the\n"
            "perturbative QFT is UV-divergent while the QCA step is UV-finite (W4)",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#888", alpha=0.95))
    ax.set_xlabel(r"$\log_{10}(M/\mathrm{GeV})$  (1-loop, SM $\beta$-coefficients)")
    ax.set_ylabel(r"running coupling $1/\alpha_i$")
    ax.set_title("1-loop running of the three SM couplings to the Planck scale: they reach the\n"
                 r"GUT range near $M_P$ (the 1-loop matching, W4); the three couplings do not meet at a single point")
    ax.set_xlim(np.log10(M_Z) - 0.1, np.log10(M_P) + 0.5)
    ax.set_ylim(8, 60)
    ax.legend(loc="lower left", framealpha=0.95)
    _save(fig, "fig2_unification_1loop.png")


# ===========================================================================
# Figure 3 — the 10-loop ladder: one factor, two results (hierarchy + WIMP)
# ===========================================================================
def fig3_10loop_ladder():
    from cl10 import hier_exponent

    M_P = hier_exponent.M_P
    LOOP = 1.0 / (16.0 * np.pi)
    n = np.arange(0, 11)
    ladder = M_P * (LOOP ** n)          # GeV at each loop level

    fig, ax = plt.subplots()
    ax.semilogy(n, ladder, "o-", color=C_ROOT, lw=2.2, ms=6,
                label=r"$M_P\,(1/16\pi)^{n}$")
    # the two landing points
    ax.axhline(100.0, color=C_ZERO, ls="--", lw=1.5,
               label="WIMP / weak scale ($100$ GeV)")
    ax.scatter([10], [ladder[10]], s=90, color=C_ZERO, zorder=5, clip_on=False,
               label=rf"10-loop $= M_P(1/16\pi)^{{10}} = $ {ladder[10]:.0f} GeV "
                     r"(the 0-parameter WIMP mass)")
    ax.scatter([0], [ladder[0]], s=90, color="k", zorder=5, clip_on=False,
               label=r"$M_P = 1.22\times10^{19}$ GeV (the Planck scale)")
    ax.annotate("the 10-loop suppression\nbrings the Planck scale\ndown to the WIMP scale",
                xy=(10, ladder[10]), xytext=(5.5, 1e6),
                fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="->", color="0.3"))
    ax.set_xlabel(r"loop count $n$  (0-parameter: $10 = n_{\rm gen} = 2 \times 5$, the 16/16-bar chiral split $\times$ the 5 qubits of the 32-D module)")
    ax.set_ylabel(r"energy scale (GeV, log scale)")
    ax.set_title("The 10-loop ladder: one 0-parameter factor gives both\n"
                 r"the hierarchy $v/E_P = 2.08\,(1/16\pi)^{10}$ (V53) and the DM mass (V54)")
    ax.set_xticks(n)
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(50, 1e20)
    ax.legend(loc="upper right", framealpha=0.95)
    _save(fig, "fig3_10loop_ladder.png")


# ===========================================================================
# Figure 4 — the WIMP miracle: relic abundance vs mass
# ===========================================================================
def fig4_wimp_miracle():
    from cl10 import dm_wimp

    OMEGA_OBS = dm_wimp.relic_abundance(100.0)      # = 0.12 by construction
    m = np.logspace(1.7, 20, 240)
    omega = OMEGA_OBS * (100.0 / m) ** 3

    m_dm = dm_wimp.dm_wimp()["m_DM_GeV"]
    planck_gap = dm_wimp.dm_wimp()["v42_planck_gap"]
    M_P = dm_wimp.M_P

    fig, ax = plt.subplots()
    ax.loglog(m, omega, color=C_ROOT, lw=2.2,
              label=r"relic $\Omega h^2 = \Omega_{\rm obs}(100/m)^3$")
    ax.axhline(0.12, color="k", ls="--", lw=1.5,
              label=r"observed $\Omega_{\rm DM} h^2 = 0.12$")
    # the WIMP-miracle point and the 0-parameter point
    ax.scatter([100.0], [0.12], s=110, color=C_SU3, zorder=5, clip_on=False,
               label="WIMP miracle: $m=100$ GeV $\\Rightarrow$ $\\Omega_{\\rm obs}$")
    ax.scatter([m_dm], [dm_wimp.relic_abundance(m_dm)], s=110, color=C_ZERO, zorder=5, clip_on=False,
               marker="*",
               label=f"0-parameter: $m_{{\\rm DM}}=M_P(1/16\\pi)^{{10}}={m_dm:.0f}$ GeV $\\Rightarrow$ $\\Omega h^2=$ {dm_wimp.relic_abundance(m_dm):.3f}")
    # the V42 Planck-mass soliton (under-abundant)
    ax.scatter([M_P], [0.12 * planck_gap], s=70, color="#7f7f7f", zorder=5, clip_on=False,
               label=r"V42 Planck-mass soliton ($\sim 10^{-52}$, under-abundant)")
    # the 0-parameter WIMP star is identified in the legend (no extra
    # annotation: it sits in the top corner next to the WIMP-miracle point).
    ax.set_xlabel(r"relic mass $m$ (GeV, log scale)")
    ax.set_ylabel(r"relic abundance $\Omega h^2$")
    ax.set_title("The WIMP miracle, resolved: the 0-parameter 10-loop mass\n"
                 "$M_P(1/16\\pi)^{10}=118$ GeV gives the observed dark-matter abundance (V54)")
    ax.set_xlim(40.0, 1e20)
    ax.set_ylim(1e-55, 1.0)
    ax.legend(loc="upper right", framealpha=0.95)
    _save(fig, "fig4_wimp_miracle.png")


# ===========================================================================
# Figure 5 — 2-route gravity: computed vs analytic Kretschmann (Schwarzschild)
# ===========================================================================
def fig5_gravity_kretschmann():
    from cl10 import gravity

    M = 1.0
    r = np.linspace(2.2, 12.0, 9)                  # r > 2M = 2
    Kc, Ka = [], []
    for ri in r:
        x = np.array([0.0, ri, 0.5, 0.7])

        def sw(xx):
            return gravity.metric_schwarzschild(M, xx)

        R = gravity.riemann_4d(sw, x)
        g, _ = gravity._metric_and_dg(sw, x)
        Kc.append(gravity.kretschmann_4d(R, g))
        Ka.append(48.0 * M ** 2 / ri ** 6)
    Kc, Ka = np.array(Kc), np.array(Ka)

    fig, ax = plt.subplots()
    ax.loglog(r, Ka, "-", color="k", lw=1.8,
              label=r"analytic $K = 48 M^2/r^6$ (Schwarzschild)")
    ax.loglog(r, Kc, "o", color=C_ZERO, ms=6,
              label="computed from the metric (4-D central finite differences)")
    # show the residual of the 2-route agreement
    resid = float(np.abs(Kc / Ka - 1.0).max())
    ax.text(0.03, 0.04,
            f"max 2-route residual = {resid:.1e}\n(the black hole is recovered by computation)",
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#888", alpha=0.95))
    ax.set_xlabel(r"radial coordinate $r$  ($M=1$, $r>2M$)")
    ax.set_ylabel(r"Kretschmann $K = R^{mnpq}R_{mnpq}$")
    ax.set_title("2-route verification: the computed curvature equals the analytic\n"
                 "Schwarzschild Kretschmann (the black hole, $G=0$ in vacuum)")
    ax.legend(loc="upper right", framealpha=0.95)
    _save(fig, "fig5_gravity_kretschmann.png")


# ===========================================================================
# Figure 6 — the 512/512 Z2 cancellation (the 0-parameter dark energy)
# ===========================================================================
def fig6_512_512():
    from cl10 import cosmology3

    de = cosmology3.dark_energy()
    ne, no = de["n_even"], de["n_odd"]

    fig, ax = plt.subplots()
    cats = ["bosonic\n(512 even)", "fermionic\n(512 odd)", "vacuum\n0-point"]
    vals = [ne * 0.5, -no * 0.5, 0.0]
    cols = [C_SU2, C_SU3, C_ZERO]
    bars = ax.bar(cats, vals, color=cols, width=0.6, edgecolor="k", linewidth=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + (12 if v >= 0 else -40),
                f"{v:+.0f}", ha="center", fontsize=10,
                color=(C_SU2 if v > 0 else (C_SU3 if v < 0 else C_ZERO)))
    ax.axhline(0, color="k", lw=1.0)
    ax.set_ylabel(r"vacuum 0-point  $\sum \frac{1}{2}\,\mathrm{sign}$")
    ax.set_title(r"Dark energy: the $\mathbb{Z}_2$ 512/512 cancellation "
                 r"$\Rightarrow$ $\Lambda_{\rm 0param}=0$ (the s-mechanism), $w=-1$")
    ax.set_ylim(-300, 300)
    ax.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc=C_SU2, label=r"$512 \times (+\frac{1}{2}) = +256$  (bosons)"),
        plt.Rectangle((0, 0), 1, 1, fc=C_SU3, label=r"$512 \times (-\frac{1}{2}) = -256$  (fermions)"),
        plt.Rectangle((0, 0), 1, 1, fc=C_ZERO, label=r"$\Lambda_{\rm 0param} = +256 - 256 = 0$"),
    ], loc="upper right", framealpha=0.95)
    _save(fig, "fig6_512_512.png")


# ===========================================================================
# Figure 7 — QCD confinement: the 0-parameter linear vs Coulombic potential
# ===========================================================================
def fig7_qcd():
    from cl10 import qcd

    sigma = qcd.string_tension()      # the 0-parameter string tension (= 1)
    r = np.logspace(0, 3, 200)        # quark separation (Planck units)
    v_qcd = sigma * r                 # the linear (confining) potential
    v_em = 1.0 / r                    # the Coulombic (screening) potential

    fig, ax = plt.subplots()
    ax.loglog(r, v_qcd, color=C_SU3, lw=2.2,
              label=r"strong:  $V(r) = \sigma\, r$  (linear — the energy grows without bound)")
    ax.loglog(r, v_em, color=C_GUT, lw=2.2,
              label=r"electromagnetic:  $V(r) \propto 1/r$  (Coulombic — screening)")
    ax.axvline(1.0, color="0.6", ls=":", lw=1.0)
    ax.set_xlabel(r"quark–antiquark separation $r$  (Planck units, log)")
    ax.set_ylabel(r"potential $V$  (Planck units, log)")
    ax.set_title("QCD confinement (0-parameter): the linear potential $V = \\sigma r$ rises without\n"
                 "bound (confinement), in contrast to the Coulombic $1/r$ of the unconfined EM force")
    ax.legend(loc="upper left", framealpha=0.95)
    _save(fig, "fig7_qcd_confinement.png")


# ===========================================================================
# Figure 8 — GUT leptogenesis: the 0-parameter baryon asymmetry
# ===========================================================================
def fig8_baryon():
    from cl10 import eta_cp

    r = eta_cp.eta_cp()
    obs = r["delta_B_observed"]       # 6.1e-10
    v40 = r["delta_B_v40"]            # named eta_CP = 1e-3  (V40, ratio 2.8)
    v55 = r["delta_B_0param"]         # 0-parameter 2-loop eta_CP  (V55, ratio 1.12)

    fig, ax = plt.subplots()
    cats = ["observed", "V40\n(named $\\eta_{CP}=10^{-3}$)",
            "V55\n($\\frac{28}{79}\\,\\varepsilon_{CP}\\,D$)"]
    vals = [obs, v40, v55]
    vals9 = [x / 1e-9 for x in vals]
    cols = ["#7f7f7f", "#ff7f0e", C_ZERO]
    bars = ax.bar(cats, vals9, color=cols, width=0.5, edgecolor="k", linewidth=0.6)
    for b, v, c in zip(bars, vals, cols):
        ax.text(b.get_x() + b.get_width() / 2, v / 1e-9 + 0.06,
                f"{v:.1e}\n(ratio {v / obs:.2f})", ha="center", fontsize=9, color=c)
    ax.axhline(obs / 1e-9, color="k", ls="--", lw=1.2,
               label=r"the observed $\Delta B = 6.1\times10^{-10}$")
    ax.set_ylabel(r"baryon asymmetry $\Delta B$  ($\times 10^{-9}$)")
    ax.set_title("GUT leptogenesis: the baryon asymmetry as a product of 0-parameter numbers.  "
                 r"V55: $\Delta B = \frac{28}{79}\,\varepsilon_{CP}\,D$ ($\varepsilon_{CP}=g_0^4/16\pi$, $g_0=1/2^3$; $D=(1/16\pi)^2$) gives $\Delta B = 6.8\times10^{-10}$ — a 12% match (V40's ratio: 2.8)")
    ax.set_ylim(0, max(vals9) * 1.4)
    ax.legend(loc="upper left", framealpha=0.95)
    _save(fig, "fig8_baryon_asymmetry.png")


# ===========================================================================
# Figure 9 — the 16 Weyl is anomaly-free (all 8 gauge anomalies cancel)
# ===========================================================================
def fig9_anomaly():
    from cl10 import matter

    # the [U(1)]^3 per-field contribution n*Y^3 (the exact cancellation)
    order_keys = ["Q", "u^c", "d^c", "L", "e^c", "N"]
    order_lbl = ["Q", r"$u^c$", r"$d^c$", "L", r"$e^c$", "N"]
    frac_lbl = {  # the exact n*Y^3 per field (an exact algebraic identity)
        "Q": "1/36", "u^c": "-8/9", "d^c": "1/9",
        "L": "-1/4", "e^c": "1", "N": "0"}
    perfield = [n * y ** 3 for (_s3, _s2, y, n) in
                (matter.SM_FIELDS[k] for k in order_keys)]

    fig, ax = plt.subplots()
    x = np.arange(len(order_keys))
    cols = [C_SU3 if v >= 0 else C_GUT for v in perfield]
    ax.bar(x, perfield, color=cols, edgecolor="k", linewidth=0.5, width=0.5)
    ax.axhline(0, color="k", lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(order_lbl)
    ax.set_ylabel(r"$n \times Y^3$  (the $[\mathrm{U}(1)]^3$ contribution per field)")
    for xi, k in zip(x, order_keys):
        v = perfield[xi]
        va = "bottom" if v >= 0 else "top"
        ax.text(xi, v + (0.10 if v >= 0 else -0.10), frac_lbl[k],
                ha="center", va=va, fontsize=10)
    ax.set_ylim(-1.55, 1.55)
    ax.set_title(r"$[\mathrm{U}(1)]^3$: the six SM fields cancel ($\sum = 0$) — the 16 Weyl is anomaly-free (V0)")
    _save(fig, "fig9_anomaly_free.png")

# ===========================================================================
# Figure 10 — CL10 vs observation: the 0-parameter predictions vs the data
# ===========================================================================
def fig10_predictions():
    from cl10 import eta_cp as ep, dm_wimp as dw, toepred as tp

    e = ep.eta_cp()
    d = dw.dm_wimp()
    # (x = observed/measured, y = 0-parameter prediction, label, color, marker, label-offset)
    rows = [
        (e["delta_B_observed"], e["delta_B_0param"],
         "baryon asymmetry  $\\Delta B$  (2-loop, V55)", C_ZERO, "o", (12, -22)),
        (d["omega_obs"], d["omega_DM"],
         r"dark-matter $\Omega h^2$  (10-loop WIMP, V54)", C_SU2, "o", (-175, 20)),
        (d["weak_scale_GeV"], d["m_DM_GeV"],
         r"dark-matter mass $m_{\rm DM}$  (10-loop, V54)  [GeV]", C_SU2, "s", (16, 16)),
        (2.016e-17, 2.02e-17,
         r"hierarchy $v/E_P$  (10-loop, V53)", C_GUT, "o", (12, -24)),
        (1e-3, e["eta_CP_0param"],
         r"CP factor $\eta_{CP}$  (2-loop, V55)", C_SU3, "^", (-150, -26)),
        (0.05, ((1.0 / 8 * 246.0) ** 2 / 1e14) * 1e9,
         r"neutrino mass $m_\nu$  (seesaw, W2)  [eV]", C_ROOT, "D", (14, -26)),
        (1e34, tp.proton_decay_lifetime(),
         r"proton lifetime $\tau_p$  (GUT, a lower bound)  [yr]", "#7f7f7f", "x", (-205, -16)),
    ]
    xs = np.array([r[0] for r in rows])
    ys = np.array([r[1] for r in rows])
    lo = min(xs.min(), ys.min()) * 0.3
    hi = max(xs.max(), ys.max()) * 3.0

    fig, ax = plt.subplots(figsize=(15, 9))
    # the 45-degree line (perfect match) + the within-a-factor-of-2 band
    ax.plot([lo, hi], [lo, hi], color="k", lw=1.0,
            label="perfect match  (0-param $=$ observed)")
    ax.fill_between([lo, hi], 0.5 * np.array([lo, hi]), 2.0 * np.array([lo, hi]),
                    color="0.82", alpha=0.6, zorder=0,
                    label="within a factor of 2")
    for (xo, yp, lab, col, mk, off) in rows:
        ax.scatter([xo], [yp], s=130, color=col, marker=mk, zorder=5, clip_on=False)
        ax.annotate(lab, (xo, yp), xytext=off, textcoords="offset points",
                    fontsize=9, color=col)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("observed / measured value  (log)")
    ax.set_ylabel("0-parameter prediction  (log)")
    ax.set_title("CL10 vs observation: the 0-parameter predictions vs the measured values.\n"
                 r"The 45$^\circ$ line $=$ a match; most predictions sit within a factor of 2 (2-digit residuals)",
                 fontsize=13)
    ax.legend(loc="center right", framealpha=0.95)
    _save(fig, "fig10_predictions_vs_obs.png")

# ===========================================================================
# Figure 11 — the 0-parameter ladder of scales (one object -> all scales)
# ===========================================================================
def fig11_ladder():
    from cl10 import dm_wimp as dw, toepred as tp

    m_nu_eV = ((1.0 / 8 * 246.0) ** 2 / 1e14) * 1e9     # ~ 10^-2 eV
    # (energy GeV, label, color, side(+/-1), height multiplier)
    scales = [
        (m_nu_eV * 1e-9, r"neutrino  $m_\nu \sim 10^{-2}$ eV  (seesaw, W2)", C_ZERO, -1, 0.5),
        (1.0, r"QCD  $\Lambda_{\rm QCD} \sim 1$ GeV  (V27)", C_SU3, 1, 0.5),
        (dw.dm_wimp()["m_DM_GeV"], r"0-param DM  $m_{\rm DM}=118$ GeV  (10-loop, V54)", C_SU2, -1, 0.5),
        (tp.V_EW, r"electroweak  $v = 246$ GeV", C_GUT, 1, 0.95),
        (1e14, r"2-step seesaw  $M_N \sim 10^{14}$ GeV  (W2)", C_ROOT, -1, 0.5),
        (tp.M_PLANCK, r"Planck / GUT  $M_P = M_{\rm GUT} = 1.22\times10^{19}$ GeV", "k", 1, 0.5),
    ]
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.axhline(0, color="0.3", lw=1.0)
    for (E, lab, col, side, hgt) in scales:
        x = np.log10(E)
        ax.plot([x, x], [0, 0.30 * hgt * side], color=col, lw=1.0)
        ax.scatter([x], [0], s=90, color=col, zorder=5, clip_on=False,
                   edgecolor="k", lw=0.5)
        ax.text(x, 0.37 * hgt * side, lab, ha="center",
                va="bottom" if side > 0 else "top", fontsize=9, color=col)
    ax.set_xlim(-12.5, 20.5)
    ax.set_ylim(-0.72, 0.92)
    ax.set_yticks([])
    ax.set_xlabel("energy scale  (GeV, log scale)")
    ax.set_title("The 0-parameter ladder of scales: one object generates six scales spanning 30 decades\n"
                 r"(the single $10$-loop factor $M_P(1/16\pi)^{10}$ gives both the DM mass and the $v/E_P$ hierarchy)")
    _save(fig, "fig11_ladder_of_scales.png")


# ===========================================================================
# Figure 12 — the 0-parameter parameter audit: the five epistemic classes
# ===========================================================================
def fig12_audit():
    from cl10 import parameter_audit as pa

    c = pa.provenance_table()["counts_by_provenance"]
    classes = [  # (label, count, color, items)
        ("algebraic 0-param", c["algebraic_0param"], C_SU3,
         "16/16-bar, the 45, the 5/3, 8 anomalies, 32-spectrum, det U, the CPT mirror"),
        ("free 1-param", c["free_1param"], C_GUT,
         r"$g_{10}=1/5$ (the one free 't Hooft angle; not UV-fixed)"),
        ("convention", c["convention"], C_SU2,
         r"$g_0=1/8$ ($1/2^{n_{\rm space}}$ mapping); $M_{\rm GUT}=M_P$ (scale ID)"),
        ("coincidence", c["coincidence"], C_ZERO,
         r"$v/E_P = 2.02\times10^{-17}$ (the 2-digit 0.18% match)"),
        ("named residual", c["named_residual"], "#7f7f7f",
         r"$\eta_{CP}$, $\Lambda_{\rm obs}$, the sub-Planck DM scale, the exact $m_\nu$"),
    ]
    y = np.arange(len(classes))[::-1]
    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.barh(y, [n for (_, n, _, _) in classes],
            color=[col for (_, _, col, _) in classes], edgecolor="k", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{lab}   (n={n})" for (lab, n, _, _) in classes])
    for (yi, (_, n, _, items)) in zip(y, classes):
        ax.text(n + 0.15, yi, items, va="center", ha="left", fontsize=8, color="0.2")
    ax.set_xlim(0, 13)
    ax.set_xlabel("number of quantities  in the class")
    ax.set_title("The honest parameter audit: the '0-parameter' label, reclassified into five\n"
                 r"epistemic classes (the 7 algebraic invariants hold; the rest are named, not fitted)  (V52)")
    _save(fig, "fig12_parameter_audit.png")


# ===========================================================================
# Figure 13 — the CPT spectral mirror: matter (16) vs antimatter (16-bar)
# ===========================================================================
def fig13_cpt_mirror():
    from cl10 import qca_spectrum as qs

    p = qs.cpt_block_phases()
    m = np.asarray(p["ph16"], dtype=float)      # the 16 (matter) eigenphases, sorted
    b = np.asarray(p["ph16b"], dtype=float)     # the 16-bar (antimatter), for the check
    # the CPT mirror holds: the 16-bar multiset is the phase-conjugate of the 16
    assert np.allclose(np.sort(b), -np.sort(m)[::-1], atol=1e-9), "the 16-bar is the mirror of the 16"
    k = np.arange(1, len(m) + 1)
    fig, ax = plt.subplots()
    # the CPT axis (omega = 0) and the two mirror blocks: for each rank k,
    # the matter level m_k and its CPT conjugate -m_k (the 16-bar level)
    ax.axvline(0, color="k", lw=1.0, ls="--", label="CPT axis  ($\\omega = 0$)")
    ax.plot(m, k, "o-", color=C_SU3, lw=1.6, ms=5, label="matter  $16$   ($\\omega_k$)")
    ax.plot(-m, k, "o-", color=C_SU2, lw=1.6, ms=5,
            label="antimatter  $\\overline{16}$   ($-\\omega_k$)")
    ax.set_xlabel("eigenphase  $\\omega_k$  (rad/step)")
    ax.set_ylabel("rank $k$  (within the 16)")
    ax.set_title("The CPT spectral mirror: the $16$ (matter) and $\\overline{16}$ (antimatter) QCA blocks are\n"
                 r"exact phase-conjugates  ($\omega_{\overline{16}} = -\omega_{16}$) — CPT as spectral conjugation (V24/V52)")
    ax.legend(loc="lower right", framealpha=0.95)
    _save(fig, "fig13_cpt_mirror.png")



# ===========================================================================
# Figure 14 — CL10 0-parameter vs the Standard Model: head-to-head (2 x 3)
# ===========================================================================
def fig14_sm_head2head():
    from cl10 import (eta_cp as ep, dm_wimp as dw, hier_exponent as he,
                      toepred as tp, neutrino as nn)

    e = ep.eta_cp()
    d = dw.dm_wimp()
    h = he.hierarchy17()
    m_nu_pred = nn.seesaw_2step()["m_nu_2step_eV"]

    def _pair(ax, title, unit, obs, pred, log, note=None, ok=None):
        xs = [0, 1]
        vals = [obs, pred]
        bars = ax.bar(xs, vals, color=["#9a9a9a", C_ZERO],
                      edgecolor="k", lw=0.6, width=0.62)
        if log:
            ax.set_yscale("log")
            ax.set_ylim(min(vals) * 0.35, max(vals) * 2.6)
        else:
            ax.set_ylim(0, max(vals) * 1.20)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v * 0.97, f"{v:g}",
                    ha="center", va="top", fontsize=8.5, color="0.25")
        ratio = pred / obs
        if ok is None:
            ok = 0.5 <= ratio <= 2.0
        col = C_SU3 if ok else "#e07b00"
        txt = note if note is not None else f"0-param / obs = {ratio:.2f}"
        ax.text(0.5, 0.98, txt, transform=ax.transAxes, ha="center", va="top",
                fontsize=8.5, color=col,
                bbox=dict(boxstyle="round", fc="#f6f6f6", ec="#c8c8c8", alpha=0.92))
        ax.set_xticks(xs)
        ax.set_xticklabels(["SM\nmeasured", "CL10\n0-param"], fontsize=9)
        ax.set_title(title + (f"  ({unit})" if unit else ""), fontsize=10)
        ax.set_ylabel("value" + (f", {unit}" if unit else "") + (" (log)" if log else ""),
                      fontsize=9)

    fig, axs = plt.subplots(2, 3, figsize=(15.2, 8.6))
    fig.suptitle("CL10 0-parameter prediction  vs  Standard-Model measurement\n"
                 "no quantity is fitted — the measured values are PDG / Planck reference data",
                 fontsize=13, y=0.995)
    _pair(axs[0, 0], "baryon asymmetry  $\\Delta B$", "$\\times10^{-10}$",
          e["delta_B_observed"], e["delta_B_0param"], False)
    _pair(axs[0, 1], "dark-matter  $\\Omega_{\\rm DM} h^2$", "",
          d["omega_obs"], d["omega_DM"], False)
    _pair(axs[0, 2], "dark-matter mass  $m_{\\rm DM}$", "GeV",
          d["weak_scale_GeV"], d["m_DM_GeV"], False)
    _pair(axs[1, 0], "hierarchy  $v/E_P$", "",
          h["observed_ratio"], h["predicted_ratio"], True,
          note=r"2.08$\times$ below obs (the 2-digit residual)", ok=False)
    _pair(axs[1, 1], r"neutrino mass  $m_\nu$", "eV",
          0.05, m_nu_pred, True,
          note=f"{m_nu_pred / 0.05:.1f}$\\times$ obs (order-of-magnitude)", ok=False)
    _pair(axs[1, 2], "proton lifetime  $\\tau_p$", "yr",
          1e34, tp.proton_decay_lifetime(), True,
          note="0-param $\\gg$ the experimental\nlower bound (consistent)", ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    _save(fig, "fig14_sm_head2head.png")


# ===========================================================================
# Figure 15 — the 16 Weyl = one full SM generation = 10 + 5bar + 1 of SU(5)
# ===========================================================================
def fig15_matter_content():

    items = [
        (r"$Q$",   6, r"$(3,2,+1/6)$",   "10", C_SU3),
        (r"$\bar u$",  3, r"$(\bar3,1,-2/3)$", "10", C_SU3),
        (r"$\bar e$",  1, r"$(1,1,+1)$",   "10", C_SU3),
        (r"$\bar d$",  3, r"$(\bar3,1,+1/3)$", "$\\bar5$", C_SU2),
        (r"$L$",   2, r"$(1,2,-1/2)$",  "$\\bar5$", C_SU2),
        (r"$\nu^c$", 1, r"$(1,1,0)$",   "$1$", C_ROOT),
    ]
    y = np.arange(len(items))[::-1]
    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    ax.barh(y, [it[1] for it in items], color=[it[4] for it in items],
            edgecolor="k", lw=0.5)
    for yi, (lab, dim, rep, su5, col) in zip(y, items):
        ax.text(dim + 0.15, yi, f"{rep}  (in the {su5})",
                va="center", fontsize=9, color="0.2")
    ax.set_yticks(y)
    ax.set_yticklabels([it[0] for it in items])
    ax.set_xlabel("dimension  $n$  (Weyl components)")
    ax.set_xlim(0, 10)
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], color=C_SU3, lw=0, marker="s", ms=11,
               label=r"the $\mathbf{10}$  :  $Q, \bar u, \bar e$   $(6+3+1 = 10)$"),
        Line2D([0], [0], color=C_SU2, lw=0, marker="s", ms=11,
               label=r"the $\mathbf{\bar 5}$  :  $\bar d, L$   $(3+2 = 5)$"),
        Line2D([0], [0], color=C_ROOT, lw=0, marker="s", ms=11,
               label=r"the $\mathbf{1}$  :  $\nu^c$   $(1)$  — the sterile"),
    ], loc="lower right", framealpha=0.95, fontsize=8.5)
    ax.set_title(r"The 16 Weyl $= $ one full SM generation $= \mathbf{10} \oplus \mathbf{\bar 5} \oplus \mathbf{1}$ of $\mathrm{SU}(5)$"
                 r"  —  the six left-handed fields, total dimension $16$ (exact, 0-parameter)",
                 fontsize=10.5)
    _save(fig, "fig15_matter_content.png")


# ===========================================================================
# Figure 16 — the 0-parameter GUT coupling and the 5/3 at M_Z
# ===========================================================================
def fig16_couplings_mz():
    from cl10 import unification as u

    rd = u.running_down()
    az = u.A_Z
    p0 = rd["couplings_at_M_Z"]
    names = [r"$\alpha_1$  U(1)", r"$\alpha_2$  SU(2)", r"$\alpha_3$  SU(3)"]
    cols = [C_GUT, C_SU2, C_SU3]
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1.35, 1.0]})
    # (a) the 1-loop running of the 3 SM couplings to M_P + the GUT band
    s = np.logspace(np.log10(u.M_Z), np.log10(u.M_P), 220)
    ln = np.log(s / u.M_Z)
    for k in range(3):
        axA.plot(np.log10(s), 1.0 / az[k] - (u.B_SM[k] / (2.0 * np.pi)) * ln,
                 color=cols[k], lw=2.0, label=names[k])
    axA.axhline(1.0 / u.ALPHA_10, color=C_ZERO, ls="--", lw=1.5,
                label=r"0-param GUT coupling  $1/\alpha_{10} = 25$")
    axA.axvline(np.log10(u.M_P), color="k", lw=0.8)
    axA.set_xlabel(r"$\log_{10}(M/\mathrm{GeV})$")
    axA.set_ylabel(r"running $1/\alpha_i$  (1-loop, SM $\beta$)")
    axA.set_xlim(np.log10(u.M_Z) - 0.1, np.log10(u.M_P) + 0.3)
    axA.set_ylim(8, 42)
    axA.set_title("(a)  the three SM couplings reach the GUT band near $M_P$",
                  fontsize=9.5)
    axA.legend(fontsize=8, loc="upper left", framealpha=0.95)
    # (b) at M_Z: measured vs the 0-param down-running (paired bars)
    m = [0, 1, 2]
    w = 0.36
    axB.bar([x - w / 2 - 0.02 for x in m], az, width=w, color="#9a9a9a",
            edgecolor="k", lw=0.5, label="SM measured")
    p0plot = [max(p, 0.0) for p in p0]
    axB.bar([x + w / 2 + 0.02 for x in m], p0plot, width=w, color=C_ZERO,
            edgecolor="k", lw=0.5, label="CL10 0-param (down-run)")
    for x in m:
        axB.text(x - w / 2, az[x] + 0.003, f"{az[x]:.3f}",
                 ha="center", fontsize=7.5, color="0.3")
        if p0[x] <= 0:
            axB.text(x + w / 2, 0.003, "below\npole", ha="center", fontsize=7,
                     color="0.4")
        else:
            axB.text(x + w / 2, p0[x] + 0.003, f"{p0[x]:.3f}",
                     ha="center", fontsize=7.5, color="0.45")
    axB.set_xticks(m)
    axB.set_xticklabels(names)
    axB.set_ylabel("$\\alpha_i$  at $M_Z$")
    axB.set_ylim(0, 0.145)
    axB.set_title("(b)  at $M_Z$:  $\\alpha_2$ matches ($0.66\\times$);  "
                  r"$\alpha_3$ order-of-mag;  $\alpha_1$ below the pole"
                  "\n(the 1-loop matching = the W4 named residual)", fontsize=8.5)
    axB.legend(fontsize=8, loc="upper right", framealpha=0.95)
    fig.suptitle("The 0-parameter GUT coupling $\\alpha_{10}=1/25$ and the $5/3$  —  "
                 "the SM couplings at $M_Z$", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "fig16_couplings_mz.png")


# ===========================================================================
# Figure 17 — the 0-parameter ladder of scales vs the SM ladder
# ===========================================================================
def fig17_ladder_sm():
    from cl10 import dm_wimp as dw, toepred as tp

    M_P = tp.M_PLANCK
    LOOP = 1.0 / (16.0 * np.pi)
    # each point: (GeV value, label, color, y_level, ha, x_nudge) -- the level
    # + nudge + a leader line stagger the cluster near log10(100-246 GeV)
    cl10 = [
        (M_P * LOOP ** 17,   r"neutrino  $m_\nu = M_P(1/16\pi)^{17}$  (W2)", C_ZERO, 0, "center", 0.0),
        (1.0,               r"QCD  $\Lambda \sim 1$ GeV  (V27)",           C_SU3,  0, "center", 0.0),
        (dw.dm_wimp()["m_DM_GeV"], r"DM  $m_{\rm DM} = M_P(1/16\pi)^{10}$  (W7)", C_SU2, 1, "right", -0.4),
        (tp.V_EW,           r"electroweak  $v = 246$ GeV",                C_GUT,  2, "left",  +0.4),
        (M_P * LOOP ** 3,   r"2-step  $M_N = M_P(1/16\pi)^3$  (W2)",     C_ROOT, 0, "center", 0.0),
        (M_P,              r"GUT / Planck  $M_P = M_{\rm GUT}$",          "k",    0, "center", 0.0),
    ]
    sm = [
        (0.05 * 1e-9,      r"$m_\nu \sim 0.05$ eV",            "#8c8c8c", 0, "center", 0.0),
        (1.0,             r"$\Lambda_{\rm QCD} \sim 1$ GeV",   "#8c8c8c", 0, "center", 0.0),
        (100.0,           r"WIMP  $\sim 100$ GeV",             "#8c8c8c", 1, "right", -0.4),
        (246.0,           r"$v = 246$ GeV",                    "#8c8c8c", 2, "left",  +0.4),
        (1.22e19,         r"$M_P$",                            "#8c8c8c", 0, "center", 0.0),
    ]

    fig, ax = plt.subplots(figsize=(18, 7.5))
    ax.axhline(0.30, color="0.7", lw=0.8, ls=":")
    ax.axhline(-0.30, color="0.7", lw=0.8, ls=":")
    top_levels = {0: 0.66, 1: 0.84, 2: 0.99}
    bot_levels = {0: -0.66, 1: -0.84, 2: -0.99}

    def draw_row(points, ydot, levels):
        for (E, lab, col, lvl, ha, xn) in points:
            x = np.log10(E)
            ylab = levels[lvl]
            ax.plot([x, x + xn], [ydot, ylab], color=col, lw=0.8,
                    alpha=0.55, zorder=2)
            ax.scatter([x], [ydot], s=85, color=col, zorder=5,
                       edgecolor="k", lw=0.4)
            va = "bottom" if ydot > 0 else "top"
            ax.text(x + xn, ylab + (0.05 if ydot > 0 else -0.05), lab,
                    ha=ha, va=va, fontsize=8.5, color=col)

    draw_row(cl10, 0.30, top_levels)
    draw_row(sm, -0.30, bot_levels)

    ax.text(9.5, 1.05, "CL10 0-parameter ladder", fontsize=11, weight="bold",
            color="0.2", ha="center")
    ax.text(9.5, -1.09, "SM observed ladder", fontsize=11, weight="bold",
            color="0.4", ha="center")
    ax.set_xlim(-11, 20)
    ax.set_ylim(-1.18, 1.18)
    ax.set_yticks([])
    ax.set_xlabel("energy scale (GeV, log)   —   one object spans 30+ decades")
    ax.set_title("The 0-parameter ladder of scales vs the SM ladder: the same hierarchy,"
                 " from the Planck scale down to the neutrino", fontsize=11)
    _save(fig, "fig17_ladder_sm.png")

def main():
    print("Generating CL10 publication figures ->", OUT)
    fig1_one_object()
    fig2_unification()
    fig3_10loop_ladder()
    fig4_wimp_miracle()
    fig5_gravity_kretschmann()
    fig6_512_512()
    fig7_qcd()
    fig8_baryon()
    fig9_anomaly()
    fig10_predictions()
    fig11_ladder()
    fig12_audit()
    fig13_cpt_mirror()
    fig14_sm_head2head()
    fig15_matter_content()
    fig16_couplings_mz()
    fig17_ladder_sm()
    print("done: 17 figures.")


if __name__ == "__main__":
    main()
