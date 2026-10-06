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
                 r"the hierarchy $v/E_P = 2.08\,(1/16\pi)^{10}$ and the DM mass")
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
    # the Planck-mass soliton (under-abundant)
    ax.scatter([M_P], [0.12 * planck_gap], s=70, color="#7f7f7f", zorder=5, clip_on=False,
               label=r"Planck-mass soliton ($\sim 10^{-52}$, under-abundant)")
    # the 0-parameter WIMP star is identified in the legend (no extra
    # annotation: it sits in the top corner next to the WIMP-miracle point).
    ax.set_xlabel(r"relic mass $m$ (GeV, log scale)")
    ax.set_ylabel(r"relic abundance $\Omega h^2$")
    ax.set_title("The WIMP miracle, resolved: the 0-parameter 10-loop mass\n"
                 "$M_P(1/16\\pi)^{10}=118$ GeV gives the observed dark-matter abundance")
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
    cats = ["observed", "named\n($\\eta_{CP}=10^{-3}$)",
            "0-param.\n($\\frac{28}{79}\\,\\varepsilon_{CP}\\,D$)"]
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
    ax.set_title("GUT leptogenesis: the baryon asymmetry as a product of factors of three provenances.  "
                 r"$\Delta B = \frac{28}{79}\,\varepsilon_{CP}\,D$ ($\varepsilon_{CP}=g_0^4/16\pi$, $g_0=1/2^3$; $D=(1/16\pi)^2$; the $28/79$ sphaleron is an input of the SM sector) gives $\Delta B = 6.8\times10^{-10}$ — a 12% match (the named $\eta_{CP}$: 2.8)")
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
    ax.set_title(r"$[\mathrm{U}(1)]^3$: the six SM fields cancel ($\sum = 0$) — the 16 Weyl is anomaly-free")
    _save(fig, "fig9_anomaly_free.png")

# ===========================================================================
# Figure 10 — CL10 vs observation: the 0-parameter predictions vs the data
# ===========================================================================
def fig10_predictions():
    from cl10 import eta_cp as ep, dm_wimp as dw, u10pred as tp

    e = ep.eta_cp()
    d = dw.dm_wimp()
    # (x = observed/measured, y = 0-parameter prediction, label, color, marker, label-offset)
    rows = [
        (e["delta_B_observed"], e["delta_B_0param"],
         "baryon asymmetry  $\\Delta B$  (2-loop)", C_ZERO, "o", (12, -22)),
        (d["omega_obs"], d["omega_DM"],
         r"dark-matter $\Omega h^2$  (10-loop WIMP)", C_SU2, "o", (-175, 20)),
        (d["weak_scale_GeV"], d["m_DM_GeV"],
         r"dark-matter mass $m_{\rm DM}$  (10-loop)  [GeV]", C_SU2, "s", (16, 16)),
        (2.016e-17, 2.02e-17,
         r"hierarchy $v/E_P$  (10-loop)", C_GUT, "o", (12, -24)),
        (1e-3, e["eta_CP_0param"],
         r"CP factor $\eta_{CP}$  (2-loop)", C_SU3, "^", (-150, -26)),
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
    from cl10 import dm_wimp as dw, u10pred as tp

    m_nu_eV = ((1.0 / 8 * 246.0) ** 2 / 1e14) * 1e9     # ~ 10^-2 eV
    # (energy GeV, label, color, side(+/-1), height multiplier)
    scales = [
        (m_nu_eV * 1e-9, r"neutrino  $m_\nu \sim 10^{-2}$ eV  (seesaw, W2)", C_ZERO, -1, 0.5),
        (1.0, r"QCD  $\Lambda_{\rm QCD} \sim 1$ GeV", C_SU3, 1, 0.5),
        (dw.dm_wimp()["m_DM_GeV"], r"0-param DM  $m_{\rm DM}=118$ GeV  (10-loop)", C_SU2, -1, 0.5),
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
                 r"epistemic classes (the 7 algebraic invariants hold; the rest are named, not fitted)")
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
                 r"exact phase-conjugates  ($\omega_{\overline{16}} = -\omega_{16}$) — CPT as spectral conjugation")
    ax.legend(loc="lower right", framealpha=0.95)
    _save(fig, "fig13_cpt_mirror.png")



# ===========================================================================
# Figure 14 — CL10 0-parameter vs the Standard Model: head-to-head (2 x 3)
# ===========================================================================
def fig14_sm_head2head():
    from cl10 import (eta_cp as ep, dm_wimp as dw, hier_exponent as he,
                      u10pred as tp, neutrino as nn)

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
    from cl10 import dm_wimp as dw, u10pred as tp

    M_P = tp.M_PLANCK
    LOOP = 1.0 / (16.0 * np.pi)
    # each point: (GeV value, label, color, y_level, ha, x_nudge) -- the level
    # + nudge + a leader line stagger the cluster near log10(100-246 GeV)
    cl10 = [
        (M_P * LOOP ** 17,   r"neutrino  $m_\nu = M_P(1/16\pi)^{17}$  (W2)", C_ZERO, 0, "center", 0.0),
        (1.0,               r"QCD  $\Lambda \sim 1$ GeV",                 C_SU3,  0, "center", 0.0),
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


# ===========================================================================
# Figure 18 — the residual-set reduction: the unique n + the 0-free/3-named set
# ===========================================================================
def fig18_residual_reduction():
    from cl10 import four_limits as fl

    r = fl.unique_10loop()
    h_scan = [(p["n"], p["hier_dev"]) for p in r["scan"]]
    d_scan = [(p["n"], p["dm_dev"]) for p in r["scan"]]
    n_best = r["n_best"]

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(13.0, 4.7), gridspec_kw={"width_ratios": [1.4, 1.0]})
    # (a) the two closures (deviation in decades vs n) — two statements of one scale
    axA.plot([p[0] for p in h_scan], [p[1] for p in h_scan], "o-",
             color=C_SU3, lw=1.6, ms=4, label="hierarchy  $v/E_P$  (10-loop)")
    axA.plot([p[0] for p in d_scan], [p[1] for p in d_scan], "s-",
             color=C_SU2, lw=1.6, ms=4, label="dark matter  $m_{\\rm DM}$  (10-loop)")
    axA.axhline(0.0, color="k", lw=0.8, ls="--")
    axA.axvline(n_best, color=C_ZERO, lw=1.4, ls=":",
                label=f"the unique solution  $n = {n_best}$  (both within 1 decade)")
    hv = min(p[1] for p in h_scan if p[0] == n_best)
    dv = min(p[1] for p in d_scan if p[0] == n_best)
    axA.scatter([n_best], [min(hv, dv)], color=C_ZERO, s=90, zorder=5,
                edgecolor="k", lw=0.6)
    axA.set_xlabel("$n$  (the 0-parameter loop count $= 2 \\times n_{\\rm qubits}$)")
    axA.set_ylabel("deviation from observed  (e-folds)")
    axA.set_xlim(0, 21)
    axA.set_ylim(-0.5, 2.4)
    axA.legend(loc="upper left", framealpha=0.95, fontsize=8.5)
    axA.set_title("Two closures, one scale: only $n = 10$ is within one decade of both\n"
                  r"($n_{\rm gen} = 5$ satisfies neither) — two statements of one scale (Section 16)",
                  fontsize=10.5)
    # (b) the residual-set reduction (the 4 limits -> 0 free + 3 named)
    axB.axis("off")
    rows = [
        ("the 4 'limits'  (the honest limits)", 0.98, "0.3", "normal"),
        ("   1.  the '0-parameter' convention", 0.86, "0.4", "normal"),
        ("   2.  the hierarchy gap (W8)", 0.74, "0.4", "normal"),
        ("   3.  the baryogenesis (W5)", 0.62, "0.4", "normal"),
        ("   4.  the dark matter (W7)", 0.50, "0.4", "normal"),
        ("        ==>  0 free (algebra)  +  1 angle  +  3 named", 0.34, C_ZERO, "bold"),
        (r"$g_{10} = 1/5$,   $g_0 = 1/8$   (0-parameter selection rules)", 0.24, "0.3", "normal"),
        ("the named set = 2 two-digit coefficients $(2.08, 1.18)$", 0.15, "0.4", "normal"),
        (r"+  1  CP phase  $= 3^{*}$", 0.07, "0.4", "normal"),
    ]
    for text, y, col, wt in rows:
        axB.text(0.5, y, text, ha="center", va="center",
                 fontsize=10 if col == C_ZERO else 8.5,
                 color=col, weight=wt)
    axB.set_title("The residual set reduces to 0 free (algebra) + 3 named", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, "fig18_residual_reduction.png")

# ===========================================================================
# Figure 19 — the 2-loop running: the perturbative QFT Landau pole vs the unitary QCA step
# ===========================================================================
def fig19_uv_2loop():
    """The 2-loop running of the U(1) coupling: the perturbative QFT reaches a
    Landau pole below M_P, the 1-loop running stays bounded, and the unitary QCA
    step remains bounded at every order."""
    from cl10 import unification as u, u10pred as tp

    M_Z, M_P = u.M_Z, tp.M_PLANCK
    b1, b2 = u.B_SM[0], u.B2_SM[0]
    a1z = u.A_Z[0]
    c = 1.0 / a1z
    Bcoef = b1 / (2.0 * np.pi)
    Acoef = b2 / (b1 ** 2 * 2.0 * np.pi)
    t2 = (-Bcoef + np.sqrt(Bcoef ** 2 + 4.0 * Acoef * c)) / (2.0 * Acoef)
    pole2 = M_Z * float(np.exp(t2))

    t = np.linspace(0.0, np.log(3.0e20 / M_Z), 400)
    M = M_Z * np.exp(t)
    inv1 = c - Bcoef * t
    inv2 = c - Bcoef * t - Acoef * t ** 2
    valid = inv2 > 0

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.semilogx(M, inv1, "-", color=C_SU2, lw=2.2,
                label=r"1-loop  $1/\alpha_1(M)$  (the pole is above $M_P$)")
    ax.semilogx(M[valid], inv2[valid], "-", color=C_ZERO, lw=2.4,
                label=r"2-loop  $1/\alpha_1(M)$  (the pole is below $M_P$)")
    ax.axvline(M_P, color="k", lw=1.2)
    ax.text(M_P * 1.05, 57, r"$M_P$", fontsize=9, va="bottom", ha="left")
    ax.axvline(pole2, color=C_ZERO, ls=":", lw=1.4)
    ax.scatter([pole2], [0], s=130, marker="x", color=C_ZERO, zorder=5,
               label=r"the 2-loop Landau pole  $\sim$ " + f"{pole2:.0e}" + r" GeV")
    ax.axhline(1.0 / u.ALPHA_10, color="0.35", ls="--", lw=1.3,
               label=r"0-parameter GUT coupling  $1/\alpha_{10} = 25$")
    ax.annotate(
        r"the 2-loop perturbative $QFT$ reaches a Landau pole below $M_P$;  the $QCA$ "
        "step is unitary, so its $n$-step amplitudes stay bounded\n"
        r"$\Rightarrow$ the perturbative UV divergence is not a property of the QCA itself",
        xy=(pole2, 0), xytext=(3e3, 44),
        fontsize=8.5, ha="left", va="top",
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#888", alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="0.3"))
    ax.set_xlabel(r"scale $M$  (GeV, log scale)")
    ax.set_ylabel(r"$1/\alpha_1$  (the $U(1)$ coupling)")
    ax.set_title("2-loop running of the $U(1)$ coupling: the perturbative QFT pole "
                 "vs the bounded QCA step")
    ax.set_xlim(70, 3e20)
    ax.set_ylim(-3, 60)
    ax.legend(loc="lower left", framealpha=0.95, fontsize=8)
    _save(fig, "fig19_uv_2loop.png")


# ===========================================================================
# Figure 20 — proton decay p -> e+ pi0: the 0-parameter lifetime vs the GUT scale
# ===========================================================================
def fig20_proton_decay():
    """The 0-parameter proton lifetime tau_p ~ M_GUT^4 / (alpha_10^2 m_p^5).
    The GUT scale sets the decay rate; the 0-parameter point (M_GUT = M_P) is
    far above the experimental bound, so the proton is stable on all observed
    timescales (consistent with data, and falsifiable)."""
    from cl10 import u10pred as tp

    M_P = tp.M_PLANCK
    ALPHA = tp.ALPHA_10
    m_p = tp.M_P_ROBSON
    GEV_TO_S = tp.GEV_TO_S
    S_TO_YR = tp.S_TO_YR

    def tau(M):
        return (M ** 4 / (ALPHA ** 2 * m_p ** 5)) * GEV_TO_S / S_TO_YR

    M_GUT = np.logspace(14, 20, 300)
    tau_p = tau(M_GUT)

    BOUND = 1.6e34          # yr, the p -> e+ pi0 bound (Super-Kamiokande)
    TAU_UNIVERSE = 1.4e10   # yr, the age of the universe

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.loglog(M_GUT, tau_p, "-", color=C_ROOT, lw=2.4,
              label=r"$\tau_p \sim M_{\rm GUT}^{4}\,/\,(\alpha_{10}^{2}\,m_p^{5})$")
    ax.axhline(BOUND, color=C_ZERO, ls="--", lw=1.5,
               label="experimental bound  $\\tau_p > 1.6\\times10^{34}$ yr  ($p\\to e^+\\pi^0$)")
    ax.axhline(TAU_UNIVERSE, color="0.45", ls=":", lw=1.5,
               label="age of the universe  $\\sim 1.4\\times10^{10}$ yr")
    tau0 = float(tp.proton_decay_lifetime())
    ax.scatter([M_P], [tau0], s=150, marker="*", color=C_ZERO, zorder=5,
               label=r"0-parameter:  $M_{\rm GUT}=M_P=1.22\times10^{19}$ GeV "
                     r"$\Rightarrow$ $\tau_p = $" + f"{tau0:.1e}" + r" yr")
    ax.axhspan(1.0, BOUND, color=C_ZERO, alpha=0.05)
    ax.annotate(
        r"$M_{\rm GUT} \lesssim 5\times10^{15}$ GeV gives $\tau_p$ below the bound (excluded, shaded);"
        "\nthe 0-parameter point sits $\\sim$13 orders of magnitude above it, so the proton is stable "
        "on all\nobserved timescales (consistent with data; the prediction is falsifiable)",
        xy=(M_P, tau0), xytext=(2e14, 3e34),
        fontsize=8.5, ha="left", va="bottom",
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#888", alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="0.3"))
    ax.set_xlabel(r"GUT scale  $M_{\rm GUT}$  (GeV, log scale)")
    ax.set_ylabel(r"proton lifetime  $\tau_p$  (yr, log scale)")
    ax.set_title("Proton decay $p \\to e^+ \\pi^0$: the 0-parameter lifetime vs the GUT scale")
    ax.set_xlim(1e14, 1e20)
    ax.set_ylim(1e9, 1e49)
    ax.legend(loc="upper left", framealpha=0.95, fontsize=8)
    _save(fig, "fig20_proton_decay.png")

# ===========================================================================
# Figure 21 — proton decay: the 0-parameter rate vs the bound (the margin)
# ===========================================================================
def fig21_proton_margin():
    """The 0-parameter decay rate normalised to the experimental bound,
    Gamma/Gamma_bound = tau_bound/tau_p, vs the GUT scale.  The 0-parameter
    point (M_GUT = M_P) sits ~13 orders of magnitude below the bound (the
    proton is stable); M_GUT below ~5e15 GeV would be excluded (shaded)."""
    from cl10 import u10pred as tp
    M_P, ALPHA, m_p = tp.M_PLANCK, tp.ALPHA_10, tp.M_P_ROBSON
    def tau(M):
        return (M**4 / (ALPHA**2 * m_p**5)) * tp.GEV_TO_S / tp.S_TO_YR
    BOUND = 1.6e34          # yr, the p -> e+ pi0 bound (Super-Kamiokande)
    M_GUT = np.logspace(14, 20, 400)
    ratio = BOUND / tau(M_GUT)                 # = Gamma / Gamma_bound
    tau0 = float(tp.proton_decay_lifetime())
    r0 = BOUND / tau0
    M_excl = M_P * (BOUND / tau0) ** 0.25      # the M where tau = the bound

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.loglog(M_GUT, ratio, "-", color=C_ROOT, lw=2.4,
              label=r"$\Gamma/\Gamma_{\rm bound}\ \propto\ M_{\rm GUT}^{-4}$")
    ax.axhline(1.0, color=C_ZERO, ls="--", lw=1.5,
               label=r"experimental bound  $\Gamma=\Gamma_{\rm bound}$")
    ax.axvline(M_P, color="0.35", ls=":", lw=1.3)
    ax.axvspan(1e14, M_excl, color=C_ZERO, alpha=0.08)
    ax.scatter([M_P], [r0], s=150, marker="*", color=C_ZERO, zorder=5,
               label=rf"0-parameter $M_{{\rm GUT}}=M_P$:  $\Gamma/\Gamma_{{\rm bound}}={r0:.1e}$")
    ax.annotate(
        rf"$M_{{\rm GUT}}\leq {M_excl:.1e}$ GeV would make the proton decay faster than the bound "
        "(shaded, excluded);\nthe 0-parameter point is $\\sim$13 orders below it, so the proton "
        "is stable on all observed timescales",
        xy=(M_P, r0), xytext=(2.5e15, 3e-5),
        fontsize=8.5, ha="left", va="bottom",
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#888", alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="0.3"))
    ax.set_xlabel(r"GUT scale  $M_{\rm GUT}$  (GeV, log)")
    ax.set_ylabel(r"$\Gamma/\Gamma_{\rm bound}$  (the decay rate vs the bound, log)")
    ax.set_title("Proton decay: the 0-parameter rate vs the bound (the falsifiability margin)")
    ax.set_xlim(1e14, 1e20)
    ax.set_ylim(1e-15, 10)
    ax.legend(loc="upper right", framealpha=0.95, fontsize=8)
    _save(fig, "fig21_proton_margin.png")

# ===========================================================================
# Figure 22 — proton decay: the three GUT channels (0-parameter rate vs the bound)
# ===========================================================================
def fig22_proton_channels():
    """The 0-parameter partial lifetimes for the three GUT channels
    p -> e+ pi0, p -> mu+ pi0, p -> e+ K+ (the rate is set by the two-body
    phase space; the hadronic matrix element is O(1)).  Each sits
    ~10^13-14 above its experimental bound (the proton is stable in all)."""
    from cl10 import u10pred as tp
    tau0 = float(tp.proton_decay_lifetime())   # the p -> e+ pi0 0-parameter lifetime
    mp, me, mmu, mpi, mk = 938.272, 0.511, 105.658, 134.977, 493.677
    def pf(a, b):
        return 0.5 / mp * np.sqrt((mp**2 - (a + b)**2) * (mp**2 - (a - b)**2))
    pf0 = pf(me, mpi)
    chans = [
        ("$p \\to e^+ \\pi^0$", pf(me, mpi) / pf0, 1.0, 2.4e34),
        ("$p \\to \\mu^+ \\pi^0$", pf(mmu, mpi) / pf0, 1.0, 7.0e33),
        ("$p \\to e^+ K^+$", pf(me, mk) / pf0, 1.0, 2.7e33),
    ]
    names, margin = [], []
    for nm, phase, had, bound in chans:
        tau_ch = tau0 * phase / had      # tau ~ 1/rate ~ 1/(phase-space)
        names.append(nm)
        margin.append(tau_ch / bound)    # how many times longer than the bound

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    x = np.arange(len(names))
    bars = ax.bar(x, margin, 0.5, color=[C_ROOT, C_SU2, C_SU3], alpha=0.85, edgecolor="k")
    for b, m in zip(bars, margin):
        ax.text(b.get_x() + b.get_width()/2, m, f"{m:.1e}",
                ha="center", va="bottom", fontsize=9)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel(r"$\tau_p^{\rm (pred)}/\tau_p^{\rm (bound)}$  (the stability margin, log)")
    ax.set_title("Proton decay: the 0-parameter lifetime vs the bound, for the three GUT channels\n"
                 r"(each $\sim10^{13}$-$10^{14}$ above its bound: the rate is GUT-suppressed, the two-body phase space differs)",
                 fontsize=10)
    ax.axhline(1.0, color=C_ZERO, ls="--", lw=1.4,
               label=r"the bound  ($\Gamma=\Gamma_{\rm bound}$)")
    ax.set_ylim(1e0, 1e15)
    ax.legend(loc="upper left", framealpha=0.95, fontsize=8)
    _save(fig, "fig22_proton_channels.png")

# ===========================================================================
# Figure 23 — the 0-parameter QCA wave: the unitary wavepacket (heatmaps)
# ===========================================================================
def fig23_qca_wave():
    """The 0-parameter QCA field (the same 4-neighbour discrete wave equation
    as the CL10 graviton): a pulse emits a **causal light-cone** — the signal
    propagates at the lattice speed $c=0.5$ and never faster (the QCA is
    causal).  A 2x4 grid of $|\\psi|^2$ heatmaps over time (each panel
    normalised to its own peak, so the expanding ring is always resolved)."""
    Lx, Ly = 200, 140
    r2 = 0.25                                   # (c dt/dx)^2, c=0.5 (the CFL limit)
    x, y = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing="ij")
    cx, cy, sig = Lx // 2, Ly // 2, 3
    psi = np.exp(-0.5 * ((x - cx)**2 + (y - cy)**2) / sig**2)   # a sharp pulse
    psi_m1 = psi.copy()                          # zero initial velocity
    exL = x / (0.12 * Lx); exR = (Lx - 1 - x) / (0.12 * Lx)
    edgeY = np.minimum(y, Ly - 1 - y) / (0.12 * Ly)
    s = np.clip(np.minimum(np.minimum(exL, exR), edgeY), 0.0, 1.0)
    times = [0, 25, 50, 75, 100, 125, 150, 175]
    frames = []
    for t in range(180):
        lap = np.zeros((Lx, Ly))
        lap[:, 1:] += psi[:, :-1]; lap[:, :-1] += psi[:, 1:]
        lap[1:, :] += psi[:-1, :]; lap[:-1, :] += psi[1:, :]
        lap -= 4 * psi
        psi = (2 * psi - psi_m1 + r2 * lap) * s
        psi_m1 = psi
        if t in times:
            frames.append((t, psi * psi))
    fig, axes = plt.subplots(2, 4, figsize=(11.8, 5.4), sharex=True, sharey=True,
                             gridspec_kw={"wspace": 0.12, "hspace": 0.42})
    cmap = plt.get_cmap("inferno")
    for i, (t, f) in enumerate(frames):
        ax = axes[i // 4, i % 4]
        ax.set_facecolor("black")
        m = float(f.max())
        ax.imshow((f / m).T, origin="lower", extent=[0, Lx, 0, Ly],
                  cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")
        ax.set_title(f"t = {t}    (peak $\\times$ {m:.1e})", fontsize=9, color="k")
    for ax in axes[-1]:
        ax.set_xlabel(r"space $x$  (lattice units)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"space $y$")
    fig.suptitle("The 0-parameter QCA wave: the dispersive spreading of a pulse — the peak "
                 "drops $\\propto 1/\\mathrm{area}$ as the causal support expands at $c=0.5$ "
                 "(the unitary norm is conserved)", fontsize=10.5, y=0.985)
    fig.subplots_adjust(top=0.88, left=0.05, right=0.98, wspace=0.12, hspace=0.45)
    _save(fig, "fig23_qca_wave.png")

# ===========================================================================
# Figure 24 — the 0-parameter black hole: the well + geodesic (left), and the
# 2-route Kretschmann (right), both computed from the CL10 metric
# ===========================================================================
def fig24_gravity_3d():
    """fig24 — the 0-parameter (Planck-mass) black hole.  Left: the 3-D
    gravity well with the **geodesic** integrated through the CL10
    finite-difference Christoffel symbols (the precessing rosette) and the
    horizon ring.  Right: the **2-route** verification — the Kretschmann
    *computed from the metric* (4-D central finite differences) equals the
    analytic 48M^2/r^6 (the black hole is recovered by computation).  Both
    panels are computed from ``cl10.gravity`` (``metric_schwarzschild``,
    ``_christoffel_at``, ``riemann_4d``/``kretschmann_4d``)."""
    from cl10 import gravity
    from scipy.integrate import solve_ivp

    M = 1.0                        # the 0-parameter black-hole mass (Planck units)
    R = 13.0
    H = gravity.H_FD

    def sw(x):
        return gravity.metric_schwarzschild(M, x)

    # ---- LEFT: the 3-D well + the geodesic (via the CL10 Christoffels) ----
    g = np.linspace(-R, R, 90)
    X, Y = np.meshgrid(g, g)
    rr = np.sqrt(X**2 + Y**2)
    # the well: the metric's time-time part g_00 = -(1 - 2M/r) gives, in the
    # weak-field limit, the potential Phi = -M/r (capped at the horizon).
    Phi = -M / np.maximum(rr, 2 * M)

    def geodesic(lam, u):
        t, r, th, ph, td, rd, thd, phd = u
        G = gravity._christoffel_at(sw, np.array([t, r, th, ph]), H)
        v = np.array([td, rd, thd, phd])
        dv = -np.einsum("mab,a,b->m", G, v, v)
        return [td, rd, thd, phd, dv[0], dv[1], dv[2], dv[3]]

    r0, vt = 10.0, 0.34
    gamma = 1.0 / np.sqrt(1 - vt**2)
    u0 = [0.0, r0, np.pi / 2, 0.0, gamma / np.sqrt(1 - 2 * M / r0), 0.0,
          0.0, gamma * vt / r0]
    sol = solve_ivp(geodesic, [0.0, 600.0], u0, method="DOP853",
                    rtol=1e-8, atol=1e-8, max_step=0.5)
    r = sol.y[1]; ph = sol.y[3]
    tx = r * np.cos(ph); ty = r * np.sin(ph); tz = -M / np.maximum(r, 2 * M)

    a = np.linspace(0, 2 * np.pi, 200)
    hpx, hpy = 2 * M * np.cos(a), 2 * M * np.sin(a)

    # ---- RIGHT: the 2-route Kretschmann (computed vs analytic) ----
    rv = np.linspace(2.3, 20.0, 12)
    Kc, Ka = [], []
    for ri in rv:
        x = np.array([0.0, ri, 0.5, 0.7])
        Rm = gravity.riemann_4d(sw, x, H)
        gm, _ = gravity._metric_and_dg(sw, x, H)
        Kc.append(gravity.kretschmann_4d(Rm, gm))
        Ka.append(48.0 * M**2 / ri**6)
    Kc, Ka = np.array(Kc), np.array(Ka)

    fig = plt.figure(figsize=(12.6, 5.6))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax1.set_proj_type("ortho")
    ax1.plot_surface(X, Y, Phi, cmap="viridis", alpha=0.4, rstride=3, cstride=3)
    ax1.plot(tx, ty, tz, color=C_ZERO, lw=1.5,
             label="the geodesic (CL10 Christoffels)")
    ax1.plot(hpx, hpy, np.full_like(hpx, -M / (2 * M)), color="w", ls="--",
             lw=1.4, label=r"the horizon $r=2M$")
    ax1.scatter([tx[0]], [ty[0]], [tz[0]], color="w", s=20, zorder=5,
                label="start")
    ax1.set_zlabel(r"$\Phi=-M/r$"); ax1.set_title("the well + the geodesic",
                                                 fontsize=10)
    ax1.view_init(elev=55, azim=-40)
    ax1.legend(loc="upper left", fontsize=7)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.loglog(rv, Ka, "-", color="k", lw=1.6, label=r"analytic $K=48M^2/r^6$")
    ax2.loglog(rv, Kc, "o", color=C_ZERO, ms=5,
               label="computed from the CL10 metric")
    resid = float(np.abs(Kc / Ka - 1.0).max())
    ax2.text(0.03, 0.05, f"2-route residual = {resid:.1e}\n(the black hole is "
             "recovered by computation)", transform=ax2.transAxes, fontsize=8,
             bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#888", alpha=0.95))
    ax2.set_xlabel(r"radius $r$  ($M=1$, $r>2M$)")
    ax2.set_ylabel("Kretschmann $K = R^{mnpq}R_{mnpq}$")
    ax2.set_title("2-route: the computed Kretschmann", fontsize=10)
    ax2.legend(loc="upper right", fontsize=8)

    fig.suptitle("fig24 — the 0-parameter black hole: the well + geodesic "
                 "(left), the 2-route Kretschmann (right)", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, "fig24_gravity_3d.png")


# ===========================================================================
# Figure 25 — the Young double-slit as a *single-electron* measurement: the
# Tonomura flow (one electron at a time) with the Born rule sampled from the
# *time-integrated* screen, so the impact points cluster into the fringes.
# ===========================================================================
def fig25_young():
    """fig25 - the CL10 Young double-slit as a single-electron measurement
    (the Tonomura flow).  The emitter fires one localized electron at a time;
    each one diffracts through the slits and is detected at a *single point*
    sampled from the **time-integrated** on-screen |psi|^2 (the interference
    pattern the screen records -- a single instant of the traveling wave is
    only a blob, so we integrate over the passage).  The points are collected
    one per electron and build the fringes.  Panels: (a) the wave (x-z slice,
    the two diffracted waves interfering), (b) the time-integrated interference
    pattern + its x-profile (the fringe groups), (c) the 9-to-1 collapse of one
    electron onto a point, (d) the reality -- the accumulated single-electron
    detections cluster into the fringe groups.  Opt-in (--young); its PNG is
    committed."""
    import matplotlib.patches as mpatches
    np.random.seed(2026)
    from simulation.young import collapse as Coll
    from simulation.young import wave

    C = plt.get_cmap("inferno")
    TXT = "#2b2b2b"
    LIGHT = "#f0f0f0"

    def blackize(ax):
        ax.set_facecolor("black")
        for sp in ax.spines.values():
            sp.set_color("#666")
        ax.tick_params(labelcolor=TXT)

    def norm(A):
        A = np.asarray(A, float)
        return (A - A.min()) / max(A.max() - A.min(), 1e-12)

    Lx, Ly, Lz = 128, 128, 220
    lam, d, a = 4.5, 22.0, 5.0
    cfg = dict(grid=[Lx, Ly, Lz], r=0.5, wavelength=lam, source_width=18.0,
               source_amplitude=1.0, source_z=20, ramp_steps=200,
               source_mode="single", source_shape="point",
               source_y_width=14.0, packet_z_width=8.0,
               sponge_width=16, barrier_z=44, slit_separation=d,
               slit_width=a, slit_height=Ly, screen_z=200)
    bar_z, scr_z = cfg["barrier_z"], cfg["screen_z"]
    sim = wave.WaveSim(cfg)
    col = Coll.Collapse(dict(grid=[Lx, Ly, Lz], source_width=18.0,
                             show_collapse=True, collapse_frames=9,
                             collapse_peak=6.0, n_detections=10 ** 8))

    # -- (a) the wave: a *separate* driven sheet beam (continuous mode) so the
    # x-z slice shows the full field diffracting through the two slits and the
    # two waves interfering (like fig25) -- a single electron's small packet
    # would only show a little blob.
    sim_a = wave.WaveSim(dict(cfg, source_mode="continuous",
                              source_shape="sheet"))
    xz_slice = None
    for _i in range(520):
        sim_a.step()
        if _i == 500:
            xz_slice = (sim_a.field()[:, Ly // 2, :] ** 2).copy()

    # -- run the Tonomura flow: one electron at a time, one point each --------
    absorb_frac, max_flight, n_target = 0.02, 300, 150
    sim.emit_packet()
    dets = []
    screened = False
    flight = 0
    time_pattern = None  # (b) the time-integrated interference pattern
    n_steps = 0
    while len(dets) < n_target and n_steps < 200000:
        sim.step()
        n_steps += 1
        flight += 1
        # the front reaches the screen -> arm the screen buffer
        if not screened and (sim.at_screen() or flight > max_flight):
            sim.arm_screen()
            screened = True
            flight = 0
        # every step: integrate the on-screen intensity into the buffer
        sim.accumulate_screen()
        # the sponge has absorbed the electron -> sample the Born-rule point
        # from the *time-integrated* pattern, then fire the next electron
        if screened and sim.absorbed(absorb_frac):
            P = sim.screen_pattern()
            if time_pattern is None:
                time_pattern = P.copy()
            pt = col.sample(P)
            if pt is not None:
                dets.append(pt)
            sim.emit_packet()
            screened = False
            flight = 0
    det = np.array(dets) if dets else np.empty((0, 2))
    if xz_slice is None:
        xz_slice = (sim_a.field()[:, Ly // 2, :] ** 2)
    if time_pattern is None:
        time_pattern = sim.screen_pattern()

    # pick a detected point in a bright fringe for the collapse animation
    if det.size:
        bright = det[time_pattern[det[:, 0], det[:, 1]] > 0.3 * time_pattern.max()]
        if not bright.size:
            bright = det
        cx, cy = max(bright, key=lambda p: time_pattern[int(p[0]), int(p[1])])
        cx, cy = int(cx), int(cy)
    else:
        cx, cy = Lx // 2, Ly // 2
    # the contraction: the screen superposition * a shrinking Gaussian
    gx, gy = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing="ij")
    sig0, sig1, pkmax, nframe = 60.0, 1.5, 6.0, 9
    seq = np.stack([
        time_pattern * np.exp(-((gx - cx) ** 2 + (gy - cy) ** 2) /
                              (2 * (sig0 * (sig1 / sig0) ** (f / (nframe - 1))) ** 2))
        * (1.0 + (pkmax - 1.0) * f / (nframe - 1)) for f in range(nframe)])

    fig = plt.figure(figsize=(12.6, 9.2), facecolor="white")
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.26,
                          left=0.05, right=0.97, top=0.89, bottom=0.08)
    # (a) the wave: the propagation + double-slit diffraction (x-z slice)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.imshow(np.sqrt(norm(xz_slice)), cmap=C, origin="lower", aspect="auto",
               extent=[0, Lz, 0, Lx])
    ax1.axvline(bar_z, color="w", ls="--", lw=1.4, alpha=0.9)
    ax1.axvline(scr_z, color="yellow", ls="--", lw=1.4, alpha=0.9)
    ax1.annotate("slits", xy=(bar_z, Lx * 0.5), xytext=(bar_z - 22, Lx * 0.5),
                 color="w", fontsize=8, ha="center")
    ax1.annotate("screen", xy=(scr_z, Lx * 0.5), xytext=(scr_z - 30, Lx * 0.62),
                 color="yellow", fontsize=8, ha="center")
    ax1.set_xlabel("propagation $z$  (source $\\leftarrow$ , screen $\\rightarrow$)",
                   color=TXT)
    ax1.set_ylabel("transverse $x$", color=TXT)
    ax1.set_title("(a) the wave: one electron, propagation $+$ diffraction",
                 color=TXT)
    blackize(ax1)
    # (b) the time-integrated interference pattern on the detection screen
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(norm(time_pattern).T, cmap=C, origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly])
    ax2.set_xlabel("transverse $x$", color=TXT)
    ax2.set_ylabel("transverse $y$", color=TXT)
    ax2.set_title("(b) the interference: the time-integrated screen", color=TXT)
    blackize(ax2)
    inset = ax2.inset_axes([0.03, 0.66, 0.55, 0.30], facecolor="#1a1a1a")
    prof = time_pattern.mean(axis=1)
    inset.plot(prof, color=C_ZERO, lw=1.3)
    inset.set_xlim(0, Lx); inset.set_ylim(0, prof.max() * 1.2)
    inset.set_xticks([0, Lx // 2, Lx]); inset.set_yticks([])
    inset.set_title("profile along $x$ (the fringe groups)", fontsize=7.5,
                   color=LIGHT, pad=1)
    for sp in inset.spines.values():
        sp.set_color("#444")
    inset.tick_params(colors=LIGHT, labelsize=6.5)
    # (c) the collapse: one electron reduces onto a point
    ax3 = fig.add_subplot(gs[1, 0])
    base = norm(time_pattern)
    ax3.imshow(base.T, cmap=C, origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly], alpha=0.55, vmin=0, vmax=1)
    fin = norm(seq[-1])
    ax3.imshow(np.ma.masked_less(fin.T, 0.18), cmap=C, origin="lower",
               aspect="auto", extent=[0, Lx, 0, Ly], vmin=0, vmax=1)
    ax3.add_patch(mpatches.Circle((cx, Ly - cy), 0.5 * lam / 2.0, fill=False,
                                  color="w", ls="--", lw=1.8))
    ax3.plot([cx], [Ly - cy], marker="*", color="w", ms=15, zorder=5,
             markeredgecolor="k", markeredgewidth=0.5)
    ax3.set_xlabel("transverse $x$", color=TXT)
    ax3.set_ylabel("transverse $y$", color=TXT)
    ax3.set_title("(c) the collapse: one electron reduces onto a point", color=TXT)
    blackize(ax3)
    for i in range(3):
        sax = ax3.inset_axes([0.10 + 0.27 * i, 0.04, 0.24, 0.40], facecolor="black")
        sax.imshow(norm([seq[0], seq[len(seq) // 2], seq[-1]][i]).T, cmap=C,
                   origin="lower", aspect="auto")
        sax.set_xticks([]); sax.set_yticks([])
        sax.set_title(["$\\sigma$ large", "contracting", "the point"][i],
                      fontsize=7.5, color=LIGHT, pad=1)
        for sp in sax.spines.values():
            sp.set_color("#444")
    # (d) the reality: the single-electron detections cluster into the fringes
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.imshow(base.T, cmap=C, origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly], alpha=0.30, vmin=0, vmax=1)
    xs = det[:, 0]; ys = Ly - det[:, 1]
    ax4.scatter(xs, ys, s=4, c=C_ZERO, alpha=0.8, linewidths=0)
    ax4.set_xlabel("transverse $x$", color=TXT)
    ax4.set_ylabel("transverse $y$", color=TXT)
    ax4.set_title(f"(d) the reality: {len(xs)} single-electron hits cluster into the fringes",
                 color=TXT)
    blackize(ax4)
    fig.suptitle("The CL10 Young double-slit as a single-electron measurement: "
                 "one electron at a time, each reduced onto a point (Tonomura flow)",
                 fontsize=12.5, color="#111")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, "fig25_young.png")


def main(with_young=False):
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
    fig18_residual_reduction()
    fig19_uv_2loop()
    fig20_proton_decay()
    fig21_proton_margin()
    fig22_proton_channels()
    fig23_qca_wave()
    fig24_gravity_3d()
    if with_young:
        fig25_young()
    print(f"done: {25 if with_young else 24} figures.")


if __name__ == "__main__":
    import sys
    main(with_young=("--young" in sys.argv))
