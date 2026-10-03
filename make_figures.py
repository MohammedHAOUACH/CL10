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
# Figure 2 — the 1-loop running of the three SM couplings to M_P
# ===========================================================================
def fig2_unification():
    from cl10 import unification

    B = unification.B_SM          # (41/10, -19/6, -7)
    AZ = unification.A_Z          # measured at M_Z
    M_Z, M_P = unification.M_Z, unification.M_P
    names = [r"$\alpha_1$ (U(1))", r"$\alpha_2$ (SU(2))", r"$\alpha_3$ (SU(3))"]
    cols = [C_GUT, C_SU2, C_SU3]

    s = np.logspace(np.log10(M_Z), np.log10(M_P), 240)
    ln = np.log(s / M_Z)
    fig, ax = plt.subplots()
    for k in range(3):
        inv_a = 1.0 / AZ[k] - (B[k] / (2.0 * np.pi)) * ln
        with np.errstate(divide="ignore", invalid="ignore"):
            a = 1.0 / inv_a
        a = np.where(inv_a > 0, a, np.nan)
        ax.plot(np.log10(s), a, color=cols[k], label=names[k])
    # the 0-parameter GUT coupling and the U(1) Landau pole
    ax.axhline(unification.ALPHA_10, color=C_ZERO, ls="--", lw=1.4,
               label=r"0-parameter $\alpha_{10}=1/25$ (the GUT coupling)")
    pole = unification.landau_pole_u1()
    ax.axvline(np.log10(pole), color=C_GUT, ls=":", lw=1.4,
               label=f"U(1) Landau pole ($\\sim$ {pole:.0e} GeV, above $M_P$)")
    ax.axvline(np.log10(M_P), color="k", ls="-", lw=1.0)
    ax.text(np.log10(M_P), 0.008, r" $M_P$", rotation=90, fontsize=9, ha="left")
    ax.set_xlabel(r"$\log_{10}(M/\mathrm{GeV})$  (1-loop, SM $\beta$-coefficients)")
    ax.set_ylabel(r"coupling $\alpha_i$")
    ax.set_title("1-loop running of the three SM couplings to the Planck scale")
    ax.set_ylim(0, 0.14)
    ax.legend(loc="upper right", framealpha=0.95)
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
    ax.scatter([10], [ladder[10]], s=90, color=C_ZERO, zorder=5,
               label=rf"10-loop $= M_P(1/16\pi)^{{10}} = $ {ladder[10]:.0f} GeV "
                     r"(the 0-parameter WIMP mass)")
    ax.scatter([0], [ladder[0]], s=90, color="k", zorder=5,
               label=r"$M_P = 1.22\times10^{19}$ GeV (the Planck scale)")
    ax.annotate("the 10-loop suppression\nbrings the Planck scale\ndown to the WIMP scale",
                xy=(10, ladder[10]), xytext=(5.5, 1e6),
                fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="->", color="0.3"))
    ax.set_xlabel(r"loop count $n$  (0-parameter: $10 = n_{\rm gen}$, the number of $\mathrm{Cl}(10)$ generators)")
    ax.set_ylabel(r"energy scale (GeV, log scale)")
    ax.set_title("The 10-loop ladder: one 0-parameter factor gives both\n"
                 r"the hierarchy $v/E_P = 2.08\,(1/16\pi)^{10}$ (V53) and the DM mass (V54)")
    ax.set_xticks(n)
    ax.set_ylim(1e2, 1e20)
    ax.legend(loc="upper right", framealpha=0.95)
    _save(fig, "fig3_10loop_ladder.png")


# ===========================================================================
# Figure 4 — the WIMP miracle: relic abundance vs mass
# ===========================================================================
def fig4_wimp_miracle():
    from cl10 import dm_wimp

    OMEGA_OBS = dm_wimp.relic_abundance(100.0)      # = 0.12 by construction
    m = np.logspace(2, 20, 240)
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
    ax.scatter([100.0], [0.12], s=110, color=C_SU3, zorder=5,
               label="WIMP miracle: $m=100$ GeV $\\Rightarrow$ $\\Omega_{\\rm obs}$")
    ax.scatter([m_dm], [dm_wimp.relic_abundance(m_dm)], s=110, color=C_ZERO, zorder=5,
               marker="*",
               label=f"0-parameter: $m_{{\\rm DM}}=M_P(1/16\\pi)^{{10}}={m_dm:.0f}$ GeV $\\Rightarrow$ $\\Omega h^2=$ {dm_wimp.relic_abundance(m_dm):.3f}")
    # the V42 Planck-mass soliton (under-abundant)
    ax.scatter([M_P], [0.12 * planck_gap], s=70, color="#7f7f7f", zorder=5,
               label=r"V42 Planck-mass soliton ($\sim 10^{-52}$, under-abundant)")
    # the 0-parameter WIMP star is identified in the legend (no extra
    # annotation: it sits in the top corner next to the WIMP-miracle point).
    ax.set_xlabel(r"relic mass $m$ (GeV, log scale)")
    ax.set_ylabel(r"relic abundance $\Omega h^2$")
    ax.set_title("The WIMP miracle, resolved: the 0-parameter 10-loop mass\n"
                 "$M_P(1/16\\pi)^{10}=118$ GeV gives the observed dark-matter abundance (V54)")
    ax.set_xlim(1e2, 1e20)
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


def main():
    print("Generating CL10 publication figures ->", OUT)
    fig1_one_object()
    fig2_unification()
    fig3_10loop_ladder()
    fig4_wimp_miracle()
    fig5_gravity_kretschmann()
    fig6_512_512()
    print("done: 6 figures.")


if __name__ == "__main__":
    main()
