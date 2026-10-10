"""Figure module: scales."""
import math
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from .style import ROOT, save
from .style import C_FIT, C_GUT, C_ROOT, C_SU2, C_SU3, C_ZERO


def fig2_wimp_miracle():
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
              label="WIMP miracle: $m=100$ GeV")
    ax.scatter([m_dm], [dm_wimp.relic_abundance(m_dm)], s=110, color=C_ZERO, zorder=5,
               clip_on=False, marker="*",
              label=f"CL10 mass: $m_{{\\rm DM}}={m_dm:.0f}$ GeV "
                    r"($\Omega h^2$ not computed)")
    # the Planck-mass soliton (under-abundant)
    ax.scatter([M_P], [0.12 * planck_gap], s=70, color="#7f7f7f", zorder=5, clip_on=False,
              label=r"Planck-mass soliton ($\sim 10^{-52}$)")
    # the 0-parameter WIMP star is identified in the legend (no extra
    # annotation: it sits in the top corner next to the WIMP-miracle point).
    ax.set_xlabel(r"relic mass $m$ (GeV, log scale)")
    ax.set_ylabel(r"relic abundance $\Omega h^2$")
    ax.set_title("The relic curve and the CL10 mass: the 118 GeV point sits on the\n"
                 "reference freeze-out curve — the abundance is not computed (W7)")
    ax.set_xlim(40.0, 1e20)
    ax.set_ylim(1e-55, 1.0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2,
              framealpha=0.95, fontsize=9.5)
    fig.subplots_adjust(bottom=0.30)
    save(fig, "fig2_wimp_miracle.png")

def fig4_qcd():
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
    ax.set_title("The area-law linear potential $V = \\sigma r$ rises without bound\n"
                 r"($\sigma = 1$ in Planck units $=M_P^2$: the QCD scale is not reproduced, Section 13.2)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30),
              framealpha=0.95, fontsize=10.5)
    fig.subplots_adjust(bottom=0.32)
    save(fig, "fig4_qcd_confinement.png")

def fig5_baryon():
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
                f"{v:.1e}\n(ratio {v / obs:.2f})", ha="center", fontsize=11, color=c)
    ax.axhline(obs / 1e-9, color="k", ls="--", lw=1.2,
               label=r"the observed $\Delta B = 6.1\times10^{-10}$")
    ax.set_ylabel(r"baryon asymmetry $\Delta B$  ($\times 10^{-9}$)")
    ax.set_title(r"GUT leptogenesis: $\Delta B = \frac{28}{79}\,\varepsilon_{CP}\,D$"
                 "\n"
                 r"(0-parameter: $6.8\times10^{-10}$, a 12% match;"
                 " the named $\\eta_{CP}$: 2.8)",
                 fontsize=12)
    ax.set_ylim(0, max(vals9) * 1.4)
    ax.legend(loc="upper left", framealpha=0.95)
    save(fig, "fig5_baryon_asymmetry.png")

def fig7_ladder():
    from cl10 import dm_wimp as dw, u10pred as tp

    LOOP = 1.0 / (16.0 * np.pi)
    MP = tp.M_PLANCK
    # (energy GeV, label, color, side(+/-1), height multiplier)
    scales = [
        (MP * LOOP ** 17, "neutrino ladder  $m_\\nu = 0.146$ eV  (W2, F3)",
         C_ZERO, -1, 0.5),
        (0.0487e-9, "$\\frac{1}{3}$ cand.  $0.0487$ eV", C_FIT, -1, 0.95),
        (0.2, r"QCD  $\Lambda_{\rm QCD} \sim 0.2$ GeV  (reference only)", C_SU3, 1, 0.5),
        (dw.dm_wimp()["m_DM_GeV"], r"DM  $118$ GeV  (W7)",
         C_SU2, -1, 0.95),
        (tp.V_EW, r"electroweak  $v = 246$ GeV  (measured)", C_GUT, 1, 0.95),
        (MP * LOOP ** 3, r"seesaw  $M_N = 9.6\times10^{13}$ GeV  ($n = 3$)",
         C_ROOT, -1, 0.5),
        (MP, r"Planck  $M_P = 1.22\times10^{19}$ GeV  ($n = 0$)", "k", 1, 0.95),
    ]
    fig, ax = plt.subplots(figsize=(11.0, 4.4))
    ax.axhline(0, color="0.3", lw=1.0)
    for (E, lab, col, side, hgt) in scales:
        x = np.log10(E)
        ax.plot([x, x], [0, 0.30 * hgt * side], color=col, lw=1.0)
        ax.scatter([x], [0], s=90, color=col, zorder=5, clip_on=False,
                   edgecolor="k", lw=0.5)
        ha = "left" if x < -9 else ("right" if x > 17 else "center")
        ax.text(x, 0.37 * hgt * side, lab, ha=ha,
                va="bottom" if side > 0 else "top", fontsize=11, color=col)
    ax.set_xlim(-13.5, 20.5)
    ax.set_ylim(-0.72, 0.92)
    ax.set_yticks([])
    ax.set_xlabel("energy scale  (GeV, log scale)")
    ax.set_title("The ladder of scales: one factor $(1/16\\pi)^n$, four scales spanning 30 decades\n"
                 r"(QCD is shown for reference only: it is not produced by the ladder)", fontsize=13)
    save(fig, "fig7_ladder_of_scales.png")

def fig9_sm_head2head():
    import json
    from cl10 import eta_cp as ep, dm_wimp as dw, hier_exponent as he
    e = ep.eta_cp()
    d = dw.dm_wimp()
    h = he.hierarchy17()
    with open(os.path.join(ROOT, "experiments", "V64_desert", "results",
                           "v64_desert.json")) as fh:
        tau_desert = json.load(fh)["data"]["selected"]["tau_p_yr"]

    def _pair(ax, title, unit, obs, pred, log, note=None, ok=None,
              left="SM"):
        xs = [0, 1]
        vals = [obs, pred]
        bars = ax.bar(xs, vals, color=["#9a9a9a", C_ZERO],
                      edgecolor="k", lw=0.6, width=0.62)
        if log:
            ax.set_yscale("log")
            ax.set_ylim(min(vals) * 0.35, max(vals) * 2.6)
            ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        else:
            ax.set_ylim(0, max(vals) * 1.20)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v * 0.97, f"{v:g}",
                    ha="center", va="top", fontsize=10.5, color="0.25")
        ratio = pred / obs
        if ok is None:
            ok = 0.5 <= ratio <= 2.0
        col = C_SU3 if ok else "#e07b00"
        txt = note if note is not None else f"CL10 / obs = {ratio:.2f}"
        ax.text(0.5, -0.28, txt, transform=ax.transAxes, ha="center", va="top",
                fontsize=10.5, color=col,
                bbox=dict(boxstyle="round", fc="#f6f6f6", ec="#c8c8c8", alpha=0.92))
        ax.set_xticks(xs)
        ax.set_xticklabels([left, "CL10"], fontsize=11.5)
        ax.set_title(title + (f"  ({unit})" if unit else ""), fontsize=12)
        ax.set_ylabel("value" + (f", {unit}" if unit else "") + (" (log)" if log else ""),
                      fontsize=11)

    fig, axs = plt.subplots(2, 3, figsize=(11.2, 6.0))
    fig.suptitle("CL10 prediction  vs  measured value\n"
                 "sources: Planck 2018 + PDG 2024, NuFIT 6.0, Super-Kamiokande "
                 "(provenance per panel, Section 14)",
                 fontsize=14, y=0.995)
    _pair(axs[0, 0], "baryon asymmetry  $\\Delta B$", "$\\times10^{-10}$",
          6.12, 6.8, False, left="Planck",
          note="coincidental product (W5);\nwashout at self-consistent $y$:\n$K = 47$, in equilibrium (V106)")
    _pair(axs[0, 1], "dark-matter mass  $m_{\\rm DM}$", "GeV",
          100.0, d["m_DM_GeV"], False, left="weak scale",
          note="mass coincidence (W7);\n$\\Omega h^2$ not computed")
    _pair(axs[0, 2], "hierarchy  $v/E_P$", "",
          246.22 / 1.2209e19, h["predicted_ratio"], True, left="PDG",
          note=r"2.08$\times$ below obs (the 2-digit residual)", ok=False)
    _pair(axs[1, 0], r"neutrino mass  $m_\nu$ (ladder)", "eV",
          0.0501, 0.146, True, left="NuFIT",
          note="2.9$\\times$ obs; excluded by\n$\\sum m_\\nu < 0.12$ eV (W2, F3)", ok=False)
    _pair(axs[1, 1], r"$\sum m_\nu$: $\frac{1}{3}$ candidate", "eV",
          0.12, 0.0487, True, left="Planck+BAO bound",
          note="inside the bound; finite window,\nstructure protected (V103)", ok=True)
    _pair(axs[1, 2], "proton lifetime  $\\tau_p$ (desert)", "yr",
          2.4e34, tau_desert, True, left="SK bound",
          note=f"{tau_desert / 2.4e34:.0f}$\\times$ above the bound\n(F5 live)", ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    save(fig, "fig9_sm_head2head.png")

def fig12_ladder_sm():
    from cl10 import dm_wimp as dw, u10pred as tp

    M_P = tp.M_PLANCK
    LOOP = 1.0 / (16.0 * np.pi)
    # each point: (GeV value, label, color, y_level, ha, x_nudge) -- the level
    # + nudge + a leader line stagger the cluster near log10(100-246 GeV)
    cl10 = [
        (M_P * LOOP ** 17,   r"neutrino  $m_\nu$  (W2)", C_ZERO, 0, "left", 0.0),
        (1.0,               r"QCD  $\Lambda \sim 1$ GeV",                 C_SU3,  0, "center", 0.0),
        (dw.dm_wimp()["m_DM_GeV"], r"DM  118 GeV  (W7)", C_SU2, 1, "right", -0.4),
        (tp.V_EW,           r"electroweak  $v = 246$ GeV",                C_GUT,  2, "left",  +0.4),
        (M_P * LOOP ** 3,   r"2-step  $M_N$  (W2)",     C_ROOT, 1, "center", 0.0),
        (M_P,              r"Planck  $M_P = 1.22\times10^{19}$ GeV",      "k",    0, "right", 0.0),
    ]
    sm = [
        (0.05 * 1e-9,      r"$m_\nu \sim 0.05$ eV",            "#8c8c8c", 0, "left", 0.0),
        (1.0,             r"$\Lambda_{\rm QCD} \sim 1$ GeV",   "#8c8c8c", 0, "center", 0.0),
        (100.0,           r"WIMP  $\sim 100$ GeV",             "#8c8c8c", 1, "right", -0.4),
        (246.0,           r"$v = 246$ GeV",                    "#8c8c8c", 2, "left",  +0.4),
        (1.22e19,         r"$M_P$",                            "#8c8c8c", 0, "center", 0.0),
    ]

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
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
                    ha=ha, va=va, fontsize=12.5, color=col)

    draw_row(cl10, 0.30, top_levels)
    draw_row(sm, -0.30, bot_levels)

    ax.text(-10.8, 1.05, "CL10 ladder", fontsize=14, weight="bold",
            color="0.2", ha="left")
    ax.text(-10.8, -1.09, "SM observed ladder", fontsize=14, weight="bold",
            color="0.4", ha="left")
    ax.set_xlim(-11, 20)
    ax.set_ylim(-1.18, 1.18)
    ax.set_yticks([])
    ax.set_xlabel("energy scale (GeV, log)   —   one object spans 30+ decades",
                  fontsize=12)
    ax.set_title("The CL10 ladder vs the SM ladder: the same hierarchy,"
                 " from the Planck scale down to the neutrino", fontsize=13)
    save(fig, "fig12_ladder_sm.png")

def fig13_residual_reduction():
    from cl10 import four_limits as fl

    r = fl.unique_10loop()
    h_scan = [(p["n"], p["hier_dev"]) for p in r["scan"]]
    d_scan = [(p["n"], p["dm_dev"]) for p in r["scan"]]
    n_best = r["n_best"]

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.4, 1.0]})
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
    axA.set_xlabel("loop count $n$  (scan $n = 1\\dots20$)")
    axA.set_ylabel("deviation from observed  (e-folds)")
    axA.set_xlim(0, 21)
    axA.set_ylim(-0.5, 2.4)
    axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.36),
               framealpha=0.95, fontsize=10.5)
    axA.set_title("One loop integer, one scale: the hierarchy and the DM mass are\n"
                  r"one equation, $v = 2.08\,m_{\rm DM}$ — not two independent confirmations",
                  fontsize=12)
    # (b) the scan n = 1..20: n = 10 nearest, n = 9 is 1.4 decades away
    worst = {n: max(h, dv) for (n, h), (nn2, dv) in zip(h_scan, d_scan)
             if n == nn2}
    ns = sorted(worst)
    axB.bar([str(n) for n in ns], [worst[n] for n in ns],
            color=[C_ZERO if n == n_best else
                   ("#e07b00" if n == n_best - 1 else "0.7") for n in ns],
            width=0.7)
    axB.set_xlabel(r"scan $n$  (distinct from the coupling scan $\alpha = 1/n^2$)")
    axB.set_ylabel(r"worst deviation of the two scales (e-folds)")
    axB.set_ylim(0, 2.75)
    axB.annotate(f"$n = {n_best}$: nearest", xy=(str(n_best), worst[n_best]),
                 xytext=(19.5, 2.5), ha="right",
                 fontsize=11.5, color=C_ZERO,
                 arrowprops=dict(arrowstyle="->", color=C_ZERO, lw=1.1))
    axB.annotate(f"$n = {n_best - 1}$: 1.4 decades away",
                 xy=(str(n_best - 1), worst[n_best - 1]), xytext=(1.2, 2.5),
                 fontsize=11.5, color="#e07b00",
                 arrowprops=dict(arrowstyle="->", color="#e07b00", lw=1.1))
    axB.set_title("The scan $n = 1\\dots20$: $n = 10$ is nearest\n"
                  "($n = 9$ is 1.4 decades away)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.subplots_adjust(bottom=0.34)
    save(fig, "fig13_residual_reduction.png")
