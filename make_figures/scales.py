"""Figure module: scales."""

import matplotlib.pyplot as plt
import numpy as np

from .style import save
from .style import C_FIT, C_GUT, C_ROOT, C_SU2, C_SU3, C_ZERO


def fig17_ladder():
    from cl10 import dm_wimp as dw, u10pred as tp

    LOOP = 1.0 / (16.0 * np.pi)
    MP = tp.M_PLANCK
    # (energy GeV, label, color, side(+/-1), height multiplier)
    scales = [
        (MP * LOOP ** 17, "neutrino ladder  $m_\\nu = 0.146$ eV  (O2, F3)",
         C_ZERO, -1, 0.5),
        (0.0487e-9, "$\\frac{1}{3}$ cand.  $0.0487$ eV", C_FIT, -1, 0.95),
        (0.2, r"QCD  $\Lambda_{\rm QCD} \sim 0.2$ GeV  (reference only)", C_SU3, 1, 0.5),
        (dw.dm_wimp()["m_DM_GeV"], r"DM  $118$ GeV  (O7)",
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
    save(fig, "fig17_ladder_of_scales.png")

def fig18_residual_reduction():
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
    axB.set_ylim(0, 3.4)
    axB.tick_params(axis="x", labelsize=8)
    axB.annotate(f"$n = {n_best}$: nearest", xy=(str(n_best), worst[n_best]),
                 xytext=(19.5, 3.15), ha="right",
                 fontsize=11.5, color=C_ZERO,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none",
                           alpha=0.92),
                 arrowprops=dict(arrowstyle="->", color=C_ZERO, lw=1.1))
    axB.annotate(f"$n = {n_best - 1}$: {worst[n_best - 1]:.1f} decades away",
                 xy=(str(n_best - 1), worst[n_best - 1]), xytext=(1.2, 3.15),
                 fontsize=11.5, color="#e07b00",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none",
                           alpha=0.92),
                 arrowprops=dict(arrowstyle="->", color="#e07b00", lw=1.1))
    axB.set_title(f"The scan $n = 1\\dots20$: $n = {n_best}$ is nearest\n"
                  f"($n = {n_best - 1}$ is {worst[n_best - 1]:.1f} decades away)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.subplots_adjust(bottom=0.34)
    save(fig, "fig18_residual_reduction.png")
