"""Figure module: proton."""
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from .style import ROOT, save
from .style import C_ROOT, C_SU2, C_SU3, C_ZERO


def fig16_proton_channels():
    """The partial lifetimes of the three GUT channels p -> e+ pi0,
    p -> mu+ pi0, p -> e+ K+ at the desert point (V64: M_X = 1.8e16 GeV,
    tau_p = 1.8e36 yr).  The rate is set by the two-body phase space
    (Gamma ~ pf, so tau ~ 1/pf; the hadronic matrix element is O(1))."""
    import json
    with open(os.path.join(ROOT, "experiments", "V64_desert", "results",
                           "v64_desert.json")) as fh:
        tau0 = json.load(fh)["data"]["selected"]["tau_p_yr"]
    mp, me, mmu, mpi, mk = 938.272, 0.511, 105.658, 134.977, 493.677
    def pf(a, b):
        return 0.5 / mp * np.sqrt((mp**2 - (a + b)**2) * (mp**2 - (a - b)**2))
    pf0 = pf(me, mpi)
    chans = [
        ("$p \\to e^+ \\pi^0$", pf(me, mpi) / pf0, 1.0, 2.4e34),
        ("$p \\to \\mu^+ \\pi^0$", pf(mmu, mpi) / pf0, 1.0, 1.6e34),
        ("$p \\to e^+ K^+$", pf(me, mk) / pf0, 1.0, 5.9e33),
    ]
    names, margin = [], []
    for nm, phase, had, bound in chans:
        tau_ch = tau0 * pf0 / phase    # tau ~ 1/rate ~ 1/(phase space)
        names.append(nm)
        margin.append(tau_ch / bound)    # how many times longer than the bound

    fig, ax = plt.subplots(figsize=(5.8, 3.7))
    x = np.arange(len(names))
    bars = ax.bar(x, margin, 0.5, color=[C_ROOT, C_SU2, C_SU3], alpha=0.85, edgecolor="k")
    for b, m in zip(bars, margin):
        ax.text(b.get_x() + b.get_width()/2, min(m, 8e3), f"{m:.0f}",
                ha="center", va="top", fontsize=12, color="w", weight="bold")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel(r"$\tau_p^{\rm (pred)}/\tau_p^{\rm (bound)}$  (the stability margin, log)")
    ax.set_title("Proton decay: the desert point (V64, $M_X = 1.8\\times10^{16}$ GeV)\n"
                 "against the three channel bounds — F5 is live",
                 fontsize=11.5)
    ax.axhline(1.0, color=C_ZERO, ls="--", lw=1.4,
               label=r"the bound  ($\Gamma=\Gamma_{\rm bound}$)")
    ax.set_ylim(1e0, 1e4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              framealpha=0.95, fontsize=10.5)
    fig.subplots_adjust(bottom=0.30)
    save(fig, "fig16_proton_channels.png")
