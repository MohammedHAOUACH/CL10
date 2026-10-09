"""Shared style: rcParams, the provenance colour system, the save helper."""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "doc", "figures")

plt.rcParams.update({
    "figure.figsize": (5.4, 3.8),
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "axes.axisbelow": True,
    "axes.unicode_minus": False,
    "lines.linewidth": 1.8,
})
os.makedirs(OUT, exist_ok=True)

# consistent 0-parameter palette
C_GUT = "#1f77b4"     # blue
C_SU2 = "#ff7f0e"     # orange
C_SU3 = "#2ca02c"     # green
C_ROOT = "#7f1c8b"    # purple (the one object / the root)
C_ZERO = "#d62728"    # red (the 0-parameter points / the measured)

# the provenance palette (Section 14)
C_ALG = "#1f4e79"     # algebraic (derived, 0-parameter)
C_CONV = "#e67e22"    # conventional (mapping choices)
C_FIT = "#c0392b"     # fitted coefficients
C_OPEN = "#7f8c8d"    # open residuals


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote", path)
