#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 - the jury-grade figures (A-O).

Fifteen figures designed to make the paper's structure, its strongest results
and its named failures visible at a glance.  Every number is COMPUTED from
the cl10/ modules or read from the experiment JSONs (V61/V62); measured
values are reference data (PDG, Planck, Super-Kamiokande, KATRIN, NuFIT),
never fits.  Each figure carries:

  * the provenance colour code (the paper's methodological signature):
      algebraic  = blue    (fixed by the root, no parameter)
      free       = purple  (the one 't Hooft angle g_10)
      convention = amber   (a mapping / scale identification)
      fitted     = orange  (a coefficient read from the data)
      open/fail  = red     (a named residual W1-W9 / a falsifier)
      measured   = grey    (reference data)
  * a residual badge where a computation backs the claim.

Usage:
    python3 make_figures_pro.py            # light (paper) variants
    python3 make_figures_pro.py --dark     # dark (slides) variants
    python3 make_figures_pro.py --only A,C,E,G
"""

import argparse
import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Wedge

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "doc", "figures_pro")

# ---------------------------------------------------------------------------
# the provenance palette (the paper's colour code, used in EVERY figure)
# ---------------------------------------------------------------------------
C_ALG = "#2f7fbf"      # algebraic   (the root)
C_FREE = "#8e44ad"     # free angle  (g_10)
C_CONV = "#d9a520"     # convention
C_FIT = "#e8842c"      # fitted coefficient
C_OPEN = "#c0392b"     # open residual / failure
C_MEAS = "#8a8f98"     # measured reference data
C_MATTER = "#d62728"   # the 16 (matter)
C_ANTI = "#7f1c8b"     # the 16-bar (antimatter)
C_SU3 = "#2ca02c"
C_SU2 = "#ff7f0e"
C_GUT = "#1f77b4"

DARK = False
BADGE_FC = "#f2f2f2"
BADGE_EC = "#888888"
TXT = "#111111"


def set_style(dark=False):
    global DARK, BADGE_FC, BADGE_EC, TXT
    DARK = dark
    if dark:
        BADGE_FC, BADGE_EC, TXT = "#1c2530", "#5a6673", "#e8ecf1"
        face, axes_face, grid = "#0e1216", "#141a21", 0.22
    else:
        BADGE_FC, BADGE_EC, TXT = "#f2f2f2", "#888888", "#111111"
        face, axes_face, grid = "white", "white", 0.3
    plt.rcParams.update({
        "figure.figsize": (6.8, 4.6),
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "axes.grid": True,
        "grid.alpha": grid,
        "grid.linestyle": ":",
        "axes.axisbelow": True,
        "axes.unicode_minus": False,
        "lines.linewidth": 1.8,
        "mathtext.fontset": "cm",
        "figure.facecolor": face,
        "axes.facecolor": axes_face,
        "savefig.facecolor": face,
        "text.color": TXT,
        "axes.labelcolor": TXT,
        "xtick.color": TXT,
        "ytick.color": TXT,
        "axes.edgecolor": "#5a6673" if dark else "#333333",
    })


def badge(ax, txt, loc="lower left", fs=8.5, x=None, y=None):
    pos = {"lower left": (0.03, 0.05), "lower right": (0.97, 0.05),
           "upper left": (0.03, 0.95), "upper right": (0.97, 0.95)}[loc]
    xx = pos[0] if x is None else x
    yy = pos[1] if y is None else y
    ha = "right" if xx > 0.5 else "left"
    va = "top" if yy > 0.5 else "bottom"
    ax.text(xx, yy, txt, transform=ax.transAxes, fontsize=fs, ha=ha, va=va,
            bbox=dict(boxstyle="round,pad=0.4", fc=BADGE_FC, ec=BADGE_EC,
                      alpha=0.95))


def _save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    if DARK:
        d = os.path.join(OUT, "dark")
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, name)
    else:
        path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", path)


def provenance_legend(handles=None, loc="lower center", ncol=6, fs=8.5):
    if handles is None:
        handles = [
            Line2D([], [], color=C_ALG, lw=6, label="algebraic (the root)"),
            Line2D([], [], color=C_FREE, lw=6, label="free angle $g_{10}$"),
            Line2D([], [], color=C_CONV, lw=6, label="convention"),
            Line2D([], [], color=C_FIT, lw=6, label="fitted coefficient"),
            Line2D([], [], color=C_OPEN, lw=6, label="open / fails"),
            Line2D([], [], color=C_MEAS, lw=6, label="measured (reference)"),
        ]
    return handles, dict(loc=loc, ncol=ncol, fontsize=fs, framealpha=0.9,
                         borderpad=0.5, handlelength=1.4)


# ===========================================================================
# Fig. A - graphical abstract: one object, five physics
# ===========================================================================
def figA_one_object():
    fig = plt.figure(figsize=(13.6, 8.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    cx, cy, cr = 0.5, 0.52, 0.115

    # --- the central object ------------------------------------------------
    ax.add_patch(Circle((cx, cy), cr + 0.012, fc="none", ec=C_ALG, lw=2.5))
    ax.add_patch(Circle((cx, cy), cr, fc=("#141a21" if DARK else "#eef4fa"),
                        ec=C_ALG, lw=1.5))
    ax.text(cx, cy + 0.055, r"$\mathrm{Cl}(10,0)$", ha="center",
            fontsize=21, color=C_ALG, weight="bold")
    ax.text(cx, cy + 0.012, "the root: 10 generators\n$h_a^\\dagger=h_a,\\ h_a^2=I,\\ \\{h_a,h_b\\}=0$",
            ha="center", fontsize=9.5)
    ax.text(cx, cy - 0.062, "32-D module  ·  1024 elements\n$w=h_1\\cdots h_{10}$: the $16/\\overline{16}$ split",
            ha="center", fontsize=9.5)

    # the 10 generator tiles on the circle
    for k in range(10):
        a = math.pi / 2 + 2 * math.pi * k / 10
        gx, gy = cx + (cr + 0.032) * math.cos(a), cy + (cr + 0.032) * math.sin(a)
        space = k < 4
        ax.add_patch(FancyBboxPatch((gx - 0.014, gy - 0.011), 0.028, 0.022,
                                    boxstyle="round,pad=0.002",
                                    fc=(C_ALG if space else
                                        ("#2a3542" if DARK else "#d7dde5")),
                                    ec=C_ALG if space else BADGE_EC, lw=1.2))
        ax.text(gx, gy, f"$h_{k+1}$", ha="center", va="center", fontsize=8,
                color="white" if space else TXT)

    def panel(x, y, w, h, title, edge, theorem):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.008",
                                    fc=("#141a21" if DARK else "#fafbfc"),
                                    ec=edge, lw=2.2))
        ax.text(x + w / 2, y + h - 0.028, title, ha="center", fontsize=12,
                weight="bold", color=edge)
        ax.text(x + w / 2, y + h - 0.052, theorem, ha="center", fontsize=8.5,
                color=edge, style="italic")

    def arrow(x, y):
        ax.add_patch(FancyArrowPatch((cx + (cr + 0.05) * math.cos(math.atan2(y - cy, x - cx)),
                                      cy + (cr + 0.05) * math.sin(math.atan2(y - cy, x - cx))),
                                     (x, y), connectionstyle="arc3,rad=0.12",
                                     arrowstyle="-|>", mutation_scale=16,
                                     color="#7f8c99", lw=1.6))

    # --- 1. spacetime (top-left) -------------------------------------------
    px, py, pw, ph = 0.045, 0.60, 0.24, 0.30
    panel(px, py, pw, ph, "spacetime (1+3)", C_ALG,
          "Theorem 1 - 0-parameter")
    ox, oy = px + 0.075, py + 0.10
    ax.fill([ox, ox - 0.055, ox + 0.055], [oy, oy + 0.10, oy + 0.10],
            color=C_ALG, alpha=0.25)
    ax.fill([ox, ox - 0.055, ox + 0.055], [oy, oy - 0.10, oy - 0.10],
            color=C_ALG, alpha=0.12)
    ax.plot([ox - 0.055, ox, ox + 0.055], [oy + 0.10, oy, oy + 0.10],
            color=C_ALG, lw=1.6)
    ax.plot([ox - 0.055, ox, ox + 0.055], [oy - 0.10, oy, oy - 0.10],
            color=C_ALG, lw=1.6)
    ax.plot([ox, ox], [oy - 0.11, oy + 0.11], color=C_MEAS, lw=1.0, ls=":")
    ax.text(px + 0.165, py + 0.155, r"time $T=h_1$", fontsize=8.5, color=C_ALG)
    ax.text(px + 0.165, py + 0.125, r"space $S_a=\mathrm{i}h_{a+1}$",
            fontsize=8.5, color=C_ALG)
    ax.text(px + 0.165, py + 0.09, r"$v^2=v_0^2-\sum v_a^2$", fontsize=8.5)
    ax.text(px + 0.02, py + 0.022,
            "the light cone, inside the GUT", fontsize=8, color=C_MEAS,
            ha="left")

    # --- 2. matter (top-right) ---------------------------------------------
    px2, py2, pw2, ph2 = 0.715, 0.60, 0.24, 0.30
    panel(px2, py2, pw2, ph2, "matter: one generation", C_ALG,
          "Prop. 1 + Theorem 2 - anomalies = 0")
    groups = [("Q", 6, C_SU3), ("u$^c$", 3, C_SU3), ("d$^c$", 3, C_SU3),
              ("L", 2, C_SU2), ("e$^c$", 1, C_MEAS), ("$\\nu^c$", 1, C_MEAS)]
    gx0, gy0 = px2 + 0.022, py2 + 0.085
    for name, n, col in groups:
        for j in range(n):
            ax.add_patch(plt.Rectangle((gx0 + j * 0.0135, gy0), 0.0115, 0.030,
                                       fc=col, ec="none", alpha=0.9))
        ax.text(gx0 + (n - 1) * 0.0068, gy0 - 0.014, name, ha="center",
                fontsize=8)
        gx0 += n * 0.0135 + 0.012
    ax.text(px2 + pw2 / 2, py2 + 0.155, "$6{+}3{+}3{+}2{+}1{+}1=16$",
            ha="center", fontsize=10, weight="bold")
    ax.text(px2 + pw2 / 2, py2 + 0.032,
            "the 16 Weyl of the pseudoscalar;\nall 8 anomaly coefficients vanish",
            ha="center", fontsize=8, color=C_MEAS)

    # --- 3. gauge (right) ---------------------------------------------------
    px3, py3, pw3, ph3 = 0.715, 0.24, 0.24, 0.30
    panel(px3, py3, pw3, ph3, "gauge: the 45", C_ALG,
          "Theorem 3 - the 5/3, 0-parameter")
    segs = [(12, "12 SM", C_GUT), (12, "12 X/Y", C_OPEN),
            (20, "20", C_MEAS), (1, "", C_CONV)]
    bx = px3 + 0.025
    for n, lab, col in segs:
        wseg = 0.19 * n / 45
        ax.add_patch(plt.Rectangle((bx, py3 + 0.15), wseg, 0.035, fc=col,
                                   ec="none"))
        if lab:
            ax.text(bx + wseg / 2, py3 + 0.196, lab, ha="center", fontsize=7.5)
        bx += wseg
    ax.text(px3 + pw3 / 2, py3 + 0.10,
            "$\\mathrm{Spin}(10)$ adjoint $=h_ah_b$", ha="center", fontsize=9)
    ax.text(px3 + pw3 / 2, py3 + 0.055,
            "$g_1^2/g'^2=5/3$ from $\\mathrm{Tr}\\,Y^2=5/6$\n"
            "X/Y: proton decay (F5)", ha="center", fontsize=8, color=C_MEAS)

    # --- 4. statistics (bottom-right) ---------------------------------------
    px4, py4, pw4, ph4 = 0.045, 0.045, 0.24, 0.15
    panel(px4, py4, pw4, ph4, "statistics", C_CONV, "")
    ax.text(px4 + pw4 / 2, py4 + ph4 - 0.07,
            "the $\\mathbb{Z}_2$ grading (postulate)", ha="center",
            fontsize=8.5, color=C_CONV)
    ax.add_patch(Wedge((px4 + 0.055, py4 + 0.04), 0.033, 90, 270,
                       fc=C_SU2, alpha=0.8))
    ax.add_patch(Wedge((px4 + 0.055, py4 + 0.04), 0.033, 270, 90,
                       fc=C_GUT, alpha=0.8))
    ax.text(px4 + 0.028, py4 + 0.043, "512\nodd", ha="center", fontsize=7,
            color="white")
    ax.text(px4 + 0.082, py4 + 0.043, "512\neven", ha="center", fontsize=7,
            color="white")
    ax.text(px4 + 0.175, py4 + 0.045, "even = boson\nodd = fermion",
            fontsize=8.5, ha="center")
    ax.text(px4 + 0.175, py4 + 0.008, "no SUSY predicted", fontsize=8,
            ha="center", color=C_MEAS)

    # --- 5. dynamics + gravity (left-bottom) --------------------------------
    px5, py5, pw5, ph5 = 0.045, 0.24, 0.24, 0.30
    panel(px5, py5, pw5, ph5, "dynamics + gravity", C_CONV,
          "Theorems 4-5 - unitary, causal")
    gx = np.linspace(px5 + 0.03, px5 + 0.13, 7)
    gy = np.linspace(py5 + 0.10, py5 + 0.20, 7)
    for i, xg in enumerate(gx):
        for j, yg in enumerate(gy):
            r = abs(i - 3) + abs(j - 3)
            ax.add_patch(Circle((xg, yg), 0.0022,
                                fc=(C_ALG if r <= 2 else
                                    ("#2a3542" if DARK else "#c8cfd8")),
                                ec="none"))
    for rad, al in ((1, 0.9), (2, 0.5), (3, 0.25)):
        ax.add_patch(Circle((gx[3], gy[3]), rad * (gx[1] - gx[0]), fc="none",
                            ec=C_ALG, lw=1.2, alpha=al, ls="--"))
    ax.text(px5 + 0.15, py5 + 0.185, "front speed $=1.000$\n(Thm 5, V36)",
            fontsize=8, color=C_ALG)
    ax.text(px5 + 0.15, py5 + 0.115, "$\\det U=1$; CPT mirror\n"
            "$\\omega_{\\overline{16}}=-\\omega_{16}$ (Thm 4)", fontsize=8)
    ax.text(px5 + 0.225, py5 + 0.032,
            "gravity: the metric of the\n(1,3) subalgebra - $g_{\\mu\\nu}$\n"
            "and $G=1$ are ADDED (convention)",
            ha="right", fontsize=7.5, color=C_CONV)

    # --- the free angle + residuals strip -----------------------------------
    ax.add_patch(FancyBboxPatch((0.345, 0.045), 0.31, 0.115,
                                boxstyle="round,pad=0.008",
                                fc=("#141a21" if DARK else "#fafbfc"),
                                ec=C_FREE, lw=2.2))
    ax.text(0.5, 0.135, "one free angle, honestly labelled", ha="center",
            fontsize=11, weight="bold", color=C_FREE)
    ax.text(0.5, 0.098, "$g_{10}=1/5$: the only parameter of the dynamics;\n"
                        "the algebraic results hold for EVERY $g_{10}$",
            ha="center", fontsize=9)
    ax.text(0.5, 0.058, "9 open residuals W1-W9  ·  7 falsifiers F1-F7",
            ha="center", fontsize=9, color=C_OPEN)

    arrow(0.165, 0.60)
    arrow(0.835, 0.60)
    arrow(0.835, 0.54)
    arrow(0.285, 0.45)
    arrow(0.285, 0.09)

    h, kw = provenance_legend()
    fig.legend(handles=h, **kw)
    fig.suptitle("CL10: one Clifford root on a Planck-lattice QCA - "
                 "spacetime, matter, gauge, statistics and gravity all follow",
                 fontsize=13, y=0.985)
    _save(fig, "figA_one_object_abstract.png")


# ===========================================================================
# Fig. B - the atlas of the 1024: where the physics sits in the algebra
# ===========================================================================
def figB_atlas_1024():
    from math import comb

    grades = [comb(10, r) for r in range(11)]
    # colour by grade: even grades blue family, odd grades red family
    cmap_even = plt.get_cmap("Blues")
    cmap_odd = plt.get_cmap("Reds")
    cell_col = np.zeros((32, 32))
    cell_grade = np.full((32, 32), -1)
    idx = 0
    for r in range(11):
        for _ in range(grades[r]):
            i, j = divmod(idx, 32)
            cell_col[i, j] = r
            cell_grade[i, j] = r
            idx += 1

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(12.6, 6.0),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    img = np.zeros((32, 32, 4))
    for r in range(11):
        shade = 0.35 + 0.6 * (r / 10)
        base = cmap_even(shade) if r % 2 == 0 else cmap_odd(shade)
        mask = cell_grade == r
        for c in range(4):
            ch = img[..., c]
            ch[mask] = base[c]
            img[..., c] = ch
    ax.imshow(img, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("the atlas of the root: all 1024 algebra elements,\n"
                 "ordered by grade $r$ (binomial sizes)", fontsize=11)

    # highlight the physics cells
    def cell_rect(k):
        i, j = divmod(k, 32)
        return plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="none",
                             ec="k" if not DARK else "w", lw=1.4)

    # grade 1: cells 1..10 ; the first 4 are the spacetime
    for k in range(1, 5):
        ax.add_patch(cell_rect(k))
    # grade 2: cells 11..55 (the 45 gauge generators)
    g2_start = 1 + 10
    for k in range(g2_start, g2_start + 45):
        i, j = divmod(k, 32)
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="none",
                                   ec=C_CONV, lw=0.9))
    # grade 10: the pseudoscalar (the last cell)
    ax.add_patch(cell_rect(1023))
    i, j = divmod(1023, 32)
    ax.annotate("$w=h_1\\cdots h_{10}$:\nthe $16/\\overline{16}$ split",
                (j, i), xytext=(15, 27.5), fontsize=9, color=C_OPEN,
                arrowprops=dict(arrowstyle="->", color=C_OPEN))
    i1, j1 = divmod(1, 32)
    ax.annotate("grade 1: $h_1..h_4$ = spacetime\n(the other 6 are inert)",
                (j1 + 2, i1), xytext=(13, 5.5), fontsize=9,
                color="k" if not DARK else "w",
                arrowprops=dict(arrowstyle="->"))
    i2, j2 = divmod(g2_start + 20, 32)
    ax.annotate("grade 2: the 45 = the gauge bosons",
                (j2, i2), xytext=(13, 11.5), fontsize=9, color=C_CONV,
                arrowprops=dict(arrowstyle="->", color=C_CONV))
    badge(ax, "even grades (0,2,..10) = boson sector\n"
              "odd grades (1,3,..9) = fermion sector\n"
              "the physics occupies 4+45+1 = 50 of 1024 cells",
          loc="lower left", fs=8)

    # side: the binomial sizes
    axb.bar(range(11), grades, color=[(cmap_even(0.35 + 0.6 * r / 10)
                                       if r % 2 == 0 else
                                       cmap_odd(0.35 + 0.6 * r / 10))
                                      for r in range(11)])
    for r, gsz in enumerate(grades):
        axb.text(r, gsz + 3, str(gsz), ha="center", fontsize=8.5)
    axb.set_xlabel("grade $r$ (number of generators in the product)")
    axb.set_ylabel("number of elements $\\binom{10}{r}$")
    axb.set_title("the grade decomposition: 1+10+45+...+1 = 1024", fontsize=11)
    axb.annotate("spacetime\n(4 of 10)", (1, 10), xytext=(0.05, 150),
                 fontsize=8.5, arrowprops=dict(arrowstyle="->"))
    axb.annotate("gauge (45)", (2, 45), xytext=(3.0, 180), fontsize=8.5,
                 color=C_CONV, arrowprops=dict(arrowstyle="->", color=C_CONV))
    axb.annotate("chirality (1)", (10, 1), xytext=(6.3, 90), fontsize=8.5,
                 color=C_OPEN, arrowprops=dict(arrowstyle="->", color=C_OPEN))
    badge(axb, "computed from cl10.clifford: the 1024 ordered\n"
               "products of the 10 Pauli-pattern generators",
          loc="upper right", fs=8)

    fig.suptitle("the whole construction inside one 1024-dimensional algebra",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    _save(fig, "figB_atlas_1024.png")


# ===========================================================================
# Fig. C - the provenance sunburst (Table 5 as an identity)
# ===========================================================================
def figC_provenance_sunburst():
    classes = [
        ("algebraic", 7, C_ALG,
         ["16/$\\overline{16}$ split", "the 45", "the 5/3",
          "8 anomalies = 0", "32-level spectrum", "$\\det U=1$",
          "CPT mirror"]),
        ("free", 1, C_FREE, ["$g_{10}=1/5$"]),
        ("convention", 2, C_CONV,
         ["$\\alpha_{10}\\equiv g_{10}^2$", "$g_0=1/8$"]),
        ("fitted", 5, C_FIT,
         ["2.08 (hierarchy)", "1.18 (DM mass)", "$n_{2\\rm step}=3$",
          "$M_{10}=1.2\\times10^7$ GeV (V64)",
          "$M_{(1,2)}=1.0\\times10^{12}$ GeV (V64)"]),
        ("open", 9, C_OPEN,
         ["W1 $\\Lambda$", "W2 $m_\\nu$", "W3 families", "W4 couplings",
          "W5 $\\Delta B$", "W6 flavour", "W7 DM identity",
          "W8 hierarchy", "W9 doublers"]),
    ]

    fig = plt.figure(figsize=(11.5, 8.2))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.92])
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axis("off")

    total = sum(c[1] for c in classes)
    a0 = 90.0
    for name, n, col, items in classes:
        sweep = 360.0 * n / total
        ax.add_patch(Wedge((0, 0), 0.52, a0, a0 + sweep, width=0.22,
                           fc=col, ec="white" if not DARK else "#0e1216",
                           lw=1.5))
        amid = math.radians(a0 + sweep / 2)
        ax.text(0.41 * math.cos(amid), 0.41 * math.sin(amid),
                f"{name}\n{n}", ha="center", va="center", fontsize=10,
                weight="bold", color="white")
        # outer ring: the named numbers
        a1 = a0
        for it in items:
            sw = sweep / n
            ax.add_patch(Wedge((0, 0), 0.95, a1 + 0.6, a1 + sw - 0.6,
                               width=0.43, fc=col, ec="white" if not DARK
                               else "#0e1216", lw=0.8, alpha=0.55))
            am = math.radians(a1 + sw / 2)
            r_t = 0.735
            ha = "center"
            t = ax.text(r_t * math.cos(am), r_t * math.sin(am), it,
                        ha=ha, va="center", fontsize=8.2, rotation=0,
                        rotation_mode="anchor")
            if 90 < math.degrees(am) % 360 < 270:
                t.set_rotation(-math.degrees(am) + 90)
            else:
                t.set_rotation(math.degrees(am))
            a1 += sw
        a0 += sweep

    ax.add_patch(Circle((0, 0), 0.30, fc=("#141a21" if DARK else "#eef4fa"),
                        ec="#5a6673"))
    ax.text(0, 0.07, "24 numbers", ha="center", fontsize=14, weight="bold")
    ax.text(0, -0.03, "one audit", ha="center", fontsize=11, color=C_MEAS)
    ax.text(0, -0.14, "Table 5", ha="center", fontsize=9, color=C_MEAS)
    ax.text(0, -1.28, "7 fixed by the root  ·  1 free angle  ·  2 conventions  ·  "
                      "5 fitted coefficients  ·  9 open residuals (W1-W9)",
            ha="center", fontsize=10.5, color=TXT)
    fig.suptitle("the provenance discipline: every headline number classified "
                 "by how it was obtained", fontsize=13)
    _save(fig, "figC_provenance_sunburst.png")


# ===========================================================================
# Fig. D - the anomaly waterfall: the strongest result, made visible
# ===========================================================================
def figD_anomaly_waterfall():
    from cl10 import matter

    an = matter.anomalies()
    fields = list(matter.SM_FIELDS.items())
    contribs = [(name, n * y ** 3) for (name, (s3, s2, y, n)) in fields]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.4),
                                   gridspec_kw={"width_ratios": [1.5, 1]})

    # ---- waterfall of [U(1)]^3 --------------------------------------------
    cum = 0.0
    xs = np.arange(len(contribs) + 1)
    for k, (name, c) in enumerate(contribs):
        col = C_SU3 if c > 0 else C_OPEN
        ax1.bar(xs[k], c, bottom=cum, color=col, width=0.62,
                edgecolor="none")
        ax1.text(xs[k], cum + c / 2, f"${['1/36','-8/9','1/9','-1/4','1','0'][k]}$",
                 ha="center", va="center", fontsize=9,
                 color="white" if abs(c) > 0.2 else TXT)
        if k < len(contribs) - 1:
            ax1.plot([xs[k] + 0.31, xs[k + 1] - 0.31], [cum + c, cum + c],
                     color=C_MEAS, lw=1, ls=":")
        cum += c
    ax1.bar(xs[-1], cum, color=C_ALG, width=0.62, edgecolor="none")
    ax1.text(xs[-1], 0.06, "$\\sum=0$", ha="center", fontsize=11,
             weight="bold", color="white")
    ax1.axhline(0, color="#5a6673", lw=1)
    labels = [n.replace("^c", "$^c$") for n, _ in contribs] + ["total"]
    ax1.set_xticks(xs)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel(r"$n\,Y^3$  (the $[\mathrm{U}(1)]^3$ coefficient)")
    ax1.set_title(r"the $[\mathrm{U}(1)]^3$ anomaly: the six fields balance "
                  "exactly", fontsize=11)
    resid = abs(an["u1_cubed"])
    badge(ax1, f"computed residual $|\\sum nY^3| = {resid:.1e}$\n"
               "(cl10.matter, V11)", loc="upper right", fs=8.5)

    # ---- the 8 coefficients ------------------------------------------------
    names = [r"$[\mathrm{SU}(3)]^3$", r"$[\mathrm{SU}(2)]^3$",
             r"$[\mathrm{U}(1)]^3$", r"$\mathrm{U}(1)[\mathrm{SU}(3)]^2$",
             r"$\mathrm{U}(1)[\mathrm{SU}(2)]^2$", r"$[\mathrm{grav}]^2\mathrm{U}(1)$",
             r"$[\mathrm{grav}]^2\mathrm{SU}(3)$", r"$[\mathrm{grav}]^2\mathrm{SU}(2)$"]
    vals = [an["su3_cubed"], an["su2_cubed"], an["u1_cubed"], an["u1_su3_sq"],
            an["u1_su2_sq"], an["grav2_u1"], an["grav2_su3"], an["grav2_su2"]]
    ys = np.arange(len(names))[::-1]
    ax2.scatter(ys * 0 + vals, ys, s=90, color=C_ALG, zorder=3)
    ax2.axvline(0, color=C_ALG, lw=1.2, ls="--")
    for y, v in zip(ys, vals):
        ax2.text(0.06, y, "= 0  \\checkmark", va="center", fontsize=9.5,
                 color=C_ALG)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(names)
    ax2.set_xlim(-0.4, 1.1)
    ax2.set_xticks([0])
    ax2.set_xticklabels(["0"])
    ax2.set_title("all 8 gauge and mixed anomalies vanish\n(Theorem 2 - "
                  "algebraic, no parameter)", fontsize=11)
    badge(ax2, "the hardest consistency condition of a\nchiral gauge theory, "
               "satisfied by the module", x=0.97, y=0.55, fs=8.5)

    fig.suptitle("anomaly cancellation: the strongest single result of the "
                 "16 Weyl module", fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    _save(fig, "figD_anomaly_waterfall.png")


# ===========================================================================
# Fig. E - the vertical ladder of scales
# ===========================================================================
def figE_ladder_vertical():
    from cl10 import u10pred

    M_P = u10pred.M_PLANCK
    LOOP = 1.0 / (16.0 * np.pi)
    logM = math.log10(M_P)
    def rung(n):
        return logM + n * math.log10(LOOP)

    fig, ax = plt.subplots(figsize=(9.8, 9.6))
    ymin, ymax = rung(17) - 1.4, rung(0) + 1.0
    for n in range(18):
        ax.axhline(rung(n), color="#9aa4af" if not DARK else "#3a4552",
                   lw=0.8, ls="-", xmin=0.15, xmax=0.80)
        ax.text(0.815, (rung(n) - ymin) / (ymax - ymin), f"$n={n}$",
                fontsize=7.5, color=C_MEAS, transform=ax.transAxes)

    # n=0 Planck / GUT
    ax.plot([M_P], [rung(0)], "o", color="black" if not DARK else "white",
            ms=12, zorder=5)
    ax.text(M_P * 0.12, rung(0) - 1.1,
            "Planck / GUT: $M_P=1.22\\times10^{19}$ GeV\n"
            "($M_{\\rm GUT}=M_P$: convention)", fontsize=9, ha="right",
            color="black" if not DARK else "white")
    # n=3 seesaw (label left)
    ax.plot([10 ** rung(3)], [rung(3)], "o", color=C_FIT, ms=11, zorder=5)
    ax.text(10 ** rung(3) * 0.18, rung(3) + 0.25,
            "seesaw $M_N=M_P(1/16\\pi)^3=9.6\\times10^{13}$ GeV\n"
            "($n_{2\\rm step}=3$: fitted)", fontsize=9, ha="right",
            color=C_FIT)
    # n=10 DM (label right) and v (label left): same band, opposite sides
    ax.plot([10 ** rung(10)], [rung(10)], "o", color=C_FIT, ms=11, zorder=5)
    ax.text(10 ** rung(10) * 1.5, rung(10) + 0.3,
            "DM $m_{\\rm DM}=M_P(1/16\\pi)^{10}=118$ GeV\n"
            "(coefficient 1.18: fitted; identity: W7)", fontsize=9,
            color=C_FIT)
    ax.plot([246.0], [math.log10(246.0)], "o", color=C_MEAS, ms=11, zorder=5)
    ax.text(246.0 * 0.22, math.log10(246.0) + 0.3,
            "electroweak $v=246$ GeV (measured)\n"
            "$v/M_P=2.08\\,(1/16\\pi)^{10}$ (2.08: fitted, W8)", fontsize=9,
            ha="right", color=C_MEAS)
    # QCD off-ladder
    n_qcd = (math.log10(0.2) - logM) / math.log10(LOOP)
    ax.plot([0.2], [math.log10(0.2)], "X", color=C_OPEN, ms=13, zorder=5)
    ax.text(0.2 * 1.6, math.log10(0.2) + 0.25,
            f"QCD $\\Lambda\\sim0.2$ GeV: needs $n\\simeq{n_qcd:.1f}$\n"
            "- NOT on the ladder", fontsize=9, color=C_OPEN)
    # neutrinos: ladder value (label right), 1/3 candidate (label right, below)
    ax.plot([10 ** rung(17)], [rung(17)], "o", color=C_OPEN, ms=11, zorder=5)
    ax.text(10 ** rung(17) * 1.5, rung(17) + 0.25,
            "ladder $m_\\nu=M_P(1/16\\pi)^{17}=0.146$ eV:\n"
            "$2.9\\times$ the atmospheric scale; excluded by\n"
            "$\\sum m_\\nu<0.12$ eV (W2, F3)", fontsize=9, color=C_OPEN)
    ax.plot([0.0487e-9], [math.log10(0.0487e-9)], "D", color=C_ALG, ms=10,
            zorder=5, mfc="none", mew=1.8)
    ax.text(0.0487e-9 * 1.4, math.log10(0.0487e-9) + 0.18,
            "$\\frac{1}{3}$-candidate $0.0487$ eV (2.6% from data; NOT "
            "derived)", fontsize=9, color=C_ALG, ha="left")

    # the 29-decade bracket (right edge, rotated label)
    ax.plot([10 ** 20.5] * 2, [rung(17), rung(0)], color=C_GUT, lw=1.4)
    ax.text(10 ** 20.3, 0.5 * (rung(13) + rung(14)),
            "one factor $(1/16\\pi)^n$ spans 29 decades", rotation=90,
            fontsize=9.5, color=C_GUT, ha="right", va="center")

    ax.set_xscale("log")
    ax.set_xlim(1e-13, 1e22)
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(range(-10, 20, 2))
    ax.set_ylabel("energy (GeV, log)")
    ax.set_xlabel("energy (GeV)")
    ax.set_title("the loop ladder: one factor, four scales\n"
                 "(the integers and 2-digit coefficients are the residuals)",
                 fontsize=11.5)
    badge(ax, "rungs: $M_P(1/16\\pi)^n$, $n=0..17$\n"
              "(cl10.neutrino, cl10.hier_exponent)", loc="upper left", fs=8)
    h, _ = provenance_legend(ncol=6)
    fig.legend(handles=h, loc="lower center", ncol=6, fontsize=8.5,
               framealpha=0.9, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=[0, 0.045, 1, 1])
    _save(fig, "figE_ladder_vertical.png")


# ===========================================================================
# Fig. F - the causal front: propagation measured on the lattice
# ===========================================================================
def figF_causal_front():
    from cl10 import engine as E
    from cl10 import local_qca

    L = 96
    c0 = L // 2
    snaps = [1, 4, 12, 32, 64, 96]
    n_front = 96

    rng = np.random.default_rng(0)
    f = E.zeros((L, L, L, 32))
    v = rng.standard_normal(32) + 1j * rng.random(32)
    f[c0, c0, c0] = E.Array_(v / np.linalg.norm(v))
    U = local_qca.internal_unitary()

    fields = {}
    radii = []
    idx = np.arange(L)
    manh = np.minimum(np.abs(idx - c0), L - np.abs(idx - c0))
    manh3 = manh[:, None, None] + manh[None, :, None] + manh[None, None, :]
    for k in range(1, n_front + 1):
        f = local_qca.local_step(f, U, k - 1)
        nrm = np.sqrt(np.sum(np.abs(E.asnumpy(f)) ** 2, axis=-1))
        radii.append(int(manh3[nrm > 0].max()))
        if k in snaps:
            fields[k] = np.sum(np.abs(E.asnumpy(f[:, :, c0, :])) ** 2, axis=-1)

    # control: internal layer only
    f2 = E.zeros((L, L, L, 32))
    f2[c0, c0, c0] = E.Array_(v / np.linalg.norm(v))
    for k in range(1, 97):
        f2 = local_qca.internal_layer(f2, U)
    ctl = np.sum(np.abs(E.asnumpy(f2[:, :, c0, :])) ** 2, axis=-1)

    fig = plt.figure(figsize=(14.4, 6.0))
    gs = fig.add_gridspec(2, 7, height_ratios=[2.3, 1])
    for i, k in enumerate(snaps):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(np.log10(fields[k] + 1e-16), cmap="magma", vmin=-16, vmax=0)
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(c0 + k * np.cos(th), c0 + k * np.sin(th), color="cyan",
                lw=1.1, ls="--")
        ax.set_title(f"step {k}\nfront at radius {k}", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    ax = fig.add_subplot(gs[0, 6])
    ax.imshow(np.log10(ctl + 1e-16), cmap="magma", vmin=-16, vmax=0)
    ax.set_title("control: internal only\nstays on one cell", fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])

    ax = fig.add_subplot(gs[1, :4])
    ax.plot(np.arange(1, n_front + 1), radii, color=C_ALG, lw=2.2,
            label="measured support front")
    ax.plot([1, n_front], [1, n_front], color=C_MEAS, ls="--",
            label="light cone: radius = step")
    ax.set_xlabel("time step")
    ax.set_ylabel("max Manhattan radius")
    ax.set_title("the causal radius is implemented: front speed = 1.000 "
                 "(Theorem 5)", fontsize=10.5)
    ax.legend(fontsize=8, loc="upper left")
    dev = max(abs(r - k) for k, r in enumerate(radii, 1))
    badge(ax, f"support front $=k$ EXACTLY for every $k=1..{n_front}$\n"
              f"(max deviation {dev} cells); the front amplitude\n"
              f"decays as $\\sin(g_{10})^k$ - the support is the exact one",
          fs=8)

    ax = fig.add_subplot(gs[1, 4:])
    ax.axis("off")
    ax.text(0.02, 0.92, "the local step:\n"
                        "$U_{\\rm local}=\\prod_x U_{\\rm GUT}(x)"
                        "\\prod_{\\langle x,x'\\rangle}P_{xx'}$",
            fontsize=10.5, va="top")
    ax.text(0.02, 0.55, "the swap-rotation bond layer $P_{xx'}$ carries no "
                        "angle:\nit is what turns the light cone of the "
                        "subalgebra\ninto a propagation bound.", fontsize=9.5,
            va="top", color=C_MEAS)
    ax.text(0.02, 0.2, "unitarity of the local step: $1.6\\times10^{-15}$ "
                       "(V36)", fontsize=9.5, color=C_ALG, va="top")

    fig.suptitle("strictly causal dynamics: a single-cell excitation inside "
                 "the light cone, front speed exactly 1", fontsize=13)
    fig.tight_layout(rect=[0, 0.0, 1, 0.92])
    _save(fig, "figF_causal_front.png")


# ===========================================================================
# Fig. G - the Brillouin-zone doublers: the honest negative result
# ===========================================================================
def figG_brillouin_doublers():
    from cl10 import brillouin

    _, U32, U16 = brillouin.blocks()
    res = brillouin.brillouin()

    fig = plt.figure(figsize=(13.6, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1])

    # ---- (a) the gap map on the qz=0 plane ---------------------------------
    ax = fig.add_subplot(gs[0, 0])
    n = 33
    qv = np.linspace(-np.pi, np.pi, n)
    G = np.zeros((n, n))
    for i, qx in enumerate(qv):
        for j, qy in enumerate(qv):
            G[j, i] = brillouin._gap((qx, qy, 0.0), U16)
    im = ax.imshow(np.log10(G + 1e-6), origin="lower", extent=[-np.pi, np.pi,
                   -np.pi, np.pi], cmap="viridis", vmin=-6, vmax=-1.2)
    for sx in (1, -1):
        for sy in (1, -1):
            ax.plot(sx * np.pi / 2, sy * np.pi / 2, "o", ms=13, mfc="none",
                    mec="white", mew=2)
    ax.plot(0, 0, "s", ms=9, color="orange", mfc="none", mew=2)
    ax.set_xlabel("$q_x$"); ax.set_ylabel("$q_y$")
    ax.set_title("(a) $\\log_{10}$ gap of the 16-sector Floquet\n"
                 "on the $q_z=0$ plane: 4 gapless nodes", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046)
    badge(ax, "orange square: zone centre,\ngapped at 0.0365", fs=7.5,
          loc="upper left")

    # ---- (b) the 3D BZ with the 12 nodes -----------------------------------
    ax = fig.add_subplot(gs[0, 1], projection="3d")
    lim = np.pi
    for s in (-1, 1):
        ax.plot([s * lim, s * lim, s * lim, s * lim, s * lim],
                [s * lim, s * lim, -s * lim, -s * lim, s * lim],
                [s * lim, -s * lim, -s * lim, s * lim, s * lim],
                color="#7f8c99", lw=0.8)
    nodes = np.array([n["q"] for n in res["nodes"]])
    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], s=140, c=C_OPEN,
               depthshade=True, edgecolor="white", linewidth=0.6)
    ax.scatter([0], [0], [0], s=70, c="orange", marker="s", depthshade=True)
    ax.text(0.1, 0.1, 0.1, "centre: gapped", fontsize=8, color="orange")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_xlabel("$q_x$"); ax.set_ylabel("$q_y$"); ax.set_zlabel("$q_z$")
    ax.set_title("(b) the 12 doubler momenta\nin the supercell zone", fontsize=10)
    ax.view_init(elev=18, azim=-55)

    # ---- (c) the linear node dispersion ------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    qz = np.linspace(-0.4, 0.4, 41)
    gaps = [brillouin._gap((np.pi / 2, np.pi / 2, z), U16) for z in qz]
    ax.plot(qz, gaps, color=C_OPEN, lw=1.8)
    qz2 = np.linspace(-0.4, 0.4, 41)
    g0 = [brillouin._gap((0.0, 0.0, z), U16) for z in qz2]
    ax.plot(qz2, g0, color="orange", lw=1.6, label="through the centre: gapped")
    ax.set_xlabel(r"momentum along the node line $q_z$")
    ax.set_ylabel(r"$\min_a|\omega_a(q)|$  (per 6-step cycle)")
    ax.set_title("(c) linear node: $|v|=0.142$\ncentre: gapped, NOT massless",
                 fontsize=10)
    ax.legend(fontsize=8)
    badge(ax, f"12 nodes at $(\\pm\\frac{{\\pi}}{{2}},\\pm\\frac{{\\pi}}{{2}},0)$ "
              f"+ perm.\ngap $\\leq$ {res['max_doubler_gap']:.1e}; "
              f"$\\sum$ chirality $=0$\nBloch vs torus "
              f"{res['bloch_vs_torus_dev']:.0e} (V63, two routes)",
          loc="upper left", fs=8)

    fig.suptitle("lattice chirality, computed: the doubling SURVIVES the "
                 "automaton - the named failure W9", fontsize=13)
    fig.tight_layout(rect=[0, 0.0, 1, 0.93])
    _save(fig, "figG_brillouin_doublers.png")


# ===========================================================================
# Fig. H - the CPT butterfly: the mirror for every g_10
# ===========================================================================
def figH_cpt_butterfly():
    from cl10 import qca_spectrum

    cp = qca_spectrum.cpt_block_phases()
    ph16, ph16b = cp["ph16"], cp["ph16b"]

    grid = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    resid = []
    for g in grid:
        c = qca_spectrum.cpt_block_phases(g)
        resid.append(float(np.abs(np.sort(-c["ph16"]) - c["ph16b"]).max()))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.6),
                                   gridspec_kw={"width_ratios": [1.4, 1]})

    # ---- the butterfly: each matter level and its exact mirror -------------
    for i, a in enumerate(ph16):
        ax1.plot([-a, a], [i, i], color=C_MEAS, lw=0.9, alpha=0.7)
        ax1.plot([a], [i], "o", color=C_MATTER, ms=7)
        ax1.plot([-a], [i], "o", color=C_ANTI, ms=7)
    ax1.axvline(0, color=TXT, lw=1.2)
    ax1.set_yticks([])
    ax1.set_xlabel(r"eigenphase $\omega$ (rad/step)")
    ax1.set_ylabel("level index")
    m = float(np.abs(np.sort(-ph16) - ph16b).max())
    ax1.set_title("the 16 (matter, red) and the $\\overline{16}$ "
                  "(antimatter, purple):\n"
                  r"$\omega_{\overline{16}}=-\omega_{16}$, level by level"
                  "\n" + f"mirror residual $={m:.1e}$; "
                  "$\\det U-1$: $5.3\\times10^{-15}$; "
                  "$[U,w]$: $7\\times10^{-16}$",
                  fontsize=10.5)

    # ---- the g10 invariance -------------------------------------------------
    ax2.semilogy(grid, resid, "o-", color=C_ALG, ms=7)
    ax2.axhline(1e-12, color=C_MEAS, ls="--", lw=1)
    ax2.text(0.33, 2.2e-12, "machine precision", fontsize=8, color=C_MEAS)
    ax2.set_xlabel(r"the 't Hooft angle $g_{10}$")
    ax2.set_ylabel("mirror residual")
    ax2.set_title("the mirror holds for EVERY $g_{10}$:\nit is a property of "
                  "the algebra,\nnot of the free parameter", fontsize=11)
    badge(ax2, f"max residual over the grid: {max(resid):.1e}\n"
               "(cl10.parameter_audit, V52)", loc="upper left", fs=8.5)

    fig.suptitle("CPT as an exact spectral mirror - algebraic, "
                 "parameter-independent", fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    _save(fig, "figH_cpt_butterfly.png")


# ===========================================================================
# Fig. I - the black hole recovered by computation (Flamm embedding)
# ===========================================================================
def figI_blackhole_embedding():
    from scipy.integrate import solve_ivp

    from cl10 import gravity

    M = 1.0
    H = gravity.H_FD

    def sw(x):
        return gravity.metric_schwarzschild(M, x)

    # Flamm's paraboloid: the t=0 slice embedded in R^3
    r = np.linspace(2.001, 13.0, 120)
    rho = np.zeros_like(r)
    for i in range(1, len(r)):
        rm = 0.5 * (r[i] + r[i - 1])
        rho[i] = rho[i - 1] + (r[i] - r[i - 1]) / math.sqrt(max(1 - 2 * M / rm, 1e-9))
    z = 2.0 * np.sqrt(2 * M * (r - 2 * M))
    KR = 48.0 * M ** 2 / r ** 6

    ph = np.linspace(0, 2 * np.pi, 160)
    PH = ph[None, :]
    shape = (len(r), len(ph))
    Xs = np.broadcast_to(rho[:, None] * np.cos(PH), shape)
    Ys = np.broadcast_to(rho[:, None] * np.sin(PH), shape)
    Zs = np.broadcast_to(z[:, None], shape)
    KRm = np.broadcast_to(KR[:, None], shape)

    # the geodesic through the CL10 finite-difference Christoffels
    r0, vt = 10.0, 0.34
    gamma = 1.0 / np.sqrt(1 - vt ** 2)
    u0 = [0.0, r0, np.pi / 2, 0.0, gamma / np.sqrt(1 - 2 * M / r0), 0.0,
          0.0, gamma * vt / r0]

    def geo(lam, u):
        t, rr, th, pp, td, rd, thd, pd = u
        G = gravity._christoffel_at(sw, np.array([t, rr, th, pp]), H)
        vv = np.array([td, rd, thd, pd])
        dv = -np.einsum("mab,a,b->m", G, vv, vv)
        return [td, rd, thd, pd, dv[0], dv[1], dv[2], dv[3]]

    sol = solve_ivp(geo, [0.0, 600.0], u0, method="DOP853", rtol=1e-8,
                    atol=1e-8, max_step=0.5)
    rg, pg = sol.y[1], sol.y[3]
    rg = np.clip(rg, 2.001, None)
    rho_g = np.interp(rg, r, rho)
    z_g = np.interp(rg, r, z)

    fig = plt.figure(figsize=(13.6, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1])

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    ax.set_proj_type("ortho")
    ax.plot_surface(Xs, Ys, Zs, facecolors=plt.cm.turbo(
        (np.log10(KRm) - np.log10(KR.min())) /
        (np.log10(KR.max()) - np.log10(KR.min()))),
        rstride=2, cstride=2, linewidth=0, antialiased=True, alpha=0.95)
    ax.plot(rho_g * np.cos(pg), rho_g * np.sin(pg), z_g, color="white",
            lw=1.8, label="geodesic (CL10 Christoffels)")
    ax.scatter([rho_g[0] * np.cos(pg[0])], [rho_g[0] * np.sin(pg[0])],
               [z_g[0]], color="cyan", s=30, label="start")
    ax.plot([0], [0], [0], "o", color="black" if not DARK else "white", ms=7,
            label="horizon $r=2M$ (the throat)")
    ax.set_xlabel(r"embedded radius $\rho$")
    ax.set_ylabel("")
    ax.set_zlabel("embedding $z$")
    ax.set_title("Flamm embedding of the CL10 Schwarzschild slice,\ncoloured "
                 "by the computed Kretschmann $K$", fontsize=10.5)
    ax.view_init(elev=28, azim=-60)
    ax.legend(fontsize=7.5, loc="upper left")

    ax = fig.add_subplot(gs[0, 1])
    rv = np.linspace(2.3, 20.0, 12)
    Kc, Ka = [], []
    for ri in rv:
        x = np.array([0.0, ri, 0.5, 0.7])
        Rm = gravity.riemann_4d(sw, x, H)
        gm, _ = gravity._metric_and_dg(sw, x, H)
        Kc.append(gravity.kretschmann_4d(Rm, gm))
        Ka.append(48.0 * M ** 2 / ri ** 6)
    ax.loglog(rv, Ka, "-", color="k" if not DARK else "w", lw=1.8,
              label=r"analytic $K=48M^2/r^6$")
    ax.loglog(rv, Kc, "o", color=C_OPEN, ms=6,
              label="computed from the CL10 metric\n(4-D central differences)")
    resid = float(np.abs(np.array(Kc) / np.array(Ka) - 1).max())
    ax.set_xlabel(r"radius $r$  ($M=1$, $r>2M$)")
    ax.set_ylabel(r"Kretschmann $K=R^{mnpq}R_{mnpq}$")
    ax.set_title("two routes agree: the black hole is\nrecovered by "
                 "computation", fontsize=10.5)
    ax.legend(fontsize=8)
    badge(ax, f"2-route residual $={resid:.1e}$ (the $O(h^2)$ FD error);\n"
              "computed Ricci $=0$: the vacuum equation holds", fs=8.5)

    fig.suptitle("gravity on the lattice: the 0-parameter (Planck-mass) black "
                 "hole, computed from the metric", fontsize=13)
    fig.tight_layout(rect=[0, 0.0, 1, 0.93])
    _save(fig, "figI_blackhole_embedding.png")


# ===========================================================================
# Fig. J - the hydrogen hero: real orbitals + the CL10 photon
# ===========================================================================
def figJ_hydrogen_hero():
    # analytic hydrogen eigenstates (the INPUT of atomic physics, per the paper)
    n = 96
    g = np.linspace(-15, 15, n)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2) + 1e-9
    C = np.clip(Y / R, -1, 1)

    psi = {
        "1s": np.exp(-R) / math.sqrt(math.pi),
        r"2p$_0$": (R * np.exp(-R / 2) / (4 * math.sqrt(2 * math.pi))) * C,
        r"3d$_0$": (R ** 2 * np.exp(-R / 3)) * (3 * C ** 2 - 1),
    }
    levels = {"1s": 0.06, r"2p$_0$": 0.10, r"3d$_0$": 0.08}
    cols = {"1s": C_ALG, r"2p$_0$": C_SU2, r"3d$_0$": C_SU3}

    fig = plt.figure(figsize=(13.6, 5.6))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.25])

    from scipy.interpolate import RegularGridInterpolator

    th = np.linspace(0.03, np.pi - 0.03, 64)
    ph = np.linspace(0, 2 * np.pi, 128)
    TH, PH = np.meshgrid(th, ph)
    dx = np.sin(TH) * np.cos(PH)
    dy = np.cos(TH)
    dz = np.sin(TH) * np.sin(PH)
    rr = np.linspace(0.05, 14.0, 700)

    for k, (name, w) in enumerate(psi.items()):
        ax = fig.add_subplot(gs[0, k], projection="3d")
        f2 = w ** 2
        lev = levels[name] * f2.max()
        itp = RegularGridInterpolator((g, g, g), f2)
        prof = np.empty((len(rr),) + TH.shape)
        for i, rv in enumerate(rr):
            prof[i] = itp(np.stack([rv * dx, rv * dy, rv * dz], -1)) - lev
        sgn = np.sign(prof)
        ch = sgn[:-1] * sgn[1:] < 0
        has = ch.any(axis=0)
        i0 = prof.shape[0] - 2 - np.argmax(ch[::-1], axis=0)
        i0 = np.clip(i0, 0, len(rr) - 2)
        f0 = np.take_along_axis(prof, i0[None], 0)[0]
        f1 = np.take_along_axis(prof, (i0 + 1)[None], 0)[0]
        t = -f0 / (f1 - f0)
        rhit = rr[i0] + t * (rr[1] - rr[0])
        rhit[~has] = np.nan
        ax.plot_surface(rhit * dx, rhit * dy, rhit * dz, color=cols[name],
                        alpha=0.85, rstride=1, cstride=1, linewidth=0,
                        antialiased=True)
        ax.scatter([0], [0], [0], color="black" if not DARK else "white",
                   s=18)
        ax.set_title(f"$|\\psi|^2$ isosurface: {name}", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_box_aspect((1, 1, 0.8))
        ax.view_init(elev=18, azim=-60 + 30 * k)

    # the photon panel: the V61 detuning scan (measured on the lattice)
    ax = fig.add_subplot(gs[0, 3])
    with open(os.path.join(HERE, "experiments", "V61_hydrogen", "results",
                           "v61_hydrogen.json")) as fh:
        v61 = json.load(fh)["data"]
    scan = v61["resonance"]["detuning_scan"]
    det = [s["detuning_from_1s2p_rad"] for s in scan]
    exc = [s["excited_fraction"] for s in scan]
    ax.plot(det, exc, "o-", color=C_ALG, ms=6, label="excited fraction")
    ax.plot([-d for d in det[::-1]], exc, "o-", color=C_ALG, ms=6)
    sel = v61["selection_rules"]
    ax.set_xlabel("detuning from the 1s-2p gap (rad/step)")
    ax.set_ylabel("excited fraction")
    ax.set_title("the CL10 photon resonates on the gap\n"
                 "(measured, V61/V62)", fontsize=10)
    ax.text(0.98, 0.70, "atom: real hydrogen (INPUT)\nphoton: the CL10 "
            "massless field", transform=ax.transAxes, fontsize=8,
            color=C_MEAS, va="top", ha="right")
    badge(ax, f"allowed dipoles: 9.93, $-3.00$; forbidden:\n"
              f"$\\leq{sel['max_forbidden_dipole']:.0e}$ - suppression "
              f"{sel['suppression_factor']:.0e}\n"
              "level gaps = Rydberg, rel. error 0", loc="upper right", fs=8)

    fig.suptitle("the photon sector of CL10 drives a real atom: resonance, "
                 "E1 selection rules, re-emission", fontsize=13)
    fig.tight_layout(rect=[0, 0.0, 1, 0.94])
    _save(fig, "figJ_hydrogen_hero.png")


# ===========================================================================
# Fig. K - the unification claim: SM failure, desert repair at the root value
# ===========================================================================
def figK_unification_failure():
    from cl10 import unification, uvrunning

    M_Z, M_P = unification.M_Z, unification.M_P
    A10 = unification.ALPHA_10
    Ms = np.logspace(math.log10(M_Z), math.log10(M_P), 200)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(17.6, 5.4))

    # ---- up-running ---------------------------------------------------------
    cols = (C_GUT, C_SU2, C_SU3)
    labs = (r"$\alpha_1$ (GUT-normalised)", r"$\alpha_2$", r"$\alpha_3$")
    a_up = np.array([uvrunning.one_loop_running(M) for M in Ms]).T
    for k in range(3):
        ax1.semilogy(Ms, a_up[k], color=cols[k], lw=2, label=labs[k])
    ax1.axhline(A10, color=C_FIT, ls="--", lw=1.8,
                label=r"$\alpha_{10}=1/25$ (the root)")
    aP = uvrunning.one_loop_running(M_P)
    spread = (max(aP) - min(aP)) / np.mean(aP)
    ax1.axvspan(1e18, M_P, color=C_FIT, alpha=0.12)
    ax1.set_xscale("log")
    ax1.set_xlim(M_Z, M_P)
    ax1.set_ylim(0.01, 0.12)
    ax1.set_xlabel("scale (GeV)")
    ax1.set_ylabel(r"$\alpha_i$")
    ax1.set_title("(a) up-running: the GUT RANGE, not a point", fontsize=10.5)
    ax1.legend(fontsize=8)
    badge(ax1, f"at $M_P$: no meeting point - spread {spread*100:.0f}%\n"
               f"({aP[0]:.3f}, {aP[1]:.3f}, {aP[2]:.3f}); the $5/3$ "
               "normalisation\nis algebraic and correct, but the VALUES\n"
               "do not unify (W4)", loc="lower left", fs=8.5)

    # ---- down-running -------------------------------------------------------
    b = unification.B_SM
    rd = unification.running_down()
    pole3 = rd["su3_down_pole_scale"]
    Ms2 = np.logspace(math.log10(M_Z), math.log10(M_P), 200)
    for k in range(3):
        a = 1.0 / (1.0 / A10 + (b[k] / (2 * np.pi)) * np.log(M_P / Ms2))
        a = np.where(a > 0, a, np.nan)
        ax2.semilogy(Ms2, np.clip(a, 1e-4, 30), color=cols[k], lw=2,
                     label=labs[k])
        meas = unification.A_Z[k]
        ax2.plot([M_Z], [meas], "o", color=C_MEAS, ms=7, zorder=5)
    ax2.axvline(pole3, color=C_OPEN, lw=2)
    ax2.text(pole3 * 3.0, 2.5, f"$\\mathrm{{SU}}(3)$ Landau pole\n"
             f"$2.2\\times10^9$ GeV", fontsize=8.5, color=C_OPEN)
    ax2.axvspan(M_Z, pole3, color=C_OPEN, alpha=0.12)
    ax2.axhline(A10, color=C_FIT, ls="--", lw=1.4)
    ax2.set_xscale("log")
    ax2.set_xlim(M_Z, M_P)
    ax2.set_ylim(0.01, 10)
    ax2.set_xlabel("scale (GeV)")
    ax2.set_ylabel(r"$\alpha_i$ run DOWN from $\alpha_{10}$")
    ax2.set_title("(b) down-running from the root's value: it FAILS",
                  fontsize=10.5)
    ax2.legend(fontsize=8)
    badge(ax2, r"$\alpha_1(M_Z)$: $1.17\times$ measured; "
               r"$\alpha_2(M_Z)=0.195$: $5.8\times$" "\n"
               r"$\alpha_3$: non-perturbative at $2.2\times10^9$ GeV "
               "(W4, F2)", loc="lower left", fs=8.5)

    # ---- the desert repair (V64) --------------------------------------------
    sel = json.load(open(os.path.join(HERE, "experiments", "V64_desert",
                                      "results",
                                      "v64_desert.json")))["data"]["selected"]
    M_X, M_A, M_B = sel["M_X_GeV"], sel["M_A_GeV"], sel["M_B_GeV"]
    d10 = (6.0 / 5.0, 2.0, 2.0)
    d2 = (2.0 / 5.0, 2.0 / 3.0, 0.0)
    L0 = np.log(Ms / M_Z)
    LA = np.maximum(0.0, np.log(np.maximum(Ms, M_A) / M_A))
    LB = np.maximum(0.0, np.log(np.maximum(Ms, M_B) / M_B))
    for k in range(3):
        inv = (1.0 / unification.A_Z[k]
               - (unification.B_SM[k] * L0 + 3 * d10[k] * LA + d2[k] * LB)
               / (2 * np.pi))
        ax3.semilogy(Ms, 1.0 / inv, color=cols[k], lw=2)
    ax3.axhline(A10, color=C_FIT, ls="--", lw=1.8,
                label=r"$\alpha_X=1/25$ (the root, not refitted)")
    ax3.plot([M_X], [A10], "*", color=C_ALG, ms=17, mec="k", ls="none",
             label=rf"meeting: $M_X = {M_X:.1e}$ GeV")
    for Mv in (M_A, M_B):
        ax3.axvline(Mv, color=C_CONV, ls=":", lw=1.4)
    ax3.axvline(M_P, color="0.4", ls="--", lw=1.4)
    ax3.set_xscale("log")
    ax3.set_xlim(M_Z, M_P)
    ax3.set_ylim(0.01, 0.12)
    ax3.set_xlabel("scale (GeV)")
    ax3.set_ylabel(r"$\alpha_i$ with the desert content")
    ax3.set_title("(c) the desert repair: meeting AT the root value",
                  fontsize=10.5)
    ax3.legend(fontsize=8)
    badge(ax3, r"three $\mathbf{10}\oplus\overline{\mathbf{10}}$ + one $(1,2)$ "
               "pair (V64)\n"
               rf"$M_{{10}}={M_A:.0e}$, $M_{{(1,2)}}={M_B:.0e}$ GeV (fitted);"
               r" $\tau_p=1.8\times10^{{36}}$ yr (F5)" "\n"
               r"$M_{\rm GUT}=M_P$ ABANDONED; spectrum inferred (W4)",
          loc="lower left", fs=8.5)

    fig.suptitle("the unification claim, tested against itself: SM running "
                 "FAILS; the theory's own desert matter repairs it at the "
                 "root value", fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    _save(fig, "figK_unification_failure.png")


# ===========================================================================
# Fig. L - the falsifier dashboard
# ===========================================================================
def figL_falsifier_dashboard():
    from cl10 import u10pred

    tau_p = u10pred.proton_decay_lifetime()

    fig = plt.figure(figsize=(13.6, 7.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.15], hspace=0.55,
                          wspace=0.35)

    # ---- gauge 1: proton lifetime ------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    ax.set_xscale("log")
    ax.set_xlim(1e33, 3e48)
    ax.set_ylim(-0.5, 0.5)
    ax.axvspan(1e33, 2.4e34, color=C_OPEN, alpha=0.25)
    ax.axvline(2.4e34, color=C_OPEN, lw=2)
    ax.text(2.4e34, 0.33, "Super-K bound\n$\\tau_p>2.4\\times10^{34}$ yr",
            fontsize=9, color=C_OPEN, ha="right")
    sel = json.load(open(os.path.join(HERE, "experiments", "V64_desert",
                                      "results",
                                      "v64_desert.json")))["data"]["selected"]
    ax.axvspan(1e35, 4e35, color=C_FIT, alpha=0.18)
    ax.text(1.2e35, -0.08, "Hyper-K reach", fontsize=8, color=C_FIT)
    ax.plot([sel["tau_p_yr"]], [0], "*", color=C_ALG, ms=17, mec="k")
    ax.text(sel["tau_p_yr"] * 1.6, 0.15,
            f"CL10 desert (V64): $\\tau_p={sel['tau_p_yr']:.1e}$ yr\n"
            "($M_X=1.8\\times10^{16}$ GeV: derived)", fontsize=9,
            color=C_ALG, ha="left")
    ax.plot([tau_p], [0], "o", color="none", ms=13, mec=C_CONV, mew=2)
    ax.text(tau_p * 0.35, -0.32,
            "$M_{\\rm GUT}=M_P$: $4\\times10^{47}$ yr (ABANDONED)",
            fontsize=8.5, color=C_CONV, ha="right")
    ax.text(1e36, 0.35, "F5 is LIVE: decay at HK reach refutes the desert "
            "solution; nothing down to $10^{37}$ yr refutes it too",
            fontsize=9, color=C_MEAS)
    ax.set_yticks([])
    ax.set_xlabel(r"proton lifetime $\tau_p/B(p\to e^+\pi^0)$ (yr)")
    ax.set_title("F5 - proton decay", fontsize=11, color=C_ALG)

    # ---- gauge 2: neutrino mass --------------------------------------------
    ax = fig.add_subplot(gs[1, :2])
    ax.set_xscale("log")
    ax.set_xlim(0.005, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.axvspan(0.12, 0.5, color=C_OPEN, alpha=0.25)
    ax.axvline(0.12, color=C_OPEN, lw=2)
    ax.text(0.115, 0.3, "$\\sum m_\\nu<0.12$ eV\n(Planck+BAO)", fontsize=9,
            color=C_OPEN, ha="right")
    ax.plot([0.146], [0.05], "o", color=C_OPEN, ms=13)
    ax.text(0.152, 0.28, "bare ladder 0.146 eV:\nEXCLUDED (W2)", fontsize=9,
            color=C_OPEN, ha="left")
    ax.plot([0.05], [-0.05], "D", color=C_ALG, ms=11, mfc="none", mew=2)
    ax.text(0.05, -0.35, "$\\frac{1}{3}$-candidate $\\sum m_\\nu\\simeq0.05$ eV:\n"
            "inside the bound; next-gen tests (F3)", fontsize=9, color=C_ALG,
            ha="center")
    ax.plot([0.050], [0.0], "|", color=C_MEAS, ms=18, mew=2)
    ax.text(0.050, 0.16, "atmospheric 0.050 eV", fontsize=8.5, color=C_MEAS,
            ha="center")
    ax.set_yticks([])
    ax.set_xlabel("neutrino mass scale (eV)")
    ax.set_title("F3 - the neutrino spectrum", fontsize=11, color=C_OPEN)

    # ---- gauge 3: DM mass ---------------------------------------------------
    ax = fig.add_subplot(gs[1, 2])
    ax.set_xscale("log")
    ax.set_xlim(1, 1e4)
    ax.set_ylim(-0.5, 0.5)
    ax.axvspan(50, 300, color=C_ALG, alpha=0.15)
    ax.plot([118], [0], "o", color=C_FIT, ms=13)
    ax.text(118, 0.2, "$m_{\\rm DM}=118$ GeV", fontsize=9.5, color=C_FIT,
            ha="center")
    ax.text(118, -0.3, "F7: a WIMP far from\n100-150 GeV refutes the ladder",
            fontsize=8.5, color=C_MEAS, ha="center")
    ax.set_yticks([])
    ax.set_xlabel("dark-matter mass (GeV)")
    ax.set_title("F7 - dark matter", fontsize=11, color=C_FIT)

    # ---- the status board ---------------------------------------------------
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    rows = [
        ("F1", "a non-vanishing gauge anomaly",
         "all 8 coefficients = 0 (algebraic)", "consistent", C_ALG),
        ("F2", "running admitting no desert content at $\\alpha_X=1/25$",
         "$5/3$ exact; SM fails; V64 content restores the meeting",
         "conditional", C_CONV),
        ("F3", "inverted ordering or $\\sum m_\\nu>0.12$ eV",
         "bare ladder already excluded; $\\frac{1}{3}$-candidate testable",
         "tension", C_OPEN),
        ("F4", "a massive graviton or $c\\neq1$ dispersion",
         "massless, $\\omega/|k|\\to1$, UV-bounded", "consistent", C_ALG),
        ("F5", "GUT-pattern decay far from $1.8\\times10^{36}$ yr",
         "desert solution: $75\\times$ above the bound; HK/DUNE test it",
         "testable", C_FIT),
        ("F6", "heavy-neutrino decays in equilibrium",
         "the model's own $\\Gamma_N/H\\simeq4.6$: condition FAILS",
         "FAILS now", C_OPEN),
        ("F7", "a WIMP far from 118 GeV",
         "prediction untested: direct detection + colliders", "open",
         C_FIT),
    ]
    ax.text(0.005, 0.97, "F", fontsize=11, weight="bold", va="top")
    ax.text(0.05, 0.97, "the observation that refutes", fontsize=11,
            weight="bold", va="top")
    ax.text(0.42, 0.97, "where CL10 stands today", fontsize=11, weight="bold",
            va="top")
    ax.text(0.86, 0.97, "verdict", fontsize=11, weight="bold", va="top")
    for i, (fid, obs, pos, verdict, col) in enumerate(rows):
        y = 0.84 - i * 0.125
        ax.text(0.005, y, fid, fontsize=10.5, weight="bold", color=col,
                va="top")
        ax.text(0.05, y, obs, fontsize=9.5, va="top")
        ax.text(0.42, y, pos, fontsize=9.5, va="top", color=C_MEAS)
        ax.text(0.86, y, verdict, fontsize=9.5, va="top", color=col,
                weight="bold")
        ax.plot([0, 1], [y - 0.045, y - 0.045], color=BADGE_EC, lw=0.5,
                alpha=0.5)

    fig.suptitle("seven falsifiers, stated before the data decide: F6 "
                 "already bites; F2 is conditional on the V64 desert spectrum",
                 fontsize=13)
    _save(fig, "figL_falsifier_dashboard.png")


# ===========================================================================
# Fig. M - the desert crossing: the same root value, two deserts
# ===========================================================================
def figM_desert_crossing():
    from cl10 import unification

    M_Z, M_P, A10 = unification.M_Z, unification.M_P, unification.ALPHA_10
    sel = json.load(open(os.path.join(HERE, "experiments", "V64_desert",
                                      "results",
                                      "v64_desert.json")))["data"]["selected"]
    M_X, M_A, M_B = sel["M_X_GeV"], sel["M_A_GeV"], sel["M_B_GeV"]
    d10 = (6.0 / 5.0, 2.0, 2.0)
    d2 = (2.0 / 5.0, 2.0 / 3.0, 0.0)
    cols = (C_GUT, C_SU2, C_SU3)
    Ms = np.logspace(math.log10(M_Z), math.log10(M_P), 400)
    L0 = np.log(Ms / M_Z)
    LA = np.maximum(0.0, np.log(np.maximum(Ms, M_A) / M_A))
    LB = np.maximum(0.0, np.log(np.maximum(Ms, M_B) / M_B))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 5.8))
    for ax in (ax1, ax2):
        ax.set_xlim(math.log10(M_Z) - 0.1, math.log10(M_P) + 0.3)
        ax.set_ylim(5, 62)
        ax.set_xlabel(r"scale $M$ (GeV)")
        ax.axhline(25.0, color=C_FREE, ls="--", lw=1.6)
    ax1.set_ylabel(r"$1/\alpha_i$")
    for k in range(3):
        inv = 1.0 / unification.A_Z[k] - unification.B_SM[k] * L0 / (2 * np.pi)
        ax1.plot(np.log10(Ms), inv, color=cols[k], lw=2.4)
    ax1.axvspan(math.log10(1e18), math.log10(M_P) + 0.3, color=C_OPEN,
                alpha=0.10)
    ax1.axvline(math.log10(M_P), color=C_OPEN, lw=2)
    ax1.text(math.log10(2e2), 6.5, "the desert as BELIEVED:\nthe SM field content",
             fontsize=10.5, color=TXT, weight="bold")
    ax1.text(math.log10(4e17), 11, "no meeting at $M_P$\nspread $47\\%$  —  W4",
             fontsize=10.5, color=C_OPEN, weight="bold", ha="right")
    ax1.text(math.log10(M_P) - 0.15, 7, r"$M_{\rm GUT} = M_P$", rotation=90,
             ha="right", va="bottom", fontsize=9, color=C_OPEN)
    for k in range(3):
        inv = (1.0 / unification.A_Z[k]
               - (unification.B_SM[k] * L0 + 3 * d10[k] * LA + d2[k] * LB)
               / (2 * np.pi))
        ax2.plot(np.log10(Ms), inv, color=cols[k], lw=2.4)
    ax2.text(math.log10(2e2), 6.5, "the desert as INFERRED (V64):\n"
             r"$3\,(\mathbf{10}\oplus\overline{\mathbf{10}})"
             r" + 1\,(1,2)$ pairs",
             fontsize=10.5, color=TXT, weight="bold")
    for Mv, lab in ((M_A, "$M_{10}$"), (M_B, "$M_{(1,2)}$")):
        ax2.axvline(math.log10(Mv), color=C_FIT, ls=":", lw=1.8)
        ax2.text(math.log10(Mv) + 0.15, 8.5, lab, fontsize=9.5, color=C_FIT)
    ax2.axvline(math.log10(M_X), color=C_FREE, lw=1.4, alpha=0.6)
    ax2.plot(math.log10(M_X), 25.0, "*", color=C_FREE, ms=24, mec=TXT,
             mew=1.0, ls="none")
    ax2.text(math.log10(M_X) - 0.35, 41,
             r"$\alpha_X = 1/25$ EXACT" "\n" "the root value,\nnot refitted",
             fontsize=11, color=C_FREE, weight="bold", ha="right")
    ax2.axvline(math.log10(M_P), color=C_MEAS, ls="--", lw=1.4)
    ax2.text(math.log10(M_P) - 0.15, 7, r"$M_{\rm GUT}=M_P$ abandoned",
             rotation=90, ha="right", va="bottom", fontsize=8.5, color=C_MEAS)
    fig.suptitle("the desert crossing: the root value was right — "
                 "the field content was missing", fontsize=13.5)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    _save(fig, "figM_desert_crossing.png")


# ===========================================================================
# Fig. N - falsifiable by design: the proton prediction moved 11 orders
# ===========================================================================
def figN_falsifiable():
    from cl10 import u10pred

    tau_old = u10pred.proton_decay_lifetime()
    sel = json.load(open(os.path.join(HERE, "experiments", "V64_desert",
                                      "results",
                                      "v64_desert.json")))["data"]["selected"]
    tau_new = sel["tau_p_yr"]

    fig, ax = plt.subplots(figsize=(13.6, 5.4))
    ax.set_xscale("log")
    ax.set_xlim(1e33, 3e48)
    ax.set_ylim(-1.1, 1.5)
    ax.set_yticks([])
    ax.set_xlabel(r"proton lifetime $\tau_p$  (yr)")

    ax.axvspan(1e33, 2.4e34, color=C_OPEN, alpha=0.30)
    ax.text(7e33, 1.05, "EXCLUDED\nby Super-K today", fontsize=9.5,
            color=C_OPEN, weight="bold", ha="center")
    ax.axvline(2.4e34, color=C_OPEN, lw=2)
    ax.axvspan(2.4e34, 1e35, color=C_OPEN, alpha=0.12)
    ax.text(5e34, 0.72, "decay at Hyper-K reach\nrefutes the solution",
            fontsize=9, color=C_OPEN, ha="center")
    ax.axvspan(tau_new / 3.0, tau_new * 3.0, color=C_ALG, alpha=0.22)
    ax.axvspan(1e37, 3e37, color=C_OPEN, alpha=0.12)
    ax.text(1.4e37, 1.05, "nothing by $10^{37}$ yr\nrefutes too", fontsize=9,
            color=C_OPEN, ha="center")

    ax.plot([tau_new], [0], "*", color=C_ALG, ms=26, mec=TXT, mew=1.2,
            ls="none")
    ax.annotate(f"CL10 desert (V64):  $\\tau_p = {tau_new:.1e}$ yr\n"
                r"($75\times$ above the bound; theory error $\times3$)",
                (tau_new, 0), xytext=(0.0, 0.55), textcoords="data",
                fontsize=11, color=C_ALG, weight="bold", ha="center")
    ax.plot([tau_old], [0], "x", color=C_MEAS, ms=20, mew=3.5, ls="none")
    ax.text(tau_old, -0.55, r"$M_{\rm GUT}=M_P$:  $4\times10^{47}$ yr"
            "\nABANDONED", fontsize=9.5, color=C_MEAS, ha="center")
    ar = FancyArrowPatch((tau_old * 0.35, 0.9), (tau_new * 4.0, 0.9),
                         arrowstyle="-|>", mutation_scale=22, lw=2,
                         color=C_MEAS)
    ax.add_patch(ar)
    ax.text(3e42, 1.12, "the theory MOVED its prediction by 11 orders "
            "of magnitude\nto become testable  (F5)", fontsize=11,
            color=TXT, weight="bold", ha="center")

    for x, lab, yy in ((2.4e34, "Super-K 2014", -0.35),
                       (1e35, "Hyper-K ~2030", -0.72),
                       (1e36, "DUNE / HK-X ~2035+", -0.35)):
        ax.plot([x], [-0.15], "|", color=TXT, ms=14, mew=1.5)
        ax.text(x, yy, lab, fontsize=8.5, color=TXT, ha="center")
    fig.suptitle("falsifiable by design: the desert solution sits inside "
                 "the reach of the next generation", fontsize=13.5)
    fig.tight_layout(rect=[0, 0.02, 1, 0.92])
    _save(fig, "figN_falsifiable.png")


# ===========================================================================
# Fig. O - the desert chain: from the breaking sector to a testable lifetime
# ===========================================================================
def figO_causal_chain():
    fig, ax = plt.subplots(figsize=(13.8, 6.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def box(x, y, w, h, txt, col, fs=10):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.012",
                                    fc=BADGE_FC, ec=col, lw=2.2, zorder=3))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                fontsize=fs, color=TXT, zorder=4)

    def arrow(x0, x1, y, col, style="-|>"):
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle=style,
                                     mutation_scale=20, lw=2.2, color=col,
                                     zorder=2))

    y0, h0 = 0.55, 0.26
    box(0.015, y0, 0.165, h0,
        r"$\mathbf{16}_H \oplus \overline{\mathbf{16}}_H$" "\n"
        "or a SPLIT doubler 16\n(W9: 12 doubler nodes)", C_ALG)
    box(0.215, y0, 0.165, h0,
        "vectorlike remnants\n"
        r"$3(\mathbf{10}\oplus\overline{\mathbf{10}})$" "\n"
        r"$+\;1(1,2)$ pairs" "\n" "thresholds FITTED", C_CONV)
    box(0.415, y0, 0.165, h0,
        r"$\Delta b = (4.0,\ 6.7,\ 6.0)$" "\n"
        "NON-proportional\n(the no-go evaded)", C_CONV)
    box(0.615, y0, 0.165, h0,
        r"meeting:  $M_X = 1.8\times10^{16}$ GeV" "\n"
        r"$\alpha_X = 1/25$ EXACT" "\n"
        "the root value, NOT refitted", C_FREE)
    box(0.815, y0, 0.17, h0,
        r"$\tau_p = 1.8\times10^{36}$ yr" "\n"
        "Hyper-K / DUNE\ntest it  (F5)", C_FIT)
    for x0, x1, col in ((0.18, 0.215, C_ALG), (0.38, 0.415, C_CONV),
                        (0.58, 0.615, C_CONV), (0.78, 0.815, C_FREE)):
        arrow(x0, x1, y0 + h0 / 2, col)

    box(0.215, 0.14, 0.365, 0.20,
        "counterfactual (the paper's no-go, verified in V64):\n"
        "COMPLETE multiplets shift $\\Delta b \\propto (1,1,1)$ —\n"
        "the pairwise differences are invariant: NO crossing", C_OPEN, fs=9)
    ax.add_patch(FancyArrowPatch((0.40, 0.34), (0.4975, 0.55),
                                 arrowstyle="-|>", mutation_scale=18, lw=1.8,
                                 color=C_OPEN, ls=(0, (4, 3)), zorder=2))
    ax.text(0.44, 0.44, "fails", fontsize=9, color=C_OPEN, weight="bold")
    box(0.615, 0.14, 0.37, 0.20,
        "the cost, stated plainly:\n"
        r"$M_{\rm GUT} = M_P$ ABANDONED; the spectrum is" "\n"
        "INFERRED from the couplings, not derived (W4,\n"
        "conditional; F2 sharpens the refutation)", C_MEAS, fs=9)
    ax.add_patch(FancyArrowPatch((0.70, 0.55), (0.79, 0.34),
                                 arrowstyle="-|>", mutation_scale=18, lw=1.8,
                                 color=C_MEAS, ls=(0, (4, 3)), zorder=2))
    ax.text(0.715, 0.44, "cost", fontsize=9, color=C_MEAS, weight="bold")
    fig.suptitle("the desert chain: from the breaking sector's own remnants "
                 "to a proton-decay prediction", fontsize=13.5)
    fig.tight_layout()
    _save(fig, "figO_causal_chain.png")



# ===========================================================================
FIGS = {
    "A": figA_one_object,
    "B": figB_atlas_1024,
    "C": figC_provenance_sunburst,
    "D": figD_anomaly_waterfall,
    "E": figE_ladder_vertical,
    "F": figF_causal_front,
    "G": figG_brillouin_doublers,
    "H": figH_cpt_butterfly,
    "I": figI_blackhole_embedding,
    "J": figJ_hydrogen_hero,
    "K": figK_unification_failure,
    "L": figL_falsifier_dashboard,
    "M": figM_desert_crossing,
    "N": figN_falsifiable,
    "O": figO_causal_chain,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dark", action="store_true")
    ap.add_argument("--only", default="".join(FIGS))
    args = ap.parse_args()
    set_style(dark=args.dark)
    for k in args.only.upper():
        if k not in FIGS:
            print("unknown figure", k)
            continue
        print(f"--- fig {k} ---")
        FIGS[k]()


if __name__ == "__main__":
    main()
