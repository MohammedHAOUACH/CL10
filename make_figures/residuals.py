"""Figure module: residuals."""
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from .style import ROOT, save
from .style import C_ALG, C_CONV, C_FIT, C_GUT, C_OPEN, C_ROOT, C_SU2, C_SU3, C_ZERO


def _treemap(items, x, y, w, h, out):
    """Binary treemap: items [(idx, size)] laid out in the rect (x, y, w, h)."""
    if not items:
        return
    if len(items) == 1:
        out.append((items[0][0], x, y, w, h))
        return
    tot = sum(s for _, s in items)
    best, best_ratio, acc = 1, None, 0
    for k in range(1, len(items)):
        acc += items[k - 1][1]
        r = max(acc / (tot - acc), (tot - acc) / acc)
        if best_ratio is None or r < best_ratio:
            best_ratio, best = r, k
    frac = sum(s for _, s in items[:best]) / tot
    if w >= h:
        _treemap(items[:best], x, y, w * frac, h, out)
        _treemap(items[best:], x + w * frac, y, w * (1 - frac), h, out)
    else:
        _treemap(items[:best], x, y, w, h * frac, out)
        _treemap(items[best:], x, y + h * frac, w, h * (1 - frac), out)


def fig8_audit():
    """The parameter audit as a treemap: every headline number of Table 5 is
    one equal-area cell; cells are grouped and coloured by provenance."""
    classes = [  # (label, color, the numbers themselves)
        ("algebraic — 0-parameter", C_ALG,
         [r"$16/\overline{16}$", "the $45$", "the $5/3$", "8 anomalies",
          "32-spectrum", r"$\det U = 1$", "CPT mirror"]),
        ("free parameter", C_ROOT, [r"$g_{10} = 1/5$"]),
        ("convention", C_CONV,
         [r"$\alpha_{10} \equiv g_{10}^2$", r"$g_0 = 1/2^{n_{\rm space}}$"]),
        ("fitted coefficient", C_FIT,
         ["2.08 (hierarchy)", "1.18 (DM mass)", r"$n_{2\text{step}} = 3$",
          r"$M_{10}$", r"$M_{(1,2)}$", r"$\lambda_* = 4{\times}10^{-4}$",
          r"$\kappa_1$", r"$\kappa_2$", r"$\kappa_3$", r"$\kappa_4$",
          r"$\kappa_5$", r"$\kappa_6$"]),
        ("open residual", C_OPEN,
         ["W1 $\\Lambda$", "W2 $m_\\nu$", "W3 $N_{\\rm gen}$", "W4 unification",
          "W5 $\\Delta B$", "W6 masses", "W7 DM", "W8 hierarchy", "W9 quantum"]),
    ]
    items = [(ci, len(nums)) for ci, (_, _, nums) in enumerate(classes)]
    rects = []
    _treemap(items, 0, 0, 12.0, 5.6, rects)

    fig, ax = plt.subplots(figsize=(12.2, 6.1))
    ax.set_xlim(0, 12.0)
    ax.set_ylim(0, 5.6)
    ax.axis("off")
    ax.set_aspect("equal")

    for (ci, x, y, w, h) in rects:
        lab, col, nums = classes[ci]
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=col, ec="white",
                                   lw=2.0, alpha=0.92, zorder=2))
        hdr = 0.30
        fs_h = 12 if w >= 2.2 else 8
        ax.text(x + 0.08, y + h - hdr / 2 - 0.03, f"{lab}  ({len(nums)})",
                fontsize=fs_h, weight="bold", color="white",
                ha="left", va="center", zorder=4)
        # equal-area cells inside the block: pick the grid that wastes the
        # fewest slots and keeps cells closest to square
        n = len(nums)
        iw, ih = w - 0.16, h - hdr - 0.16
        cols = min(range(1, n + 1),
                   key=lambda c: abs(math.log((iw / c) / (ih / math.ceil(n / c))))
                                 + 0.15 * (math.ceil(n / c) * c - n))
        rows = math.ceil(n / cols)
        cw = (w - 0.16) / cols
        ch = (h - hdr - 0.16) / rows
        for i, name in enumerate(nums):
            r_, c_ = divmod(i, cols)
            cx = x + 0.08 + c_ * cw
            cy = y + h - hdr - 0.08 - (r_ + 1) * ch
            ax.add_patch(plt.Rectangle((cx + 0.02, cy + 0.02),
                                       cw - 0.04, ch - 0.04,
                                       fc="white", ec=col, lw=1.0,
                                       alpha=0.96, zorder=3))
            fs = 11 if min(cw, ch) > 0.7 else 9.5
            ax.text(cx + cw / 2, cy + ch / 2, name, fontsize=fs,
                    ha="center", va="center", color="0.15", zorder=4)
    ax.set_title("The parameter audit, by area: 31 headline numbers, five provenance classes\n"
                 "(Table 5 — the algebraic block is 0-parameter and $g_{10}$-independent; "
                 "everything else is named)", fontsize=13.5)
    fig.tight_layout(pad=0.5)
    save(fig, "fig8_parameter_audit.png")

def fig23_gauged_layer():
    import json
    with open(os.path.join(ROOT, "experiments", "V79_gauged_bond_layer",
                           "results", "v79_gauged_bond_layer.json")) as f:
        d = json.load(f)["data"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    rows = [
        (r"link algebra $[L_A,L_B] = if_{AB}^{C}L_C$",
         d["link_algebra"]["L_commutator_err"]),
        ("Jacobi identity", d["link_algebra"]["jacobi_err"]),
        (r"adjoint Casimir $= 64$", d["link_algebra"]["adjoint_casimir_err"]),
        (r"spinor Casimir $= 45$", d["spinor_casimir"]["err"]),
        (r"Gauss algebra $[Q_A,Q_B] = if_{AB}^{C}Q_C$",
         d["gauss_algebra"]["err"]),
        ("step unitarity", d["K"]["unitarity_probe_err"]),
        ("bond charge covariance", d["K_gauss"]["Ad_K_bond_charge_err"]),
    ]
    ypos = np.arange(len(rows))[::-1]
    errs = [max(r[1], 2e-17) for r in rows]
    axA.barh(ypos, errs, color=C_SU3, alpha=0.85, height=0.6)
    axA.set_xscale("log")
    axA.set_xlim(1e-17, 1e-12)
    axA.axvline(1e-15, color="0.35", lw=1.4, ls="--",
                label="machine precision $10^{-15}$")
    axA.set_yticks(ypos)
    axA.set_yticklabels([r[0] for r in rows], fontsize=10.5)
    axA.set_xlabel(r"absolute error (exact zeros floored at $2\times10^{-17}$)")
    axA.set_title("(a)  the exactness ledger: every identity of the\n"
                  "     gauged bond layer closes at machine precision",
                  fontsize=12)
    axA.legend(fontsize=10, loc="upper right", framealpha=0.95)
    leak = d["K_gauss"]["single_site_Gauss_leak"]
    axA.text(0.44, 0.46,
             "control: a single-site (non-local)\n"
             "step leaks the Gauss charge at "
             f"{leak:.2f}\n- off-scale by 15 orders",
             transform=axA.transAxes, fontsize=10.5, color=C_ZERO)

    lv = d["content"]["casimir_levels"]
    mu = d["content"]["multiplicities"]
    axB.bar([str(int(v)) for v in lv], mu, color=C_ROOT, alpha=0.85,
            width=0.55)
    axB.set_yscale("log")
    axB.set_ylim(1, 3000)
    for x, m in zip(range(len(lv)), mu):
        axB.text(x, m * 1.15, str(m), ha="center", fontsize=11.5)
    axB.set_xlabel(r"quadratic Casimir $C_2$ of the invariant sector")
    axB.set_ylabel("multiplicity")
    axB.set_title("(b)  gauge-invariant content: $32 \\oplus 288 \\oplus 1120$"
                  f"\n     = {sum(mu)} states; no singlet (triality)",
                  fontsize=12)
    fig.suptitle("The gauged bond layer (V79): the Gauss constraint closes "
                 "on the root module", fontsize=13.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig23_gauged_bond_layer.png")

def fig27_w4_closure():
    import json
    def load(e, f):
        with open(os.path.join(ROOT, "experiments", e, "results", f)) as fh:
            return json.load(fh)["data"]
    d74 = load("V74_multiscale_nodes", "v74_multiscale_nodes.json")
    d81 = load("V81_link_desert_content", "v81_link_desert_content.json")
    d82 = load("V82_composite_desert_exclusion",
               "v82_composite_desert_exclusion.json")
    d86 = load("V86_junction_desert", "v86_junction_desert.json")
    d87 = load("V87_general_junction", "v87_general_junction.json")
    mx74 = d74["M_X_node_GeV"]
    res82 = d82["A2_full_weight"]["min_max_residual"]
    n81 = d81["meeting"]["n_solutions"]
    b86 = d86["best"]
    b87 = d87["general_best"]

    fig, ax = plt.subplots(figsize=(12.8, 6.6))
    ax.set_xlim(0, 12.8)
    ax.set_ylim(0, 7.4)
    ax.axis("off")

    def box(x, y, w, h, txt, ec, fc, lw=1.6, fs=12):
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=True, fc=fc,
                                   ec=ec, lw=lw, zorder=2))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                fontsize=fs, zorder=3)

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.4),
                    zorder=1)

    box(3.4, 6.5, 6.0, 0.7,
        "W4: which content makes the measured couplings meet\n"
        "at $\\alpha_X = 1/25$ with desert thresholds?  (forced: V73)",
        C_ROOT, "#f3ecf7", lw=2.0, fs=11)

    left = [
        ("three full $\\mathbf{16}$ nodes,\nany strengths  (V74)",
         f"rank-1 Jacobian: meeting FORCED\n"
         f"$M_X = {mx74:.1e}$ GeV - excluded by F5"),
        ("gauged bond-layer composites\n$16 \\otimes 45$  (V81)",
         f"{n81} solutions over the full scan\nof the breaking direction"),
        ("mass-ansatz families,\nrelaxed $\\alpha_X$  (V82)",
         f"min residual $\\geq {res82:.2f}$ (target 0);\n0 physical solutions"),
    ]
    right = [
        ("any COMPLETE multiplet content\n(V84, theorem)",
         "equal shift $\\sum db_1 = \\sum db_2$\nover every complete multiplet"),
        ("junction-split nodes,\ndemocratic  (V86)",
         f"pattern FABRICATED (split 55x);\nresidual {b86['residual']:.2f}, "
         f"$\\alpha_X$ {b86['implied_alpha_X']:.3f}"),
        ("junction-split nodes,\ngeneral ratios  (V87)",
         f"floor UNCHANGED at {b87['residual']:.2f};\n"
         f"$\\alpha_X$ {b87['implied_alpha_X']:.3f} (data 0.04)"),
    ]
    for k, (q, r) in enumerate(left):
        y = 4.85 - k * 1.55
        box(0.5, y + 0.85, 5.4, 0.62, q, "0.25", "#eef3f8")
        box(0.5, y, 5.4, 0.72, r, C_ZERO, "#fbeeee")
        ax.text(6.15, y + 0.36, "$\\mathbf{\\times}$", fontsize=19,
                color=C_ZERO, ha="center", va="center", weight="bold")
        arrow(3.8, 6.5, 3.2, y + 1.47)
    for k, (q, r) in enumerate(right):
        y = 4.85 - k * 1.55
        box(7.0, y + 0.85, 5.4, 0.62, q, "0.25", "#eef3f8")
        box(7.0, y, 5.4, 0.72, r, C_ZERO, "#fbeeee")
        ax.text(6.75, y + 0.36, "$\\mathbf{\\times}$", fontsize=19,
                color=C_ZERO, ha="center", va="center", weight="bold")
        arrow(9.0, 6.5, 9.7, y + 1.47)

    box(1.9, 0.30, 9.0, 0.95,
        "CLOSED CONDITIONALLY BY V96: a flavour-dependent bond layer "
        "($\\kappa_s$ per class)\nmeets exactly at $M_X = 1.8\\times10^{16}$ "
        "GeV for all $n = 3\\dots12$ - six named inputs (Fig. 33)",
        "#2c7a3f", "#eaf6ec", lw=2.0, fs=11.5)
    for x in (3.2, 9.7):
        arrow(x, 1.75, 6.4, 1.2)

    ax.set_title("W4: every exit of the BASE automaton closed - the "
                 "fragment content is feasible (V92) and the meeting "
                 "closes conditionally (V96)", fontsize=13.5)
    fig.tight_layout()
    save(fig, "fig27_w4_closure.png")

def fig28_w7_portal():
    import json
    with open(os.path.join(ROOT, "experiments", "V88_mirror_scalar",
                           "results", "v88_mirror_scalar.json")) as f:
        d = json.load(f)["data"]
    lam2, sig2 = d["V2_lam_star"], d["V2_sigma_SI_cm2"]
    lam1, sig1 = d["V1_lam_for_Omega"], d["V1_sigma_SI_cm2"]
    lz = 1.5e-47
    gap = d["ss_hh_gap_GeV"]
    cloud = d["V1_virtual_Pmax"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.25, 1.0]})

    lam = np.logspace(-5, 0, 400)
    f_N, m_N, M_H, M_DM = 0.30e-3, 0.939, 125.1, 118.0
    mu = M_DM * m_N / (M_DM + m_N)
    sig = lam ** 2 * f_N ** 2 * mu ** 2 / (np.pi * M_H ** 4) * 0.3894e-27
    axA.loglog(lam, sig, color=C_ROOT, lw=2.2,
               label=r"$\sigma_{\rm SI} = \lambda^2 f_N^2 \mu^2 / "
                     r"(\pi m_h^4)$")
    axA.axhline(lz, color=C_ZERO, ls="--", lw=1.8,
                label=r"LZ bound (4.2 t yr): $2.2\times10^{-48}$ min "
                      r"(40 GeV), $\approx1.5\times10^{-47}$ at 118 GeV")
    axA.fill_between([1e-5, 1], 1e-58, lz, color=C_ZERO, alpha=0.07)
    axA.plot([lam2], [sig2], "o", ms=11, color=C_SU2, zorder=5)
    axA.annotate(f"pair portal (Z2 BROKEN)\n$\\lambda_* = {lam2:.1e}$: "
                 f"$\\Omega h^2 = 0.120$\n$\\sigma_{{\\rm SI}} = "
                 f"{sig2:.0e}$ - below reach",
                 xy=(lam2, sig2), xytext=(3e-4, 8e-43), fontsize=11,
                 color=C_SU2, va="top",
                 arrowprops=dict(arrowstyle="->", color=C_SU2, lw=1.2))
    axA.plot([lam1], [sig1], "X", ms=13, color=C_ZERO, zorder=5)
    axA.annotate(f"Higgs portal (Z2 EXACT)\n$\\lambda = {lam1:.2f}$: "
                 f"$\\Omega h^2 = 0.120$\n$\\sigma_{{\\rm SI}} = "
                 f"{sig1:.0e}$ - EXCLUDED",
                 xy=(lam1, sig1), xytext=(0.1, 2e-53), fontsize=11,
                 color=C_ZERO,
                 arrowprops=dict(arrowstyle="->", color=C_ZERO, lw=1.2))
    axA.set_xlabel(r"portal coupling $\lambda$")
    axA.set_ylabel(r"direct detection $\sigma_{\rm SI}$ [cm$^2$]")
    axA.set_title("(a)  both portals reproduce the abundance - one is "
                  "excluded,\n     the other breaks the exact mirror "
                  "$\\mathbb{Z}_2$", fontsize=12.5)
    axA.set_ylim(1e-58, 1e-42)
    axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24),
               fontsize=9.5, framealpha=0.95)

    rows = [
        ("carrier stability", "EXACT (both portals)", "#2c7a3f"),
        ("stabilizer", "scalar parity $(-1)^{N_s}$", "#2c7a3f"),
        ("$s\\overline{s} \\to hh$", f"CLOSED (gap {gap:.1f} GeV;\n"
         f"virtual cloud {cloud:.0e})", C_SU2),
        ("$h^* \\to f\\overline{f}$", "FORBIDDEN if Z2 exact\n"
         "(pairs are mirror-odd)", C_ZERO),
        ("open channel", "$s\\overline{s} \\to hh^*$ only (Z2 exact)",
         C_SU2),
        ("mirror relic $\\overline{Q}$", "stable under both (V66 intact)",
         "#2c7a3f"),
    ]
    axB.axis("off")
    axB.set_title("(b)  the audit verdicts (V88)", fontsize=12.5)
    y = 0.98
    for (k, v, c) in rows:
        axB.text(0.02, y, k, fontsize=10.5, weight="bold", va="top")
        axB.text(0.52, y, v, fontsize=10.5, color=c, va="top")
        y -= (max(k.count("\n"), v.count("\n")) + 1.35) * 0.068
    axB.text(0.02, 0.01,
             "W7: the carrier EXISTS and is exactly stable; the observed\n"
             "abundance and the exact mirror $\\mathbb{Z}_2$ cannot both "
             "hold (F7).",
             fontsize=10.5, color=C_ROOT, weight="bold", va="bottom")

    fig.suptitle("The built dark-matter fix, audited against its own "
                 "symmetry (V88)", fontsize=13.5)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.34)
    save(fig, "fig28_w7_portal.png")

def fig29_loop_selection():
    import json
    with open(os.path.join(ROOT, "experiments", "V89_loop_selection",
                           "results", "v89_loop_selection.json")) as f:
        d = json.load(f)["data"]
    root = d["root"]
    Fn = root["F_n"]
    ns = np.arange(2, 31)
    loop = d["loop_law"]
    nstar = root["n_star_below_loop_factor"]
    fit = root["geometric_fit"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.3, 1.0]})

    axA.semilogy(ns, Fn, "o-", color=C_GUT, ms=5, lw=1.6,
                 label=r"Cartan ratio $F_n$ (root coupling)")
    nn = np.arange(10, 171)
    axA.semilogy(nn, fit["A"] * fit["r"] ** nn, "--", color=C_ROOT, lw=1.4,
                 label=fr"geometric law $A r^n$, $r = {fit['r']:.3f}$")
    axA.axhline(loop, color=C_ZERO, ls="--", lw=1.8,
                label=r"loop factor $1/(16\pi) = 0.0199$")
    axA.axvline(10, color=C_ZERO, lw=1.2, alpha=0.6)
    axA.plot([10], [Fn[8]], "*", ms=17, color=C_ZERO, zorder=6)
    axA.plot([11], [Fn[9]], "o", ms=11, mfc="none", mec=C_SU2, mew=2.0,
             ls="none", zorder=6)
    axA.plot([170], [root["F_170_extrapolated"]], "D", ms=9, color=C_ROOT,
             ls="none", zorder=6)
    axA.set_xlabel("closed-loop length $n$")
    axA.set_ylabel(r"Cartan-content ratio $F_n$")
    axA.set_title("(a)  the ladder length 10 is the CROSSING of the\n"
                  "     closed-loop suppression through $1/(16\\pi)$",
                  fontsize=12.5)
    axA.legend(loc="lower left", fontsize=9.5, framealpha=0.95)

    axB.axis("off")
    axB.set_title("(b)  the selection rule, measured to $n = 30$ (V89)",
                  fontsize=12.5)
    rows = [
        ("crossing $n_*$", f"EXACTLY {nstar} (ladder length)", "#2c7a3f"),
        ("spectral radius", f"$n_\\rho = {root['n_rho_crossing']}$ "
         "(different object)", C_SU2),
        ("mod-4 rule", "ABSENT (min per class "
         f"{min(v['min'] for v in root['mod4_pattern_n8_30'].values()):.0e}"
         " in every class)", C_ZERO),
        ("geometric ratio", f"$r = {fit['r']:.3f}$ (fit $n = 10..30$)",
         C_GUT),
        ("17-factor comp.", f"misses by x"
         f"{root['extrapolation_ratio']:.0e} - named residual", C_ZERO),
    ]
    for i, (k, v, c) in enumerate(rows):
        y = 0.95 - i * 0.115
        axB.text(0.02, y, k, fontsize=11, weight="bold", va="top")
        axB.text(0.50, y, v, fontsize=11, color=c, va="top")
    axB.text(0.02, 0.01,
             "W8: the LENGTH-10 factor is selected by the crossing;\n"
             "the COMPOSITION of 17 factors is not supplied by one\n"
             "geometric law - the residual is now precisely named.",
             fontsize=10, color=C_ROOT, weight="bold", va="bottom")

    fig.text(0.03, 0.015,
             f"star: crossing at EXACTLY $n = {nstar}$ "
             f"($F_{{10}} = {Fn[8]:.4f}$ vs 0.0199);   "
             f"open circle: $F_{{11}} = {Fn[9]:.2f}$ spike - no mod-4 rule;\n"
             f"diamond: $F_{{170}} = {root['F_170_extrapolated']:.0e}$ vs "
             f"$(1/16\\pi)^{{17}} = {root['ladder_17_factor']:.0e}$ - "
             f"composition misses by x{root['extrapolation_ratio']:.0e}",
             fontsize=10.5, color="0.25", va="bottom")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.30)
    save(fig, "fig29_loop_selection.png")

def fig30_w2_protection():
    import json
    with open(os.path.join(ROOT, "experiments", "V90_w2_protection",
                           "results", "v90_w2_protection.json")) as f:
        d = json.load(f)["data"]
    S17, bound = 0.146, 0.12

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.3, 1.0]})

    r = np.linspace(1.0, 3.0, 4001)
    with np.errstate(divide="ignore"):
        s = S17 / (1 + r) + 2 * S17 / np.abs(2 - r)
    axA.semilogy(r, s, color=C_GUT, lw=2.0,
                 label=r"$\sum m_\nu(r) = \frac{S_{17}}{1+r} + "
                       r"\frac{2 S_{17}}{|2-r|}$")
    axA.axhline(bound, color=C_ZERO, ls="--", lw=1.8,
                label=r"cosmological bound $\sum m_\nu < 0.12$ eV")
    axA.axvline(2.0, color="#2c7a3f", lw=1.6)
    axA.plot([2.0], [S17 / 3], "*", ms=17, color="#2c7a3f", zorder=6)
    axA.annotate(f"the ONLY passing point:\n$r = 2$, "
                 f"$\\sum m_\\nu = {S17/3:.4f}$ eV",
                 xy=(2.0, S17 / 3), xytext=(1.95, 0.02), fontsize=11.5,
                 color="#2c7a3f", weight="bold", ha="right",
                 arrowprops=dict(arrowstyle="->", color="#2c7a3f", lw=1.2))
    axA.annotate("any deviation makes the\nmassless states heavy "
                 "($1/|2-r|$)", xy=(2.05, 0.6), xytext=(2.3, 2.0),
                 fontsize=11, color=C_SU2,
                 arrowprops=dict(arrowstyle="->", color=C_SU2, lw=1.0))
    axA.set_xlabel(r"junction ratio $r = \kappa / w_0$")
    axA.set_ylabel(r"$\sum m_\nu$ [eV]")
    axA.set_ylim(0.012, 2e3)
    axA.set_title("(a)  the window is EXACTLY the point $r = 2$\n"
                  "     (1 passing value in 20001 scanned)", fontsize=12.5)
    axA.legend(loc="upper left", fontsize=10.5)

    axB.axis("off")
    axB.set_title("(b)  is the point protected? (V90)", fontsize=12.5)
    rows = [
        ("cubic action", "all of $S_3$ (6 permutations)",
         C_GUT),
        ("commutant", f"dim {d['commutant_nullity']} = "
         r"$\operatorname{span}\{I, J\}$ only", C_ZERO),
        ("quantization", f"one period "
         f"{d['w_period_2pi_over_maxspacing']:.2f} for $w_0$ AND $\\kappa$",
         C_ZERO),
        ("window width", f"{d['window_width']} (measure zero)", C_ZERO),
        ("mass regression", f"{d['eigenvalue_regression']:.0e} vs V76",
         "#2c7a3f"),
    ]
    for i, (k, v, c) in enumerate(rows):
        y = 0.92 - i * 0.14
        axB.text(0.02, y, k, fontsize=11.5, weight="bold", va="top")
        axB.text(0.52, y, v, fontsize=11.5, color=c, va="top")
    axB.text(0.02, 0.02,
             "W2: no symmetry enlarges the stabilizer of $\\kappa = 2w_0$,\n"
             "no quantization fixes the ratio, and the window is a\n"
             "point - the tuning is ONE exact number, unprotected.",
             fontsize=11.5, color=C_ROOT, weight="bold", va="bottom")

    fig.suptitle("The democratic neutrino point: exact, unprotected, "
                 "measure-zero (V90)", fontsize=13.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig30_w2_protection.png")

def fig31_w3_forcing():
    import json
    with open(os.path.join(ROOT, "experiments", "V91_w3_forcing",
                           "results", "v91_w3_forcing.json")) as f:
        d = json.load(f)["data"]
    mx = d["fixed_meeting_scale"]
    obs = 2.0e16

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.15, 1.0]})

    # (a) the 12 doubler nodes projected on (qx, qy), marker by zero axis
    nodes = [(1.5708, 1.5708, "z"), (1.5708, 4.7124, "z"),
             (4.7124, 1.5708, "z"), (4.7124, 4.7124, "z"),
             (1.5708, 1.5708, "y"), (1.5708, 4.7124, "y"),
             (4.7124, 1.5708, "y"), (4.7124, 4.7124, "y"),
             (1.5708, 0.0, "x"), (4.7124, 0.0, "x"),
             (0.0, 1.5708, "x"), (0.0, 4.7124, "x")]
    seen = set()
    for qx, qy, ax_ in nodes:
        mk = {"x": "o", "y": "s", "z": "^"}[ax_]
        col = {"x": C_SU3, "y": C_SU2, "z": C_GUT}[ax_]
        lab = f"zero axis {ax_}" if ax_ not in seen else None
        seen.add(ax_)
        jx = 0.10 if ax_ == "y" else (-0.10 if ax_ == "z" else 0.0)
        jy = 0.10 if ax_ == "x" else 0.0
        axA.plot(qx + jx, qy + jy, mk, ms=10, color=col, label=lab)
    for qx, qy in [(1.5708, 1.5708), (1.5708, 4.7124), (4.7124, 1.5708)]:
        axA.plot(qx, qy, "*", ms=19, color=C_ZERO, zorder=6)
    axA.annotate("the 3-node triangle (V76): realizable,\n"
                 "but its cubic orbit has 8 triangles",
                 xy=(1.5708, 1.5708), xytext=(0.3, 2.3), fontsize=11,
                 color=C_ZERO, weight="bold",
                 arrowprops=dict(arrowstyle="->", color=C_ZERO, lw=1.2))
    axA.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    axA.set_xticklabels(["0", r"$\pi/2$", r"$\pi$",
                         r"$3\pi/2$", r"$2\pi$"])
    axA.set_yticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    axA.set_yticklabels(["0", r"$\pi/2$", r"$\pi$",
                         r"$3\pi/2$", r"$2\pi$"])
    axA.set_xlabel(r"$q_x$")
    axA.set_ylabel(r"$q_y$")
    axA.set_title("(a)  the 12 doubler nodes: ONE cubic orbit of 12\n"
                  "     (the axis families are not orbits)", fontsize=12.5)
    axA.legend(loc="upper center", fontsize=10.5)

    # (b) fixed meeting scale vs size + verdict
    ns = [3, 4, 8, 12]
    vals = [mx[f"n_{n}"]["M_X_GeV"] for n in ns]
    cols = ["#2c7a3f" if v and abs(np.log10(v / obs)) <= 0.5 else C_ZERO
            for v in vals]
    axB.bar([str(n) for n in ns], vals, color=cols, width=0.6)
    axB.axhline(obs, color=C_ZERO, ls="--", lw=1.8,
                label=r"observed unification $2\times10^{16}$ GeV (F5)")
    for n, v in zip(ns, vals):
        axB.text(str(n), v * 1.15, f"{v:.1e}", ha="center", fontsize=11)
    axB.set_yscale("log")
    axB.set_ylim(1e14, 1e18)
    axB.set_xlabel("split-set size $n$ (cubic-symmetric: only 12)")
    axB.set_ylabel(r"fixed meeting $M_X(n)$ [GeV]")
    axB.set_title("(b)  every size misses the observed scale;\n"
                  "     the rank-1 rigidity holds for $n = 4, 8, 12$",
                  fontsize=12.5)
    axB.legend(loc="upper left", fontsize=10.5)
    fig.text(0.75, 0.015,
             "W3: cubic-symmetric split (12 nodes) EXCLUDED; every proper\n"
             "split BREAKS the cubic symmetry - the count 3 is an input.",
             fontsize=11, color=C_ROOT, weight="bold", ha="center",
             va="bottom")

    fig.suptitle("Is the family count 3 forced? No - the symmetric "
                 "alternative is excluded (V91)", fontsize=13.5)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig31_w3_forcing.png")

def fig33_w4_cone():
    import json
    with open(os.path.join(ROOT, "experiments", "V92_w4_feasibility",
                           "results", "v92_w4_feasibility.json")) as fh:
        d92 = json.load(fh)["data"]
    G = np.array(d92["db_generators"])          # (6, 3) fragment classes
    Tp = np.array(d92["target_Tp"])             # data-preferred meeting
    w = d92["witness_barycentric"]              # keys index {0, f_0..f_5}

    # project the 7 vertices {0, f_0..f_5} onto their plane
    V = np.vstack([np.zeros((1, 3)), G])
    cen = V.mean(axis=0)
    _u, sv, Vt = np.linalg.svd(V - cen, full_matrices=False)
    P = Vt[:2].T
    v2 = (V - cen) @ P
    tp2 = (Tp - cen) @ P

    def hull(pts):
        pts = sorted(map(tuple, pts))
        def half(Pts):
            H = []
            for p in Pts:
                while len(H) >= 2 and \
                        ((H[-1][0] - H[-2][0]) * (p[1] - H[-2][1])
                         - (H[-1][1] - H[-2][1]) * (p[0] - H[-2][0])) <= 0:
                    H.pop()
                H.append(p)
            return H
        return half(pts) + half(reversed(pts))[1:-1]

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    h = hull(v2)
    ax.add_patch(plt.Polygon(h, closed=True, fc="#dce9f5", ec=C_ALG,
                             lw=1.8, zorder=1))
    ctr = v2.mean(axis=0)
    for i, p in enumerate(v2[1:], start=0):
        ax.plot([v2[0][0], p[0]], [v2[0][1], p[1]], color=C_ALG, lw=1.1,
                zorder=2)
        ax.plot(*p, "o", color=C_ALG, ms=8, zorder=3)
        d = p - ctr
        d = d / max(np.linalg.norm(d), 1e-9)
        off = (0.15, -0.38) if i == 0 else (0.25 * d[0], 0.25 * d[1])
        ax.annotate(f"$f_{i}$", p, xytext=(p[0] + off[0],
                    p[1] + off[1]), fontsize=13, color=C_ALG)
    ax.plot(*v2[0], "o", color="0.3", ms=8, zorder=3)
    ax.annotate("$0$ (no extra content)", v2[0],
                xytext=(v2[0][0] + 0.12, v2[0][1] - 0.12), fontsize=12,
                color="0.3", ha="left")
    ax.plot(*tp2, marker="*", ms=24, color=C_FIT, markeredgecolor="black",
            lw=0.7, zorder=4)
    lam = ", ".join(
        f"$\\lambda_0$={v:.2f}" if k == "0"
        else f"$\\lambda_{{f_{int(k) - 1}}}$={v:.2f}"
        for k, v in w.items())
    ax.annotate("data-preferred meeting (V64):\n"
                f"$M_X=1.8\\times10^{{16}}$ GeV, $\\alpha_X=1/25$ -\n"
                "INSIDE the hull.\n"
                f"witness: {lam}",
                tp2, xytext=(-2.0, -1.55),
                fontsize=12, color=C_FIT, ha="left",
                arrowprops=dict(arrowstyle="->", color=C_FIT))
    ax.set_title("V92: the desert content the data require lies INSIDE "
                 "the hull of the node-$\\mathbf{16}$ fragments\n"
                 "($\\beta$-contribution plane; feasibility $\\Rightarrow$ "
                 "W4 is a realizability problem, not a no-go)",
                 fontsize=13.5)
    ax.set_xlabel("$\\beta_1,\\beta_2,\\beta_3$ plane (PCA)")
    ax.grid(True, alpha=0.3, ls=":")
    ax.set_xlim(v2[:, 0].min() - 0.35, v2[:, 0].max() + 1.85)
    ax.set_ylim(v2[:, 1].min() - 1.45, v2[:, 1].max() + 0.30)
    fig.tight_layout()
    save(fig, "fig33_w4_cone.png")

def fig36_vertex_ground_state():
    import json
    p100 = os.path.join(ROOT, "experiments", "V100_vertex_correction",
                        "results", "v100_vertex_correction.json")
    p101 = os.path.join(ROOT, "experiments",
                        "V101_bond_current_ground_state",
                        "results", "v101_bond_current_ground_state.json")
    d100 = json.load(open(p100))["data"]
    d101 = json.load(open(p101))["data"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # (a) the entropy of one Planck bond, from naive counting to Bekenstein
    I6 = d100["I"]["6"][0]
    Scap = d100["S_capacity"]
    labels = ["$2\\ln 45$\nnaive",
              "$S_{\\rm cap}$\nkernel cap.",
              "$S_{\\rm cap}/2$\nper bond",
              "$(1/3)\\ln I_6$\nvertex coeff.",
              "$1/4$\nBekenstein"]
    vals = [2 * math.log(45), Scap, Scap / 2,
            d100["vertex_coefficient"], 0.25]
    cols = [C_OPEN, C_ALG, C_ALG, C_FIT, C_ZERO]
    bars = ax1.bar(range(5), vals, color=cols, width=0.62)
    for k, v in enumerate(vals):
        ax1.text(k, v + 0.15, f"{v:.2f}", ha="center", fontsize=12)
    ax1.axhline(0.25, color=C_ZERO, ls="--", lw=1.2)
    ax1.set_xticks(range(5))
    ax1.set_xticklabels(labels, fontsize=10.5)
    ax1.set_ylabel("nats per Planck bond")
    ax1.set_ylim(0, 8.6)
    ax1.set_title("(a) the vertex correction: deficit factor "
                  "$15.2 \\to 6.6$", fontsize=13)
    ax1.text(0.30, 0.97,
             "$I_n=\\dim(\\mathbf{45}^{\\otimes n})^G$: "
             f"{round(d100['I']['2'][0])}, {round(d100['I']['3'][0])}, "
             f"{round(d100['I']['4'][0])}, {round(d100['I']['5'][0])}, "
             f"{round(I6)}",
             transform=ax1.transAxes, fontsize=11.5, va="top",
             color=C_ALG)

    # (b) the 24-level spectrum of the bond-current Hamiltonian
    spec = d101["spectrum"]
    E0 = d101["E0"]
    for k, e in enumerate(spec):
        ax2.plot([k - 0.35, k + 0.35], [e, e], color=C_ALG, lw=2.2)
    ax2.axhline(E0, color=C_FIT, ls="--", lw=1.1)
    ax2.set_xlim(-0.8, len(spec) - 0.2)
    ax2.set_ylim(E0 - 60, E0 + 180)
    ax2.annotate(f"$E_0 = {E0:.1f}\\,\\eta$, deg. {d101['degeneracy']}",
                 xy=(0.2, E0 + 8), xytext=(7.5, E0 + 12),
                 fontsize=11.5, color=C_FIT)
    ax2.annotate(f"gap $\\Delta = {d101['gap']:.1f}\\,\\eta$",
                 xy=(2.6, E0 + d101["gap"] / 2), xytext=(0.0, E0 + 150),
                 fontsize=11.5, color=C_ALG,
                 arrowprops=dict(arrowstyle="<->", color=C_ALG, lw=1))
    ax2.text(0.02, 0.06,
             f"$S_{{\\rm GS}} = {d101['S_ground']:.4f} = \\ln 45$ "
             "(middle cut)",
             transform=ax2.transAxes, fontsize=11.5, va="bottom",
             color=C_ALG)
    ax2.set_xlabel("kernel state (Gram-orthonormalised)")
    ax2.set_ylabel("$E$ (units of $\\eta$)")
    ax2.set_title("(b) the first ground state: $J\\!\\cdot\\!L$ splits the "
                  "kernel, $S_{\\rm GS} = \\ln 45$ (V101)", fontsize=13)
    fig.suptitle("The gauged area law corrected by vertices, and the "
                 "ground state of the gauged layer (V100, V101)",
                 fontsize=13.5)
    fig.tight_layout()
    save(fig, "fig36_vertex_ground_state.png")
