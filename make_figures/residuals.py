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
         ["2.08 (hierarchy; V108: $3F_{10}^{10}$, conditional)", "1.18 (DM mass)", r"$n_{2\text{step}} = 3$",
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
        d90 = json.load(f)["data"]
    with open(os.path.join(ROOT, "experiments", "V103_w2_democratic",
                           "results", "v103_w2_democratic.json")) as f:
        d103 = json.load(f)["data"]
    S17, bound = 0.146, 0.12

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.3, 1.0]})

    r = np.linspace(0.05, 3.0, 4001)
    with np.errstate(divide="ignore"):
        s_id = S17 / (1 + r) + 2 * S17 / np.abs(2 - r)
    s_demo = S17 / (1 + r)
    axA.semilogy(r, s_id, color=C_GUT, lw=2.0,
                 label=r"identity $m_D$ (V90): $\frac{S_{17}}{1+r} + "
                       r"\frac{2S_{17}}{|2-r|}$")
    axA.semilogy(r, s_demo, color="#2c7a3f", lw=2.4,
                 label=r"democratic $m_D$ (V103): $m_3 = "
                       r"\frac{S_{17}}{1+r}$, two massless for every $r$")
    axA.axhline(bound, color=C_ZERO, ls="--", lw=1.8,
                label=r"cosmological bound $\sum m_\nu < 0.12$ eV")
    r_min = S17 / bound - 1.0
    axA.axvspan(r_min, 3.0, color="#2c7a3f", alpha=0.07)
    axA.axvline(2.0, color="#2c7a3f", lw=1.2, alpha=0.5)
    axA.plot([2.0], [S17 / 3], "*", ms=17, color="#2c7a3f", zorder=6)
    r_atm = d103["r_atm"]
    axA.plot([r_atm], [bound * 0.417], "o", ms=10, color=C_SU2, zorder=6)
    axA.annotate(f"atmospheric scale selects\n$r = {r_atm:.2f}$ "
                 f"($m_3 = {S17/(1+r_atm):.4f}$ eV)",
                 xy=(r_atm, S17 / (1 + r_atm)), xytext=(0.35, 0.4),
                 fontsize=11, color=C_SU2,
                 arrowprops=dict(arrowstyle="->", color=C_SU2, lw=1.0))
    axA.annotate(f"window $r > {r_min:.2f}$ (width "
                 f"{d103['window_width_democratic']:.2f})",
                 (2.6, 0.016), fontsize=11.5, color="#2c7a3f",
                 weight="bold", ha="center")
    axA.set_xlabel(r"junction ratio $r = \kappa / w_0$")
    axA.set_ylabel(r"$\sum m_\nu$ [eV]")
    axA.set_ylim(0.012, 2e3)
    axA.set_title("(a)  V90 scanned $m_D \\propto I$ (divergent, one point);\n"
                  "     the implemented democratic $m_D$ gives a finite window",
                  fontsize=12.5)
    axA.legend(loc="upper left", fontsize=9.5)

    axB.axis("off")
    axB.set_title("(b)  what protects what (V90 + V103)", fontsize=12.5)
    rows = [
        ("cubic action", "all of $S_3$ (6 permutations)", C_GUT),
        ("commutant", f"dim {d90['commutant_nullity']} = "
         r"$\operatorname{span}\{I, J\}$ only", C_ZERO),
        ("$m_D$ (V28/V45)", "democratic, rank 1 - 0-parameter", "#2c7a3f"),
        ("massless pair", "protected for EVERY $r$ (rank-1 $m_D$)",
         "#2c7a3f"),
        ("window", f"$r > {r_min:.2f}$, width "
         f"{d103['window_width_democratic']:.2f} (was: point)", "#2c7a3f"),
        ("$r$ from data", f"atmospheric: $r = {r_atm:.2f}$; "
         "$r = 2$: $S_{17}/3$", C_SU2),
    ]
    for i, (k, v, c) in enumerate(rows):
        y = 0.94 - i * 0.108
        axB.text(0.02, y, k, fontsize=11, weight="bold", va="top")
        axB.text(0.42, y, v, fontsize=11, color=c, va="top")
    axB.text(0.02, 0.02,
             "W2 re-audited: the measure-zero tuning was an artifact of\n"
             "the identity-$m_D$ assumption. The structure (2 massless +\n"
             "1 massive, NO) is protected by the same family symmetry\n"
             "that derives $\\sin^2\\theta_{12} = 1/3$; the value sits in a\n"
             "finite window. Solar mass still needs W6.",
             fontsize=10, color=C_ROOT, weight="bold", va="bottom")

    fig.suptitle("The democratic neutrino point: the V90 tuning is an "
                 "artifact of the wrong Dirac assumption (V103)",
                 fontsize=13.5)
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


def fig37_residual_attacks():
    """Attacking the residuals head-on: V109 (W5 CP/mirror), V110 (W4 bond
    off-diagonal), V111 (W6 generation dressing), V112 (W7 dark U(1)),
    V113 (W8 axis factorization), V116 (W3 triangle orbit-8)."""
    import json

    def _load(v):
        p = os.path.join(ROOT, "experiments", v, "results",
                         v.lower() + ".json")
        return json.load(open(p))["data"]

    d109 = _load("V109_w5_cp_search")
    d110 = _load("V110_w4_offdiagonal_bond")
    d111 = _load("V111_w6_generation_bond")
    d112 = _load("V112_w7_dark_u1_portal")
    d113 = _load("V113_w8_axis_factorization")
    d116 = _load("V116_w3_triangle_orbit8")

    fig, axs = plt.subplots(2, 3, figsize=(14.4, 8.8))
    (ax1, ax2, ax3), (ax4, ax5, ax6) = axs

    # (a) V109: the mirror is the CP-forbidding structure
    cases = d109["per_case_worst"]
    names = ["Wilson\n(V71)", r"$\kappa_s$ anti", r"$\kappa_s$ sym",
             "dressing\nanti", "dressing\nsym", "combined"]
    keys = ["kappa_antisym", "kappa_sym", "dressing_antisym",
            "dressing_sym", "both"]
    vals = [d109["base_pairing_dev"]] + [cases[k] for k in keys]
    cols = [C_ALG] + [C_OPEN] * 5
    bars = ax1.bar(range(6), vals, color=cols, width=0.62)
    ax1.set_yscale("log")
    ax1.set_ylim(1e-17, 60)
    for j, (b, v) in enumerate(zip(bars, vals)):
        if j == 0:
            ax1.text(b.get_x() + b.get_width() / 2, v * 2.0, f"{v:.1e}",
                     ha="center", va="bottom", fontsize=10)
        else:
            ax1.text(b.get_x() + b.get_width() / 2, v * 0.35, f"{v:.1e}",
                     ha="center", va="top", fontsize=10, color="white",
                     rotation=90)
    ax1.axhline(1e-10, color=C_ZERO, ls="--", lw=1.2)
    ax1.text(0.02, 4e-10, "mirror preserved", fontsize=10, color=C_ZERO)
    ax1.set_xticks(range(6))
    ax1.set_xticklabels(names, fontsize=9, rotation=28, ha="left")
    ax1.set_ylabel(r"$\pm$-pairing deviation of the 32 spectrum")
    ax1.set_title("(a) W5: the mirror forbids CP — only Cartan-type\n"
                  "extensions preserve it; class-dependent dressings\n"
                  "break Theorem 4 (V109)", fontsize=12)

    # (b) V110: the off-diagonal bond floor
    per_n = d110["per_n"]
    ns = sorted(int(k) for k in per_n)
    res = [per_n[str(n)]["residual"] for n in ns]
    ax2.plot(ns, res, "o-", color=C_OPEN, lw=1.8)
    ax2.axhline(0.0, color=C_ALG, lw=1.6)
    ax2.annotate(f"floor {min(res):.1f} at n={ns[res.index(min(res))]}\n"
                 "the meeting does NOT close",
                 xy=(ns[res.index(min(res))], min(res)), xytext=(6.5, 12),
                 fontsize=10.5, color=C_OPEN,
                 arrowprops=dict(arrowstyle="->", color=C_OPEN, lw=1))
    nforced = sum(1 for f in d110["forced_signs"]
                  if f["reachable_sign"] != "mixed")
    ax2.text(0.02, 0.97, f"level repulsion forces {nforced}/6 signs;\n"
             r"the six $\kappa_s$ remain a named input",
             transform=ax2.transAxes, ha="left", va="top", fontsize=10.5,
             color=C_FIT)
    ax2.set_xlabel("desert pairs $n$")
    ax2.set_ylabel("meeting residual")
    ax2.set_title("(b) W4: one Hermitian bond operator with 12\n"
                  "parameters does not close the meeting (V110)",
                  fontsize=12)

    # (c) V111: named-parameter cost of reaching the measured mixing
    labels = ["scalar dressing\n(exactly blind)",
              "common family\ndressing",
              "chiral-asymmetric\ndressing"]
    costs = [float("nan"), 7, 13]
    ax3.bar([1, 2], [7, 13], color=[C_FIT, C_FIT], width=0.55)
    ax3.bar([0], [1], color=C_ALG, width=0.55, hatch="///")
    ax3.text(0, 8.4, "never reached\n(0 parameters help)", ha="center",
             fontsize=10, color=C_ALG)
    ax3.text(1, 7.3, "7 parameters\nREACHED", ha="center", fontsize=10.5)
    ax3.text(2, 13.3, "13 parameters\nREACHED as a fit", ha="center",
             fontsize=10.5)
    ax3.axhline(3, color=C_ZERO, ls="--", lw=1.3)
    ax3.text(0.52, 3.45, "3 constraints\n(measured angles)", fontsize=9,
             color=C_ZERO, ha="left", va="bottom")
    ax3.set_xticks(range(3))
    ax3.set_xticklabels(labels, fontsize=10)
    ax3.set_ylim(0, 17.5)
    ax3.set_ylabel("named parameters needed")
    ax3.set_title("(c) W6: generation-space dressing moves the\n"
                  "mixing — but every entry is a named input (V111)",
                  fontsize=12)

    # (d) V112: the dark-U(1) window squeeze
    e2 = 4 * math.pi / 127.9
    mu_N, ms = 0.94, 118.0
    target_sv = 2.2e-26 / 1.168e-17
    dd_lim = 2.0e-47 / 0.3894e-27
    mzp = np.geomspace(20.0, 5000.0, 200)
    gD = 1.0
    eps_relic = np.sqrt(target_sv * math.pi
                        * (mzp ** 2 + ms ** 2) ** 2
                        / (gD ** 2 * e2 * 3.0 * 4 * ms ** 2))
    eps_dd = np.sqrt(dd_lim * math.pi * mzp ** 4
                     / (gD ** 2 * e2 * mu_N ** 2))
    ax4.plot(mzp, eps_relic, color=C_FIT, lw=1.8,
             label=r"relic $\Omega h^2$: $\varepsilon$ needed")
    ax4.plot(mzp, eps_dd, color=C_OPEN, lw=1.8,
             label=r"DD bound $\sigma_{\rm SI} < 2\times10^{-47}$: ceiling")
    ax4.fill_between(mzp, eps_dd, eps_relic, facecolor="none",
                     edgecolor=C_OPEN, hatch="///", alpha=0.6,
                     label="no overlap: window empty")
    ax4.axvspan(20, 2900, ymin=0.93, color=C_CONV, alpha=0.3)
    ax4.text(0.03, 0.90, "LEP-class $Z'$ reach", fontsize=9.5,
             color=C_CONV, transform=ax4.transAxes, va="top")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_ylim(1e-8, 1e-1)
    ax4.set_xlabel(r"$m_{Z'}$ [GeV]  ($g_D = 1$)")
    ax4.set_ylabel("kinetic mixing $\\varepsilon$")
    ax4.legend(fontsize=9.5, loc="lower left")
    ax4.set_title("(d) W7: dark U(1) portal — relic needs\n"
                  r"$g_D\varepsilon \gtrsim 0.02$, DD allows "
                  r"$\lesssim 10^{-3}$:" "\nwindow EMPTY (V112)",
                  fontsize=12)

    # (e) V113: the coefficient 3 is the node-set choice
    sets = ["all 12\nnodes", "one split node\nper family", "single\nnode"]
    coef = [d113["coef_all12"], d113["coef_one_per_family"],
            d113["coef_single"]]
    obs = d113["coef_observed"]
    cols = [C_OPEN, C_ALG, C_OPEN]
    bars = ax5.bar(range(3), coef, color=cols, width=0.58)
    ax5.axhline(obs, color="#9a9a9a", ls="--", lw=1.6)
    ax5.text(2.42, obs * 1.03, f"observed {obs:.4f}", fontsize=10.5,
             color="#666", ha="right", va="bottom")
    for b, v in zip(bars, coef):
        ax5.text(b.get_x() + b.get_width() / 2, v + 0.12, f"{v:.4f}",
                 ha="center", fontsize=10.5)
    ax5.set_xticks(range(3))
    ax5.set_xticklabels(sets, fontsize=10)
    ax5.set_ylim(0, 9.6)
    ax5.set_ylabel(r"hierarchy coefficient $\sum$-nodes")
    ax5.set_title("(e) W8: node contributions EXACTLY equal across\n"
                  "families (0.00e+00); the coefficient 3 IS the choice\n"
                  "one-per-family — 99.4% agreement (V113)", fontsize=12)

    # (f) V116: the triangle and its orbit-8
    nodes = [(x, y) for x in (0, 1, -1) for y in (0, 1, -1)
             if (x != 0) + (y != 0) == 2]
    # plot the 12 nodes projected: (qx,qy) plane shows 8; mark families
    h = math.pi / 2
    plane = [(a * h, b * h) for a in (1, -1) for b in (1, -1)]
    outz = [(0.0, 0.0)]
    ax6.scatter([p[0] for p in plane], [p[1] for p in plane], s=90,
                color="#9a9a9a", zorder=3)
    tri = [(h, h), (h, 0.0), (0.0, h)]
    ax6.scatter([p[0] for p in tri], [p[1] for p in tri], s=170,
                marker="*", color=C_ALG, zorder=4)
    ax6.plot([p[0] for p in tri] + [tri[0][0]],
             [p[1] for p in tri] + [tri[0][1]], color=C_ALG, lw=1.4,
             zorder=2)
    for p in [(-h, h), (h, -h), (-h, -h)]:
        ax6.scatter([p[0]], [p[1]], s=90, color="#9a9a9a", zorder=3)
    ax6.annotate("orbit of 8 triangles\nstabilizer $S_3$ (order 6)",
                 xy=(h, h), xytext=(-h * 0.9, h * 0.55), fontsize=10.5,
                 color=C_ALG,
                 arrowprops=dict(arrowstyle="->", color=C_ALG, lw=1))
    ax6.set_xlim(-2.2, 2.2)
    ax6.set_ylim(-2.2, 2.2)
    ax6.set_xlabel(r"$q_x$")
    ax6.set_ylabel(r"$q_y$")
    ax6.set_title("(f) W3: 64 triangles → 4 cubic orbits "
                  "(8, 8, 24, 24);\n"
                  r"each choice breaks $O_h \to S_3$ — 8 distinct "
                  "vacua (V116)", fontsize=12)

    fig.suptitle("Attacking the residuals head-on: what V109-V116 close, "
                 "what they name", fontsize=14.5)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    save(fig, "fig37_residual_attacks.png")


def fig38_doubler_cones():
    """V102: the twelve doubler cones are exactly linear but never round."""
    import json
    d = json.load(open(os.path.join(
        ROOT, "experiments", "V102_lorentz_isotropy", "results",
        "v102_lorentz_isotropy.json")))["data"]
    nodes = d["per_node"]

    fig = plt.figure(figsize=(13.6, 4.5))
    gs = fig.add_gridspec(1, 3, wspace=0.30)

    # (a) one node: the velocity surface v(n) over 100 directions,
    #     inside the causal sphere |v| = 1
    ax = fig.add_subplot(gs[0], projection="3d")
    vd = np.array(nodes[0]["v_dirs"])
    pts = vd[:, :3] * vd[:, 3:4]
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=vd[:, 3],
               cmap="viridis", s=16, depthshade=False)
    u = np.linspace(0, 2 * np.pi, 25)
    v = np.linspace(0, np.pi, 13)
    sx, sy, sz = (np.outer(np.cos(u), np.sin(v)),
                  np.outer(np.sin(u), np.sin(v)),
                  np.outer(np.ones_like(u), np.cos(v)))
    ax.plot_wireframe(sx, sy, sz, color=C_OPEN, lw=0.4, alpha=0.35)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
    ax.set_zlim(-1.05, 1.05)
    ax.set_xlabel(r"$\hat n_x$"); ax.set_ylabel(r"$\hat n_y$")
    ax.set_zlabel(r"$\hat n_z$")
    ax.view_init(elev=14, azim=32)
    ax.set_title("(a) one node: raw $|v(\\hat n)|$ (min eigenphase)\n"
                 "over 100 directions, inside the causal sphere\n"
                 "(grey, speed 1) - never round (V102)", fontsize=12)

    # (b) all 12 nodes: the 1200 measured velocities
    ax2 = fig.add_subplot(gs[1])
    for i, n in enumerate(nodes):
        vv = np.array(n["v_dirs"])[:, 3]
        ax2.scatter(np.full_like(vv, i) + 0.06 * ((i % 2) - 0.5), vv,
                    s=10, color=C_GUT, alpha=0.55)
    vp = np.array(nodes[0]["v_principal"])
    ax2.axhspan(float(vp[0]), float(vp[2]), color=C_GUT, alpha=0.15)
    ax2.annotate("tensor fit: principal band\n0.20 / 0.20 / 0.28",
                 xy=(5.5, 0.26), xytext=(5.5, 0.62), ha="center",
                 fontsize=10.5, color=C_GUT,
                 arrowprops=dict(arrowstyle="->", color=C_GUT, lw=1))
    ax2.axhline(1.0, color=C_ZERO, ls="--", lw=1.6)
    ax2.text(11.4, 1.02, "causal speed 1", ha="right", fontsize=10.5,
             color=C_ZERO)
    ax2.axhline(0.0, color="k", lw=0.8)
    ax2.set_xticks(range(12))
    ax2.set_xticklabels([f"{i}" for i in range(12)], fontsize=10)
    ax2.set_xlabel("doubler node (one $O_h$ orbit)")
    ax2.set_ylabel(r"$|v(\hat n)|$  (lattice units)")
    ax2.set_ylim(0, 1.08)
    ax2.set_title("(b) 12 nodes x 100 directions: 1200/1200 linear;\n"
                  "the velocity tensor fit gives principal cone\n"
                  "0.20/0.20/0.28 - a factor ~4 below causal (V102)",
                  fontsize=12)

    # (c) anisotropy per node: the orbit is uniform
    ax3 = fig.add_subplot(gs[2])
    an = [n["anisotropy"] for n in nodes]
    ax3.bar(range(12), an, color=C_FIT, width=0.62)
    ax3.axhline(float(np.mean(an)), color=C_ALG, ls="--", lw=1.5)
    ax3.text(11.4, np.mean(an) + 0.012, f"mean {np.mean(an):.3f}",
             ha="right", fontsize=10.5, color=C_ALG)
    ax3.set_ylim(0, 0.9)
    ax3.set_xticks(range(12))
    ax3.set_xticklabels([f"{i}" for i in range(12)], fontsize=10)
    ax3.set_xlabel("doubler node")
    ax3.set_ylabel("cone anisotropy (eigen-spread / mean)")
    ax3.set_title("(c) every node of the orbit carries the SAME\n"
                  "anisotropic cone: 0.748-0.758 - the lattice\n"
                  "signature is universal, not a defect (V102)",
                  fontsize=12)

    fig.suptitle("The twelve doubler cones: exactly linear, never Lorentz-"
                 "round (V102)", fontsize=14.5)
    fig.subplots_adjust(left=0.03, right=0.99, bottom=0.12, top=0.79)
    save(fig, "fig38_doubler_cones.png")


def fig39_isotropy_search():
    """V115: no isotropy point on the whole V65 Wilson family."""
    import json
    d = json.load(open(os.path.join(
        ROOT, "experiments", "V115_w9_cone_isotropy", "results",
        "v115_w9_cone_isotropy.json")))["data"]
    pts = d["gapless_points"]
    A = np.array([p[0] for p in pts])
    th = np.array([p[2] for p in pts])
    w = np.array([p[3] for p in pts])
    n = len(A)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.6, 4.5))

    # (a) the search space: every gapless point, and the empty target
    ax1.axhspan(0.0, 1e-6, color=C_SU3, alpha=0.9)
    ax1.axhspan(1e-6, 0.05, color=C_SU3, alpha=0.15)
    sc = ax1.scatter(w, A, c=th, cmap="plasma", s=9, alpha=0.8)
    fig.colorbar(sc, ax=ax1, pad=0.02, label=r"Wilson angle $\theta$")
    ax1.annotate("Lorentz target $A<10^{-6}$:\n0 / %d points" % n,
                 xy=(1.5, 0.035), xytext=(1.5, 0.30), fontsize=11.5,
                 color=C_SU3, ha="center",
                 arrowprops=dict(arrowstyle="->", color=C_SU3, lw=1.4))
    amin = float(d["min_anisotropy"])
    ax1.annotate(f"closest approach $A$ = {amin:.3f}", xy=(w[np.argmin(A)],
                 amin), xytext=(0.35, 0.72), fontsize=11.5, color=C_ZERO,
                 arrowprops=dict(arrowstyle="->", color=C_ZERO, lw=1.4))
    ax1.set_xlabel(r"Wilson strength $w$")
    ax1.set_ylabel(r"cone anisotropy $A=\frac{v_{\max}-v_{\min}}{v_{\max}+v_{\min}}$")
    ax1.set_ylim(0, 1.03)
    ax1.set_title(f"(a) {n} gapless points of the family\n"
                  "(2 reps x 60 angles, all zeros located):\n"
                  "the isotropy band is empty (V115)", fontsize=12)

    # (b) the distribution of A: pinned at full anisotropy
    ax2.hist(A, bins=np.arange(0.90, 1.005, 0.01), color=C_FIT)
    ax2.axvline(amin, color=C_ZERO, ls="--", lw=1.5)
    ax2.text(amin - 0.004, ax2.get_ylim()[1] * 0.94,
             f"min {amin:.3f}", ha="right", va="top", fontsize=11,
             color=C_ZERO)
    ax2.set_xlabel(r"cone anisotropy $A$")
    ax2.set_ylabel("gapless points")
    ax2.set_xlim(0.895, 1.005)
    ax2.set_title("(b) every point sits at $A \\gtrsim 0.90$;\n"
                  "isotropy ($A\\to 0$) is not approached\n"
                  "anywhere on the family (V115)", fontsize=12)

    # (c) no angle helps: A vs theta
    ax3.scatter(th, A, c=w, cmap="viridis", s=9, alpha=0.8)
    ax3.axhline(float(d["bare_cone_max_anisotropy"]), color=C_OPEN,
                ls=":", lw=1.5)
    ax3.text(6.15, float(d["bare_cone_max_anisotropy"]) + 0.004,
             "bare cone (w = 0)", ha="right", fontsize=10.5, color=C_OPEN)
    ax3.set_xlabel(r"Wilson angle $\theta$ (rad)")
    ax3.set_ylabel(r"cone anisotropy $A$")
    ax3.set_ylim(0.895, 1.005)
    ax3.set_title("(c) all 60 angles x all strengths: the Wilson\n"
                  "factor moves the zeros, never the cone's\n"
                  "shape - W9 stays open (V115)", fontsize=12)

    fig.suptitle("Is continuous Lorentz invariance emergent anywhere on "
                 "the Wilson family? No (V115)", fontsize=14.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig39_isotropy_search.png")
