"""Figure module: couplings."""
import math

import matplotlib.pyplot as plt
import numpy as np

from .style import exp, save
from .style import C_GUT, C_ROOT, C_SU2, C_SU3, C_ZERO
from cl10 import config as _config      # the shared inputs (no engine import here)


def fig19_couplings_mz():
    from cl10 import unification as u

    rd = u.running_down()
    az = u.A_Z
    p0 = rd["couplings_at_M_Z"]
    names = [r"$\alpha_1$  U(1)", r"$\alpha_2$  SU(2)", r"$\alpha_3$  SU(3)"]
    fig, axB = plt.subplots(figsize=(5.8, 3.6))
    # at M_Z: measured vs the down-run from the reference reading g = 1/5 (paired bars)
    m = [0, 1, 2]
    w = 0.36
    axB.bar([x - w / 2 - 0.02 for x in m], az, width=w, color="#9a9a9a",
            edgecolor="k", lw=0.5, label="measured (PDG 2024)")
    p0plot = [max(p, 0.0) for p in p0]
    axB.bar([x + w / 2 + 0.02 for x in m], p0plot, width=w, color=C_ZERO,
            edgecolor="k", lw=0.5, label="reference reading $g = 1/5$: $\\alpha_{10} \\equiv g^2 = 1/25$")
    for x in m:
        axB.text(x - w / 2, az[x] + 0.003, f"{az[x]:.3f}",
                 ha="center", fontsize=11.5, color="0.3")
        if p0[x] <= 0:
            axB.text(x + w / 2 + 0.05, 0.003, "below\npole", ha="left", fontsize=11,
                     color="0.4")
        else:
            axB.text(x + w / 2, p0[x] + 0.003, f"{p0[x]:.3f}",
                     ha="center", fontsize=11.5, color="0.45")
    axB.set_xticks(m)
    axB.set_xticklabels(names)
    axB.set_ylabel("$\\alpha_i$  at $M_Z$")
    axB.set_ylim(0, 0.26)
    axB.set_xlim(-0.6, 3.35)
    axB.set_title(r"At $M_Z$: $\alpha_1$ matches ($1.17\times$); $\alpha_2$ $5.8\times$ too large;"
                  "\n" r"$\alpha_3$ negative — the 1-loop matching is the O4 named residual", fontsize=12)
    axB.legend(fontsize=9.5, loc="upper right", framealpha=0.95)
    fig.suptitle(r"Couplings at $M_Z$: the down-run uses the reference reading $\alpha_{10} \equiv g^2 = 1/25$;"
                 "\n" r"the closed clock of Construction 1 is $\alpha = g^2/4\pi = \pi/64 = 0.0491$ (Theorem 6)",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig19_couplings_mz.png")

def fig4_uv_2loop():
    """The 2-loop running of the U(1) coupling: the perturbative QFT reaches a
    Landau pole below M_P, the 1-loop running stays bounded, and the unitary QCA
    step remains bounded at every order."""
    from cl10 import unification as u, u10pred as tp, uvrunning as uv

    M_Z, M_P = u.M_Z, tp.M_PLANCK
    b1, b2 = u.B_SM[0], u.B2_SM[0]
    a1z = u.A_Z[0]
    c = 1.0 / a1z
    Bcoef = b1 / (2.0 * np.pi)
    Acoef = b2 / (b1 ** 2 * 2.0 * np.pi)
    pole2 = uv.two_loop_pole_u1()   # the library value, not a local copy

    t = np.linspace(0.0, np.log(3.0e20 / M_Z), 400)
    M = M_Z * np.exp(t)
    inv1 = c - Bcoef * t
    inv2 = c - Bcoef * t - Acoef * t ** 2
    valid = inv2 > 0

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.semilogx(M, inv1, "-", color=C_SU2, lw=2.2,
                label=r"1-loop $1/\alpha_1(M)$ (pole above $M_P$)")
    ax.semilogx(M[valid], inv2[valid], "-", color=C_ZERO, lw=2.4,
                label=r"2-loop $1/\alpha_1(M)$ (pole below $M_P$)")
    ax.axvline(M_P, color="k", lw=1.2)
    ax.text(M_P * 1.05, 57, r"$M_P$", fontsize=11.5, va="bottom", ha="left")
    ax.axvline(pole2, color=C_ZERO, ls=":", lw=1.4)
    ax.scatter([pole2], [0], s=130, marker="x", color=C_ZERO, zorder=5,
               label=r"2-loop Landau pole $\sim$ " + f"{pole2:.2e}" + " GeV")
    ax.axhline(1.0 / u.ALPHA_10, color="0.35", ls="--", lw=1.3,
               label=r"reference reading  $1/\alpha_{10} = 25$  ($g = 1/5$, a convention)")
    fig.text(0.5, 0.012,
             r"the 2-loop perturbative $QFT$ reaches a Landau pole below $M_P$;  the $QCA$ step is"
             " unitary,\nso its $n$-step amplitudes stay bounded — the perturbative UV divergence"
             " is not a property of the QCA itself",
             ha="center", va="bottom", fontsize=10,
             bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#888", alpha=0.95))
    ax.set_xlabel(r"scale $M$  (GeV, log scale)")
    ax.set_ylabel(r"$1/\alpha_1$  (the $U(1)$ coupling)")
    ax.set_title("2-loop running of the $U(1)$ coupling: the perturbative QFT pole "
                 "vs the bounded QCA step")
    ax.set_xlim(70, 3e20)
    ax.set_ylim(-3, 60)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
              columnspacing=1.2, handlelength=1.6,
              framealpha=0.95, fontsize=9.5)
    fig.subplots_adjust(bottom=0.38)
    save(fig, "fig4_uv_2loop.png")

def fig20_desert():
    import json
    from cl10 import unification as u

    with open(exp("V64_desert",
                           "results", "v64_desert.json")) as f:
        sel = json.load(f)["data"]["selected"]
    MZ, MP, A10 = u.M_Z, u.M_P, u.ALPHA_10
    M_X, M_A, M_B = sel["M_X_GeV"], sel["M_A_GeV"], sel["M_B_GeV"]
    d10 = (6.0 / 5.0, 2.0, 2.0)          # one Dirac 10+10bar pair
    d2 = (2.0 / 5.0, 2.0 / 3.0, 0.0)     # one Dirac (1,2)+(3bar,1)... (1,2) pair
    n10 = int(sel["content"].split("x10")[0])   # the desert content, from V64
    cols = [C_GUT, C_SU2, C_SU3]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.5),
                                   gridspec_kw={"width_ratios": [1.45, 1.0]})
    # (a) piecewise 1-loop running with the desert content vs SM-only
    s = np.logspace(np.log10(MZ), np.log10(MP), 400)
    L0 = np.log(np.maximum(s, MZ) / MZ)
    LA = np.maximum(0.0, np.log(np.maximum(s, M_A) / M_A))
    LB = np.maximum(0.0, np.log(np.maximum(s, M_B) / M_B))
    for k in range(3):
        sm = 1.0 / u.A_Z[k] - (u.B_SM[k] / (2.0 * np.pi)) * L0
        axA.plot(np.log10(s), sm, color=cols[k], lw=1.0, ls="--", alpha=0.5)
        ds = 1.0 / u.A_Z[k] - (1.0 / (2.0 * np.pi)) * (
            u.B_SM[k] * L0 + n10 * d10[k] * LA + d2[k] * LB)
        axA.plot(np.log10(s), ds, color=cols[k], lw=2.2)
    axA.axhline(1.0 / A10, color=C_ZERO, ls="--", lw=1.4,
                label=r"reference reading  $g^2 = 1/25$  (a convention, not refitted)")
    axA.plot([np.log10(M_X)], [25.0], marker="*", ms=16, color=C_ROOT,
             mec="k", mew=0.8, ls="none",
             label=rf"meeting:  $M_X = {M_X:.1e}$ GeV")
    for Mv, lab in ((M_A, rf"$M_{{10}} = {M_A:.1e}$"),
                    (M_B, rf"$M_{{(1,2)}} = {M_B:.1e}$")):
        axA.axvline(np.log10(Mv), color="0.45", ls=":", lw=1.2)
        axA.text(np.log10(Mv) + 0.12, 55.5, lab, fontsize=10.5, color="0.35",
                 rotation=90, va="top")
    axA.axvline(np.log10(M_X), color=C_ROOT, ls="-", lw=1.0, alpha=0.6)
    axA.axvline(np.log10(MP), color="0.2", ls="--", lw=1.4)
    axA.text(np.log10(MP) - 0.15, 12.5, r"$M_{\rm GUT} = M_P$  abandoned",
             fontsize=10.5, color="0.25", rotation=90, ha="right", va="bottom")
    axA.set_xlabel(r"$\log_{10}(M/\mathrm{GeV})$")
    axA.set_ylabel(r"$1/\alpha_i$  (1-loop, piecewise $\beta$)")
    axA.set_xlim(np.log10(MZ) - 0.1, np.log10(MP) + 0.3)
    axA.set_ylim(5, 62)
    axA.set_title("(a)  desert content (solid) meets at $\\alpha_{10} \\equiv g^2 = 1/25$\n"
                  r"($\alpha_X = 0.0400$ standard); the closed clock reads $\pi/64 = 0.0491$ — $1.226\times$ above",
                  fontsize=11)
    axA.legend(fontsize=10, loc="lower right", framealpha=0.95)
    # (b) proton lifetime: the desert point vs the bound and the abandoned
    Ms = np.logspace(14, 19, 200)
    ltau = np.log10(sel["tau_p_yr"]) + 4.0 * (np.log10(Ms) - np.log10(M_X))
    axB.semilogy(Ms, 10.0 ** ltau, color=C_ROOT, lw=2.0,
                 label=r"$\tau_p \propto M_X^4/\alpha_X^2$")
    axB.axhline(2.4e34, color="k", lw=1.2, ls="--",
                label="Super-K bound  $2.4\\times10^{34}$ yr")
    axB.axhspan(1e35, 4e35, color=C_SU2, alpha=0.18)
    axB.text(1.1e14, 1.4e35, "Hyper-K reach", fontsize=10.5, color="0.45")
    axB.plot([M_X], [sel["tau_p_yr"]], marker="*", ms=16, color=C_ROOT,
             mec="k", mew=0.8, ls="none")
    axB.annotate(rf"CL10 desert:  $\tau_p = {sel['tau_p_yr']:.1e}$ yr",
                 (M_X, sel["tau_p_yr"]), xytext=(14, 8),
                 textcoords="offset points", fontsize=10.5, color=C_ROOT,
                 ha="left")
    axB.plot([MP], [4e47], marker="o", ms=7, mfc="none", mec="0.4", ls="none",
             label=r"$M_{\rm GUT} = M_P$  (abandoned)")
    axB.set_xscale("log")
    axB.set_xlabel(r"$M_X$  (GeV)")
    axB.set_ylabel(r"$\tau_p$  (yr)")
    axB.set_ylim(1e33, 3e48)
    axB.set_title(f"(b)  the desert point sits "
                  f"{sel['tau_p_yr'] / 2.4e34:.0f}$\\times$ above the bound — "
                  "F5 is live", fontsize=11.5)
    axB.legend(fontsize=10, loc="upper left", framealpha=0.95)
    fig.suptitle("The desert solution:  three $\\mathbf{10}\\oplus"
                 "\\overline{\\mathbf{10}}$ pairs + one $(1,2)$ pair unify at "
                 "the reference reading $\\alpha_{10} \\equiv g^2 = 1/25$ — a convention of the step, not the clock",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.875])
    save(fig, "fig20_desert_unification.png")

def _frag_label(k):
    """fragment key -> label with the Cartan charges as exact fractions"""
    from fractions import Fraction
    if k.startswith("(") and k.endswith(")"):
        a, b = (p.strip() for p in k[1:-1].split(","))
        fr = Fraction(float(b)).limit_denominator(12)
        bs = (str(abs(fr.numerator)) if fr.denominator == 1
                  else r"\frac{%d}{%d}" % (abs(fr.numerator), fr.denominator))
        return r"$f = (%s,\,%s%s)$" % (a, "-" if fr < 0 else "", bs)
    return r"$f = %s$" % k

def fig6_wilson():
    import json
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from cl10 import supercell as SC

    with open(exp("V65_wilson",
                           "results", "v65_wilson.json")) as f:
        v65 = json.load(f)["data"]
    gg = next(r for r in v65["part_B"].values()
              if r["su5"] == "Georgi-Glashow")
    root = next(r for r in gg["selective_roots"] if r["physical"])

    # sized so the full-width float stays under the page budget: three figure
    # floats share one page in the two-column layout
    fig = plt.figure(figsize=(11.8, 4.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.30], wspace=0.16)
    axA = fig.add_subplot(gs[0], projection="3d")
    fam_col = {"x": C_SU3, "y": C_SU2, "z": C_GUT}
    # the selective reading rings three of the twelve nodes, one per axis
    # family.  which three is an input (Section 8.1), so the ring marks one
    # signed representative per family -- an abs() membership test would match
    # all four members of every family and ring all 12
    split = {(np.pi / 2, np.pi / 2, 0.0), (np.pi / 2, 0.0, np.pi / 2),
             (0.0, np.pi / 2, np.pi / 2)}
    # the ringed nodes are drawn last: in projection they overlap other nodes
    # and their black ring must never be hidden behind one
    keyed = [(q in split, q) for q in SC.M_POINTS]
    for is_split, q in sorted(keyed, key=lambda k: k[0]):
        axis = next(j for j in range(3) if abs(q[j]) < 1e-9)
        col = fam_col[["x", "y", "z"][axis]]
        axA.scatter(*q, s=230 if is_split else 85, c=[col],
                    alpha=1.0 if is_split else 0.85, zorder=5 if is_split else 3,
                    edgecolors="black" if is_split else "white",
                    linewidths=1.9 if is_split else 0.5, depthshade=False)
    # the block zone is [-pi,pi]^3 and the twelve M points are the edge
    # centres of the inner half-zone cube: draw both, so no node can sit on
    # or outside a pane boundary
    h = np.pi / 2
    for a in (-h, h):
        for b in (-h, h):
            axA.plot([a, a], [b, b], [-h, h], color="0.55", lw=0.9, alpha=0.8)
            axA.plot([a, a], [-h, h], [b, b], color="0.55", lw=0.9, alpha=0.8)
            axA.plot([-h, h], [a, a], [b, b], color="0.55", lw=0.9, alpha=0.8)
    axA.set_xlim(-np.pi, np.pi)
    axA.set_ylim(-np.pi, np.pi)
    axA.set_zlim(-np.pi, np.pi)
    from matplotlib.ticker import FixedFormatter, FixedLocator
    for axis in (axA.xaxis, axA.yaxis, axA.zaxis):
        axis.set_major_locator(FixedLocator([-h, 0.0, h]))
        axis.set_major_formatter(FixedFormatter([r"$-\pi/2$", "0",
                                                 r"$\pi/2$"]))
        axis.set_tick_params(labelsize=10)
    axA.set_xlabel("$q_x$", fontsize=11, labelpad=1)
    axA.set_ylabel("$q_y$", fontsize=11, labelpad=1)
    axA.set_zlabel("$q_z$", fontsize=11, labelpad=1)
    axA.view_init(elev=16, azim=-58)
    axA.set_box_aspect((1, 1, 1))
    for axis in (axA.xaxis, axA.yaxis, axA.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("0.78")
        axis.pane.set_linewidth(0.8)
        axis.line.set_color("0.6")      # the box edges, not black rules
    axA.grid(False)
    axA.legend(
        handles=[plt.Line2D([], [], marker="o", ls="", ms=9, color=c,
                            label=f"{n} family")
                 for c, n in ((C_SU3, "$x$"), (C_SU2, "$y$"), (C_GUT, "$z$"))]
        + [plt.Line2D([], [], marker="o", ls="", ms=11, mfc="0.8",
                      mec="black", mew=1.8,
                      label="ringed: selective node (one per family)")],
        fontsize=9.5, ncol=4, loc="upper center", frameon=False,
        bbox_to_anchor=(0.5, -0.13))
    axA.set_title("(a)  12 M points: the edge centres of the half-zone\n"
                  r"cube $[-\pi/2,\pi/2]^3$ inside the zone $[-\pi,\pi]^3$"
                  "\n(which three split is an input, not derived)",
                  fontsize=11.5)

    axB = fig.add_subplot(gs[1])
    ths = root["thresholds_GeV"]
    labels = list(ths)
    vals = [ths[k] for k in labels]
    ypos = np.arange(len(labels))
    axB.barh(ypos, vals, color=C_ROOT, alpha=0.85, height=0.6)
    axB.set_yticks(ypos)
    axB.set_yticklabels([_frag_label(k) for k in labels], fontsize=10.5)
    axB.set_xscale("log")
    axB.set_xlim(1e6, 1e19)
    axB.axvline(root["M_X_GeV"], color=C_ZERO, lw=2.0)
    axB.axvline(float(_config.CFG.desert.m_x_gev), color="0.35", lw=1.6, ls="--")
    axB.text(root["M_X_GeV"] * 1.5, 2.5,
             f"Wilson root meeting\n$M_X$ = {root['M_X_GeV']:.1e} GeV\n"
             "(excluded, F5)", color=C_ZERO, fontsize=9.5, ha="left",
             va="center", rotation=90, linespacing=1.3)
    mx = float(_config.CFG.desert.m_x_gev)
    axB.text(mx * 1.5, 2.5,
             f"data-preferred\n$M_X = {mx:.1e}$ GeV",
             color="0.3", fontsize=9.5, ha="left", va="center", rotation=90,
             linespacing=1.3)
    axB.set_xlabel("fragment mass [GeV]")
    axB.set_title("(b)  the six thresholds of the 3 split nodes\n"
                  f"     ($\\theta$ = {root['theta_rad']:.3f}, "
                  f"$w$ = {root['w']:.2e}): scale excluded", fontsize=11.5)
    fig.suptitle("The internal Wilson factor: pattern predicted, "
                 "scale falsified", fontsize=13)
    # tight_layout cannot see the 3D axes: the geometry is set by hand so the
    # suptitle clears the panel titles and the family legend clears the labels
    fig.subplots_adjust(left=0.045, right=0.985, top=0.795, bottom=0.205,
                        wspace=0.24)
    save(fig, "fig6_wilson_mechanism.png")

def fig22_composite():
    import json
    import sys
    sys.path.insert(0, exp("V65_wilson"))
    sys.path.insert(0, exp("V81_link_desert_content"))
    from V65_wilson import cartan_weights, identify_XY, physical_su5
    from V81_link_desert_content import flux_weights, sm_roots, c2_labels
    from cl10 import unification as u

    with open(exp("V81_link_desert_content",
                           "results", "v81_link_desert_content.json")) as f:
        d = json.load(f)["data"]

    lam, _ = cartan_weights()
    xs, ys = identify_XY(lam)
    yc = ys[0]
    gg = next(x for x in xs
              if physical_su5(lam @ x, lam @ yc) == "Georgi-Glashow")
    al = flux_weights()
    comp = (lam[:, None, :] + al[None, :, :]).reshape(-1, 5)
    bL, b1, b2 = sm_roots(al, lam, gg, yc)
    bL2 = float(bL @ bL)
    Ys = comp @ yc
    m2 = (comp @ bL) / bL2
    c2 = c2_labels(comp, b1, b2)
    db_state = np.stack([(2 / 3) * (3 / 5) * Ys ** 2,
                         (2 / 3) * m2 ** 2,
                         (2 / 3) * c2 / 8.0], axis=1)
    zero_j = next(j for j in range(45)
                  if np.abs(al[j]).sum() < 1e-12)
    sm_idx = np.arange(16) * 45 + zero_j
    dbs = db_state.copy()
    dbs[sm_idx] *= 0.0
    db_ex = dbs.sum(axis=0)

    x = math.log(float(_config.CFG.desert.m_x_gev) / u.M_Z)
    INV = 1 / u.ALPHA_10
    TWO_PI = 2 * math.pi
    need = [TWO_PI * (1 / u.A_Z[i] - INV - u.B_SM[i] * x / TWO_PI) / db_ex[i]
            for i in range(3)]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    names = [r"$\mathrm{U}(1)$", r"$\mathrm{SU}(2)$", r"$\mathrm{SU}(3)$"]
    cols = [C_GUT, C_SU2, C_SU3]
    ypos = np.arange(3)[::-1]
    axA.barh(ypos, need, color=cols, alpha=0.85, height=0.55)
    for y, v in zip(ypos, need):
        axA.text(v + 0.02, y, f"{v:.2f}", va="center", fontsize=11.5)
    axA.set_yticks(ypos)
    axA.set_yticklabels(names, fontsize=12)
    axA.set_xlabel(r"required threshold sum $\sum_s (x - z_s)$ "
                   rf"[$\ln M$ units], $M_X = {float(_config.CFG.desert.m_x_gev):.1e}$ GeV")
    axA.set_title("(a)  the three meeting equations demand three\n"
                  "     different sums - no single spectrum fits",
                  fontsize=11.5)
    axA.axvspan(min(need), max(need), color=C_ZERO, alpha=0.08)
    axA.text(0.68, 0.90, f"spread factor {max(need)/min(need):.1f}",
             transform=axA.transAxes, fontsize=11, color=C_ZERO)

    with open(exp("V65_wilson", "results", "v65_wilson.json")) as f65:
        th_root = (json.load(f65)["data"]["part_B"]["X_1"]
                   ["selective_roots"][0]["theta_rad"])
    v = math.cos(th_root) * gg + math.sin(th_root) * yc
    c_al = np.round(np.abs(al @ v), 6)
    vals, cnts = np.unique(c_al, return_counts=True)
    axB.bar(vals, cnts, width=0.045, color=C_ROOT, alpha=0.85)
    axB.set_ylim(0, 28)
    axB.set_xlabel(r"flux mass pattern $|\alpha \cdot v|$ "
                   rf"($\theta$ = {th_root:.3f}, the root angle)")
    axB.set_ylabel("count")
    axB.set_title("(b)  the forced pattern $m = w(c_0 + |\\alpha\\cdot v|)$:\n"
                  f"     {d['meeting']['n_solutions']} physical solutions "
                  "over the full scan", fontsize=11.5)
    fig.suptitle("The gauged layer's composite content "
                 "($16 \\otimes 45 = 720$ formal tensor-product states; the "
                 "gauge-invariant kernel is empty and does not deliver "
                 "the meeting", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig22_composite_content.png")

def fig13_closed_loops():
    import json
    with open(exp("V203_closed_branch_coefficient", "results",
                  "v203_closed_branch_coefficient.json")) as f:
        cl = json.load(f)["data"]["closed"]
    Fn = {int(k): v for k, v in cl["F_n_2_30"].items()}
    loop = float(_config.CFG.physics.loop)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ev = sorted(n for n in Fn if n % 2 == 0)
    fourk = [n for n in ev if n % 4 == 0]
    fourk2 = [n for n in ev if n % 4]
    axA.semilogy(fourk, [Fn[n] for n in fourk], "o-", color=C_ZERO, lw=1.8,
                 ms=6, label=r"$n = 4k$: the measured branch (removable)")
    axA.semilogy(fourk2, [Fn[n] for n in fourk2], "o", mfc="none", mew=1.4,
                 color="0.55", ms=7,
                 label=r"$n = 4k+2$: $0/0$, no limit at the root clock")
    odd = [n for n in sorted(Fn) if n % 2]
    axA.semilogy(odd, [1e-16] * len(odd), "x", color=C_ROOT, ms=7,
                 label=r"odd lengths: $F_n = 0$ exactly (selection rule)")
    axA.axhline(loop, color="0.35", ls=":", lw=1.5)
    axA.text(2.1, loop * 1.4, rf"$1/(16\pi) = {loop:.4f}$", fontsize=11,
             color="0.3")
    axA.set_xlabel("closed-loop length $n$ (hops on $n$ distinct links)")
    axA.set_ylabel(r"Cartan-content ratio $F_n = \alpha((MG)^n)/\alpha(G^n)$")
    axA.set_ylim(1e-17, 200)
    axA.legend(fontsize=9.5, loc="center right")
    axA.set_title("(a)  the clock filter splits the even lengths: the $4k$ class\n"
                  rf"     carries the ladder ($F_4 = {Fn[4]:.2f}$, $F_{{12}} ="
                  rf" {Fn[12]:.2f}$); the $4k + 2$ points are floor ratios",
                  fontsize=11.5)

    labels = [r"$G$", r"$MG$" + "\n(root clock)",
              r"$MG$" + "\n(reference $1/5$)"]
    vals = [1.0, cl["rho_MG"], cl["rho_MG_implemented_V85"]]
    cols = ["0.6", C_ZERO, C_ROOT]
    axB.bar(range(3), vals, color=cols, alpha=0.85, width=0.6)
    for i, v in enumerate(vals):
        axB.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=11)
    axB.axhline(1.0, color="0.35", ls=":", lw=1.2)
    axB.set_xticks(range(3))
    axB.set_xticklabels(labels, fontsize=11)
    axB.set_ylim(0, 1.15)
    axB.set_ylabel("spectral radius")
    axB.set_title("(b)  the traced hop CONTRACTS the propagation:\n"
                  f"     $\\rho(MG) = {cl['rho_MG']:.6f}$ at the root clock",
                  fontsize=11.5)

    fig.suptitle("Closed traced loops at the root clock: the first crossing of "
                 rf"$1/(16\pi)$ is $n^* = {cl['first_crossing_n']}$", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig13_closed_loops.png")

def fig21_equal_shift():
    import json
    with open(exp("V84_equal_shift_theorem",
                           "results", "v84_equal_shift_theorem.json")) as f:
        d = json.load(f)["data"]
    r16 = d["over_16"]["sum_db1"] / d["over_16"]["sum_db2"]
    r45 = d["over_flux45"]["sum_db1"] / d["over_flux45"]["sum_db2"]
    rv10 = d["vector10"]["ratio_db1_over_db2"]
    v64 = d["V64_content"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    names = [r"matter $\mathbf{16}$", r"flux $\mathbf{45}$",
             r"vector $\mathbf{10}$", r"SU(5) $\mathbf{10}$",
             "(1,2) pair\n(fragment!)"]
    ratios = [r16, r45, rv10,
              v64["ratio_db1_over_db2_for_10"],
              v64["ratio_db1_over_db2_for_(1,2)"]]
    cols = [C_ROOT, C_ROOT, C_ROOT, C_ROOT, C_ZERO]
    axA.bar(range(5), ratios, color=cols, alpha=0.85, width=0.6)
    for i, v in enumerate(ratios):
        axA.text(i, v + 0.03, f"{v:.2f}", ha="center", fontsize=11.5)
    axA.axhline(1.0, color="0.35", ls=":", lw=1.4)
    axA.text(4.35, 1.18, "equal shift: $\\sum db_1 = \\sum db_2$", fontsize=11,
             color="0.3", ha="right")
    axA.set_xticks(range(5))
    axA.set_xticklabels(names, fontsize=10.5)
    axA.set_ylim(0, 1.3)
    axA.set_ylabel(r"ratio $\sum db_1 / \sum db_2$")
    axA.set_title("(a)  every COMPLETE multiplet shifts U(1) and SU(2)\n"
                  "     equally - only the fragment evades the theorem",
                  fontsize=11.5)

    parts = [v64["mismatch_from_10s"], v64["mismatch_from_(1,2)"]]
    axB.bar([0, 1], parts, color=[C_GUT, C_ZERO], alpha=0.85, width=0.5)
    axB.text(0, 0.01, "exactly 0", ha="center", fontsize=11.5)
    axB.text(1, parts[1] - 0.03, f"{parts[1]:.3f}", ha="center", fontsize=11.5,
             color=C_ZERO)
    axB.set_xticks([0, 1])
    axB.set_xticklabels([r"$3\times(\mathbf{10} \oplus \overline{\mathbf{10}})$",
                         r"$(1,2)$ pair"], fontsize=11)
    axB.axhline(0, color="0.3", lw=1)
    axB.set_ylim(-0.35, 0.12)
    axB.set_ylabel(r"mismatch $\sum db_1 - \sum db_2$")
    axB.set_title("(b)  the desert repair's whole mismatch is carried by\n"
                  "     its fragment pair - the automaton has none",
                  fontsize=11.5)

    fig.suptitle("The equal-shift theorem: the desert content must be "
                 "$\\mathrm{SU}(5)$-fragment content", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig21_equal_shift.png")
