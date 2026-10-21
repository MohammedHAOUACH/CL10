"""Figure module: residuals."""
import math

import matplotlib.pyplot as plt
import numpy as np

from .style import exp, save
from .style import C_ALG, C_CONV, C_FIT, C_GUT, C_OPEN, C_ROOT, C_SU2, C_SU3, C_ZERO
from cl10 import config as _config      # the shared inputs (no engine import here)


def fig10_gauged_layer():
    import json
    with open(exp("V79_gauged_bond_layer",
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
    fig.suptitle("The gauged bond layer: the Gauss constraint closes "
                 "on the root module", fontsize=13.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig10_gauged_bond_layer.png")

def fig15_w7_portal():
    import json
    with open(exp("V88_mirror_scalar",
                           "results", "v88_mirror_scalar.json")) as f:
        d = json.load(f)["data"]
    lam2, sig2 = d["V2_lam_star"], d["V2_sigma_SI_cm2"]
    lam1, sig1 = d["V1_lam_for_Omega"], d["V1_sigma_SI_cm2"]
    lz = float(_config.CFG.bounds.lz_sigma_si_cm2)
    gap = d["ss_hh_gap_GeV"]
    cloud = d["V1_virtual_Pmax"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.25, 1.0]})

    lam = np.logspace(-5, 0, 400)
    f_N = float(_config.CFG.dm.f_n_nucleon)
    m_N = float(_config.CFG.physics.proton_mass_gev)
    M_H = float(_config.CFG.sm.higgs_mass_gev)
    M_DM = (float(_config.CFG.physics.m_planck_gev)
            * float(_config.CFG.physics.loop) ** 10)
    mu = M_DM * m_N / (M_DM + m_N)
    sig = (lam ** 2 * f_N ** 2 * mu ** 2 / (np.pi * M_H ** 4)
           * float(_config.CFG.units.gev2_to_cm2))
    axA.loglog(lam, sig, color=C_ROOT, lw=2.2,
               label=r"$\sigma_{\rm SI} = \lambda^2 f_N^2 \mu^2 / "
                     r"(\pi m_h^4)$")
    axA.axhline(lz, color=C_ZERO, ls="--", lw=1.8,
                label=r"LZ bound (4.2 t yr): $2.2\times10^{-48}$ min "
                      rf"(40 GeV), $\approx{lz:.1e}$ at {M_DM:.0f} GeV")
    axA.fill_between([1e-5, 1], 1e-58, lz, color=C_ZERO, alpha=0.07)
    axA.plot([lam2], [sig2], "o", ms=11, color=C_SU2, zorder=5)
    axA.annotate(f"pair portal (Z2 BROKEN)\n$\\lambda_* = {lam2:.1e}$: "
                 f"$\\Omega h^2 = {float(_config.CFG.physics.omega_dm_h2):.3f}$\n$\\sigma_{{\\rm SI}} = "
                 f"{sig2:.0e}$ - below reach",
                 xy=(lam2, sig2), xytext=(3e-4, 8e-43), fontsize=11,
                 color=C_SU2, va="top",
                 arrowprops=dict(arrowstyle="->", color=C_SU2, lw=1.2))
    axA.plot([lam1], [sig1], "X", ms=13, color=C_ZERO, zorder=5)
    axA.annotate(f"Higgs portal (Z2 EXACT)\n$\\lambda = {lam1:.2f}$: "
                 f"$\\Omega h^2 = {float(_config.CFG.physics.omega_dm_h2):.3f}$\n$\\sigma_{{\\rm SI}} = "
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
        ("mirror relic $\\overline{Q}$", "stable under both portals",
         "#2c7a3f"),
    ]
    axB.axis("off")
    axB.set_title("(b)  the audit verdicts", fontsize=12.5)
    y = 0.98
    for (k, v, c) in rows:
        axB.text(0.02, y, k, fontsize=10.5, weight="bold", va="top")
        axB.text(0.52, y, v, fontsize=10.5, color=c, va="top")
        y -= (max(k.count("\n"), v.count("\n")) + 1.35) * 0.068
    axB.text(0.02, 0.01,
             "O7: the built portal carrier is exactly stable, but the observed\n"
             "abundance and the exact mirror $\\mathbb{Z}_2$ cannot both hold (F7);\n"
             "no state of the spinor module carries the 118 GeV scale.",
             fontsize=10.5, color=C_ROOT, weight="bold", va="bottom")

    fig.suptitle("The built dark-matter fix, audited against its own "
                 "symmetry", fontsize=13.5)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.34)
    save(fig, "fig15_w7_portal.png")

def fig14_loop_selection():
    import json
    with open(exp("V203_closed_branch_coefficient",
                  "results", "v203_closed_branch_coefficient.json")) as f:
        cl = json.load(f)["data"]["closed"]
    Fn = {int(k): v for k, v in cl["F_n_2_30"].items()}
    Fev = {int(k): v for k, v in cl["F_even_2_260"].items()}
    loop = float(_config.CFG.physics.loop)
    nstar = cl["first_crossing_n"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.3, 1.0]})

    ev = sorted(n for n in Fn if n % 2 == 0)
    odd = sorted(n for n in Fn if n % 2)
    axA.semilogy(ev, [Fn[n] for n in ev], "o-", color=C_GUT, ms=5, lw=1.6,
                 label=r"even lengths (root clock $\pi/4$)")
    axA.semilogy(odd, [1e-16] * len(odd), "x", color=C_ROOT, ms=7,
                 label=r"odd lengths: identically zero (selection rule)")
    axA.axhline(loop, color=C_ZERO, ls="--", lw=1.8,
                label=rf"loop factor $1/(16\pi) = {loop:.4f}$")
    axA.set_xlabel("closed-loop length $n$")
    axA.set_ylabel(r"Cartan-content ratio $F_n$")
    axA.set_ylim(1e-17, 200)
    axA.set_title("(a)  the exact odd-$n$ selection rule and the even ratios\n"
                  "     at the root clock, measured to $n = 30$",
                  fontsize=12.5)
    axA.legend(loc="lower left", fontsize=9.5, framealpha=0.95)

    nn = sorted(Fev)
    axB.loglog(nn, [Fev[n] for n in nn], "-", color=C_GUT, lw=1.5)
    axB.axhline(loop, color=C_ZERO, ls="--", lw=1.8)
    axB.plot([nstar], [Fev[nstar]], "*", ms=18, color=C_ZERO, zorder=6)
    axB.annotate(f"first crossing\n$n^* = {nstar}$, $F = {Fev[nstar]:.5f}$\n"
                 f"deficit $\\times{cl['deficit_vs_loop']:.4f}$"
                 if "deficit_vs_loop" in cl else
                 f"first crossing\n$n^* = {nstar}$",
                 xy=(nstar, Fev[nstar]), xytext=(40, 0.006), fontsize=10.5,
                 color=C_ZERO,
                 arrowprops=dict(arrowstyle="->", color=C_ZERO, lw=1.2))
    axB.set_xlabel("closed-loop length $n$ (even lengths, log scale)")
    axB.set_ylabel(r"even-length ratio $F_n$")
    axB.set_title("(b)  the even ratios to the crossing:\n"
                  rf"     the ladder length at the root clock IS $n^* = {nstar}$",
                  fontsize=12.5)

    fig.text(0.03, 0.015,
             f"odd lengths cancel identically (below $4\\times10^{{-15}}$, "
             f"exact selection rule of the root clock); the even ratios stay "
             f"above $1/(16\\pi)$ until $n^* = {nstar}$,\n"
             f"where the closed-loop suppression first crosses the loop "
             f"factor - rung 208 is the ladder length, and the "
             f"$n = 10$ reading of the reference angle is a compensation, "
             f"not the crossing.",
             fontsize=10.5, color="0.25", va="bottom")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.24)
    save(fig, "fig14_loop_selection.png")

def fig16_w2_protection():
    import json
    with open(exp("V90_w2_protection",
                           "results", "v90_w2_protection.json")) as f:
        d90 = json.load(f)["data"]
    with open(exp("V103_w2_democratic",
                           "results", "v103_w2_democratic.json")) as f:
        d103 = json.load(f)["data"]
    S17 = (float(_config.CFG.physics.m_planck_gev)
           * float(_config.CFG.physics.loop) ** 17 * 1e9)   # eV, the ladder S(17)
    bound = float(_config.CFG.neutrino.sum_mnu_max_ev)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                   gridspec_kw={"width_ratios": [1.3, 1.0]})

    r = np.linspace(0.05, 3.0, 4001)
    with np.errstate(divide="ignore"):
        s_id = S17 / (1 + r) + 2 * S17 / np.abs(2 - r)
    s_demo = S17 / (1 + r)
    axA.semilogy(r, s_id, color=C_GUT, lw=2.0,
                 label=r"identity $m_D$: $\frac{S_{17}}{1+r} + "
                       r"\frac{2S_{17}}{|2-r|}$")
    axA.semilogy(r, s_demo, color="#2c7a3f", lw=2.4,
                 label=r"democratic $m_D$: $m_3 = "
                       r"\frac{S_{17}}{1+r}$, two massless for every $r$")
    axA.axhline(bound, color=C_ZERO, ls="--", lw=1.8,
                label=rf"cosmological bound $\sum m_\nu < {bound}$ eV")
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
    axA.set_title("(a)  the identity $m_D$ scan (divergent, one point);\n"
                  "     the democratic $m_D$ gives a finite window",
                  fontsize=12.5)
    axA.legend(loc="upper left", fontsize=9.5)

    axB.axis("off")
    axB.set_title("(b)  what protects what", fontsize=12.5)
    rows = [
        ("cubic action", "all of $S_3$ (6 permutations)", C_GUT),
        ("commutant", f"dim {d90['commutant_nullity']} = "
         r"$\operatorname{span}\{I, J\}$ only", C_ZERO),
        ("$m_D$", "democratic, rank 1 - 0-parameter", "#2c7a3f"),
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
             "O2 re-audited: the measure-zero tuning was an artifact of\n"
             "the identity-$m_D$ assumption. The structure (2 massless +\n"
             "1 massive, NO) is protected by the same family symmetry\n"
             "that derives $\\sin^2\\theta_{12} = 1/3$; the value sits in a\n"
             "finite window. Solar mass still needs O6.",
             fontsize=10, color=C_ROOT, weight="bold", va="bottom")

    fig.suptitle("The democratic neutrino point: the identity-$m_D$ tuning is an "
                 "artifact of the wrong Dirac assumption",
                 fontsize=13.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig16_w2_protection.png")

def fig24_w3_forcing():
    import json
    with open(exp("V91_w3_forcing",
                           "results", "v91_w3_forcing.json")) as f:
        d = json.load(f)["data"]
    mx = d["fixed_meeting_scale"]
    obs = float(_config.CFG.sm.m_unification_obs_gev)

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
    axA.annotate("the 3-node triangle: realizable,\n"
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
                label=r"data-preferred meeting scale $2\times10^{16}$ GeV (F5)")
    for n, v in zip(ns, vals):
        axB.text(str(n), v * 1.15, f"{v:.1e}", ha="center", fontsize=11)
    axB.set_yscale("log")
    axB.set_ylim(1e14, 1e18)
    axB.set_xlabel("split-set size $n$ (cubic-symmetric: only 12)")
    axB.set_ylabel(r"fixed meeting $M_X(n)$ [GeV]")
    axB.set_title("(b)  every size misses the data-preferred scale;\n"
                  "     the rank-1 rigidity holds for $n = 4, 8, 12$",
                  fontsize=12.5)
    axB.legend(loc="upper left", fontsize=10.5)
    fig.text(0.75, 0.015,
             "O3: cubic-symmetric split (12 nodes) EXCLUDED; every proper\n"
             "split BREAKS the cubic symmetry - the count 3 is an input.",
             fontsize=11, color=C_ROOT, weight="bold", ha="center",
             va="bottom")

    fig.suptitle("Is the family count 3 forced? No - the symmetric "
                 "alternative is excluded", fontsize=13.5)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig24_w3_forcing.png")

def fig23_w4_cone():
    import json
    with open(exp("V92_w4_feasibility",
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
    ax.annotate("data-preferred meeting:\n"
                f"$M_X=1.8\\times10^{{16}}$ GeV, $\\alpha_X=1/25$ -\n"
                "INSIDE the hull.\n"
                f"witness: {lam}",
                tp2, xytext=(-2.0, -1.55),
                fontsize=12, color=C_FIT, ha="left",
                arrowprops=dict(arrowstyle="->", color=C_FIT))
    ax.set_title("The desert content the data require lies INSIDE "
                 "the hull of the node-$\\mathbf{16}$ fragments\n"
                 "($\\beta$-contribution plane; feasibility $\\Rightarrow$ "
                 "O4 is a realizability problem, not a no-go)",
                 fontsize=13.5)
    ax.set_xlabel("$\\beta_1,\\beta_2,\\beta_3$ plane (PCA)")
    ax.grid(True, alpha=0.3, ls=":")
    ax.set_xlim(v2[:, 0].min() - 0.35, v2[:, 0].max() + 1.85)
    ax.set_ylim(v2[:, 1].min() - 1.45, v2[:, 1].max() + 0.30)
    fig.tight_layout()
    save(fig, "fig23_w4_cone.png")

def fig11_vertex_ground_state():
    import json
    p100 = exp("V100_vertex_correction",
                        "results", "v100_vertex_correction.json")
    p101 = exp("V101_bond_current_ground_state",
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
    ax1.bar(range(5), vals, color=cols, width=0.62)
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
                  "kernel, $S_{\\rm GS} = \\ln 45$", fontsize=13)
    fig.suptitle("The gauged area law corrected by vertices, and the "
                 "ground state of the gauged layer",
                 fontsize=13.5)
    fig.tight_layout()
    save(fig, "fig11_vertex_ground_state.png")


def fig25_residual_attacks():
    """Attacking the residuals head-on: V109 (O5 CP/mirror), V110 (O4 bond
    off-diagonal), V111 (O6 generation dressing), V112 (O7 dark U(1)),
    V113 (O8 axis factorization), V116 (O3 triangle orbit-8)."""
    import json

    def _load(v):
        p = exp(v, "results",
                         v.lower() + ".json")
        return json.load(open(p))["data"]

    d109 = _load("V109_w5_cp_search")
    d110 = _load("V110_w4_offdiagonal_bond")
    d113 = _load("V113_w8_axis_factorization")

    fig, axs = plt.subplots(2, 3, figsize=(14.4, 8.8))
    (ax1, ax2, ax3), (ax4, ax5, ax6) = axs

    # (a) V109: the mirror is the CP-forbidding structure
    cases = d109["per_case_worst"]
    names = ["Wilson", r"$\kappa_s$ anti", r"$\kappa_s$ sym",
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
    ax1.set_title("(a) O5: the mirror forbids CP — only Cartan-type\n"
                  "extensions preserve it; class-dependent dressings\n"
                  "break Theorem 4", fontsize=12)

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
    ax2.set_title("(b) O4: one Hermitian bond operator with 12\n"
                  "parameters does not close the meeting",
                  fontsize=12)

    # (c) V111: named-parameter cost of reaching the measured mixing
    labels = ["scalar dressing\n(exactly blind)",
              "common family\ndressing",
              "chiral-asymmetric\ndressing"]
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
    ax3.set_title("(c) O6: generation-space dressing moves the\n"
                  "mixing — but every entry is a named input",
                  fontsize=12)

    # (d) V112: the dark-U(1) window squeeze
    e2 = 4 * math.pi / float(_config.CFG.sm.alpha_em_inv)
    mu_N = float(_config.CFG.physics.proton_mass_gev)
    ms = (float(_config.CFG.physics.m_planck_gev)
          * float(_config.CFG.physics.loop) ** 10)
    target_sv = (float(_config.CFG.dm.sigma_relic_target_cm3s)
                 / (float(_config.CFG.units.gev2_to_cm2)
                    * float(_config.CFG.units.c_cm_s)))
    dd_lim = (float(_config.CFG.bounds.dd_sigma_si_cm2)
              / float(_config.CFG.units.gev2_to_cm2))
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
    ax4.set_title("(d) O7: dark U(1) portal — relic needs\n"
                  r"$g_D\varepsilon \gtrsim 0.02$, DD allows "
                  r"$\lesssim 10^{-3}$:" "\nwindow EMPTY",
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
    ax5.set_title("(e) O8: node contributions EXACTLY equal across\n"
                  "families (0.00e+00); the coefficient 3 IS the choice\n"
                  "one-per-family — 99.4% agreement", fontsize=12)

    # (f) V116: the triangle and its orbit-8
    # plot the 12 nodes projected: (qx,qy) plane shows 8; mark families
    h = math.pi / 2
    plane = [(a * h, b * h) for a in (1, -1) for b in (1, -1)]
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
    ax6.set_title("(f) O3: 64 triangles → 4 cubic orbits "
                  "(8, 8, 24, 24);\n"
                  r"each choice breaks $O_h \to S_3$ — 8 distinct "
                  "vacua", fontsize=12)

    fig.suptitle("Attacking the residuals head-on: what the named tests close, "
                 "what they name", fontsize=14.5)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    save(fig, "fig25_residual_attacks.png")


def fig7_doubler_cones():
    """V102: the twelve doubler cones are exactly linear but never round."""
    import json
    d = json.load(open(exp("V102_lorentz_isotropy", "results",
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
                 "(grey, speed 1) - never round", fontsize=12)

    # (b) all 12 nodes: the 1200 measured velocities
    ax2 = fig.add_subplot(gs[1])
    for i, n in enumerate(nodes):
        vv = np.array(n["v_dirs"])[:, 3]
        ax2.scatter(np.full_like(vv, i) + 0.06 * ((i % 2) - 0.5), vv,
                    s=10, color=C_GUT, alpha=0.55)
    vp = np.array(nodes[0]["v_principal"])
    ax2.axhspan(float(vp[0]), float(vp[2]), color=C_GUT, alpha=0.15)
    ax2.annotate(f"tensor fit: principal band\n{vp[0]:.2f} / {vp[1]:.2f} / {vp[2]:.2f}",
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
    ax2.set_title("(b) 12 nodes x 100 directions: "
                  f"{d['n_linear']}/{d['n_linear'] + d['n_quadratic'] + d['n_mixed']} linear;\n"
                  "the velocity tensor fit gives principal cone\n"
                  f"{vp[0]:.2f}/{vp[1]:.2f}/{vp[2]:.2f} - a factor ~{1.0 / float(vp[2]):.0f}"
                  " below causal",
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
                  f"anisotropic cone: {min(an):.3f}-{max(an):.3f} - the lattice\n"
                  "signature is universal, not a defect",
                  fontsize=12)

    fig.suptitle("The twelve doubler cones: exactly linear, never Lorentz-"
                 "round", fontsize=14.5)
    fig.subplots_adjust(left=0.03, right=0.99, bottom=0.12, top=0.79)
    save(fig, "fig7_doubler_cones.png", tight=False)


def fig8_isotropy_search():
    """V115: no isotropy point on the whole V65 Wilson family."""
    import json
    d = json.load(open(exp("V115_w9_cone_isotropy", "results",
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
                  "the isotropy band is empty", fontsize=12)

    # (b) the distribution of A: pinned at full anisotropy
    ax2.hist(A, bins=np.arange(0.996, 1.0001, 0.0005), color=C_FIT)
    ax2.axvline(amin, color=C_ZERO, ls="--", lw=1.5)
    ax2.text(amin + 0.0006, ax2.get_ylim()[1] * 0.45,
             f"min {amin:.3f}", ha="left", va="center", fontsize=11,
             color=C_ZERO)
    ax2.set_xlabel(r"cone anisotropy $A$")
    ax2.set_ylabel("gapless points")
    ax2.set_xlim(0.9955, 1.0005)
    ax2.set_title(f"(b) every point sits at $A \\gtrsim {amin:.3f}$;\n"
                  "isotropy ($A\\to 0$) is not approached\n"
                  "anywhere on the family", fontsize=12)

    # (c) no angle helps: A vs theta
    ax3.scatter(th, A, c=w, cmap="viridis", s=9, alpha=0.8)
    ax3.axhline(float(d["bare_cone_max_anisotropy"]), color=C_OPEN,
                ls=":", lw=1.5)
    ax3.text(6.15, float(d["bare_cone_max_anisotropy"]) + 0.004,
             "bare cone (w = 0)", ha="right", fontsize=10.5, color=C_OPEN)
    ax3.set_xlabel(r"Wilson angle $\theta$ (rad)")
    ax3.set_ylabel(r"cone anisotropy $A$")
    ax3.set_ylim(0.9955, 1.0025)
    ax3.set_title("(c) all 60 angles x all strengths: the Wilson\n"
                  "factor moves the zeros, never the cone's\n"
                  "shape - O9 stays open", fontsize=12)

    fig.suptitle("Is continuous Lorentz invariance emergent anywhere on "
                 "the Wilson family? No", fontsize=14.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig8_isotropy_search.png")


def fig9_tangent_algebra():
    """V118/V120: the massless corner node is three nodal planes, not a cone."""
    import json
    with open(exp("V118_w9_little_group",
                           "results", "v118_w9_little_group.json")) as f:
        d118 = json.load(f)["data"]
    with open(exp("V120_w9_entangling_bond",
                           "results", "v120_w9_entangling_bond.json")) as f:
        d120 = json.load(f)["data"]
    t = d118["T4"][0]
    vecs = np.array(t["branch_velocity_vectors"])
    cols = [C_SU3, C_FIT, C_GUT]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.6, 4.6))

    # (a) the node in the kz = 0 slice: three velocity vectors, three nodal
    # lines. the three velocity labels go in the empty right-hand corner: at
    # the arrow heads they sit on top of one another
    ax1.add_patch(plt.Circle((0, 0), 1.0, fill=False, color="0.6", lw=0.8))
    for a, v in enumerate(vecs):
        nrm = math.hypot(v[0], v[1])
        ax1.annotate("", xy=(v[0], v[1]), xytext=(0, 0),
                     arrowprops=dict(arrowstyle="-|>", color=cols[a], lw=1.8))
        if nrm > 1e-9:
            d = np.array([-v[1], v[0]]) / nrm
            ax1.plot([-d[0], d[0]], [-d[1], d[1]], color=cols[a], lw=1.4,
                     ls="--")
        ax1.text(2.02, -0.30 - 0.32 * a,
                 r"$v_%d=(%.2f,\,%.2f)$" % (a + 1, v[0], v[1]),
                 fontsize=10, color=cols[a], ha="right", va="center")
    ax1.text(2.02, 1.12, r"dashed: nodal line $v_\alpha\cdot q=0$",
             fontsize=9.5, color="0.35", ha="right", va="center")
    ax1.set_xlim(-1.3, 2.05)
    ax1.set_ylim(-1.25, 1.25)
    ax1.set_aspect("equal")
    ax1.set_xlabel(r"$k_x$ (lattice units)")
    ax1.set_ylabel(r"$k_y$")
    ax1.set_title("(a) the gapless corner triplet, $k_z=0$ slice:\n"
                  "three planes through the node, not a cone", fontsize=11.5)

    # (b) the tangent velocity algebra: commuting, not Clifford-type
    comm = t["commutator_norms"]
    ax2.bar(range(3), comm, color=C_OPEN)
    ax2.axhline(t["sum_Sj2_scalar"], color=C_ALG, ls="--", lw=1.6)
    ax2.text(2.45, t["sum_Sj2_scalar"] * 0.5,
             r"Clifford-type scale $\sum_j\|S_j\|^2 = %.2f$"
             % t["sum_Sj2_scalar"], ha="right", va="top", fontsize=10,
             color=C_ALG)
    ax2.text(0.03, 0.34,
             "$\\sum_j S_j^2$ is a scalar;\neach branch planar to %.1e"
             % t["planar_branch_err"], transform=ax2.transAxes, fontsize=10,
             ha="left", va="bottom")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-14, 3.0)
    ax2.set_xticks(range(3))
    ax2.set_xticklabels([r"$[S_1,S_2]$", r"$[S_2,S_3]$", r"$[S_3,S_1]$"])
    ax2.set_ylabel(r"$\|[S_i,S_j]\|$ (lattice units)")
    ax2.set_title("(b) the tangent velocity operators COMMUTE:\n"
                  "max $\\|[S_i,S_j]\\|$ = %.1e, %.0f orders below the cone scale"
                  % (max(comm), -math.log10(max(comm) / t["sum_Sj2_scalar"])),
                  fontsize=11.5)

    # (c) the entangling bond layer: the right direction, nowhere near enough
    rows = d120["scan"]
    ph = [r["phi"] for r in rows]
    pl = [max(r["planar_err"], 1e-16) for r in rows]
    cm = [r["commutator_over_norm"] for r in rows]
    ax3.semilogy(ph, pl, "s-", color=C_FIT,
                 label=r"planarity defect $\max_{d,\alpha}|s_\alpha(d)-|v_\alpha\cdot d||$")
    ax3.axhline(1.0, color=C_SU3, ls=":", lw=1.4)
    ax3.text(0.02, 0.45, "Clifford-type target", fontsize=9.5, color=C_SU3,
             ha="left", va="top")
    ax3.set_ylim(1e-16, 3)
    ax3.set_xlabel(r"entangling bond phase $\phi$")
    ax3.set_ylabel("planarity defect (log)")
    ax3b = ax3.twinx()
    ax3b.plot(ph, cm, "o-", color=C_OPEN,
              label=r"$\max\|[S_i,S_j]\|/\|S\|$ (a cone needs $\sim1$)")
    ax3b.set_ylim(0, 0.3)
    ax3b.set_ylabel(r"$\max\|[S_i,S_j]\|/\|S\|$", color=C_OPEN)
    ax3b.tick_params(axis="y", colors=C_OPEN)
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3b.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, fontsize=9, loc="center right",
               framealpha=0.95)
    ax3.set_title("(c) an entangling bond layer: the right mechanism,\n"
                  "not enough of it (defect %.0e, commutator %.2f)"
                  % (max(pl), max(cm)), fontsize=11.5)

    fig.suptitle("Why the doubler node is never a Weyl cone: the tangent "
                 "velocities commute", fontsize=14.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig9_tangent_algebra.png")
