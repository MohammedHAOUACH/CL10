"""Figure module: algebra."""
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from .style import ROOT, save
from .style import C_ALG, C_FIT, C_GUT, C_ROOT, C_SU2, C_SU3, C_ZERO


def fig1_one_object():
    from cl10 import qca_spectrum

    lv = np.asarray(qca_spectrum.levels(), dtype=float)          # 32 phases, sorted
    n = len(lv)
    zp = float(lv.sum())

    fig, ax = plt.subplots(figsize=(5.6, 4.6))
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
    fig.text(0.5, 0.02,
             r"$\sum_k \omega_k = $" + f"{zp:.1e}" + r"  (0-point, W1)      "
             + r"$\det U = $" + f"{1 - sm['det_U_minus_1']:.2e}"
             + r"  (SL(32), W1)" + "\n"
             + r"CPT mirror residual $= $" + f"{sm['mirror_resid']:.1e}"
             + r"      smallest level  $\omega_{\min} = $" + f"{sm['omega_min']:.4f}",
             ha="center", va="bottom", fontsize=10.5,
             bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#888", alpha=0.95))
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], color=C_ROOT, lw=2, label=r"$-\omega_k$  ($\overline{16}$, antimatter)"),
        Line2D([0], [0], color=C_ZERO, lw=2, label=r"$+\omega_k$  ($16$, matter)"),
    ], loc="lower right", fontsize=10.5, framealpha=0.95)
    ax.set_xlim(-1, n)
    fig.subplots_adjust(top=0.86, bottom=0.34, left=0.10, right=0.97)
    save(fig, "fig1_one_object_32spectrum.png")

def fig3_512_512():
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
                f"{v:+.0f}", ha="center", fontsize=12.5,
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
    save(fig, "fig3_512_512.png")

def fig6_anomaly():
    from cl10 import matter

    # the [U(1)]^3 per-field contribution n*Y^3 (the exact cancellation)
    order_keys = ["Q", "u^c", "d^c", "L", "e^c", "N"]
    order_lbl = ["Q", r"$u^c$", r"$d^c$", "L", r"$e^c$", "N"]
    frac_lbl = {  # the exact n*Y^3 per field (an exact algebraic identity)
        "Q": "1/36", "u^c": "-8/9", "d^c": "1/9",
        "L": "-1/4", "e^c": "1", "N": "0"}
    perfield = [n * y ** 3 for (_s3, _s2, y, n) in
                (matter.SM_FIELDS[k] for k in order_keys)]

    fig, ax = plt.subplots()
    x = np.arange(len(order_keys))
    cols = [C_SU3 if v >= 0 else C_GUT for v in perfield]
    ax.bar(x, perfield, color=cols, edgecolor="k", linewidth=0.5, width=0.5)
    ax.axhline(0, color="k", lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(order_lbl)
    ax.set_ylabel(r"$n \times Y^3$  (the $[\mathrm{U}(1)]^3$ contribution per field)")
    for xi, k in zip(x, order_keys):
        v = perfield[xi]
        va = "bottom" if v >= 0 else "top"
        ax.text(xi, v + (0.10 if v >= 0 else -0.10), frac_lbl[k],
                ha="center", va=va, fontsize=12.5)
    ax.set_ylim(-1.55, 1.55)
    ax.set_title(r"The $[\mathrm{U}(1)]^3$ anomaly: the six SM fields"
                 "\ncancel ($\\sum = 0$) — the 16 Weyl is anomaly-free")
    save(fig, "fig6_anomaly_free.png")

def fig9_cpt_mirror():
    from cl10 import qca_spectrum as qs

    p = qs.cpt_block_phases()
    m = np.asarray(p["ph16"], dtype=float)      # the 16 (matter) eigenphases, sorted
    b = np.asarray(p["ph16b"], dtype=float)     # the 16-bar (antimatter), for the check
    # the CPT mirror holds: the 16-bar multiset is the phase-conjugate of the 16
    assert np.allclose(np.sort(b), -np.sort(m)[::-1], atol=1e-9), "the 16-bar is the mirror of the 16"
    k = np.arange(1, len(m) + 1)
    fig, ax = plt.subplots()
    # the CPT axis (omega = 0) and the two mirror blocks: for each rank k,
    # the matter level m_k and its CPT conjugate -m_k (the 16-bar level)
    ax.axvline(0, color="k", lw=1.0, ls="--", label="CPT axis  ($\\omega = 0$)")
    ax.plot(m, k, "o-", color=C_SU3, lw=1.6, ms=5, label="matter  $16$   ($\\omega_k$)")
    ax.plot(-m, k, "o-", color=C_SU2, lw=1.6, ms=5,
            label="antimatter  $\\overline{16}$   ($-\\omega_k$)")
    ax.set_xlabel("eigenphase  $\\omega_k$  (rad/step)")
    ax.set_ylabel("rank $k$  (within the 16)")
    ax.set_title("The CPT spectral mirror: the $16$ (matter) and $\\overline{16}$ (antimatter) QCA blocks are\n"
                 r"exact phase-conjugates  ($\omega_{\overline{16}} = -\omega_{16}$) — CPT as spectral conjugation")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=3,
              framealpha=0.95, fontsize=10.5)
    fig.subplots_adjust(bottom=0.24)
    save(fig, "fig9_cpt_mirror.png")

def fig11_matter_content():

    items = [
        (r"$Q$",   6, r"$(3,2,+1/6)$",   "10", C_SU3),
        (r"$\bar u$",  3, r"$(\bar3,1,-2/3)$", "10", C_SU3),
        (r"$\bar e$",  1, r"$(1,1,+1)$",   "10", C_SU3),
        (r"$\bar d$",  3, r"$(\bar3,1,+1/3)$", "$\\bar5$", C_SU2),
        (r"$L$",   2, r"$(1,2,-1/2)$",  "$\\bar5$", C_SU2),
        (r"$\nu^c$", 1, r"$(1,1,0)$",   "$1$", C_ROOT),
    ]
    y = np.arange(len(items))[::-1]
    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    ax.barh(y, [it[1] for it in items], color=[it[4] for it in items],
            edgecolor="k", lw=0.5)
    for yi, (lab, dim, rep, su5, col) in zip(y, items):
        ax.text(dim + 0.15, yi, f"{rep}  (in the {su5})",
                va="center", fontsize=11.5, color="0.2")
    ax.set_yticks(y)
    ax.set_yticklabels([it[0] for it in items])
    ax.set_xlabel("dimension  $n$  (Weyl components)")
    ax.set_xlim(0, 10)
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], color=C_SU3, lw=0, marker="s", ms=11,
               label=r"the $\mathbf{10}$  :  $Q, \bar u, \bar e$   $(6+3+1 = 10)$"),
        Line2D([0], [0], color=C_SU2, lw=0, marker="s", ms=11,
               label=r"the $\mathbf{\bar 5}$  :  $\bar d, L$   $(3+2 = 5)$"),
        Line2D([0], [0], color=C_ROOT, lw=0, marker="s", ms=11,
               label=r"the $\mathbf{1}$  :  $\nu^c$   $(1)$  — the sterile"),
    ], loc="lower right", framealpha=0.95, fontsize=10.5)
    ax.set_title(r"The 16 Weyl $= $ one full SM generation $= \mathbf{10} \oplus \mathbf{\bar 5} \oplus \mathbf{1}$ of $\mathrm{SU}(5)$"
                 r"  —  the six left-handed fields, total dimension $16$ (exact, 0-parameter)",
                 fontsize=12)
    save(fig, "fig11_matter_content.png")

def fig32_graphical_abstract():
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(10.6, 6.0))
    ax.set_xlim(-5.95, 5.95)
    ax.set_ylim(-3.1, 3.1)
    ax.axis("off")
    ax.set_aspect("equal")

    # the root: Cl(10,0) with its ten generators (spokes never cross text)
    R = 1.05
    ax.add_patch(plt.Circle((0, 0), R, fc="#f3ecf7", ec=C_ROOT,
                            lw=2.4, zorder=3))
    for k in range(10):
        a = math.tau * k / 10 + math.tau / 20   # 18 deg + 36 deg k
        ax.plot([0, 0.97 * R * math.cos(a)], [0, 0.97 * R * math.sin(a)],
                color=C_ROOT, lw=2.0, zorder=4)
        c, s = math.cos(a), math.sin(a)
        ax.text(1.24 * R * c, 1.24 * R * s, f"$h_{{{k + 1}}}$", fontsize=11,
                color=C_ROOT, zorder=5,
                ha="left" if c > 0.3 else "right" if c < -0.3 else "center",
                va="bottom" if s > 0.3 else "top" if s < -0.3 else "center")
    ax.text(0, R + 0.34, "$\\mathrm{Cl}(10,0)$", ha="center", va="center",
            fontsize=15, weight="bold", color=C_ROOT, zorder=5)
    ax.text(0, -(R + 0.66), "$h_a^2 = I$,   $h_a h_b = -h_b h_a$",
            ha="center", va="center", fontsize=11.5, color="0.35", zorder=5)

    # five theorem boxes, compact, at pentagon positions
    HW, HH = 1.42, 0.55
    branches = [
        (-3.75, 2.30, "(1+3) Minkowski\ncausal front $c = 1.0$\nTheorems 1, 5"),
        (3.75, 2.30, "$\\mathbf{16}$ Weyl = ONE generation\n8 anomalies $=0$\nTheorems 2, 3"),
        (4.15, 0.0, "$\\mathbf{45}$ = Spin(10) adjoint\n$5/3$ normalisation\nTheorem 4"),
        (3.75, -2.30, "gravity on the lattice\n$M_P = 1.22\\times10^{19}$ GeV\nSection 9"),
        (-3.75, -2.30, "$\\mathbb{Z}_2$ grading\nstatistics from the root\nSection 7"),
    ]
    for (x, y, txt) in branches:
        ax.add_patch(FancyBboxPatch((x - HW, y - HH), 2 * HW, 2 * HH,
                                    boxstyle="round,pad=0.06",
                                    fc="#eef3f8", ec=C_ALG, lw=1.5,
                                    zorder=3))
        ax.text(x, y, txt, ha="center", va="center", fontsize=11.5,
                linespacing=1.35, zorder=4)
        # arrow: circle rim (outside the label ring) to the box rim
        r = math.hypot(x, y)
        lam = min(HW / max(abs(x), 1e-9), HH / max(abs(y), 1e-9))
        ax.annotate("", xy=(x * (1 - lam), y * (1 - lam)),
                    xytext=(1.42 * x / r, 1.42 * y / r),
                    arrowprops=dict(arrowstyle="-|>", color="0.35",
                                    lw=1.4, shrinkA=0, shrinkB=0), zorder=2)
    ax.set_title("One Clifford root, five theorems — every entry is "
                 "algebraic and 0-parameter", fontsize=14)
    fig.tight_layout(pad=0.4)
    save(fig, "fig32_graphical_abstract.png")

def fig34_brillouin_doublers():
    from cl10 import brillouin
    r = brillouin.brillouin()
    nodes = r["nodes"]
    fig = plt.figure(figsize=(6.2, 5.2))
    ax = fig.add_subplot(projection="3d")
    pi = math.pi
    for s, e in [((-pi, -pi, -pi), (pi, -pi, -pi)),
                 ((pi, -pi, -pi), (pi, pi, -pi)),
                 ((pi, pi, -pi), (-pi, pi, -pi)),
                 ((-pi, pi, -pi), (-pi, -pi, -pi)),
                 ((-pi, -pi, pi), (pi, -pi, pi)),
                 ((pi, -pi, pi), (pi, pi, pi)),
                 ((pi, pi, pi), (-pi, pi, pi)),
                 ((-pi, pi, pi), (-pi, -pi, pi)),
                 ((-pi, -pi, -pi), (-pi, -pi, pi)),
                 ((pi, -pi, -pi), (pi, -pi, pi)),
                 ((pi, pi, -pi), (pi, pi, pi)),
                 ((-pi, pi, -pi), (-pi, pi, pi))]:
        ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]],
                color="0.6", lw=0.8)
    for n in nodes:
        c = C_FIT if n["chirality"] > 0 else C_ALG
        ax.scatter(*n["q"], color=c, s=90, edgecolor="black", lw=0.6,
                   depthshade=False)
    ax.scatter(0, 0, 0, color="#d4a017", s=130, marker="D",
               edgecolor="black", lw=0.6, depthshade=False)
    ax.text(-0.2 * pi, -0.2 * pi, 0.78 * pi,
            f"$\\Gamma$ gap {r['gap_q0_16']:.4f}", fontsize=12)
    ax.set_xlim(-pi, pi)
    ax.set_ylim(-pi, pi)
    ax.set_zlim(-pi, pi)
    ax.set_xlabel("$q_x$")
    ax.set_ylabel("$q_y$")
    ax.set_zlabel("$q_z$")
    ax.set_title(f"{r['n_doubler_nodes']} gapless doubler nodes on the "
                 f"zone edges ($\\sum$ chirality = "
                 f"{r['chirality_sum']:.0f}) - the doubling theorem, "
                 "displayed, not hidden", fontsize=13)
    ax.view_init(elev=18, azim=-58)
    fig.tight_layout()
    save(fig, "fig34_brillouin_doublers.png")
