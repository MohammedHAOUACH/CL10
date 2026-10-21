"""Figure module: algebra."""
import math

import matplotlib.pyplot as plt
import numpy as np

from .style import save
from .style import C_ALG, C_FIT, C_GUT, C_ROOT, C_ZERO


def fig3_one_object():
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
                 r"$H=-i\log U$,  $U=\prod_{a<b}(\cos g_{10}\,I-i\sin g_{10}\,T_{ab})$ at the root clock $g_{10}=\pi/4$")
    # the 9-to-1 invariants, computed
    from cl10 import parameter_audit
    sm = parameter_audit.spectral_mirror()
    fig.text(0.5, 0.02,
             r"$\sum_k \omega_k = $" + f"{zp:.1e}" + r"  (0-point, O1)      "
             + r"$\det U = $" + f"{1 - sm['det_U_minus_1']:.2e}"
             + r"  (SL(32), O1)" + "\n"
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
    save(fig, "fig3_one_object_32spectrum.png")

def fig5_brillouin_doublers():
    from cl10 import supercell as SC
    _, _, U16 = SC.sectors()
    gap0 = SC.gap(SC.GAMMA, U16)
    nodes = []
    for q in SC.M_POINTS:
        w = np.angle(np.linalg.eigvals(SC.bloch(tuple(q), U16)))
        pinned = (int(np.sum(np.abs(w) < 1e-9))
                  + int(np.sum(np.abs(np.abs(w) - math.pi) < 1e-9)))
        nodes.append((q, pinned))
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
    t = np.linspace(-pi / 2, pi / 2, 25)
    for zero_axis in range(3):
        for s1 in (1.0, -1.0):
            a, b = [j for j in range(3) if j != zero_axis]
            pts = []
            for u in t:
                q = [0.0, 0.0, 0.0]
                q[a], q[b] = s1 * u, u
                pts.append(q)
            ax.plot([p[0] for p in pts], [p[1] for p in pts],
                    [p[2] for p in pts], color=C_ALG, lw=1.6, alpha=0.85)
    for q, pinned in nodes:
        assert pinned == 64
        ax.scatter(*q, color=C_FIT, s=90, edgecolor="black", lw=0.6,
                   depthshade=False)
    for q in SC.R_POINTS:
        ax.scatter(*q, color="none", s=150, marker="^", edgecolor=C_GUT,
                   lw=1.2, depthshade=False)
    ax.scatter(0, 0, 0, color="#d4a017", s=130, marker="D",
               edgecolor="black", lw=0.6, depthshade=False)
    ax.text(-0.2 * pi, -0.2 * pi, 0.78 * pi,
            f"$\\Gamma$ pinned (gap {gap0:.1e})", fontsize=12)
    ax.set_xlim(-pi, pi)
    ax.set_ylim(-pi, pi)
    ax.set_zlim(-pi, pi)
    ax.set_xlabel("$q_x$")
    ax.set_ylabel("$q_y$")
    ax.set_zlabel("$q_z$")
    ax.set_title("At the root clock the exactly gapless set is $\\Gamma$ plus the six\n"
                 "face-diagonal lines $\\{(\\pm t,\\pm t,0)\\}$ and axis permutations: the 12\n"
                 "M endpoints are pinned to $\\{0,\\pi\\}$ (64 of 128 levels, gap exactly 0),\n"
                 "the 8 R corners carry exactly 3-fold multiplets", fontsize=13)
    ax.view_init(elev=18, azim=-58)
    fig.tight_layout()
    save(fig, "fig5_brillouin_doublers.png")
