"""Figure module: dynamics."""

import matplotlib.pyplot as plt
import numpy as np

from .style import save
from .style import C_ZERO


def fig1_young():
    """The Young double-slit as a single-electron measurement
    (the Tonomura flow).  The emitter fires one localized electron at a time;
    each one diffracts through the slits and is detected at a *single point*
    sampled from the **time-integrated** on-screen |psi|^2 (the interference
    pattern the screen records -- a single instant of the traveling wave is
    only a blob, so we integrate over the passage).  The points are collected
    one per electron and build the fringes.  Panels: (a) the wave (x-z slice,
    the two diffracted waves interfering), (b) the time-integrated interference
    pattern + its x-profile (the fringe groups), (c) the 10-to-1 collapse of one
    electron onto a point, (d) the reality -- the accumulated single-electron
    detections cluster into the fringe groups."""
    import matplotlib.patches as mpatches
    np.random.seed(2026)
    from simulation.young import collapse as Coll
    from simulation.young import wave

    C = plt.get_cmap("inferno")
    TXT = "#2b2b2b"
    LIGHT = "#f0f0f0"

    def blackize(ax):
        ax.set_facecolor("black")
        for sp in ax.spines.values():
            sp.set_color("#666")
        ax.tick_params(labelcolor=TXT)

    def norm(A):
        A = np.asarray(A, float)
        return (A - A.min()) / max(A.max() - A.min(), 1e-12)

    Lx, Ly, Lz = 128, 128, 220
    lam, d, a = 4.5, 22.0, 5.0
    cfg = dict(grid=[Lx, Ly, Lz], r=0.5, wavelength=lam, source_width=18.0,
               source_amplitude=1.0, source_z=20, ramp_steps=200,
               source_mode="single", source_shape="point",
               source_y_width=14.0, packet_z_width=8.0,
               sponge_width=16, barrier_z=44, slit_separation=d,
               slit_width=a, slit_height=Ly, screen_z=200)
    bar_z, scr_z = cfg["barrier_z"], cfg["screen_z"]
    sim = wave.WaveSim(cfg)
    col = Coll.Collapse(dict(grid=[Lx, Ly, Lz], source_width=18.0,
                             show_collapse=True, collapse_frames=9,
                             collapse_peak=6.0, n_detections=10 ** 8))

    # -- (a) the wave: a *separate* driven sheet beam (continuous mode) so the
    # x-z slice shows the full field diffracting through the two slits and the
    # two waves interfering -- a single electron's small packet
    # would only show a little blob.
    sim_a = wave.WaveSim(dict(cfg, source_mode="continuous",
                              source_shape="sheet"))
    xz_slice = None
    for _i in range(520):
        sim_a.step()
        if _i == 500:
            xz_slice = (sim_a.field()[:, Ly // 2, :] ** 2).copy()

    # -- run the Tonomura flow: one electron at a time, one point each --------
    absorb_frac, max_flight, n_target = 0.02, 300, 150
    sim.emit_packet()
    dets = []
    screened = False
    flight = 0
    time_pattern = None  # (b) the time-integrated interference pattern
    n_steps = 0
    while len(dets) < n_target and n_steps < 200000:
        sim.step()
        n_steps += 1
        flight += 1
        # the front reaches the screen -> arm the screen buffer
        if not screened and (sim.at_screen() or flight > max_flight):
            sim.arm_screen()
            screened = True
            flight = 0
        # every step: integrate the on-screen intensity into the buffer
        sim.accumulate_screen()
        # the sponge has absorbed the electron -> sample the Born-rule point
        # from the *time-integrated* pattern, then fire the next electron
        if screened and sim.absorbed(absorb_frac):
            P = sim.screen_pattern()
            if time_pattern is None:
                time_pattern = P.copy()
            pt = col.sample(P)
            if pt is not None:
                dets.append(pt)
            sim.emit_packet()
            screened = False
            flight = 0
    det = np.array(dets) if dets else np.empty((0, 2))
    if xz_slice is None:
        xz_slice = (sim_a.field()[:, Ly // 2, :] ** 2)
    if time_pattern is None:
        time_pattern = sim.screen_pattern()

    # pick a detected point in a bright fringe for the collapse animation
    if det.size:
        bright = det[time_pattern[det[:, 0], det[:, 1]] > 0.3 * time_pattern.max()]
        if not bright.size:
            bright = det
        cx, cy = max(bright, key=lambda p: time_pattern[int(p[0]), int(p[1])])
        cx, cy = int(cx), int(cy)
    else:
        cx, cy = Lx // 2, Ly // 2
    # the contraction: the screen superposition * a shrinking Gaussian
    gx, gy = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing="ij")
    sig0, sig1, pkmax, nframe = 60.0, 1.5, 6.0, 9
    seq = np.stack([
        time_pattern * np.exp(-((gx - cx) ** 2 + (gy - cy) ** 2) /
                              (2 * (sig0 * (sig1 / sig0) ** (f / (nframe - 1))) ** 2))
        * (1.0 + (pkmax - 1.0) * f / (nframe - 1)) for f in range(nframe)])

    fig = plt.figure(figsize=(12.6, 6.6), facecolor="white")
    # three rows: the profile strip under panel (b) gets its own row instead of
    # being pasted over that panel's tick labels; panel (a) spans rows 0-1 so it
    # keeps its height.
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.22, 1.0],
                          hspace=0.70, wspace=0.26,
                          left=0.05, right=0.97, top=0.89, bottom=0.08)
    # (a) the wave: the propagation + double-slit diffraction (x-z slice)
    ax1 = fig.add_subplot(gs[0:2, 0])
    ax1.imshow(np.sqrt(norm(xz_slice)), cmap=C, origin="lower", aspect="auto",
               extent=[0, Lz, 0, Lx])
    ax1.axvline(bar_z, color="w", ls="--", lw=1.4, alpha=0.9)
    ax1.axvline(scr_z, color="yellow", ls="--", lw=1.4, alpha=0.9)
    ax1.annotate("slits", xy=(bar_z, Lx * 0.5), xytext=(bar_z - 22, Lx * 0.5),
                 color="w", fontsize=12, ha="center")
    ax1.annotate("screen", xy=(scr_z, Lx * 0.5), xytext=(scr_z - 30, Lx * 0.62),
                 color="yellow", fontsize=12, ha="center")
    ax1.set_xlabel("propagation $z$  (source $\\leftarrow$ , screen $\\rightarrow$)",
                   color=TXT, fontsize=11, labelpad=1)
    ax1.set_ylabel("transverse $x$", color=TXT)
    ax1.set_title("(a) the wave: one electron, propagation $+$ diffraction",
                 color=TXT)
    blackize(ax1)
    # (b) the time-integrated interference pattern on the detection screen
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(norm(time_pattern).T, cmap=C, origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly])
    ax2.set_xlabel("transverse $x$", color=TXT)
    ax2.set_ylabel("transverse $y$", color=TXT)
    ax2.set_title("(b) the interference: the time-integrated screen", color=TXT)
    blackize(ax2)
    prof = time_pattern.mean(axis=1)
    axp = fig.add_subplot(gs[1, 1], facecolor="#1a1a1a")
    axp.plot(prof, color=C_ZERO, lw=1.3)
    axp.set_xlim(0, Lx)
    axp.set_ylim(0, prof.max() * 1.2)
    axp.set_xticks([0, Lx // 2, Lx])
    axp.set_yticks([])
    axp.set_title("profile along $x$ (the fringe groups)", fontsize=11,
                  color=LIGHT, pad=2)
    for sp in axp.spines.values():
        sp.set_color("#444")
    axp.tick_params(colors=LIGHT, labelsize=8)
    # (c) the collapse: one electron reduces onto a point
    ax3 = fig.add_subplot(gs[2, 0])
    base = norm(time_pattern)
    ax3.imshow(base.T, cmap=C, origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly], alpha=0.55, vmin=0, vmax=1)
    fin = norm(seq[-1])
    ax3.imshow(np.ma.masked_less(fin.T, 0.18), cmap=C, origin="lower",
               aspect="auto", extent=[0, Lx, 0, Ly], vmin=0, vmax=1)
    ax3.add_patch(mpatches.Circle((cx, Ly - cy), 0.5 * lam / 2.0, fill=False,
                                  color="w", ls="--", lw=1.8))
    ax3.plot([cx], [Ly - cy], marker="*", color="w", ms=15, zorder=5,
             markeredgecolor="k", markeredgewidth=0.5)
    ax3.set_xlabel("transverse $x$", color=TXT)
    ax3.set_ylabel("transverse $y$", color=TXT)
    ax3.set_title("(c) the collapse: one electron reduces onto a point", color=TXT)
    blackize(ax3)
    for i in range(3):
        sax = ax3.inset_axes([0.10 + 0.27 * i, 0.02, 0.24, 0.26], facecolor="black")
        sax.imshow(norm([seq[0], seq[len(seq) // 2], seq[-1]][i]).T, cmap=C,
                   origin="lower", aspect="auto")
        sax.set_xticks([]); sax.set_yticks([])
        sax.set_title(["$\\sigma$ large", "contracting", "the point"][i],
                      fontsize=11, color=LIGHT, pad=1)
        for sp in sax.spines.values():
            sp.set_color("#444")
    # (d) the reality: the single-electron detections cluster into the fringes
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.imshow(base.T, cmap=C, origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly], alpha=0.30, vmin=0, vmax=1)
    xs = det[:, 0]; ys = Ly - det[:, 1]
    ax4.scatter(xs, ys, s=4, c=C_ZERO, alpha=0.8, linewidths=0)
    ax4.set_xlabel("transverse $x$", color=TXT)
    ax4.set_ylabel("transverse $y$", color=TXT)
    ax4.set_title(f"(d) the reality: {len(xs)} single-electron hits cluster into the fringes",
                 color=TXT)
    blackize(ax4)
    fig.suptitle("The CL10 Young double-slit as a single-electron measurement: "
                 "one electron at a time, each reduced onto a point (Tonomura flow)",
                 fontsize=14, color="#111")
    # the gridspec carries its own spacing (tight_layout would push panel (a)'s
    # xlabel onto panel (c)'s title); keep the manual geometry.
    fig.subplots_adjust(left=0.05, right=0.97, top=0.89, bottom=0.08,
                        hspace=0.70, wspace=0.26)
    save(fig, "fig1_young.png")

def fig12_gravity_3d():
    """The Planck-mass black hole, recovered by computation.  Left: the 3-D
    gravity well with the **geodesic** integrated through the CL10
    finite-difference Christoffel symbols (the precessing rosette) and the
    horizon ring.  Right: the **2-route** verification — the Kretschmann
    *computed from the metric* (4-D central finite differences) equals the
    analytic 48M^2/r^6 (the black hole is recovered by computation).  Both
    panels are computed from ``cl10.gravity`` (``metric_schwarzschild``,
    ``_christoffel_at``, ``riemann_4d``/``kretschmann_4d``)."""
    from cl10 import gravity
    from scipy.integrate import solve_ivp

    M = 1.0                        # the 0-parameter black-hole mass (Planck units)
    R = 13.0
    H = gravity.H_FD

    def sw(x):
        return gravity.metric_schwarzschild(M, x)

    # ---- LEFT: the 3-D well + the geodesic (via the CL10 Christoffels) ----
    g = np.linspace(-R, R, 90)
    X, Y = np.meshgrid(g, g)
    rr = np.sqrt(X**2 + Y**2)
    # the well: the metric's time-time part g_00 = -(1 - 2M/r) gives, in the
    # weak-field limit, the potential Phi = -M/r (capped at the horizon).
    Phi = -M / np.maximum(rr, 2 * M)

    def geodesic(lam, u):
        t, r, th, ph, td, rd, thd, phd = u
        G = gravity._christoffel_at(sw, np.array([t, r, th, ph]), H)
        v = np.array([td, rd, thd, phd])
        dv = -np.einsum("mab,a,b->m", G, v, v)
        return [td, rd, thd, phd, dv[0], dv[1], dv[2], dv[3]]

    r0, vt = 10.0, 0.34
    gamma = 1.0 / np.sqrt(1 - vt**2)
    u0 = [0.0, r0, np.pi / 2, 0.0, gamma / np.sqrt(1 - 2 * M / r0), 0.0,
          0.0, gamma * vt / r0]
    sol = solve_ivp(geodesic, [0.0, 600.0], u0, method="DOP853",
                    rtol=1e-8, atol=1e-8, max_step=0.5)
    r = sol.y[1]; ph = sol.y[3]
    tx = r * np.cos(ph); ty = r * np.sin(ph); tz = -M / np.maximum(r, 2 * M)

    a = np.linspace(0, 2 * np.pi, 200)
    hpx, hpy = 2 * M * np.cos(a), 2 * M * np.sin(a)

    # ---- RIGHT: the 2-route Kretschmann (computed vs analytic) ----
    scan = gravity.verify_kretschmann_scan(H)
    rv = np.asarray(scan["radii"])
    Kc, Ka = np.asarray(scan["computed"]), np.asarray(scan["expected"])

    fig = plt.figure(figsize=(11.5, 5.2))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax1.set_proj_type("ortho")
    ax1.plot_surface(X, Y, Phi, cmap="viridis", alpha=0.4, rstride=3, cstride=3)
    ax1.plot(tx, ty, tz, color=C_ZERO, lw=1.5,
             label="the geodesic (CL10 Christoffels)")
    ax1.plot(hpx, hpy, np.full_like(hpx, -M / (2 * M)), color="w", ls="--",
             lw=1.4, label=r"the horizon $r=2M$")
    ax1.scatter([tx[0]], [ty[0]], [tz[0]], color="w", s=20, zorder=5,
                label="start")
    ax1.set_zlabel(r"$\Phi=-M/r$"); ax1.set_title("the well + the geodesic",
                                                 fontsize=12)
    ax1.view_init(elev=55, azim=-40)
    ax1.legend(loc="upper left", fontsize=10)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.loglog(rv, Ka, "-", color="k", lw=1.6, label=r"analytic $K=48M^2/r^6$")
    ax2.loglog(rv, Kc, "o", color=C_ZERO, ms=5,
               label="computed from the CL10 metric")
    resid = float(np.abs(Kc / Ka - 1.0).max())
    ax2.text(0.03, 0.05, f"2-route residual = {resid:.1e}\n(the black hole is "
             "recovered by computation)", transform=ax2.transAxes, fontsize=11,
             bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#888", alpha=0.95))
    ax2.set_xlabel(r"radius $r$  ($M=1$, $r>2M$)")
    ax2.set_ylabel("Kretschmann $K = R^{mnpq}R_{mnpq}$")
    ax2.set_title("2-route: the computed Kretschmann", fontsize=12)
    ax2.legend(loc="upper right", fontsize=10.5)

    fig.suptitle("The black hole recovered by computation: the well + geodesic "
                 "(left), the 2-route Kretschmann (right)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig12_gravity_3d.png")

