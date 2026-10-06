#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The runner for the 3-D black hole (CL10).

``python -m simulation.blackhole``            -> **headless**: run the
    geodesics and render the 3-D scene off-screen (PyVista) as a video
    (``video/blackhole.mp4`` + stills).  The scripted sequence of
    ``config.yaml`` fires the particles automatically; each frame advances
    every live particle by ``steps_per_frame`` RK4 substeps, so the video
    shows the geodesics **being constructed** (the rosette precessing, the
    plunge crossing the horizon, the fly-by deflecting).

``python -m simulation.blackhole --interactive`` -> the **live VisPy
    window**: the well, the horizon ring, and the **keyboard shortcuts**
    that fire the particles (you choose what to send and watch the geodesic
    build in real time).

The physics is fig24's: the CL10 Schwarzschild metric, the CL10 finite-
difference Christoffels, the geodesic ODE.  M = 1, nothing fitted.
"""

import os
import sys
import time

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))


def load_config():
    with open(os.path.join(HERE, "config.yaml")) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# headless: the off-screen PyVista video
# ---------------------------------------------------------------------------
def headless(cfg):
    import imageio
    import pyvista as pv

    from . import geodesic as G

    pv.OFF_SCREEN = True
    bh = G.BlackHole(cfg)
    R = float(cfg.get("well_radius", 16.0))
    spf = int(cfg.get("steps_per_frame", 12))
    n_frames = int(cfg.get("n_video_frames", 240))
    fps = int(cfg.get("fps", 30))
    script = cfg.get("script", [])

    plotter = pv.Plotter(off_screen=True, border=False,
                         window_size=tuple(cfg.get("window_size", (1280, 800))))
    plotter.set_background("#05050f")

    X, Y, Phi = bh.well(R, n=90)
    grid = pv.StructuredGrid(X, Y, Phi)
    # shade the well by depth so the funnel reads even at a shallow angle
    grid["depth"] = Phi.ravel()
    plotter.add_mesh(grid, scalars="depth", cmap="coolwarm",
                     opacity=0.5, show_edges=False, lighting=True,
                     show_scalar_bar=False)

    hx, hy, hz = bh.horizon_ring()
    ring = pv.MultipleLines(np.column_stack([hx, hy, hz]))
    plotter.add_mesh(ring, color="white", line_width=2)

    lines = {}                       # particle id -> pv actor

    def refresh(p):
        pts = np.array(p.trail)
        if len(pts) < 2:
            return
        if id(p) in lines:
            lines[id(p)].mapper.SetInputData(pv.MultipleLines(pts))
        else:
            actor = plotter.add_mesh(pv.MultipleLines(pts), color=p.color,
                                     line_width=2.5)
            lines[id(p)] = actor

    hud = plotter.add_text("", position="upper_left", font_size=11,
                           color="white")

    plotter.camera_position = [(0.9 * R, -1.9 * R, 1.15 * R), (0, 0, -6.0),
                               (0, 0, 1)]

    out_video = os.path.join(HERE, cfg.get("out_video", "video/blackhole.mp4"))
    frames_dir = os.path.join(HERE, cfg.get("out_frames", "video/frames"))
    os.makedirs(frames_dir, exist_ok=True)
    os.makedirs(os.path.dirname(out_video), exist_ok=True)
    writer = imageio.get_writer(out_video, fps=fps)

    script = [(int(f), str(k)) for f, k in script]
    t0 = time.time()
    stills = {0, n_frames // 4, n_frames // 2, 3 * n_frames // 4, n_frames - 1}
    for frame in range(n_frames):
        for f, k in script:
            if f == frame:
                bh.fire(k)
        bh.advance(spf)
        for p in bh.particles:
            refresh(p)
        c = bh.counts()
        hud.set_text("upper_left",
                     "CL10 black hole 3-D   M = 1 (0-parameter)\n"
                     "geodesics via the CL10 Christoffels (fig24 physics)\n"
                     "particles: %d  (%d live, %d captured, %d escaped)"
                     % (len(bh.particles), c["orbiting"], c["captured"],
                        c["escaped"]))
        img = plotter.screenshot(frames_dir + "/frame_%04d.png" % frame)
        writer.append_data(img)
        if frame % 40 == 0:
            print("  frame %3d/%d" % (frame, n_frames), flush=True)
    writer.close()
    plotter.close()
    print("[blackhole] frames=%d  particles=%d  %.1fs"
          % (n_frames, len(bh.particles), time.time() - t0), flush=True)
    print("[blackhole] video: %s" % out_video, flush=True)
    return out_video


# ---------------------------------------------------------------------------
def main(argv):
    cfg = load_config()
    if "--interactive" in argv:
        from . import vispy_bh
        vispy_bh.run(cfg)   # opens the window; the canvas is built there
        return
    headless(cfg)


if __name__ == "__main__":
    main(sys.argv[1:])
