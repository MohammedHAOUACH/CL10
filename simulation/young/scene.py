#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The **3-D** PyVista scene -- a large scene of the Young experiment.

What the screen shows (per frame):

- the **wave** as a big **3-D point cloud**: the field intensity ``psi**2``
  on a subsampled grid, coloured by the colormap (a glowing 3-D volume of the
  diffracting beam);
- the **double-slit barrier** as solid semi-transparent **boxes** with the two
  **slits** as the gaps between them (you see the wave pass *through* them);
- the **detection screen**: a faint plane at ``screen_z`` with the **detected
  points** (the accumulated *reality*) as bright dots building the fringes,
  plus the instantaneous on-screen intensity;
- the **collapse**: when a detection happens, the wave is shown
  **contracting** onto the detected point (a bright marker) -- the
  wave-packet reduction made visible;
- the **parameter panel**: a corner text block with the live values.

The heavy field is kept small enough to draw fast (a stride + a threshold), so
the scene stays interactive.
"""

import numpy as np
import pyvista as pv


class Scene:
    def __init__(self, cfg, off_screen=True):
        self.cfg = cfg
        Lx, Ly, Lz = (int(v) for v in cfg["grid"])
        self.Lx, self.Ly, self.Lz = Lx, Ly, Lz
        self.sub = int(cfg.get("render_subsample", 2))
        self.cmap = cfg.get("colormap", "inferno")
        self.w, self.h = (int(v) for v in cfg.get("render_size", [1280, 720]))

        self.plotter = pv.Plotter(
            off_screen=off_screen,
            window_size=(self.w, self.h),
            border=False,
        )
        self.plotter.set_background("black")

        # a wide camera: pull back so the whole big box is in frame
        cam = list(cfg.get("camera", [0.0, 220.0, 780.0]))
        focal = [Lx / 2.0, Ly / 2.0, Lz / 2.0]
        self.plotter.camera_position = [cam, focal, [0.0, 1.0, 0.0]]
        self.plotter.camera.parallel_projection = False
        self.plotter.camera.view_angle = float(cfg.get("camera_fov", 50.0))

        self._static_added = False
        self._panel = None

    # -- the static geometry (added once) -------------------------------------
    def _add_static(self, sim):
        # the double-slit barrier: the (Lx x Ly) plane at barrier_z, solid
        # everywhere except the two slits (decomposed into 5 boxes).
        bz = sim.barrier_z
        Lx, Ly, Lz = self.Lx, self.Ly, self.Lz
        slits = sorted(sim.slit_rects, key=lambda r: r[0])
        y_top = max(slits[0][3], slits[1][3])
        y_bot = min(slits[0][2], slits[1][2])
        x1l, x1r = slits[0][0], slits[0][1]
        x2l, x2r = slits[1][0], slits[1][1]

        def box(x0, x1, y0, y1, c="gray", o=0.5):
            b = pv.Box((x0, max(x0 + 1, x1), y0, max(y0 + 1, y1), bz - 1, bz + 1))
            self.plotter.add_mesh(b, color=c, opacity=o, show_edges=False)

        # top and bottom strips (full width)
        box(0, Lx, y_top, Ly)
        box(0, Lx, 0, y_bot)
        # the three middle segments between / around the slits
        box(0, x1l, y_bot, y_top)
        box(x1r, x2l, y_bot, y_top)
        box(x2r, Lx, y_bot, y_top)

        # the detection screen: a faint plane at screen_z
        sz = sim.screen_z
        plane = pv.Plane(
            center=((Lx - 1) / 2.0, (Ly - 1) / 2.0, sz),
            direction=(0, 0, 1),
            i_size=Lx - 1, j_size=Ly - 1,
            i_resolution=Lx - 1, j_resolution=Ly - 1,
        )
        self.plotter.add_mesh(
            plane, color="white", opacity=0.06,
            scalars=None,
        )

        # a thin frame around the whole box (so the scene reads as a big 3-D
        # volume, "grande taille")
        frame = pv.Box((0, Lx - 1, 0, Ly - 1, 0, Lz - 1))
        self.plotter.add_mesh(
            frame, color="white", opacity=0.0, show_edges=True,
            edge_color="white", edge_opacity=0.4, line_width=1.5,
        )
        self._static_added = True

    # -- the dynamic field (the wave / the contracted packet) ------------------
    def _add_field(self, intensity):
        """Draw the (subsampled, thresholded) 3-D field as a coloured cloud."""
        I = intensity
        sub = self.sub
        Is = I[::sub, ::sub, ::sub]  # (Lx', Ly', Lz')
        # a low threshold so only the "shape" of the wave is drawn (fast)
        Imax = Is.max() if Is.size else 0.0
        if Imax <= 0:
            return
        mask = Is > (0.02 * Imax)
        idx = np.argwhere(mask)  # (n, 3) -> (x', y', z')
        if idx.size == 0:
            return
        pts = (idx * sub).astype(np.float32)
        val = Is[tuple(idx.T)]
        pd = pv.PolyData(pts)
        pd.point_data["i"] = val / Imax  # normalised 0..1
        self.plotter.add_mesh(
            pd,
            scalars="i",
            cmap=self.cmap,
            render_points_as_spheres=False,
            point_size=2.0,
            opacity=0.95,
            clim=(0.0, 1.0),
            reset_camera=False,
            render=False,
        )

    # -- the detection screen content ------------------------------------------
    def _add_screen(self, sim, screen_intensity, detections):
        sz = sim.screen_z
        Lx, Ly = self.Lx, self.Ly
        # the instantaneous on-screen intensity (a faint heatmap)
        I = np.asarray(screen_intensity, dtype=np.float32)
        Imax = I.max() if I.size else 0.0
        if Imax > 0:
            grid = pv.ImageData(
                dimensions=(Lx, Ly, 1),
                spacing=(1.0, 1.0, 1.0),
                origin=(0.0, 0.0, sz),
            )
            grid["s"] = (I / Imax).ravel(order="F")
            self.plotter.add_mesh(
                grid,
                scalars="s",
                cmap=self.cmap,
                opacity=0.35,
                show_scalar_bar=False,
                render=False,
            )
        # the detected points (the accumulated reality) as bright dots
        if detections:
            d = np.asarray(detections, dtype=np.float32)  # (n, 2) = (x, y)
            pts = np.zeros((d.shape[0], 3), dtype=np.float32)
            pts[:, 0] = d[:, 0]
            pts[:, 1] = d[:, 1]
            pts[:, 2] = sz
            pd = pv.PolyData(pts)
            pd.point_data["b"] = 1.0
            self.plotter.add_mesh(
                pd,
                scalars="b",
                cmap=["black", "yellow", "white"],
                render_points_as_spheres=True,
                point_size=3.0,
                opacity=1.0,
                show_scalar_bar=False,
                reset_camera=False,
                render=False,
            )

    # -- the collapse marker (the detected point, the flash) -------------------
    def _add_collapse_marker(self, x0, y0, sz):
        c = pv.Sphere(
            center=(x0, y0, sz), radius=3.0, theta_resolution=16,
            phi_resolution=16,
        )
        self.plotter.add_mesh(
            c, color="red", opacity=1.0, render=False,
        )

    # -- the parameter panel (a corner text block) -----------------------------
    def set_panel(self, text):
        if self._panel is not None:
            try:
                self.plotter.remove_actor(self._panel)
            except Exception:  # noqa: BLE001
                pass
        self._panel = self.plotter.add_text(
            text,
            position="upper_left",
            font_size=13,
            color="white",
        )

    # -- a full frame ----------------------------------------------------------
    def draw_frame(self, sim, intensity, screen_intensity, detections,
                   collapse_xy=None, render=True):
        self.plotter.clear_actors()
        self._add_static(sim)
        self._add_field(intensity)
        self._add_screen(sim, screen_intensity, detections)
        if collapse_xy is not None:
            self._add_collapse_marker(collapse_xy[0], collapse_xy[1], sim.screen_z)
        if render:
            self.plotter.render()

    def save_frame(self, path):
        self.plotter.screenshot(path)

    def open_writer(self, path, fps):
        self.plotter.open_movie(path, framerate=int(fps))
        return self.plotter.mwriter

    def write_frame(self):
        self.plotter.write_frame()

    def close(self):
        try:
            self.plotter.close()
        except Exception:  # noqa: BLE001
            pass
