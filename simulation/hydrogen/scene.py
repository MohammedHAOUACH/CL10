#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The **3-D** scene -- the headless (off-screen) renderer of the hydrogen
video.

The look (the old-project style, made to read clearly):

- the **electron** cloud ``|psi|**2`` as a **bright solid isosurface**
  (a teal surface at a fixed fraction of the peak density) + a faint
  glow **MIP volume** behind it.  The isosurface is a *sphere* in the 1s,
  a *two-lobe dumbbell* in the 2p, a *doughnut* in the 3d -- the excitation
  **deforms the shape** of the cloud, the way it does in the real atom;
- the **photon** ``(|F|**2 + |F_out|**2)`` as a faint MIP volume -- the packet
  that flies in, is **absorbed**, and is **re-radiated**;
- the **nucleus** as a small octahedron at the centre;
- the **parameter panel** (a corner text block).

The 2-D mid-plane cut is rendered **separately** (see the runner) and pasted
beside this 3-D view, giving the "3-D on the left, 2-D cut on the right"
layout of the live window.

The heavy fields (the FDTD + the Rabi) run on the **GPU** (CuPy); only the
small subsampled field is brought to the host for the render.
"""

import numpy as np
import pyvista as pv

from . import backend as B


class Scene:
    def __init__(self, cfg, off_screen=True):
        self.cfg = cfg
        self.Lx, self.Ly, self.Lz = (int(v) for v in cfg["grid"])
        self.L = max(self.Lx, self.Ly, self.Lz)
        self.sub = int(cfg["volume_stride"])
        self.log_gain = float(cfg["log_gain"])
        self.iso_frac = float(cfg["iso_fraction"])
        self.elec_color = "#35e0c8"  # a bright teal (the electron surface)
        self.win = tuple(cfg["render_size"])
        self.plotter = pv.Plotter(
            off_screen=off_screen, border=False, window_size=self.win,
        )
        self.plotter.set_background("#05050f")
        # the camera looks along +y from the -y side, so the z-oriented
        # 2p dumbbell shows two separated lobes (the clearest view)
        self.plotter.camera.position = (
            0.5 * self.L, -1.25 * self.L, 0.7 * self.L,
        )
        self.plotter.camera.focal_point = (
            self.L / 2, self.L / 2, self.L / 2,
        )
        self.plotter.camera.view_up = (0, 0, 1)
        self._nucleus = None

    # ------------------------------------------------------------------
    # the nucleus (a small octahedron at the centre, built once)
    # ------------------------------------------------------------------
    def _add_nucleus(self):
        if self._nucleus is None:
            r = max(2.0, self.L * 0.02)
            c = np.array([self.L / 2] * 3, dtype=np.float32)
            v = np.array(
                [
                    [r, 0, 0], [-r, 0, 0],
                    [0, r, 0], [0, -r, 0],
                    [0, 0, r], [0, 0, -r],
                ],
                dtype=np.float32,
            ) + c
            tri = [
                (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
                (0, 3, 5), (2, 0, 5), (1, 2, 5), (3, 1, 5),
            ]
            cells = np.array([[3, a, b, c3] for (a, b, c3) in tri])
            self._nucleus = pv.PolyData(v, cells)
        self.plotter.add_mesh(
            self._nucleus, color=(1.0, 0.30, 0.12), style="wireframe",
            line_width=1.5, opacity=1.0, render=False,
        )

    # ------------------------------------------------------------------
    # a (subsampled) 3-D field to a host array
    # ------------------------------------------------------------------
    def _numpy(self, f_gpu):
        a = f_gpu[:: self.sub, :: self.sub, :: self.sub]
        return B.asnumpy(a)

    def _point_image(self, a):
        nx, ny, nz = a.shape
        img = pv.ImageData(
            dimensions=(nx, ny, nz),
            spacing=(self.sub, self.sub, self.sub), origin=(0, 0, 0),
        )
        img.point_data["d"] = np.ascontiguousarray(
            a.transpose(2, 1, 0),
        ).ravel()
        img.set_active_scalars("d")
        return img

    def _cell_image(self, a01):
        nx, ny, nz = a01.shape
        img = pv.ImageData(
            dimensions=(nx + 1, ny + 1, nz + 1),
            spacing=(self.sub, self.sub, self.sub), origin=(0, 0, 0),
        )
        img.cell_data["d"] = np.ascontiguousarray(
            a01.transpose(2, 1, 0),
        ).ravel()
        img.set_active_scalars("d")
        return img

    def _glow_data(self, a):
        """log1p-compress a density to a 0..1 array (for a MIP volume)."""
        lg = np.log1p(np.maximum(a, 0.0) * self.log_gain)
        peak = max(float(lg.max()), 1e-9)
        return np.clip(lg / peak, 0.0, 1.0)

    # ------------------------------------------------------------------
    # a full 3-D frame (the electron isosurface + glow, the photon, nucleus)
    # ------------------------------------------------------------------
    def draw_frame(self, sim, panel_text, render=True):
        self.plotter.clear_actors()

        # -- the electron: a bright isosurface + a faint glow volume --
        de = self._numpy(sim.electron_density())
        if float(de.max()) > 1e-12:
            img = self._point_image(de)
            iso = img.contour(self.iso_frac * float(de.max()))
            self.plotter.add_mesh(
                iso, color=self.elec_color, render=False,
            )
            self.plotter.add_volume(
                self._cell_image(self._glow_data(de)),
                cmap="viridis", opacity=0.12,
                clim=(0.0, 1.0), reset_camera=False, render=False,
            )

        # -- the photon: a faint MIP volume (the packet + the radiation) --
        dp = self._numpy(sim.photon_density())
        if float(dp.max()) > 1e-9:
            self.plotter.add_volume(
                self._cell_image(self._glow_data(dp)),
                cmap="inferno", opacity=0.35,
                clim=(0.0, 1.0), reset_camera=False, render=False,
            )

        self._add_nucleus()
        self.set_panel(panel_text)
        if render:
            self.plotter.render()

    def set_panel(self, text):
        self.plotter.add_text(
            text, position="upper_right", color="white",
            font_size=10, render=False,
        )

    def screenshot(self, path=None):
        """Render and return the 3-D frame as an (H, W, 3) uint8 array (and
        optionally save it to ``path``)."""
        self.plotter.render()
        if path:
            return self.plotter.screenshot(path)
        return self.plotter.screenshot()

    def close(self):
        try:
            self.plotter.close()
        except Exception:
            pass
