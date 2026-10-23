#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The 3-D black hole -- the **live** VisPy visualizer (CL10).

The scene (the fig17 left panel, made live):

  * the **embedding well** Phi = -M/r (the capped Schwarzschild time-time
    part) as a translucent surface;
  * the **horizon ring** r = 2M as a dashed white circle at the well floor;
  * every **fired particle** as a growing 3-D line: the geodesic is
    integrated incrementally (fixed-step RK4 through the CL10 finite-
    difference Christoffels), so you **watch the geodesic being built** --
    the precessing rosette, the plunge through the horizon, the deflected
    fly-by -- instead of seeing a finished curve;
  * a HUD: the key list + the live per-particle status.

Keys
----
- **1** : fire the **rosette** (fig17's orbit: slightly bound, it precesses).
- **2** : fire the **tight rosette** (closer: the precession per turn is
  much larger -- the rosette visibly rotates).
- **3** : fire the **plunge** (below the circular speed: it spirals in).
- **4** : fire the **radial infall** (zero angular momentum: straight down).
- **5** : fire the **fly-by** (a hyperbolic scattering off the hole).
- **p** : fire the **near-photon** particle (v ~ c near the photon sphere:
  the extreme bending, then capture or escape).
- **g** : clear all particles (the hole stays).
- **space**: pause / resume.   **+ / -**: the speed (RK4 steps per frame).
- **q / Esc** : quit.

Nothing is fitted: M = 1 (the 0-parameter Planck-mass hole, G = c = 1), the
Christoffels are the CL10 central finite differences (H_FD = 1e-3), the
geodesic is the CL10 gravity ODE -- the same code path as fig17.
"""

import numpy as np

from vispy import app, scene

from . import geodesic as G


class BHView(scene.SceneCanvas):
    def __init__(self, cfg):
        scene.SceneCanvas.__init__(self, keys="interactive", size=(1280, 800),
                                   title="CL10 — the 3-D black hole: fire a "
                                         "particle, watch the geodesic build")
        self.unfreeze()
        self.events.key_press.connect(self.on_key_press)
        self.bh = G.BlackHole(cfg)
        self.spf = int(cfg["steps_per_frame"])
        self.paused = False
        self.frame = 0

        self.view = self.central_widget.add_view()
        self.view.camera = scene.TurntableCamera(fov=45, elevation=55,
                                                 azimuth=-40, distance=95.0)
        R = float(cfg["well_radius"])
        X, Y, Phi = self.bh.well(R, n=90)
        self._add_well(X, Y, Phi)
        self._add_horizon()

        self.hud = scene.visuals.Text("", parent=self.scene, color="white",
                                      font_size=11, anchor_x="left",
                                      anchor_y="top")
        self.hud.pos = (14, self.size[1] - 14)
        self._update_hud("")

        self._timer = app.Timer(interval=1.0 / 30.0, connect=self._tick,
                                start=True)
        self._lines = {}     # particle id -> LineVisual

    # ------------------------------------------------------------------
    def _add_well(self, X, Y, Phi):
        n = X.shape[0]
        verts = np.column_stack([X.ravel(), Y.ravel(), Phi.ravel()])
        faces = []
        for i in range(n - 1):
            for j in range(n - 1):
                a, b = i * n + j, i * n + j + 1
                c, d = (i + 1) * n + j, (i + 1) * n + j + 1
                faces += [(a, b, d), (a, d, c)]
        mesh = scene.visuals.Mesh(vertices=verts, faces=np.array(faces),
                                  color=(0.25, 0.5, 0.8, 0.35),
                                  parent=self.view.scene)
        mesh.order = 0

    def _add_horizon(self):
        hx, hy, hz = self.bh.horizon_ring()
        pts = np.column_stack([hx, hy, hz])
        scene.visuals.Line(pos=np.vstack([pts, pts[:1]]),
                           color="white", width=2.0, connect="strip",
                           parent=self.view.scene)

    def _line_for(self, p):
        if id(p) not in self._lines:
            col = np.array(
                [int(p.color[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]
                + [1.0])
            ln = scene.visuals.Line(pos=np.zeros((1, 3)), color=col,
                                    width=2.2, connect="strip",
                                    parent=self.view.scene)
            ln.order = 10
            self._lines[id(p)] = ln
        return self._lines[id(p)]

    # ------------------------------------------------------------------
    def _tick(self, _event):
        if self.paused:
            return
        self.bh.advance(self.spf)
        for p in self.bh.particles:
            self._line_for(p).set_data(np.array(p.trail), connect="strip")
        self.frame += 1
        self._update_hud()

    def _update_hud(self, extra=""):
        c = self.bh.counts()
        lines = [
            "CL10 black hole 3-D   M = 1 (0-parameter), frame %d" % self.frame,
            "keys: 1 rosette   2 tight   3 plunge   4 radial   5 fly-by   p v~c",
            "      g clear    space pause    +/- speed    q quit",
            "particles: {}  ({} live, {} captured, {} escaped)".format(
                len(self.bh.particles), c["orbiting"], c["captured"],
                c["escaped"]),
        ]
        for p in self.bh.particles[-8:]:
            lines.append("  %-22s r=%6.2f  %s"
                         % (p.name, p.u[1], p.status))
        if extra:
            lines.append(extra)
        self.hud.text = "\n".join(lines)

    # ------------------------------------------------------------------
    def on_key_press(self, event):
        k = event.key.name.lower() if hasattr(event.key, "name") \
            else str(event.key).lower()
        # VisPy names the digits "One".."Nine", "Zero" -> map to the presets
        k = {"one": "1", "two": "2", "three": "3", "four": "4",
             "five": "5", "six": "6", "seven": "7", "eight": "8",
             "nine": "9", "zero": "0", "escape": "q"}.get(k, k)
        if k in ("escape", "q"):
            self.close()
            return
        if k == " ":
            self.paused = not self.paused
            self._update_hud("paused" if self.paused else "running")
            return
        if k in ("+", "="):
            self.spf = min(400, self.spf * 2)
            self._update_hud("speed x2 -> %d steps/frame" % self.spf)
            return
        if k == "-":
            self.spf = max(1, self.spf // 2)
            self._update_hud("speed /2 -> %d steps/frame" % self.spf)
            return
        if k == "g":
            for ln in self._lines.values():
                ln.visible = False
            self._lines.clear()
            self.bh.particles.clear()
            self._update_hud("cleared")
            return
        p = self.bh.fire(k)
        if p is not None:
            self._update_hud("fired: %s (r0=%.1f)" % (p.name, p.u[1]))
        else:
            self._update_hud("unknown key: %r" % k)


def run(cfg):
    canvas = BHView(cfg)
    canvas.show()
    app.run()
