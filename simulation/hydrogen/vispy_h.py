#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The 3-D hydrogen atom -- the **live** VisPy visualizer (CL10, GPU).

The look matches the headless video (the old-project style):

  * the **electron cloud** ``|psi|**2`` as a **bright isosurface** (a teal
    surface at a fixed fraction of the peak) + a faint glow **MIP volume**
    -- a *sphere* in the 1s, a *two-lobe dumbbell* in the 2p, a *doughnut*
    in the 3d.  The excitation **deforms the shape** of the cloud;
  * the **CL10 photon** ``(|F|**2 + |F_out|**2)`` as a faint MIP volume -- the
    packet that flies in, is **absorbed**, and is **re-radiated**;
  * the **nucleus** as a small octahedron at the centre;
  * a **2-D mid-plane cut** (the shape, front-on) + a HUD on the right.

**The excitation is the dynamics, not a re-assignment** ("sans triche"):
pressing ``e`` sends a resonant CL10 photon (a free 3-D FDTD field) and the
Rabi atom absorbs it (the Fermi-golden-rule rate is non-zero only when the
photon's frequency matches a dipole-allowed gap).  The cloud then deforms
1s -> 2p in real time; nothing is imposed on the state.

**GPU-first.**  The FDTD, the Rabi, the densities and the log1p compression
run on CuPy (the GPU); only the small subsampled 0..1 field is copied to the
host for the VisPy GL render.

Keys
----
- **e**  : fire a **Lyman-alpha** photon (1s -> 2p, 10.2 eV) -- the main
  excitation (the atom is reset to the 1s first, then the photon is fired).
- **l**  : fire a **1s -> 3p** photon (12.1 eV).
- **b**  : fire a **2p -> 3d** photon (the atom is set to the 2p first).
- **o**  : fire an **off-resonant** photon (between the gaps) -- the control
  shot: it is *not* absorbed, the atom stays in the 1s.
- **g**  : reset the atom to the 1s and clear the photons.
- **1/2/3**: the 2-D **view mode** (electron / photon / potential) -- NOT a
  state re-assignment.
- **space**: pause / resume.   **+ / -**: the speed (QCA steps per frame).
- **q / Esc** : quit.
"""

import numpy as np

from vispy import app, scene
from vispy.color import Colormap

# the colormaps: purple-free (the standard "inferno" / "viridis" have a
# magenta/purple low end that reads as a "mask" over the black background).
FIRE = Colormap(
    colors=["#000000", "#3a0000", "#cc2200", "#ff8a00", "#ffeb4d"],
    controls=[0.0, 0.2, 0.5, 0.8, 1.0],
)  # the photon (black -> dark red -> red -> orange -> yellow)
DENSITY = Colormap(
    colors=["#000000", "#12237a", "#1f9aa0", "#43e08a", "#ffe94d"],
    controls=[0.0, 0.3, 0.55, 0.8, 1.0],
)  # the electron density (black -> blue -> teal -> green -> yellow)

from . import backend as B
from . import qca as QCA


# ---------------------------------------------------------------------------
# the nucleus (a small octahedron, 6 vertices / 8 faces, centred on the origin)
# ---------------------------------------------------------------------------
def _nucleus_mesh(radius):
    v = np.array(
        [
            [radius, 0, 0], [-radius, 0, 0],
            [0, radius, 0], [0, -radius, 0],
            [0, 0, radius], [0, 0, -radius],
        ],
        dtype=np.float32,
    )
    tri = [
        (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
        (0, 3, 5), (2, 0, 5), (1, 2, 5), (3, 1, 5),
    ]
    cells = np.array([[3, a, b, c] for (a, b, c) in tri], dtype=np.uint32)
    return v, cells


# ---------------------------------------------------------------------------
# the AGC (auto-gain): a flicker-free peak normalisation of a log1p field
# ---------------------------------------------------------------------------
class AGC:
    """Fast attack / slow release on a field's (log1p) peak, so a volume does
    not pulse on / off as a photon packet passes through the box."""

    def __init__(self, release=300.0, floor=1e-6):
        self.release = float(release)
        self.floor = float(floor)
        self.peak = 1.0

    def update(self, peak):
        if peak < self.floor:
            return self.peak
        if peak > self.peak:
            self.peak = peak  # fast attack
        else:
            self.peak = (
                peak / self.release + self.peak * (1.0 - 1.0 / self.release)
            )
        return self.peak


# ---------------------------------------------------------------------------
# the main window
# ---------------------------------------------------------------------------
class HydroVisPy(scene.SceneCanvas):
    """The live VisPy window over the CL10 hydrogen engine."""

    VIEW_MODES = ["electron", "photon", "potential"]
    MODE_TITRES = {
        "electron": "Nuage |psi|²",
        "photon": "Photon |F|²",
        "potential": "Potentiel W",
    }

    def __init__(self, cfg, show=True):
        self.cfg = cfg
        self.Lx, self.Ly, self.Lz = (int(v) for v in cfg["grid"])
        self.L = max(self.Lx, self.Ly, self.Lz)
        win = tuple(int(v) for v in cfg.get("window_size", (1200, 800)))

        scene.SceneCanvas.__init__(
            self,
            title="CL10 — Atome d'hydrogène — 1s→2p par absorption de photon (GPU)",
            size=win, keys="interactive", show=show, bgcolor="#050505",
        )
        self.events.key_press.connect(self.on_key_press)

        # -- the physics (the CL10 Rabi atom + the FDTD photon, on the GPU) --
        self.unfreeze()
        self.sim = QCA.AtomQCA(cfg)
        self.sim.set_ground()

        # -- the 3-D view (left) ---------------------------------------------
        self.view_3d = scene.ViewBox(parent=self.central_widget)
        self.view_3d.camera = "arcball"
        self.view_3d.camera.center = [self.L / 2.0] * 3
        self.view_3d.camera.distance = self.L * 1.9
        # a front-ish 3/4 view (the z-oriented 2p dumbbell reads in 3-D)
        self.view_3d.camera.elevation = -14
        self.view_3d.camera.azimuth = -72

        self.stride = int(cfg.get("volume_stride", 2))
        self.log_gain = float(cfg.get("log_gain", 1.0e4))
        self.gamma = float(cfg.get("alpha_gamma", 1.0))
        self.iso_frac = float(cfg.get("iso_fraction", 0.05))
        rstep = float(cfg.get("relative_step_size", 0.8))
        method = cfg.get("method", "mip")

        vs = (
            max(1, self.Lx // self.stride),
            max(1, self.Ly // self.stride),
            max(1, self.Lz // self.stride),
        )
        # the electron cloud: a **bright isosurface** (a clean surface, no
        # low-value background) -- a sphere in the 1s, a two-lobe dumbbell in
        # the 2p, a doughnut in the 3d (the excitation deforms the shape)
        self.elec_data = np.zeros(vs, dtype=np.float32)
        self.elec = scene.visuals.Isosurface(
            self.elec_data, level=self.iso_frac, parent=self.view_3d.scene,
        )
        self.elec.cmap = DENSITY
        self.elec.color = "#35e0c8"  # the bright teal surface (the headless look)
        self.elec.set_gl_state("translucent", blend=True, depth_test=True)
        self.elec.transform = scene.transforms.STTransform(
            scale=(self.stride, self.stride, self.stride),
        )
        # the photon: a faint MIP volume (the packet + the radiation)
        self.phot_data = np.zeros(vs, dtype=np.float32)
        self.phot = scene.visuals.Volume(
            self.phot_data, parent=self.view_3d.scene, method=method,
            clim=(0.0, 1.0), relative_step_size=rstep,
        )
        self.phot.cmap = FIRE
        self.phot.set_gl_state("translucent", blend=True, depth_test=True)
        self.phot.transform = scene.transforms.STTransform(
            scale=(self.stride, self.stride, self.stride),
        )

        # the nucleus (a small octahedron at the centre)
        v, f = _nucleus_mesh(max(2.0, self.L * 0.02))
        v = v + (self.L / 2.0, self.L / 2.0, self.L / 2.0)
        self.nucleus = scene.visuals.Mesh(
            v, f, color=(1.0, 0.30, 0.12, 1.0),
            parent=self.view_3d.scene,
        )
        self.nucleus.set_gl_state(
            "translucent", blend=True, depth_test=True, cull_face=False,
        )

        # -- the 2-D screen (right): the mid-plane cut -----------------------
        self.view_ecran = scene.ViewBox(
            parent=self.central_widget, camera="panzoom",
        )
        self.view_ecran.bgcolor = "#050505"
        self.view_ecran.camera.set_range(x=(0, self.L), y=(0, self.L))
        self._sl_shape = (
            max(1, self.Lz // self.stride), max(1, self.Ly // self.stride),
        )
        self.ecran_img = scene.visuals.Image(
            np.zeros(self._sl_shape, dtype=np.float32),
            parent=self.view_ecran.scene, cmap=DENSITY,
            clim=[0.0, 1.0], interpolation="bilinear",
        )

        # -- the HUD (the panel, top-right of the 3-D view) ------------------
        self.info_text = scene.visuals.Text(
            "", parent=self.scene, color="white", font_size=10,
            anchor_x="right", anchor_y="top",
        )
        self.info_text.pos = (0.70 * self.size[0], 8)

        # -- the AGC (for the glow + photon volumes) -------------------------
        rel = float(cfg.get("agc_release", 300.0))
        self._agc_g = AGC(rel)
        self._agc_p = AGC(rel)

        self._maj_layout()

        # -- state ------------------------------------------------------------
        self.vis_mode = "electron"
        self.paused = False
        self.steps_per_frame = int(cfg.get("steps_per_frame", 4))
        self._last_photon_eV = None

        # prime the visuals
        self._maj_volume()
        self._maj_ecran()
        self._maj_hud()

        if show:
            self.timer = app.Timer(
                interval=0.016, connect=self.on_timer, start=True,
            )
            self.show()

    # ------------------------------------------------------------------
    # layout: the 3-D view on the left, the 2-D screen on the right
    # ------------------------------------------------------------------
    def _maj_layout(self):
        W, H = self.size
        w2 = min(int(W * 0.30), H)
        self.view_ecran.pos = (W - w2, (H - w2) // 2)
        self.view_ecran.size = (w2, w2)
        self.view_3d.pos = (0, 0)
        self.view_3d.size = (W - w2, H)
        # keep the HUD at the top-right of the 3-D view
        self.info_text.pos = (self.view_3d.pos[0] + self.view_3d.size[0] - 8, 8)

    def on_resize(self, event):
        super().on_resize(event)
        self._maj_layout()

    # ------------------------------------------------------------------
    # the GPU -> host compression helpers
    # ------------------------------------------------------------------

    def _normalize(self, f_gpu):
        """The subsampled density normalised to 0..1 (for the isosurface
        level, which tracks a fixed fraction of the current peak)."""
        a = f_gpu[:: self.stride, :: self.stride, :: self.stride]
        peak = max(float(B.asnumpy(a.max())), 1e-12)
        return (B.asnumpy(a) / peak).astype(np.float32)

    def _compress(self, f_gpu, agc):
        """log1p-compress ``f_gpu`` to a subsampled 0..1 host array (for a
        MIP volume), normalised by the AGC-smoothed peak."""
        lg = B.log1p(B.maximum(f_gpu, B.zeros_like(f_gpu)) * self.log_gain)
        lg_s = lg[:: self.stride, :: self.stride, :: self.stride]
        peak = agc.update(float(B.asnumpy(lg_s.max())))
        data = B.clip(lg_s / max(peak, 1e-9), 0.0, 1.0) ** self.gamma
        return B.asnumpy(data)

    def _compress_slab(self, f_gpu):
        """The (2-D) mid-plane cut of ``f_gpu`` (a host 0..1 array): the
        plane at the centre, transposed so the z-oriented dumbbell is
        vertical on screen."""
        cx = self.Lx // 2
        sl = f_gpu[cx, :, :][:: self.stride, :: self.stride].T
        if self.vis_mode == "potential":
            lo = float(B.asnumpy(B.minimum(sl)))
            hi = float(B.asnumpy(B.maximum(sl)))
            img = B.clip((sl - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
        else:
            lg = B.log1p(B.maximum(sl, B.zeros_like(sl)) * self.log_gain)
            peak = max(float(B.asnumpy(lg.max())), 1e-9)
            img = B.clip(lg / peak, 0.0, 1.0) ** self.gamma
        return B.asnumpy(img)

    # ------------------------------------------------------------------
    # the periodic updates
    # ------------------------------------------------------------------
    def on_timer(self, event):
        if not self.paused:
            for _ in range(self.steps_per_frame):
                self.sim.step()
        self._maj_volume()
        self._maj_ecran()
        self._maj_hud()
        self.update()

    def _maj_volume(self):
        de = self.sim.electron_density()
        self.elec_data[:] = self._normalize(de)
        self.elec.set_data(self.elec_data)
        dp = self.sim.photon_density()
        self.phot_data[:] = self._compress(dp, self._agc_p)
        self.phot.set_data(self.phot_data)

    def _maj_ecran(self):
        if self.vis_mode == "electron":
            f = self.sim.electron_density()
        elif self.vis_mode == "photon":
            f = self.sim.photon_density()
        else:
            f = self.sim.W
        img = self._compress_slab(f)
        self.ecran_img.set_data(img)

    def _maj_hud(self):
        p = self.sim.populations()
        E = self.sim.atom_energy()
        om = self.sim.omega_ph
        lab = self.MODE_TITRES[self.vis_mode]
        n = self.VIEW_MODES.index(self.vis_mode) + 1
        lines = [
            "CL10 HYDROGEN (GPU)  [%s]" % ("PAUSE" if self.paused else "RUN"),
            "t=%-5d  E_atom=%6.2f eV  excited=%4.1f%%"
            % (self.sim.t, E, 100.0 * self.sim.excited_fraction()),
        ]
        for label, pop in zip(self.sim.labels, p):
            lines.append("  %-3s %6.3f" % (label, pop))
        if om > 0:
            eV = self._last_photon_eV
            extra = " (%.2f eV)" % eV if eV else ""
            lines.append("photon  %6.3f rad/step%s" % (om, extra))
        else:
            lines.append("photon  (none)")
        lines.append("mode[%d] %s" % (n, lab))
        lines.append("[e] Lyman-a [l] 1s->3p [b] 2p->3d [o] off [g] 1s")
        lines.append("[space] pause [+/-] speed [1/2/3] mode [q]")
        self.info_text.text = "\n".join(lines)

    # ------------------------------------------------------------------
    # the keyboard shortcuts
    # ------------------------------------------------------------------
    def _fire_gap(self, i, j):
        """Reset the atom to level ``i``, then fire a photon at the (i, j)
        gap -- the resonant shot for that transition."""
        self.sim.set_level(i)
        amp = float(self.cfg.get("photon_amp", 0.7))
        width = float(self.cfg.get("photon_width", 0.0)) or None
        self.sim.fire_transition(i, j, amp=amp, width=width)
        self._last_photon_eV = self.sim.gaps_eV.get((i, j))
        self._maj_hud()

    def on_key_press(self, event):
        k = event.key
        if k in ("q", "escape"):
            self.close()
        elif k == "space":
            self.paused = not self.paused
        elif k == "add":
            self.steps_per_frame = min(128, max(1, self.steps_per_frame * 2))
        elif k == "subtract":
            self.steps_per_frame = max(1, self.steps_per_frame // 2)
        elif k == "e":
            self._fire_gap(0, 1)  # Lyman-alpha: 1s -> 2p
        elif k == "l":
            self._fire_gap(0, 2)  # 1s -> 3p
        elif k == "b":
            self._fire_gap(1, 3)  # 2p -> 3d
        elif k == "o":
            # the off-resonant control shot (between the gaps)
            self.sim.set_ground()
            self.sim.fire_photon(
                0.13, direction="x",
                amp=float(self.cfg.get("photon_amp", 0.7)),
            )
            self._last_photon_eV = None
        elif k in ("g", "r"):
            self.sim.set_ground()
            self._last_photon_eV = None
        elif k in ("1", "2", "3"):
            self.vis_mode = self.VIEW_MODES[int(k) - 1]
        self._maj_hud()

    # ------------------------------------------------------------------
    def run(self):
        app.run()

    def close(self):
        try:
            self.timer.stop()
        except Exception:
            pass
        scene.SceneCanvas.close(self)


def run(cfg, show=True):
    """Build and return the live VisPy window."""
    win = HydroVisPy(cfg, show=show)
    if show:
        win.run()
    return win
