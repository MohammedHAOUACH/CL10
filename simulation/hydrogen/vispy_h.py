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
  * a large **2-D mid-plane cut** (the shape, front-on) on the right, plus two
    text panels: a **top-left** title + key list and a **bottom-left** live
    state with a one-line explanation that adapts to every key press;

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
- **1/2/3**: the **view mode** (electron / photon / potential) -- applies to
  BOTH the 3-D and the 2-D views (mode 2 shows the CL10 photon field).  NOT a
  state re-assignment.
- **space**: pause / resume.   **+ / -**: the speed (QCA steps per frame).
- **q / Esc** : quit.
"""

import numpy as np
import math

from vispy import app, scene
from vispy.color import Colormap
from vispy.util.quaternion import Quaternion

# the colour stops (RGB hex + controls): the RGB versions drive the **2-D
# screen** (fully opaque); the RGBA versions (the ``_make_rgba`` ramp below)
# drive the **3-D volumes** with an alpha channel, so the faint outer shell is
# transparent and you can see *through* it into the core (the nucleus, the
# inner lobes) -- "see inside the cloud".
FIRE_STOPS = (
    ["#000000", "#3a0000", "#cc2200", "#ff8a00", "#ffeb4d"],
    [0.0, 0.2, 0.5, 0.8, 1.0],
)  # the photon (black -> dark red -> red -> orange -> yellow)
DENSITY_STOPS = (
    ["#000000", "#12237a", "#1f9aa0", "#43e08a", "#ffe94d"],
    [0.0, 0.3, 0.55, 0.8, 1.0],
)  # the electron density (black -> blue -> teal -> green -> yellow)
COOL_STOPS = (
    ["#000000", "#003a5a", "#0088cc", "#a0f0ff"],
    [0.0, 0.3, 0.6, 1.0],
)  # the potential well (black -> dark blue -> blue -> light cyan)

FIRE = Colormap(*FIRE_STOPS)
DENSITY = Colormap(*DENSITY_STOPS)
COOL = Colormap(*COOL_STOPS)


def _make_rgba(stops, alpha_min):
    """Build an **RGBA** Colormap from the RGB ``stops``: the same RGB colours,
    but with an alpha channel that ramps from ``alpha_min`` (the faint outer
    shell, value 0) up to **1.0** (the bright core, value 1).

    The bright core stays fully opaque (so it keeps its colour and brightness
    -- a pure ``alpha < 1`` ramp would *dim* the whole cloud), while only the
    dilute halo is made see-through.  That is what "see inside" needs: you look
    *through* the faint shell into the dense core (the nucleus, the far lobe)
    instead of the near side hiding it.  ``alpha_min = 1.0`` recovers the old
    solid look; ``alpha_min = 0`` makes the faint shell fully transparent."""
    def _rgb(h):
        h = h.lstrip("#")
        return (
            int(h[0:2], 16) / 255.0,
            int(h[2:4], 16) / 255.0,
            int(h[4:6], 16) / 255.0,
        )
    colors, controls = stops
    rgba = [(*_rgb(c), alpha_min + (1.0 - alpha_min) * k)
            for c, k in zip(colors, controls)]
    return Colormap(colors=rgba, controls=list(controls))

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
        # the atom's length scale (Bohr radius, in cells).  The camera zooms to
        # THIS (not the whole box): with a large box and a small atom, showing
        # the box would shrink the atom to a speck.
        self.a0 = float(cfg["bohr_radius"])
        win = tuple(int(v) for v in cfg["window_size"])

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
        self.view_3d = scene.ViewBox(
            parent=self.central_widget, camera="arcball",
        )
        self.view_3d.camera.center = [self.L / 2.0] * 3
        # zoom to the atom's size (a0), not the whole box: the atom fills the
        # 3-D view even when the lattice is much larger than the cloud.
        self.view_3d.camera.distance = self.a0 * 7.0
        # a 3/4 tilt (rotate about X) so the z-oriented 2p dumbbell reads as a
        # dumbbell, not an end-on circle.  Set via the quaternion + an explicit
        # view_changed -- the same reliable pattern the Young window uses.
        _tilt = math.radians(-28.0)
        self.view_3d.camera._quaternion = Quaternion(
            *[math.cos(_tilt / 2.0), math.sin(_tilt / 2.0), 0.0, 0.0],
        )
        self.view_3d.camera.view_changed()

        self.stride = int(cfg["volume_stride"])
        self.log_gain = float(cfg["log_gain"])
        self.gamma = float(cfg["alpha_gamma"])
        # the 3-D volume opacity: the (max) alpha of the bright core.  Below
        # 1 the cloud is see-through -- the faint shell is transparent, so the
        # nucleus and the far lobe read through the near side ("see inside").
        self.vol_alpha = float(cfg["volume_alpha"])
        vs = (
            max(1, self.Lx // self.stride),
            max(1, self.Ly // self.stride),
            max(1, self.Lz // self.stride),
        )
        # the electron cloud: an **additive glow** volume (the same reliable
        # method as the Young wave) -- a sphere in the 1s, a two-lobe dumbbell
        # in the 2p, a doughnut in the 3d (the excitation deforms the shape).
        # (an isosurface here draws as flat coloured quads -- the "planes" --
        # and is flaky under some GL drivers, so a smooth glow is used instead.)
        self.elec_data = np.zeros(vs, dtype=np.float32)
        self.elec = scene.visuals.Volume(
            self.elec_data, parent=self.view_3d.scene, method="additive",
            clim=(0.0, 1.0),
        )
        self.elec.cmap = _make_rgba(DENSITY_STOPS, self.vol_alpha)
        self.elec.transform = scene.transforms.STTransform(
            scale=(self.stride, self.stride, self.stride),
        )
        # the photon: an additive glow volume (the packet + the radiation)
        self.phot_data = np.zeros(vs, dtype=np.float32)
        self.phot = scene.visuals.Volume(
            self.phot_data, parent=self.view_3d.scene, method="additive",
            clim=(0.0, 1.0),
        )
        self.phot.cmap = _make_rgba(FIRE_STOPS, self.vol_alpha)
        self.phot.transform = scene.transforms.STTransform(
            scale=(self.stride, self.stride, self.stride),
        )
        # the potential: an additive glow volume (the Coulomb well, -1/r)
        self.pot_data = np.zeros(vs, dtype=np.float32)
        self.pot = scene.visuals.Volume(
            self.pot_data, parent=self.view_3d.scene, method="additive",
            clim=(0.0, 1.0),
        )
        self.pot.cmap = _make_rgba(COOL_STOPS, self.vol_alpha)
        self.pot.transform = scene.transforms.STTransform(
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

        # -- the texts --------------------------------------------------------
        # top-left: the title + the keyboard shortcuts (mostly static).
        # NOTE: VisPy's multi-line Text anchors are inverted (font y-up vs
        # canvas y-down): anchor_y="bottom" puts the anchor at the text's TOP
        # so the block grows *down* from pos.y -- which is what we want here.
        self.info_text = scene.visuals.Text(
            "", parent=self.scene, color="white", font_size=11,
            anchor_x="left", anchor_y="bottom",
        )
        self.info_text.pos = (12, 0)
        # bottom-left: the live state + a one-line explanation (adapts to
        # every key press).  anchor_y="top" -> the anchor is at the text's
        # *bottom*, so the block grows *up* from pos.y (the bottom corner).
        self.state_text = scene.visuals.Text(
            "", parent=self.scene, color="#8fe8ff", font_size=11,
            anchor_x="left", anchor_y="top",
        )
        self.state_text.pos = (12, 12)

        # -- the AGC (for the glow + photon volumes) -------------------------
        rel = float(cfg["agc_release"])
        self._agc_p = AGC(rel)

        self._maj_layout()

        # -- state ------------------------------------------------------------
        self.vis_mode = "electron"
        self.paused = False
        self.steps_per_frame = int(cfg["steps_per_frame"])
        self._last_photon_eV = None
        self._msg = "pret: tire [e] Lyman-alpha pour exciter l'atome (1s->2p)"

        # prime the visuals
        self._maj_volume()
        self._maj_ecran()
        self._maj_text()


        if show:
            self.timer = app.Timer(
                interval=0.016, connect=self.on_timer, start=True,
            )
            self.show()

    # ------------------------------------------------------------------
    # layout: the 3-D view on the left, the 2-D screen on the right
    # ------------------------------------------------------------------
    def _maj_layout(self):
        if not hasattr(self, "view_ecran"):
            return  # the initial resize fires before the views are built
        W, H = self.size
        # the 2-D cut: a large square (bigger than before), centred vertically
        # in the right half of the window.
        w2 = min(int(W * 0.42), H)
        self.view_ecran.pos = (W - w2, (H - w2) // 2)
        self.view_ecran.size = (w2, w2)
        self.view_3d.pos = (0, 0)
        self.view_3d.size = (W - w2, H)
        # the 2-D camera: match the image extent so the cut fills (and centres)
        # the view -- no black border, the slab is centred in the panel.
        sh_y, sh_x = self._sl_shape
        self.view_ecran.camera.set_range(x=(0, sh_x), y=(0, sh_y))
        # the two texts: top-left (title + keys; anchor "bottom" grows *down*
        # from y=0) and bottom-left (state + explanation; anchor "top" grows
        # *up* from y=H) -- both anchored so the block is inside the window.
        self.info_text.pos = (12, 0)
        self.state_text.pos = (12, H)

    def on_resize(self, event):
        super().on_resize(event)
        self._maj_layout()

    # ------------------------------------------------------------------
    # the GPU -> host compression helpers
    # ------------------------------------------------------------------

    def _compress(self, f_gpu, agc):
        """log1p-compress ``f_gpu`` to a subsampled 0..1 host array (for an
        additive glow volume), normalised by the AGC-smoothed peak."""
        lg = B.log1p(B.maximum(f_gpu, B.zeros_like(f_gpu)) * self.log_gain)
        lg_s = lg[:: self.stride, :: self.stride, :: self.stride]
        peak = agc.update(float(B.asnumpy(lg_s.max())))
        data = B.clip(lg_s / max(peak, 1e-9), 0.0, 1.0) ** self.gamma
        return B.asnumpy(data)
    def _compress_electron(self, f_gpu):
        """log1p-compress the (bounded) electron density to a subsampled 0..1
        array, normalised by the *current* peak (no AGC): the brightest part
        is always at full intensity, so the glow stays bright and stable as
        the atom moves between eigenstates (the photon AGC would dim it)."""
        lg = B.log1p(B.maximum(f_gpu, B.zeros_like(f_gpu)) * self.log_gain)
        lg_s = lg[:: self.stride, :: self.stride, :: self.stride]
        peak = max(float(B.asnumpy(lg_s.max())), 1e-9)
        data = B.clip(lg_s / peak, 0.0, 1.0) ** self.gamma
        return B.asnumpy(data)

    def _compress_slab(self, f_gpu):
        """The (2-D) mid-plane cut of ``f_gpu`` (a host 0..1 array): the
        plane at the centre, transposed so the z-oriented dumbbell is
        vertical on screen."""
        cx = self.Lx // 2
        sl = f_gpu[cx, :, :][:: self.stride, :: self.stride].T
        if self.vis_mode == "potential":
            # the well's *magnitude* (deepest = brightest), so the core is
            # bright and the far field dark (normalising the raw negative well
            # would invert it: core dark, edge bright).
            m = B.maximum(-sl, B.zeros_like(sl))
            peak = max(float(B.asnumpy(m.max())), 1e-9)
            img = B.clip(m / peak, 0.0, 1.0) ** 2.0
        else:
            lg = B.log1p(B.maximum(sl, B.zeros_like(sl)) * self.log_gain)
            peak = max(float(B.asnumpy(lg.max())), 1e-9)
            img = B.clip(lg / peak, 0.0, 1.0) ** self.gamma
        return B.asnumpy(img)

    def _compress_potential(self, f_gpu):
        """Compress the (negative) Coulomb well to a subsampled 0..1 array: the
        *depth* of the well is the intensity (deepest = brightest).  A power law
        (not log) keeps the far field dark, so the well reads as a bright core
        instead of a uniform haze (a log keeps a smooth 1/r well bright
        everywhere, which renders the whole box as a white glow)."""
        m = B.maximum(-f_gpu, B.zeros_like(f_gpu))
        m_s = m[:: self.stride, :: self.stride, :: self.stride]
        peak = max(float(B.asnumpy(m_s.max())), 1e-9)
        data = B.clip(m_s / peak, 0.0, 1.0) ** 2.0
        return B.asnumpy(data)

    # ------------------------------------------------------------------
    # the periodic updates
    # ------------------------------------------------------------------
    def on_timer(self, event):
        if not self.paused:
            for _ in range(self.steps_per_frame):
                self.sim.step()
        self._maj_volume()
        self._maj_ecran()
        self._maj_text()
        self.update()

    def _maj_volume(self):
        de = self.sim.electron_density()
        self.elec_data[:] = self._compress_electron(de)
        self.elec.set_data(self.elec_data)
        dp = self.sim.photon_density()
        self.phot_data[:] = self._compress(dp, self._agc_p)
        self.phot.set_data(self.phot_data)
        dw = self.sim.W
        self.pot_data[:] = self._compress_potential(dw)
        self.pot.set_data(self.pot_data)
        # the modes apply to the 3-D too: show only the current mode's volume
        self.elec.visible = self.vis_mode == "electron"
        self.phot.visible = self.vis_mode == "photon"
        self.pot.visible = self.vis_mode == "potential"

    def _maj_ecran(self):
        if self.vis_mode == "electron":
            f = self.sim.electron_density()
        elif self.vis_mode == "photon":
            f = self.sim.photon_density()
        else:
            f = self.sim.W
        img = self._compress_slab(f)
        self.ecran_img.set_data(img)

    def _maj_text(self):
        # top-left: the title + the keyboard shortcuts (mostly static)
        self.info_text.text = (
            "CL10 HYDROGEN [GPU]  (%s)\n"
            "------------------------------\n"
            "[e] Lyman-a 1s->2p   [l] 1s->3p\n"
            "[b] 2p->3d          [o] hors-res.\n"
            "[1/2/3] mode 3D+2D: 1 elec / 2 photon / 3 potent.\n"
            "[space] pause  [+/-] vitesse  [g] 1s  [q] quit"
        ) % ("PAUSE" if self.paused else "RUN")
        # bottom-left: the live state + the one-line explanation
        p = self.sim.populations()
        E = self.sim.atom_energy()
        om = self.sim.omega_ph
        lab = self.MODE_TITRES[self.vis_mode]
        n = self.VIEW_MODES.index(self.vis_mode) + 1
        lines = [
            "t=%-5d  E_atom=%7.2f eV  exc=%5.1f%%"
            % (self.sim.t, E, 100.0 * self.sim.excited_fraction()),
        ]
        for label, pop in zip(self.sim.labels, p):
            lines.append("  %-3s %7.3f" % (label, pop))
        if om > 0:
            eV = self._last_photon_eV
            extra = "  (%.2f eV)" % eV if eV else ""
            lines.append("photon  %7.3f rad/step%s" % (om, extra))
        else:
            lines.append("photon  (aucun)")
        lines.append("mode[%d] %s" % (n, lab))
        lines.append(">>> %s" % self._msg)
        self.state_text.text = "\n".join(lines)

    # ------------------------------------------------------------------
    # the keyboard shortcuts
    # ------------------------------------------------------------------
    def _fire_gap(self, i, j):
        """Reset the atom to level ``i``, then fire a photon at the (i, j)
        gap -- the resonant shot for that transition."""
        self.sim.set_level(i)
        amp = float(self.cfg["photon_amp"])
        width = float(self.cfg["photon_width"]) or None
        self.sim.fire_transition(i, j, amp=amp, width=width)
        self._last_photon_eV = self.sim.gaps_eV.get((i, j))
        self._maj_text()

    def on_key_press(self, event):
        k = event.key
        # VisPy passes a Key object (it compares to its string names, but has
        # no __int__). Normalise to its primary string name (lowercase) so the
        # digit modes and the +/- speed keys work on every backend.
        if hasattr(k, "name"):
            k = k.name
        k = str(k).lower()
        if k in ("q", "escape"):
            self.close()
            return
        elif k == "space":
            self.paused = not self.paused
            self._msg = "pause" if self.paused else "resume"
        elif k in ("+", "=", "add"):
            self.steps_per_frame = min(128, max(1, self.steps_per_frame * 2))
            self._msg = "vitesse: %d steps/frame" % self.steps_per_frame
        elif k in ("-", "_", "subtract"):
            self.steps_per_frame = max(1, self.steps_per_frame // 2)
            self._msg = "vitesse: %d steps/frame" % self.steps_per_frame
        elif k == "e":
            self._fire_gap(0, 1)  # Lyman-alpha: 1s -> 2p
            self._msg = ("Lyman-alpha tire: l'atome absorbe (1s->2p), "
                         "le photon s'efface")
        elif k == "l":
            self._fire_gap(0, 2)  # 1s -> 3p
            self._msg = "1s->3p tire: l'atome absorbe le photon"
        elif k == "b":
            self._fire_gap(1, 3)  # 2p -> 3d
            self._msg = "2p->3d tire (l'atome est mis en 2p d'abord)"
        elif k == "o":
            # the off-resonant control shot (between the gaps)
            self.sim.set_ground()
            self.sim.fire_photon(
                0.13, direction="x",
                amp=float(self.cfg["photon_amp"]),
            )
            self._last_photon_eV = None
            self._msg = "hors-resonance: le photon traverse, AUCUNE absorption"
        elif k in ("g", "r"):
            self.sim.set_ground()
            self._last_photon_eV = None
            self._msg = "reset: l'atome retourne en 1s, les photons sont effaces"
        elif k in ("1", "2", "3"):
            self.vis_mode = self.VIEW_MODES[int(k) - 1]
            self._msg = "mode %d: le 3D et le 2D montrent le(s) %s" % (
                int(k), self.MODE_TITRES[self.vis_mode])
        self._maj_text()

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
