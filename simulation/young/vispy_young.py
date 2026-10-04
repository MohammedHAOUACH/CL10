#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The **3-D Young experiment** -- live VisPy visualizer (CL10).

Inspired by the ACT double-slit visualizer (``ACI/visualizers/young``): the same
library (VisPy), the same layout and interactions, adapted to the CL10 wave.

Layout
------
- **left**  : a 3-D ``ViewBox`` (arcball camera) with
  - the **wave** as an **additive** ``Volume`` of ``|psi|**2`` (a glowing 3-D
    cloud that propagates, diffracts through the two slits, and interferes),
  - the **double-slit wall** (a dark plane at ``barrier_z`` + the two bright
    slit openings),
  - the **detection screen** (a faint plane at ``screen_z``).
- **right** : a 2-D ``ViewBox`` (panzoom) with the **detection screen** as a
  ``hot`` ``Image`` (the on-screen intensity = the interference fringes) + the
  **accumulated detections** (the *reality*) as bright dots.

Interaction (keyboard)
----------------------
- **Space**  : pause / resume
- **+ / -**  : speed (sim steps per frame)
- **r**      : reset the wave + the detections
- **c**      : clear the detection screen (the accumulated reality)
- **v**      : toggle the 3-D camera -- *along the beam* (the screen faces you,
  you read the fringes)  <->  *perpendicular* (the beam is horizontal, you see
  it travel through the slits)
- **q / Esc** : quit

Flicker-free
------------
- the ``Volume`` transform is set **once** (never re-assigned per frame -- that
  re-attaches the render node and causes a periodic "heartbeat");
- the brightness is stabilised with an **AGC** (auto-gain): the 99th-percentile
  of the density is the target, with a fast attack but a slow exponential
  release, so the volume does not pulse on/off with each wave salve.
"""

import numpy as np

from vispy import app, scene

from . import backend as B
from .collapse import Collapse
from .wave import WaveSim


# ---------------------------------------------------------------------------
# 3-D static geometry (the box, the planes, the slit openings)
# ---------------------------------------------------------------------------
def make_box(parent, Lx, Ly, Lz):
    """The bounding box.  VisPy maps Sim (X, Y, Z) -> (depth, height, width),
    so the box is drawn with view X = Sim Z (the beam, horizontal),
    view Y = Sim Y (height), view Z = Sim X (depth)."""
    vx, vy, vz = float(Lz), float(Ly), float(Lx)
    aretes = np.array([
        [0, 0, 0], [vx, 0, 0], [vx, 0, 0], [vx, vy, 0],
        [vx, vy, 0], [0, vy, 0], [0, vy, 0], [0, 0, 0],
        [0, 0, vz], [vx, 0, vz], [vx, 0, vz], [vx, vy, vz],
        [vx, vy, vz], [0, vy, vz], [0, vy, vz], [0, 0, vz],
        [0, 0, 0], [0, 0, vz], [vx, 0, 0], [vx, 0, vz],
        [vx, vy, 0], [vx, vy, vz], [0, vy, 0], [0, vy, vz],
    ], dtype=np.float32)
    line = scene.visuals.Line(
        pos=aretes, connect="segments",
        color=(0.3, 0.6, 1.0, 0.7), width=1, parent=parent,
    )
    line.set_gl_state("translucent", blend=True, depth_test=False)


def make_plane(parent, z_sim, Lx, Ly, face, border):
    """A faint plane perpendicular to Sim Z (placed at z_sim)."""
    px = float(z_sim)
    sommets = np.array([
        [px, 0, 0], [px, Ly, 0], [px, Ly, Lx], [px, 0, Lx],
    ], dtype=np.float32)
    faces = np.array([[0, 1, 2], [0, 2, 3], [0, 2, 1], [0, 3, 2]], dtype=np.uint32)
    mesh = scene.visuals.Mesh(
        vertices=sommets, faces=faces,
        face_colors=np.array([face] * 4), parent=parent,
    )
    mesh.set_gl_state("translucent", blend=True, depth_test=True, cull_face=False)
    contour = scene.visuals.Line(
        pos=np.array([
            [px, 0, 0], [px, Ly, 0], [px, Ly, Lx], [px, 0, Lx], [px, 0, 0],
        ], dtype=np.float32),
        connect="strip", color=border, width=2, parent=parent,
    )
    contour.set_gl_state("translucent", blend=True, depth_test=False)


def make_slit_outlines(parent, sim):
    """The two slit openings as bright **yellow rectangular outlines** on the
    barrier plane (matching the reference visualizer).  Each rectangle lies in
    the (view Y = Sim Y, view Z = Sim X) plane at view X = ``barrier_z``
    (perpendicular to the beam).  Returns the created line nodes (so the
    outlines can be rebuilt when the slits are changed live)."""
    px = float(sim.barrier_z)
    nodes = []
    for (x0, x1, y0, y1) in sim.slit_rects:
        x0, x1, y0, y1 = float(x0), float(x1), float(y0), float(y1)
        line = scene.visuals.Line(
            pos=np.array([
                [px, y0, x0], [px, y1, x0], [px, y1, x1],
                [px, y0, x1], [px, y0, x0],
            ], dtype=np.float32),
            connect="strip", color=(1.0, 1.0, 0.0, 1.0), width=3,
            parent=parent,
        )
        line.set_gl_state("translucent", blend=True, depth_test=False)
        nodes.append(line)
    return nodes


# ---------------------------------------------------------------------------
# the main window
# ---------------------------------------------------------------------------
class YoungVisPy(scene.SceneCanvas):
    def __init__(self, cfg):
        self.cfg = cfg
        vis = cfg.get("visualization", {})
        self.Lx, self.Ly, self.Lz = (int(v) for v in cfg["grid"])
        win_size = tuple(int(v) for v in vis.get("window_size", (1200, 800)))

        scene.SceneCanvas.__init__(
            self, title="CL10 Young 3-D (VisPy)",
            size=win_size, keys="interactive", show=True,
        )
        self.unfreeze()
        self.central_widget.bgcolor = "#050505"

        # -- the 3-D view (left) ---------------------------------------------
        self.view_3d = scene.ViewBox(parent=self.central_widget, camera="arcball")
        self.view_3d.camera.center = [self.Lz / 2.0, self.Ly / 2.0, self.Lx / 2.0]
        self.view_3d.camera.distance = max(self.Lx, self.Ly, self.Lz) * 2.5
        # two preset orientations:
        #   _cam_side (identity)      -> the beam (Sim Z) is HORIZONTAL
        #   _cam_beam (90 deg about X) -> the screen faces the camera
        from vispy.util.quaternion import Quaternion
        self._cam_beam = Quaternion(*[
            np.cos(np.pi / 4), 0.0, 0.0, np.sin(np.pi / 4),
        ])
        self._cam_side = Quaternion()
        self.view_3d.camera._quaternion = self._cam_side
        self.view_3d.camera.view_changed()
        self._cam_is_beam = False

        # -- the 2-D detection screen (right) --------------------------------
        self.view_ecran = scene.ViewBox(parent=self.central_widget, camera="panzoom")
        self.view_ecran.bgcolor = "#050505"
        self.view_ecran.camera.set_range(x=(0, self.Lx), y=(0, self.Ly))

        self._maj_layout()

        # -- the physics ------------------------------------------------------
        self.sim = WaveSim(cfg)
        self.collapse = Collapse(cfg)
        self.collapse_every = int(cfg.get("collapse_every", 2))
        self._det_timer = 0

        # -- 3-D visuals ------------------------------------------------------
        stride = int(vis.get("volume_stride", cfg.get("render_subsample", 2)))
        self.volume_stride = stride
        vs = (self.Lx // stride, self.Ly // stride, self.Lz // stride)
        self.volume = scene.visuals.Volume(
            np.zeros(vs, dtype=np.float32),
            parent=self.view_3d.scene,
            method="additive",
            clim=tuple(vis.get("clim", (0.0, 1.0))),
        )
        self.volume.cmap = vis.get("cmap", "hot")
        # set ONCE (re-assigning per frame re-attaches the node = a flicker)
        self.volume.transform = scene.transforms.STTransform(
            scale=(stride, stride, stride),
        )

        make_box(self.view_3d.scene, self.Lx, self.Ly, self.Lz)
        # the slit wall: a faint BLUE translucent plane (the barrier) with the
        # two slit openings marked as bright yellow outlines.
        make_plane(
            self.view_3d.scene, self.sim.barrier_z, self.Lx, self.Ly,
            (0.2, 0.4, 1.0, 0.07), (0.3, 0.6, 1.0, 0.9),
        )
        self._slit_nodes = make_slit_outlines(self.view_3d.scene, self.sim)
        # the detection screen: the "projection screen" (a faint orange plane).
        make_plane(
            self.view_3d.scene, self.sim.screen_z, self.Lx, self.Ly,
            (1.0, 0.6, 0.0, 0.07), (1.0, 0.85, 0.0, 0.9),
        )

        # -- 2-D visuals ------------------------------------------------------
        self.image_2d = scene.visuals.Image(
            np.zeros((self.Ly, self.Lx), dtype=np.float32),
            parent=self.view_ecran.scene,
            cmap="hot", clim=(0.0, 1.0), interpolation="bilinear",
        )
        # the accumulated detections (the reality) as bright dots
        self.det_2d = scene.visuals.Markers(
            pos=np.zeros((1, 2)),
            parent=self.view_ecran.scene,
            symbol="disc", size=3,
            face_color=(1.0, 1.0, 0.0, 1.0),
            edge_color=(1.0, 1.0, 0.0, 1.0),
        )
        self.det_2d.set_gl_state("translucent", blend=True, depth_test=False)

        # -- HUD --------------------------------------------------------------
        self.info_text = scene.visuals.Text(
            "", parent=self.scene, color="white",
            font_size=11, anchor_x="left", anchor_y="bottom",
        )
        self.info_text.pos = (14, 14)

        # -- state -------------------------------------------------------------
        self.paused = False
        self.steps_per_frame = int(vis.get("steps_per_frame", 4))
        self.volume_refresh = int(vis.get("volume_refresh_rate", 2))
        self.show_detections = bool(vis.get("show_detections", True))
        self._step = 0
        # the single-electron lifecycle: detect when the front reaches the
        # screen, then wait for the **sponge** to absorb the wave (total
        # energy decayed) before firing the next electron -- no reset in
        # between.  _flight_timer is None until the first packet is fired;
        # _detected marks the on-screen detection; _absorb_frac / _max_flight
        # tune the (short) absorb phase after the detection.
        self._absorb_frac = float(cfg.get("absorb_frac", 0.02))
        self._max_flight = int(cfg.get("max_flight_steps", 300))
        self._flight_timer = None
        self._detected = False
        self._frame = 0
        self._vol_ref = None    # the AGC reference (volume)
        self._scr_ref = None    # the AGC reference (2-D screen)
        self._release = float(vis.get("agc_release", 600.0))
        self._gain_max = float(vis.get("agc_gain_max", 8.0))
        # the volume translucency: the density is scaled by this before the
        # additive blend, so < 1.0 makes the cloud read as a semi-transparent
        # "alpha" volume (like the reference).
        self.volume_alpha = float(vis.get("volume_alpha", 0.7))

        # prime the visuals
        self._maj_volume()
        self._maj_screen()
        self._maj_hud()

        self.timer = app.Timer(interval=0.016, connect=self.on_timer, start=True)
        self.show()

    # ------------------------------------------------------------------
    # layout: the 3-D view takes the left, the 2-D screen the right with the
    # exact (Lx x Ly) proportions of the simulation slice.
    # ------------------------------------------------------------------
    def _maj_layout(self):
        W, H = self.size
        w2 = min(int(W * 0.25), H)
        h2 = min(int(round(w2 * self.Ly / self.Lx)), H)
        self.view_3d.pos = (0, 0)
        self.view_3d.size = (W - w2, H)
        self.view_ecran.pos = (W - w2, (H - h2) // 2)
        self.view_ecran.size = (w2, h2)

    def on_resize(self, event):
        super().on_resize(event)
        if hasattr(self, "view_3d") and hasattr(self, "view_ecran"):
            self._maj_layout()

    # ------------------------------------------------------------------
    # the animation loop
    # ------------------------------------------------------------------
    def on_timer(self, _event):
        if self.paused:
            return
        single = self.sim.source_mode == "single"
        # fire the very first electron before the loop starts
        if single and self._flight_timer is None:
            self.sim.emit_packet()
            self._flight_timer = 0
        for _ in range(self.steps_per_frame):
            self.sim.step()
            self._step += 1
            if single:
                # the single-electron lifecycle: when the front reaches the
                # screen, arm the screen buffer and integrate the on-screen
                # |psi|^2 over the passage (a single instant is a blob -- the
                # fringes build up only as the diffracted wave crosses the
                # slits and interferes).  Let the wave ride until the **sponge**
                # absorbs it (the buffer now holds the time-averaged
                # interference pattern), take the ONE Born-rule detection from
                # it, and fire the next electron -- no reset in between.
                self._flight_timer += 1
                if not self._detected and (
                    self.sim.at_screen() or self._flight_timer > self._max_flight
                ):
                    self.sim.arm_screen()
                    self._detected = True
                # every step: integrate the on-screen intensity into the buffer
                self.sim.accumulate_screen()
                if self._detected and self.sim.absorbed(self._absorb_frac):
                    if self.collapse.n_detected() < self.collapse.n_detections:
                        I = self.sim.screen_pattern()
                        if I.sum() > 0:
                            pt = self.collapse.sample(I)
                            if pt is not None:
                                self.collapse.add(*pt)
                                self.collapse.set_event(*pt)
                    self.sim.emit_packet()  # fire the next electron (no reset)
                    self._detected = False
                    self._flight_timer = 0
            else:
                # the old continuous-beam mode: sample from the on-screen |psi|**2
                self._det_timer += 1
                if (
                    self._det_timer >= self.collapse_every
                    and self.collapse.n_detected() < self.collapse.n_detections
                ):
                    I = self.sim.screen_intensity()
                    if I.sum() > 0:
                        pt = self.collapse.sample(I)
                        if pt is not None:
                            self.collapse.add(*pt)
                            self._det_timer = 0

        self._maj_screen()
        if self._frame % self.volume_refresh == 0:
            self._maj_volume()
        self._frame += 1
        self._maj_hud()
        self.update()

    # ------------------------------------------------------------------
    # AGC (auto-gain control): a stable auto-exposure.  The 99th-percentile of
    # the data is the target; fast attack (follow the peak immediately), slow
    # exponential release (the gain can only rise at a bounded rate) + a hard
    # cap.  This keeps the brightness steady -- no "heartbeat" flicker.
    # ------------------------------------------------------------------
    def _agc(self, data, ref, top):
        target = float(np.percentile(data, 99))
        if target > 1e-9:
            if ref is None or target >= ref:
                ref = target                       # fast attack
            else:
                ref = ref * (1.0 - 1.0 / self._release)   # slow release
                ref = max(ref, target * 0.5)
        if ref is None or ref <= 1e-9:
            return data, None
        gain = min(self._gain_max, top / ref)
        return data * gain, ref

    def _maj_volume(self):
        stride = self.volume_stride
        # |psi|**2, subsampled on the GPU, then to the host (one transfer)
        psi = self.sim.psi[::stride, ::stride, ::stride]
        dens = B.asnumpy(psi * psi)
        dens = np.clip(dens, 0.0, 1e6)
        # a gentle log so the faint tail stays visible while the core glows
        dens = np.log1p(dens * 2000.0) / 10.0
        top = float(self.volume.clim[1]) if hasattr(self.volume.clim, "__getitem__") \
            else 1.0
        dens, self._vol_ref = self._agc(dens, self._vol_ref, top)
        if self.volume_alpha != 1.0:
            dens *= self.volume_alpha
        self.volume.set_data(dens.astype(np.float32))

    def _maj_screen(self):
        I = self.sim.screen_intensity()          # (Lx, Ly)
        Imax = float(I.max()) if I.size else 0.0
        # normalise by the frame's own peak so the (faint) fringes are
        # always readable, whatever the absolute scale of the field.
        if Imax > 1e-30:
            data = (I / Imax).T                  # (Ly, Lx) = (rows, cols)
        else:
            data = np.zeros((self.Ly, self.Lx), dtype=np.float32)
        self.image_2d.set_data(data.astype(np.float32))
        self.image_2d.update()
        # the accumulated detections (the reality) on the same screen
        if self.show_detections:
            det = self.collapse.detections
            if det:
                d = np.asarray(det, dtype=np.float32)   # (n, 2) = (x, y)
                self.det_2d.set_data(pos=d, size=3,
                                     face_color=(1, 1, 0, 1), symbol="disc")
            else:
                self.det_2d.set_data(pos=np.zeros((1, 2), dtype=np.float32),
                                     size=3, face_color=(1, 1, 0, 1),
                                     symbol="disc")
            self.det_2d.update()

    def _maj_hud(self):
        status = "PAUSE" if self.paused else "RUN"
        self.info_text.text = (
            f"CL10 Young 3-D  (VisPy, {B.backend_name()})\n"
            f"step    : {self._step}\n"
            f"wavelength : {self.sim.wavelength:.1f}\n"
            f"slits   : sep={self.sim.slit_separation:.1f} w={self.sim.slit_width:.1f}\n"
            f"detections : {self.collapse.n_detected()}\n"
            f"speed   : {self.steps_per_frame}x   {status}\n"
            f"[Space] pause | [+/-] speed | [s/a] sep | [w/x] width\n"
            f"[e/d] wavelength | [R] reset | [C] clear | [V] camera | [Q] quit"
        )
        self.info_text.update()

    def _rebuild_slit_wall(self):
        """Detach the old slit outlines and rebuild them from the current
        (live-tuned) slits."""
        for n in self._slit_nodes:
            n.parent = None
        self._slit_nodes = make_slit_outlines(self.view_3d.scene, self.sim)

    # ------------------------------------------------------------------
    # keyboard
    # ------------------------------------------------------------------
    def on_key_press(self, event):
        k = event.text
        if k == " ":
            self.paused = not self.paused
            self._maj_hud()
        elif k == "+":
            self.steps_per_frame = min(20, self.steps_per_frame + 1)
            self._maj_hud()
        elif k == "-":
            self.steps_per_frame = max(1, self.steps_per_frame - 1)
            self._maj_hud()
        elif k.lower() == "r":
            # reset: a fresh wave + a fresh screen
            self.sim = WaveSim(self.cfg)
            self.collapse = Collapse(self.cfg)
            self._step = 0
            self._det_timer = 0
            self._flight_timer = None
            self._detected = False
            self._vol_ref = None
            self._scr_ref = None
            self._maj_screen()
            self._maj_volume()
            self._maj_hud()
        elif k.lower() == "c":
            self.collapse.detections.clear()
            self.det_2d.set_data(pos=np.zeros((1, 2), dtype=np.float32),
                                 size=3, face_color=(1, 1, 0, 1),
                                 symbol="disc")
            self.det_2d.update()
            self._maj_hud()
        elif k.lower() == "s":
            self.sim.slit_separation = min(120, self.sim.slit_separation + 2)
            self.sim._build_barrier()
            self._rebuild_slit_wall()
            self._maj_hud()
        elif k.lower() == "a":
            self.sim.slit_separation = max(6, self.sim.slit_separation - 2)
            self.sim._build_barrier()
            self._rebuild_slit_wall()
            self._maj_hud()
        elif k.lower() == "w":
            self.sim.slit_width = min(30, self.sim.slit_width + 1)
            self.sim._build_barrier()
            self._rebuild_slit_wall()
            self._maj_hud()
        elif k.lower() == "x":
            self.sim.slit_width = max(2, self.sim.slit_width - 1)
            self.sim._build_barrier()
            self._rebuild_slit_wall()
            self._maj_hud()
        elif k.lower() == "e":
            self.sim.wavelength = min(40, self.sim.wavelength + 1)
            self._maj_hud()
        elif k.lower() == "d":
            self.sim.wavelength = max(4, self.sim.wavelength - 1)
            self._maj_hud()
        elif k.lower() == "v":
            self._cam_is_beam = not self._cam_is_beam
            self.view_3d.camera._quaternion = (
                self._cam_beam if self._cam_is_beam else self._cam_side
            )
            self.view_3d.camera.view_changed()
        elif k.lower() == "q" or k == "\x1b":
            self.close()


def run(cfg, show=True):
    """Build and return the VisPy window (and run the vispy event loop)."""
    cfg = dict(cfg)
    win = YoungVisPy(cfg)
    if show:
        app.run()
    return win
