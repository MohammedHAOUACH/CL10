#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The runner.

``python -m simulation.young``            -> **headless**: run the wave on the
    pinned GPU, render the big 3-D scene off-screen, and write the video
    (the wave passing the slits + the collapse events + the screen filling).
``python -m simulation.young --interactive`` -> the same, but with a live
    window and a **parameter panel** (keyboard shortcuts to tune the slits,
    the wavelength, the speed, the collapse, ...).

The **collapse** (the wave-packet reduction) is shown on the screen exactly as
the user asked: the field is seen **contracting** onto the detected point, the
particle **lands** on the detection screen, and the accumulated detections build
the interference pattern (the *reality* emerging from the *reductions*).
"""

import os
import sys
import time

import yaml

from . import backend as B
from . import collapse as Collapse
from . import scene as Scene
from . import wave as WaveSim

HERE = os.path.dirname(os.path.abspath(__file__))


def load_config():
    with open(os.path.join(HERE, "config.yaml")) as f:
        return yaml.safe_load(f)


def make_panel(cfg, sim, collapse, t):
    """The text of the parameter panel (the live corner block)."""
    return (
        f"CL10 Young 3-D   [{B.backend_name().upper()}]\n"
        f"GPU id      : {B.gpu_id()}\n"
        f"step        : {t}\n"
        f"slit sep    : {sim.slit_separation:6.1f}\n"
        f"slit width  : {sim.slit_width:6.1f}\n"
        f"slit height : {sim.slit_height:6.1f}\n"
        f"wavelength  : {sim.wavelength:6.1f}\n"
        f"r (CFL)     : {sim.r:6.2f}\n"
        f"detections  : {collapse.n_detected()}\n"
        f"collapse    : {'ON' if collapse.show else 'off'}"
    )


# ---------------------------------------------------------------------------
# headless: the video (the verified deliverable)
# ---------------------------------------------------------------------------
def run_headless(cfg):
    t0 = time.time()
    print(f"[young] backend={B.backend_name()}  {B.gpu_info()}")

    sim = WaveSim.WaveSim(cfg)
    collapse = Collapse.Collapse(cfg)
    scene = Scene.Scene(cfg, off_screen=True)

    n_steps = int(cfg["n_steps"])
    n_frames = int(cfg["n_video_frames"])
    collapse_every = int(cfg["collapse_every"])
    frame_every = max(1, n_steps // n_frames)
    # the *displayed* collapse is shown slowly (once per display_every steps) so
    # it spans ``collapse_frames`` rendered frames and completes before the next;
    # *all* detections (the screen / the reality) fill at the faster collapse_every.
    display_every = max(
        collapse_every, frame_every * int(cfg.get("collapse_frames", 3))
    )

    out_video = os.path.join(HERE, cfg.get("out_video", "video/young.mp4"))
    frames_dir = os.path.join(HERE, cfg.get("out_frames", "video/frames"))
    os.makedirs(os.path.dirname(out_video), exist_ok=True)
    os.makedirs(frames_dir, exist_ok=True)

    writer = scene.open_writer(out_video, fps=int(cfg.get("fps", 30)))

    t = 0
    frame_idx = 0
    collapse_timer = 0
    display_timer = 0
    next_frame_at = 0
    stills = {0, n_frames // 4, n_frames // 2, 3 * n_frames // 4, n_frames - 1}

    while t < n_steps:
        # -- the wave step (GPU) --
        sim.step()
        t += 1

        # -- every detection (fills the screen = the "reality") --
        collapse_timer += 1
        if (
            collapse_timer >= collapse_every
            and collapse.n_detected() < collapse.n_detections
        ):
            I = sim.screen_intensity()
            if I.sum() > 0:
                pt = collapse.sample(I)
                if pt is not None:
                    collapse.add(*pt)
                    collapse_timer = 0

        # -- the *displayed* collapse (the wave-packet reduction) --
        display_timer += 1
        if display_timer >= display_every and not collapse.event_active():
            if collapse.detections:
                x, y = collapse.detections[-1]
                collapse.set_event(x, y)
            display_timer = 0

        # -- draw a frame when due --
        if frame_idx < n_frames and t >= next_frame_at:
            field = sim.field()
            if collapse.event_active():
                field = collapse.contract(field)
                collapse.advance_event()  # one advance per rendered frame
            intensity = field * field
            cp = (collapse.event[0], collapse.event[1]) if collapse.event_active() else None
            scene.draw_frame(sim, intensity, sim.screen_intensity(),
                             collapse.detections, cp, render=False)
            scene.set_panel(make_panel(cfg, sim, collapse, t))
            scene.plotter.render()
            scene.write_frame()
            if frame_idx in stills:
                scene.save_frame(os.path.join(frames_dir, f"frame_{frame_idx:04d}.png"))
            frame_idx += 1
            next_frame_at += frame_every
            if frame_idx % 25 == 0:
                print(f"  frame {frame_idx:3d}/{n_frames}  (t={t}, "
                      f"detections={collapse.n_detected()})")

    writer.close()
    scene.close()
    dt = time.time() - t0
    print(f"[young] done in {dt:.1f}s -> {out_video}")
    print(f"[young] frames={frame_idx}  detections={collapse.n_detected()}  "
          f"steps={t}")
    print(f"[young] stills in {frames_dir}")
    return out_video, frames_dir


# ---------------------------------------------------------------------------
# interactive: the live window + the parameter panel (keyboard)
# ---------------------------------------------------------------------------
def run_interactive(cfg):
    sim = WaveSim.WaveSim(cfg)
    collapse = Collapse.Collapse(cfg)
    scene = Scene.Scene(cfg, off_screen=False)
    p = scene.plotter

    n_steps = int(cfg["n_steps"])
    collapse_every = int(cfg["collapse_every"])
    step_ms = 50  # ~20 fps; the actual pace is bounded by the render time
    # mutable state shared by the key handler and the timer callback
    st = {"running": True, "collapse_timer": 0}


    def key_q():
        st["running"] = False
        p.iren.interactor.TerminateApp()

    def key_s():
        sim.slit_separation += 2
        sim._build_barrier()

    def key_a():
        sim.slit_separation = max(6, sim.slit_separation - 2)
        sim._build_barrier()

    def key_w():
        sim.slit_width = min(30, sim.slit_width + 1)
        sim._build_barrier()

    def key_x():
        sim.slit_width = max(2, sim.slit_width - 1)
        sim._build_barrier()

    def key_e():
        sim.wavelength = min(40, sim.wavelength + 1)

    def key_d():
        sim.wavelength = max(4, sim.wavelength - 1)

    def key_o():
        collapse.show = not collapse.show

    def key_p():
        st["running"] = not st["running"]

    for key, cb in (
        ("q", key_q),
        ("s", key_s),
        ("a", key_a),
        ("w", key_w),
        ("x", key_x),
        ("e", key_e),
        ("d", key_d),
        ("o", key_o),
        ("p", key_p),
        (" ", key_p),
    ):
        p.add_key_event(key, cb)

    def tick(step):
        if not st["running"] or sim.t >= n_steps:
            return
        sim.step()
        # a detection (fills the screen)
        st["collapse_timer"] += 1
        if (
            st["collapse_timer"] >= collapse_every
            and collapse.n_detected() < collapse.n_detections
        ):
            I = sim.screen_intensity()
            if I.sum() > 0:
                pt = collapse.sample(I)
                if pt is not None:
                    collapse.add(*pt)
                    collapse.set_event(*pt)
                    st["collapse_timer"] = 0
        collapse.advance_event()
        # the field (contracted onto the detected point while the event is active)
        field = sim.field()
        if collapse.event_active():
            field = collapse.contract(field)
        cp = (collapse.event[0], collapse.event[1]) if collapse.event_active() else None
        # rebuild all actors (render=False) then set the panel; Timer.execute
        # does the single render after this callback, so the panel is included
        # in one frame (no double render, no flicker).
        scene.draw_frame(sim, field * field, sim.screen_intensity(),
                         collapse.detections, cp, render=False)
        scene.set_panel(make_panel(cfg, sim, collapse, sim.t))

    # duration is the per-tick interval in ms (a repeating VTK timer), NOT the
    # total: n_steps ticks at step_ms each.
    p.add_timer_event(max_steps=n_steps, duration=step_ms, callback=tick)
    print("[young] interactive window. keys: s/a slit sep, w/x slit width, "
          "e/d wavelength, o collapse on/off, p pause, q quit")
    # an initial draw so the window is not black before the first tick
    field = sim.field()
    scene.draw_frame(sim, field * field, sim.screen_intensity(),
                     collapse.detections, None, render=False)
    scene.set_panel(make_panel(cfg, sim, collapse, sim.t))
    scene.plotter.render()
    p.show()
    scene.close()


def main(argv):
    """The dispatcher.  The GPU is pinned in ``__main__.py`` (before the CuPy
    import) from ``config.yaml`` (default **1**); ``--interactive`` opens the
    live window + the parameter panel instead of writing the video."""
    interactive = "--interactive" in argv
    cfg = load_config()
    if interactive:
        run_interactive(cfg)
    else:
        run_headless(cfg)


if __name__ == "__main__":
    main(sys.argv[1:])
