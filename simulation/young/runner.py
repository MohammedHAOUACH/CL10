#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The runner.

``python -m simulation.young``            -> **headless**: run the wave on the
    pinned GPU, render the big 3-D scene off-screen (PyVista), and write the
    video (the wave passing the slits + the collapse events + the screen
    filling).
``python -m simulation.young --interactive`` -> the **live VisPy window**:
    the wave as a glowing 3-D volume (additive, anti-flicker AGC) passing the
    double-slit wall, the detection screen as the hot 2-D panel (the
    interference fringes + the accumulated detections = the *reality*), with a
    HUD and keyboard controls (pause, speed, camera, reset, clear, live
    slit/wavelength tuning).

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
    from . import scene as Scene  # lazy: PyVista is only needed for the headless video (not the interactive VisPy window)
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

    # the single-electron lifecycle (the real Tonomura flow): fire one packet;
    # when its front reaches the screen, **arm the screen buffer** and start
    # integrating the on-screen |psi|^2 over the passage (a single instant is a
    # blob -- the fringes build up only as the diffracted wave crosses the
    # slits and interferes).  Let the wave ride until the **sponge** absorbs
    # it (the buffer now holds the time-averaged interference pattern), take
    # the ONE Born-rule detection from it, and fire the next electron.  No
    # reset in between -- the boundary takes each electron.
    single = sim.source_mode == "single"
    absorb_frac = float(cfg.get("absorb_frac", 0.02))
    max_flight = int(cfg.get("max_flight_steps", 300))  # safety, if the sponge is slow
    if single:
        sim.emit_packet()
    screened = False
    flight_timer = 0

    while t < n_steps:
        # -- the wave step (GPU) --
        sim.step()
        t += 1

        if single:
            flight_timer += 1
            # 1) the front reaches the screen -> arm the screen buffer
            if not screened and (sim.at_screen() or flight_timer > max_flight):
                sim.arm_screen()
                screened = True
            # every step: integrate the on-screen intensity into the buffer
            sim.accumulate_screen()
            # 2) the sponge has absorbed the wave -> the buffer holds the
            #    time-averaged interference pattern; take the ONE Born-rule
            #    detection from it, then fire the next electron (no reset)
            if screened and sim.absorbed(absorb_frac):
                I = sim.screen_pattern()
                if I.sum() > 0 and collapse.n_detected() < collapse.n_detections:
                    pt = collapse.sample(I)
                    if pt is not None:
                        collapse.add(*pt)
                        collapse.set_event(*pt)
                        display_timer = 0
                sim.emit_packet()
                screened = False
                flight_timer = 0
        else:
            # the old continuous-beam mode: sample from the on-screen |psi|**2
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
# interactive: the live VisPy window (3-D volume + 2-D detection screen)
# ---------------------------------------------------------------------------
def run_interactive(cfg):
    """The live window (VisPy): the wave as a glowing 3-D additive volume
    passing the double-slit wall, the detection screen as the hot 2-D panel
    (the interference fringes + the accumulated detections = the *reality*),
    with a HUD and keyboard controls (pause, speed, camera, reset, clear, and
    live slit / wavelength tuning)."""
    from . import vispy_young

    # a real WINDOW needs a Qt GUI backend (PySide6/PyQt).  Without one VisPy
    # silently falls back to its **off-screen EGL** backend: the process runs
    # (this line prints) but **no window appears**.  Detect that up front.
    import importlib.util
    qt = next(
        (n for n in ("PySide6", "PyQt6", "PySide2", "PyQt5")
         if importlib.util.find_spec(n) is not None),
        None,
    )
    if qt is None:
        print("\n[young] WARNING: no Qt GUI backend found (PySide6/PyQt) -- VisPy will "
              "use its off-screen EGL backend and NO window will open.\n"
              "        install one to get the live 3-D window:\n"
              "            pip install PySide6        (or:  pip install cl10[young])\n")
    print(f"[young] interactive window (VisPy, backend={qt or 'EGL/off-screen'}). "
          "keys: Space pause | +/- speed | "
          "s/a slit sep | w/x slit width | e/d wavelength | R reset | "
          "C clear | V camera | Q quit")
    vispy_young.run(cfg)


def main(argv):
    """The dispatcher.  The GPU is pinned in ``__main__.py`` (before the CuPy
    import) from ``config.yaml`` (default **1**); ``--interactive`` opens the
    live VisPy window instead of writing the video."""
    interactive = "--interactive" in argv
    cfg = load_config()
    if interactive:
        run_interactive(cfg)
    else:
        run_headless(cfg)


if __name__ == "__main__":
    main(sys.argv[1:])
