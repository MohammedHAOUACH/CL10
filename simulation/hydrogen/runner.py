#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The runner for the 3-D hydrogen atom (CL10).

``python -m simulation.hydrogen``            -> **headless**: run the atom on
    the pinned GPU and render it off-screen (PyVista) as a video.  The
    layout is the live-window's: a **3-D MIP volume** of the electron cloud
    ``|psi|**2`` (a sphere in the 1s, a dumbbell in the 2p, a doughnut in
    the 3d -- the excitation deforms the shape) + the fainter photon volume
    (the packet that flies in, is absorbed, and is re-radiated) + the
    nucleus, beside a **2-D mid-plane cut** of the cloud.  The script fires
    **resonant** photon packets at the hydrogen gaps (the Lyman-alpha
    1s->2p, the 1s->3p, a 2p->3d, and one **off-resonant** control).
    Writes the video (``video/hydrogen.mp4``) + stills.

``python -m simulation.hydrogen --interactive`` -> the **live VisPy window**:
    the atom as a **MIP volume** of ``|psi|**2``, the photon fields as the
    faint overlay, the nucleus, a 2-D mid-plane cut, a HUD, and the
    **keyboard shortcuts** that fire the photons (the excitation is
    triggered by the user, the way a laser excites the atom in the lab).

The **levels** are the real hydrogen eigenstates (Rydberg gaps -- 0-parameter,
never fitted), and the photon frequencies are exactly those gaps.
"""

import os
import sys
import time

import numpy as np

import yaml

from . import backend as B
from . import qca as QCA

HERE = os.path.dirname(os.path.abspath(__file__))


def load_config():
    with open(os.path.join(HERE, "config.yaml")) as f:
        return yaml.safe_load(f)


def make_panel(cfg, sim, t):
    """The text of the parameter panel (the live corner block)."""
    p = sim.populations()
    lines = [
        "CL10 Hydrogen 3-D   [%s]" % B.backend_name().upper(),
        "step         : %d" % t,
        "populations  : " + "  ".join(
            ["%s=%.2f" % (sim.labels[k], p[k]) for k in range(len(p))]
        ),
    ]
    om = sim.omega_ph
    if om > 0:
        lines.append("photon       : %.3f rad/step" % om)
    else:
        lines.append("photon       : (none)")
    lines += [
        "E atom     : %.2f eV" % sim.atom_energy(),
        "dipole D   : %.3g" % sim.D,
        "excited    : %.1f %%" % (100.0 * sim.excited_fraction()),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# the scripted photon sequence of the headless video: (fire-step, (i, j))
# ---------------------------------------------------------------------------
def build_script(cfg):
    n_steps = int(cfg["n_steps"])
    # resonant shots at the hydrogen gaps, spread over the timeline: the
    # Lyman-alpha (1s->2p) twice (excite, decay, re-excite), the 1s->3p, a
    # 2p->3d (from the excited level), and one off-resonant control shot.
    return [
        (int(0.08 * n_steps), (0, 1)),
        (int(0.40 * n_steps), (0, 1)),
        (int(0.62 * n_steps), (0, 2)),
        (int(0.80 * n_steps), (1, 3)),
        (int(0.92 * n_steps), (-1, -1)),  # the off-resonant control
    ]


# ---------------------------------------------------------------------------
# the 2-D mid-plane cut (rendered beside the 3-D volume, the live look)
# ---------------------------------------------------------------------------
def _render_2d(f_gpu, log_gain, target_w, target_h):
    """The mid-plane cut of ``f_gpu`` (the plane at the centre, along the
    3-D axis) as a viridis RGB array of shape (target_h, target_w, 3)."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import cm
    from PIL import Image

    Lx = int(f_gpu.shape[0])
    sl = B.asnumpy(f_gpu[Lx // 2, :, :])  # (Ly, Lz)
    sl = sl.T  # (Lz, Ly): the z-oriented dumbbell reads vertically
    lg = np.log1p(np.maximum(sl, 0.0) * log_gain)
    img = np.clip(lg / max(float(lg.max()), 1e-9), 0.0, 1.0)
    rgb = (cm.viridis(img)[..., :3] * 255).astype(np.uint8)
    im = Image.fromarray(rgb).resize((target_w, target_h), Image.BILINEAR)
    return np.asarray(im)


# ---------------------------------------------------------------------------
# headless: the video (the verified deliverable)
# ---------------------------------------------------------------------------
def run_headless(cfg):
    from . import scene as Scene  # lazy: PyVista only for the headless video
    import imageio

    t0 = time.time()
    print("[hydrogen] backend=%s  %s" % (B.backend_name(), B.gpu_info()),
          flush=True)

    sim = QCA.AtomQCA(cfg)
    sim.set_ground()
    print("[hydrogen] levels (rad/step): %s"
          % ["%s=%.4f" % (l, o) for l, o in zip(sim.labels, sim.omega)],
          flush=True)
    print("[hydrogen] gaps (eV):  1s->2p=%.2f  1s->3p=%.2f  2p->3d=%.2f"
          % (sim.gaps_eV[(0, 1)], sim.gaps_eV[(0, 2)], sim.gaps_eV[(1, 3)]),
          flush=True)

    script = build_script(cfg)
    print("[hydrogen] photon script (t, (i,j)): %s" % script, flush=True)

    n_steps = int(cfg["n_steps"])
    n_frames = int(cfg["n_video_frames"])
    frame_every = max(1, n_steps // n_frames)

    out_video = os.path.join(HERE, cfg["out_video"])
    frames_dir = os.path.join(HERE, cfg["out_frames"])
    os.makedirs(os.path.dirname(out_video), exist_ok=True)
    os.makedirs(frames_dir, exist_ok=True)

    from PIL import Image

    sc = Scene.Scene(cfg, off_screen=True)
    w3, h3 = (int(v) for v in cfg["render_size"])
    w2 = int(w3 * 0.5)  # the 2-D cut is half as wide, same height
    log_gain = float(cfg["log_gain"])

    writer = imageio.get_writer(
        out_video, fps=int(cfg["fps"]), codec="libx264",
        quality=8, macro_block_size=1, pixelformat="yuv420p",
    )

    amp = float(cfg["photon_amp"])
    fired = 0
    t = 0
    frame_idx = 0
    stills = {0, n_frames // 4, n_frames // 2, 3 * n_frames // 4, n_frames - 1}

    while t < n_steps:
        sim.step()
        t += 1
        # fire the next scripted photon
        if fired < len(script) and t >= script[fired][0]:
            i, j = script[fired][1]
            if i < 0:
                # the off-resonant control shot (between the gaps)
                sim.set_ground()
                sim.fire_photon(0.13, direction="x", amp=amp)
                print("  [t=%d] fired OFF-RESONANT photon (0.13)"
                      % t, flush=True)
            else:
                sim.set_level(i)
                sim.fire_transition(i, j, amp=amp)
                print("  [t=%d] fired %s -> %s (resonant)"
                      % (t, sim.labels[i], sim.labels[j]), flush=True)
            fired += 1

        if frame_idx < n_frames and t >= (frame_idx + 1) * frame_every:
            sc.draw_frame(sim, make_panel(cfg, sim, t), render=False)
            rgb3d = sc.screenshot()
            rgb2d = _render_2d(sim.electron_density(), log_gain, w2, h3)
            frame = np.concatenate([rgb3d, rgb2d], axis=1)
            writer.append_data(frame)
            if frame_idx in stills:
                Image.fromarray(frame).save(
                    os.path.join(frames_dir, "frame_%04d.png" % frame_idx),
                )
            frame_idx += 1
            if frame_idx % 25 == 0:
                print("  frame %3d/%d  (t=%d)" % (frame_idx, n_frames, t),
                      flush=True)

    writer.close()
    sc.close()
    dt = time.time() - t0
    print("[hydrogen] done in %.1fs -> %s" % (dt, out_video), flush=True)
    print("[hydrogen] frames=%d  photons=%d  steps=%d"
          % (frame_idx, fired, t), flush=True)
    print("[hydrogen] stills in %s" % frames_dir, flush=True)
    return out_video, frames_dir


# ---------------------------------------------------------------------------
# interactive: the live VisPy window
# ---------------------------------------------------------------------------
def run_interactive(cfg):
    """The live window (VisPy) with the keyboard shortcuts."""
    from . import vispy_h
    vispy_h.run(cfg)


def main(argv):
    """The dispatcher.  The GPU is pinned in ``__main__.py`` (before the
    CuPy import) to the ``gpu`` id from the config."""
    cfg = load_config()
    if "--interactive" in argv:
        run_interactive(cfg)
    else:
        run_headless(cfg)


if __name__ == "__main__":
    main(sys.argv[1:])
