# CL10 / simulation/young — the 3-D Young (double-slit) simulation

A **3-D** visualization of the **Young double-slit experiment** as a *measurement
process*.  A Gaussian beam (a 0-parameter wave) propagates on a big 3-D grid,
**diffracts through the two slits**, and the screen shows the **wave-packet
reduction (collapse)**: the packet **contracts onto a detected point** (sampled
from the on-screen `psi**2`, the Born rule), the point **lands** on the
detection screen, and the accumulated detections build the **interference
pattern** — the *reality* emerging from the *reductions*.

This is the CL10 **9-to-1 / Z2-projection** collapse (see `doc/publication.md`)
made visible.  **Nothing is fitted**: the wave equation, the discrete Laplacian,
and the Born-rule sampling are all 0-parameter.

## Run

```bash
cd CL10
source .venv/bin/activate
# headless: renders the big 3-D scene off-screen on the pinned GPU and writes
# video/young.mp4 + stills in video/frames/
python -m simulation.young
# live window with the interactive parameter panel
# (keys: s/a slit sep, w/x slit width, e/d wavelength, o collapse on/off,
#  p pause, q quit)
python -m simulation.young --interactive
```

The **GPU** is pinned in `config.yaml` (`gpu: 1`) *before* the CuPy import
(`__main__.py`), so the whole sim runs on that card; without a GPU it degrades
to NumPy (CPU).

## The parameter panel

Every parameter lives in **`config.yaml`** (the single source of truth) and is
also shown live in the window corner:

| key | meaning |
|---|---|
| `gpu` | which GPU to pin (default **1** — the free card) |
| `grid` | the `(Lx, Ly, Lz)` box — a big 3-D volume |
| `r` | the leapfrog coefficient `c dt/dx` (CFL: `r <= 1/sqrt(3)`) |
| `wavelength` | the source wavelength (sets the fringe spacing) |
| `source_width` / `source_amplitude` / `ramp_steps` | the Gaussian beam |
| `sponge_width` | the absorbing boundary (kills reflections) |
| `barrier_z` / `slit_separation` / `slit_width` / `slit_height` | the two slits |
| `screen_z` | the detection screen |
| `n_detections` / `collapse_every` | the sampling rate (the *reality*) |
| `show_collapse` / `collapse_frames` / `collapse_peak` | the wave-packet reduction |
| `render_size` / `n_video_frames` / `fps` / `n_steps` / `camera` | the rendering |
| `colormap` | the field colormap |

## Files

- `config.yaml` — the parameter panel (edit here).
- `backend.py` — the GPU core (CuPy -> NumPy).
- `wave.py` — the 3-D FDTD wave equation (the beam, the sponge, the slits, the source).
- `collapse.py` — the wave-packet reduction (the Born-rule sampling + the contraction).
- `scene.py` — the PyVista 3-D scene (the field cloud, the barrier, the screen, the panel).
- `runner.py` / `__main__.py` — the headless video / the live window.
