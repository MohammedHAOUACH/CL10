# CL10 / simulation/young — the 3-D Young (double-slit) simulation

A **3-D** visualization of the **Young double-slit experiment** as a *measurement
process* — the **Tonomura single-electron** version.  The source fires **one
electron at a time** (a localized 0-parameter wave packet): each packet
**diffracts through the two slits**, flies to the detection screen, and is
**reduced onto a single detected point** (sampled from the on-screen `|psi|**2`,
the Born rule).  The field then **resets** and the **next electron** is fired
(electron after the other).  The **accumulated single detections build the
interference pattern** — the *reality* emerging from the *reductions*.

This is the CL10 **9-to-1 / Z2-projection** collapse (see `doc/publication.md`)
made visible.  **Nothing is fitted**: the wave equation, the discrete Laplacian,
and the Born-rule sampling are all 0-parameter.

## Run

```bash
cd CL10
source .venv/bin/activate
# headless: renders the big 3-D scene off-screen (PyVista) on the pinned GPU
# and writes video/young.mp4 + stills in video/frames/
python -m simulation.young
# live window (VisPy): the wave as a glowing 3-D volume + the detection
# screen as the hot 2-D panel (the fringes + the accumulated detections).
# keys: Space pause | +/- speed | s/a slit sep | w/x slit width |
#       e/d wavelength | R reset | C clear | V camera | Q quit
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
| `grid` | the `(Lx, Ly, Lz)` box — a wide + flat 3-D volume (`256×256×512`, the long `z` axis) |
| `r` | the leapfrog coefficient `c dt/dx` (CFL: `r <= 1/sqrt(3)`) |
| `wavelength` | the source wavelength (sets the fringe spacing) |
| `source_width` / `source_amplitude` / `source_z` / `ramp_steps` | the Gaussian beam (and its `z` position) |
| `source_mode` / `source_shape` / `source_y_width` / `packet_z_width` / `flight_steps` | the **single-electron** source: `single` fires one packet at a time (the Tonomura experiment), `point`/`sheet` shape, and the flight time to the screen |
| `sponge_width` | the absorbing boundary (kills reflections) |
| `barrier_z` / `slit_separation` / `slit_width` / `slit_height` | the two slits |
| `screen_z` | the detection screen |
| `n_detections` / `collapse_every` | the sampling rate (the *reality*) |
| `show_collapse` / `collapse_frames` / `collapse_peak` | the wave-packet reduction |
| `render_size` / `n_video_frames` / `fps` / `n_steps` / `camera` | the rendering |
| `colormap` | the field colormap (headless video) |
| `visualization.{cmap,clim,volume_stride,volume_refresh_rate,steps_per_frame,show_detections,agc_release,agc_gain_max}` | the live VisPy window |

## The live window (VisPy)

`--interactive` opens a live **VisPy** window with two panels:

  `|psi|**2` (a bright cloud that propagates, diffracts through the two slits,
  and interferes), the **double-slit wall** (a dark plane + the two bright slit
  openings), the bounding box, and the detection-screen plane.
  intensity (the **interference fringes**) with the **accumulated detections**
  (the *reality*) drawn as bright dots on top.

**Flicker-free.** The volume's transform is set once (never re-attached per
frame, which is what causes a periodic "heartbeat"), and the brightness is
stabilised by an **AGC** (auto-gain): the 99th-percentile of the density is the
exposure target, with a fast attack and a slow exponential release, so the
volume glows steadily instead of pulsing on/off with each wave salve.

**Keys.** `Space` pause · `+`/`-` speed · `s`/`a` slit separation ·
`w`/`x` slit width · `e`/`d` wavelength · `V` toggle the camera (along the beam
vs perpendicular) · `R` reset the wave + screen · `C` clear the detections ·
`Q`/`Esc` quit.  Slit/wavelength changes are applied **live** (the barrier and
the slit wall rebuild on the fly).

## Files

- `config.yaml` — the parameter panel (edit here).
- `backend.py` — the GPU core (CuPy -> NumPy).
- `wave.py` — the 3-D FDTD wave equation (the beam, the sponge, the slits, the source).
- `collapse.py` — the wave-packet reduction (the Born-rule sampling + the contraction).
- `scene.py` — the PyVista 3-D scene (used by the **headless video**: the field
  cloud, the barrier, the screen, the corner panel).
- `vispy_young.py` — the **live VisPy window**: the 3-D additive volume of
  `|psi|**2` (with the anti-flicker AGC), the double-slit wall, the detection
  screen, the hot 2-D fringe panel with the accumulated detections, the HUD,
  and the keyboard controls.
- `runner.py` / `__main__.py` — the headless video / the live window dispatch.
