# CL10 / simulation/hydrogen — the 3-D hydrogen atom

A **3-D** visualization of the **hydrogen atom** as a CL10 QCA, in which the
atom can be **excited like in reality — by absorbing a resonant photon**.
The **photon** is the free CL10 massless 3-D field (the same 6-neighbour
leapfrog engine as `simulation/young`); the **atom** is the **semiclassical
Rabi / Fermi-golden-rule** model whose levels and wavefunctions are the
**real hydrogen** eigenstates.

The **excitation is not imposed**: a photon packet is fired at the atom, and
the dynamics **decides** whether it is absorbed.  The (Fermi-golden-rule)
dipole drive transfers population **only** when the packet frequency matches
a **transition (gap)** of the atom — the Rabi / selection rule.  A photon at
the **Lyman-alpha** gap (1s→2p, 10.2 eV) is absorbed and the electron cloud
deforms from a *sphere* into a *two-lobe dumbbell*; a photon **between** the
gaps is not absorbed and the atom stays in the 1s.  The excited atom then
**decays by radiating** a photon (the Einstein-A process: the transition
dipole sources an outgoing wave at the gap frequency — the *natural line*).

**Nothing is fitted.**  The level energies are the **Rydberg** formula
`E_n = -13.606/n²` eV, the shapes are the hydrogen eigenfunctions
`R_nl(r) Y_lm`, and the photon is the CL10 engine.  The only free knobs are
the light-matter coupling `A`, the radiated source strength `h`, and the
radiative lifetime `gamma` — physical parameters of the atom, not data fits.

The **whole computation runs on the GPU** (CuPy): the FDTD photon fields,
the dipole drive and the emission source are all on the device; only the
final fields are copied to the host for rendering.

## Run

```bash
cd CL10
source .venv/bin/activate
# headless: runs the atom on the pinned GPU and renders the big 3-D scene
# off-screen (PyVista) — fires resonant photons at the hydrogen gaps and
# writes video/hydrogen.mp4 + stills in video/frames/
python -m simulation.hydrogen
# live window (VisPy): the electron cloud as a glowing 3-D additive volume,
# the photon fields as the hot overlay, a 2-D energy-level panel, and the
# keyboard shortcuts that fire the photons (you trigger the excitation).
python -m simulation.hydrogen --interactive
```

The **GPU** is pinned in `config.yaml` (`gpu: 1`) *before* the CuPy import
(`__main__.py`), so the whole sim runs on that card; without a GPU it
degrades to NumPy (CPU).

## The live window (VisPy) — the keyboard shortcuts

`--interactive` opens a live **VisPy** window with two panels:

- **left** — a 3-D `ViewBox` (arcball camera) with the **electron cloud**
  `|psi|²` as a **bright teal isosurface** (a clean surface, no low-value
  background: a sphere in the 1s, a two-lobe dumbbell in the 2p, a doughnut in
  the 3d — the excitation deforms the shape), the **photon** `|F|² + |F_out|²`
  as a MIP volume with a **purple-free "fire" colormap** (black→red→orange→
  yellow: the incoming packet flying in, the emitted radiation flying out), and
  the **nucleus** (an orange octahedron);
- **right** — a 2-D `ViewBox` (panzoom) with the **mid-plane cut** of the
  current mode (the electron density `|psi|²`, the photon `|F|²`, or the
  Coulomb potential `W`) on a **purple-free "density" colormap** (black→blue→
  teal→green→yellow): the 1s reads as a round dot, the 2p as a **two-lobe
  dumbbell** (two lobes + a node in the middle), the 3d as a doughnut;
- a **HUD** (top-right of the 3-D view): the live populations, the atom energy
  in eV, the excited fraction, the fired photon's frequency, and the mode.

**Keys** (fire the photon — the excitation is triggered by you, like a laser
in the lab):

| key | action |
|---|---|
| `e` | fire a **Lyman-alpha** photon (1s→2p, 10.2 eV, 121.6 nm) — the main excitation |
| `l` | fire a **1s→3p** photon (12.1 eV, 102.6 nm) |
| `b` | fire a **2p→3d** photon (1.89 eV, 656 nm) — the atom is set to the 2p first |
| `o` | fire an **off-resonant** photon (between the gaps) — the control shot: not absorbed |
| `1` / `2` / `3` | the 2-D cut mode: the electron `|psi|²` / the photon `|F|²` / the potential `W` |
| `Space` | pause / resume |
| `+` / `-` | the speed (QCA steps per frame) |
| `g` / `r` | reset the atom to the 1s and clear the photon fields |
| `q` / `Esc` | quit |

**Flicker-free.**  Each volume's transform is set once (never re-attached per
frame, which is what causes a periodic "heartbeat"), and the brightness is
stabilised by an **AGC** (auto-gain): the (log1p) peak of the density is the
exposure target, with a fast attack and a slow exponential release, so the
cloud glows steadily instead of pulsing on/off as a photon packet passes.

## The parameter panel

Every parameter lives in **`config.yaml`** (the single source of truth) and is
shown live in the window corner:

| key | meaning |
|---|---|
| `gpu` | which GPU to pin (default **1** — the free card) |
| `field_dtype` | `float32` (fast) or `float64` (robust for long runs) |
| `grid` | the `(Lx, Ly, Lz)` box of Planck cells |
| `r` | the leapfrog coefficient `c dt/dx` (CFL: `r <= 1/sqrt(3)`) |
| `sponge_width` | the absorbing boundary (kills reflections) |
| `bohr_radius` | the Bohr radius `a0` in cells (the scale of the drawn eigenfunctions) |
| `coupling_absorb` (`A`) | the strength of the (Fermi-golden-rule) dipole drive |
| `coupling_emit` (`h`) | the strength of the dipole source that radiates |
| `reaction_gamma_max` (`gamma`) | the Einstein-A spontaneous-decay rate (per step) |
| `line_width` | the natural linewidth of the absorption line (0 = the lifetime) |
| `photon_width` / `photon_amp` | the fired packet (Gaussian width, amplitude) |
| `n_steps` / `n_video_frames` / `fps` / `out_video` / `out_frames` | the headless video |
| `window_size` / `volume_stride` / `volume_alpha` / `steps_per_frame` / `volume_refresh_rate` / `agc_release` / `agc_gain_max` | the live VisPy window |

## Files

- `config.yaml` — the parameter panel (edit here).
- `backend.py` — the GPU core (CuPy → NumPy).
- `qca.py` — the atom (the real hydrogen levels + eigenfunctions) coupled to
  the CL10 photon: the (Fermi-golden-rule) dipole drive, the Einstein-A
  emission source, the 3-D FDTD photon fields, and the observables.
- `scene.py` — the PyVista 3-D scene (used by the **headless video**: the
  electron cloud, the incoming + emitted photon clouds, the Coulomb well, the
  corner panel).
- `vispy_h.py` — the **live VisPy window**: the 3-D additive volumes (with the
  anti-flicker AGC), the 2-D energy-level panel, the HUD, and the keyboard
  shortcuts that fire the photons.
- `runner.py` / `__main__.py` — the headless video / the live window dispatch.
