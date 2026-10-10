# CL10 / simulation/blackhole — the 3-D black hole: fire particles, watch the geodesic build

The black-hole twin of `simulation/hydrogen`: a **3-D live simulation** in
which you **send particles with the keyboard** and watch the **geodesic
construct itself in real time** on the embedding well.

The physics is **exactly fig17's left panel** (`make_figures.fig17_gravity_3d`):
the Schwarzschild metric of `cl10.gravity` (`metric_schwarzschild`, M = 1 —
the 0-parameter Planck-mass hole, G = c = 1), the Christoffel symbols by
**central finite difference** (`_christoffel_at`, H_FD = 1e-3), and the
geodesic equation `d²x/dλ² = −Γ dx/dλ dx/dλ`.  The difference is the
integration: fig17 hands the ODE to `solve_ivp` and draws the finished
curve; here the same ODE is integrated **incrementally (fixed-step RK4)**,
one chunk per rendered frame, so every fired particle leaves a **trail that
grows live** — you see the rosette precess, the plunge cross the horizon,
the fly-by deflect.

**Verified against fig17**: with fig17's initial conditions (r₀ = 10,
v_t = 0.34), the fixed-step RK4 reproduces the DOP853 reference to
r(λ=600) = 7.170386 on both sides (relative agreement ~1e-4 along the way).

Nothing is fitted: the constants of motion (E, L) emerge from the geometry —
the initial 4-velocity is built from a local orthonormal (v_r, v_t) exactly
as in fig17.

## Run

```bash
cd CL10
source .venv/bin/activate
# live window (VisPy): the well + horizon ring + the fire shortcuts
python -m simulation.blackhole --interactive
# headless: the scripted sequence rendered off-screen (PyVista) ->
# video/blackhole.mp4 + stills in video/frames/
python -m simulation.blackhole
```

## The live window — the keyboard shortcuts

| key | fires | what you watch |
|-----|-------|----------------|
| `1` | rosette (fig17's orbit, r₀=10, v=0.34) | the **precessing rosette** building turn by turn |
| `2` | tight rosette (r₀=8, v=0.42) | a much larger **advance per turn** (GR strong field) |
| `3` | plunge (r₀=8, v=0.18) | the **spiral in** through the horizon ring |
| `4` | radial infall (r₀=20, v_r only) | zero angular momentum: **straight down** the well |
| `5` | fly-by (r₀=40, hyperbolic) | the **deflected** scattering trajectory |
| `p` | near-photon (v ≈ c at r₀=8) | extreme bending near the photon sphere |
| `g` | — | clear all trails |
| `space` | — | pause / resume |
| `+` / `-` | — | speed (RK4 steps per frame) |
| `q` / `Esc` | — | quit |

The HUD lists every particle with its live radius and status
(`orbiting` / `captured` past the horizon / `escaped` / `stopped`).

## Files

- `config.yaml` — the parameter panel: M, the presets (r₀, v_r, v_t, colour),
  the RK4 step, the scripted headless sequence.
- `geodesic.py` — the engine: `Particle` (incremental RK4 through the CL10
  Christoffels, growing trail, horizon/escape detection) and `BlackHole`
  (the fired-particle list, the well, the horizon ring).
- `vispy_bh.py` — the live VisPy window (well mesh, horizon ring, growing
  `LineVisual` per particle, HUD, the shortcuts).
- `runner.py` — `--interactive` (VisPy) or headless (PyVista off-screen ->
  `video/blackhole.mp4`).

## Numerical note

The RK4 stages are guarded at the horizon: the Christoffels blow up at
r ≤ 2M, so a stage sampling inside the horizon aborts the step and marks
the particle `captured` (this is what keeps the plunge trails physical
instead of a numerical runaway).
