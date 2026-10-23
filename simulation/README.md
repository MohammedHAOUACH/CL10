# simulation — three live 3-D simulations

Interactive visualizations of the theory's own objects: each is a runnable
package, not a movie. They share one computation core (`cl10/engine.py`: CuPy →
JAX → NumPy) and one parameter-panel convention (`config.yaml`). Back to the
[CL10 README](../README.md).

| package | what it shows | what is CL10, what is input | paper artifact |
|---|---|---|---|
| [`young/`](young/README.md) | the Young double-slit as a **measurement process**: one electron at a time, diffracting, reduced to one detected point; the accumulated detections build the pattern | the discrete wave equation, the Laplacian and the Born-rule sampling are all 0-parameter | **Fig. 1** (`make_figures` runs it as a figure) |
| [`hydrogen/`](hydrogen/README.md) | the 3-D hydrogen atom driven by a photon wave packet: resonant absorption, the $1s\to2p$ deformation, re-emission | the **atom is an input** (the real Rydberg eigenstates); the **photon is the CL10** free massless 3-D field | **Fig. 2** — the panels are the product of `experiments/phenomenology/V62_hydrogen_production`; this package makes the video and stills |
| [`blackhole/`](blackhole/README.md) | fire particles with the keyboard and watch the geodesic draw itself on the embedding well | the Schwarzschild metric of `cl10.gravity` ($M = 1$, $G = c = 1$), the Christoffel symbols by central finite difference, the geodesic equation; the integration is fixed-step RK4 (Fig. 12 uses `solve_ivp`) | none (the visualization twin of Fig. 12's left panel) |

## Run

```bash
python3 -m simulation.young        # the double-slit (the Fig. 1 scene)
python3 -m simulation.hydrogen     # the hydrogen atom + the photon (the Fig. 2 scene)
python3 -m simulation.blackhole    # the black-hole geodesics
```

Each package also has a `--help`/non-interactive mode for batch rendering, and
its own README for the keyboard controls and the parameter panel.

## Conventions

- **`config.yaml` is the parameter panel.** Every physical setting lives there;
  nothing is fitted. The runner mirrors the same values in a live panel in
  interactive mode.
- **The GPU pin is shared.** `gpu:` in `simulation/young/config.yaml` is the
  card `make_figures` and `run_figures.sh` pin too — one place decides which
  card the repository uses.
- **Generated output is git-ignored**: `video/` in each package
  (`simulation/blackhole/video/frames/` ignored, its `blackhole.mp4` committed
  because it is small). Regenerate with the command above.
- **The physics is the same code as the paper's figures**: the wave equation,
  the lattice photon and the geodesic machinery come from `cl10/`, so a
  simulation and its figure cannot disagree.
