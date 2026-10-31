# make_figures — the 25 publication figures, one command

Every figure in `doc/final_publication.md` is **computed from the construction**;
none contains hand-entered data. This package is the single entry point.

The figures are cited by number in [doc/final_publication.md](../doc/final_publication.md);
after regenerating, rebuild the PDFs with
`review/tools/build_pdfs.sh` or the REVTeX build
([publication/README.md](../publication/README.md)). Back to the
[CL10 README](../README.md).

```bash
python3 -m make_figures                # all 25 figures -> doc/figures/
python3 -m make_figures --only 7 12    # a subset, by paper number
python3 -m make_figures --fast         # skip the two slow GPU figures (1, 2)
python3 -m make_figures --cpu          # force the numpy backend
```

## What one run does

1. **`setup_gpu()`** pins the GPU *before* anything imports CuPy: the card id
   comes from `gpu:` in `simulation/young/config.yaml` (exported as
   `CUDA_VISIBLE_DEVICES`), and `CL10_BACKEND=cupy` is set. If cupy is
   missing, or the pinned card has less than 2 GiB free (checked via
   `nvidia-smi`, since creating a CUDA context can itself OOM), the run
   falls back to numpy with a printed notice — never a crash.
2. **The registry** (`__init__.py`, `FIGURES`) runs the 24 matplotlib
   figures in paper order, grouped by subject:

   | module | figures |
   |---|---|
   | `algebra.py` | 3, 5 — the Cl(10) root: the 32-level spectrum, the Brillouin-zone doubler census |
   | `couplings.py` | 4, 6, 13, 19, 20, 21, 22 — 2-loop UV running, the Wilson mechanism, closed loops, couplings at M_Z, the desert solution, the equal-shift theorem, composite content |
   | `dynamics.py` | 1, 12 — the Young double-slit (the Tonomura flow), the 3+1 gravity sector: geodesic + 2-route Kretschmann |
   | `residuals.py` | 7–11, 14–16, 23–25 — doubler cones, isotropy search, tangent velocity algebra, gauged bond layer, vertex correction + ground state, loop selection, the portal dilemma, neutrino protection, the fragment hull, the family count, the named-gap tests |
   | `scales.py` | 17, 18 — the ladder of scales, the loop-integer scan |

3. **Fig. 1 (young)** runs the `simulation/young` QCA (wave + collapse,
   GPU through the same backend switch). It is the slowest figure (~15 min
   on an RTX 3090); `--fast` skips it.
4. **Fig. 2 (hydrogen)** is *not* a matplotlib figure of this package: it
   is the product of `experiments/phenomenology/V62_hydrogen_production`
   (production 296³ GPU run), invoked here as a subprocess so it re-pins its
   own card and writes `doc/figures/fig2_hydrogen.png` +
   `experiments/phenomenology/V62_hydrogen_production/results/fig2_hydrogen.meta.json`
   itself. A nonzero V62 exit code fails the whole run.

## Conventions

- Each function is named `fig<N>_<slug>` and `save()`s exactly
  `doc/figures/fig<N>_<slug>.png`, the number used in the paper. The
  registry numbers must match the `![...](figures/fig<N>_<slug>.png)`
  placements of `doc/final_publication.md`; each figure appears once.
- Figure text carries physics only: no experiment ids, no module names.
- Measured values appearing in figures (PDG, Planck, Super-Kamiokande,
  KATRIN, NuFIT) are reference data, not fits.
- Style and the output directory live in `style.py` (`OUT`).

## After regenerating

Rebuild the PDFs: `./review/tools/build_pdfs.sh` (`doc/final_publication.pdf`),
or `bash publication/build.sh` for the REVTeX build
of the paper (it copies `doc/figures/` into `publication/build/`).
