# make_figures — the 37 publication figures, one command

Every figure in `doc/publication.md` is **computed from the construction**;
none contains hand-entered data. This package is the single entry point.

```bash
python3 -m make_figures                # all 37 figures -> doc/figures/
python3 -m make_figures --only 7 12    # a subset, by paper number
python3 -m make_figures --fast         # skip the two slow GPU figures (19, 20)
python3 -m make_figures --cpu          # force the numpy backend
```

## What one run does

1. **`setup_gpu()`** pins the GPU *before* anything imports CuPy: the card id
   comes from `gpu:` in `simulation/young/config.yaml` (exported as
   `CUDA_VISIBLE_DEVICES`), and `CL10_BACKEND=cupy` is set. If cupy is
   missing, or the pinned card has less than 2 GiB free (checked via
   `nvidia-smi`, since creating a CUDA context can itself OOM), the run
   falls back to numpy with a printed notice — never a crash.
2. **The registry** (`__init__.py`, `FIGURES`) runs the 36 matplotlib
   figures in paper order, grouped by subject:

   | module | figures |
   |---|---|
   | `algebra.py` | 1, 3, 6, 9, 11, 32, 34 — the Cl(10) root: spectrum, grading, anomalies, CPT mirror, matter content, graphical abstract, doublers |
   | `scales.py` | 2, 4, 5, 7, 10, 13, 14 — WIMP curve, area law, baryon product, ladder, head-to-head, residual reduction |
   | `couplings.py` | 12, 15, 21, 22, 24, 25, 26 — running, 2-loop UV, desert, Wilson, composites, closed loops, equal shift |
   | `proton.py` | 16 — the three decay channels at the desert point |
   | `dynamics.py` | 17, 18, 19, 35 — QCA causal front, 3+1 gravity, **young double-slit**, black-hole frames |
   | `residuals.py` | 8, 23, 27-31, 33, 36, 37 — parameter audit, gauged layer, W2-W9 closures, fragment hull, vertex ground state, named-gap tests |

3. **Fig. 19 (young)** runs the `simulation/young` QCA (wave + collapse,
   GPU through the same backend switch). It is the slowest figure (~15 min
   on an RTX 3090); `--fast` skips it.
4. **Fig. 20 (hydrogen)** is *not* a matplotlib figure of this package: it
   is the product of `experiments/V62_hydrogen_production` (production
   296³ GPU run), invoked here as a subprocess so it re-pins its own card
   and writes `doc/figures/fig20_hydrogen.png` +
   `experiments/V62_hydrogen_production/results/fig20_hydrogen.meta.json`
   itself. A nonzero V62 exit code fails the whole run.

## Conventions

- Each function is named `fig<N>_<slug>` and `save()`s exactly
  `doc/figures/fig<N>_<slug>.png`, the number used in the paper. The
  registry numbers must match the `**Fig. N —**` blocks of
  `doc/publication.md`; the paper references every figure exactly once.
- Measured values appearing in figures (PDG, Planck, Super-Kamiokande,
  KATRIN, NuFIT) are reference data, not fits.
- Style and the output directory live in `style.py` (`OUT`).

## After regenerating

Rebuild the PDF: `bash publication/build.sh` (it copies `doc/figures/`
into `publication/build/`).
