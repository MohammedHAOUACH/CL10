# results — the run artifacts

The verdicts and the logs. The **verdicts are tracked** (they are the evidence
the paper cites); the **logs are not** (they are re-runnable output). Back to
the [CL10 README](../README.md).

## Tracked

| path | what it is |
|---|---|
| `validation_summary.json` | the full-suite rollup written by `run_validation.py`: `n`, `n_pass`, `pass`, `hard_tier` + `hard_tier_criterion`, `backend_parent` + `backend_requested` + `n_backend_reported` + `n_backend_gpu` (the GPU coverage of the run), `n_named_residuals`, `named_residuals` (the open problems O1–O9), and one entry per experiment — `id`, `dir`, `tier`, `property`, `pass`, `status`, `residual`, `wall_time`, `returncode`, `crashed`, `backend` (where the experiment says which core computed it), and the experiment's own `verdict` object embedded whole (`data`, `tol`, …). |
| `../experiments/<category>/<Vxx>/results/*.json` | each experiment's own verdict, written next to its code. The rollup folds them, so a subset run still produces a complete summary. |

## Ignored (regenerable)

| pattern | produced by |
|---|---|
| `experiments_run_<date>.log` | `./run_experiments.sh` — the full suite, one line per experiment (the tqdm bar is stderr: it stays on the terminal, out of the log) |
| `figures_run_<date>.log` | `./run_figures.sh` — the 25 figures |
| `legacy/` | run logs kept from before the logging convention |
| `experiments/**/results/kernel.npz` | the Gauss-Laplacian kernel caches the gauged-layer experiments rebuild on demand |

## Regenerate

```bash
python3 run_validation.py --report-only   # re-roll the verdicts already on disk, run nothing
python3 run_validation.py --jobs 8        # the whole suite
./run_experiments.sh                          # the whole suite pinned to cupy (aborts if the card cannot host the workers)
./run_experiments.sh V43_qca_spectrum         # one experiment by exact directory name
```

A verdict is a claim about a run: it records the backend that computed it where
the experiment reports one, the tolerance used, and the wall time. If a number
in the paper moves, the matching verdict file is where it moved first.

A verdict is reported **STALE** (counted as a failure) when
`crashes/<Vxx>.log` is newer than it: the experiment crashed after writing that
verdict, and a subset run would otherwise fold the older PASS into the rollup.
Re-running the experiment replaces the verdict and clears the flag.
