# doc — the paper and its references

One manuscript, three companions, two references, 25 computed figures. The
manuscript is the **only** source: the companion papers are generated from it,
the figures are computed from the construction, and the PDFs are built from the
Markdown. Nothing here is hand-copied from another file. Back to the
[CL10 README](../README.md).

## The files

| file | what it is |
|---|---|
| `final_publication.md` | **the master**: the whole theory in one manuscript (982 lines, 56 pp compiled, 25 figures, 10 tables, 55 references). Every claim, number and provenance label lives here first. |
| `paper1_core.md` | the core: the root, the automaton, the spectrum, the S-matrix (432 lines, 22 pp). |
| `paper2_gauge_gravity.md` | the gauge and gravity sectors (188 lines, 9 pp). |
| `paper3_phenomenology.md` | the couplings, the flavour sector, the cosmology, the falsifiers (431 lines, 28 pp). |
| `reference.md` | the technical reference: the exact constructions, the numerical protocol, the reproducibility statement, the reading map of the code. |
| `api.md` | the API reference of the 42 `cl10` modules — signatures, returned keys, provenance vocabulary, the runner contract. Regenerated from the live modules, never from memory. |
| `figures/` | the 25 PNGs, every one computed by [`make_figures`](../make_figures/README.md). |
| `*.pdf` | the compiled artifacts (the master + the three companions). |
| `*.tex`, `*.aux`, `*.log`, `*.out` | LaTeX intermediates of the build, **git-ignored**, safe to delete. |

## The dependency direction

```
doc/final_publication.md
   ├─ review/tools/build_series.py ──> doc/paper{1,2,3}_*.md      (headers: review/tools/headers/)
   ├─ review/tools/build_pdfs.sh ────> doc/paper{1,2,3}_*.pdf, doc/final_publication.pdf
   └─ publication/build.sh ──────────> publication/publication.pdf (the REVTeX 4.2 two-column build)

cl10/ + experiments/ ── make_figures ──> doc/figures/*.png ──(embedded by the builds above)
cl10/ ── introspection ───────────────> doc/api.md
```

Edit the master, then regenerate. Editing a companion paper by hand is a bug:
the next `build_series.py` run overwrites it.

## Commands

```bash
python3 review/tools/build_series.py   # master -> the three companion papers
./review/tools/build_pdfs.sh           # fit_tables -> pandoc -> fix_tables -> xelatex x2 -> doc/*.pdf
bash publication/build.sh              # the REVTeX build of the master -> publication/publication.pdf
python3 -m make_figures                # the 25 figures -> doc/figures/  (./run_figures.sh adds the card check)
```

`build_pdfs.sh` needs `pandoc >= 3` and `texlive` (`xelatex`, the pandoc
filters). After any figure or manuscript change, rebuild the PDFs — a PDF older
than its figures is a stale artifact, not a preview.

## Citation

The four compiled PDFs are attached to the
[v1.0.3 release](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0.3)
and archived on [Zenodo](https://zenodo.org/records/23321228); cite
[10.5281/zenodo.23321106](https://doi.org/10.5281/zenodo.23321106).
