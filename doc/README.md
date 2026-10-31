# doc — the paper and its references

One consolidated manuscript, two references, 25 computed figures. The
manuscript is the **only** paper file: the figures are computed from the
construction, the PDFs are built from the Markdown, and the journal layout is
generated from the same source. Nothing here is hand-copied from another file.
Back to the [CL10 README](../README.md).

## The files

| file | what it is |
|---|---|
| `final_publication.md` | **the consolidated paper**: the whole theory in one manuscript (1223 lines, 90 pp compiled, 75 pp in the journal layout, 25 figures, 8 tables, 69 references). Every claim, number and provenance label lives here. |
| `reference.md` | the technical reference: the exact constructions, the numerical protocol, the reproducibility statement, the reading map of the code. |
| `api.md` | the API reference of the 44 `cl10` modules — signatures, returned keys, provenance vocabulary, the runner contract. Regenerated from the live modules, never from memory. |
| `figures/` | the 25 PNGs, every one computed by [`make_figures`](../make_figures/README.md). |
| `*.pdf` | the compiled paper. |
| `*.tex`, `*.aux`, `*.log`, `*.out` | LaTeX intermediates of the build, **git-ignored**, safe to delete. |

## The dependency direction

```
doc/final_publication.md
   ├─ review/tools/build_pdfs.sh ────> doc/final_publication.pdf
   └─ publication/build.sh ──────────> publication/publication.pdf (the REVTeX 4.2 two-column build)

cl10/ + experiments/ ── make_figures ──> doc/figures/*.png ──(embedded by the builds above)
cl10/ ── introspection ───────────────> doc/api.md
```

Edit the manuscript, then regenerate.

## Commands

```bash
./review/tools/build_pdfs.sh           # fit_tables -> pandoc -> fix_tables -> xelatex x2 -> doc/final_publication.pdf
bash publication/build.sh              # the REVTeX build of the manuscript -> publication/publication.pdf
python3 -m make_figures                # the 25 figures -> doc/figures/  (./run_figures.sh adds the card check)
```

`build_pdfs.sh` needs `pandoc >= 3` and `texlive` (`xelatex`, the pandoc
filters). After any figure or manuscript change, rebuild the PDFs — a PDF older
than its figures is a stale artifact, not a preview.

## Citation

The compiled PDFs are the assets of the
[v1.2 GitHub release](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.2).
The archived [Zenodo snapshot](https://zenodo.org/records/23321228) is from
the v1.0.2 era (four then-current paper PDFs + source archive); the concept DOI
[10.5281/zenodo.23321106](https://doi.org/10.5281/zenodo.23321106) always
resolves to the latest archived version.
