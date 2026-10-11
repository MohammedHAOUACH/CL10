# publication — the journal build

One manuscript source, one command, one PDF. `doc/final_publication.md` is the
**only** source of the paper; nothing in this directory edits it.

```bash
bash publication/build.sh        # -> publication/publication.pdf
```

Current artifact: 40 pp., 25 figures, 8 tables, 53 references, REVTeX 4.2
`aps,prd,twocolumn`. The build is clean: 0 errors, 0 undefined references,
0 stuck floats; the remaining overfull hboxes are the dense display-math paragraphs.

## The pipeline

```
doc/final_publication.md
  └─ preprocess.py     manuscript -> build/publication_gen.md
                       (Divs around figures / tables / the bibliography)
  └─ pandoc            build/publication_gen.md
       --lua-filter=filter.lua    Divs -> \figure / \figure* / \table* /
                                  thebibliography; injects the verified DOIs
                                  from doi.json as clickable doi.org links
       --template=template.tex    revtex4-2, aps/prd, twocolumn
  └─ pdflatex x3       build/publication.tex -> build/publication.pdf
  └─ cp                -> publication/publication.pdf   (the committed artifact)
```

## Files

| file | role |
|---|---|
| `build.sh` | the one command (copies `../doc/figures/` into `build/` first) |
| `preprocess.py` | manuscript -> pandoc-friendly markdown (Divs for figures/tables/refs) |
| `filter.lua` | Divs -> figure / figure* / table* / thebibliography + DOI links |
| `template.tex` | the REVTeX 4.2 document template |
| `fetch_dois.py` | resolves + verifies every reference against Crossref (title similarity, year, volume/page) -> `doi.json`; re-run it when the bibliography changes |
| `doi.json` | the verified `{"<ref number>": "<doi>"}` map consumed by `filter.lua` |
| `publication.pdf` | the committed compiled paper |
| `build/` | git-ignored scratch: intermediates + the figure copies; safe to delete, `build.sh` regenerates everything |

## Prerequisites

`pandoc >= 3`, `texlive` with `pdflatex`, `texlive-publishers` (revtex4-2),
`poppler-utils` (`pdfinfo`, used by the build's final report).

## Figures

The figures live in `doc/figures/` and are generated from the construction
by `python -m make_figures` (see `make_figures/README.md`). Regenerate the
figures first, then rebuild the PDF — `build.sh` copies them into `build/`.

## Checks after a build

`build.sh` prints pages + page size. For the full clean-build claim, grep
`build/publication.log` for `!`, `undefined`, `Overfull`.
