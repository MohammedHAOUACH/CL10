# publication — the journal build

One manuscript source, one command, one PDF. `doc/final_publication.md` is the
**only** source of the paper; nothing in this directory edits it.

The unified master builds here; the three companion papers build through
`review/tools/build_pdfs.sh`. All four compiled PDFs are attached to the
[v1.0 release](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0)
and archived on [Zenodo](https://zenodo.org/records/23321228)
([10.5281/zenodo.23321106](https://doi.org/10.5281/zenodo.23321106)). Back to the
[CL10 README](../README.md).

```bash
bash publication/build.sh        # -> publication/publication.pdf
```

Current artifact: 71 pp., 25 figures, 18 table divs, 69 references, REVTeX 4.2
`aps,prd,twocolumn`. The build is clean: 0 errors, 0 undefined references, all
floats placed in ascending order within two pages of their first textual
reference; the remaining overfull hboxes are the dense display-math paragraphs.
The three companion papers (`doc/paper1_core.md`, `paper2_gauge_gravity.md`,
`paper3_phenomenology.md`) are built by the other pipeline,
`review/tools/build_pdfs.sh` (pandoc + xelatex); this directory builds only
the unified master.

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

In the manuscript each figure is an image line followed by its italic caption
line, and every figure is cited in the prose (`(Fig. 7)` / `(Figs. 13, 14)`).
`preprocess.py` binds each pair into a `{.figure}` / `{.figurewide}` Div (the
wide class is chosen from the PNG aspect ratio, > 1.8:1) and moves it after
the paragraph of its first `\ref{fig:N}`; `filter.lua` emits
`\begin{figure}[!ht]` / `\begin{figure*}[!t]` with `\caption` + `\label`.
Adding a figure therefore needs no LaTeX in the manuscript: image, caption,
and one in-text reference.

Tables are sized the same way, from an estimate of the wrapped height (real
word wrapping per cell, TeX kept verbatim the way pandoc's `stringify` does):
a one-column `{.table}` float when the table fits a column, a full-width
`{.tablewide}` float otherwise. A table taller than one float becomes a
**page-reserved group**: `{.tablepage}` + `{.tablecont .tablepage}` floats
emitted with `[p]` placement (their own pages, no interleaved text, no blank
page from a `\clearpage` flushing unrelated figures), set `\scriptsize` with
`\arraystretch` 0.90, each repeating the header row, the continuation carrying
an italic "Table N (continued)" line without consuming a table number. The
number of floats `K` is the smallest that fits the page budget and the row
boundary is the DP-optimal one minimising the tallest float: a greedy fill
lands a point over budget and buys a page the table does not need, and a
`table*` taller than the page is **clipped** — rows silently vanish, which is
the bug this group exists to prevent.

## Checks after a build

`build.sh` prints pages + page size. For the full clean-build claim, grep
`build/publication.log` for `!`, `undefined`, `Overfull`. Then
`python3 review/tools/number_audit.py` re-traces every decimal the manuscript
prints to the artefact that produces it — a build is only clean if that passes.
