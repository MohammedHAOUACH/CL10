# review — the manuscript pipeline

This directory holds the **tooling** that turns the manuscript into its PDFs,
plus the provenance guard the paper claims. Nothing here edits the physics:
the manuscript is `doc/final_publication.md`. Back to the
[CL10 README](../README.md).

## tools/ — the manuscript pipeline

| file | role |
|---|---|
| `fit_tables.py` | pre-pass: measures each Markdown table's wrapped height and marks the ones that need a smaller font or a split, so LaTeX never clips rows. |
| `fix_tables.py` | post-pass on the generated LaTeX: applies the sizing decisions and repairs the float placement. |
| `build_pdfs.sh` | the one command: `fit_tables` → pandoc → `fix_tables` → `xelatex` ×2 → `doc/*.pdf`. Needs `pandoc >= 3` and `texlive`. |
| `number_audit.py` | the provenance check the paper claims: every decimal the manuscripts print must be a number an experiment artefact prints (leaf or status text), a power-of-ten multiple of one, a bibliography identifier, a section cross-reference, or arithmetic declared in its `DERIVED` table. Exits non-zero on a number nothing in the repository produced. |

```bash
./review/tools/build_pdfs.sh           # -> doc/final_publication.pdf
bash publication/build.sh              # the REVTeX build of the master (see publication/README.md)
python3 review/tools/number_audit.py   # every printed decimal traced to an artefact
```

The bibliography resolver is **not** here: `publication/fetch_dois.py` resolves
and verifies every entry against Crossref and writes the `doi.json` map that the
LaTeX filter turns into clickable links. Re-run it whenever the reference list
changes.
