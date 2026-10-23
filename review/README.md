# review — the audit trail of the current paper, and the build tooling

Two things live here: the **evidence** about the current manuscript (a referee
simulation, its claim-by-claim verification, the response, the bibliography
audit), and the **tooling** that turns the master into the four PDFs. Nothing
here edits the physics: the manuscript is `doc/final_publication.md`. Back to
the [CL10 README](../README.md).

## The audit documents

| file | what it is |
|---|---|
| `glmreview.md` | a simulated PRD referee report on the master as submitted (934 lines, ≈29 000 words, 25 figures, 55 references): the major and minor concerns, in the order a referee would raise them. |
| `glmreview_verification.md` | every claim of that report checked against the manuscript, the repository code, and new computations run during the check. Verdicts: **CONFIRMED** (the report is right, the manuscript changed), **CONFIRMED-WRONG-DETAIL** (the concern is real, the report's arithmetic or citation is not), **REFUTED** (the claim does not survive checking). |
| `response_letter.md` | the point-by-point response, and the decision to resubmit as the three-paper series (`doc/paper1_core.md`, `paper2_gauge_gravity.md`, `paper3_phenomenology.md`). |
| `refs_audit.json` | the bibliography audit of the **current** numbering: 55 entries, contiguity, duplicate titles, uncited entries, dangling citations, and the per-entry Crossref verdict (`verdict`, `score`, `crossref_title`). |

## tools/ — the manuscript pipeline

| file | role |
|---|---|
| `build_series.py` | the master → the three companion papers: it splits `doc/final_publication.md` by section, renumbers the sections and figures per paper, rewrites the cross-references, and prepends/appends the per-paper headers. |
| `headers/` | `paper{1,2,3}_head.md` and `paper{1,2,3}_foot.md`: the title block, the author line, the citation and availability footers of each companion. |
| `fit_tables.py` | pre-pass: measures each Markdown table's wrapped height and marks the ones that need a smaller font or a split, so LaTeX never clips rows. |
| `fix_tables.py` | post-pass on the generated LaTeX: applies the sizing decisions and repairs the float placement. |
| `build_pdfs.sh` | the one command: `fit_tables` → pandoc → `fix_tables` → `xelatex` ×2 → `doc/*.pdf`. Needs `pandoc >= 3` and `texlive`. |
| `number_audit.py` | the provenance check the paper claims: every decimal the manuscripts print must be a number an experiment artefact prints (leaf or status text), a power-of-ten multiple of one, a bibliography identifier, a section cross-reference, or arithmetic declared in its `DERIVED` table. Exits non-zero on a number nothing in the repository produced. |

```bash
python3 review/tools/build_series.py   # master -> doc/paper{1,2,3}_*.md
./review/tools/build_pdfs.sh           # -> doc/final_publication.pdf + the three companion PDFs
bash publication/build.sh              # the REVTeX build of the master (see publication/README.md)
python3 review/tools/number_audit.py   # every printed decimal traced to an artefact
```

The bibliography resolver is **not** here: `publication/fetch_dois.py` resolves
and verifies every entry against Crossref and writes the `doi.json` map that the
LaTeX filter turns into clickable links. Re-run it whenever the reference list
changes.
