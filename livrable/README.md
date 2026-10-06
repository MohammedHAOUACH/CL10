# livrable/ — journal deliverable

Built from `doc/publication.md`, which the build **reads but never writes**.

## Build

```bash
python3 livrable/make_paper.py --no-compile     # -> livrable/publication.tex + figures/ copy
bash ~/.claude/skills/latex-document-skill/scripts/compile_latex.sh \
     livrable/publication.tex --engine xelatex --auto-fix --preview
```

`make_paper.py` also compiles by itself (`python3 livrable/make_paper.py`,
two `xelatex` passes); the skill's script is the preferred path because it
adds `--auto-fix` (float placement, `microtype`) and page previews.

REVTeX 4.2 is not in the system TeX Live on this machine. Installed once into
the user tree (`~/texmf/tex/latex/revtex/`) from CTAN:

```bash
curl -L -o /tmp/revtex.zip https://ctan.math.illinois.edu/macros/latex/contrib/revtex.zip
unzip -o /tmp/revtex.zip -d /tmp/revtex_pkg
cd /tmp/revtex_pkg/revtex && for d in ltxutil ltxgrid ltxfront ltxdocext revtex4-2 aip4-2; do tex $d.ins; done
mkdir -p ~/texmf/tex/latex/revtex && cp *.cls *.sty *.rtx ~/texmf/tex/latex/revtex/ && mktexlsr
```

## Files

| file | role |
|---|---|
| `make_paper.py` | Markdown → REVTeX 4.2 converter (stdlib only) |
| `publication.tex` | generated manuscript (`revtex4-2`, `aps,prd,preprint`) |
| `publication.pdf` | compiled paper, 26 figures, 32 references |
| `doi.json` | DOIs verified against Crossref, appended to the bibliography |
| `orcid.txt` | optional: the author ORCID, appended to the affiliation |
| `cover_letter.md` | submission letter, journal order, arXiv listing, referees |
| `figures/` | copies of `doc/figures/*.png` used by the manuscript |

The manuscript text lives in `doc/publication.md` and the build writes nothing
there: the abstract, the code-availability statement and the corrected
bibliography entries are in the source itself.

## Converter constraints found the hard way

- `\begin{document}` must precede `\title` with REVTeX 4-2.
- The `array` package conflicts with REVTeX 4-2 (`Extra \or` in every table).
- REVTeX 4-2 rejects `p{}` columns inside `longtable` (`Incomplete \iffalse`)
  → long columns are `\parbox`ed there, and `\caption` inside `longtable` is
  overridden by REVTeX → the Markdown "Table N." label becomes a header row.
- Section numbering is forced to arabic because the manuscript cross-references
  "Section 8.4" in prose.
- Figure captions drop the Markdown "Fig. N —" label: REVTeX prints `FIG. n.`.
- `p{}` widths are computed per table from the printed content
  (`column_widths`); a fixed width for every long column stacked the
  math-heavy cells into unreadable towers.
- `xurl` (not `url`) so the repository URL breaks; `\allowbreak{}` after `/`
  and `_` inside `\texttt{}` and after commas inside `$...$`, otherwise a long
  `\texttt{path}` or `\mathrm{diag}(...)` pushes a line past the margin.
- `floatfix` + `\topfraction`/`\totalnumber`: 26 figures in one section strand
  under the default TeX float limits.
- `$$ ... \\ ... $$` becomes `\[ \begin{gathered} ... \end{gathered} \]`;
  `\\` inside `\[ \]` alone is invalid.
- In a Python replacement string, `\a` is the BEL character: a LaTeX
  replacement must be written `\\allowbreak{}` (raw string is not enough).
- Nothing may sit between `\begin{longtable}{...}` and the `\multicolumn`
  caption row: `\footnotesize` there starts the first row, the `\multicolumn`
  becomes a **misplaced `\omit`**, and longtable then loops forever (a 600 s
  hang). Font and `\tabcolsep` go in a group *around* the table.
- Dense tables switch to `\scriptsize` when the widest cell needs more than two
  wrapped lines, and a prose column is never squeezed under 3.6 cm (mid-word
  hyphenation); the width is taken from the widest math column.
