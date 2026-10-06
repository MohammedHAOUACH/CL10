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
| `abstract.tex` | **override**: the 250-word journal abstract replaces the 410-word one in `doc/publication.md` |
| `code_availability.tex` | **override**: repository URL + reproduction commands, appended to Appendix B |
| `doi.json` | DOIs verified against Crossref, appended to the bibliography |
| `figures/` | copies of `doc/figures/*.png` used by the manuscript |

## Converter constraints found the hard way

- `\begin{document}` must precede `\title` with REVTeX 4-2.
- The `array` package conflicts with REVTeX 4-2 (`Extra \or` in every table).
- REVTeX 4-2 rejects `p{}` columns inside `longtable` (`Incomplete \iffalse`)
  → long columns are `\parbox`ed there, and `\caption` inside `longtable` is
  overridden by REVTeX → the Markdown "Table N." label becomes a header row.
- Section numbering is forced to arabic because the manuscript cross-references
  "Section 8.4" in prose.
- Figure captions drop the Markdown "Fig. N —" label: REVTeX prints `FIG. n.`.
