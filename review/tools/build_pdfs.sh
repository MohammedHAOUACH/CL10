#!/usr/bin/env bash
# Build the submission PDFs: pandoc -> standalone LaTeX -> atomic-row table
# fix (fix_tables.py) -> xelatex x2. Run from the repository root.
set -euo pipefail
cd "$(dirname "$0")/../../doc"
python3 ../review/tools/fit_tables.py paper1_core.md paper2_gauge_gravity.md \
  paper3_phenomenology.md final_publication.md
for f in paper1_core paper2_gauge_gravity paper3_phenomenology final_publication; do
  pandoc "$f.md" -o "$f.tex" --standalone --pdf-engine=xelatex \
    -V geometry:margin=0.9in -V fontsize=10pt -V colorlinks -V linkcolor=blue -V lang=en-US
  python3 ../review/tools/fix_tables.py "$f.tex"
  xelatex -interaction=nonstopmode "$f.tex" >/dev/null
  xelatex -interaction=nonstopmode "$f.tex" >/dev/null
  echo "$f: $(pdfinfo "$f.pdf" | awk '/^Pages/{print $2}') pages"
done
