#!/usr/bin/env bash
# Build the journal PDF from doc/publication.md (read-only) using pandoc + REVTeX 4.2.
# Requires: pandoc >= 3, texlive (pdflatex), texlive-publishers (revtex4-2).
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p build
cp -r ../doc/figures build/ 2>/dev/null || true

python3 preprocess.py
pandoc build/publication_gen.md \
  --lua-filter=filter.lua \
  --template=template.tex \
  -s -o build/publication.tex

cd build
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode publication.tex >/dev/null
done
cd ..
cp build/publication.pdf publication.pdf

echo "OK: $(pwd)/publication.pdf"
pdfinfo publication.pdf | grep -E 'Pages|Page size'
