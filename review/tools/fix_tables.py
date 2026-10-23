#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Wrap longtable cells in unbreakable \parbox so rows page-break atomically
instead of splitting mid-sentence across pages.

Usage: fix_tables.py file.tex   (edits in place)

Pandoc's longtable rows end with '\\' + newline but are soft-wrapped across
physical lines, and a cell boundary '&' can land at a line start. So: split
the header (between \toprule\noalign{} and \midrule) and the body (after
\endlastfoot) into logical rows at the row terminator, flatten each row to
one line (LaTeX treats newlines as spaces), split cells on ' & ', and wrap
each cell in \parbox[t]{\linewidth}{...}. Rows whose cells already use
minipage machinery pass through untouched.
"""
import re
import sys

LT = re.compile(r"(\\begin\{longtable\}.*?\\end\{longtable\})", re.S)
SKIPCELL = re.compile(r"minipage|parbox|multicolumn")
ROWEND = re.compile(r"\\\\[ \t]*\n")


def wrap_cells(row):
    if "&" not in row or SKIPCELL.search(row):
        return row
    cells = row.split(" & ")
    out = []
    for c in cells:
        s = c.strip()
        out.append("\\parbox[t]{\\linewidth}{" + s + "}" if s else s)
    return " & ".join(out)


def process(text):
    rows = ROWEND.split(text)
    out = []
    for r in rows:
        if not r.strip() or SKIPCELL.search(r):
            out.append(r)
        else:
            out.append(wrap_cells(" ".join(r.split("\n"))))
    return "\\\\\n".join(out)


def fix_block(m):
    block = m.group(1)
    block = re.sub(r"(?s)(?<=\\toprule\\noalign\{\}\n)(.*?)(?=\\midrule)",
                   lambda mm: process(mm.group(1)) + "\n", block)
    parts = block.split("\\endlastfoot", 1)
    if len(parts) == 2:
        body = parts[1]
        end = body.rindex("\\end{longtable}")
        body = "\\endlastfoot\n" + process(body[:end].lstrip("\n")) \
            + body[end:]
        block = parts[0] + body
    return block


def main(path):
    src = open(path, encoding="utf-8").read()
    open(path, "w", encoding="utf-8").write(LT.sub(fix_block, src))


if __name__ == "__main__":
    main(sys.argv[1])
