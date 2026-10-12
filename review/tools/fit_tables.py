#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Fit pipe-table column widths to their content (pandoc honours the
relative widths of the separator row for wide tables).

Usage: fit_tables.py file.md [file.md ...]   (edits in place)

For each column j, on markup-stripped text: L_j = longest body cell,
wl_j = longest word. Target percent p_j = max(sqrt(L_j), (wl_j + 4)/0.9)
— sqrt(L) balances row heights across columns, the second term keeps the
longest word on one line (~90 chars per full text width at 10pt).
Normalised to 100 dashes in the separator row (alignment colons kept).
"""
import re
import sys

SEP = re.compile(r"^\|[\s:|-]+\|[-\s:|]*$")


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]
def visible(s):
    s = re.sub(r"\\(?:mathrm|mathbf|text|operatorname|tfrac|bar|hat)\{", "{", s)
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    s = re.sub(r"\*\*|\*|`", "", s)
    return re.sub(r"[$&^_{}\\]", "", s)


def fit_block(block):
    header, sep, *body = block
    ncol = len(cells(header))
    L, wl = [], []
    for j in range(ncol):
        lens, words = [0], [1]
        for row in body:
            cs = cells(row)
            if j < len(cs):
                v = visible(cs[j])
                lens.append(len(v))
                words += [len(w) for w in v.split() if w]
        L.append(max(lens))
        wl.append(max(words))
    w = [L[j] ** 0.5 for j in range(ncol)]
    floor = [(wl[j] + 4) / 0.9 for j in range(ncol)]
    for _ in range(ncol):
        total = sum(w)
        w = [100 * x / total for x in w]
        if all(w[j] >= floor[j] for j in range(ncol)):
            break
        w = [max(w[j], floor[j]) for j in range(ncol)]
    dashes = [max(3, round(x)) for x in w]
    align = cells(sep)
    parts = []
    for j in range(ncol):
        a = align[j] if j < len(align) else "---"
        left = ":" if a.startswith(":") else ""
        right = ":" if a.endswith(":") else ""
        parts.append(left + "-" * dashes[j] + right)
    newsep = "|" + "|".join(parts) + "|"
    return [header, newsep] + body


def main(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    out, i, n = [], 0, len(lines)
    while i < n:
        if (lines[i].startswith("|") and i + 1 < n and SEP.match(lines[i + 1])
                and set(lines[i + 1]) <= set("|-: ")):
            block = [lines[i], lines[i + 1]]
            j = i + 2
            while j < n and lines[j].startswith("|"):
                block.append(lines[j])
                j += 1
            out.extend(fit_block(block))
            i = j
        else:
            out.append(lines[i])
            i += 1
    open(path, "w", encoding="utf-8").write("\n".join(out))


if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
        print("fitted:", p)
