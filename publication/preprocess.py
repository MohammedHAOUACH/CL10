#!/usr/bin/env python3
r"""Preprocess doc/final_publication.md into pandoc-friendly markdown for REVTeX.

Read-only on doc/. Writes build/publication_gen.md:
  - YAML metadata (title, abstract, keywords)
  - [n], [n,m], [n-m] citations -> \cite{refNN} raw spans (math masked)
  - Fig. N / Figs. N-M / Table N refs -> \ref{...} raw spans
  - figure blocks (**Fig. N — title.** + image + caption) -> .figure/.figurewide Divs
  - **Table N. title.** + pipe table -> .table Div
  - References ordered list -> .bibliography Div
"""
import math
import re
import struct
import sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "../doc/final_publication.md"
OUT = sys.argv[2] if len(sys.argv) > 2 else "build/publication_gen.md"

text = open(SRC, encoding="utf-8").read()

# ---------------------------------------------------------------- helpers
MATH_RE = re.compile(r"\$\$.*?\$\$|\$[^$\n]*\$", re.S)


def mask(s):
    store = []

    def r(m):
        store.append(m.group(0))
        return "\x00%d\x00" % (len(store) - 1)

    return MATH_RE.sub(r, s), store


def unmask(s, store):
    return re.sub(r"\x00(\d+)\x00", lambda m: store[int(m.group(1))], s)


CITE_RE = re.compile(r"\[(\d{1,2}(?:\s*(?:,|\u2013|-){1,2}\s*\d{1,2})*)\]")


def cite_repl(m):
    keys = []
    for part in m.group(1).split(","):
        part = part.strip().replace("\u2013", "-").replace("--", "-")
        if "-" in part:
            a, b = part.split("-")
            keys += ["ref%02d" % i for i in range(int(a), int(b) + 1)]
        else:
            keys.append("ref%02d" % int(part))
    return "``\\cite{%s}``" % ",".join(keys)


def figrefs(s):
    """Fig. 9 / Fig. 15 / Figs. 1-24 -> Fig.~\ref{fig:9} raw spans."""
    out = []
    i = 0
    pat = re.compile(r"\b(Figs?)\.(?=\s*\d)")
    while i < len(s):
        m = pat.search(s, i)
        if not m:
            out.append(s[i:])
            break
        out.append(s[i:m.start()])
        plural = m.group(1) == "Figs"
        j = m.end()
        refs = []
        while j < len(s):
            while j < len(s) and s[j] == " ":
                j += 1
            dm = re.match(r"(\d{1,2})\s*(?:\u2013|-)\s*(?=\d)", s[j:])
            if dm:
                j += dm.end()
                m2 = re.match(r"(\d{1,2})", s[j:])
                if not m2:
                    break
                a, b = int(dm.group(1)), int(m2.group(1))
                j += m2.end()
                refs += ["fig:%d" % k for k in range(a, b + 1)]
            else:
                nm = re.match(r"(\d{1,2})", s[j:])
                if not nm:
                    break
                refs.append("fig:%s" % nm.group(1))
                j += nm.end()
            # continue past ", " only if another number follows
            k = j
            while k < len(s) and s[k] == " ":
                k += 1
            if k < len(s) and s[k] == ",":
                k += 1
                while k < len(s) and s[k] == " ":
                    k += 1
                if re.match(r"\d", s[k:]):
                    j = k
                    continue
            break
        if not refs or len(refs) > 6:
            # too many (e.g. 'Figs. 1-24' in prose): leave the text as-is
            out.append(s[m.start():j])
            i = j
            continue
        label = "Figs.~" if plural and len(refs) > 1 else "Fig.~"
        out.append("``%s%s``" % (label, ", ".join("\\ref{%s}" % r for r in refs)))
        i = j
    return "".join(out)


def tablerefs(s):
    """Table 5 -> Table~\ref{tab:5}; '**Table N.' caption lines untouched."""
    return re.sub(r"(?<!\*\*)\bTable\s+(\d{1,2})\b",
                  lambda m: "``Table~\\ref{tab:%s}``" % m.group(1), s)


def xform(s):
    masked, store = mask(s)
    masked = CITE_RE.sub(cite_repl, masked)
    masked = figrefs(masked)
    masked = tablerefs(masked)
    s = unmask(masked, store)
    # long display equations -> equation + gathered (centered, one number).
    # multline was wrong here: it left-aligns the first line and
    # right-aligns the last, which reads as mis-positioned math.
    def multline(m):
        body = m.group(1).strip()
        if len(body) < 80 or "\\begin{" in body:
            # already laid out by the manuscript: do not re-wrap it
            return m.group(0)
        if re.search(r"\\oplus\s", body):
            # direct-sum chain: break after each \oplus (operator at EOL)
            chunks = [c.strip() for c in re.split(r"\s*\\oplus\s*", body)]
            if len(chunks) < 2:
                return m.group(0)
            lines = [c + r" {} \oplus" for c in chunks[:-1]] + [chunks[-1]]
        else:
            # system of conditions: break at each \qquad, centered
            chunks = [c.strip() for c in re.split(r"\\qquad\s*", body)]
            if len(chunks) < 2:
                return m.group(0)
            lines = chunks
        return ("\\begin{equation}\n\\begin{gathered}\n    "
                + "\n    \\\\\n".join(lines)
                + "\n\\end{gathered}\n\\end{equation}")

    s = re.sub(r"\$\$(.*?)\$\$", multline, s, flags=re.S)
    return s


def png_ratio(path):
    try:
        with open("../doc/" + path, "rb") as fh:
            fh.read(16)
            w, h = struct.unpack(">II", fh.read(8))
        return w / h
    except OSError:
        return 1.0


# ---------------------------------------------------------------- metadata
lines = text.split("\n")
title = None
abstract_lines, keyword_line = [], None
i = 0
while i < len(lines):
    ln = lines[i]
    if title is None and ln.startswith("# "):
        title = ln[2:].strip()
    elif ln == "## Abstract":
        i += 1
        while i < len(lines) and not lines[i].startswith("**Keywords.**"):
            if lines[i].strip():
                abstract_lines.append(lines[i].strip())
            i += 1
        continue
    elif ln.startswith("**Keywords.**"):
        keyword_line = ln[len("**Keywords.**"):].strip()
    i += 1
assert title and abstract_lines and keyword_line, "metadata extraction failed"

# ---------------------------------------------------------------- body
body = []
in_refs = False
bib_open = False
in_abstract = False
for ln in lines:
    if ln.startswith("# "):
        continue
    if ln == "## Abstract":
        in_abstract = True
        continue
    if in_abstract:
        if ln.startswith("**Keywords.**"):
            in_abstract = False
        continue
    if ln.startswith("**Mohammed Haouach**") or ln.startswith("*Independent researcher"):
        continue
    if ln.startswith("**Keywords.**"):
        continue
    if ln.strip() == "---":
        continue
    if ln == "## Figures":
        continue  # figures become floats; the section heading is dropped
    if ln.startswith("All figures are computed from the construction"):
        continue
    if ln == "## References":
        in_refs = True
        # thebibliography supplies its own heading; no markdown heading
        body.append("")
        body.append("::: {.bibliography}")
        body.append("")
        bib_open = True
        continue
    if in_refs and bib_open and ln.strip() and not re.match(r"^\d{1,2}\. ", ln):
        body.append("")
        body.append(":::")
        body.append("")
        bib_open = False
    if ln.startswith("### "):
        sub = re.sub(r"^\d+\.\d+\s+", "", ln[4:])
        body.append("## " + sub)  # ### -> ##, numbers auto by LaTeX
        continue
    if ln.startswith("## "):
        title2 = ln[3:]
        if title2.startswith("Appendix"):
            body.append("\\appendix")
            body.append("")
            title2 = re.sub(r"^Appendix [A-Z]\.\s+", "", title2)
            body.append("# " + title2)
            continue
        title2 = re.sub(r"^\d+\.\s+", "", title2)
        body.append("# " + title2)  # ## -> #  (pandoc: # -> \section)
        continue
    body.append(ln)
assert in_refs, "references section not found"
if bib_open:
    body.append("")
    body.append(":::")

body_text = "\n".join(body)

# ---------------------------------------------------------------- figures
# the paper places each figure as an image line immediately followed by its
# italic caption line; bind the pair into a float div. images wider than
# 1.8:1 span both columns (figure*), the rest take one column.
FIGPAIR = re.compile(
    r"^!\[[^\]]*\]\((figures/fig(\d+)_\w+\.png)\)\n\*([^\n]+)\*[ \t]*$", re.M)


def _ratio(img):
    # png_ratio resolves against ../doc/ and falls back to 1.0
    return png_ratio(img)


def _figdiv(m):
    img, num, cap = m.group(1), m.group(2), m.group(3).strip()
    cls = "figurewide" if _ratio(img) > 1.8 else "figure"
    return ("::: {.%s id=fig:%s src=%s}\n\n%s\n\n:::" % (cls, num, img, cap))


body_text = FIGPAIR.sub(_figdiv, body_text)

# ---------------------------------------------------------------- tables
# wrap every bare pipe-table block first, then merge "**Table N. title.**"
# captions into the wrapper div (avoids nesting divs)
_auto = [0]


def bare_table(m):
    _auto[0] += 1
    return ("::: {.table id=tabx%d}\n\n%s\n\n:::\n"
            % (_auto[0], m.group(0).rstrip("\n")))


body_text = re.sub(
    r"(?:^\|[^\n]*\n){2,}",
    bare_table, body_text, flags=re.M)


def table_block(m):
    num, ttitle, inner = m.group(1), m.group(2), m.group(3)
    # LaTeX supplies the number; keep only the descriptive title
    return ("::: {.table id=tab:%s}\n\n**%s**\n\n%s\n\n:::\n"
            % (num, ttitle, inner.strip()))


# the manuscript writes "**Table N.** title" (the number is LaTeX's job)
body_text = re.sub(
    r"^\*\*Table (\d+)\.\*\* (.+?)\n\n::: \{\.table id=tabx\d+\}\n\n((?:\|[^\n]*\n)+)\n:::\n",
    table_block, body_text, flags=re.M)
# the older form, caption inside the bold
body_text = re.sub(
    r"^\*\*Table (\d+)\. (.+?)\.\*\*\n\n::: \{\.table id=tabx\d+\}\n\n((?:\|[^\n]*\n)+)\n:::\n",
    table_block, body_text, flags=re.M)

# ---------------------------------------------------------------- table sizing
# one-column float when the table fits a column, a full-width float when it
# does not, and a full-width float plus a continuation float when it is
# taller than a page (a table* taller than \textheight spills past the
# bottom margin onto the next page)
CHAR, LINE, PAD = 4.1, 11.4, 6.0      # footnotesize pt/char, pt/line, pt/cell
COLW, TEXTW = 246.0, 510.0            # revtex PRD column and page width
MAXCOL, MAXFULL = 400.0, 560.0        # tallest one-column / full-width float

# the page-reserved group is set scriptsize with arraystretch 0.92; these are
# the measured scale (rendered rule-to-rule heights of a nine-row table) and
# the usable body height of one full-width float page: text area 672pt less
# caption, repeated header row and rules
CHART, LINET = 3.6, 7.0
MAXPAGE = 600.0


def _cells(row):
    return [c.strip() for c in row.strip().strip("|").split("|")]


def _plain(s):
    # mimic pandoc's stringify: markdown emphasis goes, TeX stays verbatim
    # (a formula is one long unbreakable word, exactly as LaTeX sees it)
    s = re.sub(r"\*\*|`|(?<!\\)_", "", s)
    return re.sub(r"[{}&\\]", "", s)


def _lines(text, wchars):
    """Lines a p{} column of wchars characters needs: real word wrapping,
    because a long unbreakable formula costs a whole line of its own."""
    n = 0
    for para in text.split("\n"):
        words = para.split()
        if not words:
            n += 1
            continue
        line, cur = 1, 0
        for w in words:
            lw = len(w)
            if lw > wchars:
                line += 2 if cur else 1
                cur = 0
            elif cur + lw > wchars:
                line += 1
                cur = lw + 1
            else:
                cur += lw + 1
        n += line
    return n


def _fracs(rows, unit, char):
    """Column fractions the filter will emit, or None if a word cannot fit."""
    ncol = len(_cells(rows[0]))
    maxlen = [0] * ncol
    maxword = [0] * ncol
    for r in rows:
        for j, c in enumerate(_cells(r)[:ncol]):
            p = _plain(c)
            maxlen[j] = max(maxlen[j], len(p))
            maxword[j] = max(maxword[j], max((len(w) for w in p.split()),
                                             default=0))
    fmin = [(maxword[j] * char + PAD) / unit for j in range(ncol)]
    if sum(fmin) > 0.97:
        return None
    total = sum(maxlen) or 1
    f = [max(maxlen[j] / total, fmin[j]) for j in range(ncol)]
    s = sum(f)
    if s > 0.97:
        flex = sum(max(0.0, f[j] - fmin[j]) for j in range(ncol)) or 1.0
        excess = s - 0.97
        f = [max(fmin[j], f[j] - excess * max(0.0, f[j] - fmin[j]) / flex)
             for j in range(ncol)]
    return f


def _rowh(r, f, unit, char, line):
    lines = 1
    for j, c in enumerate(_cells(r)[:len(f)]):
        avail = max(f[j] * unit - PAD, 20.0)
        lines = max(lines, _lines(_plain(c), max(int(avail / char), 4)))
    return lines * line + 4.0


def _height(rows, unit, char=CHAR, line=LINE):
    """Estimated float height in pt for the pipe rows at the given width."""
    f = _fracs(rows, unit, char)
    if f is None:
        return 1e9
    return sum(_rowh(r, f, unit, char, line) for r in rows)


TABLEDIV = re.compile(r"::: \{\.table id=([^}\s]+)\}\n\n((?:.|\n)*?)\n:::\n")


def size_table(m):
    tid, inner = m.group(1), m.group(2)
    cap = ""
    mm = re.match(r"(\*\*.*?\*\*)\n\n", inner)
    if mm:
        cap, inner = mm.group(1), inner[mm.end():]
    rows = [r for r in inner.split("\n") if r.strip()]
    if len(rows) < 3:
        return m.group(0)
    head, sep, body = rows[0], rows[1], rows[2:]
    if _height(rows, COLW) <= MAXCOL:
        cls = "table"
    elif _height(rows, TEXTW) <= MAXFULL:
        cls = "tablewide"
    else:
        # taller than one float: a page-reserved group. the body is split
        # into K full-width floats, each repeating the header row. K is the
        # smallest count that fits the page budget and the split is the one
        # that minimises the tallest float (DP over row boundaries): greedy
        # filling lands on a split a point over budget and buys a page the
        # table does not need, and a table* taller than the page is clipped
        # -- rows simply vanish
        f = _fracs(rows, TEXTW, CHART)
        hs = [_rowh(r, f, TEXTW, CHART, LINET) for r in body]
        n = len(hs)
        K = min(n, max(2, math.ceil(sum(hs) / MAXPAGE)))
        while True:
            inf = float("inf")
            dp = [[inf] * (n + 1) for _ in range(K + 1)]
            cut = [[n] * (n + 1) for _ in range(K + 1)]
            for i in range(n):
                dp[1][i] = sum(hs[i:])
            for k in range(2, K + 1):
                for i in range(n - 1, -1, -1):
                    acc = 0.0
                    for j in range(i, n - 1):
                        acc += hs[j]
                        cand = max(acc, dp[k - 1][j + 1])
                        if cand < dp[k][i]:
                            dp[k][i] = cand
                            cut[k][i] = j + 1
            if K >= n or dp[K][0] <= MAXPAGE:
                break
            K += 1
        chunks, at, k = [], 0, K
        while k > 1 and at < n:
            chunks.append(body[at:cut[k][at]])
            at = cut[k][at]
            k -= 1
        chunks.append(body[at:])
        out = ("::: {.tablepage id=%s}\n\n%s\n\n%s\n\n:::\n"
               % (tid, cap, "\n".join([head, sep] + chunks[0])))
        for k, ch in enumerate(chunks[1:], 1):
            out += ("\n::: {.tablecont .tablepage id=%s_cont%d}\n\n%s\n\n:::\n"
                    % (tid, k, "\n".join([head, sep] + ch)))
        return out
    return "::: {.%s id=%s}\n\n%s\n\n%s\n\n:::\n" % (cls, tid, cap,
                                                     "\n".join(rows))


body_text = TABLEDIV.sub(size_table, body_text)

# bare URLs -> \url{} so hyperref can break them (avoids overfull lines)
body_text = re.sub(
    r"(?<!\S)(https?://[^\s)]+)",
    lambda m: "``\\url{%s}``" % m.group(1), body_text)

# ---------------------------------------------------------------- citations/refs
body_text = xform(body_text)

# ---------------------------------------------------------------- place figures
# move each figure div right after the paragraph of its first textual
# reference (Fig.~\ref{fig:N}); a figure with no textual reference stays
# where the paper put it
figdiv_re = re.compile(
    r"::: \{\.figure(?:wide)? id=fig:(\d+) src=[^\n]*\}\n\n(?:.|\n)*?\n:::\n")
divs = {}


def grab(m):
    divs[m.group(1)] = (m.group(0), m.start())
    return ""


body_text = figdiv_re.sub(grab, body_text)

# descending: figures sharing one anchor paragraph end up in ascending
# order at that anchor (each insert lands before the previous one)
for num in sorted(divs, key=int, reverse=True):
    div, orig = divs[num]
    # the raw span may be singular ("``Fig.~\ref{fig:9}``") or plural
    # ("``Figs.~\ref{fig:7}, \ref{fig:8}``"): match the \ref itself
    pos = body_text.find("\\ref{fig:%s}" % num)
    if pos == -1:
        pos = min(orig, len(body_text))
    end = body_text.find("\n\n", pos)
    if end == -1:
        end = len(body_text)
    body_text = body_text[:end] + "\n\n" + div + body_text[end:]

# ---------------------------------------------------------------- emit
# a ": " in the title becomes a hard line break (backslash-newline in the
# literal block, pandoc's markdown hard break); titleplain keeps the one-line
# form for the PDF metadata
title_block = title
if ": " in title:
    head, rest = title.split(": ", 1)
    title_block = head + ":\\\n  " + rest
yaml = ["---", "title: |", "  " + title_block,
        'titleplain: "%s"' % title,
        "abstract: |", "  " + "\n\n  ".join(abstract_lines),
        "keywords: |", "  " + keyword_line, "---", ""]

open(OUT, "w", encoding="utf-8").write("\n".join(yaml) + "\n" + body_text + "\n")

print("wrote", OUT, len(body_text), "chars")
print("cite spans:", body_text.count("\\cite{"))
print("fig ref spans:", body_text.count("\\ref{fig:"))
print("tab ref spans:", body_text.count("\\ref{tab:"))
print("figure divs:", body_text.count("{.figure"))
print("table divs:", body_text.count("{.table"))
print("bibliography div:", body_text.count("{.bibliography"))
