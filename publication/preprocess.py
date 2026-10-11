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
        if len(body) < 80:
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
        continue
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


body_text = re.sub(
    r"^\*\*Table (\d+)\. (.+?)\.\*\*\n\n::: \{\.table id=tabx\d+\}\n\n((?:\|[^\n]*\n)+)\n:::\n",
    table_block, body_text, flags=re.M)

# the one caption-less table (the seven criteria) gets a real caption
body_text = body_text.replace(
    "::: {.table id=tabx1}\n\n",
    "::: {.table id=tabx1}\n\n**The seven criteria against which the paper "
    "judges itself.**\n\n", 1)

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
yaml = ["---", "title: |", "  " + title,
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
