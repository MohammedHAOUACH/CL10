#!/usr/bin/env python3
r"""Preprocess doc/publication.md into pandoc-friendly markdown for REVTeX.

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

SRC = sys.argv[1] if len(sys.argv) > 1 else "../doc/publication.md"
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
FIG_TITLE = re.compile(r"^\*\*Fig\. (\d+) \u2014 (.+)\.\*\*$", re.M)
IMG = re.compile(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$", re.M)

# split the figure region on Fig. headers
parts = FIG_TITLE.split(body_text)
# parts: [pre, num, title, mid, num, title, mid, ...]
if len(parts) > 1:
    rebuilt = [parts[0]]
    for k in range(1, len(parts), 3):
        num, ftitle, mid = parts[k], parts[k + 1], parts[k + 2]
        mi = IMG.search(mid)
        if not mi:
            rebuilt.append("**Fig. %s \u2014 %s.**\n%s" % (num, ftitle, mid))
            continue
        img = mi.group(2)
        cap_lines = []
        leftover = []
        seen_img = False
        for r in mid.split("\n"):
            if not seen_img:
                if r.strip() == mi.group(0).strip():
                    seen_img = True
                    continue
            if seen_img:
                if (r.strip().startswith("---") or r.strip().startswith("#")
                        or r.strip().startswith("\\appendix")):
                    leftover.append(r)
                    continue
                if leftover:
                    leftover.append(r)
                    continue
            if r.strip():
                cap_lines.append(r.strip())
        wide = png_ratio(img) > 1.8
        cls = "figurewide" if wide else "figure"
        cap = "**%s.** %s" % (ftitle, "\n\n".join(cap_lines))
        rebuilt.append("::: {.%s id=fig:%s src=%s}\n\n%s\n\n:::\n%s\n"
                       % (cls, num, img, cap, "\n".join(leftover)))
    body_text = "".join(rebuilt)

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

# two figures are never cited in the body (journal rule: every figure must
# be referenced); add the citations at their natural anchor sentences
body_text = body_text.replace(
    "\u2014 a ratio $1.12$.",
    "\u2014 a ratio $1.12$ (Fig. 8).", 1)
body_text = body_text.replace(
    "It should not be read as a reduction of the theory's difficulty.",
    "It should not be read as a reduction of the theory's difficulty "
    "(``Fig.~\\ref{fig:18}``).", 1)

# bare URLs -> \url{} so hyperref can break them (avoids overfull lines)
body_text = re.sub(
    r"(?<!\S)(https?://[^\s)]+)",
    lambda m: "``\\url{%s}``" % m.group(1), body_text)

# ---------------------------------------------------------------- citations/refs
body_text = xform(body_text)

# ---------------------------------------------------------------- place figures
# move each figure div right after the paragraph of its first textual
# reference (Fig.~\ref{fig:N}); unreferenced figures stay at the end
figdiv_re = re.compile(
    r"::: \{\.figure(?:wide)? id=fig:(\d+) src=[^\n]*\}\n\n(?:.|\n)*?\n:::\n")
divs = {}
def grab(m):
    divs[m.group(1)] = m.group(0)
    return ""
body_text = figdiv_re.sub(grab, body_text)

for num in sorted(divs, key=int):
    # fig 26: its caption block swallowed the "Provenance" paragraph, making
    # a ~515pt float that fits nowhere; keep only the figure description in
    # the caption, return the provenance paragraph to the appendix prose,
    # and shrink the image so the float lands on the page with free space.
    if num == "26" and "Provenance of" in divs[num]:
        i = divs[num].find("**Provenance of")
        j = divs[num].find("\n:::", i)
        prov = divs[num][i:j].strip()
        divs[num] = (divs[num][:i].rstrip() + "\n\n:::\n")
        divs[num] = divs[num].replace(
            "{.figure id=fig:26", '{.figure id=fig:26 width="0.85\\columnwidth"', 1)
        divs[num + "_prov"] = prov + "\n"
    marker = "``Fig.~\\ref{fig:%s}``" % num
    if marker not in body_text:
        marker = "``Figs.~"  # range reference: place at first range mention
    # fig 26 is cited in the code-availability paragraph at the very end
    # of Appendix B: declaring it there forces the float past the
    # bibliography. Anchor it at the start of Appendix B instead so it
    # floats within the appendix pages, before REFERENCES.
    if num == "26":
        anchor = body_text.find("\\appendix\n\n# Numerical protocol")
        if anchor != -1:
            end = body_text.find("\n\n", anchor)
            body_text = (body_text[:end] + "\n\n" + divs[num]
                         + ("\n" + divs.get("26_prov", ""))
                         + body_text[end:])
            continue
    pos = body_text.find(marker)
    if pos == -1:
        body_text += "\n" + divs[num]
        continue
    # end of the paragraph containing the reference
    end = body_text.find("\n\n", pos)
    if end == -1:
        end = len(body_text)
    body_text = body_text[:end] + "\n\n" + divs[num] + body_text[end:]

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
