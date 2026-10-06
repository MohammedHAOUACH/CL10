#!/usr/bin/env python3
r"""Build the journal deliverable from the Markdown manuscript.

Reads ``doc/publication.md`` (READ-ONLY: this script never writes to ``doc/``)
and produces, in ``livrable/``:

  publication.tex   REVTeX 4.2 source (``revtex4-2``, APS/PRD options)
  refs.bib          optional (only with --bib); default uses \thebibliography
  figures/*.png     copies of the manuscript figures
  publication.pdf   compiled with xelatex (twice, for the cross-references)

Editorial inputs kept here, not in ``doc/``:

  livrable/doi.json       verified DOIs, keyed by reference number; appended
                          to the bibliography entries as \url{https://doi.org/..}
  livrable/orcid.txt      if present, the author ORCID, appended to the affiliation

Usage:
  python3 livrable/make_paper.py [--twocolumn] [--no-compile]
"""
import json
import os
import re
import shutil
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC_MD = os.path.join(ROOT, "doc", "publication.md")
FIG_SRC = os.path.join(ROOT, "doc", "figures")

REF_KEYS = {n: f"ref{n:02d}" for n in range(1, 100)}

# text-mode characters that LaTeX needs escaped
ESCAPES = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "#": r"\#",
           "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
           "^": r"\textasciicircum{}"}

CITE_RE = re.compile(r"\[(\d{1,2}(?:\s*(?:,|\u2013|-){1,2}\s*\d{1,2})*)\]")
UNICODE = {"\u2014": "---", "\u2013": "--", "\u2018": "`", "\u2019": "'",
           "\u201c": "``", "\u201d": "''", "\u00a7": r"\S", "\u00d7": r"$\times$",
           "\u2192": r"$\rightarrow$", "\u2248": r"$\approx$", "\u2265": r"$\ge$",
           "\u2264": r"$\le$", "\u03a3": r"$\Sigma$"}


def escape_text(s):
    out = []
    for ch in s:
        out.append(ESCAPES.get(ch, UNICODE.get(ch, ch)))
    return "".join(out)




def citations(s):
    """[12], [25,26] and [1-5] -> \\cite{ref12}, ... (ranges expanded)."""
    def repl(m):
        keys = []
        for part in m.group(1).split(","):
            part = part.strip()
            part = re.sub(r"[\u2013-]{1,2}", "-", part)
            if "-" in part:
                a, b = part.split("-")
                keys += [REF_KEYS[int(i)] for i in range(int(a), int(b) + 1)]
            else:
                keys.append(REF_KEYS[int(part)])
        return "\\cite{" + ",".join(keys) + "}"
    return CITE_RE.sub(repl, s)


def inline(s, mbox_math=False):
    """Convert Markdown text to LaTeX; $...$ math is stashed, never rewritten.

    mbox_math=True additionally shields each inline formula in \\mbox{} so it
    cannot be split across lines (table cells: no orphaned symbols).
    """
    maths = []

    def stash(m):
        # short formulas are atomic: \mbox keeps them on one line (no orphaned
        # symbols in table cells); long ones keep the comma-break behaviour
        if mbox_math and rendered_len(m.group(0)) <= 14:
            maths.append("\\mbox{" + m.group(0) + "}")
            return "\x00%d\x00" % (len(maths) - 1)
        # a long unbreakable inline expression (diag(...), charge lists) is what
        # pushes a paragraph past the margin: permit a break after each comma
        maths.append(re.sub(r",(?=\\|\s|\+|-|[A-Za-z0-9(])", r",\\allowbreak{}", m.group(0)))
        return "\x00%d\x00" % (len(maths) - 1)

    if mbox_math:
        # a number glued to its unit ($X$ eV, $4\times10^{47}$ yr) must not
        # wrap across lines: stash the whole pair as one raw-LaTeX \mbox
        # (inserting \mbox directly would be mangled by escape_text below)
        def unit_pair(m):
            maths.append("\\mbox{" + m.group(1) + " " + m.group(2) + "}")
            return "\x00%d\x00" % (len(maths) - 1)
        s = re.sub(r"(\$[^$]{1,16}\$)\s(GeV|MeV|TeV|keV|eV|yr)\b",
                   unit_pair, s)

    s = re.sub(r"\$\$[^$]*\$\$|\$[^$]*\$", stash, s)
    s = escape_text(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"\\emph{\1}", s)
    def code(m):
        body = re.sub(r"(/|\\_|_)", r"\1\\allowbreak{}", m.group(1))
        return "\\texttt{" + body + "}"
    s = re.sub(r"`([^`]+)`", code, s)
    s = re.sub(r"https?://[^\s)\]]+", r"\\url{\g<0>}", s)
    s = citations(s)
    return re.sub(r"\x00(\d+)\x00", lambda m: maths[int(m.group(1))], s)


def col_spec(rows):
    """l for short columns, p{3.4cm} where the text is long."""
    ncol = max(len(r) for r in rows)
    spec = []
    for c in range(ncol):
        longest = max((len(r[c]) for r in rows if c < len(r)), default=0)
        spec.append("p{3.4cm}" if longest > 42 else "l")
    return "".join(spec)


def split_row(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def is_sep(line):
    return bool(re.match(r"^\|[\s:|-]+\|$", line.strip()))


def rendered_len(cell):
    """Approximate the printed width of a Markdown cell (markup costs nothing)."""
    s = re.sub(r"\\[a-zA-Z]+", "", cell)          # \mathrm, \tfrac, \alpha -> printed glyphs only
    s = re.sub(r"[\$&{}\\]", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    return len(s)


CM_PER_PT = 0.0351455
CHAR_CM = 0.19          # average glyph width at \footnotesize


def prose_share(cell):
    """Fraction of the printed characters that are text rather than math."""
    tot = rendered_len(cell)
    if not tot:
        return 1.0
    maths = sum(rendered_len(m) for m in re.findall(r"\$[^$]*\$", cell))
    return 1.0 - maths / tot


def column_widths(rows, ncol):
    """Natural width per column, then a share of \\textwidth for the wide ones.

    The share is proportional to the content, except that a prose column is
    never squeezed under 3 cm: mid-word hyphenation in a 2 cm column is what
    made the residual tables unreadable. The deficit is taken from the widest
    column, which is a math column and wraps badly at any width anyway.
    """
    avail = LAYOUT[0] * CM_PER_PT - 0.45 * ncol   # \\tabcolsep: 6pt on each side
    longest = [max((r[c] for r in rows if c < len(r)), key=rendered_len, default="")
               for c in range(ncol)]
    nat = [rendered_len(longest[c]) * CHAR_CM + 0.5 for c in range(ncol)]
    if sum(nat) <= avail:
        return [None] * ncol, max(nat)             # everything fits unwrapped
    narrow = [c for c in range(ncol) if nat[c] < 3.0]
    wide = [c for c in range(ncol) if nat[c] >= 3.0]
    if not wide:
        wide, narrow = list(range(ncol)), []
    budget = max(1.5, avail - sum(nat[c] for c in narrow))
    share = sum(nat[c] for c in wide)
    out = {c: max(2.2, budget * nat[c] / share) for c in wide}
    total = sum(out.values())
    if total > budget:                             # scale the proportional share
        for c in wide:
            out[c] = out[c] * budget / total
    # a prose column under ~3.6 cm hyphenates mid-word; take the width from the
    # widest column, which is a math column and wraps badly at any width.
    # Only columns whose own content is genuinely long prose qualify: a
    # single-word label column must keep its natural width.
    for c in wide:
        if out[c] < 3.6 and nat[c] >= 4.5 and prose_share(longest[c]) > 0.6:
            others = [d for d in wide if d != c]
            donors = [d for d in others if prose_share(longest[d]) < 0.6] or others
            donor = max(donors, key=lambda d: out[d])
            need = min(3.6 - out[c], out[donor] - 2.2)
            out[c] += need
            out[donor] -= need
    # no unbreakable unit may overflow its column: floor each width at the
    # longest single word OR short inline formula (formulas <= 14 glyphs are
    # \mbox'ed in cells, so they are unbreakable too), scanned over ALL cells
    for c in wide:
        units = []
        for r in rows:
            cell = r[c] if c < len(r) else ""
            units += [rendered_len(m) for m in re.findall(r"\$[^$]*\$", cell)
                      if rendered_len(m) <= 14]
            units += [rendered_len(w) for w in
                      re.sub(r"\$[^$]*\$", " ", cell).split()]
        word_min = (max(units, default=0) * CHAR_CM + 0.35)
        if out[c] < word_min:
            donors = [d for d in wide if d != c and out[d] > 2.6]
            if donors:
                donor = max(donors, key=lambda d: out[d])
                out[donor] -= word_min - out[c]
                out[c] = word_min
    return [None if c not in out else round(out[c], 2) for c in range(ncol)], max(nat)


def table_block(lines, caption=None):
    """Short tables float as a tabular; long ones become a longtable.

    Column widths are proportional to the printed content: a fixed 3.4 cm for
    every long column stacked the math-heavy cells into unreadable towers.
    REVTeX 4-2 rejects p{} columns inside longtable, so there the wrapped
    columns are \\parbox{\\dots} of the same computed width.
    """
    rows = [split_row(l) for l in lines if not is_sep(l)]
    ncol = max(len(r) for r in rows)
    widths, tallest = column_widths(rows, ncol)
    # a cell needing more than two wrapped lines at \footnotesize is more
    # readable at \scriptsize than as a six-line tower
    font = ("scriptsize" if tallest > 2 * (LAYOUT[0] * CM_PER_PT - 0.45 * ncol)
            else "footnotesize")
    pad = ("\\setlength{\\tabcolsep}{4pt}\n" if ncol >= 4 else "")

    def row(r, wrap, bold=False):
        cells = []
        for i, c in enumerate(r):
            if wrap and i < ncol and widths[i]:
                # \\raggedright: no hyphenation and no stretched spaces inside
                # cells; words stay whole, the right edge stays ragged.
                # [c]: centre short cells vertically so row rules read cleanly.
                # \\vspace top/bottom: the text must breathe inside the cell —
                # never flush against the row rules above or below.
                t = inline(c, mbox_math=True)
                if bold:
                    t = r"\textbf{" + t + "}"
                t = (r"\parbox[c]{" + f"{widths[i]}cm" + "}{\\vspace{2.5pt}"
                     "\\raggedright \\hyphenpenalty=10000000 "
                     + t + "\\vspace{2.5pt}}")
            else:
                t = inline(c)
                if bold:
                    t = r"\textbf{" + t + "}"
            cells.append(t)
        return " & ".join(cells) + r" \\"

    def caption_text():
        # APS style: bold "Table N." label, plain caption text
        if not caption:
            return ""
        m = re.match(r"^(Table \d+\.)(.*)$", caption.strip())
        return (r"\textbf{" + m.group(1) + "}" + inline(m.group(2))
                if m else inline(caption))

    def caption_row():
        return ("\\multicolumn{" + str(ncol) + "}{l}{" + caption_text() + "} \\\\\n"
                + "\\addlinespace[2pt]\n")

    cap_tab = (caption_text() + r" \\[4pt]" + "\n") if caption else ""
    stretch = "\\renewcommand{\\arraystretch}{1.45}\n"
    # light rule between body rows: multi-line parbox cells otherwise read as
    # one grey block; a faint rule keeps each row visually distinct
    sep = "\n\\arrayrulecolor{black!40}\\hline\\arrayrulecolor{black}\n"

    if len(rows) > 8 or any(w for w in widths):
        # the Markdown caption already carries the "Table N." label, so it is
        # typeset as a header row: REVTeX overrides longtable's \caption.
        # The caption sits in the FIRST head only (\endfirsthead); continuation
        # pages repeat just the bold column headers.
        head = row(rows[0], True, bold=True)
        rest = sep.join(row(r, True) for r in rows[1:])
        # \footnotesize / \arraystretch must NOT sit between \begin{longtable}
        # and the \multicolumn caption row: that starts the first row and the
        # \multicolumn then becomes a misplaced \omit (and longtable loops).
        return ("{\\" + font + "\n" + pad + stretch + "\n"
                + "\\begin{longtable}{" + "l" * ncol + "}" + "\n" + caption_row()
                + "\\toprule" + "\n" + head + "\n" + "\\midrule" + "\n"
                + "\\endfirsthead" + "\n"
                + "\\toprule" + "\n" + head + "\n" + "\\midrule" + "\n"
                + "\\endhead" + "\n" + "\\midrule" + "\n" + "\\endfoot" + "\n"
                + "\\bottomrule" + "\n" + "\\endlastfoot" + "\n"
                + rest + "\n" + "\\end{longtable}" + "\n}\n")

    head = row(rows[0], False, bold=True)
    rest = sep.join(row(r, False) for r in rows[1:])
    spec = "".join("p{" + f"{widths[c]}cm" + "}" if widths[c] else "l"
                   for c in range(ncol))
    return (r"\begin{table}[htb]" + "\n" + cap_tab
            + r"\begin{center}" + "\\" + font + "\n" + pad + stretch
            + r"\begin{tabular}{" + spec + "}" + "\n"
            + r"\toprule" + "\n" + head + "\n" + r"\midrule" + "\n"
            + rest + "\n" + r"\bottomrule" + "\n"
            + r"\end{tabular}" + "\n" + r"\end{center}" + "\n"
            + r"\end{table}" + "\n")


def png_size(path):
    with open(path, "rb") as f:
        head = f.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return tuple(struct.unpack(">II", head[16:24]))


def measure_layout(opts):
    """Ask the document class for its own text width/height (no magic numbers)."""
    probe, log = os.path.join(HERE, ".layout.tex"), os.path.join(HERE, ".layout.log")
    open(probe, "w", encoding="utf-8").write(
        "\\documentclass[" + opts + "]{revtex4-2}\n\\begin{document}\nx\n"
        "\\typeout{MEASURE tw=\\the\\textwidth th=\\the\\textheight}\n"
        "\\end{document}\n")
    subprocess.run(["xelatex", "-interaction=nonstopmode", ".layout.tex"],
                   cwd=HERE, capture_output=True, text=True)
    txt = open(log, encoding="utf-8", errors="replace").read()
    for f in (probe, log, os.path.join(HERE, ".layout.aux")):
        if os.path.exists(f):
            os.remove(f)
    m = re.search(r"MEASURE tw=([\d.]+)pt th=([\d.]+)pt", txt)
    return (float(m.group(1)), float(m.group(2))) if m else (468.0, 665.5)


LAYOUT = (468.0, 665.5)          # overwritten by main() with the class's own values


def figure_block(path, caption):
    name = os.path.basename(path)
    # REVTeX prints "FIG. n." itself: drop the Markdown "Fig. N - " label
    if caption:
        caption = list(caption)
        caption[0] = re.sub(r"^Fig\.\s*\d+\s*[\u2014\u2013-]{1,2}\s*", "", caption[0])
    cap = " ".join(inline(c) for c in caption if c)
    tw, th = LAYOUT
    dims = png_size(os.path.join(ROOT, "doc", path))
    width = tw if not dims else min(tw, 0.72 * th * dims[0] / dims[1])
    return (r"\begin{figure}[htbp]" + "\n"
            + f"\\includegraphics[width={width:.0f}pt]{{{name}}}" + "\n"
            + r"\caption{" + cap + "}" + "\n"
            + r"\label{fig:" + name + "}" + "\n"
            + r"\end{figure}" + "\n")


def strip_number(title):
    return re.sub(r"^(?:Appendix\s+[A-Z]\.?\s+|\d+(?:\.\d+)*\.?\s+)", "", title).strip()


def build_body(lines, opts):
    out, i, n = [], 0, len(lines)
    pending_caption = None          # **Fig. N — ..** / **Table N. ..** line
    in_appendix = False
    while i < n:
        raw = lines[i].rstrip()
        s = raw.strip()

        if not s or s == "---":
            i += 1
            continue

        m = re.match(r"^(#{2,4})\s+(.*)$", s)
        if m:
            depth = len(m.group(1))
            title = inline(strip_number(m.group(2)))
            if m.group(2).startswith("Appendix") and depth == 2 and not in_appendix:
                out.append("\\appendix")
                in_appendix = True
            cmd = {2: "section", 3: "subsection", 4: "subsubsection"}[depth]
            out.append(f"\\{cmd}{{{title}}}\n")
            i += 1
            continue

        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)$", s)
        if m:
            path = m.group(2)
            cap = [pending_caption] if pending_caption else []
            # the caption paragraph follows the image
            j = i + 1
            while j < n and not lines[j].strip() == "":
                cap.append(lines[j].strip())
                j += 1
            out.append(figure_block(path, [c for c in cap if c]))
            pending_caption = None
            i = j
            continue

        if s.startswith("**") and s.endswith("**") and re.match(r"^\*\*(Fig\.|Table)", s):
            pending_caption = s.strip("* ")
            i += 1
            # a table caption is attached to the table that follows
            j = i
            while j < n and not lines[j].strip():
                j += 1
            if j < n and lines[j].strip().startswith("|"):
                tbl = []
                while j < n and lines[j].strip().startswith("|"):
                    tbl.append(lines[j])
                    j += 1
                out.append(table_block(tbl, pending_caption))
                pending_caption = None
                i = j
                continue
            continue

        if s.startswith("|"):
            tbl = []
            while i < n and lines[i].strip().startswith("|"):
                tbl.append(lines[i])
                i += 1
            cap = (r"\caption{" + inline(pending_caption) + "}\n"
                   if pending_caption else "")
            pending_caption = None
            out.append(table_block(tbl, pending_caption))
            continue

        if s.startswith("$$") and s.endswith("$$"):
            body = s[2:-2].strip()
            if "\\\\" in body:      # a multi-line display needs an alignment env
                out.append("\\[\n\\begin{gathered}\n" + body + "\n\\end{gathered}\n\\]\n")
            else:
                out.append("\\[" + body + "\\]\n")
            i += 1
            continue

        # paragraph: consecutive non-blank, non-special lines
        para = []
        while i < n:
            t = lines[i].strip()
            if (not t or t == "---" or t.startswith("#") or t.startswith("|")
                    or t.startswith("![") or t.startswith("$$")
                    or (t.startswith("**") and t.endswith("**")
                        and re.match(r"^\*\*(Fig\.|Table)", t))):
                break
            para.append(t)
            i += 1
        if para:
            out.append(inline(" ".join(para)) + "\n")
            pending_caption = None

    return "\n".join(out)


def parse_references(block):
    """The numbered Markdown list -> \\bibitem entries (faithful, no re-parsing)."""
    dois = {}
    p = os.path.join(HERE, "doi.json")
    if os.path.exists(p):
        for k, v in json.load(open(p)).items():
            if v and k.isdigit():
                dois[int(k)] = v
    items = []
    for line in block:
        m = re.match(r"^(\d{1,2})\.\s+(.*)$", line.strip())
        if not m:
            if items:
                items[-1][1] += " " + line.strip()      # wrapped entry
            continue
        items.append([int(m.group(1)), m.group(2)])
    out = [r"\begin{thebibliography}{99}" + "\n"]
    for n, text in items:
        doi = dois.get(n)
        tail = (" \\url{https://doi.org/" + doi + "}" if doi else "")
        out.append("\\bibitem{" + REF_KEYS[n] + "} " + inline(text) + tail + "\n")
    out.append(r"\end{thebibliography}" + "\n")
    return "".join(out)


def main():
    argv = sys.argv[1:]
    twocolumn = "--twocolumn" in argv
    compile = "--no-compile" not in argv

    md = open(SRC_MD, encoding="utf-8").read().split("\n")

    title = md[0].lstrip("# ").strip()
    author = next(l.strip().strip("*") for l in md[1:8] if l.strip().startswith("**"))
    affil = next(l.strip().strip("*") for l in md[1:8] if l.strip().startswith("*")
                 and not l.strip().startswith("**"))

    # split: abstract block, keywords line, body, references block
    def idx(pat):
        return next(i for i, l in enumerate(md) if re.match(pat, l))
    i_abs, i_ref = idx(r"^## Abstract"), idx(r"^## References")
    i_kw = next(i for i in range(i_abs, i_ref) if md[i].startswith("**Keywords"))
    abstract_md = [l for l in md[i_abs + 1:i_kw] if l.strip()]
    keywords = md[i_kw]
    body_lines = md[i_kw + 2:i_ref]
    refs_block = md[i_ref + 1:]

    abstract = "\n\n".join(inline(" ".join(p.split("\n")))
                           for p in " ".join(abstract_md).split("\n\n"))
    abstract_src = "doc/publication.md"

    doc = []
    opts = "aps,prd,amsmath,amssymb,superscriptaddress,longbibliography,floatfix"
    opts += ",twocolumn" if twocolumn else ",preprint"
    global LAYOUT
    LAYOUT = measure_layout(opts)
    doc.append("\\documentclass[" + opts + "]{revtex4-2}\n")
    doc.append("\\usepackage{graphicx}\n\\usepackage{booktabs}\n"
               "\\usepackage{longtable}\n\\usepackage{microtype}\n"
               "\\usepackage{xcolor}\n\\usepackage{colortbl}\n"
               "\\usepackage{ragged2e}\n")   # faint grey row rules; ragged cells
    doc.append("\\usepackage{xurl}\n")   # no 'array': REVTeX 4-2 redefines the column machinery
    doc.append("\\graphicspath{{figures/}}\n\n")
    # 26 figures in one section: the default TeX float limits strand them
    doc.append("\\renewcommand{\\topfraction}{0.95}\n"
               "\\renewcommand{\\bottomfraction}{0.85}\n"
               "\\renewcommand{\\textfraction}{0.05}\n"
               "\\renewcommand{\\floatpagefraction}{0.7}\n"
               "\\setcounter{totalnumber}{4}\n\n")
    # the manuscript cross-references sections by arabic number ("Section 8.4")
    doc.append("\\renewcommand{\\thesection}{\\arabic{section}}\n"
               "\\renewcommand{\\thesubsection}{\\arabic{section}.\\arabic{subsection}}\n"
               "\\renewcommand{\\thesubsubsection}{\\arabic{section}.\\arabic{subsection}"
               ".\\arabic{subsubsection}}\n\n")
    doc.append("\\begin{document}\n\n")
    doc.append("\\title{" + inline(title) + "}\n")
    doc.append("\\author{" + inline(author) + "}\n")
    orcid = ""
    op = os.path.join(HERE, "orcid.txt")
    if os.path.exists(op):
        orcid = open(op, encoding="utf-8").read().strip()
    doc.append("\\affiliation{" + inline(affil)
               + ("; ORCID: " + orcid if orcid else "") + "}\n\n")
    doc.append("\\begin{abstract}\n" + abstract + "\n")
    if keywords:
        doc.append("\n\\medskip\n" + inline(keywords.replace("**Keywords.**",
                   "*Keywords.*")) + "\n")
    doc.append("\\end{abstract}\n\n\\maketitle\n\n")
    doc.append(build_body(body_lines, opts))
    doc.append("\n\\clearpage\n\\bibliographystyle{apsrev4-2}\n")
    doc.append(parse_references(refs_block))
    doc.append("\n\\end{document}\n")

    tex = "".join(doc)
    open(os.path.join(HERE, "publication.tex"), "w", encoding="utf-8").write(tex)

    figdir = os.path.join(HERE, "figures")
    os.makedirs(figdir, exist_ok=True)
    used = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", "\n".join(md))
    for u in used:
        src = os.path.join(ROOT, "doc", u)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(figdir, os.path.basename(u)))

    print(f"publication.tex: {len(tex.splitlines())} lines, "
          f"{len(used)} figures, abstract from {abstract_src}")

    if not compile:
        return 0
    for _ in range(2):
        r = subprocess.run(["xelatex", "-interaction=nonstopmode",
                            "publication.tex"], cwd=HERE,
                           capture_output=True, text=True)
        if r.returncode != 0:
            log = open(os.path.join(HERE, "publication.log"), encoding="utf-8",
                       errors="replace").read()
            errs = [l for l in log.split("\n") if l.startswith("!")][:8]
            print("LATEX FAILED:\n" + "\n".join(errs))
            return 1
    pdf = os.path.join(HERE, "publication.pdf")
    pages = re.search(r"Output written on publication.pdf \((\d+)",
                      open(os.path.join(HERE, "publication.log"),
                          encoding="utf-8", errors="replace").read())
    print(f"publication.pdf: {os.path.getsize(pdf)/1e6:.2f} MB, "
          f"{pages.group(1) if pages else '?'} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
