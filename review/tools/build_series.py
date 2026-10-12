#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the three companion papers from the corrected master manuscript.

Mechanical assembly: header prose + verbatim body ranges + footer prose,
with per-paper remapping of section cross-references, figure numbers,
table numbers and reference numbers. The master doc/final_publication.md
is read-only input.
"""
import re
MASTER = "doc/final_publication.md"
HDR = "review/tools/headers"

lines = open(MASTER, encoding="utf-8").read().split("\n")


def rng(a, b):
    """1-based inclusive line range."""
    return "\n".join(lines[a - 1:b])


# ------------------------------------------------- section boundaries (dynamic)
def sect_line(title):
    return next(i + 1 for i, l in enumerate(lines) if l.strip() == title)


def body_end(heading_line):
    """last line before the '---' separator preceding this heading"""
    sep = next(i for i in range(heading_line - 1, 0, -1)
               if lines[i - 1].strip() == "---")
    return sep - 1


H10 = sect_line("## 10. The gauged bond layer")
H12 = sect_line("## 12. The loop factor and the scale ladder")
H18 = sect_line("## 18. Discussion")
P1 = (24, body_end(H10))
P2 = (H10, body_end(H12))
P3 = (H12, body_end(H18))
HA = sect_line("## Appendix A. The explicit 32-dimensional representation")
REF = sect_line("## References")
HA_END = next(i for i in range(HA + 1, REF)
              if lines[i - 1].strip() == "---") - 1


# ---------------------------------------------------------------- reference list
ref_start = next(i for i, l in enumerate(lines) if l.strip() == "## References")
refs = {}
for l in lines[ref_start:]:
    m = re.match(r"^(\d+)\.\s", l)
    if m:
        refs[int(m.group(1))] = l

avail = next(l for l in lines if l.startswith("**Code and data availability.**"))
coi = next(l for l in lines if l.startswith("**Competing interests.**"))

# ---------------------------------------------------------------- helpers
CITE = re.compile(r"\[(\d+(?:\s*[–,]\s*\d+)*)\]")


def expand(token):
    out = []
    for part in re.split(r"\s*[,]\s*", token):
        if "–" in part:
            a, b = (int(x) for x in part.split("–"))
            out.extend(range(a, b + 1))
        else:
            out.append(int(part))
    return out


def render(nums):
    nums = sorted(set(nums))
    if len(nums) == 1:
        return str(nums[0])
    if nums == list(range(nums[0], nums[-1] + 1)):
        return f"{nums[0]}–{nums[-1]}"
    return ",".join(str(n) for n in nums)


def remap_refs(text, order):
    def sub(m):
        old = expand(m.group(1))
        for n in old:
            if n not in order:
                order[n] = len(order) + 1
        return "[" + ",".join(str(order[n]) for n in
                              sorted(set(old), key=lambda n: order[n])) + "]"
    return CITE.sub(sub, text)


SECT = re.compile(
    r"\b(Sections?)\s+(\d+(?:\.\d+)?"
    r"(?:\s*(?:[–,]|and|&)\s*\d+(?:\.\d+)?)*)")
NUM = re.compile(r"\d+(?:\.\d+)?")


def remap_sections(text, smap):
    """smap: master section int -> (paper_label or None, new number or label)."""
    def token(p):
        key = int(p.split(".")[0])
        lab, new = smap[key]
        if isinstance(new, str) and not NUM.fullmatch(new):
            return f"[{lab} {new}]" if lab else new
        return f"{new}{p[len(str(key)):]}"  # keep sub-number like ".4"

    def sub(m):
        word, phrase = m.group(1), m.group(2)
        parts = re.split(r"(\s*(?:[–,]|and|&)\s*)", phrase)
        papers = [smap[int(p.split('.')[0])][0] for p in parts if NUM.fullmatch(p)]
        if all(x is None for x in papers):
            return (f"{word} " + "".join(
                token(p) if NUM.fullmatch(p) else p for p in parts))
        uniq = list(dict.fromkeys(papers))
        if len(uniq) == 1:
            lab = uniq[0]
            joined = "".join(token(p) if NUM.fullmatch(p) else p for p in parts)
            if all(isinstance(smap[int(p.split('.')[0])][1], str)
                   and not NUM.fullmatch(smap[int(p.split('.')[0])][1])
                   for p in parts if NUM.fullmatch(p)):
                return f"[{lab} {joined}]"
            return f"[{lab} {word} {joined}]"
        # mixed internal/external: label each external token individually
        out2 = [token(p) if NUM.fullmatch(p) else p for p in parts]
        head = f"{word} " if papers and papers[0] is None else ""
        return head + "".join(out2)
    return SECT.sub(sub, text)



FIG = re.compile(r"\bFigs?\.\s*(\d+(?:[a-z])?(?:\s*[–,]\s*\d+(?:[a-z])?)*)")

def remap_figs(text, fmap):
    def sub(m):
        phrase = m.group(1)
        parts = re.split(r"(\s*[–,]\s*)", phrase)
        out = []
        for p in parts:
            mm = re.fullmatch(r"(\d+)([a-z]?)", p)
            if mm:
                out.append(str(fmap[int(mm.group(1))]) + mm.group(2))
            else:
                out.append(p)
        word = "Figs." if ("," in phrase or "–" in phrase) else "Fig."
        return word + " " + "".join(out)
    return FIG.sub(sub, text)


TBL = re.compile(r"\bTable\s+(\d+)")


def remap_tables(text, tmap):
    def sub(m):
        return str(tmap[int(m.group(1))])
    return TBL.sub(sub, text)


HEAD = re.compile(r"(?m)^(##+ )(\d+)((?:\.\d+)?)\.")


def renumber_headings(body, smap):
    def sub(m):
        key = int(m.group(2))
        lab, new = smap[key]
        assert lab is None, f"external heading {key}"
        return f"{m.group(1)}{new}{m.group(3)}."
    return HEAD.sub(sub, body)

def build(out_path, head, body_ranges, foot, smap, fmap, tmap, with_appendix):
    order = {}
    body = [rng(a, b) for a, b in body_ranges]
    if with_appendix:
        body.append(rng(HA, HA_END))
    body = "\n\n".join(body)
    body = renumber_headings(body, smap)
    body = remap_sections(body, smap)
    body = remap_figs(body, fmap)
    body = remap_tables(body, tmap)
    text = "\n\n".join(
        [open(head, encoding="utf-8").read().rstrip(), body,
         open(foot, encoding="utf-8").read().rstrip()])
    text = remap_refs(text, order)
    ref_lines = ["## References", ""]
    for old, new in sorted(order.items(), key=lambda kv: kv[1]):
        ref_lines.append(re.sub(r"^\d+\.", f"{new}.", refs[old]))
    text += "\n\n" + "\n".join(ref_lines) + "\n\n" + avail + "\n\n" + coi + "\n"
    if not with_appendix:
        text = text.replace("Appendix A", "[Paper I] Appendix A")
    open(out_path, "w", encoding="utf-8").write(text)
    print(f"{out_path}: {len(order)} references, "
          f"{text.count(chr(10)) + 1} lines")


I, II, III = "Paper I", "Paper II", "Paper III"

# Paper I: body = Sections 1-9 (24-341) + Appendix A
s1 = {n: (None, n) for n in range(1, 10)}
s1.update({10: (II, 2), 11: (II, 3), 12: (III, 2), 13: (III, 3),
           14: (III, 4), 15: (III, 5), 16: (III, 6), 17: (III, 7),
           18: (II, "Discussion"), 19: (III, "Conclusion")})
f1 = {n: n for n in range(1, 26)}
t1 = {1: "Table 1", 2: "Table 2", 3: f"[{II}, Table 1]",
      4: f"[{III}, Table 1]", 5: f"[{III}, Table 2]",
      6: f"[{III}, Table 3]", 7: f"[{III}, Table 4]",
      8: f"[{III}, Table 5]"}
build("doc/paper1_core.md", f"{HDR}/paper1_head.md", [P1],
      f"{HDR}/paper1_foot.md", s1, f1, t1, True)

# Paper II: body = Sections 10-11 (343-482)
s2 = {n: (I, n) for n in range(1, 10)}
s2.update({10: (None, 2), 11: (None, 3)})
s2.update({12: (III, 2), 13: (III, 3), 14: (III, 4), 15: (III, 5),
           16: (III, 6), 17: (III, 7), 18: (III, "Discussion"),
           19: (III, "Conclusion")})
f2 = {n: n for n in range(1, 10)}
f2.update({10: 1, 11: 2, 12: 3, 13: 1, 14: 2, 15: 3, 16: 4, 17: 5,
           18: 6, 19: 7, 20: 8, 21: 9, 22: 10, 23: 11, 24: 12, 25: 13})
t2 = {1: f"[{I}, Table 1]", 2: f"[{I}, Table 2]", 3: "Table 1",
      4: f"[{III}, Table 1]", 5: f"[{III}, Table 2]",
      6: f"[{III}, Table 3]", 7: f"[{III}, Table 4]",
      8: f"[{III}, Table 5]"}
build("doc/paper2_gauge_gravity.md", f"{HDR}/paper2_head.md", [P2],
      f"{HDR}/paper2_foot.md", s2, f2, t2, False)

# Paper III: body = Sections 12-17 (484-849)
s3 = {n: (I, n) for n in range(1, 10)}
s3.update({10: (II, 2), 11: (II, 3)})
s3.update({12: (None, 2), 13: (None, 3), 14: (None, 4), 15: (None, 5),
           16: (None, 6), 17: (None, 7), 18: (III, "Discussion"),
           19: (III, "Conclusion")})
f3 = {n: n for n in range(1, 10)}
f3.update({13: 1, 14: 2, 15: 3, 16: 4, 17: 5, 18: 6, 19: 7, 20: 8,
           21: 9, 22: 10, 23: 11, 24: 12, 25: 13})
t3 = {1: f"[{I}, Table 1]", 2: f"[{I}, Table 2]",
      3: f"[{II}, Table 1]", 4: "Table 1", 5: "Table 2", 6: "Table 3",
      7: "Table 4", 8: "Table 5"}
build("doc/paper3_phenomenology.md", f"{HDR}/paper3_head.md", [P3],
      f"{HDR}/paper3_foot.md", s3, f3, t3, False)
