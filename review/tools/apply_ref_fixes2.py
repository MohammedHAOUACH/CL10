#!/usr/bin/env python3
r"""Remove the duplicate bibliography entry (refs 2 and 30 are both
K. G. Wilson, "Confinement of quarks", found by verify_refs7.py).

Entry 30 is deleted; entries 31-56 shift to 30-55. Body citations remap:
[30] -> [2] (same paper), and every cited number n in 31..56 -> n-1.
Ranges are mapped endpoint-wise; a range straddling 30 is split.
"""
import os
import re
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   os.pardir, os.pardir, "doc", "final_publication.md")
text = open(SRC, encoding="utf-8").read()
body, refs = text.split("## References", 1)

CITE_RE = re.compile(r"\[(\d{1,2}(?:\s*(?:,|\u2013|-){1,2}\s*\d{1,2})*)\]")


def mapnum(n):
    if n == 30:
        return 2
    return n - 1 if 31 <= n <= 56 else n


def cite_repl(m):
    parts = []
    for part in m.group(1).split(","):
        part = part.strip().replace("\u2013", "-").replace("--", "-")
        if "-" in part:
            a, b = (int(x) for x in part.split("-"))
            ma, mb = mapnum(a), mapnum(b)
            if a <= 30 <= b and not (a == b == 30):
                # straddles the deleted entry: emit explicit list
                nums = sorted({mapnum(n) for n in range(a, b + 1)})
                parts.append(",".join(str(n) for n in nums))
            else:
                parts.append("%d\u2013%d" % (ma, mb) if ma < mb else str(ma))
        else:
            parts.append(str(mapnum(int(part))))
    return "[%s]" % ",".join(parts)


body_new = CITE_RE.sub(cite_repl, body)

# references: drop entry 30, renumber 31..56 -> 30..55
lines = refs.splitlines()
out = []
for line in lines:
    m = re.match(r"^(\d{1,2})\.\s", line)
    if m:
        n = int(m.group(1))
        if n == 30:
            continue
        if n > 30:
            line = re.sub(r"^\d{1,2}\.", "%d." % (n - 1), line, count=1)
    out.append(line)
refs_new = "\n".join(out)

new = body_new + "## References" + refs_new
open(SRC, "w", encoding="utf-8").write(new)

# report
cited = set()
for m in CITE_RE.finditer(body_new):
    for part in m.group(1).split(","):
        part = part.strip().replace("\u2013", "-").replace("--", "-")
        if "-" in part:
            a, b = (int(x) for x in part.split("-"))
            cited.update(range(a, b + 1))
        else:
            cited.add(int(part))
nums = [int(m.group(1)) for m in
        re.finditer(r"^(\d{1,2})\.\s", refs_new, re.M)]
print("entries:", len(nums), "contiguous:", nums == list(range(1, len(nums) + 1)))
print("uncited:", [n for n in nums if n not in cited])
print("dangling:", sorted(c for c in cited if c not in nums))
