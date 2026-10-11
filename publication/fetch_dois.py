#!/usr/bin/env python3
"""Resolve and verify DOIs for the 53 references in doc/final_publication.md.

For each numbered reference: search Crossref by the quoted title (math
stripped), verify title similarity + year + (when present) volume/page.
Multi-paper entries ("A ... (1987); B ... (2003)") are resolved on the
FIRST paper. Only verified DOIs are written.

Output: publication/doi.json  {"<ref number>": "<doi>", ...}
"""
import difflib
import json
import re
import time
import urllib.parse
import urllib.request

SRC = "../doc/final_publication.md"
OUT = "doi.json"

text = open(SRC, encoding="utf-8").read()
m = re.search(r"^## References\s*$", text, re.M)
assert m, "References section not found"
refs_block = text[m.end():]
entries = re.findall(r"^(\d{1,2})\. (.*?)(?=^\d{1,2}\. |\Z)",
                     refs_block, re.S | re.M)
print("references parsed:", len(entries))


def get(url):
    req = urllib.request.Request(
        url, headers={"User-Agent": "CL10-refs/1.0 (mailto:haouach.mohammed@gmail.com)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


MATH_WORDS = {
    r"p\s*\\to\s*e\^\+\\pi\^0": "p to e pi0",
    r"p\s*\\to\s*\\mu\^\+\\pi\^0": "p to mu pi0",
    r"p\s*\\to\s*\\nu\s*K\^+": "p to nu K",
    r"\\nu\s*K\^+": "nu K",
}


def norm(s):
    for k, v in MATH_WORDS.items():
        s = re.sub(k, v, s)
    s = re.sub(r"\$[^$]*\$", " ", s)
    s = re.sub(r"[*_\\{}^+⋅]", " ", s)
    return re.sub(r"\s+", " ", s.lower()).strip()


def first_paper(entry):
    """'A ... (1987); B ... (2003)' -> first paper only."""
    parts = re.split(r";\s*(?=[A-Z][^;]*?\(\d{4}\))", entry)
    return parts[0]


def title_of(entry):
    m = re.search(r'"([^"]+)"', entry)
    t = m.group(1) if m else None
    return t.rstrip(",;") if t else None


def year_of(entry):
    ys = re.findall(r"\((19[5-9]\d|20[0-2]\d)\)", entry)
    return int(ys[0]) if ys else None


def vol_page(entry):
    m = re.search(r"\*\*(\d[\w()\-]*)\*\*,?\s*([A-Za-z]?\d[\w\-]*)", entry)
    return (m.group(1), m.group(2)) if m else (None, None)


def authors_of(entry):
    head = entry.split('"')[0]
    m = re.match(r"\s*(?:[A-Z]\.\s*)*([A-Z][\w'\-]+)", head)
    return m.group(1) if m else None


results = {}
for num, entry in entries:
    entry1 = first_paper(entry)
    title = title_of(entry1)
    if not title:
        print(f"[{num:>2}] no quoted title -> skip")
        continue
    year = year_of(entry1)
    vol, page = vol_page(entry1)
    author = authors_of(entry1)
    q = urllib.parse.quote(norm(title))
    aq = f"&query.author={urllib.parse.quote(author)}" if author else ""
    try:
        data = get(f"https://api.crossref.org/works?query.bibliographic={q}{aq}"
                   f"&rows=5&select=DOI,title,issued,volume,page,container-title")
    except Exception as e:
        print(f"[{num:>2}] crossref error: {e}")
        time.sleep(1)
        continue
    best, best_score = None, 0.0
    for it in data.get("message", {}).get("items", []):
        t = (it.get("title") or [""])[0]
        score = difflib.SequenceMatcher(None, norm(title), norm(t)).ratio()
        iy = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        if score > 0.80:
            if year and iy and abs(iy - year) > 1:
                score -= 0.10
            if vol and it.get("volume") and str(it["volume"]) == str(vol):
                score += 0.05
            if page and it.get("page") and str(it["page"]).startswith(str(page)):
                score += 0.05
        if score > best_score:
            best, best_score = it, score
    if best and best_score >= 0.85:
        doi = best["DOI"]
        iy = (best.get("issued", {}).get("date-parts") or [[None]])[0][0]
        ok_year = (not year) or (iy and abs(iy - year) <= 1)
        if ok_year:
            results[num] = doi
            print(f"[{num:>2}] OK  {doi}  ({iy}) {best['title'][0][:55]}")
        else:
            print(f"[{num:>2}] YEAR MISMATCH {iy} vs {year}: {doi}")
    else:
        bt = (best.get("title") or [""])[0] if best else ""
        print(f"[{num:>2}] NO MATCH (best {best_score:.2f}) {title[:45]} | got: {bt[:45]}")
    time.sleep(0.3)

# DOIs verified by direct resolution (doi.org / Crossref works API / INSPIRE)
# on 2026-10-10; the fuzzy title search above cannot reach these entries
# (generic titles, translations, or multi-paper entries).
MANUAL_DOIS = {
    5: "10.1007/s11467-016-0616-z",
    13: "10.1103/PhysRevD.110.030001",
    14: "10.1051/0004-6361/201833910",
    15: "10.1103/PhysRevD.102.112011",
    16: "10.1126/science.adq9592",
    29: "10.1016/0550-3213(85)90418-3",
    49: "10.1016/0550-3213(83)90610-7",
    52: "10.1103/PhysRevD.90.072005",
    53: "10.1103/4dyc-z8zf",
}
for _n, _d in MANUAL_DOIS.items():
    results.setdefault(_n, _d)

json.dump(results, open(OUT, "w"), indent=2)
print(f"\nwrote {OUT}: {len(results)} verified DOIs of {len(entries)}")
