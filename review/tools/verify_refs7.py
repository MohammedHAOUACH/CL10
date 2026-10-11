#!/usr/bin/env python3
"""Verify the CURRENT bibliography of doc/final_publication.md against Crossref.

Successor to verify_refs.py: the claims are PARSED from the manuscript instead
of hardcoded, so renumbering cannot stale the audit. For each numbered entry
" N. Authors, "Title," *Journal* **vol**, page (year)." the FIRST paper of a
multi-paper entry is claimed (as in verify_refs.py); entries with no volume
(books) are not Crossref-checked. Also checks numbering contiguity, duplicate
entries, and that every entry is cited in the body and every body citation
resolves. Writes review/refs_audit.json.
"""
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify_refs import check  # noqa: E402

SRC = os.path.join(HERE, os.pardir, os.pardir, "doc", "final_publication.md")
OUT = os.path.join(HERE, os.pardir, "refs_audit.json")

ENTRY_RE = re.compile(r"^(\d{1,2})\.\s+(.*)$")
VOL_RE = re.compile(r"\*\*(?P<vol>[^*]+?)\*\*")
SEG_RE = re.compile(
    r'^(?P<authors>[^"]+?),\s*"(?P<title>[^"]+?),?"\s*,?\s*'
    r"\*(?P<journal>[^*]+?)\*\*?\s*(?P<rest>.*)$")
PAGE_RE = re.compile(r"\*\*[^*]+?\*\*,?\s*(?P<page>[A-Za-z]?\d[\w.\-]*)")
YEAR_RE = re.compile(r"\((\d{4})\)")


def parse(text):
    """(number, authors, title, journal, volume, page, year) per entry."""
    entries = []
    in_refs = False
    for line in text.splitlines():
        if line.strip() == "## References":
            in_refs = True
            continue
        if in_refs and line.startswith("## "):
            break
        if not in_refs:
            continue
        m = ENTRY_RE.match(line)
        if not m:
            continue
        n, body = int(m.group(1)), m.group(2)
        seg = body.split(";")[0].strip()          # first paper of the entry
        sm = SEG_RE.match(seg)
        if not sm:
            if '"' not in seg:                   # book: italic title, no quote
                ym = YEAR_RE.search(seg)
                entries.append((n, seg.split(",")[0], seg, "book", None, None,
                                int(ym.group(1)) if ym else None))
            else:
                entries.append((n, seg, None, None, None, None, None))
            continue
        rest = sm.group("rest")
        vm = VOL_RE.search(rest)
        pm = PAGE_RE.search(rest)
        ym = YEAR_RE.search(rest)
        entries.append((
            int(n), sm.group("authors").strip(), sm.group("title").strip(),
            sm.group("journal").strip(),
            vm.group("vol").strip() if vm else None,
            pm.group("page").strip() if pm else None,
            int(ym.group(1)) if ym else None))
    return entries


def body_citations(text):
    """Set of cited reference numbers in the body (before ## References)."""
    body = text.split("## References")[0]
    cited = set()
    for m in re.finditer(r"\[(\d{1,2}(?:\s*(?:,|\u2013|-){1,2}\s*\d{1,2})*)\]",
                         body):
        for part in m.group(1).split(","):
            part = part.strip().replace("\u2013", "-").replace("--", "-")
            if "-" in part:
                a, b = part.split("-")
                cited.update(range(int(a), int(b) + 1))
            else:
                cited.add(int(part))
    return cited


def main():
    text = open(SRC, encoding="utf-8").read()
    entries = parse(text)
    numbers = [e[0] for e in entries]
    contiguous = numbers == list(range(1, len(numbers) + 1))

    # duplicate titles
    by_title = {}
    for e in entries:
        if e[2]:
            by_title.setdefault(e[2].lower(), []).append(e[0])
    dupes = {t: ns for t, ns in by_title.items() if len(ns) > 1}

    cited = body_citations(text)
    uncited = [n for n in numbers if n not in cited]
    dangling = sorted(c for c in cited if c not in numbers)
    results = []
    for n, authors, title, journal, vol, page, year in entries:
        if title is None:
            results.append({"n": n, "verdict": "PARSE_FAIL", "claim": authors})
            print(f"[{n:2d}] !! PARSE_FAIL {authors[:70]}")
            continue
        claim = f"{authors}, {title}, {journal} {vol}, {page} ({year})"
        if vol is None:
            results.append({"n": n, "verdict": "BOOK (not Crossref-checked)",
                            "claim": claim})
            print(f"[{n:2d}] BOOK  {title[:70]}")
            continue
        surname = re.search(r"\(([A-Z][\w'\-]+)", authors)
        lead = (surname.group(1) if surname
                else authors.split(",")[0].split()[-1])
        pg = re.search(r"\d+", page) if page else None
        r = check(n, title, lead, journal, vol,
                  int(pg.group()) if pg else None, year)
        r["claim"] = claim
        results.append(r)
        flag = "OK " if r["verdict"] == "OK" else "!! "
        print(f"[{n:2d}] {flag}{r['verdict']:9s} {r.get('crossref_title', '')[:70]}")
        for p in r.get("problems", []):
            print(f"      - {p}")
        time.sleep(0.4)

    audit = {
        "source": "doc/final_publication.md (parsed, current numbering)",
        "n_entries": len(entries),
        "contiguous": contiguous,
        "duplicate_titles": dupes,
        "uncited_entries": uncited,
        "dangling_citations": dangling,
        "entries": results,
    }
    with open(OUT, "w") as f:
        json.dump(audit, f, indent=2, ensure_ascii=False)

    ok = sum(1 for r in results if r["verdict"] == "OK")
    books = sum(1 for r in results if r["verdict"].startswith("BOOK"))
    bad = [r for r in results
           if r["verdict"] not in ("OK",) and not r["verdict"].startswith("BOOK")]
    print(f"\nentries: {len(entries)}  contiguous: {contiguous}  "
          f"OK: {ok}  books: {books}  flagged: {len(bad)}")
    if dupes:
        print("duplicate titles: " + "; ".join(
            f"{t[:40]} -> {ns}" for t, ns in dupes.items()))
    if uncited:
        print(f"uncited entries: {uncited}")
    if dangling:
        print(f"dangling citations: {dangling}")
    for r in bad:
        print(f"  ref {r['n']}: {r['verdict']} {r.get('problems') or r.get('error', '')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
