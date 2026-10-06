#!/usr/bin/env python3
"""Fourth pass: arXiv API verification for entries Crossref indexes poorly."""
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

NS = {"a": "http://www.w3.org/2005/Atom"}


def arxiv(q, max_results=6):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({
        "search_query": q, "max_results": max_results,
        "sortBy": "relevance"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return ET.fromstring(r.read())
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                print(f"   ERR {e}")
                return None
            time.sleep(3 * (attempt + 1))
    return None


def show(root):
    if root is None:
        return
    for e in root.findall("a:entry", NS):
        t = re.sub(r"\s+", " ", (e.findtext("a:title", "", NS) or "")).strip()
        aid = (e.findtext("a:id", "", NS) or "").split("/abs/")[-1]
        pub = (e.findtext("a:published", "", NS) or "")[:10]
        jm = e.findtext("a:journal_ref", "", NS) or e.findtext("a:arxiv_comment", "", NS) or ""
        au = ", ".join((x.findtext("a:name", "", NS) or "") for x in e.findall("a:author", NS))
        print(f"   {aid} [{pub}] {t[:78]}")
        print(f"      authors: {au[:110]}")
        if jm:
            print(f"      journal_ref: {jm[:110]}")


QUERIES = [
    ("4  Arrighi Nesme Scherrer QCA intro", 'ti:"introduction to quantum cellular automata" AND au:Arrighi'),
    ("5  D'Ariano Perinotti Bahreyni Ford", 'all:"equivalence of quantum cellular automata and circuits"'),
    ("5b Bahreyni Ford QCA", 'au: Bahreyni AND au: Ford AND all:"quantum cellular automata"'),
    ("13 PDG 2024 review", 'ti:"Review of Particle Physics" AND cat:hep-ex'),
    ("15 SK proton decay 0.37 Mt yr", 'all:"0.37 megaton-years exposure of Super-Kamiokande"'),
    ("19b N&N proof CMP", 'ti:"proof of the Nielsen-Ninomiya theorem"'),
    ("19c N&N no-go 1981", 'au:Ninomiya AND ti:"no-go theorem"'),
    ("29 t Hooft deterministic QM field", 'ti:"Deterministic quantum mechanics and field theories"'),
    ("31 quasi-degenerate leptogenesis PRL", 'all:"leptogenesis with quasi-degenerate neutrinos"'),
    ("31b Blanchet Di Bari", 'au:Blanchet AND au:"Di Bari" AND all:leptogenesis'),
    ("32 Barbieri Hall Strumia leptogenesis", 'au:Barbieri AND au:Strumia AND all:"minimal thermal leptogenesis"'),
    ("12c Distler Fraser heterotic", 'ti:"Heterotic fermions in noncommutative geometry"'),
    ("24 Sakharov", 'ti:"Violation of CP invariance" AND au:Sakharov'),
]
for label, q in QUERIES:
    print(f"=== {label} ===")
    show(arxiv(q))
    time.sleep(1.2)
