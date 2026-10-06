#!/usr/bin/env python3
"""Sixth pass: INSPIRE HEP (the authoritative HEP bibliography) for the
entries Crossref could not settle."""
import json
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "CL10-bibliography-audit/1.0 (mailto:haouach.mohammed@gmail.com)"}


def insp(q, size=5):
    url = "https://inspirehep.net/api/literature?" + urllib.parse.urlencode(
        {"q": q, "size": size, "fields": "titles,authors,arxiv_eprints,dois,publication_info,citation_count"})
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.load(r).get("hits", {}).get("hits", [])
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                print(f"   ERR {e}")
                return []
            time.sleep(3 * (attempt + 1))
    return []


def show(hits):
    for h in hits:
        md = h.get("metadata", {})
        title = (md.get("titles") or [{}])[0].get("title", "?")
        pi = (md.get("publication_info") or [{}])[0]
        journal = pi.get("journal_title", "")
        vol = pi.get("journal_volume", "")
        pg = pi.get("artid") or pi.get("page_start") or ""
        yr = pi.get("year") or (md.get("imprints") or [{}])[0].get("date", "")
        dois = ", ".join(md.get("dois", [])[:2])
        arx = ", ".join(e.get("value", "") for e in md.get("arxiv_eprints", [])[:1])
        cc = md.get("citation_count", 0)
        au = ", ".join(a.get("full_name", "") for a in md.get("authors", [])[:5])
        print(f"   [{cc:5d} cit] {title[:76]}")
        print(f"       {journal} {vol} ({yr}) art {pg} | DOI {dois} | arXiv {arx}")
        print(f"       {au[:120]}")


QUERIES = [
    ("4  Arrighi Nesme Scherrer", 'a:"Arrighi, Pablo" AND a:"Nesme, Timothee" AND a:"Scherrer, Renaud" AND title:"quantum cellular automata"'),
    ("5  D'Ariano Perinotti Bahreyni Ford", 'title:"equivalence of quantum cellular automata"'),
    ("5b Bahreyni Ford", 'a:"Bahreyni, Davood" AND a:"Ford, G.W."'),
    ("13 PDG 2024", 'title:"Review of Particle Physics" AND a:"Navas, S."'),
    ("15 SK 0.37 Mt yr", 'title:"Search for proton decay" AND a:"Abe, K." AND reportnum:"Super-Kamiokande"'),
    ("19a N&N no-go", 'a:"Nielsen, H.B." AND a:"Ninomiya, M." AND 1981'),
    ("19b proof NN", 'title:"proof of the Nielsen-Ninomiya theorem"'),
    ("23 John blow-up", 'a:"John, F." AND title:"blow-up"'),
    ("24 Sakharov 1967", 'title:"Violation of CP invariance"'),
    ("29 tHooft deterministic", 'a:"t Hooft, G." AND title:"deterministic"'),
    ("31 quasi-degenerate leptogenesis", 'title:"Leptogenesis with quasi-degenerate neutrinos"'),
    ("31b Alloul Ariga", 'a:"Alloul, B." OR a:"Ariga, A."'),
    ("32 Barbieri Hall Strumia", 'a:"Barbieri, Riccardo" AND a:"Hall, Lawrence J." AND a:"Strumia, Alessandro"'),
    ("12c Distler Fraser", 'a:"Distler, Johanna" AND a:"Fraser, H."'),
    ("12c2 heterotic fermions NCG", 'title:"Heterotic fermions in noncommutative geometry"'),
]
for label, q in QUERIES:
    print(f"=== {label} ===")
    show(insp(q))
    time.sleep(1.0)
