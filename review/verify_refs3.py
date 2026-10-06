#!/usr/bin/env python3
"""Third pass: direct DOI resolution for candidate DOIs + author/container
queries for the entries Crossref fuzzy search missed."""
import json
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "CL10-bibliography-audit/1.0 (mailto:haouach.mohammed@gmail.com)"}


def get(url):
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                return {"error": str(e)}
            time.sleep(2 * (attempt + 1))
    return {"error": "unreachable"}


def fmt(item):
    if "error" in item:
        return f"NOT FOUND ({item['error'][:60]})"
    m = item.get("message", item)
    return (f"{(m.get('title') or ['?'])[0][:75]} | "
            f"{(m.get('container-title') or ['?'])[0][:35]} | "
            f"vol {m.get('volume')} | pg {m.get('page')} | "
            f"{(m.get('issued',{}).get('date-parts') or [[None]])[0][0]} | "
            f"{', '.join(a.get('family','') for a in m.get('author',[])[:5])}")


print("=== C. direct DOI resolution (candidate DOIs) ===")
CAND = {
    "13 PDG2024": "10.1103/PhysRevD.110.070001",
    "14 Planck2018VI": "10.1051/0004-6361/201833910",
    "15 SK 0.37Mtyr": "10.1103/PhysRevD.102.012006",
    "19a NN no-go PLB105": "10.1016/0370-2693(81)91221-9",
    "19b NN proof CMP78": "10.1007/BF01208280",
    "29 tHooft PLB364": "10.1016/0370-2693(95)01296-5",
    "31 Blanchet PRL110": "10.1103/PhysRevLett.110.191601",
    "12c DistlerFraser": "10.1007/BF02101611",
    "24 Sakharov (doi guess)": "10.1070/JETP.1967v005p0024",
}
for label, doi in CAND.items():
    print(f"[{label}] {fmt(get('https://api.crossref.org/works/' + urllib.parse.quote(doi)))}")
    time.sleep(0.3)

print("\n=== D. author+container targeted queries ===")
QUERIES = [
    ("4 Arrighi NatComput", "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.author": "Arrighi Nesme Scherrer", "query.container-title": "Natural Computing",
         "rows": 6, "select": "title,author,container-title,volume,page,DOI,issued"})),
    ("5 D'Ariano JCP", "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.author": "D'Ariano Perinotti Bahreyni Ford", "rows": 6,
         "query": "quantum cellular automata circuits",
         "select": "title,author,container-title,volume,page,DOI,issued"})),
    ("7 Fritzsch Minkowski", "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.author": "Fritzsch Minkowski", "rows": 6,
         "query.bibliographic": "Unified interactions of strong and electroweak forces of quarks and leptons",
         "select": "title,author,container-title,volume,page,DOI,issued"})),
    ("23 John CPAM", "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.author": "Fritz John", "rows": 8,
         "query.bibliographic": "Blow-up of solutions of nonlinear wave equations in three space dimensions",
         "select": "title,author,container-title,volume,page,DOI,issued"})),
    ("32 Barbieri Hall Strumia", "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.author": "Barbieri Hall Strumia", "rows": 6,
         "query.bibliographic": "Minimal thermal leptogenesis neutrino mass spectrum",
         "select": "title,author,container-title,volume,page,DOI,issued"})),
    ("12b Castro AAC 2015", "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.author": "Castro", "rows": 6,
         "query.bibliographic": "Note on a Clifford algebra based grand unification program of gravity and the Standard Model",
         "select": "title,author,container-title,volume,page,DOI,issued"})),
]
for label, url in QUERIES:
    d = get(url)
    items = d.get("message", {}).get("items", []) if "message" in d else []
    print(f"[{label}]")
    for it in items:
        print(f"   {fmt(it)}  DOI={it.get('DOI')}")
    if not items:
        print("   (no hits)")
    time.sleep(0.3)
