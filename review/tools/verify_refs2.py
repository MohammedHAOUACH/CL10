#!/usr/bin/env python3
"""Second pass: verify flagged entries by DIRECT DOI resolution and exact-title
query (query.title + container filter), instead of fuzzy bibliographic search."""
import json
import os
import re
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
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
        return f"ERR {item['error']}"
    m = item.get("message", item)
    return (f"{(m.get('title') or ['?'])[0][:80]} | "
            f"{(m.get('container-title') or ['?'])[0][:38]} | "
            f"vol {m.get('volume')} | pg {m.get('page')} | "
            f"{(m.get('issued',{}).get('date-parts') or [[None]])[0][0]} | "
            f"DOI {m.get('DOI')} | "
            f"{', '.join(a.get('family','') for a in m.get('author',[])[:5])}")


# --- A. claimed DOIs (formerly livrable/doi.json, now removed): resolve each directly
CLAIMED_DOIS = {
    1: "10.1103/PhysRevLett.59.521", 2: "10.1103/PhysRevD.14.3432",
    3: "10.1007/BF02199356", 6: "10.1103/PhysRevLett.32.438",
    8: "10.1016/0370-1573(81)90092-2", 10: "10.1016/0040-9383(64)90003-5",
    12: "10.1139/cjp-2013-0686", 16: "10.1126/science.adq9592",
    17: "10.1007/JHEP09(2020)178", 20: "10.1103/PhysRev.138.B988",
    21: "10.1016/0370-2693(80)90212-9", 22: "10.1103/PhysRevLett.70.9",
    25: "10.1103/PhysRevD.42.3344", 26: "10.1016/0550-3213(88)90133-2",
    27: "10.1016/S0370-2693(02)01336-9", 28: "10.1103/PhysRevD.10.2445",
}
print("=== A. claimed DOIs resolved directly ===")
for n, doi in CLAIMED_DOIS.items():
    print(f"[{n:2d}] {fmt(get('https://api.crossref.org/works/' + urllib.parse.quote(doi)))}")
    time.sleep(0.25)

# --- B. exact-title queries for the flagged ones
EXACT = [
    (4,  "An introduction to quantum cellular automata", "Arrighi"),
    (5,  "Equivalence of quantum cellular automata and circuits", "D'Ariano"),
    (7,  "Unified interactions of strong and electroweak forces of quarks and leptons", "Fritzsch"),
    (13, "Review of Particle Physics", "Navas"),
    (14, "Planck 2018 results. VI. Cosmological parameters", "Aghanim"),
    (15, "Search for proton decay via p→e+π0 and p→μ+π0 in 0.37 megaton-years exposure of Super-Kamiokande", "Abe"),
    (16, "Direct neutrino-mass measurement based on 259 days of KATRIN data", "Aker"),
    (19, "A no-go theorem to avoid the fermion doubling in lattice gauge theory", "Nielsen"),
    (19, "A proof of the Nielsen-Ninomiya theorem", "Nielsen"),
    (23, "Blow-up of solutions of nonlinear wave equations in three space dimensions", "John"),
    (24, "Violation of CP Invariance, C Asymmetry, and Baryon Asymmetry of the Universe", "Sakharov"),
    (29, "Deterministic quantum mechanics and field theories", "t Hooft"),
    (31, "Leptogenesis with quasi-degenerate neutrinos", "Blanchet"),
    (32, "Minimal thermal leptogenesis: The neutrino mass spectrum", "Barbieri"),
]
print("\n=== B. exact-title queries (query.title) ===")
for n, title, author in EXACT:
    q = urllib.parse.urlencode({"query.title": title, "rows": 5,
                                "select": "title,author,container-title,volume,page,DOI,issued"})
    d = get("https://api.crossref.org/works?" + q)
    items = d.get("message", {}).get("items", []) if "message" in d else []
    print(f"[{n:2d}] {title[:55]}")
    for it in items:
        print(f"      {fmt(it)}")
    if not items:
        print("      (no hits)")
    time.sleep(0.3)
