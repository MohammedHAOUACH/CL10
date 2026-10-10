#!/usr/bin/env python3
"""Verify every bibliography entry of doc/publication.md against Crossref.

For each claimed entry we query Crossref by title (+ author) and compare the
top hit's title / container / volume / first page / year against the claim.
Writes review/refs_audit.json with a per-entry verdict.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, "refs_audit.json")
UA = {"User-Agent": "CL10-bibliography-audit/1.0 (mailto:haouach.mohammed@gmail.com)"}

# The 53 entries exactly as printed in doc/publication.md (2026-10-10 rewrite).
# Multi-paper entries are claimed on their FIRST paper; the second paper is
# checked by verify_refs3/5 where needed. Books are not Crossref-checked.
CLAIMS = [
    (1,  "Space-time as a causal set", "Bombelli", "Phys. Rev. Lett.", 59, 521, 1987),
    (2,  "Computation of the quantum effects due to a four-dimensional pseudoparticle", "'t Hooft", "Phys. Rev. D", 14, 3432, 1976),
    (3,  "From quantum cellular automata to quantum lattice gases", "Meyer", "J. Stat. Phys.", 85, 551, 1996),
    (4,  "An overview of quantum cellular automata", "Arrighi", "Nat. Comput.", 18, 885, 2019),
    (5,  "Quantum cellular automata and free quantum field theory", "D'Ariano", "Front. Phys.", 12, 120301, 2016),
    (6,  "Unity of all elementary-particle forces", "Georgi", "Phys. Rev. Lett.", 32, 438, 1974),
    (7,  "Unified interactions of leptons and hadrons", "Fritzsch", "Ann. Phys.", 93, 193, 1975),
    (8,  "Group theory for unified model building", "Slansky", "Phys. Rep.", 79, 1, 1981),
    (9,  "Lie Algebras in Particle Physics", "Georgi", "book", None, None, 1999),
    (10, "Clifford modules", "Atiyah", "Topology", 3, 3, 1964),
    (11, "Spinors and Space-Time", "Penrose", "book", None, None, 1984),
    (12, "A Clifford algebra-based grand unification program of gravity and the Standard Model: a review study", "Castro", "Can. J. Phys.", 92, 1501, 2014),
    (13, "Review of Particle Physics", "Navas", "Phys. Rev. D", 110, 30001, 2024),
    (14, "Planck 2018 results. VI. Cosmological parameters", "Aghanim", "Astron. Astrophys.", 641, 6, 2020),
    (15, "Search for proton decay via p to e+ pi0 and p to mu+ pi0 with an enlarged fiducial volume in 0.37 megaton-years exposure of Super-Kamiokande", "Takenaka", "Phys. Rev. D", 102, 112011, 2020),
    (16, "Direct neutrino-mass measurement based on 259 days of KATRIN data", "Aker", "Science", 388, 180, 2025),
    (17, "The fate of hints: updated global analysis of three-flavor neutrino oscillations", "Esteban", "JHEP", 2020, 178, 2020),
    (18, "The Quantum Theory of Fields", "Weinberg", "book", None, None, 1995),
    (19, "A no-go theorem for regularizing chiral fermions", "Nielsen", "Phys. Lett. B", 105, 219, 1981),
    (20, "Photons and gravitons in perturbation theory: Derivation of Maxwell's and Einstein's equations", "Weinberg", "Phys. Rev.", 138, 988, 1965),
    (21, "Limits on massless particles", "Weinberg", "Phys. Lett. B", 96, 59, 1980),
    (22, "Universality and scaling in gravitational collapse of a massless scalar field", "Choptuik", "Phys. Rev. Lett.", 70, 9, 1993),
    (23, "Blow-up of solutions of nonlinear wave equations in three space dimensions", "John", "Manuscripta Math.", 28, 235, 1979),
    (24, "Violation of CP invariance, C asymmetry, and baryon asymmetry of the universe", "Sakharov", "JETP Lett.", 5, 24, 1967),
    (25, "Cosmological baryon and lepton number in the presence of electroweak fermion interactions", "Harvey", "Phys. Rev. D", 42, 3344, 1990),
    (26, "The statistical theory of anomalous fermion number nonconservation", "Khlebnikov", "Nucl. Phys. B", 308, 885, 1988),
    (27, "Tri-bimaximal mixing and the neutrino oscillation data", "Harrison", "Phys. Lett. B", 530, 167, 2002),
    (28, "Confinement of quarks", "Wilson", "Phys. Rev. D", 10, 2445, 1974),
    (29, "On the quantum structure of a black hole", "'t Hooft", "Nucl. Phys. B", 256, 727, 1985),
    (30, "The Cellular Automaton Interpretation of Quantum Mechanics", "'t Hooft", "book", None, None, 2016),
    (31, "The minimal scenario of leptogenesis", "Blanchet", "New J. Phys.", 14, 125012, 2012),
    (32, "A lower bound on the right-handed neutrino mass from leptogenesis", "Davidson", "Phys. Lett. B", 535, 25, 2002),
    (33, "General relativity without coordinates", "Regge", "Il Nuovo Cimento", 19, 558, 1961),
    (34, "Constraint damping in the Z4 formulation and harmonic gauge", "Gundlach", "Class. Quantum Grav.", 22, 3767, 2005),
    (35, "Black holes and entropy", "Bekenstein", "Phys. Rev. D", 7, 2333, 1973),
    (36, "The world as a hologram", "Susskind", "J. Math. Phys.", 36, 6377, 1995),
    (37, "The cosmological constant problem", "Weinberg", "Rev. Mod. Phys.", 61, 1, 1989),
    (38, "The glueball spectrum from an anisotropic lattice study", "Morningstar", "Phys. Rev. D", 60, 34509, 1999),
    (39, "Baryogenesis through leptogenesis", "Barbieri", "Nucl. Phys. B", 575, 61, 2000),
    (40, "Resonant leptogenesis", "Pilaftsis", "Nucl. Phys. B", 692, 303, 2004),
    (41, "The Early Universe", "Kolb", "book", None, None, 1990),
    (42, "Freeze-in production of FIMP dark matter", "Hall", "JHEP", 2010, 80, 2010),
    (43, "Weak-interaction mixing in the six-quark theory", "Fritzsch", "Phys. Lett. B", 73, 317, 1978),
    (44, "Parametrization of the Kobayashi-Maskawa matrix", "Wolfenstein", "Phys. Rev. Lett.", 51, 1945, 1983),
    (45, "Supersymmetric dark matter", "Jungman", "Phys. Rep.", 267, 195, 1996),
    (46, "Two-loop renormalization factors of dimension-six proton decay operators in the minimal supersymmetric SU(5) gauge theory", "Hisano", "Phys. Lett. B", 724, 283, 2013),
    (47, "Demonstration of single-electron buildup of an interference pattern", "Tonomura", "Am. J. Phys.", 57, 117, 1989),
    (48, "A method for simulating chiral fermions on the lattice", "Kaplan", "Phys. Lett. B", 288, 342, 1992),
    (49, "Two-loop renormalization group equations in a general quantum field theory. I. Proof of the equivalence theorem", "Machacek", "Nucl. Phys. B", 222, 83, 1983),
    (50, "NuFIT 6.0: Updated global analysis of three-flavor neutrino oscillations", "Esteban", "JHEP", 2024, 216, 2024),
    (51, "Proton stability in grand unified theories, in strings and in branes", "Nath", "Phys. Rep.", 441, 191, 2007),
    (52, "Search for proton decay via p to nu K+ using 260 kiloton year data of Super-Kamiokande", "Abe", "Phys. Rev. D", 90, 72005, 2014),
    (53, "Dark Matter Search Results from 4.2 Tonne-Years of Exposure of the LUX-ZEPLIN (LZ) Experiment", "Aalbers", "Phys. Rev. Lett.", 135, 11802, 2025),
]


def crossref(query, rows=4):
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.bibliographic": query, "rows": rows,
         "select": "title,author,container-title,volume,page,DOI,issued,type"})
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)["message"]["items"]
        except Exception as e:  # noqa: BLE001
            if attempt == 2:
                return {"error": str(e)}
            time.sleep(2 * (attempt + 1))
    return {"error": "unreachable"}


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def first_page(p):
    if not p:
        return None
    m = re.search(r"\d+", p)
    return int(m.group()) if m else None


def check(n, title, author, journal, volume, page, year):
    items = crossref(f"{title} {author}")
    if isinstance(items, dict):
        return {"n": n, "verdict": "API_ERROR", "error": items["error"]}
    best, score = None, 0.0
    nt, na = norm(title), norm(author)
    for it in items:
        t = norm((it.get("title") or [""])[0])
        if not t:
            continue
        # substring overlap of the claimed title against the hit title
        overlap = 0
        for w in re.findall(r"[a-z0-9]+", title.lower()):
            if len(w) > 3 and w in t:
                overlap += 1
        words = max(1, len([w for w in re.findall(r"[a-z0-9]+", title.lower()) if len(w) > 3]))
        s = overlap / words
        auths = " ".join(norm(a.get("family", "")) for a in it.get("author", []))
        if na and na in auths:
            s += 0.15
        yr = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        if yr == year:
            s += 0.10
        if s > score:
            score, best = s, it
    if best is None:
        return {"n": n, "verdict": "NO_HIT", "claim": title}
    hit_year = (best.get("issued", {}).get("date-parts") or [[None]])[0][0]
    hit_vol = best.get("volume")
    hit_pg = first_page(best.get("page"))
    hit_j = (best.get("container-title") or [""])[0]
    problems = []
    if score < 0.6:
        problems.append(f"title match weak ({score:.2f})")
    if volume is not None and str(hit_vol) != str(volume):
        problems.append(f"volume: claimed {volume}, Crossref {hit_vol}")
    if page is not None and hit_pg is not None and hit_pg != page:
        problems.append(f"page: claimed {page}, Crossref {hit_pg}")
    if hit_year != year:
        problems.append(f"year: claimed {year}, Crossref {hit_year}")
    return {
        "n": n,
        "verdict": "OK" if not problems else "MISMATCH",
        "problems": problems,
        "score": round(score, 2),
        "crossref_title": (best.get("title") or [""])[0],
        "crossref_journal": hit_j,
        "crossref_volume": hit_vol,
        "crossref_page": best.get("page"),
        "crossref_year": hit_year,
        "crossref_doi": best.get("DOI"),
        "crossref_authors": ", ".join(
            f"{a.get('given','')} {a.get('family','')}".strip() for a in best.get("author", [])[:6]),
    }


def main():
    results = []
    for c in CLAIMS:
        n, title, author, journal, volume, page, year = c
        if journal == "book":
            results.append({"n": n, "verdict": "BOOK (not Crossref-checked)",
                            "claim": f"{author}, {title} ({year})"})
            print(f"[{n:2d}] BOOK  {title}")
            continue
        r = check(*c)
        r["claim"] = f"{author}, {title}, {journal} {volume}, {page} ({year})"
        results.append(r)
        flag = "OK " if r["verdict"] == "OK" else "!! "
        print(f"[{n:2d}] {flag}{r['verdict']:9s} {r.get('crossref_title','')[:70]}")
        for p in r.get("problems", []):
            print(f"      - {p}")
        time.sleep(0.4)
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    ok = sum(1 for r in results if r["verdict"] == "OK")
    bad = [r for r in results if r["verdict"] not in ("OK", "BOOK (not Crossref-checked)")]
    print(f"\nOK: {ok}/{len(CLAIMS)}   flagged: {len(bad)}")
    for r in bad:
        print(f"  ref {r['n']}: {r['verdict']} {r.get('problems') or r.get('error','')}")


if __name__ == "__main__":
    main()
