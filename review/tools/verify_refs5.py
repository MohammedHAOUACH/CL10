#!/usr/bin/env python3
"""Fifth pass: doi.org resolution + Semantic Scholar title search for the
entries Crossref/arXiv fuzzy search could not settle."""
import json
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "CL10-bibliography-audit/1.0 (mailto:haouach.mohammed@gmail.com)"}


def jget(url, headers=None):
    h = dict(UA)
    if headers:
        h.update(headers)
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=h)
            with urllib.request.urlopen(req, timeout=35) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                return {"error": str(e)}
            time.sleep(2 * (attempt + 1))
    return {"error": "unreachable"}


def citeproc(doi):
    """doi.org content-negotiation -> CSL JSON (works for old/Elsevier items)."""
    d = jget("https://api.crossref.org/works/" + urllib.parse.quote(doi))
    if "message" in d:
        m = d["message"]
        return (f"{(m.get('title') or ['?'])[0][:78]} | {(m.get('container-title') or ['?'])[0][:34]} | "
                f"vol {m.get('volume')} | pg {m.get('page')} | "
                f"{(m.get('issued',{}).get('date-parts') or [[None]])[0][0]} | "
                f"{', '.join(a.get('family','') for a in m.get('author',[])[:6])}")
    return f"NOT FOUND ({d.get('error','')[:50]})"


print("=== E. candidate DOIs via Crossref ===")
CAND = {
    "13 PDG2024": "10.1103/PhysRevD.110.070001",
    "13b PDG2024 alt": "10.1103/PhysRevD.110.070001v1",
    "15 SK epi0 0.37mtyr": "10.1103/PhysRevD.102.012006",
    "15b SK alt 012006?": "10.1103/PhysRevD.102.012005",
    "19a NN PLB105-219": "10.1016/0370-2693(81)91221-9",
    "19a2 NN PLB105-219b": "10.1016/0370-2693(81)91222-0",
    "19b NN CMP78-589": "10.1007/BF01208281",
    "23 John CPAM23-377": "10.1002/cpa.3160230306",
    "24 Sakharov": "10.1070/JETP.1967v005p0024",
    "24b Sakharov alt": "10.1016/0370-2693(67)90091-2",
    "29 tHooft PLB364-357": "10.1016/0370-2693(95)01296-5",
    "29b tHooft PLB364-357b": "10.1016/0370-2693(95)01396-8",
    "31 Blanchet PRL110-191601": "10.1103/PhysRevLett.110.191601",
    "31b Alloul PRL": "10.1103/PhysRevLett.110.191601",
    "32 BHS NPB776-77": "10.1016/j.nuclphysb.2007.04.011",
    "12c DistlerFraser CMP161-195": "10.1007/BF02101611",
    "12c2 DistlerFraser": "10.1007/BF02099391",
    "7 FritzschMinkowski": "10.1016/0003-4916(75)90211-0",
}
for k, doi in CAND.items():
    print(f"[{k}] {doi}\n     {citeproc(doi)}")
    time.sleep(0.3)

print("\n=== F. Semantic Scholar title search ===")
TITLES = [
    (4, "An introduction to quantum cellular automata Arrighi Nesme Scherrer"),
    (5, "Equivalence of quantum cellular automata and circuits D'Ariano Perinotti Bahreyni Ford"),
    (13, "Review of Particle Physics 2024 Navas Physical Review D 110 070001"),
    (15, "Search for proton decay p to e+ pi0 p to mu+ pi0 0.37 megaton-years Super-Kamiokande"),
    (19, "A no-go theorem to avoid the fermion doubling in lattice gauge theory Nielsen Ninomiya"),
    (19, "A proof of the Nielsen-Ninomiya theorem Communications Mathematical Physics 1981"),
    (23, "Blow-up of solutions of nonlinear wave equations in three space dimensions John 1970 Communications Pure Applied Mathematics"),
    (24, "Violation of CP invariance C asymmetry baryogenesis of the universe Sakharov 1967 JETP Letters"),
    (29, "Deterministic quantum mechanics and field theories t Hooft Physics Letters B 1996"),
    (31, "Leptogenesis with quasi-degenerate neutrinos Alloul Ariga Blanchet Di Bari"),
    (32, "Minimal thermal leptogenesis the neutrino mass spectrum Barbieri Hall Strumia"),
    (12, "Heterotic fermions in noncommutative geometry Distler Fraser Communications Mathematical Physics 1994"),
]
for n, t in TITLES:
    d = jget("https://api.semanticscholar.org/graph/v1/paper/search?" + urllib.parse.urlencode(
        {"query": t, "limit": 4, "fields": "title,year,venue,externalIds,authors,publicationTypes"}))
    print(f"[{n}] {t[:60]}")
    for p in (d.get("data") or []):
        ids = p.get("externalIds") or {}
        doi = ids.get("DOI") or ids.get("arXiv") or "-"
        au = ", ".join(a.get("name", "") for a in (p.get("authors") or [])[:5])
        print(f"     {p.get('year')} | {doi} | {(p.get('venue') or '')[:45]} | {p.get('title','')[:70]}")
        print(f"        {au[:110]}")
    if not d.get("data"):
        print(f"     (none / {d.get('error','')[:40]})")
    time.sleep(1.5)
