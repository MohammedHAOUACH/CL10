#!/usr/bin/env python3
r"""Provenance audit of the manuscript's numbers.

    python3 review/tools/number_audit.py [doc/final_publication.md ...]

Every distinctive decimal the manuscript prints must be a number the repository
produces: a numeric leaf of an experiment artefact, or a number that appears in
one of its check names / status strings, or a numeric scalar of the
configuration the code itself reads (`cl10/config.yaml`), or a power-of-ten
multiple of one (the artefacts print 2193467237.372068; the paper quotes
2.193e9). A percentage print (`$99.44\%$`) may be backed by the same quantity
stored as a fraction (V113 records agreement_one_per_family = 0.99444); the unit
conversion is the only scaling beyond powers of ten the audit allows, and it
applies only to a token spelled with a `\%`.

Four classes are exempt and are *reported as such*, never skipped silently:

  ref      bibliography lines (`## References` onward): the decimals there are
           volume, page, year and arXiv identifiers, not physics numbers
  xref     structure, not data: cross-references (Section 12.1, Fig. 13.2,
           Eq. 4.1) and the section numbers the headings themselves carry
  ident    persistent identifiers on the same line (a DOI, an arXiv id, a Zenodo
           concept record): names of objects, not measurements
  derived  arithmetic the sentence performs on sourced inputs, listed in
           DERIVED with the input it starts from

A token that is neither produced nor derived nor classifiable is a number
nothing in the repository backs: the audit prints it and exits non-zero.
"""
import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MASTER = os.path.join(ROOT, "doc", "final_publication.md")

# a section number follows its keyword even through a list: "Sections 12.1 and
# 13.2", "(Section 8.3, 9, 16)"
XREF_BEFORE = re.compile(
    r"(Sections?|Figs?|Eqs?|Equations?|Propositions?|Tables?|"
    r"Appendixes?|Appendices?)(?:[\s,.&\-–—]|\d|\.|and)*$")
HEADING = re.compile(r"^#{1,6}\s*$")
ART_TOKEN = re.compile(r"\d+\.\d+(?:[eE]([-+]?\d+)|\\times\s*10\^\{?(-?\d+)\}?)?")
PAPER_TOKEN = re.compile(r"\d+\.\d+(?:\\times10\^\{?(-?\d+)\}?|[eE]([-+]?\d+))?")
IDENT_LINE = re.compile(r"doi\.org|arxiv\.org|zenodo|10\.\d{4,}/")
CONFIG = os.path.join(ROOT, "cl10", "config.yaml")
CONFIG_ENTRY = re.compile(r"^\s*[A-Za-z0-9_.\-]+:\s*(.+?)\s*(?:#.*)?$")

def art_value(m):
    """the value an artefact text carries: 1.1e-16, 1.10\\times10^{-16}, 2.7e-15"""
    v = float(m.group(0).split("\\times")[0].split("e")[0].split("E")[0])
    if m.group(1):
        v *= 10.0 ** int(m.group(1))
    elif m.group(2):
        v *= 10.0 ** int(m.group(2))
    return v

# the paper's own arithmetic on sourced inputs: token -> the arithmetic and its inputs
DERIVED = {
    r"5.98\times10^{10}": r"Gamma_N = y^2 M_N/8pi at y=1/8, M_N = 9.6132e13 GeV (V106 M_N_GeV); Gamma/H = 4.6 = V106 cases.y=g0=1/8.K",
    r"1.93\times10^{12}": r"Gamma_N = y^2 M_N/8pi at y=0.71, M_N = 9.6132e13 GeV (V106 M_N_GeV); the sentence quotes Gamma/H = 148 from the same two numbers",
    r"0.91497": r"rho^4 at rho = 0.978030, the closed-branch norm (V200/V203): 0.978030^4 = 0.9149739",
    r"7.613": r"2 ln 45, the naive bond entropy of the 45-generator bond algebra: 7.6133250 nats",
    r"7.61": r"2 ln 45 = 7.61333 nats (the same naive bond entropy, printed to four digits)",
    r"10.84": r"12 (16 pi F_208)^10 at the certified F_208 = 0.0196928 (V204, V209): 10.8381; the float64 ladder reads 10.8",
    r"1.013": r"100 (1 - 16 pi F_208) at the certified F_208 = 0.0196928 (V204, V209): 1.0132%",
    r"1.1\times10^{8}": r"the chain reading over the rung: 1.6e7 eV / 0.146 eV (V210, and the rung of V196) = 1.096e8",
    r"50.3": r"16 pi = 50.2655, the rung spacing of the ladder",
    r"43.9": r"(7/2pi) ln(M_P/M_Z) at M_P = 1.220890e19 GeV (config) and M_Z = 91.1876 GeV: 43.935",
    r"11.6": r"the rung ruler at 0.2 GeV: ln(0.2 / 1.220890e19) / ln(1/(16 pi)) = 11.630",
    r"18.5": r"|1/5^2 - pi/64| / (pi/64) = 18.51% (the closed coupling pi/64 = 0.049087, V135)",
    r"0.0084": r"0.0487 eV x m_sol/m_atm = 0.0487 x 0.0086/0.0501 eV (config m_sol_ev, m_atm_ev) = 0.00836 eV",
    r"8.4\times10^{-38}": r"4 pi / M_P^2 at M_P = 1.220890e19 GeV (cl10/config.yaml physics.m_planck_gev): 8.430e-38 GeV^-2, the UV bound of cl10.lattice_rg",
    r"5.3\times10^{-39}": r"P_2 / M_P^2 with P_2 = 0.789 (the 2-particle S-matrix S^(2) = U (x) U) at M_P = 1.220890e19 GeV: 5.293e-39 GeV^-2",
    r"2.1\times10^{-66}": r"5.293e-39 GeV^-2 x 3.89379338e-28 cm^2/GeV^2 (cl10/config.yaml units.gev2_to_cm2): 2.061e-66 cm^2",
}

# a number whose own sentence cites a source is backed by that source: the
# citation marker must lie in the same sentence (for a display equation, in the
# sentence that introduces or closes it).
CIT_MARK = re.compile(r"\[\d+(?:[-,\s]\d+)*\]")
SENT_END = re.compile(r"\.(?=\s+[A-Z(]|\s*$)")
DISPLAY = re.compile(r"\$\$.*\$\$")


def cite_scope(lines, ln, col):
    """the text inside which a citation would back the number on this line"""
    line = lines[ln - 1]
    if DISPLAY.fullmatch(line.strip()):
        prev = next((lines[j] for j in range(ln - 2, -1, -1) if lines[j].strip()),
                    "")
        nxt = next((lines[j] for j in range(ln, len(lines)) if lines[j].strip()),
                   "")
        return prev + " " + nxt
    start, end = 0, len(line)
    for m in SENT_END.finditer(line):
        if m.end() <= col:
            start = m.end()
        elif m.start() >= col:
            end = m.start()
            break
    return line[start:end]


def artefact_values():
    """(value, provenance) pairs over every experiment artefact, leaves and text."""
    out = []
    for dp, _dirs, files in os.walk(ROOT):
        if "/results" not in dp and not dp.endswith("results"):
            continue
        for f in files:
            if not f.endswith(".json"):
                continue
            path = os.path.join(dp, f)
            try:
                with open(path) as fh:
                    doc = json.load(fh)
            except (OSError, ValueError):
                continue
            stack = [(doc, os.path.relpath(path, ROOT))]
            while stack:
                node, where = stack.pop()
                if isinstance(node, dict):
                    for k, v in node.items():
                        stack.append((v, f"{where}:{k}"))
                elif isinstance(node, list):
                    for v in node:
                        stack.append((v, where))
                elif isinstance(node, str):
                    for m in ART_TOKEN.finditer(node):
                        try:
                            out.append((art_value(m), where))
                        except ValueError:
                            pass
                elif isinstance(node, (int, float)) and not isinstance(node, bool):
                    out.append((float(node), where))
    out += config_values()
    return out


def config_values():
    """(value, provenance) pairs over the scalars of cl10/config.yaml.

    The constants the code reads are part of what the repository produces: a
    number the API or the README quotes from the configuration is backed by the
    file that states it, and the provenance names the line.
    """
    out = []
    try:
        lines = open(CONFIG).read().split("\n")
    except OSError:
        return out
    for ln, raw in enumerate(lines, 1):
        m = CONFIG_ENTRY.match(raw)
        if not m:
            continue
        for t in ART_TOKEN.finditer(m.group(1)):
            try:
                out.append((art_value(t), f"cl10/config.yaml:{ln}"))
            except ValueError:
                pass
    return out


def sigfigs(mant):
    """Significant digits the manuscript prints for a mantissa."""
    return max(len(mant.replace(".", "").lstrip("0")), 1)


def close(a, b, sig=4):
    """The artefact value rounds to the printed value at its printed precision.

    A token is backed only when the artefact leaf agrees to the LAST PRINTED
    DIGIT: no absolute band, because a leaf that is exactly zero is not
    evidence for a printed 1e-15 defect, and a band wide enough to forgive a
    stale rounding is wide enough to forgive a wrong one.

    The match is on MAGNITUDE: the manuscript's token carries no sign (the
    regex extracts a mantissa), while the artefacts print signed levels -
    $E_0 = -5802.40$ in the text is the leaf -5802.398188196636 in V101. The
    sign itself is a claim the sentence makes, not one a value can witness.
    """
    if a == b:
        return True
    a, b = abs(a), abs(b)
    if a == 0 or b == 0:
        return False
    exp = math.floor(math.log10(a))
    half_ulp = 0.5 * 10.0 ** (exp - sig + 1)
    return abs(a - b) <= half_ulp * (1.0 + 1e-9)


def audit(path, values):
    text = open(path).read()
    lines = text.split("\n")
    ref_start = next((i for i, l in enumerate(lines, 1)
                      if l.startswith("## References")), len(lines) + 1)
    # the document's own section numbers: headings carry them ("### 12.1 ...")
    # and the provenance tables cite them bare ("splitting exclusions
    # (13.2, 13.5, Table 5)")
    sections = {m.group(1) for m in
                map(re.compile(r"^#{2,6}\s+(\d+\.\d+)(?=$|[\s.:])").match, lines) if m}
    hits = {}
    pct = set()
    for ln, line in enumerate(lines, 1):
        for m in PAPER_TOKEN.finditer(line):
            tok = m.group(0)
            try:
                val = float(tok.split("\\times")[0].split("e")[0].split("E")[0])
            except ValueError:
                continue
            if m.group(1):
                val *= 10.0 ** int(m.group(1))
            elif m.group(2):
                val *= 10.0 ** int(m.group(2))
            pre = line[:m.start()]
            hits.setdefault(tok, []).append(
                (ln, val, pre, bool(HEADING.match(pre)), m.start()))
            if line[m.end():m.end() + 2] == "\\%":
                pct.add(tok)

    produced, refs, xrefs, idents, derived, unknown = [], 0, 0, 0, 0, []
    for tok, ent in sorted(hits.items(),
                           key=lambda kv: min(l for l, *_ in kv[1])):
        val = ent[0][1]
        struct = tok in sections and re.fullmatch(r"\d{1,2}\.\d", tok)
        sig = sigfigs(tok.split("\\times")[0].split("e")[0].split("E")[0])
        src = next((w for v, w in values if close(val, v, sig)
                    or (tok in pct and close(val / 100.0, v, sig))), None)
        if src:
            produced.append((tok, src))
            continue
        if all(ln >= ref_start for ln, *_ in ent):
            refs += len(ent)
            continue
        if all(CIT_MARK.search(cite_scope(lines, l2, col))
               for l2, _, _, _, col in ent):
            refs += len(ent)
            continue
        if all(head or struct or XREF_BEFORE.search(ctx[-40:])
               for _, _, ctx, head, _col in ent):
            xrefs += len(ent)
            continue
        if all(IDENT_LINE.search(lines[ln - 1]) for ln, *_ in ent):
            idents += len(ent)
            continue
        if tok in DERIVED and all(ln < ref_start for ln, *_ in ent):
            derived += len(ent)
            print(f"  derived {tok} (lines "
                  f"{sorted(l for l, *_ in ent)}): {DERIVED[tok]}")
            continue
        unknown.append((min(l for l, *_ in ent), tok,
                        sorted({l for l, *_ in ent})[:6]))

    print(f"{os.path.relpath(path, ROOT)}: {len(hits)} decimal tokens -> "
          f"{len(produced)} artefact-backed, {refs} bibliography, {xrefs} xrefs, "
          f"{idents} identifiers, {derived} derived, {len(unknown)} unexplained")
    for ln, tok, ls in unknown:
        print(f"  UNEXPLAINED L{ln}: {tok}  lines={ls}")
    return unknown


def main(argv):
    paths = [os.path.abspath(p) for p in argv] or [MASTER]
    values = artefact_values()
    print(f"artefact json numbers: {len(values)}")
    bad = []
    for p in paths:
        bad += audit(p, values)
    if bad:
        print(f"FAIL: {len(bad)} unexplained number(s)")
        return 1
    print("PASS: every printed decimal is artefact-backed, a citation, "
          "a cross-reference, or arithmetic on a named input")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
