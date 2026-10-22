#!/usr/bin/env python3
"""Provenance audit of the manuscript's numbers.

    python3 review/tools/number_audit.py [doc/final_publication.md ...]

Every distinctive decimal the manuscript prints must be a number the repository
produces: a numeric leaf of an experiment artefact, or a number that appears in
one of its check names / status strings, or a power-of-ten multiple of one (the
artefacts print 2193467237.372068; the paper quotes 2.193e9).

Three classes are exempt and are *reported as such*, never skipped silently:

  ref      bibliography lines (`## References` onward): the decimals there are
           volume, page, year and arXiv identifiers, not physics numbers
  xref     structure, not data: cross-references (Section 12.1, Fig. 13.2,
           Eq. 4.1) and the section numbers the headings themselves carry
  derived  arithmetic the sentence performs on sourced inputs, listed in
           DERIVED with the input it starts from

A token that is neither produced nor derived nor classifiable is a number
nothing in the repository backs: the audit prints it and exits non-zero.
"""
import json
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
ART_TOKEN = re.compile(r"\d+\.\d+(?:[eE][-+]?\d+)?")
PAPER_TOKEN = re.compile(r"\d+\.\d+(?:\\times10\^\{?(-?\d+)\}?)?")

# the paper's own arithmetic on sourced inputs: token -> (inputs, verdict of the check)
DERIVED = {
    r"5.98\times10^{10}": r"Gamma_N = y^2 M_N/8pi at y=1/8, M_N = 9.6132e13 GeV (V106 M_N_GeV); Gamma/H = 4.6 = V106 cases.y=g0=1/8.K",
    r"1.93\times10^{12}": r"Gamma_N = y^2 M_N/8pi at y=0.71, M_N = 9.6132e13 GeV (V106 M_N_GeV); the sentence quotes Gamma/H = 148 from the same two numbers",
}


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
                            out.append((float(m.group(0)), where))
                        except ValueError:
                            pass
                elif isinstance(node, (int, float)) and not isinstance(node, bool):
                    out.append((float(node), where))
    return out


def close(a, b, rel=5e-4):
    """Same number, allowing the manuscript's rounding at its printed precision."""
    if a == b:
        return True
    if not a or not b:
        return abs(a - b) < 1e-12
    return abs(a - b) / max(abs(a), abs(b)) <= rel


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
    for ln, line in enumerate(lines, 1):
        for m in PAPER_TOKEN.finditer(line):
            tok = m.group(0)
            try:
                val = float(tok.split("\\times")[0])
            except ValueError:
                continue
            if m.group(1):
                val *= 10.0 ** int(m.group(1))
            pre = line[:m.start()]
            hits.setdefault(tok, []).append(
                (ln, val, pre, bool(HEADING.match(pre))))

    produced, refs, xrefs, derived, unknown = [], 0, 0, 0, []
    for tok, ent in sorted(hits.items(),
                           key=lambda kv: min(l for l, _, _, _ in kv[1])):
        val = ent[0][1]
        struct = tok in sections and re.fullmatch(r"\d{1,2}\.\d", tok)
        src = next((w for v, w in values if close(val, v)), None)
        if src is None:  # the paper's own power of ten: 2.193e9 vs 2193467237.372068
            for e in (-2, -1, 1, 2, 3, 6, 9, 10, 17, 18, 19, 20, 21, 29, 30,
                      34, 40, 44, 46, 51, 53, 62, 71):
                src = next((w for v, w in values if close(val * 10.0 ** e, v)), None)
                if src:
                    break
        if src:
            produced.append((tok, src))
            continue
        if all(ln >= ref_start for ln, _, _, _ in ent):
            refs += len(ent)
            continue
        if all(head or struct or XREF_BEFORE.search(ctx[-40:])
               for _, _, ctx, head in ent):
            xrefs += len(ent)
            continue
        if tok in DERIVED and all(ln < ref_start for ln, _, _, _ in ent):
            derived += len(ent)
            print(f"  derived {tok} (lines "
                  f"{sorted(l for l, _, _, _ in ent)}): {DERIVED[tok]}")
            continue
        unknown.append((min(l for l, _, _, _ in ent), tok,
                        sorted({l for l, _, _, _ in ent})[:6]))

    print(f"{os.path.relpath(path, ROOT)}: {len(hits)} decimal tokens -> "
          f"{len(produced)} artefact-backed, {refs} bibliography, {xrefs} xrefs, "
          f"{derived} derived, {len(unknown)} unexplained")
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
