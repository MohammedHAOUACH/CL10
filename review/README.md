# review — audit trail of the publication

What was checked, when, by what method, and with which result. Nothing here
is part of the paper; everything here is reproducible evidence about it.

## Active documents

| file | what it is |
|---|---|
| `prepublication.md` | the current pre-publication audit (2026-10-06/07): every headline number re-derived from the code before compiling, the build checks (errors / undefined refs / floats / overfull boxes), the DOI re-verification. Supersedes the review rounds for everything bibliographic. |
| `refs_audit.md` | the bibliography audit report (52 entries): real errors found and the corrections applied, each with its Crossref / INSPIRE / arXiv proof. The definitive record for the references. |
| `refs_audit.json` | the machine-readable per-entry verdict produced by `tools/verify_refs.py`. |
| `measured_values_audit.md` | the 2026-10-09 audit of every *measured* value in Fig. 10 and Table 4 against its primary source (PDG 2024, Planck 2018, NuFIT 6.0, Super-Kamiokande): verdicts, the corrections applied to code and manuscript, and the one flagged inconsistency (the $6.6\times10^{33}$ yr $K^+$-mode number in the Fig. 16 caption). |

## tools/ — the bibliography verifiers

Six one-shot scripts that check the bibliography against external sources
(Crossref DOI + title/author/volume/page matching, INSPIRE HEP, the arXiv
API). They are kept because the audit is re-runnable:

```bash
python3 review/tools/verify_refs.py    # full pass -> review/refs_audit.json
python3 review/tools/verify_refs2.py   # the 46-DOI re-check cited by the audits
python3 review/tools/verify_refs3..6.py  # targeted re-checks (the entries the first pass flagged)
```

`verify_refs.py` writes `review/refs_audit.json` (one level up from `tools/`).
The scripts encode the claims **as printed in `doc/publication.md`**; if the
bibliography changes, update the claim lists before re-running.

## archive/ — closed review rounds

Sequential review cycles, kept for the trace, not for reference:

| file | round |
|---|---|
| `review2.md` | the critique |
| `review2_response.md` | the response |
| `review3.md` | the independent verification |
| `review4.md` | the final audit of that cycle (2026-10-06, 62/62 PASS) — superseded by `prepublication.md`; it predates the 102-experiment suite and the 36-figure renumbering, so its figure numbers and test counts are historical. |
