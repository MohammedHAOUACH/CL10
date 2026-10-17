# archive/review — the audit trail of the *previous* manuscript

Evidence about `archive/doc/publication.md`, the superseded manuscript: what was
checked, when, by what method, with which result. Nothing here is part of the
current paper — the live audit trail is
[`review/`](../../review/README.md), and the build tooling is
[`review/tools/`](../../review/README.md).

| file | what it is |
|---|---|
| `prepublication.md` | the pre-publication audit (2026-10-06/07): every headline number re-derived from the code before compiling, the build checks (errors / undefined refs / floats / overfull boxes), the DOI re-verification. |
| `refs_audit.md` | the bibliography audit report (52 entries): the real errors found and the corrections applied, each with its Crossref / INSPIRE / arXiv proof. |
| `refs_audit.json` | the machine-readable per-entry verdict of that audit. |
| `refs_audit_direct.json` | the direct-lookup pass over the same entries (no cached metadata). |
| `measured_values_audit.md` | the 2026-10-09 audit of every *measured* value in Fig. 10 and Table 4 against its primary source (PDG 2024, Planck 2018, NuFIT 6.0, Super-Kamiokande): verdicts, the corrections applied to code and manuscript, and the one flagged inconsistency (the $6.6\times10^{33}$ yr $K^+$-mode number in the Fig. 16 caption). |

The figure numbers, experiment counts and section numbers inside these files are
those of the **superseded** manuscript: `doc/final_publication.md` renumbered the
figures, consolidated the sections and relabelled the residuals W1–W9 → O1–O9.

The one-shot bibliography verifiers (`verify_refs*.py`) that produced
`refs_audit.json` were retired with that manuscript. The current bibliography is
kept honest at build time by `publication/fetch_dois.py`, which resolves and
verifies every entry against Crossref and writes the `doi.json` map the LaTeX
filter consumes — see [`publication/README.md`](../../publication/README.md).
