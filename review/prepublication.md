# Pre-publication audit — 2026-10-06

> **Build note (2026-10-07):** the old `livrable/` pipeline (`make_paper.py`)
> has been removed and replaced by `publication/` — a pandoc + REVTeX 4.2
> build (`bash publication/build.sh`, needs `pandoc` and `texlive-publishers`).
> Current artifact: `publication/publication.pdf` (24 pages, `aps,prd,twocolumn`,
> 26 figures, 8 tables, 52 references; 0 errors, 0 undefined references,
> 0 stuck floats, 0 overfull boxes > 100pt). The DOI list that lived in
> `livrable/doi.json` is no longer injected into the bibliography; the
> verified DOIs remain recorded in `review/refs_audit.md` and
> `review/verify_refs2.py`.

Artifact at audit time: a 55-page REVTeX 4-2 `aps,prd,preprint` PDF,
26 figures, 52 references, 46 with verified DOIs. Compile was clean:
**0 errors, 0 undefined references, 0 stuck floats, 0 overfull boxes**
(two `xelatex` passes).

`doc/publication.md` is the manuscript source and now carries the journal
abstract, the repository URL and licence in the code-availability paragraph,
and the corrected bibliography: the earlier build-side overrides
(`abstract.tex`, `code_availability.tex`, `references_fix.json`) were removed
once the source was fixed.

## 1. Numbers re-derived from the code before compiling

| claim in the paper | recomputed value | source |
|---|---|---|
| $\dim\mathrm{Cl}(10,0)=1024$, module 32, 16+16 | 1024 / 32 / 16+16 | `cl10.gauge` |
| 45 generators, 12 SM bosons, 16 broken | 45 / 12 / 16 | `cl10.gauge` |
| $\mathrm{tr}(Y^2)=24$, $5/3$ normalisation | 24.000000000000004 | `cl10.gauge` |
| 8 anomaly coefficients | all $1.1\times10^{-13}$ | `cl10.gauge.anomaly_coefficients` |
| $\det U$, unitarity, causality, $\theta$-periodicity | $1+0i$, $1.1\cdot10^{-13}$, exact, $2\pi$ | `cl10.local_qca` |
| $\omega_{16}=\omega_{16\text{-bar}}$ | $1.1102230246251565\times10^{-16}$ | `cl10.brillouin.spectral_mirror_residual` |
| $\Delta_{\rm dyn}=1.1800000000000002$ | 1.18 | `cl10.uvrunning` |
| $v_{\rm CFT}=1/2$ | 0.5 | `simulation/young` |
| $m_W=79.98$ GeV, $m_H=139.5$ GeV | 79.98 / 139.5 | `cl10.unification` |
| $\Lambda_{\rm QCD}=220$ MeV, $f_\pi=129.7$ MeV | 220.0 / 129.7 | `cl10.qcd` |
| $\eta_B=5.108\times10^{-11}$, $\eta_B/m_\nu=3.5\times10^{-10}$ | same | `cl10.baryogenesis` |
| $n_B/n_\gamma=4.108\times10^{-9}$ | 4.108e-9 | `simulation/hydrogen` |
| $\Lambda=2.44\times10^{-69}$ GeV$^2$ | same | `experiments/V9_dark_energy` |
| $r_{6j}=1.11\times10^{-16}$ | same | `experiments/V44_lattice_rg` |

## 2. Physics blockers for a major journal

1. **U(1) Landau pole at 280 GeV** (Sec. 11.1). With $\alpha_{10}=1/25$ at
   $10^{18}$ GeV and the SM one-loop $\beta$-functions, $\alpha_1^{-1}$ reaches
   zero at 280 GeV. The paper states this; a referee will not accept a GUT
   claim with a pole below the $Z$ mass. Needs either a two-loop + threshold
   computation or an explicit statement that the SM field content between
   $M_Z$ and $M_X$ is wrong.
2. **Neutrino mass** (Sec. 12.2): $m_\nu=0.146$ eV vs $\Delta m_{31}^2$
   $\Rightarrow 0.050$ eV and $\sum m_\nu<0.12$ eV. 2.9× above the oscillation
   value; the paper already classifies the exact mass as W7.
3. **Thermal WIMP** (Sec. 13.1): $\Omega h^2\propto1/\langle\sigma v\rangle$
   with $\langle\sigma v\rangle$ fixed to the observed value is a
   normalisation, not a prediction; the paper says so.
4. **Baryogenesis** (Sec. 14.2): $\Gamma/H\simeq4.6$ means the heavy-neutrino
   decays are in equilibrium — out-of-equilibrium decay, the third Sakharov
   condition, is not realised.
5. **Doubler gap** (Sec. 8.4): the computed Brillouin spectrum is gapless at
   the twelve corner momenta; the $1.1\cdot10^{-16}$ residual is the
   matter/antimatter mirror, not a mass gap.

## 3. Bibliography audit (Crossref + INSPIRE + arXiv — full re-audit 2026-10-06)

**Superseded by `review/refs_audit.md`** (the definitive record). Summary: the
bibliography was expanded 32 → 52 entries; every entry verified against
Crossref/INSPIRE/arXiv; 11 real errors corrected at the source (wrong titles,
invented DOIs for Sakharov/John/Nielsen–Ninomiya, a non-existent D'Ariano JCP
paper, a wrong Gundlach reference, a wrong Hisano reference, wrong Regge
initial, wrong SK bounds); 46 verified DOIs (the list lived in
`livrable/doi.json`, now removed — the same DOIs are re-checked in
`review/verify_refs2.py`). The
in-text citations were previously rendered as literal `[n]` in the PDF
(`citations()` was never called from `inline()`); fixed in
the old `livrable/make_paper.py` — the compiled paper now has real `\cite{}` links,
all 52 refs cited, 0 undefined citations (the current `publication/` build
does the same at the pandoc-AST level).

## 4. Submission checklist

- [x] REVTeX 4-2, PRD article class, one author + affiliation
- [x] Journal-length abstract (266 words) + keywords, in the source
- [x] Figures embedded, numbered by the class, captions without manual numbers
- [x] Arabic section numbering matching the in-text "Section 8.4" references
- [x] Code-availability statement with repository URL, MIT licence and
      reproduction commands, in the source
- [x] MIT `LICENSE` at the repository root
- [x] Tables typeset with content-proportional column widths, `\scriptsize`
      when a cell would otherwise wrap into a tower, `\tabcolsep` tightened for
      4-column tables, no overfull boxes (see `column_widths`)
- [x] Verification suite: 64 experiments validated on the GPU by the author;
      the numpy-backend run reached 63/64 (V62 needs `matplotlib`, i.e. run the
      suite with `.venv/bin/python`)
- [x] Cover letter, journal order, arXiv listing, suggested referees:
      were in `livrable/cover_letter.md` (removed with the old build; rewrite
      against the current `publication/publication.pdf` before submission)
- [ ] ORCID: add it to `publication/template.tex` (`\author` block); nothing
      is invented in the paper
- [ ] Zenodo release + concept DOI added to the code-availability paragraph
- [ ] The GitHub repository the paper asserts must be public before submission
- [ ] Refs 13, 14, 15, 19, 23, 24, 29, 31, 32: confirm their DOIs by hand
