# Pre-publication audit — 2026-10-06

Artifact: `livrable/publication.pdf` (57 pages, REVTeX 4-2 `aps,prd,preprint`,
26 figures, 32 references, 16 with verified DOIs).
Build: `python3 livrable/make_paper.py --no-compile` then
`bash ~/.claude/skills/latex-document-skill/scripts/compile_latex.sh livrable/publication.tex --engine xelatex --auto-fix`.

`doc/publication.md` is the manuscript source and was **not modified**: the
journal abstract, the code-availability block and the corrected bibliography
entries are overrides in `livrable/` (`abstract.tex`, `code_availability.tex`,
`references_fix.json`, `doi.json`).

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

## 3. Bibliography audit (Crossref, title+volume+page+year)

**Verified and DOIs now printed in the paper** (16): 1, 2, 3, 6, 8, 10, 12,
16, 17, 20, 21, 22, 25, 26, 27, 28.

**Corrected in the deliverable** (`livrable/references_fix.json`; the source
text in `doc/publication.md` still carries the wrong values and must be fixed
there too):

| # | as written in doc/ | registered record |
|---|---|---|
| 3 | J. Stat. Phys. 85, **167** (1996) | 85, **551–574** (1996) |
| 6 | 'Unity of all elementary particles' | 'Unity of all elementary-**particle forces**' |
| 10 | Bull. Amer. Math. Soc. **69, 115 (1963)** | **Topology 3, 3–38 (1964)** |
| 12 | Can. J. Phys. **92, 1989 (2014)** | 92, **1501–1527** (2014) |
| 16 | Science 388, **1159** (2025) | 388, **180–185** (2025) |
| 25 | 'Covariant baryon number of the electroweak instanton/sphaleron processes' | the record at PRD **42, 3344** is 'Cosmological baryon and lepton number in the presence of electroweak fermion interactions' |
| 26 | Nucl. Phys. B **374, 37 (1992)** | **B 308, 885–912 (1988)** |

**Not verifiable through Crossref** (kept without DOI): 4, 5, 7, 9, 11, 13,
14, 15, 18, 19, 23, 24, 29, 30, 31, 32 — books, arXiv-only items, and older
or non-indexed journal records. Refs 13 (PDG), 14 (Planck), 15 (Super-K),
19 (Nielsen–Ninomiya), 24 (Sakharov), 29 ('t Hooft 1996), 31 (Blanchet),
32 (Barbieri–Hall–Strumia) are the ones a referee will check: each needs its
DOI confirmed by hand before submission.

## 4. Submission checklist

- [x] REVTeX 4-2, PRD article class, one author + affiliation
- [x] Journal-length abstract (250 words) + keywords
- [x] Figures embedded, numbered by the class, captions without manual numbers
- [x] Arabic section numbering matching the in-text "Section 8.4" references
- [x] Code-availability statement with repository URL and reproduction commands
- [x] MIT `LICENSE` at the repository root
- [ ] Repository URL is asserted in the paper but the GitHub repository is not
      yet public/archived — register a Zenodo release and put the concept DOI
      in `code_availability.tex`
- [ ] ORCID for the author (`\author[0000-0000-0000-0000]{...}`)
- [ ] Cover letter + suggested referees
- [ ] arXiv category (hep-ph primary; cross-list hep-th, gr-qc)
- [ ] The GPU suite (V52, V53, V55, V57, V58, V60, V61, V62) has not been
      re-run since the `run_validation.py` fix; their `result.json` verdicts
      are currently unverified and Appendix B quotes them
