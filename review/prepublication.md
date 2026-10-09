# Pre-publication audit — 2026-10-06

> **Build note (2026-10-07):** the old `livrable/` pipeline (`make_paper.py`)
> has been removed and replaced by `publication/` — a pandoc + REVTeX 4.2
> build (`bash publication/build.sh`, needs `pandoc` and `texlive-publishers`).
> Current artifact: `publication/publication.pdf` (24 pages, `aps,prd,twocolumn`,
> 26 figures, 8 tables, 52 references; 0 errors, 0 undefined references,
> 0 stuck floats, 0 overfull boxes > 100pt). The DOI list that lived in
> `livrable/doi.json` is regenerated and re-verified against Crossref by
> `publication/fetch_dois.py` (-> `publication/doi.json`) and re-injected
> into the compiled bibliography as clickable https://doi.org/ links.

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

1. **Coupling values fail the down-running** (Sec. 11.1). With
   $\alpha_{10}=1/25$ at $M_P$ and the SM one-loop $\beta$-functions,
   $\alpha_2(M_Z)$ lands $5.8\times$ above its measured value and
   $\alpha_3$ hits a Landau pole at $2.2\times10^9$ GeV. **RESOLVED
   conditionally (2026-10-07, V64):** the SM field content between $M_Z$ and
   $M_X$ IS wrong — with the theory's own vectorlike remnant matter (three
   $10\oplus\overline{10}$ pairs at $1.2\times10^7$ GeV + one $(1,2)$ pair at
   $1.0\times10^{12}$ GeV) the measured couplings meet at $\alpha_X=1/25$
   **exactly** at $M_X=1.8\times10^{16}$ GeV; the spectrum is inferred, not
   derived, and $M_{\rm GUT}=M_P$ is abandoned ($\tau_p=1.8\times10^{36}$ yr,
   F5). The paper now states the repair and its cost (Sec. 11.1, W4, F2/F5).
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
`review/tools/verify_refs2.py`). The
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
- [x] Verification suite: 65 experiments validated on the GPU by the author
      (65/65 PASS, 2026-10-07, cupy backend, incl. V62 production 296³ and
      V64 desert); the numpy-backend run passes 64/65 (V62 needs the GPU)
- [x] Cover letter, journal order, arXiv listing, suggested referees:
      were in `livrable/cover_letter.md` (removed with the old build; rewrite
      against the current `publication/publication.pdf` before submission)
- [ ] ORCID: add it to `publication/template.tex` (`\author` block); nothing
      is invented in the paper
- [ ] Zenodo release + concept DOI added to the code-availability paragraph
- [ ] The GitHub repository the paper asserts must be public before submission
- [ ] Refs 13, 14, 15, 19, 23, 24, 29, 31, 32: confirm their DOIs by hand

## 5. Addendum 2026-10-09 — W2 re-audit (V103)

The V90 claim "the democratic neutrino point is an exact, unprotected,
measure-zero tuning" was re-derived and found to rest on an assumption
V90 never justified: a Dirac matrix proportional to the identity. The
implemented 0-parameter structure (V28/V45, `cl10/yukawa.py`) is the
democratic rank-1 $m_D = m_0 J$. Re-running the window scan with it
(`experiments/V103_w2_democratic`, PASS):

- $m_D^{\mathsf T} M^+ m_D = 3m_0^2 J/(2(w_0+\kappa))$ for every junction
  ratio $r = \kappa/w_0$ (dev $<10^{-12}$, rank-deficient $r = 2$ via the
  Moore–Penrose pseudo-inverse);
- the two massless neutrinos are protected by the rank-1 family symmetry
  for every $r$ — not only at $r = 2$;
- $m_3(r) = S_{17}/(1+r)$ is finite everywhere; the cosmological window
  is $r > 0.217$ (width $2.78$), and the atmospheric scale selects
  $r = 1.91$; $r = 2$ ($S_{17}/3 = 0.0487$ eV) is a junction-amplitude
  convention, not a tuning.

What does NOT change: the candidate still requires a democratic $m_D$
(else gapless coloured doublers, V76); the solar mass still needs the W6
splitting; the ladder value $0.146$ eV is still in tension.

Changes: `experiments/V103_w2_democratic` (new, PASS),
`make_figures/residuals.py` (fig30 redrawn with both Dirac assumptions),
`doc/publication.md` (Section 10.3, Table 4 $m_\nu$ row, W2 row, Fig. 30
caption, audit-trail paragraph, code-availability line, check count
102 → 104 — the count was already stale at 102 when V102 landed).

## 6. Addendum 2026-10-09 — W6 vev placements (V104)

The three mixing deviations ($\sin^2\theta_{12} = 1/3$ vs $0.308$,
$\theta_{13} \neq 0$, CKM = identity vs $\lambda = 0.225$) were tested
as ONE effect: the named W6 10-Higgs vev. `experiments/V104_w6_mixing_manifold`
(PASS, negative result):

- **Branch A (charged-lepton vev)**: democratic base + rank-$\le2$
  correction at the PDG spectrum. The determinant lemma (bilinear
  form, rank-1 update) reduces the constraint to
  $\sum_i (p_i^*)^2 e^{-it_i}/D_i = 1/\sum D$ with $p = W^T v$; the
  lightest family decouples from the democratic direction, and that
  direction IS the TBM solar column: $\sin^2\theta_{12} \le 0.005$,
  excluded at $23\sigma$.
- **Branch B (free-direction neutrino rank-1 vev)**: 4 parameters vs
  3 angles — covers the whole angle cube, predicts nothing.
- **Branch C (cubic-axis neutrino vev, the symmetry-natural choice)**:
  0-parameter prediction; 12 branches, best $20.3\sigma$ from NuFIT.

Net: the W6 residual is narrowed by computation — the vev must live
in the neutrino sector, with a non-cubic direction, and $\ge 2$ vev
terms in generic relative orientation are needed. No natural form of
the single named 10-Higgs input explains the deviations; they remain
named residuals, now with their mechanism-space bounded.

Changes: `experiments/V104_w6_mixing_manifold` (new, PASS),
`doc/publication.md` (W6 row, Table 4 $\sin^2\theta_{12}$ row,
audit-trail paragraph, check count 104 → 105).
