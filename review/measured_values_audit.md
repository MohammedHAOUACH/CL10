# Measured-values audit — 2026-10-09

Scope: every **measured** value that enters Fig. 10 (`fig10_sm_head2head`)
and Table 4 (Section 11.6, "the scoreboard"). Method: each value was
re-checked on 2026-10-09 against the primary experimental source
(PDG 2024, Planck 2018, NuFIT 6.0, Super-Kamiokande) via web lookup of the
published paper / PDG listing, not from memory. The CL10 side of every
comparison is untouched; only the reference values and their citations
were audited.

## 1. Verdict per quantity

| quantity | value used | primary source | verdict |
|---|---|---|---|
| baryon asymmetry $\Delta B$ | $(6.12\pm0.04)\times10^{-10}$ | Planck 2018 VI, A&A **641**, A6 (2020) [14]; PDG 2024 $Y_B = (8.59\pm0.04)\times10^{-11}$, i.e. $\eta = 6.12\times10^{-10}$ | **corrected** (was $6.1\pm0.04$, unrounded) |
| hierarchy $v/M_P$ | $246.22 / 1.2209\times10^{19} = 2.0167\times10^{-17}$ GeV | PDG 2024 [13]: $v = (246.21967\pm0.00026)$ GeV, $M_P = 1.2209\times10^{19}$ GeV | **corrected** (was $2.016\times10^{-17}$ from the paper's own rounded $246/1.22\times10^{19}$; the internal convention $246/1.22\times10^{19}$ is kept in Section 6 where it defines the CL10 ratio) |
| neutrino mass scale $m_\nu$ | $\sqrt{|\Delta m^2_{3\ell}|} = 0.0501$ eV (normal ordering) | NuFIT 6.0, arXiv:2410.05380 [50]: $|\Delta m^2_{3\ell}| = 2.512\times10^{-3}$ eV$^2$ | **corrected** (was $0.050$) |
| neutrino mass bound $\sum m_\nu$ | $< 0.12$ eV (95% CL) | Planck 2018 + BAO [14] | confirmed, unchanged (DESI-era combinations are tighter but model-dependent; the paper keeps the Planck+BAO number) |
| proton lifetime bound $\tau_p/B$ | $> 2.4\times10^{34}$ yr ($p\to e^+\pi^0$) | Super-Kamiokande, PRD **102**, 112011 (2020) [15] | confirmed, unchanged; it is the strongest of the three bounds and matches the CL10-dominant channel |
| proton bound, $K^+$ mode | $> 5.9\times10^{33}$ yr ($p\to\nu K^+$) | Super-Kamiokande, PRD **90**, 072005 (2014) [52] | **corrected**: the Fig. 16 caption and the figure code carried the older SK-I number $6.6\times10^{33}$; now $5.9\times10^{33}$ |
| $\Omega_{\rm DM}h^2$ | $0.120$ | Planck 2018 VI [14] | confirmed, unchanged |
| $\alpha_i(M_Z)$ | $0.0169,\ 0.0338,\ 0.118$ | PDG 2024 [13] | confirmed, unchanged |
| dark-matter mass | no measured value exists | — | the grey bar in Fig. 10 is the **weak-scale reference** (100 GeV), now labelled "weak scale", not "SM" |

## 2. Changes applied

- `cl10/eta_cp.py`: `DELTA_B_OBS = 6.1e-10` -> `6.12e-10` (Planck 2018).
  The CL10 product is unchanged ($6.8\times10^{-10}$); the ratio moves
  $1.115 \to 1.111$, so the paper text now reads "a ratio $1.11$" and
  "the $11\%$ agreement".
- `doc/publication.md`:
  - Section 11.1 and Table 4: $(6.12\pm0.04)\times10^{-10}$ [14];
    $0.0501$ eV [50]; $\sum m_\nu < 0.12$ eV [14];
    $2.0167\times10^{-17}$ [13]; $>2.4\times10^{34}$ yr [15].
  - Fig. 10 caption: every measured value now carries its source and
    reference number.
- `make_figures/scales.py` (`fig10_sm_head2head`): the grey bars use the
  audited numbers; each panel's left bar is labelled with its source
  (Planck / weak scale / PDG / NuFIT / Planck+BAO bound / SK bound);
  the suptitle names the four data sources.

## 3. Notes

1. The 2014 SK paper [52] sets $5.9\times10^{33}$ yr for $p\to\nu K^+$
   (260 kt·yr). The Fig. 16 caption and the `fig16_proton_channels` code
   carried the older SK-I number $6.6\times10^{33}$; both corrected to
   $5.9\times10^{33}$ in this pass.
2. No measured value exists for $m_{\rm DM}$; the panel compares CL10's
   118 GeV against the weak-scale reference only, and the note keeps
   "$\Omega h^2$ not computed (W7)".
3. The internal CL10 convention $v/M_P = 246/1.22\times10^{19}$
   (Section 6, Eq. before the ladder) is deliberately **not** changed:
   it is the theory's own rounded input, distinct from the measured
   $2.0167\times10^{-17}$ it is compared against.

## 4. Fig. 12 — the gauge couplings at M_Z (same-day addendum)

`cl10/unification.py` carries `A_Z = (0.0169, 0.0338, 0.118)`. Re-checked
2026-10-09 against the PDG 2024 constants table
(`rpp2024-rev-phys-constants.pdf`):

| coupling | value used | PDG 2024 input | derived | verdict |
|---|---|---|---|---|
| $\alpha_1$ (GUT-normalised) | $0.0169$ | $\alpha^{-1}(M_Z)\simeq 127.95$, $\sin^2\hat\theta(M_Z) = 0.23129(4)$ ($\overline{\rm MS}$) | $\tfrac53\,\alpha/c^2 = 0.01695$, $\alpha_1^{-1} = 59.0$ | confirmed |
| $\alpha_2$ | $0.0338$ | same | $\alpha/s^2 = 0.03379$ | confirmed |
| $\alpha_3$ | $0.118$ | $\alpha_s(m_Z) = 0.1180(9)$ | — | confirmed |
| $M_Z$ | was $91.2$ GeV | $m_Z = 91.1880(20)$ GeV | — | **corrected** to $91.1876$ GeV (world average; effect on the running $\sim 10^{-4}$, plots unchanged) |

Changes: `cl10/unification.py` `M_Z = 91.1876`; `fig12_couplings_mz`
legend now reads "measured (PDG 2024)"; the Fig. 12 caption and the
Table 4 $\alpha_i(M_Z)$ row now cite [13] with the $\overline{\rm MS}$
inputs. No measured value in Fig. 12 needed correction beyond $M_Z$.

## 5. Fig. 12 CL10 side — hand cross-check of the down-run (no error)

The CL10 bars are not a bug. Hand 1-loop from $\alpha_{10} = 1/25$ at
$M_P$, $\ln(M_P/M_Z) = 39.435$:
$1/\alpha_i(M_Z) = 25 + (b_i/2\pi)\cdot 39.435$ gives
$(50.73,\ 5.13,\ -18.93)$, i.e. $(0.0197,\ 0.1951,\ -0.0528)$ — identical
to `running_down` to 4 digits. The reverse run of the measured values to
$M_P$ gives $(33.4,\ 49.5,\ 52.4)$: the SM couplings do not meet at
$M_P$ (the well-known 1-loop SM non-unification). The gap IS the named
residual W4; the existing experiments cover it: V38 (quantifies it),
V56 ($\alpha_{10} = 1/25$ is a root number, within 1.73 of the 1-loop
GUT value), V58 (UV-finite discrete RG), V64 (with $\alpha_X = 1/25$
fixed and the CL10-realizable desert content
$3\times\mathbf{10}+\overline{\mathbf{10}} + 1\times(1,2)$,
$\Delta b = (4.0, 6.67, 6.0)$, the measured couplings meet at
$M_X = 1.78\times10^{16}$ GeV with $\alpha_X = 0.040 = \alpha_{10}$
exactly; $\tau_p = 1.8\times10^{36}$ yr), V81/V82 (the composite content
is excluded as the desert supplier), V86 (the meeting is not yet derived
from the root). Fig. 12 shows the 0-parameter SM-content failure;
Fig. 21 shows the desert resolution.
