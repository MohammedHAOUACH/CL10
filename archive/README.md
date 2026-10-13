# Archive

Material the final paper (`doc/final_publication.md`) no longer leans on, kept so
that every historical number stays traceable. Nothing here is run by
`run_validation.py`; the archived scripts still use the pre-categorisation path
conventions (`ROOT` two levels up, sibling experiments by flat name), so they are
kept for reading, not for re-execution.

What supersedes this material lives in [../experiments/](../experiments/)
(catalogue: [experiments/README.md](../experiments/README.md)) and
[../doc/final_publication.md](../doc/final_publication.md). Back to the
[CL10 README](../README.md).

## `experiments/` — 35 superseded experiments

Each row names the experiment that now carries the result.

| id | what it established | superseded by |
|---|---|---|
| `V5_hierarchy` | the 0-parameter hierarchy statement | `V53` (the exponent), `V108` (the coefficient) |
| `V6_gravity` | the $(1{+}3)$ metric dynamics of the step | `V37`, `V50` |
| `V8_dark_matter` | the Planck-scale stable soliton as the carrier | `V54` (the ladder mass), `V66` (no carrier in the spectrum) |
| `V9_dark_energy` | the vacuum sum reduced by the $\mathbb{Z}_2$ cancellation | `V78` ($\Lambda = 0$ protected by $\det U = 1$) |
| `V10_antimatter` | the $16 \oplus \overline{16}$ CPT pair | `V24` |
| `V12_unification` | the $5/3$ normalisation | `V30`, then `V38`/`V64` |
| `V13_gauge_gravity` | matter, gauge and gravity in one object | `V50` |
| `V14_proton` | the $\tau_p$ prediction at $M_{\rm GUT} = M_P$ | `V64`, `V133` — that reading is abandoned in the final paper |
| `V15_seesaw` | the two-step seesaw | `V57` |
| `V16_families` | $3 = n_{\text{space}}$ as the family count | `V26` |
| `V17_w1_screening` | the screening mechanism for the vacuum energy | `V78` |
| `V18_w2_seesaw2step` | the seesaw as the W2 mechanism | `V57` |
| `V19_w3_families3` | the family count as the W3 mechanism | `V91` (the count is an input) |
| `V20_massgen` | the $\mathbf{10}_H$ mass term | `V28`, `V45` |
| `V21_running` | the 0-parameter GUT coupling and one-loop structure | `V38`, `V64` |
| `V22_baryon` | the step satisfies no Sakharov condition | `V55`, `V106`, `V117` |
| `V25_holography` | an area-scaling entropy claim from raw counting | `V97`–`V100` — raw counting is volume-law; the area law is a property of the gauge-invariant sector |
| `V27_qcd_confinement` | the linear potential and area law of the QCD sector | `V130` (the gauged layer's tension); Fig. 4 is computed live from `cl10.qcd` |
| `V29_yukawa_hierarchy` | $\sin^2\theta_{12} = 1/3$ from the democratic base | `V45` |
| `V30_unification` | the GUT-range running at $M_P$ | `V38`, `V64` |
| `V31_derive_gr` | Einstein–Hilbert, curvature, the black hole | `V50` |
| `V33_baryogenesis` | the named GUT ingredients of the asymmetry | `V55`, `V106` |
| `V34_cosmology` | the vacuum cancellation and the hierarchy value | `V78`, `V108` |
| `V35_cosmology_unified` | dark matter, dark energy and antimatter from one object | `V54`, `V78`, `V24` |
| `V39_families_qca` | the three families inside the QCA dynamics | `V26`, `V91` |
| `V40_leptogenesis` | the asymmetry at the right order | `V55`, `V106` (the mechanism fails) |
| `V42_dm_scale` | the Planck-scale candidate and its abundance gap | `V54`, `V123` |
| `V44_lattice_rg` | the UV-finite GUT cross-section | `V72`, `V73` |
| `V46_nonlinear_gravity` | the 4-D nonlinear gravity and self-gravity | `V50` |
| `V47_2loop_uvfinite` | UV finiteness of the two-step and the two-particle sector | `V48` |
| `V58_uvr` | the discrete RG resolution of the one-loop non-unification | `V64` |
| `V60_tbm_regression` | the TBM regression of the democratic base | `V45` |
| `V61_hydrogen` | the $120^3$ hydrogen resonance check | `V62` (the $296^3$ production run) |
| `V71_w5w1_stability` | the mirror under the Wilson extension | `V109` |
| `V119_w9_double_multiplet_cone` | the double-multiplet cone search at the cubic corner | `V115`, `V118` — its numbers are not quoted in the final paper |

## `doc/publication.md`

The previous manuscript: the same theory told through its addenda, with the
experiment identifiers (V60–V137) in the body, the nine residuals labelled
W1–W9, and a separate code/protocol appendix. `doc/final_publication.md`
supersedes it: consolidated sections, residuals relabelled O1–O9, no experiment
identifiers in the text, figures placed where they are used.

## `publication/publication.pdf`

The compiled PDF of that manuscript (REVTeX 4.2, 45 pp.). The build pipeline
itself — `publication/build.sh`, `preprocess.py`, `filter.lua`, `template.tex`,
`fetch_dois.py` — stayed in `publication/` and now points at the final manuscript.

## `review/`

The audit trail of the previous manuscript: `prepublication.md` (the
pre-publication number re-derivation and build checks), `refs_audit.md` /
`refs_audit.json` / `refs_audit_direct.json` (the bibliography audit), and
`measured_values_audit.md` (the measured values against PDG/Planck/NuFIT/SK).
The re-runnable verifiers stayed in `review/tools/`.
