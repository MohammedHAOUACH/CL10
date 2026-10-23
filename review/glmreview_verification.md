# Verification of the referee report (`review/glmreview.md`)

Every claim of the report checked against the manuscript (`doc/final_publication.md`), the repository code, and new computations run during this verification. Verdicts: **CONFIRMED** (the report is right and the manuscript changed), **CONFIRMED-WRONG-DETAIL** (the concern is real, the report's arithmetic or citation is not), **REFUTED** (the report's claim does not survive checking). Actions reference the corrected master and the three companion papers (`doc/paper1_core.md`, `doc/paper2_gauge_gravity.md`, `doc/paper3_phenomenology.md`) built by `review/tools/build_series.py`.

---

## Major concerns

**M1 — title overstates the gravity claim. CONFIRMED.**
Table 1 itself records "one object: not met for the gravitational sector" ($g_{\mu\nu}$ and $G = 1$ are inputs; spin-2 declared; no massless spin-2 candidate in the implemented spectrum, Section 8.3). The old title asserted the opposite.
*Action:* master retitled "Causal Clifford Theory (CL10): a quantum cellular automaton on Cl(10,0) — Standard Model kinematics, exactly unitary causal dynamics, and a lattice gravitational sector"; the abstract's first tension sentence states "the gravitational sector adds $g_{\mu\nu}$ and $G = 1$ by hand". In the series the gravity word appears only in Paper II's title, with no "from a single algebra" claim.

**M2 — Lorentz failure at the nodes, softened in the abstract. CONFIRMED.**
Reproduced from the repository: V102 gives principal velocities $(0.200, 0.200, 0.284)$ per unit block momentum, anisotropy $0.748$–$0.758$ at all twelve nodes, $1200/1200$ linear dispersions; V115's exhaustive Wilson-family scan finds $0$ of $2726$ gapless points with $A < 10^{-6}$; the corner multiplet has commuting tangent velocity operators — three nodal planes.
*New computation (V142, `experiments/dynamics_spectrum/V142_w9_cone_isotropy` family):* frame averaging over the internal frames — the escape the report asks about in Question 4 — is blind: rebuilding the step in four random $\mathrm{SO}(10)$ frames leaves the node gap $2.003\times10^{-4}$, the principal velocities and the anisotropy $0.7545$ unchanged to $3\times10^{-15}$. The momentum dependence lives entirely in the sub-site factor of the Kronecker bond layer, so no internal-frame choice can move the cone.
*Action:* the abstract now names the anisotropy and the failed remedies (Wilson family, two bond layers, frame averaging); O9 and falsifier C3 (Paper I) state it as a named tension.

**M3 — unification conditional on six fitted inputs, abstract omits it. CONFIRMED.**
Verified against Section 13.1: SM-only running does not meet (spread $47\%$); the Wilson mechanism's own meeting is rigid at $M_X = 8.9\times10^{14}$ GeV (rank-1 Jacobian, split-invariance exact to $10^{-14}$), excluded by Super-Kamiokande ($M_X > 6\times10^{15}$ GeV); the $\alpha_X = 1/25$ solution needs the flavour-dependent bond layer with six $\kappa_s$, proven not to come from any on-site algebra element or class invariant.
*Action:* abstract states "holds only for a desert spectrum with two fitted thresholds and a flavour-dependent bond layer with six fitted couplings".

**M4 — coupling-convention inconsistency between Sections 9 and 13. CONFIRMED (the report's arithmetic partly wrong, the inconsistency real).**
Recomputed the desert inverse problem in the standard normalisation $\alpha = g^2/4\pi$ from the repository scripts: the three measured couplings $\alpha_i(M_Z) = (0.0169, 0.0338, 0.118)$ meet at $\alpha_X = 0.04005$, $M_X = 1.8\times10^{16}$ GeV — the meeting value is standard-normalised. The implemented angle gives $\alpha = (1/5)^2/4\pi = 0.00318$, a factor $12.6$ **below** the meeting value; the closed angle gives $\alpha = \pi^2/16/4\pi = \pi/64 = 0.0491$, a factor $1.226$ **above**. The old text's "exact match at exactly the root coupling" was the $\alpha_{10} \equiv g^2$ convention ($g^2 = 1/25 = 0.04 \approx 0.04005$) — a coincidence of conventions, not a physical match. The report's alternative ("the closed value is $g^2 = \pi^2/16$, a factor 15.4 larger") is its own convention slip: in the $g^2$ convention the meeting value is also $g^2$-normalised, and the consistent comparison is $0.617$ vs $0.405$, the same factor $1.54$; in the standard convention it is $1.226$. Either way the manuscript mixed the two, and the "factor 1.54" and "2.41×" statements inherited the mix.
*Action:* one convention stated once in Section 1.4 ($\alpha$ always $g^2/4\pi$; $\alpha_{10} \equiv g^2$ a named factor-$4\pi$ convention used only where explicit); Section 9.2 rewritten with the single-normalisation accounting ($1.226$ above / $12.6$ below); the convention-dependent "1.54/2.41" sentences deleted; F2 restated in standard normalisation with the closed-branch prediction $\alpha_X = \pi/64 = 0.0491$.

**M5 — "zero continuous unfixed numbers" contradicted by continuous inputs. CONFIRMED.**
The closure axiom quantizes only the on-site rotation angle. The construction retains: the Wilson strength $\theta = 1.731$ (nothing in the root fixes it — a fit selecting the split), the junction amplitude $r$ (window $r > 0.217$, atmospheric scale selects $r = 1.91$), the ordering of the $45$ plane-rotation factors (five shuffled orderings break the order-8 closure, deviation $\sim 1$), the split node set, and the six $\kappa_s$.
*Action:* the claim is now scoped — "zero continuous unfixed numbers **in the on-site update**" (abstract, Section 9 verdict); the audit row reads "free parameter 1 → 0 … the junction amplitude, the node set and the $\kappa_s$ remain inputs of the extended construction"; the Table 1 parameters row states "two conventions, eleven fitted coefficients, one conditional coefficient".

**M6 — closure is a postulate and worsens empirical contacts. CONFIRMED as description; the framing is a matter of presentation.**
The costs are exactly as the report lists: crossing $n^* = 10 \to 208$, desert destroyed ($\tau_p = 2.8\times10^{39}$ yr), coupling $1.226$ off the meeting, nodal lines replacing isolated nodes. The axiom is presented as an axiom (Section 9: "Call this the closure axiom"), not derived. Its non-trivial content, now stated: it produces new exact structure — the selection rule $F_n = 0$ for every odd $n$ (below $7\times10^{-15}$, no analogue at $g = 1/5$), a fundamental step of order $8$, and the finite algebraic spectrum $\{\pm\pi/4\,(10), \pm3\pi/4\,(6)\}$ — evidence the quantized branch is a different theory, not a reparametrisation.
*Action:* Section 9 states the axiom's status and both branches' numbers; the series separates the branches (Question 3 below).

**M7 — dark matter: mass without carrier. CONFIRMED.**
Section 12.2's no-go reproduced: $0$ survivors over $720$ Wilson angles; mirror-sector lightest states $1.3\times10^{8}$, $2\times10^{10}$ GeV; exactly $2$ of $1024$ algebra products commute with the step (the needed parity does not exist in the root); relic abundances $0.798$, $0.0246$, $2.3\times10^{-5}$–$7.6\times10^{-4}$ — none is $0.120$.
*Action:* abstract now says "no state of the implemented spectrum can be the dark matter" and calls $118$ GeV a mass *scale*; the one-equation status of $v = 2.08\, m_{\rm DM}$ is stated in Section 12.2.

**M8 — neutrino result is a candidate with fitted exponent. CONFIRMED.**
$n_{\text{nodes}} = 3$ is proven not to be a link-word length; $17 = 2\times10 - 3$ is an integer identity the algebra does not supply; $r = 1.91$ is selected by the atmospheric scale. The report's "a convention chosen to land on the measured value is a fit by another name" — accepted: the $r = 2$ democratic point is now labelled a candidate convention, and the audit counts the junction amplitude as an input.
*Action:* abstract wording "a candidate: the exponent $17$ is fitted and the junction amplitude is a free input inside a finite window".

**M9 — Sakharov-3 conflict belongs in the abstract. CONFIRMED.**
$\Gamma_N/H = 47$ at $y = 0.40$ vs required $y < 0.058$ reproduced from the repository.
*Action:* the tension list in the abstract (and Paper III's) now includes "the ladder's seesaw scale is incompatible with out-of-equilibrium decay".

**M10 — "machine precision" contradicts Table 3. CONFIRMED.**
Table 3 residuals $6.5\times10^{-7}$, $5.5\times10^{-5}$, $10^{-5}$ are $O(h^2)$ truncation errors; Section 1.4 defines machine precision as $\le 10^{-14}$.
*Action:* abstract now reads "reproduces general relativity to $O(h^2)$, not to machine precision", with the Kretschmann ratio $1.000009$ and maximum relative deviation $5.5\times10^{-5}$ quoted.

**M11 — determinism/Born rule. CONFIRMED.**
*Action:* abstract and Summary of results: "the projection and the Born rule are inputs of the interpretation (O9), not derived from the unitary step".

**M12 — spin-statistics postulate must be flagged. CONFIRMED.**
*Action:* abstract: "the Born rule, the spin–statistics identification and the split-node set are declared inputs"; Section 6 flags the identification as a postulate of the interpretation.

**M13 — format and scope: split into three papers. CONFIRMED.**
The master carried 52 lines over 1000 characters and table rows up to 2525.
*Action:* the series implements exactly the report's disposition: Paper I (Sections 1–9 + Appendix A), Paper II (Sections 10–11), Paper III (Sections 12–17 + audit + falsifiers), assembled mechanically by `review/tools/build_series.py` (verbatim body ranges, renumbered sections/figures/tables/references, cross-paper pointers). No line in any paper exceeds 2000 characters; the O2, O4, O8 and O9 rows, which the report measured at up to 2525 characters, are all below 2000 now.

**M14 — reproducibility: no independent verification channel. CONFIRMED as a limitation; partially answered by existing material.**
The Appendix A double construction (numeric $32\times32$ vs symbolic Clifford package, agreement as identities) and the exact $3^n$-word loop census ($55326$ of $3^{10}$) are the analytic cross-checks the report asks for on the loop computation; they were already in the manuscript and are now cited as such. The physical claim "the phase supplies the factor" rests on the Wilson phase value — now labelled an input in the audit (the $\theta = 1.731$ fit, M5).
*Action:* the availability statement says explicitly: numbers verified internally, pipeline not independently executed by a second party, independent re-execution invited.

---

## Minor points

1. **REFUTED-WRONG-DETAIL.** The report asks to correct "T. D'Ariano" to "A. D'Ariano". The Frontiers paper (Front. Phys. **12**, 120301 (2016)) is by **A. Bisio, G. M. D'Ariano, P. Perinotti** — the initial G. M. is right, the report's "A. D'Ariano" conflates Bisio's initial with D'Ariano's. The master's entry was nonetheless incomplete (Bisio missing): fixed to "A. Bisio, G. M. D'Ariano, and P. Perinotti".
2. **REFUTED.** The report doubts the Hisano author list. Verified against arXiv:1302.2194 and Phys. Lett. B **724**, 283 (2013): **J. Hisano, D. Kobayashi, Y. Muramatsu, N. Nagata** — the master is correct as written.
3. **VERIFIED-CORRECT.** Reference 13: Carlos Castro (Castro Perelman), *Can. J. Phys.* **92**, 1501 (2014) — "C. Castro" and the title stand.
4. **CONFIRMED.** Abstract was ~450 words in three paragraphs plus keywords. Split: a single ~250-word abstract + a "Summary of results" section; keywords moved to the end of the abstract block only in the master (dropped in the papers).
5. **CONFIRMED.** Three symbols for the angle/coupling. Fixed in Section 1.4: $g$ the rotation angle, $\alpha = g^2/4\pi$ the only coupling normalisation, $\alpha_{10} \equiv g^2$ a named convention.
6. **CONFIRMED.** The "nine-to-one" paragraph is now an explicit remark with the disclaimer that it is bookkeeping, not a reduction of difficulty.
7. **CONFIRMED.** Fig. 2 caption compressed to one input sentence.
8. **CONFIRMED.** The Section 12.2 paragraph and the O7 row were over 2000 characters; both now below. (Master line 855, the Discussion paragraph, remains 2103 characters but is not part of any submitted paper.)

9. **CONFIRMED.** The underscore leak is gone; "exactly 1" used throughout.
10. **CONFIRMED.** The "$16$-sector block-Brillouin spectrum" row compressed to a status line; the physics moved to text below the table.
11. **CONFIRMED.** No electroweak-mixing prediction is claimed in any abstract; $\sin^2\theta_W = 3/8$ appears only as a conventional tree-level relation in Section 13.2 with the honest running statement kept.
12. **CONFIRMED.** The abstract now says the gravitational sector "adds $g_{\mu\nu}$ and $G = 1$ by hand".
13. **CONFIRMED.** Availability statement now records the verification status exactly as the report phrases the limitation.

---

## Questions

1. **Normalisation of $\alpha = 1/25$:** standard, $\alpha = g^2/4\pi$ (the desert solution meets at $\alpha_X = 0.04005$). The old "exact match" was the $\alpha_{10} \equiv g^2$ coincidence. Recomputed consistently: closed angle $1.226\times$ above the meeting, implemented angle $12.6\times$ below. Fixed in Sections 1.4, 9.2, 13.1 and F2.
2. **$\theta = 1.731$ and $r$:** not fixed by the root. $\theta$ is a fit (it selects the split and the breaking chain); $r = 1.91$ is fixed by the atmospheric scale (a fit), $r = 2$ is a candidate convention, not established. Both are now counted as inputs in the audit and named in the abstract's scope sentence.
3. **Which theory:** both branches are theories of the same root, and the series separates them: Paper I carries the core and the closure axiom with both branches' measured costs; Paper III's phenomenology is computed at the implemented $g = 1/5$, with the closed branch's distinct predictions stated where they differ (F2: $\alpha_X = \pi/64$; F5: $\tau_p = 2.8\times10^{39}$ yr). No abstract mixes the branches.
4. **Frame averaging:** no. New experiment V142: averaging over four random $\mathrm{SO}(10)$ frames of the ten generators leaves gap, velocities and anisotropy unchanged to $3\times10^{-15}$ — the cone is fixed by the sub-site factor of the bond layer, not by the internal frame. The failure sits in the tension list (abstract, O9, C3).
5. **The 24-point orbit:** these are gapless residual doublers of the Wilson mechanism at the implemented angle, not the matter of the low-energy universe; the one-generation statement is a statement about the internal module (the $\mathbf{16}$ and its charges), and every abstract now displays the doubling next to it. The orbit is a named defect of the implemented angle (O9); at the closed angle the gapless set becomes six nodal lines and the per-node chirality bookkeeping must be redone (Section 9.2).
6. **Confinement:** the gauged layer is **not** the QCD confining sector — its string tension is $5.4\times10^{40}$ times QCD, it has no light vector sector (gap $64\,M_P$), and $\Lambda_{\rm QCD}$ is off the ladder ($n \simeq 11.6$, not an integer). The honest answer, now stated in Paper II Section 2.2 and Paper III: the theory does not currently confine quarks; this is a named tension, not a claim.
7. **Three families:** the identification $N_{\rm gen} = 3 = n_{\rm space}$ plus the three axis families of four nodes is an index, not a derivation — the count $3$ is a symmetry-breaking input (O3). The distinguishing observation exists: if the family number were the orbit count, the cubic-symmetric split (all $12$ nodes) would be forced, and its rigid meeting scale $M_X(12) = 2.9\times10^{15}$ GeV is excluded by proton decay (F5). So the observed three families are already distinguished from the twelve-node orbit reading by Super-Kamiokande; what is not distinguished is three families from three accidental families, and the paper says so.

---

## Summary of the disposition

The report's rejection rests on four fixable defects — the title (M1), the abstract/body gap (M2, M3, M7, M10, M11, M12), one inconsistent convention (M4), and the format (M13) — plus the scope judgement that the work is a research programme. All four are fixed, and the series realisation follows the report's own disposition. Two of the report's own corrections (minor 1, minor 2) do not survive checking against the cited papers; the manuscript's entries were right, with one author list completed.
