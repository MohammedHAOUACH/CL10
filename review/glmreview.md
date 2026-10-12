# Referee Report

**Journal:** Physical Review D (referee simulation, "grand journal physique")
**Manuscript:** *"The Universe as a Quantum Cellular Automaton: the Standard Model and Gravity from a Single Clifford Algebra"* (`doc/final_publication.md`, file reviewed in full)
**Author:** Mohammed Haouach
**Referee:** GLM (report recorded as `review/glmreview.md`)
**Length of manuscript reviewed:** 934 lines, ≈ 29,000 words, 25 figures, 8 tables, 55 references.

---

## Recommendation

**Not suitable for publication in its present form.** The manuscript contains genuinely original, carefully audited work, but as a submission to a major physics journal it fails on scope, on the central claim (a derivation of the Standard Model *and* gravity from one algebraic object), and on a number of internal inconsistencies that the abstract hides or softens relative to the body. I recommend **rejection with encouragement to restructure and resubmit** (see Section "Disposition" at the end).

---

## Summary of the manuscript

The paper proposes that the universe is a quantum cellular automaton whose cell is a 2×2×2 block of Planck sites, each of the eight sub-sites carrying the Clifford algebra Cl(10,0); the six-stage update is a product of Spin(10) plane rotations and swap-rotations. From this single object the author derives: a (1+3) Minkowski subalgebra; one SM generation as the 16-dimensional Weyl module with all eight anomaly coefficients vanishing; the 45 adjoint and the 5/3 hypercharge normalisation; an exactly unitary, strictly causal dynamics; an exact matter/antimatter spectral mirror; a gauged bond layer with an exact Gauss constraint and an area law; a lattice gravitational sector that numerically reproduces GR; a loop-factor ladder (1/16π)^n organising four scales; and a "closure axiom" that quantizes the single continuous coupling to g ∈ {π/2, π/4}.

The manuscript is unusually honest: it displays the lattice doubling rather than hiding it, reports seven tensions with observation, classifies every headline number by provenance (algebraic / conventional / fitted / open), lists nine open problems pinned to named quantities and seven falsifiers. This provenance discipline is, in my view, the manuscript's most distinctive virtue.

---

## Strengths

1. **The algebraic core is correct and cleanly presented.** Theorem 1 (Minkowski subalgebra), Theorem 2 (the eight anomaly coefficients — I verified the U(1)^3, [grav]²U(1) and mixed traces by hand), and Theorem 3 (the 5/3 normalisation) are standard but correct, and the paper correctly attributes them to the 1974–1981 GUT literature rather than claiming novelty.
2. **The dynamical results are sharp.** Exact unitarity with det U = 1 for every angle, the ±-paired spectrum and spectral mirror holding for every coupling, and the measured front speed of exactly one Planck cell per stage are clean, checkable statements, and the Lieb–Robinson replacement of the light cone after the Wilson deformation is properly cited [6].
3. **The honesty of the audit.** Seven tensions are reported rather than hidden (SU(2) overshoot 5.8×, SU(3) Landau pole at 2.2×10⁹ GeV, neutrino ladder excluded by Σm_ν < 0.12 eV, no dark-matter state in the implemented spectrum, computed relic abundance ≠ 0.120, string tension 5.4×10⁴⁰ × QCD, no massless spin-2 candidate). Very few manuscripts in this genre display their failures with this discipline.
4. **Numerical internal consistency is excellent.** I spot-checked ~20 headline numbers and they reproduce to the quoted precision: the hierarchy coefficient (3F₁₀¹⁰ = 2.04×10⁻¹⁷ vs v/M_P = 2.017×10⁻¹⁷), the ladder scales (M_P(1/16π)³ = 9.6×10¹³ GeV; M_P(1/16π)¹⁷ = 0.146 eV), the proton lifetime τ_p = M_X⁴/(α_X² m_p⁵) = 1.8×10³⁶ yr at M_X = 1.8×10¹⁶ GeV, the string-tension ratio 5.4×10⁴⁰, the exclusion scale M_X > 6×10¹⁵ GeV from the Super-Kamiokande bound, and the closure-angle bookkeeping (π/64 vs 1/25, factor 1.226; relative gap 1.128 vs 0.734). The claim in the availability statement that every number is reproducible from public scripts is credible, and the referee commendation for this is real.
5. **The falsifier list is concrete** (F3 neutrino pattern, F5 proton decay at Hyper-Kamiokande reach) and correctly identifies the most accessible tests.

---

## Major concerns

**M1 — The central claim of the title is not delivered, and the title overstates it.**
"A Standard Model and Gravity from a Single Clifford Algebra" is falsified by the paper's own Table 1: the gravitational sector adds a symmetric tensor field g_μν(x) and the value G = 1 by hand; diffeomorphism invariance is imposed, not derived; spin-2 is "declared, not derived"; and the one-particle spectrum provably contains no massless spin-2 candidate (Section 8.3: the bond factor carries T₂g but never E_g at any coupling). The parameter criterion is "not met" (two conventions, eleven fitted coefficients, one conditional coefficient). A paper that fails its own parameter criterion cannot carry this title without qualification. The abstract should state in its first sentence that gravity enters as a declared field on a subalgebra of the module, not as a derivation.

**M2 — The Lorentz problem is fatal for the continuum limit, and the abstract softens it.**
The computed low-energy spectrum is linear in all 1200 sampled directions but strongly anisotropic (principal velocities (0.200, 0.200, 0.284) per unit block momentum; anisotropy 0.748–0.758), the exhaustive scan over the Wilson family finds 0 of 2726 gapless points with A < 10⁻⁶, and at the exactly gapless corner multiplet the three tangent velocity operators commute — three nodal planes, not a cone. The two built remedies (entangling Peierls bond, three-axis covariant bond) each pay a measured price and neither produces a cone. The paper calls this "structural, not a finite-radius artefact", and I agree — but then the theory does not recover special relativity at its nodes, while its own Section 3 derives a Minkowski light cone and Theorem 5 measures front speed exactly 1. These are two different statements about the same lattice, and the manuscript's own framing ("the light cone and the causal radius are distinct statements") does not resolve the contradiction: a universe whose matter propagates at (0.100, 0.100, 0.142) Planck-cell units with anisotropy ~1 in every direction is observationally excluded at the 10⁻²² level by gamma-ray timing and by the isotropy of the cosmic-ray sky. The paper must either (a) elevate this to tension number eight, stated in the abstract, or (b) show a mechanism (renormalisation of the cone, a continuum limit, an averaging over frames of the 10-dimensional generator space) that restores isotropy. Neither is present. This is, in my judgement, the single deepest physics problem of the construction.

**M3 — The unification result of the abstract is conditional on an extension with six fitted inputs; the abstract omits this.**
The abstract states "the three measured couplings meet at M_X = 1.8×10¹⁶ GeV at exactly the root coupling α_X = 1/25, predicting τ_p = 1.8×10³⁶ yr". Section 13.1 shows that (i) under Standard-Model running the couplings do not meet (spread 47%; the SU(2) coupling overshoots 5.8×, SU(3) hits a Landau pole); (ii) the automaton's *own* Wilson mechanism forces M_X = 8.9×10¹⁴ GeV, excluded by Super-Kamiokande, and this is proven rigid (rank-1 Jacobian, split-invariance exact); (iii) the surviving solution at α_X = 1/25 exists only for a **flavour-dependent bond layer** — an *extension* of the automaton with six named inputs κ_s that are "a solution family, not a unique prediction", and the paper proves these six are not obtainable from any on-site algebra element, no class invariant, nor the most general off-diagonal bond operator (floor residual 30.4). The correct abstract statement is therefore: "unification is achieved *conditionally* on a named extension whose six couplings are fitted and not supplied by the root". As written, the abstract claims a derivation the body explicitly refuses to claim. The abstract of a paper must never be less honest than its Section 13.

**M4 — Coupling-convention inconsistency between Sections 9 and 13.**
Section 13.1 defines α₁₀ ≡ g² = 1/25 at g = 1/5 — explicitly flagged as "a convention". Section 9.2 compares the closed angle against "the 1/25 = 0.04 at which the desert solution meets" using α = g²/4π = π/64 = 0.0491. These are two different normalisations of the same rotation angle: g² = π²/16 = 0.617 at g = π/4, so the comparison 0.0491 vs 0.04 mixes α = g²/4π with α = g². Either the 1/25 meeting claim is in the g² convention (in which case the closed value is g² = π²/16, a factor 15.4 larger, not 1.226), or it is in the 4π-normalisation (in which case the implemented value is g²/4π = 0.0032, not 0.04). The "factor 1.54 worsened" and "2.41× above the data-optimal" statements inherit whichever convention is chosen. This must be fixed; it touches a headline falsifier (F2).

**M5 — "Zero continuous unfixed numbers" is contradicted by the paper's own continuous inputs.**
The closure axiom quantizes the *rotation angle*, and the paper counts this as the removal of the one continuous parameter. But the construction retains, by the paper's own statements: the Wilson angle θ = 1.731 at "the physical root" (a real number selecting which three of twelve nodes split and at which scales — nothing in the root fixes it); the junction amplitude r, explicitly "not quantized... a free junction amplitude inside a finite window" (Section 12.3); the ordering of the 45 plane-rotation factors, an input on which the order-8 closure depends (five shuffled orderings break it, deviation ~1); the split node set (O9); and the six fitted κ_s. The parameter audit counts only the on-site dynamics. The abstract's "zero continuous unfixed numbers" is therefore true in a narrow sense and misleading in the broad sense in which a reader will take it. It should be rephrased: "zero continuous unfixed numbers *in the fundamental update*; the phenomenological sector carries a continuous Wilson angle, a continuous junction amplitude, an ordering convention and six fitted couplings."

**M6 — The closure axiom is a new postulate, not a consequence, and it worsens every empirical contact.**
The demand that the fundamental update "be an element of the root" (an algebraic, not analytic, map) is a plausible ontological axiom but it is introduced to remove a parameter the data already disfavoured (g = 1/5 rejected against the best-fit 1/n², n = 7). Its price is measured and heavy: the hierarchy crossing moves from n* = 10 to n* = 208, the desert solution is destroyed (τ_p jumps to 2.8×10³⁹ yr, (1,2) threshold above M_X), the coupling moves 1.54× further from the one-loop target, and the doubler census becomes a nodal-line network on which the Nielsen–Ninomiya counting must be redone. A referee must ask: why is this axiom compelling? 't Hooft's determinism does not require it; unitarity does not require it. As presented it is a postulate whose sole empirical support is aesthetic. The paper's own framing ("the closure repairs the parameter criterion algebraically and worsens it phenomenologically") is honest, but the abstract's "zero continuous unfixed numbers" headline then advertises as a *victory* a move that the body reports as a net empirical loss, and the implemented g = 1/5 is retained as "the fallback theory". Two theories now coexist in one manuscript, sharing one algebra and disagreeing on the angle; the manuscript never decides between them, and the abstract of neither branch is stated separately. This is unacceptable for a journal reader: the paper must either (i) declare g = π/4 the theory and g = 1/5 the falsified fallback, with the abstract rewritten accordingly, or (ii) declare g = 1/5 the implemented theory and the closure axiom a proposed refinement, with both abstracts stated. As written, the reader cannot tell which theory is being claimed.

**M7 — The dark-matter prediction is a mass without a carrier, and the abstract's "predicting a 118 GeV dark-matter mass" is not a prediction of a particle.**
Section 12.2 proves that *no state of the implemented spectrum* can be the dark matter (0 survivors over 720 Wilson angles; the mirror sector's lightest states at 1.3×10⁸ and 2×10¹⁰ GeV; the required parity "does not exist inside the Clifford algebra of the root"), the computed relic abundance is not 0.120 for any coupling the root supplies (0.798, 0.0246, or 2.3×10⁻⁵–7.6×10⁻⁴ with the nodal-plane density D = 1485), and the minimal fix (a mirror-odd ladder scalar) is excluded by direct detection in its Z₂-preserving form and breaks the exact mirror Z₂ in its surviving form. What remains is a *mass coincidence*: one equation, v = 2.08 m_DM, correctly identified by the paper itself as "one equation, not two independent confirmations" (Fig. 18a). The abstract should read: "the ladder predicts a mass scale 118 GeV at which no state of the implemented spectrum can be stable, and the computed abundance fails for every supplied coupling." F7 is then not a falsifier of a prediction but a test of a coincidence.

**M8 — The neutrino result is a candidate, not a result, and the exponent is fitted.**
The bare ladder value 0.146 eV is excluded by Σm_ν < 0.12 eV — the paper says so — and the ⅓ candidate at 0.0487 eV requires (i) a junction structure built by the author (the cubic-symmetric zero-mean Wilson-strength junction), (ii) a democratic m_D "else gapless coloured doublers", and (iii) a value r = 1.91 selected by the atmospheric scale (or r = 2 "a junction-amplitude convention, not a tuning" — a statement the referee does not find established: a convention chosen to land on the measured value is a fit by another name). The exponent identity 2×10 − 3 = 17 is exact, but n_2step = 3 is proven *not* to be a link-word length, so 17 = 20 − 3 is a fitted integer that the Clifford algebra demonstrably does not supply (a genuinely nice negative result, stated honestly). The abstract's "n = 17 the neutrino mass" must carry the word "fitted".

**M9 — Baryogenesis fails, as the paper states, and the abstract correctly reports it; but then Sakharov-3 is in permanent conflict.**
The computed theorem (all CP-odd invariants of the step vanish; none of the three Sakharov conditions met; stable under the Wilson extension) is a strong negative result. The Peierls-flux route supplies a CP phase but not the out-of-equilibrium condition at the Yukawa fixed by the atmospheric mass (Γ_N/H = 47 at y = 0.40, vs the required y < 0.058). The incompatibility of m_ν ≃ 0.05 eV with out-of-equilibrium decay at M_N = 9.6×10¹³ GeV is a real theorem of this framework and belongs in the abstract's tension list (it is currently only in Section 13.3 and O5).

**M10 — "Reproduces general relativity to machine precision" (abstract) contradicts Table 3.**
The residuals of Table 3 are finite-difference truncation errors: 6.5×10⁻⁷ (de Sitter Riemann), 5.5×10⁻⁵ relative (Kretschmann), 10⁻⁵ (2-sphere), stated by the paper itself as O(h²). Section 1.4 correctly defines a machine-precision identity as ≤ 10⁻¹⁴. The abstract's "machine precision" therefore violates the paper's own tolerance discipline within two pages of stating it. The honest sentence is: "reproduces GR on de Sitter, Schwarzschild and 3+1 evolutions to O(h²), with the Kretschmann invariant agreeing with the analytic form to 5.5×10⁻⁵ relative."

**M11 — The measurement postulate and the Born rule are inputs; the determinism claim must be softened correspondingly.**
The paper honestly declares the projection as an input (O9). But then the abstract's "'t Hooft's cellular-automaton determinism exact at the fundamental update, with randomness entering only through the measurement projection" is a restatement of 't Hooft's interpretation, not a result of this construction: any unitary operator admits a deterministic re-reading of the evolution of states; what 't Hooft's programme requires — and what this manuscript does not supply — is a derivation of the Born rule from the automaton. The manuscript should state in the abstract: "the Born rule is not derived (O9)".

**M12 — Spin-statistics is postulated, and the CPT claim is properly weakened but only deep in Section 7.2.**
The Z₂ grading identifying even elements with bosons is "a postulate of the interpretation". No spin-statistics theorem is derived, no antiunitary T is constructed on the many-cell Hilbert space. The abstract's "exact matter/antimatter spectral mirror" is accurate; but the word "CPT" appears in Section 7.2 in a context the paper itself flags as overstating. Fine — but then Section 6's "the even sector is identified with bosons and the odd sector with fermions" is doing load-bearing work in the ontology and must be flagged in the abstract's list of inputs alongside O9.

**M13 — Format and scope.**
The manuscript is ≈ 29,000 words with 25 figures and 52 lines exceeding 1000 characters (table rows up to 2525 characters, e.g. Table 7 rows O2, O4, O8, O9). A journal referee cannot referee a 9-table audit in which individual table rows exceed the length of an abstract, and a reader cannot extract the physics from tables written as prose. PRD accepts long papers, but not papers in which the audit tables carry the narrative. The manuscript must be restructured: (i) a paper on the algebra, the step, unitarity, causality and the spectral mirror (Sections 2–9), with the doubling and its Wilson remedy; (ii) a paper on the gauged bond layer, the area law and gravity (Sections 10–11); (iii) a paper on the phenomenology with the provenance audit (Sections 12–17). Each is publishable in roughly its present content; together they are not.

**M14 — Reproducibility, and what a referee can and cannot verify.**
The claim that all numbers are produced by public scripts is the right practice, and my spot-checks all pass. However, the manuscript presents no independent verification channel (no co-verification of the 32-dimensional representation by a second party, no cross-check of the closed-loop F_n computation against an analytic limit where one exists). The symbolic cross-check of Appendix A (a second construction of Cl(10,0), agreement as identities) is good practice and I commend it; the same discipline should be extended to the traced-loop computation, where an analytic limit (the 3^n-word enumeration) does exist and is quoted — good — but where the *physical* claim ("the phase supplies the factor, not the topology") rests on a Wilson phase whose value is an input.

---

## Minor points and errors

1. **Reference 5**: "T. D'Ariano" — the author of the QCA free-field papers is **A. D'Ariano** (with P. Perinotti). Correct the initial.
2. **Reference 47**: the author list "Y. Hisano, M. Kobayashi, K. Muramatsu, and N. Nagata" does not match the standard citation for the two-loop proton-decay renormalisation factors (Phys. Lett. B **724**, 283 (2013)); verify the authors against the DOI.
3. **Reference 13**: "C. I. Castro" — verify initials; the Clifford-unification review author is Carlos Castro (sometimes Castro Perelman), and the two 2014/2015 citations should be checked for the exact titles.
4. **Abstract length**: the abstract is three paragraphs and ~450 words with a keyword list; journals cap it. Split into a 200-word abstract and the body of the summary.
5. **Notation**: three symbols for the same angle/coupling — the rotation angle g, the coupling g₁₀ = 1/5 of Section 10.4, and α₁₀ ≡ g². Introduce them in one place, in one convention (see M4).
6. **Section 15, "The nine-to-one observation"**: this paragraph states it is "a bookkeeping observation... not a derivation" — then why is it in the parameter audit at all? Either make it a remark or delete it.
7. **Fig. 2 caption**: the hydrogen atom is an input, correctly stated twice in two consecutive sentences; compress.
8. **Section 12.2, line-length**: the O7 row of Table 7 and the paragraph of Section 12.2 exceed 2000 characters; a journal cannot typeset them.
9. **English/typos**: e.g. Section 7.3 "the front *saturates* the cone" (underscore formatting leaked from a draft); "front speed 1.0" vs "exactly 1" — pick one form.
10. **Table 5, row "16-sector block-Brillouin spectrum"**: this is a physics paragraph in a table; move to text.
11. **Section 13.2**: the sentence "running which, as just shown, the theory does not currently perform consistently" is a good honest statement; but then the claim "sin²θ_W = 3/8... close to the measured 0.23129" in the abstract's neighbourhood invites misreading. State in the abstract that no electroweak mixing prediction is claimed.
12. **Falsifier F4**: correctly states that the Section 10.4 quintet channel is *not* claimed to be the graviton; but the abstract says "the lattice gravitational sector reproduces general relativity" — the reader will conflate the two. Add "(the g_μν field is an input; see Section 11.1)" to the abstract.
13. **Data availability**: the GitHub link is given; the referee could not execute the pipeline within the review and relied on internal-consistency spot-checks. State this limitation explicitly in the paper ("numbers verified internally; pipeline reproducible from the repository").

---

## Questions to the authors (answer before resubmission)

1. In which normalisation is α = 1/25 (g² or g²/4π)? Reconcile Section 9.2 with Section 13.1 and recompute the "factor 1.226" and "1.54×" statements (M4).
2. Is θ = 1.731 fixed by anything in the root? If not, is the "physical root" a fit? Same question for r = 1.91 vs the "convention" r = 2 (M5, M8).
3. Which theory is being claimed: g = 1/5 (implemented, Lorentz-failing at the nodes, α = 1/25 meeting with a fitted desert spectrum) or g = π/4 (closed, order 8, nodal lines, α = π/64, no desert)? Give one abstract per branch or pick one (M6).
4. Does any mechanism (frame-averaging over the 10 generators, a continuum limit, a deformation) restore cone isotropy at any gapless point? If not, the theory predicts observable direction-dependent propagation — where does it sit in the tension list (M2)?
5. What is the exact status of the 16-sector gapless orbit of 24 points (Section 8.3): are these the matter degrees of freedom of the low-energy universe? If yes, the one-generation statement of the abstract must be qualified: the lattice carries 12 doubler families plus an orbit, not one generation.
6. The 5.4×10⁴⁰ string tension of the gauged layer: since no light vector sector is supplied (gap 64 M_P) and Λ_QCD is not on the ladder, is the layer intended as the QCD confining sector or not? Section 10.2 says confinement is "absent from the implemented step"; then what confines quarks in this theory?
7. How does the automaton produce the *three* light families given that the implemented spectrum has twelve gapless nodes organised in three families of four — is the observed family number identified with the axis-family label, and what observation would distinguish this from three accidental families (O3)?

---

## Assessment against the journal's criteria

| criterion | verdict |
|---|---|
| originality | **high** — the coupling-invariance of the algebraic invariants, the closure axiom and the provenance audit are, to this referee's knowledge, new in this combination |
| validity of the results claimed | **partial** — the algebraic and dynamical results (unitarity, causality, anomalies, 5/3, spectral mirror) are sound; the unification, dark-matter, neutrino and baryogenesis claims are conditional, fitted or negative, and are stated with varying degrees of candour between abstract and body |
| clarity | **insufficient** in present form — audit tables with rows exceeding 2000 characters, three competing theories in one manuscript, one inconsistent coupling convention |
| suitability for this journal | **no** — the manuscript is a research programme; it should be a series of papers |
| reproducibility | **credibly claimed, not independently verified by this referee; internal spot-checks all pass |

---

## Disposition

**Reject**, with an explicit encouragement: the algebraic/dynamical core (Sections 2–9) is, in this referee's judgement, publishable after a focused rewrite that (a) resolves M4 (coupling convention), (b) states M2 (Lorentz failure) as a named tension in the abstract, (c) separates the closed-branch and implemented-branch phenomenologies, and (d) restores journal formatting. The phenomenology (Sections 12–17) should be a second paper whose abstract contains the seven tensions verbatim and none of whose headline claims outruns its own Section 15 audit. The gravity and gauged-layer material (Sections 10–11) is a third paper; its "one object" claim must be removed from the title.

The manuscript's most valuable feature — the refusal to hide a single failure — is rare, and the referee wishes to state that plainly: a work that reports "the theory fails the parameter criterion" against its own seven-criteria table, displays its twelve gapless doublers, and computes its own relic abundance only to report that it is wrong, is closer in spirit to good physics than most unification preprints. But good physics also requires that the abstract never say more than Section 15 does, and in the present version it does.

*Referee report filed as* `review/glmreview.md`.
