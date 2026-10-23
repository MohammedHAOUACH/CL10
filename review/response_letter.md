# Response to the referee report — resubmission as a three-paper series

**Manuscript reviewed:** "The Universe as a Quantum Cellular Automaton: the Standard Model and Gravity from a Single Clifford Algebra" (`doc/final_publication.md`)
**Referee report:** `review/glmreview.md`
**Point-by-point verification:** `review/glmreview_verification.md` (every claim of the report checked against the manuscript, the repository code, and new computations)

---

Dear Editor,

The referee recommended rejection with an explicit disposition: split the work into three papers, fix the coupling convention, state the Lorentz failure and the conditionality of the unification claim in the abstracts, and restore journal formatting. We have followed that disposition exactly. The single 29,000-word manuscript is withdrawn and replaced by three companion papers built mechanically from the corrected master (verbatim body ranges, renumbered sections, figures, tables and references, explicit cross-paper pointers):

- **Paper I — "Causal Clifford Theory (CL10): Minkowski kinematics, one anomaly-free generation, and exactly unitary causal dynamics from a single Clifford algebra"** (`doc/paper1_core.md`): the algebraic and dynamical core, the doubling and its Wilson gapping, the Lorentz structure at the nodes with its measured anisotropy, and the closure axiom with its measured costs. This is the paper the referee judged "publishable after a focused rewrite".
- **Paper II — "The gauged bond layer of Causal Clifford Theory: exact Gauss constraint, bounded area law, string tension, and the lattice gravitational sector"** (`doc/paper2_gauge_gravity.md`): the gauged layer, the area law, and the gravitational sector — whose title makes no "single algebra" claim, as required.
- **Paper III — "The loop-factor ladder and the phenomenology of Causal Clifford Theory: four scales, coupling unification, dark matter, neutrinos, and the provenance audit"** (`doc/paper3_phenomenology.md`): the ladder, the couplings, the cosmology, the provenance audit and the falsifiers; its abstract carries the nine tensions verbatim and outruns no statement of its own audit section.

## The four defects the report identified, and their resolution

1. **The coupling convention (M4).** The referee is right that Sections 9 and 13 used two normalisations of the same angle. We re-ran the desert inverse problem in the standard normalisation $\alpha = g^2/4\pi$: the meeting value is $\alpha_X = 0.04005$ (standard), the implemented angle gives $\alpha = 0.00318$ (a factor 12.6 below), the closed angle gives $\alpha = \pi/64 = 0.0491$ (a factor 1.226 above). The old "exact match at exactly the root coupling" was the $\alpha_{10} \equiv g^2$ convention — a coincidence of conventions, not a physical match, and we now say so. One convention is stated once (Section 1.4 of Paper I) and used everywhere; the convention-dependent "factor 1.54 / 2.41" statements are deleted; F2 states the closed-branch prediction $\alpha_X = \pi/64$ in the same normalisation.
2. **The Lorentz failure (M2).** Now named in the abstract of Paper I with its numbers (anisotropy 0.748–0.758; 0 of 2726 gapless points of the Wilson family isotropic; three nodal planes at the corner multiplet). The referee's Question 4 asked whether frame averaging over the ten generators could restore a cone. We ran the experiment (V142 in the repository): it cannot — averaging over four random SO(10) frames leaves the gap, the principal velocities and the anisotropy unchanged to $3\times10^{-15}$, because the momentum dependence lives entirely in the sub-site factor of the bond layer. The failure is a named tension (abstract, O9, falsifier C3), not a softened footnote.
3. **The abstract/body gap (M1, M3, M5, M7–M12).** Each abstract now says no more than its audit section: gravity adds $g_{\mu\nu}$ and $G=1$ by hand and reproduces GR to $O(h^2)$, not to machine precision; unification holds only for a desert spectrum with two fitted thresholds and a flavour-dependent bond layer with six fitted couplings; "zero continuous unfixed numbers" is scoped to the on-site update, with the Wilson strength, the junction amplitude, the node set and the $\kappa_s$ counted as inputs; the dark-matter statement is a mass scale with no carrier in the implemented spectrum; the Born rule, the spin–statistics identification and the split-node set are declared inputs.
4. **Format (M13).** The series realises the referee's split exactly. No line in any paper exceeds 2000 characters; the four table rows the report measured at up to 2525 characters are all below 2000; the abstracts are single paragraphs.

## Answers to the seven questions

The full answers are in `review/glmreview_verification.md`; in brief: (1) $\alpha = 1/25$ is standard-normalised, and the "exact match" was a convention coincidence — recomputed consistently above. (2) $\theta = 1.731$ and $r = 1.91$ are fits, not numbers of the root; both are now counted as inputs. (3) Both branches ($g = 1/5$ implemented, $g = \pi/4$ closed) are theories of the same root; the series separates them and no abstract mixes them. (4) Frame averaging is blind — new experiment V142, numbers above. (5) The 24-point orbit is a gapless doubler residual of the Wilson mechanism, not the matter of the low-energy universe; the one-generation statement is about the internal module and every abstract displays the doubling beside it. (6) The gauged layer is not the QCD confining sector (tension $5.4\times10^{40}$, no light vector sector, $\Lambda_{\rm QCD}$ off the ladder); the theory does not currently confine quarks, and this is stated as a named tension. (7) The family count 3 is a symmetry-breaking input; it is already distinguished from the twelve-node orbit reading by proton decay ($M_X(12) = 2.9\times10^{15}$ GeV excluded), but not from three accidental families, and the paper says so.

## Two corrections in the report that do not survive checking

Minor 1: the Frontiers 2016 paper is by A. Bisio, **G. M.** D'Ariano and P. Perinotti — the master's initial was right; we completed the author list rather than introduce the report's "A. D'Ariano". Minor 2: the Hisano reference is verified against arXiv:1302.2194 and Phys. Lett. B 724, 283 (2013): J. Hisano, **D.** Kobayashi, **Y.** Muramatsu, N. Nagata — the master is correct as written.

## What we thank the referee for

The report's judgement that the provenance discipline is the manuscript's distinctive virtue, and its spot-checks (all of which reproduce), encouraged us to make the discipline load-bearing rather than decorative: every headline number now carries its class in the same normalisation, every failure named in an abstract carries its number, and the repository carries the new V142 experiment answering Question 4 by computation rather than assertion.

We believe the series meets the journal's standard for honesty and precision, and we are ready to answer further questions.

Sincerely,
Mohammed Haouach
