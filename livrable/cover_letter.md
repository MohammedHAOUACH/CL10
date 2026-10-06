# Cover letter — submission

**Manuscript:** *Causal Clifford Theory (CL10): the Standard Model and gravity
from a single Clifford root on a Planck-lattice quantum cellular automaton*
(`livrable/publication.pdf`, 26 figures, 32 references).

**Primary target:** *Physical Review D* (lattice formulation + GUT
phenomenology + explicit falsifiers).
**Fallbacks, in order:** *European Physical Journal C*, *Journal of High Energy
Physics*, *Foundations of Physics* (if the referees judge it a
quantum-foundations proposal rather than a phenomenological model).

**arXiv:** primary `hep-ph`; cross-list `hep-th`, `gr-qc`, `quant-ph`.

---

Dear Editor,

We submit the enclosed manuscript, which reports a construction rather than a
model-building claim: a unitary quantum cellular automaton on a Planck lattice
whose on-site internal space is the 32-dimensional module of
$\mathrm{Cl}(10,0)$, and the statement of exactly what that one algebraic
object does and does not fix.

Demonstrated in the paper: the $(1{+}3)$ Minkowski subalgebra inside
$\mathrm{Cl}(10,0)$; the 16-component Weyl module whose charges are one
Standard Model generation with all eight gauge and mixed anomalies cancelling;
the 45-dimensional $\mathrm{Spin}(10)$ adjoint containing the twelve gauge
bosons; the $5/3$ hypercharge normalisation; an exactly unitary, strictly
causal update ($\det U = 1$) with an exact matter–antimatter spectral mirror;
and a lattice theory of a massless spin-2 field whose discrete
Einstein–Hilbert action reproduces general relativity to machine precision
(two independent numerical routes, including a 3+1 BSSN/Z4 evolution).

The paper's methodological content is a **provenance discipline**: every
number is labelled *algebraic*, *conventional*, *fitted* or *open*, and the
theory is judged against seven criteria it states for itself. It fails its own
parameter criterion and says so. Five tensions with observation are reported
rather than repaired: a $U(1)$ Landau pole at $\sim 280$ GeV when
$\alpha_{10} = 1/25$ is run down with the Standard Model $\beta$-functions; a
seesaw neutrino mass $0.146$ eV above the oscillation scale and the
$\sum m_\nu < 0.12$ eV bound; a relic abundance that is normalised rather than
computed; heavy-neutrino decays that are *in* equilibrium
($\Gamma/H \simeq 4.6$), so the third Sakharov condition fails; and a
computed Brillouin-zone spectrum that is gapless at twelve doubler momenta.
Nine open residuals (W1–W9) and seven falsifiers (F1–F7) are listed with the
size of each gap quantified.

We believe this honesty is the paper's contribution: a discrete unification
is easy to overstate, and here the boundary between what follows from the
algebra and what does not is explicit, machine-checkable and falsifiable.

**Reproducibility.** The implementation, the 64-experiment verification suite
and the figure generators are public (MIT licence) at
https://github.com/MohammedHAOUACH/CL10; every figure is produced by that
code. No measured value in the paper is produced by the code — the comparison
values are the published data of PDG 2024, Planck 2018, Super-Kamiokande,
KATRIN and NuFIT.

**Declarations.** No competing interests. No funding. Single author.

**Suggested referees** (each active in a literature the paper is judged
against; none has a co-authorship or institutional link with the author):

1. **Pawel Arrighi** (Luminy, Marseille) — quantum cellular automata, the
   index-theoretic classification of QCA steps.
2. **Renate Loll** (Radboud, Nijmegen) — discrete quantum gravity / causal
   dynamical triangulations, for the lattice-gravity claims.
3. **H. B. Nielsen** (Niels Bohr Institute) — lattice fermions and the
   doubling theorem the paper's W9 residual confronts.
4. **M. S. Turner** (Chicago / Fermilab) — baryogenesis and the
   sphaleron/leptogenesis literature the paper's W5 residual cites.
5. **G. 't Hooft** (Utrecht) — deterministic cellular-automaton
   interpretations of quantum mechanics; the paper cites his construction
   heavily, which the editor may prefer to avoid.

Sincerely,
Mohammed Haouach
haouach.mohammed@gmail.com

---

## Before submitting (author-only steps)

1. **ORCID.** Put the number in `livrable/orcid.txt`; the build appends it to
   the affiliation. Without the file the compiled paper has no ORCID.
2. **Zenodo.** Link the GitHub repository to Zenodo (`zenodo.org/github`),
   publish a release, and add the concept DOI to the code-availability
   paragraph of `doc/publication.md` so the repository citation is versioned.
3. **Repository.** The paper asserts https://github.com/MohammedHAOUACH/CL10;
   make sure that URL is public before submission.
4. **Uncited DOIs.** Refs 13, 14, 15, 19, 23, 24, 29, 31, 32 carry no DOI in
   the compiled bibliography (not resolvable through Crossref from this
   machine); confirm them by hand.
