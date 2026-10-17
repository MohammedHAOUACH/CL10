# The alternative-cell and signature exploration

**Programme:** `alt_cell` (merged into `master`).
**Experiments:** `experiments/alt_cell/V149`-`V160` (12, tier 4).
**Status:** exploration only. No publication claim changes; the 115 experiments
of the publication programme are untouched and their verdicts are unchanged.

---

## 1. The question

The master construction puts matter on a **2x2x2 block** (8 sub-sites),
advances it through **six bond stages** (three axes, two parities), and
takes the internal algebra to be **Cl(10,0)**. The published node
census (V63/V115) finds 12 gapless half-zone nodes at the M momenta
with anisotropic cones (A = 0.998). Two independent suspicions recur
in any reading of the construction:

1. **Signature.** Is Cl(10,0) forced, or do other Cl(p,q), p+q = 10,
   give the same physics with different node geometry?
2. **Cell.** Is the node doubling an artifact of the 2x2x2 block and
   the six-stage cycle - would a different cell (3x3x3, an FCC/diamond
   1-factorization) or a different stage order remove it?

This branch answers both with a pre-registered census. The gates were
fixed **before** the candidate scans ran:

- **G1 unitarity:** the 45-plane product is exactly unitary at the
  coupling (order <= 64 at one of the pre-registered angles
  {pi/2, pi/4, g_10 = 1/5}).
- **G2 mirror:** the step commutes with the pseudoscalar w
  (16/16-bar leakage < 1e-12).
- **G3 frame:** an exact (1+3) frame exists (four +I generators,
  residual 0).
- **G4 charges:** the SM charge table is reproduced (anomaly-free,
  5/3 hypercharge).
- **G5 geometry:** gapless nodes survive, one O_h orbit, cone
  isotropy A < 1e-6, chirality sum 0, front speed exactly 1 with zero
  outside-light-cone amplitude.

Every candidate is checked on **two routes**: the exact Bloch matrix
of the cell and the real-space step on a small torus (the V63/V154
method). A candidate must pass both.

## 2. The signature census (V149-V152, V159)

All eleven signatures Cl(p,q), p+q = 10, were run through the full
step machinery.

| q | signature | Hermitian planes | unitary at g_10 | (1+3) frame | SM charges | nodes (27 pts) | best A |
|---|---|---|---|---|---|---|---|
| 0 | Cl(10,0) | 45/45 | yes (9.9e-32) | yes (res 0) | exact | 12 (M) | 0.998 |
| 1 | Cl(9,1) | 36 | no (res 1.0 at pi/4) | yes | broken | - | - |
| 2 | Cl(8,2) | 29 | no | yes | multiset ok, multiplets broken | - | - |
| 3 | Cl(7,3) | 24 | no | yes | broken | - | - |
| 4 | Cl(6,4) | 21 | no | yes | permuted (anomaly-free) | - | - |
| 5 | Cl(5,5) | 20 | no | yes | broken | - | - |
| 6 | Cl(4,6) | 21 | no | yes | multiset ok, multiplets broken | - | - |
| 7 | Cl(3,7) | 24 | no | **no** (res 4) | broken | - | - |
| 8 | Cl(2,8) | 29 | no | no | multiset ok, multiplets broken | - | - |
| 9 | Cl(1,9) | 36 | no | no | broken | - | - |
| 10 | Cl(0,10) | 45/45 | yes (9.9e-32) | **no** (res 4) | conjugated (anomaly-free) | 0 (gapped) | - |
(V149: plane hermiticity, unitarity residual, closure order; V150:
w^2, 16/16 split, frame realizability; V151: charge-table match;
V152: full node/cone/chirality/torus census.)

**Findings.**

- **Only q = 0 and q = 10 keep the 45-plane product exactly unitary**
  (V149). Mixed signatures make some planes anti-Hermitian (T^2 = -I);
  the plane product then has unitarity residual O(1) at every
  pre-registered angle. The GUT step machinery exists only at the two
  Euclidean/pseudo-Euclidean ends.
- **The 16/16-bar split survives every signature** (w^2 = -1 for q
  even; V150), but an exact (1+3) frame needs four +I generators:
  q <= 6. q = 10 has no frame (residual 4).
- **Only q = 0 reproduces the SM charge table** (V151): q = 4 permutes
  the fields, q = 10 conjugates all hypercharges, odd q break the
  multiplets (complex hypercharges).
- **q = 10 gaps the sampled M nodes** (M gap 0.0310 at all 27
  high-symmetry momenta, V152) - the signature *does* move the node
  spectrum. But q = 10 has no Minkowski frame and conjugated charges:
  it is not our world. **No signature with a (1+3) frame removes the
  12 nodes or isotropizes a cone** (best A = 0.998).
- **CP (V159):** q = 0 is real (Majorana) type - the V117 no-phase
  theorem. q = 10 is **quaternionic** (p - q = -10 = 6 mod 8) and its
  step is unitary: a CP phase is *not* forbidden there. It is the one
  interesting corner of the family - and it has no frame and
  conjugated charges. The CP window and the Minkowski frame are
  mutually exclusive in the p+q = 10 family: O5 stays closed for
  every viable candidate.

**Verdict (signature):** Cl(10,0) is the unique signature with
unitary GUT step + (1+3) frame + SM charges + anomaly freedom. The
signature is not a free knob; the node doubling is not caused by it
in any usable way.

## 3. The cell and cycle census (V153-V156)

| census | space | result |
|---|---|---|
| V153 stage cycles | all 720 orderings of the six stages | **720/720 gapless** at the 27 high-symmetry momenta, 720 single-orbit node sets, 0 beat the baseline A = 0.998, ZERO reach A < 1e-6: the order is not the cause |
| V154 blocks | {2,3}^3, implemented + canonical cycle | 8 shapes: front speed exactly 1, zero outside-light-cone amplitude everywhere; (2,2,2) keeps 12 M nodes under both cycles; every 3-division block gapped at the 27 sampled momenta (sampling statement); no block reaches A < 1e-6 |
| V155 FCC/diamond | 524 1-factorizations (24 diamond, 500 FCC cap) | **zero gapless** at the FCC high-symmetry points: the graphene/doubling mechanism does not transfer |
| V156 node forcing | 1254 candidates from V152-V155 | 722 gapless; representatives recomputed: 4 single-orbit (split forced: only empty or full - all restatements of the root cell), 0 several-orbit; V153 independently: 720/720 cycles single-orbit; none combines a forced split with an isotropic cone |

**Findings.**

- The causal structure is cell-independent: **every** block shape
  keeps front speed exactly 1 with zero outside-light-cone amplitude
  (V154). The cell size is not the cause of anything causal.
- Blocks with a 3-division are gapped *at the sampled momenta* - but
  the (2,2,2) anchor keeps its 12 M nodes under both the implemented
  and the canonical stage cycle, and the stagger only moves nodes off
  the sampled set or gaps them: a sampling statement, not a removal
  theorem.
- The FCC/diamond route (the "graphene would fix it" intuition) is
  dead in the opposite direction: those lattices are **gapped** at
  their high-symmetry points with this step - the doubling mechanism
  does not transfer.
- **The stage order is irrelevant:** every one of the 720 orderings
  of the six stages is gapless with a single-orbit node set, and none
  improves the cone anisotropy at all (V153). The doubling is carried
  by the stage *set* and the cell, not by their order.
- **The split stays an input:** across all 1254 census candidates the
  single-orbit (forced-split) sets are all restatements of the root
  cell, and none pairs a forced split with an isotropic cone (V156):
  O3/O9 stay inputs for the whole family.

## 4. The hard tier and unification (V157-V158)

- **V157 (gate table):** zero candidates survive the census gates.
  All 722 algebraically eligible near-misses - the 720 stage cycles
  and the root's two forms - pass all five hard gates (unitarity,
  front speed 1, anomalies, 5/3, closure at a pre-registered angle):
  the hard tier is jointly satisfiable but does **not** discriminate.
  The discriminating gate is the isotropy gate (A = 0.998 for every
  candidate - the root's own published W9 residual). Every
  elimination is recorded; no candidate is repaired.
- **V158 (unification-forcing):** the node content of every gapless
  candidate is a **complete 16 + 16-bar per node** - a proportional
  shift of the beta-function coefficients. V64's proportional-shift
  no-go therefore kills the gauge coupling meeting for every
  candidate: 0 meetings found, none within a factor 2 of 2e16 GeV.
  **No cell or signature forces the meeting; F5 stays open for the
  whole family.**

## 5. The dark-matter carrier (V160)

The lightest excitation of every gapless candidate is a **gapless
node** (eigenphase 2.0e-4), never an isolated neutral massive mode;
the mirror Z2 pairs every mode with its 16-bar partner. 0 carriers
found. The 118 GeV ladder mass (V122) is not reproducible as a
geometry artifact of the cell: **W7 is not a cell artifact.**

## 6. What this branch settles

1. **Cl(10,0) is forced** by the conjunction of the four algebraic
   gates (unitary 45-step, (1+3) frame, SM charges, anomalies). The
   only unitary alternative, Cl(0,10), has no frame and conjugated
   charges.
2. **The node doubling is not a cell artifact.** It survives the
   canonical stage cycle and all 720 stage orderings (720/720
   gapless, none isotropic), every block shape's causal structure is
   clean, the FCC/diamond alternative is gapped. The doubling is a
   property of the six-stage fermion-doubling structure itself (the
   Wilson term of V65 is the published answer).
3. **No candidate forces gauge unification** (V158): the node content
   is always a complete generation pair, so V64's no-go applies
   uniformly. F5/O1/O9 stay open exactly as published.
4. **The CP window (O5) and the Minkowski frame are mutually
   exclusive** across the whole signature family (V159).
5. **W7 (the 118 GeV carrier) is not a geometry artifact** (V160).

Every number above traces to a verdict JSON under
`experiments/alt_cell/*/results/`.
