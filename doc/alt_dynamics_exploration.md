# The alternative-dynamics exploration

**Programme:** `alt_dynamics` (merged into `master`).
**Experiments:** `experiments/alt_dynamics/V161`-`V170` (10, tier 4).
**Status:** exploration only. No publication claim changes; the 115 experiments
of the publication programme are untouched and their verdicts are unchanged.

---

## 1. The question

The master construction fixes the **motion** of the 2x2x2 cell: one
time-independent step, six bond stages at the fixed angle $g_{10}=\pi/12$,
the internal layer the ordered 45-plane Wilson product, the same at every
stage. The published node census (V63/V115) finds 12 gapless half-zone
nodes with anisotropic cones ($A=0.998$). The companion branch
([`alt-cell-exploration`](alt_cell_exploration.md), V149-V160) varied the **cell** and the
**signature**; neither removed the doubling. This branch varies the third
free axis: the **step rule itself**. Five suspicions recur in any reading
of the construction:

1. **Per-stage angles.** Is the single angle $g_{10}$ forced, or do
   per-stage angles (e.g. $\pi/4$ on three stages) isotropize the cones?
2. **Internal factorization.** Is the ordered 45-plane Wilson product
   forced, or does another factorization of the same internal layer move
   the nodes off the half-zone?
3. **Mass stage.** Could one extra internal (mass) stage per cycle open
   the node gaps without breaking the mirror parity?
4. **Bond matching.** Is the axis/parity 1-factorization of the cubic
   bond set forced, or does another matching (graphene-like) remove the
   doubling?
5. **Time dependence.** Is the time-independent step forced - does a
   Floquet splitting or a quantized Wilson layer change the node
   structure?

This branch answers all five with a pre-registered census. The gates were
fixed **before** the candidate scans ran, and are the same gates as the
alt_cell branch:

- **G1 unitarity:** the implemented step is exactly unitary
  ($\|UU^\dagger-I\| < 10^{-12}$).
- **G2 front speed:** the light-cone front moves exactly one Planck cell
  per stage (measured by the real-space support growth of V154's
  `front_speed`).
- **G3 anomalies:** the charge table reproduces the SM anomaly
  cancellation.
- **G4 five-thirds:** the $\alpha_3$ channel closes at the root angle.
- **G5 closure:** the internal layer is the Wilson product (or a
  conjugate/Cartan-preserving restatement of it).

The candidate families are enumerated in `cl10/altcell.py` (the
`bloch_stages`, `step_stages`, `front_speed_stages` helpers added on this
branch); every candidate is a full step (bond stages + internal layer)
diagonalized on the 16-sector Bloch at the 27 high-symmetry momenta of
V63. The node/isotropy measurement is identical to V153-V157: a node is a
high-symmetry point with spectral gap $<10^{-3}$; $A$ is the max/min
eigenvalue ratio of the cone metric $M=\partial^2\omega/\partial q^2$ at
the node; a cone is *true* only if $M$ is positive definite and
$A<10^{-6}$.

---

## 2. The motion census (V161-V166)

| family | experiment | candidates | gapless | nodes | best $A$ | isotropic | closures |
|---|---|---|---|---|---|---|---|
| per-stage angles $\{\pi/2,\pi/4\}^6$ | V161 | 729 | 9 | 8-12 | 0.949 | 0 | 0 |
| internal factorizations | V162 | 12 | 2 | 12-13 | 0.998 | 0 | 0 |
| mass-stage insertions | V163 | 90 | 0 | - | - | 0 | 0 |
| cubic 1-factorizations | V164 | 500 (cap) | 487 | 2-12 | 0.998 | 0 | 0 |
| Floquet splittings | V165 | 3 | 0 | - | - | 0 | 0 |
| quantized Wilson layers | V166 | 8 | 1 | 8 | 0.993 | 0 | 0 |

**V161 - per-stage angles.** All $3^6=729$ assignments of
$g\in\{\pi/2,\pi/4\}$ to the six stages (the anchor $g_{10}$ row added).
720 are fully gapped at the 27 high-symmetry points. The 9 gapless rows
are the anchor (12 M nodes, $A=0.998$) and 8 rows with three $\pi/4$
stages and three $\pi/2$ stages (8 nodes, $A=0.949$): the $\pi/2$ stages
are pure lattice translations (V134), so the gapless rows are the
$\pi/4$-support of the root motion, not new physics. Zero isotropic, zero
closures. The two route orderings agree to $1.2\times10^{-14}$.

**V162 - internal factorizations.** 12 factorizations of the internal
layer: the Wilson product, its reverse, two Cartan-subalgebra
reorderings, and eight non-Cartan regroupings. Only the two
Cartan-preserving ones stay gapless - the baseline (12 nodes, $A=0.998$)
and `cartan5` (13 nodes, $M$ gap exactly 0, $A=1.0$, two orbits). Every
non-Cartan factorization gaps the entire zone. The nodes are
**Cartan-protected**: the gapless set is carried by the Cartan
subalgebra of the Wilson product, not by its ordering. The mirror
commutator $\|[U,w]\|_{\max}=7.8\times10^{-16}$ throughout.

**V163 - mass stages.** 90 insertions of one extra internal (mass) stage
per cycle: all 45 planes $\times$ $\{g_{10},\pi/4\}$ plus 0-parameter
representatives. 81 rows are fully gapped at the high-symmetry points;
the 9 that retain tolerance points have $M$ gaps $\ge5\times10^{-4}$ with
no positive-definite $M$ (the two $A<10^{-6}$ rows are massive zone
centers, cone growth $\sim0.9$, not cones). The mirror Z2 is never
broken: leakage 0 for every candidate. A mass stage gaps the nodes
without touching the mirror parity - but it never produces a true cone.

**V164 - bond matchings.** 500 cubic 1-factorizations (cap; 5 skipped by
the cap), the canonical axis/parity factorization included. 487 are
gapless somewhere; 31 have a single node orbit (all 12 nodes, best
$A=0.998$). Zero isotropic, zero closures. The node doubling is not an
artifact of the axis/parity matching: generic matchings keep it, and none
isotropizes it.

**V165 - Floquet.** Three time-splittings of the root step (AA, AB, BA:
internal/bond interleavings over a doubled period). All three open a gap
at every high-symmetry point ($2\times10^{-3}$ to $1.8\times10^{-2}$);
only the time-independent anchor stays gapless (12 nodes, $2\times$ the
anchor gap). Floquet time-dependence destroys the gaplessness rather
than protecting it.

**V166 - quantized Wilson layers.** The 8 quantized Wilson variants
built from the published coefficients ($X=(1,1,1,1,1)$,
$Y=(-\tfrac16,\tfrac16,\tfrac16,\tfrac14,\tfrac14)$, $w\in\{\pi/2,\pi/4\}$,
$\pm$ branches). 7/8 are fully gapped; the one gapless variant
($w=\pi/2$, $X_-$) keeps 8 nodes with $A=0.993$. The left/right shell
ratios of the gapped variants are $\sim1.0$ (quasi-local): the quantized
layers do not restore the half-zone nodes. The published Wilson
coefficients are not a free choice that happens to give nodes - the
node-carrying layer is the Wilson product itself.

---

## 3. The hard tier on motion (V167)

The five hard gates of V157, re-applied to the **motion** candidates
(the cell fixed at the root 2x2x2 block, the motion varied):

| gate | outcome |
|---|---|
| eligible (unitary + anomalies + five-thirds) | 509 |
| front speed exactly 1 Planck cell/stage | eliminates the mass rows (0.9167) |
| closure (Wilson product or conjugate) | eliminates every non-Wilson candidate |
| pass all five | **0** |

Failure patterns: 221 candidates fail the closure gate alone; 288 fail
front speed and closure together (the mass-stage rows: the extra internal
stage carries no bond transport, so the front moves 11 cells in 12
stages, 0.9167). The Wilson layer is the unique causal, closure-
consistent motion; the mass stages are subluminal. Zero survivors.

---

## 4. The downstream censuses on the candidate spectra (V168-V170)

The three downstream gates of the alt_cell branch (V158 unification,
V159 CP, V160 DM), re-applied to the **dynamics** candidate spectra:
| gate | candidates | survivors | reading |
|---|---|---|---|
| V168 forced meeting | 509 | 0 | one mode per node, not a complete $16+\overline{16}$ per node at a single momentum |
| V169 CP phase | 1336 | 0 | every effective internal unitary is real type ($J^2=+1$) |
| V170 DM carrier | 509 | 0 | lightest mode is a gapless node; lightest non-node mode at the massive scale |

**V168 - forced meeting.** For every gapless candidate the node content
is one mode per node (506/509 exactly; the three exceptions are the
degenerate Cartan row and two mass rows with 2-3 light modes per node),
not a complete $16+\overline{16}$ at a single
momentum. V64's proportional-shift no-go (the node eigenspace shifts
proportionally under any candidate coupling, so no meeting is forced)
kills the meeting for every motion. 0 meetings, 0 forced.

**V169 - CP phase.** 1336 candidates (all motion families): the
Frobenius-Schur indicator of every effective internal unitary is real
type ($J^2=+1$); no candidate admits a CP-odd phase (max CP-odd
$1.9\times10^{-16}$, max mirror commutator $1.3\times10^{-15}$). The
$q=0$ CP obstruction of V159 is motion-independent: changing the step
rule never makes the internal layer quaternionic.

**V170 - DM carrier.** 509 gapless candidate spectra: the mirror Z2
survives every motion (max $w$-commutator $1.3\times10^{-15}$), so every
mode is mirror-stable; but the lightest excitation of every candidate is
a gapless node (paired mirror content). A carrier would be a NON-node
mode that is light (eigenphase $<0.01$, far below the root GAMMA gap
$0.0365$) **and** separated from the rest of its spectrum (second
eigenphase $>0.02$) - a single flat mode, not a cone branch passing near
a node. The lightest non-node excitation of every candidate is such a
cone branch (min $10^{-3}$ step-phase, never isolated): no candidate
supplies an isolated neutral light mode, and the 118 GeV carrier is not
a motion artifact.

---

## 5. Verdict

The motion is as rigid as the cell and the signature. Of the five
suspicions, none survives:

- per-stage angles gap the zone or restate the root ($A\ge0.949$);
- the nodes are Cartan-protected, not ordering artifacts;
- a mass stage gaps the nodes but never yields a true cone and breaks
  the causal front budget;
- generic bond matchings keep the doubling;
- Floquet splittings and quantized Wilson layers destroy the gaplessness.

The hard tier (V167) eliminates every non-Wilson motion on closure and
front speed alone. The downstream gates (V168-V170) confirm that the
unification meeting, the CP phase, and the DM carrier behave on every
candidate spectrum exactly as they do on the root: the meeting is not
forced, the internal layer is real type, the lightest mode is a node.

**Boundary lines.** This branch sharpens two open questions of
`doc/final_publication.md`:

- **F7 (why this motion).** The Wilson six-stage cycle is not one option
  among many: it is the unique motion that is unitary, causal (front
  speed exactly 1), anomaly-free, five-thirds-closing, and closed under
  its own Cartan structure. Every alternative either gaps the fermions
  or breaks the light cone. F7 remains open as to *why* the root is
  chosen, but it is now bounded: the motion is rigid, not a free
  parameter.
- **O3 (the node doubling / isotropy).** The doubling is not a motion
  artifact. It survives all 729 angle assignments, all 12
  factorizations, all 500 matchings, every Floquet splitting, and the
  quantized Wilson family; it is Cartan-protected. O3 stands as an open
  problem of the *representation*, not of the *dynamics*.

The master construction is unchanged: the published node census
(12 nodes, $A=0.998$) is the unique gapless motion of the root cell.
