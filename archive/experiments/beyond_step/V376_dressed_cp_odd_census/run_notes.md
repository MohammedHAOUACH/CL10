# V376 - CP-odd census of the dressed class: run notes

Verdict: **CP-EVEN-DRESSED-CLASS-CLOSED** (tier 2, theorem
extension). Wall 40.6 s, CPU (mpmath 40-digit verification tier).
Suite 289.

## Class members (gated identifications, declared)
- W U: W is the published dressed hop of V251, W = v101_gap x G10
  (gate: V251's omega rows are +-W = gap x 0.2, the then-clock
  G10 = 0.2 recorded; current constants give W = 33.057904663).
  The step is the same 45 plane rotations at the dressed angle;
  unitarity defect 3.1e-15.
- e^{i eta J.L} U: on the certified internal 32-dim space the
  bond-current operator is the MEASURED scalar sum_a T_a.T_a =
  C2 I, C2 = 45.000000 (dev 0.0e+00) - so the dressing is a
  global phase and the 45 rotation factors are unchanged. That
  is a measurement, not an assumption: without the scalar gate
  the question would be ill-posed on this space (a non-scalar
  dressing would leave the booked rotation-factor census without
  a factorization).
- V373 full generator: its internal-space action reduces to the
  same measured scalar dressing (the string Casimir and H_hop
  act on flux/lattice sectors outside the internal invariant
  space) - declared reduction.

## Result
Every booked invariant (2- and 3-rotation products full, the
n4_sample = 200000 seed-0 4-rotation sample, every B-violating
h_i R_ab h_j R_cd sandwich full, the democratic Yukawa) is <=
1e-12 at every eta in {0, 0.02, ..., 1} for all three members,
with the mpmath-40-digit tier agreeing with float64 to 2.9e-16.
The eta = 0 column reproduces the published V52 census exactly
(0.0e+00 on every entry - the re-use requirement). The dressed
class carries no CP-odd content the bare QCA lacks: the
algebraic cancellations (odd numbers of symmetric traceless
plane generators in every imaginary part) are angle-blind, so
the dressed angle W and the scalar dressings cannot break them.

Recorded as data, explicitly NOT booked invariants: the certified
channel-space spectra (V371/V373 dk) carry eta-dependent
relative phases (max |Im Tr| over the channel record 22.81), but
channel-space traces are not invariants of the booked internal
class and CP is not defined on the channel labels by the booked
machine. The named eta_CP residual of Table 6 stands. The worst
booked invariant over the whole grid is 2.7e-15 (roundoff scale)
- the closed branch is not marginal.

## Machinery lessons
1. The booked CP-operator gate "w* = -w" is a basis phase
   convention of the PRE-rebuild Clifford representation. The
   current published 32x32 Pauli source carries a REAL volume
   element (w* = +w, recorded as data). The basis-invariant
   content - w^2 = -1 and exact anticommutation with all 10
   Cl(10,0) generators, both at 0.0e+00 - is what gates; a
   convention string must never be gated as physics.
2. Verification tiers: evaluate float64 and mpmath on the SAME
   random products (not two independent samplings), or the
   "deviation bound" measures sampler luck instead of precision.
3. mpmath matrix API: rows attribute is the dimension (int),
   there is no .shape/.trace() on this mpmath version.
4. The factor sets are eta-independent because the dressings are
   measured scalars - the per-eta booking is honored by recording
   the (declared) once-per-set verification tier at every eta;
   the float census itself runs per eta.
