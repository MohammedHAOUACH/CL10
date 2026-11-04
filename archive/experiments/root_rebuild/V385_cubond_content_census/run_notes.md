# V385 - content census of the Clifford-bond cone: run notes

Verdict: **CUBOND-CONTENT-VERBATIM** (tier 2). Wall 0.1 s, host
numpy, 256-dim. Suite 298.

## The measurement

The booking asked: which particles live in the cone of the Clifford
bond? The machine is the V383 engine's 256 states, bond(8) x
internal(32), read through the five 32-dim Cartan generators
C_k = i h_a h_{a+1} (a in {0,2,4,6,8}, clifford.generators), common
eigenbasis via the sqrt(1,2,3,5,7) combination - the V65
cartan_weights pattern lifted from the 16-dim Weyl sector to the full
32-dim module.

- **T0** the gamma algebra re-certifies exactly (gamma^2 = I and
  anticommutators 0.0e+00); the certified 16-dim reference form
  (V65.cartan_weights) reproduces its own 16-weight multiset,
  off-diagonal residual 2.3e-13 (<= 1e-9 booked).
- **T1** the 32 internal weights are EXACTLY the 16 certified
  weights plus their mirror: multiset equality with matched-component
  deviation 3.7e-15, nothing missing, nothing extra (<= 1e-9 booked).
  The weights are the 32 sign vectors of {+-1}^5 - the 16 of the
  Weyl sector (odd sign product) and the 16 of its conjugate - the
  V52 mirror of the internal cell, verified here on the entangling
  engine for the first time. The 32-dim Cartan diagonalizes to
  1.1e-13 off-diagonal.
- **T2** at Gamma the engine is the identity: U(0) = I exactly
  (0.0e+00), all 256 states at phase 0 (gap 0.0e+00); the 256-weight
  census is eight copies of the same multiset, equal to
  8 x (lam16 + (-lam16)) exactly.

## Reading

The cone of the Clifford bond carries exactly the certified cell
inventory: no new particle, no lost one. The existence theorem
(V383) changes the inventory in no way - the content question of the
booking is answered by census equality, not by construction. The
mirror-pair structure (16 + 16-conjugate) is the V52 matter/antimatter
statement in the weight language, and it is the full internal content
of the 256-state bond engine.

## Declared implementation notes (post-first-run declared, physics gates untouched)

- lam16 recomputed in-script via V65.cartan_weights() and
  auto-recorded: the V65 artifact carries no "weights" key (checked
  by grep), the fallback named in the booking.
- Weight components are exact sign vectors in {+-1}^5 (deviation
  from +-1 at most 2.0e-15), so multiset equality is taken after
  rounding components to integers; the booked 1e-9 gate is applied
  to the matched-component deviation between the unrounded multisets.
