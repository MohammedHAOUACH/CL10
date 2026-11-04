# V387 - Wilson masses for the parasites of the Clifford bond: run notes

Verdict: **CUBOND-DOUBLER-MASS-SELECTIVE** (tier 2). Suite 300.

## The measurement

What happens to the masses of the sixteen parasites when the
certified Wilson factor (V65) acts on the Clifford bond?
U'(q) = (I_8 x W(q)) U_gamma(q) with W(q) =
diag(exp(-1j w phi(q) d32)), phi(q) = sum_j (1 - cos q_j) - the
exact shape of V65.step_wilson lifted from the 16-dim sector to the
32-dim module - and d32 the charges of
D(theta) = cos(theta) X + sin(theta) Y at the booked physical
junction angle theta = 1.7308.

## Declared implementation notes (post-first-run declared, physics gates untouched)

- **theta provenance.** The booking says: grep "theta" in the V86
  artifact first, use the booked constant if absent. The grep ran:
  the V86 artifact's "theta" field is its scan best 1.56 (not the
  physical junction angle); the values near 1.73 in that file are
  F12/F23 RESIDUALS, not angles. The supporting artifact readings
  of the physical angle are 1.73078 (V96/V187 solutions, key
  "theta") and 1.730784 (V187 theta_check/theta_star). The booked
  constant 1.7308 stands as the machine input, recorded in
  data.theta_used with data.theta_support = 1.73078.
- **lam16 and weights read from the V385 artifact** (booking:
  "poids V385 (artefact lu)"); the 32-dim Cartan eigenbasis is
  reconstructed in-script (same construction as V385) for the W
  factor, and the multiset is checked equal to the artifact's
  (0.0e+00 deviation).
- **N_pred reading.** The booked prediction "nombre de poids 32 de
  charge nulle ... x 8 facteurs de liaison" is implemented as
  N_pred = (# weights of |lambda . c5| <= 1e-9, c5 =
  cos(theta) ux[0] + sin(theta) ys[0]) x 8, computed and printed
  BEFORE the scan: **N_pred = 0** - at theta = 1.7308 NO weight of
  the 32 has zero D-charge; the smallest nonzero charge is
  **0.005216**. The gate "gapless == N_pred" is then evaluated at
  the state level: for every w >= 0.5, every one of the 15
  nonzero-momentum gapless sites of the family holds exactly 0
  states at phase {0, pi} within 1e-9. (The 16th gapless site is
  Gamma itself - the cone, measured separately in T2.)
- **T2 slope metric.** The booked quantity ||S'_j - gamma_j|| is
  read through the Clifford projection of the finite-difference
  slope S'_j (M_jk = 1/256 Re tr(S'_j gamma_k), the same exact
  projection the V386 ledger uses): the projection kills the
  off-Clifford term sin(q_j) (I x d32) - it is orthogonal to the
  Clifford basis by tr(P_a) = 0 - and the ungated tilt norm
  (max 1.0e-04 at q = 1e-04, consistent with sin(q) x max|charge|
  0.19) is recorded in data.T2_off_clifford_tilt_max. The gated
  statement is the booked parenthetical: the W factor is O(q^2) at
  the origin, so the cone slopes stay (1,1,1): max |M - I_3| =
  **5.0e-09 <= 1e-5** at every w.

## Results

- **T0** the V65 certification replays on the lifted basis: the
  16-dim X multiset equals {1x10, -3x5, 5x1} to 7.1e-15; the
  32-weight multiset matches the V385 artifact exactly; W(w=0) = I
  reproduces U_gamma to 7.8e-16.
- **T1 survivors.** N_pred = 0: at the physical angle every
  nonzero-momentum gapless site loses every exactly-gapless state
  at every w >= 0.5 - all 15 parasite sites acquire a mass.
- **T2 cone.** The cone at Gamma survives: gap 4.4e-16 at every
  w; slopes (1,1,1) to 5.0e-09.
- **T3 masses recorded** at w = 1.0, per site: gap floor
  **0.010432**, max **0.031296** (per-site values in
  data.t3_masses_w1). The pattern is the charge structure: the
  smallest mass is w phi(q*) x 0.005216 at the NN sites
  (phi = 2), the largest at phi = 6 corners.

## Reading

The certified Wilson mechanism solves the sixteen-parasite problem
on the Clifford bond: the parasites acquire mass at every strength
w >= 0.5 while the Gamma cone survives intact with velocities
(1,1,1). What survives exactly is the zero-charge content, and at
the physical junction angle there is none: the smallest charge
0.005216 is the clock tension written into the mass floor - the
physical reading 1.7308 does not close the neutral direction, and
the gap floor inherits that miss linearly.
