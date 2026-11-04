# V378 - the crystalline dilemma, witnessed: run notes

Verdict: **CRYSTALLINE-DILEMMA-WITNESSED** (tier 2, witness census of
the V377 named surviving structure). Wall 18.8 s, host numpy (V337
ran the same machinery on cupy GPU 1: the 9 overlap rows reproduce at
max |d gapM| 0.0e+00 - bit-exact across backends). Suite 291.

## The class (all booked pre-run)
The 2I engine (cl10/supercell_2i, the repository's only clock-free
closure engine) over omega = pi*j/5, j = 1..9 (the 5-fold angles,
the binary-icosahedral element at finite order) x beta in {G10,
pi/5, 2pi/5, 1/5 rad}: 36 points. The conjunction C: unitarity
<= 1e-11 + gapless M (min |phase| <= 1e-9 over all M-points) + the
published R cluster signature [32, 48] + internal mirror <= 1e-10 +
closure at GAMMA, order <= 240.

## The witness (measured)
- The internal 2I layer closes at EVERY class point, order <= 10
  (36/36): the clock is born free in the internal layer at every
  5-fold angle. The full step closes at 27/36 - exactly the
  rational-multiple-of-pi beta columns; beta = 1/5 rad NEVER closes
  (0.2 rad/stage vs 2pi: incommensurate, closure_k None to cap 240).
- gapM = 0 fires only at (omega/pi, beta) in {(0.4, 2pi/5),
  (0.8, pi/5), (1.2, pi/5), (1.6, 2pi/5)} - everywhere the R
  signature FAILS.
- Where the R signature holds (the published beta class), gapM >=
  0.104720 (the floor equals the V337 row value to backend roundoff).
- Conjunction C: 0 of 36. The DILEMMA-BROKEN branch - the only
  module-side solution carrier booked - did not fire.

## Reading
On the only engine the repository has built where closure costs no
clock, the gapless M-points and the published R multiplets never
co-occur. Together with V377 (no built member moves a contact) the
module-side of the ten-to-one is CLOSED AS A SEARCH: a module-side
solution, if it exists, must sit outside everything built - it must
carry a clock-free closure core AND a crystallographic walk ladder
AND the published R multiplets at once, and every partial carrier
built so far (2I, star, pentagon, icosahedral cell) fails at least
one leg by measurement.

## Machinery lessons
- Closure must be stated at the level it is measured: internal
  layer vs full step gave 36/36 vs 27/36 - the booking conflated
  them (post-run correction booked same day).
- Transport phases at rational RADIANS (0.2) are incommensurate
  with 2pi; only beta rational in PI closes the staggered swap
  layer. A "fitted pin" beta is not a clock candidate at all.
