# V379 - the hierarchy-coefficient shadow tribunal: run notes

Verdict: **E-SHADOW-COEFFICIENT-NO-WEB-CONTACT** (tier 2, O8-facing;
the E-WEB-CONTACT-FOUND branch did not fire). Wall 0.5 s (mpmath
40 dp + small numpy). Suite 292.

## The trap, recorded (T1)
3*exp(-1/e) = 2.0766018826660391... vs the observed hierarchy
coefficient 2.076564111173832 (V233 field): deviation 3.777149e-05,
inside the booked +-1e-04 observation band, ZERO fitted parameters.
As a numerical fact this is the cleanest hit ever found against the
theory's measured numbers, and it is recorded as such.

## The tribunal (T2/T3)
- Required-exponent census: with the V312 p0 F_cross reads, the
  exponent the observed coefficient REQUIRES, n_req = ln(obs/3) /
  ln(16 pi F), reads 35.9408 (pi*1/4), 35.6340 (pi*3/4), 13.8336
  (pi*1/2) - the nearest web integer for the pi/4 read is 36 (the
  module dimension) at distance 0.0592; the n uncertainty (obs
  band + F-cross precision) is ~0.01, so the 36 reading misses by
  ~6 sigma. No rung, no vertex invariant, no walk length reads it.
- Web membership: exp(+1/e) and exp(-1/e) tested against the full
  certified web basis (22 numbers: 16pi, 1/pi, phi, loop, the three
  F_cross reads, rho, e, 1/e, and the structural integers
  1,6,10,17,23,30,32,36,45,64,141,208): pslq at maxcoeff 200 with
  40-dp verification and degeneracy guards: ZERO relations.
  e-type transcendentals are off the web - the hit has no
  mechanism: SHADOW.

## The near-miss map (T4)
At the four 2I gapless points the R cluster signature is [16, 48]
(published: [32, 48]): the 48-cluster survives intact at every
point; the 32-cluster splits exactly dyadically 16+16. The shape
constraint on any three-leg carrier: the clock-free closure core
must act trivially enough at the R corners to preserve the
32-cluster - its defect at the best built points is precisely one
dyadic split.

## Machinery lessons (two pslq degeneracies, both real traps)
- pslq([1, X, w]) with INTEGER w fires trivially (-w*1 + 1*w);
  membership requires the target to participate (nonzero target
  coefficients).
- The 4-column form fires on collinear columns X and X*w: the
  grouped coefficients (a + c w) and (b + d w) must not both
  vanish. Booked guards now: target participation + grouped
  non-degeneracy + 40-dp exact verification.

## Reading
The two searches that the closed module side leaves - find a
number-relation for the open values (V379 kills the cleanest
candidate), or build the three-leg carrier (V377/V378 map its
shape) - now have certified negative/positive boundary data:
the open numbers are not reachable by any bounded relation to
the constructed web, and the carrier, if it exists, must preserve
the R 32-cluster.
