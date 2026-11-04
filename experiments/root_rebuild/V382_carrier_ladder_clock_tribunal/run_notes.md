# V382 - the carrier-ladder clock tribunal: run notes

Verdict: **CARRIER-LADDER-NO-CLOCK-RESCUE** (tier 2, O8; the
booked branch held - F2 stands COMPLETELY unchanged by the
carrier). Wall 26 s. Suite 295.

## The tribunal

The certified V206 hierarchy functional, imported whole (its
1/5-rad reading reproduced to 1.4e-13, bit-level; the three
closure-clock readings to max 3.6e-03 - the float64 ladder
drift V209 certifies against), with the clock twist REPLACED by
the carrier internal u_2i(pi/2, e_z) - hop phases, projectors,
Cartan norm, ladder walk all untouched.

Carrier-ladder reads at the five carrier-closed clocks:

| clock | n* | coef | miss vs 2.0765641 |
|---|---|---|---|
| pi/4 | 193 | 2.3708397 | +2.94e-01 |
| pi/3 | 15 | 0.0049570 | -2.07 |
| pi/5 | 33 | 0.2733869 | -1.80 |
| 2pi/5 | 9 | 0.6437110 | -1.43 |
| pi/6 | 11 | 2.7153494 | +6.39e-01 |

Dense scan, 320 points on [0.05, 1.60]: ZERO band-hits (1e-4
constants band) anywhere - closed or not. The rescue set is
empty.

## What it means

- The carrier does change the ladder (the pi/4 read moves from
  the published 2.7002078 to 2.3708397 - the 30.0% miss narrows
  to 14.2%), but never within 2900x the band, and no clock on
  the whole scanned interval hits the observed value.
- The clock arbitration of V233 (the 2x2 score: closure picks
  pi/4, data picks 1/5 rad, no closed clock reads the
  coefficient) survives the carrier discovery intact. The
  carrier is a STRUCTURAL discovery - the three legs co-occur
  (V380) - with NO clock consequence.
- Recorded for the record: the carrier ladder is the first
  certified machinery that moves the pi/4 coefficient read
  meaningfully (2.70 -> 2.37): the direction the carrier twists
  the loop matters, and the residue to 2.0766 is a measured
  0.294, not a fitted one.

## Declared correction (same day)

T0's booked replica threshold (1e-9) was mis-set against the
machinery's own certified lesson (float64 closure-clock reads
detune ~1e-3 at n* = 208 first-crossings - V209's raison
d'etre); corrected pre-reissue to: 1/5-rad bit-level + closure
reads < 5e-3, declared in the pre-registration file. The status
branch was also fixed (an empty rescue set must not print
CLOCK-HIT). T1/T2 tolerances untouched; their readings identical
before and after.
