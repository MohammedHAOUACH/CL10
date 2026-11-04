# V384 - the cubond ladder functional: run notes

Verdict: **CUBOND-LADDER-NO-CONTACT** (tier 2, O8/F2). Wall 352 s
(256-dim, 101 loop evaluations at nmax 600). Suite 297.

## The design, and how it was tested

The certified hierarchy machinery (V206) lifted exactly to the
entangling engine. The lift is FORCED by linearity: the sector
hop equals ma_lambda = A_lambda C with the content operator
A_lambda = (1/45) sum_j P_lambda^{jj} (verified against V206's
own sector_hops at 1e-16 before booking anything), so the
dressing is sum_lambda phase_lambda (I_8 x A_lambda) G_engine
and the loop is B_n = (M G)^n exactly as certified.

Two design bugs died exactly where the booking planned they
would - the T0 anchor gate:
1. the loop was built as (M) x G once (MG = L@G = M alone):
   the anchor 1/5 read n*=10 coef 0.687 - wrong, caught at
   1e+00 detune against the certified 2.08811;
2. the booked norm-lift (Cartan content of the BOND TRACE) is
   IDENTICALLY ZERO on any entangling engine - the Pauli bond
   coefficients are traceless, and the engine's bond trace is
   only cos^3(ng) x I - the division by zero IS the
   measurement of that structural fact; the corrected lift
   (Cartan content of every bond block) equals the published
   cartan_norm term by term on the anchor (8/256 = 1/32).

After the fixes: the anchor reproduces the certified readings -
1/5-rad to 2.8e-13 (n* = 10, the rung itself), pi/2 to 5.8e-05,
pi/4 to 1.6e-03, 3pi/4 to 3.6e-03 (all within the booked
float64 ladder-drift convention of V382).

## The measurement

| clock | n* | coef | miss vs 2.0765641 |
|---|---|---|---|
| pi/4 | 206 | 2.3828879 | +3.06e-01 |
| pi/3 | 11 | 0.0115086 | -2.07 |
| pi/5 | 19 | 0.0000046 | -2.08 |
| 2pi/5 | 8 | 0.0161871 | -2.06 |
| pi/6 | 11 | 0.0000132 | -2.08 |

Dense scan, 96 points on [0.05, 1.60]: ZERO band-hits, no VOIDs.

Reading: the cubond engine moves the pi/4 ladder read the
furthest ever (2.7002 published -> 2.3708 carrier -> 2.3829
cubond, n* 208 -> 193 -> 206) but every entangling clock reads
either +0.31 above or catastrophically below; the coefficient
is NOT carried by bond geometry. The design space of the ladder
mechanism is closed one family wider: published engine (V206 /
F2 FIRED), module replacements (V377), the carrier family
(V382), and now the entangling Clifford bond (V384).

## What the cubond family still owns

V383's theorem stands untouched: the cone algebra exists on the
certified algebra. V384 only says its ladder does not read the
hierarchy. The two facts together sharpen the residue: the cone
and the coefficient are carried by DIFFERENT circuit classes in
everything built - which is, as of today, a measured statement,
not a guess.
