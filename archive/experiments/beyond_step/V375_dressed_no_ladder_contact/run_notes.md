# V375 - dressed spectra vs the empirical ladder: run notes

Verdict: **DRESSED-NO-LADDER-CONTACT** (tier 2, honest negative).
Wall 1.2 s, CPU. Suite 288.

## Machine (exactly as booked)
Contact harvest MECHANICAL from the V52 provenance table: the
filter "quantity or note contains rung/ladder" returns exactly the
five booked rows (harvest gate): 0.019894 (ladder base, declared),
0.1999980 (g_obs rung-10 clock), 1.18 (DM-mass coefficient), 3
(n_2step seesaw exponent), 118 (sub-Planck DM scale, rung-10) -
value extraction is the first numeric token of the row value
string; the extracted set equals the booked set (gate).
Dressed eigenphase differences: V373 kernel spectra (dk_factored x
4 families + dk_fiber_ell2: 328 pairwise differences) + V374
lattice-sum per-row g = E_c - E0 (867 rows): 1195 differences x 5
contacts = T = 5975 machine-exact Decimal tests at relative
1e-12. Adjacency retest per booking for V374-family hits (same
(N,K) at eta +- 0.02, independent re-diagonalization); for the
eta-INDEPENDENT V373 kernel spectra the retest is adapted,
DECLARED: reproduction in the independent second measurement of
the same channels (factored route vs fiber route).

## Result
ZERO hits at relative 1e-12 (zero even at the pre-filter). The
dressing machinery is blind to the empirical ladder: the dressed
eigenphase differences {-96,-64,-32,... integer/algebraic values
and 1224 lattice-sum reals} share no contact with the five
published ladder numbers. Booked false-hit chance T x 4e-24 =
2.39e-20 << 1e-16. The V364 double-naming cross-check path was
never entered (nonzero-hit branch).

## Why this is a real negative, not a null test
The contacts span 0.02 to 118 (4.5 decades); a dressed difference
within relative 1e-12 of ANY of them would have been reported with
its family and adjacency evidence. Closest misses (far outside the
1e-12 band): 0.019894 vs 0.020699 (rel 4.1e-2), 0.1999980 vs
0.192971 (rel 3.5e-2), 1.18 vs 0.887356 (rel 0.25), 3 vs 0.887356,
118 vs 96. The V373 spectra are exact integers (differences 0, 32,
64, 96, ...); the V374 g values are generic reals.

## Machinery lesson
Decimal exactness at relative 1e-12 over artifact floats: use
Decimal(repr(stored_float)) so the test compares the artifact
digits, not binary artifacts of binary artifacts.
