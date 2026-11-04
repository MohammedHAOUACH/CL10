# V386 - chiral ledger and the Nielsen-Ninomiya book: run notes

Verdict: **CUBOND-CHIRAL-LEDGER-NN-CLOSED** (tier 2). Wall 2.7 s,
host numpy, 256-dim. Suite 299.

## The measurement

Is the cone chiral like the Standard Model? Chirality operator on
the V383 engine: G5 = gamma_1 gamma_2 gamma_3. Tangent derivatives
exact: the factor identity F_j' = i gamma_j F_j gives
dU/dq_j = U (i gamma_j), so S_j = i dU/dq_j = -U gamma_j at every
site, no finite differences anywhere. The 16 gapless sites are
enumerated explicitly: the 8 of {0,pi}^3 and the 8 cubic corners
(+-pi/2,+-pi/2,+-pi/2).

- **T1** G5 is anti-Hermitian and G5^2 = -I exactly (0.0e+00); its
  spectrum is pure +-i with multiplicities **128/128** - the 256
  states split into equal helicity sectors: the cone is vectorlike.
- **T2** exact Clifford projection M_jk = (1/256) tr(S_j gamma_k)
  (the gammas are an orthogonal Clifford basis): scalar components
  max 1.2e-16; chi = sign det M on the {0,pi}^3 octet reads
  **(-1, +1, +1, -1, +1, -1, -1, +1)** - the sum over the 8
  Nielsen-Ninomiya doublers is **0 exactly**: the chirality book
  closes. At the 8 corner doublers the projection is degenerate
  (det M = 0, chi recorded as 0): there S_j = +-i G5 gamma_j is a
  product of two gammas, the non-Clifford residual max 1.0 recorded,
  not corrected.
- **T3** at the 8 NN sites the S operators close the Clifford
  algebra: min ||[S_i,S_j]|| = 2.000000 >= 1 at 8/8 sites.

## Reading

The content of the Clifford-bond family is NOT chiral like the
Standard Model - it is vectorlike: equal 128+128 helicities and an
exactly vanishing Nielsen-Ninomiya sum over the doubling octet. This
is booked as a measurement: the family carries the chirality book
closed with zero net chirality, so the chiral asymmetry of the
Standard Model must come from the selection, not from the cone -
which is exactly the question V387/V388 put to the certified Wilson
mechanism.
