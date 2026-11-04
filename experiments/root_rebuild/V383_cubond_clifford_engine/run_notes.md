# V383 - the Clifford-bond (cubond) engine: run notes

Verdict: **CUBOND-EXIST-DOUBLERS-PROLIFER** (tier 2, O9). Wall
4.6 s, host numpy, 256-dim. Suite 296.

## The theorem (existence, certified exact)

On the certified algebra alone - the V337 spatial-plane spin
triple (S_aS_b + S_bS_a = 0.0e+00, (2S_a)^2 = I.0e+00) and the
three bond-link Paulis of V138's 8-dimensional bond factor -
the operators

    gamma_j = P_j (x) 2 S_j,   P = (X_1, Y_2, Z_3)

form an EXACT Clifford set on 256: gamma^2 = I at 0.0e+00,
anticommutators 0.0e+00, Hermiticity 0.0e+00. The engine
U(q) = prod_j (cos q_j + i sin q_j gamma_j) is unitary (2.2e-16),
its Fourier support is exactly the 8 corner coefficients of one
lattice step, and its velocity algebra - the thing 285
experiments of Kronecker bonds measured at exactly 0.0e+00 and
V185 proved necessary for a cone - measures ||[gamma_i,gamma_j]||
= 2.0 at all three pairs. The spectra of the gamma_j are exactly
{+-1} x 128+128: principal velocities (1,1,1) - against the
published node's (0.00000, 0.14189, 0.14189). Leading cones are
exactly isotropic at speed 1 (8 random directions at a = 1e-3:
every phase within 1.9e-07 of +-*|q|, under the booked 3a^2
band-bend bound).

## The honest ledger

- The carrier (V380): (I + 2iS_z)/sqrt2 = u_2i(-pi/2, e_z) at
  1.1e-16 - the carrier is a half-angle factor of the family,
  Pauli-dressed. The sign is booked (-pi/2, not +pi/2).
- Doublers PROLIFERATE: the {0,pi}^3 octet (8, counting the
  origin cone site) plus the (pi/2,pi/2,pi/2) corner - exactly 9
  gapless sites on the 5^3 census grid, 16 over the full BZ by
  cubic symmetry. Nielsen-Ninomiya is not evaded; it is
  over-satisfied. The pi/2-type face/edge sites are gapped with
  floor pi/2 (1.571).
- Front radius sqrt(3) (corner hops) against the published 1 -
  recorded, not gated.
- The ladder functional on an entangled engine is the NAMED open
  design problem: the V206/V204/V209 machinery is built on
  internal-only 32-dimensional twist steps; the cubond engine is
  256-dimensional bond-entangled and was NEVER run against the
  published anchors because it is not the published engine.
  F2/O8 stand untouched by V383.

## Method note (declared)

The construction was derived on paper (symmetry reasoning: Pauli
commuting (x) spin anticommuting = Clifford), verified by one
user-authorized 9-second probe, then booked in full and run
formally - the probe and the formal run agree on every number.
The cone-metric bug of the probe (comparing negative phases
against +|q|) was caught before booking and the corrected metric
is the one gated here.
