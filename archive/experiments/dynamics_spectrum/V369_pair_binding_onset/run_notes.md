# V369 run notes - pair-binding onset at the certified layer interaction

Run: GPU 1 pinned, 526 s, exit 1, tier 2, residual True.
VERDICT (as booked, kept verbatim): **V369-ANCHOR-DRIFT** - T0 and T1
fail, T3 slope 0. The machine is sound; the BOOKED MATHEMATICS was
wrong. Honest negative, committed as run.

## What failed and why (the algebra, on the record)

Booked identity T0: J.L = 1/2 (C2_tot - C2_1 - C2_2 - C2_L). FALSE:
    1/2 (C2_tot - C2_1 - C2_2 - C2_L) = J_1.J_2 + J_1.L + J_2.L,
the TOTAL pair of cross-contractions - it contains the current-
current term J_1.J_2 which the published bond-current Hamiltonian
(H = sum_b J_b.L_b, V101) does NOT contain. Measured rel defect 1.1
(T0) and ||(J.L+77)W|| = 1.9 (T1): the kernel shift is NOT -77 for
the bond-current operator; -77 belongs to the operator
(C2_tot - 154)/2 which is a DIFFERENT interaction (it includes the
particle-particle channel). The booked 18-eta Casimir-slope mechanism
therefore does not exist for the published interaction: measured
slope -0.000.

Second, independent formal defect (found by the machine's flat
eta-profile, gap ~1e-5 at EVERY eta): the additive phase bookkeeping
is invalid - e^{i eta J.L} does not commute with the step's internal
layer (both carry the link algebra), so spectra are not sums of
shifts. The channel bands (40/149/108 internal phases x 1536 bond
points) cover the circle densely at the sampled grid; additive-shift
comparisons measure only grid spacing.

## What the machine certifies (all green, reusable)

- Channel tower of E4 under C2_tot (k=400, cupy matvec):
  {36 x 40, 64 x 149, 80 x 108, 84 x 101} - exact (C2 resid <= 5e-12);
- channel projectors orthonormal to 3e-15 after the V366-precedent
  Gram -1/2 (the Lanczos degenerate-cluster defect appeared again at
  2.9e-01 / 5.8e-01 and the booked hygiene fixed it);
- V368 eta=0 anchor 2.017083e-05 reproduced to 1e-9 by the SAME
  machinery; V101 published E0 = -5802.398188 eta re-read.
- ARPACK cache results/v369_sm400.npz (~300 MB, gitignored): reruns
  of the analysis now cost seconds.

## What this says about door n°2 (careful, no overreading)

The naive Casimir shortcut to a threshold is DEAD: the published
interaction (bond-current x link) does not shift singlet and channel
sectors by Casimir arithmetic, and the phase-additive binding test is
formally invalid for it. The binding question in its correct form is
a non-commuting-product unitary spectral problem (eigenphases of
e^{i eta J.L} . U_step) - a bigger machine, and NOT claimed today.
What is immediately decidable and booked as V370: the exact
first-order coupling of the interaction to the singlet subspace -
the leakage ||(1 - P_ker) J.L W|| and its channel decomposition.
If the bond-current interaction preserves neutrality (leakage 0)
the door stays shut differently; if it leaks, the interaction does
not bind the pair - it mixes the neutral state into charged
channels, which is also a closure of the naive route.
