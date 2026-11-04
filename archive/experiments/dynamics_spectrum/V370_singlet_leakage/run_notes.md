# V370 run notes - singlet leakage under the certified bond-current operator

Run 1 (cache path): GPU 1 pinned, 19.5 s, exit 0, tier 2, residual
True; total leakage 7.321011e-14.
Run 2 (fallback branch, cache deleted - fresh-checkout path): GPU 1,
644 s, exit 0, same checks green; total leakage 8.780863e-14. The
7.3/8.8e-14 band is the floor of the cached eigenbasis decomposition
(the degenerate kernel re-rotates between rebuilds), not a physical
drift; J.L = -64.000000 on both bases. Prose in catalogue, paper and
register quotes the band.
VERDICT: **LEAKAGE-ZERO**. T0 (corrected V369 decomposition) fires at
4.6e-15; kernel gate green (dim 2, ||C2 W|| 2.4e-15).

## Facts

1. **J.L W = -64 W exactly**: the bond-current interaction acts on
   the singlet kernel as MINUS the adjoint Casimir, scalar on both
   intertwiner columns - no multiplicity splitting at first order
   (B = W^dag J.L W eigenvalues -64.000000 / -64.000000). Via Gauss
   (L_a W = -(T1_a + T2_a) W) this is the statement that the total
   current squared of the neutral pair is the Casimir 64.
2. **Total leakage l = 7.3e-14**: the kernel is an invariant subspace
   of the certified bond-current operator to machine precision. The
   neutral pair is FIRST-ORDER ROBUST: e^{i eta J.L} preserves it and
   acts there as e^{-64 i eta}.
3. Channel weights: NOT RESOLVED and MOOT in this branch - a
   projector-free bound gives weight_c <= ||eps|| = 7e-14 for every
   channel. The first pass's l_36 = 90.51 was exposed as Lanczos
   cross-cluster contamination (90.51 = 64*sqrt(2) = ||J.L W||
   exactly: the contaminated cached cluster contained a full kernel
   direction); the booked correction resolves channel weights only
   in the LEAKING branch, where a clean channel basis would be built
   from scratch (the k=400 SM cache has O(1) cross-cluster overlap -
   recorded hygiene fact about which='SM' at large k on degenerate
   towers).

## Machine hygiene record (third appearance of the Lanczos-degeneracy
family in this programme)

V366 (intra-cluster Gram -1/2), V369 (same, three clusters), V370
(CROSS-cluster O(1) overlap at k=400 SM - the eigenvalues remain
trustworthy, the eigenvector spans across clusters are not). The
script now branches: zero-leakage branch uses only W (projector-
free); the leaking branch cross-orthogonalizes clusters sequentially
with a rank assert.

## What it means for door n°2

The bond-current interaction does NOT destroy the candidate (neutrality
robust, no splitting) and does NOT bind it at this cut (the -64 shift
is common to the whole 2-dim kernel; whether it opens a threshold
against the hybridized channel bands is the exact-spectral question:
the deformed step e^{i eta J.L} U_step restricted consistently needs
the lattice sum H = sum_b J_b.L_b on the torus - a genuine
interacting-lattice spectral programme, booked in pre-registrations
as V371 SHAPE without a run, since today's pieces (cache contamination,
single-bond restriction) do not bound it honestly).
