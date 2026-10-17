# Alt-nonperturbative exploration: exact flux dressing and the 4-site gauged sector against W9

Exploration programme `alt_nonpert` (merged into `master`), following
`alt_gauged`
([alt_gauged_exploration.md](alt_gauged_exploration.md)),
`alt-interactions-exploration`
([alt_interactions_exploration.md](alt_interactions_exploration.md)),
`alt-dynamics-exploration` ([alt_dynamics_exploration.md](alt_dynamics_exploration.md))
and `alt-cell-exploration` ([alt_cell_exploration.md](alt_cell_exploration.md)).
Experiments V180-V183, category `experiments/alt_nonpert/`, tier 4,
0 new continuous parameters (eta = g_10 = 1/5, the root angle).

## The question

Five routes have failed to isotropize the 12 M-node cones of the master
automaton (A = 0.998): signature, cell, motion, the minimal collision
(V171-V174), and the gauged bond-current interaction at one loop / under
coarse-graining (V175-V179, zero candidates). The alt-gauged report
names the one route never measured: the **non-perturbative gauged
sector** - the exact flux dressing of the cone, and the two-fermion
spectrum one volume up (the IR-decoupling question is a VOLUME
question). This branch computes both exactly, no truncation.

## V180 - the exact cone dressing by the flux operators (PASS)

The cone node (pi/2, pi/2, 0) dressed by the 24 exact flux operators of
the V101 kernel: the anchor sum of flux projectors is 64 I (err
2.26e-14), the sum rule over the 24 bonds is 6.000000, and the dressed
cone speed is measured against the free one on the V115 sphere:

| dq | A_free | A_np | dA |
|---|---|---|---|
| 1e-2 | 0.846307 | 0.846307 | +2.1e-12 |
| 1e-3 | 0.997985 | 0.997985 | +1.3e-11 |
| 1e-4 | 0.100185 | 0.100185 | +9.1e-12 |

dA is at machine precision (2-route 8.5e-13; eta -> 0 gives dev 0;
eps anisotropy 3-6e-15). **The exact one-order flux dressing does not
move the W9 anisotropy at all.**

## V181 - the 4-site gauged kernels and S-matrices (PASS)

The V99/V101/V175 machinery generalized one volume up: the 4-site OPEN
chain (Is - Ts - Ts - Is, 144 states) and the 4-site RING (four Ts,
cyclic link identification, 1296 states). New sandwich tensors UQ
(Ts-Ts current, keeps the bond slot); every ring bond uses UQ on the
J-site tensor and V2 on the other, convention fixed by the ring Gram.

| system | states | Gram rank | E0 (eta) | deg | gap (eta) | H 2-route | S unitarity | phase dev / tol |
|---|---|---|---|---|---|---|---|---|
| N=3 chain (V101 anchor) | 24 | 24 | -5802.3982 | 6 | 42.0906 | - | - | - |
| N=4 open chain | 144 | 144 | -1664197.02 | 1 | 12787.60 | 2.1e-14 | 2.7e-15 | 2.2e-8 / 9.2e-7 |
| N=4 ring | 1296 | 976 | -6.076e8 | 1 | 2.751e8 | 1.3e-15 | 4.2e-15 | 7.0e-4 / 1.5e-1 |

Anchors: N=3 reproduces the V101 JSON (E0 dev 4.2e-11, gap dev
9.1e-13); Gauss residual 5.0e-13 on every basis tensor. The S-matrix is
built from the whitened H by exponentiation (unitary by construction,
verified < 1e-12); the 2-route phase comparison uses two genuinely
different whitening factorizations (Cholesky/QR-rectangular vs
eigh-rectangular) and the tolerance is the pre-registered
whitening-noise bound 100 * eta * ||H|| * eps * kappa(N) - the ring
Gram is ill-conditioned (kappa 5.6e4, rank 976/1296: the cyclic
identification makes 320 ring states dependent), so its bound is
looser; the deviations sit inside it.

Structural facts recorded: fusion channels are PRESERVED on the chain
(mixing 3.4e-11, as at N=3) but strongly MIXED on the ring (~8): the
cyclic link identification couples the r = 1/45/210 sectors. The ring
spectrum is NOT isospectral with the chain: the eta-units of V101 are
volume-dependent (below).

## V182 - IR decoupling across volumes (INCONCLUSIVE, recorded honestly)

Pre-registered clauses: LIGHT if some two-fermion gap < 1e-2 eta;
MASSIVE if gap >= 42 eta at all volumes AND the volume drift < 1 eta.
Measured: gap 42.09 (N=3) -> 12787.60 (chain) -> 2.75e8 (ring) eta;
drift >> 1. Neither clause fires: **INCONCLUSIVE**. What IS decisive:
the LIGHT branch is false (no light two-fermion mode exists at any
volume), and the drift clause fails for a structural reason - the V101
eta-units are NOT extensive: the per-link Gram shrinks ~30x per added
link, energies grow superlinearly, and gap/|E0| stays constant
(0.0073 -> 0.0077). The flux sector is massive at every volume probed;
no two-fermion excitation can dress the IR cone.

(The rollup labels V182 `CRASH (rc=1)` because the repo convention ties
the exit code to the verdict - this is the honest FAIL verdict written
in the JSON, not a crash.)

## V183 - the non-perturbative isotropy gate (PASS, 0 candidates)

9 census rows (V180 dressing + V181/V182 spectra), candidacy = A_np <
1e-6 or |dA| > 0.5 or a light gap < 1e-2: **zero candidates**. Hard
tier untouched (front speed 1, unitarity, closure). No
non-perturbative gauged-sector mechanism isotropizes the W9 cone.

## The engine-bug hypothesis, refuted

Could a wrong dispersion hide an isotropic physics? The baseline
dispersion is exactly O_h-invariant: max |E(Pq) - E(q)| = 2.0e-15 over
all 48 signed permutations of the momentum (generic q and the cone
node), on top of the existing anchors (Weyl algebra 8.9e-16, Casimir
45 I, plane_product(T16, g_10) = U16 to 5.6e-16, front speed exactly 1,
anomalies 0, hypercharge 5/3). The A = 0.998 anisotropy is exactly
O_h-covariant - not an artifact.

## Verdict

The non-perturbative gauged sector is now measured exactly: the flux
dressing leaves the cone anisotropy untouched (dA ~ 1e-12), the
two-fermion flux spectrum is massive at every volume (gap/|E0|
constant), and the ring's channel mixing does not produce a light
isotropic mode. W9 (cone isotropy) remains open within the gauged
sector; the remaining variable is the SPECIFICATION of the step and
cell itself - see `alt-step-search`
([alt_step_search_exploration.md](alt_step_search_exploration.md)).
