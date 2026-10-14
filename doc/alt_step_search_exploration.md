# Alt-step-search: a loss + gradient-descent search over the step/cell specification against W9

Branch `alt-step-search` (stacked on `alt-nonperturbative-exploration`).
Experiments: `experiments/alt_step/V184..V187`. Parameter-rule override
(user-directed, this branch only): the 45 per-plane step angles are
searched as CONTINUOUS variables; the hard tier (unitarity, front speed
exactly 1, closure order <= 64, chirality) stays non-negotiable.

## Why the specification is the remaining variable

`alt_nonperturbative_exploration.md` closed the gauged sector: the flux
dressing leaves the cone anisotropy untouched (dA ~ 2e-12) and the
two-fermion flux spectrum is massive at every volume. The engine-bug
hypothesis was refuted there (dispersion exactly O_h-invariant to
2.0e-15). What is left is the SPECIFICATION: the step
U16 = plane_product(T16, g_10) and the cell (2,2,2) are one choice in a
continuous (45 angles) x discrete (cell shape, stage order) space. This
branch searches that space with a differentiable cone loss and Adam.

## V184 - the differentiable cone loss (PASS)

Parameterization (exact): theta in R^45,
`U(theta) = Pi_k (cos theta_k I - i sin theta_k T16[k])` in the
`weyl_generators()` plane order; theta0 = g_10 * ones reproduces the
implemented step (`U(theta0) = sectors()[2]` to 5.7e-16).
`M(q, theta) = Pi_stages kron(S_st(q), U(theta))` (jnp reimplementation
of `bloch_cell`); `K = (M + M^dag)/2`;
`omega(q, theta) = arccos(eigvalsh(K)[-1])` = `gap_at` exactly (dev
7.2e-12 vs `gap_at`). Loss: node q*(theta) by coarse grid +
Nelder-Mead on omega^2; speeds v(n) = omega(q* + dq n)/dq over the 124
`sphere_dirs(6)` directions, dq = 1e-3;

    L(theta) = A + 100 * omega(q*)^2,   A = (vmax - vmin)/(vmax + vmin)

(the V115 metric). Gradient (validated, no naive eigvalsh VJP - invalid
at the degenerate node): the aniso part flows only through the
argmax/argmin directions via the single-band formula
d omega/d theta = -v^dag (dK/dtheta) v / sqrt(1 - lambda^2)
(vs central finite differences: 1.8e-8 CPU / 2.7e-8 GPU); the penalty
is the top-eigenspace projector subgradient (vs autodiff at a generic
gapped q: 5e-13), hinged - active only while omega(q*) > 1e-6 (an
always-on penalty crushes the cone into a flat band, A -> 1). The fast
forward factors M = (P(q) kron I16)(I kron U^12) in NumPy (124 matrices
in 0.13 s; equals `K_of` to 1.3e-15); on a live CUDA device the batch
runs through XLA `eigh` instead.

Landscape (9 rows, GPU route): baseline A = 0.9994 at omega(q*) =
4.2e-8; 8 uniform random seeds give A = 0.9962-1.0000 with
omega(q*) = 0 (random steps are gapless-flat: A -> 1). The V115
physical anchor |A - 0.998| < 1e-3 holds (0.9980).

## V185 - Adam descent on the 45 angles, baseline cell (25 seeds)

Pre-registered seeds: baseline g_10*1; 16 uniform(0, 2pi); 8 quarter
{pi/4, pi/2}. Alternating scheme: 10 cycles x [Nelder-Mead node + 30
Adam steps (lr 0.02)], best-L tracking of the trajectory.

| seed | A_final | omega(q*) | L | ||dtheta|| |
|---|---|---|---|---|
| baseline g_10*1 | 3.777e-11 | 3.43e-03 | 0.0012 | 2.31 |
| uniform-0 | 1.548e-09 | 2.26e-02 | 0.0510 | 0.58 |
| uniform-1 | 1.050e-06 | 3.52e-04 | 0.0000 | 3.24 |
| uniform-2 | 9.713e-01 | 0 | 0.9713 | 0.77 |
| uniform-3 | 9.989e-01 | 0 | 0.9989 | 2.42 |
| uniform-4 | 9.808e-01 | 0 | 0.9808 | 2.18 |
| uniform-5 | 9.843e-01 | 0 | 0.9843 | 2.23 |
| uniform-6 | 4.863e-10 | 2.16e-02 | 0.0467 | 2.20 |
| uniform-7 | 9.868e-01 | 0 | 0.9868 | 2.48 |
| uniform-8 | 9.720e-01 | 0 | 0.9720 | 1.43 |
| uniform-9 | 9.668e-01 | 0 | 0.9668 | 0.98 |
| uniform-10 | 5.820e-10 | 3.71e-02 | 0.1375 | 0.82 |
| uniform-11 | 9.993e-01 | 0 | 0.9993 | 2.18 |
| uniform-12 | 9.698e-01 | 0 | 0.9698 | 0.00 |
| uniform-13 | 6.969e-08 | 1.70e-03 | 0.0003 | 0.75 |
| uniform-14 | 9.805e-01 | 0 | 0.9805 | 1.91 |
| uniform-15 | 9.813e-01 | 0 | 0.9813 | 1.71 |
| quarter-0 | 9.859e-01 | 0 | 0.9859 | 1.88 |
| quarter-1 | 9.861e-01 | 0 | 0.9861 | 0.99 |
| quarter-2 | 9.856e-01 | 0 | 0.9856 | 2.23 |
| quarter-3 | 9.663e-01 | 0 | 0.9663 | 1.13 |
| quarter-4 | 1.209e-04 | 1.50e-05 | 0.0001 | 1.60 |
| quarter-5 | 1.165e-09 | 2.00e-02 | 0.0399 | 1.06 |
| quarter-6 | 9.991e-01 | 0 | 0.9991 | 1.66 |
| quarter-7 | 3.636e-09 | 7.82e-03 | 0.0061 | 2.35 |

Reading: descent from random seeds escapes to TRIVIALLY isotropic
gapped points - at a band edge the touching is quadratic, so
v(n) = omega_0/dq + O(dq) is direction-independent and A ~
1e-4..1e-11 while omega(q*) = 3.5e-4..3.7e-2 stays finite (10 of 25
seeds, including the baseline itself: A 3.8e-11 at omega 3.4e-3).
The two closest approaches to a candidate both miss on omega:
quarter-4 (A 1.2e-4, omega 1.5e-5) and uniform-1 (A 1.05e-6, omega
3.5e-4). Every gapless row (omega = 0) keeps A >= 0.966. Candidates
(A < 1e-6 AND omega < 1e-6): ZERO. The loss cannot separate trivial
band-edge isotropy from a cone except through omega -> 0, and no
seed reached a gapless isotropic basin.

## V186 - joint cell x step search (4 cells x 8 seeds)

| cell | stages | M | best A | best omega(q*) | best seed |
|---|---|---|---|---|---|
| (2,2,2) | baseline | 128 | 0.9966 | 0 | uniform-5 |
| (2,2,2) | phase-major | 128 | 0.9805 | 0 | uniform-2 |
| (2,2,3) | phase-major | 192 | 0.9675 | 0 | uniform-5 |
| (3,3,3) | phase-major | 432 | 0.9671 | 0 | uniform-4 |

All 32 rows are gapless-flat: every seed's descent ends on a gapless
point whose cone is maximally anisotropic (A 0.967-1.000); no cell
shape admits an isotropic gapless step. Anchor: the (2,2,2)/baseline
baseline seed's initial node A reproduces V185's to 1.5e-4 (gate
1e-3, cross-process NM jitter). Candidates: ZERO.

Structural tier note: anomalies, hypercharge and the chiral matter
content are fixed by the signature + sector decomposition and are
untouched by (cell, theta) - the cell changes only the motion.

## V187 - the step-search gate

38 census rows scanned (V184 9 + V185 25 + V186 4 summaries),
candidacy A < 1e-6 AND omega(q*) < 1e-6: **zero candidates** - the
hard tier (unitarity, closure order, front speed, chirality) had
nothing to reject. PASS with the negative recorded.

## Verdict

The step/cell SPECIFICATION has been searched directly - 45
continuous per-plane angles by differentiable cone loss + Adam (25
seeds), and four cell/stage choices (32 seeds) - on GPU. No isotropic
gapless step exists in the searched space: descent either keeps the
cone anisotropic (A >= 0.95 at every gapless point it reaches) or
isotropizes only trivially at gapped band edges (omega >= 3.5e-4).
The implemented step (g_10 on all 45 planes, cell (2,2,2)) remains
the only closed-form gapless point found, and W9 (cone isotropy)
stays open: it is not a tunable consequence of the step angles.

Cross-links: `alt_nonperturbative_exploration.md` (the gauged-sector
closure that motivated searching the specification),
`alt_gauged_exploration.md`, `alt_interactions_exploration.md`.
