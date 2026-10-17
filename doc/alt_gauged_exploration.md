# Alt-gauged exploration: the two-fermion S-matrix of the gauged bond layer

Exploration programme `alt_gauged` (merged into `master`), following
`alt_interactions`
([alt_interactions_exploration.md](alt_interactions_exploration.md)),
`alt-dynamics-exploration` ([alt_dynamics_exploration.md](alt_dynamics_exploration.md))
and `alt-cell-exploration` ([alt_cell_exploration.md](alt_cell_exploration.md)).
Experiments V175-V179, category `experiments/alt_gauged/`, tier 4,
0 new continuous parameters (the gauged coupling is eta = g_10 = 1/5,
the root angle - the only strength the construction admits).

## The question

Four exhaustive censuses have failed to isotropize the 12 M-node cones
of the master automaton (A = 0.998; V118: one tangent axis exactly
flat): signature (V149-V152), cell (V153-V157), motion (V161-V167), and
the minimal admissible collision (V171-V174: the block-contact vertex is
an O_h scalar, its one-loop self-energy is marginal and not
grid-convergent, A(b) = 0.9980 exactly under coarse-graining, zero
candidates). The alt-interactions report names the only surviving route:
the interacting **gauged sector** (V79-V101, Spin(10) adjoint quantum
links), whose two-fermion S-matrix is not built in the repo. This branch
builds it and measures what it can do to the cone.

## What the gauged sector already was (ground facts, cited)

| experiment | established fact |
|---|---|
| V79 | adjoint Spin(10) quantum links per bond (C^45); Gauss law C^A = T^A + L_left - L_right; bond gate K = exp(i g_10 J.L); content theorem 32 (x) 45 = 32 + 288 + 400 + 720, NO singlet (triality): the physical fermion is the DRESSED 32 |
| V98 | NO one-body physical states (center obstruction: every site carries an even number of spinors) - two fermions are the minimal physical system |
| V99 | three-site, two-link chain with paired matter: physical kernel dim ker = 24, fusion channels r = 1, 45, 210 |
| V101 | bond-current H = sum_b J.L diagonalized exactly on the 24-state kernel (Gram rank 24/24): E0 = -5802.3982 eta, degeneracy 6, gap 42.0906 eta |
| V129/V130 | the gauged layer is MASSIVE (Gauss-Laplacian gap 64.0000 = C2(45)) and CONFINING (string tension sigma = C2(45) at the Planck scale) |
| V140 | directed-link Gauss closes (conjugate head convention); KS electric gate exactly invariant; dressed hopping obstructed on the fixed adjoint link; the quintet directional layer is ISOTROPIC (spread 4e-8) but lives in the directional link layer, not the Weyl cone |
| V141 | the regular link of the root's finite group closes hopping but reduces the gauge group to a finite monomial group (order 1024, Spin(10) centralizer 0) and gaps the Eg half of the quintet at lambda^2 - computed no-go |

The missing object was the S-matrix. `cl10/smatrix.py` builds only the
free one (S = U (x) U, V48 Kronecker property).

## V175 - the S-matrix exists (the named missing object, built)

On the minimal physical two-fermion system (V99 kernel, 24 states, Gram
rank 24/24) the bond-current interaction exponentiates exactly:

```
S = exp(-i g_10 L^-1 H L^-dagger),   N = L L^dagger (exact Gram)
```

| quantity | value |
|---|---|
| unitarity ‖S S† - I‖ | 2.4e-13 (route 2: 2.3e-13) |
| spectrum vs V101 JSON | 7.8e-16 |
| 2-route (einsum/per-generator vs explicit-sum/per-bond; Cholesky vs orthogonal whitening) | H dev 2.7e-14, eigenphase-set dev 5.8e-11 |
| E0 / degeneracy / gap | -5802.3982 eta / 6 / 42.0906 eta |
| ground scattering phase (6-fold) | -1.909644 rad (= -g_10 E0 mod 2pi) |
| fusion-channel mixing (r = 1, 45, 210) | 2.8e-30, 6.1e-29, 6.1e-29 - the interaction PRESERVES the channels |
| eta -> 0 anchor | S = I exactly |

The S-matrix is a pure phase diagonal in the V101 energy eigenbasis:
the bond-current interaction is elastic and channel-blind on this
system (H is a gauge scalar). This is the exact two-fermion S-matrix of
CL10's only interaction.

## V176 - the vertex is NOT a scalar (the discriminating measurement)

The momentum-space image of the interaction is the Spin(10) current
Gamma^A = I_8 (x) T^A_16 on the 128-D Bloch space (T^A_16 the Weyl
generators, anchored: [T,T] = i f T to 8.9e-16, Casimir 45 I, and
plane_product(T16, g_10) = U16 exactly - the generators exponentiate to
the implemented step). The node-mode current form factor
F(kappa n) = sum_A |<u(q*)| Gamma^A |u(q* + kappa n)>|^2 over the 26 O_h
directions:

| \|kappa\| | A1g | Eg | T2g | F_mean |
|---|---|---|---|---|
| 0.05 | 0.8949 | 0.0072 | 0.0435 | 4.422 |
| 0.10 | 0.6618 | 0.0294 | 0.0588 | 3.267 |
| 0.20 | 0.6618 | 0.0294 | 0.0588 | 3.261 |
| 0.40 | 0.1667 | 0.3333 | 0.5000 | 0.381 |

Contrast with V171: the contact vertex is A1g = 1.0000 (scalar) with
anisotropic weight ~ k^2.00 (irrelevant). The gauged vertex carries the
spin-2 pieces (Eg + T2g = 5/6 at |kappa| = 0.4 - exactly the adjoint
flux's sym^2 content A1g + Eg + T2g measured in V140), and its
anisotropic part is scale-INDEPENDENT (alpha_aniso = 0.01 ~ 0): the
gauged vertex is MARGINAL on the cone, not irrelevant. The on-node
current expectation is F(0) = 5.000000 (exact). Anchors: Casimir
sum_A Gamma^A Gamma^A = 45 I_128 to 7.1e-15; 2-route (16-D vs 32-D
module) dev 8.9e-16.

## V177 - one-loop flow with the gauged vertex

The V172 protocol (static Fock self-energy, normal-ordered tadpole,
all 12 M-node singularities excluded and ray-subtracted) with
`current_vertex=True`, U_c = g_10, grids 16 and 24:

| grid | A_dressed(1e-2) | A_dressed(1e-3) | A_dressed(1e-4) | alpha_aniso(Sigma) |
|---|---|---|---|---|
| 16 | 0.9975 | 0.9969 | 0.9998 | 1.00 |
| 24 | 0.9985 | 0.9748 | 0.9656 | 1.00 |

The gauged self-energy is NOT zero on the cone (max |Sigma - tadpole| =
0.199 - the vertex reaches the cone, discriminating check), and its
anisotropic part scales as |k|^1.00: exactly the MARGINAL power of the
contact self-energy (V172), even though the VERTEX itself is marginal
(alpha = 0, V176) - the 1/omega propagator, not the vertex, sets the
power. The dressed anisotropy is NOT grid-convergent (drift 2.2e-2 on
A_dressed(1e-3)): the flat tangent axis (V118) makes the static
propagator non-integrable, as for the contact. The grid-24 drift toward
0.9656 at dq = 1e-4 is an unconverged measurement and cannot nominate
(V179 exclusion). Bare anchor A(1e-3) = 0.9980 (V63/V115).

## V178 - block RG with the gauged insertion

Same protocol as V173 (energy-window downfold of H(q) = i log U(q) onto
the node subspace, cuts pi/2, pi/4, pi/8 = b = 1, 2, 4, dq = 1e-3), with
the gauged self-energy (grid 16, normal-ordered) in the cone mode:

| b | n_low | window margin | A_bare | A_gauged |
|---|---|---|---|---|
| 1 | 71 | 3.9e-2 | 0.997985 | 0.997985 |
| 2 | 39 | 1.0e-4 | 0.997985 | 0.997985 |
| 4 | 17 | 2.7e-3 | 0.997985 | 0.997985 |

A(b) is EXACTLY unchanged by the gauged insertion at every scale (the
self-energy shifts the cone speeds, not their anisotropy ratio): the
cubic anisotropy stays marginal under coarse-graining with the gauged
vertex too. 2-route (b=1 downfold vs full gap_at) dev 1.3e-12; the b=1
bare row reproduces A = 0.998 (anchor).

## V179 - the gate

Scans the 12 census rows of V175-V178. The 2 V177 rows are excluded
from candidacy (grid drift 2.2e-2 >= 1e-6: an unconverged measurement
cannot nominate); minimum anisotropy over the remaining rows 0.9656;
ZERO candidates (no A < 1e-6, no converged flow toward isotropy > 0.5).
No hard-tier run was triggered; the gate's defensive branch (bond-layer
S unitarity 2.4e-13 + free front speed exactly 1) stands ready.

## Verdict

The two-fermion S-matrix of the gauged layer NOW EXISTS: exact, 24x24,
unitary to 1e-13, on the minimal physical system (V175). Its
momentum-space vertex is not the scalar the contact census used: it
carries the Eg + T2g (spin-2) pieces that can rotate the cone, and it is
marginal (alpha = 0) rather than irrelevant (alpha = 2) on the cone
(V176). Yet at one loop and under block RG it does not move the W9
anisotropy: A(b) = 0.997985 exactly with the insertion (V177/V178), and
the gate finds zero candidates (V179). This is consistent with what the
gauged sector is known to be: the flux sector is MASSIVE (gap 64 eta =
C2(45), V129) and CONFINING (V130), so it decouples from the IR cone as
a light mode - the interacting fixed point, if it exists, is not
reachable at one loop or at tree-level coarse-graining. The W9 isotropy
question remains open, now with the S-matrix built and the last
perturbative route measured.
