# Alt-interactions exploration: do collisions wash out the cone anisotropy?

Exploration programme `alt_interactions` (merged into `master`),
following `alt-dynamics-exploration`
([alt_dynamics_exploration.md](alt_dynamics_exploration.md)) and
`alt-cell-exploration` ([alt_cell_exploration.md](alt_cell_exploration.md)).
Experiments V171-V174, category `experiments/alt_interactions/`, tier 4,
0 new continuous parameters (collision strength from the pre-registered
set {g_10, pi/4, pi/2}).

## The question

The node census (V63) and every free generalization - signature
(V149-V152), cell (V153-V157), motion (V161-V167), the Wilson family
(V115), bond-layer generalizations (V120/V121), frame averaging (V142) -
leave the 12 M-node cones anisotropic (A = 0.998; V118: one tangent axis
exactly flat). The one route none of those censuses touched is the
standard lattice-statistics miracle: at an INTERACTING fixed point the
cubic anisotropy is irrelevant and rotational invariance is restored in
the continuum (the 3D Ising / lattice-phi^4 case; Aharony: irrelevant for
N < N_c ~ 3). This branch tests that route computationally.

## The free-theory fact

The master automaton has NO collision vertex. `cl10/smatrix.py` builds
the 2-particle S-matrix as the Kronecker product ``U (x) U``, and V48
confirms every n-particle amplitude is a Kronecker product ("the
Kronecker product property"). The bond-layer swap-rotation entangles two
fermions sharing a cell (V97, S_A up to 0.783) but a one-body gate never
scatters - it is not an interaction. The interacting sector of CL10 is
the GAUGED bond layer (V79-V101: Spin(10) adjoint quantum links, Gauss
kernel, matter-flux entangled physical sector); its two-fermion S-matrix
is not built in the repo. The branch therefore uses the minimal
admissible collision, clearly labelled exploration-only: the lattice-gas
BLOCK CONTACT ``C = exp(i U_c N_block)`` (a pair of fermions in the same
2x2x2 block picks up the phase ``exp(i U_c)``), ``U_c`` from the
pre-registered set; ``U_c = 0`` is the free anchor.

## The cone subspace

At an M node the gapless branch is a SINGLE mode (|theta| = 2.0029e-4,
next band 0.0147 - V63), so the cone projector is rank-1 and
well-isolated; all V171-V173 kernels use it. The 4-mode window first
tried is NOT isolated (bands at 0.0717/0.0747 are quasi-degenerate) and
was rejected: the rank-1 projector agrees between independent solvers
(eigh of i log U vs eig of U + QR) to 1.3e-14 (V171 2-route).

## The collision census

### V171 - the vertex (experiments/alt_interactions/V171_collision_vertex)

The contact vertex is O_h-scalar (weight = U_c, k-independent). The
node-mode form factor ``f(n) = |<u_1(q* + k n)|u_1(q*)>|^2`` - the only
place the lattice geometry enters the vertex - is A1g-dominant
(w_A1g = 1.0000 at k = 1e-3) and in the IR regime (k below the 0.0147
band gap) both its anisotropic weight and its isotropic leakage scale as

    alpha_aniso = alpha_iso = 2.00   (all three node representatives)

The collision vertex is IRRELEVANT on the cone: the interaction adds no
new anisotropy at low energy - but a scalar vertex also cannot remove
the anisotropy that lives in the propagator. 12 rows, free anchor
A = 0.998, 2-route 1.3e-14.

### V172 - the one-loop flow (V172_self_energy_flow)

The static Fock (exchange) self-energy of the contact on the cone,
normal-ordered (the k-independent tadpole absorbed in the node) with all
12 M-node singularities excluded and ray-subtracted (the pre-registered
contingency for the flat-axis IR sensitivity). Result: the anisotropic
part of Sigma scales as |k|^1.00 - the SAME power as the bare cone, so
the contact self-energy is MARGINAL, not irrelevant. The dressed
anisotropy is NOT grid-convergent (grids 24 vs 32 disagree by
0.019-0.66 on A_dressed(1e-3)): the exactly flat tangent axis (V118)
makes the static propagator non-integrable and the one-loop coefficient
grows with resolution. The one-loop flow is INCONCLUSIVE at this order -
and the flat direction is IR-singular against the collision. 3 rows,
anchor A_bare = 0.9980.

### V173 - the block RG (V173_block_rg)

Loop-free coarse-graining: downfold H(q) = i log U(q) onto the node
subspace at energy windows cut in {pi/2, pi/4, pi/8} (block scales
b = 1, 2, 4; subspace sizes n_low = 70, 38, 18):

| b | cut | A_bare | A_interacting (U_c = g_10) |
|---|-----|--------|----------------------------|
| 1 | pi/2 | 0.9980 | 0.9980 |
| 2 | pi/4 | 0.9980 | 0.9980 |
| 4 | pi/8 | 0.9980 | 0.9980 |

The cubic anisotropy is MARGINAL under coarse-graining - it neither
grows nor shrinks - and the self-energy insertion does not change that
at this order. The downfold is exact on the low band (b = 1 reproduces
the full gap_at cone to 3.5e-13), so this is a statement about the
theory, not the truncation. (Window margin at b = 2 is 7.1e-5: a band
sits near the cut; recorded in the verdict, the cone branch itself is
deep inside the subspace.)

### V174 - the gate (V174_isotropy_gate)

18 census rows scanned; the 3 unconverged V172 rows excluded from
candidacy (an unconverged measurement cannot nominate); minimum
converged anisotropy 0.9980; ZERO candidates. No hard-tier run needed:
the minimal admissible collision does not isotropize the W9 cone.

## Verdict

At every level computable with the minimal admissible collision:

1. the collision VERTEX is irrelevant (alpha = 2, V171) - the premise
   of the miracle holds: the interaction does not add anisotropy;
2. the one-loop self-energy is MARGINAL (alpha = 1, V172) and its
   coefficient is not grid-convergent - the flat tangent axis is
   IR-singular against the collision;
3. the cone anisotropy itself is MARGINAL under coarse-graining
   (A(b) = 0.9980 exactly, V173) - the free value is a fixed point of
   the loop-free flow.

The Ising/Wilson-Fisher miracle is a statement about the INTERACTING
fixed point - a non-perturbative object. One-loop irrelevance of the
VERTEX (measured) is NOT restoration of the CONE (not achieved): the
restoration would need the fixed point of the actual CL10 interaction,
i.e. the gauged bond layer's 2-fermion S-matrix (V79-V101), which is
not built. W9's anisotropy survives every interaction test admissible
without new continuous parameters; the surviving route to isotropy is
the non-perturbative gauged-layer fixed point, not the automaton-level
collision.

## GPU route

The heavy grid pass (batched 128x128 Bloch eigendecomposition) has a
numba CUDA implementation in `cl10/altgpu.py` (batched cyclic Jacobi
for normal matrices: each 2x2 core is diagonalized by the unitary
polar factor G = X(X^dagger X)^{-1/2} of its eigenvector matrix,
applied as U <- G^dagger U G; gtc2017-numba style, one thread block
per matrix). `run_altgpu.sh` verifies it: CONFIRMED - eigenphase error
1.2e-11 and projector error 2.8e-15 vs numpy eig on random unitaries,
and the Fock self-energy on the grid-16 pass agrees GPU vs numpy on
the actual Bloch matrices to 2.9e-13. The verification exposed a
structural fact: the automaton's anti-unitary symmetry pairs the
spectrum in +-theta, so at 1315/4095 grid points the two smallest
|theta| coincide exactly and the single-mode projector of V172 is an
arbitrary basis choice there (the subspace projector is canonical);
the comparison uses the subspace convention at degenerate points. It
is a verification route, not a verdict; the CPU verdicts above are
the record.

## Boundary

- F7 (why this g_10 and this root): untouched - the collision census
  changes no coupling.
- O3 (why three families): untouched.
- W9 (cone isotropy): the interaction route is now also censused and
  closed at every perturbative order computable at the automaton level;
  only the non-perturbative gauged-layer fixed point remains open.
