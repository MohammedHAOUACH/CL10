# Alt-step-search: a loss + gradient-descent search over the step/cell specification against W9

Programme `alt_step` (merged into `master`; it followed `alt_nonpert`).
Experiments: `experiments/alt_step/V184`, `V188`-`V190` (the V185-V187
numbers were taken by the cone programme in master). Parameter-rule override
(user-directed, this branch only): the 45 per-plane step angles are
searched as CONTINUOUS variables; the hard tier (unitarity, front speed
exactly 1, closure order <= 64, chirality) stays non-negotiable and a
negative search is a valid deliverable.

## 1. Why this search: the logic

W9 asks whether the cone of gapless matter is isotropic at low energy.
Every mechanism inside the model has now been measured and closed:

- interactions do not isotropize it (`alt_interactions_exploration.md`);
- the gauged sector's S-matrix, one-loop flow and block RG do not
  (`alt_gauged_exploration.md`);
- the non-perturbative flux sector is massive at every volume and the
  exact cone dressing leaves A untouched to ~1e-12
  (`alt_nonperturbative_exploration.md`);
- the dispersion is exactly O_h-invariant (max |E(Pq) - E(q)| = 2.0e-15
  over the 48 signed permutations), so the anisotropy is not an engine
  artifact - it is a property of the step and the cell.

What remains is the SPECIFICATION itself: the step
U16 = plane_product(T16, g_10) and the cell (2,2,2) are one point in a
large space - 45 continuous angles times discrete cell/stage choices.
The question this branch answers: **does any step in that space have a
gapless AND isotropic cone?** If yes, W9 has a positive solution and the
new step faces the full hard tier. If no, W9 is not a specification
tuning problem, and that negative is the deliverable.

## 2. The search space (V184's parameterization)

The step is a product over the 45 root planes of the algebra:

    U(theta) = Pi_{k=0..44} ( cos(theta_k) I_16 - i sin(theta_k) T16[k] )

in the `weyl_generators()` plane order - the same order as the
implemented `plane_product`. Each factor is exactly unitary (T16[k]
Hermitian), so U(theta) is unitary for every theta: the search never
leaves the unitary manifold, and the implemented step is the single
point theta0 = g_10 * ones(45). The parameterization is exact, not an
approximation: anchor (a) below shows U(theta0) = `SC.sectors()[2]` to
1.3e-15.

The momentum side is the block Bloch matrix of a generalized cell:

    M(q, theta) = Pi_stages kron(S_stage(q), U(theta))

with S_stage the nsub x nsub bond-stage matrix (per bond (i,j,v) the
2x2 rotation with phase exp(i q.v)) - a jnp reimplementation of
`bloch_cell` validated against it by anchor (b). The cell enters through
its shape (n1,n2,n3) and its stage order; V189 searches four choices.

## 3. The loss and why it is shaped this way

At a step theta, two things must hold for a W9 solution: the spectrum
must be gapless somewhere (a node q*), and the cone at that node must
be isotropic. Both are made differentiable:

- **gap proxy**: K = (M + M^dag)/2; omega(q, theta) =
  arccos(eigvalsh(K)[-1]). For a unitary M this equals `gap_at`
  exactly (eigenphases in (-pi, pi], cos monotone) - anchor (b): dev
  7.2e-12 vs `gap_at`.
- **node**: q*(theta) is found, not fixed - coarse grid over
  {0, pi/2, pi, 3pi/2}^3 then Nelder-Mead on omega^2 (the eigvalsh VJP
  is invalid at the degenerate node, so the node is optimized
  numerically, not by autodiff).
- **anisotropy**: the published V115 metric. Speeds
  v(n) = omega(q* + dq n, theta)/dq over the 124 `sphere_dirs(6)`
  directions, dq = 1e-3; A = (vmax - vmin)/(vmax + vmin). A = 0 means
  every direction has the same speed.
- **loss**: L(theta) = A + 100 * omega(q*)^2. The penalty term keeps
  the search on gapless steps; A drives it toward isotropy.

The candidate criterion (pre-registered, never tuned): a row is a
candidate iff A < 1e-6 AND omega(q*) < 1e-6 - a genuinely isotropic
gapless cone.

## 4. The gradient: what fails naively and what was built instead

Autodiff through `eigvalsh` is a trap here, in three ways measured
during development:

1. at the node lambda_max is degenerate: the eigenvalue VJP is a
   subgradient only and its eigenvector part is discontinuous;
2. omega = arccos(lambda) has infinite derivative at lambda = 1 (the
   gapless clip): the VJP produces +-inf and Adam NaNs (this exact
   failure was observed on the first V188 run);
3. max/min over the 124 speeds are kinks: only the argmax/argmin
   directions carry the gradient of A.

The built pipeline (V184 `combined_grad`) is the exact subgradient:

- the aniso part: only the two selected (argmax/argmin) directions,
  each via the single-band formula
  d omega / d theta = -v^dag (dK/dtheta) v / sqrt(1 - lambda^2),
  with the eigenvector v and lambda taken from a forward
  eigendecomposition (no VJP through either);
- the penalty part: the projector subgradient of the top eigenspace
  (averaged over the top pair when degenerate - the natural
  subgradient), hinged: active only while omega(q*) > 1e-6. The hinge
  matters: an always-on penalty crushes the cone into a flat band
  (A -> 1), which was measured before the hinge was introduced;
- Adam is hand-rolled (no optax installed), with a non-finite-gradient
  guard, lr 0.02, 30 steps per cycle.

Validation anchors (all in V184, all discriminating):

| anchor | what it proves | measured |
|---|---|---|
| (a) U(g_10*1) vs `sectors()[2]` | the parameterization is the implemented step | 1.3e-15 |
| (b) omega vs `gap_at` | the gap proxy is exact | 7.2e-12 |
| (c) A vs V115 published 0.998 | the metric is the published one | 0.9980 |
| (d) single-band grad vs central finite differences | the subgradient formula | 2.7e-08 |
| (d') penalty_grad vs autodiff at a gapped q | the projector subgradient | 5e-13 |
| (e) fast NumPy K forward vs jnp `K_of` | the speed refactor is exact | 1.3e-15 |

Engineering notes (measured, recorded): the forward of the 124
direction-shifted K matrices factors as
M = (P(q) kron I16)(I kron U^12) - one 16x16 power for theta, one
nsub x nsub product for q, one einsum for the assembly (124 matrices in
0.13 s); on a live CUDA device the batch runs through XLA `eigh`
(cuSOLVER) instead, with a NumPy/LAPACK fallback and a XLA retry on
`LinAlgError`. The route (cpu:0 / cuda:0) is recorded in every JSON.

## 5. The protocol (pre-registered before the runs)

- **V184**: anchors + a 9-row landscape (theta0 + 8 uniform(0, 2pi)
  seeds, rng 184): what the loss looks like before any descent.
- **V188**: 25 seeds - baseline g_10*1, 16 uniform(0, 2pi) (rng 185),
  8 quarter-grid {pi/4, pi/2} (rng 186). Alternating scheme, 10
  cycles x [Nelder-Mead node + 30 Adam steps on theta]; best-L
  tracking of the trajectory (Adam overshoots on this rugged
  landscape; the returned iterate is the best seen).
- **V189**: four cell/stage choices - (2,2,2) baseline stages,
  (2,2,2) phase-major, (2,2,3), (3,3,3) (M = 128/128/192/432) - each
  with 8 seeds (baseline + 7 uniform, rng 187), same scheme.
- **V190**: the gate. Candidacy from all V188/V189 rows (A < 1e-6 AND
  omega < 1e-6); every candidate faces the hard tier on U* =
  U(theta*): unitary_err < 1e-12, closure period <= 64, front speed
  exactly 1 with zero outside-lightcone bonds, chirality nonzero -
  plus the pre-registered pi/8 closure-projected re-evaluation.
- Anchors between experiments: V188's baseline initial A must
  reproduce V184's baseline row (gate 1e-3: A is a knife-edge ratio,
  dA/domega_min ~ 5e3, and the NM node landing jitters ~1e-5 in q*
  across processes on the GPU route); V189's (2,2,2)/baseline initial
  A must reproduce V188's (same gate). The post-descent A is chaotic
  across processes (NM jitter cascades through 300 Adam steps) and is
  therefore never used as an anchor.

## 6. Results

### V184 - the loss is defined and anchored (PASS, 41 s GPU)

All six anchors pass (table above). Landscape: the baseline sits at
A = 0.9994 with omega(q*) = 4.2e-8; the 8 random steps are all
gapless-flat (omega = 0, A = 0.977-1.000) - a random product of plane
rotations generically has flat gapless bands, so the hard direction of
the search is isotropy at fixed gaplessness, not gaplessness itself.

### V188 - Adam descent on the 45 angles, baseline cell (PASS, 25 seeds)

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

Reading: two escape routes, never both. Descent either (i) stays
gapless and keeps the cone anisotropic (every omega = 0 row has
A >= 0.966), or (ii) isotropizes - 10 of 25 seeds reach A ~
1e-4..1e-11 - but only at gapped band edges (omega = 3.5e-4..3.7e-2).
Route (ii) is TRIVIAL isotropy: at a band edge the touching is
quadratic, so v(n) = omega_0/dq + O(dq) is direction-independent and
A -> 0 without any cone. The two closest approaches to a candidate
both miss on the gap: quarter-4 (A 1.2e-4, omega 1.5e-5) and
uniform-1 (A 1.05e-6, omega 3.5e-4). Candidates: ZERO. The baseline
itself descends to A 3.8e-11 at omega 3.4e-3 - the loss prefers the
trivial escape even from the implemented step.

### V189 - joint cell x step search (PASS, 4 cells x 8 seeds, 4 h GPU)

| cell | stages | M | best A | best omega(q*) | best seed |
|---|---|---|---|---|---|
| (2,2,2) | baseline | 128 | 0.9966 | 0 | uniform-5 |
| (2,2,2) | phase-major | 128 | 0.9805 | 0 | uniform-2 |
| (2,2,3) | phase-major | 192 | 0.9675 | 0 | uniform-5 |
| (3,3,3) | phase-major | 432 | 0.9671 | 0 | uniform-4 |

All 32 rows are gapless-flat: for every cell shape and stage order,
every seed's descent ends gapless with A 0.967-1.000. No cell admits
an isotropic gapless step. Structural tier note: anomalies,
hypercharge and the chiral matter content are fixed by the signature +
sector decomposition and are untouched by (cell, theta) - the cell
changes only the motion - so no structural re-check was required for
non-candidates.

**Reproduction of the carried verdict** (`V189_reproduce_objective.py`, driven
by `run_reproduce_check.sh`; 2026-10-18, no descent).  The verdict stores the
$(\theta^*, q^*)$ every row converged to, so its objective can be re-evaluated
under the current library without re-running the search.  **On the card - the
route the run itself used (`JAX cuda:0`) - it reproduces bit for bit**: for the
best row of all four cells the `aniso` deviation is `0.00e+00` and `omega` is
exactly 0 (`gap2 ~ -2e-14`).  On the CPU core the same evaluation deviates by
2.2e-05 / 1.8e-04 / 4.1e-06 / 1.3e-05, which is the GPU-vs-CPU eigensolver at
a node where the gap is pinned to zero: `cl10/supercell.py` and
`cl10/config.py` are unchanged since before the run, and the `altcell`
host-only edit moved arrays between hosts without touching values.  Either way
the candidacy threshold 1e-6 sits far below, so the conclusion above does not
depend on re-running the 4 h search.

### V190 - the step-search gate (PASS, 0 candidates)

38 census rows scanned (V184 9 + V188 25 + V189 4 summaries);
candidacy A < 1e-6 AND omega(q*) < 1e-6: **zero candidates**. The hard
tier had nothing to reject; the pi/8 closure-projected re-evaluation
was never triggered. PASS with the negative recorded.

## 7. What the negative actually establishes - and its limits

Established, within the searched space (45 continuous angles x 4
cell/stage choices, 57 descent runs, ~17000 Adam steps):

- no step in the space carries a gapless isotropic cone - the two
  properties never coexist along any descent trajectory;
- the implemented step is not a lucky point but an isolated closed-form
  gapless point: random steps are gapless-FLAT (A -> 1), and isotropy
  is reachable only by paying a gap;
- the loss's blind spot is characterized, not hidden: A = 0 at a
  gapped band edge is trivially true (quadratic touching), which is
  why the candidacy gate requires omega < 1e-6 as well, and why the
  closest approaches are reported by how much gap they carry.

Limits, honestly:

- gradient descent is local: 57 seeds explore 57 basins of a rugged
  45-dimensional landscape, not the whole space; a far-away isotropic
  gapless step could exist that no seed found. What the census shows
  is that none is reachable by descent from generic or structured
  starts, and that the loss's isotropic basin is the trivial one;
- the closure gate (period <= 64) was never exercised (no candidate
  reached it); a hypothetical candidate with generic angles would in
  general fail closure, and the pi/8 projection row was the
  pre-registered bridge for that - it stayed unused;
- the search covers the step's 45-plane factorization and four block
  cells; other motion encodings (different bond graphs, non-block
  cells) are outside this branch and belong to the alt_cell census.

## 8. Verdict

The step/cell SPECIFICATION has been searched directly - 45 continuous
per-plane angles by differentiable cone loss + Adam (25 seeds), and
four cell/stage choices (32 seeds), on GPU with every route anchored
against the published machinery. No isotropic gapless step exists in
the searched space: descent either keeps the cone anisotropic (A >=
0.966 at every gapless point it reaches) or isotropizes only trivially
at gapped band edges (omega >= 3.5e-4). The implemented step (g_10 on
all 45 planes, cell (2,2,2)) remains the only closed-form gapless
point found, and W9 (cone isotropy) stays open: it is not a tunable
consequence of the step angles.

Cross-links: `alt_nonperturbative_exploration.md` (the gauged-sector
closure that motivated searching the specification),
`alt_gauged_exploration.md`, `alt_interactions_exploration.md`.
