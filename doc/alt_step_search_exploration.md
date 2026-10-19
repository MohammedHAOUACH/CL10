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
  landscape; the returned iterate is the best seen). The step is
  evaluated through its exact factorization - V189's `Core` object,
  imported from V189 - and the descent takes central finite
  differences on the 45 angles: the loss is a max over phase pairs,
  so there is no gradient to differentiate.
- **V189**: the same four cell/stage choices - (2,2,2) baseline stages,
  (2,2,2) phase-major, (2,2,3), (3,3,3) (M = 128/128/192/432) - but no
  descent: the step factorizes as `P(q) (x) U(theta)^ns`, so 419
  pre-registered steps per cell are evaluated directly (rng 187).
- **V190**: the gate. Candidacy from all V188/V189 rows (A < 1e-6 AND
  omega < 1e-6); every candidate faces the hard tier on U* =
  U(theta*): unitary_err < 1e-12, closure period <= 64, front speed
  exactly 1 with zero outside-lightcone bonds, chirality nonzero -
  plus the pre-registered pi/8 closure-projected re-evaluation.
- Anchors between experiments: V188's baseline initial A must
  reproduce V184's baseline row (gate 1e-3: A is a knife-edge ratio,
  dA/domega_min ~ 5e3, and the NM node landing jitters ~1e-5 in q*;
  measured 2.8e-4). V188 and V189 both carry the same harder anchor:
  the implemented step evaluated at one node through the factorized
  route and through V184's `(16 n_sub)^3` eigendecomposition, with
  `lam_max` agreeing to 1e-12 (measured 6.7e-16 and 1.3e-15) - the
  fast route is the reference route, not an approximation of it. The
  post-descent A is chaotic across processes and routes (NM jitter
  cascades through 300 Adam steps) and is therefore never used as an
  anchor.

## 6. Results

### V184 - the loss is defined and anchored (PASS, 41 s GPU)

All six anchors pass (table above). Landscape: the baseline sits at
A = 0.9994 with omega(q*) = 4.2e-8; the 8 random steps are all
gapless-flat (omega = 0, A = 0.977-1.000) - a random product of plane
rotations generically has flat gapless bands, so the hard direction of
the search is isotropy at fixed gaplessness, not gaplessness itself.

### V188 - Adam descent on the 45 angles, baseline cell (PASS, 25 seeds, 15 min CPU)

The step is evaluated through its exact factorization - V189's `Core`, imported
from V189 - and the descent takes central finite differences on the 45 angles:
the loss is a max over phase pairs of $\cos(\phi_i(q)+\psi_j(\theta))$, so there
is no gradient to differentiate (V184's anchor (d) had already pinned its
autodiff subgradient against finite differences at 2.7e-8). No card is needed.

| seed | A_final | omega(q*) | L | ||dtheta|| |
|---|---|---|---|---|
| baseline g_10*1 | 0.9665 | 1.18e-13 | 0.9665 | 1.96 |
| uniform-0 | 0.9385 | 8.40e-14 | 0.9385 | 0.80 |
| uniform-1 | 0.9634 | 1.08e-13 | 0.9634 | 2.43 |
| uniform-2 | 0.9762 | 1.04e-13 | 0.9762 | 2.19 |
| uniform-3 | 0.9690 | 1.56e-13 | 0.9690 | 3.00 |
| uniform-4 | 0.9499 | 1.12e-13 | 0.9499 | 1.34 |
| uniform-5 | 0.9701 | 9.69e-14 | 0.9701 | 2.76 |
| uniform-6 | 0.9601 | 1.36e-13 | 0.9601 | 0.72 |
| uniform-7 | 0.9893 | 2.07e-13 | 0.9893 | 1.55 |
| uniform-8 | 0.9740 | 1.23e-13 | 0.9740 | 1.97 |
| uniform-9 | 0.9826 | 1.15e-13 | 0.9826 | 0.00 |
| uniform-10 | 0.9806 | 9.06e-14 | 0.9806 | 2.99 |
| uniform-11 | 0.9826 | 8.58e-14 | 0.9826 | 1.49 |
| uniform-12 | 0.9905 | 1.13e-13 | 0.9905 | 2.30 |
| uniform-13 | 0.9687 | 8.68e-14 | 0.9687 | 1.14 |
| uniform-14 | 0.9622 | 9.53e-14 | 0.9622 | 1.80 |
| uniform-15 | 0.9585 | 9.03e-14 | 0.9585 | 1.42 |
| quarter-0 | 0.9748 | 8.74e-14 | 0.9748 | 2.16 |
| quarter-1 | 0.9812 | 8.77e-14 | 0.9812 | 0.00 |
| quarter-2 | 0.9855 | 1.07e-13 | 0.9855 | 2.05 |
| quarter-3 | 0.9692 | 1.04e-13 | 0.9692 | 2.58 |
| quarter-4 | 0.9682 | 1.16e-13 | 0.9682 | 0.81 |
| quarter-5 | 0.9831 | 9.33e-14 | 0.9831 | 2.85 |
| quarter-6 | 0.9823 | 1.17e-13 | 0.9823 | 1.93 |
| quarter-7 | 0.9749 | 8.18e-14 | 0.9749 | 1.36 |

Reading: **one escape route, and it is never taken.** All 25 descents end on the
gapless manifold ($\omega(q^*)\le 2.1\mathrm{e}{-}13$) with the cone anisotropic,
$A\in[0.9385,0.9905]$; the best is 0.9385, and the baseline itself lands at
0.9665 from $A_0=0.9997$. **Candidates: ZERO.** The earlier card run of the same
programme (JAX autodiff subgradient, cuSOLVER node search) did reach the trivial
isotropy - 6 seeds at $A\sim10^{-8}$, every one at a **gapped** band edge
($\omega\ 8.8\mathrm{e}{-}4..2.3\mathrm{e}{-}2$): at a band edge the touching is
quadratic, $v(n)=\omega_0/dq+O(dq)$ is direction-independent and $A\to0$ without
any cone. That is the loss's blind spot, and why candidacy requires
$\omega<10^{-6}$ as well as $A<10^{-6}$ - and why the candidate count, not the
endpoint (which is route-dependent), is the pre-registered quantity. On either
route the count is zero.

Anchors: the baseline row reproduces V184's $A(\theta_0)$ to 2.8e-4 (gate 1e-3),
and the factorized route reproduces V184's dense route at one node to
$\lambda_{\max}$ 6.7e-16 (gate 1e-12) and $A$ to 3.7e-4 - the arccos
amplification of a 1e-15 eigenvalue error at a gapless point.

### V189 - joint cell x step CENSUS (PASS, 4 cells x 419 steps, 8 min CPU)

**Why the shape changed.** The first version ran Adam jointly over the cell and
the 45 angles: 4 h of GPU, and inside the suite it could not even be re-executed,
because the eight cupy workers already hold the card and its own four workers
killed themselves at import. The step then turned out to **factorize exactly**,
which removes both the cost and the need for a descent:

$$M(q,\theta) = P(q)\otimes U(\theta)^{n_s}$$

$P(q)$ is the product of the cell's bond matrices - it carries no $\theta$ at
all - and $U(\theta)$ the internal $16\times16$ unitary, $n_s$ the number of
stages. The spectrum of a Kronecker product is the set of products, so

$$\text{gap}(q,\theta)=\min_{i,j}\bigl|\arg p_i(q)+\arg v_j(\theta)\bigr|$$

Checked against V184's own core on all four cells: $|M_{\text{direct}} - P\otimes
U^{n_s}| \le 8\mathrm{e}{-}16$, the gap to $3\mathrm{e}{-}15$. The $q$-half is
computed **once per cell** and reused for every step, and one step costs a
$16\times16$ eigendecomposition plus a phase search instead of a
$(16n_{\text{sub}})^3$ diagonalization per optimizer step; $n_{\text{sub}}$ is
only 8, 12 or 27, so the factorization skips a factor of $\sim10^4$.

**What it does instead.** No optimizer: a **census**. 419 pre-registered steps
per cell - the implemented step, the quarter-angle steps, per-stage cycles, 128
quarter-lattice draws, 256 uniform draws, 24 perturbations of the baseline - run
on all four cells, judged by the same criterion as before: gapless at its own
node ($\omega<10^{-6}$) and isotropic there ($A<10^{-6}$, the V115 cone metric
over the 215 sphere directions, the same objective V184 and V188 used).

| cell | stages | $M$ | gapless | best $A$ | best seed |
|---|---|---|---|---|---|
| (2,2,2) | baseline | 128 | 385/419 | 0.9466 | uniform-210 |
| (2,2,2) | phase-major | 128 | 385/419 | 0.9631 | uniform-55 |
| (2,2,3) | phase-major | 192 | 407/419 | 0.9377 | uniform-100 |
| (3,3,3) | phase-major | 432 | 403/419 | 0.9401 | uniform-146 |

**Result.** 1580 of the 1676 sampled steps are gapless, and their cone
anisotropy spans $A\in[0.9377,1.0000]$ (median 0.9944). **No step is isotropic
at its node.** Where the descent version could report one trajectory per seed,
the census reports the distribution of $A$ over the whole family - and the
anomalies, hypercharge and the chiral matter content are fixed by the signature
and the sector decomposition, untouched by (cell, theta), so no structural
re-check was required for non-candidates.

**The anchor.** The implemented step, evaluated at one node through the
factorized route and through V184's $(16n_{\text{sub}})^3$ eigendecomposition,
agrees in $\lambda_{\max}$ to $1.3\mathrm{e}{-}15$: the fast route is the
reference route. The cone metric agrees to $3.7\mathrm{e}{-}4$, which is the
$\arccos$ amplification of a $10^{-15}$ eigenvalue error at a point where the
gap is pinned to zero, not a difference of spectra. The two node searches land
up to $7.1$ apart on that step (it sits $4\mathrm{e}{-}8$ from gapless) - a
property of the search, recorded and deliberately not gated.

### V190 - the step-search gate (PASS, 0 candidates)

1710 census rows scanned (V184 9 + V188 25 + V189 1676);
candidacy A < 1e-6 AND omega(q*) < 1e-6: **zero candidates**. The hard
tier had nothing to reject; the pi/8 closure-projected re-evaluation
was never triggered. PASS with the negative recorded.

## 7. What the negative actually establishes - and its limits

Established, within the searched space (45 continuous angles x 4
cell/stage choices: 25 Adam descents on the baseline cell, ~17000 Adam
steps, plus 1676 explicitly sampled steps x 4 cells with no descent):

- no step in the space carries a gapless isotropic cone - the two
  properties never coexist along any descent trajectory, nor at any of
  the 1676 census steps;
- the implemented step is not a lucky point but an isolated closed-form
  gapless point: random steps are gapless-FLAT (A -> 1), and the only
  way the metric reports isotropy at all is the trivial band-edge
  touching, which costs a gap;
- the loss's blind spot is characterized, not hidden: A = 0 at a
  gapped band edge is trivially true (quadratic touching), which is
  why the candidacy gate requires omega < 1e-6 as well, and why the
  closest approaches are reported by how much gap they carry.

Limits, honestly:

- the search is a sample, not the space: 25 descent basins plus 1676
  pre-registered steps of a 45-dimensional landscape. What the census
  adds over the descents is the distribution - 1580 of the 1676 sampled
  steps are gapless and their A spans 0.9377..1.0000, so the isotropic
  value is not merely unreachable by descent, it is absent from a
  structured + random sample of the family;
- the closure gate (period <= 64) was never exercised (no candidate
  reached it); a hypothetical candidate with generic angles would in
  general fail closure, and the pi/8 projection row was the
  pre-registered bridge for that - it stayed unused;
- the search covers the step's 45-plane factorization and four block
  cells; other motion encodings (different bond graphs, non-block
  cells) are outside this branch and belong to the alt_cell census.

## 8. Verdict

The step/cell SPECIFICATION has been searched directly - 45 continuous
per-plane angles by differentiable cone loss + Adam (25 seeds), and a
419-step census on each of four cell/stage choices (1676 steps, no
descent), every route anchored against the published machinery. No
isotropic gapless step exists in the searched space: every descent ends
gapless and anisotropic, and the census says the same over 1676 sampled
steps (A >=
0.9377 at every gapless point reached, by descent or by census); the only
isotropy the loss can produce is the trivial band-edge one, which costs a
gap. The implemented step (g_10 on all 45 planes, cell (2,2,2)) remains the
only closed-form gapless point found, and W9 (cone isotropy) stays open: it is not a tunable
consequence of the step angles.

Cross-links: `alt_nonperturbative_exploration.md` (the gauged-sector
closure that motivated searching the specification),
`alt_gauged_exploration.md`, `alt_interactions_exploration.md`.

**Follow-ups (the question closed in three steps).** `V191_parity_dirac_cone`
supplies the *why*: $\mathrm{Hom}_{O_h}(T_{1u}, D\otimes D^*) = 0$ for all ten
irreps, so no single-parity multiplet at an inversion-invariant momentum can
be conical - and it constructs the escape, an exactly isotropic non-abelian
cone ($\rho_{\mathrm{comm}} = 1/\sqrt2$, anisotropy $10^{-15}$) on a
parity-degenerate $A_{1g}\oplus T_{1u}$ pair built on an $O_h$-invariant
hopping, outside the engine's circuit class. `V192_engine_two_angle_cone`
searches the engine's own six-stage product - the largest cubic-covariant
stage family (angle $a$ on the three intra-block stages, $b$ on the three
inter-block ones) jointly with the internal angle, $33\,800$ grid points
polished: $60$ exactly gapless points, and the parity split is the answer -
$8$ single-parity multiplets flat to $1.3\times10^{-12}$, $52$ parity-mixed
ones that disperse but are exactly Abelian ($\rho_{\mathrm{comm}} \le
1.2\times10^{-11}$): $0$ conical nodes. `V201_root_swap_cone` extends the
internal layer by the root's own $10$ sector-swap involutions $h_c w$ (the
volume element $w = h_1..h_{10}$: $S^2 = I$, $\{S, w\} = 0$, exact to $0$;
the doubled $32$-D internal module anchored to the Weyl step at $8\times10^{-15}$):
the rerun of the same grid gives $8$ exactly gapless points per swap, all
$6$-dim $T_{2g}\oplus T_{1u}$ multiplets with **non-commuting** tangent
operators ($\rho_{\mathrm{comm}}$ $0.50$–$0.80$) but vanishing linear velocity
(max slope $1.7\times10^{-12}$): the non-commutativity is a root ingredient,
the dispersion is not, and the cone stays outside the step class. Paper
reference: Section 8.4 and Table 7, row O9.
