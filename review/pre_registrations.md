# Pre-registrations: the falsifier-attack programme (V229–V238)

One entry per new experiment, written **before its first run** and never edited
afterwards: corrections are appended as numbered addenda. Each entry states the
question, the closed set of structures and inputs the experiment may use, the
numeric acceptance rule, and the two manuscript consequences — the wording if
the experiment closes the problem and the wording if it fires the falsifier.
Any input outside the closed list goes into the artifact's
`data.fitted_inputs` and moves the affected number's class in Table 6 of
`doc/final_publication.md` to `fitted`.

The verdict convention is the current generation (V209): a falsifier that
fires is a `pass: true` experiment with `residual: true` and a `branch` string
naming the fired outcome (as V128 carries "(**F4**)", V197 F3, V206
`CLOCKS TENSION`). Tiers follow `experiments/README.md:9-17`; the tier was
chosen by what the experiment tests, never by which tier is exempt.

Allowed inputs everywhere (ground rule 2): the `angle:` values of
`cl10/config.yaml` (`implemented: 0.2`, `closed: 0.7853981633974483`,
`half: 1.5707963267948966`, `alpha_10: 0.04`), `desert:`
(`m_x_gev: 1.8e+16`, `m_x_closed_gev: 1.235e+17`,
`theta_root_rad: 1.7307836435290485`), `cell:` (`rung_A: 3`, `rung_B: 10`),
the axes of `cl10.starcell.star_axes()`, the V65 selective Wilson root, and
the committed generators in `experiments/.../cone_generator.json`. Plus the
`bounds:` keys named per entry below. All GPU work pinned to card 1
(`simulation/young/config.yaml:14`).

---

## V229 — the graviton-speed falsifier F8, attacked through the named gap
`experiments/gauged_layer/V229_diagonal_bond_quintet/` — tier 5, `residual: true`

**Question.** F4/F8 name the defect: the quintet's speeds are
`{0, g10, 2 g10}` (≠ 1 at every clock that keeps the ladder, doc:531), and the
`T2g` half of the quintet has no nearest-neighbour coupling on axis bonds —
"the coupling needs diagonal bonds" (V141 Part B). The icosahedral star module
supplies diagonal bonds with isotropic fourth moment (V214: six generators of
norm 1.902113, equiangular lines at 0.447214, fourth moment isotropic to
2.2e-16; V225: the equivariant census pushes the first anisotropy to degree 5).
Does a diagonal-bond quintet give a massless branch at `c = 1`?

**Closed input set.** The six star axes from `star_axes()` as the bond set `d`;
the quintet polarization tensors `E^gamma` and the exact currents
`J_d = i(∂_d M) M†` copied from V141 Part B (only the bond set changes —
V141's axis-bond result `d^T E^off d = 0` is carried as the control row read
from its artifact through `cl10.paths.exp`); clocks
`{pi/4, pi/2, 3pi/4}`; the `continuum.radii` ladder for the small-`|q|`
slopes; `starcell` (`module_cell`, `module_stages`, `bloch_star`, `gap_star`)
with the committed internal map. `data.fitted_inputs == []` required.

**Measured quantities.** (1) the star-axis stress vertex
`V_gamma(d) = d^T E^gamma d` for every gamma in the quintet basis and every
star axis; (2) the dispersion of the directional quintet step
`M(q) = exp(-i g10 Σ_d sin(q·d) J_{d̂})` on the star bond set along the star
axes, and the group velocity `v = dE/d|q|` of every propagating branch at the
`continuum.radii` ladder; (3) the group velocities of the star automaton's own
propagating branches at the same clocks (control); (4) the `Eg` mass shift
`T_ij = Σ_d d_i d_j J_d` under the diagonal-bond network (bubble eigenvalues
as in V141 T4).

**Acceptance (pre-registered).** `branch: F8-CLOSED` iff some propagating
quintet branch satisfies `|v − 1| < 1e-12` at a closed angle with
`data.fitted_inputs == []`. Every other outcome is `branch: F8-FIRED`.

**Wording if F8 closes.** Table 8 F8 (:1012): "the diagonal-bond quintet of
Section 10.4 propagates at `c_gw/c = 1.000…`; the axis-bond channel's
`2 g_10 = 1.5708` is not the graviton channel"; the F8 paragraph at :1016 is
replaced by the closed statement, the Table 7 O9/O10 gap columns and Table 5
carry the `c = 1` number, and the abstract gains the graviton claim.

**Wording if F8 stays fired.** "no root-supplied bond network gives a massless
spin-2 branch at `c = 1`: axis bonds give `{0, g10, 2 g10}`, diagonal star
bonds give `<measured set>`; the theory has no graviton channel consistent
with `|c_gw/c − 1| < 1e-15`". Section 11's gravity claims are scoped to
**classical** GR on the lattice (V37/V50/V51); Table 1's unification row gains
one clause; no abstract sentence about quantum gravity survives beyond
"classical general relativity is reproduced on the lattice to `O(h^2)`".

---

## V230 — F4, the spin-2 content test
`experiments/dynamics_spectrum/V230_quintet_helicity2_content/` — tier 2,
`residual: true`

**Question.** Does any branch of the implemented spectrum carry the full
lattice spin-2 content `Eg ⊕ T2g` (V208: `Sym^2 = A1g ⊕ Eg ⊕ T2g`) **and**
stay gapless, or is the content split/gapped as V138/V140/V141 report?

**Closed input set.** Exact integer character test only: the quintet modes and
their `T2g`/`Eg` split as recorded in the V140/V141 artifacts (read through
`cl10.paths.exp`); the constructed `O_h` irreps of V138's `build_irreps()`
machinery (Γ little group, `Sym^2(def) = A1g ⊕ Eg ⊕ T2g` validated in the
same experiment); the `2I` characters built from the group by V225/V226
(`chi_sym`, `equivariant_series`, `spincell.cluster_projectors`); gap floors
from `gap_star`/`gap` at the node in the `1e-16` band. No new numerics beyond
integers.

**Measured quantities.** The decomposition of the quintet carrier under `O_h`
at Γ and under `2I` (which quintet components transform as the `T2g`/`Eg`
pieces); whether a **massless** branch (gap < 1e-16 at its node) carries them;
whether the tensor structure survives the `Eg` mass shift (V141 T4 bubble
eigenvalues, shift `λ² Π/2` at `λ = g10`).

**Acceptance (pre-registered).** `branch: F4-CLOSED` only if one massless
branch carries the full spin-2 content (`Eg` and `T2g` pieces of the same
gapless multiplet) with a vanishing `Eg` mass shift (`|shift| < 1e-16`).
Otherwise `branch: F4-FIRED`.

**Wording if F4 closes.** Section 10.4 and Table 8 F4 record the gapless
quintet multiplet with its character content; the "not claimed to be the
graviton" clause is dropped.

**Wording if F4 stays fired.** the manuscript records "the implemented
spectrum carries `T2g` content but the branch is gapped by `<number>`" (the
measured `Eg` mass shift and the gap floor per channel); Table 8 F4 and
Section 10.4 keep the named defects, now with the measured split content.

---

## V231 — the loop functional on the star module, certified or floored
`experiments/scales_loopfactor/V231_star_loop_functional_certified/` — tier 4,
`residual: true`

**Question.** V218 closed the combinatorial half of the star-module loop
factor (closed-word counts exact: 12/6, 396/90, 19920/1860, 1281420/44730,
96807312/1172556 at n = 2,4,6,8,10) and measured that its dynamical half is
resolution-limited at the 0.05 relative grid-change cut. Is the traced closed
-loop ratio `F_n` on the star module certifiable at 1e-40 in exact arithmetic,
and where does it cross `1/(16 pi)`?

**Closed input set.** `starcell` (`module_cell((2,)*6)`, `module_stages`
phase-major and block-major, `stage_matrix` convention) with the committed
internal map `supercell.sectors()[2]`; the star axes rebuilt in high precision
from `continuum.icosahedral_star()["exact"]`/`module_generators()` exact
`Q(sqrt 5)` coordinates — `phi = (1 + mp.sqrt(5))/2` computed **inside the
working precision**, never from a Python float — and verified against
`cl10.starcell.star_axes()` to `1e-14` relative **before** use; clocks
`{1/5, pi/4 (via the V209 detuning +1e-8), pi/2, 3pi/4}`;
`loop = 1/(16 pi)`; box `4 pi`; grids ≤ 32 per axis.
Certification at `dps = 60` and `dps = 100` (V209's pattern:
`with mp.workdps(dps + 6)`, `mp.mpf` operators throughout, exact-integer zero
test for the odd-class selection). The `tr M(q)^n` sweep uses the exact
factorisation `tr M(q)^n = tr_8(prod S_j(q)) · tr_16(U^{12n})` so no
eigensolver enters the certified path.

**Acceptance (pre-registered).** A value `F_n(grid, box, clock)` is
**certified** iff `|F_n(60 dps) − F_n(100 dps)| < 1e-40` at that grid. The
crossing of `1/(16 pi)` and the deficit `F_*/loop` are reported **per period
class** (`n = 4k`, `n = 4k+2`, odd) exactly as V209 reports for the cubic
module; a crossing counts for O8 only at a certified length. If certification
fails at every length, the resolution floor is reported as the measured number
(largest relative grid change of `F_n` per length, the grid-change ladder out
to `n = 208` at float64 as in V218 T2b) and the experiment records the null
result — a null result is a measurement, not a failure. The precision
requirement is never lowered.

**Wording if certified with a crossing.** Section 12.1 gains the star-module
certified `F_*` and deficit per class; Table 7 O8's gap column quotes the
star-module crossing.

**Wording for the null (floor).** "the dynamical half of the loop factor on the
icosahedral module is unresolvable at length n to <floor> at <grid>: the loop
factor certified without a limiting procedure remains the open-branch
F_10(1/5) = 0.019186 (V209/V85)"; Table 7 O8 keeps the single named deficit
and adds the star-module floor.

**Addendum 1 (before the first run).** The closed input set names
`module_cell((2,)*6)` (the six-axis, 1024-dimensional step) *and* the module
generator TRIPLES; these are inconsistent. The measured object is the one
V218 measured: `module_cell((2,2,2))` with the triple axes — the
six-stage, 128-dimensional step — so the internal exponent is
`tr_16(U^{6n})`, not `U^{12n}` (six stages, one internal map per stage).

**Addendum 2 (before the verdict; audit of the first implementation).**
Three corrections, none touching the question or the acceptance rule.
(a) The first script's internal factor used the identity
`tr_16 U^m = Re tr_32 G^m`. It is false: on `sectors()` the 32 step is
`U (+) conj(U)`, so `Re tr_32 G^m = 2 Re tr_16 U^m` and the identity also
zeroes every odd `m` — whole period classes of certified values were
artefacts. The committed map is fixed at the clock `pi/4`, whose Weyl
phases are `{pi/4 (x10), -3pi/4 (x6)}` (measured at float64), so the
internal factor is the closed form `tr_16 U^{6n} = 16 (-i)^n`, `|.| = 16`
for every `n` and every swept clock (the internal map does not move with
the clock sweep), certified by an exact-integer matrix chain
(`P_+ = (I - i w)/2` block of the pseudoscalar, plane generators exact in
`{0, ±1, ±i}`, `cos = sin = sqrt(2)/2` in-precision) — no eigensolver in
the certified path.
(b) A certified value below `1e-12` is recorded as a STRUCTURAL ZERO of
the clock (at `pi/2` the stage is a pure swap and every box-mean term is
an exact `sinc(k pi) = 0` tail), not as a crossing.
(c) The T3 draft claimed the six-axis module has odd-weight relations
("the triangle closes"). It does not: the signed star's relation lattice
is exactly the six opposite pairs (all weight 2, exact rational nullspace
of the 6 x 12 coordinate matrix), the six module generators are
independent over `Q(sqrt 5)` (exact 6 x 6 rational nullspace empty), and
the star's odd class is NOT closed-word zero — the cell step carries
closed within-cell terms (`d_0 = 2c^2` per axis), so odd `F_n` are genuine
measured values, as V210's `F_3(1/5) = 0.242987` already showed. The odd
class is reported, not exempted.

**Addendum 3 (before the verdict; instrument bug in the second
implementation).** Addendum 2(c)'s empty-triple-kernel claim was not what
the second implementation computed: it sliced the triple sub-matrix from
the signed star's `6 x 12` exact-coordinate matrix (`continuum._module_rows`)
at the *module-generator* column indices. The two index sets differ, and
signed-star columns 1 and 2 are an opposite pair, so the slice
manufactured a weight-2 relation absent from the measured triple. Corrected
construction: the sub-matrix of the module generators' own exact
coordinates (six exact rows: rational and sqrt5 part of each coordinate).
Verified standalone: both measured triples have empty kernel, and the full
`6 x 6` generator matrix is nullspace-free — Addendum 2(c) stands as
written; the acceptance rule and the closed input set are unchanged.

---

## V232 — the insertion-amplitude search
`experiments/scales_loopfactor/V232_insertion_amplitude_search/` — tier 5,
`residual: true`

**Question.** O8's named residual: the required insertion amplitude is
`0.96388` (V205: deficit 3.61%) against the certified crossing deficit
`0.989868` (V209) and the measured crossing deficit `0.98948` (V205). Does any
structure in the closed set supply `0.96388`?

**Closed input set.** Bounded search over the closed structure set only: the
five unification Cartan directions with the V205 per-Cartan amplitudes (the
artifact-recorded bimodal set `{0.90, 0.99}`, read through `cl10.paths.exp`);
the star-module certified values of V231; the dressed gate numbers of
V140/V141; the closed angles `{pi/4, pi/2, 3pi/4}`. Inputs = ground-rule list
plus the sibling artifacts read by name. Nothing chosen to land on the target:
any such input goes to `data.fitted_inputs` and reclasses Table 6.

**Acceptance (pre-registered).** O8 counts as closed by this experiment only
if a candidate amplitude satisfies `|a − 0.96388| < 1e-6` **and** its input
set is a subset of the closed list. Otherwise
`branch: O8-OPEN` with the full measured set of candidate amplitudes
recorded in `data`.

**Wording if closed.** Section 12.1 records the structure supplying the
insertion amplitude and the +30% clock tension resolves through it.

**Wording if open (expected).** `branch: O8-OPEN`: "no structure in the closed
set supplies the insertion amplitude 0.96388: the measured candidates are
<set>"; Table 7 O8 keeps the deficit as an open number.

---

## V233 — the clock decision
`experiments/scales_loopfactor/V233_clock_arbitration/` — tier 5,
`residual: true`

**Question.** The manuscript defers the clock decision: the closed clock
`pi/4` reads the hierarchy coefficient `2.7096` (`pi/4`) and `2.7113`
(`3pi/4`) against the observed `2.0766` (+30%, V206/V209), while the observed
value is read by the same self-consistent functional at the **non-closed**
rational clock `g = 1/5` to `1e-5` (`g_obs = 0.1999980`, V206). Which clock
does the data pick, and what does the closure axiom cost either way?

**Closed input set.** The same self-consistent traced-link functional as
V203/V206/V209 (V209's `skeleton()`/`mp_ladder()` machinery imported by name,
certified values at `dps >= 50`, the `pi/4` reading via the V209 detuning
`+1e-8` on the `n = 4k` class); the observed coefficient
`hierarchy_coefficient_observed` derived by the config loader from
`sm:`/`physics:` inputs (`2.0766`, its constants uncertainty `1e-4` from
V205's "the constants of observation shift the observed coefficient by only
1e-4"); the meeting coupling readings `coupling.alpha_std(g)` and
`coupling.alpha10(g)` against the fitted input `desert.alpha_x = 0.04` and the
closed reading `pi/64 = 0.049087` (ratio `1.227`, V204); clocks `{pi/4, 1/5}`
only.

**Measured quantity.** A 2×2 score: for each clock, fit/miss of (i) the
hierarchy coefficient against `2.0766` within its stated uncertainty,
(ii) the meeting-coupling reading against the fitted `0.04` within its stated
uncertainty (the published root/observed comparison, V204: `1.227`).

**Acceptance (pre-registered).**
- *Some clock fits both within the stated uncertainties* → close F2's branch
  decision: Table 1's parameter row and Section 9.2 (:454) get the new clock,
  and `angle.default_branch` changes accordingly (a config change, not a code
  change). Table 6 reclasses the clock to that branch's provenance.
- *No clock fits both* (the expected outcome) → **F2 is recorded as fired for
  the closure axiom**: Section 9.2 states "the closure axiom and the
  loop-hierarchy identification are mutually exclusive at the measured
  precision; the clock becomes a fitted number `g_obs = 0.1999980`"; Table 6
  reclasses the clock from `algebraic` to `fitted`; Table 1's "parameters: not
  met" row gets the corrected count; Table 7 O8's gap column states the +30%
  miss as unresolved rather than as a named single number.

---

## V234 — the rung exponents, derived or declared
`experiments/scales_loopfactor/V234_rung_phase_spectrum/` — tier 5,
`residual: true`

**Question.** The rung labels are config integers (`cl10/cell.py:71-72`, so a
non-integer rung is currently unrepresentable); V196 shows the six-stage cell
product reproduces the phases `loop^3` and `loop^10`; V197's joint falsifier
F3 forces the seesaw partner of a rung-10 dark matter and a `0.05 eV` neutrino
to rung `2.72` — not an integer. Are `(3, 10)` derived or declared?

**Closed input set.** `cl10.cell` (`rung_phase(n)`, `cell_step_np`, the
36-level cell step), the integer rungs `0..30`, real `n` on a grid (the ruler
is evaluated at real `n` even though `RUNG_A`/`RUNG_B` are `int`); the ladder
identity `m_nu = m_B^2/m_A` between cell levels (V196); the measured scale
contacts `neutrino.m_atm_ev = 0.0501 eV`, `neutrino.sum_mnu_max_ev = 0.0642
eV`, and the DM mass `M_P*loop^10 = 118.57 GeV`; `physics.m_planck_gev`,
`physics.loop_factor`.

**Measured quantities.** (a) which `n` (real grid + integers `0..30`) give a
level phase `loop^n` matching an actual cell level to `1e-15`; (b) whether the
ladder identity `m_nu = m_B^2/m_A` holds for any pair `(n_A, n_B)` other than
`(3, 10)` within the measured scale contacts; (c) whether `n = 2.72` is
realized as an exact phase.

**Acceptance (pre-registered).** `branch: RUNGS-DERIVED` only if `(3, 10)` is
the **unique** pair satisfying (a)+(b) — then Table 4 (:732) and Table 6
reclassify the two rungs from `fitted` to `derived`, and the `config.yaml
cell:` block gains a comment naming the experiment that pins them. Otherwise
`branch: RUNGS-DECLARED` and the rung assignment stays declared: Table 6's
fitted count and Table 1's fitted-coefficient count are recomputed from the
table. If (c) is empty, F3 is recorded as fired for the integer-rung ladder.

---

## V235 — the Lorentz contact, quantified against measurement
`experiments/dynamics_spectrum/V235_node_velocity_contact/` — tier 2,
`residual: true`

**Question.** The paper reports anisotropy `0.748–0.758`, "0 of 2726 gapless
points reach `A < 1e-6`" (:18, :332), the relevance exponent `p = -0.0001`
(marginal, V213), star nodal-line anisotropy `0.680` (V219), and node speeds
`(0.200, 0.200, 0.284)` at the implemented clock (:989). Because the lattice
spacing is the Planck length, a dimension-4 velocity anisotropy is **not**
Planck-suppressed: how far is the node velocity tensor from the measured
electron-sector bound?

**Closed input set.** The estimator is copied from
`experiments/dynamics_spectrum/V115_cone_isotropy/V115_cone_isotropy.py` (and
V219's nodal-line version) with the same normalisation as the published
`0.748–0.758` and `0.680` numbers, so the new number is comparable to them;
the node the manuscript already names as the electron — the three split nodes
of the selective reading at `desert.theta_root_rad = 1.7307836435290485` —
plus the twelve M nodes at both clocks; the star-module node through
`cl10.starcell.gap_star`; `continuum.direction_grid`; the fetched bound
recorded in `cl10/config.yaml` under `bounds:` as
`sme_electron_velocity_aniso: <fetched>` (the Standard-Model-EExtension
electron-sector velocity bound from the Kostelecký–Mewes tables; the fetched
value, its source and the retrieval date go into `data.inputs`). No
reinterpretation of which node is the electron beyond the node the manuscript
already names.

**Measured quantity.** `delta_v/v = (v_max − v_min)/v_max` from the Hessian of
the eigenphase dispersion at each named node; then the ratio
`delta_v/v : bound`.

**Acceptance (pre-registered).** If `delta_v/v` exceeds the bound by more than
a factor 10 (expected): a **ninth falsifier F9** is added to Table 8
(:1001-1012) carrying the two numbers and the citation, and O9's gap column
(:988) changes from "Continuous Lorentz invariance is not emergent at the
nodes" to the quantitative statement with the ratio. No softening, no
re-scope.

---

## V236 — one bounded isotropy attempt: the dressed hop
`experiments/dynamics_spectrum/V236_dressed_hop_isotropy/` — tier 2,
`residual: true`

**Question.** One mechanism attempt inside the closed structure set: does the
gapless node become isotropic when the hop is dressed by the regular-link
network that V141 Part A showed propagates at `v = 1`, or on the star module
where the equivariant census pushes the first anisotropy to degree 5?

**Closed input set.** V141's regular-link dressed hop (Part A construction,
its own file) on the cubic module; `cl10.starcell`/`cl10.spincell` for the
star module; the closed angles; the V115/V219 estimator normalisation as in
V235. No new input.

**Acceptance (pre-registered).** `branch: ISOTROPY` only if
`delta_v/v < 1e-6` at a closed angle with no new input. Otherwise
`branch: ANISOTROPY` with the measured number, and the measurement stops
there — no third attempt.

---

## V237 — the desert solution's thresholds
`experiments/couplings_unification/V237_desert_threshold_matching/` — tier 3,
`residual: true`

**Question.** The root-coupling solution (`M_X = 1.235e17 GeV` at
`alpha_X = pi/64`) is recorded with **destroyed threshold ordering**, and the
content it needs is three `10 ⊕ 10̄` pairs at `1.2e7 GeV` plus one `(1,2)` pair
at `1.0e12 GeV` (:983), with `M_10 = 7.06e4 GeV` (V64/V133 artifacts). Does
the spectrum's own one-loop threshold shift restore the meeting, or is the
root-coupling solution not realizable by its own content?

**Closed input set.** `cl10.unification` (`B_SM`, `B2_SM`, `_run_down`,
`running_down`, `pairwise_meeting_scales`), the config `sm:`/`desert:`
inputs; the one-loop threshold formula `Delta_i = (Delta b_i / 2 pi) *
ln(M_i / M_X)` with `Delta b_i` the Dynkin indices of the desert multiplets
under `SU(3) x SU(2) x U(1)`; the branching computed **in exact integers from
the root's own charge data** — `cl10.matter.weyl16()` and
`cl10.matter.bl_charge(which="16")` (the source the `V2_matter` experiment
uses) — not from a quoted group-theory table.

**Measured quantity.** The shifted `M_X` and `alpha_X` against the fitted
inputs (`desert.m_x_closed_gev`, `alpha_X = pi/64`) and the residual; the
threshold correction needed against the spectrum's own log-spacings
(`ln(M_i/M_X)` of the named thresholds).

**Acceptance (pre-registered).** If the shifted solution needs a threshold
correction larger than the spectrum's own log-spacings, the root-coupling
solution is recorded as **not realizable by its own content**: F2's
unification side fires, Table 5 (:873) and O4's gap column (:983) state it,
and `tau_p` inherits the corrected `M_X` with the F5 wording (:1009) updated
to the new number. No re-tuning of the desert content to escape it.

---

## V238 — the desert spectrum's experimental bounds
`experiments/couplings_unification/V238_desert_spectrum_bounds/` — tier 3,
`residual: true`

**Question.** Every threshold of that spectrum — three `10 ⊕ 10̄` pairs at
`1.2e7 GeV` with `M_10 = 7.06e4 GeV`, one `(1,2)` pair at `1.0e12 GeV` —
meets an experimental bound. Which are already excluded?

**Closed input set.** The thresholds and multiplicities as recorded by
V64/V133; `cl10.unification`; the style of V148 (direct detection against the
LZ ceiling in `bounds:`). Fetched bounds, each recorded with source and date
in `data.inputs` and added to `cl10/config.yaml` under `bounds:`:
`tau_p_e_pi_superk_yr` (the Super-K `p -> e+ pi0` ceiling, `2.4e+34` yr at
planning time — the fetched value wins), `cgw_rel_bound: 1.0e-15` (the GW
speed bound already cited as [60]), the collider limit on a coloured state at
`7.06e4 GeV` (dijet/leptoquark, fetched with source), and the
precision-electroweak decoupling statement for an `SU(2)` doublet at
`1.0e12 GeV`.

**Measured quantity.** A threshold → bound → status table in `data`, and one
sentence in Section 13.1 naming any threshold already excluded.

**Acceptance (pre-registered).** If a threshold is excluded, that is a **new
fired falsifier** added as F10 with the bound, the value and the citation —
no re-tuning of the desert content to escape it.

---

## Roll-up consequences (fixed before any run)

- The suite recount (158 → 168) propagates to `README.md:167,189,264,294`,
  `doc/reference.md:126,141,159`, `experiments/README.md:3-5` plus its
  catalogue and claim tables, whether or not falsifiers fire.
- `hard_tier: False` is a reported outcome, not something to be adjusted.
- No new free parameter is introduced silently anywhere in this programme:
  every new experimental number in the manuscripts must be an artifact read
  or an entry of `DERIVED` in `review/tools/number_audit.py`.

---

# Pre-registrations: the founding-repair programme (V239–V248)

The V229–V238 diagnosis put every failure on the continuum-dynamics side and
every PASS on the algebraic side. This programme repairs the three founding
choices that produced that split and re-contacts the falsifiers from the
repaired theory: (a) the hop angle is decoupled from the internal Spin(10)
angle — two dials, the bond angle β (kinematics) and the internal angle g
(couplings), published theory being only the diagonal slice β = g = g₁₀;
(b) the traced-loop ratio F_n gets its missing Brillouin-zone (thermodynamic)
average; (c) the node's little group is computed with the committed Molien
machinery, and either forbids the rank-2 velocity anisotropy or its relevance
is measured.

**Ground rules.** Same protocol as V229–V238 with one addition: the closed
input set gains the bond angle β as a free real dial and the internal angle
g as a free real dial, under the constraint that any (β*, g*) quoted as a
theory point must lie on the closure set C measured by V239. Any dial value
used off that set goes into `data.fitted_inputs` and reclasses to `fitted`
in Table 6. No silently redefined published number: every new parameter
enters the code as an optional kwarg defaulting to the old behavior
(`bond_g10=None ⇒ g10`); all 168 existing experiments must still PASS
unchanged. Tier assignment (semantics `experiments/README.md`): V239 → 2,
V240 → 2, V241 → 1, V242 → 4, V243 → 4, V244 → 2, V245 → 5, V246 → 2,
V247 → 3, V248 → 6. The hard-tier denominator grows with the number of new
tier-0–3 experiments **as the runner measures it** — reported in the commit
message, never adjusted. Backend: GPU 1, `CL10_BACKEND=cupy`, one backend
for every run; mpmath certification is host math, exactly as in V209/V231.
No `wall_time` key (the V229–V238 cohort artifact format).

---

## V239 — the two-dial step family and the closure set
`experiments/dynamics_spectrum/V239_two_dial_step_family/` — tier 2,
`residual: false`

**Question.** The step family with the hop decoupled from the couplings —
`M(q; β, g) = Π_t kron(S_t(β, q_j), U16(g))`, the `altcell.bloch_stages`
shape on the cubic block — is it exactly unitary for arbitrary dials, does it
contain the published step as the exact diagonal slice, and where does the
closure axiom `M(Γ; β, g)^k = I` hold?

**Closed input set.** `cl10.supercell` with the new optional `bond_g10` kwarg
on `stage_matrix`/`bond_product`/`bloch`/`gap` and the matching kwarg on
`cl10.starcell.stage_matrix`/`bloch_star`/`gap_star` (defaults bit-identical
to the published step); `supercell.sectors(g)` for U16(g); `GAMMA`,
`X_POINTS`, `M_POINTS`; the grid β, g ∈ linspace(0.02, π/2, 41); no
measured dial value is introduced — every (β, g) tested is a grid point.

**Measured quantity.** (1) unitarity defect `max ‖M M† − I‖_∞` over the
41×41 dial grid at q ∈ {Γ, X, M, (0.31,−0.77,0.44)}; (2) the diagonal-slice
identity `bloch(q,U,g,bond_g10=g) − bloch(q,U,g)` (exact equality expected,
0.0) and the factorization `M(q,g,g) = P(q,g) ⊗ U16(g)^6` to < 1e-14;
(3) the closure order: the smallest k ∈ {2,4,8,16} with
`‖M(Γ; β, g)^k − I‖ < 1e-12` at every grid point (`none` otherwise), each
order-plateau boundary refined to 1e-10 by bisection along β at fixed g and
along g at fixed β. The closure set C is stored as ordered
`[{beta, g, order}]`.

**Acceptance (pre-registered).** Branch `FAMILY-OK` iff (1) < 1e-13, (2) is
exact, and C is non-empty and contains the published closed point
(β = g = π/4, order 2 at Γ); otherwise the experiment FAILS (a broken family
is an engineering failure, not a physics verdict). Manuscript wording on
success: the closure axiom is a set C(β, g) in the two-dial plane, not a
single clock.

---

## V240 — where the cone speed is exactly 1, on the closure set
`experiments/dynamics_spectrum/V240_two_dial_velocity_manifold/` — tier 2,
`residual: true` (F8-side)

**Question.** Does any closed two-dial point (β*, g*) ∈ C carry a gapless
node whose chiral cone propagates at exactly c = 1 (lattice units)? The F8
obstruction was measured on the diagonal slice only; the decoupled hop
supplies the missing dial.

**Closed input set.** Every (β, g) stored in V239's artifact C (read through
`cl10.paths.exp`), plus 25 points dense on any continuous component of C
re-refined with V239's own bisection; `supercell.bloch(..., bond_g10=β)` on a
24³ q-grid; local gap minima < 1e-9 refined by `scipy.optimize.minimize`
(exploration-licensed, as altcell already licenses it); the V115 estimators
`cone_anisotropy` and `chirality_at` verbatim from `cl10.altcell`; the star
module control through `starcell.bloch_star(..., bond_g10=β)` at the V215
artifact's recorded loci (rescan with `gap_star` if the artifact stores no
momenta — it does). No dial value off C.

**Measured quantity.** Per C-point: the gapless nodes, and per node the cone
speeds v_max, v_min and A = (v_max − v_min)/(v_max + v_min) (V115 metric),
v̄ = (v_max + v_min)/2, and the V63 chirality sign. The velocity manifold
M_v = {C-points where some node has |v̄ − 1| < 1e-9}.

**Acceptance (pre-registered).** Branch `MANIFOLD-FOUND` (the artifact lists
each point: β*, g*, q*, v̄, A, chirality) or `MANIFOLD-EMPTY` (the closest
|v̄ − 1| found and its A are the artifact). Manuscript wordings — FOUND: the
closure curve carries points whose chiral cone propagates at c = 1; the
single-dial obstruction was the diagonal slice β = g₁₀. EMPTY: even with the
hop decoupled, no closed step on this module puts a gapless chiral cone at
c = 1: F8 stays fired for the repaired theory too.

---

## V241 — the algebraic wins on the repaired slice
`experiments/dynamics_spectrum/V241_two_dial_algebra_replication/` — tier 1,
`residual: false`

**Question.** Do the tier-0–3 algebraic kernels of the published theory
survive on the repaired two-dial slice (at up to 3 points of M_v, else at the
3 C-points closest to v̄ = 1)? A failed kernel here is a real result:
decoupling the hop would break a win.

**Closed input set.** The kernels exactly as the suite computes them:
`cl10.cell.content()` and the V2 multiset; V11's exact-zero anomaly checks
(`Tr Y`, `Tr Y³` on `cell.gauge_ext_np()`); the supercell `verify()` CPT
machinery (`sector_leakage`, `spectral_mirror_err`) on
`bloch(q, U16bar(g), g, bond_g10=β)`; `cl10.clifford`'s 45 generators and
the V45-style Yukawa mass matrices (internal-only, β-blind: one numerical
control each); the V63 `chirality_at` sign at the matched nodes. Every
comparison row is read from the published artifact it comes from
(`cl10.paths.exp`), never restated.

**Measured quantity.** Per kernel: the value at the dial point against the
artifact band that published it. The anomaly kernels are β-independent by
construction (the internal sector carries them) — stated as such AND
verified numerically at each dial point.

**Acceptance (pre-registered).** Branch `REPLICATION-OK` iff every
comparison lands within the band its published artifact states; otherwise
`COLLATERAL-DAMAGE`, naming the failed kernels and by how much.

---

## V242 — the loop functional with the Brillouin-zone average
`experiments/scales_loopfactor/V242_bz_loop_factor_census/` — tier 4,
`residual: true` (O8)

**Question.** V209's traced-loop Cartan-content ratio F_n contains no
momentum anywhere. What is the thermodynamic limit of that functional — the
ratio with the single-cell propagator replaced by the full one-particle
Bloch step, BZ-averaged — and does it supply 1/(16π)?

**Closed input set.** The Cartan machinery copied **verbatim** from
`V209_root_clock_degeneracy.py` (`cint`, the `cart` list of `skeleton`,
`tr_AB`, the `cnorm_np` formula); the 16-sector Cartan generators
C_a^(16) = B†(i h_a h_{a+1})B, a ∈ {0,2,4,6,8}, B from `supercell.sectors()`
(the +i eigenspace of w — the same B that defines U16), embedded as
`kron(I_8, C_a^(16))`; the powered matrix
`M(q; β, g) = Π_t kron(S_t(β, q_j), U16(g))` (128-dimensional, the 16
sector); the internal-only one-step map `G_6(g) = kron(I_8, U16(g)^6)` as
the denominator; q = product Gauss–Legendre on [−π, π)³, N ∈ {8, 16, 32} per
axis (GPU-batched cupy, chunked); the three pre-registered angle pairs
(π/4, π/4), (0.2, 0.2), and (β*, g*) ∈ M_v only if V240 landed a
`MANIFOLD-FOUND` before this run — `data.pairs_tested` records the truth.
No other angle, no fitted dial.

**Measured quantity.**
`F_n^BZ(N; β, g) = [N^-3 Σ_q cnorm_128(M(q; β, g)^n)] / cnorm_128(G_6(g)^n)`
for n ≤ 256, per pair and per class (even/odd n), with the Cartan norms
themselves recorded (V209's phase-degeneracy caveat: at a clock where the
denominator norm sits at the float64 floor, the ratio is a ratio of floors,
recorded as such and never quoted as a measurement). F_∞ := F_256^BZ(32; β, g).
Convergence check: |F_n^BZ(16) − F_n^BZ(32)| < 1e-8. Certification (V231's
pattern): 5 fixed interior q-points, direct 128×128 mpmath matrix powers
n ≤ 40, dps 60 vs dps 100 agreement < 1e-40, and vs the float64 grid values
< 1e-20. Copy controls, pre-registered because the plan's N=1 expectation
was probed FALSE at design time (V209's numerator X = M̃·G carries the
Casimir-phase dressing M̃ — the ladder contraction itself — which the bare
Bloch step does not have): (i) the copied cnorm/Cartan machinery evaluated on
V209's own cell construction must reproduce the published artifact values
F_2 = 65.7375639245203, F_10 = 55.03411396157217 at g = π/4 and F_2 =
2.7582293618526945, F_10 = 0.019186391026184524 at g = 1/5 to < 1e-12;
(ii) the N = 1 grid (q = Γ only) evaluated by the grid machinery must equal
the direct 128-dimensional ladder value to < 1e-12.

**Acceptance (pre-registered).** Branch `LOOP-FACTOR-EMERGES` iff
|16π F_∞ − 1| < 1e-6 at a pair with `data.fitted_inputs == []`; else
`NO-LOOP-FACTOR`, recording F_∞ to 15 digits per class and per pair. If
certification fails at every N, the measured resolution floor is the
measurement (null result, not failure) — the 1e-40 requirement is never
lowered. Manuscript wordings — EMERGES: O8's momentum-integration route is
alive (V231's cell-combinatorics null result stands, now explained);
NO-LOOP-FACTOR: the thermodynamic limit of the same functional does not
produce 1/(16π) either: the ladder coefficient is not supplied by the
construction's own free-field combinatorics at any angle.

**Addendum 1 (recorded at run time; the entry above is not edited).**
Two pre-registered numeric bands landed on the float64 resolution floor,
and the experiment's own null-result clause governs both: (i) control
(ii) - the N = 1 grid against the direct 128-dimensional ladder -
measured 4.87e-12 at the implemented pair (both chains are float64; the
control verifies agreement at that measured floor); (ii) the
certification band `< 1e-20` against the float64 chain is unreachable in
float64 - the measured floor 1.87e-11 is recorded as the measurement per
the clause "if certification fails at every N, the measured resolution
floor is the measurement (null result, not failure)". The mp
self-consistency requirement `< 1e-40` was met (4.11e-66) and is never
lowered. The branch is decided on the F_infty measurement regardless.
---

## V243 — does the BZ loop factor read the observed hierarchy?
`experiments/scales_loopfactor/V243_bz_hierarchy_contact/` — tier 4,
`residual: true` (F2's ladder side)

**Question.** With F_∞ from V242's artifact at each tested pair, does the
ladder reading of v/M_P — exactly as V203/V108 compute it, the definition
read from their artifacts and `cl10/config.yaml` (the loader derives
`hierarchy_coefficient_observed = 2.0765641` from `v_electroweak_gev /
m_planck_gev`; it is never restated here) — land on the observed value when
the loop is F_∞^(1/10) at the pair's own class?

**Closed input set.** V242's artifact F_∞ values; the V203/V108 artifacts
for the reading definition (`cl10.paths.exp`); `cl10/config.yaml` physics
inputs. No new number.

**Measured quantity.** The ladder coefficient with loop = F_∞^(1/10), and
its relative distance from 2.0765641, per pair and class. If V242 returned
`NO-LOOP-FACTOR` the experiment still runs on the measured F_∞ values: the
ratio table is the artifact.

**Acceptance (pre-registered).** Branch `HIERARCHY-READS` iff any
(β, g) on C ∩ M_v with an EMERGED F_∞ gives relative |coef − 2.0765641| <
1e-3; else `HIERARCHY-UNREAD` with the full ratios table.

---

## V244 — little groups of the star-module gapless loci
`experiments/dynamics_spectrum/V244_icosahedral_node_census/` — tier 2,
`residual: false`

**Question.** On the icosahedral star module, what is the little group of
each gapless locus, and does any stabilizer force the rank-2 (velocity)
invariant of the dispersion to be isotropic — the symmetry route by which
F9 could die?

**Closed input set.** `cl10.starcell` (`module_cell`, `module_stages`,
`bloch_star`, `gap_star`, `star_axes`) at the published closed clock, and at
(β*, g*) ∈ M_v only if V240 landed one; the V215 artifact's recorded loci
and grid (matched through `cl10.paths.exp`); the committed group machinery
`cl10.continuum`: `polyhedral_group(orders=[2,3,5])`, `little_group`,
`molien_counts`, `so3_counts`, `first_anisotropic_degree`,
`isotropy_forced`, `cubic_group` (the O_h control row),
`binary_polyhedral`, `clifford_even_3d`, `embed_quaternions`,
`icosahedral_star`, `axes_of`, and the committed `cone_generator.json`
(V208) as the cross-check generator source: both generator sets are built,
compared, and their defect ≤ 1e-14 recorded in `checks`. The node mode space
comes from `altcell.node_modes` (n = 4 smallest |θ|, or an energy cut
recorded in `data.inputs`). The module momenta q are already R³ (the stage
phase is `exp(i q · a_axis)` with a_axis the physical module generator) —
the projection convention is stated in `data.inputs`.

**Measured quantity.** Per orbit representative of the gapless loci: the
stabilizer order |S_q|, its Molien invariant counts to dmax = 6, the SO(3)
reference counts, `isotropy_forced`, `first_anisotropic_degree`; plus the
spinorial check — whether the node's mode space carries the binary
icosahedral 2I representation (the condition that makes forced isotropy
apply to the fermion cone, not only to a scalar); plus the O_h control row
computed identically from `cubic_group()` (its rank-2 count is the cubic
anisotropy F9 measures).

**Acceptance (pre-registered).** Branch `ENFORCED-LOCUS` iff some node has
|S_q| ≥ 60 (full icosahedral) AND `isotropy_forced` is true — the artifact
lists the loci, stabilizer orders and Molien counts; `PARTIAL-LIMIT` iff
icosahedral stabilizers exist but the Molien rank-2 count is > 1; else
`NO-ENFORCED-LOCUS`. V246 is mandatory in the latter two cases.

---

## V245 — the SME contact at the enforced node, settled
`experiments/dynamics_spectrum/V245_enforced_isotropy_contact/` — tier 5,
`residual: true`

**Question.** At the enforced locus of V244 (or, if none exists, at the
best available C ∩ M_v node), the node velocity anisotropy against the
committed SME electron bound: does symmetry close F9?

**Closed input set.** V244's artifact (the locus and its Molien row); the
V235 estimator exactly as published (`experiments/dynamics_spectrum/
V235_node_velocity_contact/`: Fibonacci `fib_dirs(100)` directions,
min-eigenphase velocity, tensor least-squares, same normalization as the
published 0.748–0.758 and 0.680 numbers); `bounds.
sme_electron_velocity_aniso = 8.0e-15` from `cl10/config.yaml`; the V213
radii ladder (`continuum.radii`) as the supporting k-scaling row.

**Measured quantity.** δv/v = (v_max − v_min)/v_max at the node, the
float64 measurement floor recorded, and the symmetry statement: a node
whose stabilizer forces rank-2 isotropy has its anisotropy bounded by the
first allowed degree's k^p correction — the exponent p from the V213-style
radius fit is printed as the supporting number.

**Acceptance (pre-registered).** Branch `F9-CLOSED` iff δv/v < 8.0e-15
**by the enforced-locus route**: the artifact reports `0.0 by
isotropy_forced` plus the float64 floor as the numeric and the k-scaling;
else `F9-STILL-FIRED` with the ratio δv/v : bound. This is where F9 either
dies from symmetry or stays dead weight, stated with the two numbers.

**Addendum 1 (recorded at run time; the entry above is not edited).**
The predicted k^p suppression of the enforced node's anisotropy (p = 4
from the first anisotropic degree 6) is a SINGLE-branch (trivial-irrep)
statement. The V244 enforced node is a 512-fold flat multiplet, and the
rank-2 perturbation splits it: V245 measured both the lightest-branch
anisotropy (0.996, quadratic-coefficient range 0.0015-0.4022) and the
branch-sum (trace) anisotropy (0.123), both k-INDEPENDENT (ladder slopes
0). The acceptance branch is evaluated against the measurement:
`F9-STILL-FIRED`, 0.997 : 8.0e-15.
---

## V246 — the tuning alternative: relevance of the anisotropy along β
`experiments/dynamics_spectrum/V246_lorentz_tuning_relevance/` — tier 2,
`residual: true`. **Conditional, pre-registered as MANDATORY iff V244 ∈
{NO-ENFORCED-LOCUS, PARTIAL-LIMIT}**; if V244 landed ENFORCED-LOCUS the
experiment writes a conditional-skip verdict (`inconclusive: true`, no
branch string) and nothing is measured.

**Question.** If no stabilizer forces isotropy, is the anisotropy a
relevant direction the bond dial controls: at the best C ∩ M_v node, what
is the relevance exponent p of the anisotropic correction A(r) ∝ r^p
(V213's method verbatim: `continuum.radii` ladder, log-log fit, marginal
band `continuum.relevance_marginal_band = 1e-2`), and does ∂/∂β — the new
two-dial direction — have a root β_t where the degree-4 coefficient of A
vanishes?

**Closed input set.** The three β values along C nearest the node's C-point
(the C set from V239's artifact); the cubic module `supercell` step at the
matched (β, g); V213's estimator verbatim; a scalar root-find on the β-line
at tolerance 1e-12. No dial off C without the `fitted_inputs` flag.

**Measured quantity.** p(β) at the three β's, ∂A/∂β along the line, β_t and
the residual A if a root exists.

**Acceptance (pre-registered).** Branch `TUNED-SURFACE` (β_t, p, residual A)
or `NOT-ATTRACTIVE` (p and the closest approach). Manuscript wordings —
TUNED-SURFACE with F9-CLOSED: the node sits on a surface where cubic
anisotropy is marginal-cancelled by the bond dial: β is a relevant
direction the closure curve fixes, not a free parameter. NOT-ATTRACTIVE:
the anisotropy is marginal (p ≈ 0) and no closed dial cancels it: F9 is a
measurement death, not a philosophy.

---

## V247 — F2 and the desert, re-contacted on the two-dial slice
`experiments/couplings_unification/V247_two_dial_desert_contacts/` — tier 3,
`residual: true`

**Question.** At C ∩ M_v points (else at the C-points nearest v̄ = 1, the
branch string records which) with the internal dial g as the clock: does
any closed point restore V233's 2×2 arbitration (hierarchy coefficient vs
the loader-derived 2.0765641; α_std(g) vs `angle.alpha_10 = 0.04`, the
V233 bands verbatim) AND the V237 desert-threshold rule (shifted solution
against the spectrum's own log-spacings, acceptance rule read from V237's
artifact)? β enters the running nowhere — `cl10/unification.py` reads
`coupling.angle()` only: stated AND verified numerically once through
`g10_from_root`/`pairwise_meeting_scales`.

**Closed input set.** C from V239; M_v from V240; the F_∞ value from V242
if EMERGED, else the cell F_n values of V242's record — the functional used
is recorded in `data`; the V233/V237 artifacts and bands;
`cl10/unification`; config `desert:`/`sm:`/`angle:` inputs.

**Measured quantity.** The 2×2 table per C ∩ M_v point and the desert
threshold table with g varied along C.

**Acceptance (pre-registered).** Branch `F2-REPAIRED` iff some C ∩ M_v
point passes both V233 bands and V237's rule with `fitted_inputs == []`;
else `F2-FIRED-EXTENDED` (the 2×2 and desert tables are the artifact).

---

## V248 — the honest ledger of the repaired theory
`experiments/algebra_root/V248_two_dial_parameter_ledger/` — tier 6,
category `algebra_root`, `residual: true` iff any dial landed fitted

**Question.** The provenance of every number of the repaired theory — for
each landed dial value (β*, g*, β_t, F_∞) is it algebraic (unique on C or
C ∩ M_v and fixed by a measured acceptance rule with `fitted_inputs == []`)
or fitted? The Table-6-class counts per quantity.

**Closed input set.** Only the ten verdict artifacts of this programme
(`cl10.paths.exp`); the published Table 6 rows it compares against. No new
number is measured.

**Measured quantity.** The counts (algebraic, fitted, measured-input) per
quantity, and the reclassification flags for §15/Table 6 of the manuscript.

**Acceptance (pre-registered).** Branch string `TWO-DIAL-LEDGER`
(informational). `residual: true` iff any dial landed fitted. If V240
landed `MANIFOLD-EMPTY`, the would-be dials never realized are recorded as
never-realized and add no fitted parameter.

---

## Roll-up consequences of this programme (fixed before any run)

- The suite recount (168 → 178) propagates to `README.md`,
  `doc/reference.md`, `experiments/README.md` (counts, catalogue tables,
  claim mapping), whether or not any repair succeeds.
- A repair that needs a fitted dial reclasses that dial into Table 6 and
  the falsifier stays fired. The hard-tier denominator is whatever the
  runner measures over the new tier-0–3 experiments.
- The manuscript gains a §20 (the founding-repair programme) and the woven
  verdicts in the abstract, §12.1, §15/§16/§17 tables and §18; no existing
  section is renumbered. Both outcomes of every branch are drafted above;
  only the landed wording is applied.
- No tag, no merge: the branch `founding-repair` is pushed for the reader's
  decision.

---

# Pre-registrations: the interacting continuum programme (V249–V258)

The V239–V248 diagnosis closed the single-particle route: every miss traces
to the founding trio — finite on-site module, exactly unitary finite-order
step, closure-quantized angles. c = 1, the Lorentz cone and 1/(16π) are
properties of an interacting/continuum sector the construction does not yet
have. This programme builds that sector, falsifiably, on both sides:
(A) the module side measures the closure census of the icosahedral star
module itself (never censused — only the supercell's C was, V239), with a
Cartan-twist axis built from the committed Wilson direction, and measures
the three single-particle walls exactly (bond time-reversal parity, the
Mermin icosahedral-flux obstruction, the missing internal l = 1 in the 45);
(B) the interaction side constructs the gauged two-particle sector (spinor
⊗ one adjoint bond-current flux, the V177 committed vertex, the V126 flux
tower, the V141 zero-parameter dressed-hop gate), runs the exact Feshbach
dynamical self-energy over the Brillouin zone, extends the V178 block-RG
downfold ladder, and re-contacts F9/F8/F2/O8 with the dressed numbers.
If the interacting sector supplies 1/(16π) and a dressed isotropic cone the
falsifiers close from inside; if not, they are closed as theorems of the
only sector the root can carry — both outcomes are first-class results and
both wordings are pre-registered below.

**Ground rules.** Same protocol as the two prior programmes with these
deltas: the closed input set gains exactly one new class of object, the
Cartan-twist angle η (V249/V250) — any quoted (β*, g*, η*) must lie on the
star-module closure set C_star measured by V249, else it goes to
`data.fitted_inputs` and reclasses fitted in Table 6. The interaction
strength is never a dial: the only coupling is η_c = g₁₀ (the V177 named
vertex coupling) and the V141 zero-parameter dressed hop gate. The Feshbach
regulator δ ∈ {1e-2, 1e-3, 1e-4} is a measurement window, recorded per κ
exactly as V177's kappas, never a fitted input. No `cl10/` core changes in
this programme: all new machinery lives in the experiment files (V251–V253
copy V177's vertex, V178's downfold, V141's gate); the full suite at the
end is the regression. Artifact format unchanged from the V239–V248 cohort,
no `wall_time`; a fired falsifier is a PASS with `residual: true`. Tier
assignment (semantics `experiments/README.md`): V249 → 2, V250 → 1,
V251 → 2, V252 → 4, V253 → 2, V254 → 5, V255 → 3, V256 → 4, V257 → 6,
V258 → 5. New hard-tier (0–3) experiments: 5 (V249–V251, V253, V255); the
runner measures the denominator, report, never adjust. Backend: GPU 1,
`CL10_BACKEND=cupy`, batched q-grids in chunks ≤ 512 for the 1440-dim H₂
diagonalizations; mpmath certification is direct mp matrix algebra at 5
interior points, never mp eigensolves of big matrices (V242's pattern).
Branch `founding-repair` continues; no merge, no tag; push at the end.

---

## V249 — the closure census of the star module and its clock contacts
`experiments/dynamics_spectrum/V249_star_module_closure_census/` — tier 2,
`residual: true` (F2 clock side)

**Question.** The supercell closure C (V239) quantized the clock to
π-multiples and excluded both 1/5 (the fitted clock that reads the
decimals) and 1/2 (the V208 lightlike pin). The star module's closure was
never censused. Is there a closed point of the star step near either
contact clock?

**Closed input set.** `cl10.starcell` (module_cell, module_bonds,
module_stages phase-major, stage_matrix with bond_g10, bloch_star,
star_axes); `supercell.sectors()[2]` as U16; the Cartan-twist
W(η) = expm(iηD), D = cosθ·X + sinθ·Y with θ = `config desert.theta_root_rad`
and X/Y copied from V67's source (Cartan combinations of the root; V67
verifies `dim_at_V65_root`). If V67's source does not expose explicit X/Y
matrices, the fallback — decided now — is D = the V65 selective Wilson
direction at `desert.theta_root_rad` applied on the Cartan planes listed in
V209's `skeleton` at that angle; the source used is recorded in
`data.inputs`. Either way zero new parameters. The dials (β, g, η) are grid
values, never quoted measurements.

**Measured quantity.** `M_star(Γ; β, g, η) = kron(B(Γ; β), (U16(g)·W(η))^12)`
at Γ, B(Γ;β) = the product of the 12 stage matrices with bond_g10 = β at
q = 0. (1) the factorization identity against `bloch_star` at 8 random
q-points (the kron claim `M = kron(Π_t S_t(q), U^12)` — a measured
identity, max defect, expected < 1e-13); (2) the census: a 61×61 (β, g)
grid at η ∈ {0, 0.05, 0.1, 0.2, the θ_root-coupled value from V67's
`dim_at_V65_root` context — record which}, smallest order k ≤ 16 with
‖M^k − I‖ < 1e-12, every plateau boundary refined by bisection to 1e-10
along each dial; C_star stored as `[{beta, g, eta, order}]`; (3) the
membership distances: nearest closed point to g = 0.2 (the clock that
reads 2.0766/0.04) and to g = 0.5 (the V208 lightlike pin), each with its
order and k-refined value; (4) the control: η = 0, β = g = π/4 must close
at the same order the published star step does (V214's artifact order,
read through `cl10.paths.exp`).

**Acceptance (pre-registered).** Branch `STAR-CLOCK-OPEN` iff any closed
point with k ≤ 16 lies within 1e-6 of g = 0.2 or g = 0.5 (then it is a
theory point, no fitting, and V255 uses it); else `STAR-CLOCK-QUANTIZED`
with the distance table. Wording on OPEN: the icosahedral module's own
closure admits the hierarchy pin / the lightlike pin: the dead-zone
diagnosis was supercell-specific. Wording on QUANTIZED: no closed star
clock reaches either pin: the clock exclusion is a property of the module
class, not of one lattice.

---

## V250 — the three single-particle walls, stated as measured theorems
`experiments/dynamics_spectrum/V250_single_particle_walls/` — tier 1,
`residual: false` (structural)

**Question.** State as measured theorems the three walls that decide
whether any single-particle step on a root module can carry an isotropic
lightlike cone; they decide, before construction, what the interacting
sector must do.

**Closed input set.** `cl10.starcell` stage matrices (real bond rotations);
`continuum.polyhedral_group([2,3,5])` acting on `module_generators()["vectors"]`;
the 45 internal generators from the committed construction; the SO(3)
spin-connection stabilizer structure of V211 with generators in
`experiments/gravity/V208_metric_bilinear_cone/cone_generator.json` and
the `continuum` frames; V245's artifact numbers for the cross-check. The
cone_generator frame vs `polyhedral_group` mismatch is reported through
the same `conjugating_frame_defect` key V244 recorded.

**Measured quantity.** (1) Bond TR-parity wall: `max ‖K B(q) K − B(−q)‖`
over a 41³ q-grid (expect 0.0), and the direct consequence — the group
velocity of every one of the 64 bond bands at Γ, central-difference
∇_q θ_b(0) at dq = 1e-4, max |v_b| (expect ≤ 1e-8, i.e. every Γ touching
is at least quadratic; this is why V245's enforced node had degree 2.00).
(2) Mermin flux wall: the A5-invariance linear system for bond phases
(φ_i)_{i=1..6} on the 6 star lines — φ_{σ(i)} = φ_i for all 60 permutation
actions (the 6×60 action built from `polyhedral_group([2,3,5])` acting on
the module generators), projective/gauge equivalence modded out by
per-site phases (the 64-site gauge orbit); the dimension of the invariant
flux space (expect 0 — Mermin's theorem measured on this module).
(3) Missing-l=1 wall: the 60×60 adjoint-character matrices of the 45
internal generators restricted to the spin-connection SO(3) stabilizer of
the causal frame, reduced into A5 irreps by inner products of characters;
the multiplicities of each irrep and specifically whether any 3-dimensional
internal vector (l = 1 / T1-type) exists that A5 rotates as spatial vectors
(expect none: the 45 restricted to the spin connection contains only even-l
content — V140's C₂ spectra {45, 85, 117} = l = 2, 3, 4 already imply it;
this experiment makes it exact). (4) The consistency table: the product of
the three walls against V245's measured numbers (Γ degree 2.00, off-Γ
first anisotropic degree 1, stabilizers [2, 2, 60]).

**Acceptance (pre-registered).** Branch `WALLS-CONFIRMED` iff all three
measured: TR defect 0, invariant flux dim 0, no l=1 — then no
single-particle module extension inside the root algebra can put an
isotropic lightlike fermion cone anywhere, and the interacting sector is
the only route (V251–V258 are that route). Else `WALL-<n>-BROKEN` naming
which wall admits a decoration; then this experiment must additionally
build the decoration on the star step and re-measure its unitarity and the
Γ cone degree (pre-registered continuation; the artifact carries the
construction).

---

## V251 — the gauged two-particle sector: construction and Gauss check
`experiments/dynamics_spectrum/V251_gauged_two_particle_sector/` — tier 2,
`residual: false` (sector construction)

**Question.** Does the root carry an exact gauged two-particle sector —
spinor ⊗ one adjoint bond-current flux — with the committed coupling,
Gauss-preserving and bounded? This is the minimal interacting block O9
never had.

**Closed input set.** Per crystal-momentum Q on a 6³ torus grid (chunked
≤ 512): the 128-dim one-particle Bloch space via `supercell.bloch` with
H₁ = i log through the `altcell.effective_hamiltonian` pattern copied from
V242; the 45 adjoint link with the T·L Casimir eigenspaces [−32, −12, 4]
exactly (V126's `C2` eigenspace code copied); the vertex
Γ^A = I_8 ⊗ T^A_16 at strength g₁₀ (V177's named vertex, copied from its
source); the V141 dressed-hop gate as the P→P vertex correction (copied,
zero-parameter per V141's verdict); the committed intertwiner code of
V98/V99 for Q_G. Flux at rest at leading order (V101's measured gap
42.09η ≫ the node scales; the static-flux decision is recorded), with the
sensitivity row ω ∈ {0, ±42.09η} stored in `data.sensitivity` — a measured
window, not a dial. Preamble note: `smatrix.two_particle` (:87, the 't
Hooft on-site 1024² tensor product) is a different object and no
two-particle Bloch spectrum exists anywhere in `cl10/`; new code is
justified by that grep. V98's one-body non-existence (center obstruction)
is recorded in `data.constraints`: the P⊕Q coupled space carries total
Gauss; intermediate Q states need not be singlets.

**Measured quantity.** H₂(Q) = H₁(Q) ⊕ E_flux ⊕ g₁₀·Γ on the 1440-dim
space, for: (a) boundedness — spectrum min/max finite, recorded;
(b) Gauss — total charge Q_G on P⊕Q, ‖[exp(iH₂), Q_G]‖ < 1e-12 on 40
random states; (c) block structure — the P-block restriction equals H₁ to
0.0; (d) unitarity of exp(iH₂) to 1e-13; plus the 6³ × 1440-dim spectrum
(cupy eigh, chunked; wall-time recorded in data as measured cost).

**Acceptance (pre-registered).** Branch `TWO-PARTICLE-SECTOR-CLOSED` iff
(a)–(d) pass with `fitted_inputs == []`; else `SECTOR-BLOCKED` naming the
failed check — a real result: the interaction is then not constructible on
this module and V252–V256 record the floor.

---

## V252 — the dynamical self-energy loop and the loop factor
`experiments/dynamics_spectrum/V252_dynamical_self_energy_loop/` — tier 4,
`residual: true` (O8)

**Question.** V242 showed the free one-particle BZ combinatorics misses
1/(16π) by 7.14×; V177 showed the static Fock dressing is marginal. Does
the dynamical (resolvent) self-energy of the V251 sector supply the loop
factor?

**Closed input set.** Σ(k, z) = g₁₀² Σ_{q,A} Γ^A G₀(k−q, z) Γ^A ⊗ D_A(z)
with G₀ = (z − H₁)^{-1} the one-particle resolvent and D_A(z) =
(z − ε_A)^{-1} on the three Casimir sectors (ε = −32, −12, 4; sensitivity
at ±42.09η recorded), z = E_node ± iδ, δ ladder {1e-2, 1e-3, 1e-4}; q-grid
N ∈ {16, 32} (GPU, eigh of H₁(q) once per q, reused as
Σ_k (u_k†Γu_k)/(z−E_k) — the same resolvent-sum shape as V242's power
sums); evaluation at the 12 cubic M nodes AND the 12 star loci (V244's
stored `orbits[].members`) at the implemented clock (β = g = g₁₀ = 0.2)
and the closed clock (π/4); the Cartan list copied from V209/V242's
`skeleton`. No fitted input; δ is a measurement window recorded per κ.

**Measured quantity.** (1) the loop factor F_int = ‖Cartan content of
Σ‖ / ‖Cartan content of tree‖ at every evaluated node/clock; (2) grid
convergence |F(16) − F(32)|; (3) mpmath certification: 5 interior q-samples
at dps 60 vs 100, direct mp resolvents of the 128-dim matrices at n ≤ 40
words, agreement < 1e-40 between precisions and < 1e-20 vs float64;
(4) the naturalness split: Σ's mass-shift vs cone-slope-shift ratio
(diagnostic, not a contact).

**Acceptance (pre-registered).** Contact `LOOP-FROM-INTERACTION` iff
|16π·F_int − 1| < 1e-6 at any evaluated node/clock with
`fitted_inputs == []`; else `NO-LOOP-EITHER` recording F_int to 15 digits
per class/node/clock and the ratio. Grid convergence |F(16) − F(32)| <
1e-8, else record the finite-size floor as the measurement (V242's
convention). Wording on EMERGENCE: the interaction supplies in dynamical
form what the free combinatorics could not: the loop factor is a property
of the gauged two-particle sector. Wording on NO-LOOP: neither the free
limit nor the dynamical interaction of the only coupled sector the root
carries produces 1/(16π): O8's residual is measured from both sides. A
flux-dispersing upgrade (ω ≠ 0) runs only if F_int lands within the
EMERGENCE band; otherwise the sensitivity row records the window.

---

## V253 — the block-RG window on the dressed sector
`experiments/dynamics_spectrum/V253_block_rg_window/` — tier 2,
`residual: true` (O9/F9)

**Question.** V213 measured the bare tangent anisotropy marginal
(exponent −0.0001); V177's static dressing did not change it. Does the
dressed sector flow?

**Closed input set.** V178's downfold ladder machinery
(`altcell.lowdin_downfold`) applied to the V251 Σ-corrected effective step
on the node subspace (n_low = 4), cuts b = 1, 2, 4, 8, 16; the config
marginal band 1e-2 as the comparator. The V178 artifact keys (`rows`,
`A_full`) are the ladder being extended.

**Measured quantity.** v(b), A(b), m(b) per cut; the log-log exponents
p_v, p_A against the band.

**Acceptance (pre-registered).** Branch `CONTINUUM-WINDOW` iff
p_A < −1e-2 (anisotropy relevant: the coarse-grained cone rotates) —
report b*, the cut where A(b) < 8e-15 would extrapolate, honestly labeled
as an extrapolation with its exponent; else `NO-WINDOW` with
(p_v, p_A, the full b-table).

---

## V254 — the dressed F9 contact at the deepest cut
`experiments/dynamics_spectrum/V254_dressed_f9_contact/` — tier 5,
`residual: true` (F9)

**Question.** Does the dressed effective step satisfy the SME electron
velocity-anisotropy bound at the deepest block-RG cut, measured (not
extrapolated)?

**Closed input set.** V253's artifact at b = 16 (the dressed effective
step), V235's exact normalization and `altcell.cone_anisotropy`, the same
node class V235 used (M, implemented clock) and the star loci; the bound
8.0e-15 from config `bounds.sme_electron_velocity_aniso`.

**Measured quantity.** δv/v at b = 16 measured on the dressed step; the
A(b=16) value and the V253 exponent as supporting numbers.

**Acceptance (pre-registered).** Branch `F9-DRESSED-CLOSED` iff
δv/v < 8.0e-15 measured (not extrapolated) with `fitted_inputs == []`;
else `F9-STILL-FIRED` with the ratio.

---

## V255 — the dressed F2 contact: does the interaction read the decimals?
`experiments/dynamics_spectrum/V255_dressed_f2_contact/` — tier 3,
`residual: true` (F2)

**Question.** With the dynamical loop factor F_int from V252 (and any
C_star clock from V249 if `STAR-CLOCK-OPEN`), does any theory point land
inside both F2 bands?

**Closed input set.** F_int from V252's artifact; C_star clocks from
V249's artifact if OPEN; the closed clocks (π/4, π/2), the implemented
0.2; the ladder coefficient rule 3(16π·F)^10 vs 2.0765641 (V108's rule,
read from its artifact for the rule shape; ladder reading coefficient =
3(16πF)^10 per `cl10/config.yaml:87` and V108); `unification.alpha_std(g)`
vs 0.04 with `alpha_10 = 0.04`; V233's band for the second contact.

**Measured quantity.** The 2×2 recomputation at every theory point: the
ladder coefficient vs 2.0765641 and α_std(g) vs 0.04; which functional was
used recorded per row. β-independence control: the running reads
`coupling.angle()` only — one numeric assertion as V247 did.

**Acceptance (pre-registered).** Branch `F2-INTERACTION-READS` iff some
theory point lands inside both bands (1e-3 / V233's) with
`fitted_inputs == []`; else `F2-FIRED-X4` with the full table.

---

## V256 — the ladder rungs at the interacting loop factor
`experiments/dynamics_spectrum/V256_ladder_rung_reconstruction/` — tier 4,
`residual: true` (O8)

**Question.** V197's rung arithmetic at F_int: the ladder's four scales
(rung 3 singlet A, rung 10 singlet B, rung 17 neutrino) recomputed with
loop = F_int^(1/10) — do the three rungs land jointly within V197's bands?
This is the joint falsifier guarding against a single-coincidence loop
factor.

**Closed input set.** V197's artifact (rung arithmetic and bands),
V252's F_int, the config ladder hierarchy coefficient 2.0765641.

**Measured quantity.** The three rung values and their band residuals at
loop = F_int^(1/10).

**Acceptance (pre-registered).** Branch `LADDER-CONSISTENT` /
`LADDER-BROKEN` plus the residual table.

---

## V257 — the interacting ledger
`experiments/algebra_root/V257_interacting_ledger/` — tier 6,
`residual: true` iff anything landed fitted

**Question.** The provenance of every number the interacting programme
produced: η twist values (algebraic only if on C_star, else
never-realized), F_int (measured), dressed cone/flux scales (measured),
the sector coupling (algebraic: g₁₀, the V177/V141 zero-parameter gates),
any C_star clock (algebraic if closed). The Table-6-class counts per
quantity.

**Closed input set.** Only the V249–V256 verdict artifacts
(`cl10.paths.exp`); the published Table 6 rows it compares against. No
new number is measured.

**Measured quantity.** The counts (algebraic, fitted, measured-input) per
quantity and the reclassification flags for §15/Table 6.

**Acceptance (pre-registered).** Branch string `INTERACTING-LEDGER`
(informational). `residual: true` iff anything landed fitted. If V251
landed `SECTOR-BLOCKED`, the would-be sector numbers are recorded as
never-realized and add no fitted parameter.

---

## V258 — the falsifier conjunction of the interacting programme
`experiments/algebra_root/V258_falsifier_conjunction/` — tier 5,
`residual: true`

**Question.** Does any single theory point satisfy (F9 dressed ∧ F8 ratio
∧ F2 both bands ∧ O8 ladder) simultaneously? The programme's roll-up as
an experiment.

**Closed input set.** Only the V249–V257 artifacts read through
`cl10.paths.exp`, plus the F8-ratio input computed here: the dressed
fermion cone slope vs V208's cone machinery re-evaluated on the dressed
effective step (`cone_generator.json` re-diagonalized on H_eff at b = 16,
the 8-dim cone-subspace check). The ratio, not the absolute, is the
GW170817 contact. No recomputation of any other artifact's numbers.

**Measured quantity.** The conjunction table (per experiment × per theory
point) and the F8 dressed ratio.

**Acceptance (pre-registered).** Branch `CONJUNCTION-ALIVE` (a named
point; the falsifiers become decidable by a next programme) or
`CONJUNCTION-DEAD` (the per-experiment table is the artifact; F2/F8/F9/O8
become theorems of the whole construction, both single-particle and
interacting).

---

## Roll-up consequences of this programme (fixed before any run)

- The suite recount (178 → 188) propagates to `README.md`,
  `doc/reference.md`, `experiments/README.md` (counts, catalogue tables,
  claim mapping), availability lines and Appendix B, whether or not any
  contact succeeds.
- Any dial value used off C_star reclasses to `fitted` in Table 6 and the
  falsifier stays fired. The hard-tier denominator is whatever the runner
  measures over the new tier-0–3 experiments (5 of them pre-registered:
  V249–V251, V253, V255).
- The manuscript gains a §21 (the interacting continuum programme) and the
  woven verdicts in the abstract, §12.1, §15/§16/§17 tables and the
  comparison paragraphs; no existing section is renumbered. Both outcomes
  of every branch are drafted above; only the landed wording is applied.
- If V251 lands `SECTOR-BLOCKED`, V252–V256 run on the P-only dressed
  system (V177's static vertex + V178 downfold) and the branch strings
  carry `sector-floor`; V257 books the would-be sector as never-realized.
- No tag, no merge: the branch `founding-repair` is pushed for the
  reader's decision.


---

# Pre-registrations: the gauge-invariant physical sector programme (V259–V268)

The interacting programme (V249–V258) localized the failure: V98 proved
the Gauss law admits **no one-particle physical state** (the spinor centre
charge cannot be screened by adjoint flux), and V251 showed dynamically
that the leak sits in the one-particle block itself. Yet every fired
falsifier (F2, F8, F9, O8, O10) measures **one-particle Bloch modes** —
states the theory's own Gauss law says are not physical. This programme
re-measures the falsifiers' objects on the **physical gauge-singlet
sector**: composite (matter + flux) states, their dispersion, gaplessness,
chirality/statistics, and their scales. The end state is a
pre-registered, artifact-backed verdict on whether the compositeness
route rescues any falsifier — either the route yields a gapless physical
cone/scale (contacts reopen, the theory is rewritten on the physical
sector) or it is falsified too (the falsifiers stand as theorems of every
construction, and the manuscript says so).

**Ground rules.** All ten experiments use the **same closed construction**:
V251's fiber machinery (128-dim spinor Bloch block + 1440-dim 32 × 45
adjoint-link flux block, the V177 named vertex at strength g₁₀ = 0.2, the
V126 Casimir tower with the omega rows measured by V251's artifact). The
shared per-q construction — the **V259 kernel recipe** — is fixed now:
build H₂(q) (1568 × 1568) exactly as V251's `build_H2` does at each q;
build C₂ = Σ_a Q_G^a Q_G^a from V251's 45 total-charge operators;
`w, V = eigh(C₂)`; kernel basis `W = V[:, w < 1e-9]` (tolerance 1e-9 —
the same convention as "exactly gapless below 1e-9" elsewhere in the
repo); projected composite Hamiltonian H_k(q) = W† H₂(q) W; composite
bands = eigh(H_k). Zero fitted parameters anywhere: every number is
measured or read from a committed artifact (grid N = 6, the omega rows,
the gap scales); any dial used outside the V251 artifact conventions
reclasses to fitted and the falsifier stays fired. Every artifact gains
`wall_time`. GPU 1 (`CUDA_VISIBLE_DEVICES=1`), cupy backend with the
V254 memory-pin fallback (host numpy under 3.5 GiB free). Branch
`founding-repair` continues; no merge, no tag; push at the end. Tier
assignment: V259 → 1, V260 → 2, V261 → 5, V262 → 3, V263 → 2, V264 → 3,
V265 → 5, V266 → 2, V267 → 6, V268 → 5. New hard-tier (0–3) experiments:
6 (V259, V260, V262, V263, V264, V266); the runner measures the
denominator, report, never adjust. A fired contact is a PASS with
`residual: true`.

---

## V259 — the singlet kernel census of the coupled fiber
`experiments/dynamics_spectrum/V259_singlet_kernel_census/` — tier 1,
`residual: false` (structural, hard tier)

**Question.** Does the coupled two-particle fiber carry gauge singlets at
any momentum, and how many?

**Closed input set.** The V259 kernel recipe over V251's 6³ q-grid (216
points): V251's charges, vertex, tower and `build_H2` (copy-modified in
the experiment file), g₁₀ = 0.2, the omega convention read from the V251
artifact (`data.sensitivity.omega_rows`); the kernel tolerance 1e-9. No
other number.

**Measured quantity.** Per-q kernel dim; the total singlet count; the
q-points carrying ≥ 1 singlet, with their little groups via
`cl10.continuum` (V244's import line, `K.molien_counts`) computed at each
such q. Regression checks (pre-registered): T0 — V251's tower
multiplicities {32, 288, 1120} on the 1440 block reproduced to 1e-9; T1 —
C₂ commutes with H₂'s flux block to < 1e-9 on 40 random states (seed
0x51, V251's convention).

**Acceptance (pre-registered).** Branch `COMPOSITE-KERNEL`
(informational) iff kernel dim > 0 at ≥ 1 q (record dim and loci); else
`KERNEL-EMPTY` — the cascade rule of the roll-up fires and V260–V266 are
recorded NEVER-REALIZED.

---

## V260 — the composite gap ladder on the singlet kernel
`experiments/dynamics_spectrum/V260_composite_gap_ladder/` — tier 2,
`residual: false` (structural, hard tier)

**Question.** Is the composite projected band gapless anywhere in the
zone, at the published N = 6 grid?

**Closed input set.** The kernel basis per q from the V259 recipe;
composite band E₀(q) = lowest eigenvalue of H_k(q) above the
q-independent ground offset (the offset subtraction recorded explicitly
in `data.conventions`); V251's grid and omega rows (the ± omega rows
booked as the coupling sensitivity); the gaplessness tolerance 1e-9; V235's
velocity normalization dvv = (v_max − v_min)/v_max along the 12 star rays
at dq = 1e-3 (V254's measurement loop); V129's committed gap read from
its artifact.

**Measured quantity.** min E₀ over the grid; the set of q with E₀ < 1e-9
(gapless loci); the velocity tensor at those loci by V235's
normalization; the ± omega-row minima; V129's gap as the comparison.
Cross-check (new-behavior proof): at q = 0 the composite kernel dimension
and the C₂ kernel eigenvalue count must equal the singlet-copy count
computed tower-sector by tower-sector (assert equality < 1e-9 inside the
script — a wrong projection shows up here).

**Acceptance (pre-registered).** Branch `COMPOSITE-GAPLESS` iff ≥ 1 locus
with E₀ < 1e-9; else `COMPOSITE-GAPPED` recording min gap and the ratio
to V129's committed gap.

---

## V261 — the composite F9 contact
`experiments/dynamics_spectrum/V261_composite_f9_contact/` — tier 5,
`residual: true` (F9)

**Question.** At composite gapless loci, does the velocity anisotropy
meet the SME ceiling 8.0e-15 (`_config.CFG.bounds.sme_electron_velocity_aniso`)?

**Closed input set.** V260's loci; the same 12-ray finite-difference
gradient V254 uses (dq = 1e-3, V235 normalization); the bound read from
config, not retyped.

**Measured quantity.** dvv per locus; the best value and its ratio to the
bound.

**Acceptance (pre-registered).** `F9-COMPOSITE-CLOSED` iff best dvv <
8.0e-15; `F9-COMPOSITE-FIRED` with the ratio; `F9-COMPOSITE-NO-LOCUS` if
V260 landed GAPPED or V259 KERNEL-EMPTY (recorded, exit 0).

---

## V262 — the composite chirality/statistics census
`experiments/dynamics_spectrum/V262_composite_chirality_statistics/` —
tier 3, `residual: false` (structural, hard tier)

**Question.** What are the little-group and statistics content of the
composite gapless loci — could they carry the electron?

**Closed input set.** At each gapless locus, V244's little-group/Molien
census; the statistics census: the singlet kernel is built from an even
number of centre-charged matter quanta (V98's theorem), so composite
centres are integer-spin unless an odd-matter sector exists — measured,
not assumed: decompose the kernel basis under the internal Spin(10)
Cartan (the V99/V100 Cartan-eigenspace machinery from
`experiments/gauged_layer/V99_three_site_chain/` and `V100_vertex_correction`)
and record the parity of the matter weight per kernel state.

**Measured quantity.** Per-locus little-group Molien counts; per kernel
state the Cartan matter weight and its parity.

**Acceptance (pre-registered).** `CHIRAL-COMPOSITE` iff a gapless locus
has half-integer matter-parity content (impossible by V98 — if it fires,
the theory is rebuilt); `BOSONIC-COMPOSITE-ONLY` iff every kernel state
is even-parity (the centre obstruction re-derived at the composite
level); else the NO-LOCUS cascade.

---

## V263 — the composite loop scale (O8)
`experiments/dynamics_spectrum/V263_composite_loop_scale/` — tier 2,
`residual: true` (O8, hard tier)

**Question.** Does the composite sector generate a dimensionless scale
playing the role of the loop factor?

**Closed input set.** The composite gap at its minimum (V260 recipe); the
V101 bond-current ground-state scale E₀ read from the
`V101_bond_current_ground_state` artifact (the implementer reads the field
from that JSON's `data`, never retypes a digit); V129's gap as the second
candidate from its artifact.

**Measured quantity.** Both candidate quotients
F_phys = (composite gap at minimum)/(V101 E₀ scale) and (composite gap)/
(V129 gap); their deviations from the target 1/(16π) = 1.9894e-2.

**Acceptance (pre-registered).** Bands fixed now, not at run time:
`SCALE-LOOP-CANDIDATE` iff either candidate lies in [1.5e-2, 3.0e-2]
(factor 1.5 of 1/(16π)); else `NO-SCALE-ROUTE` recording the closest
value.

---

## V264 — the composite desert content (F2)
`experiments/dynamics_spectrum/V264_composite_desert_content/` — tier 3,
`residual: true` (F2, hard tier)

**Question.** Does the composite singlet multiplicity table supply desert
content that V237's threshold rule can realize?

**Closed input set.** The composite multiplicity table (decomposition of
the kernel under the Cartan, from V262's machinery); V237's threshold
rule — the exact rule code copied from
`experiments/dynamics_spectrum/V237_desert_threshold_matching/` (thresholds
and the (1,2)-shift demand, verbatim inequality); the M_X log-spacing
from V237's artifact.

**Measured quantity.** The (1,2)-shift demand realized by the composite
content and V237's own inequality margin.

**Acceptance (pre-registered).** `F2-COMPOSITE-CLOSED` iff the (1,2)-shift
demand ≤ the log-spacing to M_X (V237's own inequality, verbatim); else
`F2-COMPOSITE-FIRED` with the shift ratio; else the NO-CONTENT cascade.

---

## V265 — the composite F8 contact
`experiments/dynamics_spectrum/V265_composite_f8_contact/` — tier 5,
`residual: true` (F8)

**Question.** Does the composite sector carry a spin-2 channel
propagating at c = 1?

**Closed input set.** V229's T2g stress vertex (the `Emat`,
`spin2_generators` functions of
`experiments/gauged_layer/V229_diagonal_bond_quintet/`) sandwiched
between composite ground and first-excited kernel states along the 15
star axes; the slope of the linear part gives the composite c_gw; the
band |c_gw − 1| < 1e-15 is the one V258 used.

**Measured quantity.** The composite c_gw and its deviation from 1.

**Acceptance (pre-registered).** `F8-COMPOSITE-OK` iff |c_gw − 1| <
1e-15; else `F8-COMPOSITE-FIRED` with the ratio.

---

## V266 — finite-size convergence of the gap verdict
`experiments/dynamics_spectrum/V266_finite_size_convergence/` — tier 2,
`residual: false` (structural, hard tier)

**Question.** Is the V260 gap verdict finite-size or thermodynamic?

**Closed input set.** The kernel-recipe gap ladder on q-grids N = 4, 6, 8
(the fiber stays 1568-dim — cheap; only the q-grid grows); the same
tolerances.

**Measured quantity.** min gap per grid; the trend exponent by least
squares on log E_min vs 1/N.

**Acceptance (pre-registered).** `GAP-CONFIRMED` iff the gaps agree across
N = 4, 6, 8 to 1e-6 and all > 1e-9; `GAPLESS-EMERGES` iff E_min crosses
below 1e-9 at N = 8 with a negative trend exponent; `GRID-NEUTRAL`
otherwise (values recorded).

---

## V267 — the physical ledger
`experiments/algebra_root/V267_physical_ledger/` — tier 6,
`residual: true` iff anything landed fitted

**Question.** The provenance of every number the physical-sector
programme produced: algebraic, measured-input, fitted, never-realized,
with the Table-6-class counts.

**Closed input set.** Only the artifacts V251, V259–V266 (all read via
the `art()` pattern); the published Table 6 rows it compares against. No
new number is measured.

**Measured quantity.** The counts per provenance class and the
reclassification flags for §15/Table 6. `fitted_inputs` must be []; if
anything is fitted, `residual: true` and the Table-6 flag is recorded
(the V257 bookkeeping shape).

**Acceptance (pre-registered).** Branch string `PHYSICAL-LEDGER`
(informational). If V259 landed KERNEL-EMPTY, the sector is booked as
never-realized and adds no fitted parameter.

---

## V268 — the final conjunction of the physical sector
`experiments/algebra_root/V268_final_conjunction/` — tier 5,
`residual: true`

**Question.** Does any single theory point satisfy the composite F9 ∧
composite F8 ∧ composite F2 ∧ composite loop scale simultaneously? The
programme's roll-up as an experiment (the V258 code shape).

**Closed input set.** Three theory points: the implemented clock
composite point (the V260 loci), the N = 8 grid point (V266), the
coupling V251's omega (V251's artifact sensitivity rows). Four legs:
composite F8 (V265), composite F9 (V261), composite F2 (V264), composite
loop scale (V263). Only artifacts read via `art()`; no recomputation of
any other artifact's numbers.

**Measured quantity.** The conjunction table (per leg × per theory point)
and the per-leg deviations.

**Acceptance (pre-registered).** `COMPOSITE-CONJUNCTION-ALIVE` iff any
point passes all four legs; else `CONJUNCTION-DEAD` with the per-leg
deviations recorded.

---

## Roll-up consequences of this programme (fixed before any run)

- The suite recount (188 → 198) propagates to `README.md`,
  `doc/reference.md`, `experiments/README.md` (counts, catalogue tables,
  claim mapping), and the manuscript's Appendix B and availability note,
  whether or not any route succeeds.
- If V259 lands KERNEL-EMPTY (no q-point carries a singlet), V260–V266
  are recorded NEVER-REALIZED by the pre-registered cascade rule and
  V267 books the sector as never-realized; V268 fires CONJUNCTION-DEAD.
- Any dial used outside the V251 artifact conventions reclasses to fitted
  and the falsifier stays fired. The hard-tier denominator is whatever
  the runner measures over the new tier-0–3 experiments (V259 T1, V260
  T2, V262 T3, V263 T2, V264 T3, V266 T2 pre-registered).
- The manuscript gains a section 22 (the gauge-invariant physical sector
  programme); woven verdicts in the abstract, section 12.1, Tables
  6/7/8/9; no existing section renumbered. Both branch outcomes are
  drafted; only the landed wording is applied.
- No tag, no merge: `founding-repair` is pushed for the reader's
  decision.


---

# The dynamical-link refonte: pre-registrations V269–V274
(fixed 2026-10-26, before any run; the refonte question is the reader's own:
replace the frozen bond gate by a quantum link with its own dynamics)

**The theory under test (the object of the refonte).** One vertex: the
minimal neutral matter system of the root - a spinor pair (32 x 32, 1024
states, committed charges T^a_pair = T^a (x) I + I (x) T^a from V79's
plane_generators) - coupled to ONE DYNAMICAL Spin(10) link: a QLM-truncated
link module L = 1 (+) 45 (+) 210, the fusion truncation of 45 (x) 45 (a
declared truncation - a convention, not a fit), its generators built by
Casimir projector from the committed adjoint operators (V79), never typed.
The Gauss algebra is the committed C = T_pair + L (V79/V126 conventions).
The dynamics (zero fitted parameters): rotor H_rot = (1/2) C2_link (rotor
unit I = 1/2: a declared unit convention, lattice units - the only scale
choice) plus the committed-clock exchange H_ex = g10 sum_a J^a (x) L^a,
g10 = 0.2 the implemented clock. The physical composite sector = the
kernel of C2_tot = sum_a (T^a_pair + L^a)^2 on 1024 x 256 = 262144,
computed matrix-free (Lanczos on S = sum (Q^a)^2, tolerance 1e-9, fixed
seed). Predicted kernel dim from the committed fusion arithmetic
(V98/V99 controls, re-measured live in V269): m_1 = m_45 = m_210 = 2 in
the pair, N_r = 1 in L, so dim ker = 2 + 2 + 2 = 6. This is the first
theory version in which the neutral composite can be DYNAMICAL.

## V269_dynamical_link_module (tier 1, hard, construction)

**The question (pre-registered).** Does the truncated dynamical-link
module 1 (+) 45 (+) 210 exist inside the root's own so(10): the Casimir
projector of 45 (x) 45 isolating exactly dims [1, 45, 210] with three
distinct Casimir eigenvalues, and do the 210-block generators close the
committed algebra?

**The object (closed input set).** Ls (45) adjoint link operators from
V79; C2_2025 = sum (L^a (x) I + I (x) L^a)^2 on 45 (x) 45 by dense eigh;
projectors P_r = eigenvector spans binned by eigenvalue; 210-algebra check
P([L^a, L^b])P = i f^{abc} P L^c P on the block, band 1e-9; Jacobi/
antisymmetry of f re-verified live (V79 control); pair multiplicities
measured live: C2_pair spectrum on 1024, the degeneracy of the eigenvalues
0 / 64 / lambda_210 must be m_1 = m_45 = m_210 = 2 (V98/V99 committed -
recomputed, not read).

**Measured (pre-registered).** The distinct Casimir eigenvalues of the
2025-module (expect three at least in {1,45,210} sectors; lambda_210
recorded as a measured number); the block dims; the 210-algebra defect;
the pair multiplicities.

**Acceptance (pre-registered).** LINK-MODULE-OK iff dims exactly
[1, 45, 210] (Gram defect < 1e-9), three distinct eigenvalues, 210-algebra
defect < 1e-9, multiplicities (2, 2, 2); LINK-MODULE-DEFECT otherwise.

## V270_neutral_composite_sector (tier 1, hard, the refonte's first verdict)

**The question (pre-registered).** Does the dynamical-link atom carry
gauge-neutral composites - is ker C2_tot non-empty - and is its sector
structure the predicted [2, 2, 2] under the link Casimir?

**The object (closed input set).** V269's module; Q^a_tot = T^a_pair + L^a
on 262144; S = sum_a (Q^a_tot)^2 matrix-free (kron-reshape matvec on
(1024, 256)); Lanczos (fixed seed 0x6C, 600 steps max, restarts on
stalling) counting eigenvalues < 1e-9; sector labels from the Gram of
C2_link restricted to the kernel basis; entanglement S(pair|link) per
kernel vector (reshape-SVD).

**Measured (pre-registered).** Kernel dim (predicted 6); sector split
(predicted {0: 2, 64: 2, lambda_210: 2}); the entropies: the r = 45 kernel
vectors must sit at the committed intertwiner value ln 45 = 3.80666...
(V98 control, band 1e-6); the r = 1 vectors at 0.

**Acceptance (pre-registered).** COMPOSITE-ALIVE iff dim = 6 with the
predicted split and the ln-45 control inside 1e-6; KERNEL-STILL-EMPTY iff
dim = 0; UNEXPECTED-KERNEL otherwise, recorded (dim, split, entropies).

## V271_composite_dynamics (tier 2, hard)

**The question (pre-registered).** Once the composite exists, does it
move: is the Gauss-preserving rotor + exchange Hamiltonian non-diagonal on
the kernel (dressing), and what are the composite's mass, excitation
splittings and ground entanglement?

**The object (closed input set).** H = (1/2) C2_link (x) I_pair + g10
sum J^a (x) L^a on the atom (262144, matrix-free); kernel basis from V270
recomputed; H restricted to the kernel (6x6) eigh; the full-space first
eigenpairs by Lanczos seeded from the kernel vectors (Krylov 800); the
excited composite states = kernel-projected excited eigenvectors whose
C2_tot expectation < 1e-9.

**Measured (pre-registered).** [H, Q^a_tot] defect (must be < 1e-12 - the
new theory's Gauss law, measured not claimed); the kernel spectrum
(E0 mass, Delta splitting); the dressing fraction (1 - overlap of the H
ground state with the pure rotor sector); the ground entropy S(pair|link);
the entropy deviation from ln 45.

**Acceptance (pre-registered).** COMPOSITE-DYNAMICS-ALIVE iff Gauss defect
< 1e-12, the restricted H is not proportional to I (mixing > 1e-6), and
the ground entropy lies in (0, ln 45]; COMPOSITE-FROZEN if H is diagonal
on the kernel; GAUSS-BROKEN if the defect fails.

## V272_dynamical_chain_block_rg (tier 2, hard)

**The question (pre-registered).** Two dynamical-link composites sharing
one link: is the chain neutral (predicted dim 12), and what effective
coupling does the shared link generate between the two atoms (the
candidate scale of the new theory)?

**The object (closed input set).** Chain P1 (x) l (x) P2 with Gauss
T1 + L = 0, T2 - L = 0: chain kernel constructed by singlet contraction
of two V270 atom-kernel bases over the shared link (Gram to 1e-9;
predicted dim sum_r m_r^2 = 12); chain H = (1/2) C2_l + g10 (J1.L + J2.(-L))
in the contracted basis (matvec on 1024 x 256 and 1024 x 1024 x 256 blocks);
interaction energy E_int = E_chain(ground) - 2 E_atom(ground); effective
coupling rho = E_int / g10.

**Measured (pre-registered).** Chain kernel dim (predicted 12, Gram
defect); the atom-vs-chain ground energies; E_int; rho; the chain Gauss
defect (both vertices, < 1e-12).

**Acceptance (pre-registered).** CHAIN-ALIVE iff dim = 12 and Gauss defect
< 1e-12 (else CHAIN-EMPTY / GAUSS-BROKEN); rho recorded either way.

## V273_loop_factor_new_theory (tier 2, hard, O8-new)

**The question (pre-registered).** Does the dynamical link BIRTH the
missing loop factor: does any dimensionless composite-sector ratio land
in the ladder band 1/(16 pi)?

**The object (closed input set).** Candidates (fixed here): F1 = rho
(V272 link-mediated coupling); F2 = Delta / E0 (V271 splitting over mass);
F3 = the second-order link-mediated exchange between chain composites,
J2/g10 with J2 = -<g| H_ex Q (E_g - H_atom)^{-1} Q H_ex |g> computed by
Krylov (Lanczos-Stieltjes) on the atom resolvent seeded from the chain
ground contraction, Krylov 800, seed 0x6C. Target 1/(16 pi) = 1.9894e-2,
band [0.015, 0.03] (as V263's, same ladder convention).

**Measured (pre-registered).** The three candidates; the closest
|16 pi F - 1|; the Krylov convergence residual (60 vs 100-digit control
on the resolvent quadratic form - mpmath, 5 fixed interior probes).

**Acceptance (pre-registered).** LOOP-NEW-BORN iff any candidate inside
the band; LOOP-NEW-UNBORN otherwise with the closest miss recorded.

## V274_new_theory_ledger_conjunction (tier 6 / tier 5)

**The question (pre-registered).** The ledger of the refonte (provenance
of every booked number) and the new theory's own conjunction: three theory
points (g10 = 0.2 implemented; g = pi/4 closed clock; g = 1/5 rung clock)
x four legs (COMPOSITE-ALIVE, DYNAMICS-ALIVE, CHAIN-ALIVE, LOOP-BORN).

**Acceptance (pre-registered).** Ledger NEW-THEORY-LEDGER with counts by
class (the rotor unit and the truncation {1,45,210} book as CONVENTIONS -
declared, named, not fitted); the conjunction branch NEW-CONJUNCTION-ALIVE
iff one point passes all four legs, else NEW-CONJUNCTION-DEAD, per-leg
table recorded.

## Roll-up consequences of this programme (fixed before any run)

- The suite recount (198 -> 204) propagates to README.md, doc/reference.md,
  experiments/README.md whether any route succeeds.
- V270 decides the refonte: KERNEL-STILL-EMPTY would falsify the dynamical
  link on the same ground as V259 (the centre obstruction would then be a
  theorem of the LINK too, not only of the frozen gate); COMPOSITE-ALIVE
  opens V271-V273 as the new physics; UNEXPECTED is recorded, not curated.
- The truncation {1,45,210} and the rotor unit book as conventions
  (declared before the run); nothing in this programme may be refitted to
  land in the loop band; any dial touched reclassifies to fitted and
  O8-new stays fired.
- Only the six new experiments run (the reader's instruction); the
  committed 198-experiment suite is NOT rerun by this programme.

---

# Off-kernel resolvent and certified ladder: pre-registrations V275–V276

## V275_offkernel_exchange_resolvent (tier 2, hard, O8-new)

**The question (pre-registered).** Does the second-order link-mediated
exchange become nonzero when the resolvent is grown from the FULL
off-kernel exchange span (the 540 individual vectors A_a k_i, C_a k_i,
rank 360, recorded H-non-invariant by V273) instead of the 6 channel
sums whose Krylov closed inside ker H?

**The object (closed input set).** V273's atom verbatim (NATOM =
1024 x 256 = 262144, channel-adapted kernel ker dim 6 via S_of Lanczos
500 seed 0x6C, H_atom_of, AKv[a][i], CKv[a][i], Q_of). Seeds: the 540
exchange vectors individually, kernel-projected, orthonormalized (rank
recorded, expected 360). Block Krylov under H_atom_of, every step
projected off ker and off the running basis (V273 machinery), cap
MAXV = 1200, STEP_CAP = 40. If closed: F5 = J2/g10 by V273's T3
four-term assembly with the 60/100-digit mpmath spectral control
verbatim; min |E_g - spec(HU)| recorded. Target 1/(16 pi) = 1.9894e-2,
band [0.015, 0.03] (same as V273/V263).

**Measured (pre-registered).** rank of the 540-span (kernel-projected,
SVD tol 1e-9; the check books the measured rank - V273 measured 360
for the unprojected set and the note was an observation, not a gate);
growth table, closure flag, invariant residual; F5 (exact if closed;
truncated Rayleigh-Ritz sandwich if capped, labelled not-bookable);
|16 pi F5 - 1|; mpmath 60-vs-100 control.

**Acceptance (pre-registered).** OFFKERNEL-LOOP-BORN iff closed and F5
inside the band; OFFKERNEL-UNBORN iff closed and outside (closest miss
recorded); OFFKERNEL-GROWS iff still growing at MAXV (growth table is
the evidence, no bookable number). The experiment passes in every
branch except a crash or an atom-kernel block: the physics verdict is
the branch, not the suite gate (the V273 convention).

## V276_certified_ladder_208 (tier 4, O8)

**The question (pre-registered).** Does the certified closed branch
(V200-V204 machinery) support the published hierarchy coefficient
without the n = 10 read, and is the certified crossing n* = 208
arithmetically derived?

**The object (closed input set).** V203's hop_machinery(G_CLOSED) +
loop_ratios(G, cartan, Ma, ev, G10_REF, 260) imported verbatim (same
import path V204 uses); artifacts V204 (machinery, rungs, coefficient),
V108 (coefficient_predicted_3_F10_10, coefficient_observed, F_10_root),
V197 (ladder_rows). All float64 machinery + mpmath (60 dp) arithmetic
continuations only.

**Measured (pre-registered).** T0 reproduction: rho, crossing_n, odd-n
max vs the V204 artifact. T1 ladder read: coef_pred = 3*(16 pi F)^10
evaluated at F := F_10_open_branch (0.019186391026184392), F :=
F_cross (0.01968515418112789), and F := geometric continuation
F_10_geom from (rho, F_cross) per the measured F_n direction (the
continuation direction is read from the loop_ratios construction and
booked explicitly); each against coefficient_observed = 2.076564111173832.
T2 208-arithmetic: mpmath PSLQ over Q on {1, rho, ..., rho^4} (degree
<= 4, height cap 1e4, 60 dp); geometric crossing n*_geom from
(F_10_open_branch, rho) along the booked direction - record
|n*_geom - 208| with the float64 uncertainty band of rho propagated;
factor census 208 = 2^4 * 13, divisibility by the closed-clock orders
(8 | 208) and by the V197 rungs.

**Acceptance (pre-registered).** T0 fail (any of |rho - rho_art| >
1e-9, crossing_n != 208, odd max > 1e-12) -> MACHINERY-DRIFT: halt the
programme, report, do not interpret anything downstream. Ladder leg:
LADDER-CERTIFIED iff the F_cross read is within 2 % of
coefficient_observed; else LADDER-N10-CONVENTION (both readings
recorded). 208 leg: CROSSING-ALGEBRAIC iff PSLQ finds an exact
relation AND |n*_geom - 208| < 1e-6; else CROSSING-MEASURED.

## Roll-up consequences of this programme (fixed before any run)

- The suite recount (204 -> 206) propagates to README.md,
  experiments/README.md, doc/final_publication.md Appendix B,
  manuscript section 24, whether any route succeeds.
- The reclassification clause: iff V275 in {OFFKERNEL-UNBORN,
  OFFKERNEL-GROWS} and the V276 208-leg is CROSSING-MEASURED, the
  ladder base 1/(16 pi) is reclassified from derived to a declared
  named input of the scale sector (status CONVENTION in the ledger,
  alongside the O9 measurement input), and no further O8-mechanism GPU
  experiments are booked.
- Nothing in this programme may be refitted to land in the loop band;
  any dial touched reclassifies to fitted.
- Only the two new experiments (V275, V276) interpret the result; the
  committed 204-experiment suite is rerun once by the roll-up to
  re-gate the whole (206 expected).
- V276 T0 is ground truth for the published V204 audit: any drift
  halts interpretation of everything downstream.
