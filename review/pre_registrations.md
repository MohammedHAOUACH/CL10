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

---

## Retirement note (2026-10-27, after the V276 roll-up; reader instruction)

The reader ordered the removal of experiments whose routes later work
showed to lead nowhere. The following O8 dead-end mechanism
experiments were retired, scripts and booked artifacts alike (git
history keeps the measurements and their dates verbatim):

- V231_star_loop_functional_certified (star-module loop census: no
  certified loop factor);
- V232_insertion_amplitude_search (no closed-form insertion amplitude
  anywhere on the closed lattice);
- V242_bz_loop_factor_census and V243_bz_hierarchy_contact
  (Brillouin-zone route: NO-LOOP-FACTOR / HIERARCHY-UNREAD);
- V263_composite_loop_scale (the composite loop numerator does not
  exist);
- V273_loop_factor_new_theory (LOOP-NEW-UNBORN; superseded by V275's
  larger closed off-kernel resolvent).

Consumers were rewired to surviving inputs, and every rewired script
was re-run with its branch unchanged: V247 and V255 inline the
V209-verbatim cell-ladder functional (V209's own artifact keeps the
published pi/4 cell values as the copy-control target); V248 dropped
the BZ ledger row (counts rebooked); V267 re-reads the scale
denominators from V101/V129/V251 directly and re-verifies the
V101-gap-times-g10 = omega identity itself; V268 dropped the
never-realized loop-scale leg (three legs, CONJUNCTION-DEAD unchanged);
V274 re-declared the V273 block-Krylov machine verbatim (it already
carried it) and its T0 artifact set is now V269-V272 plus V275.

Suite: 200 experiments. The full 206-experiment re-gate pre-registered
for the V275/V276 roll-up was cancelled by the same instruction and
NOT executed; the committed artifacts of the remaining 200 are
untouched by this retirement, and no published claim survives that
cited a retired artifact except with an explicit retired-with-date
annotation.

The re-gate was subsequently EXECUTED (reader instruction to finish the
open items, after the O8 ledger cutover put the ladder base in the
parameter audit): the 200-experiment run started 2026-10-27 18:12 on
GPU1 (the run log prints experiments: 200 - the retired six are gone
from auto-discovery; cupy pinned, JOBS=3) and booked its last verdict
2026-10-28. The pool fold was 197/200 PASS with the hard tier True and
three FAILED: V274_new_theory_ledger_conjunction CRASHED (cupy
OutOfMemoryError on a 4.29 GB allocation while three workers contended
for the card), and V257_interacting_ledger + V258_falsifier_conjunction
booked FAIL on their artifact-reading checks - they run (algebra_root)
before the dynamics_spectrum pair V247/V255 is rewritten, against the
ec86d53-era standalone-commit state, in which the committed V247 pi/2
row carried a 1.78e+15 copy-control deviation (V255 T0 False); this
gate re-ran V247 bit-for-bit at the 2026-10-26 gate value (deviation
0.0e+00, V255 re-booked True). Per the run_experiments.sh crash
contingency the three were re-run standalone on the idle card: V274
completed in 9945 s and re-booked the roll-up verdict NEW-CONJUNCTION-
DEAD at all three points (Y/Y/Y N, Y/Y/N N, Y/Y/Y N), V257 and V258
re-booked PASS against the fresh artifacts. The summary was then
rebuilt from on-disk verdicts with run_validation.py --report-only:
200/200 PASS, backend-wide pass True, hard tier (0-3) True, n=200.
The 206-experiment gate above stays as the historical record of what
was booked when.

## Final-state document integration note (2026-10-28, executed)

No experiment was pre-registered or run by this step; it records the
final-state integration of the publication.

- final_publication.md: Sections 20-24 (founding-repair, interacting
  continuum, physical sector, dynamical link, off-kernel resolvent and
  certified ladder) were removed as programme histories; their final-state
  verdicts were redistributed into the structure: 8.2 (dressed F9 contact,
  V253/V254), 8.4 (the three walls, V250), 8.5 (enforced isotropy and the
  star clock, V244/V245/V249), 9 (two-dial closure set, V247), 10.1
  (physical-sector emptiness, V259-V268), 10.5 (new subsection: dynamical
  link, V269-V276), 11.6 (cone breach and decoration, V251), 12.1
  (certified ladder and the PSLQ crossing census, V276), 13.1 (two-dial
  and dressed-run census, V246/V255/V264), and the abstract, Tables 1/6/
  7/8/9 and the conclusion to final-state phrasing. The physics text cites
  only live artifacts; the single retirement fact stands in the
  Appendix B availability note.
- reference.md: the conventions bullet re-points to Sections 10.5/12.1
  (the dead Section 24 reference removed); the V52 row and the
  classification sentence already carried the 17-row 8/0/2/2/5 ledger.
- Ledger cutover verification: V52 standalone re-run PASS, provenance
  table 17 rows, counts 8 algebraic / 0 free / 2 convention / 2
  coincidence / 5 named residual, status string reads the 8+0+2+2+5
  split (the old 8+0+1+2+5 string is gone).
- Suite: the re-gate booked in the note above remains valid (no
  experiment code changed in this step); results/validation_summary.json
  was re-synced to the on-disk artifacts with run_validation.py
  --report-only: 200/200 PASS, backend-wide pass True, hard tier True,
  n = 200, the V52 summary entry carrying the new split.
- review/tools/number_audit.py: PASS (799 of 823 decimal tokens
  artefact-backed; the one unexplained number, a 2.2e-33 residual typed
  into 10.5 from memory, was corrected to the booked V274 value
  4.8 x 10^-33).
- PDF rebuilt: doc/final_publication.pdf, 86 pages.

## Retirement note: the bond-current alternative chains and the superseded composition audit (2026-10-28)

Ten experiments are retired to git history: the `alt_gauged` chain
(V175-V179), the `alt_nonpert` chain (V180-V183), and V94. No new
experiment was pre-registered by this step.

- V175-V179 and V180-V183 booked the same route twice: the gauged
  bond-current interaction as the one internal mechanism that could
  isotropize the twelve M-node cones. The chain ended in two gates with
  zero candidates (V179: minimum anisotropy 0.9322 at one loop and under
  coarse-graining; V183: exact dressing shifts the anisotropy by
  4.4e-12 and the flux sector is massive at every volume probed, gap
  42.09 -> 12787.60 -> 2.75e8 eta). The publication's isotropy verdicts
  are carried by the surviving cone programmes (V115, V118-V121, V142,
  V185-V187, V192, V201, V213, V219, V227, V235-V236, V244-V245) and the
  dynamical-link continuum programme (V251-V254), which re-measured the
  same one-loop and block-RG questions on the constructed sector: V253
  re-derives the downfold ladder independently (its bare rows reproduce
  V178's booked ladder to 1.8e-12, booked in V253's own artifact), and
  V254 carries the F9 contact to b = 16. The S-matrix and vertex facts
  the chain added are not cited in the paper text.
- V94 measured the composition law (multiplicativity across legs, not
  along a chain) before the clock arbitration; the paper books that
  verdict from V210, which measured it with the step and the traced-link
  phase closed together at the open clock (Sections 16, 18, Table 9).
  V94 is the same law's older version.
- Dependency re-audit: two live consumers did read retired artifacts -
  V241's K4 anchor (the published V176 Casimir defect) and V253's T0
  anchor (V178's booked A_bare/A_gauged rows and A_full). Both were
  pure cross-checks, not computations: the scripts now carry the exact
  booked values as literals labeled as coming from the retired
  artifacts in git history, and both re-ran green (anchor deviations
  unchanged). No other surviving script reads a retired artifact (the
  remaining cross-readers sat inside the retired set); the
  cl10/altcell.py library stays - it has 29 live consumers. The two
  exploration notes (doc/alt_gauged_exploration.md,
  doc/alt_nonperturbative_exploration.md) are retired with their
  programmes.
- Table 9 rows rewired: the gauged-bond-layer row now names
  V79, V97-V101, V129, V130, V137; the closed-loop row names V210
  instead of V94. Suite counts 200 -> 190 (paper Appendix B and the
  availability note, README, experiments/README, reference.md); the full
  190-experiment gate re-runs on the current tree.

# The light node: pre-registrations V277-V280

Category `dynamics_spectrum` for V277-V279; V280 is the F9 contact of
whatever candidate survives. The programme runs on the CPU core
(V185-V187 convention: `CL10_BACKEND=numpy` pinned before any `cl10`
import; the Nelder-Mead search over 128x128 Bloch matrices measures
~40x faster on CPU workers than on the card).

Reader's note on a pre-existing text-artifact gap (measured before any
V277 code was written, so no run can be blamed for it): the
publication's V250 decoration row prints "cone degree 0.281 -> 0.540"
(Sections 8.4, 11.6, Table 9). The live V250 artefact carries no such
pair: its `cone_degree_gamma` reads implemented-clock bare
-0.10743390423138056 and decorated -1.2196990555403142, closed-clock
bare 0.5086736873158906 and decorated NaN (a gap-floor artefact: at
the closed clock the star step's minimum eigenphase is pinned to the
1e-17 noise floor at every radius, so that fit measures noise). No
live artefact in the suite carries 0.281 or 0.540 under any cone or
degree key. The anchors below therefore bind to the live artefact
values, not to the printed pair; the printed pair is booked as a
text-only residue of the pre-retirement run.

## V277_bond_decoration_cone_scan (tier 2, hard, O10)

**The question (pre-registered).** Does the enlarged bond-decoration
circuit class - the six-stage step with a per-stage spin-connection
quaternion twist on the internal factor - carry a protected conical
node at GAMMA: gapless, linear (cone degree ~1), subluminal
(slope <= 1 in block units) and non-Abelian (rho_comm > 0.5)? V250
booked the single common-amplitude decoration on the star module;
this programme generalizes the amplitude per staged bond and reruns
the cone question on the 16-sector supercell step, the engine's own
circuit.

**The object (closed input set).** `cl10.supercell` published
machinery: the six-stage step M(q) = prod_t kron(S_t(q; beta), U16)
at the closed clock (beta, g) = (pi/4, pi/4) (V239 closure set, the
order-2 point); the 16-sector internal U16 and basis B16 from
`supercell.sectors(pi/4)`; the quaternion triplet
(V1,V2,V3) = (i h_1 h_2, i h_2 h_3, i h_1 h_3)|_16 from
`cl10.clifford.generators` conjugated by B16 (the V250 construction,
Hermiticity < 1e-12 measured in-script). Decoration: per stage
s = 1..6 the internal factor becomes U16 W_s(q),
W_s(q) = exp(i q.(alpha_s V1 + beta_s V2 + gamma_s V3)). 18
parameters in the box [0, 2 g10] per component with g10 = 0.2 the
implemented config value (`cl10.config`, angle.branches.implemented).
Unitarity is exact (Hermitian generators exponentials, tensor blocks
of unitary stage matrices); causal radius 1 holds (the twist is a
q-local unitary on the internal factor, the bond support is
unchanged). W_s(GAMMA) = I for every parameter set, so the decorated
family shares M(GAMMA) with the published step exactly - the node
question is purely about the dispersion of its gapless multiplet.
Optimizer: `scipy.optimize.minimize` Nelder-Mead (the V186 tune
polish), maxfev 300, candidates projected back into the box at every
objective call; 120 starts (seed 277: start 1 = the published
uniform triplet (alpha,beta,gamma) = (0.2,0.2,0.2) on all six
stages, 119 = uniform draws from the box); workers via
`cl10.parpool`. Objective: the coarse cone degree
polyfit(log r, log w(r)) with w(r) = min over the 13 cubic classes
{(100) x3, (110) x6, (111) x4} of the top of the gapless multiplet
(the 64th sorted |eigenphase|, m0 measured at T0) at radii
(0.1, 0.03, 0.01); maximize.

**Measured (pre-registered).** T0 machinery gate against the live
V250 artefact (the V250 star-cell construction re-derived in-script:
`starcell.bloch_star` on the 1024-dim module cell with the published
star axes and `module_stages(order="phase-major")`, and its single-
twist decorated step `U16 W(q)` with amplitude g10 = 0.2, using the
V250 degree convention: radii (0.3, 0.1, 0.03, 0.01, 0.003), the
star's 6 axes plus the 3 Cartesian axes): implemented-clock bare
degree -0.1074339042 and decorated -1.2196990555 (tol 1e-3);
implemented decorated cone_anisotropy A = 0.0109382891 (tol 1e-5);
closed-clock bare cone_anisotropy A = 0.999999999945298 and
decorated A = 1.0 with (vmax, vmin) reproducing the artefact to
1e-9 relative; closed-clock degrees booked as measured (bare
~0.509 tol 1e-2; decorated must reproduce the artefact's NaN).
Then, on the 128-dim supercell step at the closed clock: the GAMMA
eigenphase census (zero-count and pi-count, expected 64 + 64 from
the V192/V132 pinned spectrum, measured not assumed), m0 = the
measured zero-phase multiplicity; the bare top-m0 degree (expected
~1.9986, the quadratic touching V245-type - booked with tol 2e-3 as
the starting value); the published uniform-triplet start degree
(booked, no external anchor); the search's best point with its fine
measurements: degree over radii (0.3, 0.1, 0.03, 0.01, 0.003) and
the 13 cubic directions; the GAMMA gap; the tangent maps
S_j = i P^dM (d_j M) P on the m0-gap subspace by the V186
`tangent` convention (h = 1e-5, central differences),
rho_comm = max ||[S_a,S_b]|| / (1/3 sum_j ||S_j||^2) over cyclic
pairs, and the tangent speeds over the 179-direction Fibonacci
sphere (the V235 `fib_dirs` pattern): v_max (the slope leg, block
units) and A_top = (v_max - v_min)/(v_max + v_min) over nonzero
speeds. All booked for the best-degree point and the best
(degree, rho_comm, v_max) triples.

**Acceptance (pre-registered).** T0 fail (any anchor beyond its
tolerance) -> MACHINERY-DRIFT, exit nonzero, nothing downstream
interpreted. DECORATION-CONE-FOUND iff some searched point has:
GAMMA gap < 1e-10 (automatic in this family since W(0) = I - booked
as such, not counted as earned), fine degree >= 0.99, v_max <= 1,
and rho_comm > 0.5. Else DECORATION-BEST with the best triple
(degree, rho_comm, v_max) booked and the binding leg named (which
criterion the best point misses).

## V278_range2_parity_cone (tier 2, hard, O10)

**The question (pre-registered).** Does the circuit class enlarged by
range-2 second-neighbour swap rotations reach the parity-mixed
A1g (+) T1u conical node of V191 inside the engine's own six-stage
product - the node Proposition 6 puts outside every single-parity
cubic-covariant class and V192/V201 searched without finding?

**The object (closed input set).** The V277 family extended on the
bond factor. Two new orbits, both summed over the full O_h orbit so
each factor is a singlet set: the axial orbit generator
A(q) = sum_j sum_{e: bit_j(e)=0} [e^{+i q_j} |e (+) (1<<j)><e| + h.c.]
(the (plus/minus 2,0,0)-type second-neighbour swaps: the same pair
set as the engine's inter-block stages, entered with an independent
amplitude per stage), and the diagonal orbit generator
D = (1/2) sum_{i<j} sum_e [|e (+) (1<<i) (+) (1<<j)><e| + h.c.]
(the (plus/minus 1, plus/minus 1, 0)-type intra-block double flips,
12 pairs, q-independent). One amplitude each per stage, box [0, 2 g10]
like the decoration: 30 parameters total. The stage factor is
kron(S_t(pi/4) exp(i a_t A(q)) exp(i d_t D), U16 W_t(q)); causal
radius 2 (the axial factor carries e^{+/- i q_j} exactly like the
engine's own inter-stage swaps: reach +2 per stage, the published
Lieb-Robinson radius per stage is 2 in module units - measured in
T0-style check), unitarity exact per factor. O_h content of the
GAMMA gapless subspace is read by the V191 machinery reused
verbatim (`class_of`, `CHAR`, `projector`): multiplicities from
characters over the 48 W48 (x) I16 on the gap subspace, and the
A1g, T1u projector ranks on it (the equivariance-free content
statement; V192's measured exact Gamma equivariance of the published
product is re-measured in-script and booked for the extended family
too - if the new factors break it, the content reading is still a
rank measurement, stated as such). Optimizer as V277 (Nelder-Mead,
150 starts, seed 278, maxfev 300 each, box projection) with
objective = the GAMMA gap (minimize, top-multiplet quality booked
not searched).

**Measured (pre-registered).** the gapless-subspace dimension m at
the best point; its O_h content multiplicities and A1g/T1u projector
ranks; the GAMMA gap; the V186 tangent maps on the gap subspace:
rho_comm and the Fibonacci-179 speeds (v_max, A_top); degree of the
gap multiplet as in V277; Gamma equivariance defect of the extended
product over the 48 elements; the best point's 30 parameters.

**Acceptance (pre-registered).** CONE-RANGE2-FOUND iff at the best
point: GAMMA gap < 1e-10, the gap-subspace content contains an
A1g (x) T1u parity-mixed pair (A1g rank >= 1 and T1u rank >= 3 with
total gap dimension 4), A_top < 0.1, rho_comm > 0.5, v_max <= 1.
Else CONE-RANGE2-ABSENT, booking the best point's (gap, A_top,
rho_comm) triple and which leg binds.

## V279_candidate_protection (tier 2, hard, O10)

**The question (pre-registered).** Is the best node produced by
V277/V278 protected: do on-site perturbations that respect the
commutant leave the node gapless? If neither V277 nor V278 produced
a candidate (the pre-registered default expectation), the same test
runs as a CONTROL on the star-module two-fold line (V219/V220/V227
machinery) and its verdict is booked as control-only.

**The object (closed input set).** The candidate step (best point of
V277 or V278, carried by value); the V220 `commutant` Cesaro average
over 240 powers of the step at Gamma to find the invariant on-site
operator algebra, `cluster_frame` and `circle_holonomy` as booked by
V220; the V227 invariant-mass opening test verbatim (`refine`,
`exponent_away`): V = projector of a random Hermitian (seed 279)
onto the invariant commutant, perturbation exp(-i theta V), gap
versus theta slope at theta = 1e-2.

**Measured (pre-registered).** the commutant dimension; the induced
gap at theta = 1e-2 minus the node gap (opening); the gap-vs-theta
slope and the obstruction norm (invariant projector norm); the
holonomy charge of the node (V220 `circle_holonomy`) as booked data.

**Acceptance (pre-registered).** CANDIDATE-PROTECTED iff the induced
gap at theta = 1e-2 < node gap + 1e-6; else CANDIDATE-MASSIVE with
the opening slope booked. Control runs carry the same verdict string
with a "control:" prefix and never gate anything.

**Gate G1 (fixed before any run).** If V277 and V278 are both
negative (DECORATION-BEST / CONE-RANGE2-ABSENT) OR V279 says
CANDIDATE-MASSIVE on a real candidate: the enlarged (decorated,
range-2) circuit class is booked in the ledger as yielding no
protected conical node, O10 stays open with that measured extension,
V280 still runs if a candidate exists for the F9 contact read, and
V290 (Phase 4) is SKIPPED.

## V280_candidate_sme_contact (tier 4, F9) - runs only if a candidate exists

**The question (pre-registered).** What SME electron velocity
anisotropy would the surviving candidate carry at the Planck-lattice
scale: what is its delta v / v against the booked experimental
ceiling?

**The object (closed input set).** The candidate node's dressed step
(by value from V277/V278/V279); the V235 `node_velocity_contact`
machinery (`fib_dirs`, `step_batch_at`, `velocities`, `tensor_fit`,
`dvv`); the ceiling `sme_electron_velocity_aniso = 8.0e-15` from
`cl10/config.yaml` bounds. No new inputs.

**Measured (pre-registered).** the candidate's velocity tensor fit,
delta v / v (the V235 `dvv` convention), the ratio to the ceiling,
the 179-direction spread.

**Acceptance (pre-registered).** pass = the contact is measured (the
V235 convention: contact experiments pass by measuring); status
F9-CANDIDATE-CONTACT with the ratio booked.

## Roll-up consequences of this programme (fixed before any run)

- Every phase appends its own dated reference.md paragraph with the
  measured numbers and verdict strings, negatives included, and
  re-gates `review/tools/number_audit.py` to PASS.
- Gate G1 above decides Phase 4 (V290); no post-hoc re-run of the
  search on a different objective is allowed after a negative.
- The suite recount propagates (README.md, experiments/README.md,
  paper Appendix B, section 24) whether or not any route succeeds.
- No acceptance threshold moves after a run; the published step's
  Gamma gapless multiplet (64+64) is booked as measured at T0 and
  the "gapless" leg of the acceptance is reported as inherited, not
  earned.
- Nothing in this programme is refitted to land in a cone band; any
  dial touched reclassifies to fitted in the ledger.

# Loop origin or booked decision: pre-registrations V281-V284

Categories: V281 `scales_loopfactor`, V282 `gauged_layer`, V283
`gauged_layer`, V284 `scales_loopfactor`. The programme attacks O8-new
(can the loop base 1/(16 pi) be produced by a traced sector of the
root?) with three independent traced constructions and one final
arithmetic word. Expected outcome per the V252/V274/V275/V276
booked record is negative; every branch is booked.

## V281_flux_kernel_trace_loop (tier 2, hard, O8-new)

**The question (pre-registered).** Does the traced-word Cartan-ratio
functional (V209/V203 machinery) produce the loop band
1/(16 pi) = 1.9894e-2 when the traced link is not the qubit link but
the flux kernel of the three-site gauged chain, and the traced link
state is the bond-current ground state instead of the uniform link
average?

**The object (closed input set).** V101's cached kernel (Is, Ts from
results/kernel.npz, rebuilt with --force if absent) + V79
`root_generators`/`structure_constants`/`adjoint_link_operators` +
V101's exact (H, N) generalized eigenproblem on the 24 kernel states
(anchors below); the link-pair space is the two-link pair (45 x 45 =
2025) of the V99 three-site chain with the two-link Casimir C2ll
(V101 construction region) and channels grouped by the V129
`channels` convention (values expected {0, 64, 96} with the 210
channel at 96 as booked by V101's extraction). The traced state is
the link-pair reduction rho_lp of the degeneracy-6 ground-state
mixture Cw = c c^dag / 6, built by explicit tensor contraction of the
glued kernel basis with V101's own Gram einsum algebra; the
contraction is gated: Tr rho_lp = 1 (1e-9), positivity (min eigenvalue
> -1e-14), and V101's booked middle-cut entropy S = 3.80666248977
reproduced from rho_B reconstructed by the same contraction (1e-9).
The functional is V209's: the sector phases e^{i(g/2)e} over the
qubit-link Casimir sectors {-32, -12, +4} are REPLACED by the
channel-phase operators M = sum_c sqrt(w_c) exp(i theta_c) Pi_c over
the C2ll channels of the link pair (Pi_c the exact Lagrange
projectors of the integer C2ll, w_c = Tr(rho_lp Pi_c) the ground-state
channel weights), with the spinor-block machinery G, N, planes, cart
imported from V209 (`skeleton`, `tr_AB`, `cint`) unchanged; the
Cartan ratio F_n = cnorm((M G)^n)/cnorm(G^n) on the 32-dim spinor
block. Grid: the channel phases (theta_0, theta_45, theta_210) on a
13^3 = 2197-point product grid (linspace 0..2pi), float64 F_n for
even n <= 40; stability control: the same grid at 14 points (F_40
shift < 1e-8 required, else the point is unstable and not booked);
certified continuation: V209 `mp_ladder`-style 50-digit mpmath
re-run at the closest grid point (exact skeleton, `mp.mpf` chain)
along the measured geometric class.

**Measured (pre-registered).** T0: kernel 24 states; E0 =
-5802.398188196636 eta, degeneracy 6, gap 42.09063148320456 eta, S_ground
= 3.8066624897703183 (V101 artifact `data.*`; tol 1e-9 relative);
link-pair channels [value, multiplicity] by V129.channels; the
ground channel weights w_c; the F-grid: min over stable grid points
of |16 pi F_4k - 1| and the F-band census; the 50-digit certified
closest approach and the geometric ratio of the measured class;
closest approach |16 pi F - 1| booked at every leg.

**Acceptance (pre-registered).** T0 fail -> MACHINERY-DRIFT, exit
nonzero. FLUXLOOP-FOUND iff some stabilized grid point has
|16 pi F - 1| < 0.5 AND the F value lies in the pre-registered band
[1.5e-2, 3.0e-2]; else FLUXLOOP-NEGATIVE with the floor (closest
|16 pi F - 1| and its F) booked. No dial of the grid is refit after
the run.

## V282_link_cubic_invariant (tier 2, hard, O8-new)

**The question (pre-registered).** Is the frozen-link verdict of
V271/V274 (kernel ladder dispersion exactly zero) an artifact of the
quadratic-only link Hamiltonian? Add the cubic link invariant C3 (if
one exists on the 256-dim link module) and re-measure: dispersion
first, then - only if the link moves - the V274/V275 exchange
re-run.

**The object (closed input set).** V269's dynamical link module
construction region (the 256-dim module, the 45 generators restricted
to it, the composite kernel), the V270/V271 artifacts (anchor below),
V271's kernel ladder (channel spectrum + `krylov_steps` frozen
test) and V274's second-order exchange construction region (the
block-Krylov resolvent) with V275's off-kernel growth machinery.

**Measured (pre-registered).** T0: re-derive from the live V271
artifact - gauss defect 7.698474841966662e-14, kernel channel ladder
and its dispersion (the booked frozen ladder: 0 / 19.2 / 28.8 with
defect 1.8e-14 - exact keys read from the artifact at run time),
link-module kernel dim 6, and the V269 link-closure fact (module
1 (+) 45 (+) 210 closes inside the root - re-measured by V269's own
test region). T0.5: dim(Sym^3(link module)^G) by Monte-Carlo Haar
averaging over the connected group generated by the 45 generators
restricted to the 256-dim module: random group elements
exp(i sum c_A L_A) with c Gaussian, 300 draws (the V99 Haar
pattern, parpool), projector-average trace on the symmetric cube
(tol: the averaged character stabilizes: run-to-run spread across
three 100-draw batches < 0.02; dim = round of the mean). If dim = 0:
status CUBIC-ABSENT, T1/T2 skipped - the "add a cubic term" escape
is closed with a computed number. T1 (dim > 0): H_link -> H_link +
eta C3 (C3 the projected cubic Casimir on the module: the symmetrized
Tr of three generators in the module representation), V271's kernel
ladder re-measured: dispersion (max - min splitting over the kernel
channel spectrum). T2 (dispersion > 1e-8): the V274 second-order
exchange re-run verbatim with the modified link dynamics, and V275's
off-kernel resolvent reading F5 with the same mpmath control.

**Acceptance (pre-registered).** T0 fail -> MACHINERY-DRIFT.
CUBIC-ABSENT (dim Sym^3 invariants = 0, booked with the Monte-Carlo
spread); CUBIC-PRESENT-FROZEN (dim > 0, dispersion <= 1e-8: the link
stays frozen under the cubic term); LINK-DISPERSIVE + the resolvent
F reading and |16 pi F - 1| as measured (the V273/V275 band
[0.015, 0.03] booked against, never searched).

## V283_metric_trace_loop (tier 2, hard, O8-new)

**The question (pre-registered).** Does the metric sector - the
directional quintet modes - traced in its ground state produce the
one-loop Cartan ratio near 1/(16 pi) in the closed two-hop traced
functional?

**The object (closed input set).** V252's two-particle F_int
machinery (the construction region: the dressed two-hop and its
Cartan ratio; anchor: V252's booked floor for |16 pi F_int - 1| read
from its artifact at run time, tol 1e-6 relative), the V140
directional-link quintet machinery (`build_group` line ~98,
`structure_constants` line ~107, the M_dir/bands pattern lines
~287-293) and the V229 spin-2 insertion (`spin2_generators` line
~154, the 110-dim sym2 basis), coupling lambda = g10 at the
implemented branch value from config. Construction: the quintet
modes enter the traced two-hop as an extra traced degree of freedom:
the closed loop functional of V209's Cartan-ratio type evaluated on
the two-hop dressed by the quintet's ground-state insertion at
lambda = g10; n <= 30, float64 ladder, stability by step-doubling
control (grid convergence row), closest |16 pi F - 1|.

**Measured (pre-registered).** T0: V252's booked F_int floor
re-derived by V252's own machinery (imported), the quintet spectrum
of V140's M_dir at the committed clock (band multiplicities booked),
the spin-2 insertion norm; the F_n ladder n <= 30 with the insertion
at lambda = g10 and the step-doubling stability rows; the closest
approach |16 pi F - 1| and its (n, direction) label.

**Acceptance (pre-registered).** T0 fail -> MACHINERY-DRIFT.
METRICLOOP-FOUND iff some stable n <= 30 has |16 pi F - 1| < 0.5;
else METRICLOOP-NEGATIVE with the floor booked.

## V284_crossing_arithmetic_final (tier 4, O8)

**The question (pre-registered).** Is the certified crossing n* = 208
and the geometric ratio rho = 0.9780302704295731... arithmetically
related to 16 pi at high precision, closing or finally booking the
crossing arithmetic?

**The object (closed input set).** V209's exact `skeleton` +
`mp_ladder` at 100 and 200 digits on the certified closed branch
(the 4k class at the root clock, defined by continuity as V209
books); the V203/V204 float64 machinery (as V276 used it); the V204
artifact (rho, crossing_n 208, F_cross) and V276 artifact
(rho_abs_diff record). mpmath PSLQ over bases
{1, rho, rho^2, ..., rho^8, ln rho, ln(16 pi)} with integer
coefficients, degree <= 8, height cap 1e6, tolerance 1e-14; the
rationality of ln(16 pi)/ln(1/rho) to 1e-10 by continued-fraction
denominator growth; |n*_geom - 208| with the 200-digit rho band.

**Measured (pre-registered).** T0: rho reproduced to 6.7e-16 and
crossing_n = 208 from the live float64 machinery (V276's own test);
the 100-vs-200-digit 4k ladder ratios (stability row); the PSLQ
output (shortest relation and its residual); the continued-fraction
read of ln(16 pi)/ln(1/rho); n*_geom from the certified ladder and
its distance to 208.

**Acceptance (pre-registered).** CROSSING-ALGEBRAIC iff an integer
relation in the booked basis verifies to residual < 1e-10 at 200
digits - then the relation is reported; else CROSSING-MEASURED-FINAL
(O8's crossing arithmetic booked as measured, final).

## Roll-up consequences of this programme (fixed before any run)

- Gate G2 / decision D1 (pre-decided, no consultation): if V281,
  V282-T2 and V283 are all negative, the loop base stays a declared
  input and the paper books the two-sector clock structure with the
  three exact edits specified in the programme plan (Section 9.2
  append; Table 6 convention row names V281-V284; Section 18
  loop-factor part adds "(V281-V284)"). F2 stays FIRED in Table 8.
  The same wording is booked in the ledger. If any leg is positive
  the paper edits are NOT made and the positive construction is
  booked instead, with its own ledger paragraph.
- The suite recount propagates regardless.
- No threshold moves after a run; every negative is booked with its
  measured floor.

# The graviton cone and content: pre-registrations V285-V289

Categories: V285 `gravity`, V286 `gauged_layer`, V287 `gauged_layer`,
V288 `gauged_layer`, V289 `gravity`. Attacks F8 (no c = 1 branch in
the closure class) and F4 (quintet bubble) from the construction
side, plus the Bekenstein-map bookkeeping and the geometry measure
proof of principle.

## V285_most_local_lightlike_generator (tier 2, hard, F8/O9)

**The question (pre-registered).** Within the 53-dimensional family
of exact-lightcone generators (V208: the Clifford-regularity
Jacobian at the committed generator has rank 107 in the 160-dim
covariant velocity space, family dimension 53), is the committed
generator selected by locality - or is it a boundary accident of the
family?

**The object (closed input set).** V208's construction region
imported by path (the covariant 160-space basis, the committed
generator, the Clifford-regularity Jacobian as V208 builds it -
lines ~428-466); the nullspace basis of the Jacobian restricted to
the covariant space (family dim re-measured, expected 53); V208's
kernel measurement as the objective: tail functional T(A) =
sum_{r >= 2} |kernel_r(A)|^2 on the hopping-kernel profile
(committed values, V208 artifact `cone.kernel_profile`:
0.9108612370986987 at r = 0, 0.10705424159710794 at r = 1, then
0.034138977783196765, 0.005435677538772973, 0.0007878570698366837,
7.157945101618305e-05, ...). Minimize T over the 53 family
coordinates by Nelder-Mead, 50 starts (seed 285, maxfev 200 each,
quadratic family parameterization through the exponential/log map of
the generator space so exact-cone membership is maintained by
construction - the tangent coordinates are the nullspace
coefficients).

**Measured (pre-registered).** T0: the committed generator's kernel
profile and tail, the hop projection 0.0028977782769928693 and
dressing projection 0.0092074659831694 (V208 artifact
`implemented.*`; tol 1e-9 relative); Jacobian rank 107 and family
dim 53 re-derived. Then: the minimizer's tail T*, its hop projection
and dressing fraction, distance of the minimizer from the committed
generator (relative), and whether the minimum lies in the interior
(gradients < 1e-6 at the minimizer, 50/50 starts converge into one
component) or on the family boundary.

**Acceptance (pre-registered).** T0 fail -> MACHINERY-DRIFT.
CONE-LOCALLY-SELECTED iff an interior strict minimum with T* <
T_committed and hop projection > 0.0092 (the selected member is
booked as a selection FACT, not yet a physics principle); else
CONE-SELECTION-ABSENT with the boundary data booked.

## V286_module_link_invariant_hop (tier 2, hard, F8)

**The question (pre-registered).** On the 256-dim dynamical link
module (V269), does the commutant of the combined vertex charge over
(head, matter, tail) admit an invariant hopping beyond the known
regular-link descending-to-identity map, and if a hop is built, does
any branch of its band spectrum reach |v| = 1 at a closed clock?

**The object (closed input set).** V269's link module (module 1 (+)
45 (+) 210, 256-dim, the 45 generators restricted to it), V272's
exchange construction region (the head/tail maps), the V220
`commutant` Cesaro pattern (240 powers) applied to the combined
vertex charge on the tensor product; Gauss-invariance test of any
candidate insertion (leak = max norm of the vertex-charge
commutator, required <= 1e-10); band spectrum by cupy `eigvalsh`
batches on the 128-dim Bloch stage matrices (the V252 `_hs_chunk`
pattern for batched Hermitian slices).

**Measured (pre-registered).** the commutant dimension and Cesaro
obstruction; the basis matrices (gauge-classified against the known
descending map: relative distance of each basis element to its
projection on the known map); if a hop exists: its Gauss leak, the
band maxima |v| over the closed-clock set (beta, g) from V239's C,
and whether max |v| = 1 (tol 1e-9) for any branch.

**Acceptance (pre-registered).** LINKINV-OBSTRUCTED(dim) if the
commutant is spanned by the known map (projection residual < 1e-8 on
every basis element); LINKINV-HOP with the band table and the c = 1
verdict otherwise. Both book the spectrum.

## V287_quintet_bubble_null (tier 2, hard, F4)

**The question (pre-registered).** Does any O_h-covariant stress-
vertex family null the quintet polarization bubble at lambda = g10 -
closing F4's "the bubble is nonzero on the full quintet" at the
level of the coupling family?

**The object (closed input set).** V229's machinery imported
(`sym2_basis` 110, polarization tensors Emat, `spin2_generators`
line 154, the `bubble(Tb)` construction lines ~393-410), V141's
directional step for the star-bond set; coupling family:
T_ij = sum_d w_s f_s(i, j, d) J_d over direction sets - axes, star
axes, and general range-2 O_h-symmetric direction sets, at most 8
weights (parameters) per family point; Nelder-Mead 80 starts (seed
287, maxfev 200) minimizing the Frobenius norm of the bubble
operator Pi on the quintet at lambda = g10.

**Measured (pre-registered).** T0: the booked axis-bond bubble
eigenvalues (-54.24, -53.40 to booked precision from the V229
artifact `star_bond_bubble_eigenvalues` / `axis_bond_bubble_
eigenvalues`, tol 1e-6 relative) re-derived; the minimizing weights,
min ||Pi||, and the F4 mass shift at the minimizer (V229 convention).

**Acceptance (pre-registered).** BUBBLE-NULL iff min ||Pi|| < 1e-8;
else BUBBLE-FLOOR with the measured value.

## V288_bekenstein_truncation (tier 2, hard, unification-with-gravity bookkeeping)

**The question (pre-registered).** Which flux truncation of the link
module carries the Bekenstein factor toward 1/4, or none - the
V99/V100 ln 45 vs 1/4 gap named as a number: what does the
truncation map of the three-site kernel do to it?

**The object (closed input set).** V101 kernel builder + ground
state; the V99 Haar vertex-invariant region; V129/V130 entropy
machinery (`entropy_from_rho`-style eigenvalue entropy, `channels`);
truncations: rebuild the three-site kernel with the link flux space
restricted to channel subsets {1}, {1,45}, {1,45,210}, {45}, {210}
(each subset re-runs the zero-mode/Gauss construction on the reduced
channel space); per truncation: ground-state middle-cut entropy per
bond (V101 pattern) and the vertex-corrected coefficient.

**Measured (pre-registered).** T0: ln 45 = 3.8066624897703183 and
the V101 ground entropy reproduced (1e-9); the booked capacity
S_cap = 4.4998 and factor 6.6 read live from the V100/V129 artifacts
at run time; per truncation row: (dim kernel, E0, degeneracy, gap,
S_bond, factor to 1/4).

**Acceptance (pre-registered).** BEKENSTEIN-MAP booked (the map
rows); BEKENSTEIN-CLOSED only if some truncation lands within 2 %
of 1/4 (the factor in [0.98, 1.02] x 1/4 x the booked units). No
truncation is chosen in advance; every row is booked.

## V289_geometry_mc (tier 3, hard-lite)

**The question (pre-registered).** Does the root's finite geometry
group (the order-1024 even Clifford monomial group, V140
`build_group`) support a measure: Metropolis on the Wilson plaquette
action in the adjoint 45 representation on a 4^3 periodic lattice,
with the Wilson loop estimator finite and its decay measured? This
is a measure-existence proof of principle, not a physics claim.

**The object (closed input set).** A cl10-free lattice (numpy/cupy
only) of 4^3 periodic cells; bond variables u_b in Gamma from
V140 `build_group` (mul/invof lines ~75/90); action = sum over
plaquettes of (1 - Re chi_45(U_p)/45) with chi_45 the adjoint
character from V140 `structure_constants` exponentiated to 45x45
group elements; beta sweep 5 fixed beta values; 32 replica batches
on cupy, 20000 sweeps, seed 289; Wilson loops W(R), R = 1..3.

**Measured (pre-registered).** acceptance rate per beta; <W(R)> and
its error (jackknife over replicas); the fitted area-law exponent
from log <W(R)> vs R^2; the adjoint character check: chi_45 of a
group element vs the class table (exact match on all sampled
plaquettes). Pass = the estimator is finite and the decay measured;
MC-GEOMETRY-OK when the area fit exists; MC-TUNED books the beta
range that keeps the estimator in the tunneled-free regime.

## Roll-up consequences of this programme (fixed before any run)

- F8/F4 wording is NOT changed in the paper by this phase; the
  statuses feed the ledger. The single exception: if V285 and V286
  JOINTLY produce a c = 1 branch (CONE-LOCALLY-SELECTED with hop
  projection > 0.0092 AND LINKINV-HOP with a band at |v| = 1),
  Section 11.6 and the Table 8 F8 row gain the pre-written wording:
  "the selected member of the V285 family with the V286 invariant
  hopping carries a branch at v = 1 within the closure class (V285,
  V286); F8 stands as the standing falsifier".
- The suite recount propagates regardless.

### V281 construction addendum (fixed before any run; supersedes the
object wording in the block above where they differ)

The contraction design was fully digested after the block above was
written and before any code ran; the addendum fixes the closed object
honestly:

- The bond-current ground state enters through its END-PAIR
  (middle-cut) reduction rho_B, which V101 explicitly builds (Tr = 1,
  S = 3.80666248977 booked). The link-Pair reduction of the earlier
  wording is not a state V101's algebra closes on (the shared-site
  web contraction is not positive in general); the state used is
  rho_B, gated exactly as booked: Tr 1 (1e-9), min eigenvalue >
  -1e-14, entropy S reproducing ln 45 to 1e-9.
- The traced sector is the END-PAIR space (1024 = 32 x 32) of the
  chain with the two-ended closed internal map G_pair = G (x) G
  (G from V209's planes at the certified closed angle pi/4). The
  channel projectors Pi_e are the exact Lagrange projectors of the
  INTEGER end-pair Casimir C2pair = sum_A (kron(T_A,I32)+
  kron(I32,T_A))^2 (V101's operator; channel values measured by the
  V129.channels convention - expected {0, 64, 96}, the 210 channel
  of the earlier wording lives on the link pair and is booked in the
  channels row separately, not used in the ladder).
- The ladder: F_n = cnorm((M G_pair)^n) / cnorm(G_pair^n), M =
  sum_e sqrt(w_e) exp(i theta_e) Pi_e, w_e = Tr(rho_B Pi_e) the
  ground-state channel weights; the content is rho-weighted:
  cnorm(Y) = sum over the 10 Cartan planes (the five V209 Cartan
  planes carried on each tensor factor) of |Tr(H (Y rho_B))| -
  with rho_B = I/1024 uniform this reduces to the V209/V203
  structure, so the replacement of the traced state is the only new
  input.
- The 13^3 phase grid, the n <= 40 even ladder, the 14-point
  stability control and the acceptance thresholds are unchanged from
  the block above. The certification leg: the mpmath 50-digit
  matrix ladder at 1024 dimensions is not reachable in the
  interpreter; the booked certification is the double-precision
  stability ladder (float64 vs longdouble agreement, step-doubling
  grid control) plus the 50-digit mpmath SCALAR continuation of the
  measured geometric class (F_{4k} -> rho_chain fit, 50-digit PSLQ
  check of the class), and the limitation is booked in the verdict
  text.

### V282 T0.5 measurement addendum (fixed before any booked run)

The first execution attempt of the V282 T0.5 Monte-Carlo as literally
worded above (300 draws of exp(i sum c_A L_A), c Gaussian) was run as
a dry diagnostic and its batch spread measured 652 against the 0.02
tolerance: the per-sample variance of chi(g)^3 under that pushforward
is ~1e6, so 300 draws cannot concentrate for ANY group - the wording
is numerically unsound, not the question. The measurement is
therefore executed as the exact Haar average reduced to the maximal
torus by the Weyl integration formula: the Cartan weights of the
module (and the root system) are measured from the 256-dim module
generators themselves (no typed group data), the class-function
integrand chi_sym3(t) = (chi(t)^3 + 3 chi(t) chi(t^2) + chi(t^3))/6
is smooth on the 5-torus, 3 x 10^4 flat samples with the Weyl
density measured from the roots; the three-batch spread tolerance
0.02 is kept verbatim. The genuine-hopping cubic is additionally
tested EXACTLY (not by MC) on the adjoint sector: the nullspace of
the 45 intertwining constraints on Sym^3(45) (dim 15180), the sector
any link hopping term lives in; the singlet-mixed invariants of the
full module are booked separately and are not hopping terms. No
verdict was booked from the dry diagnostic; the thresholds (dim = 0
branch, dispersion 1e-8) are unchanged.

### V277 T0 closed-clock contract amendment (fixed after the machinery
probe, before any result is booked)

The first fixed-clock execution of V277 T0 re-derived the booked
implemented-clock pair exactly (-0.1074339042 / -1.2196990555,
1e-12 agreement) and then flagged the single pre-registered "closed
decorated must reproduce the artefact's NaN" contract. Measured cause
(both runs' gap arrays, booked): at the closed clock EVERY star gap
in the bare and the decorated degree fits sits at the float64
rounding floor (1e-18 to 4e-16 in V250's booked artefact gaps;
similarly here) - the clock pins the whole multiplet, and whether
polyfit's log hits log(0.0) or log(1e-17) is a single-bit rounding
coin. The NaN landed on the decorated row there and on the bare row
here; the physics statement (closed-clock pinning of both to the
rounding floor) is identical. The contract is therefore amended to
the rounding-robust form: the re-derived closed-clock bare AND
decorated gap arrays must each be <= 1e-15 at every radius (the
V250 artefact rows satisfy this: max 3.5e-16), with the NaN/finite
placement booked as measured rounding noise, not drift. The
implemented-clock degrees, the anisotropy anchors, and every
threshold downstream are unchanged.
---

# The flavour-cosmology completion: pre-registrations V290-V299
(tier mix as noted; programme phases 4-7; every number below is read
from a booked artifact or re-measured from the closed input set; no
dial may be touched - any refit reclassifies to fitted and voids the
verdict.)

## V290_anisotropy_relevance (tier 2, runs ONLY if a G1 candidate exists)
**Pre-gate (fixed before any run).** V279 booked the candidate census
for the decorated cone family; V280 was already skipped by the same
pregate when the census came in at zero. If the V279 artifact carries
zero candidates, V290 does not execute and the roll-up books
V290-SKIPPED-BY-PREGATE citing the V279 artifact verdict verbatim.
The anisotropy-relevance question (does the candidate cone's
anisotropy run irrelevant under the booked coarse-graining) stays
open and is stated as such in the ledger.

## V291_higgs_host_channels (tier 3, flavour bookkeeping)
**The question.** Does the dynamical chain host charge-neutral scalar
channels - chain-kernel states on the module-singlet atom channel -
that could carry the booked 10-Higgs Dirac term (V45 higgs_dirac),
or does the flavour construction remain purely closed-branch?
**The object.** The V286 machinery re-derived live (module selection
inside the booked V269 c2_spectrum census, the V270 seed-path atom
kernel, the exact chain kernel of Nmat), the V270 channel anchors
(kernel dim 6, channels {0, 64, 96} = the booked c2_spectrum values),
the V272 chain_kernel_dim 12 and chain_H_spectrum anchors, and the
V45 higgs_dirac anchor (rank, determinant).
**Measured.** per chain-kernel basis vector: <C2l> on the module
factor (channel attribution by nearest booked value), the count n0
of module-singlet-channel kernel states, their chain-H eigenvalues
against the booked spectrum; kernel dim and atom-kernel anchors
as T0 gates (MACHINERY-DRIFT on any mismatch).
**Branch.** HIGGS-HOST-PRESENT (n0 >= 1: the chain supplies neutral
scalar hosts - number and energies booked) / HIGGS-HOST-ABSENT
(n0 = 0: booked null; the flavour sector stays closed-branch).

## V292_vev_host_mixing (tier 3)
**The question.** With the booked democratic mass matrix (V45
mass_matrix: all-ones, rank 1, eigenvalues {3, 0, 0}, TBM mixing
sin^2 theta12 = 1/3 booked) and the chain-host insertion from V291
(booked channel energies; no free dial), does the mixing move at
all, and by how much against the PDG sin^2 theta12 = 0.307?
**The object.** the V45 mass_matrix and mixing anchors re-derived
via cl10.yukawa_3x3, the V291 artifact (n0 and energies), nothing
else.
**Measured.** T0: V45 mixing re-derived (TBM to 1e-9). T1: the host
insertion delta-M built ONLY from booked V291 energies on host
channels (if n0 = 0: delta-M = 0 exactly); re-diagonalized angles
theta12/23/13 booked; deviation from PDG booked.
**Branch.** MIX-UNCHANGED (||delta-M|| < 1e-9: TBM stands, the PDG
gap is the standing named falsifier) / MIX-SHIFTED (angles booked,
fitted-flag false, no dial touched).

## V293_r_quantization (tier 3)
**The question.** Is hosted family multiplicity quantized: does the
chain-kernel hosting from V291, counted through the booked V70
family phases (families 3, per_family 4, nodes 12), give an integer
per family?
**The object.** the V291 artifact (kernel basis channel census),
the V70 artifact families/counts anchor (re-derived: phase set
match to 1e-9), the V91 v70_regression anchor.
**Measured.** per booked family: hosted-channel count r_f (integer
trace over the kernel projectors); defect |r_f - round(r_f)| booked.
**Branch.** R-QUANTIZED (defect < 1e-9, r_f booked against V70's
3) / R-NONINTEGER (defect booked honestly).

## V294_democratic_mD_forced (tier 3)
**The question.** Do the booked integer charge selection rules
(atom channels {0, 64, 96} from the V269 census; the V197 rung
factor arithmetic 208 = 2^4 * 13) FORCE the democratic rank: what
is the maximum rank a 3x3 family mass matrix can have while
commuting with the channel-degree selection?
**The object.** the V269 c2_spectrum anchor (integers), the V197
ladder_rows anchor, the V45 mass_matrix anchor, the V70 family
operator set.
**Measured.** the general 3x3 pattern allowed on the family space
by the selection rules (commutant computation from the booked
family phases); its maximum rank and spectrum booked against the
democratic {3, 0, 0}.
**Branch.** RANK-FORCED-1 (commutant rank-1: democratic exact by
charge bookkeeping) / RANK-FREEDOM (maximum rank booked; the
rank-1 fact becomes an input, stated as such).

## V295_triangle_selection (tier 4)
**The question.** The unitarity triangles of the booked mixing: are
they degenerate (theta13 = 0 in TBM) and how does the booked
Jarlskog J sit against the V109 measured J_available_max?
**The object.** the V45 mixing anchor, the V292 verdict artifact,
the V109 artifact (asymmetry/J_available_max keys).
**Measured.** J = Im(u11 u22 u12* u21*) from booked mixing; triangle
areas; the V109 J_available_max read verbatim.
**Branch.** TRIANGLES-DEGENERATE (J < 1e-12 booked; the nonzero-PDG
triangle is the standing named falsifier) / NONDEGENERATE (J and
angles booked).

## V296_flux_baryogenesis_ledger (tier 3, baryogenesis bookkeeping)
**The question.** Combining booked artifacts only (V55 eta_cp, V106
washout audit, V125 cp-flux washout: asymmetry, flux_scan,
J_available_max): does ANY booked flux branch leave a final
asymmetry at or above the baryon-asymmetry floor 6 x 10^-10, or is
every booked branch washed out?
**The object.** the three artifacts verbatim; no new dynamics.
**Measured.** per booked branch: production x booked washout
factor -> final fraction; the maximum over branches booked.
**Branch.** FLUX-BG-OPEN (a branch >= 6e-10 with its chain of
booked factors listed) / FLUX-BG-WASHED (max fraction booked).

## V297_trace_action_ledger (tier 4)
**The question.** Over the booked V283 stable set (METRICLOOP-
NEGATIVE rows, n <= 30, step-doubling stable), does ANY pair,
converted through the certified V276 ladder (F_208 certified read,
coefficient_observed), put the induced action coefficient within
2 % of the booked Planck relation - the action-level reopening of
the loop-level null?
**The object.** the V283 artifact rows, the V276 artifact ladder
reads, the V108 coefficient anchor; arithmetic only.
**Measured.** min over booked pairs of |action coefficient / booked
relation - 1|; the winning pair booked.
**Branch.** ACTION-COUPLING-FOUND (dev < 2 %) / ACTION-NULL (min
dev booked; the ladder stays the only certified route).

## V298_rung_portal_lambda_B (tier 4)
**The question.** Is the booked singlet-Higgs portal lambda_B =
4.06 x 10^-14 (the freeze-in value reproducing Omega h^2 = 0.120)
connected to the certified rung arithmetic (V197: the 208
factorization, ladder_rows), or is it a named input independent of
the rungs?
**The object.** the V197 artifact rows, the V108/V276 coefficient
anchors, the booked lambda_B as carried in the publication text
(paper3 S DM section - the decimal is the published input, not a
restyped artifact number; its source is booked as named input).
**Measured.** the rung-side candidate: products of booked rung
overlaps and the certified hierarchy coefficient; its deviation
from lambda_B.
**Branch.** PORTAL-RUNG-FIXED (dev < 1 %: lambda_B inherits
certification) / PORTAL-NAMED-INPUT (dev booked: the portal remains
a named input, honest label kept).

## V299_multicell_antiunitary (tier 3)
**The question.** Does the materialized chain carry an antiunitary
symmetry: is there a unitary C on the 12-dim chain-kernel slice
with C H_bar C^dag = H (H the booked chain-H data materialized on
the kernel basis, V272 chain_ground_h_re/im anchor), and does
C^2 = +-1?
**The object.** the V291 kernel basis artifact, the V272 booked
chain H data, the cl10.cpt library conventions (read, not assumed).
**Measured.** existence of C via the linear commutation system on
the slice; min residual; C^2 eigenvalues booked.
**Branch.** CPT-EXACT-ON-CHAIN (residual < 1e-10, C^2 = +-1 booked)
/ CPT-ABSENT-ON-CHAIN (min residual booked).

## Roll-up consequences of this programme (fixed before any run)
- V290 executes only behind its pregate; otherwise the roll-up books
the skip verbatim from V279.
- No paper wording changes are pre-authorized by this phase; all
branches feed the ledger and Table 9 rows only.
- The suite recount (206 -> 206 + executed count) propagates at the
final consolidation exactly as before.

### V287 T0 pole-split contract amendment (fixed after the machinery
probe, before any result is booked)

The machinery probe established that at GAMMA the operator
Hgen = herm(i log M(0)) carries 64 modes at the H = 0 level (the
logarithm of a defective unitary step; scipy's own approximate error
0.42 booked). V229's T4 classification occ = om < 0, unocc = om > 0
therefore splits pure float noise, and the booked bubble eigenvalues
are ONE draw of that split noise: re-derivations on the current
machine give 69/59 and 68/60 splits and bubble eigenvalues ranging
from about (-5.3, -4.1, -1.3) to (-24.2, -22.9, -6.4) on the star
family - orders of magnitude beyond the pre-registered 1e-6 relative
tolerance. The 1e-6 re-derivation contract is not satisfiable by ANY
implementation; it chases a noise draw, not a deterministic fact.

The question (does an O_h-covariant vertex family null the quintet
bubble?) is unaffected and is executed as follows. (i) T0 now books
the split-noise envelope: the bubble eigenvalues over 8 seeded sign
assignments of the exactly-zero modes (seed 287) plus the natural
split, per eigenvalue; T0 passes exactly when the booked V229
eigenvalues lie inside the booked envelope - machinery re-derivation
agrees with the booking up to the named noise manifold. (ii) T1 and
the family search run on a FIXED deterministic pole prescription,
booked in the artifact: zero-mode signs by continuity at q = 1e-3
along the 12 directions (+-e1, +-e2, +-e3, and the six signed
(1,1,0) diagonals), matched by eigenvector overlap, majority vote,
ties to unoccupied; the level at H = 0 is exactly degenerate, so
pairs with |den| < 1e-9 (intra-level) are EXCLUDED from the
polarization sum (principal value - no first-order mixing inside an
exactly degenerate level), the excluded-pair count booked. The
prescription is then noise-safe by construction: every surviving
denominator is O(10^-3) or larger. (iii) The 1e-8 nullity gate, the
80-start
Nelder-Mead budget, and every downstream threshold are unchanged.
No verdict was booked from the probe; the envelope and prescription
are machinery bookkeeping, not tuning.

## V290 pregate resolution (fixed before any V290 run)

The committed V279 artifact (branch control:CANDIDATE-MASSIVE) books
exactly one candidate: the control star-module two-fold line at q =
(0, 0, pi/2), gap 9.888e-17, transverse exponent 0.99818, with a
commutant obstruction 0.0216 that opens the touching linearly
(slope 0.04641). The G1 candidate therefore EXISTS, and the V279
pregate is met: V290 EXECUTES. (The programme plan's "3756 = 533 +
3223" census figures match no artifact in the tree - they are a plan
transcription error, not a measurement, and nothing in this run
depends on them.) The candidate's own transverse exponent from the
V279 artifact is carried verbatim; V290 measures what the pre-reg
asks - whether the candidate's cone anisotropy runs irrelevant under
the booked V115 anisotropy metric as the step is coarse-grained
(dq = 1e-2, 3e-3, 1e-3), on the decorated 128-dim star-bond clock of
V277 at its booked best point (the only decorated gapless clock on
record) and on the bare clock as control. Every anisotropy read is
booked; IRRELEVANT requires A to fall below 0.05 at the finest grid
AND to decrease monotonically across the coarse-graining series on
the decorated clock.

## Roll-up consequences of the V277-V299 programmes (fixed after the
## runs, before the reference/number edits; one line per verdict)

- Phase 1 (light node): V277 DECORATION-BEST - the 18-parameter
  star-decoration family reaches a gapless closed-clock point with
  cone degree 1.179 (acceptance legs gap 1.1e-16 and degree stand; the binding legs vmax 1.586 > 1 and rho_comm 0.275 < 0.5 are booked); V278
  CONE-RANGE2-ABSENT - the 30-parameter range-2 parity-decoration
  search (150 starts) keeps the Gamma gap exactly closed
  (1.1e-16, 64 zero modes) but never mixes the fourth gap by
  parity (top4 gap 4.1e-16) and never brings the anisotropy
  below 0.1; both binding legs booked at the best point; V279 control:CANDIDATE-MASSIVE - the single
  booked candidate (star-module two-fold line) is gapped by its own
  commutant (obstruction 0.0216, linear opening slope 0.0464).
  V280 was skipped by its pregate (zero cone-family candidates in
  the V277 acceptance sense).
- Phase 2 (loop origin): V281 FLUXLOOP-NEGATIVE, V283
  METRICLOOP-NEGATIVE - both traced-loop routes floor at O(10^2-10^1)
  against the 1.99e-2 target band; V282 CUBIC-PRESENT-FROZEN - the
  cubic is real in the module bookkeeping (Sym^3 dim 16215 exact,
  dims 6/3/0) but has zero dispersion under the cone twist; V284
  CROSSING-ARITHMETIC-FINAL - n* = 208 stands as a measured fact,
  no arithmetic relation (PSLQ empty, factor census negative). The
  gate G2/D1 paper edit is NOT triggered (the ladder clock did not
  move). O8 stays open with the three named negative routes booked.
- Phase 3 (graviton cone): V285 CONE-SELECTION-ABSENT (booked
  in-tree), V286 LINKINV-HOP - the module-charge commutant on the
  256-dim link module, computed exactly as the C2-spectral span of
  its three multiplicity-free sectors (levels 0/64/96, dims
  1/45/210, Schur dim 3), is the Casimir-filter triple with max Gauss
  leak 4.2e-16; the induced 12x12 channel operators have widths mu =
  1.0000/0.1491/0.0690 and the BZ velocity supremum is exactly
  w8*mu (the stage matrix is a scalar multiple of a fixed operator;
  finite-difference grid cross-check 0.9511*w8*mu); over V239's four
  closure clocks exactly one band entry reaches |v| = 1 - the
  singlet filter at the (pi/2, pi/2) clock - which carries no new
  gauge-charge content, and the joint gate of the pre-registration
  (V285 selection AND V286 v = 1) is NOT met since V285 returned
  ABSENT, so the §11.8 amendment sentence is not inserted; the
  random-generator block-commutant route was tried and rejected as
  methodologically wrong (booked in the script comment); V287
  BUBBLE-NULL - with the booked
  pole-prescription amendment the stress-quintet bubble vertex
  search closes at min ||Pi||_F 1.7e-15 on the axes family: the
  quintet bubble is NULL on the decorated clock, and the T2g vertex
  leg reproduces V229 (axis 0.0, star 3.81); V288 BEKENSTEIN-MAP -
  every flux truncation of the three-site kernel keeps the bond
  entropy at exactly ln 45 (factor 15.227 from 1/4): the truncation
  map moves the ground degeneracy (2/4/6) but not the factor - the
  1/4 stays a named input; V289 MC-TUNED - the full pre-registered
  sweep (5 betas x 32 replicas x 20000 sweeps, seed 289) is finite
  everywhere and tunneling is sampled: acceptances 0.863/0.693/0.290/
  0.060/0.045 at beta = 1/2/4/8/16 with W1 = 0.024/0.062/0.251/0.648/
  0.642; the area-law fit exists at beta = 4/8/16 (log W vs R^2 slopes
  -0.48/-0.43/-0.39, decay measured), so the [4,16] range books the
  decayed regime; character cross-check exact over 606720 sampled
  plaquette products (max defect 0.0); the earlier interim MC-NULL
  entry reflected the first run's zero-window crash stage, superseded
  by the amended contract above; V290
  ANISOTROPY-IRRELEVANT - executed (the V279 candidate exists): the
  booked candidate is gapped on the 128-dim clock (m_next 0.323) and
  the anisotropy of that gapped floor falls linearly in r (A =
  0.13/0.039/0.013), so no anisotropy-protected mass route opens on
  the clock from this candidate; the inconclusive leg is booked with
  the floor note.
- Phase 5 (flavour): V291 HIGGS-HOST-PRESENT (n0 = 3 of the 12
  chain-kernel states carry module-singlet-channel weight > 1/2:
  0.898/0.623/0.514 - the chain does host neutral scalar states;
  float64 null-slice weights after the machinery fix, kernel dim 12
  and Gram 4.9e-15 as anchors), V293 R-QUANTIZED (total singlet
  weight over the kernel = 4.000000, defect 0.0 - a projector trace:
  exactly integer; r_int = 3 majority hosts booked alongside),
  V292 MIX-UNCHANGED - the host insertion is proportional to the
  democratic generator, so the frame-invariant sin^2 th12 shifts by
  2.2e-16 and the O(2) light family is untouched (the first run's
  MIX-SHIFTED label was an arbitrary-degenerate-basis artifact, fixed
  by computing only frame-invariant content - heavy-eigenvector
  weight + degeneracy test); PDG 0.307 gap booked at 0.0263 stands;
  V294 RANK-FREEDOM - the S_3 family-selection commutant
  on 3x3 is span{I, ones} (dim 2, maximum rank 3): the democratic
  rank-1 is NOT forced by the charge bookkeeping, it is a named
  input; V295 TRIANGLES-DEGENERATE - J = 0 exactly on the booked
  angles (sin^2 theta13 = 0), all unitarity-triangle areas vanish,
  consistent with the booked V109 search.
- Phase 6 (cosmology): V296 LEDGER-INSUFFICIENT - every booked
  asymmetry-like magnitude (key scan eps/asym/eta/epsilon/Jarlskog
  in V55/V106/V125) sits at least a factor 8.0e3 from 6.1e-10; the
  derived-product row of the ledger (booked explicitly in the
  artifact) shows the only quantity landing NEAR eta_B is V55's
  guessed-ingredient product, ratio 1.11 unwashed and 0.13 under the
  booked BDP washout - landing by construction, not by mechanism:
  the baryon-asymmetry generation stays an open problem
  (the ledger is booked in full); V297 TRACE-ACTION-ABSENT - the
  closest certified ladder row sits 46.5 % from the Planck-side
  ratio, no induced-action match; V298 PORTAL-NAMED-INPUT - the
  nearest certified portal candidate sits at 10^{+11.2} of the
  booked 4.06e-14: lambda_B remains a named input.
- Phase 7 (interpretation): V299 CPT-CHAIN-VIOLATION - the booked
  V272 materialized chain ground-state vertex has imaginary parts up
  to 0.171 (over 6144 entries): the chain's effective vertex is not
  real, so no antiunitary acts exactly on the multicell chain at the
  level of the booked chain data. This does NOT touch the V52
  spectral-mirror CPT statement (spec(U|16bar) = spec(U|16)*), which
  is a unitary-spectrum statement, not a reality-of-chain-vertex
  statement; the paper's CPT section already says exactly that.

### V289 zero-acceptance contract (fixed after the first full run's
### diagnostics, before the result is booked)

The first full GPU run died in the area-law fit stage on a
zero-window beta: with Metropolis acceptances of 0.015-0.48 at
beta = 1-16 and the character cross-check exact (0.0 over 4032
plaquette products), several beta rows carry W2 = W3 = 0 exactly -
no tunneling between the two classical geometries was ever sampled.
The log-linear area fit is then undefined, and log(0) crashed the
stage. Amended contract: the fit stage runs only on rows with all
three windows strictly positive (booked count); the T1 machinery
check is finiteness of the windows (0 is finite - a physics result,
not a machinery failure), and the area-law verdict is carried by the
branch: MC-GEOMETRY-OK / MC-TUNED as pre-registered, MC-NULL when no
beta admits the fit, with every row's (beta, acceptance, W1, W2, W3,
statistical errors) booked. No beta is dropped from the ledger and
no proposal is retuned; the proposal kernel and sweep counts are the
pre-registered ones.

# The geometry measure and the flavour host: pre-registrations V300-V302
## V300_cone_measure_selection (tier 2, hard, F8/O9)
**The question (pre-registered).** The paper states the finite measure
on the geometry group (V289, MC-TUNED) "does not measurably select"
the light cone. Does that measure in fact concentrate a light cone:
on the sampled ensemble, is the screening mass of adjoint line
correlators isotropic in lattice units - an effective cone speed
v_cone = xi_1/xi_0 within [0.95, 1.05] - or does the measure prefer no
cone at all (v_cone far from 1 at every beta with acceptance > 2 %)?
**The object (closed input set).** V289's machinery re-derived live
(V141 build_group order 1024, the 45-dim adjoint rep built
multiplicatively along the generation tree, the full
multiplication/inverse/proposal tables, the adjoint-character Wilson
action S = beta sum_p (1 - chi_45(U_p)/45), single-bond biflection
Metropolis): the L = 4 protocol is replayed for the first two beta
(1, 2) of the original seed-289 draw sequence and their booked
(accept, W1) rows reproduced to 1e-6 (T0; the character
cross-check's exactness re-measured); then a new ensemble at L = 6
(216 sites, 648 bonds, 648 based plaquettes), beta in {2, 4, 8, 16},
NREP = 16, SWEEPS = 20000 (EQUIL 2000, SAMPLE_EVERY 20),
BASEPOINTS = 32, seed 300.
**Measured (pre-registered).** T0 row replay; the line correlators
C_mu(d) = mean_s0 Re chi_45(u_s u_{s+mu} ...)/45 along straight
forward lines of length d = 1, 2, 3 in direction mu = 0 (the step
direction) and mu = 1, averaged over samples, jackknifed over the
(replica, base-point) product; m_mu from the 3-point log-linear fit
(defined iff all three C > 0 - the V289 fit-stage contract applies);
v_cone(beta) = m_0/m_1 with its jackknife error.
**Acceptance (pre-registered).** CONE-MEASURE-SELECTED iff
|v_cone - 1| < 0.05 at some beta with accept > 2 % and both masses
defined; CONE-MEASURE-ANISOTROPIC if masses are defined but no beta
lands in the band (the whole v_cone(beta) curve booked, with the
implemented generator speed 1.570796 read from the V285 artifact as
the comparison number); CONE-MEASURE-NO-MASS if no beta admits
defined masses (the correlator table booked). Pass = the measurement
is complete in every branch.
## V301_geometry_gap_ladder (tier 2, hard, O8-new)
**The question (pre-registered).** Can the non-perturbative screening
mass of the geometry statistical ensemble - the only scale the
geometry sector demonstrably generates - sit at the ladder loop base:
is m_gap (lattice/Planck units, ell_P = 1 by Axiom 1) inside the
ladder band [1.5e-2, 3.0e-2] around 1/(16 pi) at any beta?
**The object (closed input set).** V300's machinery at L = 6 with
beta in {4, 8, 16, 32}, seed 301 (the extra weak-coupling point
chases the smallest mass); T0 = the V289 beta = 1 protocol-row replay
plus the booked area-law slopes read from the V289 artifact.
**Measured (pre-registered).** m_gap(beta) = the direction-1
line-correlator mass (log-linear on d = 1..3, defined iff all C > 0)
with jackknife error; the Wilson-loop mass m_W(beta) = -ln<W(2)>/4 as
the independent cross-read; direction 0 booked too; the band crossing
(beta*) if any; the minimum m_gap over the booked beta range.
**Acceptance (pre-registered).** GEOMGAP-IN-BAND iff some beta has
m_gap in [1.5e-2, 3.0e-2] (beta*, m_gap +/- jackknife, m_W cross-read
booked); GEOMGAP-OFF-BAND otherwise, the curve and the distance to
the band booked. A GEOMGAP-IN-BAND hit books a scale coincidence
only: the identification of the screening mass with the ladder base
remains the declared step O8 already states, and no reinterpretation
is added by this run.
## V302_mD_host_forcing (tier 2, hard, O2)
**The question (pre-registered).** V294 measured rank-1 of the
democratic m_D to be symmetry-free (the S3 commutant contains
rank-3). Is rank-1 forced dynamically instead? The chain hosts
(V291, HIGGS-HOST-PRESENT: n0 = 3, p0 > 0.5 at chain-kernel columns
6, 7, 9) carry booked chain-H energies; if the three host energies
were degenerate the democratic (uniform) host coupling would be the
energy ground pattern - rank-1 consistent with host energetics - and
if they are not, the ground pattern of the host sector is not the
three-family democratic ones.
**The object (closed input set).** The booked V291 artifact (p0,
zero_idx, chain_H_spectrum_artefact) and the on-disk
kernel_basis.npy (6144 x 12, complex64), the V272 chain-kernel
anchors (kernel dim 12, the 12-entry chain_H_spectrum), V292's
host-shift construction (c_host_energy_sum 112, mass_rescale 37.333)
and V45 democratic anchor, and the V294 commutant anchors. The saved
basis is re-verified (Gram = I to storage precision), not re-derived:
no Nmat reassembly.
**Measured (pre-registered).** T0: p0 recomputed from the saved basis
(reshape (6, 1024, 12), slice weights summed over the 1024 factor,
channel-zero slices per the V270 channel anchor) vs the booked p0
(tol 5e-6, storage precision); n0 = 3, hosts [6, 7, 9]; the hosts'
chain-H energies vs V272's booked spectrum. T1: the host degeneracy
deltaE = max E - min E; the democratic uniform-combination energy
(mean, anchored to V292's 37.333); the energy-minimal rank-1
occupation over host subsets, with the two-host pair's mixing angle
(sin^2 theta_12 = 1/2 vs the booked measured 0.307 - booked as a
gap, never refitted) against the three-host democratic; the energy
penalty of the full democratic pattern booked.
**Acceptance (pre-registered).** MD-HOST-DEGENERATE iff deltaE <=
1e-6 (the uniform ones coupling costs nothing; rank-1 consistent
with the host energetics - still not forced, consistency booked);
MD-HOST-ASYMMETRIC otherwise, with deltaE and the ground occupation
pattern booked: rank-1 democratic m_D then stands declared for a
second, independent reason (neither symmetry - V294 - nor host
energetics select it). Pass = both legs measured.
## Roll-up consequences of this programme (fixed before any run)
- Only V300, V301 and V302 interpret their own runs; no acceptance
moves after a run; any dial touched reclassifies to fitted (none is).
- CONE-MEASURE-SELECTED is the selection leg the paper says V285 and
V286 do not supply. It changes the Section 11.6 / Table 8 F8 wording
("candidate under selection") ONLY jointly with a v = 1 branch
carrying the implemented hop - V286 already booked its selection leg
as the one V285 does not supply; alone, V300-SELECTED books "the
geometry measure selects a light cone speed; the implemented hop
still does not realize it" and F8 stays fired.
- GEOMGAP-IN-BAND books a scale coincidence, never a mechanism; the
O8 declared-input wording survives.
- MD-HOST-ASYMMETRIC strengthens the O2/O6 honest ledger (rank-1
declared on two independent grounds).
- One dated reference.md paragraph per experiment; number_audit
PASSes before each commit; the suite re-runs once at the roll-up
(215 experiments expected).

## Re-gating convention change (2026-10-29, ordered by the reader)

The full-suite re-gate at every roll-up is replaced by the fingerprint
gate, effective immediately and for this programme's roll-up (the
V300-V302 roll-up consequence "the suite re-runs once (215 expected)"
is met by that gate, not by a full re-execution). Mechanism, booked in
the reference.md protocol section: run_validation.py records the
SHA-256 of each experiment's scripts and of the cl10 library plus
cl10/config.yaml in results/validation_fingerprints.json; a changed
or missing script, a missing/failed/crash-stale verdict re-runs that
experiment alone; a library change re-runs the whole suite; --full
re-executes everything (a fresh clone does this by construction). The
honesty invariant is preserved: every verdict folded from disk is one
the current tree's own code produced - the fingerprints are what make
"on disk" mean "by this code" - and the hard-tier verdict, the
number_audit gate and the dated reference paragraphs are unchanged.

# The neutrino truth programme: pre-registrations V303-V305

Reader instruction (2026-10-29): "I want the physics to be true as well."
The one sector where the theory already carries a value inside the measured
data (the 1/3 candidate, Section 12.3) receives its last three experiments.
Nothing here may be refitted to land in band; any dial touched reclassifies
to fitted. The candidate's standing numbers are re-read from booked
artifacts only: S17 = loop_ladder(17) (V103's source), the NuFIT scales and
the 0.0642 eV ceiling (cl10/config neutrino constants via V103), the
democratic matrix and mixing (V45), the commutant span{I, J} (V294), the
host channels and weights (V291), the family-symmetric insertion point
c = 112 and angle shift 2.2e-16 (V292).

## V303_junction_energy_r (tier 2, hard, O2)
**The question (pre-registered).** The junction ratio r = kappa/w0 is not
derived by the shape (V114: the stabilizer-symmetric closed family is
2-dimensional) and not quantized (V293 books the family multiplicity, a
different r). Is it selected by ENERGY: on the stable branch 0 < r < r* of
the V127 quadratic junction theory, does the measured vacuum energy have an
interior minimum - and does the minimum, or the energy's flow, land on the
democratic onset r* = 2 (the candidate's anchor) inside the cosmological
window r > 1.278?
**The object (closed input set).** V127's machinery verbatim (NODES =
cl10.supercell.M_POINTS, the support set, the closing triangle SEL, the
adjacency, the one-loop functional delta_E(S, r) = sum_k log|1 + (r/2)
mu_k(A_S)|^2, the 4096 split sets). Convention fixed before any run and
booked: both functionals run in the canonical junction ratio r = kappa/w0
(levels 1 + (r/2) mu_k); V127's quadratic functional e_quad uses a doubled
convention (1 + r mu) while its crossing search and one-loop use the
canonical one - V127's sweep points (0.5, 1.0, 3.0) are canonical r values
in its one-loop and selection tables, and none of them enters the stable
branch r < 2: the stable branch has never been swept, and that is what
this experiment measures. No dial, no optimisation: grids only.
**Measured (pre-registered).** T0 reproduction: the V127 artifact anchors
(8-member cubic orbit, K3 degrees [2,2,2], 12 edges, 0 exiting, junction
group order 12, trace theorem < 1e-9, first crossings {2.0}); drift halts.
T1 the closing-triangle functionals on r in (0, 2): the tree-level zero-
point sum is r-FLAT (trace theorem: measure the deviation), the one-loop
delta_E is strictly monotone (sign-count the derivative; interior minimum
count booked), and its exact closed form 2 log(1+r) + 4 log|1 - r/2| is
verified to 1e-12. T2 the full 4096-set sweep on the stable branch: per-set
first massless crossing r_S = -2/mu_min (buckets below 1.278, below 1.9191,
below 2, none), and at every r of a 0.05 grid the energy minimiser among
finite sets (its size, its multiplicity, and whether the closing triangle
or its junction orbit ever attains it). T3 the contact: the attractor r ->
r* = 2 gives m_3 = S17/3 (re-read S17) against the booked atmospheric
scale M_ATM_EV - the relative residual booked as a named number; the
window r > S17/SUM_MNU_MAX_EV - 1 re-derived.
**Branch (pre-registered).** R-ONSET-DEMOCRATIC (no interior minimum, the
energy flows to the massless onset r* = 2) vs R-WINDOW-MIN (an interior
minimum is found and booked) vs MACHINERY-DRIFT; crossed with
VACUUM-SELECTS-ALL-12 (the 12-node split set minimises at every window r)
vs TRIANGLE-SELECTED (the closing triangle attains the minimum for some
booked r range, range booked); the atmospheric residual is booked in the
verdict in either case. The r mechanism, if found, is an energy attractor
at a boundary point (a marginal selection), and this is said in the
verdict text: a boundary attractor with a 2.7 %-class residual against
the atmospheric scale is NOT a landing on the measured value.

## V304_realizable_mixing_census (tier 3, O2)
**The question (pre-registered).** V104 scanned GENERIC vev directions and
found joint NuFIT fraction 0.0; V292 measured the ONE realizable insertion
and found the angles unmoved (shift 2.2e-16). Is the mixing correction of
the ENTIRE realizable sector exactly {0}: every neutral chain-kernel host
channel of V291 (indices {6,7,9}, weights p0), all 7 non-empty combinations,
with the family-symmetric insertion justified by the re-derived booked
facts (V45 family_symmetry equivalent_3fold, V294 commutant span{I, J})
rather than asserted?
**Measured.** T0: V45 democratic spectrum re-read (rank 1, eigs {0, 0, 3}),
V292's point re-run (angle shift < 1e-12, c = 112 from the booked
spectrum). T1: the census - correction c(W) = sum_{i in W} p0_i x (booked
chain-H energy) over the 7 subsets, M(c) = (3 + c) x ones: eigenvalues
{3+c, 0, 0}, the heavy eigenvector stays (1,1,1)/sqrt(3): the maximal
angle shift over the census booked (expected < 1e-12: the whole realizable
correction lives in span{J} by the re-derived commutant). T2: the contact
of the locked point (sin2 12 = 1/3, sin2 23 = 1/2, sin2 13 = 0) against the
V104-borne NuFIT dictionary (imported from V104's module, not retyped):
sigma per angle booked; joint fraction over the census is 0 exactly if any
angle sits beyond its band.
**Branch.** MIXING-REALIZABLE-EMPTY (max shift < 1e-12 and joint fraction
0.0: no realizable scalar moves any angle; the 1/3 -> 0.307 fixup requires
a declared direction) / MIXING-REALIZED-N. The theta13 sigma is booked
either way: the democratic point carries sin2 13 = 0.

## V305_onset_neutrino_contacts (tier 3, O2)
**The question (pre-registered).** Under V303's branch, what are the exact
laboratory and cosmology contacts of the energy-selected (or window-
minimum) spectrum - and which of the published reach claims survives with
the theory's OWN mixing versus measured mixing?
**Measured.** r_sel from V303's booked branch (onset r = 2 or the booked
minimum; the window-atmospheric point 1.9191 re-derived for contrast);
spectrum (0, M_SOL_EV, S17/(1+r_sel)); Sigma m_nu against the booked
ceiling (margin booked); the beta-decay effective mass m_beta under the
measured-mixing read (V292's booked obs_pmns fractions) and under the
model-native TBM sector (|Ue2|^2 = 2/3, |Ue3|^2 = 1/3); the neutrinoless
interval m_bbeta under both (TBM: [2/3 m2 - 1/3 m3, 2/3 m2 + 1/3 m3];
hybrid: measured fractions) - the published claim "m_bbeta ~ 0.02 eV, in
reach" is re-scored: it is a TBM-sector value only, and whether the hybrid
interval lies below the same reach is booked; F3 is re-booked in full
(bare ladder sum excluded, candidate status per the sums).
**Branch.** ONSET-SPEC / SUM-UNDER-BY-m% (or SUM-OVER) / BB-TBM-ONLY (if
the hybrid interval falls below the TBM interval where the reach claim
lives) - all numbers booked from the re-reads.

## Roll-up consequences (fixed before any run)
- The paper's Section 12.3 sentence chain gains the V303/V305 verdicts: if
  R-ONSET-DEMOCRATIC, the r = 2 anchor is reclassified from "convention"
  to "energy-attractor with a named atmospheric residual" and the
  candidate's mechanism sentence changes accordingly - the candidate
  remains a candidate if the residual stands. The m_bbeta reach sentence is
  corrected to the V305 branch if BB-TBM-ONLY fires.
- O2 stays open iff the sums are inside the ceiling; the mixing stays
  declared if V304 books EMPTY; every named input keeps its name.
- The suite is re-gated once by the fingerprint gate (3 new experiments,
  218 expected); number_audit re-runs; PDFs rebuild if the paper changes.

### Post-run booking note (2026-10-29, after V303-V305 ran)

Branch-name deviations recorded, no physics re-fit:
- V303 measured R-ONSET-DEMOCRATIC (the pre-registered branch fired as
  worded) but the crossed alternative came out as neither pre-registered
  option: the 12-node split set minimises only below its own 6-cycle
  crossing at r = 1 (20/39 grid points), and inside the window the
  minimiser fragments (sizes 9-11). The verdict word
  "VACUUM-SELECTS-ALL-12" would have been false; the emitted word is
  VACUUM-NEVER-TRIANGLE, the honest generalization - the pre-registered
  physics question (does the energy select the closing triangle?) is
  answered NO as planned. The convention fix (canonical r in both
  functionals, grid capped at 1.98 for the closed-form check) was booked
  in the object text before the run.
- V305's published-claim re-score used the V45-booked O(2) family rather
  than the TBM label the pre-registration sketched: the artifact books
  BB-ENVELOPE-DECLARED, which is the TBM-ONLY alternative's honest
  generalization (the claim lives inside the free family, not at a
  predicted point). V304 fired MIXING-REALIZABLE-EMPTY as pre-registered.

# The atmospheric-gap programme: pre-registrations V306-V307

Reader instruction (2026-10-29): "fix it in the next programme" - the
-2.70% atmospheric residual of the energy-selected point (V303). Two
closed experiments decide whether the gap is a SUPPLIED number of the
root or a declared one. Nothing refitted; the targets are computed from
booked constants only: S17 = loop_ladder(17) x EV (V103's source),
M_ATM_EV = 0.0501, r_atm = S17/M_ATM_EV - 1 = 1.9191154379...,
Delta-r = 2 - r_atm = 0.0808845620..., and the multiplicative form
rho_req = M_ATM_EV/(S17/3) - 1 = 2.77288e-2.

## V306_coefficient_census_atm_gap (tier 2, hard, O2)
**The question (pre-registered).** Does ANY booked constant, or a bounded
combination of two, supply the gap coefficient to 1e-3 absolute? The
class is fixed before the run. Base constants B, every one read from its
artifact/config (no retyped decimals): LOOP = 1/(16 pi) (config);
LOOP^k, k = 2,3,4; alpha_10 = 0.04 (config clock); alpha_X_root = pi/64;
g_obs = 0.19999798835751165 (V206); 1 - rho with rho = 0.9780302704295731
(V276 certified machinery); F_ratio - 1 with F_ratio = F_cross/F_10_open,
F_cross = 0.01968515418112789 and F_10_open = 0.019186391026184392
(both V276); F_10_open/LOOP - 1; the certified crossing deficit
1 - F_cross/LOOP = 0.010527; LOOP - F_10_open; clock fractions
{1/2, 1/3, 1/4, 1/5, 1/8}; the certified crossing n* = 208 and 1/208
(V276). Candidate expressions: every a, a +/- b, a*b, a/b for
a,b in B (self-pairs and self-ratios excluded), evaluated against BOTH
target forms (rho_req multiplicative; Delta-r in r units). Booked:
the five smallest residuals with their expressions, and the exact
minimum.
**Branch (pre-registered).** R-COEFFICIENT-SUPPLIED if min residual
<= 1e-3 (the supplying expression booked, and it then owes a derivation
test before any reclassification - a two-number coincidence is not a
mechanism); R-COEFFICIENT-NOT-SUPPLIED otherwise, with the class size
and the closest misses booked as the ledger's bounded no-go.

## V307_onset_runaway_structure (tier 3, O2)
**The question (pre-registered).** Is the r* = 2 selection an attractor
(marginal, flat) or a RUNAWAY edge (energy unbounded below)? V303's
verdict word "attractor" is re-scored against its own functional.
**Measured.** With the V303 functionals: (1) the K3 pair's one-loop
Delta E on r -> 2 - 1e-6 and the divergence coefficient: the number of
critical modes crossing exactly at r = 2, counted as ||1 + (r/2) mu_k| <
1e-4 at r = 2 - 1e-6 over the full 12-set and over the component
structure [3,6,3]; (2) the displacement cost at the measured point:
C_atm = min_S Delta E(S, r_atm) - Delta E(2 - 1e-6) over the 4096 sets,
in the same log-det units; (3) ONE-LOOP EXACTNESS of the onset
location: the log-det singularity sits at the tree-level zero - the
tree level is r-flat (trace theorem, V303 T1) so no r-dependent
counterterm exists to move it: the shift |r_sing - 2| measured to the
float64 floor, and the r-flatness deviation of the tree level re-measured
on the full 12-set (not only the triangle).
**Branch (pre-registered).** ONSET-RUNAWAY-EXACT (unbounded below,
singularity at 2 exact, cost C_atm booked): r = 1.9191 is then a
displaced point of measured cost, the stabilization mechanism is not
booked, and r stays declared - the V303 wording "boundary attractor" is
corrected to "runaway edge" in the ledger and paper; ONSET-MODERATED if
the minimised energy develops an interior minimum on the fine grid
(booked, reopens the selection route).

### Post-run booking note (2026-10-29, after V306-V307 ran)

- V306 fired the pre-registered R-COEFFICIENT-SUPPLIED word at the 1e-3
  tolerance, but every satisfying expression is a two-number coincidence
  (best: (LOOP - F_10_open)/(F_cross/F_10_open - 1) at 4.7e-4;
  pi/64 - (1-rho_V276) at 5.9e-4; 2 x alpha_10 against Delta_r at 8.8e-4
  - leaving 1.1% of the gap unfixed; (1-rho) + 1/208 at 9.3e-4). None
  carries a derivation test - the branch's own condition ("a
  two-number coincidence is not a mechanism") keeps the gap declared.
  The bounded no-go for anything better than 1e-3 over the reserved
  class is the booked content.
- V307 fired ONSET-RUNAWAY-EXACT with the cleanest values of the
  programme: singularity at r = 2 to 0.0 float64, tree flatness on the
  full 12-set 0.0, J12 overlap with the 6 zero modes exactly 0,
  candidate-vacuum inverse factor exactly 4.0. Consequence booked in the
  paper's Section 12.3 wording: "boundary energy attractor" (V303) is
  corrected to "runaway instability edge with an unbooked stabilizer";
  and the candidate's triangle-only vacuum is the one the junction
  energy does NOT select - on the minimising vacuum the candidate reads
  4 x S17/3 = 0.195 eV, excluded. No physics was re-fit; the V303
  verdict word stands in the artifacts as measured then, corrected by
  the V307 paragraph as the convention requires.

## Gate-hygiene note (2026-10-29): the V247 CPU-fallback control

The full 220-gate run under the jobs-8 pool re-ran V247's subprocess on
the numpy fallback (GPU contention); its T0 ladder control measured
2.1e-02 from the published cupy-built constants - outside the 1e-3
cross-process GPU floor. This is booked as a measured CPU-vs-GPU
sensitivity of the float32-cast ladder, not physics: V247's T0 band is
now backend-aware (1e-3 on cupy, 3e-2 labelled on numpy fallback) and
the verdict JSON carries the backend, so the gate's OFF-CORE line
surfaces any such run. Standalone GPU re-run: 2.9e-04, inside the GPU
band unchanged.

# The pillar-surgery programme: pre-registrations V308-V314

Reader instruction (2026-10-30): the defects are in the pillars themselves,
not in the details. Four pillar choices go under bounded surgery: the
clock/closure (V308-V310), the root signature (V311-V312), the lightcone
dictionary (V313), the scale-formula ledger (V314). Success = a
zero-fitted-parameter formula inside the booked bands (provenance
*algebraic*); honest bounded no-gos are equally bookable outcomes.

**Programme rules (fixed before any run).**
(i) Nothing may be refitted to land in a band.
(ii) Every constant is booked (a `cl10/config.yaml` path or a booked
artifact via the `art()` pattern) - a new dial reclassifies the result to
fitted.
(iii) The gate re-runs incrementally only (fingerprint gate; never
`--full`). The bands live in the `pillar_surgery:` config section;
changing a band mid-run is a reclassification to fitted and must be said
in the booking.

## V308_coefficient_variants_fifth_clock (tier 3, O8)
**The question (pre-registered).** the rung-10 clock 1/5 reads the
coefficient at 2.088106865906787 (V206 `coef_at_fifth`) vs the observed
2.0765641 - a gap of -0.55%. Does any *variant reading of the loop
machinery at that clock* close it?
**Object.** import from
`experiments/scales_loopfactor/V203_closed_branch_coefficient/V203_closed_branch_coefficient.py`:
`hop_machinery`, `loop_ratios` (both clock-parameterized:
`hop_machinery(g10)`). Gate T0: reproduce V206's `data_clock` values to
1e-9 by artifact read (no literals). Variants, each scored against
`hierarchy_coefficient_observed` (loader-derived, read from config tree)
with band `coef_band_abs`:
- (a) plain rung-10 read at g = 0.2 (must reproduce coef_at_fifth);
- (b) certified-continuation read at g = 0.2: geometric continuation
  F_n = (rho, F_cross) recipe exactly as V276 T1 (import recipe by
  reading `V276_certified_ladder_208.py` lines 124-240 and copying the
  continuation construction), measured (rho, F_cross) re-derived at 0.2
  inside the experiment;
- (c) icosahedral-module closed-word ladder at g = 0.2 (copy the V218
  construction: read
  `experiments/scales_loopfactor/V218_closed_word_census_module/V218_closed_word_census_module.py`
  and reuse its ladder function; rung-10 coefficient from its F_10).
**Branch (pre-registered).** `COEF-VARIANT-CLOSED` (some variant within
`coef_band_abs`) / `COEF-GAP-PERSISTS` (book the five smallest variant
residuals with their provenance).

## V309_clock_algebraicity (tier 3, O8/O4)
**The question (pre-registered).** is g_obs = 0.19999798835751165
(artifact read V206) an algebraic number of small height, i.e. could a
root rule generate it?
**Measured.** `mpmath.findpoly(g_obs, d)` for every d <=
`pslq_max_degree` (mpmath pattern as V276); plus the rational census:
every a/b, b <= 5000 inside `clock_band_abs` of g_obs (count + closest
listed); plus the booked difference g_cross_rung10 - g_obs (both from
V206 artifact: 0.20013316730041675 - 0.19999798835751165) reported to
full precision as a named gap 1.35e-4 for V310's census.
**Branch (pre-registered).** `CLOCK-ALGEBRAIC` (relation found,
polynomial booked - then the V314 ledger must use it) /
`CLOCK-NO-RELATION-DECLARED`.

## V310_coefficient_expression_census (tier 2, hard, O8)
**The question (pre-registered).** does a bounded expression over BOOKED
constants supply the -0.55% coefficient correction (multiplicative
target: expr = 2.0765641/2.0881069 ratio form AND the direct form
= 2.0765641 +- 1e-4) or the clock gap 1.35e-4?
**Object.** the expression generator pattern from
`experiments/neutrino/V306_coefficient_census_atm_gap/V306_coefficient_census_atm_gap.py`
(`itertools.permutations` over a booked dict, pairs; read its lines
100-130), extended to depth 3 (`expr_depth_max`): pairs, then one more
booked constant with ops {+, -, x, /, pow^+-1}. Base constants (each
read via config path or `art()`, no literals): LOOP
(u10pred.LOOP), `angle.alpha_10`, the three `angle.branches` values,
`physics.m_planck_gev`, `physics.v_electroweak_gev`, loader-derived
`hierarchy_coefficient_observed`, `sm.m_z_gev`, sm.alpha1_mz/alpha2_mz/
alpha3_mz, `sm.m_unification_obs_gev`, `sm.delta_b_obs`,
neutrino.m_atm_ev/m_sol_ev/sum_mnu_max_ev, `cosmology.omega_dm_h2`,
V206 artifact `g_obs/g_cross_rung10/coef_at_fifth`, V276 artifact
`rho/F_cross`, V206 closure-class `coef(pi/4) = 2.7002077805447326` via
art. Report the five smallest residuals per target with their expression
strings; supplied word only at the booked bands.
**Branch (pre-registered).** `COEF-EXPRESSION-SUPPLIED` /
`COEF-CENSUS-NEGATIVE` (bounded no-go booked).

## V311_quaternionic_root_step (tier 3, O5/O8)
**The question (pre-registered).** the real form Cl(10,0) carries no
CP-odd invariant (V117). The quaternionic real forms of the SAME complex
algebra Cl_10(C) - p in `quaternionic_p` - do. Does the automaton's step
survive a signature swap?
**Object.** generators by exact construction: from
`cl10.clifford.generators()` the 10 kron matrices Gamma_j; set
Gamma'_j = Gamma_j for j <= p, Gamma'_j = i Gamma_j for j > p
(anticommutation preserved, squares flip). Plane generators and step:
copy the 4-line pattern of `cl10/sm_qca.py:48-73` (T'_ab = i
Gamma'_a Gamma'_b, `plane_rotation`, `gut_qca_step`) into the experiment
as `signature_step(p, g)`. Measured per p: FS indicator of the module
via import `fs_indicator` from
`experiments/baryogenesis/V117_real_form/V117_real_form.py` (line 83) -
expect -1 (T0 gate); the CP-odd invariant test of V117 T3 pattern (the
J U J^-1 vs U^-1 deviation, computed in-experiment on the step at the
closed clock `angle.branches.closed`); unitarity residual
||U'U'^dagger - I||_max (V208 `unitary_err` pattern, computed
in-experiment); the closure census: for every g = pi a/b with 0 < a < b
<= `closure_den_max`, smallest k <= `closure_k_max` with
||U'(g)^k -/+ I||_max < 1e-12 - book the found closure clocks per
signature. The quaternionic Gamma construction (Gamma_j -> i Gamma_j for
j > p) is a choice of real-form embedding inside the fixed complex rep;
the step's unitarity verdict is with respect to this embedding - stated
here so the no-go, if it fires, is class-bounded and honest.
**Branch (pre-registered).** `STEP-UNITARY` per signature (some p keeps
exact unitarity - record its closure clocks) /
`STEP-HYPERBOLIC-DECLARED` (all quaternionic steps fail unitarity - the
signature surgery dies at the step level; this is a theorem-grade
booking).

## V312_quaternionic_loop_ladder (tier 3, O8; conditional)
**Contingency (pre-decided).** If V311 books `STEP-HYPERBOLIC-DECLARED`
for ALL p, DO NOT run the ladder; the experiment file still exists, its
main() reads V311's artifact, books the skip explicitly
(`branch: "SKIPPED-STEP-HYPERBOLIC"`), passes its single T0 check, and
exits 0. If at least one signature is `STEP-UNITARY`, run for those only:
**Object.** rebuild the V203 machinery over the new root in-experiment:
Ts = the 45 plane generators T'_ab; structure constants by calling
cl10's `structure_constants` (import site per V203's imports) on the new
Ts; adjoint links via the V79 `adjoint_link_operators` (import as V203
does); then the ladder F_n, crossing n*, deficit and coefficient at
V311's closure clocks for that signature.
**The question (pre-registered).** does any quaternionic-root closure
clock read the coefficient inside `coef_band_abs` of the observed, or
reproduce the odd-n exact selection rule with a crossing whose deficit
matches the observed coefficient?
**Branch (pre-registered).** `SIGNATURE-LADDER-MATCH` /
`SIGNATURE-LADDER-NEGATIVE` (closest reading booked).

## V313_lightcone_dictionary_conjunction (tier 3, F8/O10)
**The question (pre-registered).** the exactly lightlike graviton cone
exists but only in a 53-parameter family nothing selects, whose known
repair clock g = 1/2 collapses the loop ladder (paper Section 11.6,
V208/V285/V229). Re-parameterize the family and test the conjunction
member-wise.
**Object.** reuse V285's machinery wholesale - read
`experiments/gravity/V285_most_local_lightlike_generator/V285_most_local_lightlike_generator.py`
(module header lines 1-90: `build_Oh`, `module_reps`, `velocity_ops`,
covariant 160-space basis, Clifford-regularity Jacobian rank 107 /
nullity 53, tangent coordinates delta in R^53, Nelder-Mead locality
functional T with `nelder_starts` starts). Re-measure the
Jacobian/nullity at every clock in `scan_clocks` (0.125, 0.2, 0.5,
pi/4, pi/2): report the family dimension at each. For the clocks with
nonzero family: minimize T; for the minimizer member evaluate (i) ladder
intactness: the V206 closure-deficit at that clock from the V203 ladder
(deficit band `ladder_intact_band` vs the pi/4 value 0.98953), (ii) the
selector tests: locality-minimum uniqueness (V285's T* vs
T_committed), bond-entropy weight G = 1/(4 ln 45) as a weight only
(report, never refit). T0 gate: the Jacobian rank/nullity at pi/4 must
reproduce V285's booked 107/53 before any scan result is trusted; if it
does not, STOP and report the drift as `MACHINERY-DRIFT`.
**Branch (pre-registered).** `LIGHTCONE-SELECTED` (a member at a
ladder-intact clock at which a booked selector functional is uniquely
minimized) / `LIGHTCONE-CLOCK-CONFLICT` (families empty at every
ladder-intact clock, or selectors never unique - bounded no-go with the
per-clock family dimensions booked).

## V314_four_scale_ledger (tier 2, hard, O8/O10)
**The question (pre-registered).** assemble every candidate supplied by
V308-V313 into compound scale-theories and score them against the four
measured scale contacts. The four contacts and their current standing
(all re-derived in-experiment from config/artifacts, gate
`ledger_scale_rel_tol`): seesaw M_P LOOP^3 = 9.6e13 GeV (algebraic, no
band issue), DM M_P LOOP^10 = 118 GeV vs the observed window
`physics.wimp_window_gev`, neutrino ladder M_P LOOP^17 and the candidate
S17/3 vs `neutrino.m_atm_ev` (the -2.7% booked gap), and the hard one:
rung-10 coefficient vs `hierarchy_coefficient_observed` at the best
clock. A compound theory = (clock rule from V309/V311/V313 or the fitted
1/5, coefficient variant from V308 or census hit from V310); enumerate
all combinations (bounded: the variant/candidate lists from the
artifacts, cross product < 50 rows), score every row on all four
contacts, print the full ledger table sorted by worst-contact residual.
**Branch (pre-registered).** `LEDGER-CONSISTENT-THEORY-CANDIDATE` (some
row: all four contacts inside their booked bands AND every input
provenance algebraic/conventional/measured - zero fitted - then this is
the true formula and the roll-up reclassifies the paper) /
`LEDGER-BEST-DEFICIT-DECLARED` (book the closest row's worst residual
with the full table; this bounds the current programme's best compound
theory).

### Post-run booking note (2026-10-30, after V308-V314 ran)

- V308 fired the pre-registered COEF-GAP-PERSISTS. T0 reproduced V206's
  coef_at_fifth 2.088106865906787 to <1e-9. The certified-continuation
  variant is DEGENERATE at g = 1/5 (the crossing at that clock is
  n* = 10, so the V276 continuation is the identity); the icosahedral
  return-amplitude ladder reads 2.5e26 (the box mean of tr M(q)^10 does
  not decay to the loop scale, V218's T4 as measured). Five smallest
  residuals booked; plain read closest at |d| 1.15e-2.
- V309 fired CLOCK-NO-RELATION-DECLARED exactly as the verification
  predicted: findpoly empty at every degree 1..6 (height cap 20), zero
  rationals a/b with b <= 5000 inside 1e-6 of g_obs (closest 1/5 at
  2.012e-6), gap g_cross_rung10 - g_obs = 1.3517894290510091e-4
  booked as V310's target.
- V310 fired COEF-EXPRESSION-SUPPLIED at the clock-gap target only:
  ((alpha_10*g_impl)*alpha1_mz) = 1.352e-4 vs the booked 1.35178943e-4,
  residual 2.106e-8 inside the 1e-6 band. Class discipline deviation
  booked before the word was accepted: the depth-3 op list was scored
  in two layers - CLASS-PRE (the pre-registered {+,-,x,/ both orders}
  plus the pow^+-1 reading p^-1 = 1/p) decides the word; CLASS-POW
  (arbitrary booked exponents) is scored and booked as annotation only
  (44 in-band power coincidences, V306's lesson kept). The tautology
  exclusion (a candidate may not use its target's defining constants)
  is booked in the artifact. Qualifier: alpha_1(M_Z) = 0.0169 is a
  three-digit measured input - the match's 1.6e-5 relative rate is
  inside the published uncertainty of that input, so the expression
  cannot be finer than its coarsest ingredient.
- V311 fired STEP-UNITARY on p = 0 (Cl(0,10)) only (unitarity residual
  1.443290e-15 at the closed clock; closure clocks pi*1/2 k = 2,
  pi*1/4 and pi*3*4 k = 4, deviations 4.9e-16 to 2.6e-15); p = 8, 7,
  4, 3 fail unitarity (residuals 2.0, 2.0, 4.0, 2.0 - the mixed-sign
  planes make the rotation factors singular at pi/4). T0 FS = -1
  re-derived at every booked root. The verdict is class-bounded to the
  booked embedding (Gamma_j -> i Gamma_j for the last 10-p), exactly
  as pre-registered.
- V312 ran (the SKIPPED-STEP-HYPERBOLIC contingency was NOT needed
  because V311 fired a unitary signature). SIGNATURE-LADDER-NEGATIVE:
  the exact odd-n selection rule survives the Cl(0,10) embedding
  (odd_max < 1e-12 at all three closure clocks), no closure clock
  reads the coefficient in band; closest pi*1/2 at 2.2995015
  (|d| 2.23e-01).
- V313 fired LIGHTCONE-CLOCK-CONFLICT. Machinery note: the first run
  misread the ladder leg - on the exact-rule clocks a vanishing odd
  length is not a crossing (it took n = 3 with deficit 5.5e-14 at
  pi/4). Fixed to the V276/V204 scan convention (skip vanishing odds
  where the rule holds) with a pi/4 deficit anchor to the booked V206
  value 0.989526873, on the BACKEND-AWARE band of the V247 precedent
  (the float64 root-clock deficit is a floor artefact, V209: 1e-3 on
  cupy, 3e-2 labelled on the numpy fallback). Re-run: backend cupy,
  anchor |d| 8.79e-05 inside the 1e-3 hard band; intact set {1/5, pi/4,
  pi/2}; 1/2 collapses at 1.44e-3 and 1/8 at 0.4836. The branch word is
  unchanged between the two runs; the first-run artifact was superseded,
  not booked.
- V314 fired LEDGER-BEST-DEFICIT-DECLARED: 11 compound rows; contacts
  1-2 exact (V197 re-derivation gates 1e-14), the neutrino contact
  S17/3 vs m_atm at -2.70% gives every row normalized 0.539; the
  fitted g_obs row passes all four contacts (coefficient contact
  normalized 1.94e-09 - the inversion is the definition) but carries
  fitted provenance; the best ZERO-FITTED row - the V310 corrected
  clock - reads 2.0764436008045477 vs 2.076564111173832: |d|
  1.205e-4 vs band 1e-4, worst normalized 1.2051. The programme's
  bound on its best compound theory, booked with the full table.
- Contingencies used: V312 skip NOT used; MACHINERY-DRIFT stops NOT
  triggered (the V313 first-run defect was corrected machinery, the
  T0 gates themselves passed); no band widened; the pillar_surgery
  config section was never edited after the pre-registration commit.
- Paper-editing note: the plan's conditional text assumed all seven
  words negative. Two fired locally positive (V310 gap-supplied, V311
  p = 0 unitary) while the arbiter V314 declared the ledger deficit,
  so the Table-9 row and the Section 17 sentence were kept as
  bounded no-gos with the two surviving facts named inside them, and
  the O8/O10 rows carry the no-go words with the ledger deficit
  1.205e-4 as the tightest booked residual. No published claim was
  re-provenanced algebraic: the ledger rejected every zero-fitted
  compound theory, so no claim met the plan's success condition
  (LEDGER-CONSISTENT-THEORY-CANDIDATE with zero fitted inputs).

## Gate-hygiene note (2026-10-30): the pillar-surgery full re-run and five repairs

The new pillar_surgery config section changes the library fingerprint, so
the gate re-ran all 227 (GPU 1, --jobs 3, cupy): 222 PASS, 5 FAILED, none
a physics regression.

- V311/V312 crashed on the cupy core: swapped_generators did
  np.asarray(engine array) - implicit cupy-to-numpy conversion is illegal
  on the routed engine. Fixed with the documented host transfer
  cl10.engine.asnumpy. The re-run on cupy is BYTE-IDENTICAL to the booked
  artifacts (p=0 unitarity 1.443290e-15; crossings 2.2995015/2.7040813/
  2.7159955): the original standalone bookings had silently run on the
  numpy fallback. Physics unchanged on both cores; the fix is the honest
  backend label.
- V272/V274 crashed cupy OOM (4.3 GB allocation) under the jobs-3
  concurrency. Solo reruns on a quiet GPU: V272 CHAIN-ALIVE (kernel dim
  12, gap 3.600e+01), V274 NEW-CONJUNCTION-DEAD - values as booked.
- V278 reported MACHINERY-DRIFT as an ORDERING artefact of the single
  pass: it ran at t+6824 s comparing its host recomputation (1.998504)
  against the still-master-era V277 artifact (1.998493); V277's own
  rerun - later in the same pass - reads 1.9985044777. The 1.1e-5 shift
  on the top-64 degree fit is environment-level (both cores move
  together). Solo rerun against the settled V277: T0 inside 1e-6, full
  scan PASS, CONE-RANGE2-ABSENT (gap 5.55e-17).
- The fingerprint ledger was seeded (--seed-fingerprints) from the fresh
  green verdicts of this session - every one of the five was re-run on
  GPU after its repair - and the default gate then reports
  "0 re-run, 227 folded, library unchanged": 227/227 PASS, hard tier MET
  over 87, 51/51 GPU-backend reports, no numpy-fallback label anywhere.

V313 ledger-state note: the gate re-run supersedes the bg_13 rerun as
the booked artifact. Deficits move inside the documented bands (1/2 at
0.0014400, pi/4 anchor |d| 4.3e-05 cupy); the intact set {1/5, pi/4,
pi/2} and both selector verdicts (interior False, unique False) are
invariant; T* moves 1.195211e-3 -> 1.195142e-3 - boundary-kink
sensitivity at the fifth significant figure, booked, no verdict effect.

# The theory-hunt programme: pre-registrations V315-V318 (branch theory-hunt, 2026-10-30)

Programme rules (fixed for all four experiments): (i) nothing is refitted to
land in a band; (ii) every input carries a provenance tag - algebraic /
config-declared / measured / combinatorial / fitted; a dial tuned to a target
is `fitted` and disqualifies only the zero-fitted success tier, never the
ranking; (iii) the classes below are closed sets fixed before the run - any
class widening mid-run is a reclassification and must be said in the booking;
(iv) the config is untouched this programme, so the gate re-runs
incrementally only (library fingerprint unchanged). No config.yaml edit; all
bands come from existing keys: pillar_surgery.{coef_band_abs, clock_band_abs,
ledger_scale_rel_tol, closure_den_max, closure_k_max, scan_clocks,
nelder_starts}, tolerances.strict, angle.branches.closed, neutrino.m_atm_ev,
physics.wimp_window_gev, angle.alpha_10, desert.{alpha_x, m_x_closed_gev},
continuum.{radii}.

Implementer-pinned artifact keys (the plan's "read the artifact first and
pick the key, both named"): V315 closure anchor
`V311_quaternionic_root_step data.signatures."0".closure_clocks` (the p=0
real-form row: pi/2 k=2, pi/4 k=4, 3pi/4 k=4, 32-dim signature step). The
V206_clock_tension artifact carries NO step-closure-table key (its
`data.closure_class` is the coefficient-crossing class, and
`data.data_clock` the inversion data) - so the V311 row is the booked clock
anchor and the 128-dim six-stage step's pi/4 closure is booked against
V203_closed_branch_coefficient T2's fact M(0)^2 = I at pi/4 (k = 2, dev
< CLOSE_EPS). NOTE booked before the run: the plan's token "{pi/4: k = 8}"
matches no booked artifact (V311 books k=4 on the 32-dim signature step,
V203 M(0)^2 = I on the 128-dim step); the booking uses the artifacts.
V316 W1 eta*: `V202_loop_angle_selection data.selection.nearest_crossing`
(the crossing selected against g_10 = 1/5 under the booked ANGLE-OUTSIDE
verdict). V316 rung-17 anchor: `V314_four_scale_ledger
data.rederived.S17_over_3_eV`. V318 T0: `V314_four_scale_ledger
data.best_zero_fitted_row.worst_normalized` (the current gated leaf
1.2051036935645243 - the plan's literal 1.2051036928450998 is the superseded
pre-gate rerun artifact per the booked gate-hygiene note; the artifact wins).
V318 C5: `V237_desert_threshold_matching data.matching.shift_over_spacing_12`
(booked 1.1171513611688386), `data.matching.dz12_efolds` (booked
13.097478575754906) and `data.matching.alpha_X`; the pass criterion is
V247_two_dial_desert_contacts' booked T0 desert-reproduction rule
|dz12_efolds(g) - booked| < 1e-9 (host float64 rule, zero matrix noise).

## V315_closure_quantizer_census (tier 2, hard; O4/O8; target A)

**Question.** Is there a booked-class deformation of the six-stage step whose
closure class contains a RATIONAL clock pi*a/b whose order-k closure coexists
with the exact odd-n selection rule and whose rung-10 read lands on the
observed hierarchy coefficient inside pillar_surgery.coef_band_abs? 1/5 =
pi*(1/(5 pi)) is not rational in pi - the search is for the closest rational
closure AND the coefficient contact at it, not for 1/5; the distance
|b - 10| is booked as the clock-order distance.

**Class (closed, fixed pre-run).** Two families over the committed 6-stage
product, machinery from cl10/supercell.py (STAGES, stage_matrix(even, j, q,
g10, bond_g10)):
- D-family (decorations): per-stage internal twist on/off at the booked
  amplitude (alpha_10, alpha_10, alpha_10), alpha_10 = angle.alpha_10 - 2^6 =
  64 assignments. dec_step(q, g, amps) is V277_bond_decoration_cone_scan's
  bloch_dec copied with the fixed GCLOSED replaced by the scan clock g in
  both stage_matrix slots; twist(q, amp) copied from V277 lines 94-98
  (exact 16x16 eigh route). BOOKED BEFORE THE RUN: V277 books W(0) = I
  (inherited), so every D-variant has the identity step at q = 0 - the
  closure census of the D-family is degenerate with the identity variant by
  construction, and the rung-10 ladder leg (V308's traced_ladder body at the
  clock) carries no stage-wise U16 link at all: its link layer is the 32-dim
  GUT hop and the adjoint link operators, so the plan's U16 -> U16 x twist
  substitution is inapplicable to it and the ladder leg IS the plain traced
  ladder at the variant clock. This is the class, stated pre-run, not a
  post-run reclassification.
- P-family (stage orderings): all 720 permutations of SC.STAGES fed to the
  same copy with zero decorations (the product is non-commutative;
  SC.bond_product's own stages= parameter is the engine's blessing of this
  class).
Total: 64 + 720 = 784 variants.

**Measurement per variant.** Closure census exactly in V311's pattern
(closure_census): scan g = pi*a/b, 0 < a < b <= pillar_surgery.closure_den_max
(79 distinct), smallest k <= pillar_surgery.closure_k_max with
||U(0)^k -/+ I||_max < CLOSE_EPS (CLOSE_EPS = 1e-12, V311's own definition
line copied verbatim) on U = dec_step(q = 0, g, amps). For every variant with
>= 1 closure clock: ladder leg at that clock, traced_ladder_dec(g, amps,
nmax=12) = V308's traced_ladder body at g (see the booked D-family note);
odd-rule check odd_max = max |F_n|, n odd <= 11, against
tolerances.strict (V312's test); rung-10 coefficient via V308's imported
read_coef(F10, LOOP) = 3.0*(F10/LOOP)**10; contact against V308's COEF_OBS
formula read (u10pred.V_EW/u10pred.M_PLANCK)/LOOP**10.
**T0 gates.** (1) identity decoration + default order: the closure clock SET
{pi/2, pi/4, 3pi/4} closes with dev < CLOSE_EPS (anchor: V311 p=0 artifact
row keys verbatim; the booked 128-dim fact M(0)^2 = I at pi/4 from V203 T2
sets the expected k = 2 there - the k difference vs V311's 32-dim k=4 is
booked, not hidden); (2) V308's coef_at_fifth reproduced to < 1e-9 (V308's
own T0 pattern, via art). On failure: MACHINERY-DRIFT, nothing interpreted.
**Branch (pre-registered words).** QUANTIZER-CLOSURE-FOUND (a variant +
closure clock with odd-n exact AND |coef - COEF_OBS| <= coef_band_abs) /
QUANTIZER-NEGATIVE (book the full table: every variant with any closure
clock, its clocks, order k, odd_max, rung-10 coef, |coef - COEF_OBS|;
closest row named; variants with no closure clock counted per family).
Contingency pre-registered: if NO variant had any rational closure, the
search widens only by the engine's second dial bond_g10 != g10 over the
booked closed set {pi/4, pi/2} - same 64+720 census, booked as an ADDITION.
(The identity variant is a class member with a booked pi/4 closure, so this
contingency is unreachable unless T0 itself drifts.)

## V316_weight_double_contact (tier 2, hard; O8 + neutrino; target B)

**Question.** Does a booked-class reweighting of the traced-word sum move
BOTH rung contacts simultaneously inside their bands - rung-10:
|coef(g=1/5) - COEF_OBS| <= coef_band_abs; rung-17: |S17(variant)/3 -
neutrino.m_atm_ev|/m_atm <= ledger_scale_rel_tol - WITHOUT touching the
mechanism or the clock (g = 0.2 fixed)?

**Object (closed class, mechanism and clock g = angle.branches.implemented
fixed).** Three weight families:
- W1 eta-law: V200's F10_of(eta) (its line 157 body lifted verbatim,
  NMAX = 10 as V200 books) at the booked eta set {angle.branches.closed
  members = {pi/4}, V202's selected eta* =
  data.selection.nearest_crossing} - 2 values.
- W2 phase-class renormalization: the traced words reweighted by V202's
  phase_classes(B0, cartan, half_ev, n) class index k (the even sum of
  half-eigenvalues); booked multiplier class {1, cos(eta_class),
  cos(2 eta_class)} with eta_class = k * g at the clock g = 1/5 - three
  variants. F_n^m = sum_j |sum_k m_k e^{i k g} C_k^{(n)}[j]| / Gnorm_n with
  the class vectors C_k^{(n)} built by V202's prefix construction per length
  n (n = 1..20; the rung-10 leg reads n = 10). Consistency gate booked with
  the class: the m = 1 member reproduces loop_ratios' F_n at the clock to
  < 1e-9 for all n (V202's "the phase classes ARE the census" equality).
- W3 node-set channel factors: the 12 M-channel node-set factors at g = 0.2,
  computed by V205's machinery pattern - its block-3 k=0 propagator-trace
  filter lifted verbatim (w = SC.phases(SC.bloch(q, U16(g), g10 = g)) at
  each q in SC.M_POINTS, g = 1/5; f_s = |sum_j exp(i Re w_j)|/128); all
  C(12,3) = 220 triples. The coefficient is read with the triple's combined
  factor per V205's own composition (V203/V204/V205's booked
  required-sum-factor relation): coef(T) = (sum_{s in T} f_s/f_bar) x
  (F_10(1/5)/LOOP)**10, f_bar = the 12-factor mean at the clock - the unique
  normalization that reproduces the booked sum-factor-3 plain reading when
  the factors are equal, stated pre-run (the plan's contingency anticipated
  the alternative multiplicative form; the sum-factor form is the one
  V205's main() books and is used verbatim). Node-set choice is
  combinatorial provenance, never fitted.
Rung-17 leg per variant: S17^var = nv.loop_ladder(17) * nv.EV_PER_GEV x
(|F_17^var| / |F_17^plain|) - V314's re-derivation block recompute with the
variant's F-ratios at rung 17 (the plain member has ratio exactly 1 and
reproduces the booked leaf; the ladder runs nmax = 20 - F_n for n <= 12 is
bit-identical, loop_ratios is incremental). W3 does not touch the F
machinery: its ratio is 1.
**T0.** The plain read reproduces V206 coef_at_fifth (leaf 2.088106865906787,
via art, < 1e-9) AND V314's booked S17/3 (leaf 0.048749228..., via art,
< 1e-9). The W3 leg runs exactly 220 triples (count booked in the artifact).
**Branch.** WEIGHT-DOUBLE-CONTACT (a variant passes both bands; its full
reweight rule booked) / WEIGHT-NO-DOUBLE (book the full table ranked by the
max of the two normalized residuals, top 20 rows + per-family best).

## V317_selection_conjunction (tier 3; F8/O9/O4; target C)

**Question.** Does ONE booked selector functional uniquely select at >= 2 of
the three unselected sites - S1: a member of the 53-parameter lightcone
family at the clock 1/5; S2: the six desert couplings kappa_s (V105/V122/
V145/V146 - cross-site consistency, not re-supply); S3: the 3-of-12 channel
triple (the 12 chain-kernel channel indices in V272's booked
chain_H_spectrum order; the triples are the object).

**Measurement.** Functional x site grid, 4 x 3 = 12 cells; uniqueness per
cell fixed pre-run; cells not listed below are booked `N/A` pre-run, never
invented later:
- F1 locality T (V285 machinery via V313's rebuilt block; T_committed =
  sum(profile(H)^2, r = 2..7) as V313 books): S1 EVALUATED at V313's
  artifact members - data.search.x_best (the T* minimizer) and the
  committed member x = 0; V313's top5 = [0.0011951..., 1000 x 4] books only
  ONE non-penalty member - booked here so the member count is not invented
  later. Unique iff T(member) < T_committed strictly, interior by V313's own
  rule (53-grad norm < 1e-6 x max(1, T_committed) and no box-face
  coordinate), and one distinct minimum. S2: N/A. S3: N/A.
- F2 bond-entropy weight G = 1/(4 ln 45) (V313's weight block): S1
  weight-only, never a selector (V313's own wording booked); S2 N/A; S3
  EVALUATED: G-weighted triple score G x sum_{s in T} E_s over all 220
  triples (E = V272 chain_H_spectrum), unique iff the gap to the runner-up
  triple > 1e-6 (note booked: a uniform positive weight cannot reorder the
  triples - the F2/F3 S3 cells share their uniqueness gap by construction).
- F3 host energetics: S1 N/A; S2 EVALUATED via V247's v237_rule(AX) +
  meet_solve(AX) (its lines 113-127) over the candidate kappa_s families
  {V105_cartan_kappas, V122_class_invariant_kappas, V145_superexchange_kappas,
  V146_two_number_kappas}: family F CLOSES the meeting iff the (1,2)-shift
  solution from meet_solve at the V237-booked alpha_X
  (data.matching.alpha_X) is consistent with the family's own booked
  desert-closure flag: V105 models.*.<n>.desert_ok any true;
  V122 embeddings.<e>.free_kappa.rows any closed==true with its booked
  [z_min, z_max] inside the desert e-fold window (0 < z < ln(desert.
  m_x_closed_gev / M_Z)); V145 meeting rows thresholds_in_desert any true;
  V146 joint-row desert flag (key named in the check string after reading
  the artifact - the flag exists in data.joint). Unique iff exactly one
  family closes. S3 EVALUATED: V302_mD_host_forcing's energetics lifted
  verbatim (ones_shift = sum of the chosen channels' chain_H_spectrum
  entries, V302 T1 convention); unique iff gap to runner-up > 1e-6.
- F4 relevance exponent (V213's p_A pattern): S1 N/A (V213 measures
  anisotropy at a fixed member, not members); S2 N/A; S3 EVALUATED: p_A of
  the triple's splitting = mean over the triple's nodes of V213's
  power_ladder slope of A_transverse over continuum.radii (V213's
  tangent_field/lowest machinery verbatim, its M-node convention generalized
  to the triple's M nodes) at clock 1/5; unique iff gap to runner-up > 1e-6.
**T0.** V313's T_committed = 0.0011956422840809166 reproduced < 1e-9 from
the same machinery; V247's booked (1,2)-shift row reproduced < 1e-9 from
v237_rule at alpha_std(pi/4) (dz12 leaf 13.097478575754906).
**Branch.** SELECTION-UNIFIED (one functional unique-selects at >= 3 sites)
/ SELECTION-PARTIAL (exactly 2) / SELECTION-DISPERSED (every functional at
<= 1); the full 12-cell table booked in all cases with the selecting
member/triple/family values when a cell is unique.

## V318_closest_theory_ledger (tier 1, hard; the arbiter and the deliverable)

**Question.** Among ALL compound theories = (clock from the booked set) x
(weight variant) x (selection), ranked on FIVE contacts, which comes closest
to the measured values, and does any zero-fitted row pass them all?

**Rows.** Clocks: 1/5 (fitted-inverted reference), the V309+V310 corrected
clock g_cross_rung10 - alpha_10*(1/5)*alpha1_mz (V314's row-(d) construction
copied verbatim from the V310 supplied-hit value), every V315 closure clock
(with its odd-rule status), V311 p=0 clocks pi/2, pi/4, 3pi/4 (deduplicated
against the V315 set). Weights: plain, V316 double-CONTACT winners if any
else top-3 by the double ranking, W1 eta values. Selection: none, V317 S3
unique triple if a cell fires, V317 S1 member if a cell fires. Cross product
cap 200: if over, keep top-3 per dimension and book the pruning rule. Every
weight functional is read at the ROW's clock (V314's per-row clock
convention, copied unchanged); rung-17 legs use the variant F-ratio at the
row clock as in V316.
**Contacts** (V314's scoring block copied verbatim; normalized-residual
convention unchanged): C1 seesaw M_P*LOOP^3 (algebraic, normalized 0);
C2 DM window (physics.wimp_window_gev, band ledger_scale_rel_tol); C3
rung-17 neutrino (same band); C4 rung-10 coefficient (coef_band_abs); C5 NEW
desert-coupling contact: V247's booked criterion - residual
|dz12_efolds(alpha_std(clock)) - V237-booked| , band 1e-9 (V247's T0 desert
reproduction rule), the row's shift_over_spacing_12 booked against the
booked 1.1171513611688386.
**Provenance per input:** clock/weight/selection tagged {algebraic, config,
measured, combinatorial, fitted}; zero-fitted = no fitted input anywhere in
the row.
**T0.** V314's best_zero_fitted_row.worst_normalized reproduced < 1e-9 by
re-running V314's 11 rows inside V318 (then extended; the artifact leaf wins
over the plan's superseded literal, per the pinned-key note above).
**Branch.** LEDGER-CONSISTENT-THEORY-CANDIDATE (a zero-fitted row passes all
five contacts) / LEDGER-CLOSEST-THEORY-DECLARED (no such row: book the FULL
ranked table and name THE closest row explicitly - its five normalized
residuals, its provenance, and its distance = worst normalized; this row is
the programme's closest-theory deliverable regardless of provenance; closest
is not true and the booking language will not blur that).

### Post-run booking note (2026-10-30): the theory-hunt programme V315-V318

Fired words (all pre-registered verbatim): V315 `QUANTIZER-NEGATIVE`,
V316 `WEIGHT-NO-DOUBLE`, V317 `SELECTION-DISPERSED`, V318
`LEDGER-CLOSEST-THEORY-DECLARED`. No zero-fitted candidate passed the five
contacts, so no paper claim was re-provenanced algebraic; the programme's
deliverable is the named closest row with its distance booked below.

V315: all 784 booked-class variants (64 D-decorations - identity-degenerate
at q = 0 by V277's booked W(0) = I - and 720 P stage-orderings) have rational
closures: 43 of the 79 pi*a/b clocks (b <= 16) close the 128-dim six-stage
step at Gamma within CLOSE_EPS (max dev 5.1e-14); the identity row closes
pi/2 (k = 2), pi/4 (k = 2 - V203's M(0)^2 = I, dev 7.5e-15) and 3pi/4
(k = 2), and the P-family adds no new closure clock. The odd-n rule is exact
only at the order-8 clocks (pi/4, 3pi/4), where the rung-10 read at nmax = 12
is the V203-booked 55.03 -> 3(F/LOOP)^10 = 7.86e+34 against COEF_OBS 2.0766
(band 1e-4): the closest odd-exact closure misses by 35 orders, clock-order
distance |b - 10| = 6. The second-dial contingency (bond_g10) did not fire.
Contingency used as pinned: V206 carries no step-closure-table key, so the
V311 p=0 row anchored the clock SET and V203 booked k = 2 at pi/4; the plan's
token "{pi/4: k = 8}" matches no booked artifact and was not used.

V316: 225 variants, zero double-band hits. W3: the 12 M-node V205
propagator-trace factors at the clock are EXACTLY equal (f_bar = 0.0838951,
max spread 1.3e-15) - all C(12,3) = 220 triples reproduce the booked
sum-factor-3 reading 2.0881 (channel-blindness, combinatorial restatement);
the contingency anticipated a multiplicative alternative - the sum-factor
form is what V205 books, applied verbatim with the f_bar normalization that
fixes the plain member. W1: eta* = 0.1981004295 (V202
data.selection.nearest_crossing) reads the rung-10 coefficient exactly 3.000
(+0.92 off, 9234x band) but moves the rung-17 leg to S17/3 = 0.051847 eV
(+3.49%, INSIDE the 5% ledger band) - the only single-rung pass booked in
the programme; eta = pi/4 reads 2.44e10. W2 multipliers cos(k g)/cos(2 k g)
read 2.19e8 / 0.003281 (rung-17 rel 710.7x / 1.101x - both out). Machinery
note (booked, not a class change): the first run reported MACHINERY-DRIFT
because the class index int() truncated the noisy half-eigenvalues; V202's
own even-eigenvalue guard with int(round()) fixed it - the class-polynomial
route then reproduces loop_ratios to 2.6e-14 (n = 1..20).

V317: SELECTION-DISPERSED - sites per functional {F1: 0, F2: 0, F3: 1,
F4: 0}; the 12-cell table booked with the 7 pre-registered N/A cells. T0s
bit-exact (T_committed |d| = 0.0; dz12 |d| = 0.0). F1 x S1 not unique: the
T* member (V313 x_best) reads T = 0.00119559 strictly below T_committed but
the 53-gradient is 1.03 (V313's boundary-kink pattern), interior False.
F3 x S2 UNIQUE: exactly one kappa_s family closes the (1,2)-meeting
cross-site - V105_cartan_kappas, model B_general_cartan (rows n = 3..10
book desert_ok True, meeting residual ~1e-14, z0 22.0-32.7 inside the
desert e-fold window (0, 34.842), six kappa_s supplied) - while V122's
free-kappa windows (z_max 32.92) and V145's threshold windows miss the
meet_solve solved positions at the V237-booked alpha_X (zA 6.652, zB 36.213
e-folds - zB above M_X) and V146 books n_desert = 0. S3 never unique under
any functional: the chain-H triple-energy gap 5.86e-15 and the V213 p_A
triple gap 1.42e-11 - the 12 channels sit in three isospectral p_A classes
(0.69830 / 0.70134 / 0.70484), so no 3-of-12 triple separates at 1e-6.

V318: 48 compound rows; the raw cross product 45 closure clocks x 6 weights
x 1 selection = 270 exceeded the 200 cap, so the pre-registered pruning
fired as booked: mandatory clocks {1/5, corrected, pi/2, pi/4, 3pi/4} +
top-3 V315 closure clocks by plain contact (pi*15/16, pi*1/16, pi*4/7);
weight dimension kept its booked composition; the selection dimension is
{none} (no V317 selection-eligible cell fired). T0b: V314's 11 rows re-run
inside reproduce the gated leaf 1.2051036935645243 to |d| 1.4e-10; the
plan's literal 1.2051036928450998 is the superseded pre-gate rerun artifact
(pinned in the pre-reg; the artifact wins). 42 of 48 rows are zero-fitted.
C5 passes 48/48: V237's dz12 = 13.0974785757549 is reproduced at every
row clock to <= 1.35e-13 under alpha_std(g) - the desert (1,2)-shift does
not discriminate clocks, booked rather than hidden.

THE closest row (named explicitly, distance booked):

| input | value | provenance |
|---|---|---|
| clock | g_cross_rung10 - ((alpha_10 g_impl) alpha1_mz) = 0.199997967300417 | algebraic/conventional/measured (V310 census hit) |
| weight | W3 node triple [2, 10, 11], combined sum factor 3.000000 | combinatorial (node set) - degenerate with the plain read to 1e-14 (the equal node factors) |
| selection | none | trivial |
| C1 seesaw | normalized 0.000e+00 | algebraic |
| C2 DM window | normalized 0.000e+00 (118.5709 GeV in window) | config |
| C3 rung-17 | 0.048749228 eV vs 0.0501 eV, normalized 0.539230 | algebraic |
| C4 rung-10 coef | 2.0764436 vs 2.0765641, normalized 1.205104 (BINDING) | algebraic machinery, corrected clock |
| C5 desert shift | |d| 1.35e-13 vs booked 13.0974785757549, normalized 1.35e-04 | measured (V237 booking) |

Distance = worst normalized = 1.205103693702192; zero-fitted: yes. Closest
is not true: the row inherits its clock from the V310 census hit on the
V206-inverted data clock (measured/conventional inputs), and the ranking
does not promote any claim to algebraic. Gate hygiene: the plan's expected
"223 folded" is a miscount - the runner reports the actual fold count
(4 re-run, 227 folded, library unchanged); no config edit was made, so the
library fingerprint is unchanged and the incremental gate re-ran only the
four new experiments.

# The rung-weight programme: pre-registration V319 (2026-10-30)

**Context (from booked artifacts, no new inputs).** V318 named the closest
zero-fitted row at worst-normalized 1.205103693702192; inside that row the
ONLY failing contact is C4 (coefficient 2.076443600804462 at the V310
corrected clock vs COEF_OBS 2.076564111173832, band 1e-4); C3 passes at
0.539230 (the neutrino -2.6962% floor), C5 at 1.35e-04. The required rung-10
multiplier window is exact: W(10)^10 in
[(COEF_OBS - coef_band_abs)/C_CORR, (COEF_OBS + coef_band_abs)/C_CORR].

**Question (pre-registered).** Does a rung-dependent weight built ONLY of
named constants (LOOP from u10pred - the booked ladder base, V276;
angle.alpha_10 - config conventional per the V52 audit), no dial tuned to
any target, move the rung-10 contact into its band AT THE ZERO-FITTED
CORRECTED CLOCK while leaving C3 inside its 5% band - and is the census
then distinguishable from density?

**Class (closed, fixed now, BEFORE the run; V306 rule: the target-defining
constants COEF_OBS, C_15, C_CORR, S17, m_atm, g_cand NEVER appear inside
any W):**
- A loop-order shifts: W(n) = 1 +/- LOOP^k, k = 1..6. 12 members.
- B rung-linear drift: W(n) = exp(+/- LOOP^k * n/17), k = 1..6 (anchor rung
  17 = the ledger's rung B, config cell.rung_B; no free parameter). 12.
- C coupling-mixed shifts: W(n) = 1 +/- alpha_10^a * LOOP^b, a in {1,2},
  b in {1,2}. 8.
- D fitted reference (EXCLUDED from the success tier, ranking only):
  W(n) = exp(lambda*(n-10)) with lambda solved to land rung-10 exactly at
  COEF_OBS; provenance fitted.
32 parameter-free members total.

**Reading conventions (booked, inherited).** coef' = 3(|F_10| W(10)/LOOP)^10
= C * W(10)^10; rung-17 leg' = S17 * W(17) (odd rungs have F = 0 exactly by
the certified odd rule, so the ledger leg carries ratio 1 and the weight
multiplies the leg directly - V316's rung-17 convention); C1 stays
algebraic-normalized (rung-A identity, no measured target; booked); C2 =
DM rung-10 mass times W(10) vs physics.wimp_window_gev; C5 = V247
v237_rule at alpha_std(row clock), weight-free.

**T0 (pre-registered).** C_15 = plain read at g = 1/5 reproduces V206
data.data_clock.coef_at_fifth = 2.088106865906787 to < 1e-9; C_CORR = plain
read at g_cand = g_cross_rung10 - V310-hit reproduces the V318 closest-row
coef leaf 2.076443600804462 to < 1e-9; S17 reproduces the V303 contact leaf
rel < ledger_scale_rel_tol.

**Coincidence discipline (V306's lesson, quantitative, fixed now).** For
each family the artifact books raw hits/total AND the log-uniform expected
random hits = n_candidates x window-log-width / family-log-range. A member
is SUPPLIER-GRADE only if its family's expected random hits < 0.5;
otherwise every in-band hit is booked as a COINCIDENCE (the V306 lesson:
a rich enough in-band class supplies its target by density) and does not
fire the success tier. Consistency leg booked (not branch-deciding): for
every in-band hit and the D reference, also record the certified-branch
crossing residual |F_208 W(208)/(1/16pi) - 1| against the plain booked
1.013% deficit (V209/V276) - a real supplier should not widen it silently.

**Branch (pre-registered words).**
- LEDGER-CONSISTENT-THEORY-CANDIDATE iff >= 1 parameter-free member is
  supplier-grade (density gate above) AND its corrected-clock row passes
  all five contacts zero-fitted; the full five-contact row booked.
- RUNG-WEIGHT-NO-SUPPLIER otherwise: all 32 members ranked by
  max(norm_rung10, norm_rung17), every in-band hit booked with its family
  density and the exact required window; the candidate formulas are booked
  explicitly, with the density verdict stated in the same breath.
Closest is not true: a census hit that lands in a 1e-4-wide band from a
family with O(1) expected random hits is booked as supplied-by-density,
never as derived. No paper claim re-provenances algebraic unless the
density gate fires.

### Post-run booking note (2026-10-30): the rung-weight census V319

Fired word: `RUNG-WEIGHT-NO-SUPPLIER`. T0 exact: V206 coef_at_fifth
reproduced |d| 2.4e-13, the V318 closest-row coef leaf |d| 1.0e-14, the
scale legs 0.0. The required window (corrected clock) is lnW(10) in
[9.878e-07, 1.062e-05]; two parameter-free members land in both bands:

| member | W(10) | coef | C4 norm | C3 norm | certified F_208 deficit at pi/4 |
|---|---|---|---|---|---|
| A `1+LOOP^3` | 1.000007874 | 2.0766071 | 0.430 | 0.539 | 1.055e-2 (plain 1.056e-2) |
| B `exp(+LOOP^3*n/17)` | 1.000004632 | 2.0765398 | 0.243 | 0.539 | 1.047e-2 |

Both would pass the full five-contact ledger at the corrected clock
(zero-fitted provenance, C2 in window, C5 unchanged, distance = C3 =
0.539 - the neutrino floor binds, not the coefficient). They do NOT fire
the success tier: the pre-registered chance-density rule (supplier-grade
iff family expected log-uniform random hits < 0.5) gives E = 1.45 (A),
1.46 (B), 2.66 (C) - the LOOP^3 scale sits inside the decade the window
samples, so the hits are indistinguishable from density and are booked
SUPPLIED-BY-DENSITY, exactly the V306 lesson one level down. The fitted
reference reaches worst-normalized 0.539038: the ceiling of the whole
five-contact ledger is the neutrino -2.6962% contact, not the coefficient.

Booked corrections and facts (not class widenings): (i) the pre-registered
D anchor form exp(lam*(n-10)) is degenerate at rung 10 (W(10) = 1
identically); corrected to exp(lam*n) with lam = ln(COEF_OBS/C_CORR)/100 =
5.803522543772898e-07 - tier-irrelevant, D is excluded by rule. (ii) The
certified F_208 leg is evaluated at the certified ROOT clock (pi/4,
V209/V276); at g_cand the 208-rung ratio underflows float64
(1.4e-23) - booked as fact. (iii) Neither in-band weight closes the
certified crossing deficit (1.056e-2 -> 1.047e-2 at best): the
LOOP^3-scale weight and the certified n* = 208 deficit are separate
residues; no class member touches both.

Meaning for the hunt: a zero-fitted five-contact PASSING theory EXISTS in
the booked formula space - the corrected clock with a LOOP^3-scale rung
weight, distance 0.539 - but the census cannot separate it from
coincidence at band 1e-4. What would break the tie, booked: an
independent consistency the formula must pass by derivation, the named
candidates being the certified n* = 208 deficit (fails), the -2.6962%
neutrino floor (untouched by any rung weight at rung 10's decade), and
the rung-3 seesaw row (weight-free by convention - a derived weight would
have to say what it does there).

# The rung-sector programme: pre-registration V320 (2026-10-30)

**Question (pre-registered).** The V319 census found the required rung
weight (1+LOOP^3-scale) but booked it supplied-by-density. Is that weight
a REAL contribution when the booked twisted sectors of V315 are SUMMED -
treated as internal degrees of freedom, uniform measure, no tunable
parameter? Predict W_eff(10) = (sector-mean rung factor)/(plain rung
factor) and compare to the V319 booked window.

**Object (closed class, fixed now).** The 64 D-family decorations of V315
(amps[t] in {(0,0,0), (G10,G10,G10)} per stage, G10 = angle.alpha_10)
over the engine step dec_step(q, g, amps, SC.STAGES) (V315's copied
machinery). Rung leg per sector via the Floquet phase sum on the BOOKED
12-point M mesh (SC.M_POINTS, the same mesh V205/V318-W3 use):
rho_n^{(a)}(g) = (1/12) sum_q (1/128) tr M_a(q;g)^n, exact via eigen-
values. W_eff(n) = mean_a rho_n^{(a)} / rho_n^{(0)}; secondary measure
booked: average over the CLOSING sectors only (closure_row(M_a(0,
pi/4)) closes within CLOSE_EPS - V315's row). No new q-mesh, no weights,
no amplitudes: everything booked.

**Guard (fixed now, not an after-run escape).** The ratio is meaningful
only when the plain denominator is not phase-cancelled: a branch requires
|rho_n^{(0)}| > 1e-9 for the n compared; if the guard fails at some n the
row is booked GUARD-FAILED at that n and cannot contribute at it.

**Targets (artifact-read, no literals).** The V319 window leaves: at the
corrected clock lnW(10) in [9.877594588304732e-07, 1.0619053724142347e-05];
at 1/5 lnW(10) in [-5.591e-4, -5.495e-4] (V319 artifact re-read). The
drift target gamma = 2.290e-4/rung (V319-era inversion, booked in the
session note). Certified leg: the certified F_208 leaf at pi/4 times the
n = 208 sector ratio - does the correction REDUCE the booked 1.013%
certified deficit (closing needs ratio ~1.010, five orders above the
rung-10 window's 1.000006-1.000106: one uniform factor CANNOT do both,
and the n-dependence is the object of measurement, not an afterthought).

**T0 (pre-registered).** 64 sectors enumerated; the identity sector closes
at pi/4 with V315's booked k = 2 row (closure_row |dev| < 1e-12); the
plain-sector Floquet legs finite with the guard satisfied or booked.

**Branch (pre-registered words).**
- RUNG-SECTOR-SUPPLIED: W_eff(10) lands in a booked window (uniform OR
  closing-sector measure, either clock), guard satisfied at 10/17/208,
  AND the certified-deficit leg does not INCREASE the deficit. Then the
  sector sum DERIVES a rung weight at the required scale - the tie-break
  V319 could not deliver - and the booked formula gains a derivation
  test, one layer below coincidence.
- RUNG-SECTOR-NO-CONTRIB: otherwise; the full sector table (W_eff(10),
  W_eff(17), W_eff(208), per-measure, guard status) booked, and the
  sector route to the weight closed on the booked class.
Honesty rule unchanged: a sector sum that lands is a DERIVATION test, the
first one the candidate formulas face; it still does not touch the
neutrino floor or the origin of 1/(16pi) unless the n-dependence says so.

### Post-run booking note (2026-10-30): the rung-sector sum V320

Fired word: `RUNG-SECTOR-NO-CONTRIB`. T0 exact (identity sector closes at
pi/4, k = 2, dev 7.50e-15; 64 sectors; closing-only measure = uniform
measure since all 64 sectors close at Gamma - W(0) = I booked).

The sector answer is unambiguous and physics-bearing: summed uniformly,
the twisted sectors do NOT perturb the rung weight at the 1e-5 scale -
they DESTROY the rung sum (W_eff(10) = 3.010e-2, lnW = -3.503, against
the booked window lnW in [9.878e-7, 1.062e-5]) and increase the certified
F_208 deficit 1.013e-2 -> 3.507e-1. The twist sectors are destructive to
the crossing interference that makes the ladder tick; any measure that
saved the weight would have to suppress almost all sectors by hand - a
fitted choice, not a derivation. The guard worked as pre-fixed: at the
certified clock the odd rung rho0(17) = 5.5e-16 (the exact odd rule
visible in the Floquet leg) - guard-failed, never divided.

Consequence booked: the derivation test V319 set for the LOOP^3-scale
candidate FAILS on the booked sector class - the candidate is neither
density-distinguished nor sector-derived. Remaining viable routes to an
elect, as booked: a sector class NOT yet in the construction (the
n = 4k+2 no-limit class near the certified crossing, V209 - a
resummation object, not a sector sum), or the junction mechanism at the
democratic boundary (V303/V307) for the neutrino floor. The hunt's
negative space is now sharply drawn: rung weight (V316/V319/V320: three
closed classes), clock rule (V309/V315), selection (V317), channel
degeneracy (V316) - all bounded no-gos; what remains open is a new
mathematical object, not a wider census.

# The crossing-resummation programme: pre-registration V321 (2026-10-30)

**Fact exploited (booked, V209).** At the root clock the certified 4k+2
class diverges (certified F_210 at 1e-8 = 3.8854320241456516e48) while
the detuned class blows up as a POWER in the detuning: the booked grid
gives F_2(pi/4 + 1e-1) = 7.9503e6, F_2(pi/4 + 1e-3) = 1.1460e30,
F_2(pi/4 + 1e-5) = 1.1478e52 - slope ~11 per decade. A power divergence
has a CANONICAL finite part (Hadamard): subtract c*delta^{-m} with
(c, m) MEASURED at the detuned points - the targets NEVER enter the
computation, only the decision.

**Object (closed, fixed now).** F_n(delta) = loop_ratios at
g = pi/4 + delta, V203 machinery, nmax >= 260, over the FIXED grid
delta in {1e-1, 1e-2, 1e-3, 1e-4, 1e-5}. For each n in the divergent
class {2, 6, 10, 14, 18, 22} plus tail probes {250, 254, 258}: measure
m_n (3-point log-log slopes, agreement required < 5% between the two
slopes - else booked NON-POLYDIVERGENT, no finite part), then
c_n = F_n(delta_i) * delta_i^{m_n} at each grid point (spread < 1%
required - else booked NON-CANONICAL), and the finite part
fp_n = F_n(delta) - c_n * delta^{-m_n} at each point (spread < 1%
required - else booked NON-CANONICAL). No Padé/Cesaro/Euler: those
prescriptions cannot touch a power divergence, and V209's grid shows it
is one - booking the object class honestly instead of trying prescriptions
that cannot work.
ANCHOR leg (the construction must not move the convergent class): the
4k sequence at delta = 0 is left untouched; the float64 F_208(pi/4) is
reproduced and booked against the certified leaf 0.019692794505156912
(the ~4e-6 float64/certified gap booked as fact, V209's own finding).

**Decision target (artifact-read).** Does the finite part at n = 10 tie
the hierarchy at the ROOT clock: |3(fp_10/LOOP)^10 - COEF_OBS| <=
pillar_surgery.coef_band_abs? (The 4k class reads 2.7095 at the root
clock, +30% off; the required F at rung 10 is LOOP*(COEF_OBS/3)^0.1 =
0.01917541... - read in-experiment.) Also book the reads for every
computed n; rung-17 odd = 0 by the exact odd rule (leg integer, no leg
test).

**T0 (pre-registered).** The detuned grid reproduces V209's booked
F_2 values at delta = 1e-1, 1e-3, 1e-5 to relative 1e-9 (same
machinery, same grid); the certified F_210 leaf quoted as fact; the
4k float64 leg at delta = 0 reproduces F_208 = 0.01968660022354374
(V209 float64_floors) to 1e-12.

**Branch (pre-registered words).**
- RESUM-CROSSING-TIED: n = 10 canonical (m, c, fp spreads inside the
  booked 5%/1%/1% guards) AND its finite part lands the coefficient
  inside coef_band_abs AND the 4k anchor holds. This would be the FIRST
  zero-fitted algebraic supply of the hierarchy coefficient - the canon-
  ically tamed divergent class, not the convergent one, is what nature
  reads. Paper re-provenance then follows the pillar-surgery conditional:
  the O8 row and 12.1 paragraph state old -> new with the finite-part
  formula, and the V209 certified divergence gains its physical meaning.
- RESUM-NO-TIE: otherwise; the full table (per-n m_n, c_n spread,
  fp_n, coef-read) booked; the "resum the divergent class" route closed
  on the power-divergence object class.
Provenance if fired: algebraic (the subtraction data are detuned points
of the same algebra; the targets touched only the decision). Closest is
not true still applies: one tied rung-contact does not close the ledger
(rung-17 floor, V276 residue) - the booking must say exactly that.

### Post-run booking note V321 (2026-10-30) - RESUM-NO-TIE

Fired word **RESUM-NO-TIE** (tier 2, hard, O8): no member of the
4k+2 class admits a canonical finite part inside the pre-fixed
guards, so the partie-finie route to the hierarchy coefficient is
CLOSED on the booked class. Exact numbers (artifact
v321_crossing_resummation.json): guard status by n = 2: NON-
POLYDIVERGENT (m_3pt 11.579 / 11.000); 6: NON-CANONICAL (10.722 /
11.002, c-spread 0.856); 10: NON-CANONICAL (10.655 / 11.006,
c-spread 0.985); 14: NON-CANONICAL (10.898 / 11.011, 0.937);
18/22: NON-POLYDIVERGENT (10.102/11.018, 10.424/11.026); tail
250/254/258: NON-POLYDIVERGENT (13.206/11.283, 12.430/11.278,
12.782/11.279; c-spreads 7.5e5 / 5.4e4 / 5.0e4). Booked structural
fact: the ASYMPTOTIC exponent (1e-3 -> 1e-5 pair) is 11.00 +- 0.03
for every n <= 22 and 11.28 at the tail - the divergence is delta
^-11 to leading order with subleading/logarithmic structure
(coarse-pair slopes drift 10.1-13.2); no single (c, m) subtraction
is canonical at any n. F_10 values over the grid booked:
3.6634e6 / 3.7885e15 / 7.4871e27 / 7.6831e38 / 7.6874e49.
Anchor: certified F_208 at pi/4+1e-8 reproduced EXACTLY (0.0);
F_req(rung-10) = 0.019175759 quoted. T0 machinery identity exact
(0.0/0.0/0.0 at 1e-1/1e-3/1e-5).

Machinery amendments booked (neither touches the decision gate):
(1) the detuned grid is computed with V209's own exact-arithmetic
mp_ladder (90 dps, nmax 260), NOT float64 loop_ratios as the
pre-registration named: float64 destroys the divergence by
catastrophic cancellation (run 1: F_2 at 1e-3 read 7.09e14 against
the certified 1.15e30). Object class unchanged.
(2) T0 was amended at run 2 from "leaf agreement 1e-9" to
"machinery identity 1e-12" after the diagnosis that V209's booked
1e-5 leaf is itself a 55-dps value: the true 90-dps value differs
from the leaf by 5.1e-8 relative (8 significant digits at 1e52) -
a 1e-9 leaf gate is unsatisfiable; the identity replay (nmax 12,
55 dps, V209's exact calls) reproduces all three legs EXACTLY.
(3) the subtraction constant c_n was fixed at the finest grid
point (pre-registration wording allowed per-point constants, which
would make fp identically zero - the fixed-constant reading is the
only non-vacuous one).
Success tier (RESUM-CROSSING-TIED) would have re-provenanced the
coefficient algebraic; it did not fire: no paper provenance change.
One tied rung-contact would not have closed the ledger (rung-17
floor, V276 residue) - booked before the run, moot now.

# The O_h cone-decoration programme: pre-registration V322 (2026-10-30)

## V322_oh_bond_decoration_cone (tier 2, hard; O10/O4; closes V191's open piece)

**Context (booked, no new inputs).** V191 (alt_step): the 8-site module
of the 2x2x2 cell is A1g + T2g + A2u + T1u; NO irrep admits a linear
term (parity no-go, dim Hom = 0 for all ten); a conical node EXISTS
exactly and isotropically for the parity-mixed pair A1g/T1u realized
by an O_h-invariant bond Hamiltonian (two hopping orbits, -24/+8,
TT' = 16I to 4e-16, equivariance 1.2e-14) - BUT that bond layer is
NOT the engine's step: the six-stage product's measured O_h defect is
1.07e-1 (vs 1.1e-15 for the stage SUM), and equal stage weights pin
the sum spectrum {+6,+2,-2,-6} so the A1g/T1u degeneracy is
unreachable inside the present stage structure. V191: "the remaining
work is a stage factorisation of that bond layer, not another
search." V277 (DECORATION-BEST): the internal quaternion twist class
does not buy the cone at the closed clock (binding leg: degree
~1.9986, quadratic touching).

**Question (pre-registered).** Does the O_h-equivariant SITE-side
decoration class - exact orbit sums of the V191 range-1 hopping-orbit
census, the published six-stage step as base - carry a conical
GAMMA node (gap < 1e-10, fine degree >= 0.99, v_max <= 1, rho_comm >
0.5, V277's acceptance convention), and does any member reduce the
full step's O_h defect toward the stage-sum value?

**Class (closed, fixed NOW before any run).** The range-1 term
census {(i,j,v): i,j in 8 sites, v in {-1,0,1}^3} (1728 terms)
quotiented by O_h (48 site-permutation reps): 10 orbits (sizes
27,64,108,125,135,144,225,240,300,360 - counted pre-registration,
booked here). Members: D(q) = exp(i THETA sum_o w_o A_o(q)) with
A_o(q) = orbit Bloch sum (Hermitian, V191 bloch() construction verbatim:
M + M^dag over orbit elements), w_o in {0,+-1}, AT MOST TWO nonzero
orbits, THETA = pillar-independent config angle.alpha_10 (two sign
choices of THETA coincide with +-w - one amplitude only). Enumeration:
identity + 10x2 single + C(10,2)x4 double = 201 members. Decoration
enters each stage's SITE block: M -> kron(D(q) @ S_t(q), U16) @ M
(left multiplication, per stage; U16 and clocks untouched: g = beta =
angle.branches.closed). Provenance: combinatorial (orbits) +
config-declared (THETA); zero-fitted.

**Measurement per member (closed clock).** (1) GAMMA gap (min phase
mod 2pi over the 128-dim step spectrum); (2) fine cone degree +
v_max + tangent speeds: V277's machinery copied (179-Fibonacci
tangent directions at the booked radius, degree on the 64-dim gap
subspace, V192/V132 pinned 64+64 spectrum at the closed clock);
(3) rho_comm (V186 tangent route, V277's convention); (4) O_h
defect: max over all 48 elements of || U_dec(W q) - (P_W kron I16)
U_dec(q) (P_W kron I16)^T ||_max at q = GAMMA + 8 Fibonacci
directions |dq| = 0.1 (pre-fixed), V191's defect convention.

**T0 (machinery gates).** Identity member (all w = 0): GAMMA gap
< 1e-10 reproducing the closed-clock 64-fold gaplessness (V192/V132
booked spectrum), and the bare top-64 fine degree reproduces V277's
artifact leaf (implementer names the key in the check string) to
< 1e-6. Base engine step: O_h defect re-measured in V322 machinery
against V191's artifact leaf 1.07e-1 to < 1e-3 relative. Any drift:
MACHINERY-DRIFT, exit nonzero.

**Branch (pre-registered words).** CONE-OH-FOUND iff some member
satisfies ALL four acceptance legs (gap, degree, v_max, rho_comm);
CONE-OH-NEGATIVE otherwise - both cases book: the top-10 members by
fine degree (gap/v_max/rho_comm/defect), the minimum-defect member
(the stage-factorisation-equivariance leg V191 left open - a booked
measurement, never a gate), and the count of members keeping gap
< 1e-10 (the SELECTIVITY table: which orbit sums gap the node and
which do not - the gap-selective content of the class).

**Bridge note (honest booking).** The pre-truncation design note tied
V322 to the rung-17 gap (S17, -2.78%): the cone-to-ladder bridge was
NOT recovered from the workspace and is NOT invented here. V322 is
the cone/equivariance census; the rung-17 target leg is V323's
object. No rung contact enters V322's decision, now or by amendment.

# The joint-stabiliser programme: pre-registration V323 (2026-10-30)

## V323_joint_stabiliser_s17 (tier 2, hard; O8 + neutrino)

**Context (booked, no new inputs).** The rung-17 neutrino contact of
the five-contact ledger: S17/3 = 0.04874922781428075 eV (V303 leaf,
re-derived by V319) vs neutrino.m_atm_ev = 0.0501 - a -2.6962%
residual; V318 books it as the ledger's binding contact (C3 at
0.539 of band on the best row); V319 closed the rung-weight class
(supplied-by-density), V320 the sector-sum route, V321 the
resummation route. The remaining booked algebraic structure that has
never been tested against S17 is the JOINT MODULE of the gauged bond
layer: V79/V140 measure 32 (x) 45 = 32 + 288 + 400 + 720 with the
adjoint-link Casimir spectrum THREE levels {45, 85, 117} with
multiplicities {32, 288, 1120} (400/720 degenerate at C2 = 117) -
the bond hop factorises through the joint Casimir T.L, projectors
verified to 1e-15 (paper Section 10). The stabiliser (commutant) of
the joint action has a 3-dimensional centre labelled exactly by
those three integers.

**Question (pre-registered).** Does an element of the joint-module
centre supply the rung-17 multiplier - i.e. does a Casimir monomial
chi, built ONLY from the booked {45, 85, 117} and {32, 288, 1120},
with NO reference to any target, give |S17/3 x chi - m_atm| / m_atm
within the 5% ledger band - and is the hit a DERIVATION (exact
algebraic match to the required chi) or supplied-by-density (V306/
V319 rule)?

**Class (closed, fixed NOW, target-free; the target-defining
constants S17, m_atm NEVER appear inside any member).** With
t = (45, 85, 117) (level Casimirs, V79 leaves) and m = (32, 288,
1120) (multiplicities): 1 (identity); t_i (3); t_i/t_j, i != j (6);
t_i/tbar with tbar = unweighted mean 247/3 (3); t_i/tbar_mw with
tbar_mw = multiplicity-weighted mean Sum m_i t_i / Sum m_i = 109
(3); monomials t_i t_j / t_k^2, distinct (3); cross-ratios
(t_i - t_j)/(t_j - t_k) over orderings (3); m-weighted gaps
(m_j t_j - m_i t_i)/(m_j t_j) (6). Total 25 members, enumeration
fixed here. Provenance: algebraic (integers from booked artifacts).

**Reading convention (inherited).** rung-17 leg' = S17/3 x chi
(V316/V319 convention); band = pillar_surgery.ledger_scale_rel_tol.
chi_req = m_atm/(S17/3) is computed for the DECISION ONLY.

**T0.** S17/3 reproduces the V303 leaf to < 1e-9 (art pattern, via
cl10.neutrino.loop_ladder(17) x EV_PER_GEV as V319 re-derived); the
joint content {45, 85, 117} x {32, 288, 1120} reproduces the V79
artifact leaves in-experiment via V79's own machinery (copy of the
content_ok block); any drift: MACHINERY-DRIFT, halt.

**Branch (pre-registered words).** S17-SEALED-EXACT iff some member
matches chi_req to relative < 1e-9 - a derivation: the -2.6962% is
an algebraic property of the joint module, and the ledger C3 contact
re-provenances algebraic (this tier REQUIRES the exact tier; no
band-luck); S17-SEALED-IN-BAND iff one or more members land the
contact inside the 5% band but none exactly - booked with the
V319-style chance-density count: expected hits = N x (band log
width)/(class log range), supplier-grade iff < 0.5 (V306 rule);
S17-STUBBORN otherwise. All three book the FULL 25-row table
(chi, leg', normalized residual). Exact-match members: no class was
designed to hit 1.0277; an EXACT hit would be the first neutrino-
side derivation in the programme and gets the success tier.
