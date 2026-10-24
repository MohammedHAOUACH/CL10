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
