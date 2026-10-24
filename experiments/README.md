# The verification program

168 experiments (the 156 of the publication programme plus the 12 exploration
experiments of the `alt_*` programmes: 5 `alt_gauged`, 4 `alt_nonpert`,
3 `alt_step`), grouped by
the part of `doc/final_publication.md` they validate.  Each experiment
recomputes a published number from the construction, states the provenance of
every input it uses, and writes a JSON verdict to `<experiment>/results/`:
`pass`, `status`, and the tier of the claim it tests - `0` the root, `1` the
object and its kinematics, `2` the dynamics, `3` the anomaly-free matter and
the couplings, `4` exploration, `5` the falsifiers, `6` the provenance and
audit checks.  **84 of the 156 publication experiments still write the older
verdict (`pass`, `status`, `data`) with no `tier`, `property` or `wall_time`.**
They count in PASS/FAILED but sit outside the tier accounting, and the runner
prints that denominator on the hard-tier line instead of hiding it.  `python3
run_validation.py` runs all of them and rolls the verdicts up into
`results/validation_summary.json`.

The papers these experiments validate: [Paper I](../doc/paper1_core.md) ·
[Paper II](../doc/paper2_gauge_gravity.md) · [Paper III](../doc/paper3_phenomenology.md) ·
[the unified master](../doc/final_publication.md) — compiled PDFs in `../doc/`, on the
[v1.0 release](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0) and on
[Zenodo](https://zenodo.org/records/23321228); cite:
[10.5281/zenodo.23321106](https://doi.org/10.5281/zenodo.23321106).
Back to the [CL10 README](../README.md).

A script never hard-codes its own category: sibling experiments and their
verdicts are addressed by name through `cl10.paths.exp`, so moving an
experiment between categories moves a directory and nothing else.

Each experiment names the open problem its verdict speaks to (O1–O10 of the
paper's Section 16) inside its own verdict, and rewrites its verdict file
whenever it is re-run.

| category | what it establishes | paper |
|---|---|---|
| `algebra_root` | the root itself, the closure of the coupling, the provenance audit | §2, §9, §16 |
| `spacetime_matter` | the $(1{+}3)$ subalgebra, one anomaly-free generation, the $45$, the statistics, the spectral mirror, the flavour structure | §3–§6, §13.5 |
| `dynamics_spectrum` | unitarity, the causal radius, the $n$-particle sector, the doubling and its Wilson gap, the node geometry, the cone-isotropy no-gos, the graviton content | §7, §8, §9 |
| `gauged_layer` | the adjoint bond layer: Gauss constraint, centre obstruction, area law, vertex correction, ground state, confinement | §10 |
| `gravity` | the discrete Einstein–Hilbert action, curvature from an arbitrary metric, $3{+}1$ evolution | §11 |
| `scales_loopfactor` | where a loop factor can come from, the closed-loop factor, the ladder length, the composition law, the rung arithmetic of every scale, the traced-loop weight functional | §12.1 |
| `couplings_unification` | $\alpha_{10} = 1/25$, the desert solution, the proton lifetime, every proposed supplier of the desert spectrum | §13.1, §13.4 |
| `neutrino` | the seesaw ladder, the $\tfrac13$ candidate, the protection of the two massless states, the $\Delta L = 2$ no-go | §12.3 |
| `cosmology_darkmatter` | the $118$ GeV ladder mass, the carrier no-go, the built scalar, the relic abundance, the portal window | §12.2, §14.2 |
| `baryogenesis` | the absence of every Sakharov condition, the asymmetry product, the washout, the real-form obstruction | §13.3 |
| `family_number` | why the count $3$ is an input and what its symmetric alternative costs | §4, §16 (O3) |
| `phenomenology` | the lattice photon driven through hydrogen (Fig. 2) | §6 |
| `alt_gauged` | the exploration programme: the two-fermion S-matrix of the gauged bond layer (V79-V101) built and measured - scattering phases, vertex O_h content, one-loop flow, block RG - against the cone-isotropy gate (tier 4) | §8.4, §16 (O9) |
| `alt_nonpert` | the exploration programme: non-perturbative gauged-sector cone probes - exact cone dressing by the flux operators, 4-site chain/ring kernels and S-matrices, IR decoupling across volumes (tier 4) | §8.4, §16 (O9) |
| `alt_step` | the exploration programme: the parity no-go for a conical node and its reruns inside the engine's own six-stage circuit class and with the root's sector-swap involutions - the cone is outside the class, searches included (tier 4) | §8.4, §16 (O9) |

The `alt_*` programmes run on the **numpy core**: each script pins
`CL10_BACKEND=numpy` before importing `cl10`, and `cl10/altcell.py` is
host-only by contract (every call into an engine-routed module is moved back
to the host where it enters). Their linear algebra is host-side, so their
verdicts do not depend on what else holds the card; the runner lists them
OFF-CORE. The GPU engine stays with the publication programme.

The same self-containment rule holds for the clock. The closure discovery
(V203-V206) made the CLOSED branch the config default, and the census and
gate programmes are stated at the closed root. An experiment whose own
anchors are 1/5-reference-branch numbers - the word-census family
(V126/V200/V202), the relic DM chain (V123/V147), the gauged cone-flow and
block-RG chain (V175-V179), the cone-flux dressing chain (V180-V183) - pins
`CL10_ANGLE = 0.2` at the top of the script, before `cl10` latches the
configuration; closed-branch experiments take their angles explicitly
(`pi/4`, `pi/2`). The pin direction is one-way: a 1/5 anchor never silently
reads the closed clock.

## Which experiment backs which numbered claim

The paper's numbered claims are not named inside the scripts. This is the
map. A claim counts as **backed** when an experiment measures the quantity
the claim asserts, not when a script merely mentions the word.

| claim | measured by | what is measured |
|---|---|---|
| Theorem 1 - the charges of the `16` are exactly one SM generation | `V2_matter`, `V23_full_sm_qca`, `V26_family`, `V70_family_index`, `V107_mixing_from_kappas` | the charge table of the six fields, the family index |
| Theorem 1 - the `5/3` hypercharge normalisation | `V32_qca_qft`, `V52_parameter_audit` | the normalization of the diagonal generator |
| Theorem 2 - all eight gauge and mixed anomalies cancel | `V11_anomaly`, `V120_entangling_bond`, `V121_three_axis_bond`, `V52_parameter_audit`, `V78_lambda_protection` | the eight cubic and mixed coefficients of one `16` |
| Theorem 3 - the `45` adjoint and its twelve gauge bosons | `V131_algebraic_closure`, `V79_gauged_bond_layer`, `V129_gauged_layer_gap`, `V130_string_tension`, `V175_gauged_smatrix`, `V176_gauged_vertex_content` | the adjoint content of the bond layer and its vertex |
| Theorem 6 - the closure axiom quantizes `g`, step order 8 | `V131_algebraic_closure`, `V134_pi2_branch`, `V135_coupling_selection_closed`, `V133_desert_pi64`, `V132_gapless_census_pi4` | the order of the on-site update, the finite internal spectrum, the branch census |
| Proposition 1 - the (1+3) Minkowski subalgebra | `V1_spacetime`, `V93_cartan_junction`, `V137_broken_flux` | the signature of the ten-generator frame, the Cartan junction |
| Proposition 2 - lattice chirality and the doubler census | `V63_bz_chirality`, `V132_gapless_census_pi4`, `V139_nodal_plane_crossing` | the `16`-sector closure, the gap at Gamma, the 12 M nodes and their charges |
| Proposition 3 - the momentum derivative is internal-blind; the Wilson factor is a spectator of the tangent algebra | `V191_parity_dirac_cone` | the tangent operators of the joint step, and the stage-sum levels |
| Proposition 4 - no cone over the 880-point entangling-bond space | `V186_cone_inverse`, `V120_entangling_bond`, `V121_three_axis_bond`, `V115_cone_isotropy` | the tuned bond census as (gap, anisotropy, rho_comm) triples |
| Proposition 5 - the rigidity wall at gap 8e-3 | `V185_separation_lemma`, `V187_rigidity_frontier`, `V186_cone_inverse` | the (gap, rho_comm) frontier over 5292 points |
| Proposition 6 - the parity no-go, and the cone outside the class | `V191_parity_dirac_cone`, `V192_engine_two_angle_cone`, `V201_root_swap_cone` | `Hom_O_h(T1u, D (x) D*) = 0` for all ten irreps; the constructed node; the search for the same node inside the six-stage product, which does not find it, and its rerun on the Dirac module with the 10 root sector-swap involutions `h_a w`: the gapless points turn non-Abelian (rho to 0.80) but stay flat - not found either |
| `det U = 1`, exact unitarity, the +/- phase pairing | `V43_qca_spectrum`, `V131_algebraic_closure`, `V52_parameter_audit`, `V78_lambda_protection`, `V117_real_form`, `V125_cp_flux_washout` | the determinant of the step and its eigenphase pairing |
| front speed exactly one Planck cell per stage | 21 scripts across 8 programmes, plus `cl10/supercell.verify()` invariant 2 | the support growth of the real-space step |
| Nielsen-Ninomiya bookkeeping | `V63_bz_chirality`, `V132_gapless_census_pi4` | the node charges and their sum |
| the cubic equivariance of the bond product (Section 8.4) | `cl10/supercell.verify()` invariant 6, `V191_parity_dirac_cone`, `V192_engine_two_angle_cone` | `W P(q) W^T = P(Rq)` per momentum: 1.1e-16 at Gamma, 1.07e-1 generic, 1.40e-1 M, 1.46e-1 R |
| the cell module of Axiom 1 and its B-parity (Sections 2, 4, 12.2) | `cl10/cell.py`, `cl10/supercell.verify()` invariant 7, `V194_carrier_in_cell`, `V195_cell_redefinition`, `V196_cell_definition` | 36-level step unitary, spectrum +-paired, `det = 1`, quasi-energy sum zero on the full cell; `[G,U] = [G,45] = 0`; the ladder identity `m_nu = m_B^2/m_A` exact between cell levels; the rung phases `loop^3`, `loop^10` reproduced by the six-stage product; `M(q) = kron(P(q), U^6)` |
| the rung arithmetic and its joint falsifier (Table 4, F3) | `V197_rung_arithmetic` | rung exponents of every quoted scale: the four ladder rows integer-exact, every other scale off-integer by 0.06–0.42 rungs, and the [rung-10 DM, 0.05 eV neutrino] pair forcing the seesaw partner to rung 2.72 |
| the traced-loop weight functional (Section 12.1, O8) | `V200_loop_weight_functional` | `F_n = N_n r_n I_n` exact to 2.3e-16 for n = 2..10 on the word census; word norms angle-independent; `F_10` crosses 1/(16 pi) at 8 angles, deficit factor 0.9644 at the root angle |
| the loop angle-selection census (Section 12.1, O8) | `V202_loop_angle_selection` | `F_10(eta)` is an exact real cosine functional: even Casimir eigenvalues, 63 even phase classes, period pi, eta -> pi - eta pairing, coefficient mirrors to 2e-17; the four crossings of 1/(16 pi) in (0, pi/2] contain no closed or rational angle of the root - nearest 0.198100, g_10 = 1/5 short by 1.9e-3 (branch ANGLE-OUTSIDE) |
| the closed-branch self-consistency (Section 12.1, O8) | `V203_closed_branch_coefficient` | step and phase both at pi/4: odd-n selection rule exact (3.3e-15), rho(MG) = 0.978030 angle-independent, the first crossing stays at n* = 208 (deficit 0.9895) - the published ladder length survives closure - while F_10(pi/4) = 55.03 (the eta-law 0.1949 is a step@1/5 hybrid); the 12 M channels exactly gapless with identical pinned spectra (64/128), gapless set Gamma + six face-diagonal lines, gapped bulk; required node-set factor 2.308 not natural |
| channel-blindness of the O8 deficit (Section 12.1, O8) | `V205_channel_filtered_deficit` | required insertion amplitude $0.96388$ against measured $0.98948$ at the crossing, blind to channel: sector phases fourth roots $\{1,-1,-1\}$ at $\pi/2$ and $\{1,i,i\}$ at $\pi/4$, $\mathbb{Z}_3$-filtered traces empty-kernel to 3e-14 on the pinned channel spectra, per-Cartan amplitudes bimodal $\{0.90, 0.99\}$, none at $0.964$; observed coefficient with current constants $2.07656$; adjacent lengths $6.5$/$1.1$ bracket, none match |
| the clock tension of the hierarchy deficit (Sections 12.1, 13.1; O4 and O8) | `V206_clock_tension` | closure-class readings of the coefficient: 2.7002 (pi/4, n* = 208), 2.2994 (pi/2, n* = 102), 2.7084 (3pi/4) - the closest closed reading is +10.7% off the observed 2.0766; the observed value is read by the SAME self-consistent functional at rung 10 at g_obs = 0.1999980, the rational clock 1/5 to relative 1e-5 (rung-10 crossing at 0.200133); the exact odd-length cancellation exists only on the order-8 clock family (odd loops 1e17 at pi/2); dead-zone collapse to n = 3 between 1/3 and 1/2; no principled link mixture or sector reading bridges it |
| the A->BB channel bound on O5 (Section 13.3) | `V199_a_to_bb_channel` | `Gamma_B/H` from 4.7 to 1.2e5 over the 13 root bond strengths; 9/13 drive `Br(A->LH)` below the 6.4% CP-strength ratio; survivable below kappa 3.04, negligible below 0.080 |
| the spin-2 content census (Table 8 F4) | `V230_quintet_helicity2_content` | the quintet as one irreducible 5-dim icosahedral representation (exact Q(sqrt5) characters, commutant 2/1/1), the massless-branch census per bond channel, the closed-clock gaps 6.6047 (diagonal) and 16.7287 (axis, Eg half, T2g bubble exactly 0) |
| the diagonal-bond graviton channel (Table 8 F8) | `V229_diagonal_bond_quintet` | the T2g stress vertex on the fifteen star axes (0.6325 on 6/15 entries, nonzero on all fifteen), the exact leading-order isotropy identity (9e-16), the doubled ladder {0, 2 g10, 4 g10} x |d| |
| the Lorentz contact (Table 8 F9) | `V235_node_velocity_contact`, `V236_dressed_hop_isotropy` | node velocity tensors 0.296/1/0.295/0.9999 against the electron-sector SME ceiling 8.0e-15 (factor 3.7e13 past it); the closed-clock dressing is degenerate and cannot move the cone (delta_v/v = 1.0000) |
| the clock arbitration (O8, Table 8 F2 clock side) | `V233_clock_arbitration` | 2x2 score over {pi/4, 1/5}: neither clock reads the observed hierarchy coefficient 2.076564 +- 1e-4 (2.7095 certified, +30.5%, 6330x the uncertainty; 2.0881, +0.55%, 115x); the 1/5 agreement is a fit at g_obs = 0.1999980; the clock books as FITTED |
| the certified star-module loop ladder (O8) | `V231_star_loop_functional_certified` | 416 exact box means certified at 1e-40 for $n \le 12$; no certified crossing on the closed clocks (closest approach 7.75x loop), a one-length undershoot at $n = 6$ on $1/5$, grid floor 0.24 out to $n = 208$; the star module does not reproduce the loop factor |
| the bounded insertion-amplitude search (O8) | `V232_insertion_amplitude_search` | 60 candidates inside the closed structure set (14 named structures plus the 46 V231-certified star-module deficits); none within 1e-6 of the required 0.96388, none within 1e-3 among the closed-input candidates (closest closed-input 9.8e-3; closest overall 5.3e-4, off-closure); O8 stays open |
| the rung phases and neutrino contacts (Table 4, F3) | `V234_rung_phase_spectrum` | only loop^3 numerically realized on the 36-level step (loop^10 = 9.71e-18 below float64 phase resolution); rung_B = 10 pinned by the DM window admits no rung_A landing inside the neutrino contacts; (3, 10) is a config declaration |
| the desert threshold matching and bounds (O4, F2, F10) | `V237_desert_threshold_matching`, `V238_desert_spectrum_bounds` | the (1,2) shift demand +13.1 e-folds = 1.12x the log-spacing to M_X - not realizable by the spectrum's own content; every fitted threshold 14x-1e9x beyond its fetched experimental ceiling - F10 open |

## The experiments

### `algebra_root`
| id | verdict |
|---|---|
| `V0_root` | $\mathrm{Cl}(10,0)$ is the root: the ten generators, $w^2=-1$, the algebra basis |
| `V52_parameter_audit` | the provenance audit: algebraic / conventional / fitted / open, and the CP-odd invariants that vanish |
| `V59_four_limits` | the four limits reduce to the seven algebra invariants plus one angle |
| `V78_lambda_protection` | $\Lambda = 0$ is protected by $\det U = 1$ under arbitrary traceless perturbations |
| `V131_algebraic_closure` | the root quantizes its own step: $g \in \{\pi/2,\pi/4\}$, the order-$8$ spectrum |
| `V134_pi2_branch` | the $\pi/2$ branch is a pure lattice translation |
| `V135_coupling_selection_closed` | the closed coupling $\pi/64$ against the one-loop data |
| `V196_cell_definition` | the 36-level cell of the engine (`cl10/cell.py`) as a root fact: the ladder identity $m_\nu = m_B^2/m_A$ holds exactly between cell levels, the B-parity is exact against the step and the $45$, the rungs $\mathrm{loop}^3$, $\mathrm{loop}^{10}$ are reproduced as phases by the six-stage product, transport factorises as $\mathrm{kron}(P,U^6)$, anomaly bookkeeping intact; +- pairing, `det U = 1` and vanishing quasi-energy sum on all 36 levels |

### `spacetime_matter`
| id | verdict |
|---|---|
| `V1_spacetime` | the $(1{+}3)$ Minkowski subalgebra inside the root |
| `V2_matter` | the $16$ Weyl module is exactly one generation |
| `V3_gauge` | the $45$ is the $\mathrm{Spin}(10)$ adjoint |
| `V4_statistics` | the $\mathbb{Z}_2$ grading is the spin-statistics sign |
| `V11_anomaly` | all eight gauge and mixed anomaly coefficients vanish |
| `V24_cpt` | the $16$ and $\overline{16}$ are CPT conjugates; the step is CPT symmetric |
| `V7_qca` | the step is exactly unitary |
| `V23_full_sm_qca` | the full Standard Model is a 0-parameter unitary QCA on this module |
| `V32_qca_qft` | the step gives the SM effective action with 0-parameter electroweak breaking |
| `V43_qca_spectrum` | one 32-level spectrum, ten named projections |
| `V28_yukawa` | the democratic rank-1 Dirac mass matrix |
| `V45_yukawa_matrix` | the $3\times3$ Yukawa matrix, the TBM mixing, $\sin^2\theta_{12} = 1/3$ |
| `V49_quark_ckm` | the computed CKM is the identity |
| `V67_breaking_chain` | the Wilson direction preserves the Model algebra and breaks to $\mathrm{SU}(3)\times\mathrm{SU}(2)\times\mathrm{U}(1)^3$ |
| `V104_mixing_manifold` | the natural vev placements are excluded or vacuous |
| `V107_mixing_from_kappas` | gauge-class bond couplings are generation-blind: they cannot supply mixing |
| `V111_generation_bond` | the parameter cost of a generation-space dressing (7 / 13) |

### `dynamics_spectrum`
| id | verdict |
|---|---|
| `V36_local_qca` | the causal radius is 1 cell per stage; the light cone is realized |
| `V41_smatrix` | the step is the $S$-matrix; the first processes and the cross-section bound |
| `V48_3particle` | the $3$- and $4$-particle amplitudes are UV finite; bounded for all $n$ |
| `V63_bz_chirality` | chirality is a dynamical symmetry; the 12 gapless node momenta, $\sum$ chirality $=0$ |
| `V65_wilson` | an internal Wilson factor of the same step gaps all 12 doublers |
| `V102_lorentz_isotropy` | the node dispersion: 1200 sampled directions, anisotropic linear cone |
| `V115_cone_isotropy` | no isotropy point anywhere on the Wilson family |
| `V118_little_group` | tangent velocities and node little groups: three nodal planes, not a cone |
| `V120_entangling_bond` | the entangling bond layer reaches $0.071$ of the required commutator |
| `V121_three_axis_bond` | the three-axis covariant bond reaches $0.658$ and costs the causal front |
| `V128_f4_graviton_content` | no massless spin-2 content in the implemented spectrum (**F4**) |
| `V132_gapless_census_pi4` | at $\pi/4$: $\Gamma$ plus six face-diagonal nodal lines |
| `V138_gamma_spin2` | the zone centre carries $T_{2g}$ but never $E_g$: no lattice spin-2 in the qubit link, the quintet lives in the directional link layer |
| `V139_nodal_plane_crossing` | the grid census was too coarse: an exactly gapless orbit of $24$ points on the nodal planes at $g_{10} = 1/5$, little group of order $2$ |
| `V142_frame_averaging` | frame averaging over the ten generators is exactly blind to the cone anisotropy: in four random SO(10) frames the node gap moves by 1.8e-15, the anisotropy by 0, the tangent commutators stay at the noise floor — all momentum dependence lives in the sub-site factor of the Kronecker bond layer |
| `V185_separation_lemma` | the separation lemma: every tangent operator is exactly $A_k(q) \otimes I_{16}$ (defect 4.7e-16 over 40 momenta), so the whole low-energy velocity structure lives in the 8×8 sub-site factor |
| `V186_cone_inverse` | the cone inverse problem over the axis-covariant entangling family (880 bonds, 220 triples): no cone in the family; the implemented corner triplet is exactly radius-1 with the Wilson factor off (leakage 8.7e-17) |
| `V187_rigidity_frontier` | 5292 sweep points around the top-commuting bonds: gaplessness and non-commuting tangent algebras are mutually exclusive, the frontier floor is a wall at gap 8.0e-3, no power law links the regimes |
| `V207_corner_pair_zero` | **the corner-pair measurement at the two-triplet objective**: over the same Wilson family both lowest corner triplets go to zero simultaneously — six states with all six phases below 2.5e-16 at $w = 0.241747$, $\theta = 0.913137$ ($g = 1/5$ member), isolated from the next level by 8.7e-3, six such points in one scan of the plane; still Abelian (max commutator 5.3e-12 against $\sum_j\lVert S_j\rVert^2 = 0.82$) and exactly planar to 4.8e-12 over 2196 directions — six nodal planes, not a cone, and $\sum_j S_j^2$ is not a scalar on six states ($7.2\times10^{-2}$ about the mean 0.137, exactly the branch splitting). The objective has to be taken on the two triplets at once: the corner spectrum is particle-hole paired, so its two lowest phases sit in one triplet, and minimising the wrong objective lands at $w = 2.8894$, $\theta = 1.5717$, where the six phases are in fact $1.4\times10^{-4}$ and $2.6\times10^{-4}$ — excluded |
| `V212_continuum_census` | **the continuum census, locus by locus**: little groups read from the measured point group ($\\Gamma$ 48, $X$ 8, $M$ 4, $R$ 6) with their Molien counts; isotropy is forced at exactly one locus, $\\Gamma$ — which is quadratic there (exponent 1.9999998, gap 1.4e-16) — while the loci that touch linearly ($M$, exponent 1.0004, and the whole $\\Gamma$–$M$ segment with max gap 6.4e-16) carry 3 quadratic invariants, so anisotropic Fermi velocities there are required by symmetry rather than tuning: both jaws of the trap measured. The cubic group's first anisotropic invariant is degree 4, the icosahedral group's (order 60, from the $[2,3,5]$ presentation) degree 6 |
| `V213_relevance_exponent` | **the quartic anisotropy is marginal, not irrelevant**: over 179 directions the node's leading exponent is $\\alpha = 0.9998$, and the anisotropy saturates at 1.0000 because 4 of 179 directions are exactly nodal at every radius of the ladder — the saturated ratio therefore carries no exponent ($p = -0.0001$ over 7 radii), which is what the fixed-probe reading of that same ratio reports. The block ladder with window $\\pi/(2b)$ reproduces the exponent to 0.000 from an independent ladder; the node's little group (order 4, 3 quadratic invariants) already allows an anisotropic invariant at degree 1, predicting $p = 0$. The fixed-probe convention (V178's windows, probe held at $10^{-3}$) returns slope 0.0000 because it measures the downfolded subspace, not the flow |
| `V214_star_module_step` | the shift module replaced by the icosahedral star, measured: 6 independent generators of norm 1.902113 (the Lieb–Robinson radius per step), rank 6 in three dimensions with relation rank 6 (exactly the six opposite pairs), dense (shortest element 1.0705 → 0.2527 between word bounds 1 and 5); six equiangular lines at 0.447214 with fourth moment isotropic to 2.2e-16 where the three cubic axes are off by 0.400. The price of the move is exact: the cubic axis the engine hopped on is not a module point |
| `V215_star_gapless_locus` | the new module's gapless locus: 125 momenta out to $|q| = 1.6516$, smallest gap 8.9e-18, 13 exactly gapless grid points and no whole ray gapless to the zone edge; 12 of 12 rays leave with exponent below 1.4 (0.996–0.999 along the two-fold directions), so a linear cone candidate exists; the Clifford test at $\\Gamma$ has 512 of 1024 phases zero with first-order splitting 1.7e-13, 2.0e-13, 1.3e-13 — no linear term, $\\Gamma$ stays quadratic on the new module too |
| `V216_exclusivity_boundary` | the V187 exclusivity boundary measured off-axis: $\\kappa = \\max_{ab}\\lVert[V_a,V_b]\\rVert/\\max_a\\lVert V_a\\rVert^2$ on the whole gap subspace over 33 cubic and 34 star momenta. Away from the gapless set $\\kappa \\sim \\mathrm{gap}^\\alpha$ with $\\alpha = 12.987$ (cubic, correlation 0.656) and 8.860 (star, 0.566) — a law of the gapped region only, it does not extrapolate onto the touching. On the gapless set the cubic $\\Gamma$–$M$ line is Abelian ($\\kappa \\le 1.3\\times10^{-11}$, multiplicity 32/64) while all six star axis momenta are non-Abelian, $\\kappa$ to 23.062 at gap 0 with multiplicity 128; at the most non-Abelian star momentum the normalised squares of the tangent generators are $[57.111, 0, 16.704]$ and the largest anticommutator 84.2: a linear, non-Abelian, strongly anisotropic touching, not a cone |
| `V217_internal_sector_survival` | the internal sector is blind to which shift module carries it: $w^2 = -I$, $g_1^2 = I$ and $\\{g_1,w\\} = 0$ all to exactly 0, sector projectors of rank 16/16, the 32-internal step commuting with $1\\otimes w$ to 5.0e-16 (cubic) and 6.7e-16 (star), the step block diagonal in the $16/\\overline{16}$ split to 4.2e-16 on the 2048-dimensional star step, its 16-block equal to the committed step to 3.6e-15, matter and antimatter isospectral to 6.2e-15 |
| `V219_nodal_line_anisotropy` | the anisotropy of the new module's nodal lines: the gap stays on the 1e-16 floor along the line and closes transversely with exponent 0.995; the transverse velocity profile over 96 azimuths at radius $3\\times10^{-3}$ has leading azimuthal degree 2 with odd-harmonic power 1.2e-12 (the period $\\pi$ of the two-fold little group) and total anisotropy 0.680 of its mean against 0.500 for the cubic $\\Gamma$–$M$ line; the line's own two-fold invariance is exact in the linear-response limit (defect slope 1.0000 in the radius); the two-band Clifford test fails at its maximum, defect 2.000 with generator norms 2.7e-6, 0.177, 1.318 |
| `V220_nodal_line_charge` | **no charge on the new module's nodal lines**: at the line the touching is a cluster of 128 of 1024 eigenphases at phase zero (cubic $\\Gamma$–$M$ control: 64), largest inside 2.3e-14 and first outside 3.2e-01, and consecutive eigenvectors of the single closest phase overlap at 1e-14 worst case, so no band-resolved Berry phase is defined on it (inside a degenerate cluster the eigensolver's basis is non-orthogonal - condition up to 4.4e+02 - so the frame is re-orthonormalised before transport). On that frame the cluster's holonomy is unitary to 4.1e-04 and the cluster is isolated along every circle by 5.2 to 1479, so its determinant phase is an invariant of the step restricted to the touching - and it depends on the circle: -0.0998 / -1.3300 / -2.7335 rad at radii 3e-3, 1e-2, 3e-2 in the plane across the line, -2.2902 / -1.3300 / -0.3790 tilted 45°, spread 2.6337, where a monodromy is constant; the cubic control gives exactly 0.00000 at the two smallest radii. Curvature flux, not a quantized charge. The mass test against a *computed* symmetry operator (commutant of the internal map by Cesàro over 240 powers, obstruction 4.2e-15) opens a gap at the line growing linearly in the coupling: 1.277e-04 at $\\theta = 10^{-3}$, 1.277e-03 at $10^{-2}$, ratio 10.00; the little group of the line is 2 of the point group's 60 and those elements fix the spectrum at the line to exactly 0; the cluster's closing exponent is 1.006 across the line, 1.008 tilted, 1.000 for the cubic line |
| `V225_covariant_clifford_cell` | **the cell the exclusivity boundary points to, built and measured**: bond operators $U_d=-\\hat v_d\\cdot\\vec\\sigma$ intertwine $2I$ exactly (3.5e-15 over the 120 lifts, $U_d^2=1$ to 2.2e-16, $U_{-d}=-U_d$ to exactly 0); the quaternionic assignment intertwines just as exactly but squares to $-I$ and its cell has no node at all (gap 1.571), the scalar assignment is trivially covariant ($U_{-d}=U_d$, defect 2.0); the step is $2I$-invariant at the zone centre (2.5e-15) and inversion is broken (1.414 = the defect that lets the linear term in), with the equivariant census pushing the first anisotropy to degree 5 while admitting a mass (degree-0 count 1); the node closes exactly when the internal map rides at exponent 0 or 4 mod 8 ($U_{16}^8-I=4.9\\times10^{-15}$ against $U_{16}^4-I=2.0$), 192 zero eigenphases of 384, the gap cluster splitting as $\\mathrm{sym}^1\\oplus\\mathrm{sym}^3\\oplus\\mathrm{sym}^5$ with multiplicity 16 each by projectors built from the group; the tangent algebra on the whole cluster is Abelian (one norm 10.7600 to 3.6e-15, commutators 2.1e-14) and the Clifford identity is restored on the 32-dimensional spinor piece ($V_j^2=0.201002\\,\\mathbb 1$ to 5.8e-16, anticommutators 8.0e-16, one speed 0.448332) against the 2.000 and 0.444 V215 measured on the committed cells; 19 of 20 sampled rays carry exponent 1.0000 at both edges of the phase spread, the twentieth being a direction of the gapless set (64 of 384 bands pinned at phase zero at every radius) |
| `V226_node_charge_protection` | **no Berry charge at the zone-centre touching, and the invariant class masses it linearly**: the piercing is measured before any flux is read. On the mesh's 92 icosahedral directions no sampled point of the three spheres about the touching is gapless (their largest gaps 1.8e-03, 6.0e-03, 2.4e-02), while on the cell's six signed coordinate axes 4, 5 and 6 of 6 come back exactly gapless at those radii (worst line gap 2.6e-08) - the mesh is not sampling the lines the touching sits on. So every sphere is cut: 1 of 9 keeps full link rank and the refinement ladder reads 5, 12, 15, 16, 25, 34 at the same radii, a spread of 29, which is a mesh artefact and not a number. The instrument is right: $C=-1$ exactly on a two-level Weyl step, exactly 0 on spheres about a gapped point of this cell (defects 1.7e-15, 9.6e-16, 2.6e-15, all three radii trusted). The $2I$-invariant operator class is non-empty (exact group average over the 120 lifts, obstruction 1.4e-15, norm 7.31) and its on-site switch-on $S(q)\\mapsto e^{-i\\theta V}S(q)$ opens the touching linearly - 8.05e-06, 8.05e-05, 8.02e-04 at $\\theta = 10^{-3}, 10^{-2}, 10^{-1}$, gap over coupling constant to 0.3% - with the minimum relocating off the zone centre (3.6e-06 at radius 5e-03); the matched scalar-bond control is opened by the same class as well (2.03e-06, 2.14e-05, 3.24e-04) but at a rate that grows with the coupling rather than a constant one |
| `V227_locus_doubling_charges` | **the new cell has a gapless set, not a gapless point, no Berry charge on it, and no Brillouin zone to double over**: the grid over $\\lvert q\\rvert\\le\\pi/\\ell = 1.6516$ ($\\ell = 1.902113$, V214's generator length; 1791 points, spacing 0.22) puts its smallest gap at exactly 0 and counts 40 points below $10^{-8}$ (43 at $10^{-6}$, $10^{-3}$, $10^{-2}$); 12 components resolved at the grid's own resolution give **0 isolated point nodes and 12 points pierced by gapless rays** (10 to 18 of the mesh's 42 directions), transverse exponent 1.0000000; the gapless directions collect into **3 lines**, the coordinate axes, gapless along them to $10^{-12}$ (largest probe gap 2.6e-08). **0 of 12** enclosing spheres keep full link rank (smallest overlaps 2e-17 to 6e-16) and their readings pin at integer or half-integer defects, summed trusted charge 0; of the four spheres about the zone centre exactly one is trusted (flux 18 at $r=1.2387$, defect 1.4e-14), two pin at a half-integer and the zone-edge sphere refuses because the mesh crosses a gap closing - one trusted radius cannot establish radius independence, and the run says so. The doubling theorem's hypothesis is measured absent: the twelve directions in the basis of three give a $12\\times3$ matrix with entries in $\\mathbb Q(\\sqrt5)$, integrality search over $\\lvert n\\rvert\\le3$ bottoming out at defect 0.1459; the module's shortest rational hop is norm 2.0 along $(0,0,2)$ in 17 steps. The spinless control calibrates the instrument: 37 bands at phase zero at $\\Gamma$, its sphere at $r=0.040$ keeping rank (overlap 0.9998) and returning exactly 0 (defect 8.0e-16) |
| `V228_internal_sector_loop_return` | **the internal sector survives the spinor bond layer whole; the CPT operator and the loop return do not**: with the committed $32$-internal map the Clifford step is an internal scalar ($\\max\\lvert[S_{32}(q),\\mathbb 1\\otimes w]\\rvert = 8.4\\times10^{-31}$; $w^2=-1$ and $g_1wg_1^{-1}+w$ to exactly 0, $U_{32}^8-I=4.9\\times10^{-15}$); it is block diagonal in the pseudoscalar basis (off-block 2.9e-16 on 768 dimensions split $384+384$), its 16-block is the step built directly on $U_{16}$ (4.1e-15, 3.0e-15 - the reuse is forced, not chosen) and matter and antimatter are isospectral to 3.8e-15, 5.4e-15; the internal data are untouched ($5/3 = 1.6666666667$, the eight anomaly coefficients zero, $32 = 2\\times16$). CPT at the operator level splits on the bond layer against a matched control (same clock, layers, exponent, scalar bonds): $\\mathrm P g_1 S(q) g_1^{-1}(\\mathrm P g_1)^{-1}=S(-q)$ holds there at 5.0e-16, while on the Clifford cell **no candidate holds** - the three built without the inversion sit at 1.1348, the three with it at exactly $\\sqrt2 = 1.414214$, the very number V225 T1 measured as the broken inversion. The loop return $F_n=\\langle\\mathrm{tr}\\,S(q)^n\\rangle/d$ over the ball (grids 5, 9, 17; 2109 points; dimension 384) **splits by parity**: every odd $n\\le12$ vanishes in modulus below 9e-18, while the even ones alternate 0.460, -0.334, -0.294, 0.212, 0.209, -0.181 and stay resolution-limited under refinement (3.1e-02 rising to 0.23), against a scalar-bond control nonzero at every length - so the bond layer is what kills the odd returns, and the loop factor of Section 12 has to be rebuilt on this cell rather than carried over |
| `V230_quintet_helicity2_content` | **F4 content census, F4 stays fired** (zero fitted inputs): under the cubic group the quintet carrier $\mathrm{Sym}^2$ is $A_{1g} \oplus E_g \oplus T_{2g}$ (V208 control), and the qubit bond factor carries $T_{2g}$ but never $E_g$; under the icosahedral group of order 60 built from the presentation the SAME space is $A \oplus W$ with $W$ ONE irreducible 5-dimensional piece - exact $\mathbb{Q}(\sqrt5)$ character inner products, commutant ranks 2/1/1, the $E_g/T_{2g}$ split a cubic-subgroup artifact; the only MASSLESS LINEAR quintet branches live on the directional layer (axis bonds [0, 0, pi/4, pi/4, pi/2, pi/2] V140; diagonal star bonds [0, pi/2, pi] at pi/4, gap 0 at Gamma V229); at the closed clock the diagonal-bond shift gaps the quintet by up to 6.6047 and the axis-bond shift the $E_g$ half by 16.7287 with the $T_{2g}$ bubble exactly 0 |
| `V235_node_velocity_contact` | **F9 fires at the SME contact** (zero fitted inputs): node velocity tensors measured against the electron-sector ceiling $\delta v/v \le 8.0\times10^{-15}$ (Kostelecky-Mewes Table D6, fetched 2026-10-23, registered in the release `bounds`): $0.296$ at the implemented clock $1/5$, $1$ at the closed clock $\pi/4$ (one axis exactly flat), $0.295$ on the selective Wilson cones at $\theta^* = 1.7308$ (the node the manuscript names as the electron's), $0.9999$ on the star nodal line; the best of the named set, $0.292$, is a factor $3.7\times10^{13}$ past the ceiling - the constructed gapless nodes are refuted wherever an electron candidate is named |
| `V236_dressed_hop_isotropy` | **the dressing cannot restore isotropy at the closed clock** (zero fitted inputs): at $\pi/4$ the V141 dressing is degenerate as a mechanism - $A(h) = \rho(h^2) = -I$ per bond, the six insertions multiply to $+I$, dressed and bare spectra identical to exactly 0, the electric gate a scalar ($5\pi/4$, err 3e-14) - so it cannot move the cone: the dressed network's gapless points carry $\delta v/v = 1.0000$ (the $v = 1$ of the regular link is the front-speed maximum, measured to 1e-9); the star module's equivariant census starts at degree 5 - no closed-clock dressing moves the anisotropy |

### `gauged_layer`
| id | verdict |
|---|---|
| `V79_gauged_bond_layer` | the adjoint link layer exists and its Gauss constraint closes exactly |
| `V97_entanglement_bond` | the entanglement is the bond layer, not the shared cell |
| `V98_gauged_invariant_sector` | the centre obstruction: no one-body physical states |
| `V99_three_site_chain` | the 24-state kernel and the named Bekenstein gap |
| `V100_vertex_correction` | the vertex-corrected area-law coefficient ($I_6 = 141$) |
| `V101_bond_current_ground_state` | the first ground state of the gauged layer |
| `V129_gauged_layer_gap` | the layer is massive and its bond entropy is exactly saturated |
| `V130_string_tension` | the layer confines, at the Planck scale |
| `V140_dressed_gate_quintet` | the dressed gate on the directed-link network: the KS electric gate is exactly invariant, the dressed hopping is obstructed; the quintet disperses massless, linear, leading-order isotropic |
| `V141_regular_link_quintet_coupling` | the regular link of the root's finite group closes the hopping obstruction exactly and the dressed hop runs at v = 1; the quintet–matter coupling gaps the Eg half; the velocity trilemma is computed: no O_h-covariant step on the 6-bond module is isotropic at v = 1 |
| `V229_diagonal_bond_quintet` | **the named diagonal-bond gap closes and F8 stays fired** (zero fitted inputs): on the fifteen icosahedral star axes the $T_{2g}$ stress vertex that vanishes identically on axis bonds (V141 control 0.0) acquires 0.6325 on 6/15 per-axis entries and is nonzero on all fifteen axes; the quintet step on the star bonds is massless and linear with six gapless bands at Gamma, leading-order isotropy EXACT by the measured identity $\sum_a (\mathbf{n}\cdot\mathbf{d}_a)\hat{\mathbf{J}}_{\hat d_a} = 2\lvert d\rvert (\mathbf{n}\cdot\mathbf{J})$ to 9e-16; the helicity ladder doubles to $\{0, 2g_{10}, 4g_{10}\} \times \lvert d\rvert$, $\lvert d\rvert = 1.9021$, and no closed clock places any branch at $c = 1$ |

### `gravity`
| id | verdict |
|---|---|
| `V37_eh_action` | the discrete Einstein–Hilbert action derives the metric dynamics |
| `V50_complete_gravity` | Riemann from an arbitrary metric, two-route verification, the black hole by computation |
| `V51_numrel` | $3{+}1$ BSSN/Z4 evolution: the constraint rows agree to $O(h^2)$ |
| `V208_metric_bilinear_cone` | the metric is a bilinear of the bond layer: bond(6 oriented) $= A_{1g}\oplus E_g\oplus T_{1u}$ lacks the quintet's $T_{2g}$ and $\mathrm{Sym}^2 = A_{1g}\oplus E_g\oplus T_{2g}$ lacks the $T_{1u}$, so one bond amplitude cannot be the metric (exact characters); the commutant of the $36$-state bilinear has dimension $34$ — $34$ allowed self-couplings, none chosen; parity wall: the hop velocity is polar and its block inside the even (metric) sector vanishes ($2.7\times10^{-16}$), the axial dressing survives with spectrum $g\{0,\pm1,\pm2\}$ in every direction; the committed triple `cone_generator.json` is Clifford-regular (defect $1.4\times10^{-13}$) with an exactly direction-independent $36$-branch spectrum ($5.2\times10^{-15}$ over $44$ directions), a cone subspace of dim $8$ at slope $1$ with content $E\oplus T_1\oplus T_2$ (quintet weight $5.0000$), $4$ branches at slope $0$, $\mathrm{M}(0) = I$ exactly — no mass term writable in the class — and an exponentially local kernel; the covariant velocity space is $160$-dimensional and the exact-cone family through the generator $53$-dimensional (Jacobian rank $107$), and neither implemented triple lies in it (hop $0.0029$, dressing $0.0092$): the layer gives $c_{\rm gw}/c = 1.5708 = 2\pi/4$, excluded by GW170817 (**F8**), and the clock $g = 1/2$ that would fix it is in V206's dead zone ($n^* = 3$, deficit $1.4\times10^{-3}$); $G = 1/(4\ln 45) = 6.567\times10^{-2}$, cell $2\sqrt{\ln 45} = 3.902\,\ell_{\rm P}$, Bekenstein gap $4\ln 45 = 15.227$ |
| `V211_frame_stabiliser_gravity_gauge` | **what the root forces in the gravity/gauge sector, measured on the committed generators**: the $45$ bivectors close on $\mathfrak{so}(10)$ (defect $1.8\times10^{-15}$, Killing $64\,\mathbb 1$) with one trace form $\mathrm{Tr}(T_aT_b) = -32\delta_{ab}$; the background metric is the Clifford pairing of the causal span, $\mathrm{diag}(1,-1,-1,-1)$ exactly; the stabiliser of that plane is $\mathrm{spin}(4)\oplus\mathrm{spin}(6)$ ($21 = 6+15$) and the other $24$ generators form a $\mathbb Z_2$-graded coset ($10^{-15}$ both ways); $\mathfrak{su}(2)_L\oplus\mathfrak{su}(2)_R$ closes with residual exactly $0$; $\mathfrak{u}(3) = \mathfrak{su}(3)\oplus\mathfrak u(1)$ sits inside $\mathrm{spin}(6)$; the $B-L$ spectrum of the $\mathbf{16}$ is $(+\tfrac13)^6(-1)^2(-\tfrac13)^6(+1)^2$ and $\mathrm{Tr}\,Y^2/\mathrm{Tr}\,T_{3L}^2 = 5/3$ as a trace ratio. The boundary is recorded at V208's committed values (hop $0.29\%$, dressing $0.92\%$, a $53$-dimensional exact-cone family, $c_{\rm gw} = 1.570796$): content and couplings forced, the selected cone not |

### `scales_loopfactor`
| id | verdict |
|---|---|
| `V53_hier_exponent` | the ladder exponent $10$ and its algebraic reading |
| `V69_link_division_tower` | the coarse-graining tower is exactly kinematic |
| `V72_block_rg` | the exact real-space block-spin decimation of the bond layer |
| `V73_bond_gate_loophole` | the gate loophole closed: the factorized flow is the classical shadow |
| `V77_link_division` | one traced qubit link divides by exactly $\cos 2\theta$ |
| `V80_adjoint_link_division` | the adjoint link divides, but not to the loop factor (factor $26$) |
| `V83_link_rg_levels` | the traced-link division does not compound to the loop law |
| `V85_closed_loop_factor` | the closed-loop factor $F_{10} = 0.0192$ |
| `V89_loop_selection` | the selection rule measured to length $30$: the crossing is at $n = 10$ |
| `V94_composition_law` | the composition law: multiplicativity across legs, not along a chain |
| `V108_hierarchy_multiplicity` | the hierarchy coefficient as a product of measured root numbers |
| `V113_axis_factorization` | the axis-family factorization audit of that coefficient |
| `V126_loop_enumeration` | which loops supply the ladder factor |
| `V197_rung_arithmetic` | the rung ruler on every quoted scale: ladder rows integer-exact, every other scale off-integer by $0.06$–$0.42$ rungs, and the [$m_B$ at rung 10, $\sum m_\nu$ at $0.05$ eV] pair forces the seesaw partner to rung $2.72$ — mutually exclusive (F3) |
| `V200_loop_weight_functional` | the weight of the traced loop census is exactly $F_n = N_n\, r_n\, I_n$ (multiplicity, mean-word Cartan amplitude, interference), word norms angle-independent, and the O8 residual is the single deficit factor $0.9644$ at the root angle |
| `V202_loop_angle_selection` | the weight functional as an exact real cosine series of the angle (63 even phase classes, period $\pi$, $\eta \to \pi - \eta$ pairing, coefficient mirrors): the $1/(16\pi)$ crossings $0.1981004$, $0.5778509$, $1.0048608$, $1.3473566$ in $(0,\pi/2]$ - no root angle among them, $g_{10} = 1/5$ short by $1.9\times10^{-3}$ (branch ANGLE-OUTSIDE) |
| `V203_closed_branch_coefficient` | **the closed branch measured self-consistently** (step and traced-link phase both at $\pi/4$): the odd-$n$ selection rule is exact ($3.3\times10^{-15}$), the word norms stay angle-independent, the first crossing of $1/(16\pi)$ stays at $n^* = 208$ with deficit $0.9895$ - the published ladder length survives closure - and $F_{10}(\pi/4) = 55.03$, not the hybrid eta-law $0.1949$; at $\pi/4$ the 12 M channels are EXACTLY gapless with identical spectra $\{0^{32}, \pm\pi/3, \pm2\pi/3, \pi^{32}\}$ ($64/128$ pinned, $3\times10^{-15}$ across the 12), $\Gamma$ fully pinned $64+64$, and the exactly gapless set is $\Gamma$ plus the six face-diagonal lines with a gapped bulk ($\ge 1.2\times10^{-2}$) - the V132 line network (V65 re-gaps the M nodes under the Wilson factor); the reference branch's $99.44\%$ is identified exactly as a compensation $2.0881 = 3 \times 0.9644^{10}$; the per-channel factor at $n^*$ is $0.8997$ and the node-set readings $0.900\,/\,2.699\,/\,10.8$ for $\{1,3,12\}$ channels miss the required sum factor $2.308$ (closest $+30.0\%$) - the mechanism survives with zero free parameters, the precision is angle-generic |
| `V204_closed_ladder_limits_audit` | the closed-ladder concepts audited end-to-end at $\pi/4$: per-channel $(16\pi F_*)^{10} = 0.8997$, readings $0.900/2.699/10.8$ vs observed $2.0766$ and required $2.308$; rung 3 $= 9.6\times10^{13}$ GeV and rung 10 $= 118.57$ GeV (the $m_B$ the cell carries); bare rung-17 neutrino $0.146$ eV (the $1/3$ candidate $0.0487$); $\alpha_X$ root/observed $= 1.227$, $\tau_p = 2.769\times10^{39}$ yr; the $f_6$ out-of-equilibrium clause fails ($2.8\times$ observed rate, washout survives) |
| `V205_channel_filtered_deficit` | the O8 residual measured channel-blind: the required insertion amplitude $0.96388$ (deficit $3.61\%$) against the measured crossing deficit $0.98948$ ($1.05\%$) — the traced-word sector phases at the closed clocks are fourth roots of unity, every filtered propagator trace on the pinned channel spectrum is an empty kernel ($3\times10^{-14}$), the per-Cartan insertion amplitudes are bimodal $\{0.90, 0.99\}$; the constants of observation shift the observed coefficient by only $10^{-4}$ ($2.07656$); the required deficit matches the discarded reference branch's interference deficit to $0.055\%$ |
| `V206_clock_tension` | **O4 and O8 are one clock tension, measured**: the closure class reads the hierarchy coefficient at $2.7002$ ($\pi/4$) and $2.7084$ ($3\pi/4$) — both crossings on the $n = 4k$ class, $30.0\%$ / $30.4\%$ off the observed $2.0766$ — while the $\pi/2$ reading $2.2994$ ($10.7\%$ off) is floor-conditioned ($n^* = 102$ is a $4k{+}2$ length with no limit at the closure angle, V209) and not a stable number; the observed $2.0766$ is read by the same self-consistent functional at $g_{\rm obs} = 0.1999980$, the rational clock $1/5$ to $10^{-5}$, with the rung-10 crossing at $0.20013$; nature's clock is not closed, the exact ladder (odd-loop cancellation below $2\times10^{-15}$) lives only on the order-$8$ closed clock family, and F2 (deficit $0.9896$ predicted against $0.96388$ required) is the decision |
| `V209_root_clock_degeneracy` | **the root clock is phase-degenerate: the loop ratio there is exactly $0/0$** — $8g \in 2\pi\mathbb Z$ makes the $\mathbf{288}$ and $\mathbf{1120}$ traced-link phases coincide, and exact integer arithmetic at $\pi/4$ gives $\mathrm{Tr}(C_a N^n) = 0$ for every even $n \le 60$ and $\mathrm{Tr}(C_a W^n) = 0$ for every $n \le 60$; the period-4 classes resolve it differently — $n = 4k$ is removable (both norms vanish linearly in $g - \pi/4$; the float64 ladder is the limit to $4\times10^{-4}$, certified at 50 digits: $F_4 = 1.82995$, $F_8 = 1.67436$, $F_{208} = 0.0196928$, deficit $0.989868$, coefficient $2.7096$ at $\pi/4$ and $2.7113$ at $3\pi/4$, $n^* = 208$ at both detunings $10^{-6}$ and $10^{-8}$), while $n = 4k+2$ has no limit (the denominator vanishes to far higher order, $F_{210} = 3.7\times10^{7}$ at $\pi/4 + 10^{-8}$, poles accumulating at the closure angle) — so $F_2 = 65.7$, $F_6 = 60.1$, $F_{10} = 55.0$, $F_{170} = 1.57$, $(F_{10})^{17}$, $R_{17} = 4.0\times10^{-30}$, the "$33\times$ branch separation" and the $\pi/2$ crossing are ratios of round-off floors there, not measurements; the open branch needs no limit: $F_{10}(1/5) = 0.019186391$ stable at 15/30 digits, coefficient $2.0881$ |
| `V210_open_branch_composition` | **the composition law, measured where the functional is non-degenerate**: at $g = 1/5$ the whole ladder $F_n$, $n \le 172$, is stable between 30- and 50-digit evaluation; one chain of 170 hops on one line gives $F_{170} = 1.3089\times10^{-21}$ while 17 insertions on different lines multiply to $(F_{10})^{17} = 6.4697\times10^{-30}$ — a chain suppresses $2.02\times10^{8}$ times *less* than the product of its pieces, so the ladder's 17 factors cannot come from one geometric chain; the rung-17 suppression the ladder carries is $0.5401 \times (1/16\pi)^{17}$ (the clock deficit $0.964413$ to the 17th power); the odd-$n$ cancellation is a closure-class property, not a law ($F_3(1/5) = 0.242987$, matching V206's record) |
| `V218_closed_word_census_module` | **the loop factor on the icosahedral module, exactly where it can be computed**: closed words of length $n$ on the star's six axes against the cubic three - 12/6, 396/90, 19920/1860, 1281420/44730, 96807312/1172556 at $n = 2,4,6,8,10$, brute-force walk counts agreeing at every checked length - the count ratio running 10.7, 28.6, 82.6, 251.6 against $2^n = 64, 256, 1024$, the shortfall being the return-law prefactor. The dynamical half is resolution-limited and that is measured, not assumed: on the cubic step the trapezoid is the exact Fourier projector (grids 12, 16, 24, 32 agreeing to 4.8e-14, $F_n$ below the loop scale 0.019894 from $n = 3$), while on the star's almost-periodic box mean no length resolves at the 0.05 cut - relative grid drift 0.24 at $n = 2$ rising to 1.99 at $n = 10$ - because the six axes have irrational images in $\mathbb Q(\sqrt5)$, so no finite grid is a projector. The ratio $F_n$(star)/$F_n$(cubic) at the finest grid swings by a relative 1.21 over the unresolved lengths and is not quotable; resolving the box mean at $n = 208$ would need $832^3 = 5.8\times10^8$ diagonalisations of a 128-dimensional step |
| `V231_star_loop_functional_certified` | **the star-module box mean certifies - and supplies no loop factor** (zero fitted parameters, 1119 s): 416 exact closed-form box means $F_n(B)$ agree between dps 60 and 100 below 1e-40 at every $n \le 12$, both triples, boxes $\pi$-$64\pi$, four clocks, with the internal factor $\mathrm{tr}_{16} U^{6n} = 16(-i)^n$ certified by the exact-integer chain to 0e+00 (the first implementation's $\mathrm{tr}_{16} U^m = \mathrm{Re}\,\mathrm{tr}_{32} G^m$ was a factor-2 error, retired); NO certified crossing on the closed clocks - every $n \ge 2$ value sits above $1/(16\pi)$, closest approach 7.75x at $n = 12$; at $1/5$ a single length, $n = 6$ ($4k{+}2$), dips below, deficit ratios 1.4e-3-0.84, nearest value above 1.62x - a one-length dip, not a crossing; the grid ladder to $n = 208$ floors at 0.238 and certifies nothing beyond the closed form; $\pi/2$ structurally degenerate (46 values below 1e-12); exact-integer structure: signed-star relation lattice = six weight-2 opposite pairs (the draft's closing-triangle claim retracted), triple kernels empty, and the star's odd class NOT closed-word zero ($d_0 = (2c^2)^3$; min odd $\lvert F_n\rvert = 0.087$, $\lvert F_1(1/5,4\pi)\rvert = 113.4$) - the cubic odd-class zero is sublattice grading, not module combinatorics; the certified O8 number remains $F_{10}(1/5) = 0.019186$ (V209/V85) |
| `V232_insertion_amplitude_search` | **O8 stays open** (zero fitted inputs): the required insertion amplitude $0.96388$ is not supplied by any of 60 candidates inside the closed structure set only - the five unification Cartan directions (V205 bimodal 0.94-1.02), the certified crossing deficits 0.989868 ($\pi/4$) and 0.989933 ($3\pi/4$), the dressed electric-gate trace amplitudes at the closed clocks, the regular-link insertion (1 exactly) and the dressed hop, the 46 V231-certified star-module deficits, and the flagged reference-branch deficit 0.9644132 at $1/5$ (not a closed angle); closest closed-input candidate 9.8e-3 off, closest overall 5.3e-4 off and off-closure |
| `V233_clock_arbitration` | **the clock decision, F2 clock side FIRED** (zero fitted inputs): the 2x2 score over the clocks $\{\pi/4, 1/5\}$ against (i) the observed hierarchy coefficient $2.076564 \pm 10^{-4}$ and (ii) the fitted desert coupling $0.04$: $\pi/4$ reads 2.7095 certified (+30.5%, 6330x the uncertainty - the V209 number), $1/5$ reads 2.0881 (+0.55%, 115x) - NEITHER fits, and the claimed $1/5$ agreement is a FIT at $g_{\rm obs} = 0.1999980$ ($2.0\times10^{-6}$ off the fifth); $\pi/4$'s standard coupling reading $\pi/64 = 0.049087$ is $1.227\times$ the fitted value and $1/5$'s $0.0400$ matches only by the fit; the clock moves from the algebraic to the fitted class, and no closed clock passes either empirical contact |
| `V234_rung_phase_spectrum` | **the rungs are declared, not solved** (zero fitted inputs): on the 36-level cell step only $\mathrm{loop}^3$ is numerically realized (matched to 1e-16); $\mathrm{loop}^{10} = 9.71\times10^{-18}$ sits below the float64 phase resolution and is certified only by the exact ladder construction; no other integer $0..30$ matches any level phase; the DM window pins $\mathrm{rung}_B = 10$ uniquely, and rung 10 admits NO rung_A whose identity $m_B^2/m_A$ lands inside the measured neutrino contacts (rung_A 3 gives 0.146 eV - the published EXCLUDED ladder; rung_A 2 gives 0.0029 eV, and 2.7265 demanded by the atmospheric scale is not realizable by integer rungs); the pair $(3, 10)$ is a config declaration; F3 recorded FIRED |

### `couplings_unification`
| id | verdict |
|---|---|
| `V38_couplings_root` | the readings $g_{10} = 1/5$ and $g_0 = 1/8$ are root-data conventions (the clock closed at $\pi/4$); the running is the named residual |
| `V56_g10_selection` | the one-loop data check of the $1/5$ reading (right order, factor 1.73) |
| `V64_desert` | the desert solution: $\alpha_X = 1/25$, $M_X = 1.8\times10^{16}$ GeV, $\tau_p$ |
| `V74_multiscale_nodes` | multi-scale node thresholds cannot move $M_X$ (rank-1 Jacobian) |
| `V81_link_desert_content` | the composite spectrum does not deliver the meeting |
| `V82_composite_desert_exclusion` | the composite content excluded under every natural mass ansatz |
| `V84_equal_shift_theorem` | the equal-shift theorem over complete multiplets |
| `V86_junction_desert` | the fragment pattern is fabricated, the meeting is not |
| `V87_general_junction` | the general non-democratic junction scanned as the supplier |
| `V92_unification_feasibility` | the data-preferred meeting lies inside the node-$\mathbf{16}$ fragment hull |
| `V93_cartan_junction` | a general Cartan direction realizes the meeting but breaks the Model algebra |
| `V95_junction_split` | no exact meeting for a universal junction strength |
| `V96_kappa_s_bond` | the flavour-dependent bond layer closes the meeting exactly, at six named inputs |
| `V105_cartan_kappas` | the six $\kappa_s$ are not obtainable from any on-site algebra element |
| `V110_offdiagonal_bond` | off-diagonal bond realizability: the meeting floor |
| `V122_class_invariant_kappas` | class invariants as the source of the $\kappa_s$: excluded |
| `V133_desert_pi64` | the desert re-closure at $\alpha_X = \pi/64$ leaves the desert |
| `V137_broken_flux` | the broken-generator flux is a per-class bond number, but not the one that closes |
| `V145_superexchange_kappas` | the second-order bond-current (superexchange) shift is a class-constant $\tilde\kappa_s$ but misses the meeting: residual 0.577, 0/10 desert closures |
| `V146_two_number_kappas` | superexchange + broken flux + ground-state Cartan shift solved jointly over all node counts: residual 3.46, $R^2 \le 0.30$ — the supplier space is exhausted |
| `V237_desert_threshold_matching` | **the root-coupling solution is not realizable by its own content, F2 unification side fires** (zero fitted inputs): at $\alpha_X = \pi/64$, $M_X = 1.235\times10^{17}$ GeV the three $\mathbf{10} \oplus \overline{\mathbf{10}}$ thresholds sit at $7.052\times10^4$ GeV - below the $10^6$ GeV desert floor - and the $(1,2)$ threshold at $4.877\times10^{17}$ GeV, above $M_X$ itself; restoring the ordering moves the $(1,2)$ threshold $+13.1$ e-folds, $1.12\times$ the full log-spacing from $1.0\times10^{12}$ GeV to $M_X$ - no one-loop threshold within the spectrum's own spacings restores the meeting; the driving beta-content is strictly integer-theoretic ($k_1 = 3/5$, $\Delta b = (6/5, 2, 2)$) |
| `V238_desert_spectrum_bounds` | **no desert threshold is excluded today - F10 open** (zero fitted inputs): the colour pair at $M_{10} = 7.060\times10^4$ GeV sits $14.1\times$ above the strongest LHC colour probe (5 TeV, CMS $t$-channel dimuon + dielectron, fetched 2026-10-23) and pair production is kinematically closed at the LHC ($2M_{10} = 141.2$ TeV > 13.6 TeV); the $(1,2)$ doublet at $1.0\times10^{12}$ GeV is $10^9\times$ above the electroweak-precision doublet ceiling (oblique $STU$, PDG 2025) and fully decoupled ($(v/M)^2 = 6\times10^{-20}$); $\tau_p = 2.76\times10^{39}$ yr is $1.2\times10^5$ past the Super-Kamiokande bound - the spectrum's failure (V237) is internal ordering, not exclusion |

### `neutrino`
| id | verdict |
|---|---|
| `V57_neutrino` | the two-step seesaw on the ladder |
| `V68_neutrino_node_audit` | the $n_{\text{nodes}} \ge 3$ bound and the democratic node-mixing rival |
| `V76_mixing_pattern` | the cubic-symmetric junction that closes the node triangle |
| `V90_neutrino_protection` | the protection audit: the cubic group induces all of $S_3$, $r$ is not quantized |
| `V103_democratic_point` | with the implemented democratic Dirac matrix the two massless states are protected |
| `V114_junction_r` | the junction stabilizer: $r$ is a convention, not a tuning |
| `V136_majorana_word` | the $\Delta L = 2$ operator is not a gauge-link composite |

### `cosmology_darkmatter`
| id | verdict |
|---|---|
| `V54_dm_wimp` | the $118$ GeV ladder mass, the same factor as the hierarchy |
| `V66_dm_identity` | the carrier scan: no viable state of the spinor module (the carrier is the cell singlet B, V196) |
| `V75_dm_nogo_extension` | the no-go completed, plus the minimal fix |
| `V88_mirror_scalar` | the ladder scalar is stable under either portal |
| `V112_dark_u1_portal` | the dark-$\mathrm{U}(1)$ portal window is empty |
| `V123_relic_abundance` | the relic abundance computed from the couplings the root supplies |
| `V124_z2_gradings` | the exact $\mathbb{Z}_2$ gradings of the step: the mirror is the pseudoscalar |
| `V147_relic_from_bond` | every bond number of the scalar class as annihilation strength: $\Omega h^2$ lands at $10^{-12}$–$10^{-3}$ with $D$, never $0.120$ |
| `V148_direct_detection` | the same numbers as Higgs-portal strength: $\sigma_{\rm SI} = 2.5\times10^{-45}$–$6.3\times10^{-41}$ cm$^2$, all above the LZ-class bound — the bond-strength portal is excluded |
| `V193_freeze_in` | the production side computed — freeze-in and freeze-out with exact thermal convolution, the regime decided at the self-consistency root: the exact-$\mathbb{Z}_2$ portal's only solution inside the LZ-class ceiling is freeze-in at $\lambda = 4.1\times10^{-14}$ ($\sigma_{\rm SI} = 6\times10^{-71}$ cm$^2$, DD-blind); its thermal solution $\lambda = 0.31$ is excluded as published; the pair portal reproduces $0.120$ at its computed coupling $\lambda = 1.3\times10^{-3}$ ($\sigma_{\rm SI} = 6.7\times10^{-50}$) |
| `V194_carrier_in_cell` | the carrier census of the built cell: the exact grading group of step+gauge is $\mathbb{Z}_2\times\mathbb{Z}_2$ (global $\times$ mirror, no fragment-selective parity), and no built root places a neutral spinor-module fragment in the WIMP window (best $2.0\times10^{10}$ GeV) — the certificate the engine cell meets |
| `V195_cell_redefinition` | the certificate realized on the engine's cell object (`cl10/cell.py`): the cell content is the 36-level module — two neutral singlet sectors (four levels), $A$ at the $n=3$ seesaw floor ($m_\nu = m_B^2/M_A$, the F3 ladder intact) and $B$ odd at the $n=10$ floor — the carrier is a state of the automaton: exactly stable, neutral, at the DERIVED $118.57$ GeV, its $0.120$ by freeze-in at $\lambda = 4.1\times10^{-14}$ (DD-blind); the residual is the portal strength |
| `V199_a_to_bb_channel` | the $A \to BB$ channel with the $13$ bond numbers of the singlet class as coupling candidates: $9$ of $13$ drive $\mathrm{Br}(A \to LH)$ below $6.4\%$ and kill the resonant leptogenesis route; the route survives only for $\kappa < 3.04$, the channel is negligible below $0.080$ (branch SPLIT), and no value produces the relic |

### `baryogenesis`
| id | verdict |
|---|---|
| `V55_eta_cp` | the asymmetry as a product of factors of three provenances |
| `V106_washout_audit` | at the neutrino Yukawa the decay is in equilibrium: the mechanism fails |
| `V109_cp_phase_search` | mirror-preserving extensions supply no CP phase; class-dependent dressings break the mirror |
| `V117_real_form` | the module is of real (Majorana) type: that is what forbids the phase |
| `V125_cp_flux_washout` | the CP flux spent on the asymmetry |

### `family_number`
| id | verdict |
|---|---|
| `V26_family` | the three generations as a 0-parameter $3$, three ways |
| `V70_family_index` | the 12 nodes split into three axis families of 4 |
| `V91_count3_forcing` | the count $3$ is a symmetry-breaking input, not a forced number |
| `V116_triangle_orbit8` | the cubic orbits of node triangles; the 3-node set breaks the cubic symmetry |
| `V127_triangle_energies` | the energies of the eight triangle vacua: the degeneracy is not the junction's |

### `phenomenology`
| id | verdict |
|---|---|
| `V62_hydrogen_production` | the $296^3$ hydrogen run: resonant absorption, the $1s \to 2p$ deformation (Fig. 2) |

### `alt_gauged` (exploration programme — the two-fermion S-matrix of the gauged bond layer)
| id | verdict |
|---|---|
| `V175_gauged_smatrix` | THE TWO-FERMION S-MATRIX OF THE GAUGED LAYER EXISTS: S = exp(-i g_10 L^-1 H L^-dagger) on the 24-state three-site kernel, unitary to 2.2e-13, spectrum = V101 exactly (0.0e+00), fusion channels preserved (mixing up to 6.1e-29), eta->0 gives S = I |
| `V176_gauged_vertex_content` | the gauged current vertex is NOT a scalar: Eg + T2g weight 5/6 at |k| = 0.4 (vs contact A1g = 1.0000), anisotropic part scale-independent (alpha = 0 - MARGINAL, vs contact alpha = 2); on-node current F(0) = 5.000000 |
| `V177_gauged_self_energy_flow` | one-loop Fock flow with the gauged vertex: Sigma nonzero (0.204), anisotropic part ~ |k|^1.00 = MARGINAL (propagator sets the power, not the vertex); dressed A not grid-convergent (drift 1.2e-2, flat-axis IR — the contact baseline of the same run has the same non-convergence) |
| `V178_gauged_block_rg` | block RG b = 1, 2, 4 with the gauged insertion: A(b) = 0.997985 EXACTLY unchanged (bare = gauged at every scale) - anisotropy stays marginal (2-route 1.3e-12) |
| `V179_gauged_isotropy_gate` | 12 census rows, 0 candidates (2 unconverged V177 rows excluded, min A 0.9322): the gauged bond-current interaction does not isotropize the cone at one loop / under coarse-graining - the flux sector is massive (gap 64 eta, V129) and decouples from the IR cone |

### `alt_nonpert` (exploration programme — exact cone flux dressing, ring kernels, IR decoupling, non-perturbative isotropy gate)
| id | verdict |
|---|---|
| `V180_cone_flux_dressing` | exact cone dressing by the 24 gauged flux operators: A_np - A_free = +4.43e-12 (2-route eps dev 1.7e-16, eta->0 dev 0.0) - the flux dressing does NOT move the cone anisotropy at one dressing order |
| `V181_ring_kernel_smatrix` | the 4-site gauged kernels exist: chain 144 (full rank, E0 -1664197.02, gap 12787.60 eta) and ring 1296 (Gram rank 976, E0 -6.08e8, gap 2.75e8 eta); S unitary to 2.7e-15/4.2e-15, 2-route phase dev within the whitening-noise bound; chain channels preserved (mixing 3e-11), ring channels mixed (~8) |
| `V182_ir_decoupling` | INCONCLUSIVE under the pre-registered clauses: gap grows 42.09 -> 12787.60 (chain) -> 2.75e8 (ring) eta, drift >> 1 - the LIGHT branch is decisively false; V101 eta-units are not extensive (Gram shrinks ~30x/link), gap/\|E0\| ~ 0.0073-0.0077 at every volume |
| `V183_nonpert_isotropy_gate` | 9 census rows, 0 candidates: no non-perturbative gauged-sector mechanism isotropizes the cone - the flux sector is massive at every volume probed; the cone anisotropy (O9) stays open |

### `alt_step` (exploration programme)
| id | verdict |
|---|---|
| `V191_parity_dirac_cone` | **why** every node found so far is flat, and the node that is not: `Hom_O_h(T1u, D (x) D*) = 0` for all ten irreps (and dim 0 by direct solve on the T1u block), so no single-parity multiplet at an inversion-invariant momentum can be conical; the parity-degenerate `A1g+T1u` pair escapes it - constructed on the O_h-invariant hopping `-24 O14 + 8 O50`, gap 2.9e-11 at Gamma, multiplicity 4, cone anisotropy 1.0e-15, rho_comm = 1/sqrt 2, `TT^dag = 16 I` to 5.6e-16; the implemented 6-stage product is exactly O_h-equivariant at Gamma (1.1e-16) and breaks it at every other momentum (1.07e-1 generic, 1.40e-1 M, 1.46e-1 R), where its stage sum is invariant to 1.1e-15 |
| `V192_engine_two_angle_cone` | **the search for the cone inside the engine's own circuit - negative**: over the two-angle axis-symmetric stage family (`a` on the three intra-block stages, `b` on the three inter-block ones - the largest six-stage family keeping Gamma exactly covariant: 1.1e-16 at Gamma, 9.0e-1 at M, 6.3e-1 at generic q; a generic six-angle product misses Gamma covariance by 6.0e-2 best of 200 draws, a per-axis generalisation by 0.90), searched jointly with the internal angle over 33800 points of (a, b, theta) and polished to residuals below 4.3e-14, 60 points carry on-site states at phase zero, the full 128x128 step gap at most 4.4e-15. They split by parity of the gapless eigenspace, and the split is the answer: 8 points have a 4-dimensional single-parity gapless space (`A1g (+) T2g` or `A2u (+) T1u`, content read off the eigenspace by character inner product) and are flat to 1.3e-12 over 179 sampled directions, which is what Proposition 6 forces; the other 52 have the 8-dimensional parity-mixed space `A1g (+) T2g (+) A2u (+) T1u`, disperse (slope 0.09-1.73 per unit block momentum, 0.03-0.58 Planck cells/stage, inside the front speed 0.67, no amplitude outside the light cone in 12 stages) but are exactly Abelian, rho_comm <= 1.2e-11 at anisotropy 0.423 - nodal planes, not a cone. Conical nodes: 0. On the diagonal a = b, at pi/4 and 3pi/4: 64 gapless states carrying 16 x (A1g (+) T2g), gap at most 8e-17, slope at most 3e-14 - the flat zone-centre band of Section 8.3. The gapless subspace is taken from the eigenvectors of the full step, the internal phases from $U^6$ itself, and the $O_h$ content from the measured character table - a root-finder on a wrapped phase difference, an assumed two-level subspace, a small-angle extrapolation of the internal phases or an unverified character table each put a node where there is none |
| `V201_root_swap_cone` | **the cone search with the root's own sector swaps - negative**: the internal layer interleaved with each of the 10 involutions $S_c = h_c w$ of the Clifford module ($w = h_1..h_{10}$: $S^2 = I$, $SS^\dagger = I$, $\{S,w\} = 0$, $S = S^\dagger$, all exact to $0$, exchanging $\mathbf{16}$ and $\overline{\mathbf{16}}$ on the doubled module, anchored to the V192 Weyl step by $B^\dagger U_{32}^6 B = U_{16}^6$ to $8\times10^{-15}$, Gamma equivariance intact to $6.2\times10^{-17}$); the V192 grid rerun per swap (33800 points): $8$ exactly polished four-cancellation points each, all a $6$-dim $T_{2g}\oplus T_{1u}$ gapless multiplet with non-commuting tangent operators ($\rho_{\mathrm{comm}}$ 0.50–0.80 - the V192 Abelianness is removed) and vanishing linear velocity (max slope $1.7\times10^{-12}$): the non-commutativity is bought from the root, the dispersion is not, conical nodes 0 |

