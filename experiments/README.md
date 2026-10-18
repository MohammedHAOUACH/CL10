# The verification program

154 experiments (the 115 of the publication programme plus the 39 exploration
experiments of the `alt_*` programmes: 12 `alt_cell`, 10 `alt_dynamics`,
4 `alt_interactions`, 5 `alt_gauged`, 4 `alt_nonpert`, 4 `alt_step`), grouped by
the part of `doc/final_publication.md` they validate.  Each experiment
recomputes a published number from the construction, states the provenance of
every input it uses, and writes a JSON verdict (`pass`, `status`, `tier`,
`property`, `residual`, `wall_time`) to `<experiment>/results/`.  `python3
run_validation.py` runs all of them and rolls the verdicts up into
`results/validation_summary.json`.

The papers these experiments validate: [Paper I](../doc/paper1_core.md) ·
[Paper II](../doc/paper2_gauge_gravity.md) · [Paper III](../doc/paper3_phenomenology.md) ·
[the unified master](../doc/final_publication.md) — compiled PDFs in `../doc/`, on the
[v1.0.2 release](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0.2) and on
[Zenodo](https://zenodo.org/records/23321228); cite:
[10.5281/zenodo.23321106](https://doi.org/10.5281/zenodo.23321106).
Back to the [CL10 README](../README.md).

A script never hard-codes its own category: sibling experiments and their
verdicts are addressed by name through `cl10.paths.exp`, so moving an
experiment between categories moves a directory and nothing else.

The old wall labels map to the paper's open problems as W1→O1, W2→O2, W3→O3,
W4→O4, W5→O5, W6→O6, W7→O7, W8→O8, W9→O9; every verdict in this tree is
written in the current wording, and each experiment rewrites its own verdict
file whenever it is re-run.

| category | what it establishes | paper |
|---|---|---|
| `algebra_root` | the root itself, the closure of the coupling, the provenance audit | §2, §9, §16 |
| `spacetime_matter` | the $(1{+}3)$ subalgebra, one anomaly-free generation, the $45$, the statistics, the spectral mirror, the flavour structure | §3–§6, §13.5 |
| `dynamics_spectrum` | unitarity, the causal radius, the $n$-particle sector, the doubling and its Wilson gap, the node geometry, the cone-isotropy no-gos, the graviton content | §7, §8, §9 |
| `gauged_layer` | the adjoint bond layer: Gauss constraint, centre obstruction, area law, vertex correction, ground state, confinement | §10 |
| `gravity` | the discrete Einstein–Hilbert action, curvature from an arbitrary metric, $3{+}1$ evolution | §11 |
| `scales_loopfactor` | where a loop factor can come from, the closed-loop factor, the ladder length, the composition law | §12.1 |
| `couplings_unification` | $\alpha_{10} = 1/25$, the desert solution, the proton lifetime, every proposed supplier of the desert spectrum | §13.1, §13.4 |
| `neutrino` | the seesaw ladder, the $\tfrac13$ candidate, the protection of the two massless states, the $\Delta L = 2$ no-go | §12.3 |
| `cosmology_darkmatter` | the $118$ GeV ladder mass, the carrier no-go, the built scalar, the relic abundance, the portal window | §12.2, §14.2 |
| `baryogenesis` | the absence of every Sakharov condition, the asymmetry product, the washout, the real-form obstruction | §13.3 |
| `family_number` | why the count $3$ is an input and what its symmetric alternative costs | §4, §17 (O3) |
| `phenomenology` | the lattice photon driven through hydrogen (Fig. 2) | §6 |
| `alt_cell` | the exploration programme: the Cl(p,q) signature census and the alternative-cell census against the node/isotropy gates (tier 4) | §17 (O1, O9) |
| `alt_dynamics` | the exploration programme: the alternative step-rule (motion) census - per-stage angles, internal factorizations, mass stages, bond matchings, Floquet, quantized Wilson - against the same gates (tier 4) | §17 (O1, O5, O9, W4, W7) |
| `alt_interactions` | the exploration programme: do interactions (the minimal admissible block-contact collision) restore cone isotropy at low energy - vertex irrelevance, one-loop flow, block RG (tier 4) | §17 (W9) |
| `alt_gauged` | the exploration programme: the two-fermion S-matrix of the gauged bond layer (V79-V101) built and measured - scattering phases, vertex O_h content, one-loop flow, block RG - against the W9 isotropy gate (tier 4) | §17 (W9) |
| `alt_nonpert` | the exploration programme: non-perturbative gauged-sector probes of W9 - exact cone dressing by the flux operators, 4-site chain/ring kernels and S-matrices, IR decoupling across volumes (tier 4) | §17 (W9) |
| `alt_step` | the exploration programme: the step/cell SPECIFICATION searched directly - a differentiable cone loss over the 45 continuous per-plane step angles (Adam, JAX) and a joint cell-shape census - against the W9 isotropy gate and the hard tier (tier 4) | §17 (W9) |

The `alt_*` programmes run on the **numpy core**: each script pins
`CL10_BACKEND=numpy` before importing `cl10`, and `cl10/altcell.py` is
host-only by contract (every call into an engine-routed module is moved back
to the host where it enters). Their linear algebra is host-side, so their
verdicts do not depend on what else holds the card; the runner lists them
OFF-CORE. The GPU engine stays with the publication programme.

## The experiments

### `algebra_root`
| id | verdict |
|---|---|
| `V0_root` | $\mathrm{Cl}(10,0)$ is the root: the ten generators, $w^2=-1$, the algebra basis |
| `V52_parameter_audit` | the provenance audit: algebraic / conventional / fitted / open, and the CP-odd invariants that vanish |
| `V59_four_limits` | the four limits reduce to the seven algebra invariants plus one angle |
| `V78_w1_lambda` | $\Lambda = 0$ is protected by $\det U = 1$ under arbitrary traceless perturbations |
| `V131_algebraic_closure` | the root quantizes its own step: $g \in \{\pi/2,\pi/4\}$, the order-$8$ spectrum |
| `V134_pi2_branch` | the $\pi/2$ branch is a pure lattice translation |
| `V135_coupling_selection_closed` | the closed coupling $\pi/64$ against the one-loop data |

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
| `V43_qca_spectrum` | one 32-level spectrum, nine named projections |
| `V28_yukawa` | the democratic rank-1 Dirac mass matrix |
| `V45_yukawa_matrix` | the $3\times3$ Yukawa matrix, the TBM mixing, $\sin^2\theta_{12} = 1/3$ |
| `V49_quark_ckm` | the computed CKM is the identity |
| `V67_breaking_chain` | the Wilson direction preserves the Model algebra and breaks to $\mathrm{SU}(3)\times\mathrm{SU}(2)\times\mathrm{U}(1)^3$ |
| `V104_w6_mixing_manifold` | the natural vev placements are excluded or vacuous |
| `V107_mixing_from_kappas` | gauge-class bond couplings are generation-blind: they cannot supply mixing |
| `V111_w6_generation_bond` | the parameter cost of a generation-space dressing (7 / 13) |

### `dynamics_spectrum`
| id | verdict |
|---|---|
| `V36_local_qca` | the causal radius is 1 cell per stage; the light cone is realized |
| `V41_smatrix` | the step is the $S$-matrix; the first processes and the cross-section bound |
| `V48_3particle` | the $3$- and $4$-particle amplitudes are UV finite; bounded for all $n$ |
| `V63_bz_chirality` | chirality is a dynamical symmetry; the 12 gapless node momenta, $\sum$ chirality $=0$ |
| `V65_wilson` | an internal Wilson factor of the same step gaps all 12 doublers |
| `V102_lorentz_isotropy` | the node dispersion: 1200 sampled directions, anisotropic linear cone |
| `V115_w9_cone_isotropy` | no isotropy point anywhere on the Wilson family |
| `V118_w9_little_group` | tangent velocities and node little groups: three nodal planes, not a cone |
| `V120_w9_entangling_bond` | the entangling bond layer reaches $0.071$ of the required commutator |
| `V121_w9_three_axis_bond` | the three-axis covariant bond reaches $0.658$ and costs the causal front |
| `V128_f4_graviton_content` | no massless spin-2 content in the implemented spectrum (**F4**) |
| `V132_gapless_census_pi4` | at $\pi/4$: $\Gamma$ plus six face-diagonal nodal lines |
| `V138_gamma_spin2` | the zone centre carries $T_{2g}$ but never $E_g$: no lattice spin-2 in the qubit link, the quintet lives in the directional link layer |
| `V139_nodal_plane_crossing` | the grid census was too coarse: an exactly gapless orbit of $24$ points on the nodal planes at $g_{10} = 1/5$, little group of order $2$ |
| `V142_w9_frame_averaging` | frame averaging over the ten generators is exactly blind to the cone anisotropy: in four random SO(10) frames the node gap moves by 1.8e-15, the anisotropy by 0, the tangent commutators stay at the noise floor — all momentum dependence lives in the sub-site factor of the Kronecker bond layer |
| `V185_separation_lemma` | the separation lemma: every tangent operator is exactly $A_k(q) \otimes I_{16}$ (defect 4.5e-16 over 40 momenta), so the whole low-energy velocity structure lives in the 8×8 sub-site factor |
| `V186_cone_inverse` | the cone inverse problem over the axis-covariant entangling family (880 bonds, 220 triples): no cone in the family; the implemented corner triplet is exactly radius-1 with the Wilson factor off (leakage 8.7e-17) |
| `V187_rigidity_frontier` | 5292 sweep points around the top-commuting bonds: gaplessness and non-commuting tangent algebras are mutually exclusive, the frontier floor is a wall at gap 8.0e-3, no power law links the regimes |

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

### `gravity`
| id | verdict |
|---|---|
| `V37_eh_action` | the discrete Einstein–Hilbert action derives the metric dynamics |
| `V50_complete_gravity` | Riemann from an arbitrary metric, two-route verification, the black hole by computation |
| `V51_numrel` | $3{+}1$ BSSN/Z4 evolution: the constraint rows agree to $O(h^2)$ |

### `scales_loopfactor`
| id | verdict |
|---|---|
| `V53_hier_exponent` | the ladder exponent $10$ and its algebraic reading |
| `V69_w8_tower` | the coarse-graining tower is exactly kinematic |
| `V72_block_rg` | the exact real-space block-spin decimation of the bond layer |
| `V73_bond_gate_loophole` | the gate loophole closed: the factorized flow is the classical shadow |
| `V77_link_division` | one traced qubit link divides by exactly $\cos 2\theta$ |
| `V80_adjoint_link_division` | the adjoint link divides, but not to the loop factor (factor $26$) |
| `V83_link_rg_levels` | the traced-link division does not compound to the loop law |
| `V85_closed_loop_factor` | the closed-loop factor $F_{10} = 0.0192$ |
| `V89_loop_selection` | the selection rule measured to length $30$: the crossing is at $n = 10$ |
| `V94_composition_law` | the composition law: multiplicativity across legs, not along a chain |
| `V108_hierarchy_multiplicity` | the hierarchy coefficient as a product of measured root numbers |
| `V113_w8_axis_factorization` | the axis-family factorization audit of that coefficient |
| `V126_w8_loop_enumeration` | which loops supply the ladder factor |

### `couplings_unification`
| id | verdict |
|---|---|
| `V38_couplings_root` | $g_{10} = 1/5$ and $g_0$ are numbers of the root; the running is the named residual |
| `V56_g10_selection` | the one-loop data check of $g_{10} = 1/5$ |
| `V64_desert` | the desert solution: $\alpha_X = 1/25$, $M_X = 1.8\times10^{16}$ GeV, $\tau_p$ |
| `V74_multiscale_nodes` | multi-scale node thresholds cannot move $M_X$ (rank-1 Jacobian) |
| `V81_link_desert_content` | the composite spectrum does not deliver the meeting |
| `V82_composite_desert_exclusion` | the composite content excluded under every natural mass ansatz |
| `V84_equal_shift_theorem` | the equal-shift theorem over complete multiplets |
| `V86_junction_desert` | the fragment pattern is fabricated, the meeting is not |
| `V87_general_junction` | the general non-democratic junction scanned as the supplier |
| `V92_w4_feasibility` | the data-preferred meeting lies inside the node-$\mathbf{16}$ fragment hull |
| `V93_cartan_junction` | a general Cartan direction realizes the meeting but breaks the Model algebra |
| `V95_junction_split` | no exact meeting for a universal junction strength |
| `V96_kappa_s_bond` | the flavour-dependent bond layer closes the meeting exactly, at six named inputs |
| `V105_w4_cartan_kappas` | the six $\kappa_s$ are not obtainable from any on-site algebra element |
| `V110_w4_offdiagonal_bond` | off-diagonal bond realizability: the meeting floor |
| `V122_w4_class_invariant_kappas` | class invariants as the source of the $\kappa_s$: excluded |
| `V133_desert_pi64` | the desert re-closure at $\alpha_X = \pi/64$ leaves the desert |
| `V137_w4_broken_flux` | the broken-generator flux is a per-class bond number, but not the one that closes |
| `V145_w4_superexchange_kappas` | the second-order bond-current (superexchange) shift is a class-constant $\tilde\kappa_s$ but misses the meeting: residual 0.577, 0/10 desert closures |
| `V146_w4_two_number_kappas` | superexchange + broken flux + ground-state Cartan shift solved jointly over all node counts: residual 3.46, $R^2 \le 0.30$ — the supplier space is exhausted |

### `neutrino`
| id | verdict |
|---|---|
| `V57_neutrino` | the two-step seesaw on the ladder |
| `V68_w2_audit` | the $n_{\text{nodes}} \ge 3$ bound and the democratic node-mixing rival |
| `V76_w2_mixing` | the cubic-symmetric junction that closes the node triangle |
| `V90_w2_protection` | the protection audit: the cubic group induces all of $S_3$, $r$ is not quantized |
| `V103_w2_democratic` | with the implemented democratic Dirac matrix the two massless states are protected |
| `V114_w2_junction_r` | the junction stabilizer: $r$ is a convention, not a tuning |
| `V136_w2_majorana_word` | the $\Delta L = 2$ operator is not a gauge-link composite |

### `cosmology_darkmatter`
| id | verdict |
|---|---|
| `V54_dm_wimp` | the $118$ GeV ladder mass, the same factor as the hierarchy |
| `V66_dm_identity` | the carrier scan: no viable state in the implemented spectrum |
| `V75_dm_nogo_extension` | the no-go completed, plus the minimal fix |
| `V88_mirror_scalar` | the ladder scalar is stable under either portal |
| `V112_w7_dark_u1_portal` | the dark-$\mathrm{U}(1)$ portal window is empty |
| `V123_w7_relic_abundance` | the relic abundance computed from the couplings the root supplies |
| `V124_w7_z2_gradings` | the exact $\mathbb{Z}_2$ gradings of the step: the mirror is the pseudoscalar |
| `V147_w7_relic_from_bond` | every bond number of the scalar class as annihilation strength: $\Omega h^2$ lands at $10^{-12}$–$10^{-3}$ with $D$, never $0.120$ |
| `V148_w7_direct_detection` | the same numbers as Higgs-portal strength: $\sigma_{\rm SI} = 2.5\times10^{-45}$–$6.3\times10^{-41}$ cm$^2$, all above the LZ-class bound — the bond-strength portal is excluded |

### `baryogenesis`
| id | verdict |
|---|---|
| `V55_eta_cp` | the asymmetry as a product of factors of three provenances |
| `V106_w5_washout_audit` | at the neutrino Yukawa the decay is in equilibrium: the mechanism fails |
| `V109_w5_cp_search` | mirror-preserving extensions supply no CP phase; class-dependent dressings break the mirror |
| `V117_w5_real_form` | the module is of real (Majorana) type: that is what forbids the phase |
| `V125_w5_cp_flux_washout` | the CP flux spent on the asymmetry |

### `family_number`
| id | verdict |
|---|---|
| `V26_family` | the three generations as a 0-parameter $3$, three ways |
| `V70_w3_index` | the 12 nodes split into three axis families of 4 |
| `V91_w3_forcing` | the count $3$ is a symmetry-breaking input, not a forced number |
| `V116_w3_triangle_orbit8` | the cubic orbits of node triangles; the 3-node set breaks the cubic symmetry |
| `V127_w3_triangle_energies` | the energies of the eight triangle vacua: the degeneracy is not the junction's |

### `phenomenology`
| id | verdict |
|---|---|
| `V62_hydrogen_production` | the $296^3$ hydrogen run: resonant absorption, the $1s \to 2p$ deformation (Fig. 2) |

### `alt_cell` (exploration programme — the signature census and the cell/cycle census; not one of the publication's 115)
| id | verdict |
|---|---|
| `V149_signature_algebra` | only Cl(10,0) and Cl(0,10) keep the 45-plane product exactly unitary; mixed signatures give anti-Hermitian planes and an O(1) unitarity residual |
| `V150_signature_frame` | the 16/16 split survives every signature; the exact (1+3) frame needs four +I generators (q <= 6): with V149 the Euclidean root is the only signature with both |
| `V151_signature_anomaly` | only q=0 reproduces the SM charge table; q=4 permutes the fields, q=10 conjugates the charges, odd q break the multiplets (complex hypercharges) |
| `V152_signature_step` | q=10 gaps the sampled M nodes but has no (1+3) frame; no signature with a frame removes the 12 nodes or isotropizes a cone |
| `V153_stage_cycle_census` | all 720 stage orderings: none removes the half-zone nodes, none reaches A < 1e-6 - the order is not the cause |
| `V154_block_census` | blocks {2,3}^3: front speed exactly 1 everywhere, no isotropy; the (2,2,2) cell keeps its 12 M nodes under both cycles |
| `V155_fcc_diamond` | 524 FCC/diamond 1-factorizations: zero gapless at the FCC high-symmetry points - the graphene mechanism does not transfer |
| `V156_node_forcing` | the 12-node set is one orbit (split forced) for the implemented cell; no candidate combines a forced split with an isotropic cone |
| `V157_hardtier` | zero census survivors; all 722 algebraically eligible near-misses pass the five hard gates (every one a restatement of the root cell) - the hard tier is jointly satisfiable but the isotropy gate eliminates everything |
| `V158_unification_forced` | the node content is a complete 16+16-bar per node: V64's proportional-shift no-go kills the meeting for every candidate |
| `V159_cp_census` | real type at q=0; q=10 is quaternionic AND unitary (a CP window) but has no frame: the CP phase and the Minkowski frame are mutually exclusive |
| `V160_dm_spectrum` | the lightest mode of every candidate is a gapless node with paired mirror content: no 118 GeV carrier from geometry |

### `alt_dynamics` (exploration programme — the step-rule / motion census: per-stage angles, factorizations, Floquet, quantized Wilson layers)
| id | verdict |
|---|---|
| `V161_angle_per_stage` | 729 angle assignments: 720 fully gapped at the high-symmetry points; the 9 gapless rows keep the M nodes (8-12) with A >= 0.949 - no per-stage angle set removes the doubling |
| `V162_internal_factorization` | 12 factorizations of the internal layer: only the two Cartan-preserving ones stay gapless (12/13 nodes, A = 1.0 at M); every non-Cartan factorization gaps the whole zone - the nodes are Cartan-protected |
| `V163_mass_stage` | 90 mass-stage insertions: the M nodes survive every plane/gauge choice (gap >= 5e-4, no true cone); the mirror Z2 is never broken (leakage 0) |
| `V164_bond_matchings` | 500 cubic 1-factorizations: 487 gapless, 31 single-orbit (12 nodes, best A 0.998); zero isotropic, zero closures - the matching is not the cause |
| `V165_floquet` | 3 Floquet time-splittings of the root step: all three open a gap (2e-3 to 1.8e-2) at every high-symmetry point; only the time-independent anchor stays gapless |
| `V166_wilson_quantized` | the quantized Wilson family (X: 1,1,1,1,1; Y: -1/6,1/6,1/6,1/4,1/4): 7/8 variants fully gapped; the one gapless variant keeps 8 nodes with A = 0.993 - the published coefficients are not a free choice |
| `V167_hardtier_dynamics` | zero dynamics survivors: the closure gate eliminates every non-Wilson candidate, the front gate eliminates the mass stages - the Wilson layer is the unique causal/closure-consistent motion |
| `V168_unification_dynamics` | 509 candidate spectra: zero forced meetings (the node content is one mode per node, not a complete 16+16-bar per node at a single momentum) - V64's no-go stands for every motion |
| `V169_cp_dynamics` | 1336 candidates: every effective internal unitary is real type (J^2 = +1); no motion admits a CP phase - the q=0 CP obstruction is motion-independent |
| `V170_dm_dynamics` | 509 candidate spectra: the mirror Z2 survives every motion; the lightest mode of every candidate is a gapless node and no light mode is isolated from its spectrum - no 118 GeV carrier from any motion |

### `alt_interactions` (exploration programme — the collision vertex, self-energy flow, block RG and isotropy gate)
| id | verdict |
|---|---|
| `V171_collision_vertex` | block-contact collision U_c in {g_10, pi/4, pi/2}: the vertex is O_h-scalar and its node-mode form factor is irrelevant on the cone (anisotropic weight ~ k^2.00, A1g = 1.0000) - the interaction adds no anisotropy but a scalar vertex cannot remove the propagator's |
| `V172_self_energy_flow` | one-loop Fock self-energy of the contact: anisotropic part ~ |k|^1.00 = MARGINAL (same power as the cone); dressed anisotropy not grid-convergent (dev 0.019-0.66) - the flat tangent axis (V118) is IR-singular against the collision; one-loop flow inconclusive |
| `V173_block_rg` | energy-window downfold at b = 1, 2, 4: A(b) = 0.9980 exactly, bare and with the self-energy insertion - the cubic anisotropy is MARGINAL under coarse-graining (2-route 3.5e-13) |
| `V174_isotropy_gate` | 18 census rows, 0 candidates (3 unconverged rows excluded): the minimal admissible collision does not isotropize the W9 cone - only the non-perturbative gauged-layer fixed point (V79-V101) remains |

### `alt_gauged` (exploration programme — the two-fermion S-matrix of the gauged bond layer)
| id | verdict |
|---|---|
| `V175_gauged_smatrix` | THE TWO-FERMION S-MATRIX OF THE GAUGED LAYER EXISTS: S = exp(-i g_10 L^-1 H L^-dagger) on the 24-state three-site kernel, unitary to 2.4e-13, spectrum = V101 to 7.8e-16, fusion channels preserved (mixing ~1e-29), eta->0 gives S = I |
| `V176_gauged_vertex_content` | the gauged current vertex is NOT a scalar: Eg + T2g weight 5/6 at |k| = 0.4 (vs contact A1g = 1.0000), anisotropic part scale-independent (alpha = 0 - MARGINAL, vs contact alpha = 2); on-node current F(0) = 5.000000 |
| `V177_gauged_self_energy_flow` | one-loop Fock flow with the gauged vertex: Sigma nonzero (0.199), anisotropic part ~ |k|^1.00 = MARGINAL (propagator sets the power, not the vertex); dressed A not grid-convergent (drift 2.2e-2, flat-axis IR, V172 precedent) |
| `V178_gauged_block_rg` | block RG b = 1, 2, 4 with the gauged insertion: A(b) = 0.997985 EXACTLY unchanged (bare = gauged at every scale) - anisotropy stays marginal (2-route 1.3e-12) |
| `V179_gauged_isotropy_gate` | 12 census rows, 0 candidates (2 unconverged V177 rows excluded, min A 0.9656): the gauged bond-current interaction does not isotropize the W9 cone at one loop / under coarse-graining - the flux sector is massive (gap 64 eta, V129) and decouples from the IR cone |

### `alt_nonpert` (exploration programme — exact cone flux dressing, ring kernels, IR decoupling, non-perturbative isotropy gate)
| id | verdict |
|---|---|
| `V180_cone_flux_dressing` | exact cone dressing by the 24 gauged flux operators: A_np - A_free = +2.1e-12 (2-route 8.5e-13, eta->0 dev 0) - the flux dressing does NOT move the W9 anisotropy at one dressing order |
| `V181_ring_kernel_smatrix` | the 4-site gauged kernels exist: chain 144 (full rank, E0 -1664197.02, gap 12787.60 eta) and ring 1296 (Gram rank 976, E0 -6.08e8, gap 2.75e8 eta); S unitary to 2.7e-15/4.2e-15, 2-route phase dev within the whitening-noise bound; chain channels preserved (mixing 3e-11), ring channels mixed (~8) |
| `V182_ir_decoupling` | INCONCLUSIVE under the pre-registered clauses: gap grows 42.09 -> 12787.60 (chain) -> 2.75e8 (ring) eta, drift >> 1 - the LIGHT branch is decisively false; V101 eta-units are not extensive (Gram shrinks ~30x/link), gap/\|E0\| ~ 0.0073-0.0077 at every volume |
| `V183_nonpert_isotropy_gate` | 9 census rows, 0 candidates: no non-perturbative gauged-sector mechanism isotropizes the W9 cone - the flux sector is massive at every volume probed; W9 stays open |

### `alt_step` (exploration programme)
| id | verdict |
|---|---|
| `V184_step_loss_landscape` | the differentiable cone loss over the 45-plane step space: U(g_10*1) = the implemented step to 1.3e-15, omega = gap_at to 7.2e-12, subgradient pipeline vs finite differences to 2.7e-8, K forward routes to 1.3e-15; 9-row landscape |
| `V188_step_gradient_search` | Adam descent, 25 pre-registered seeds, baseline cell: descent escapes the anisotropy only by closing the gap, or stays gapless at A ~ 1 |
| `V189_cell_step_joint_search` | joint cell x step **census** - the step factorizes as `P(q) x U(theta)^ns`, so no descent: 419 pre-registered steps x 4 cells ((2,2,2) baseline / (2,2,2) phase-major / (2,2,3) / (3,3,3)); 1580/1676 gapless, A 0.9377..1.0000, no isotropic candidate |
| `V190_step_search_gate` | the step-search gate: 0 candidates from 1710 rows scanned (V184 9 + V188 25 + V189 1676); the baseline step remains the only closed-form gapless point reached; W9 stays open |

`V189` runs in 8 min on the CPU and is re-executed by the driver like any
other experiment. The factorization `M(q,theta) = P(q) (x) U(theta)^ns` (bond
layer times internal layer, exact to 8e-16) makes the `q`-half a per-cell
precomputation and a step a 16x16 eigendecomposition, instead of a `(16 n_sub)^3`
diagonalization per optimizer step. Its internal anchor - the implemented step
evaluated through the factorized route and through V184's `(16 n_sub)^3`
eigendecomposition at one node, agreeing in `lam_max` to 1.3e-15 - is what proves
the fast route is the reference route.
