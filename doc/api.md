# CL10 — the API reference

The Python interface of Causal Clifford Theory. One section per module of the
`cl10` package (45 modules), in the order the theory is built: the core, the
algebra root, the automaton, the spectrum and the S-matrix, gravity, the
couplings, the flavour sector, the cosmology, the audits.

Every entry is the real symbol: the signature is the real `def`, the *returns*
column lists the keys the dict actually contains. A symbol that no longer exists
is not listed.

## The provenance vocabulary

The publication (`doc/final_publication.md`, Section 15) classifies every
number. The same classes are used here:

| class | meaning | example |
|---|---|---|
| **exact / algebraic** | an identity of $\mathrm{Cl}(10,0)$ or of the step | the 1024 basis, $5/3$, $\det U = 1$, the 16/16-bar mirror |
| **convention** | a chosen reading, not a derivation | $g_{10} = 1/5$, $\alpha_{10} \equiv g_{10}^2 = 1/25$, $g_0 = 1/8$ |
| **measured** | computed on the automaton | the closed-loop factor $F_{10} = 0.019186$ |
| **fitted** | a threshold or scale read in | $M_N$, the two desert thresholds, the six $\kappa_s$ |
| **conditional** | holds only under a stated choice | $3F_{10}^{10}$ under the one-node-per-family set |
| **open** | one of the ten open problems O1–O10 | the family count, the DM carrier, the couplings meeting |

The old wall labels W1–W9 map one-to-one onto O1–O9 (Section 16); `cl10.walls`
records the mapping.

## 1. `cl10.engine` — the single computation core

The backend is chosen once at import, in the order `backend.priority` of
`cl10/config.yaml`: **CuPy** (GPU, gated on free memory and probed through
cuBLAS + cuSOLVER) → **JAX** (XLA, float64 enabled) → **NumPy** (CPU). Every
module and experiment goes through these primitives.

**The primitives work where the data lives.** A device array is never pulled to
the host by a host call, and a host array is never pushed to the device: that
is what used to run the real-space automaton on the CPU whatever the backend
was. The linear algebra is routed by the size thresholds of `routing:` — above
`eig_gpu_min_dim` an eigendecomposition runs on the device, below it on the
host, because that is measured to be faster here (CuPy ~50 ms vs NumPy ~8 ms on
a 128×128 matrix, and CuPy has no batched eigendecomposition).

| symbol | signature | returns |
|---|---|---|
| `Array` | (type) | the array type of the active backend (`cupy.ndarray` / `jax.Array` / `numpy.ndarray`) |
| `Array_` | `(*args, **kwargs)` | an array on the active backend |
| `eye` | `(n, dtype=None)` | the identity |
| `zeros` | `(shape, dtype=None)` | a zero array (`shape` an int or a tuple) |
| `kron` / `matmul` | `(a, b)` | the product, on the backend where the operands live |
| `einsum` | `(subscripts, *operands)` | `einsum` where the operands live |
| `roll` | `(a, shift, axis)` | the neighbour shift used by the bond layer, where the field lives |
| `real` / `imag` | `(x)` | the real / imaginary part |
| `absmax` | `(x)` | $\max\lvert x\rvert$ as a host float |
| `modulus` | `(a)` | $\lvert a\rvert$ elementwise, where the data lives |
| `summ` / `mean` | `(a)` | the total sum / mean as a host float, reduced where the data lives |
| `sum_axis` | `(a, axis=None)` | a partial sum, where the data lives |
| `where` | `(cond, a, b)` | the elementwise choice, where the data lives |
| `field_norm` | `(a)` | the L2 (Frobenius) norm of a field, reduced on the device when the field is on it |
| `norm` | `(a, *args, **kwargs)` | `linalg.norm`, reduced where the data lives |
| `eigvals` / `eigvalsh` / `eig` / `eigh` | `(a)` | routed by `routing.eig_gpu_min_dim` (host below, device above); a device spectrum from `eigvals`/`eig` is kept only if it satisfies the two trace identities $\sum_a\lambda_a=\operatorname{tr}A$ and $\sum_a\lambda_a^2=\operatorname{tr}A^2$ to `EIG_VERIFY_TOL` ($10^{-6}$ relative), otherwise the whole call is redone on the host — the guard that catches the non-Hermitian device solver folding the spectrum of a large unitary onto $\pm1$ |
| `eig_verify_count` | `()` | `{"calls", "rejected"}` — how many device eigendecompositions were taken and how many failed those identities; carried in the `engine.eig_verification` block of every verdict |
| `det` / `inv` / `qr` / `matrix_rank` | `(a, ...)` | routed by `routing.device_min_dim` (and residency) |
| `to_device` | `(x)` | move a host array to the active backend; an array already on the device comes back untouched (cupy refuses the implicit host conversion `np.asarray` would attempt, so re-moving a device array is a crash unless it is decided here) |
| `asnumpy` | `(x)` | the host copy (printing, JSON) |
| `einops_rearrange` | `(x, pattern_in, pattern_out)` | `einops.rearrange` on the active backend when supported |
| `laplacian_3d` | `(u)` | the 6-neighbour discrete Laplacian on the last three (spatial) axes, where the field lives |
| `spectral_laplacian` | `(k, L)` | the momentum-space Laplacian eigenvalue on an $L^3$ torus |
| `backend_info` | `()` | `{backend, einops, kingdon, array, gpu_min_free_mib, config, ...route_report()}` |
| `route_report` | `()` | where each class of work runs under the current configuration |

### `cl10.config` — the single configuration source

`cl10/config.yaml` holds the numbers the code used to carry as literals;
`cl10/config.py` loads it, applies the environment over it, and exposes the
singleton `CFG`. Precedence: environment → YAML → the built-in defaults.

| symbol | signature | returns |
|---|---|---|
| `CFG` | (object) | attribute access over the effective configuration: `CFG.backend.gpu_min_free_mib`, `CFG.routing.eig_gpu_min_dim`, `CFG.physics.m_planck_gev`, `CFG.physics.loop`, `CFG.physics.n_gen`, `CFG.cell.rung_A` / `CFG.cell.rung_B`, `CFG.sm.m_z_gev`, `CFG.angle.g10`, `CFG.desert.m_x_gev`, `CFG.desert.theta_root_rad`, `CFG.bounds.sigma_si_lz_cm2`, `CFG.cosmology.g_star` / `CFG.cosmology.hubble_coeff`, `CFG.tolerances.strict`, `CFG.pool.njobs_effective` |
| `CFG.get` | `(dotted_key, default=None)` | dotted lookup, `default` when the key is absent |
| `reload` | `()` | re-read the file and the environment into the same `CFG` |
| `describe` | `()` | the effective configuration as a flat dict (for a verdict file) |
| `smoke` | `()` | `True` when `CL10_SMOKE` asks for the shrunken grids |
| `PATH` | (path) | the YAML file being read |

## 2. `cl10.parpool` — the CPU process pool

| symbol | signature | returns |
|---|---|---|
| `njobs_default` | `()` | `int` — `pool.njobs_effective` of `cl10/config.yaml` (`min(njobs_cap, cpus//2)`), overridable by `CL10_NJOBS` |
| `pmap` | `(fn, items, njobs=None, chunk=None, progress=None)` | the results in input order, computed in a forked process pool (the experiments' batch worker) |

## 3. `cl10.paths` — addressing experiments by name

| symbol | signature | returns |
|---|---|---|
| `exp` | `(*parts)` | a path inside `experiments/<category>/<Vxx>/...`, resolved by experiment name (the category is never hard-coded) |

## 4. `cl10.coupling` — the one parameter of the step

The single source of the rotation angle and of the two coupling readings. The
clock is the order-8 root $\pi/4$ (closure axiom, Theorem 6); reference angles
stay addressable. The *values* live in `cl10/config.yaml` (`angle:`,
`default_branch: closed`); the branch is selected by `angle.default_branch` or
by the environment `CL10_ANGLE=implemented|closed|half` (the aliases `pi/4`,
`pi/2`, `1/5`, `impl`, `quantized` are resolved by `cl10/config.py`).

| symbol | signature | returns |
|---|---|---|
| `G_IMPLEMENTED` | (constant) | `0.2` — the reference spectrum $g_{10} = 1/n_{\rm qubits}$ (Wilson-family member): a **convention**, not the clock |
| `G_CLOSED` | (constant) | $\pi/4$ — the closure root: **algebraic** (V131, order 8) |
| `G_HALF` | (constant) | $\pi/2$ — the degenerate second root: **algebraic** (V134) |
| `G0` | (constant) | `1/8` — the $g_0 = 1/2^{n_{\rm space}}$ reading: a **convention** |
| `ALPHA_10` | (constant) | `1/25` — the rotation-squared reading $\alpha_{10} = g^2$: a **convention** |
| `BRANCHES` | (constant) | name → `(angle, provenance, the experiment that fixes it)` |
| `branch` | `()` | the active branch name, from `CL10_ANGLE` (default `closed` — the root clock) |
| `angle` | `(name=None)` | `float` — the angle of the branch (default: the active one) |
| `provenance` | `(name=None)` | `str` — the provenance class of that angle |
| `alpha_std` | `(g=None)` | `float` — the standard-normalised coupling $\alpha = g^2/4\pi$ |
| `alpha10` | `(g=None)` | `float` — the rotation-squared reading $\alpha_{10} = g^2$ |

## 5. `cl10.clifford` — the root $\mathrm{Cl}(10,0)$

| symbol | signature | returns |
|---|---|---|
| `generators` | `()` | the 10 generators as a list of 10 complex `32×32` matrices (exact: $h_a^2 = I$, $\{h_a,h_b\} = 0$) |
| `basis` | `()` | the 1024 ordered products (the algebra basis, exact) |
| `pseudoscalar` | `()` | the volume element $w = h_1\cdots h_{10}$ (anti-Hermitian, $w^2 = -I$): the 16/16-bar split |
| `weyl_spinor` | `()` | the 16-D and 16-bar projectors from the $\pm i$ eigenspaces of $w$ |
| `verify` | `(tol=1e-12)` | `{n, square_to_I, anticommute, ok}` — the algebra relations (exact) |
| `main` | `()` | runs the self-check, exit code = verdict |

## 6. `cl10.spacetime` — the $(1,3)$ inside the GUT

| symbol | signature | returns |
|---|---|---|
| `causal_generators` | `()` | `(T, S)`: time $= h_1$ ($+$), space $= i h_{2..4}$ ($-$) — four of the ten generators, so the Minkowski relations are identities |
| `light_cone` | `(T, S)` | the null vectors $v = v_0 T + v_a S_a$ with $v^2 = 0$ |
| `verify` | `(tol=1e-12)` | `{t2, s2, ts, ss, ok}` — the sub-algebra check (exact) |
| `main` | `()` | runs the self-check |

## 7. `cl10.matter` — the 16 Weyl spinor

| symbol | signature | returns |
|---|---|---|
| `weyl16` | `()` | the 16 and 16-bar eigenspaces of the pseudoscalar |
| `dimension_check` | `()` | the six SM fields sum to Weyl dimension 16 (algebraic) |
| `bl_charge` | `(which='16')` | $B-L$ of the module: 16 → $-1$, 16-bar → $+1$ (algebraic) |
| `anomalies` | `()` | `{su3_cubed, su2_cubed, u1_cubed, grav2_su3, grav2_su2, grav2_u1, u1_su3_sq, u1_su2_sq}` — the eight coefficients, all zero (algebraic) |
| `mass_generation` | `()` | `{n16, n16bar, off_diagonal, diag_16, diag_16bar, A_16_16bar, B_16bar_16, mass_term, sterile_N_in_16, residuals}` — the Dirac structure $16\,\mathbf{10}_H\,\overline{16}$ |
| `main` | `()` | runs the self-check |

## 8. `cl10.cpt` — the CPT theorem

| symbol | signature | returns |
|---|---|---|
| `pseudoscalar` | `()` | $w = g_1\cdots g_{10}$ (anti-Hermitian, $w^2 = -I$) |
| `cpt_operator` | `()` | $C = g_1$ (the time direction): a Hermitian involution |
| `projectors` | `(w)` | the 16 / 16-bar projectors |
| `cpt_theorem` | `(seed=0)` | `{C_involution_err, C_unitary_err, C_w_anticommute_err, map_16_to_16bar_err, gut_commute_w_err, norm_before, norm_after_cpt, norm_rel_err, n_matter, n_antimatter, cpt_equal_norm, c16_into_16bar, backend}` — exact |
| `main` | `()` | runs the self-check |

## 9. `cl10.gauge` — the 45 and the $5/3$

| symbol | signature | returns |
|---|---|---|
| `generators_45` | `()` | the 45 $\mathrm{Spin}(10)$ generators $h_a h_b$ ($a<b$) as `32×32` matrices |
| `rank_check` | `()` | the 45 span the 45-D adjoint (exact) |
| `five_thirds` | `()` | $g_1^2/g'^2 = 5/3$ (algebraic, from the $\mathbf 5$ normalisation) |
| `sm_subgroup` | `()` | `{n_sm, n_gluon, n_w, n_b, n_heavy}` = 12 + 33 |
| `gauge_gravity_overlap` | `()` | `{n_causal, n_gauge, n_generators, one_object}` — the $(1,3)$ and the 45 live in the same $\mathrm{Cl}(10)$ |
| `running` | `()` | `{alpha_10, g_gut, m_gut_0param, five_thirds, b1, b2, b3, unification_spread, cleanly_unifies, residual}` — $\alpha_{10}$ is a **convention**; the 1-loop non-unification is **O4** |
| `main` | `()` | runs the self-check |

## 10. `cl10.sm_qca` — the 45 plane rotations

| symbol | signature | returns |
|---|---|---|
| `plane_generators` | `()` | the 45 $T_{ab} = i h_a h_b$ (Hermitian, $T^2 = I$) |
| `plane_rotation` | `(T, theta)` | $\cos\theta\, I - i \sin\theta\, T$ (exact, since $T^2 = I$) |
| `gut_qca_step` | `(planes, theta)` | the product of the 45 plane rotations: exactly unitary, $\det U = 1$ at every angle |
| `unitary_err` | `(U, tol=1e-12)` | $\lVert U U^\dagger - I\rVert_{\max}$ |
| `norm_preservation` | `(U, n_steps, L=8, seed=0)` | the norm before/after `n_steps` steps on the active backend |
| `sm_qca` | `(n_steps=50, L=8)` | `{n_planes, g10, plane_hermitian_err, plane_square_to_I_err, step_unitary_err, n_steps, norm_before, norm_after, norm_rel_err, backend, causal_radius}` |
| `main` | `()` | runs the self-check |

## 11. `cl10.cell` — the cell module of Axiom 1

The $36$-level internal module $M = S \oplus R_A \oplus R_B$: the $32$-D
Clifford module $S$ and the two neutral singlet sectors $R_A$ (the seesaw
singlet $A$) and $R_B$ (the dark matter $B$), each a conjugate pair of levels
carrying opposite rung phases. The rung integers have one source,
`CFG.cell.rung_A` and `CFG.cell.rung_B` (`verify` pins `rung_B` equal to
`CFG.physics.n_gen`); V196 verifies the cell as a root fact.

| symbol | signature | returns |
|---|---|---|
| `DIM` | (constant) | `36` |
| `RUNG_A` / `RUNG_B` | (constants) | the rung integers $3$ and $10$ from `CFG.cell` |
| `rung_phase` | `(n)` | $\mathrm{loop}^n$ — the phase a rung-$n$ level carries per stage ($= m/M_P$ at the elementary time $t_P$) |
| `phase_A` / `phase_B` | `()` | $\mathrm{loop}^{\mathrm{rung}_A}$, $\mathrm{loop}^{\mathrm{rung}_B}$ |
| `cell_step_np` | `(g10=None)` | the $36\times36$ internal layer of one stage: $\operatorname{blkdiag}(U(g_{10}), e^{+i\phi_A}, e^{-i\phi_A}, e^{+i\phi_B}, e^{-i\phi_B})$, host numpy |
| `parity_np` | `()` | the B-parity $G = \operatorname{diag}(I_{32}, +1, +1, -1, -1)$ |
| `gauge_ext_np` | `()` | the $45$ internal-plane generators $h_a h_b$ extended by $0$ on the singlet sectors |
| `masses` | `()` | `{m_A_GeV, m_B_GeV, m_nu_eV_identity}` — the two rung masses and the exact seesaw identity $m_\nu = m_B^2/m_A$ |
| `content` | `()` | hypercharge content as $(Y, \text{multiplicity})$ classes — the six fragment classes and their conjugates plus $(0, 4)$; $\operatorname{Tr}Y = \operatorname{Tr}Y^3 = 0$ |
| `verify` | `()` | `{ok, unitarity, pairing, det_phase, vacuum_sum, phase_A, phase_B, parity_step_comm, parity_gauge_max, G_squared, trY, trY3, m_A_eq_M_N, m_B_eq_dm_ladder, m_nu_identity_rel, m_A_GeV, m_B_GeV, rung_A, rung_B, rung_B_eq_n_gen}` — $\pm$ pairing to $2\times10^{-15}$, $\det$ phase to $3\times10^{-16}$, quasi-energy sum to $2\times10^{-16}$ on the full cell |
| `main` | `()` | runs the self-check (VERDICT line) |

## 12. `cl10.supercell` — the automaton

The cell is the $2×2×2$ block: 8 sub-sites × the $36$-level cell module of
Axiom 1 (`cl10.cell`) = on-site $\mathbb{C}^{288}$. One time-independent step:
six stages, each an internal layer plus a bond stage on one parity of one
axis. The Bloch and real-space machinery below acts on the $32$-D spinor block
(the physical sector the spectrum and light cone are measured on):
$\mathbb{C}^{256}$, $128\times128$ on the $16$; the full-cell transport is
$\operatorname{kron}(P(q), U_{\rm cell}^6)$, $288\times288$ (V196).

| symbol | signature | returns |
|---|---|---|
| `internal` | `(g10=0.2)` | the $32×32$ internal unitary on the active backend |
| `sectors` | `(g10=0.2)` | the 16 / 32 internal unitaries and the 16 basis |
| `stage_matrix` | `(even, j, q, g10=0.2)` | the $8×8$ sub-site matrix of one bond stage at block momentum `q` |
| `bloch` | `(q, U, g10=0.2, stages=DEFAULT)` | the exact one-step Bloch matrix $M(q)$ on the block zone |
| `phases` | `(M)` | the eigenphases of a Bloch matrix, in $(-\pi,\pi]$ |
| `gap` | `(q, U, g10=0.2, stages=DEFAULT)` | $\min_a \lvert\omega_a(q)\rvert$ |
| `degeneracy` | `(q, U, g10=0.2, stages=DEFAULT, tol=1e-9)` | the eigenphase multiplets `[(gap, multiplicity)]` |
| `star` | `(q)` | the $O_h$ orbit of `q` |
| `stabilizer_order` | `(q)` | the order of the little group ($48/\lvert\text{star}\rvert$) |
| `internal_layer` | `(fields, U)` | `U` on every sub-site of every block |
| `bond_layer` | `(fields, even, j, theta=0.2)` | one bond stage: the swap-rotation on the bonds of parity `even` along axis `j` |
| `stage` | `(fields, U, t, theta=0.2)` | one elementary stage (internal layer + bond stage) |
| `step` | `(fields, U, theta=0.2, stages=DEFAULT)` | one automaton step (six stages) |
| `run` | `(L=8, n_steps=4, x0=(4,4,4), e0=0, theta=0.2, seed=0, tol=1e-12, sector=True)` | `{backend, L, n_steps, theta, g10, sector, steps, final_norm}` — the real-space light cone from one sub-site |
| `verify` | `(tol=1e-9)` | the seven invariants on the block: unitarity, the causal radius, the 16-sector closure, the spectral mirror, the gap structure, the zone-centre closure and cubic equivariance, the stars and little groups, and the cell module of Axiom 1 (invariant 7: the 36-level step unitary, the B-parity commuting with the step and with all $45$ generators) |
| `main` | `()` | prints the invariants and the verdict |

## 13. `cl10.qca` — the metric field and the on-site factor

The metric is a **field added to the algebra**, not a number of it: the sector
assumes a symmetric rank-2 field and $G = 1$ (Section 11.1). The 3+1 evolution
lives in `cl10.einstein` / `cl10.gravity` / `cl10.numrel`.

| symbol | signature | returns |
|---|---|---|
| `METRIC` | (constant) | the $(1,3)$ vacuum metric $\mathrm{diag}(1,-1,-1,-1)$ |
| `metric_field` | `(L, h=0.0)` | the `(4,4,L,L,L)` metric field on an $L^3$ torus |
| `qca_step_internal` | `(gates, fields)` | the on-site internal factor on the last axis of `(L,L,L,n)` fields |
| `unitary_check` | `(U, tol=1e-12)` | `bool` — whether $UU^\dagger = I$ (backend-correct) |
| `main` | `()` | runs the self-check |

## 14. `cl10.qca_spectrum` — the one object: the 32-level spectrum

| symbol | signature | returns |
|---|---|---|
| `levels` | `(g10=None)` | the 32 eigenphases of the step, in $(-\pi,\pi]$ |
| `cpt_block_phases` | `(g10=None)` | `{ph16, ph16b}` — the two 16-element phase blocks |
| `qca_spectrum` | `(g10=None)` | `{g10, n_levels, spectrum, zero_point, det_U_minus_1, det_U_arg, sl32, pairing_resid, max_trace_T45, cpt_0point_16, cpt_0point_16bar, commutes_with_pseudoscalar, max_level}` — the 0-point, the $\pm$ pairing and $\det U = 1$ are exact for **every** angle |
| `verify` | `(tol=1e-10)` | `{spectrum, ok}` |
| `main` | `()` | runs the self-check |

## 15. `cl10.smatrix` — the S-matrix of the step

| symbol | signature | returns |
|---|---|---|
| `qca_step` | `()` | the $32×32$ GUT step at the active angle |
| `one_particle` | `(U)` | `{n, unitary_err, colsum_err, colsum_err_col, raw_colsum_dev, max_offdiag, nontrivial_offdiag}` |
| `two_particle` | `(U)` | `{shape, unitary_err, colsum_err, max_abs, max_offdiag, nontrivial}` — $S^{(2)} = U \otimes U$ |
| `n_particle` | `()` | `{backend, n1_unitary_residual, n2_unitary_residual, n3_unitary_residual, n4_unitary_residual, n1_bounded, uv_finite, residual}` — bounded at 3 and 4 particles |
| `smatrix` | `()` | `{backend, g10, n_planes, n_internal, one_particle, two_particle, only_g10, residual}` |
| `main` | `()` | runs the self-check |

## 16. `cl10.lattice_rg` — cross-sections and the UV bound

| symbol | signature | returns |
|---|---|---|
| `qca_step` | `()` | the $32×32$ step at the active angle |
| `one_particle_inelastic` | `(U)` | `{n, dbar, P1}` — the 1-particle inelastic probability per step |
| `two_particle` | `(U)` | `{shape, unitary_err, P2_direct, P2_closed, P2, closed_form_match}` |
| `uv_bound` | `()` | `{M_P_GeV, inv_M_P2_GeV2, sigma_max_GeV2, sigma_max_cm2}` — the unitarity bound set by the cutoff |
| `dispersion` | `()` | `{n_levels, max_omega, omega_max_over_pi, uv_bounded, zero_point, det_U_minus_1}` |
| `cross_section` | `()` | `{g10, M_P_GeV, P1, P2, elastic, optical_theorem, sigma_0_GeV2, sigma_0_cm2, sigma_pert_GeV2, sigma_pert_cm2, R_qca_over_pert, sigma_max_GeV2, uv_finite}` |
| `froissart` | `()` | `{m_min_GeV, s0_GeV2, sigma_qca_GeV2, sigma_max_GeV2, constant_uv_finite, below_froissart, strongest_uv_finite, rows}` |
| `lattice_rg` | `()` | the full verdict: `{backend, g10, M_P_GeV, one_particle, two_particle, uv_bound, dispersion, cross_section, froissart, uv_finite, optical_theorem, s2_unitary_err, closed_form_match, only_g10, residual}` |
| `main` | `()` | runs the self-check |

## 17. `cl10.einstein` — the linearised graviton and the EH action

| symbol | signature | returns |
|---|---|---|
| `light_cone` | `()` | the $(1,3)$ sub-algebra check |
| `linearized_graviton` | `(L=20)` | `{backend, L, c, massless, helicity, n_polarizations, spin2_declared, light_cone_ok, continuum_ratio, omega_max, uv_bounded, stable}` — the massless lattice dispersion at $c = 1$ |
| `riemann_scalar_2d` | `(metric, x0, x1, h=1e-3)` | the 2-D Ricci scalar of a prescribed metric |
| `einstein_hilbert` | `(r=2.0, tol=1e-4)` | `{backend, G, coupling, R_sphere, R_sphere_expected, R_flat, quadratic, cubic, quartic, n_deriv_cubic, n_deriv_quartic, einstein_hilbert_ok}` — $G = 1$ is an **input** of the sector |
| `eh_quadratic_action` | `(u)` | the discrete quadratic EH action on a lattice field |
| `einstein_hilbert_from_dynamics` | `(L=8, T=24, n_probe=16, seed=0, tol=1e-6)` | `{backend, L, T, n_probe, G, coupler, elm_on_qca_trajectory, stationary_ok, elm_off_dynamics, selects_ok, coeff_ok, k_values, continuum_ratio, continuum_ok, ok}` — the action selects the QCA trajectory |
| `nonlinear_soliton` | `()` | `{soliton, ode_residual, mass, nonlinear_ok}` — the static kink of $u'' = u - u^3$ (a planar domain wall, not a relic) |
| `schwarzschild` | `(M=1.0, r=3.0)` | `{M, horizon, Kretschmann, R_vacuum}` |
| `verify` | `(tol=1e-9)` | `{linearized, einstein_hilbert, eh_from_dynamics, light_cone, soliton, black_hole, ok}` |
| `main` | `()` | runs the self-check |

## 18. `cl10.gravity` — curvature computed from any metric

| symbol | signature | returns |
|---|---|---|
| `christoffel_4d` | `(g, dg)` | the 4-D Christoffel symbols |
| `riemann_4d` | `(mfn, x, h=1e-3)` | the full $R^\rho{}_{\sigma\mu\nu}$ of an arbitrary metric (finite differences) |
| `ricci_4d` / `ricci_scalar_4d` / `einstein_tensor_4d` / `kretschmann_4d` | `(R, …)` | the contractions of the computed Riemann tensor |
| `gamma_gamma_4d` | `(mfn, x, h=1e-3)` | the non-linear $\Gamma\Gamma$ part of the Ricci tensor |
| `metric_flat` / `metric_de_sitter` / `metric_schwarzschild` / `metric_2sphere` | `(…)` | the test metrics |
| `discrete_eh_action` | `(mfn, grid, Lambda, h=1e-3)` | the complete discrete Einstein-Hilbert action |
| `dynamics_from_action` | `(mfn, Lambda, x0, h=1e-3, eps=1e-4)` | the functional derivative of the action (dynamics ← action) |
| `nonlinear_qca` | `(L=16, n_steps=12, x0=None, seed=0, lam=0.2, tol=1e-12)` | `{backend, L, n_steps, g10, lam, norm_before, norm_after, norm_rel_err, unitary_ok, superpose_err, nonlinear_ok, front_speed, causal_ok, uv_finite, residual, ok}` — the state-dependent metric evolution |
| `verify_analytic_vs_computed` | `(h=1e-3)` | `{H, max_abs_diff, h2, ok}` — the computed curvature against analytic de Sitter |
| `verify_blackhole` | `(h=1e-3)` | `{M, r, ricci_max, ricci_scalar, einstein_max, kretschmann, kretschmann_expected, kretschmann_ratio, vacuum_ok, kretschmann_ok}` |
| `verify_kretschmann_scan` | `(h=1e-3, radii=None, m=1.0)` | `{h, radii, computed, expected, max_rel_dev, max_rel_dev_at, ok}` — the 2-route Kretschmann scan along a radial line (the right panel of Fig. 12; the maximum is carried in V50's verdict) |
| `verify_flat` / `verify_nonlinear` | `(h=1e-3)` | the flat and the $\Gamma\Gamma$ checks |
| `verify` | `(tol=1e-9)` | `{backend, analytic_vs_computed, black_hole, flat, nonlinear_gamma_gamma, nonlinear_qca, residual, ok}` |
| `main` | `()` | runs the self-check |

## 19. `cl10.nonlinear_gravity` — the analytic de Sitter reference

| symbol | signature | returns |
|---|---|---|
| `de_sitter_metric` | `(H, x0=0.0)` | the flat-FRW de Sitter metric |
| `riemann_4d` | `(H, x0=0.0)` | its analytic Riemann tensor |
| `nonlinear_gravity` | `(H=None)` | `{backend, H, einstein_space_err, einstein_space_ok, g_mu_nu_minus_3h2_g, g_ok, ricci_scalar_minus_12h2, scalar_ok, gamma_gamma_max, nonlinear_present, lambda_0param_GeV4, lambda_obs_GeV4, gap, residual, ok}` |
| `main` | `()` | runs the self-check |

## 20. `cl10.numrel` — the 3+1 evolution

| symbol | signature | returns |
|---|---|---|
| `gamma` / `K` / `gamma_inv` / `K_contra` / `K_trace` / `K_K` / `K_frob` | `(g4fn, x, …)` or `(gamma_, K_)` | the 3+1 geometric variables and their contractions |
| `G3` / `D3_divergence` / `R3` / `R3_ricci` / `R3_scalar` | `(g4fn, x, h=1e-3)` | the 3-D Christoffel, covariant divergence, Riemann, Ricci and scalar curvature |
| `ricci_0row_3p1` / `ricci_0row_4d` | `(g4fn, x, h=1e-3)` | the 4-D Ricci 0-row by the two routes |
| `verify_2route` | `(g4fn, n_pts=24, h=1e-3, tol=1e-4, seed=11)` | `{max_r00, max_r0i, h2, ok}` |
| `constraints` | `(g4fn, x, h=1e-3)` | the Hamiltonian and momentum constraints |
| `nonlinear_superpose` | `(a_fn, b_fn, x, h=1e-3, tol=1e-8)` | `{superpose_err, ok}` — non-additivity of the non-linear theory |
| `z4_damping` | `(g4fn, x, n_steps=20, h=1e-3, lam=1.0)` | `{H0, damped_final, undamped_final, damped_ratio, undamped_ratio, ok}` |
| `uv_finite` | `(g4fn, x, h=1e-3)` | `{max_abs, bounded, ok}` |
| `metric_31` | `(x)` | the 3+1 test metric |
| `verify` | `(tol=1e-4, seed=11)` | `{backend, two_route_0row, r00_3p1, r00_4d, r0i_3p1, r0i_4d, nonlinear, z4, uv_finite, damp, ok, residual}` |
| `main` | `()` | runs the self-check |

## 21. `cl10.holo` — the area law

| symbol | signature | returns |
|---|---|---|
| `area_law` | `(L_values=None, per_plaquette=0.25)` | `{L_values, per_plaquette, S, R, area_law_exponent, holographic_ratio_exponent, backend}` — $S = A/4$ at $G = 1$ (the value of $G$ is an input) |
| `main` | `()` | runs the self-check |

## 22. `cl10.qcd` — the area law and the linear potential

| symbol | signature | returns |
|---|---|---|
| `color_count` | `()` | the 8 gluons as $N_C^2 - 1$ with $N_C = 3$ |
| `string_tension` | `()` | the flux-tube energy per unit length: $\sigma = 1$ in Planck units, an **input** |
| `confining_potential` | `(sigma, r_values=None)` | the linear potential against the Coulombic one |
| `area_law_coeff` | `(rows, sigma)` | the slope of (loop action / area) |
| `main` | `()` | runs the self-check |

## 23. `cl10.unification` — the running and the $5/3$

| symbol | signature | returns |
|---|---|---|
| `M_Z` / `M_P` / `B_SM` / `A_Z` / `B2_SM` | (constants) | 91.1876 GeV; the cutoff $M_P$; the 1-loop and 2-loop SM $\beta$-function content and the measured $\alpha_i(M_Z)$ |
| `sm_running_to` | `(M)` | `(a1, a2, a3, spread)` — the 1-loop SM running to `M` |
| `landau_pole_u1` | `()` | the SM $\mathrm U(1)$ Landau-pole scale |
| `g10_from_root` | `()` | `{backend, n_module, n_qubits, g10, alpha_10, g10_matches_sm_qca, alpha_10_matches_u10pred, n_space, g0, g0_matches_coupling}` — the root numbers; the readings are **conventions** |
| `running_down` | `()` | `{backend, alpha_10_at_M_P, couplings_at_M_Z, measured_at_M_Z, ratio_to_measured, a1_within_factor_2, u1_up_pole_scale, u1_pole_above_cutoff, su3_down_pole_scale, su3_pole_inside_lattice, residual}` — the mismatch is **O4** |
| `unification` | `()` | `{alpha_10, five_thirds, M_GUT_0param, couplings_at_M_P, spread_at_M_P, mssm_unif_scale, mssm_spread, u1_landau_pole, uv_bounded, backend}` — $5/3$ algebraic; $M_{\rm GUT} = M_P$ is **not** assumed |
| `running_2loop` | `()` | `{backend, qca_2step_unitary_residual, qca_2step_finite, qca_2particle_bounded, qca_2particle_max_offdiag, pert_1loop_landau_pole, pert_2loop_landau_pole, pert_2loop_divergent, uv_finite, residual, ok}` |
| `main` | `()` | runs the self-check |

## 24. `cl10.uvrunning` — the discrete RG

| symbol | signature | returns |
|---|---|---|
| `one_loop_running` | `(M)` | the 1-loop SM couplings at `M` |
| `landau_poles` | `(M)` | the 1-loop pole scale of each coupling |
| `two_loop_pole_u1` | `()` | the 2-loop $\mathrm U(1)$ pole (below the cutoff: the perturbative divergence) |
| `discrete_rg` | `(max_power=8)` | `{backend, n, rows, all_bounded, two_step_unitary_err, two_particle_max_offdiag, two_particle_bounded}` — the coarse-graining tower $U^{2^p}$ is bounded (exact) |
| `uv_verdict` | `()` | `{backend, alpha_10, M_P_GeV, couplings_at_M_P, spread_at_M_P, one_loop_landau_poles, two_loop_u1_pole, u1_1loop_above_cutoff, u1_2loop_below_cutoff, discrete_rg, discrete_rg_bounded, uv_bounded, residual, ok}` — boundedness is exact, the coupling values are **O4** |
| `verify` | `(tol=1e-9)` | `{ok, checks, verdict}` |
| `main` | `()` | runs the self-check |

## 25. `cl10.g10_selection` — the data check of the angle reading

| symbol | signature | returns |
|---|---|---|
| `gut_1loop` | `()` | `{a_at_M_P, mean_GUT, spread}` — the 1-loop GUT-range mean and the non-meeting spread (**O4**) |
| `g10_selection` | `()` | `{backend, n_qubits, g10, alpha_10_0param, g10_matches_sm_qca, gut_1loop, ratio_to_1loop, best_n_1loop, n_closest_to_1loop, candidate_scan_gap, gap_root_n5, gap_best_n, fit_prefers_root, data_check_selects_n, one_line}` — the scan prefers $n = 7$, so $1/n_{\rm qubits}$ is a **convention** |
| `verify` | `(tol=1e-12)` | `{r, ok}` |
| `main` | `()` | runs the self-check |

## 26. `cl10.electroweak` — the electroweak sector

| symbol | signature | returns |
|---|---|---|
| `qca_step_spectrum` | `()` | `{backend, theta, n_levels, spectrum_spread, unitary_err, effective_action}` — the step as the QFT time evolution ('t Hooft's program) |
| `electroweak_breaking` | `()` | `{backend, g10, g2, gprime, vev_v, m_W, m_Z, sin2_w, sin2_w_0p, mw_over_mz, mw_over_mz_0p, electroweak_breaking_ok}` — $\sin^2\theta_W = 3/8$ is the **GUT boundary condition** from the $5/3$, not a prediction; the running down to $M_Z$ is **O4** |
| `verify` | `(tol=1e-9)` | `{spectrum, electroweak, ok}` |
| `main` | `()` | runs the self-check |

## 27. `cl10.family` — the family index

| symbol | signature | returns |
|---|---|---|
| `family_construct` | `()` | `{n_space, n_planes_C32, n_rank, n_family_candidate, n_families_observed, match, all_three_ways_3, n_weyl_1gen, n_weyl_3gen, residual}` — three routes to the number 3; the count itself is **O3** |
| `tribimaximal_3x3` | `()` | the Harrison-Perkins-Scott unitary |
| `family_qca` | `()` | `{backend, g10, n_planes, n_module_1gen, n_module_3gen, n_weyl_1gen, n_weyl_3gen, n_cpt_3gen, unitary_err, family_symmetry_off_block, family_symmetry_block_diff, max_45_pseudoscalar_commutator, step_32_pseudoscalar_commutator, step_96_pseudoscalar_commutator, tribimaximal_unitary_err, tribimaximal_angles_ok, tribimaximal_diagonalizes_democratic, tribimaximal_sin2, residual}` |
| `main` | `()` | runs the self-check |

## 28. `cl10.yukawa` — the democratic base

| symbol | signature | returns |
|---|---|---|
| `yukawa_matrix` | `()` | `{n_gen, n_dim, rank, spectrum, n_distinct_levels, n_heavy, n_light_2fold, democratic, singular, off_diagonal_Dirac, diag_zero, mass_term, residual}` |
| `democratic_matrix` | `()` | the all-ones rank-1 matrix |
| `tribimaximal_mixing` | `()` | `{sin2_12, sin2_23, sin2_13, obs, backend, diagonalizes_democratic, sin2_12_is_0param_derived, sin2_23_is_0param_derived, sin2_13_is_0param_derived, sin2_23_O2_range, sin2_13_O2_range, mixing_status}` — only $\sin^2\theta_{12} = 1/3$ is fixed |
| `hierarchy_residual` | `()` | `{base_spectrum, degenerate, residual}` — the 2-level base, the 3-level hierarchy is **O6** |
| `main` | `()` | runs the self-check |

## 29. `cl10.yukawa_3x3` — the three-generation matrices and mixings

| symbol | signature | returns |
|---|---|---|
| `one_object_spectrum` | `(g10=None)` | `{g10, n_levels_1gen, spectrum_1gen, spectrum_3gen, n_levels_3gen, zero_point_1gen, det_U_minus_1, family_triplication}` |
| `higgs_mass_16` | `()` | `{shape, n16, rank, off_diagonal, diag_16, non_singular, dirac}` — the $16\,\mathbf{10}_H\,\overline{16}$ Dirac block |
| `family_symmetry` | `()` | `{n_space, n_families_C32, S_a2_err, anticommutator_Sa_Sb, spatial_norms, norm_spread, equivalent_3fold, democratic_rank, democratic_spectrum, two_level}` |
| `mass_matrix_3x3` | `()` | `{n, matrix, eigenvalues, rank, n_heavy, n_light_2fold, n_distinct_levels, democratic}` |
| `mixing_from_democratic` | `()` | the honest analysis: `{spectrum_fixed, heavy_eigenvector, heavy_eigenvector_unique, light_degeneracy, diagonalizer_is_unique, sin2_12_fixed, sin2_12_value, sin2_12_range, sin2_23_range, sin2_13_range, sin2_12_vs_observed, max_diagonalization_err, tbm_diagonalizes_democratic, tbm_in_the_O2_family, tbm_angles, s3_transpositions_equivalent, ckm_is_identity, verdict_0param}` |
| `mixing_3x3` | `()` | `{unitary, diagonalizes_democratic, sin2_12, sin2_23, sin2_13, obs_pmns, obs_ckm, dev_pmns, dev_ckm, sin2_12_is_0param_derived, sin2_23_is_0param_derived, sin2_13_is_0param_derived, sin2_23_O2_range, sin2_13_O2_range, tbm_is_the_23_adapted_choice, mixing_status}` |
| `hierarchy_3level` | `()` | `{n_levels_0param, n_levels_observed, n_heavy, n_light_2fold, structure, degenerate_light, obs_ratio_order, residual}` |
| `quark_ckm` | `()` | `{backend, up_eigenvalues, down_eigenvalues, two_level, ckm_unitary, ckm_is_identity, ckm_computed_angles, sin2_12, sin2_23, sin2_13, obs_ckm, dev_ckm, residual, ok}` — the computed CKM **is** the identity (**O6**) |
| `w6_verdict` | `()` | the assembled verdict of the above (the old W6 label = **O6**) |
| `main` | `()` | runs the self-check |

## 30. `cl10.neutrino` — the seesaw on the ladder

| symbol | signature | returns |
|---|---|---|
| `loop_ladder` | `(n)` | the ladder scale $S(n) = M_P (1/16\pi)^n$ |
| `dirac_mass` | `()` | $m_D = S(10)$ |
| `twostep_scale` | `()` | $M_N = S(3)$ — the ladder choice, **fitted** against the desert scale |
| `twostep_unitary` | `(tol=1e-9)` | `{backend, n, unitary_err_1step, unitary_err_2step, two_step_unitary, uv_finite}` |
| `seesaw_2step` | `()` | `{backend, M_P_GeV, n_gen, n_dirac, n_2step, exp_2step, m_D_GeV, M_N_GeV, m_nu_2step_eV, m_nu_1step_eV, m_nu_ladder_eV, m_nu_obs_eV, off_2step, off_1step, right_order, ladder_match, residual}` — the ladder value 0.146 eV is **excluded** by $\sum m_\nu < 0.12$ eV (F3); **O2** |
| `verify` | `(tol=1e-12)` | `{ok, checks, seesaw, twostep}` |
| `main` | `()` | runs the self-check |

## 31. `cl10.u10pred` — the headline numbers

| symbol | signature | returns |
|---|---|---|
| `M_PLANCK` / `V_EW` / `ALPHA_10` / `M_P_ROBSON` | (constants) | 1.220890e19 GeV (the cutoff, not $M_{\rm GUT}$); 246.22 GeV (measured); the $\alpha_{10}$ reading; 0.93827208943 GeV (the proton mass) |
| `LOOP` / `RUNG_A` / `RUNG_B` | (constants) | the ladder factor $1/(16\pi)$; the rung integers $3$ and $10$ from `CFG.cell` (the same source `cl10.cell` reads) |
| `proton_lifetime` | `(M_X, alpha_X)` | $\tau_p \sim M_X^4/(\alpha_X^2 m_p^5)$ in years (a dimensional estimate) |
| `proton_lifetime_desert` | `()` | the lifetime in force: $M_X = 1.8\times10^{16}$ GeV, $\alpha_X = 1/40$ — **conditional** on two fitted thresholds (V64) |
| `proton_lifetime_closed` | `()` | the same content at the closed angle $\alpha_X = \pi/64$ (V133): the threshold ordering is destroyed |
| `ladder_scale` | `(n)` | $M_P (1/16\pi)^n$ in GeV |
| `dm_ladder_mass` | `()` | the mass **scale** $M_P(1/16\pi)^{\mathrm{RUNG\_B}} = 118$ GeV (rung from `CFG.cell`); the carrier is the cell singlet $B$ (**O7**) |
| `neutrino_ladder` | `()` | `{m_D_GeV, M_N_GeV, m_nu_eV_ladder, m_nu_eV_third, sum_m_nu_candidate_eV, status_ladder, status_third, provenance_exponent}` |
| `hierarchy` | `(F10=0.019186)` | `{F10, v_over_MPlanck, coefficient, observed_v_over_MPlanck, observed_coefficient}` — $v/M_P = 3F_{10}^{10}$, **conditional** on the one-node-per-family set (**O8**) |
| `family_number` | `()` | `{value, identification, provenance}` — provenance **open** (**O3**) |
| `verify` | `(tol=0.02)` | `{values, quoted, relative_deviation, ok}` — code-and-paper agreement, nothing more |
| `main` | `()` | prints the numbers with their provenance |

## 32. `cl10.walls` — the three named walls

The module keeps the historical W labels and records the mapping W1→O1,
W2→O2, W3→O3.

| symbol | signature | returns |
|---|---|---|
| `w1_cosmological_constant` | `()` | `{n_even, n_odd, cancellation, lambda_0param, lambda_obs, residual}` — the 512/512 cancellation is exact; the observed value is **O1** |
| `w2_seesaw` | `()` | `{g0, m_D_GeV, M_N_GeV, m_nu_2step_eV, m_nu_1step_eV, m_nu_obs_eV, off_1step, off_2step, right_order, residual}` — $M_N$ is the ladder 2-step scale `neutrino.twostep_scale` (the loop integer fitted); the right order only; **O2** |
| `w3_families` | `()` | `{n_time, n_space, n_fam_candidate, n_observed, match, derived, residual}` — `derived` is **False**: the index labels the copies, it does not force the count (**O3**) |
| `baryon_asymmetry` | `()` | `{baryon_number_violated, c_or_cp_violated, out_of_equilibrium, all_sakharov, generated, residual}` — the unitary step meets no Sakharov condition (**O5**) |
| `main` | `()` | runs the self-check |

## 33. `cl10.cosmology` — the vacuum and the hierarchy

| symbol | signature | returns |
|---|---|---|
| `cosmological_constant` | `()` | `{backend, n_even, n_odd, cancellation, lambda_0param, lambda_0param_natural, lambda_obs, w, residual}` — the cancellation is exact, $w = -1$; the observed $10^{-122}$ is **O1** |
| `hierarchy` | `()` | `{backend, E_P, v, v_over_EP_0param, observed, rel_err, vast_gap, residual}` — the two scales, the gap unexplained (**O8**) |
| `verify` | `(tol=1e-9)` | `{cc, hierarchy, ok}` |
| `main` | `()` | runs the self-check |

## 34. `cl10.cosmology3` — the dark sector

| symbol | signature | returns |
|---|---|---|
| `dark_matter` | `()` | `{backend, ode_residual, mass_energy_E, non_dispersive, gauge_charge, stable_massive_non_dispersive_neutral}` — the exact kink of $u'' = u - u^3$: a planar domain wall, **not** a carrier |
| `dm_mass_scale` | `()` | `{backend, ode_residual, mass_energy_E, non_dispersive, gauge_charge, mass_gev, sigma_gev2, omega_obs, omega_0param, abundance_gap, under_abundant, residual}` — the V42 quantification of **O7** |
| `dark_energy` | `()` | `{backend, n_even, n_odd, cancellation, lambda_0param, w, accelerated, lambda_obs}` |
| `antimatter` | `()` | `{backend, n16, n16bar, cpt_anticommute, cpt_to_bar, ortho_16_16bar, bl16, bl16bar, bl_cancel, cpt_ok}` — exact |
| `verify` | `(tol=1e-9)` | `{dark_matter, dark_energy, antimatter, ok}` |
| `main` | `()` | runs the self-check |

## 35. `cl10.dm_wimp` — the ladder WIMP scale

| symbol | signature | returns |
|---|---|---|
| `F10` | `()` | $(1/16\pi)^{10}$ |
| `relic_abundance` | `(m_GeV)` | the V42-convention contrast $\Omega_{\rm obs}(100\,\mathrm{GeV}/m)^3$ — normalised to the observed value at the reference scale **by construction**, so it is a contrast device, not an abundance prediction |
| `dm_wimp` | `()` | `{backend, F10, n_gen, M_P_GeV, m_DM_GeV, weak_scale_GeV, m_over_weak, omega_DM, omega_obs, omega_ratio_to_observed, v42_planck_omega, v42_planck_gap, v42_quantified_gap, carrier_status, relic_convention, one_line}` — a scale, not a carrier (**O7**) |
| `verify` | `(tol=1e-12)` | `{r, ok}` — `ok` is the ladder arithmetic, the WIMP-window placement, the relic-convention identity and the cross-module agreement with `cosmology3`; it does not test the carrier |
| `main` | `()` | runs the self-check |

## 36. `cl10.hier_exponent` — the exponent of the hierarchy

| symbol | signature | returns |
|---|---|---|
| `loop_factor` | `(n=10)` | $(1/16\pi)^n$ |
| `loop_ladder` | `()` | `{F10, n_gen, n_qubits, g10_as_2x5, hierarchy, dm}` — one ladder length, two scales |
| `nstep_spectrum` | `(max_n=100, target=2.02e-17)` | `{n_levels_1step, min_level_1step, min_level_small_n, max_n_scanned, any_level_at_target, n_diophantine_est, diophantine_is_0param}` — the one-step spectrum does **not** carry the exponent |
| `hierarchy17` | `()` | `{backend, loop_factor, n_gen, predicted_ratio, predicted_exponent, observed_ratio, observed_exponent, exponent_abs_diff, ratio_predicted_over_observed, scan, n_best_scan, nstep, loop_ladder, one_line}` — the crossing at $n = 10$ is selected by the data (**O8**) |
| `verify` | `(tol=1e-12)` | `{r, ok}` |
| `main` | `()` | runs the self-check |

## 37. `cl10.baryogenesis` — the Sakharov audit

| symbol | signature | returns |
|---|---|---|
| `sakharov_conditions` | `()` | `{backend, b_violation_met, max_45_pseudoscalar_commutator, cpt_invariant, in_equilibrium, unitary_err, any_sakharov_met}` — none met (a computed theorem) |
| `gut_ingredients` | `()` | `{backend, xy_process, proton_lifetime_yr, proton_lifetime_provenance, seesaw_scale_GeV, hubble_at_seesaw_GeV, equilibrium_by_yukawa, out_of_equilibrium_yukawa_ceiling, cp_phase_in_the_step, cp_phase_carrier}` |
| `K_STRONG_WASHOUT_MIN` | `3.9345` | the branch point of the BDP washout solution: below it the strong-washout equation has no real root |
| `bdp_washout_efficiency` | `(K)` | `(kappa_f, z_B, in_domain)` — the Buchmüller-Di Bari-Plumacher washout efficiency (hep-ph/0401240) at $K = \Gamma_N/H$; out of domain the efficiency is the weak-washout bound $1$ and $z_B$ is `None` |
| `leptogenesis_estimate` | `()` | `{backend, M_N_GeV, y, gamma_over_h, out_of_equilibrium, out_of_equilibrium_yukawa_ceiling, g0, eps_loop, loop_factor_0_parameter, eta_cp, sphaleron_factor, washout_efficiency, washout_z_B, washout_strong_domain, delta_b_unwashed, unwashed_ratio_to_observed, delta_b, observed_delta_b, ratio_to_observed, order_ok, agreement_survives_washout, 0param_qca_satisfies_no_sakharov, residual}` — the right order with named inputs, **with the physical washout** (**O5**): the decay at the model's own Yukawa is in equilibrium, and the unwashed product is reported only to expose the coincidence |
| `verify` | `(tol=1e-9)` | `{sakharov, ingredients, ok}` — the negative result, reproduced |
| `main` | `()` | runs the self-check |

## 38. `cl10.eta_cp` — the asymmetry as a product of factors

| symbol | signature | returns |
|---|---|---|
| `cp_odd_loop_factor` | `()` | $g_0^4/(16\pi)$ with the $g_0 = 1/8$ reading |
| `eta_cp` | `()` | `{backend, loop_factor, structure, eta_CP_0param, eta_CP_observed, eta_CP_over_observed, g0, cp_odd_loop_factor, sphaleron, delta_B_0param, delta_B_observed, delta_B_ratio, delta_B_v40, delta_B_ratio_v40, one_line}` — three provenances: the loop factor (convention), the sphaleron $28/79$ (an SM input), $\eta_{CP}$ (named) |
| `verify` | `(tol=1e-12)` | `{r, ok}` |
| `main` | `()` | runs the self-check |

## 39. `cl10.four_limits` — the bookkeeping of the residuals

| symbol | signature | returns |
|---|---|---|
| `invariance_g10` | `()` | `{backend, grid, max_mirror_residual_over_grid, max_abs_detU_dev_over_grid, g10_invariant, n_algebraic_invariants, free_angles, verdict}` — the invariants hold for every angle (exact) |
| `unique_10loop` | `(max_n=20)` | `{backend, n_gen, n_gen_is_2x5, 2_is_chiral_split, 5_is_qubit_count, n_best, unique_n_within_1decade_of_both, n_candidates_within_1decade, factor2_test_of_both_empty, coeff_hierarchy, coeff_dm, n_best_combined_dev, n_best_2nd_dev, scan, verdict}` |
| `cp_phase_reduction` | `()` | `{backend, rank1_cp_phases, rank1_cpodd, two_N_2flavor_cp_phases, two_N_3flavor_cp_phases, cpodd_at_delta_0, cpodd_at_delta_pi4, cpodd_at_delta_pi2, minimal_cp_phases_required, delta_B_0param, delta_B_observed, delta_B_ratio, verdict}` |
| `coupling_provenance` | `()` | `{backend, n_qubits, n_space, g10_root_reading, g10_from_alpha10, g10_is_root_reading, g0_root_reading, g0_is_root_reading, qca_cp_even, delta_independent_of_g10, n_free_parameters_in_algebra, n_free_parameters_in_dynamics, mapping_conventions, rejected_identifications, verdict}` |
| `four_limits` | `()` | the assembly: `{backend, unifying_invariant, L1_conventionality, L2_hierarchy_plus_L4_DM, L3_baryogenesis, coupling_provenance, residual_set, total_named_residuals, mapping_conventions, all_resolved, verdict}` |
| `verify` | `(tol=1e-12)` | `{r, ok}` |
| `main` | `()` | runs the self-check |

## 40. `cl10.parameter_audit` — the audit of the label

| symbol | signature | returns |
|---|---|---|
| `kingdon_2route` | `()` | `{backend, kingdon_available, kingdon_sq_minus_1, kingdon_anticommute, kingdon_w2_scalar, pauli_sq_minus_I, pauli_anticommute, pauli_w2_plus_I, two_route_agree}` — the root built twice, by two independent routes |
| `spectral_mirror` | `()` | `{backend, g10, n_levels, omega_min, omega_min_mass_GeV, omega_min_16weyl, levels_16weyl, w_pure_imaginary, commute_Uw, preserves_chiral_split, mirror_resid, mirror_ok, mirror_resid_all_g10, mirror_holds_for_all_g10, g10_grid, zp_16, zp_16bar, det_U_minus_1, sl32}` — the 16/16-bar mirror is exact at every angle |
| `cp_odd_theorem` | `(n4_sample=200000, seed=0)` | `{backend, max_cpodd_2rot, max_cpodd_3rot, max_cpodd_4rot_sample, max_cpodd_b_sandwich, cpodd_democratic_yukawa, n4_sample, all_cp_odd_zero, sakharov, verdict}` — every CP-odd invariant of the on-site module vanishes |
| `g10_status` | `()` | `{backend, g10_0param, n_qubits, mapping_convention, uv_finiteness_fixes_g10, pole_2loop_below_MP_for_all_g10, g10_is_free_t_hooft_angle, grid, verdict}` |
| `hierarchy_status` | `()` | `{backend, v_over_EP_observed, v_over_EP_claimed_0param, match_rel_err, match_is_2digit_coincidence, omega_min, omega_min_mass_GeV, smallest_sub_planck_scale_is, spectrum_exponent_matches_observed, hierarchy_exists_0param, residual}` |
| `dm_status` | `()` | `{backend, soliton_mass_GeV, abundance_gap_0param_over_obs, omega_min_mass_GeV, dm_observed_range_GeV, smallest_sub_planck_scale_above_observed_DM_orders, observed_DM_not_0param, residual}` |
| `provenance_table` | `()` | `{backend, rows, counts_by_provenance}` — every headline number classified |
| `verify` | `(tol=1e-9)` | `{backend, two_route, spectral_mirror, cp_odd_theorem, g10_status, hierarchy_status, dm_status, provenance_table, invariants_ok, limits_ok, ok}` |
| `main` | `()` | runs the audit |


## 41. `cl10.altcell` — the exploration generalizations of the signature and the cell

The machinery the `alt_*` exploration programmes - the W9 interaction, gauged and non-perturbative chains, the cone searches - need, kept
out of `clifford` / `sm_qca` / `supercell` so the published root is never
touched. Two axes. **The signature:** `signature_generators(q)` is the
$\mathrm{Cl}(10-q,q)$ set — the ten qubit generators with the first $q$
multiplied by $i$ — and everything downstream (the plane generators, the
45-plane product, the pseudoscalar, the $16/16$-bar split) is recomputed per
signature. **The cell:** `cell(shape)` is the mixed-radix block, the bit trick
of `supercell` being its $(2,2,2)$ special case (kept as the regression
anchor); `bloch_stages` / `step_stages` / `front_speed_stages` carry the
step-rule axis, and `cone_anisotropy`, `chirality_at`, `fock_self_energy`,
`lowdin_downfold`, `oh_group` support the censuses. `FactorizedCore` is the
`alt_step` cone programmes' joint step - the theta-free bond half `P(q)` and
the Weyl half `U(theta)^ns`, whose commuting sum spectrum replaces the
`(16 n_sub)^3` eigendecomposition. Exploration code: no published number
reads it. The module is **host-only**: the `clifford` and
`supercell` calls it makes are moved back to the host where they enter, so it
returns numpy arrays under either backend. `python -m cl10.altcell` runs the
self-check.

## 42. `cl10.altgpu` — the batched CUDA eigensolver of the interaction census

The heavy loop of V172/V173 is the batched eigendecomposition of the
$128\times128$ Bloch step over a $13824$/$32768$-point momentum grid, and numba
CUDA has no LAPACK. `build_kernels` compiles a cyclic **Jacobi** solver for
**normal** matrices (the Fasi/HSVD variant: each $2\times2$ core is
diagonalized by the unitary polar factor $G = X(X^\dagger X)^{-1/2}$ of its
eigenvector matrix, applied as the similarity $U \leftarrow G^\dagger UG$, with
$V$ accumulating the eigenvectors), one thread block per matrix; `gpu_modes`
runs the grid on the card. Measured against numpy: eigenphase
$1.2\times10^{-11}$, Fock $\Sigma$ $2.9\times10^{-13}$. Exploration code;
`python -m cl10.altgpu` is the check (`./run_altgpu.sh`).

## 43. `cl10.continuum` — the continuum criterion: point groups, the star module, and the shift group of a module

The library side of Section 8.5: which symmetries a gapless locus carries, what
the icosahedral route costs, and what algebra the loop ladder of Section 12.1
lives on. Exact rational arithmetic for the star (`exact` coordinates as
$(a,b)$ meaning $a + b\sqrt5$), float matrices for the group actions.

| symbol | signature | returns |
|---|---|---|
| `cubic_group` | `()` | the $48$ signed permutation matrices of $O_h$ |
| `polyhedral_group` | `(orders=None)` | `{elements, generators, generator_orders, order, element_orders, orders, axis, product_order, trace_defect}` — the $[2,3,5]$ presentation: $60$ elements for the icosahedral group |
| `polyhedral_axes` | `(group, order)` | the group's axes of that order, one per conjugacy-pair of directions |
| `molien_counts` | `(mats, dmax=None)` | the invariant counts of a matrix group by averaging over its elements: $O_h$ gives $[0,1,0,2,0,3,0,4]$, the icosahedral group $[0,1,0,1,0,2,0,2]$ |
| `so3_counts` | `(dmax=None)` | the rotation-group reference counts $[0,1,1,1,1,\dots]$ |
| `first_anisotropic_degree` / `isotropy_forced` | `(counts)` | $4$ for $O_h$ against $6$ for the icosahedral group; whether the quadratic count is exactly $1$ |
| `little_group` | `(q, group, tol=None)` | the stabiliser of a momentum: measured $48$ at $\Gamma$, $8$ at $X$, $4$ at $M$, $6$ at $R$ |
| `locus_census` | `(loci, group=None, dmax=None)` | per locus `{locus, q, little_group_order, invariants, quadratic_invariants, isotropy_forced, first_anisotropic_degree}` — the two jaws of the trap: isotropy is forced only where the touching is quadratic |
| `root_link_group` | `()` | `{order, element_orders, generators, is_2_group}`: the $2048$-element group of the $10$ Clifford involutions, a $2$-group, so it contains no order-$5$ element and no icosahedral subgroup |
| `binary_polyhedral` | `(group_generators)` | the double cover inside the quaternions: `{order, element_orders, unit_norm_defect}` (order $120$ for the binary icosahedral group) |
| `clifford_even_3d` | `()` | `{E, identity, square_to_minus_one, quaternion_products}` — the quaternion units inside $\mathrm{Cl}^0(\mathbb R^3)$ |
| `embed_quaternions` | `(qel, even)` | the quaternion image of a rotation, in the even Clifford basis |
| `icosahedral_star` | `()` | `{directions, exact, axes, norm, phi}` — the $12$ signed directions of the six $5$-fold axes, exact coordinates in $\mathbb Q(\sqrt5)$ |
| `axes_of` / `angle_spectrum` / `frame_defect` | `(star \| vectors)` | the six axes; the multiset of pairwise angles; whether a single frame carries the star |
| `star_invariance` | `(star, group)` | the defect with which the group permutes the star (frame-dependent by nature; the cubic group fails) |
| `module_generators` | `(star=None)` | `{indices, count, exact, vectors, norms}` — the $\mathbb Z$-basis of the module: $6$ independent axes of norm $1.902113$ |
| `module_shortness` | `(star, bounds=(1,2,3,4,5))` | per bound `{bound, generators, min_norm, coefs, displacement}` — the signature of density: the shortest element keeps falling |
| `star_module` | `(star=None, group=None)` | `{direction_count, axis_count, directions, axes, axial, angle_spectrum, fourth_moment_defect, cubic_axes_fourth_moment_defect, frame_defect, module_rank, relation_rank, relations, dense, shortness}` |
| `fourth_moment_defect` | `(vectors)` | the deviation of the direction set's fourth moment from the isotropic one: $0.4$ for the three cubic axes, $2.2\times10^{-16}$ for the star's six |
| `axial_cost` | `(star, target)` | what it costs, in module steps, to move along a direction the cubic lattice offered |
| `shortest_rational_displacement` | `(star, search=1)` | `{displacement, norm, steps}` — the shortest nonzero element of the module at that word bound |

## 44. `cl10.starcell` — the automaton step on the icosahedral shift module

The committed construction (altcell/supercell) is a staggered matching step on a
mixed-radix block cell whose bonds hop along the three cubic axes; the only
cubic input is that list of displacements. This module rebuilds the same step
with the six independent axes of the icosahedral star substituted, on the exact
$6$-dimensional footing the module deserves (rank $6$, dense in $\mathbb R^3$):
$\mathbb Z^6$ is the module's word group and the physical step is its image
under the projection carrying the $i$-th basis vector to the $i$-th module
generator. Used by V214-V217 and V219-V220.

| symbol | signature | returns |
|---|---|---|
| `module_cell` | `(shape=(2,)*6)` | `{kind: "module", shape, ndim, nsub, radix, digit, site_index, sites}` — the mixed-radix periodic cell over the six axes ($64$ sub-sites at the default shape) |
| `module_bonds` | `(c, i, phase)` | the perfect matching of stage $(i,\text{phase})$ as `(a, b, axis, wraps)`; a wrap is the inter-cell hop and carries the Bloch phase |
| `module_stages` | `(shape=(2,)*6, order="phase-major")` | the stagger; `phase-major` alternates intra-/inter-cell as `supercell.STAGES` does, `block-major` runs the cell through - V214 measures both, the choice is part of the construction |
| `stage_matrix` | `(bonds, q, axes, g10, nsub)` | the matching matrix: on each bond the $2\times2$ rotation $[[c, -is\, e^{+i\mathbf q\cdot\mathbf v}], [-is\, e^{-i\mathbf q\cdot\mathbf v}, c]]$ with $\mathbf v$ the module generator, $c = \cos g_{10}$ |
| `bloch_star` | `(cell, U, g10, stages, q, axes)` | the one-step matrix $\prod \operatorname{kron}(S_{\text{stage}}, U)$, built on the device when the engine is on a GPU core ($1024$-dimensional) |
| `eigenphases` / `gap_star` | `(M)` / `(cell, U, g10, stages, q, axes)` | the eigenphases; $\min\lvert\text{eigenphase}\rvert$ |
| `parity` | `(cell)` | the digit-wise reflection $d \to -d \bmod n$: the operator the mirror $\mathbb Z_2$ is compared against |
| `lr_radius` | `(axes)` | the Lieb-Robinson radius per step, the longest module generator: $1.902113$ |
| `word_displacement` | `(coefs, axes)` | `{steps, displacement, norm}` — the physical reach of a word in module steps |
| `star_axes` | `(shape=None)` | the module generators as physical displacements (delegates to `continuum.module_generators`) |

## 45. `cl10.spincell` — the icosahedral star cell with a Clifford bond layer

The step the equivariant census asks for, on the star's own kinematics. The
committed construction (`supercell`) and `starcell` below both factorize — the
step is transport $\otimes$ internal map — and V215/V216 measured the
consequence: the touching at the zone centre has vanishing first order. The
census (`cl10.continuum`) says what a linear isotropic node needs: a little
group $2I$, for which $\mathrm{Hom}_{2I}(D, D\otimes V_1)$ is $1$-dimensional
so that $v\,\boldsymbol\sigma\cdot\mathbf q$ is forced and the first anisotropic
equivariant sits at degree $5$; no inversion symmetry; and a gap subspace on
which the point group acts through a single irrep class. This module delivers
the covariance by making the coin the twelve directed directions of the star and
the bond operator the Hermitian Clifford generator along the direction,
$U_d = -\hat v_d\cdot\vec\sigma$, so that $D(R)\,U_d\,D(R)^\dagger = U_{Rd}$
holds exactly — measured over all $120$ lifts, not assumed. The clock keeps its
committed role: the matching the icosahedral group preserves pairs each
direction with its opposite and rotates each pair by $g_{10}$. The committed
$16$-sector map rides the block layer with an exponent the node arithmetic
selects ($4$, half the order of the map); it is a spectator for the covariance
and carries the quantum numbers. Used by V225-V228. Where `starcell` is the
spinless staggered matching step on the star's shift module ($1024$-dimensional
$S$-matrix), this is the spinor coin cell ($384$ at one layer): the same star, a
different Hilbert space.

| symbol | signature | returns |
|---|---|---|
| `directions` / `opposite_pairs` | `()` | the twelve directed star directions (six axes, both signs) — the coin — and their pairing |
| `star_geometry` | `()` | the six axes as lines, $\sum_d \hat v_d\hat v_d^{\mathsf T}$ over the twelve signed directions, and its defect from $8\,\mathbb 1$: the isotropy the geometry gives for free (the cubic axes give $0.400$) |
| `spinor_lifts` | `()` | the $120$ elements of $2I$ as $2\times2$ unitaries (`order`, `su2`, the element orders, the unitarity defect) |
| `image_group` | `(su2=None)` | the $120$ elements of $2I$ inside $\mathrm{SU}(2)$, with the closure and the $\pm\mathbb 1$ pair checked |
| `axis_permutations` | `(spinor)` | for each lift, the permutation of the six axes and the matching matrix: the two pieces of the covariance statement |
| `pauli` | `()` | the three Pauli matrices the bond operators are built from |
| `bond_ops` | `(kind)` | the twelve bond operators for one of the four assignments (`clifford`, `vector`, `fivefold`, `scalar`); V225 measures which identities each satisfies |
| `shift` / `coin` | `(q, bonds, layers)` / `(g10, layers)` | the hop with Bloch phases over the coin and layer index; the committed rotation $[[c,-is],[-is,c]]$ on each direction pair |
| `internal_ride` | `(U16, legs, layers)` | the committed $16$-sector map riding the block layer at the requested exponent |
| `step` | `(q, U16, g10, kind, layers, legs)` | the one-step unitary: internal ride $\cdot\;$ `kron(coin·shift, spinor)`; $12\cdot$ layers $\cdot\,2\cdot 16$ dimensions |
| `phases` / `gap` | `(M)` / `(q, U16, **kw)` | eigenphases; $\min\lvert\text{eigenphase}\rvert$ and the count of zero eigenphases |
| `tangent` | `(q, U16, g10, kind, layers, legs)` | the three tangent generators of the step at $q$, their norms, commutators and anticommutators on the gap cluster |
| `covariance` / `invariance_defect` | `(kind, layers)` / `(U16, g10, kind, layers, legs)` | the measured hop-covariance and bond-intertwining defects over the $120$ lifts; the defect of step invariance at the zone centre |
| `cluster_frame` / `sector_frame` | `(S, window)` / `(S, sector, window)` | an orthonormal frame of the gap cluster, or of a phase window cut further by the sign of $\sin$: the Hermitian route, because the generic eigenvectors of a degenerate unitary are not stable |
| `cluster_projectors` | `(basis, n_internal, spinor, layers)` | the $2I$ character projectors onto the gap cluster: which irreps the touching carries |
| `module_length` | `()` | the longest star generator, the Lieb-Robinson radius of one step on this module |
| `symmetry_operator` | `(Pm, D, n_internal, layers)` | the permutation-on-axes operator lifted to the cell, against which inversion and rotation covariance are measured |
| `icosahedral_sphere` | `(freq=3)` | a symmetry-uniform triangulation of $S^2$ built from the star's own directions (`vertices`, `triangles`, `faces`, `solid_angle_defect`) |
| `berry_chern` | `(step_fn, q0, mesh, radius, sector, window)` | the first Chern number of the bands inside a sector over the mesh, by Fukui–Harigai–Watanabe link products — no gauge fixing, no eigenvector derivative. Reports `chern`, `integer_defect`, `bands_in_sector`, `min_link_overlap`, `link_rank_ok`, `window`; raises when the sector changes size on the sphere or when no window stays constant |

## Shared constants and switches

| name | value | provenance |
|---|---|---|
| `cl10.u10pred.M_PLANCK` | `1.220890e19` GeV | the lattice cutoff; $M_{\rm GUT} = M_P$ is **not** assumed (PDG 2024 / CODATA 2022) |
| `cl10.u10pred.V_EW` | `246.22` GeV | measured: $v = (\sqrt2 G_F)^{-1/2}$, muon decay (PDG) |
| `cl10.u10pred.M_P_ROBSON` | `0.93827208943` GeV | the proton mass (CODATA 2022) |
| `cl10.dm_wimp.OMEGA_OBS` | `0.120` | the observed cold-DM density $\Omega_c h^2$ (Planck 2018 $0.1200\pm0.0012$) |
| `cl10.eta_cp.DELTA_B_OBS` | `6.12e-10` | the observed baryon asymmetry (Planck 2018) |
| `cl10.u10pred.GEV_TO_S`, `S_TO_YR`, `GEV2_TO_CM2`, `C_CM_S` | `6.582119569e-25`, `3.15576e7`, `3.89379338e-28`, `2.99792458e10` | exact conversions ($\hbar$; the Julian year; $(\hbar c)^2$; $c$) |
| `cl10.unification.M_X_OBS` | `2.0e16` GeV | the observed 1-loop meeting scale (the MSSM reference) |
| `cl10.dm_wimp.M_H` | `125.25` GeV | the Higgs pole mass (PDG 2024) |
| `cl10.dm_wimp.LZ_BOUND` | `1.5e-47` cm$^2$ | the LZ spin-independent direct-detection limit at ~118 GeV |
| `cl10.dm_wimp.F_N`, `SIGMA_FIDUCIAL` | `3.0e-4`, `3.0e-26` cm$^3$/s | the nucleon Higgs coupling (named hadronic input); the thermal-relic reference |
| `cl10.cosmology.G_STAR`, `HUBBLE_COEFF`, `S0_CM3`, `RHO_C_OVER_H2` | `106.75`, `1.66`, `2891.2` cm$^{-3}$, `1.0537e-5` GeV cm$^{-3}$ | the standard-cosmology inputs of the relic / washout checks |
| `cl10.neutrino.M_ATM_EV`, `M_SOL_EV`, `SUM_MNU_MAX_EV` | `0.0501`, `0.0086`, `0.0642` eV | the measured neutrino splittings (NuFIT 6.0, NO) and the cosmological $\sum m_\nu$ ceiling — DESI DR2 BAO + Planck PR4 + CMB lensing, $\Lambda$CDM, 95% CL (the F3 test; the same data gives `0.163` eV in $w_0w_a$CDM, and the Feldman–Cousins boundary limit `0.053` eV) |
| `cl10.dm_wimp.SIGMA_RELIC_TARGET` | `2.2e-26` cm$^3$/s | the relic-closing $\langle\sigma v\rangle$ of the dark-$U(1)$ portal scan (V112); the direct-detection ceiling it is compared against is the single `dm_wimp.LZ_BOUND` above — one observable, one limit |
| `cl10.u10pred.ALPHA_EM_INV` | `127.9` | the electromagnetic coupling at $M_Z$ (PDG) |
| `CFG.desert.theta_root_rad` | `1.7307836435290485` | the root direction of the Wilson term (V65, the $\mathrm{X}_1$ selective root), the angle the ladder is read at; written at the precision it was measured at, because the regression checks test it to $10^{-5}$ |
| `CFG.physics.wimp_window_gev` | `[100, 1000]` GeV | the observed WIMP window the 0-parameter DM scale is compared against (`parameter_audit`, point 4) |
| `CFG.physics.f_10_measured` | `0.019186` | the measured hierarchy input of Section 12.1: $F_{10}$ on the closed traced word census at the rung-10 clock $g = 1/5$ — the clock the data reads, to relative $10^{-5}$ (V85/V113, V206) |
| `CFG.physics.v_over_mplanck_observed`, `hierarchy_coefficient_observed` | derived, never written | `2.0167255e-17` $=$ `physics.v_electroweak_gev / physics.m_planck_gev` and `2.0765641` $=$ that quotient in units of $loop^{10}$: the loader derives both from the two measured constants, so neither ratio has a second copy (V59/V108/V203 compute the same value; V113 reads it). The closed clock reads the coefficient as $2.710$ certified (float64 ladder $2.699$–$2.700$) — $+30\%$, the clock tension — read by `u10pred.hierarchy()` |
| `CFG.units.ev_per_gev` | `1e9` | the eV↔GeV conversion (an exact decimal prefix): the one conversion every eV quote of the ladder uses |
| `CFG.neutrino.obs_sin2_theta12`, `_theta23`, `_theta13` | `0.307`, `0.579`, `0.022` | the observed PMNS mixings the flavour comparisons run against (NuFIT 6.0, NO; the second octant — the octant is unresolved, Section 13.5) |
| `CFG.ckm.obs_sin2_theta12`, `_theta23`, `_theta13` | `0.050`, `0.949`, `0.020` | the measured CKM-side numbers of the same mixing comparisons |
| `cl10.u10pred` (atomic) | `13.605693` eV | the hydrogen Rydberg energy (CODATA), the input of `simulation/hydrogen` |
| `cl10.coupling.G10` | the active branch angle | `closed` $\pi/4$ (default, the root clock), reference `implemented` $1/5$ (convention), `half` $\pi/2$ |
| `cl10.coupling.ALPHA_10` | `1/25` | the rotation-squared reading (convention) |
| `cl10.coupling.G0` | `1/8` | the $g_0 = 1/2^{n_{\rm space}}$ reading (convention) |
| `cl10.unification.M_Z` | `91.1876` GeV | measured (PDG) |
| `cl10.unification.B_SM`, `A_Z`, `B2_SM` | the SM $\beta$-function content and the measured $\alpha_i(M_Z)$ | Standard Model inputs |
| `cl10.qcd` $\sigma$ | `1` in Planck units | an input of the confinement check |
| `cl10.einstein` / `gravity` $G$ | `1` | an input of the gravity sector |

Every number above has exactly one written source: `cl10/config.yaml`, read by
`cl10/config.py` (`CFG`). The environment overrides the file; each variable is
applied once, in the loader, so every consumer sees the same value. The
experiments that need these inputs read the library names for them
(`u10pred.M_PLANCK`, `u10pred.M_X_DESERT`, `u10pred.V_EW`, `u10pred.LOOP`,
`u10pred.ALPHA_EM_INV`, `dm_wimp.OMEGA_OBS`, `dm_wimp.LZ_BOUND`,
`neutrino.M_ATM_EV`, `neutrino.SUM_MNU_MAX_EV`, `eta_cp.DELTA_B_OBS`,
`eta_cp.SPHALERON`, `coupling.G10`),
as do the entry points (the runner's `--jobs`, the `JOBS` and card gate of
`run_experiments.sh`, the gate of `run_figures.sh` and of `make_figures`) and
the figure modules (`make_figures` reads `CFG` directly, never a copy).
Measured constants are written ONCE, at the best published precision: the
file's value is the PDG/CODATA/Planck value, so no experiment carries a second
precision of the same quantity: the fewer-digit duplicates of the Planck
scale, the electroweak vev, the Higgs and proton masses, the neutrino mass
splittings, the entropy density and the relativistic degrees of freedom, the
unit conversions and the coupling values, in the DM, baryogenesis, neutrino
and coupling experiments, and the cross-section conversion of
`cl10.lattice_rg`, are gone. The shared tolerance bands follow the same rule
inside the library: the modules compare against `CFG.tolerances.<band>` (bound
once per module as `_TOL`), never against a literal; an experiment that
measures at its own band states that band itself.

Each **simulation** owns
its scene panel (`simulation/*/config.yaml`), which the code now reads with no
fallback default: a key the code requires must be in the panel.

| variable | overrides | effect |
|---|---|---|
| env `CL10_BACKEND` | `backend.force` | `numpy` / `cpu` forces the CPU core; `cupy` / `jax` pins a GPU core — a pinned GPU core skips the gate and never falls back, so a run that cannot get the card fails loudly |
| env `CL10_GPU_MIN_FREE_MIB` | `backend.gpu_min_free_mib` | the free-memory gate for the GPU in auto mode (default 2048) |
| env `CL10_ANGLE` | `angle.branch` | `closed` (default) / `implemented` / `half` (aliases `pi/4`, `pi/2`, `1/5`) |
| env `CL10_NJOBS` | `pool.njobs` | the worker count: `parpool`, the runner's `--jobs`, `run_experiments.sh` |
| env `CL10_SMOKE` | `smoke` | `1` shrinks the grids, pools and bond counts (`cl10.config.smoke()`) |

## The runner contract

Every experiment writes `results/<name>.json` with at least:

| key | meaning |
|---|---|
| json `pass` | the verdict |
| json `status` | one-line statement of what was measured |
| json `tier` | 0–3 the hard tier (the root, the kinematics, the dynamics, the matter, the $5/3$); 4–6 supporting |
| json `property` | the paper property it validates |
| json `residual` | `true` when the experiment reports an open problem |
| json `wall_time` | seconds |
| json `backend` | optional: the core that computed the verdict (`numpy` / `cupy` / `jax`); the rollup counts it |

`run_validation.py` aggregates them into `results/validation_summary.json` and
prints the hard-tier verdict, the backend coverage of the run (`backend_parent`,
`backend_requested`, `n_backend_gpu` of `n_backend_reported`) and the named
residuals. A subset run folds the untouched experiments in from the verdicts on
disk, so the rollup stays complete; a folded verdict is reported **STALE** (and
counted as a failure) when `results/crashes/<Vxx>.log` is newer than it, and the
OFF-CORE line lists the experiments this invocation ran on a different core than
the parent. V186/V187 pin `CL10_BACKEND=numpy` themselves: the cone programme is
CPU work (independent 128×128 Floquet eigendecompositions, ~40× faster on 12 CPU
workers than on a 3090), so they appear OFF-CORE by construction in a GPU run.
