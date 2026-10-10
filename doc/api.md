# CL10 — API Reference

The complete **Python** API of the `cl10` package and the `run_validation.py`
runner.  Every symbol is **0-parameter** (no fit) and **backend-agnostic**
(the active backend is selected once at import; the physics is identical on
CuPy / JAX / NumPy).

```
import cl10

cl10.engine        # the single computation core (Section 1)
cl10.clifford      # the Cl(10) root (Section 2)
cl10.spacetime     # the (1,3) inside the GUT (Section 3)
cl10.matter        # the 16 Weyl + the anomalies (Section 4)
cl10.gauge         # the 45 + the 5/3 + the GUT coupling/running + gauge-gravity (Section 5)
cl10.qca           # the metric + the QCA (Section 6)
cl10.u10pred       # the headline numbers (Section 7)
cl10.walls         # the 3 reduced walls W1/W2/W3 + the baryon asymmetry (Section 8)
```

`cl10.__version__` is the package version.

---

## 1. `cl10.engine` — the single computation core

The **one** backend-aware core.  All other modules and all experiments go
through it, so the whole theory runs on the GPU out of the box.

### backend

| symbol | meaning |
|---|---|
| `engine._BACKEND` | the active backend: `"cupy"` \| `"jax"` \| `"numpy"` (chosen at import: CuPy → JAX → NumPy) |
| `engine.Array` | the array type: `cp.ndarray` \| `jax.Array` \| `np.ndarray` |
| `engine._HAS_EINOPS` | `True` if `einops` is importable (used for reshapes) |
| `engine._HAS_KINGDON` | `True` if `kingdon` is importable (GPU physics utilities) |
| `engine.backend_info()` | `{backend, einops, kingdon, array}` (a dict) |

### array primitives (all backend-aware)

| symbol | signature | returns |
|---|---|---|
| `Array_` | `(*args, **kwargs)` | a new array on the active backend (like `np.array`) |
| `zeros` | `(shape, dtype=None)` | a zero array; `shape` is a tuple (e.g. `(L,L,L)`) or an int (default `dtype=complex`) |
| `kron` | `(a, b)` | the Kronecker product |
| `matmul` | `(a, b)` | the matrix product |
| `real` / `imag` | `(x)` | the real / imaginary part (backend-correct) |
| `absmax` | `(x)` | `float` — the max absolute value |
| `asnumpy` | `(x)` | the array on the **host** (NumPy) — for printing / JSON |
| `einops_rearrange` | `(x, pattern_in, pattern_out)` | an `einops` rearrange (NumPy fallback if available) |

### lattice / QCA primitives

| symbol | signature | returns |
|---|---|---|
| `laplacian_3d` | `(u)` | the 6-neighbour discrete Laplacian of a `(L,L,L)` field (the metric operator) |
| `qca_step` | `(fields, unitaries)` | one QCA step: apply the product of local unitaries on the last (internal) axis |
| `spectral_laplacian` | `(k, L)` | the spectral (momentum-space) Laplacian eigenvalue on an `L`-torus |

**Note.**  `asnumpy` is required before printing or JSON-serializing a GPU
array (CuPy forbids implicit NumPy conversion).

---

## 2. `cl10.clifford` — the $\mathrm{Cl}(10)$ root

| symbol | signature | returns |
|---|---|---|
| `generators` | `()` | a list of **10** `(32,32)` complex matrices — the 10 mutually-anticommuting $h_a$ (the symplectic Pauli pattern) |
| `verify` | `(tol=1e-12)` | `{n, square_to_I, anticommute, ok}` — the 0-parameter check that the 10 anticommute and square to $I$ |
| `basis` | `()` | a list of **1024** `(32,32)` matrices — the full algebra basis (all ordered products) |
| `pseudoscalar` | `()` | the `(32,32)` **anti-Hermitian** volume element $w = h_1\cdots h_{10}$ (the 16/16-bar split) |
| `weyl_spinor` | `()` | `(plus, minus)` — the `(32,16)` and `(32,16)` **16 Weyl / 16-bar** eigenvectors of $w$ |
| `main` | `()` | runs the V0 self-check (prints the verdict) |

---

## 3. `cl10.spacetime` — the $(1,3)$ inside the GUT

| symbol | signature | returns |
|---|---|---|
| `causal_generators` | `()` | `(T, [S1,S2,S3])` — the 4 $(1,3)$ generators **inside** the GUT: $T = h_1$ ($+$), $S_a = i h_{a+1}$ ($-$) |
| `verify` | `(tol=1e-12)` | `{t2, s2, ts, ss, ok}` — the 0-parameter Minkowski check (all four relations $< \varepsilon$) |
| `light_cone` | `(T, S)` | a callable `f(v0, v) -> v0**2 - sum(v_a**2)` (the Minkowski norm) |
| `main` | `()` | runs the V1 self-check |

---

## 4. `cl10.matter` — the 16 Weyl + the anomalies

| symbol | signature | returns |
|---|---|---|
| `SM_FIELDS` | (module constant) | `{name: (su3_type, su2_type, Y, n_weyl)}` — the **0-parameter** charges of the 6 SM fields (the input to the anomalies) |
| `weyl16` | `()` | `(plus, minus)` — the 16 Weyl / 16-bar (delegates to `clifford.weyl_spinor`) |
| `dimension_check` | `()` | `int` — the total Weyl dimension (= 16) |
| `bl_charge` | `(which="16")` | `int` — the 0-parameter **B–L** charge: $-1$ on the 16, $+1$ on the 16-bar |
| `anomalies` | `()` | `{su3_cubed, su2_cubed, u1_cubed, grav2_su3, grav2_su2, grav2_u1, u1_su3_sq, u1_su2_sq}` — the **8** SM gauge-anomaly coefficients (all **0** = consistent) |
| `mass_generation` | `()` | `{off_diagonal, diag_16, diag_16bar, A_16_16bar, B_16bar_16, ...}` — the **0-parameter 10-Higgs**: the 10 (grade-1) is **off-diagonal** (16$\leftrightarrow$16-bar) $\Rightarrow$ the Dirac mass term $16\,10\,16\bar\to 1$ (the exact values are named) |
| `main` | `()` | runs the V2 self-check |

**`anomalies`** is the **hard consistency requirement**: a chiral gauge theory is
inconsistent unless all eight coefficients cancel.  `V11_anomaly` checks them.

---

## 5. `cl10.gauge` — the 45 + the $5/3$ + gauge-gravity

| symbol | signature | returns |
|---|---|---|
| `M_GUT` | (module constant) | `float` — the GUT scale $= 1.22\times10^{19}$ GeV (the 0-parameter Planck scale) |
| `generators_45` | `()` | a list of **45** `(32,32)` matrices — the $h_a h_b$ for $a<b$ (the $\mathrm{Spin}(10)$ adjoint) |
| `rank_check` | `()` | `int` — the rank of the 45 (the 45-D adjoint) |
| `five_thirds` | `()` | `float` — the **0-parameter** $\mathrm{U}(1)$ normalization $g_1^2/g'^2 = 5/3$ (from $\mathrm{Tr}(\mathcal{T}^2) = 5/24$) |
| `sm_subgroup` | `()` | `{n_sm, n_gluon, n_w, n_b, n_heavy}` — the 12 SM gauge bosons inside the 45 (8 gluon + 3 W + 1 B) |
| `gauge_gravity_overlap` | `()` | `{n_causal, n_gauge, n_generators, one_object}` — the $(1,3)$ (4) and the 45 in the **one** $\mathrm{Cl}(10)$ |
| `running` | `()` | `{alpha_10, g_gut, m_gut_0param, five_thirds, b1, b2, b3, unification_spread, cleanly_unifies, residual}` — the **0-parameter** GUT coupling $\alpha_{10}=1/25$ + the $5/3$ + the 1-loop $\beta$; the **1-loop SM matching** (the SM does not cleanly unify) is the **named** residual |
| `main` | `()` | runs the V3/V12/V13 self-check |

---

## 6. `cl10.qca` — the metric + the QCA

| symbol | signature | returns |
|---|---|---|
| `METRIC` | (module constant) | the $(1,3)$ metric $\mathrm{diag}(1,-1,-1,-1)$ |
| `metric_field` | `(L, h=0.0)` | the `(4,4,L,L,L)` metric field on an $L^3$ torus (the perturbation $h$ on flat) |
| `wave_step` | `(u, dt=0.5)` | one step of the free wave equation (the linearised **graviton**) |
| `qca_step_internal` | `(gates, fields)` | one QCA step on the internal (Weyl) axis of `(L,L,L,n)` fields |
| `unitary_check` | `(U, tol=1e-12)` | `bool` — whether a gate $U$ is unitary ($UU^\dagger = I$) |
| `main` | `()` | runs the V6 self-check |

---

## 7. `cl10.u10pred` — the headline numbers

| symbol | signature | returns |
|---|---|---|
| `proton_decay_lifetime` | `()` | `float` — the **0-parameter** proton lifetime $\tau_p \sim M_{\rm GUT}^4/(\alpha_{10}^2 m_p^5)$ in **years** ($\approx 4\times10^{47}$; the proton appears stable) |
| `seesaw_mass` | `()` | `float` — the **1-step** seesaw mass $m_\nu \sim m_D^2/M_{\rm GUT}$ in **eV** with $M_N = M_{\rm GUT} = M_P$ ($\approx 1.2\times10^{-6}$; the ladder version of V57, $M_N = M_P(1/16\pi)^3$, gives $0.146$ eV — the two differ by the choice of $M_N$, which is not fixed) |
| `family_candidate` | `()` | `{candidate, from_spatial, from_internal_half, derived}` — the 0-parameter **3-family** candidate (`derived=False`: a named wall) |
| `main` | `()` | runs the V14–V16 self-check |

The 0-parameter constants: `M_PLANCK = 1.22e19` GeV, `M_P_ROBSON = 0.938`
GeV, `ALPHA_10 = 1/25`, `V_EW = 246.0` GeV.

---

## 8. `cl10.walls` — the three reduced walls (W1/W2/W3) + the baryon asymmetry

The mechanisms that reduce three of the open residuals (see
`doc/publication.md`, Section 15, Table 6):

| symbol | signature | returns |
|---|---|---|
| `w1_cosmological_constant` | `()` | `{n_even, n_odd, cancellation, lambda_0param, lambda_obs, residual}` — the **W1**: the $\mathbb{Z}_2$ $512/512$ **cancellation** $\Rightarrow$ $\Lambda = 0$ at 0-parameter (the **s-mechanism**); the exact $10^{-122}$ the named residual |
| `w2_seesaw` | `()` | `{g0, m_D_GeV, M_N_GeV, m_nu_2step_eV, m_nu_1step_eV, m_nu_obs_eV, off_1step, off_2step, right_order, residual}` — the **W2**: the **2-step** seesaw ($M_N$ intermediate, $g_0 v$ Dirac) $\Rightarrow$ $m_\nu \sim 10^{-2}$ eV (the right order) |
| `w3_families` | `()` | `{n_time, n_space, n_fam_candidate, n_observed, match, derived, residual}` — the **W3**: the **3** families identified with the **3** spatial directions of the $(1,3)$ (one number stated three ways; a hypothesis, not a derivation) |
| `baryon_asymmetry` | `()` | `{baryon_number_violated, c_or_cp_violated, out_of_equilibrium, all_sakharov, generated, residual}` — the 0-parameter QCA satisfies **none** of the 3 Sakharov conditions $\Rightarrow$ the baryon asymmetry is **not** generated (the **named** residual, V22) |
| `main` | `()` | runs the three wall mechanisms (V17–V19) and prints the verdict |

The 0-parameter constants: `M_PLANCK = 1.22e19` GeV, `V_EW = 246.0` GeV,
`G0 = 1/8` (the GUT coupling), `M_N = 1e14` GeV (the 2-step intermediate
scale).

---
## 9. `cl10.sm_qca` — the full SM QCA (the 0-parameter unitary QCA)

The **full Standard Model as a 0-parameter unitary QCA** (V23): the QCA step
is the product of the **45 Spin(10) plane rotations** $T_{ab} = i\,h_a h_b$
(Hermitian, $T_{ab}^2 = I$) at the **0-parameter** angle $g_{10} = 1/5$,
**exactly unitary** (no matrix exponential, no fit), run on the **GPU** engine.

| symbol | signature | returns |
|---|---|---|
| `plane_generators` | `()` | the 45 Spin(10) plane generators $T_{ab} = i h_a h_b$ (a list of 45 $(32\times32)$, Hermitian, $T^2=I$) |
| `plane_rotation` | `(T, theta)` | $\cos\theta\, I - i\sin\theta\, T$ (the exact plane rotation, $T^2=I$) |
| `gut_qca_step` | `(planes, theta)` | the GUT QCA step $U = \prod_{a<b}[\cos\theta\, I - i\sin\theta\, T_{ab}]$ (exactly unitary) |
| `unitary_err` | `(U)` | the unitarity residual $\|U U^\dagger - I\|_{\max}$ |
| `norm_preservation` | `(U, n_steps, L)` | $(n_0, n_1, \text{rel\_err})$ — the norm before/after $n$ steps on a random $32$-D field (GPU) |
| `sm_qca` | `(n_steps, L)` | `{n_planes, g10, plane_square_to_I_err, step_unitary_err, norm_before, norm_after, norm_rel_err, backend, causal_radius}` — the full SM QCA verdict (V23) |
| `main` | `()` | runs the V23 self-check |

---

## 10. `cl10.cpt` — the CPT theorem (the 16/16-bar are CPT conjugates)

The **CPT theorem** (V24): the CPT operator $C = g_1$ (the time direction) is a
**Hermitian involution** ($C^2=I$) that **anticommutes** with the pseudoscalar
$w$ (so it maps the 16 to the 16-bar); the 45 GUT generators **commute** with
$w$ (the chiral CPT split is preserved).  Verified on the **GPU**.

| symbol | signature | returns |
|---|---|---|
| `pseudoscalar` | `()` | the volume element $w = g_1 \\cdots g_{10}$ (anti-Hermitian, $w^2 = -I$) |
| `cpt_operator` | `()` | the CPT operator $C = g_1$ (the time direction, a Hermitian involution) |
| `projectors` | `(w)` | $(P_{16}, P_{16\\text{bar}}) = (\\tfrac{1}{2}(I - i w), \\tfrac{1}{2}(I + i w))$ (the 16/16-bar projectors) |
| `cpt_theorem` | `(seed)` | `{C_involution_err, C_unitary_err, C_w_anticommute_err, map_16_to_16bar_err, gut_commute_w_err, n_matter, n_antimatter, cpt_equal_norm, c16_into_16bar, backend}` — the CPT relations + the GPU in-action check (V24) |
| `main` | `()` | runs the V24 self-check |

---


## 11. The four deepening modules (V25–V28)

The four **deepening** modules extend the construction to four further
hallmarks expected of a unification (all run on the **GPU** engine).

**`cl10.holo` — the holographic area law (V25).**  The QCA vacuum entropy of
a region of side $L$ scales as $S \propto L^2$ (the **surface area** $6L^2$),
**not** $L^3$ (the volume); the entropy per unit area is the argument
`per_plaquette = 0.25`, i.e. the Bekenstein–Hawking $S = A/4$ at $G=1$ — an **input** of the model, not a derived density; the holographic ratio $S/V \propto 1/L$.

| symbol | signature | returns |
|---|---|---|
| `area_law` | `(L_values, per_plaquette=0.25)` | `{L_values, per_plaquette, S, R, area_law_exponent, holographic_ratio_exponent, backend}` — the 0-parameter area law (V25) |

**`cl10.family` — the 3-family construct (V26).**  The **3** observed families
are identified with one number appearing in **three** ways (the 3 space
directions; the 3 two-planes $\binom{3}{2}$; the rank-3 spatial metric) — three
names for one number, not three derivations (W3); the 16 Weyl is one generation
$\to$ 3 generations $= 16 \times 3 = 48$ Weyl.

| symbol | signature | returns |
|---|---|---|
| `family_construct` | `()` | `{n_space, n_planes_C32, n_rank, n_family_candidate, n_families_observed, match, all_three_ways_3, n_weyl_1gen, n_weyl_3gen, residual}` — the 0-parameter 3, three ways (V26) |

**`cl10.qcd` — QCD confinement (V27).**  The colour count $N_C = 3$ (read off
the decomposition of the 45; the $8 = N_C^2 - 1$ gluons are then a consistency
check); the **area law** $\langle W(C)\rangle \propto \exp(-\sigma \cdot
\mathrm{Area})$; the string tension used is $\sigma = 1$ in Planck units
($= 1.5\times10^{38}$ GeV$^2$, against the physical $0.18$ GeV$^2$ — a 38-order
mismatch, the named residual); the **linear** potential $V(r) = \sigma r$ (in
contrast to the Coulombic $1/r$).

| symbol | signature | returns |
|---|---|---|
| `color_count` | `()` | `(n_c, n_gluon, ok)` — the 0-parameter $N_C=3$ |
| `string_tension` | `()` | the 0-parameter string tension $\sigma$ |
| `confining_potential` | `(sigma, r_values)` | `[{r, v_qcd, v_qcd_over_r, v_em}, …]` — the linear QCD potential vs the Coulombic EM (V27) |

**`cl10.yukawa` — the Yukawa / mass-hierarchy (V28).**  The **10-Higgs
16→16-bar** (the grade-1 GUT vector is **off-diagonal** in the 16/16-bar split)
gives the 0-parameter Dirac mass term; the **3×3** mass matrix (the 3
generations) is the **democratic** base $|1\rangle\langle1|$: **rank 1**,
spectrum $\{3,0,0\}$ (one heavy + a 2-fold light level).  A rank-1 3×3 is
**singular**: the three distinct masses are the **named residual** (the 10-Higgs
vev), not a 0-parameter statement.

| symbol | signature | returns |
|---|---|---|
| `yukawa_matrix` | `()` | `{n_gen, n_dim, rank, spectrum, n_distinct_levels, n_heavy, n_light_2fold, democratic, singular, off_diagonal_Dirac, diag_zero, mass_term, residual}` — the 0-parameter **democratic rank-1** 3×3 Dirac mass matrix, spectrum {3,0,0} (V28; delegates to `yukawa_3x3.mass_matrix_3x3`) |


## 12. The construction modules (V36–V51)

These modules build the objects the paper's Sections 8–9 use: the local causal
QCA, the discrete Einstein–Hilbert action, the coupling readings, the
3-family step, the Yukawa base, leptogenesis, the S-matrix, the lattice RG and
the DM scale.  Each function is **backend-agnostic** (the GPU engine), takes no
fit unless a coefficient is explicitly named, and the matching experiment
writes a JSON verdict.

**`cl10.local_qca` — the local QCA (V36).**  The **local** QCA step on
the $(L,L,L)\times32$ lattice: the 45 Spin(10) plane rotations **per cell**
plus the **0-parameter exchange (swap) gates** between **adjacent** cells
(the 't Hooft swap-gate structure).  The **causal radius 1 is
implemented**, not an axiom: a single-cell excitation after $n$ steps is
supported **only** inside the $n$-cell light cone (the front speed is
**exactly 1**).

| symbol | signature | returns |
|---|---|---|
| `local_qca` | `(L, n_steps, seed=0)` | `{unitary_residual, support, front_speed, structural_ok, causal_radius_1, …}` — the local step + the realized light cone (V36) |

**`cl10.einstein` — gravity.**  The metric dynamics of the $(1,3)$
sub-algebra: the massless spin-2 (2 polarizations, $c=1$, UV-bounded), the
Einstein–Hilbert action (the computed Riemann curvature, the non-linear
soliton, the Schwarzschild black hole), and — **V37** — the
**derivation** in the direction **dynamics ← action**:

| symbol | signature | returns |
|---|---|---|
| `einstein_hilbert` | `(r=2.0, tol=1e-4)` | `{R_sphere, R_flat, G, coupling=1/16π, n_deriv, einstein_hilbert_ok, …}` — the 0-parameter EH action on prescribed metrics (V31) |
| `eh_quadratic_action` | `(u)` | the 0-parameter **discrete** quadratic EH action on the Planck lattice ($G=1$, $c=1$, $dt^2=1/4$) |
| `einstein_hilbert_from_dynamics` | `(L=8, T=24, n_probe=16, seed=0, tol=1e-6)` | `{elm_on_qca_trajectory, stationary_ok, elm_off_dynamics, selects_ok, coeff_ok, continuum_ratio, continuum_ok, ok, …}` — the **Euler–Lagrange equation** of the 0-parameter action **is** the QCA metric dynamics (V37): the numerical functional derivative vanishes on the QCA trajectory, an off-dynamics field has a non-zero ELM, the coefficient is $1/16\pi$ ($G=1$), and the discrete dispersion gives the massless $c=1$ continuum |

**`cl10.unification` — the coupling readings from the root (V38).**

| symbol | signature | returns |
|---|---|---|
| `g10_from_root` | `()` | `{n_module=32, n_qubits=5, g10=1/5, alpha_10=1/25, n_space=3, g0=1/8, …}` — the proposed readings $g_{10}=1/n_{\rm qubits}$ (the 5 qubits of the 32-D module) and $g_0 = 1/2^{n_{\rm space}}$ (the 3-space).  Both are **conventions** (Table 5 of the paper): V56 shows the data prefer $n=6$ over the root $n=5$ |
| `running_down` | `()` | `{alpha_10_at_M_P, couplings_at_M_Z, measured_at_M_Z, ratio_to_measured, u1_up_pole_scale, su3_down_pole_scale, residual, …}` — the 1-loop running **down** from $\alpha_{10}=1/25$ at the cutoff to $M_Z$: $\alpha_1(M_Z) = 0.0197$ ($1.17\times$), $\alpha_2(M_Z) = 0.195$ ($5.8\times$), and an $\mathrm{SU}(3)$ Landau pole at $2.2\times10^9$ GeV.  This is a **failure of the coupling values** (W4, falsifier F2), not a loose match |

**`cl10.family` — the 3 families inside the dynamics (V39).**

| symbol | signature | returns |
|---|---|---|
| `family_qca` | `()` | `{n_module_3gen=96, n_weyl_3gen=48, n_cpt_3gen=48, unitary_err, family_symmetry_off_block, …, tribimaximal_sin2, residual, …}` — the **3-generation QCA step**: the 45 act **identically on each family** (the family symmetry; the 96×96 is exactly unitary, the per-family 16/16-bar split is preserved) + the tribimaximal mixing (V29; only $\sin^2\theta_{12}=1/3$ is a derived number) |

**`cl10.yukawa_3x3` — the 3-gen Yukawa mass matrix: the democratic base (V45, W6).**

| symbol | signature | returns |
|---|---|---|
| `w6_verdict` | `()` | `{g10, one_object (32-level spectrum + the 96-level 3-family), family_symmetry (the 3-fold symmetry: the 3 spatial directions S_a are equivalent -> the democratic base is DERIVED, not asserted; S_a^2=-I, {S_a,S_b}=0, equal norm), mass_matrix (the democratic 3x3: eigenvalues {3,0,0}, 1 heavy + 2-fold light, rank 1), mixing (the TBM, which **diagonalizes** |1><1|; honest split: sin^2_12 = 1/3 is derived 0-parameter from the unique heavy eigenvector, while sin^2_23 and sin^2_13 sweep [0, 2/3] under the O(2) freedom of the 2-fold light subspace), hierarchy (the 0-parameter 2-level structure; the 3-level is the named residual), higgs_dirac (the computed 16->16-bar Dirac term, rank-16, off-diagonal), …}` — the **W6** projection of the one object: the **3-gen Yukawa mass matrix** (the **derived** democratic base) + the **TBM mixing** (only sin^2_12 = 1/3 is 0-parameter as a number) + the **2-level structure** + the **0-parameter rank-16 10-Higgs Dirac term**; the 3-gen hierarchy (the 10-Higgs vev), sin^2_23, sin^2_13 (the O(2) choice) and the non-trivial CKM are the named residual |

**`cl10.baryogenesis` — GUT leptogenesis (V40, W5).**

| symbol | signature | returns |
|---|---|---|
| `leptogenesis_estimate` | `()` | `{M_N_GeV=1e14, gamma_over_h, g0, eps_loop=g0^4/16π, eta_cp (named), sphaleron_factor=28/79, delta_b, observed_delta_b, ratio_to_observed, …}` — the baryon asymmetry at the **right order** as a product of an SM input, a guessed loop rate and a loop factor.  **Honest:** this V40 estimate hard-codes $M_N = 10^{14}$ GeV and its out-of-equilibrium kinematics is superseded — at the ladder value $M_N = 9.6\times10^{13}$ GeV (V57) the decay is **in** equilibrium ($\Gamma_N/H \simeq 4.6$ at $y=1/8$), so the Sakharov condition fails at the model's own couplings |

**`cl10.smatrix` — the S-matrix of the QCA (V41).**

| symbol | signature | returns |
|---|---|---|
| `qca_smatrix` | `()` | `{g10, n_planes=45, s1 (32×32), s2 (1024×1024), unitary_residual, probability_conservation, max_offdiag, …}` — the **'t Hooft program** completed: the QCA step **is** the 0-parameter S-matrix (the 1- and 2-particle transition amplitudes; the only input $g_{10}=0.2$) |

**`cl10.lattice_rg` — the QCA lattice RG: the S-matrix to the cross-sections (V44, W9).**

| symbol | signature | returns |
|---|---|---|
| `lattice_rg` | `()` | `{g10, one_particle (P1), two_particle (P2, the 1024x1024 S^(2) = U (x) U), uv_bound (sigma_max = 4 pi / M_P^2), cross_section (sigma_0 = P2 / M_P^2, sigma_pert, R, elastic, optical_theorem), dispersion (the 32-level spectrum H = -i log U, UV-bounded, zero_point), froissart (the 0-parameter Froissart-Martin bound, m_min = M_P: the QCA constant cross-section is below the ln^2 growth -> strongest UV-finiteness), uv_finite, …}` — the **'t Hooft program** at the cross-section level: the QCA is **UV-finite** (the Planck cutoff bounds every cross-section, `sigma_0 < sigma_max`), the **0-parameter** 2->2 cross-section `sigma_0 = P2 / M_P^2`, the **0-parameter dispersion** (the 32-level spectrum, UV-bounded `|omega| <= pi`, no 1-loop pairing), the **0-parameter optical theorem** (the inelastic `P2` **is** the cross-section), the **0-parameter Froissart-Martin bound** (the QCA's *constant* cross-section is **below** the `ln^2` growth — the **strongest** UV-finiteness), and the **0-parameter** QCA/perturbative ratio `R = 4 pi P2 / g_10^4`; the **named residual** (W9) is the lattice RG matching to the measured cross-sections |

**`cl10.cosmology3` — the unified cosmology + the DM scale (V42, W7).**

| symbol | signature | returns |
|---|---|---|
| `dm_mass_scale` | `()` | `{E (the soliton energy, Planck units), m_DM_GeV, sigma, Omega_0param, Omega_obs, abundance_gap, residual, …}` — the **quantified** dark-matter residual: the Planck-mass soliton (stable, massive, neutral) and the computed 0-parameter abundance gap (the sub-Planck scale is the named residual) |

**`cl10.qca_spectrum` — the one object and its invariants (V43).**

| symbol | signature | returns |
|---|---|---|
| `levels` | `(g10=None)` | the $32$ sorted phases of $U$ (the levels of $H = -i\log U$), in $(-\pi,\pi]$ |
| `qca_spectrum` | `(g10=None)` | `{g10, n_levels=32, spectrum, zero_point (=0: the ±-paired 0-point), det_U_minus_1 (=0: the SL(32) identity), sl32, pairing_resid, max_trace_T45 (=0: the 45 are traceless), cpt_0point_16, cpt_0point_16bar, commutes_with_pseudoscalar, max_level, …}` — the **one** object (the 32-level QCA spectrum) and its invariants (the 0-point $=0$, $\det U=1$ for every $g_{10}$, the chiral split) |
| `verify` | `(tol=1e-10)` | `{spectrum, ok}` — the checks: the 0-point $=0$, the SL(32) identity, the $\pm$-pairing, the chiral split |

**The four remaining-gaps modules (V46–V49).**  Four completions of the
remaining gaps (W1, W4, W9, W6): the **4-D non-linear
gravity**, the **2-loop UV-finite** running (W4), the **3- and
4-particle** amplitudes (W9), and the **quark CKM on the democratic
base** (W6 — the 2-level mass structure is fixed by the democratic base; the computed CKM is the
identity, see `quark_ckm` below).

**`cl10.nonlinear_gravity` — the 4-D non-linear gravity (V46, W1).**

| symbol | signature | returns |
|---|---|---|
| `nonlinear_gravity` | `(H=None)` | `{backend, H, einstein_space_err, g_mu_nu_minus_3h2_g, ricci_scalar_minus_12h2, gamma_gamma_max, nonlinear_present, lambda_0param_GeV4, lambda_obs_GeV4, gap, residual, ok}` — the **0-parameter 4-D non-linear gravity**: the full 4-D Riemann + Einstein tensor of the 0-parameter de Sitter metric (the Einstein-space identity $R_{\mu\nu\rho\sigma}=H^2(g g - g g)$, $G_{\mu\nu}=-3H^2 g_{\mu\nu}$, the non-zero $\Gamma\Gamma$ self-gravity, the 0-parameter vacuum energy $\Lambda\sim M_P^4$; the exact $10^{-122}$ is the W1 named residual) |
| `riemann_4d` | `(H, x0=0)` | the full 4-D Riemann tensor $R^\mu_{\nu\rho\sigma}$ of the de Sitter metric (analytic) |

**`cl10.gravity` — the COMPLETE 4-D non-linear gravity (V50, C2, W9).**

The 4-D non-linear gravity **completed**: the general Riemann / Ricci /
Einstein / Kretschmann are **computed by central finite differences from an
*arbitrary* metric** (no hard-coded Christoffels, not restricted to one
solution), **2-route verified** (the computed de Sitter Riemann matches the
independent analytic one — `cl10.nonlinear_gravity` — to $O(h^2)$; the
Schwarzschild Ricci is $0$ and the **computed** Kretschmann is $48M^2/r^6$ —
the **black hole recovered by computation**; the flat metric gives
$R=G=K=0$), the **complete** 0-parameter discrete Einstein–Hilbert action
$S = \frac{1}{16\pi}\sum_x \sqrt{-\det g}\,(R - 2\Lambda)$ ($G=1$) whose
functional derivative **is** the non-linear Einstein equation
(dynamics $\leftarrow$ action), and the **non-linear** (state-dependent)
QCA metric evolution (the 't Hooft non-linear completion: unitary/UV-finite,
causal $c=1$, superposition fails).

| symbol | signature | returns |
|---|---|---|
| `H_FD` / `G10` | (module constants) | the 0-parameter finite-difference step ($10^{-3}$, a *numerical* convergence parameter) / the GUT angle ($g_{10} = 1/5$) |
| `christoffel_4d` | `(g, dg)` | $\Gamma^\rho_{\mu\sigma}$ — the 4-D Christoffel symbols from a 4×4 metric and its first derivatives (a 4×4×4 array) |
| `riemann_4d` | `(mfn, x, h=H_FD)` | the **full 4-D Riemann** $R^\rho{}_{\sigma\mu\nu}$ (4×4×4×4) **computed by central finite differences from an *arbitrary* metric** `mfn(x) -> 4x4` |
| `ricci_4d` | `(R)` | $R_{\mu\nu} = R^\rho{}_{\mu\rho\nu}$ |
| `ricci_scalar_4d` | `(Ric, gin)` | $R = g^{\mu\nu}R_{\mu\nu}$ |
| `einstein_tensor_4d` | `(Ric, g, scalar)` | $G_{\mu\nu} = R_{\mu\nu} - \tfrac12 g_{\mu\nu} R$ |
| `kretschmann_4d` | `(R, g)` | $K = R^{mnpq}R_{mnpq}$ (the coordinate-independent curvature invariant) |
| `gamma_gamma_4d` | `(mfn, x, h=H_FD)` | $\max\lvert\Gamma\Gamma\rvert$ — the non-linear (cubic) self-gravity part of the Riemann (0 for flat, $>0$ for a curved metric) |
| `metric_flat` / `metric_de_sitter` / `metric_schwarzschild` / `metric_2sphere` | `(x)` / `(H, x)` / `(M, x)` / `(r, x)` | the 0-parameter metric factories (the known solutions, each a 4×4) |
| `discrete_eh_action` | `(mfn, grid, Lambda, h=H_FD)` | $S = \frac{1}{16\pi}\sum_x \sqrt{-\det g(x)}\,(R(x)-2\Lambda)$ (the complete 0-parameter discrete EH action, $G=1$) |
| `dynamics_from_action` | `(mfn, Lambda, x0, h=H_FD, eps=1e-4)` | the **functional derivative** $\delta S/\delta g$ (the numerical Euler–Lagrange) — the non-linear Einstein equation: vanishes on a solution (flat), nonzero on a curved metric (the self-gravity source) |
| `nonlinear_qca` | `(L=16, n_steps=12, x0=None, seed=0, lam=G10, tol=1e-12)` | the **non-linear** (state-dependent) QCA metric evolution → `{norm_rel_err, superpose_err, unitary_ok, nonlinear_ok, causal_ok, front_speed, g10, lam, residual, ok}` (norm preserved to machine precision; the superposition fails) |
| `verify_analytic_vs_computed` / `verify_blackhole` / `verify_flat` / `verify_nonlinear` | `(h=H_FD)` | the 2-route checks: the computed vs analytic de Sitter Riemann; the Schwarzschild vacuum ($G=0$ by computation) + the computed Kretschmann; the flat vacuum; the ΓΓ self-gravity (absent for flat, present for de Sitter) |
| `verify` | `(tol=1e-9)` | the complete 4-D non-linear gravity verdict (all the checks above + the non-linear QCA) |

**`cl10.numrel` — the 3+1 BSSN/Z4 numerical relativity (V51, C2, W9).**

The **complete** 3+1 numerical relativity (the BSSN/Z4 completion of the W9
gap): the 3+1 Ricci **0-row** $(R_{00}, R_{0i})$ — the **constraints** $H, M$
of the 4-D Einstein tensor — is **2-route verified** (the 4-D Riemann,
`cl10.gravity`, vs the 3+1 $(\gamma, K)$ decomposition — $K_{ij} =
-\tfrac12\partial_0\gamma_{ij}$, the 3+1 Christoffels $\Gamma^{(3)}$ —
agreeing to $O(h^2)$, residual $\sim 10^{-13}$); the **Z4**
constraint-damping (a perturbed constraint decays to 0; the undamped case
does not); the **non-linearity** (the Ricci 0-row is quadratic in $K$, so
superposition fails); and the **UV-finite** property (the 3+1 quantities are
bounded by the lattice cutoff).  The **open** part of W9 is now only the
**quantum** completion.

| symbol | signature | returns |
|---|---|---|
| `H_FD` / `DAMP` | (module constants) | the 0-parameter finite-difference step ($10^{-3}$, a *numerical* convergence parameter) / the Z4 damping rate ($1$) |
| `gamma` / `K` | `(g4fn, x)` / `(g4fn, x, h=H_FD)` | the spatial 3×3 metric $\gamma_{ij} = g_{ij}$ / the extrinsic curvature $K_{ij} = -\tfrac12\partial_0\gamma_{ij}$ (the 0-parameter gauge $\alpha=1$, $\beta=0$) |
| `G3` / `R3` | `(g4fn, x, h=H_FD)` | the 3+1 Christoffels $\Gamma^{(3)a}{}_{ij}$ / the **3-D Riemann** $R^{(3)a}{}_{bcd}$ (the BSSN $\mathscr{R}$) |
| `ricci_0row_4d` / `ricci_0row_3p1` | `(g4fn, x, h=H_FD)` | **Route 1** (the 4-D Riemann, `cl10.gravity`) / **Route 2** (the 3+1 $(\gamma, K)$) 0-row $(R_{00}, R_{0i})$ — the 2-route verification agrees to $O(h^2)$ |
| `constraints` | `(g4fn, x, h=H_FD)` | the 3+1 **constraints**: the **Hamiltonian** $H = {}^{(3)}\!R + K^2 - K_{ij}K^{ij}$ and the **momentum** $M_i = D_j K^j{}_i - D_i K$ |
| `verify_2route` | `(g4fn, n_pts=24, h=H_FD, tol=1e-4, seed=11)` | the **2-route** check over `n_pts` random points → `{max_r00, max_r0i, h2, ok}` |
| `nonlinear_superpose` | `(a_fn, b_fn, x, h=H_FD, tol=1e-8)` | the **non-linearity**: the 0-row of the sum metric $(a+b)$ vs the sum of the 0-rows (fails) |
| `z4_damping` | `(g4fn, x, n_steps=20, h=H_FD, lam=DAMP)` | the **Z4** damping → `{H0, damped_final, undamped_final, damped_ratio, undamped_ratio, ok}` (the damped $H$ decays, the undamped does not) |
| `uv_finite` | `(g4fn, x, h=H_FD)` | the **UV-finite** (bounded) check → `{max_abs, bounded, ok}` |
| `metric_31` | `(x)` | the 3+1 **test** metric (a spatial + time-dependent $\gamma_{ij}$ with an off-diagonal term) |
| `verify` | `(tol=1e-4, seed=11)` | the complete 3+1 NR verdict (the 2-route 0-row + the Z4 damping + the non-linearity + the UV-finite property) |

**`cl10.unification` — the 2-loop running (V47, W4).**

| symbol | signature | returns |
|---|---|---|
| `running_2loop` | `()` | `{backend, qca_2step_unitary_residual, qca_2particle_bounded, qca_2particle_max_offdiag, pert_1loop_landau_pole, pert_2loop_landau_pole, pert_2loop_divergent, uv_finite, residual}` — the **2-loop** running: the **UV-finite** QCA (the 2-step is unitary, the 2-particle is bounded) vs the **UV-divergent** perturbative QFT (the 2-loop Landau pole is below M_P) |

**`cl10.smatrix` — the n-particle amplitudes (V48, W9).**

| symbol | signature | returns |
|---|---|---|
| `n_particle` | `()` | `{backend, n1_unitary_residual, n2_unitary_residual, n3_unitary_residual, n4_unitary_residual, n1_bounded, uv_finite, residual}` — the **3- and 4-particle** amplitudes: the n-particle S-matrix is the Kronecker product of the 1-particle, so it is unitary (0) for all n (the **UV-finite** property, the n-step is bounded for all n) |

**`cl10.yukawa_3x3` — the quark CKM, honestly (V49, W6).**

| symbol | signature | returns |
|---|---|---|
| `quark_ckm` | `()` | `{backend, up_eigenvalues, down_eigenvalues, two_level, ckm_unitary, ckm_is_identity, ckm_computed_angles, sin2_12, sin2_23, sin2_13, obs_ckm, dev_ckm, residual}` — the quark mass matrices on the democratic base (the 2-level structure, 0-parameter). **Honest:** with H_u = H_d both democratic the computed CKM **is the identity** (`ckm_is_identity: true`), not the TBM; a non-trivial CKM needs an up/down misalignment the democratic base does not supply — that is the W6 named residual |

---

## 13. `cl10.parameter_audit` — the honest parameter audit (V52)

The **solution to the four limits of the theory**.  All numerics on
the **CuPy** GPU engine; the root cross-checked by a **second route** (the
symbolic `kingdon` $\mathrm{Cl}(10,0)$).

| symbol | signature | returns |
|---|---|---|
| `kingdon_2route` | `()` | `{backend, kingdon_available, kingdon_sq_minus_1, kingdon_anticommute, kingdon_w2_scalar, pauli_sq_minus_I, pauli_anticommute, pauli_w2_plus_I, two_route_agree}` — the **2-route root**: the symbolic, exact `kingdon` $\mathrm{Cl}(10,0)$ ($e_a^2=1$, $\{e_a,e_b\}=0$, $w^2=-1$) vs the numeric $32\times32$ symplectic-Pauli on the GPU |
| `spectral_mirror` | `()` | `{backend, g10, n_levels, omega_min, omega_min_mass_GeV, omega_min_16weyl, levels_16weyl, w_pure_imaginary, commute_Uw, preserves_chiral_split, mirror_resid, mirror_ok, mirror_resid_all_g10, mirror_holds_for_all_g10, g10_grid, zp_16, zp_16bar, det_U_minus_1, sl32}` — the **16/16-bar spectral mirror = the spectral form of CPT**: $U$ commutes with the pseudoscalar (`[U,w]=0`, preserves the chiral split), and the **antimatter** spectrum is the **mirror** of the **matter** one, $\omega_{16\bar}=-\omega_{16}$, i.e. $\mathrm{spec}(U|_{16\bar})=\mathrm{spec}((U|_{16})^\dagger)=\overline{\mathrm{spec}(U|_{16})}$ (C = matter↔antimatter since $w$ is purely imaginary; T = $U\to U^\dagger$). **Holds for every 't Hooft angle** (`mirror_resid_all_g10`, `g10_grid`) — a 0-parameter property of the algebra, refining V24's CPT theorem |
| `cp_odd_theorem` | `(n4_sample=200000, seed=0)` | `{backend, max_cpodd_2rot, max_cpodd_3rot, max_cpodd_4rot_sample, max_cpodd_b_sandwich, cpodd_democratic_yukawa, n4_sample, all_cp_odd_zero, sakharov, verdict}` — the **0-parameter CP-odd theorem**: every CP-odd (T-odd) invariant vanishes (2- and 3-rotation full, 4-rotation sample, every B-violating $h_i R_{ab} h_j R_{cd}$ sandwich full, rank-1 real democratic Yukawa) $\Rightarrow$ the QCA is CP-conserving, B-conserving, in equilibrium $\Rightarrow$ **none** of the 3 Sakharov conditions $\Rightarrow$ **no** baryon asymmetry (a theorem) |
| `g10_status` | `()` | `{backend, g10_0param, n_qubits, mapping_convention, uv_finiteness_fixes_g10, pole_2loop_below_MP_for_all_g10, g10_is_free_t_hooft_angle, grid, verdict}` — $g_{10}$ is **not** fixed by the UV-finiteness (2-loop Landau-pole) condition (the pole is below $M_P$ for every $g_{10}$) $\Rightarrow$ the **one free** 't Hooft angle |
| `hierarchy_status` | `()` | `{backend, v_over_EP_observed, v_over_EP_claimed_0param, match_rel_err, match_is_2digit_coincidence, omega_min, omega_min_mass_GeV, smallest_sub_planck_scale_is, spectrum_exponent_matches_observed, hierarchy_exists_0param, residual}` — the hierarchy **exponent** $-17$ is **not** in the single-step spectrum (the smallest level $0.0115$); the $0.18\%$ match of $2.02\times10^{-17}$ is a **2-digit coincidence**; the 0-parameter content is the **existence** of the hierarchy (W8) |
| `dm_status` | `()` | `{backend, soliton_mass_GeV, abundance_gap_0param_over_obs, omega_min_mass_GeV, dm_observed_range_GeV, smallest_sub_planck_scale_above_observed_DM_orders, observed_DM_not_0param, residual}` — the **only** sub-Planck 0-parameter scale (the smallest spectral gap $\sim10^{17}$ GeV) is ~14 orders above the observed DM; the 0-parameter abundance gap $\sim10^{-52}$ (W7, quantified) |
| `provenance_table` | `()` | `{backend, rows, counts_by_provenance}` — every "0-parameter" number **reclassified** by its true provenance: **algebraic_0param** (7), **free_1param** (1: the free 't Hooft angle $g_{10}$), **convention** (2: $g_0$, $M_{\rm GUT}=M_P$), **coincidence** (1: the $2.02\times10^{-17}$ match), **named_residual** (4) — the direct answer to "0-parameter is partly conventional" |
| `verify` | `(tol=1e-9)` | `{backend, two_route, spectral_mirror, cp_odd_theorem, g10_status, hierarchy_status, dm_status, provenance_table, invariants_ok, limits_ok, ok}` — the full audit; `ok` when the 0-parameter invariants **hold** and the four limits are **quantified** |
| `main` | `()` | runs the audit and prints the verdict |

**What V59 adds (Section 13.7).**  The `g10_status` / `provenance_table`
above are the audit's **baseline**: $g_{10}$ is the **one free** parameter of
the dynamics.  `cl10.four_limits` (V59) adds the bookkeeping observation that
the **algebra** itself has no free parameter — $g_{10}=1/5$ and $g_0=1/2^3$ are
readings of the root.  They are **conventions**, not derivations: V56 shows the
data prefer $n = 7$ over $n = 5$.  The named set is then 2 two-digit
coefficients read from the data + 1 CP phase, plus the one 't Hooft angle of
the dynamics (paper Section 14).

---

## 13.5. The 10-loop / 2-loop / QCA-step refinement modules (V53–V58)

The **loop-factor representations** of four named residuals.  The loop integer
is **represented** (chosen so the 2-digit coefficient lands on the data), not
computed: no loop expansion of the QCA exists (paper §10.1, W8).  All numerics
on the **CuPy** GPU engine.

**`cl10.hier_exponent` — the hierarchy exponent (V53, W8).**

| symbol | signature | returns |
|---|---|---|
| `loop_factor` | `(n=N_GEN)` | the $n$-loop factor $(1/16\pi)^n$; the **10** $= n_{\rm gen} = 2\times5$ (the 16/16-bar **chiral split** × the **5 qubits** of the 32-D module, $2^5 = 32$) — an identification, not a loop calculation |
| `loop_ladder` | `()` | `{F10, n_gen, n_qubits, g10_as_2x5, hierarchy, dm}` — the **same** 10-loop factor $(1/16\pi)^{10} \approx 9.7\times10^{-18}$ gives **two** results: the **hierarchy** ($v/E_P = 2.08(1/16\pi)^{10}$) and the **DM mass** ($M_P(1/16\pi)^{10} = 118$ GeV) — **one factor, two results** |
| `hierarchy17` | `()` | `{backend, loop_factor, n_gen, predicted_ratio, predicted_exponent, observed_ratio, observed_exponent, exponent_abs_diff, ratio_predicted_over_observed, scan, n_best_scan, nstep, loop_ladder, one_line}` — $v/E_P = 2.02\times10^{-17}$ **is** the 10-loop factor $\times 2.08$; the exponent $-17$ is obtained by **selecting** $n = 10$ (the $0.18\%$ match fixes the coefficient, not the exponent); the **2-digit** coefficient is **fitted** |
| `verify` | `(tol=1e-9)` | `{r, ok}` — `ok` when the 2-digit ratio is in $[1/3,3]$ and the unique best loop count is $n_{\rm gen}=10$ |

**`cl10.dm_wimp` — the sub-Planck DM scale (V54, W7).**

| symbol | signature | returns |
|---|---|---|
| `F10` | `()` | the 0-parameter 10-loop factor (same as V53) |
| `relic_abundance` | `(m_GeV)` | the **normalised** relic abundance $\Omega h^2 = \Omega_{\rm obs}(100/m)^3$ — a normalisation anchored on the observed value at 100 GeV, **not** a computed freeze-out (the required $\langle\sigma v\rangle$ is not derived) |
| `dm_wimp` | `()` | `{backend, F10, n_gen, M_P_GeV, m_DM_GeV, weak_scale_GeV, m_over_weak, omega_DM, omega_obs, omega_ratio_to_observed, v42_planck_omega, v42_planck_gap, v42_quantified_gap, one_line}` — $m_{\rm DM} = 1.18\,M_P(1/16\pi)^{10} = 118$ GeV (the WIMP/weak scale, the $1.18$ **fitted**); the **V42** Planck-soliton gap (the $10^{-52}$) is **quantified**, not resolved |
| `verify` | `(tol=1e-12)` | `{r, ok}` — `ok` when $10 \le m_{\rm DM} \le 1000$ GeV, $m_{\rm DM}/100 \in [1/3,3]$, and $0.02 \le \Omega_{\rm DM} \le 0.5$ |

**`cl10.eta_cp` — the baryon asymmetry as a product of factors of three provenances (V55, W5).**

| symbol | signature | returns |
|---|---|---|
| `cp_odd_loop_factor` | `()` | the 0-parameter CP-odd **rate** factor $g_0^4/16\pi \approx 4.9\times10^{-6}$ — the **1-loop CP-odd** decay asymmetry $\varepsilon_{CP}$ ($g_0 = 1/8 = 1/2^3$, **3** = the spatial dimension) |
| `eta_cp` | `()` | `{backend, loop_factor, structure, eta_CP_0param, eta_CP_observed, eta_CP_over_observed, g0, cp_odd_loop_factor, sphaleron, delta_B_0param, delta_B_observed, delta_B_ratio, delta_B_v40, delta_B_ratio_v40, one_line}` — the baryon asymmetry as a **product of factors of three provenances** $\Delta B = (28/79)\,\varepsilon_{CP}\,D$ (the SM sphaleron $28/79$ = an **input of the SM sector**, the LH d.o.f. count; the 1-loop CP-odd $\varepsilon_{CP} = g_0^4/16\pi$, a dimensional guess; the 2-loop dilution $D = (1/16\pi)^2$); $\Delta B = 6.8\times10^{-10}$ (**ratio 1.12**, improved from V40's 2.8); `structure` = `{sphaleron, eps_CP, D_dilution, g0, g0_is_1_over_2_cubed, loop_count}`.  **Not** a computed vertex, and the out-of-equilibrium condition fails at the model's couplings (§11.3) |
| `verify` | `(tol=1e-12)` | `{r, ok}` — `ok` when the 2-loop factor is within a factor 3 of the V40 value and the $\Delta B$ ratio is in $[1/3,3]$ (and improved on V40) |

**`cl10.g10_selection` — the 1-loop $g_{10}$ data check (V56, W4).**

| symbol | signature | returns |
|---|---|---|
| `gut_1loop` | `()` | `{a_at_M_P, mean_GUT, spread}` — the 1-loop running **up** to $M_P$: the mean GUT coupling ($\sim0.023$) and the (non-unification) spread (the W4 named residual) |
| `g10_selection` | `()` | `{backend, n_qubits, g10, alpha_10_0param, g10_matches_sm_qca, gut_1loop, ratio_to_1loop, best_n_1loop, n_closest_to_1loop, candidate_scan_gap, gap_root_n5, gap_best_n, fit_prefers_root, data_check_selects_n, one_line}` — the 0-parameter $g_{10}=1/5$ (a number of the root) gives $\alpha_{10}=1/25=0.04$, **within a factor 1.73** of the 1-loop GUT coupling; **honest scan**: the data-optimal candidate is $n=7$ (gap 0.12), **not** the root $n=5$ (gap 0.73) — `data_check_selects_n = False`, so the reading $1/n_{\rm qubits}$ is a **mapping convention**, not a fit |
| `verify` | `(tol=1e-12)` | `{r, ok}` — `ok` when the 0-parameter $g_{10}$ matches the SM QCA and the ratio to the 1-loop GUT coupling is in $[1/3,3]$ |

**`cl10.neutrino` — the 2-step seesaw on the 10-loop ladder (V57, W2).**

| symbol | signature | returns |
|---|---|---|
| `loop_ladder` | `(n)` | the 0-parameter $n$-loop scale $M_P(1/16\pi)^n$ (the "one object") |
| `dirac_mass` | `()` | the 0-parameter Dirac mass $m_D = M_P(1/16\pi)^{10} \approx 118$ GeV (the 10-loop / WIMP scale) |
| `twostep_scale` | `()` | the 0-parameter 2-step (Majorana) scale $M_N = M_P(1/16\pi)^3 \approx 9.6\times10^{13}$ GeV (the $10^{14}$ scale, **replacing** the hard-coded $10^{14}$) |
| `seesaw_2step` | `()` | `{backend, m_D_GeV, M_N_GeV, m_nu_2step_eV, m_nu_1step_eV, m_nu_ladder_eV, m_nu_obs_eV, off_2step, right_order, ladder_match, exp_2step, residual}` — the 0-parameter 2-step seesaw $m_\nu = m_D^2/M_N = M_P(1/16\pi)^{2\times10-3} = M_P(1/16\pi)^{17} \approx 0.15$ eV (**2.9x** the observed; `exp_2step = 17` the 0-parameter loop identity) |
| `twostep_unitary` | `(tol=1e-9)` | `{backend, n, unitary_err_1step, unitary_err_2step, two_step_unitary, uv_finite}` — the 2-step $U^2$ of the GUT QCA is a 0-parameter **unitary** (the 2-loop), computed on the GPU |
| `verify` | `(tol=1e-12)` | `{ok, checks, seesaw, twostep}` — `ok` when the 10-loop Dirac = the DM/WIMP scale, the 2-step is unitary, the loop identity holds, and the result is the right order |

**`cl10.uvrunning` — the UV-finite discrete RG (V58, W4).**

| symbol | signature | returns |
|---|---|---|
| `one_loop_running` | `(M)` | the 1-loop SM couplings at scale $M$ (the non-SUSY SM) |
| `landau_poles` | `(M)` | the 1-loop Landau-pole scale of each coupling ($U(1)$ a real pole; $SU(2)/SU(3)$ asymptotically free) |
| `two_loop_pole_u1` | `()` | the 2-loop $U(1)$ Landau-pole scale ($\sim10^{12}$ GeV, **below** the cutoff — the perturbative UV divergence) |
| `discrete_rg` | `(max_power=8)` | `{backend, n, rows, all_bounded, two_step_unitary_err, two_particle_max_offdiag, two_particle_bounded}` — the **coarse-graining tower** $U^{2^p}$ ($p=1\ldots8$): each step is **unitary** (bounded) (the GPU) $\Rightarrow$ the **discrete RG is UV-finite** (no Landau pole); the 2-particle S-matrix is bounded |
| `w4_verdict` | `()` | `{backend, alpha_10, M_P_GeV, couplings_at_M_P, spread_at_M_P, one_loop_landau_poles, two_loop_u1_pole, u1_1loop_above_cutoff, u1_2loop_below_cutoff, discrete_rg, w4_resolved, residual, ok}` — the full W4: the 1-loop non-unification (the $U(1)$ pole) **resolved** by the UV-finite discrete RG |
| `verify` | `(tol=1e-9)` | `{ok, checks, verdict}` — `ok` when the discrete RG tower is all bounded, the 2-step is unitary, and the 2-loop $U(1)$ pole is below the cutoff (the W4 resolution) |


## 13.7. `cl10.four_limits` — the four limits and their common object (V59)

The four "limits" of the theory are treated as the **four projections of one
invariant**, $10 = 2\times5$ (the chiral $16/16$-bar split $\times$ the 5 qubits
of the 32-D module), **plus one irreducible CP phase**.  Computed on the GPU.
This is **bookkeeping**: it groups the residuals, it does not derive any of
their values (paper Section 14).

| symbol | signature | returns |
|---|---|---|
| `invariance_g10` | `()` | **L1** — `{backend, grid, max_mirror_residual_over_grid, max_abs_detU_dev_over_grid, g10_invariant, n_algebraic_invariants, free_angles, verdict}` — the 7 algebraic invariants are **$g_{10}$-invariant**: the 16/16-bar CPT-mirror residual stays $\sim1.8\times10^{-15}$ for $g_{10}=0.1..0.5$ (a property of the $\mathrm{Cl}(10)$ **algebra**, not of the free angle) |
| `unique_10loop` | `(max_n=20)` | **L2+L4** — `{..., n_best, 5_is_qubit_count, unique_n_within_1decade_of_both, n_candidates_within_1decade, factor2_test_of_both_empty, coeff_hierarchy, coeff_dm, scan, verdict}` — the **same** $n=10$ is the **unique** $n\in1..20$ within **one decade** of **both** closures: the hierarchy ($2.08(1/16\pi)^{10}$) and the WIMP mass ($M_P(1/16\pi)^{10}=118$ GeV) — **two statements of one scale** (the strict factor-2 test of both is **empty**); the 2-digit coefficients are **named residuals** (computed from the data, not derived) |
| `cp_phase_reduction` | `()` | **L3** — `{..., rank1_cp_phases, two_N_2flavor_cp_phases, cpodd_at_delta_0, cpodd_at_delta_pi4, cpodd_at_delta_pi2, minimal_cp_phases_required, delta_B_0param, delta_B_ratio, verdict}` — rank-1 (1 $N$) has **0** CP phases (no baryon); the 2-step (2 $N$) has **1** ($\mathrm{Im}\det(Y)=|\sin\delta|$) → baryogenesis needs **one** phase, $\Delta B = 6.8\times10^{-10}$ (ratio 1.12) |
| `four_limits` | `()` | the **assembly** — `{backend, unifying_invariant, L1_conventionality, L2_hierarchy_plus_L4_DM, L3_baryogenesis, coupling_provenance, residual_set, total_named_residuals, mapping_conventions, all_resolved, verdict}` — the **algebra** has no free parameter, while $g_{10}=1/5$ and $g_0=1/2^3$ are **conventions**; the named residual set is **2 two-digit coefficients** (2.08, 1.18, read from the data) + **1 CP phase** + the mapping conventions + the one 't Hooft angle of the dynamics; the $\delta\leftrightarrow g_{10}$ identification is **false** (the QCA $U(g_{10})$ is real / CP-even) |
| `verify` | `(tol=1e-12)` | `{r, ok}` — `ok` when L1 is $g_{10}$-invariant, the unique $n=10$ is within one decade of both closures, and the CP-phase count is 0 (rank-1) / 1 (2-step) |
| `main` | `()` | runs the module, prints the L1 / L2+L4 / L3 results + the reduced residual set + the verdict, returns `0` if `ok` |

---

## 13.8. `cl10.cosmology` and `cl10.electroweak` (V34, V32)

**`cl10.cosmology` — the cosmological constant and the hierarchy (V34).**

| symbol | signature | returns |
|---|---|---|
| `cosmological_constant` | `()` | `{backend, n_even=512, n_odd=512, cancellation, lambda_0param=0.0, lambda_0param_natural=1.0, lambda_obs=1e-122, w=-1.0, residual=1e-122}` — the $\mathbb{Z}_2$ grading is **512 even / 512 odd**, each mode carrying a $0$-point $\pm\tfrac12$, so the tree-level vacuum energy cancels **exactly**; the residual is the observed $10^{-122}M_P^4$ (W1) |
| `hierarchy` | `()` | `{backend, E_P=1.22e19, v=246.0, v_over_EP_0param=2.02e-17, observed=2.016e-17, rel_err=1.8e-3, vast_gap, residual}` — the $v/E_P$ match; the residual is the **derivation of the exponent** (W8) |
| `verify` | `(tol=1e-9)` | the V34 verdict |

**`cl10.electroweak` — the QCA$\to$QFT correspondence and the electroweak sector (V32).**

| symbol | signature | returns |
|---|---|---|
| `qca_step_spectrum` | `(theta=0.2)` | `{backend, theta, n_levels=32, spectrum_spread, unitary_err, effective_action}` — the QCA step as the QFT time-evolution (the 't Hooft dictionary) |
| `electroweak_breaking` | `()` | `{backend, g10, g2, gprime, vev_v=246.0, m_W=87.2, m_Z=110.3, sin2_w=0.375, mw_over_mz=0.791, electroweak_breaking_ok}` — the breaking **at the GUT angle**: $\sin^2\theta_W = 3/8$ is the unification-scale boundary condition, and $m_W, m_Z$ here are the tree-level values at that scale, **not** the measured $80.4/91.2$ GeV.  The low-energy values require running, which the theory does not currently perform consistently (§11.1) |
| `verify` | `(tol=1e-9)` | the V32 verdict |

---

## 13.9. `cl10.brillouin` — the Brillouin-zone chirality test (V63, W9)

The lattice-chirality check the suite had never run: does the **implemented**
local step carry a single chiral 16 sector, or does the Nielsen–Ninomiya
doubling survive?  It does not evade the theorem by hypothesis (no lattice
Dirac operator), so the question is answered on the spectrum.  Protocol and
cross-validation: `doc/reference.md` §5.

| symbol | signature | returns |
|---|---|---|
| `blocks` | `(g10=G10)` | `(B, U32, U16)` — the pseudoscalar $+i$ basis $B$ (the 16), the internal 't Hooft step on the 32, and $U_{16} = B^\dagger U_{32}B$ |
| `floquet` | `(q, U, order=ORDER, g10=G10)` | the exact $128\times128$ Bloch Floquet $U_{\mathrm{F}}(q) = \prod_{t=0..5}\big(S_{j(t),p(t)}(q)\otimes U\big)$ over the 6-step cycle, on the $2\times2\times2$ supercell (8 sites) at supercell momentum $q \in [-\pi,\pi]^3$ |
| `sector_closure` | `(n_steps=6, L=8)` | the $16\to\overline{16}$ leakage of the **implemented** evolution (internal layer + staggered bond layer, `cl10.local_qca`): $9.7\times10^{-16}$ — chirality is a dynamical symmetry |
| `torus_validation` | `()` | the maximum phase-spectrum deviation between the Bloch Floquet and the real-space 6-step Floquet on a $4^3$ torus, over the 8 supercell momenta $\{0,\pi\}^3$: $4.4\times10^{-15}$ |
| `brillouin` | `()` | `{backend, g10, floquet_cycle=6, sector_leakage_16_to_16bar, bloch_vs_torus_dev, gap_q0_16=0.036499, gap_q0_32=0.036499, n_doubler_nodes=12, max_doubler_gap=2.0e-4, nodes[{q, gap, slope, chirality}], chirality_sum=0.0, v_typical=0.1419, g10_scan[{g10, gap_q0, gap_doubler}]}` — the zone centre is **gapped**, the 12 momenta $(\pm\pi/2,\pm\pi/2,0)$ are **gapless** with paired chiralities |
| `verify` | `(tol=1e-9)` | the verdict dict with `ok`: leakage $<\varepsilon$, Bloch-vs-torus $<10^{-9}$, $\mathrm{gap}(0) > 10^{-3}$, $\max\mathrm{gap}(\text{doubler}) < 10^{-3}$, $\sum\text{chirality} = 0$, 12 nodes.  `ok` is a statement about the **computation**; its physics content is a **negative** result for the theory (the doublers are present, W9) |

---


## 14. `run_validation.py` — the verification runner

| symbol | signature | returns |
|---|---|---|
| `find_experiments` | `()` | a list of the experiment directories under `experiments/` |
| `run_one` | `(exp_dir)` | `{id, dir, tier, property, residual, pass, status, returncode, wall_time, verdict}` — runs one experiment and reads its JSON verdict |
| `main` | `()` | runs all experiments (or a subset from `argv`), writes `results/validation_summary.json`, prints the hard-tier verdict, returns `0` if all pass |

**The summary** (`results/validation_summary.json`) has:
`{n, n_pass, pass, u10, u10_criterion, named_residuals, reduced_walls, experiments}`
where `u10` is the runner's **hard-tier verdict** (the hard tiers 0–3 all pass;
the paper does not use the `U10` label) and
`named_residuals` lists **all** the named/quantified residuals (the exact
values the theory does not derive — W1–W9: the cosmological constant, the seesaw, the families, the running, the Yukawa, the baryon asymmetry, the dark-matter scale, the hierarchy, the quantum sector).

**Usage:**

```bash
python3 run_validation.py              # all 117 experiments (V0-V116), in parallel
python3 run_validation.py V11 V12      # a subset (--jobs N caps the pool)
```

---

## 15. The simulation packages (V61, V62)

**`simulation.young` — the double slit as a single-electron measurement (Fig. 18).**

| symbol | signature | role |
|---|---|---|
| `wave.WaveSim` | `(config)` | the 3-D FDTD leapfrog on the 6-neighbour discrete Laplacian (`r = 0.5` CFL), the source matched to the FDTD dispersion, the quadratic absorbing sponge, the double-slit barrier |
| `collapse.Collapse` | `()` | `sample` (Born-rule inverse-CDF draw of one detection point from the time-integrated on-screen $\psi^2$), `contract` (the non-unitary wave-packet reduction), `add` (the accumulated detections) |
| `runner` / `__main__` | — | the headless run (PyVista off-screen render) and the live VisPy window |

```bash
python -m simulation.young                 # headless -> simulation/young/video/
python -m simulation.young --interactive   # live VisPy window
```

**`simulation.hydrogen` — the hydrogen atom driven by the lattice photon (V61, V62; Fig. 19).**

| symbol | signature | role |
|---|---|---|
| `qca` | — | the free CL10 massless 3-D field (the 6-neighbour leapfrog), the discrete Laplacian and the field observables |
| the atom model | — | the real Rydberg levels $E_n = -13.606/n^2$ and eigenfunctions $R_{nl}Y_{lm}$ on a $296^3$ grid, the Fermi-golden-rule dipole drive with a Lorentzian line shape, the Einstein-$A$ re-emission |
| observables | — | `populations`, `atom_energy`, `excited_fraction`, `electron_density` $=|\psi|^2$, `photon_density` $=|F|^2+|F_{\rm out}|^2$ |

```bash
python -m simulation.hydrogen     # -> simulation/hydrogen/video/ (stills + mp4; git-ignored)
python3 experiments/V61_hydrogen/V61_hydrogen.py        # the fast 120^3 check (part of the suite)
python3 experiments/V62_hydrogen_production/V62_hydrogen_production.py   # the production run -> Fig. 19
```

The atom's knobs (`A`, `h`, `gamma`) are physical parameters of the atom, not
fits to data.  The atomic spectrum is an **input of atomic physics**: the
experiment tests that the CL10 photon couples to it correctly (the gaps, the
E1 selection rules, the resonant absorption and the re-emission), not that the
theory derives hydrogen.
