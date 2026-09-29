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
cl10.toepred       # the TOE predictions (Section 7)
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

**`anomalies`** is the **hard TOE requirement**: a chiral gauge theory is
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

## 7. `cl10.toepred` — the TOE predictions

| symbol | signature | returns |
|---|---|---|
| `proton_decay_lifetime` | `()` | `float` — the **0-parameter** proton lifetime $\tau_p \sim M_{\rm GUT}^4/(\alpha_{10}^2 m_p^5)$ in **years** ($\approx 4\times10^{47}$; the proton appears stable) |
| `seesaw_mass` | `()` | `float` — the **0-parameter** light-neutrino mass $m_\nu \sim m_D^2/M_{\rm GUT}$ in **eV** ($\approx 1.2\times10^{-6}$; light, the exact value a wall) |
| `family_candidate` | `()` | `{candidate, from_spatial, from_internal_half, derived}` — the 0-parameter **3-family** candidate (`derived=False`: a named wall) |
| `main` | `()` | runs the V14–V16 self-check |

The 0-parameter constants: `M_PLANCK = 1.22e19` GeV, `M_P_ROBSON = 0.938`
GeV, `ALPHA_10 = 1/25`, `V_EW = 246.0` GeV.

---

## 8. `cl10.walls` — the three reduced walls (W1/W2/W3) + the baryon asymmetry

The **0-parameter mechanisms** that reduce the three walls (see
`doc/publication.md`, Section 13, Table 2):

| symbol | signature | returns |
|---|---|---|
| `w1_cosmological_constant` | `()` | `{n_even, n_odd, cancellation, lambda_0param, lambda_obs, residual}` — the **W1**: the $\mathbb{Z}_2$ $512/512$ **cancellation** $\Rightarrow$ $\Lambda = 0$ at 0-parameter (the **s-mechanism**); the exact $10^{-122}$ the named residual |
| `w2_seesaw` | `()` | `{g0, m_D_GeV, M_N_GeV, m_nu_2step_eV, m_nu_1step_eV, m_nu_obs_eV, off_1step, off_2step, right_order, residual}` — the **W2**: the **2-step** seesaw ($M_N$ intermediate, $g_0 v$ Dirac) $\Rightarrow$ $m_\nu \sim 10^{-2}$ eV (the right order) |
| `w3_families` | `()` | `{n_time, n_space, n_fam_candidate, n_observed, match, derived, residual}` — the **W3**: the **3** families $=$ the **3** spatial directions of the $(1,3)$ (one 0-parameter number, dual use) |
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

The four **deepening** modules extend the 0-parameter structure to new
hallmarks of a true TOE (all run on the **GPU** engine).

**`cl10.holo` — the holographic area law (V25).**  The QCA vacuum entropy of
a region of side $L$ scales as $S \propto L^2$ (the **surface area** $6L^2$),
**not** $L^3$ (the volume); the 0-parameter entropy per unit area is the
Bekenstein-Hawking $S = A/4$ ($G=1$); the holographic ratio $S/V \propto 1/L$.

| symbol | signature | returns |
|---|---|---|
| `area_law` | `(L_values, per_plaquette=0.25)` | `{L_values, per_plaquette, S, R, area_law_exponent, holographic_ratio_exponent, backend}` — the 0-parameter area law (V25) |

**`cl10.family` — the 3-family construct (V26).**  The **3** observed families
are a 0-parameter consequence of the one object in **three** ways (the 3 space
directions; the 3 two-planes $\binom{3}{2}$; the rank-3 spatial metric); the 16
Weyl is one generation $\to$ 3 generations $= 16 \times 3 = 48$ Weyl.

| symbol | signature | returns |
|---|---|---|
| `family_construct` | `()` | `{n_space, n_planes_C32, n_rank, n_family_candidate, n_families_observed, match, all_three_ways_3, n_weyl_1gen, n_weyl_3gen, residual}` — the 0-parameter 3, three ways (V26) |

**`cl10.qcd` — QCD confinement (V27).**  The 0-parameter colour count
$N_C = 3$ (the 8 gluon $= N_C^2 - 1$); the **area law**
$\langle W(C)\rangle \propto \exp(-\sigma \cdot \mathrm{Area})$; the 0-parameter
**string tension** $\sigma$; the **linear** potential $V(r) = \sigma r$ (in
contrast to the Coulombic $1/r$).

| symbol | signature | returns |
|---|---|---|
| `color_count` | `()` | `(n_c, n_gluon, ok)` — the 0-parameter $N_C=3$ |
| `string_tension` | `()` | the 0-parameter string tension $\sigma$ |
| `confining_potential` | `(sigma, r_values)` | `[{r, v_qcd, v_qcd_over_r, v_em}, …]` — the linear QCD potential vs the Coulombic EM (V27) |

**`cl10.yukawa` — the Yukawa / mass-hierarchy (V28).**  The **10-Higgs
16→16-bar** (the grade-1 GUT vector is **off-diagonal** in the 16/16-bar split)
gives the 0-parameter Dirac mass term; the **3×3** mass matrix (the 3
generations) is **non-singular** (rank 3 — all generations get mass).

| symbol | signature | returns |
|---|---|---|
| `yukawa_matrix` | `()` | `{n_gen, n_dim, rank, non_singular, off_diagonal_Dirac, diag_zero, mass_term, residual}` — the 0-parameter 3×3 non-singular Dirac mass matrix (V28) |


## 12. The completion program modules (V36–V45)

The **completion program** (C1–C7, `doc/publication.md` Section 12) closes the gap
from a consistent 0-parameter framework to a **complete** TOE.  Each module
function is **0-parameter** (no fit) and **backend-agnostic** (the GPU
engine), and the matching experiment writes a JSON verdict.

**`cl10.local_qca` — the local QCA (V36, C1).**  The **local** QCA step on
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
soliton, the Schwarzschild black hole), and — **V37 (C2)** — the
**derivation** in the direction **dynamics ← action**:

| symbol | signature | returns |
|---|---|---|
| `einstein_hilbert` | `(r=2.0, tol=1e-4)` | `{R_sphere, R_flat, G, coupling=1/16π, n_deriv, einstein_hilbert_ok, …}` — the 0-parameter EH action on prescribed metrics (V31) |
| `eh_quadratic_action` | `(u)` | the 0-parameter **discrete** quadratic EH action on the Planck lattice ($G=1$, $c=1$, $dt^2=1/4$) |
| `einstein_hilbert_from_dynamics` | `(L=8, T=24, n_probe=16, seed=0, tol=1e-6)` | `{elm_on_qca_trajectory, stationary_ok, elm_off_dynamics, selects_ok, coeff_ok, continuum_ratio, continuum_ok, ok, …}` — the **Euler–Lagrange equation** of the 0-parameter action **is** the QCA metric dynamics (V37): the numerical functional derivative vanishes on the QCA trajectory, an off-dynamics field has a non-zero ELM, the coefficient is $1/16\pi$ ($G=1$), and the discrete dispersion gives the massless $c=1$ continuum |

**`cl10.unification` — the 0-parameter numbers derived from the root (V38, C3).**

| symbol | signature | returns |
|---|---|---|
| `g10_from_root` | `()` | `{n_module=32, n_qubits=5, g10=1/5, alpha_10=1/25, n_space=3, g0=1/8, …}` — the GUT coupling $g_{10}=1/n_{\rm qubits}$ (the 5 qubits of the 32-D module) and the Dirac coupling $g_0 = 1/2^{n_{\rm space}}$ (the 3-space), derived from the root |
| `running_down` | `()` | `{alpha_10_at_M_P, couplings_at_M_Z, measured_at_M_Z, ratio_to_measured, u1_down_pole_scale, residual, …}` — the 1-loop running **down** from the 0-parameter UV input at the QCA cutoff to $M_Z$ (the 1-loop matching is the **quantified** named residual W4) |

**`cl10.family` — the 3 families inside the dynamics (V39, C4).**

| symbol | signature | returns |
|---|---|---|
| `family_qca` | `()` | `{n_module_3gen=96, n_weyl_3gen=48, n_cpt_3gen=48, unitary_err, family_symmetry_off_block, …, tribimaximal_sin2, residual, …}` — the **3-generation QCA step**: the 45 act **identically on each family** (the family symmetry; the 96×96 is exactly unitary, the per-family 16/16-bar split is preserved) + the 0-parameter tribimaximal mixing (V29) |

**`cl10.yukawa_3x3` — the 3-gen Yukawa mass matrix: the 9-to-1 projection (V45, W6, C4).**

| symbol | signature | returns |
|---|---|---|
| `w6_verdict` | `()` | `{g10, one_object (32-level spectrum + the 96-level 3-family), family_symmetry (the 3-fold symmetry: the 3 spatial directions S_a are equivalent -> the democratic base is DERIVED, not asserted; S_a^2=-I, {S_a,S_b}=0, equal norm), mass_matrix (the democratic 3x3: eigenvalues {3,0,0}, 1 heavy + 2-fold light, rank 1), mixing (the 0-parameter TBM = the 0-parameter PMNS, sin^2 = 1/3, 1/2, 0), hierarchy (the 0-parameter 2-level structure; the 3-level is the named residual), higgs_dirac (the computed 16->16-bar Dirac term, rank-16, off-diagonal), …}` — the **W6** projection of the one object: the **3-gen Yukawa mass matrix** (the **derived** democratic base) + the **0-parameter TBM mixing** + the **2-level structure** + the **0-parameter rank-16 10-Higgs Dirac term**; the **exact** 3-gen hierarchy (the 10-Higgs vev) and the **exact** CKM/PMNS are the named residual |

**`cl10.baryogenesis` — GUT leptogenesis (V40, C5).**

| symbol | signature | returns |
|---|---|---|
| `leptogenesis_estimate` | `()` | `{M_N_GeV=1e14, gamma_over_h, g0, eps_loop=g0^4/16π, eta_cp (named), sphaleron_factor=28/79, delta_b, observed_delta_b, ratio_to_observed, …}` — the baryon asymmetry at the **right order** (the mechanism 0-parameter, the exact value the named residual) |

**`cl10.smatrix` — the S-matrix of the QCA (V41, C6).**

| symbol | signature | returns |
|---|---|---|
| `qca_smatrix` | `()` | `{g10, n_planes=45, s1 (32×32), s2 (1024×1024), unitary_residual, probability_conservation, max_offdiag, …}` — the **'t Hooft program** completed: the QCA step **is** the 0-parameter S-matrix (the 1- and 2-particle transition amplitudes; the only input $g_{10}=0.2$) |

**`cl10.lattice_rg` — the QCA lattice RG: the S-matrix to the cross-sections (V44, W9, C6).**

| symbol | signature | returns |
|---|---|---|
| `lattice_rg` | `()` | `{g10, one_particle (P1), two_particle (P2, the 1024x1024 S^(2) = U (x) U), uv_bound (sigma_max = 4 pi / M_P^2), cross_section (sigma_0 = P2 / M_P^2, sigma_pert, R, elastic, optical_theorem), dispersion (the 32-level spectrum H = -i log U, UV-bounded, zero_point), froissart (the 0-parameter Froissart-Martin bound, m_min = M_P: the QCA constant cross-section is below the ln^2 growth -> strongest UV-finiteness), uv_finite, …}` — the **'t Hooft program** at the cross-section level: the QCA is **UV-finite** (the Planck cutoff bounds every cross-section, `sigma_0 < sigma_max`), the **0-parameter** 2->2 cross-section `sigma_0 = P2 / M_P^2`, the **0-parameter dispersion** (the 32-level spectrum, UV-bounded `|omega| <= pi`, no 1-loop pairing), the **0-parameter optical theorem** (the inelastic `P2` **is** the cross-section), the **0-parameter Froissart-Martin bound** (the QCA's *constant* cross-section is **below** the `ln^2` growth — the **strongest** UV-finiteness), and the **0-parameter** QCA/perturbative ratio `R = 4 pi P2 / g_10^4`; the **named residual** (W9) is the lattice RG matching to the measured cross-sections |

**`cl10.cosmology3` — the unified cosmology + the DM scale (V42, C7).**

| symbol | signature | returns |
|---|---|---|
| `dm_mass_scale` | `()` | `{E (the soliton energy, Planck units), m_DM_GeV, sigma, Omega_0param, Omega_obs, abundance_gap, residual, …}` — the **quantified** dark-matter residual: the Planck-mass soliton (stable, massive, neutral) and the computed 0-parameter abundance gap (the sub-Planck scale is the named residual) |

**`cl10.qca_spectrum` — the one 0-parameter object + the 9-to-1 reduction (V43).**

| symbol | signature | returns |
|---|---|---|
| `levels` | `(g10=None)` | the $32$ sorted phases of $U$ (the levels of $H = -i\log U$), in $(-\pi,\pi]$ |
| `qca_spectrum` | `(g10=None)` | `{g10, n_levels=32, spectrum, zero_point (=0: the ±-paired 0-point), det_U_minus_1 (=0: the SL(32) identity), sl32, pairing_resid, max_trace_T45 (=0: the 45 are traceless), cpt_0point_16, cpt_0point_16bar, commutes_with_pseudoscalar, max_level, …}` — the **one** 0-parameter object (the 32-level QCA spectrum) and its 0-parameter invariants (the 0-point $=0$, $\det U=1$ for every $g_{10}$, the CPT split) |
| `verify` | `(tol=1e-10)` | `{spectrum, ok}` — the 0-parameter checks: the 0-point $=0$, the SL(32) identity, the $\pm$-pairing, the CPT split |

---

## 13. `run_validation.py` — the TOE runner

| symbol | signature | returns |
|---|---|---|
| `find_experiments` | `()` | a list of the experiment directories under `experiments/` |
| `run_one` | `(exp_dir)` | `{id, dir, tier, property, residual, pass, status, returncode, wall_time, verdict}` — runs one experiment and reads its JSON verdict |
| `main` | `()` | runs all experiments (or a subset from `argv`), writes `results/validation_summary.json`, prints the **TOE verdict**, returns `0` if all pass |

**The summary** (`results/validation_summary.json`) has:
`{n, n_pass, pass, toe, toe_criterion, named_residuals, reduced_walls, experiments}`
where `toe` is the **TOE verdict** (the hard tiers 0–3 all pass) and
`named_residuals` lists **all** the named/quantified residuals (the exact
values the 0-parameter theory does not derive — W1, W2, W3, the running, the
Yukawa, the baryon asymmetry, the dark-matter scale, …).

**Usage:**

```bash
python3 run_validation.py              # all 46 experiments
python3 run_validation.py V11 V12      # a subset
```
