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
`doc/walls.md`):

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

## 9. `run_validation.py` — the TOE runner

| symbol | signature | returns |
|---|---|---|
| `find_experiments` | `()` | a list of the experiment directories under `experiments/` |
| `run_one` | `(exp_dir)` | `{id, dir, tier, property, residual, pass, status, returncode, wall_time, verdict}` — runs one experiment and reads its JSON verdict |
| `main` | `()` | runs all experiments (or a subset from `argv`), writes `results/validation_summary.json`, prints the **TOE verdict**, returns `0` if all pass |

**The summary** (`results/validation_summary.json`) has:
`{n, n_pass, pass, toe, toe_criterion, named_residuals, experiments}` where
`toe` is the **TOE verdict** (the hard tiers 0–3 all pass) and
`named_residuals` lists the three walls (W1, W2, W3).

**Usage:**

```bash
python3 run_validation.py              # all 23 experiments
python3 run_validation.py V11 V12      # a subset
```
