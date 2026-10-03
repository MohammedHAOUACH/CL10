# CL10 — Technical Reference

A **professional** reference for the **implementation** of Causal Clifford
Theory on the GUT root (CL10).  It fixes the **exact** constructions used
by the code, the **numerical** protocol, and the **reproducibility**
information.

The **physics** (the theorems, the claims, the verification program, the
nine residuals W1–W9, the falsifiers F1–F7) lives in
`doc/publication.md`; the **Python interface** lives in `doc/api.md`.
This reference is the companion to both: it says *how* each construction
is built, so that every number in the publication is reproducible.

**Conventions.**  Planck units $\hbar = c = \ell_P = t_P = 1$.  All
Clifford generators are **Hermitian** and **mutually anticommuting** with
$h_a^2 = I$.  The array backend is selected once at import (Section 3);
the physics is backend-agnostic.  All floating-point residuals are
reported against a global tolerance $\varepsilon = 10^{-12}$ (machine-
precision checks use $10^{-14}$).

---

## 1. The root: the exact $\mathrm{Cl}(10)$ construction

The **10** generators are built from **5 qubits** (a $2^5 = 32$-D complex
module) with the **symplectic Pauli pattern**: for $k = 0,1,\dots,4$,

$$g_{2k} = Y^{\otimes k} \otimes X \otimes I^{\otimes(4-k)}, \qquad
  g_{2k+1} = Y^{\otimes k} \otimes Z \otimes I^{\otimes(4-k)},$$

where $X, Y, Z$ are the Pauli matrices and $I$ is the $2\times2$ identity
(implementation: `cl10.clifford.generators`; the full basis:
`cl10.clifford.basis`).

**Properties (all verified to $\varepsilon$):**

| property | value | check |
|---|---|---|
| number of generators | $10$ | $\lvert\{g_a\}\rvert = 10$ |
| $g_a^2$ | $= I$ | $\max_a \lVert g_a^2 - I\rVert < \varepsilon$ |
| $\{g_a, g_b\}$ | $= 0$ for $a\ne b$ | $\max_{a<b} \lVert g_a g_b + g_b g_a\rVert < \varepsilon$ |
| algebra dimension | $2^{10} = 1024$ | the number of independent ordered products |
| module dimension | $32$ (complex) | the matrix size of each $g_a$ |

**Theorem (correctness).**  The pattern above yields exactly $10$ mutually
anticommuting Hermitian generators: within a block $k$, $X$ and $Z$
anticommute (on qubit $k$); between blocks, the $Y^{\otimes k}$ prefix makes
a generator at block $k$ anticommute with one at block $j<k$ (the $Y$
anticommutes with the $X$ or $Z$ of the later block).  This is the standard
construction of a $2n$-generator Clifford algebra on $n$ qubits.

The **basis** of the $1024$-D algebra is the set of all ordered products
$\{g_{a_1}\cdots g_{a_r}\}$ for $0 \le r \le 10$ (the identity, the 10
singles, the 45 doubles, the 120 triples, $\dots$, the volume element);
the **$\mathbb{Z}_2$ grading** (the 512 even / 512 odd sub-algebras) is
read off from the **parity** of $r$.

---

## 2. The construction map: every object of the publication, exactly

Each object of `doc/publication.md` is computed by **one** module (the
signatures are in `doc/api.md`):

| object (publication.md) | module | the exact construction in the code |
|---|---|---|
| the root $\mathrm{Cl}(10)$ (§1) | `cl10.clifford` | the 10 generators of Section 1 (the 5-qubit symplectic Pauli pattern) |
| the $(1,3)$ inside the GUT (§2) | `cl10.spacetime` | $T = h_1$ (time, $+$), $S_a = i\,h_{a+1}$ for $a = 1,2,3$ (space, $-$) — **four of the ten** generators, so the Minkowski relations are **identities** of the root |
| the 16 Weyl / 16-bar (§3) | `cl10.matter`, `cl10.cpt` | the **pseudoscalar** $w = h_1\cdots h_{10}$ (anti-Hermitian, $w^2 = -I$) diagonalized with a general `eig` (not `eigh` — $w$ is not Hermitian): the **$+i$** eigenspace = the 16, the **$-i$** = the 16-bar |
| the 45 and the $5/3$ (§4) | `cl10.gauge` | the 45 $= \{h_a h_b : a < b\}$; the hypercharge in the $\mathbf{5}$ is $\mathcal{T} = \mathrm{diag}(-\tfrac16,-\tfrac16,-\tfrac16,+\tfrac14,+\tfrac14)$ with $\mathrm{Tr}(\mathcal{T}^2) = 5/24$; GUT-normalized to $\tfrac12$, with $Y = 2\mathcal{T}$ $\Rightarrow$ $g_1^2/g'^2 = 8 \cdot 5/24 = 5/3$ |
| the 8 anomalies (§4) | `cl10.matter.anomalies` | the sums over the 0-parameter charges of the 6 fields of the 16 (Table 1 of the publication), with the Dynkin indices $T(3) = T(2) = \tfrac12$ |
| the metric + the graviton (§6) | `cl10.qca` | the **6-neighbour** discrete Laplacian on an $L^3$ torus; the linearised perturbation is a free wave (the discrete CFL $dt = h/2$ gives $c = 1$, the momentum is UV-bounded) |
| the non-linear metric, the soliton, the black hole (§6) | `cl10.einstein` | the **cubic** $\phi^3$ (Choptuik) completion: the static $u'' = u - u^3$ has the **exact** finite-energy solution $u = \sqrt{2}\,\mathrm{sech}\,x$; the Ricci scalar is **computed** on prescribed metrics (flat: $0$; the 2-sphere: $2/r^2$; the Schwarzschild: the Kretschmann $48M^2/r^6$) |
| the **complete** 4-D non-linear gravity (§6, C2, W9; V50) | `cl10.gravity` | the **general** 4-D Riemann / Ricci / Einstein / Kretschmann **computed by central finite differences from an *arbitrary* metric** $g_{\mu\nu}(x)$ (no hard-coded Christoffels, not restricted to one solution); **2-route verified**: the computed de Sitter $R^{\mu}_{\nu\rho\sigma}$ **matches** the independent analytic Riemann (`cl10.nonlinear_gravity`) to $O(h^2)$, the **Schwarzschild** Ricci is $0$ (the vacuum equation $G=0$ is satisfied **by computation**) and the **computed** Kretschmann is $48M^2/r^6$ (the **black hole recovered**); the complete 0-parameter discrete EH action $S = \frac{1}{16\pi}\sum_x \sqrt{-\det g(x)}(R(x)-2\Lambda)$ ($G=1$) whose **functional derivative** (the numerical Euler–Lagrange) **is** the non-linear Einstein equation (dynamics $\leftarrow$ action); and the **non-linear** (state-dependent) QCA metric evolution (the 't Hooft non-linear completion: the GUT angle is modulated by the local $\Gamma\Gamma$ curvature — unitary/UV-finite, causal $c=1$, superposition fails) |
| the **3+1 BSSN/Z4 numerical relativity** (§6, C2, W9; V51) | `cl10.numrel` | the **complete** 3+1 NR: the 3+1 Ricci **0-row** $(R_{00}, R_{0i})$ (the **constraints** $H, M$ of the 4-D Einstein tensor) by **2 routes** — the **4-D Riemann** (`cl10.gravity`, the general curvature) and the **3+1** $(\gamma, K)$ decomposition ($K_{ij} = -\tfrac12\partial_0\gamma_{ij}$, the 3+1 Christoffels $\Gamma^{(3)}$) — agreeing to $O(h^2)$ (residual $\sim 10^{-13}$); the **BSSN/Z4** completion of the W9 gap: the **Z4** constraint-damping (a perturbed constraint decays to 0; the undamped case does not), the **non-linearity** (the Ricci 0-row is quadratic in $K$, so superposition fails), and the **UV-finite** property (the 3+1 quantities are bounded by the lattice cutoff) — the **classical** 3+1 numerical relativity is **complete** (the **open** part of W9 is now only the **quantum** completion) |
| the 4-D non-linear gravity (analytic de Sitter, V46) | `cl10.nonlinear_gravity` | the **0-parameter** de Sitter metric (flat FRW: $g_{00}=-1$, $g_{ij}=a^2\delta_{ij}$, $a=(1/H)e^{Hx^0}$); the **full 4-D Riemann tensor** (analytic Christoffel $\Gamma^0_{ij}=Ha^2\delta_{ij}$, $\Gamma^i_{0j}=\Gamma^i_{j0}=H\delta_{ij}$ and $\partial_0\Gamma$); the **Einstein tensor** $G_{\mu\nu}=R_{\mu\nu}-\tfrac12 g_{\mu\nu}R=-3H^2 g_{\mu\nu}$ (the 0-parameter cosmological constant); the **nonzero** $\Gamma\Gamma$ (cubic) self-gravity; at the Planck Hubble the 0-parameter vacuum energy $\Lambda\sim M_P^4$ (the W1 named residual).  This is the **analytic reference** that `cl10.gravity` (V50) cross-checks against (the 2-route verification) |
| the full SM QCA step (§11, P2) | `cl10.sm_qca` | $U = \prod_{a<b}(\cos g_{10}\, I - i\sin g_{10}\, T_{ab})$ with $T_{ab} = i\,h_a h_b$ (each $T_{ab}$ Hermitian, $T_{ab}^2 = I$) — the product is **exactly** unitary; no matrix exponential |
| the local QCA (§12, C1) | `cl10.local_qca` | the **internal** layer (the 45 rotations per cell) + the **staggered** 2-cell **swap-rotation** layer $U_{\rm bond} = \cos g_{10}\, I - i\sin g_{10}\, P$ along one spatial axis per time step ($P$ = the 0-parameter swap of two adjacent cells — a permutation, no angle) |
| the EH action from the dynamics (§12, C2) | `cl10.einstein.einstein_hilbert_from_dynamics` | the **discrete** quadratic action $S_{\rm EH} = \frac{1}{16\pi}\frac12\sum_{x,t}\, u_t\big[u_{t+1} - 2u_t + u_{t-1} - dt^2\,\nabla^2 u_t\big]$ ($G = 1$, $dt^2 = \tfrac14$, $\nabla^2$ = the 6-neighbour Laplacian); the **Euler–Lagrange equation** is the **numerical** functional derivative (the central finite difference of the action) evaluated on the GPU |
| the numbers of the root (§12, C3) | `cl10.unification` | $n_{\rm qubits} = \log_2 32 = 5$ $\Rightarrow$ $g_{10} = 1/5$, $\alpha_{10} = 1/25$; $g_0 = 1/2^{n_{\rm space}} = 1/2^3 = 1/8$ ($n_{\rm space} = 3$ = the spatial dimension of the $(1,3)$); the **down-running** is the 1-loop $1/a(M) = 1/a(M_{\rm UV}) - (b/2\pi)\ln(M_{\rm UV}/M)$ |
| the 3-family step (§12, C4) | `cl10.family` | the 96-D internal space $= \bigoplus_{k=1}^3 H_{32}$; $U_{96} = \mathrm{blockdiag}(U_{32}, U_{32}, U_{32})$ at the **same** $g_{10}$ (the family symmetry); the tribimaximal PMNS (V29) is a property of the **Yukawa** sector, not of the time-evolution step |
| the 3-gen Yukawa matrix (V45, W6, C4) | `cl10.yukawa_3x3` | the **democratic** 3×3 mass matrix $M = \|1\!\rangle\!\langle1\|$ (all-ones, rank 1, eigenvalues $\{3,0,0\}$: 1 heavy + a 2-fold light level) — **derived** from the **3-fold family symmetry** (the 3 spatial directions $S_a$ are equivalent: $S_a^2=-I$, $\{S_a,S_b\}=0$, equal norm, so the 10-Higgs couples the 3 families equally); the **TBM** ($\sin^2 = 1/3, 1/2, 0$) diagonalizes it (the 0-parameter PMNS); the 3-family 96-level spectrum $= 3 \times$ (the 32 levels); the **10-Higgs** $16\to16$-bar **Dirac** term (computed, **rank-16**, off-diagonal) is the 0-parameter mass; the **3-level** hierarchy (the 10-Higgs vev) and the **exact** CKM/PMNS are the **named** residual |
| the GUT leptogenesis (§12, C5) | `cl10.baryogenesis` | $\Delta B \sim \frac{28}{79}\cdot\frac{g_0^4}{16\pi}\cdot\eta_{CP}$ with the out-of-equilibrium kinematics $\Gamma/H \sim 6\times10^{-4}$; $\eta_{CP}$ is the **named** input (taken $10^{-3}$, the observed range) |
| the S-matrix (§12, C6) | `cl10.smatrix` | $S = U$ (the 32×32) and $S^{(2)} = U\otimes U$ (the 1024×1024) — the 0-parameter transition amplitudes of the 't Hooft program |
| the lattice RG / cross-sections (V44, W9, C6) | `cl10.lattice_rg` | the **UV-finite** statement: the Planck cutoff bounds every cross-section to $\sigma_{\max} = 4\pi / M_P^2$; the 0-parameter 2→2 cross-section $\sigma_0 = P_2 / M_P^2$ ($P_2 = 1 - \bar d^2$ from $S^{(2)} = U \otimes U$); the **0-parameter dispersion** (the 32-level spectrum $H = -i\log U$, UV-bounded $|\omega| \le \pi$, no 1-loop pairing); the **0-parameter optical theorem** (the inelastic $P_2$ **is** the cross-section, elastic + inelastic = 1); the **0-parameter Froissart-Martin bound** ($m_{\min}=M_P$: the QCA's *constant* cross-section is **below** the $\ln^2$ growth — the strongest UV-finiteness); the QCA/perturbative ratio $R = 4\pi P_2 / g_{10}^4$; the **named** residual (W9) is the matching of these to the measured cross-sections (the lattice RG) |
| the dark-matter scale (§12, C7) | `cl10.cosmology3` | the **soliton energy** $E$ (in Planck units) of $u = \sqrt{2}\,\mathrm{sech}\,x$ $\Rightarrow$ $m_{\rm DM} = E \cdot M_P$; the 0-parameter **freeze-out** $\Omega h^2 \sim 1/(\sigma m)$ with $\sigma \sim 1/m^2$ (the only scale in the theory) |
| the 2-loop running (§12, W4, V47) | `cl10.unification.running_2loop` | the **UV-finite** QCA: the 2-step $U^2$ is **unitary** (residual $\sim10^{-15}$) and the 2-particle S-matrix is **bounded** ($\max|S|\le1$); **vs** the **UV-divergent** perturbative 2-loop: the $b_2$ coefficient makes the U(1) Landau pole *closer* (the 1-loop pole is above $M_P$, the 2-loop pole is below $M_P$, $\sim10^{12}$ GeV) |
| the n-particle amplitudes (§12, W9, V48) | `cl10.smatrix.n_particle` | the **n-particle** S-matrix $S^{(n)} = U\otimes\cdots\otimes U$ is the **Kronecker** product of the 1-particle, so it is **unitary iff the 1-particle is** ($(UU^\dagger)\otimes\cdots\otimes(UU^\dagger)=I$); the 1- and 2-particle are computed on the GPU (the $32\times32$ and $1024\times1024$), the 3- and 4-particle ($32^3$, $32^4$) follow by the property — the **UV-finite** property (the n-step is bounded for all n) |
| the 0-parameter quark CKM (V49, W6) | `cl10.yukawa_3x3.quark_ckm` | the **quark** mass matrices (the up-type and down-type) on the **democratic** base (the 2-level structure, eigenvalues $\{3,0,0\}$); the **CKM** (the up/down mass-matrix misalignment, the unitary 3×3) is the **TBM** ($\sin^2 = 1/3, 1/2, 0$); the **deviation** to the observed CKM $(0.050, 0.949, 0.020)$ is the W6 named residual, now quantified |
| the honest parameter audit (§15, V52) | `cl10.parameter_audit` | the **solution to the four limits of the theory**: (i) the **2-route root** — the **symbolic, exact** $\mathrm{Cl}(10,0)$ of the `kingdon` library ($e_a^2=1$, $\{e_a,e_b\}=0$, $w^2=-1$) **vs** the **numeric** $32\times32$ symplectic-Pauli (the same identities, GPU); (ii) the **16/16-bar spectral mirror** — the **spectral form of CPT**: $[U,w]=0$ (the chiral split is preserved) and the **antimatter** spectrum is the **phase-conjugate** of the **matter** one, $\omega_{16\bar}=-\omega_{16}$, $\mathrm{spec}(U|_{16\bar})=\mathrm{spec}((U|_{16})^{\dagger})$ — **holds for every** 't Hooft angle (a property of the algebra, not of $g_{10}$) — + the smallest spectral level ($0.0115$); (iii) the **0-parameter CP-odd theorem** (the baryogenesis wall W5: every CP-odd invariant vanishes — the QCA satisfies none of the 3 Sakharov conditions, the exact value is the named $\eta_{CP}$); (iv) **$g_{10}$ is not fixed** by the 2-loop Landau-pole condition (the **free** 't Hooft angle); (v) the **hierarchy exponent** $-17$ and the **sub-Planck** DM scale are **not** in the **single-step** 32-level QCA spectrum (W8/W7) |
| the hierarchy exponent (§7, W8; V53) | `cl10.hier_exponent` | the **0-parameter** value $v/E_P = 2.02\times10^{-17}$ **is** the **10-loop** factor $(1/16\pi)^{10}\times 2.08$; the **10** $= 2\times5$ (the 16/16-bar **chiral split** $\times$ the **5 qubits** of the 32-D module); the **same** factor also gives the DM mass (V54) — **one factor, two results**; the **exponent** $-17$ is 0-parameter-derived (match to 0.18%); the **2-digit** coefficient is the named residual |
| the sub-Planck DM scale (§9, W7; V54) | `cl10.dm_wimp` | the **10-loop suppression** of the Planck scale: $m_{\rm DM} = M_P (1/16\pi)^{10} = 118$ GeV (the WIMP/weak scale; the **same** 10-loop factor as V53); the **V42** Planck-mass soliton is under-abundant (the $10^{-52}$ gap, quantified); the 0-parameter **freeze-out** $\Omega h^2 \sim 0.12$ (the WIMP miracle) |
| the baryon asymmetry (§9, W5; V55) | `cl10.eta_cp` | a **product of 0-parameter numbers**: $\Delta B = (28/79)\,\varepsilon_{CP}\,D$ (the SM sphaleron $28/79$, the **1-loop CP-odd** $\varepsilon_{CP} = g_0^4/16\pi$ with $g_0 = 1/8 = 1/2^3$, and the **2-loop dilution** $D = (1/16\pi)^2$); $\Delta B = 6.8\times10^{-10}$ (**ratio 1.12** to the observed $6.1\times10^{-10}$, improved from V40's 2.8); the **CP phase** and the 2-digit ratio are the named residual |
| the $g_{10}$ 1-loop data check (§12, W4/W6; V56) | `cl10.g10_selection` | the 0-parameter $g_{10} = 1/5$ (a number of the root: the 32-D module, 5 qubits) gives $\alpha_{10} = 1/25 = 0.04$, **within a factor 1.65** of the **1-loop** GUT coupling (the mean of the three running to $M_P$, $0.024$); the 1-loop matching is the W4 named residual |
| the 2-step seesaw on the 10-loop ladder (§13, W2; V57) | `cl10.neutrino` | the **0-parameter** 2-step seesaw: both scales are numbers of the root on the **10-loop ladder** — the **Dirac** $m_D = M_P(1/16\pi)^{10}$ (the WIMP/hierarchy scale, the *same* 10-loop as V53/V54) and the **2-step (Majorana)** $M_N = M_P(1/16\pi)^{3}$ (the intermediate $10^{14}$ scale, **replacing the hard-coded $10^{14}$** of V15/V18) $\Rightarrow$ $m_\nu = M_P(1/16\pi)^{2\times10-3} = M_P(1/16\pi)^{17} \approx 0.15$ eV (**2.9x** the observed, the 2-digit coefficient + $n_{2\text{step}}=3$ the named residual); the **2-step** $U^2$ is a 0-parameter **unitary** of the GUT QCA (the 2-loop) |
| the UV-finite *discrete RG* (§12, W4; V58) | `cl10.uvrunning` | the **W4 open part** (the lattice RG) made **explicit**: the **coarse-graining tower** $U^{2^p}$ ($p=1\ldots8$) of the GUT QCA is **unitary** (bounded) at **every** step (the GPU) — so the **discrete RG is UV-finite** (no Landau pole), **in contrast to** the *continuous* perturbative RG (2-loop $U(1)$ pole $\sim10^{12}$ GeV, **below** the cutoff); the 1-loop non-unification (the $U(1)$ pole, the spread $\sim0.6$) is a **feature of the non-SUSY SM** (not a defect); the **2-step** ($U^2$, $p=1$) is the **same** object as V57 (one object, two uses); the 1-loop spread is the named residual |


---

## 3. The engine and the numerical protocol

**Backend selection** (at import, in `cl10/engine.py`): **CuPy** (GPU)
$\to$ **JAX** (XLA/GPU) $\to$ **NumPy** (CPU).  The active backend is
reported by `engine._BACKEND` and `engine.backend_info()`.  On the
reference machine (an **RTX 3090**, CUDA 13.4) the backend is **CuPy**
$14.2.0$; JAX $0.11.2$ sees $2$ CUDA devices.  All modules and all
experiments go through this **one** core, so the whole theory runs on the
GPU out of the box and degrades cleanly to CPU (the $32\times32$ matrices
are small).
**The second route (the symbolic cross-check).**  The root
$\mathrm{Cl}(10,0)$ is verified by a **second, independent** route: the
**symbolic, exact** Clifford algebra of the `kingdon` library (the
`Algebra(10,0)` frame, where $e_a^2=1$, $\{e_a,e_b\}=0$, and $w^2=-1$ hold
**exactly**, with no float), cross-checked against the **numeric**
$32\times32$ symplectic-Pauli representation of `clifford` on the GPU
(V52, `cl10.parameter_audit.kingdon_2route`).  The two agree — the root is
not an artifact of one construction.  `kingdon` is used for **symbolic**
algebraic verification only (the heavy numerics stay on the CuPy engine).


**The protocol.**  Each experiment (V0–V58 under `experiments/`):
1. constructs the relevant object from the **0-parameter** input (no fit);
2. computes the **residual** against the 0-parameter prediction;
3. writes a **JSON verdict** `{id, tier, property, pass, residual,
   status, data, wall_time}` to `experiments/<id>/results/`;
4. **passes** iff the residual is below the stated tolerance
   (`tol`; default $\varepsilon = 10^{-12}$).

`run_validation.py` runs all **59** experiments (V0–V58) in name order,
reads the JSONs, and writes `results/validation_summary.json` with the
per-experiment verdicts, the **reduced walls**, the **named residuals**,
and the **TOE verdict**:

$$\mathrm{TOE} = \bigwedge_{V \in \{0..7,\,11..13\}} \mathrm{PASS}(V)
\quad\text{(the hard tiers 0–3.)}$$

**The result (59/59, 2026-10-03).**  Every experiment **passes** on the
GPU; the **TOE verdict is positive**.  (For the *meaning* of each pass —
which walls are reduced and which residuals are named — see
`doc/publication.md` Sections 11–13.)

**The 0-parameter numbers** (no fit anywhere): the hierarchy
$2.02\times10^{-17}$, the $5/3$, the B–L $(-1/+1)$, the 16/16-bar, the
$5/24$ $\mathcal{T}$ normalization, the $G = 1$, the causal radius $1$,
the $1.22\times10^{19}$ GeV GUT scale, the $1/25$ GUT coupling, the $246$
GeV electroweak scale, and the off-diagonal **10-Higgs** mass term
$16\,10\,16\bar\to 1$ (V20).  **From the completion program (V36–V51):**
$g_{10} = 1/5$ (the **one free** 't Hooft angle) and $g_0 = 1/8$ (the
$1/2^{n_{\rm space}}$ mapping; the 3-space, V38) — the **1-parameter and
convention** status of these two couplings is **audited** in V52; the
tribimaximal angles
$\sin^2\theta_{12} = 1/3$, $\sin^2\theta_{23} = 1/2$, $\theta_{13} = 0$
$(\mathrm{V29/V39/V45})$; the SM sphaleron factor $28/79$ (V40); the 0-parameter EH
coefficient $1/16\pi$ and the lattice CFL $dt^2 = 1/4$ (V37); the **UV-finite** 2→2
cross-section $\sigma_0 = P_2 / M_P^2$, the **0-parameter dispersion** (the 32-level spectrum, UV-bounded) and the **0-parameter optical theorem** (the inelastic $P_2$ **is** the cross-section), and the QCA/perturbative ratio
$R = 4\pi P_2 / g_{10}^4$ (V44, W9); and the **complete** 4-D non-linear gravity (V50): the **general** Riemann computed from an *arbitrary* metric (the 2-route analytic/computed agreement $6.5\times10^{-7}$), the **black hole** recovered by computation (the computed Kretschmann $48M^2/r^6$, the vacuum $G=0$), and the **non-linear** (state-dependent) QCA metric evolution.
the **complete 3+1 BSSN/Z4 numerical relativity** (V51): the 2-route 0-row (the **constraints** $R_{00}, R_{0i}$ by the 4-D Riemann vs the 3+1 $(\gamma, K)$ decomposition, agreeing to $O(h^2)$), the **Z4** constraint-damping, the **non-linearity** (superposition fails), and the **UV-finite** property.

---

## 4. Reproducibility

**Environment.**  Python $\ge 3.11$; a **venv** (`.venv`) with the
`pyproject.toml` dependencies: `kingdon`, `einops`, `jax[cuda13]`,
`cupy-cuda13x`, `numpy`.  A GPU (CUDA $\ge 12$) is required for the
CuPy/JAX backends; without it the engine **degrades** to NumPy (CPU) and
all experiments still run (the $32\times32$ matrices are small).

**Run:**

```bash
cd CL10
source .venv/bin/activate
python3 run_validation.py            # -> results/validation_summary.json (59/59, TOE)
python3 run_validation.py V11 V12    # a subset
```

**Lint.**  `python3 -m ruff check cl10/ experiments/ run_validation.py`
passes (the compact numerical style allows E702/E402/E741; the real bugs
F4xx/F5xx are enforced).

**Figures.**  The nine figures of `doc/publication.md` (the "Figures"
section) are **computed, not hand-drawn**, by `make_figures.py` (at the
repo root).  It imports the `cl10` modules and writes PNGs to
`doc/figures/`:

```bash
python3 make_figures.py            # -> doc/figures/fig1..fig9_*.png
```

The nine, by module: the **32-level QCA spectrum** (the one object;
`cl10.qca_spectrum`), the **1-loop running** of the three SM couplings to
$M_P$ ($1/\alpha_i$, with the 2-loop Landau-pole contrast;
`cl10.unification`), the **10-loop ladder** (the hierarchy $+$ the DM mass;
`cl10.hier_exponent`), the **WIMP miracle** (`cl10.dm_wimp`), the **2-route
Kretschmann** (`cl10.gravity`), the **512/512 cancellation**
(`cl10.cosmology3`), **QCD confinement** (the 0-parameter linear vs
Coulombic potential; `cl10.qcd`), **GUT leptogenesis** (the baryon asymmetry
as a product of 0-parameter numbers; `cl10.eta_cp`), and **anomaly
cancellation** (the $[\mathrm{U}(1)]^3$ per-field divergence;
`cl10.matter`).  It uses the same backend selection as Section 3
(CuPy/GPU, degrading to NumPy/CPU).

---

## 5. The status of two theorems (the honest wording)

**The Nielsen–Ninomiya status.**  The lattice **chirality** is **not**
realised by a naive lattice Dirac operator (which would double); it is a
**representation** of the root: the 16 Weyl is a **module** of
$\mathrm{Cl}(10)$ (the pseudoscalar splits the 32 into 16/16-bar), and the
QCA is a 't Hooft-type **finite-state** theory (one **exactly-unitary**
per-cell step).  The **doubling** theorem (which constrains a naive lattice
Dirac) does **not** apply to this construction — the chiral
content is a 0-parameter property of the **algebra**, not of a lattice
operator.

**The CPT status.**  The CPT map (V24) is a **0-parameter consequence of
the $\mathbb{Z}_2$ grading**: the **45** GUT generators are **even** in
the grading (they **commute** with the pseudoscalar), while the **time**
direction is **odd** (it **anticommutes**).  The "CPT theorem" of the
theory is this **grading identity** — stated as the *status* of the
theorem (a real physicist would call it a grading property, not a deep
theorem).  The **spectral** form of this CPT is **computed** in V52
(`cl10.parameter_audit.spectral_mirror`): because the dynamics commute
with the pseudoscalar ($[U,w]=0$, the chiral split is preserved), the
**antimatter** ($16$-bar) QCA spectrum is the **mirror** (phase-conjugate)
of the **matter** ($16$) one — $\omega_{16\bar} = -\omega_{16}$, i.e.
$\mathrm{spec}(U|_{16\bar}) = \mathrm{spec}((U|_{16})^{\dagger})$ — and this
holds for **every** 't Hooft angle (a property of the algebra, not of the
chosen $g_{10}$): **CPT $=$ spectral conjugation**.
