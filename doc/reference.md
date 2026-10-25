# CL10 — Technical Reference

A **professional** reference for the **implementation** of Causal Clifford
Theory on the GUT root (CL10).  It fixes the **exact** constructions used
by the code, the **numerical** protocol, and the **reproducibility**
information.

The **physics** (the theorems, the claims, the verification program, the
ten open problems O1–O10, the falsifiers F1–F10) lives in
`doc/final_publication.md`; the **Python interface** lives in `doc/api.md`.
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

Each object of `doc/final_publication.md` is computed by **one** module (the
signatures are in `doc/api.md`):

| object (final_publication.md) | module | the exact construction in the code |
|---|---|---|
| the root $\mathrm{Cl}(10)$ (§2) | `cl10.clifford` | the 10 generators of Appendix A of the publication (the 5-qubit symplectic Pauli pattern) |
| the $(1,3)$ inside the GUT (§3) | `cl10.spacetime` | $T = h_1$ (time, $+$), $S_a = i\,h_{a+1}$ for $a = 1,2,3$ (space, $-$) — **four of the ten** generators, so the Minkowski relations are **identities** of the root |
| the 16 Weyl / 16-bar (§4) | `cl10.matter`, `cl10.cpt` | the **pseudoscalar** $w = h_1\cdots h_{10}$ (anti-Hermitian, $w^2 = -I$) diagonalized with a general `eig` (not `eigh` — $w$ is not Hermitian): the **$+i$** eigenspace = the 16, the **$-i$** = the 16-bar |
| the 45 and the $5/3$ (§5) | `cl10.gauge` | the 45 $= \{h_a h_b : a < b\}$; the hypercharge in the $\mathbf{5}$ is $\mathcal{T} = \mathrm{diag}(-\tfrac16,-\tfrac16,-\tfrac16,+\tfrac14,+\tfrac14)$ with $\mathrm{Tr}(\mathcal{T}^2) = 5/24$; GUT-normalized to $\tfrac12$, with $Y = 2\mathcal{T}$ $\Rightarrow$ $g_1^2/g'^2 = 8 \cdot 5/24 = 5/3$ |
| the 8 anomalies (§4, Table 1) | `cl10.matter.anomalies` | the sums over the charges of the 6 fields of the 16 (fixed by the $\mathrm{SU}(5)$ embedding), with the Dynkin indices $T(3) = T(2) = \tfrac12$ |
| the metric + the graviton (§11) | `cl10.qca` | the **6-neighbour** discrete Laplacian on an $L^3$ torus; the linearised perturbation is a free wave (the discrete CFL $dt = h/2$ gives $c = 1$, the momentum is UV-bounded) |
| the non-linear metric, the kink, the black hole (§11) | `cl10.einstein` | the **cubic** $\phi^3$ (Choptuik) completion: the static $u'' = u - u^3$ has the **exact** finite-energy solution $u = \sqrt{2}\,\mathrm{sech}\,x$; the Ricci scalar is **computed** on prescribed metrics (flat: $0$; the 2-sphere: $2/r^2$; the Schwarzschild: the Kretschmann $48M^2/r^6$) |
| the **complete** 4-D non-linear gravity (§11, O9; V50) | `cl10.gravity` | the **general** 4-D Riemann / Ricci / Einstein / Kretschmann **computed by central finite differences from an *arbitrary* metric** $g_{\mu\nu}(x)$ (no hard-coded Christoffels, not restricted to one solution); **2-route verified**: the computed de Sitter $R^{\mu}_{\nu\rho\sigma}$ **matches** the independent analytic Riemann (`cl10.nonlinear_gravity`) to $O(h^2)$, the **Schwarzschild** Ricci is $0$ (the vacuum equation $G=0$ is satisfied **by computation**) and the **computed** Kretschmann is $48M^2/r^6$ (the **black hole recovered**); the discrete EH action $S = \frac{1}{16\pi}\sum_x \sqrt{-\det g(x)}(R(x)-2\Lambda)$ ($G=1$, a **convention**: Newton's constant is an input of this sector) whose **functional derivative** (the numerical Euler–Lagrange) **is** the non-linear Einstein equation (dynamics $\leftarrow$ action); and the **non-linear** (state-dependent) QCA metric evolution (the 't Hooft non-linear completion: the GUT angle is modulated by the local $\Gamma\Gamma$ curvature — unitary/UV-finite, causal $c=1$, superposition fails) |
| the **3+1 BSSN/Z4 numerical relativity** (§11, O9; V51) | `cl10.numrel` | the **3+1** Ricci **0-row** $(R_{00}, R_{0i})$ (the **constraints** $H, M$ of the 4-D Einstein tensor) by **2 routes** — the **4-D Riemann** (`cl10.gravity`, the general curvature) and the **3+1** $(\gamma, K)$ decomposition ($K_{ij} = -\tfrac12\partial_0\gamma_{ij}$, the 3+1 Christoffels $\Gamma^{(3)}$) — agreeing to $O(h^2)$ (residual $\sim 10^{-13}$); the **Z4** constraint-damping (a perturbed constraint decays to 0; the undamped case does not), the **non-linearity** (the Ricci 0-row is quadratic in $K$, so superposition fails), and the **UV-finite** property (the 3+1 quantities are bounded by the lattice cutoff) — the **classical** 3+1 sector is therefore **implemented and cross-checked**; the open part of O9 remains the **sum over geometries**, the **lattice chirality check** and the **quantum** completion |
| the 4-D non-linear gravity (analytic de Sitter, §11, §12.1) | `cl10.nonlinear_gravity` | the de Sitter metric (flat FRW: $g_{00}=-1$, $g_{ij}=a^2\delta_{ij}$, $a=(1/H)e^{Hx^0}$); the **full 4-D Riemann tensor** (analytic Christoffel $\Gamma^0_{ij}=Ha^2\delta_{ij}$, $\Gamma^i_{0j}=\Gamma^i_{j0}=H\delta_{ij}$ and $\partial_0\Gamma$); the **Einstein tensor** $G_{\mu\nu}=R_{\mu\nu}-\tfrac12 g_{\mu\nu}R=-3H^2 g_{\mu\nu}$ (the vacuum energy **of this metric**, not the observed $\Lambda$); the **nonzero** $\Gamma\Gamma$ (cubic) self-gravity; at the Planck Hubble the vacuum energy $\Lambda\sim M_P^4$ (the O1 named residual).  This is the **analytic reference** that `cl10.gravity` (V50) cross-checks against (the 2-route verification) |
| the internal layer of the step (§8.1) | `cl10.sm_qca` | $U = \prod_{a<b}(\cos g_{10}\, I - i\sin g_{10}\, T_{ab})$ with $T_{ab} = i\,h_a h_b$ (each $T_{ab}$ Hermitian, $T_{ab}^2 = I$) — the product is **exactly** unitary; no matrix exponential |
| the automaton: the block cell, the step, the causal radius (§8.1, Theorem 5) | `cl10.supercell` | the cell is the **2×2×2 block** (8 sub-sites × the 36-level cell module of `cl10.cell` = $\mathbb{C}^{288}$); one **time-independent step** of six stages, each the internal layer on every sub-site plus one **bond stage** $\cos g_{10}\, I - i\sin g_{10}\,P$ on the four disjoint pairs of one parity along one axis (three intra-block, three inter-block); the exact Bloch matrix $M(q)$ built on the spinor block ($256\times256$; $128\times128$ on the 16) on the block zone, and the full cell transport factorizing as $\operatorname{kron}(P(q), U_{\rm cell}^6)$ ($288\times288$, V196); the measured front speed is **exactly 1 Planck cell per stage** |
| the cell module of Axiom 1 (§2, §4; V196) | `cl10.cell` | the 36-level internal module $S \oplus R_A \oplus R_B$: the 32-D Clifford module and the two neutral singlet conjugate pairs carrying the rung phases $e^{\pm i\,\mathrm{loop}^{n_A}}$, $e^{\pm i\,\mathrm{loop}^{n_B}}$ ($n_A, n_B$ from `CFG.cell`); $\pm$ pairing to $2\times10^{-15}$, $\det$ phase to $3\times10^{-16}$, quasi-energy sum to $2\times10^{-16}$ on the full 36; the B-parity and its exact commutation with step and $45$ |
| the **closure of the step** (§8.5; V131–V135) | `cl10.sm_qca`, `cl10.supercell` | the same 45-plane product at a **general** angle $g$; the closure test is $U^{k} = \pm I$ on the module (internal layer) and $M(0)^2 = I$ on the exact block Bloch matrix at the zone centre. The **only** exactly closed angles of the canonical Appendix-A ordering are $g = \pi/2$ (the step **is** the pseudoscalar $w$, ordering-independent) and $g = \pi/4$ (order 8; internal spectrum exactly $\{\pm\pi/4\,(10),\ \pm3\pi/4\,(6)\}$, leak below $5\times10^{-16}$); five shuffled factor orderings break it (deviation $\sim1$) — the ordering, already an input of the step, carries the quantization. The closed-angle consequences are measured, not assumed: the gapless census (V132, `cl10.supercell`), the $\pi/2$ branch (V134), the desert re-closure at $\alpha_X = \pi/64$ (V133, `cl10.unification` + the V64 inverse problem), the coupling contact against the one-loop mean (V135, `cl10.g10_selection`) |
| the EH action from the dynamics (§11) | `cl10.einstein.einstein_hilbert_from_dynamics` | the **discrete** quadratic action $S_{\rm EH} = \frac{1}{16\pi}\frac12\sum_{x,t}\, u_t\big[u_{t+1} - 2u_t + u_{t-1} - dt^2\,\nabla^2 u_t\big]$ ($G = 1$, $dt^2 = \tfrac14$, $\nabla^2$ = the 6-neighbour Laplacian); the **Euler–Lagrange equation** is the **numerical** functional derivative (the central finite difference of the action) evaluated on the GPU |
| the numbers of the root (§8.1, §14) | `cl10.unification` | $n_{\rm qubits} = \log_2 32 = 5$ $\Rightarrow$ $g_{10} = 1/5$, $\alpha_{10} \equiv g_{10}^2 = 1/25$; $g_0 = 1/2^{n_{\rm space}} = 1/2^3 = 1/8$ ($n_{\rm space} = 3$ = the spatial dimension of the $(1,3)$).  Both readings are **conventions** (Table 5 of the publication), not derivations; the **down-running** is the 1-loop $1/a(M) = 1/a(M_{\rm UV}) + (b/2\pi)\ln(M_{\rm UV}/M)$ |
| the 3-family step (§4, §11.5) | `cl10.family` | the 96-D internal space $= \bigoplus_{k=1}^3 H_{32}$; $U_{96} = \mathrm{blockdiag}(U_{32}, U_{32}, U_{32})$ at the **same** $g_{10}$ (the family symmetry); the tribimaximal PMNS is a property of the **Yukawa** sector, not of the time-evolution step; only $\sin^2\theta_{12}=1/3$ is a **derived number** (§11.5) |
| the 3-gen Yukawa matrix (§13, V45, O6) | `cl10.yukawa_3x3` | the **democratic** 3×3 mass matrix $M = \|1\!\rangle\!\langle1\|$ (all-ones, rank 1, eigenvalues $\{3,0,0\}$: 1 heavy + a 2-fold light level) — obtained from the **3-fold family symmetry** (the 3 spatial directions $S_a$ are equivalent: $S_a^2=-I$, $\{S_a,S_b\}=0$, equal norm, so the 10-Higgs couples the 3 families equally); the **TBM** diagonalizes it (a **regression check**): only $\sin^2\theta_{12}=1/3$ is a derived number, $\sin^2\theta_{23}=1/2$ and $\theta_{13}=0$ are the $(23)$-adapted choice in the $O(2)$ family; the 3-family 96-level spectrum $= 3 \times$ (the 32 levels); the **10-Higgs** $16\to16$-bar **Dirac** term (computed, off-diagonal, full-rank $16\times16$ block of norm 4, both diagonal blocks $<10^{-15}$) is the Dirac mass; the **3-level** hierarchy (the 10-Higgs vev) and the **exact** CKM/PMNS are the **named** residual |
| the GUT leptogenesis (§13.3, O5) | `cl10.baryogenesis` | $\Delta B \sim \frac{28}{79}\cdot\frac{g_0^4}{16\pi}\cdot\eta_{CP}$; $\eta_{CP}$ is the **named** input (taken $10^{-3}$).  **Honest:** the out-of-equilibrium Sakharov condition **fails** at the model's own couplings — at $M_N = 9.6\times10^{13}$ GeV, $\Gamma_N/H \simeq 4.6$ for $y = 1/8$ ($1.9\times10^{2}$ for $y=1$), and the condition needs $y \lesssim 0.06$ |
| the S-matrix (§7, O9) | `cl10.smatrix` | $S = U$ (the 32×32) and $S^{(2)} = U\otimes U$ (the 1024×1024) — the transition amplitudes of the 't Hooft program (no interaction vertex: $S^{(n)}$ is the free Kronecker product) |
| the lattice RG / cross-sections (§13.1, O9) | `cl10.lattice_rg` | the **UV-finite** statement: the Planck cutoff bounds every cross-section to $\sigma_{\max} = 4\pi / M_P^2 = 8.4\times10^{-38}$ GeV$^{-2}$; the 2→2 cross-section $\sigma_0 = P_2 / M_P^2 = 5.3\times10^{-39}$ GeV$^{-2} = 2.1\times10^{-66}$ cm$^2$ ($P_2 = 0.789$ from $S^{(2)} = U \otimes U$); the **dispersion** (the 32-level spectrum $H = -i\log U$, UV-bounded $|\omega| \le \pi$, the 0-point $=0$ from the $\pm$-pairing, $\det U = 1$); the **optical theorem** as implemented is the unitarity bookkeeping (elastic $1-P_2$ + inelastic $P_2 = 1$), not an independent dynamical statement; the **Froissart–Martin** comparison ($m_{\min}=M_P$: the QCA's *constant* cross-section stays **below** the $\ln^2 s$ growth); the QCA/perturbative ratio $R = 4\pi P_2 / g_{10}^4$; the **named** residual (O9) is the matching of a Planck-suppressed $\sigma_0$ to the measured cross-sections (the lattice RG) |
| the dark-matter scale (§12.2, O7) | `cl10.cosmology3` | the **kink energy** $E$ (in Planck units) of $u = \sqrt{2}\,\mathrm{sech}\,x$ (a planar domain wall: not a carrier, O7) $\Rightarrow$ $m_{\rm DM} = E \cdot M_P$; the freeze-out abundance is **not computed** — the code normalises $\Omega h^2$ to the observed value, and the required $\langle\sigma v\rangle \simeq 1.4\times10^{-8}$ GeV$^{-2}$ is not derived from the root |
| the 2-loop running (§13.1, O4) | `cl10.unification.running_2loop` | the **UV-finite** QCA: the 2-step $U^2$ is **unitary** (residual $\sim10^{-15}$) and the 2-particle S-matrix is **bounded** ($\max|S|\le1$); **vs** the **UV-divergent** perturbative 2-loop: the $b_2$ coefficient makes the U(1) Landau pole *closer* (run from the measured $M_Z$ input the 1-loop pole is at $2.196\times10^{41}$ GeV, above $M_P$; the 2-loop pole is at $2.473\times10^{12}$ GeV, below $M_P$; `V38`) |
| the n-particle amplitudes (§7, O9; V48) | `cl10.smatrix.n_particle` | the **n-particle** S-matrix $S^{(n)} = U\otimes\cdots\otimes U$ is the **Kronecker** product of the 1-particle, so it is **unitary iff the 1-particle is** ($(UU^\dagger)\otimes\cdots\otimes(UU^\dagger)=I$); the 1- and 2-particle are computed on the GPU (the $32\times32$ and $1024\times1024$), the 3- and 4-particle ($32^3$, $32^4$) follow by the property — the **UV-finite** property (the n-step is bounded for all n) |
| the quark CKM, honestly (§13, V49, O6) | `cl10.yukawa_3x3.quark_ckm` | the **quark** mass matrices (the up-type and down-type) on the **democratic** base (the 2-level structure, eigenvalues $\{3,0,0\}$); but with $H_u = H_d$ both democratic the computed **CKM is the identity**, not the TBM — a non-trivial CKM needs an up/down **misalignment** the democratic base does not supply (the observed $\lambda = 0.225$ is not reproduced); that, plus $\sin^2\theta_{23}$ and $\sin^2\theta_{13}$ (which sweep $[0,2/3]$ under the $O(2)$ freedom of the 2-fold light subspace), is the O6 named residual |
| the honest parameter audit (§15; V52) | `cl10.parameter_audit` | the **audit of the four limits of the theory**: (i) the **2-route root** — the **symbolic, exact** $\mathrm{Cl}(10,0)$ of the `kingdon` library ($e_a^2=1$, $\{e_a,e_b\}=0$, $w^2=-1$) **vs** the **numeric** $32\times32$ symplectic-Pauli (the same identities, GPU); (ii) the **16/16-bar spectral mirror** — the **spectral form of CPT**: $[U,w]=0$ (the chiral split is preserved) and the **antimatter** spectrum is the **phase-conjugate** of the **matter** one, $\omega_{16\bar}=-\omega_{16}$, $\mathrm{spec}(U|_{16\bar})=\mathrm{spec}((U|_{16})^{\dagger})$ — **holds for every** 't Hooft angle (a property of the algebra, not of $g_{10}$) — + the smallest spectral level ($0.0115$); (iii) the **CP-odd theorem** (the baryogenesis wall O5: every CP-odd invariant vanishes — the QCA satisfies none of the 3 Sakharov conditions, the exact value is the named $\eta_{CP}$); (iv) **$g_{10}$ is not fixed** by the 2-loop Landau-pole condition (the pole is below $M_P$ for every $g_{10}$ in the perturbative range) — the **one** 't Hooft angle of the dynamics, and the *selected* value $1/5$ is a **convention** (the data-optimal candidate of V56 is $n=7$, not $n=5$); (v) the **hierarchy exponent** $-17$ and the **sub-Planck** DM scale are **not** in the **single-step** 32-level QCA spectrum (O8/O7).  `provenance_table()` classifies 15 headline numbers as 7 algebraic / 1 free angle / 2 conventions / 1 coincidence / 4 named residuals |
| the hierarchy exponent (§12.1, O8; V53) | `cl10.hier_exponent` | $v/E_P = 2.02\times10^{-17}$ **is** the **10-loop** factor $(1/16\pi)^{10}\times 2.08$; the **10** $= 2\times5$ (the 16/16-bar **chiral split** $\times$ the **5 qubits** of the 32-D module) is a **representation**, not a loop calculation; the **same** factor also gives the DM mass (V54) — **one factor, two statements of one scale**; the exponent is **selected** so that the 2-digit coefficient lands on the data (the $0.18\%$ match fixes the coefficient, not the exponent); the coefficient $2.08$ is **fitted** |
| the sub-Planck DM scale (§12.2, O7; V54) | `cl10.dm_wimp` | the **10-loop suppression** of the Planck scale: $m_{\rm DM} = 1.18\,M_P (1/16\pi)^{10} = 118$ GeV (the WIMP/weak scale; the **same** 10-loop factor as V53, the $1.18$ **fitted**); the earlier **Planck-mass soliton** reading is under-abundant (the $10^{-52}$ gap, quantified); the abundance $\Omega h^2$ is **not computed** — the code normalises it to the observed $0.120$ at $100$ GeV, so there is **no WIMP-miracle prediction**, only a mass scale |
| the baryon asymmetry (§13.3, O5; V55) | `cl10.eta_cp` | a **product of factors of three provenances**: $\Delta B = (28/79)\,\varepsilon_{CP}\,D$ (the SM sphaleron $28/79$ — an **input of the SM sector**, the LH d.o.f. count; the **1-loop CP-odd** $\varepsilon_{CP} = g_0^4/16\pi$ with $g_0 = 1/8 = 1/2^3$, a dimensional guess; and the **2-loop dilution** $D = (1/16\pi)^2$, itself a **guess**: the BDP washout at the same Yukawa is $14\times$ lower, V106); $\Delta B = 6.8\times10^{-10}$ (**ratio 1.11** to the observed $6.12\times10^{-10}$) is the **unwashed** product — the physical washout $\kappa_f = 0.12$ brings it to $2.1\times10^{-10}$ (**ratio 0.35**); the **CP phase** and the 2-digit ratio are the named residual |
| the $g_{10}$ 1-loop data check (§13.1, §15; V56) | `cl10.g10_selection` | $g_{10} = 1/5$ (the 32-D module, 5 qubits) gives $\alpha_{10} = 1/25 = 0.04$, **within a factor 1.73** of the **1-loop** GUT coupling (the mean of the three running to $M_P$, $0.0231$, from the PDG inputs $\alpha_1(M_Z) = 0.0169$, $\alpha_2 = 0.0338$, $\alpha_3 = 0.118$); the 1-loop matching is the O4 named residual; **honest scan**: the data-optimal candidate is $n=7$ ($\alpha = 1/49 = 0.020$, relative deviation $0.12$), **not** the root $n=5$ (deviation $0.73$) — the data check confirms the **order only** and does **not** select $n$ |
| the 2-step seesaw on the 10-loop ladder (§12.3, O2; V57) | `cl10.neutrino` | the 2-step seesaw on the **loop ladder**: the **Dirac** $m_D = M_P(1/16\pi)^{10}$ (the WIMP/hierarchy scale, the *same* 10-loop as V53/V54) and the **2-step (Majorana)** $M_N = M_P(1/16\pi)^{3} = 9.6\times10^{13}$ GeV (**replacing the earlier hard-coded $10^{14}$**) $\Rightarrow$ $m_\nu = M_P(1/16\pi)^{2\times10-3} = M_P(1/16\pi)^{17} \approx 0.146$ eV — **2.9×** the atmospheric $0.050$ eV and **in tension** with $\sum m_\nu < 0.12$ eV; $n_{2\text{step}}=3$ and the 2-digit coefficient are **fitted**; the **2-step** $U^2$ is a **unitary** of the GUT QCA; **V136** shows the exponent is not a link-word length: the $\Delta L=2$ operator is a grade-5 (126-type) element of the root's own $1024$-D algebra with a 1-dimensional weight space, while products of the grade-2 adjoint link factors are always even-grade (the multiplication table shifts grade by $\pm2$ exactly), so no word in the 45 reaches it |
| the UV-finite *discrete RG* (§13.1, O4) | `cl10.uvrunning` | the **O4 open part** (the lattice RG) made **explicit**: the **coarse-graining tower** $U^{2^p}$ ($p=1\ldots8$) of the GUT QCA is **unitary** (bounded) at **every** step (the GPU) — so the **discrete RG is UV-finite** (no Landau pole), **in contrast to** the *continuous* perturbative RG (2-loop $U(1)$ pole $\sim10^{12}$ GeV, **below** the cutoff; and run *down* from $\alpha_{10}=1/25$ the $\mathrm{SU}(3)$ coupling hits its pole at $2.2\times10^9$ GeV, §11.1); the 1-loop non-unification (the $U(1)$ pole, the spread $\sim0.47$) is a **failure of the coupling values**, not a feature to be dismissed; the **2-step** ($U^2$, $p=1$) is the **same** object as V57 (one object, two uses); the 1-loop spread is the named residual |
| the **four limits and their common object** (§14; V59) | `cl10.four_limits` | the four "limits" of the theory are the **four projections of one invariant** $10 = 2\times5$ (the chiral $16/16$-bar split $\times$ the 5 qubits of the 32-D module) **plus one irreducible CP phase** — computed on the GPU: (L1) the 7 algebraic invariants are **$g_{10}$-invariant** (the CPT-mirror residual stays $\sim1.8\times10^{-15}$ for $g_{10}=0.1..0.5$); (L2+L4) the **same** $n=10$ is the **unique** $n\in1..20$ within **one decade** of **both** closures: the hierarchy ($2.08(1/16\pi)^{10}$) and the WIMP mass ($M_P(1/16\pi)^{10}=118$ GeV) — **two statements of one scale**; (L3) the **rank-1** (1 $N$) Yukawa has **0** CP phases (no baryon) while the **2-step** (2 $N$) has **1** irreducible CP phase ($\mathrm{Im}\det(Y)=|\sin\delta|$) — so baryogenesis needs **one** phase, $\Delta B=6.8\times10^{-10}$ (ratio 1.12).  **Bookkeeping, not a derivation:** the algebra has no free parameter, but $g_{10}=1/5$ and $g_0=1/2^3$ are **conventions**, and the named residual set is **2 two-digit coefficients** (2.08, 1.18, read from the data) **+ 1 CP phase + 1 't Hooft angle** |
| the **3-D hydrogen atom** driven by the CL10 photon (Fig. 2; V62) | `simulation.hydrogen` | the hydrogen atom (the **real** Rydberg levels $E_n=-13.606/n^2$ and the eigenfunctions $R_{nl}(r)Y_{lm}$ drawn on a $296^3$ lattice — an **input of atomic physics**, NOT derived from the root) coupled to the **CL10 free massless 3-D field** (the 6-neighbour leapfrog FDTD, the same engine as `simulation.young`): the **Fermi-golden-rule** dipole drive $A\,|F|^2\,L(\delta)$ with a **Lorentzian** line shape (a photon *at* a gap is absorbed, one far from every gap mostly is not — the **tail is real**, not zero), the **Einstein-A** emission (the transition dipole sources an outgoing wave at the gap frequency and the population relaxes), the **sponge** absorbing boundary, and the observables (`populations`, `atom_energy`, `excited_fraction`, `electron_density` $=|\psi|^2$, `photon_density` $=|F|^2+|F_{\rm out}|^2$).  The module's own `dipole_mat` measures: the gaps equal the Rydberg differences (relative error $0$), the **E1 selection rules come out of the lattice wavefunctions** — the parity-forbidden elements sit at $2.5\times10^{-18}$ ($1s|z|1s$), $4.6\times10^{-18}$ ($2p|z|2p$), $2.7\times10^{-17}$ ($3p|z|3p$) and $5.8\times10^{-19}$ ($1s|z|3d$) in cells, against $13.64$ ($1s|z|2p$), $44.60$ ($2p|z|3d$), $-57.50$ ($3p|z|3d$) and $-2.24$ ($1s|z|3p$) for the allowed ones — and the whole table is recorded, in `results/hydrogen_transition_table.json`, written by `python -m simulation.hydrogen --table-only` —, the absorption **peaks on the gap and falls with the detuning**, and the excited atom **radiates while it decays**.  The knobs $A$, $h$, $\gamma_{\rm sp}$ are physical parameters of the atom, **not fitted**.  **V62** is the production run ($296^3$, the full photon script) that regenerates the Fig. 2 figure. |
| the Young wave + the **single-electron** measurement / collapse (Fig. 1, §6) | `simulation.young` | the **3-D QCA wave** (the same discrete wave equation as the CL10 graviton, §11) run through a **double-slit barrier** as a **single-electron measurement** (the Tonomura flow): `wave.WaveSim` is the FDTD (the `r=0.5` CFL leapfrog, the source matched to the FDTD dispersion, the quadratic absorbing sponge) that fires **one localized electron at a time** — each one diffracts through the slits and is detected at a **single point** sampled from the **time-integrated** on-screen $\psi^{2}$ (the interference pattern the screen records; a single instant of the traveling wave is only a blob, so the detection is integrated over the electron's passage), with the **sponge** absorbing the wave between emissions (no reset); the **`collapse.Collapse`** class implements the **measurement** — `sample` = the **Born-rule** inverse-CDF sampling of one detection point, `contract` = the **wave-packet reduction** ($\psi\to\psi\times$ a shrinking Gaussian centred on the detected point: the **non-unitary** projection, the collapse), and `add` = the **detection accumulation** whose points cluster into the fringe groups.  The collapse is an **input of the interpretation** (O9), not a consequence of the unitary step |


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

**The protocol.**  Each experiment (198 of them under `experiments/<category>/`: the 186 of the publication programme, and the 12 of the `alt_*` exploration programmes):
1. constructs the relevant object from the **declared** inputs (no fit, unless the experiment names a fitted coefficient);
2. computes the **residual** against the stated value;
3. writes a **JSON verdict** `{id, tier, property, pass, residual,
   status, data, wall_time}` to `experiments/<category>/<id>/results/` — `tier` ($0$–$6$)
   is written by the tiered checks (the root, the kinematics, the dynamics, the
   cosmology and the falsifiability tiers); the residual-attack checks record
   their verdict against a named open problem and omit it;
4. **passes** iff the residual is below the stated tolerance
   (`tol`; default $\varepsilon = 10^{-12}$).
5. a **crash** of the experiment script (non-zero exit code) is recorded as
   a **FAIL**, any JSON verdict already on disk is **discarded** (an experiment
   never inherits the verdict of a previous run), and the **full** traceback is
   kept in `results/crashes/<id>.log`.

`run_validation.py` runs all **198** experiments in name order,
reads the JSONs, and writes `results/validation_summary.json` with the
per-experiment verdicts, the **reduced residuals**, the **named residuals**, and
the runner's hard-tier verdict.  `--report-only` skips the runs and rolls up
the verdict JSONs already on disk: a partial rerun never forces the whole
suite to be re-executed just to rewrite the summary.  The driver exports the
backend the card can host; an experiment that does its linear algebra on the
host pins it (`CL10_BACKEND=numpy` before `cl10` is imported) and the runner
lists it **OFF-CORE** - the whole `alt_*` suite and the cone programme's
inverse/rigidity searches do this, so their verdicts do not depend on what
else is using the card.  The runner labels the hard-tier verdict `U10`; the
paper does not use that label:

$$\mathrm{hard} = \bigwedge_{V:\ \mathrm{tier}(V)\in\{0,1,2,3\}} \mathrm{PASS}(V)$$

and the runner prints that denominator next to the verdict, because an
experiment that writes no tier is not folded into it.

**The result (186/186 publication, 198 PASS of 198).**  Every
experiment of the publication programme **passes**; the hard-tier verdict is
positive. The 12 exploration experiments all pass, each recording its own
verdict; the cone programmes (V175–V183, V191–V201) run on the host numpy
core, so their verdicts do not depend on what else holds the GPU.
Two cautions, stated in the paper (Appendix B): most of these checks are
**internal consistency** (they would pass for any theory built on this root),
and a pass is only as good as its criterion.  (For the *meaning* of each pass —
which residuals are reduced and which are named — see `doc/final_publication.md` Sections 17–18.)

**The provenance of the headline numbers** (the paper's Table 5;
`cl10.parameter_audit.provenance_table()` gives the same classification
over 15 rows as 7 algebraic / 1 free angle / 2 conventions / 1 coincidence
/ 4 named residuals):

- **algebraic** (fixed by the root): the $5/3$ and the $5/24$
  $\mathcal{T}$ normalisation, the $16/16\bar$ split and the singlet
  $B\!-\!L$ assignments $(-1/+1)$, the 8 anomaly coefficients, the 32-level
  spectrum and $\det U = 1$, the $16/16\bar$ spectral mirror, and the
  off-diagonal **10-Higgs** mass term $16\,10\,16\bar\to 1$ (V28, V45);
- **dynamical, measured**: the causal radius $1$ (the front speed is
  exactly $1$ lattice unit per step, V36);
- **conventions** (mapping or scale identification): $M_{\rm GUT} = M_P =
  1.22\times10^{19}$ GeV, $\alpha_{10} \equiv g_{10}^2 = 1/25$,
  $g_0 = 1/2^{n_{\rm space}} = 1/8$ (the 3-space, V38), the EH coefficient
  $1/16\pi$ ($G = 1$) and the lattice CFL $dt^2 = 1/4$ (V37);
- **fitted coefficients** (read from the data): $2.08$ (the hierarchy
  $v/M_P = 2.02\times10^{-17}$), $1.18$ (the DM mass),
  $n_{2\text{step}} = 3$ (the seesaw), with the $246$ GeV electroweak scale
  as the measured anchor;
- **inputs of other sectors**: the SM sphaleron factor $28/79$ (an
  **input of the SM sector**, not a root number: see §11.3 and Table 4 of
  `doc/final_publication.md`);
- **computed but unmatched to data**: the UV-finite 2→2 cross-section
  $\sigma_0 = P_2 / M_P^2$, the dispersion (the 32-level spectrum,
  UV-bounded), the unitarity bookkeeping behind the "optical theorem" (the
  inelastic $P_2$), and the QCA/perturbative ratio $R = 4\pi P_2 / g_{10}^4$
  (O9);
- **mixing**: only $\sin^2\theta_{12} = 1/3$ is **derived** (the unique
  heavy eigenvector); $\sin^2\theta_{23}$ and $\sin^2\theta_{13}$ sweep
  $[0, 2/3]$ under the $O(2)$ freedom of the 2-fold light subspace, and the
  tribimaximal values $(1/2, 0)$ are the $(23)$-adapted **choice** within
  that family (V45); the computed CKM is the **identity** (V49);
- **the gravity sector, implemented and cross-checked**: the **general**
  Riemann computed from an *arbitrary* metric (the 2-route analytic/computed
  agreement $6.5\times10^{-7}$), the **black hole** recovered by computation
  (the computed Kretschmann $48M^2/r^6$, the vacuum $G=0$) and the
  **non-linear** (state-dependent) QCA metric evolution (V50); the 2-route
  3+1 constraint 0-row (the **constraints** $R_{00}, R_{0i}$ by the 4-D
  Riemann vs the 3+1 $(\gamma, K)$ decomposition, agreeing to $O(h^2)$), the
  **Z4** constraint-damping, the **non-linearity** (superposition fails) and
  the **UV-finite** property (V51).

$g_{10} = 1/5$ is the theory's **one free parameter**: V52 flags it as such
and V56 shows the data prefer $n = 7$ over the root $n = 5$.  The **algebra**
has no free parameter; that is not the same as the **dynamics** being
parameter-free — until the **closure axiom** (V131) is imposed: requiring the
fundamental step to be an element of the root ($U^k = \pm I$) leaves exactly
$g = \pi/2$ and $g = \pi/4$ (the canonical Appendix-A ordering), so the angle
becomes a discrete index and the implemented $1/5$ is the **fallback** theory.
What the quantization costs is measured: the closed coupling $\alpha = \pi/64$
sits a factor $1.54$ further from the one-loop mean than $1/25$ did (V135), and
the desert meeting at the closed value leaves the desert (V133).

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
python3 run_validation.py            # -> results/validation_summary.json (162 verdicts)
python3 run_validation.py V11_anomaly V28_yukawa   # a subset (the exact directory names)
CL10_BACKEND=numpy python3 run_validation.py   # force the CPU backend (GPU busy)
python3 run_validation.py --report-only   # roll up the verdict JSONs on disk; run nothing
```

`CL10_BACKEND=numpy` overrides the GPU-first backend selection in
`cl10.engine`; the physics is unchanged, the run is slower (the hydrogen
experiment V62 dominates).

**Lint.**  `python3 -m ruff check cl10/ experiments/ run_validation.py`
passes (the compact numerical style allows E702/E402/E741; the real bugs
F4xx/F5xx are enforced).

**Figures.**  The **twenty-five** figures of `doc/final_publication.md` (placed
in the sections that use them) are **computed, not hand-drawn**.  They are generated by
 the `make_figures` **package**:

 ```bash
 python3 -m make_figures            # -> doc/figures/fig1..fig25_*.png (GPU if the pinned card is free, else numpy)
 python3 -m simulation.hydrogen     # -> simulation/hydrogen/video/ (the Fig. 2 stills + the mp4; git-ignored)
 ```

 The registry in `make_figures/__init__.py` maps every figure number to its
 generator (`algebra`, `scales`, `couplings`, `dynamics`,
 `residuals`); the per-figure provenance - which experiment or module
 computes which figure, including **Fig. 7** (V102 doubler cones),
 **Fig. 8** (V115 isotropy search) and **Fig. 9** (V118/V120 tangent
 velocity algebra) - is listed in the paper's Code
availability statement and in `make_figures/README.md`.  Two figures are
 produced by the simulation engines rather than computed inline: **Fig. 1**
 (the Young double-slit, the `simulation/young` 3-D QCA wave) and **Fig. 2**
 (hydrogen, the **V62** production $296^3$ run of the `simulation/hydrogen`
 engine on the GPU, invoked as a subprocess; its meta JSON lands in
 `experiments/phenomenology/V62_hydrogen_production/results/`, and its content is
 checked by **V62** itself).  `--fast` skips these two heavy
 figures; `--only N [N..]` regenerates a subset; `--cpu` forces numpy.

 The measured values in the comparison figures (the coupling and ladder
 panels) are **reference data**: the SM
 gauge couplings and the electroweak / Planck scales from the **Particle
 Data Group** (`pdg.lbl.gov`), the dark-matter abundance and the baryon
 asymmetry from **Planck 2018** ($\Omega_{\rm DM}h^2 = 0.120$), the
 proton-decay limits from **Super-Kamiokande**, the neutrino scales from
 **NuFIT / KATRIN**, and the direct-detection ceiling from **LZ**.  These
 are **reference data**: the coefficients the framework reads from the data
 ($2.08$, $1.18$, $n_{2\text{step}} = 3$) are **named** in Table 5 of the
 publication, not hidden in a fit.

 It uses the same backend selection as Section 3 (CuPy/GPU, degrading to
 NumPy/CPU).

---

## 5. The status of two theorems (the honest wording)

**The Nielsen–Ninomiya status (V63, `cl10.supercell`).**  The lattice
**chirality** is **not** realised by a naive lattice Dirac operator: the 16 Weyl
is a **module** of $\mathrm{Cl}(10)$ (the pseudoscalar splits the 32 into
16/16-bar), and the QCA is a 't Hooft-type **finite-state** theory (one
**exactly-unitary** per-cell step).  The doubling theorem constrains a **local,
Hermitian, translation-invariant** lattice Dirac operator, so it does not by
itself forbid this construction.  **What the theorem leaves open has now been
computed** rather than asserted.

*Protocol.*  The automaton (`cl10.supercell`) is written on the **cell = the
$2\times2\times2$ block**: eight sub-sites $e \in \{0,1\}^3$, each carrying the
32-D module.  One **step** is six stages in the fixed order
`[(intra,0),(inter,1),(intra,2),(inter,0),(intra,1),(inter,2)]`, each stage the
internal layer on every sub-site followed by one bond stage; the staggering of
the one-cell description is the sub-lattice structure of the cell, so the step
is **time-independent** and the good quantum number is the **block crystal
momentum** $q \in [-\pi,\pi]^3$.  On the eight sub-sites the bond stage along
axis $j$ is $S = \mathbb{1}$ except on the pairs $(e, e \oplus 2^j)$ selected by
the parity, where $S_{ee} = S_{oo} = \cos g_{10}$,
$S_{eo} = -i\sin g_{10}\,\phi$, $S_{oe} = -i\sin g_{10}\,\phi^*$ with $\phi = 1$
on the intra-block stage and $\phi = e^{+iq_j}$ ($e^{-iq_j}$) on the inter-block
stage according to which member of the pair crosses the cell face.  The one-step
Bloch matrix is $M(q) = \prod_{s=1}^{6}\big(S_s(q)\otimes U_{16}\big)$,
$U_{16} = B^\dagger U_{32}B$ the internal factor in the $+i$ eigenspace $B$ of
the pseudoscalar, built exactly ($8\times8$ bond $\otimes$ $16\times16$ internal
= $128\times128$; $256\times256$ on the full module).  It is cross-validated
against the **real-space** step on a $2^3$-block torus: maximum deviation
$9.8\times10^{-15}$ (numpy backend) / $7.4\times10^{-15}$ (cupy backend).

*Results (V63, tier 2, named residual).*  (1) 16-sector closure of the
step: $16\to\overline{16}$ leakage $5.7\times10^{-16}$ (numpy) / $7.0\times10^{-16}$ (cupy) — chirality is a
**dynamical** symmetry, not only an on-site one.  (2) Zone centre:
$\min_a|\omega_a(0)| = 0.036499$ per step (16 sector) and $0.036499$ (32
sector) — **gapped**.  (3) The 12 M-point momenta $(\pm\pi/2,\pm\pi/2,0)$ and their
permutations: maximum gap $2.0\times10^{-4}$, linear slopes $|v| = 0.1419$, node
chiralities $\pm1$ with sum exactly $0$.  (4) $g_{10}$ scan:
$\mathrm{gap}(0)/\mathrm{gap}(\text{doubler}) = 0.0365/0.0002$ at $g_{10}=0.2$,
$0.0715/0.0430$ at $0.15$, $2.3548/2.4722$ at $0.1$ — the zone-centre gap is not
protected.  **The doublers are present in the bare step**; the internal Wilson
factor of V65 gaps all twelve (Section 8.4), and the choice of the split node
set is what remains.  This is the lattice-chirality item of O9: the paper
reports the doubling as a **negative** result, and its removal by a mechanism of
the step itself in Section 8.4.

*The closed angle changes the census (V131, V132, V134).*  The closure axiom
(Section 8.5 of the publication) replaces the implemented $g_{10} = 1/5$ by
$g = \pi/4$, and the doubling picture changes with it.  At $\pi/4$ the exact
block step satisfies $M(0)^2 = I$ below $1\times10^{-14}$, and the **exactly
gapless set** of the block zone is $\Gamma$ together with the **six
face-diagonal lines** $\{(\pm t, \pm t, 0)\}$ and its axis permutations — the
twelve M points are only the high-symmetry points lying on them (gaps
$0$–$2.2\times10^{-16}$), X and R stay gapped at $0.5236$, and a $13^3$ grid
counts 73 gapless points against 0 at $1/5$ (V132).  The touching is
**quadratic**: through a node the gap is zero along the line, linear in one
transverse direction (slope $0.577$) and quadratic in the other, the three axis
exponents being a permutation of $(1,1,2)$ at all twelve nodes.  The
**chirality protocol above therefore does not apply**: it needs isolated zeros,
and at 12/12 nodes the gap along the line through the node stays below
$10^{-9}$ at steps $10^{-3},10^{-2},10^{-1}$, so the $\pm1$ charges whose sum
was exactly $0$ become sign noise — the naive sum is not quantized ($-2$ on the numpy backend, $-6$ on cupy).  The Nielsen–Ninomiya count
at the closed angle must be redone **per line**; the per-line programme is reported below (V219, V220), and it dissolves the count rather than completing it.  At the other
closed angle $\pi/2$ the bond layer is diagonal to $1.2\times10^{-16}$ with the
identity permutation and the eight corner shifts $\{\pm1\}^3$ — a pure lattice
translation, $M(0) = I$ on all 128 internal states, bands exactly linear with
slope $\pm1$, 513 gapless grid points: no mixing, no mass generation (V134).

*The per-line programme, run (V219, V220).*  The count is taken per
line on the icosahedral-module step with the cubic module's $\Gamma$–$M$ line as
control, and the first finding is that **there is nothing band-resolved to
count**: at the line the touching is a cluster of **128 of 1024** eigenphases at
phase zero (cubic $M$: 64), the largest phase inside the cluster at
$2.3\times10^{-14}$ with the first one outside at $3.2\times10^{-1}$, and the
cluster stays isolated along circles of every radius (smallest relative gap 5.2
at $r = 3\times10^{-2}$, rising to 1479 at $r = 3\times10^{-3}$).  Consecutive
eigenvectors of the single closest phase overlap at $10^{-14}$ in the worst
case, so a per-band Berry phase there is a choice of basis inside the
degeneracy: the $\pm1$ charges of the protocol do not exist at these lines, and
their sum can be neither quantized nor sign noise.  (Inside a 128-fold
degenerate cluster the eigensolver's basis is non-orthogonal, so the frame has
to be re-orthonormalised before transport — `V220` records its condition.)  What
is gauge invariant is the cluster's own parallel transport.  The transverse
closing exponent measures $1.006$ for a circle cut across the line and $1.008$
for one tilted $45^\circ$ about it (cubic line: $1.000$): the touching is
linear.  The symmetry test then decides stability: the commutant of the internal
map is computed (obstruction $4.2\times10^{-15}$ over 240 powers of the map), and
the unitary perturbation it supplies opens a gap at the line that grows linearly
with the coupling ($1.277\times10^{-4}$ at $\theta = 10^{-3}$,
$1.277\times10^{-3}$ at $\theta = 10^{-2}$, ratio $10.00$) — the lines are
**not** protected against a symmetry-respecting mass.  The symmetry available to
protect them is small: of the icosahedral point group's 60 elements exactly 2
preserve the line direction, and those leave the spectrum at the line fixed to
exactly $0$ (V220).

*The bond-layer programme, run (V225–V228).*  Proposition 6 said a
parity-degenerate cone cannot live in the bond class: the committed step
factorises, transport against internal layer, and a factorised step carries no
non-Abelian tangent algebra on a degenerate subspace.  The boundary it leaves
open is the **bond layer** of the shift module, and `cl10/spincell.py` builds
it: the coin stays the twelve directed axes of the icosahedral star, and what a
hop carries becomes $U_d=-\hat v_d\cdot\vec\sigma$ — Hermitian, an involution, a
vector under the binary icosahedral group, with $U_{-d}=-U_d$.  `V225` measures
the construction: covariance exact over the $120$ lifts ($3.5\times10^{-15}$),
the quaternionic assignment equally covariant but squaring to $-I$ with a step
gapped by $1.571$, the node present exactly at internal exponent $0$ or
$4\bmod 8$, the $192$-vector gap cluster splitting as
$\mathrm{sym}^1\oplus\mathrm{sym}^3\oplus\mathrm{sym}^5$ under projectors built
from the group, the tangent algebra Abelian on the whole cluster (one norm to
$3.6\times10^{-15}$, commutators $2\times10^{-14}$) and the Clifford relation
holding on the $32$-dimensional spinor piece ($V_j^2 = 0.201002\,I$ to
$5.8\times10^{-16}$, one speed $0.4483$).  `V227` classifies what is gapless:
three coordinate-axis lines, no isolated point, and no reciprocal vector at all
(integrality defect $0.1459$ over $\lvert n\rvert\le3$, the direction
coefficients living in $\mathbb Q(\sqrt5)$) — so the doubling theorem has no
hypothesis to satisfy — and no Berry charge: no sphere about any of the $12$
resolved locus points keeps its link rank.  `V226` settles the charge at the
zone-centre touching by *measuring* the piercing instead of assuming it — on the
icosahedral mesh's $92$ directions no sampled point of the three spheres is
gapless at all (their largest gaps $1.8\times10^{-3}$, $6.0\times10^{-3}$,
$2.4\times10^{-2}$), while on the cell's own six signed coordinate axes $4$, $5$
and $6$ of $6$ come back exactly gapless at those three radii (worst line gap
$2.6\times10^{-8}$): the mesh is simply not sampling the lines the touching sits
on, which is why one sphere of nine keeps rank and the readings run over
$5,12,15,16,25,34$ — and then the protection question: the
$2I$-invariant class is non-empty (obstruction $1.4\times10^{-15}$, norm $7.31$)
and its on-site switch-on opens the touching linearly, gap over coupling
constant to $0.3\%$.  `V228` measures the price: the internal sector is
untouched (chirality conserved to $8\times10^{-31}$, the $5/3$ trace ratio and
the anomaly coefficients unchanged, matter and antimatter blocks isospectral to
$3.8\times10^{-15}$), but no CPT **operator** identity survives on the cell —
the best of six candidates sits at $1.1$, and the one that holds on the
scalar-bond control to $5.0\times10^{-16}$ fails by exactly $\sqrt2$, the same
number as the cell's inversion defect — and the cell's loop return splits by
parity, every odd length vanishing in modulus below $9\times10^{-18}$ while the
even ones are not converged at the grids run ($2109$ points, changes up to
$2.3\times10^{-1}$).  The cell buys a linear node with an isotropic velocity and
a genuine spinor piece; it pays with the inversion, with the CPT operator, and
with a loop factor that has to be recomputed from scratch.

**The CPT status.**  The grading identity is exact (V24): the **45** GUT
generators are **even** in the $\mathbb{Z}_2$ grading (they **commute** with the
pseudoscalar), while the **time** direction is **odd** (it **anticommutes**).
Its **spectral** form is **computed** in V52
(`cl10.parameter_audit.spectral_mirror`): because the dynamics commute with the
pseudoscalar ($[U,w]=0$, the chiral split is preserved), the **antimatter**
($16$-bar) QCA spectrum is the **mirror** (phase-conjugate) of the **matter**
($16$) one — $\omega_{16\bar} = -\omega_{16}$, i.e.
$\mathrm{spec}(U|_{16\bar}) = \mathrm{spec}((U|_{16})^{\dagger})$ — and this
holds for **every** 't Hooft angle (a property of the algebra, not of the
chosen $g_{10}$).  This is a **spectral-mirror** statement.  It is **not** the
CPT theorem of local QFT, which ties spin to statistics through locality, the
spectrum condition and the vacuum's cyclicity; none of that is derived here,
and the paper does not claim it.
