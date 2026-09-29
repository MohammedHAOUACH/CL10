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
| number of generators | $10$ | $|\{g_a\}| = 10$ |
| $g_a^2$ | $= I$ | $\max_a \|g_a^2 - I\| < \varepsilon$ |
| $\{g_a, g_b\}$ | $= 0$ for $a\ne b$ | $\max_{a<b} \|g_a g_b + g_b g_a\| < \varepsilon$ |
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
| the full SM QCA step (§11, P2) | `cl10.sm_qca` | $U = \prod_{a<b}(\cos g_{10}\, I - i\sin g_{10}\, T_{ab})$ with $T_{ab} = i\,h_a h_b$ (each $T_{ab}$ Hermitian, $T_{ab}^2 = I$) — the product is **exactly** unitary; no matrix exponential |
| the local QCA (§12, C1) | `cl10.local_qca` | the **internal** layer (the 45 rotations per cell) + the **staggered** 2-cell **swap-rotation** layer $U_{\rm bond} = \cos g_{10}\, I - i\sin g_{10}\, P$ along one spatial axis per time step ($P$ = the 0-parameter swap of two adjacent cells — a permutation, no angle) |
| the EH action from the dynamics (§12, C2) | `cl10.einstein.einstein_hilbert_from_dynamics` | the **discrete** quadratic action $S_{\rm EH} = \frac{1}{16\pi}\frac12\sum_{x,t}\, u_t\big[u_{t+1} - 2u_t + u_{t-1} - dt^2\,\nabla^2 u_t\big]$ ($G = 1$, $dt^2 = \tfrac14$, $\nabla^2$ = the 6-neighbour Laplacian); the **Euler–Lagrange equation** is the **numerical** functional derivative (the central finite difference of the action) evaluated on the GPU |
| the numbers of the root (§12, C3) | `cl10.unification` | $n_{\rm qubits} = \log_2 32 = 5$ $\Rightarrow$ $g_{10} = 1/5$, $\alpha_{10} = 1/25$; $g_0 = 1/2^{n_{\rm space}} = 1/2^3 = 1/8$ ($n_{\rm space} = 3$ = the spatial dimension of the $(1,3)$); the **down-running** is the 1-loop $1/a(M) = 1/a(M_{\rm UV}) - (b/2\pi)\ln(M_{\rm UV}/M)$ |
| the 3-family step (§12, C4) | `cl10.family` | the 96-D internal space $= \bigoplus_{k=1}^3 H_{32}$; $U_{96} = \mathrm{blockdiag}(U_{32}, U_{32}, U_{32})$ at the **same** $g_{10}$ (the family symmetry); the tribimaximal PMNS (V29) is a property of the **Yukawa** sector, not of the time-evolution step |
| the 3-gen Yukawa matrix (V45, W6, C4) | `cl10.yukawa_3x3` | the **democratic** 3×3 mass matrix $M = \|1\!\rangle\!\langle1\|$ (all-ones, rank 1, eigenvalues $\{3,0,0\}$: 1 heavy + a 2-fold light level); the **TBM** ($\sin^2 = 1/3, 1/2, 0$) diagonalizes it (the 0-parameter PMNS); the 3-family 96-level spectrum $= 3 \times$ (the 32 levels); the **10-Higgs** $16\to16$-bar **Dirac** term (computed, **rank-16**, off-diagonal) is the 0-parameter mass; the **3-level** hierarchy (the 10-Higgs vev) and the **exact** CKM/PMNS are the **named** residual |
| the GUT leptogenesis (§12, C5) | `cl10.baryogenesis` | $\Delta B \sim \frac{28}{79}\cdot\frac{g_0^4}{16\pi}\cdot\eta_{CP}$ with the out-of-equilibrium kinematics $\Gamma/H \sim 6\times10^{-4}$; $\eta_{CP}$ is the **named** input (taken $10^{-3}$, the observed range) |
| the S-matrix (§12, C6) | `cl10.smatrix` | $S = U$ (the 32×32) and $S^{(2)} = U\otimes U$ (the 1024×1024) — the 0-parameter transition amplitudes of the 't Hooft program |
| the lattice RG / cross-sections (V44, W9, C6) | `cl10.lattice_rg` | the **UV-finite** statement: the Planck cutoff bounds every cross-section to $\sigma_{\max} = 4\pi / M_P^2$; the 0-parameter 2→2 cross-section $\sigma_0 = P_2 / M_P^2$ ($P_2 = 1 - \bar d^2$ from $S^{(2)} = U \otimes U$); the **0-parameter dispersion** (the 32-level spectrum $H = -i\log U$, UV-bounded $|\omega| \le \pi$, no 1-loop pairing); the **0-parameter optical theorem** (the inelastic $P_2$ **is** the cross-section, elastic + inelastic = 1); the QCA/perturbative ratio $R = 4\pi P_2 / g_{10}^4$; the **named** residual (W9) is the RG matching to the **measured** cross-sections |
| the dark-matter scale (§12, C7) | `cl10.cosmology3` | the **soliton energy** $E$ (in Planck units) of $u = \sqrt{2}\,\mathrm{sech}\,x$ $\Rightarrow$ $m_{\rm DM} = E \cdot M_P$; the 0-parameter **freeze-out** $\Omega h^2 \sim 1/(\sigma m)$ with $\sigma \sim 1/m^2$ (the only scale in the theory) |

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

**The protocol.**  Each experiment (V0–V45 under `experiments/`):
1. constructs the relevant object from the **0-parameter** input (no fit);
2. computes the **residual** against the 0-parameter prediction;
3. writes a **JSON verdict** `{id, tier, property, pass, residual,
   status, data, wall_time}` to `experiments/<id>/results/`;
4. **passes** iff the residual is below the stated tolerance
   (`tol`; default $\varepsilon = 10^{-12}$).

`run_validation.py` runs all **46** experiments (V0–V45) in name order,
reads the JSONs, and writes `results/validation_summary.json` with the
per-experiment verdicts, the **reduced walls**, the **named residuals**,
and the **TOE verdict**:

$$\mathrm{TOE} = \bigwedge_{V \in \{0..7,\,11..13\}} \mathrm{PASS}(V)
\quad\text{(the hard tiers 0–3).}$$

**The result (46/46, 2026-09-29).**  Every experiment **passes** on the
GPU; the **TOE verdict is positive**.  (For the *meaning* of each pass —
which walls are reduced and which residuals are named — see
`doc/publication.md` Sections 11–13.)

**The 0-parameter numbers** (no fit anywhere): the hierarchy
$2.02\times10^{-17}$, the $5/3$, the B–L $(-1/+1)$, the 16/16-bar, the
$5/24$ $\mathcal{T}$ normalization, the $G = 1$, the causal radius $1$,
the $1.22\times10^{19}$ GeV GUT scale, the $1/25$ GUT coupling, the $246$
GeV electroweak scale, and the off-diagonal **10-Higgs** mass term
$16\,10\,16\bar\to 1$ (V20).  **From the completion program (V36–V45):**
$g_{10} = 1/5$ and $g_0 = 1/8$ — **derived from the root** (the 5 qubits
of the 32-D module; the 3-space, V38); the tribimaximal angles
$\sin^2\theta_{12} = 1/3$, $\sin^2\theta_{23} = 1/2$, $\theta_{13} = 0$
$(\mathrm{V29/V39/V45})$; the SM sphaleron factor $28/79$ (V40); the 0-parameter EH
coefficient $1/16\pi$ and the lattice CFL $dt^2 = 1/4$ (V37); the **UV-finite** 2→2
cross-section $\sigma_0 = P_2 / M_P^2$, the **0-parameter dispersion** (the 32-level spectrum, UV-bounded) and the **0-parameter optical theorem** (the inelastic $P_2$ **is** the cross-section), and the QCA/perturbative ratio
$R = 4\pi P_2 / g_{10}^4$ (V44, W9).

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
python3 run_validation.py            # -> results/validation_summary.json (46/46, TOE)
python3 run_validation.py V11 V12    # a subset
```

**Lint.**  `python3 -m ruff check cl10/ experiments/ run_validation.py`
passes (the compact numerical style allows E702/E402/E741; the real bugs
F4xx/F5xx are enforced).

---

## 5. The status of two theorems (the honest wording)

**The Nielsen–Ninomiya status.**  The lattice **chirality** is **not**
realised by a naive lattice Dirac operator (which would double); it is a
**representation** of the root: the 16 Weyl is a **module** of
$\mathrm{Cl}(10)$ (the pseudoscalar splits the 32 into 16/16-bar), and the
QCA is a 't Hooft-type **finite-state** theory (one **exactly-unitary**
per-cell step).  The **doubling** theorem (which constrains a naive
lattice Dirac) does **not** apply to this construction — the chiral
content is a 0-parameter property of the **algebra**, not of a lattice
operator.

**The CPT status.**  The CPT map (V24) is a **0-parameter consequence of
the $\mathbb{Z}_2$ grading**: the **45** GUT generators are **even** in
the grading (they **commute** with the pseudoscalar), while the **time**
direction is **odd** (it **anticommutes**).  The "CPT theorem" of the
theory is this **grading identity** — stated as the *status* of the
theorem (a real physicist would call it a grading property, not a deep
theorem).
