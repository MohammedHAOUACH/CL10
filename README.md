# CL10 — The Standard Model as a Single Clifford Root

**Mohammed Haouach** — independent researcher  \haouach.mohammed@gmail.com\

**One object generates the physics.** The theory is a unitary quantum cellular
automaton (QCA) on a Planck lattice whose **single root object is the Clifford
algebra $\mathrm{Cl}(10)$** of ten mutually anticommuting generators.

The key idea of the theory: **the root is not the spacetime but the GUT.**
The $(1{+}3)$-dimensional spacetime is *inside*
the $\mathrm{Cl}(10)$ — the **time is the first GUT direction** ($h_1$) and the
**three space directions are $i\,h_2, i\,h_3, i\,h_4$** (GUT generators times
$i$). Because all four are GUT generators, they anticommute, so the **Minkowski
$(1,3)$ relations hold 0-parameter** and the **light cone is a sub-algebra of
the GUT**. The matter is the **16-D Weyl spinor** of $\mathrm{Cl}(10)$ (one SM
generation, $16 = 10 \oplus 5\bar \oplus 1 \oplus \nu^c$), the gauge is the
**45** (the $\mathrm{Spin}(10)$ adjoint), the statistics come from the
**$\mathbb{Z}_2$ grading**, and gravity is the **metric dynamics** of the
$(1,3)$ sub-algebra.

**Nothing is added.** Every structure (spacetime, matter, gauge, statistics,
gravity) is a property of the **one** object $\mathrm{Cl}(10)$; the
**algebra** has **0 free parameters** (the 't Hooft angle of the dynamics
carries a 0-parameter selection rule).

## The single object

The root $\mathrm{Cl}(10)$ has **10** mutually anticommuting generators
$h_1,\dots,h_{10}$ with $h_a^2 = I$:

| structure | source in the $\mathrm{Cl}(10)$ |
|---|---|
| the **spacetime** $(1,3)$ | the **4-generator sub-algebra**: time $=h_1$, space $=i\,h_{2,3,4}$ (the light cone is a GUT sub-algebra) |
| the **matter** (1 SM generation) | the **16-D Weyl** spinor ($16 = 10 \oplus 5\bar \oplus 1 \oplus \nu^c$) and its $\overline{16}$ |
| the **gauge** (the GUT) | the **45** (the $\mathrm{Spin}(10)$ adjoint, the $h_a h_b$ for $a<b$) |
| the **statistics** (boson/fermion) | the **$\mathbb{Z}_2$ grading** (even/odd in the 10 generators) |
| the **gravity** | the **metric** $g^{\mu\nu}$ of the $(1,3)$ sub-algebra + its dynamics (the **complete** 4-D non-linear gravity: the general Riemann computed from any metric, the 0-parameter EH action, dynamics $\leftarrow$ action, the non-linear QCA) + the **3+1 BSSN/Z4 numerical relativity** (V51) |
| the **cosmology** | dark matter (the stable **soliton**), dark energy (the **vacuum**), antimatter (the **16/16-bar** pair) |

## The layout

```
CL10/
├── pyproject.toml            # deps: kingdon, einops, jax[cuda13], cupy-cuda13x
├── doc/
│   ├── publication.md        # the publication (the theory + the U10 verification program + the 9 residuals W1-W9 + the falsifiers)
│   ├── reference.md          # the technical reference (the exact constructions, the numerical protocol, reproducibility)
│   └── api.md                # the complete CL10 API reference
├── cl10/                     # the unique engine (GPU-first, CPU-fallback)
│   ├── engine.py             # the single computation core (CuPy/JAX/NumPy)
│   ├── clifford.py           # the Cl(10) root (the 10 generators, the 1024)
│   ├── spacetime.py          # the (1,3) inside the GUT (time+space GUT dirs)
│   ├── matter.py             # the 16 Weyl + 16-bar + the anomaly cancellation
│   ├── gauge.py              # the 45 + the 5/3 unification + gauge-gravity
│   ├── qca.py                # the unitary QCA + the metric (gravity) dynamics
│   ├── u10pred.py            # the U10 predictions (proton, seesaw, families)
│   ├── walls.py              # the 3 reduced walls (W1/W2/W3) + the baryon asymmetry
│   ├── sm_qca.py             # the full SM QCA (the 45 plane rotations at g_10)
│   ├── local_qca.py          # the LOCAL QCA (V36): the causal radius implemented
│   ├── smatrix.py            # the S-matrix of the QCA (V41): the first processes
│   ├── qca_spectrum.py       # the 32-level QCA spectrum: the one 0-parameter object + the 9-to-1 reduction (V43)
│   ├── unification.py        # the 5/3, the 0-parameter couplings derived from the root (V38), the running
│   ├── electroweak.py        # the electroweak scale from the root
│   ├── family.py             # the 3-family construct + the 3-generation QCA step (V39)
│   ├── yukawa.py             # the democratic rank-1 3x3 ({3,0,0}) + the tribimaximal PMNS
│   ├── qcd.py                # the 0-parameter QCD area law + the linear potential
│   ├── holo.py               # the holographic area law (the entropy scales with area)
│   ├── baryogenesis.py       # the Sakharov conditions + GUT leptogenesis (V40)
│   ├── cosmology.py          # the cosmology (the W1 cancellation, the hierarchy)
│   ├── cosmology3.py         # the unified cosmology (dark matter + dark energy + antimatter) + the DM scale (V42)
│   ├── einstein.py           # gravity: the Riemann curvature, the EH action derived from the QCA dynamics (V37)
│   ├── nonlinear_gravity.py  # the 4-D non-linear gravity (analytic de Sitter): the full Riemann + Einstein tensor (V46) — the analytic reference
│   ├── gravity.py            # the COMPLETE 4-D non-linear gravity (V50): the general Riemann/Ricci/Einstein/Kretschmann COMPUTED from any metric (2-route verified), the black hole recovered by computation, the complete 0-parameter EH action (dynamics <- action), and the non-linear (state-dependent) QCA metric evolution
│   ├── numrel.py           # the 3+1 BSSN/Z4 numerical relativity (V51, W9): the 2-route 0-row (R_00, R_0i) by 4-D Riemann vs 3+1 gamma/K, the constraints (H, M), Z4 damping, non-linear + UV-finite
│   ├── cpt.py                # the CPT theorem (the 16/16-bar are CPT conjugates)
│   ├── lattice_rg.py         # the QCA lattice RG: the S-matrix to the cross-sections (V44, W9)
│   ├── yukawa_3x3.py         # the 3-gen Yukawa mass matrix: the 9-to-1 projection (V45, W6) + the quark CKM, honestly (V49)
│   ├── parameter_audit.py    # the honest parameter audit (V52): the solution to the four limits of the theory — the 2-route root (kingdon symbolic vs Pauli), the 16/16-bar spectral mirror, the 0-parameter CP-odd theorem, g10 the one free 't Hooft angle, the provenance table
│   ├── hier_exponent.py      # the hierarchy exponent -17 = the 10-loop factor (10 = n_gen) (V53, W8)
│   ├── dm_wimp.py            # the sub-Planck DM scale = the 10-loop WIMP scale, M_P*(1/16pi)^10 = 118 GeV (V54, W7)
│   ├── eta_cp.py             # the CP-violating factor = the 2-loop factor (1/16pi)^2 (V55, W5)
│   ├── g10_selection.py      # the 1-loop data check of the 0-parameter g_10 = 1/5 (V56, W6)
│   ├── neutrino.py           # the 2-step seesaw on the 10-loop ladder: m_D=M_P*(1/16pi)^10, M_N=M_P*(1/16pi)^3, m_nu=M_P*(1/16pi)^17 (V57, W2)
│   ├── uvrunning.py          # the UV-finite discrete RG: the coarse-graining tower U^(2^p) is bounded (no Landau pole) (V58, W4)
│   └── four_limits.py        # the UNIFIED solution to the four limits (V59, Section 16): 10 = 2 x 5 (2 = 16/16-bar chiral split, 5 = qubits) + one CP phase -> the 4 limits reduce to 6 named (1 't Hooft angle + 2 conventions + 2 2-digit coeffs + 1 CP phase)
├── experiments/              # the U10 verification program (63 experiments, V0-V62)
│   ├── V0..V10/              # T0-T1: the root, the kinematics, the dynamics (V8-V10: the cosmology)
│   ├── V11..V13/             # T3: anomaly-free, the 5/3, matter+gauge+gravity in one object
│   ├── V14..V19/             # T5: the falsifiable predictions + the 3 REDUCED walls (W1/W2/W3)
│   ├── V20..V24/             # the deepening: the 10-Higgs, the GUT running, the Sakharov, the full SM QCA, the CPT
│   ├── V25..V29/             # holography, the family construct, QCD, the Yukawa, the tribimaximal PMNS
│   ├── V30..V35/             # the GUT completion: the running at M_P, the QCA->GR, the 't Hooft program, the baryon, the cosmology capstones
│   ├── V36..V42/             # THE COMPLETION PROGRAM (C1-C7): the local QCA, the EH action from the dynamics,
│                              #   the couplings from the root, the 3 families in the dynamics, leptogenesis, the S-matrix, the DM scale
│   ├── V43/                  # the 9-to-1 reduction: the 32-level QCA spectrum is the one object; the 9 residuals its projections
│   ├── V44/                  # the QCA lattice RG: the S-matrix to the cross-sections (W9, the UV-finite 0-parameter GUT cross-section)
│   ├── V45/                  # the 3-gen Yukawa mass matrix: the 9-to-1 projection (W6, the democratic base + the TBM mixing)
│   ├── V46/                  # the 4-D non-linear gravity: the Einstein equations from the QCA metric (C2, the analytic de Sitter reference)
│   ├── V47/                  # the 2-loop running: the UV-finite QCA vs the UV-divergent perturbative (W4)
│   ├── V48/                  # the 3- and 4-particle amplitudes: the UV-finite property (W9)
│   ├── V49/                  # the quark CKM, honestly: the democratic 2-level base is 0-parameter, the computed CKM is the identity (W6)
│   ├── V50/                  # the COMPLETE 4-D non-linear gravity (W9, C2): the general Riemann COMPUTED from any metric, the 2-route verification, the black hole recovered by computation, the non-linear QCA dynamics
│   ├── V51/                  # the 3+1 BSSN/Z4 numerical relativity (W9, C2): the 2-route 0-row (R_00, R_0i) by 4-D Riemann vs 3+1 gamma/K, the constraints, Z4 damping, non-linear + UV-finite
│   ├── V52/                  # the honest parameter audit: the solution to the four limits of the theory — the 2-route root (kingdon symbolic vs Pauli), the 16/16-bar spectral mirror, the 0-parameter CP-odd theorem, g10 the free 't Hooft angle, the provenance table
│   ├── V53/                  # the hierarchy exponent -17 = the 10-loop factor (10 = n_gen) (W8)
│   ├── V54/                  # the sub-Planck DM scale = the 10-loop WIMP scale, M_P*(1/16pi)^10 = 118 GeV (W7)
│   ├── V55/                  # the CP-violating factor = the 2-loop factor (1/16pi)^2 (W5)
│   ├── V56/                # the 1-loop data check of the 0-parameter g_10 = 1/5 (W6)
│   ├── V57/                # the 2-step seesaw on the 10-loop ladder: m_nu = M_P*(1/16pi)^17 (W2)
│   ├── V58/                # the UV-finite discrete RG: the bounded tower U^(2^p) (W4)
│   ├── V59/                # the UNIFIED solution to the four limits (Section 16): 10 = 2 x 5 (2 = chiral split, 5 = qubits) + one CP phase -> the 4 limits reduce to 6 named
│   ├── V60/                # the TBM regression test: the TBM DIAGONALIZES the democratic base (the check the old suite lacked), the wrong-sign matrix is caught, sin^2_12 = 1/3 derived, sin^2_23 / sin^2_13 sweep [0, 2/3], the CKM from identical democratic matrices is the identity
│   ├── V61/                # the 3-D hydrogen atom (Fig. 26; fast 120^3 check): the gaps ARE the Rydberg differences (rel. error 0), the E1 selection rules come out of the lattice wavefunctions, the absorption peaks on the gap and falls with the detuning, the excited atom radiates while it decays.  The ATOM is an input (real hydrogen); the PHOTON is CL10
│   └── V62/                # the PRODUCTION hydrogen experiment on the GPU (296^3, the full photon script): every Fig. 26 claim measured at production resolution, the figure regenerated as a product of the experiment
├── simulation/               # the 3-D visualizations on the CL10 engine
│   ├── young/                # the Young double-slit as a single-electron measurement (Fig. 25)
│   └── hydrogen/             # the 3-D hydrogen atom + the CL10 photon (Fig. 26): the levels are the real Rydberg eigenstates (an INPUT), the photon is the CL10 free massless field; the video/stills are `python -m simulation.hydrogen`, the Fig. 26 panels are the V62 product
├── run_validation.py         # the runner: runs all 63, writes the U10 verdict
└── results/
    └── validation_summary.json   # the 63/63 PASS + the U10 verdict
```

## The engine: GPU-first, one core

`cl10/engine.py` is the **single** computation core. It picks the backend at
import time — **CuPy** (GPU) when available, else **JAX** (XLA/GPU), else
**NumPy** (CPU) — and exposes one array type `E.Array` plus the primitives every
experiment and every module uses (matmul, the spectral operators, the QCA step,
the integrator). All of `cl10/` and `experiments/` go through this one core, so
the whole theory runs on the **RTX 3090** out of the box and degrades to CPU
cleanly.

```python
import cl10
cl10.engine.BACKEND            # 'cupy' | 'jax' | 'numpy'
cl10.engine.Array              # the array type (cp.ndarray / jax.Array / np.ndarray)
h = cl10.clifford.generators() # the 10 Cl(10) generators (10 x 32x32)
```

## Run it

```bash
cd CL10
source .venv/bin/activate        # GPU (cupy) if installed, else CPU
python3 run_validation.py        # runs all 63 experiments, prints the U10 verdict
python3 make_figures.py          # regenerates the 25 matplotlib figures (fig 1-25) -> doc/figures/
python3 -m simulation.hydrogen   # regenerates the hydrogen video + stills (the Fig. 26 scene)
                                 # -> simulation/hydrogen/video/  (git-ignored)
python3 experiments/V62_hydrogen_production/V62_hydrogen_production.py
                                 # the production hydrogen experiment (296^3, GPU):
                                 # regenerates doc/figures/fig26_hydrogen.png + the meta JSON
```

## The U10 validation

Section 11 of `doc/publication.md` defines **what a 0-parameter unification
must be** (seven criteria, P0–P6), the **tiered** validation against them,
and the **falsifiers** (F1–F7).  `run_validation.py`
runs all **63** experiments (V0–V62) across the seven tiers and writes

> **U10** = all of the hard tiers (0–3: V0–V7, V11–V13) PASS — the single
> object, the kinematics, the dynamics, the **anomaly-free** matter, the
> **$5/3$ coupling unification**, and the **matter-gauge-gravity in one
> object**.  The cosmology (tier 4) and the falsifiability (tier 5) are
> supporting, with the **reduced walls** (the cosmological constant, the
> seesaw scale, the three families, the baryon asymmetry, the dark-matter
> abundance) reported, not hidden.

The **honest parameter audit** (Section 15, V52) turns the four limits
(`doc/remarques.txt`) into four computed theorems + a provenance table, and
**Section 16** (V59, `cl10.four_limits`) goes further and **unifies** them:
the four "limits" are the four projections of one 0-parameter invariant,
$10 = 2 \times 5$ (the chiral $16/16$-bar split times the 5 qubits of the
32-D module), **plus one irreducible CP phase** — and the **algebra** has
**0 free parameters** ($g_{10}=1/5$ and $g_0=1/2^3$ are 0-parameter
**selection rules** of the root; their *reading* is the audit's mapping
convention).  The named residual set is the minimum: **1 't Hooft angle of
the dynamics + 2 two-digit coefficients (2.08, 1.18 — *named*, read from
the data) + 1 CP phase = 3**, plus the mapping conventions (Point 1).

The **completion program** (V36–V42, V46, V50, V51, tier 6 — Section 12 of
`doc/publication.md`)
closes the gap from a consistent 0-parameter framework to a **complete**
U10: the **local** QCA (the causal radius *implemented*, C1), the
**dynamics derived from the 0-parameter Einstein–Hilbert action** (C2), the
the **0-parameter selection rules** of the numbers of the root ($g_{10}=1/5$
from the 5 qubits, $g_0=1/8$ from the 3-space, C3), the **three families inside the
dynamics** (the 96-D 0-parameter unitary step, C4), the **baryon asymmetry
at the right order** (the 0-parameter GUT leptogenesis mechanism, C5), the
**S-matrix of the QCA** (the 't Hooft program, C6), the **quantified
dark-matter scale** (C7) — and the **9-to-1 reduction** (V43): the
$32$-level QCA spectrum is the **one** 0-parameter object, and the nine
residuals (W1–W9) are its **nine named projections** (the 0-point is the
most invariant — the $\mathrm{SL}(32)$ identity, $=0$ for every coupling).

The **complete gravity** (V50) completes the 4-D curvature and the action part of C2: the 4-D non-linear
gravity is **completed** from the analytic de Sitter (V46) to the **general**
computation — the full Riemann / Ricci / Einstein / Kretschmann are **computed
by finite differences from an *arbitrary* metric** (not hard-coded, not
restricted to one solution), **2-route verified** (the computed de Sitter
matches the independent analytic one to $O(h^2)$; the Schwarzschild Ricci is
$0$ and the Kretschmann is $48M^2/r^6$, so the **black hole is recovered by
computation**), the **complete** 0-parameter discrete EH action
($S = \frac{1}{16\pi}\sum \sqrt{-\det g}(R-2\Lambda)$, $G=1$) whose functional
derivative **is** the non-linear Einstein equation (dynamics $\leftarrow$
action), and the **non-linear** (state-dependent) QCA metric evolution (the
't Hooft non-linear completion: unitary/UV-finite, causal, superposition
fails).  The **3+1 BSSN/Z4 numerical relativity** (V51) completes the
remaining classical piece: the 3+1 Ricci **0-row** ($R_{00}$, $R_{0i}$) is
**2-route verified** (the 4-D Riemann, V50, vs the 3+1 $(\gamma, K)$
decomposition, agreeing to $O(h^2)$), the 0-row is the 0-row of the 4-D
Einstein tensor (the **constraints** $H, M$), the evolution is **Z4-damped**
(a constraint violation decays) and **non-linear** (superposition fails),
and **UV-finite** (bounded).  The **open** part of W9 is now precisely the
**quantum** completion — not the classical gravity (V50) nor the classical
3+1 numerical relativity (V51), which are complete here.

The **10-loop / 2-loop refinement** (V53–V58, Section 15) resolves the
remaining "fixed-not-explained" limits.  **One** 0-parameter factor,
$(1/16\pi)^{10}$, gives **two** results: the **hierarchy** $v/E_P =
2.08\,(1/16\pi)^{10}$ (V53) **and** the **sub-Planck dark-matter scale**
$m_{\rm DM} = M_P(1/16\pi)^{10} = 118$ GeV (the WIMP scale, resolving the
V42 Planck-soliton gap; V54).  The **10** is **$2\times5$**: the 16/16-bar
**chiral split** $\times$ the **5 qubits** of the 32-D module ($2^5 = 32$).
The **baryon asymmetry** (V55) is a **product of factors of three provenances**
(root / QFT / SM): $\Delta B = (28/79)\,\varepsilon_{CP}\,D$ (the SM sphaleron
$28/79$ — an **input of the SM sector**, the LH d.o.f. count, a root number only
conditional on V23, the
1-loop CP-odd $\varepsilon_{CP} = g_0^4/16\pi$ with $g_0 = 1/2^3$, and the
2-loop dilution $D = (1/16\pi)^2$), giving a **12%** match to the observed
$\Delta B$.  The **2-digit** coefficients and the **CP phase** remain the
named residuals — the loop counts are numbers of the root, not fits.

Each experiment is a **0-parameter** check: it **confirms** a prediction of
the publication (or reports an honest **residual**) and writes a JSON verdict
(`tier`, `property`, `pass`, `residual`) that `run_validation.py` aggregates
into `results/validation_summary.json`.
