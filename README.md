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

**Nothing is added to the internal space.** Every structure (spacetime, matter,
gauge, statistics) is a property of the **one** object $\mathrm{Cl}(10)$. The
**algebra** has no free parameter; the **dynamics** has one free 't Hooft angle
$g_{10}$ (the proposed reading $1/n_{\rm qubits} = 1/5$ is a convention), and
three coefficients are read from the data — see the parameter audit.

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
│   ├── publication.md        # the paper (the theory + the seven criteria + the 9 open residuals W1-W9 + the 7 falsifiers)
│   ├── reference.md          # the technical reference (the exact constructions, the numerical protocol, reproducibility)
│   └── api.md                # the complete CL10 API reference
├── publication/              # the journal build: doc/publication.md -> REVTeX 4.2 PDF
│   ├── build.sh              # one command: preprocess -> pandoc -> pdflatex x3 -> publication.pdf
│   ├── preprocess.py         # manuscript -> pandoc-friendly markdown (Divs for figures/tables/refs)
│   ├── filter.lua            # Divs -> figure / figure* / table* / thebibliography
│   ├── template.tex          # revtex4-2, aps/prd, twocolumn
│   └── publication.pdf       # the compiled paper (24 pages, 26 figures, 8 tables, 52 refs)
├── cl10/                     # the unique engine (GPU-first, CPU-fallback)
│   ├── engine.py             # the single computation core (CuPy/JAX/NumPy)
│   ├── clifford.py           # the Cl(10) root (the 10 generators, the 1024)
│   ├── spacetime.py          # the (1,3) inside the GUT (time+space GUT dirs)
│   ├── matter.py             # the 16 Weyl + 16-bar + the anomaly cancellation
│   ├── gauge.py              # the 45 + the 5/3 unification + gauge-gravity
│   ├── qca.py                # the unitary QCA + the metric (gravity) dynamics
│   ├── u10pred.py            # the headline numbers (proton, seesaw, families)
│   ├── walls.py              # the 3 reduced walls (W1/W2/W3) + the baryon asymmetry
│   ├── sm_qca.py             # the full SM QCA (the 45 plane rotations at g_10)
│   ├── local_qca.py          # the LOCAL QCA (V36): the causal radius implemented
│   ├── smatrix.py            # the S-matrix of the QCA (V41): the first processes
│   ├── qca_spectrum.py       # the 32-level QCA spectrum: the one object + its invariants (V43)
│   ├── unification.py        # the 5/3, the coupling readings from the root (V38, conventions), the running
│   ├── electroweak.py        # the electroweak scale from the root
│   ├── family.py             # the 3-family construct + the 3-generation QCA step (V39)
│   ├── yukawa.py             # the democratic rank-1 3x3 ({3,0,0}) + the tribimaximal PMNS
│   ├── qcd.py                # the QCD area law + the linear potential (sigma = 1 in Planck units)
│   ├── holo.py               # the holographic area law (the entropy scales with area)
│   ├── baryogenesis.py       # the Sakharov conditions + GUT leptogenesis (V40)
│   ├── cosmology.py          # the cosmology (the W1 cancellation, the hierarchy)
│   ├── cosmology3.py         # the unified cosmology (dark matter + dark energy + antimatter) + the DM scale (V42)
│   ├── einstein.py           # gravity: the Riemann curvature, the EH action derived from the QCA dynamics (V37)
│   ├── nonlinear_gravity.py  # the 4-D non-linear gravity (analytic de Sitter): the full Riemann + Einstein tensor (V46) — the analytic reference
│   ├── gravity.py            # the 4-D non-linear gravity (V50): the general Riemann/Ricci/Einstein/Kretschmann COMPUTED from any metric (2-route verified), the black hole recovered by computation, the discrete EH action (G = 1, dynamics <- action), and the non-linear (state-dependent) QCA metric evolution
│   ├── numrel.py           # the 3+1 BSSN/Z4 numerical relativity (V51, W9): the 2-route 0-row (R_00, R_0i) by 4-D Riemann vs 3+1 gamma/K, the constraints (H, M), Z4 damping, non-linear + UV-finite
│   ├── cpt.py                # the CPT theorem (the 16/16-bar are CPT conjugates)
│   ├── lattice_rg.py         # the QCA lattice RG: the S-matrix to the cross-sections (V44, W9)
│   ├── yukawa_3x3.py         # the 3-gen Yukawa mass matrix: the democratic base (V45, W6) + the quark CKM, honestly (V49)
│   ├── parameter_audit.py    # the parameter audit (V52): the 2-route root (kingdon symbolic vs Pauli), the 16/16-bar spectral mirror, the CP-odd theorem, g10 as the one free 't Hooft angle, the provenance table
│   ├── hier_exponent.py      # the hierarchy exponent -17 = the 10-loop factor (10 = n_gen) (V53, W8)
│   ├── dm_wimp.py            # the sub-Planck DM scale = the 10-loop WIMP scale, M_P*(1/16pi)^10 = 118 GeV (V54, W7)
│   ├── eta_cp.py             # the CP-violating factor = the 2-loop factor (1/16pi)^2 (V55, W5)
│   ├── g10_selection.py      # the 1-loop data check of g_10 = 1/5 (V56, W4: the data prefer n = 7)
│   ├── neutrino.py           # the 2-step seesaw on the 10-loop ladder: m_D=M_P*(1/16pi)^10, M_N=M_P*(1/16pi)^3, m_nu=M_P*(1/16pi)^17 (V57, W2)
│   ├── uvrunning.py          # the UV-finite discrete RG: the coarse-graining tower U^(2^p) is bounded (no Landau pole) (V58, W4)
│   ├── four_limits.py        # the four limits and their common object (V59, paper Section 14): 10 = 2 x 5 (2 = 16/16-bar chiral split, 5 = qubits) + one CP phase — bookkeeping that groups the residuals, not a derivation
│   └── brillouin.py          # the 16-sector Brillouin-zone Floquet spectrum: the 2x2x2 supercell momenta, the doubler nodes, the node chiralities, the g10 scan (V63, W9)
├── experiments/              # the verification program (64 experiments, V0-V63)
│   ├── V0..V10/              # T0-T1: the root, the kinematics, the dynamics (V8-V10: the cosmology)
│   ├── V11..V13/             # T3: anomaly-free, the 5/3, matter+gauge+gravity in one object
│   ├── V14..V19/             # T5: the falsifiable predictions + the 3 REDUCED walls (W1/W2/W3)
│   ├── V20..V24/             # the deepening: the 10-Higgs, the GUT running, the Sakharov, the full SM QCA, the CPT
│   ├── V25..V29/             # holography, the family construct, QCD, the Yukawa, the tribimaximal PMNS
│   ├── V30..V35/             # the GUT completion: the running at M_P, the QCA->GR, the 't Hooft program, the baryon, the cosmology capstones
│   ├── V36..V42/             # the construction modules: the local QCA, the EH action from the dynamics,
│                             #   the coupling readings, the 3 families in the dynamics, leptogenesis, the S-matrix, the DM scale
│   ├── V43/                  # the 32-level QCA spectrum as the one object; the nine residuals are its named projections
│   ├── V44/                  # the QCA lattice RG: the S-matrix to the cross-sections (W9, the UV-finite 0-parameter GUT cross-section)
│   ├── V45/                  # the 3-gen Yukawa mass matrix (W6, the democratic base + the TBM mixing)
│   ├── V46/                  # the 4-D non-linear gravity: the analytic de Sitter reference
│   ├── V47/                  # the 2-loop running: the UV-finite QCA vs the UV-divergent perturbative (W4)
│   ├── V48/                  # the 3- and 4-particle amplitudes: the UV-finite property (W9)
│   ├── V49/                  # the quark CKM, honestly: the democratic 2-level base is fixed, the computed CKM is the identity (W6)
│   ├── V50/                  # the 4-D non-linear gravity (W9): the general Riemann COMPUTED from any metric, the 2-route verification, the black hole recovered by computation, the non-linear QCA dynamics
│   ├── V51/                  # the 3+1 BSSN/Z4 numerical relativity (W9): the 2-route 0-row (R_00, R_0i) by 4-D Riemann vs 3+1 gamma/K, the constraints, Z4 damping, non-linear + UV-finite
│   ├── V52/                  # the parameter audit: the 2-route root (kingdon symbolic vs Pauli), the 16/16-bar spectral mirror, the CP-odd theorem, g10 as the free 't Hooft angle, the provenance table
│   ├── V53/                  # the hierarchy exponent -17 = the 10-loop factor (10 = n_gen) (W8)
│   ├── V54/                  # the sub-Planck DM scale = the 10-loop WIMP scale, M_P*(1/16pi)^10 = 118 GeV (W7)
│   ├── V55/                  # the CP-violating factor = the 2-loop factor (1/16pi)^2 (W5)
│   ├── V56/                # the 1-loop data check of g_10 = 1/5 (W4)
│   ├── V57/                # the 2-step seesaw on the 10-loop ladder: m_nu = M_P*(1/16pi)^17 (W2)
│   ├── V58/                # the UV-finite discrete RG: the bounded tower U^(2^p) (W4)
│   ├── V59/                # the four limits and their common object (paper Section 14): 10 = 2 x 5 + one CP phase (bookkeeping)
│   ├── V60/                # the TBM regression test: the TBM DIAGONALIZES the democratic base (the check the old suite lacked), the wrong-sign matrix is caught, sin^2_12 = 1/3 derived, sin^2_23 / sin^2_13 sweep [0, 2/3], the CKM from identical democratic matrices is the identity
│   ├── V61/                # the 3-D hydrogen atom (Fig. 26; fast 120^3 check): the gaps ARE the Rydberg differences (rel. error 0), the E1 selection rules come out of the lattice wavefunctions, the absorption peaks on the gap and falls with the detuning, the excited atom radiates while it decays.  The ATOM is an input (real hydrogen); the PHOTON is CL10
│   ├── V62/                # the PRODUCTION hydrogen experiment on the GPU (296^3, the full photon script): every Fig. 26 claim measured at production resolution, the figure regenerated as a product of the experiment
│   └── V63/                # the Brillouin-zone chirality test (W9): the 16 sector is closed by the FULL local step (leakage 1e-15), the Bloch Floquet matches the 4^3 torus (4e-15), the zone centre is GAPPED (0.0365), the 12 doubler momenta are GAPLESS (|v| = 0.142), the node chiralities sum to 0 — the doubling survives
├── simulation/               # the 3-D visualizations on the CL10 engine
│   ├── young/                # the Young double-slit as a single-electron measurement (Fig. 25)
│   └── hydrogen/             # the 3-D hydrogen atom + the CL10 photon (Fig. 26): the levels are the real Rydberg eigenstates (an INPUT), the photon is the CL10 free massless field; the video/stills are `python -m simulation.hydrogen`, the Fig. 26 panels are the V62 product
├── run_validation.py         # the runner: runs all 64, writes the hard-tier verdict
└── results/
    └── validation_summary.json   # the 64/64 PASS + the hard-tier verdict
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
python3 run_validation.py        # runs all 64 experiments, prints the hard-tier verdict
python3 make_figures.py          # regenerates the 25 matplotlib figures (fig 1-25) -> doc/figures/
python3 -m simulation.hydrogen   # regenerates the hydrogen video + stills (the Fig. 26 scene)
                                 # -> simulation/hydrogen/video/  (git-ignored)
python3 experiments/V62_hydrogen_production/V62_hydrogen_production.py
                                 # the production hydrogen experiment (296^3, GPU):
                                 # regenerates doc/figures/fig26_hydrogen.png + the meta JSON
bash publication/build.sh        # rebuild the journal PDF from doc/publication.md
                                 # needs pandoc >= 3 and texlive-publishers (revtex4-2)
                                 # -> publication/publication.pdf
```

## The verification program

`doc/publication.md` states the **seven criteria** a unification claim must
meet (Section 1.2, with the paper's own verdict on each), the **tiered**
validation against them, the **nine open residuals** (W1–W9) and the **seven
falsifiers** (F1–F7).  `run_validation.py` runs all **64** experiments
(V0–V63) across the seven tiers and writes

> **hard-tier verdict** = all of tiers 0–3 (V0–V7, V11–V13) PASS — the single
> object, the kinematics, the dynamics, the **anomaly-free** matter, the
> **$5/3$ normalisation**, and **matter–gauge–gravity in one algebra**.  The
> cosmology (tier 4) and the falsifiability (tier 5) are supporting; the named
> residuals are reported, not hidden.

**What the 64 passes mean.**  Most of them check **internal consistency** —
algebraic identities, unitarity of the step, agreement between two numerical
routes — and would pass for any theory built on this root; they are not tests
against nature.  The six comparisons against data ($v/M_P$, $m_\nu$,
$\Delta B$, $m_{\rm DM}$, $\Omega_{\rm DM}h^2$, $\tau_p$) are the
phenomenological content, and each carries the residual named in the paper.

**The parameter audit** (Section 14, V52) classifies every headline number by
provenance: **7 algebraic** invariants, **1 free** 't Hooft angle, **3
conventions** ($\alpha_{10} \equiv g_{10}^2$, $g_0 = 1/8$, $M_{\rm GUT} = M_P$),
**3 fitted coefficients** ($2.08$, $1.18$, $n_{2\text{step}} = 3$) and **9 open
residuals**.  V59 (`cl10.four_limits`) adds the observation that the *same* loop
integer $n = 10$ organises the hierarchy, the DM mass and the seesaw ladder —
bookkeeping that groups the residuals, not a derivation of any value.

**Where the framework fails, stated plainly.**
- Running $\alpha_{10} = 1/25$ down with the SM $\beta$-functions gives
  $\alpha_1(M_Z) = -1.36$ and a $U(1)$ Landau pole at $280$ GeV (W4, F2).
- $m_\nu = 0.146$ eV is $2.9\times$ the atmospheric scale and above
  $\sum m_\nu < 0.12$ eV (W2, F3).
- The relic abundance is **normalised** to the measured value, not computed (W7).
- At the model's own couplings the heavy-neutrino decay is **in** equilibrium
  ($\Gamma_N/H \simeq 4.6$ at $y = 1/8$), so the third Sakharov condition fails
  (W5); the automaton itself satisfies **none** of the three (a computed
  theorem: every CP-odd invariant vanishes).
- The computed CKM is the **identity**; only $\sin^2\theta_{12} = 1/3$ is
  derived (W6).
- The QCD string tension used is $\sigma = 1$ in Planck units, $10^{38}\times$
  the physical value (W9).
- The 16-sector Brillouin-zone spectrum of the implemented automaton (V63) is
  gapped at $q = 0$ and gapless at **12 doubler momenta**: the fermion doubling
  the on-site construction was supposed to evade is present in the spectrum of
  the dynamics, and the theory has no Wilson / overlap / domain-wall term (W9).

**Gravity.**  V50 computes the general 4-D Riemann / Ricci / Einstein /
Kretschmann from an **arbitrary** metric by finite differences, **2-route
verified** against the analytic de Sitter ($6.5\times10^{-7}$), with the
Schwarzschild Ricci $= 0$ and the Kretschmann $48M^2/r^6$ recovered **by
computation**; V51 does the same for the 3+1 constraint 0-row with Z4 damping.
Newton's constant is an **input** ($G = 1$).  The open part of W9 is the sum
over geometries and the quantum completion — and, since V63, the **lattice
chirality** of the implemented step: the 16-sector Brillouin-zone Floquet
spectrum is gapped at the zone centre ($0.0365$) and gapless at the 12 doubler
momenta $(\pm\pi/2,\pm\pi/2,0)$ with node chiralities summing to $0$, so the
fermion doubling is **present** and nothing in the construction removes it.

**The loop ladder.**  One factor $(1/16\pi)^{10}$ **represents** the hierarchy
$v/M_P = 2.08(1/16\pi)^{10}$ and the DM mass $1.18\,M_P(1/16\pi)^{10} = 118$
GeV.  The integer $10$ is chosen so the two-digit coefficients land on the
data: no loop expansion of the QCA exists (W8).

Each experiment constructs its object from the **declared** inputs, computes the
residual against the stated value, and writes a JSON verdict (`tier`,
`property`, `pass`, `residual`) that `run_validation.py` aggregates into
`results/validation_summary.json`.
