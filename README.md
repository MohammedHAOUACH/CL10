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
│   ├── final_publication.md  # the paper (the theory + the seven criteria + the 9 open problems O1-O9 + the 7 falsifiers)
│   ├── reference.md          # the technical reference (the exact constructions, the numerical protocol, reproducibility)
│   ├── api.md                # the complete CL10 API reference
│   ├── figures/              # the 25 publication figures (make_figures/; fig 1 = young QCA, fig 2 = V62 hydrogen product)
├── publication/              # the journal build: doc/final_publication.md -> REVTeX 4.2 PDF
│   ├── build.sh              # one command: preprocess -> pandoc -> pdflatex x3 -> publication.pdf
│   ├── preprocess.py         # manuscript -> pandoc-friendly markdown (Divs for figures/tables/refs)
│   ├── filter.lua            # Divs -> figure / figure* / table* / thebibliography
│   ├── template.tex          # revtex4-2, aps/prd, twocolumn
│   ├── fetch_dois.py         # resolves + verifies the 53 refs against Crossref -> doi.json (clickable DOIs in the PDF)
│   ├── README.md             # the build pipeline, files, prerequisites
│   └── publication.pdf       # the compiled paper (25 figures, 8 tables, 53 refs)
├── make_figures/             # the figure generators: `python -m make_figures` regenerates all 25 (see make_figures/README.md)
│   ├── __init__.py           # the registry in paper order + GPU pinning + the V62 subprocess (fig 2)
│   └── algebra/scales/couplings/dynamics/residuals/style .py
├── cl10/                     # the unique engine (GPU-first, CPU-fallback)
│   ├── clifford.py           # the Cl(10) root (the 10 generators, the 1024)
│   ├── spacetime.py          # the (1,3) inside the GUT (time+space GUT dirs)
│   ├── matter.py             # the 16 Weyl + 16-bar + the anomaly cancellation
│   ├── gauge.py              # the 45 + the 5/3 unification + gauge-gravity
│   ├── qca.py                # the unitary QCA + the metric (gravity) dynamics
│   ├── u10pred.py            # the headline numbers (proton, seesaw, families)
│   ├── walls.py              # the 3 reduced walls (W1/W2/W3) + the baryon asymmetry
│   ├── sm_qca.py             # the full SM QCA (the 45 plane rotations at g_10)
│   ├── supercell.py          # the automaton: the cell is the 2x2x2 block (8 sub-sites x 32 = C^256), the one time-independent six-stage step (internal layer + staggered bond stages), the exact Bloch matrix M(q) on the block zone, the causal radius (V36)
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
│   └── paths.py              # the experiment-path resolver: an experiment is addressed by name, never by category
├── experiments/              # the verification program: 105 experiments in 12 categories
│   ├── algebra_root/         # the root, the closure of the coupling, the provenance audit
│   ├── spacetime_matter/     # the (1+3) subalgebra, one anomaly-free generation, the 45, the mirror, the flavour structure
│   ├── dynamics_spectrum/    # unitarity, the causal radius, the n-particle sector, the doubling, the node geometry
│   ├── gauged_layer/         # the adjoint bond layer: Gauss, the centre obstruction, the area law, confinement
│   ├── gravity/              # the discrete EH action, curvature from any metric, 3+1 evolution
│   ├── scales_loopfactor/    # the closed-loop factor, the ladder length, the composition law
│   ├── couplings_unification/# alpha_10 = 1/25, the desert solution, the proton lifetime, the supplier search
│   ├── neutrino/             # the seesaw ladder, the 1/3 candidate, the protection, the Delta L = 2 no-go
│   ├── cosmology_darkmatter/ # the 118 GeV ladder mass, the carrier no-go, the relic abundance, the portal
│   ├── baryogenesis/         # the Sakharov audit, the asymmetry product, the washout, the real-form obstruction
│   ├── family_number/        # why the count 3 is an input, and what its symmetric alternative costs
│   ├── phenomenology/        # the 296^3 hydrogen run (Fig. 2)
│   └── README.md             # every experiment, its verdict, and the paper section it validates
├── archive/                  # superseded material: 35 earlier experiments, the previous manuscript, its PDF
│   └── README.md             # what was superseded, and by what
├── simulation/               # the 3-D visualizations on the CL10 engine
│   ├── young/                # the Young double-slit as a single-electron measurement (Fig. 1)
│   └── hydrogen/             # the 3-D hydrogen atom + the CL10 photon (Fig. 2): the levels are the real Rydberg eigenstates (an INPUT), the photon is the CL10 free massless field; the video/stills are `python -m simulation.hydrogen`, the Fig. 2 panels are the V62 product
├── run_validation.py         # the runner: runs all 105, writes the hard-tier verdict (`./run_gpu.sh` runs the full suite on the GPU)
├── run_figures.sh            # the 25 figures on the pinned GPU card (aborts if the card is too full for cupy)
└── results/
    └── validation_summary.json   # the 105/105 PASS + the hard-tier verdict

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
python3 run_validation.py        # runs all 105 in parallel (--jobs N), prints the verdict
./run_gpu.sh                     # the full 105 on the GPU (cupy, 8 jobs; JOBS=N to change); numpy fallback per experiment if the card is full
./run_gpu.sh V131_algebraic_closure V132_gapless_census_pi4 V133_desert_pi64 V134_pi2_branch V135_coupling_selection_closed
                                   # the closed-angle batch alone (the runner takes **exact** directory names)
./run_figures.sh                   # the 25 figures on the GPU card make_figures pins; it
                                   # ABORTS if that card has < 2 GiB free, so the figures
                                   # never silently become numpy products (FORCE=1 to run
                                   # anyway, NEED=MiB to move the gate)
python3 -m make_figures            # the same generator without the card check:
                                   # --only N [N..] | --fast (skip 1+2) | --cpu
python3 -m simulation.hydrogen   # regenerates the hydrogen video + stills (the Fig. 2 scene)
                                 # -> simulation/hydrogen/video/  (git-ignored)
bash publication/build.sh        # rebuild the journal PDF from doc/final_publication.md
                                 # needs pandoc >= 3 and texlive-publishers (revtex4-2)
                                 # -> publication/publication.pdf
```

## The verification program

`doc/final_publication.md` states the **seven criteria** a unification claim must
meet (Section 1, with the paper's own verdict on each), the **nine open problems**
(O1–O9) and the **seven falsifiers** (F1–F7).  `run_validation.py` runs the **105**
experiments of `experiments/` (12 categories, catalogued in `experiments/README.md`)
and writes

> **hard-tier verdict** = every tier 0–3 check PASS — the single object, the
> kinematics, the dynamics, the **anomaly-free** matter and the **$5/3$
> normalisation**.  The cosmology and falsifiability tiers are supporting; the
> open problems are reported, not hidden.

**What the 86 passes mean.**  Most of them check **internal consistency** —
algebraic identities, unitarity of the step, agreement between two numerical
routes — and would pass for any theory built on this root; they are not tests
against nature.  The six comparisons against data ($v/M_P$, $m_\nu$,
$\Delta B$, $m_{\rm DM}$, $\Omega_{\rm DM}h^2$, $\tau_p$) are the
phenomenological content, and each carries the residual named in the paper.
**The parameter audit** (Section 14, V52) classifies every headline number by
provenance: **7 algebraic** invariants, **1 free** 't Hooft angle, **2
conventions** ($\alpha_{10} \equiv g_{10}^2$, $g_0 = 1/8$), **6 fitted
coefficients** ($2.08$, $1.18$, $n_{2\text{step}} = 3$, the two V64 desert
thresholds, and the V75 DM portal $\lambda_* = 4\times10^{-4}$) and **9 open residuals**.
V59 (`cl10.four_limits`) adds the observation that the *same* loop
integer $n = 10$ organises the hierarchy, the DM mass and the seesaw ladder —
bookkeeping that groups the residuals, not a derivation of any value.

**Where the framework fails, stated plainly.**
- Running $\alpha_{10} = 1/25$ down with the SM $\beta$-functions gives
  $\alpha_2(M_Z) = 0.195$ ($5.8\times$ the measured value) and an
  $\mathrm{SU}(3)$ Landau pole at $2.2\times10^9$ GeV (W4, F2).  V64 repairs
  this with the theory's own vectorlike remnant matter: three
  $10\oplus\overline{10}$ pairs at $1.2\times10^7$ GeV plus one $(1,2)$ pair
  at $1.0\times10^{12}$ GeV make the measured couplings meet at
  $\alpha_X = 1/25$ **exactly** at $M_X = 1.8\times10^{16}$ GeV — and V65–V68
  connect that content to the automaton (the split doublers of the internal
  Wilson factor carry exactly the vectorlike split-$\mathbf{16}$ fragments;
  the Wilson direction completes the breaking chain; meeting roots need
  $n_{\rm nodes} \ge 3$), **but** the mechanism's own single-angle root meets
  at $M_X = 8.9\times10^{14}$ GeV, excluded by Super-Kamiokande — and V74
  proves no redistribution of the three node strengths can move it (the
  meeting equations depend on them only through their sum, rank-1 Jacobian):
  the pattern is derived, the scale is forced-and-excluded, and the surviving
  solution keeps two fitted thresholds with $M_{\rm GUT} = M_P$ abandoned
  (predicted $\tau_p = 1.8\times10^{36}$ yr, F5).
- $m_\nu = 0.146$ eV is $2.9\times$ the atmospheric scale and above
  $\sum m_\nu < 0.12$ eV (W2, F3).  V76 BUILDS the node-mixing junction the
  $\tfrac13$ candidate needs (cubic-symmetric, closes the node triangle
  exactly; democratic point, eigenvalues $\{3M,0,0\}$).  V103 re-audited its
  protection with the IMPLEMENTED democratic rank-1 $m_D$: the massless pair is
  protected for every junction ratio and the candidate sits in a finite window
  ($r > 0.217$), the atmospheric scale selecting $r = 1.91$.
- The relic abundance is **normalised** to the measured value, not computed,
  and the V66 scan plus the V75 completion find **no viable carrier** at
  118 GeV anywhere in the implemented spectrum — no Wilson node can host a
  stable neutral mirror at the ladder mass (0 survivors over 720 angles);
  the minimal named fix is one mirror-odd ladder scalar, portal
  $\lambda_* = 4\times10^{-4}$, predicted $\sigma_{\rm SI} = 5\times10^{-51}$
  cm$^2$ — nearly decoupled (W7).
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
  the dynamics.  V65 supplies the remedy from the root's own algebra — an
  internal Wilson factor that gaps all 12 nodes — but the choice of the
  desert-scale node set is an input and the scale it predicts is excluded
  (W9/W4).
- The tree-level $\Lambda = 0$ is **protected**: V78 shows $\det U = 1$ pins
  the quasi-energy sum against every traceless internal deformation
  ($|S| < 4\times10^{-14}$ for random $e^{i\varepsilon A}U$, $\varepsilon \le 1$);
  the observed $\Lambda$ is then equivalent to a trace-phase breaking of
  $1.6\times10^{-120}$ — one named number, no mechanism (W1).

**Gravity.**  V50 computes the general 4-D Riemann / Ricci / Einstein /
Kretschmann from an **arbitrary** metric by finite differences, **2-route
verified** against the analytic de Sitter ($6.5\times10^{-7}$), with the
Schwarzschild Ricci $= 0$ and the Kretschmann $48M^2/r^6$ recovered **by
computation**; V51 does the same for the 3+1 constraint 0-row with Z4 damping.
Newton's constant is an **input** ($G = 1$).  The open part of W9 is the sum
over geometries and the quantum completion — and, since V63, the **lattice
chirality** of the implemented step: the 16-sector block-Brillouin one-step
spectrum is gapped at the zone centre ($0.0365$) and gapless at the 12 doubler
momenta $(\pm\pi/2,\pm\pi/2,0)$ with node chiralities summing to $0$; V65's
internal Wilson factor gaps them all, with the split node set an input.

**The loop ladder.**  One factor $(1/16\pi)^{10}$ **represents** the hierarchy
$v/M_P = 2.08(1/16\pi)^{10}$ and the DM mass $1.18\,M_P(1/16\pi)^{10} = 118$
GeV.  The integer $10$ is chosen so the two-digit coefficients land on the
data: no loop expansion of the QCA exists; V69 proves the coarse-graining
tower $U^{2^p}$ is purely kinematic; V72–V73 BUILD the genuine block-spin
RG (exact Schur decimation of the bond layer, any bond gate) and prove it
composes internal angles but can never divide them; and V77 BUILDS the
division mechanism in its smallest form — one quantum link d.o.f. per bond,
traced over, multiplies the level-1 Cartan content by exactly $\cos 2\theta$
(dev $<10^{-15}$).  The value $1/(16\pi)$ needs the dynamical Spin(10)
adjoint link with the Gauss constraint (a many-cell RG with a gauged bond
layer) (W8).

Each experiment constructs its object from the **declared** inputs, computes the
residual against the stated value, and writes a JSON verdict (`tier`,
`property`, `pass`, `residual`) that `run_validation.py` aggregates into
`results/validation_summary.json`.
