# CL10 — Causal Clifford Theory: a quantum cellular automaton on Cl(10,0)

**Mohammed Haouach** — independent researcher  \haouach.mohammed@gmail.com\

[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--2907--0530-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0006-2907-0530)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23321106.svg)](https://doi.org/10.5281/zenodo.23321106) [![Release](https://img.shields.io/github/v/release/MohammedHAOUACH/CL10?label=release)](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0.2) [![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Publication & citation

| | |
|---|---|
| **Concept DOI** — always resolves to the latest Zenodo version | [10.5281/zenodo.23321106](https://doi.org/10.5281/zenodo.23321106) |
| **Version DOI** — this snapshot (v1.0.2, Preprint, CC BY 4.0) | [10.5281/zenodo.23321228](https://doi.org/10.5281/zenodo.23321228) |
| **Zenodo record** — the four paper PDFs + the source archive | [zenodo.org/records/23321228](https://zenodo.org/records/23321228) |
| **GitHub release** — the four PDFs as assets + source archives | [v1.0.2](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0.2) |

```bibtex
@misc{haouach2026cl10,
  author = {Haouach, Mohammed},
  title  = {Causal Clifford Theory (CL10): a quantum cellular automaton on
            Cl(10,0) --- Standard Model kinematics, exactly unitary causal
            dynamics, and a lattice gravitational sector},
  year   = {2026},
  doi    = {10.5281/zenodo.23321106},
  url    = {https://doi.org/10.5281/zenodo.23321106},
  note   = {version v1.0.2: \url{https://doi.org/10.5281/zenodo.23321228}}
}
```

**Licence.** Everything in this repository — code, documentation, papers and
figures — is under the **CC BY 4.0** licence ([LICENSE](LICENSE)), the same
licence as the [Zenodo record](https://zenodo.org/records/23321228).

**The root object.** The theory is a unitary quantum cellular automaton (QCA)
on a Planck lattice whose 2x2x2 cells carry the Clifford algebra
$\mathrm{Cl}(10,0)$ of ten mutually anticommuting generators on eight sub-sites
(on-site space $\mathbb{C}^{256}$), updated by one six-stage product of
$\mathrm{Spin}(10)$ plane rotations and swap-rotations.

What the root gives **with no added field** (demonstrated, Papers I–II): the
$(1{+}3)$ Minkowski subalgebra (time $=h_1$, space $=i\,h_{2,3,4}$ — the light
cone is a GUT sub-algebra); the 16-component Weyl module = one SM generation
with all eight anomalies cancelling; the **45** adjoint containing the twelve
gauge bosons; the $5/3$ hypercharge normalisation; exactly unitary dynamics
($\det U = 1$, $\pm$-paired spectrum, exact matter/antimatter spectral mirror)
with front speed exactly one cell per stage. The **closure axiom** — the update
must be an element of the root — quantizes the rotation angle to
$g \in \{\pi/2, \pi/4\}$: the fundamental step has order 8 and the on-site
update has **zero continuous unfixed numbers**.

What is **not** derived from the root, stated plainly: the gravitational sector
adds the metric $g_{\mu\nu}$ and $G = 1$ by hand; the electroweak/Planck
hierarchy, the DM mass and the neutrino ladder are organised by one loop-factor
**ansatz** $(1/16\pi)^n$ with fitted coefficients; the Born rule, the
spin–statistics identification and the split-node set are declared inputs.
Every headline number carries a provenance class — *algebraic*, *conventional*,
*fitted*, *conditional* or *open* (the audit, Paper III Section 5).

## The papers

The work is published as **three companion papers** split by what is proven,
plus the unified master:

| file | content |
|---|---|
| [`doc/paper1_core.md`](doc/paper1_core.md) · [PDF](doc/paper1_core.pdf) | **Paper I** — the root, the (1+3) subalgebra, the 16, the anomalies, the 45, the unitary causal step, the doubling and its Wilson gap, the provenance discipline |
| [`doc/paper2_gauge_gravity.md`](doc/paper2_gauge_gravity.md) · [PDF](doc/paper2_gauge_gravity.pdf) | **Paper II** — the gauged bond layer (Gauss, area law, confinement), the gravitational sector (EH from the dynamics, 3+1 numerics), what each sector does *not* claim |
| [`doc/paper3_phenomenology.md`](doc/paper3_phenomenology.md) · [PDF](doc/paper3_phenomenology.pdf) | **Paper III** — the loop ladder, the coupling unification, the neutrino/DM/baryogenesis audits, the parameter audit, the failures, the falsifiers |
| [`doc/final_publication.md`](doc/final_publication.md) · [PDF](doc/final_publication.pdf) | the unified master (all of the above, one reference list) |

`doc/*.pdf` are the compiled versions (`review/tools/build_pdfs.sh`); the same four PDFs are attached to the [v1.0.2 release](https://github.com/MohammedHAOUACH/CL10/releases/tag/v1.0.2) and archived on [Zenodo](https://zenodo.org/records/23321228).

## Selected results

All 25 figures are **computed from the construction** — `python3 -m make_figures`,
no hand-entered data ([make_figures/README.md](make_figures/README.md)); each is
placed in the papers and validated by the named experiment
([experiments/README.md](experiments/README.md)).

| | |
|---|---|
| <img src="doc/figures/fig3_one_object_32spectrum.png" width="440"><br>**Fig. 3 — the one object**: the 32-level on-site QCA spectrum and its invariants ([`V43`](experiments/spacetime_matter/V43_qca_spectrum)). | <img src="doc/figures/fig1_young.png" width="440"><br>**Fig. 1 — the Young double-slit as a single-electron measurement** on the automaton (the Tonomura flow): the wave + collapse run of `simulation/young`. |
| <img src="doc/figures/fig2_hydrogen.png" width="440"><br>**Fig. 2 — hydrogen driven by the CL10 lattice photon**: resonant $1s\to 2p$ absorption, the $296^3$ GPU run ([`V62`](experiments/phenomenology/V62_hydrogen_production), master §6). | <img src="doc/figures/fig17_ladder_of_scales.png" width="440"><br>**Fig. 17 — the ladder of scales**: one loop factor $(1/16\pi)^n$ organises the four scales — the ansatz, its fitted coefficients named (Paper III). |
| <img src="doc/figures/fig20_desert_unification.png" width="440"><br>**Fig. 20 — the desert solution**: three $10\oplus\overline{10}$ pairs + one $(1,2)$ pair meet at $M_X = 1.8\times10^{16}$ GeV — at the cost of two fitted thresholds (Paper III). | <img src="doc/figures/fig8_isotropy_search.png" width="440"><br>**Fig. 8 — the isotropy search, a computed failure**: no isotropy point exists on the full Wilson family ([`V115`](experiments/dynamics_spectrum/V115_w9_cone_isotropy)) — reported as O9/F1, not hidden. |

## The layout

```
CL10/
├── pyproject.toml            # deps: kingdon, einops, jax[cuda13], cupy-cuda13x
├── doc/
│   ├── paper1_core.md        # Paper I (the algebraic core)
│   ├── paper2_gauge_gravity.md # Paper II (gauge + gravity sectors)
│   ├── paper3_phenomenology.md # Paper III (phenomenology + audit)
│   ├── final_publication.md  # the unified master
│   ├── *.pdf                 # the four compiled papers
│   ├── reference.md          # the technical reference (exact constructions, numerical protocol, reproducibility)
│   ├── api.md                # the complete CL10 API reference
│   └── figures/              # the 25 publication figures (make_figures/; fig 1 = young QCA, fig 2 = V62 hydrogen product)
├── review/                   # the referee round (simulated PRD)
│   ├── glmreview.md          # the referee report (14 major concerns, 13 minor, 7 questions)
│   ├── glmreview_verification.md # point-by-point verification of every concern against body + code
│   ├── response_letter.md    # the response to the referee
│   └── tools/
│       ├── build_series.py   # master -> the three companion papers (section split, per-paper reference lists)
│       ├── fit_tables.py     # per-table column widths from content (sqrt(L) balance + longest-word floor)
│       ├── fix_tables.py     # post-pandoc: atomic table cells (no mid-row page splits), repeated headers
│       ├── build_pdfs.sh     # fit_tables -> pandoc -> fix_tables -> xelatex x2 -> doc/*.pdf
│       └── headers/          # per-paper LaTeX headers (revtex-style, twocolumn)
├── publication/              # the REVTeX 4.2 journal build of the master (alt pipeline)
│   ├── build.sh              # preprocess -> pandoc -> pdflatex x3 -> publication.pdf
│   ├── preprocess.py / filter.lua / template.tex / fetch_dois.py
│   └── publication.pdf
├── make_figures/             # the figure generators: `python -m make_figures` regenerates all 25 (see make_figures/README.md)
│   ├── __init__.py           # the registry in paper order + GPU pinning + the V62 subprocess (fig 2)
│   └── algebra/scales/couplings/dynamics/residuals/style .py
├── cl10/                     # the library: 39 modules on one core (GPU-first, CPU-fallback)
│   ├── engine.py             # the single computation core: the backend choice, the array type, the primitives (matmul/einsum/roll, the spectral operators, the integrator)
│   ├── coupling.py           # the one parameter of the step: the angle branches (implemented 1/5, closed pi/4, half pi/2) and the two coupling readings (alpha_10 = g^2, g_0 = 1/8) - all conventions, Section 15
│   ├── parpool.py            # the forked process pool the experiments batch through (CL10_NJOBS)
│   ├── clifford.py           # the Cl(10) root (the 10 generators, the 1024)
│   ├── spacetime.py          # the (1,3) inside the GUT (time+space GUT dirs)
│   ├── matter.py             # the 16 Weyl + 16-bar + the anomaly cancellation
│   ├── gauge.py              # the 45 + the 5/3 unification + gauge-gravity
│   ├── qca.py                # the (1,3) metric field + the on-site QCA factor (the 3+1 evolution: einstein / gravity / numrel)
│   ├── u10pred.py            # the headline numbers (proton, seesaw, families)
│   ├── walls.py              # the three named walls and their current status (O1/O2/O3) + the Sakharov audit
│   ├── sm_qca.py             # the full SM QCA (the 45 plane rotations at g_10)
│   ├── supercell.py          # the automaton: the cell is the 2x2x2 block (8 sub-sites x 32 = C^256), the one time-independent six-stage step (internal layer + staggered bond stages), the exact Bloch matrix M(q) on the block zone, the causal radius (V36)
│   ├── smatrix.py            # the S-matrix of the QCA (V41): the first processes
│   ├── qca_spectrum.py       # the 32-level QCA spectrum: the one object + its invariants (V43)
│   ├── unification.py        # the 5/3, the coupling readings from the root (V38, conventions), the running
│   ├── electroweak.py        # the QCA step -> the SM effective action; sin^2 = 3/8 as the GUT boundary condition (not a prediction)
│   ├── family.py             # the 3-family construct + the 3-generation QCA step (V39)
│   ├── yukawa.py             # the democratic rank-1 3x3 ({3,0,0}) + the tribimaximal PMNS
│   ├── qcd.py                # the QCD area law + the linear potential (sigma = 1 in Planck units)
│   ├── holo.py               # the holographic area law (the entropy scales with area)
│   ├── baryogenesis.py       # the Sakharov audit + GUT leptogenesis (V40)
│   ├── cosmology.py          # the cosmology (the O1 vacuum cancellation, the hierarchy)
│   ├── cosmology3.py         # the unified cosmology (dark matter + dark energy + antimatter) + the DM scale (V42)
│   ├── einstein.py           # gravity: the Riemann curvature, the EH action derived from the QCA dynamics (V37)
│   ├── nonlinear_gravity.py  # the 4-D non-linear gravity (analytic de Sitter): the full Riemann + Einstein tensor (V46) — the analytic reference
│   ├── gravity.py            # the 4-D non-linear gravity (V50): the general Riemann/Ricci/Einstein/Kretschmann COMPUTED from any metric (2-route verified), the black hole recovered by computation, the discrete EH action (G = 1, dynamics <- action), and the non-linear (state-dependent) QCA metric evolution
│   ├── numrel.py             # the 3+1 BSSN/Z4 numerical relativity (V51, O9): the 2-route 0-row, the constraints, Z4 damping
│   ├── cpt.py                # the CPT theorem (the 16/16-bar are CPT conjugates)
│   ├── lattice_rg.py         # the QCA lattice RG: the S-matrix to the cross-sections (V44, O9)
│   ├── yukawa_3x3.py         # the 3-gen Yukawa mass matrix: the democratic base (V45, O6) + the quark CKM, honestly (V49)
│   ├── parameter_audit.py    # the parameter audit (V52): the 2-route root, the 16/16-bar spectral mirror, the CP-odd theorem, the provenance table
│   ├── hier_exponent.py      # the hierarchy exponent -17 = the 10-loop factor (10 = n_gen) (V53, O8)
│   ├── dm_wimp.py            # the sub-Planck DM scale = the 10-loop WIMP scale, M_P*(1/16pi)^10 = 118 GeV; the carrier is open (V54, O7)
│   ├── eta_cp.py             # the CP-violating factor = the 2-loop factor (1/16pi)^2 (V55, O5)
│   ├── g10_selection.py      # the 1-loop data check of the angle reading (V56, O4: the data prefer n = 7, so 1/n_qubits is a convention)
│   ├── neutrino.py           # the 2-step seesaw on the 10-loop ladder; the ladder value is excluded by sum m_nu (V57, O2, F3)
│   ├── uvrunning.py          # the UV-bounded discrete RG: the coarse-graining tower U^(2^p) is bounded (no Landau pole); the coupling values stay open (V58, O4)
│   ├── four_limits.py        # the four limits and their common object (V59): bookkeeping that groups the residuals, not a derivation
│   └── paths.py              # the experiment-path resolver: an experiment is addressed by name, never by category
├── experiments/              # the verification program: 115 experiments in 12 categories
│   ├── algebra_root/         # the root, the closure of the coupling, the provenance audit
│   ├── spacetime_matter/     # the (1+3) subalgebra, one anomaly-free generation, the 45, the mirror, the flavour structure
│   ├── dynamics_spectrum/    # unitarity, the causal radius, the n-particle sector, the doubling, the node geometry, the cone isotropy no-gos (V115/V118/V142), the graviton content
│   ├── gauged_layer/         # the adjoint bond layer: Gauss, the centre obstruction, the area law, confinement, the dressed-gate quintet (V140/V141)
│   ├── gravity/              # the discrete EH action, curvature from any metric, 3+1 evolution
│   ├── scales_loopfactor/    # the closed-loop factor, the ladder length, the composition law
│   ├── couplings_unification/# alpha_10 = 1/25, the desert solution, the proton lifetime, the supplier search
│   ├── neutrino/             # the seesaw ladder, the 1/3 candidate, the protection, the Delta L = 2 no-go
│   ├── cosmology_darkmatter/ # the 118 GeV ladder mass, the carrier no-go, the relic abundance, the portal
│   ├── baryogenesis/         # the Sakharov audit, the asymmetry product, the washout, the real-form obstruction
│   ├── family_number/        # why the count 3 is an input, and what its symmetric alternative costs
│   ├── phenomenology/        # the 296^3 hydrogen run (Fig. 2)
│   └── README.md             # every experiment, its verdict, and the paper section it validates
├── archive/                  # superseded material: 35 earlier experiments, the previous manuscript, its PDF, the earlier review audit trail
│   └── README.md             # what was superseded, and by what
├── simulation/               # the 3-D visualizations on the CL10 engine
│   ├── young/                # the Young double-slit as a single-electron measurement (Fig. 1)
│   ├── hydrogen/             # the 3-D hydrogen atom + the CL10 photon (Fig. 2): the levels are the real Rydberg eigenstates (an INPUT), the photon is the CL10 free massless field; the video/stills are `python -m simulation.hydrogen`, the Fig. 2 panels are the V62 product
│   └── blackhole/            # the black-hole visualization (`python -m simulation.blackhole`)
├── run_validation.py         # the runner: runs all 115, writes the hard-tier verdict (`./run_gpu.sh` runs the full suite on the GPU)
├── run_figures.sh            # the 25 figures on the pinned GPU card (aborts if the card is too full for cupy)
└── results/
    └── validation_summary.json   # the full-suite verdict + the hard-tier verdict
```

## The engine: GPU-first, one core

`cl10/engine.py` is the **single** computation core. It picks the backend at
import time — **CuPy** (GPU) when available, else **JAX** (XLA/GPU, with
float64 enabled), else **NumPy** (CPU) — and exposes one array type `E.Array`
plus the primitives every experiment and every module uses (matmul, `einsum`,
`roll`, the spectral operators, the QCA step, the integrator). All of `cl10/`
and `experiments/` go through this one core, so the whole theory runs on the
**RTX 3090** out of the box and degrades to CPU cleanly. The GPU is used for
the array work only: the control flow, the algebra construction and the verdict
logic stay on the host, and every array leaves the core through `E.asnumpy`.

Environment switches the code actually reads: `CL10_BACKEND=numpy|cpu` (force
the CPU core), `CL10_GPU_MIN_FREE_MIB` (the free-memory gate, default 2048),
`CL10_ANGLE=implemented|closed|half` (the rotation-angle branch, Section 1.4:
`implemented` $g_{10} = 1/5$ is the convention the phenomenology is computed
on, `closed` $g_{10} = \pi/4$ is the algebraic closure root, `half`
$g_{10} = \pi/2$; `cl10/coupling.py` is the single source), `CL10_NJOBS` (the
`parpool` worker count), and `CL10_SMOKE=1` (the shrunken experiment grids).

```python
import cl10
cl10.engine.backend_info()     # {'backend': 'cupy', 'einops': True, 'kingdon': True, 'array': 'ndarray'}
cl10.engine.Array              # the array type (cp.ndarray / jax.Array / np.ndarray)
h = cl10.clifford.generators() # the 10 Cl(10) generators (10 x 32x32)
```

## Run it

```bash
cd CL10
source .venv/bin/activate        # GPU (cupy) if installed, else CPU
python3 run_validation.py        # runs all 115 in parallel (--jobs N), prints the verdict
./run_gpu.sh                     # the full suite on the GPU (cupy, 8 jobs; JOBS=N to change); numpy fallback per experiment if the card is full
./run_gpu.sh V142_w9_frame_averaging V141_regular_link_quintet_coupling
                                 # single experiments by exact directory name
./run_figures.sh                 # the 25 figures on the pinned GPU card; ABORTS if the card
                                 # has < 2 GiB free (FORCE=1 to run anyway, NEED=MiB to move the gate)
python3 -m make_figures          # the same generator without the card check: --only N | --fast | --cpu
python3 -m simulation.hydrogen   # regenerates the hydrogen video + stills (the Fig. 2 scene)
python3 review/tools/build_series.py   # master -> the three companion papers (doc/paper*.md)
./review/tools/build_pdfs.sh     # fit_tables -> pandoc -> fix_tables -> xelatex x2 -> doc/*.pdf
                                 # needs pandoc >= 3 + texlive (xelatex, pandoc filters)
bash publication/build.sh        # alt: the REVTeX 4.2 two-column build of the master -> publication/publication.pdf
```

## The verification program

The master states the **seven criteria** a unification claim must meet (with the
paper's own verdict on each), the **nine open problems** (O1–O9) and the
**seven falsifiers** (F1–F7). `run_validation.py` runs the **115** experiments
of `experiments/` (12 categories, catalogued in [experiments/README.md](experiments/README.md)) and
writes

> **hard-tier verdict** = every tier 0–3 check PASS — the single object, the
> kinematics, the dynamics, the **anomaly-free** matter and the **$5/3$
> normalisation**. The cosmology and falsifiability tiers are supporting; the
> open problems are reported, not hidden.

**What the passes mean.** Most check **internal consistency** — algebraic
identities, unitarity of the step, agreement between two numerical routes — and
would pass for any theory built on this root; they are not tests against
nature. The comparisons against data ($v/M_P$, $m_\nu$, $\Delta B$,
$m_{\rm DM}$, $\Omega_{\rm DM}h^2$, $\tau_p$) are the phenomenological content,
and each carries the residual named in the paper.

**The parameter audit** (Paper III Section 5) classifies every headline number:
**7 algebraic** invariants; **free parameters 1 → 0** (the closure axiom
quantizes the 't Hooft angle to $g \in \{\pi/2, \pi/4\}$ — the physical root
$\theta = 1.731$ is one of two discrete roots fixed by the closure of the
coupling); **2 conventions** ($\alpha_{10} \equiv g_{10}^2$, $g_0 = 1/8$);
**11 fitted coefficients** ($1.18$, $n_{\rm 2step} = 3$, the two desert
thresholds, the DM portal $\lambda_*$, the six flavour-dependent bond-layer
couplings $\kappa_i$, …); **1 conditional** (the hierarchy coefficient built
from measured inputs); **9 open residuals** (O1–O9).

**Where the framework fails, stated plainly.**
- **Lorentz isotropy (O9/F1).** The doubler cones are linear in every sampled
  direction but anisotropic; no isotropy point exists on the full Wilson family
  (V115), the node little groups give three nodal planes, not a cone (V118,
  V139), the entangling bond layers reach at most $0.658$ of the required
  commutator (V120/V121), and **frame averaging over the ten generators is
  exactly blind to the anisotropy**: rebuilt in four random $\mathrm{SO}(10)$
  frames the node gap moves by $1.8\times10^{-15}$, the anisotropy by $0$, the
  tangent commutators stay at the noise floor — all momentum dependence lives
  in the sub-site factor of the Kronecker bond layer (V142).
- **Unification (O4/F5).** In the standard normalisation $\alpha = g^2/4\pi$,
  running $\alpha_{10} = 1/25$ down gives $\alpha_2(M_Z) = 0.0086$ and an
  SU(3) Landau pole; the desert solution (three $10\oplus\overline{10}$ pairs
  at $1.2\times10^7$ GeV + one $(1,2)$ pair at $1.0\times10^{12}$ GeV) meets
  exactly at $M_X = 1.8\times10^{16}$ GeV but keeps two fitted thresholds; the
  mechanism's own single-angle root meets at $8.9\times10^{14}$ GeV, excluded
  by proton decay (V64–V74); the bond-layer fix closes the meeting for all
  $n = 3..12$ at the cost of six named $\kappa_i$ (residual 0.005–0.065); the
  supplier space is exhausted — the $\kappa_i$ come from no on-site algebra
  element (V105), no class invariant (V122), the broken flux (V137, residual
  0.605), and no second-order bond number: the superexchange shift misses at
  0.577 and its combination with flux and Cartan shift at 3.46, $R^2 \le 0.30$
  (V145/V146).
- **Dark matter (O7/F6).** No state of the implemented spectrum is a viable
  carrier at 118 GeV (0 survivors over 720 angles, V66/V75); the relic
  abundance is a **computed** number: at the two supplied couplings
  $\Omega h^2 = 0.081$–$0.156$, not $0.120$ for any coupling the root supplies;
  the nodal-plane density $D = 1485$ enhances the final-state phase space by
  $(100/m)^3$-type factors the standard curve cannot see; and the bond-strength
  portal is closed — read as annihilation strength no bond number of the scalar
  class lands on $0.120$ (V147), and read as Higgs portal every one of them
  predicts $\sigma_{\rm SI} = 2.5\times10^{-45}$–$6.3\times10^{-41}$ cm$^2$,
  above the LZ-class bound (V148).
- **Baryogenesis (O5/F7).** The heavy-neutrino decay is in equilibrium
  ($\Gamma_N/H \simeq 4.6$ at $y = 1/8$); the automaton satisfies **none** of
  the three Sakharov conditions (a computed theorem: every CP-odd invariant
  vanishes); the bond-layer route needs $J = 10^{-4}$ for one extra input.
- **Doubling (O9).** The implemented 16-sector spectrum is gapless at 12 doubler
  momenta; the internal Wilson factor gaps them, but the split-node set is an
  input and the scale it predicts is excluded.
- **$\Lambda$ (O1/F4).** The tree-level $\Lambda = 0$ is protected
  ($\det U = 1$ pins the quasi-energy sum, V78); the observed $\Lambda$ is a
  named trace-phase breaking of $1.6\times10^{-120}$, no mechanism.
- The computed CKM is the identity; only $\sin^2\theta_{12} = 1/3$ is derived
  (O6). The QCD string tension used is $\sigma = 1$ in Planck units (O9).

**Gravity.** V50 computes the general 4-D Riemann/Ricci/Einstein/Kretschmann
from an **arbitrary** metric, 2-route verified against analytic de Sitter
($6.5\times10^{-7}$), Schwarzschild recovered by computation; V51 does the 3+1
constraint 0-row with Z4 damping. $G = 1$ is an **input**; the sector adds
$g_{\mu\nu}$ by hand and reproduces GR to $O(h^2)$, not machine precision. The
massless quintet of the directed-link layer disperses linearly and
leading-order isotropically (V140), the regular link closes the hopping
obstruction (V141) — but the **velocity trilemma** is computed: no
$O_h$-covariant step on the 6-bond module is isotropic at $v = 1$, and V142
closes the frame-averaging escape.

Each experiment constructs its object from the **declared** inputs, computes
the residual against the stated value, and writes a JSON verdict (`tier`,
`property`, `pass`, `residual`) that `run_validation.py` aggregates into
[`results/validation_summary.json`](results/validation_summary.json).

**Reproducibility.** Every number the master quotes is produced by a script in
this repository, and the suite is deterministic: fixed seeds, fixed grids, and
each verdict is recomputed from the construction rather than read from disk. The
numbers were verified **internally** — spot-checks against the text, the second
symbolic construction of Appendix A, the exact $3^n$-word loop enumeration — and
the pipeline has **not** been independently re-executed by a second party;
independent re-execution is invited. The published runs used one RTX 3090 with
CuPy; `CL10_BACKEND=numpy` runs the same suite on the CPU core, and
`CL10_ANGLE=closed|half` re-runs it at the other angle branches.
