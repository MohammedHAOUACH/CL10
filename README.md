# CL10 — Causal Clifford Theory on the GUT root

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
gravity) is a property of the **one** object $\mathrm{Cl}(10)$; the theory has
**0 free parameters**.

## The single object

The root $\mathrm{Cl}(10)$ has **10** mutually anticommuting generators
$h_1,\dots,h_{10}$ with $h_a^2 = I$:

| structure | source in the $\mathrm{Cl}(10)$ |
|---|---|
| the **spacetime** $(1,3)$ | the **4-generator sub-algebra**: time $=h_1$, space $=i\,h_{2,3,4}$ (the light cone is a GUT sub-algebra) |
| the **matter** (1 SM generation) | the **16-D Weyl** spinor ($16 = 10 \oplus 5\bar \oplus 1 \oplus \nu^c$) and its $\overline{16}$ |
| the **gauge** (the GUT) | the **45** (the $\mathrm{Spin}(10)$ adjoint, the $h_a h_b$ for $a<b$) |
| the **statistics** (boson/fermion) | the **$\mathbb{Z}_2$ grading** (even/odd in the 10 generators) |
| the **gravity** | the **metric** $g^{\mu\nu}$ of the $(1,3)$ sub-algebra + its dynamics |
| the **cosmology** | dark matter (the stable **soliton**), dark energy (the **vacuum**), antimatter (the **16/16-bar** pair) |

## The layout

```
CL10/
├── pyproject.toml            # deps: kingdon, einops, jax[cuda13], cupy-cuda13x
├── doc/
│   ├── publication.md        # the publication (13 sections + abstract)
│   ├── toe_plan.md           # the TOE validation plan (the 7 criteria P0-P6)
│   ├── reference.md          # the professional technical reference
│   └── api.md                # the complete CL10 API reference
├── cl10/                     # the unique engine (GPU-first, CPU-fallback)
│   ├── __init__.py
│   ├── engine.py             # the single computation core (CuPy/JAX/NumPy)
│   ├── clifford.py           # the Cl(10) root (the 10 generators, the 1024)
│   ├── spacetime.py          # the (1,3) inside the GUT (time+space GUT dirs)
│   ├── matter.py             # the 16 Weyl + 16-bar + the anomaly cancellation
│   ├── gauge.py              # the 45 + the 5/3 unification + the GUT coupling/running + gauge-gravity
│   ├── qca.py                # the unitary QCA + the metric (gravity) dynamics
│   ├── toepred.py            # the TOE predictions (proton, seesaw, families)
│   ├── walls.py              # the 3 reduced walls (W1/W2/W3) + the baryon asymmetry
│   ├── sm_qca.py             # the full SM QCA (the 0-parameter unitary QCA, GPU)
│   └── cpt.py                # the CPT theorem (the 16/16-bar are CPT conjugates)
├── experiments/              # the TOE verification program (29 experiments)
│   ├── V0_root/              # T0: the Cl(10) root
│   ├── V1..V4/               # T1: spacetime, matter, gauge, statistics
│   ├── V5..V7/               # T2: hierarchy, gravity, QCA
│   ├── V11_anomaly/          # T3: the 16 Weyl is anomaly-free (the hard one)
│   ├── V12_unification/      # T3: the 5/3 coupling unification
│   ├── V13_gauge_gravity/    # T3: matter+gauge+gravity in one object
│   ├── V8..V10/              # T4: dark matter, dark energy, antimatter
│   ├── V14..V16/             # T5: proton, seesaw, the 3 families
│   ├── V17..V19/             # the 3 REDUCED walls (W1/W2/W3 mechanisms)
│   ├── V20_massgen/          # T2: the 0-parameter 10-Higgs mass term
│   ├── V21_running/          # T3: the 0-parameter GUT coupling + the 1-loop running
│   ├── V22_baryon/           # T5: the baryon asymmetry (the Sakharov conditions)
│   ├── V23_full_sm_qca/      # T2: the full SM as a 0-parameter unitary QCA (the deepening)
│   ├── V24_cpt/              # T3: the CPT theorem (the 16/16-bar are CPT conjugates)
│   ├── V25_holography/       # T4: the holographic area law (the entropy scales with AREA)
│   ├── V26_family/           # T1: the 3-family construct (the 3 generations from the (1,3))
│   ├── V27_qcd_confinement/  # T2: QCD confinement (the 0-parameter area law + linear potential)
│   ├── V28_yukawa/           # T2: the Yukawa / mass-hierarchy 0-parameter structure
├── run_validation.py         # the runner: runs all 29, writes the TOE verdict
└── results/
    └── validation_summary.json   # the 29/29 PASS + the TOE verdict
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
python3 run_validation.py        # runs all 29 experiments, prints the TOE verdict
```

## The TOE validation

`doc/toe_plan.md` defines **what a Theory of Everything must be** (seven
criteria, P0–P6) and the **tiered** validation against it.  `run_validation.py`
runs all **29** experiments (V0–V28) across the six tiers and writes
`results/validation_summary.json`, including the **TOE verdict**:

> **TOE** = all of the hard tiers (0–3: V0–V7, V11–V13) PASS — the single
> object, the kinematics, the dynamics, the **anomaly-free** matter, the
> **$5/3$ coupling unification**, and the **matter-gauge-gravity in one
> object**.  The cosmology (tier 4) and the falsifiability (tier 5) are
> supporting, with the **three reduced walls** (the cosmological constant, the
> seesaw scale, the three families) reported, not hidden.

Each experiment is a **0-parameter** check: it **confirms** a prediction of
the publication (or reports an honest **residual**) and writes a JSON verdict
(`tier`, `property`, `pass`, `residual`) that `run_validation.py` aggregates
into `results/validation_summary.json`.
