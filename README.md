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
│   └── publication.md        # the publication (12 sections + abstract)
├── cl10/                     # the unique engine (GPU-first, CPU-fallback)
│   ├── __init__.py
│   ├── engine.py             # the single computation core (CuPy/JAX/NumPy)
│   ├── clifford.py           # the Cl(10) root (the 10 generators, the 1024 algebra)
│   ├── spacetime.py          # the (1,3) inside the GUT (time+space GUT directions)
│   ├── matter.py             # the 16 Weyl + 16-bar (one SM generation)
│   └── qca.py                # the unitary QCA + the metric (gravity) dynamics
├── experiments/              # the verification program (one experiment per claim)
│   ├── V0_root/              # the Cl(10) 10 generators, the 1024, the 16/16-bar
│   ├── V1_spacetime/         # the (1,3) Minkowski inside the GUT, the light cone
│   ├── V2_matter/            # the 16 Weyl (10+5bar+1+nuc), the 16/16-bar
│   ├── V3_gauge/             # the 45 (Spin(10) adjoint), the SM breaking
│   ├── V4_statistics/        # the Z2 grading (boson/fermion)
│   ├── V5_hierarchy/         # the electroweak / QCD hierarchy
│   ├── V6_gravity/           # the (1,3) metric dynamics (the graviton, the soliton)
│   ├── V7_qca/               # the unitary QCA (the full dynamics)
│   ├── V8_dark_matter/       # the stable soliton (the dark matter candidate)
│   ├── V9_dark_energy/       # the QCA vacuum (the dark energy / cosmological constant)
│   └── V10_antimatter/       # the 16/16-bar pair (the CPT, the asymmetry)
├── run_validation.py         # the runner: runs every experiment, writes the summary
└── results/
    └── validation_summary.json
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
python3 run_validation.py        # runs V0..V10, prints the verdicts, writes results/validation_summary.json
```

Each experiment is a **0-parameter** check: it **confirms** a prediction of the
publication (or reports an honest **residual**) and writes a JSON verdict that
`run_validation.py` aggregates into `results/validation_summary.json`.
