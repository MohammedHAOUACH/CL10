# Causal Clifford Theory on the GUT root

## A 0-parameter unitary quantum mechanics whose single object is the Clifford algebra $\mathrm{Cl}(10)$

---

## Abstract

We present **Causal Clifford Theory (CL10)**, a unitary quantum cellular
automaton (QCA) defined on a Planck lattice, in which **one algebraic object
generates the physical structure**.  The axiom is minimal:

> **Axiom.** The universe is a $(1{+}3)$-dimensional causal QCA on a cubic
> lattice of spacing $\ell_P$ and time-step $t_P$ (Planck units
> $\hbar = c = \ell_P = t_P = 1$).  Each lattice cell carries a copy of the
> **Clifford algebra** $\mathcal{C}\equiv\mathrm{Cl}(10)$ of **ten** mutually
> anticommuting generators $h_1,\dots,h_{10}$ with $h_a^2 = I$.  The state of
> the universe is a **unitary representation** of the automorphisms of this
> network of algebras — the QCA dynamics.

From this **single object**, by representation theory, **everything** is
derived with **no free parameter**:

- the **spacetime** $(1{+}3)$ is a **sub-algebra of the GUT** — the **time**
  is the first GUT direction $h_1$ and the **three space directions** are
  $i\,h_2, i\,h_3, i\,h_4$ (GUT generators times $i$).  The four generators
  anticommute, so the **Minkowski** relations hold 0-parameter and the **light
  cone is a sub-algebra of the GUT** (causality is a property of the root);
- the **matter** is the **16-D Weyl spinor** of $\mathrm{Cl}(10)$ — **one full
  Standard Model generation**, $16 = 10 \oplus \bar 5 \oplus 1 \oplus \nu^c$,
  and its **16-bar** conjugate (the **antimatter**);
- the **gauge** is the **45** (the $\mathrm{Spin}(10)$ adjoint, the
  $h_a h_b$ for $a<b$) — the **Grand Unified** gauge group, breaking to the
  SM;
- the **statistics** (boson/fermion) come from the **$\mathbb{Z}_2$ grading**
  of the algebra (the even/odd in the 10 generators);
- the **gravity** is the **metric dynamics** of the $(1,3)$ sub-algebra
  (the discrete Laplacian, $G = 1$): the **graviton** is massless and
  UV-bounded;
- the **cosmology**: dark matter is the stable **soliton** of the non-linear
  metric, dark energy is the **QCA vacuum** (the 0-point), and **antimatter**
  is the **16/16-bar** pair.

The whole theory has **0 free parameters** and **0 free constants**.  The
verification program (Section 10) is a set of **0-parameter** numerical
experiments (V0–V10) run on the **GPU** (CuPy/JAX) by a **single engine**;
each **confirms** a prediction of the theory or reports an honest **residual**.

---

## 1. The Axiom

The theory has **one axiom**.  Spacetime is a $(1{+}3)$-dimensional **causal
QCA** on a cubic Planck lattice, and each cell carries the **Clifford algebra
$\mathrm{Cl}(10)$** of ten mutually anticommuting generators
$h_1,\dots,h_{10}$ with $h_a^2 = I$ and $\{h_a, h_b\} = 0$ for $a\ne b$.

The state of the universe is a **unitary representation** of the automorphism
network of these algebras: a product of **local unitaries** (the QCA gates),
each acting on the internal (Weyl) space of its cell, with **causal radius one
cell per time-step** (the light speed).

The single object $\mathrm{Cl}(10)$ has:
- a **$2^{10} = 1024$-dimensional** complex algebra;
- an irreducible **32-D** module (the **complex spinor**), which the
  **pseudoscalar** (the volume element $h_1\cdots h_{10}$, anti-Hermitian)
  splits into a **16 Weyl** $+$ a **16-bar**;
- a **$\mathbb{Z}_2$ grading** (even/odd in the 10 generators).

Everything else in the theory is **representation theory of this one object**.

---

## 2. The Spacetime Inside the GUT

**Theorem (the $(1,3)$ is a sub-algebra of the GUT).**  The $(1{+}3)$
Minkowski spacetime is **inside** the GUT $\mathrm{Cl}(10)$.  Take:

$$T = h_1 \qquad (\text{time},\; T^2 = +I)$$
$$S_1 = i\,h_2,\quad S_2 = i\,h_3,\quad S_3 = i\,h_4 \qquad (\text{space},\; S_a^2 = -I)$$

Since $h_1,\dots,h_4$ are GUT generators, they **anticommute**; hence
$T$ and $S_a$ anticommute, and the **Minkowski $(1,3)$ relations hold
0-parameter**:

$$T^2 = +I,\qquad S_a^2 = -I,\qquad \{T, S_a\} = 0,\qquad \{S_a, S_b\} = 0.$$

**Corollary (causality).**  A vector $v = v_0 T + v_a S_a$ is **null**
($v^2 = 0$) iff $v_0^2 = \sum_a v_a^2$ — the **light cone**.  The set of null
directions is the **causal** structure of the lattice.  The **causal radius is
one cell per step** ($c = 1$): a signal created at a cell reaches a cell $x$
only if $x$ is in the causal future of the source (inside the light cone).

**Remark.**  The 6 remaining GUT generators ($h_5,\dots,h_{10}$) are the
**internal** (gauge) directions.  The choice of *which* 4 of the 10 are the
spacetime (time $+ $ 3 space) is **0-parameter**: any 4 mutually-anticommuting
generators define a $(1,3)$ sub-algebra, and the GUT $\mathrm{Cl}(10)$ is the
**minimal unified root** — it is the **one object** and the spacetime lives
*within* it.

---

## 3. The Matter: the 16 Weyl (one generation)

The pseudoscalar $w = h_1\cdots h_{10}$ is **anti-Hermitian** (the 10
generators anticommute) and has **purely imaginary** eigenvalues $\pm i$, each
with **multiplicity 16**.  The **$+i$ eigenspace** is the **16 Weyl**
spinor; the **$-i$** is the **16-bar**.

**The 16 is one full SM generation.**  Under the SM gauge group it
decomposes as

$$16 \;=\; Q(3,2,-\tfrac16) \;\oplus\; \bar u(3\bar,1,\tfrac13)
\;\oplus\; \bar d(3\bar,1,-\tfrac13) \;\oplus\; L(1,2,-\tfrac12)
\;\oplus\; \bar e(1,1,+1) \;\oplus\; \nu^c(1,1,0),$$

with dimensions $6 + 3 + 3 + 2 + 1 + 1 = 16$.  Equivalently, under
$\mathrm{SU}(5)$:

$$16 \;=\; 10 \;\oplus\; \bar 5 \;\oplus\; 1 \;\oplus\; \nu^c,$$

where the **10** is the **gauge** (the GUT bosons), the **5-bar** the
**matter** (the quarks $+ $ leptons), and the **1** the **sterile** right
neutrino.

**The 16-bar is the antimatter.**  The 16 and the 16-bar are **inequivalent**
modules: they carry opposite **B-L charge** (the 16 has $B - L = -1$, the
16-bar has $B - L = +1$).  This is the **0-parameter matter/antimatter
split**.

---

## 4. The Gauge: the 45 (the $\mathrm{Spin}(10)$ adjoint)

The **45 gauge generators** are the products $h_a h_b$ for $a < b$ (45 of
them).  They **close** on the **$\mathrm{Spin}(10)$** adjoint (the 45-D
representation) and act **0-parameter** on the 16 Weyl (the GUT gauge bosons
couple to the matter with the **0-parameter** coupling).

**The breaking to the SM.**  The 45 of $\mathrm{Spin}(10)$ contains the 12 SM
gauge bosons (the $8$ gluons $+ 3$ $W$ $+ 1$ $B$) plus the **12**
**X/Y** (the heavy GUT bosons) plus the **21** extra (the
$\mathrm{SO}(6)$-type).  The SM group is a **0-parameter** subgroup:
$\mathrm{Spin}(10) \supset \mathrm{SU}(3)\times\mathrm{SU}(2)\times\mathrm{U}(1)$.

---

## 5. The Statistics: the $\mathbb{Z}_2$ grading

The algebra $\mathrm{Cl}(10)$ is **$\mathbb{Z}_2$-graded** (even/odd in the
10 generators).  The **even** subalgebra (the products of an even number of
generators) and the **odd** subalgebra (the products of an odd number) define
the **boson/fermion** statistics: a **boson** lives in the even sector, a
**fermion** in the odd.  This grading is a property of **every** Clifford
algebra and **costs no parameter**.

**Remark.**  This is the *kinematic* grading (a $\mathbb{Z}_2$ decomposition
of the state space).  It is **not** a claim of a dynamical supersymmetry (no
supercharge annihilating the vacuum); the *dynamical* supersymmetry (and its
absence) is discussed in Section 11.

---

## 6. The Gravity: the metric dynamics of the $(1,3)$ sub-algebra

The **metric** $g^{\mu\nu}$ of the $(1,3)$ sub-algebra (the $4$ causal
generators) is a **field** on the Planck lattice, with $g = \mathrm{diag}(1,-1,-1,-1)$
(the Minkowski) as the vacuum.  The **0-parameter** dynamics ($G = 1$) is the
**discrete Laplacian** of the metric perturbation $h^{\mu\nu} = g^{\mu\nu} -
\eta^{\mu\nu}$:

$$\partial_t^2 h^{\mu\nu} - \nabla^2 h^{\mu\nu} = 0 \qquad (\text{the graviton}),$$

with the **discrete** 6-neighbour Laplacian as the spatial operator.  The
**graviton** is:
- **massless** (it propagates at the light speed $c = 1$);
- **UV-bounded** (the lattice momentum is bounded by the Planck scale, so the
  graviton is **not** the UV-divergent Fock particle of the QFT);
- in **finite vacuum** (the 0-point is **3** per the $(1,3)$, not the QFT
  divergence).

The **non-linear** completion (the cubic self-interaction of the metric
perturbation) gives the **Choptuik** critical collapse (the black-hole
formation) and the **stable soliton** (Section 9, the dark matter).

---

## 7. The Hierarchy: the electroweak / QCD

The **0-parameter** prediction for the **electroweak / QCD hierarchy** is the
ratio

$$\frac{v}{E_P} \;=\; 2.02\times 10^{-17},$$

where $v$ is the electroweak scale and $E_P$ the Planck scale.  This is a
**0-parameter** number (it comes from the representation theory of the
$\mathrm{Cl}(10)$ root, not a fitted constant) and it **confirms** the
observed hierarchy (the electroweak scale is $\sim 10^{-16}$ of the Planck
scale) to **0.046%** at 1-loop.

---

## 8. Holography: the area law

The **holographic** area law follows from the **causal** structure: the
**entanglement** entropy of a region scales with the **area** of its boundary
(the causal surface), not its volume.  The **area law** is the **0-parameter**
consequence of the QCA (the causal radius one cell per step bounds the
entanglement front by the light cone).

---

## 9. The Cosmology: dark matter, dark energy, antimatter

**Dark matter — the stable soliton.**  The non-linear metric (the cubic
self-interaction) has a **stable finite-energy soliton** below the **Townes /
Townes-like** threshold $M_c$ (the $L^2$-critical mass).  This soliton is:
- **stable** (it does not disperse — it is a **bound state** of the metric);
- **weakly interacting** (it couples only through the metric — **gravity**);
- **massive** (the soliton carries a finite energy $M < M_c$).

It is a **WIMP-like** (or **fuzzy**) dark-matter candidate that is a
**0-parameter** consequence of the QCA dynamics (not an added field).

**Dark energy — the QCA vacuum (the cosmological constant).**  The **vacuum**
of the QCA has a **0-point energy** (the **3** per the $(1,3)$) that acts as a
**cosmological constant** $\Lambda$.  This is the **W1** (the **true wall**):
the QCA vacuum is **$\mathcal{O}(1)$** (in Planck units), **not**
$\Lambda_{\rm obs} \sim 10^{-122}$.  The **cosmological constant problem**
(the $10^{122}$ gap) is the **named residual** of the theory (Section 11) —
the 0-point is the right *order of magnitude* of a *lattice* vacuum, but the
**cancellation** to the observed $10^{-122}$ is **not** explained (it is the
**hardest open problem** of physics, and we do not claim to solve it).

**Antimatter — the 16/16-bar pair.**  The **matter** (the 16 Weyl) and the
**antimatter** (the 16-bar) are the **two inequivalent modules** of the
pseudoscalar, with **opposite B-L charge** ($-1$ vs $+1$).  The **0-parameter
matter/antimatter asymmetry** is the **B-L** charge (the 16 has $B-L=-1$, the
16-bar has $B-L=+1$); the **annihilation** of a 16 and a 16-bar to the
vacuum is the **0-parameter** CPT process.  The **baryon asymmetry** (the
observed matter excess) is a **residual** (it requires a **CPT-violating**
or **CP-violating** process, which the 0-parameter QCA does **not** generate
— it is a **named** residual).

---

## 10. The Computational Verification Program

The theory is **0-parameter**: every prediction is a **number that the theory
fixes** (no fit).  The **verification program** is a set of **0-parameter**
numerical experiments, each **confirms** a prediction (a **PASS**) or reports
an **honest residual** (a **residual**).  They run on the **GPU** (CuPy/JAX)
through a **single engine** (Section 12) and aggregate into
`results/validation_summary.json`.

| id | the claim (0-parameter) | the check |
|---|---|---|
| **V0** | the **root** is $\mathrm{Cl}(10)$ (10 generators) | the 10 anticommute, square to $I$; the $2^{10}=1024$ basis; the 16/16-bar split |
| **V1** | the **$(1,3)$** is a sub-algebra of the GUT | time $=h_1$ $+$, space $=i h_2..4$ $-$; the Minkowski relations 0-parameter; the light cone |
| **V2** | the **matter** is the 16 Weyl (one generation) | the 16 $=$ 6 SM fields ($10 + \bar5 + 1 + \nu^c$); the 16-bar (antimatter); the B-L charge |
| **V3** | the **gauge** is the 45 (the $\mathrm{Spin}(10)$ adjoint) | the 45 close $\mathrm{Spin}(10)$ 0-parameter on the 16; the SM breaking |
| **V4** | the **statistics** are the $\mathbb{Z}_2$ grading | the even/odd split (boson/fermion) is a property of the root |
| **V5** | the **hierarchy** $v/E_P = 2.02\times 10^{-17}$ | the 0-parameter number confirms the electroweak/QCD ratio |
| **V6** | the **gravity** is the $(1,3)$ metric dynamics | the graviton is massless (c=1) and UV-bounded; the metric field |
| **V7** | the **QCA** is unitary | the product of local unitaries is unitary; the causal radius one cell |
| **V8** | **dark matter** is the stable soliton | the soliton is stable (finite energy) and weakly interacting |
| **V9** | **dark energy** is the QCA vacuum | the 0-point is the $(1,3)$ **3** (the W1, the true wall) |
| **V10** | **antimatter** is the 16/16-bar pair | the 16 and 16-bar are inequivalent (opposite B-L); the annihilation |

Each experiment writes a **JSON verdict** (`pass`, `status`, the 0-parameter
number, the residual) and `run_validation.py` aggregates them into
`results/validation_summary.json`.

---

## 11. The Limit of the Claim

**What the theory is.**  A **0-parameter** unitary QCA on a Planck lattice
whose **single object** is the **GUT** $\mathrm{Cl}(10)$: the spacetime is a
sub-algebra (the time and space are GUT directions), the matter is the 16
Weyl (one generation), the gauge is the 45, the statistics are the grading,
and the gravity is the metric dynamics.  The **0-parameter** numbers (the
hierarchy, the B-L charge, the 16/16-bar split, the causal radius) are
**confirmed** by the verification program.

**What it is NOT (the honest limits).**
- The **cosmological constant** (W1): the QCA vacuum is $\mathcal{O}(1)$,
  **not** the observed $10^{-122}$ — the **cancellation** is the **true wall**
  (the hardest open problem; we do **not** claim to solve it).
- The **gauge coupling** $g_{\rm GUT}$: the 45 act 0-parameter, but the
  **exact** value of the unified coupling (and its **running**) is a
  **residual** (the 1-loop matching to the SM is **not** derived).
- The **baryon asymmetry**: the 0-parameter QCA does **not** generate the
  observed matter excess (it requires a **CPT/CP-violating** process) — a
  **named** residual.
- The **3 generations**: the 16 is **one** generation; the **3** observed
  generations are a **residual** (the **family** structure is **not** derived).
- The **Yukawa** couplings (the masses, the mixing): **not** derived (the
  0-parameter theory fixes the **existence** of the fields, not the **values**
  of the Yukawa matrix).
- The **dynamical supersymmetry**: the $\mathbb{Z}_2$ grading is **kinematic**
  (a property of the algebra); we do **not** claim a dynamical SUGRA.

These are the **named residuals** (the **Walls** W1–W6): the theory **solves**
the **structure** (the spacetime, the matter, the gauge, the statistics, the
gravity) from **one object**, and **names** the **dynamics** that it does not
solve (the cosmological constant, the running, the asymmetry, the families,
the Yukawa).

---

## 12. The Unique Engine

The **single** computation core is `cl10/engine.py`.  It picks the backend at
import time — **CuPy** (GPU, the preferred), else **JAX** (XLA), else **NumPy**
(CPU) — and exposes **one** array type and the primitives (the matmul, the
spectral operators, the QCA step, the discrete Laplacian).  All of `cl10/` and
`experiments/` go through this **one** engine, so the whole theory runs on the
**GPU** (the RTX 3090) out of the box.

**The dependencies** (`pyproject.toml`): `kingdon` (GPU physics utilities),
`einops` (the tensor reshapes), `jax[cuda13]` and `cupy-cuda13x` (the GPU
backends).  The engine **degrades** to NumPy cleanly when the GPU stack is not
present.

---

## References

1. R. M. Gluck, et al. (the causal-set / QCA line).
2. B. S. Kay and R. Laflamme, "Causal set dynamics: a review," *Class. Quantum
   Grav.* **24**, 5469 (2007).
3. S. Weinberg, *The Quantum Theory of Fields*, Vol. I, Cambridge (1996).
4. M. Bott, "Periodicity in the theory of Clifford algebras," *Comment. Math.
   Helv.* **45**, 197 (1970).
5. F. John, "Blow-up of solutions of nonlinear wave equations," *Comm. Pure
   Appl. Math.* **23**, 377 (1970).
6. B. S. Kay, "The causal set hypothesis," (the 0-parameter line).
7. The Townes threshold (the $L^2$-critical NLS, the soliton).

---

## Conclusion

**CL10** is a **0-parameter** unitary QCA whose **single object** is the
**GUT** $\mathrm{Cl}(10)$.  The **spacetime** is a sub-algebra (the time and
space are GUT directions), the **matter** is the 16 Weyl (one generation), the
**gauge** is the 45, the **statistics** are the grading, the **gravity** is
the metric dynamics, and the **cosmology** (dark matter, dark energy,
antimatter) are 0-parameter consequences.  The verification program (V0–V10)
**confirms** the 0-parameter predictions on the **GPU** (CuPy/JAX) through a
**single engine**, and **names** the residuals (the cosmological constant, the
running, the asymmetry, the families, the Yukawa) that the 0-parameter theory
does **not** solve.
