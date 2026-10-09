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
│   ├── api.md                # the complete CL10 API reference
│   ├── figures/              # the 36 publication figures (make_figures/; fig 19 = young QCA, fig 20 = V62 product)
├── publication/              # the journal build: doc/publication.md -> REVTeX 4.2 PDF
│   ├── build.sh              # one command: preprocess -> pandoc -> pdflatex x3 -> publication.pdf
│   ├── preprocess.py         # manuscript -> pandoc-friendly markdown (Divs for figures/tables/refs)
│   ├── filter.lua            # Divs -> figure / figure* / table* / thebibliography
│   ├── template.tex          # revtex4-2, aps/prd, twocolumn
│   ├── fetch_dois.py         # resolves + verifies the 53 refs against Crossref -> doi.json (clickable DOIs in the PDF)
│   ├── README.md             # the build pipeline, files, prerequisites
│   └── publication.pdf       # the compiled paper (39 figures, 8 tables, 53 refs)
├── make_figures/             # the figure generators: `python -m make_figures` regenerates all 39 (see make_figures/README.md)
│   ├── __init__.py           # the registry in paper order + GPU pinning + the V62 subprocess (fig 20)
│   └── algebra/scales/couplings/proton/dynamics/residuals/style .py
├── cl10/                     # the unique engine (GPU-first, CPU-fallback)
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
├── experiments/              # the verification program (117 experiments, V0-V116)
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
│   ├── V61/                # the 3-D hydrogen atom (Fig. 20; fast 120^3 check): the gaps ARE the Rydberg differences (rel. error 0), the E1 selection rules come out of the lattice wavefunctions, the absorption peaks on the gap and falls with the detuning, the excited atom radiates while it decays.  The ATOM is an input (real hydrogen); the PHOTON is CL10
│   ├── V62/                # the PRODUCTION hydrogen experiment on the GPU (296^3, the full photon script): every Fig. 20 claim measured at production resolution, the figure regenerated as a product of the experiment
│   ├── V63/                # the Brillouin-zone chirality test (W9): the 16 sector is closed by the FULL local step (leakage 1e-15), the Bloch Floquet matches the 4^3 torus (4e-15), the zone centre is GAPPED (0.0365), the 12 doubler momenta are GAPLESS (|v| = 0.142), the node chiralities sum to 0 — the doubling survives
│   ├── V64/                # the desert inverse problem (W4, F2, F5): with alpha_X = 1/25 FIXED (the root value, not refitted), the inverse problem asks which vectorlike content makes the measured couplings meet; the required Delta b is non-universal; the lightest CL10-realizable solution (three 10+10bar + one (1,2) pair) meets at M_X = 1.8e16 GeV, cures the SU(3) down-pole, predicts tau_p = 1.8e36 yr; M_GUT = M_P is abandoned
│   ├── V65/                # the internal Wilson factor (W9 -> W4): W(q) = exp[-i w sum(1-cos q_j)(cos(th) X + sin(th) Y)] from the root's own Cartan algebra; exactly unitary, W(0) = I, quasi-local (measured Lieb-Robinson front); gaps all 12 doubler nodes; X, Y identified by charge multisets (both SU(5) embeddings); the selective 3-node reading solves the meeting equations at M_X = 8.9e14 GeV — EXCLUDED by Super-K (reported, not hidden)
│   ├── V66/                # the dark-matter identity scan (W7): the mirror Z2 is exact (leak 1e-15), the 16bar weights are exact conjugates; the only exact-stability carrier is the mirror sector, lightest state a colored Qbar at 1.3e8 GeV (excluded); NO viable 118 GeV carrier in the implemented spectrum
│   ├── V67/                # the breaking chain (W6): the Wilson direction at the V65 root leaves Spin(10) -> SU(3)xSU(2)xU(1)^3; every unbroken root is intra-fragment (the Model algebra is preserved), the 32 broken roots are exactly leptoquark-type; residual: two extra U(1)'s
│   ├── V68/                # the W2 audit: meeting roots need n_nodes >= 3 (1 and 2 excluded — a lower bound from unification); the diagonal node pattern gives only integer seesaw coefficients, the democratic one exactly 1/2 (0.073 eV, falsifiable); the exact 1/3 needs a 3x-enhanced Majorana eigenvalue no node structure supplies
│   ├── V69/                # the W8 audit: the tower U^(2^p) is exactly kinematic (phases 2^p w mod 2pi, verified p <= 10); no (16pi)^-p suppression at any level — the tower generates NO loop factor; the missing tool is named (a genuine lattice RG)
│   ├── V70/                # the W3 audit: the 12 doubler nodes form 3 families of 4 labeled by the zero axis, cubic-transitive (spectra equal to 5e-15); the index 3 EXISTS and labels the unification-required 3 split nodes, but does not force the count
│   ├── V71/                # the W5/W1 stability audit: the Wilson-extended step keeps the exact spectral mirror (1e-15), the +- pairing, and the zero point (Tr D = 0); the Sakharov no-go and the tree-level Lambda = 0 are UNCHANGED by V65
│   ├── V72/                # the missing tool BUILT: exact block-spin RG of the bond layer (Schur decimation, U_BB = 0, no approximation), 5 levels; level-1 identity onsite = 2(GR)^2 exact; no level follows (16pi)^-p for the factorized OR the non-factorized ([G,R] != 0) gate — single-particle RG composes angles, never divides them; the loop factor needs a gauged bond layer (many-cell RG)
│   ├── V73/                # the bond-gate loophole closed: 3 Haar-random Spin(10) gates + structured R = G^dagger through the V72 decimation — no realization follows the (16pi)^-p law; the Cartan content depends only on path products; W8's missing tool uniquely specified
│   ├── V74/                # the W4 structural closure: while all node thresholds sit below M_X the meeting equations depend on the three node strengths ONLY through their sum — the Jacobian has three identical columns (rank 1, split-invariance exact to 1e-14); the node mechanism's meeting point M_X = 8.9e14 GeV (tau_p 1.1e31 yr, F5-excluded) is FORCED for every split; the data-preferred spectrum is not a three-full-16 pattern: a named incompatibility
│   ├── V75/                # the W7 no-go completed: no Wilson node hosts a 118 GeV neutral mirror — for every Wilson angle the same node puts a stable colored/charged fragment below the desert (0 survivors over 720 angles); the minimal named fix: one mirror-odd ladder scalar, portal lambda_* = 4e-4 for Omega h^2 = 0.120, sigma_SI = 5e-51 cm^2 (nearly decoupled; falsified above 1e-50)
│   ├── V76/                # the W2 junction BUILT: a cubic-symmetric zero-mean Wilson-strength junction whose Fourier support closes the selected node triangle exactly (verified on the full 12x12 doubler difference table); node mass matrix c_s(2w0 I + kappa(J-I)), democratic at kappa = 2w0 ({3M,0,0}); the 1/3 candidate survives only at that measure-zero point and only with democratic m_D (else gapless colored doublers)
│   ├── V77/                # the W8 division mechanism BUILT: one quantum link d.o.f. per bond (dimer), traced over — the level-1 Cartan content is exactly cos(2*theta) times the V72 value (dev < 1e-15, tunable); a traced FLUCTUATING bond d.o.f. divides angles where every gate only composes them; the gauged Spin(10) adjoint link is now backed by a demonstrated mechanism
│   ├── V78/                # the W1 protection theorem: det U = 1 pins the quasi-energy sum against every traceless internal deformation (50 random Hermitian A, eps to 1: |S| < 4e-14); no branch crossing below eps = pi; the only internal knob is the U(1) trace phase — the observed Lambda is equivalent to a trace-phase breaking of size 1.6e-120
│   ├── V79/                # the GAUGED BOND LAYER BUILT: one Spin(10) adjoint link d.o.f. per bond + matter 16 on sites — the local Gauss constraint closes as an exact operator identity (err 0), link algebra exact, Casimirs 64/45 exact, step unitary 3.8e-16, causal radius 1; gauge-invariant content 32 + 288 + 1120 = 1440 states, no singlet (triality); the single-site non-local control leaks at 1.9
│   ├── V80/                # the W8 division on the ADJOINT link measured: the V77 qubit-link mechanism reproduces exactly (1.1e-16) and the traced two-hop contracts (max |eig| 0.59), but the level-1 factor at the root coupling is F = 1.06 and F_min = 0.52 over the anisotropy — a factor 26 above 1/(16 pi); level-N suppression is the named next test
│   ├── V81/                # the W4 content test: the gauged layer's composite fermions (16 x 45 = 720 states, five net chiral 16-copies, shifts (145.3, 145.3, 226.3), flux-pattern thresholds m = w(c0 + |alpha.v|)) do NOT deliver the data-preferred meeting at alpha_X = 1/25 — 0 physical solutions over the full scan of the breaking direction; the V64 spectrum is not the gauged layer's own content
│   ├── V82/                # the W4 content test CLOSED: full-weight / broken-flux / light-composite ansatz families and a relaxed alpha_X all fail (min residual >= 5.7, 0 physical solutions); the obstruction is the exact identity sum db1 = sum db2 = 144.0 (embedding index (3/5) sum Y^2 = sum T2 on 16 x 45) — the content shifts U(1) and SU(2) equally, the data require different shifts
│   ├── V83/                # the W8 level-N test RUN: the bond hop factorizes through the joint Casimir T.L (three irreps 32+288+1120, projectors 1e-15), the traced-link decimation iterates exactly to level 4 (V80 reproduced to 2e-16); the division DOES NOT COMPOUND: F_p = 1.06, 1.90, 0.75, 1.01 at the root (7e6 above the loop law at p=4)
│   ├── V84/                # the W4 EQUAL-SHIFT THEOREM: sum db1 = sum db2 is the index-proportionality of Tr Y^2 and Tr T3^2, verified over the 16, the flux 45, the vector 10 and the SU(5) 10 (ratios 1.0 to 1e-6); the meeting mismatch can only come from SU(5)-fragment content - the V64 repair's mismatch is carried entirely by its (1,2) pair (db1/db2 = 3/5), which the automaton does not supply
│   ├── V85/                # the W8 CLOSED-LOOP factor: the traced hop matrix contracts (rho(MG) = 0.776 vs 1) and closed n-hop loops suppress geometrically (r = 0.54); at the root coupling F_10 = 0.0192 vs 1/(16 pi) = 0.0199 - the loop factor at the ladder length 10 (3.6%); V80 regression 2e-16
│   ├── V86/                # the W4 constructive attempt: the junction-split triangle nodes as desert supplier - the Wilson direction splits the node 16 by a factor 55 (cQ +0.005, cL -0.015 vs partners ~0.82): the (3,2) and (1,2) repair states ARE automaton states (Q, L components of node 16's); the junction ratio rho breaks V74's rank-1 degeneracy; the fragment pattern is fabricated (Q, L at 4-8e15 GeV, partners above M_X) but the meeting misses: best physical residual 0.94, alpha_X = 0.030 vs 0.04 (n_nodes 3-12 scanned) - a named near-miss
│   ├── V87/                # the W4 coefficient-pattern release: the GENERAL zero-mean junction (three free branch scales a_i = 2 w0 (1 + r_i), r1 + r2 + r3 = 0) scanned over the full (r1, r2) plane + n_nodes + scale + angle: the best physical residual EQUALS the democratic floor (0.81), alpha_X = 0.090 (factor 2.3 off) - the near-miss is a property of the junction MECHANISM, not of the democratic ansatz; closure needs fragment content beyond the node 16
├── V88/                # the W7 named fix BUILT and audited: the mirror-odd ladder scalar (118 GeV) is exactly stable under either portal (scalar parity (-1)^Ns, odd sector untouched); but the exact mirror Z2 forbids h* -> f f-bar: ss -> hh closed (gap 13.2 GeV, virtual cloud 1.5e-8 = theory), only ss -> hh* open: Omega h^2 = 0.120 needs lambda = 0.13, sigma_SI = 6.8e-46 cm^2 - EXCLUDED by direct detection; the V75 numbers (4e-4, 5e-51) hold only with a portal that BREAKS the mirror Z2 at O(lam) - a named tension (F7)
├── V89/                # the W8 selection rule measured to length 30: the closed-loop ratio F_n crosses 1/(16 pi) at EXACTLY n = 10 (the ladder length is the crossing); NO mod-4 plaquette selection (every length contributes); the geometric law F_n = A r^n (r = 0.698) extrapolates to F_170 = 1.1e-26 vs (1/16 pi)^17 = 1.2e-29 - the 17-factor composition misses by x929, a named residual of its own
├── V90/                # the W2 protection test: the cubic group induces ALL of S_3 on the junction triangle and its commutant is exactly span{I, J} (nullity 2) - no symmetry protects the democratic point; w0 and kappa share ONE Wilson period (38.19) so r = kappa/w0 is NOT quantized; the cosmological window sum m_nu < 0.12 eV is EXACTLY the point r = 2 (1 passing value in [1,3] at step 1e-4) - the sharpest form of the measure-zero tuning
├── V91/                # the W3 forcing test: the cubic group acts TRANSITIVELY on the 12 doubler nodes (one orbit of 12 - the V70 axis families are not orbits), so the only cubic-symmetric split is all 12 nodes; the rank-1 meeting rigidity holds for n = 4, 8, 12; the fixed scales M_X(3..12) = 8.9e14..2.9e15 GeV all miss the observed 2e16 (F5); the 3-node triangle has a cubic orbit of 8 - the family count 3 is a symmetry-BREAKING input, its symmetric alternative excluded by computation
├── V92/                # the W4 feasibility test: the data-preferred meeting point (12.7, 21.2, 20.1) IS inside the convex hull of the six node-16 fragment beta-vectors (barycentric witness) - the required content exists; junction realizability was the remaining test
├── V93/                # the W4 Cartan-junction test: a GENERAL Cartan direction realizes the meeting exactly (6 solutions, desert thresholds) but breaks 120/120 roots including the SM ones - unphysical; the SM-preserving direction space is 2-dimensional (the class-degenerate span{X,Y} family)
├── V94/                # the W8 composition law: EXACT multiplicativity across LEGS (the x929 was a chain-vs-legs category error); on a CHAIN (MG)^{10k} is not multiplicative (R_17 = 2.0e+08) - legs and chains are different objects
├── V95/                # the W4 junction-split test: NO exact meeting with a universal kappa (uniform span{X,Y}, n = 3..12, 300 starts) - the split is invisible at first order (sum_l 2cos(2 pi l/n) - 1 = -n); the V92 witness needs INDEPENDENT per-class thresholds
├── V96/                # the W4 constructive closure: the flavour-dependent bond layer (kappa_s per fragment class, 6 fitted coefficients) supplies an EXACT meeting (residual < 1e-10) at M_X = 1.8e16 GeV, alpha_X = 1/25, for 10/10 node counts, all thresholds in the desert, tau_p = 1.9e36 yr, SM desert group preserved
├── V97/                # the entanglement structure measured: entanglement is the BOND layer, not the shared cell - same-cell fermions carry S = 0 (Slater), the internal-only step keeps S_A = 0 forever, the full step drives it to 0.783, and the entanglement front crosses radius r at exactly t = r; the implemented step is QUADRATIC (Slater stays Slater) - the named gap to the area law
├── V98/                # the PHYSICAL SECTOR of the gauged layer: CENTER OBSTRUCTION - dim ker = 0 with one 32/site (adjoint flux cannot screen a spinor; gap 146 on the two-site chain); minimal sector 2 matter/site: m_45(32x32) = 2, intertwiners dim 2, glued physical sector dim 4, all exact singlets (C2 < 1e-25); S(pair|link) = ln 45 exactly, bounded per bond and boundary-supported - the SHAPE of an area law, Bekenstein coefficient unnamed
├── V99/                # the multi-bond area law of the gauged layer: 3-site chain physical sector dim ker = 24 (flux decomposition N_1 = N_45 = N_210 = 1, m = 2 each); additivity ln 45 across a single bond EXACT; two-bond region S = ln 45 not 2 ln 45 (vertex mutual information ln 45); Bekenstein gap as a number: ln 45 vs 1/4, factor 4 ln 45 = 15.23 - named, open
├── V100/               # the vertex correction of the gauged area law: I_n = dim(45^⊗n)^G by Haar averaging = 1, 1, 6, 23, 141 (I_2 = I_3 = 1 validates); 24-state kernel orthonormalised by its Gram (rank 24); middle-cut capacity S_cap = 4.4998 nats vs naive 2 ln 45 = 7.613 (vertices subtract 3.11); per bond 2.25 nats, factor 9.0 from 1/4; cubic-lattice vertex-corrected coefficient (1/3) ln I_6 = 1.65 nats/bond, factor 6.6 - the Bekenstein deficit moved 15.2 -> 6.6
├── V101/               # the first ground state of the gauged layer: bond-current Hamiltonian H = sum J.L of V79 diagonalised exactly (generalized eigenproblem, Gram metric) on the 24-state kernel; intertwiner identity I† J_A I = L_A/45 to 3e-17; the interaction SPLITS the kernel: E0 = -5819.18 eta, degeneracy 6, gap 58.59 eta; ground-state middle-cut entanglement S = 3.8067 = ln 45 exactly - the interaction selects the flux line
├── V102/               # the node cone tensor measured: at all 12 doubler nodes and 1200 sampled directions the dispersion is LINEAR (p = 1.00) but strongly ANISOTROPIC - principal velocities (0.200, 0.200, 0.284), spread 0.75; node gap 2e-4 at the root angle (Fig. 38)
├── V103/               # the W2 re-audit with the IMPLEMENTED democratic rank-1 m_D (V28/V45): the V90 measure-zero point was an artifact of the identity-m_D assumption; the seesaw collapses to m3(r) = S17/(1+r), the massless pair is protected for EVERY r, the bound window is r > 0.217 (width 2.78), and the atmospheric scale selects r = 1.91 (m3 = 0.0501 eV)
├── V104/               # the W6 mixing manifold: the three mixing deviations tested as ONE effect over the natural U6 vev placements (12 branches x 2 orientations, democratic base + rank-1 correction): every placement excluded (best residual 0.308 vs the needed pattern) - the deviations are not a natural-vev artifact; the U6 form itself is pending
├── V105/               # the W4 Cartan test: are V96's six kappa_s images of ONE Cartan element? NEGATIVE - SM preservation forces span{X,Y} at theta* (class-degenerate, meeting residual >= 0.43) and every Cartan element splitting the six classes breaks all 8 unbroken desert roots: the flavour dependence is off-diagonal bond structure, not a Cartan dressing
├── V106/               # the W5 washout audit: at the SAME Yukawa giving m_nu = 0.0501 eV (y = 0.40) the decay is IN equilibrium (Gamma/H = 47) and the physical efficiency kappa_f = 5.4e-3 (Buchmuller-Di Bari-Plumacher) sits 14x from the guessed dilution (1/16pi)^2 - the 11% Delta B agreement is an accident of two wrong factors
├── V107/               # the W6 orthogonality test: V96's kappa_s are per-GAUGE-class (generation-blind) - PMNS/CKM exactly invariant under the whole kappa pattern (max angle shift 2e-13): the W4 bond-layer fix cannot supply the mixing deviations; W6 needs generation-space structure
├── V108/               # the hierarchy coefficient un-fitted: v/M_P = 3*F_10^10 with F_10 = 0.019186 measured (V85), 3 = the axis families (V70), 10 = the crossing (V89): 2.088 vs 2.076 observed (99.4%) - conditional on the axis-family factorization of the closed-loop sum
├── V109/               # the W5 CP search in the EXTENDED structure: the mirror is the CP-forbidding structure - the Wilson extension preserves the +- pairing (2.4e-15, V71 regression) and supplies NO phase; every class-dependent dressing (general kappa_s layer, general real-symmetric bond generator, both mirror sign conventions) BREAKS the one-cell mirror (dev ~ 1): a CP phase costs Theorem 4
├── V110/               # the W4 off-diagonal bond program V105 demanded: 80 root-vector bond channels, one Hermitian bond operator (12 parameters) - the meeting does NOT close (floor 30.4 at n=3 over n=3..12); level repulsion forces 2 of the 6 kappa_s signs: the pattern stays a named input
├── V111/               # the W6 generation-space dressing: a scalar dressing is exactly blind (1.4e-13); a COMMON non-abelian family dressing (7 parameters) reaches the measured mixing triple; chiral-asymmetric (13 parameters) reaches it as a fit - W6 is reachable only at named-input cost
├── V112/               # the W7 dark-U(1) portal: the extra U(1) of the V67 chain (E ~ (-1,-1,1,1,1)/sqrt5) has charged desert matter but is UNHIGGSABLE inside the spectrum (all E-charged 45-directions are colour triplets); the named dark-Higgs window is EMPTY (relic needs gD*eps >~ 0.02, direct detection allows <~ 1e-3)
├── V113/               # the W8 axis-family factorization MEASURED: the 12 node channels carry exactly equal loop factors (across-family spectral identity 0.00); the hierarchy coefficient is a node-set choice - one split node per family: 2.0881 vs 2.0765 observed (99.44%), all 12: 8.35, single: 0.696
├── V114/               # the W2 junction ratio r = kappa/w0 NOT derived: the V76 modulation's exact stabilizer is a 12-element subgroup (not the 24 chiral cubic); NO cubic-symmetric closed junction exists on (Z/4)^3 (closed space trivial); inside the V76 stabilizer the closed space is 2-dimensional - w0 independent of kappa: r stays free
├── V115/               # the W9 cone-isotropy scan: the node gap g(w) of the V65 Wilson family located at every zero (grid + golden-section) over 2 reps x 60 angles - 2726 gapless points, smallest cone anisotropy A = 0.997: no isotropy point on the family (W9 confirmed beyond the V65 root; Fig. 39)
├── V116/               # the W3 triangle orbit-8 audit: the 64 one-per-family triangles split into 4 cubic orbits (8, 8, 24, 24); a triangle's stabilizer is exactly S3 (order 6) - the split reduces O_h -> S3 and the 8 choices of an orbit are physically distinct vacua (8-fold degeneracy, domain walls)
├── simulation/               # the 3-D visualizations on the CL10 engine
│   ├── young/                # the Young double-slit as a single-electron measurement (Fig. 19)
│   └── hydrogen/             # the 3-D hydrogen atom + the CL10 photon (Fig. 20): the levels are the real Rydberg eigenstates (an INPUT), the photon is the CL10 free massless field; the video/stills are `python -m simulation.hydrogen`, the Fig. 20 panels are the V62 product
├── run_validation.py         # the runner: runs all 117, writes the hard-tier verdict (`./run_gpu.sh` runs the full suite on the GPU)
└── results/
    └── validation_summary.json   # the 117/117 PASS + the hard-tier verdict

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
python3 run_validation.py        # runs all 117 in parallel (--jobs N), prints the verdict
./run_gpu.sh                     # the full 117 on the GPU (cupy, 12 jobs); numpy fallback per experiment if the card is full
python3 -m make_figures          # regenerates ALL 39 figures -> doc/figures/
                                 # GPU if the pinned card is free, else numpy;
                                 # --only N [N..] | --fast (skip 19+20) | --cpu
python3 -m simulation.hydrogen   # regenerates the hydrogen video + stills (the Fig. 20 scene)
                                 # -> simulation/hydrogen/video/  (git-ignored)
bash publication/build.sh        # rebuild the journal PDF from doc/publication.md
                                 # needs pandoc >= 3 and texlive-publishers (revtex4-2)
                                 # -> publication/publication.pdf
```

## The verification program

`doc/publication.md` states the **seven criteria** a unification claim must
meet (Section 1.2, with the paper's own verdict on each), the **tiered**
validation against them, the **nine open residuals** (W1–W9) and the **seven
falsifiers** (F1–F7).  `run_validation.py` runs all **117** experiments
(V0–V116) across the seven tiers and writes

> **hard-tier verdict** = all of tiers 0–3 (V0–V7, V11–V13) PASS — the single
> object, the kinematics, the dynamics, the **anomaly-free** matter, the
> **$5/3$ normalisation**, and **matter–gauge–gravity in one algebra**.  The
> cosmology (tier 4) and the falsifiability (tier 5) are supporting; the named
> residuals are reported, not hidden.

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
chirality** of the implemented step: the 16-sector Brillouin-zone Floquet
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
