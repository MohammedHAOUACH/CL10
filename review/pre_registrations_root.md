# Root-rebuild programme: pre-registrations (branch root-rebuild, 2026-11-01)

The branch's question is the one the consolidated paper ends on: the 36-level
module of Cl(10,0) is the hypothesis of everything proved there ("the
hypothesis of everything proved below is this algebra and this module -
nothing else"), and the transcendence theorem (V334) is an arithmetic
statement about matrices of THIS module: Gaussian-integer entries, frozen
algebraic internal factor. Every number the suite fails to forge (the clock
g_obs = 0.19999798835751165, the desert crossing g_cross = 0.20013316730041675,
the ladder base LOOP = 0.019894367886486918, the coefficient
COEF_OBS = 2.076564111173832, the rung pair) lives in the declared periphery
of that axiom. A root that is to REPLACE Axiom 1 must therefore do two things
at once: force the skeleton (the Main-Theorem list as corollaries) and forge
the numbers as periods/eigenvalues of the same object.

Candidate N1, booked here as the first candidate of the programme, is the
object where the two arithmetics the construction already touches meet in one
datum:

  THE BINARY ICOSAHEDRAL GROUP 2I AND ITSMcKAY CORRESPONDENCE WITH E8.
  Facts of the candidate (all checkable, none assumed):
  - 2I = the unit quaternion group of order 120 = 2|A5|, element orders
    {1,2,3,4,5,6,10}: the digit 5 native;
  - McKay correspondence: 9 irreps with dimensions {1,2,2,3,3,4,4,5,6},
    sum d_i^2 = 120 = |2I|, tensoring by the 2-dimensional generating
    rep with adjacency = the affine E8 Dynkin diagram: the binary
    icosahedral group IS the E8 singularity;
  - E8: 240 roots, Weyl group |W| = 2^14 3^5 5^2 7 = 696729600 (5^2
    native), Coxeter number h = 30, exponents {1,7,11,13,17,19,23,29},
    Coxeter element of exact order 30; its Coxeter-plane projection of
    the 240 roots is the icosahedral 30-gon figure with Z[phi] radius
    ratios: the "8" (Bott/Clifford, the E8 lattice IS the D8 even
    unimodular lattice) and the "5" (icosahedral, phi) in ONE object;
  - the CL10 module resonance: E8's 248 under D5 x D3 branches with
    dimension identity 45 + 15 + 60 + 64 + 64, whose Spin(10) factors
    contain the Weyl pair 16 + 16' whose sum is exactly CL10's 32-level
    Dirac module S, and the 45 = the plane-generator count of CL10.
  What N1 is NOT: a proof that CL10's module IS the branching - only the
  dimension identities and the structure are certified here; the
  representation-content statements are standard branching tables, cited,
  not re-derived.

The tribunal (booked targets, read from the consolidated paper's artifacts
BEFORE either run): g_obs, g_cross, COEF_OBS as booked above; LOOP; n* = 208.

## V336_two_five_root (tier 2, hard; root-rebuild N1; the structure certification)

BOOKED CHECKS (all must pass for the positive branch word):
T0a the binary icosahedral group is built EXACTLY over the biquadratic field
    K = Q(sqrt5, s), s^2 = 10 - 2 sqrt5, as unit quaternions from the two
    generators q1 = (1+sqrt5)/4 + (s/4) i (order 10) and q2 = (1+i+j+k)/2
    (order 6): closure by BFS, count exactly 120, all elements unit-norm in
    K, order multiset {1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24};
T0b representation-ring/McKay certification by integer linear algebra:
    d = (1,2,2,3,3,4,4,5,6): sum d_i^2 = 120; the affine E8 Cartan (booked
    9x9) times d equals 0 exactly; the McKay recursion 2 d_i = sum_j N_ij d_j
    holds for every vertex with N the affine-E8 adjacency; N symmetric;
T1 E8 root system exact: the 8 simple roots (half-integer standard form),
    Gram = E8 Cartan, Cartan determinant 1 (Bareiss exact); the 240 roots
    enumerated (112 + 128) and every root's squared norm 2;
T2 Coxeter: exact Fraction matrix C (product of the 8 simple reflections in
    the standard coordinate basis), C^30 = I exactly, and C^k != I for every
    proper divisor k of 30 (so order exactly 30 = h); Weyl order via the
    exponent formula prod(e_i + 1) = 696729600 = 2^14 3^5 5^2 7 (exact
    factorization printed);
T3 Coxeter plane at 50 dp (adjugate null vector of C - lambda I for
    lambda = exp(2 pi i / 30)): the 240 projections fall on exactly 4
    distinct radii, multiplicities summing to 240; the angle set is closed
    under rotation by 2 pi/30 AND by 2 pi/10 (cluster tolerance 1e-25); the
    four radii ratios to the smallest are named by PSLQ (degree <= 6,
    coefficients <= 1e6) and the polynomials printed - booked expectation:
    the ratios lie in Z[phi] (degree <= 2 each);
T4 the volume vol(S^3 / 2I) = 8 pi^2 / 120 = pi^2 / 15 evaluated at 50 dp
    (the 5 in the denominator of a volume: number-theoretic, not geometric);
T5 the module resonance dimension identities: 45 + 15 + 60 + 64 + 64 = 248
    (E8 over D5 x D3); 16 + 16 = 32 (the Spin(10) Weyl pair against CL10's
    S = 32 as booked by Axiom 1); binom(10,2) = 45 = the plane count.
BRANCH WORDS: ROOT-2I-E8-STRUCTURE-EXACT (all green) / ROOT-STRUCTURE-FAIL
(any red; no class change is possible mid-run).
Cost booked: a few seconds (exact algebra + 50 dp 8x8 linear algebra).

## V337_root_period_census (tier 2, hard; root-rebuild N1; the numbers census)

The second half of the acceptance test: does N1 FORGE the tribunal numbers,
or does it only rhyme with them? Dictionary tiers booked BEFORE the targets
are looked at, V335's guards adopted verbatim: all comparisons in exact
50-digit mpmath, nonzero-numerator guard on vanishing products, convergents
never evaluated through float64 rounding at large magnitude.

TIER D0 (geometric periods allowed - the honest floor): the multiplicative
monoid generated by {vol S^2 = 4 pi, vol S^3 = 2 pi^2, vol S^7 = pi^4 / 6,
vol(S^3/2I) = pi^2/15, phi, sqrt5, h = 30, 240, 120, 12, |W| = 696729600,
primes^e (p <= 7, |e| <= 6), rationals p/q (q <= 1e6)} closed under
inverses. NOTE booked in advance: D0 can name LOOP = 1/(16 pi) =
(4 vol S^2)^-1 trivially once pi is admitted; a D0 hit is therefore NOT
scored as a discovery, it is scored as "the geometric tier is not selective"
and the margin of the best D0 hit against each target is reported anyway.
TIER D1 (pi-FREE structural tier - the actual test): the ring generated by
{integers <= 300, divisors of |W|, the exponents + 1, 240, 120, 60, 45, 32,
30, 16, 12, 8, 6, 4} under multiplication and +- , times phi^k (|k| <= 6),
times the pi-free cyclotomic units of Q(zeta_5)^+ and Q(zeta_6)^+ (phi,
1/phi, sqrt3: native to 2I's orders 5 and 6). No pi, no volumes: a D1 name
for the clock would be a NUMBER-THEORETIC birth, the character sheet's
requirement.
BOOKED TESTS, per target (g_obs, g_cross, COEF_OBS, LOOP):
(a) best D0 and D1 name with margin at the 1e-12 gate;
(b) PSLQ over [1, x, ..., x^6], height <= 20, 60 dp, on x = target
    (V309's pattern with 1e-30 substitution verification);
(c) the booked Coxeter probe: 30 * g_obs against 6 - is the clock 6/30 to
    structural precision? booked expectation from V206's own numbers: off by
    ~6.0e-5, five orders above the 1e-12 gate, so NAMED-NO;
(d) n* = 208 additive census: sums of +-c (c <= 4) over the structural
    integers {240,120,60,45,32,30,16,12,11,13,17,8,6,4}: every hit printed
    with the number of total hits as the density guard (a single hit among
    thousands proves nothing - V335's O(1)-per-number lesson).
BRANCH WORDS: ROOT-NUMBERS-NAMED (any target named at D1 within 1e-12) /
ROOT-NUMBERS-UNNAMED (all targets miss at D1; D0 LOOP hit reported as the
booked non-discovery).
Cost booked: under a minute.

## What the programme's first verdict will mean, stated before running

N1 passes the skeleton half of the acceptance test or it does not; either
answer is a fact about the first candidate. If both halves are green, the
branch has found an object that must be developed into a construction, not a
paper claim. If the numbers half is red (the booked expectation given V335's
experience: NO period of a known root names g_obs at 1e-12), N1 is
adjudicated as a skeleton resonance, the tribunal keeps its teeth, and the
programme books N2 with what N1 learned. No branch word will call N1 "the
new theory" on rhymes.

## V337_2i_engine (tier 2, hard; root-rebuild; THE ENGINE ITSELF CHANGED)

The user directive on this branch: do not only certify the candidate -
CHANGE THE ENGINE and use it. Booked before the run.

THE CHANGE (branch-local module cl10/supercell_2i.py, the published
supercell untouched): the internal layer of the step stops being the
45-plane product at a single angle and BECOMES A BINARY-ICOSAHEDRAL
ELEMENT: U_2I = exp(-i omega n.S) with S = (T23/2, -T13/2, T12/2) the
Spin(3) triple of the three SPATIAL planes of the Clifford generators
(pairwise anticommuting, S_a^2 = I/4 - verified in T0), and n the
icosahedral 5-fold axis n = (0, 1, phi)/sqrt(1+phi^2) in the cubic
frame - the first time the digit 5 sits INSIDE the engine's internal
algebra. Rodrigues closes exactly: (n.S)^2 = I/4, so U_2I is a true
order-10 element (U_2I^5 = -I booked). The bond layer keeps the cubic
stage algebra (the lattice is cubic; A5 has no transitive action on 8
points - 60/8 is not an integer - the cell itself cannot carry the
icosahedral group, booked as the structural reason the change is
internal-only). bloch(q, U) takes any internal unitary, so the step is
built with the PUBLISHED machinery, and bond_g10 decouples the hop
angle beta (V239 two-dial family).

GRID (booked): internal omega in {pi/5, 2pi/5, 4pi/5} (the 2I
angle family: orders 20/10/10 conjugates) x bond beta in
{pi/4, pi/5, 2pi/5} (published transport vs 5-fold transport).

T0 GATES (must all pass before any verdict is read):
T0a published-step anchor: with U = internal(pi/4), beta = pi/4 the
    module reproduces supercell's own numbers at GAMMA/M/R: exactly
    gapless M points (<= 1e-9), exactly 3-fold R multiplets, mirror
    pairing <= 1e-12, M(0)-order closure of the published step;
T0b S-triple certification: S_a^2 = I/4, anticommutators zero, the
    adjoint of U_2I on the spatial generator span an SO(3) matrix with
    trace = phi exactly (<= 1e-10) - the golden ratio ENTERING the
    engine; U_2I^5 + I <= 1e-10 for omega = 2pi/5;
T0c backend invariance (GPU 1 booked, user directive "utilise gpu 1"):
    every observable on the cupy path with CUDA_VISIBLE_DEVICES=1
    against the numpy host path, worst <= 1e-11; the sweep and the
    closure-scan matmuls are issued through cl10.engine routing.

MEASUREMENTS per grid point (the tribunal of the engine change):
(1) step unitarity max |M M^dag - I| at GAMMA, X, M, R;
(2) closure order at GAMMA: least k <= 240 with M^k = +-I (tol 1e-8),
    or NO-CLOSE-TO-240;
(3) M-point gap min |eigenphase| (is the M-node still exactly gapless
    when the internal layer breaks cubic equivariance?);
(4) R-point degeneracy: cluster of 256 eigenphases at 1e-9 - does the
    exact 3-fold multiplet survive;
(5) +- mirror pairing residual of the phase set (sorted s vs -s);
(6) 0-point sum of phases at GAMMA; det M(0).
BRANCH WORDS (booked): ENGINE2I-STRUCTURE-HOLDS (unitarity <= 1e-11 at
every q, gapless M preserved <= 1e-9, R 3-fold preserved, mirror
<= 1e-10 at every grid point); ENGINE2I-STRUCTURE-BREAKS (first
failure listed per observable in data.failures - the 8 was load-
bearing); PARTIAL (holds at some grid points, breaks at others; the
per-point table is the finding).
Cost booked: 9 configurations x (4 momenta + M/R stars + closure scan
to 240 at 256 dims) ~ minutes on GPU 1 (256-dim host path measured
faster for single matmuls; the power scan and eigens go through the
engine routing).
