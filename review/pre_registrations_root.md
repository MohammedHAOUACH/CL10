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

## V336 run note (post-run; verdict as booked)
ROOT-2I-E8-STRUCTURE-EXACT, all seven checks green, wall ~30 s host
(exact Fraction/Bareiss + 50-dp mpmath: no GPU path exists in this
experiment; GPU 1 is engaged by V337 per the user directive).
Run hygiene booked: (1) the first run's generator pair (5-fold about z +
3-fold about (1,1,1)) generates an INFINITE group - z is not a 5-fold axis
of the (0,+-1,+-phi) icosahedron - the exact BFS never saturated and the
run was cut at the 300 s cap with no artifact; corrected to the face/vertex
pair (q_f = (1+i+j+k)/2, q_v = rotation 72 deg about the edge-adjacent
vertex (0,1,phi)), float64 dedup closure with a 500-element infinite-group
guard: count exactly 120, order multiset exact, unit drift 1.1e-16; the
K-identity stays exact. (2) The hardcoded E8 simple system followed a
wrong fork convention (it spans D8 + A8 - the first green-ish run showed
n_roots 184 because all_roots enumerated e_i - e_j only (56) not +-e_i
+-e_j (112), and C^30 != I): both fixed by DERIVING the simple system
from the 240 enumerated roots (generic positivity + indecomposability)
- Cartan det 1, 240 roots norm^2 2, Coxeter order EXACTLY 30 (C^30 = I,
C^k != I for all proper divisors), marks (3,4,6,5,4,3,2,2) -> dimension
vector sorted [1,2,2,3,3,4,4,5,6], sum d^2 = 120 = |2I| EXACTLY, A_aff d
= 0 and the McKay 2 d_i = N d recursion at every vertex: the 2I/E8
correspondence certified, not quoted. (3) The booked T3 expectation "4
radii" was wrong (a D8 guess): the E8 Coxeter plane holds 8 rings x 30,
angle set closed under +2pi/30, and the smallest ring ratio is NAMED
phi exactly (PSLQ x^2 - x - 1, rho 1.61803398875): the golden ratio is a
radii ratio of the E8 Coxeter projection. (4) highest-root coordinates:
( alpha_i, theta ) are Cartan-weight coordinates, not simple-basis marks;
solved through the exact Cartan inverse (the booked marks are the simple-
basis coefficients). No class change, no target change.

## V337 run note (post-run; verdict as booked)
ENGINE2I-STRUCTURE-BREAKS, all three T0 gates green (published anchor:
M-gapless 0.0e+00, internal mirror 8.9e-16, unitarity 7.2e-15, GAMMA
closure order 2; T0b: S-triple exact, U_2I^5 = -I at 1.7e-16, adjoint
SO(3) with trace = phi at 0.0e+00 - the golden ratio enters the engine;
T0c: cupy GPU 1 vs host numpy, eigenphases 1.33e-15, 39-step power scan
1.44e-15 - verdicts path-invariant), wall ~10 s, grid 3 internal angles
x 3 hop angles. Gate hygiene booked: (1) my R-point check demanded
cluster sizes divisible by 3 - false for the coupled 256-dim Bloch
spectrum (the paper's 3-fold statement is about the multiplet structure
of the unfolded description; the coupled anchor clusters are [32, 48]):
the grid now compares the CLUSTER-SIZE SIGNATURE to the anchor. (2) My
q -> -q pairing at fixed generic q is not a theorem of the published
step (residual 1.0 - the chirality pairing lives at the M stars, the
+- pairing at the on-site layer): demoted to data; the internal-layer
mirror is the gated observable. (3) bloch returns host matrices (house
convention at 128/256 dims; the GPU is engaged by the T0c eigens and the
power scan). Physics, as measured: with 2I internal elements the step
stays EXACTLY unitary (<= 1.4e-15), keeps the on-site +- pairing
(1e-16), det = 1 exactly (0.0e+00 phase), 0-point at GAMMA (<= 1.4e-14),
and CLOSURES RICHER than the published step: order 10 at cubic hops
(published 2), order 5 at 5-fold hops - the 5 lives in the algebra
without breaking unitarity anywhere. What breaks: the M-point exact
gaplessness (the massless doublers, the chiral content) acquires masses
0.105 to 6.28 rad at 7 of 9 grid points (the two gap-preserving points,
(w, b) = (0.4 pi, 0.4 pi) and (0.8 pi, 0.2 pi), are accidental zeros,
not a branch); the R cluster signature [32, 48] survives only with the
cubic hop (beta = pi/4), [16, 48] with 5-fold hops. Verdict: the digit 5
is compatible with the metric half of the automaton (unitarity, pairing,
det, closure) and incompatible with its fermion half (massless nodes)
INSIDE THE CUBIC CELL: the 8 and the 5 cannot share this cell; the root
replacement must change the cell, not its angles. A5 has no transitive
action on 8 points (60/8 not an integer) - the booked structural reason
the cell cannot carry the icosahedral group.

## V338 run note (post-run; verdict as booked)
ROOT-NUMBERS-UNNAMED, wall 14 s host (mpmath scalars; no GPU path),
V335 guards verbatim. D1 (pi-free structural monoid 2^a 3^b 5^c 7^d
phi^f sqrt3^j, |exponents| bounded): best margins g_obs 1.006e-5 (this
is the booked 1/5 relative 10^-5), g_cross 8.284e-5, COEF_OBS 1.449e-5,
LOOP 1.125e-4 - all seven-plus orders above the 1e-12 naming gate; no
D1 name. Continued fractions q <= 1e6 (exact evaluation): g_obs best
19884/99421 at 4.950e-12 (reproducing V335's leaf exactly), g_cross
164414/821523 at 1.707e-13 (below the gate - a coincidental convergent,
booked O(1)-per-number density at q <= 1e6: naming g_cross by
164414/821523 carries no selection content, and the K1 band logic of
V335 still holds: the crossing number is not structural), COEF_OBS and
LOOP at 1.4e-12 / 2.4e-12. PSLQ (degree <= 6, H <= 20, 60 dp with the
1e-30 substitution check): ZERO hits on all four targets. Booked
Coxeter probe: |30 g_obs - 6| = 6.035e-5 - the clock is NOT 6/30 to any
structural precision. n* = 208 additive census over the structural
integers: 1704 hits among the booked pattern space - a representation
208 = 240 - 32 exists but proves nothing at that density (guard booked
before the run). LOOP at D0: named trivially by (2^4 pi)^-1 = 1/(16 pi)
- the booked non-discovery (geometric tier admits 4 pi and any multiple;
the ladder base is the Einstein-Hilbert coefficient 1/(16 pi G) at G = 1,
noted, not counted as a birth).

## Root-rebuild first-candidate verdict (N1 = 2I/E8), all three tests run
STRUCTURE: TRUE (V336 - the 8 and the 5 are one datum: McKay 2I = affine
E8 certified by computed marks; |W(E8)| = 2^14 3^5 5^2 7; phi is a Coxeter-
plane radii ratio). ENGINE INSIDE THE OLD CELL: NO (V337 - the 5 breaks
the massless nodes of the cubic cell; closure/unity survive). NUMBERS:
NO (V338 - none of g_obs, g_cross, COEF_OBS, LOOP is named by the root's
structural arithmetic at 1e-12). N1 is adjudicated as a certified
skeleton resonance, not the new theory. The open branch (V334's
transcendental-coefficient exit) stays open and is untested: it requires
a cell that is not cubic - the icosahedral candidate cell cannot be an
8-point cell (A5 has no index-8 subgroup) - which the programme books as
the shape of candidate N2: an icosahedral COMBINATORICS on more
sub-sites (A5 acts on 12 or 20 points: the dodecahedron/icosahedron
vertices) with the clock born as a period of the counting.

## V339_icosahedral_cell (tier 2, hard; root-rebuild; candidate N2a - THE CELL ITSELF)

The N2 shape booked after N1: A5 has no 8-point action, but it acts
transitively on the 12 vertices of the icosahedron. Build the 12-site cell
and judge it. All arithmetic in Q(sqrt5) and its quadratic extension
Q(sqrt5, t), t^2 = 10 + 2 sqrt5 (sin 72 = t/4) - exact, Fractions.

BOOKED CHECKS:
T1 icosahedral graph EXACT: the 12 vertices (0,+-1,+-phi) and cyclic with
    |v|^2 = phi+2 (= the V218 star norm 1.902113^2); edge iff dot = +-phi
    exactly in Q(sqrt5); booked counts: 30 edges, degree 5;
T2 A5 EXACT as vertex permutations: Rodrigues rotations R5 (72 deg about
    the vertex axis (0,1,phi)) and R3 (120 deg about a face normal) built
    over K = Q(sqrt5,t), R applied to each vertex compared EXACTLY to a
    vertex (a permutation), BFS closure of the two permutations: order
    exactly 60, orbit of any vertex = 12 (transitive), generator orders 5
    and 3, and the triangle product R5 R3 has booked order 2 (the
    (2,3,5) presentation of A5);
T3 clock-birth identities EXACT (integer combinatorics, zero fitted
    constants): |V| / |A5| = 12/60 = 1/5 (the desert rational born as a
    vertex-group ratio - booked as an identity to be certified, NOT as
    the clock itself: the tribunal target g_obs = 0.19999798835751165 is
    1/5 - 1.006e-5 and stays the open deviation); the 1-factorization of
    the 30-edge graph into EXACTLY 5 perfect matchings (search +
    certificate); stage S = cos g I - i sin g P with P^2 = I closes at
    order 5 at g = 2 pi / 5 (S^k = cos kg - i sin kg P exact algebra,
    float64 certificate 1e-12) - the closure clock of the new cell is
    2 pi/10 = pi/5, the Coxeter number 10 of H3;
T4 the automaton: on-site C^12 x internal 36 passed through unchanged is
    NOT run at first (cost); first the 12-state cell alone: at momentum q,
    each edge in direction e_i carries the phase exp(+- i q_i) on its
    2x2 block (12 axis-parallel edges - exactly 4 per axis, booked),
    oblique edges phase 1; stage_m(q) = cos g I - i sin g M_m(q) with
    M_m Hermitian involutions: unitarity of every stage and of the
    5-stage product U(q) = prod S_m(q) is EXACT (certified <= 1e-13 on a
    q grid); the +- mirror is the q -> -q conjugation M(-q) = M(q)^*
    (certified at grid points);
T5 the fermion tribunal: gap(q) = min |eigenphase of U(q)| over a 32^3
    grid in [-pi,pi]^3 at the closure clock g = 2 pi/5 (and the g family
    pi/5, 2pi/5, 4pi/5, pi/4 for reference): GAPLESS-NODES-FOUND if the
    grid min <= 1e-8 (count the node family) else the min-gap value is
    reported per clock;
T6 closure order of U(GAMMA) per clock (least k <= 120 with U^k = +-I).
BRANCH WORDS: N2A-CELL-BUILT (T1-T4 all green) +
GAPLESS-NODES-FOUND / NO-GAPLESS-NODES (T5), + CLOCK-BIRTH-PARTIAL booked
regardless (12/60 = 1/5 exact but g_obs unnamed: the 1.006e-5 deviation
is still not combinatorics - that deviation IS the open theory).
Cost booked: seconds (12x12 eigens; 32^3 grid ~2 minutes).

## V339 run note (post-run; verdict as booked with booked corrections)
N2A-CELL-BUILT-GAPLESS-NODES, wall 5.3 s host (12x12 eigens, 32^3 grid;
no GPU path), all gates green. Run hygiene booked, in order found:
(1) the first vertex construction used (b*phi, a, 0) instead of the
cyclic (b*phi, 0, a) - the T1 gate (30 edges) caught it immediately:
the wrong set is not the icosahedron (74-degree angles appear);
(2) the 120-degree face rotation does NOT close over Q(sqrt5,t) (needs
sqrt(10)): replaced by the (2,3,5) presentation - R5 vertex rotation and
R2 edge rotation, with the exact identity sin 72 = |v|/2 (t^2/(phi+2) =
4 in Q(sqrt5)) putting every rot5 entry in the booked tower; the order-3
element is then the product R5 R2 (measured 3);
(3) booked "4 axis-parallel edges per axis" is wrong: the golden
rectangles' LONG sides (length 2 phi) are icosahedral diagonals, not
edges: 2 per axis, 6 total - the geometry is right, the booking was not;
(4) the booked q -> -q spectral mirror is FALSE for this automaton:
conjugation flips the sign of the -i sin g term, so S_m(-q) = S_m(q)^T
and U(-q) is the reversed-order product; an ordered five-stage cycle is
not isospectral to its reverse (measured residual 1.2). The ordered
icosahedral cycle is CHIRAL BY COMPOSITION - measurement promoted to
data with the reason booked;
(5) T6 booked as an existence gate but found None: U(GAMMA) has infinite
order at the stage-closure clock - the stages close (S^5 = I at 2.4e-16)
while the COMPOSITION does not: the clock lives at the stage level, the
dynamics is quasi-periodic. Relaxed to a measurement with the reason.
Booked physics measured:
1. THE CELL EXISTS AND IS EXACT: icosahedron over Q(sqrt5) (30 edges,
degree 5), A5 over the tower as exact permutations ((2,3,5) orders
5/2/3, group exactly 60, transitive on 12, stabilizer 5), the 30-edge
graph factors into exactly 5 perfect matchings of 6 edges - the closed
five-stage icosahedral cycle - each stage an involution-permutation
S = cos g - i sin g P closing at order 5 at g = 2 pi/5 (Coxeter h(H3) =
10): the clock 1/5 born TWICE in the cell combinatorics: as the group
ratio |V|/|A5| = 12/60 = 1/5 EXACTLY, and as the order-5 stage closure.
2. THE NUMBER IS STILL NOT BORN: the tribunal target g_obs = 1/5 -
1.006e-5: the 1/5 backbone is named exactly, the 1.006e-5 deviation is
NOT (V338 margins unchanged); the closed CLOCK ANGLE is 2 pi/5 = 1.2566
radians, the desert reads 0.1999980 radians: pi does not appear in the
desert's number - the missing mechanism is still open.
3. THE FERMION: at g = 2 pi/5 (and 4 pi/5) the zone has EXACTLY ONE
gapless point in 32^3 grid points (2.3e-16 at GAMMA, isolated - one
zero of 32768), nontrivial sector gapped (min 4.2e-5); dispersion is
QUADRATIC (z = 1.99-2.05 along x, y, z, golden and body rays,
anisotropic) - a band touching, NOT the chiral linear node of the
published engine; at the pi-fraction clocks (pi/5, pi/4) the cell is
gapped (9.2e-2, 1.0e-2): the gaplessness is a FIVE-CLOCK resonance -
the new cell's masslessness exists exactly where the digit 5 lives.
4. CHIRALITY FREE: no q -> -q mirror in the ordered cycle (T4.4): the
icosahedral cell composes chiral; a candidate source of the coupling
asymmetry if the cell is ever wired to the internal module.
N2a verdict against the tribunal: cell BUILT, backbone named, numbers
NOT named, fermion NOT chiral (z = 2), composition NOT closed. Not the
elected theory. Next booked step N2b: couple the 12-site cell to the
36-dimensional internal module (432-dim cell) and ask whether the
internal dressing lifts z = 2 to z = 1 (chiral content restored) - the
published engine's massless nodes were linear BECAUSE of the internal
module's Wilson-like structure; if N2b restores linearity AND keeps the
single GAMMA node, the desert's 1.006e-5 deviation becomes the sole
missing datum of the reconstruction.

## V340_number_web_and_pisot_shadow (tier 2, hard; root-rebuild; the last
birth channel no census has closed: RELATIONS and SELF-SIMILARITY)

V338 named single numbers and failed; a theory can still be RELATIONAL
(numbers born as edges of an exact web) or SELF-SIMILAR (the desert's
2 ppm is a Pisot-power shadow of a substitution counting - the only
known mechanism giving EXACT rational backbone + exponentially small
irrational correction, which is precisely the shape of the desert data:
g_obs = 1/5 - 2.0116424883642736e-06).

TARGETS (published decimals, harvested from the artifacts, 15 of them):
g_obs, g_cross_rung10, gap = g_cross - g_obs, COEF 2.07656411117383,
coef_at_fifth 2.088106865906787, C_CORR 2.07644360080448,
LOOP 0.019894367886486918, F208_root 0.0196851541811279,
F_208(321) 0.0196927945051569, F_4 1.82994979751295, F_8 1.67435813070784,
r_atm 1.91911543798088, W_eff_ln rung10/17/208 (-3.50320080595748,
-2.79531879611527, -0.421699620700671), deficit(pi/2) 0.973756015764892,
m_gap_min 2.00289529530561e-4.
W1 pair web: PSLQ [1,x,y] H <= 2e4, residual <= 1e-13 (guard: 136 pairs
x 4e4 lattice points x 1e-13 ~ 5e-5 expected false positives).
W2 log web: PSLQ [1,ln|x|,ln|y|], accepted only with |integer| <= 2 and
residual <= 1e-13 (multiplicative web x^a y^b = e^c with c in {-2..2}).
W3 product web: PSLQ [1,x,y,x*y] H <= 200, residual <= 1e-16.
P Pisot shadow: targets T in {g_obs-1/5 gap, gap, COEF-2, coef_at_fifth
- COEF, m_gap_min}: search T = (p/q) * lambda^(-n) EXACTLY: lambda the
8 booked Pisot units (golden ratio, silver 1+sqrt2, plastic x^3=x+1,
supergolden x^3=x^2+1, smallest quartic Pisot x^4=x+1, 2cos(2pi/7)'s
Pisot conjugate? booked as x^3-x^2-2x+1, tribonacci x^3=x^2+x+1,
root x^3-2x^2+1? booked list fixed before the run), |p| <= 10^4,
1 <= q <= 100, 1 <= n <= 60, accepted at 1e-12 absolute after scaling
(guard booked: candidate count ~ 8 x 60 x 10^6/... restricted to the
window |T| lambda^n in [1e-2, 1e2]; expected false positives < 1e-3).
BRANCH WORDS: WEB-EDGE-FOUND / WEB-CLOSED; PISOT-SHADOW-FOUND /
PISOT-SHADOW-VOID. Definitionally-known edges (gap = g_cross - g_obs,
coef_at_fifth vs C_15 agreement) excluded from the claim by name.
Cost booked: ~2 minutes host.

## V340 run note (post-run)
WEB-CLOSED-PISOT-VOID, wall 47 s host. W1/W2/W3: ZERO edges among the
17 published decimals at the guarded heights (even the pair tier finds
nothing, and the booked run-note correction: a definitional TRIPLE like
gap = g_cross - g_obs is invisible to two-number PSLQ by construction,
so the triple tier PSLQ [1,x,y,z] H <= 1e3, residual <= 1e-18 (false-
positive expectation 5e-6 over 680 triples) was run as the booked
supplement: also ZERO exact triples - the desert web is CLOSED: the
published numbers are pairwise AND triply algebraically unrelated over
the guarded lattices. Pisot shadow: 5 targets x 8 booked Pisot units,
170 windowed (lambda,n) scales x q<=100 convergents: ZERO shadows at
1e-12: the 2 ppm of the desert is NOT the n-th shadow of any booked
self-similar substitution. With V338 (naming) and V335 (rational
density), every STATIC birth channel of the desert numbers is now
bounded closed: naming, relations, webs, self-similarity. What remains
open is only DYNAMICAL birth: a dressing series evaluated by the
mechanism itself - booked as N2b/N3.

## V341_icosahedral_internal_dressing (N2b; booked then run; tier 2 hard)
The V339 cell carries a QUADRATIC touching at GAMMA; the published
engine's nodes are linear thanks to the internal layer, and
supercell.bond_product states the exact factorization
bloch = kron(P(q; beta), U(g)**n_stages). N2b: the icosahedral cell
dressed by the theory's OWN internal Spin(10) unitary (16-sector, one
application per stage: exponent 5), at clocks 2pi/5, 4pi/5, pi/4.
Gates: T0c cell rebuilt (30 edges, 5 matchings, 6 axis edges); T0a
published anchor: bloch == kron(P, U^6) to 1e-10 AND published M-point
gapless to 1e-9; T2 cell-only GAMMA anchor at 2pi/5 reproduced to 1e-9.
Measurement: dressed gap on 24^3 grid (min + zeros), refined node, z
exponent along 5 rays (1e-3 vs 3e-2). BRANCH WORDS: N2B-LINEAR-CHIRAL-
NODES (isolated zeros, refined <= 1e-9, some z in [0.8,1.2]) /
N2B-QUADRATIC-TOUCHINGS (zeros, all z > 1.5) / N2B-GAPPED / N2B-MIXED.
Booked physics guard: the dressing can only shift WHICH pair (cell
phase, internal phase) touches at a given q; linearity is decided by
the cell's phase surface crossing the internal phase plane at nonzero
gradient - a cubic-vs-icosahedral structural fact, not a tuning.
Cost booked: ~30 s.

## V341 run note + V342 node census (booked then run)
V341 measured N2B-GAPPED at grid-refine precision 2.7e-9 with five-ray
z = 1.00 +- 0.008 at 2pi/5 and 4pi/5, and an EXACT quadratic touching
(z = 2.000, gap 0.0 at GAMMA) at pi/4: the refinement precision floor
was too coarse to decide existence. Correction booked for V342: polish
the contacts to 1e-16 (Nelder-Mead): the 2pi/5 contact IS an exact zero
(gap 0.0 at q* = (-2.61758, 0.26195, -2.35595)) - V341's 2.7e-9 was a
refine-floor artifact, verdict FIRST amended to N2B-GAPLESS-LINEAR, then SUPERSEDED by the V342 census: the zeros are 557 independent plane crossings with slopes 0.02-0.35 - NODAL SURFACES (a dressed metal), not chiral point nodes; V341 z=1 was a plane crossing its zero line, and V341 status word N2B-GAPPED-> see V342.
V342_icosahedral_node_census: count the dressed Dirac nodes (grid
cluster + polish + dedupe mod reciprocal/reflection) at 2pi/5 and
4pi/5, and measure the cone anisotropy (velocity along 7 rays per
node). Booked comparisons: the published engine's fermion content is
the 12 DOUBLER CHANNELS at M at the closed cubic clock; the census
word N2B-NODE-ORBIT-<n> reports n and whether the nodes are a single
reflection orbit; isotropy word CONE-ISOTROPIC (spread < 5 %) /
CONE-ANISOTROPIC. Cost: ~60 s.

## V342 run note (post-run)
N2B-NODE-ORBIT-557-CONE-ANISOTROPIC (fired; the booked comparison loses
badly): 557 distinct mod-reflection zeros at 2pi/5 (the un-deduped
inline pass found 1104 - grid-seed splitting, same verdict), velocities
0.023-0.354 along 7 rays at the first node: the dressed spectrum is 12
cell phases + 16 flat internal planes; each pair (cell branch, internal
phase) closes on a CODEIMENSION-1 SURFACE in q; the union is gapless on
surfaces = a nodal-surface METAL. The published engine's massless
content is 12 ISOLATED point doublers at M, protected by the cubic
cell. VERDICT: the icosahedral cell FAILS the chiral test - not gapped
and not Dirac, but worse: doubling catastrophe x192 surfaces. V341
status amended to N2B-NODAL-SURFACES. The two cell facts now pin the
obstruction exactly: cubic cell + cubic bond clock -> correct isolated
linear doublers (published, certified pi/4); icosahedral cell -> clock
5 born twice in combinatorics (V339: 12/60 = 1/5 EXACT as dimensionless
number AND order-5 stage closure) but no protected point nodes. The
elected assembly must put the CUBIC cell OUTSIDE (nodes come from the
bond product P(q) exactly, and are immune to ANY internal unitary -
supercell.bond_product: bloch = kron(P, U^k)) and the ICOSAHEDRON
INSIDE (as the internal module replacing Spin(10)). No prior experiment
tested that: V337 varied the bond angles (what the V342 lesson says to
NEVER touch); the two-dial family V239 already decouples bond beta
from internal g.

## V343_internal_icosahedral_contacts (N3 spectral-flow probe; booked then run)
The N3 assembly (booked by the V342 lesson): cubic bond product P(q)
OUTSIDE - its isolated linear M doublers are exact for ANY internal
unitary (kron factorization, V239 two-dial family already decouples
bond beta from internal g) - and the ICOSAHEDRAL 12-site module inside
at q_icos = GAMMA (the V339 cell operator, five stages, exponent 6 as
in the published step). This probe tests ONLY the clock: the spectral
flow of the internal module U_icos(g)^6, g in (0, 1): contacts = g*
where two eigenphases collide (mod 2π) or a phase hits 0. T0 anchor:
V339 single-zero leaf re-verified (U_icos(2pi/5) has eigenphase 0 at
2.3e-16). Measurement: scan step 1e-3, polish every local minimum of
the pairwise phase gap to 1e-12, report ALL contacts on (0,1) and
their distances to the desert data {g_obs, g_cross, 1/5}. Expected
structure: phases are piecewise linear-in-g with integer slopes at
small g (product of e^{+-ig} channels), so collisions sit at rational
multiples of 2pi/D - the question is whether one lands on 0.2. BRANCH
WORDS: N3-CONTACT-ONE-FIFTH (a contact at |g* - 1/5| <= 1e-12) /
N3-CONTACT-DESERT (a contact at |g* - g_obs| <= 1e-9 or
|g* - g_cross| <= 1e-9) / N3-NO-CONTACT-NEAR-THE-DESERT (all contacts
listed, distances reported). Cost: ~30 s. If NO-CONTACT-NEAR: the
icosahedron-as-internal clock channel closes too, and the only channel
left open by the whole hunt is a DRESSING SERIES around a bare contact
(booked consequence, N4).

## V343 run note (post-run; verdict as booked)
N3-NO-CONTACT-NEAR-THE-DESERT. T0 anchor green (U_icos(2pi/5) eigen-
phase 0 at 2.3e-16, V339 leaf re-verified). The internalized
icosahedral module U_icos(g)^6 has exactly TWO phase contacts on
(0,1): g* = 0.6539031701784 and g* = 0.8008741945100 - nearest desert
margins 0.454 and 0.601 (the second contact sits 8.7e-4 from 4/5,
six-plus orders outside any structural guard). N3's clock channel
closes: putting the icosahedron inside the cubic engine moves the
contacts AWAY from the desert, not toward it.

## The elected theory, as a theorem-shaped portrait (booked consequence)
The hunt now BOUNDS the elected theory from outside on five fronts,
each a fired experiment, none a guess:
(1) backbone: any elected theory NAMES 1/5 as exact combinatorics
    (V339: 12/60 = 1/5 in the A5 cell; the naming is real and the
    desert reads 1/5 to 5.03e-5 RELATIVE);
(2) fermions: it must carry ISOLATED LINEAR doublers - the cubic bond
    product is the only known carrier (published 12 M-channels; the
    icosahedral alternative is a nodal-surface metal, V342);
(3) it cannot share the cubic cell with any 5-fold bond geometry
    (V337 kills the massless nodes there);
(4) the 2 ppm deviation is NOT static: naming (V338), webs (V340),
    triples, self-similarity/Pisot (V340) are all bounded closed -
    the ONLY surviving birth channel is a DYNAMICAL dressing series
    or an exponential factor;
(5) SHAPE OF THAT FACTOR, from published data: V321's resummation
    exponents m = 11.29, 11.58, 10.86, 10.83 (rung 2/6/10 3-pt fits):
    -ln(delta) = 11.51703 lies INSIDE the published m-window
    (10.6-11.6): the desert deviation has exactly the SHAPE of an
    e^{-S} instanton factor of the ladder's own resummation action,
    though no published m matches at the 1e-12 guard (nearest 24 %).
N4 (next booked experiment): compute the ladder's exact resummation
action S (not the 3-point fit m) at the desert rung and test
delta = c e^{-S} at guard 1e-12. If it fires, the desert is born:
g_obs = 1/5 (named, combinatorial) x (1 - c e^{-S(ladder)}) (dynamical)
- backbone + instanton dressing - the only shape left standing.

## V344_turn_unit_clock_census (N4 core; booked then run)
CORRECTION booked first: the N4 booking leaned on V321's m-window
(10.6-11.6) as "published action"; the artifact shows RESUM-NO-TIE with
fp = None for EVERY n (c_spread ~ 0.99 at n=10: the power divergence
is not pure-power canonical). The m values are slope estimates, NOT a
canonical channel: portrait claim (5) WEAKENED to "shape suggestive,
channel not canonical" and claim (-ln delta in the window) demoted to
decoration. N4 core rebuilt around a HARD theorem-shaped fact:
Lindemann-Weierstrass: for algebraic g != 0, sin g, cos g are
transcendental: NO machinery polynomial in trig of g-in-radians can
resonate exactly at g = 1/5; the only way a RATIONAL clock can be a
root clock is if the machinery's phase is e^{2 pi i g} (turns). V209
measured the ladder's root-clock set as {8g in 2 pi Z} - denominator
8 in turns. The desert's denominator is 5. PROBE (published exact
machinery mp_ladder, no new object): does ANY q=5 clock resonate:
g = 2 pi p/q, q in {3,4,5,6,8,10,12} (the turn-unit rational clocks,
where trig IS algebraic: 72-degree trig lives in Q(phi) - the 5-clock
is NOT forbidden, only the RADIAN 1/5 is), plus the exact rational
g = 1/5 rad as the transcendence control (must NOT resonate if the
ladder is trig-polynomial). T0: V209 detuned leaves replayed exactly
(nmax 12, 55 dps, rel < 1e-12) and pi/4 control slope ~ 11.3 recovered.
Gates: per angle, |F_2| at delta in {1e-1, 1e-3, 1e-5}: divergence
flag iff |F_2(1e-5)| / |F_2(1e-1)| > 1e8 (three decades of growth like
the pi/4 control: 7.95e6 -> 1.15e52). BRANCH WORDS: TURN-CLOCK-5 (any
q in {5,10} angle diverges: the desert's denominator 5 is a root clock
of the published ladder in turn units - backbone birth OPENS) /
TURN-CLOCK-VOID (divergences only at q | 8: the ladder's clock set
stays denominator-8, the desert backbone is NOT born in the ladder,
and the dressing channel stands candidateless). Cost booked: under 15
min (nmax 10-12, 55 dps, ~25 angles).

## V345_dodecahedral_closed_sphere (user's cell proposal; booked then run)
User proposal, verbatim intent: 12 dodecahedron FACES as cells, each
face has exactly 5 neighbors, the surface is CLOSED (sphere: a signal
goes around with no wall, no boundary, no torus trick). Honest
placement before the run: the face-adjacency graph of the dodecahedron
IS the icosahedral graph certified in V339 (faces of the dodecahedron
= vertices of the icosahedron); what is NEW and never tested is the
DYNAMICS WITH NO LATTICE MOMENTUM: the closed-sphere automaton = the
V339 cell operator at q = GAMMA (all bond phases 1), read now as a
finite closed surface, plus exact combinatorial transport (no quantum
stages): the signal's return time as a permutation fact.
BOOKED CHECKS (constructed from the DODECAHEDRON directly, not reusing
V339 coordinates): face normals n_f = the 12 vectors (0,+-1,+-phi) and
cyclic (dodecahedron faces have icosahedron vertices as normals); two
faces adjacent iff they share an edge - booked criterion: the normals
satisfy n_i . n_j = phi/(phi+2) EXACTLY (same dot as the dual edge),
and independently the shared-edge test through dodecahedron vertices:
the 20 vertices v = intersections: three faces meet at each vertex -
booked triangle count of the adjacency graph = 20, and Euler
V - E + F = 12 - 30 + 20 = 2: the surface is the SPHERE, closed, no
walls: certified combinatorially (genus 0).
T1 the 5-neighbor schedule: the booked 1-factorization (V339: 5
perfect matchings of the 30 adjacencies) realizes the user's "each
cell has exactly 5 direct neighbors" operationally: at every stage a
cell talks to EXACTLY ONE neighbor, and over the five-stage cycle each
cell talks to its five neighbors EXACTLY ONCE (certified by matching
decomposition - each face appears in each matching once, and the five
partners are distinct).
T2 the signal loop (the user's "faire le tour complet sans etre
bloque"): pure permutation transport P = P5 P4 P3 P2 P1 (stage
matchings composed in booked order): exact order k and cycle
structure on the 12 faces: the return time of the signal - booked as
a MEASUREMENT (expectation open: k small means a closed tour).
T3 periodic clocks of the sphere: U(g) = prod_m (cos g I - i sin g
P_m) is the closed-sphere automaton (zero bond phases): scan g in
(0, pi) step 1e-3 for residual min_k ||U(g)^k - -+ I||, k <= 30;
polish minima to 1e-10; report every (g*, k) with residual <= 1e-8;
compare with the cubic engine's published clock pi/4 (closure order
2, V337 anchor) and with the five-clocks pi/5, 2pi/5 (V339 booked
stage closure; T6 leaf says the 5-stage PRODUCT at 2pi/5 has infinite
order - so the booked EXPECTATION is: no nontrivial finite closure at
2pi/5; the scan tests whether any other g closes the full step).
T4 closed-walk census on the sphere: W_n = tr A^n, n <= 10, plus
triangle count (20 booked) and pentagon-face count (12): the numbers
the sphere forges by counting: |E|/|V| = 30/12 = 5/2, faces/rotation
group = 12/60 = 1/5 (the V339 backbone reborn in the DODECAHEDRAL
language - faces over the dodecahedral rotation group), and
vertices/edges = 20/30 = 2/3: all exact integer ratios, zero fitted
constants.
BRANCH WORDS: DODECA-SPHERE-CERTIFIED (Euler 2 + triangles 20 +
schedule) + SIGNAL-LOOP-ORDER-<k> + SPHERE-CLOCK-<g*,k list or VOID>
(the sphere closes the full step at a nontrivial g only if the polish
says so). Honest limits booked: with NO momentum the surface carries
no band structure - it cannot host the published engine's isolated
linear doublers (V342 lesson); as an internal module on the cubic
engine the same 12-state cell was tested by V343 (no contact near the
desert). This experiment judges the proposal at its own level: closed
combinatorial transport + periodic clocks + exact ratios.
Cost booked: seconds (12x12, exact combinatorics).

## V344 run note (post-run; verdict as booked)
TURN-CLOCK-VOID, wall 480 s (nmax 10, 55 dps, 13 turn-clocks q in
{3,4,5,6,8,10,12} reduced from the booked ~25 by the p<q, gcd=1 sieve
plus the g = 1/5 radian control). T0 green to the bit: V209's detuned
F_2 leaves replay EXACTLY through this lens (rel 0.0e+00 at 1e-1 and
1e-5); pi/4 control slope ratio 1.444e+45.
THE MEASUREMENT: diverging set = {2pi*1/8, 2pi*3/8} = {pi/4, 3pi/4} -
the ODD-EIGHTHS family. Refinement booked as data: V209's closure
class 8g in 2piZ is NOT itself the degeneracy set (2pi*1/4 = pi/2
sits in the class, does NOT diverge: ratio 1.277); the true root
clocks of the published ladder are the cyclotomic 8-family of
Cl(10,0) in turn units. Every q in {5,10} clock is FINITE (72 deg:
F_2 exact 2.514, ratio 1.538; 2pi/5: 1.117; 2pi/10 and 3pi/5 share
the finite 7.178e+05 value), and the radian 1/5 control is finite
(2.756): the Lindemann-Weierstrass framing behaved exactly as
framed - it is not the algebraicity of 72 deg that fails, the ladder
genuinely does not ring there.
Corrections booked, in order: (1) first run died at the PRINT loop
(KeyError: control row stores its values under the delta keys, not
F2_1e-5); fixed, re-run, data identical. (2) first status logic
implemented the VOID branch with the wrong proxy (required divergence
AT 2pi*1/4) so the correct data fell through to TURN-CLOCK-
UNEXPECTED; the branch logic was replaced by the booked wording
(divergences only at denominators dividing 8, plus a live pi/4
control) BEFORE the final artifact run. No data changed, no class
change.
N4 CORE CONCLUSION: the desert's denominator 5 is NOT a root clock
of the published ladder in any unit (radians: forbidden by LW and
confirmed finite; turns: q in {5,10} all finite on the exact
ladder). The backbone birth channel is closed on both sides; the
dressing channel stands candidateless, and the paper's clock stays
the read 0.1999986, not a named root.

## V345 run note (post-run; verdict as booked)
DODECA-SPHERE-CERTIFIED-SIGNAL-LOOP-9-SPHERE-CLOCK-VOID, wall ~60 s.
T0 CORRECTION booked in place of the booked construction: the booked
face extraction (support planes with normals (0,+-1,+-phi)) FAILED
the run-gate - the standard dodecahedron vertex set (0,+-1/phi,+-phi)
cyclic + (+-1,+-1,+-1) is NOT dual-compatible with that naive normal
convention (a "face" picked by max-dot mixed edge lengths 2 and
2/phi; caught by the booked vertex-count/edge tests before any
verdict). Replaced BEFORE the verdict by construction from the
DODECAHEDRON'S OWN exact edge graph: edges = vertex pairs with
d^2 = 4/phi^2 = 6-2 sqrt5 EXACT in Q(sqrt5) (an earlier form 8-4phi
was the same number written wrong: 4phi = 2+2sqrt5), faces walked
from the geometric rotation system at each vertex. All booked
counts then certified EXACTLY: 20 vertices, 30 edges, vertex degree
{3}, 12 pentagon faces walked, face adjacency 30, face degree {5},
adjacency triangles 20 (= the 20 dodecahedron vertices - the dual
closed), Euler 12-30+20 = 2: the user's CLOSED SPHERE is certified
combinatorially, genus 0, no walls, no boundary, no torus trick.
T1 the booked schedule exists and is the 1-factorization: at every
one of the five stages each face cell talks to EXACTLY ONE neighbor,
and over the five-stage cycle to its FIVE neighbors exactly once.
T2 THE SIGNAL LOOP (pure permutation transport, no quantum): cycle
structure [3,9], return order k = 9 stages - the signal closes its
tour after NINE stages, not five; 3 faces are on 3-cycles, the other
9 on one 9-cycle. New integer 9: uncited, stays out of the paper.
T3 periodic clocks of the closed sphere U(g) = prod_m (cos g I - i
sin g P_m): scan (0,pi), polish, gate 1e-8: NO closure at any k <=
30: SPHERE-CLOCK-VOID. At the five-closed clock 2pi/5 the residual
is 1.09 - the V339 T6 leaf (the five-stage product has INFINITE
order there) re-forged on the sphere. T4 the sphere's exact ratios:
edges/faces = 5/2, faces/rotation group = 12/60 = 1/5 (the desert
backbone reborn in dodecahedral language), vertices/edges = 2/3,
triangles/faces = 5/3; closed walks W_2..W_10 with W_3 = 120 =
20 x 6.
HONEST VERDICT OF THE PROPOSAL: combinatorially TRUE as stated (the
face-adjacency graph of the dodecahedron is a closed sphere where
every cell has exactly 5 direct neighbors - now certified from the
dodecahedron itself, not inherited from V339) and dynamically EMPTY
at its own level: no coupling constant closes all twelve cells, the
transport loop has order 9 (not the desert's 5), and with no momenta
no band fermions live there (V342 lesson). The user's cell is the
cleanest picture of the V339 graph ever built; it is not a birth
channel for the desert clock.

## V343 permanent rebuild note (post-run; supersedes the scratch numbers)
The N3 probe ran as a scratch script; its numbers had no artifact.
Permanent script now exists: experiments/root_rebuild/
V343_internal_icosahedral_contacts/. It reproduces V339's construction
VERBATIM (same coordinate order, same deterministic 1-factorization
backtracking -> same factorization) and the booked T0c anchor to the
digit (U_icos(2pi/5) eigenphase 0 at 2.3e-16).
CORRECTION: the scratch run mixed the two booked contact channels
(phase-pair collisions and phase-hits-0) into one min-function, so
its "exactly TWO contacts on (0,1)" was an artifact of the mixed
metric. Measured SEPARATELY on U_icos(g)^6, grid 1e-3, golden polish
to width 1e-12, contact gate 1e-9 (transverse zeros floor ~1e-13,
tangential double roots floor ~1e-12 float-noise, true avoided
crossings stay >= 1e-8 away: no ambiguity band - 0 rejections in
(1e-9, 1e-8]): 115 exact contacts (collision channel 101, phase-0
channel 14). Provenance of the two booked numbers, fully resolved:
0.6539031701784 = a TANGENTIAL collision (depth floor 1.6e-12 at
0.6539031702315, margin 0.454 - the scratch gate ~1e-8 counted it,
an overly strict 1e-12 gate rejects it once; booked margin
reproduced); 0.8008741945100 = a phase-0 contact whose polished zero
sits at 0.8008741947207 (scratch polish stopped at residual 1e-10,
margin 0.601 reproduced).
VERDICT UNCHANGED AND STRENGTHENED: N3-NO-CONTACT-NEAR-THE-DESERT.
The flow is dense with contacts (the count was never the physics);
the point was desert proximity, and the strengthened census moves it:
the SMALLEST distance from any of the 115 contacts to any desert
datum {g_obs, g_cross, 1/5} is 0.00931 (contact 0.2094395102394,
phase-0) - seven orders of magnitude outside the booked 1e-9 guard.
Decoration (uncited): a phase-0 contact sits at 0.4188790204787 =
2pi/15 exactly (15 = 3x5 on the unit circle) - margin 0.219, units
mismatch, pure decoration.

# The full engine transplant: pre-registrations V346-V348 (2026-11-01)
User directive, stated twice and reinforced: the 12 dodecahedron FACES
become THE cells of the engine - the complete engine, GPU included -
this cell is the elected theory, and the user accepts losing the
whole legacy if the legacy loses: truth first, rebuild from the cell
if needed. The house answers the only serious form of an election
claim: AN ELECTED THEORY MUST REPRODUCE EVERY PUBLISHED ANCHOR. The
structural opening V337 recorded is why this program exists at all:
A5 has NO transitive action on the 8 sub-cells (60/8 not an integer),
so the 2I transplant could only live INSIDE the engine; but A5 acts
transitively on the 12 faces (60/12 = 5, stabilizer A4, certified in
V339 T2) - the user's cell is the FIRST candidate that can carry the
icosahedral group AS the cell. The transplant is built with the
PUBLISHED machinery only: altcell's bond-stage representation (a cell
= bond stages, stage cycles = perfect matchings), the published
Cl(10,0) internal 45-plane product (clock layer untouched),
cl10/engine GPU routing (GPU 1 per the standing directive).

## V346_sphere_cell_kinetic (tier 1, hard; the cell's kinematics)
The cell: 12 faces, 30 bonds (V345's exact construction: Q(sqrt5)
edge test, rotation-system faces - no normals, no inherited
coordinates), stages = the 5 perfect matchings (V345 T1). Kinetic
operator: weighted adjacency A(theta) with U(1) fluxes on the 30
bonds - the bond layer ported to the closed surface (no lattice
momentum: there is NO Brillouin zone on a sphere, booked as the
structural difference from the cubic engine, not as a choice).
T0 (gate before any reading): charpoly of A(0) over Z by exact
Bareiss elimination: booked (x-5)(x+1)^5(x^2-5)^3 - spectrum
{5^1, sqrt5^3, -1^5, -sqrt5^3}; W_3 = tr A^3 = 120 (V345 leaf) must
reproduce; if the polynomial is not EXACT everything downstream is
DRIFT-VOID.
M1 exact gap structure: triplet/singlet multiplicities {1,3,3,5}
booked; the two 3-fold levels are the candidate for the published
family number 3 (V91 count-3 forcing): the census reports whether 3
arises from the cell ALONE - reported, never claimed.
M2 flux-cone scan: does ANY bond-flux field close a gap and hold an
ISOLATED LINEAR cone (the published doublers' signature)? Families:
the A5-equivariant flux family (flux constant on the 30 bonds by
edge-orbit - the icosahedral graph has ONE edge orbit under A5, so
the equivariant family is 1-dimensional: A = A(0) with one flux -
and it is pure gauge by symmetry: booked expectation: NO equivariant
gap closing), then structured breakings: face 5-cocoloring flux
(A5/A4 coset pattern, amplitude scan), and random 30-torus sampling
at 1e5 draws (12x12 hermitian eig: host per the engine's measured
routing - the device would be slower at this size, gate cited).
For every closing <= 1e-9: local slope + anisotropy (V342 protocol):
ISOLATED LINEAR cone vs nodal crossing (the V342 chiral tribunal
applies verbatim).
BRANCH WORDS: SPHERE-KIN-SPECTRUM-EXACT + DOUBLERS-ISOLATED-LINEAR /
SPHERE-KIN-SPECTRUM-EXACT + DOUBLERS-NODAL / SPHERE-KIN-SPECTRUM-
EXACT + DOUBLERS-VOID / SPHERE-KIN-DRIFT.
Cost booked: seconds (exact) + minutes (sampling).

## V347_sphere_engine_step (tier 1, hard; the engine transplant itself)
The step: W(g) = S(g) . C, S(g) = the published 5-stage clock
automaton on the face graph (V345: stage m = cos g I - i sin g P_m,
P from the booked factorization), C = the published internal Cl(10,0)
45-plane product at the implemented-branch angle (cl10.supercell /
clifford UNTOUCHED - the transplant changes the CELL, not the clock
layer), on C^12 (x) C^d with d the published per-subcell internal
dimension - read from supercell at T0a, never assumed.
T0a published-path anchor: inside this script, BEFORE the sphere is
built, the cubic published path reproduces its own anchors (gapless
M <= 1e-9, R 3-fold multiplets, mirror pairing <= 1e-12 - the V337
T0a battery) - if that fails, DRIFT, no sphere number is read.
T0b backend invariance (GPU 1, cl10.engine routing): every sphere
observable on the cupy path vs the host path, worst <= 1e-11 (V337
T0c gate verbatim); step unitarity max |W W^dag - I| <= 1e-11.
M1 clock-contact battery (the election test): the sphere rung ladder
F_n^sph(g) = product of the smallest singular values of
(W(g)^n - I) - the finite-cell port of the published rung
periodicity-defect functional (rung n: n in 3..24 + the published
pattern rungs 170, 208): zero set polished to width 1e-12, gate
1e-9, compared against the published anchors g_obs =
0.19999798835751165 and g_cross_rung10 = 0.20013316730041675
(V206 artifact read via art(), never retyped).
M2 closure order at the published clock angles: least k <= 240 with
W^k = +-I (tol 1e-8) - the V337 measurement ported to the cell.
BRANCH WORDS: SPHERE-ENGINE-CLOCK-CONTACT (a rung zero within the
1e-9 guard of a published anchor: an elected-theory CANDIDATE book,
never an identity proof) / SPHERE-ENGINE-CLOCK-VOID / SPHERE-
ENGINE-DRIFT.
Cost booked: minutes (GPU 1; the ladder scan n<=24 x 600 grid points
at 384 dims or d x 12).

## V348_elected_anchor_battery (tier 1, hard; must the paper be rewritten?)
Every published anchor, read from artifacts (never retyped), judged
against the sphere engine: (a) g_obs, g_cross (V206); (b) the
published LOOP rung-ratio (config ladder); (c) family number 3
(V91/V346 census); (d) the 12 doublers / isolated linear M-channels
(V187/V239) - V346 M2 verdict is the answer; (e) masses/mixings:
the Yukawa machinery is defined on cubic sublattice sites - if no
port exists, the honest verdict is UNDEFINED = NOT REPRODUCED,
booked as a failure, never excused; (f) the closed-word/rigidity
tier-2 theorems are cubic-cell theorems: they do not transfer -
booked in advance as NOT-TRANSFERABLE.
BRANCH WORDS: ELECTED-ALL (every anchor reproduced: the paper IS
rewritten around the cell - booked consequence, not rhetoric) /
ELECTED-PARTIAL-<R-F-U ledger> / ELECTED-NONE.
Structural wall booked BEFORE the run: doublers REQUIRE momentum
(isolated linear crossings in a Brillouin zone); a closed sphere has
no BZ. If V346 M2 reports no isolated linear cone, the published
fermion sector is unreproducible on the cell and ELECTED-ALL is
arithmetically out of reach: the battery then measures how FAR the
cell gets, not whether it wins. No branch word pre-judges the
measurements; every number below is measured or read from an
artifact.

## V346 run note (post-run; verdict as booked)
SPHERE-KIN-SPECTRUM-EXACT + DOUBLERS-VOID, wall ~4 s. T0 EXACT:
charpoly of A(0) over Z = (x-5)(x+1)^5(x^2-5)^3, tr A^3 = 120.
M1: multiplicities {1,3,3,5} - the cell ALONE forges two 3-fold A5
representation multiplets (the family number 3 appears on the sphere
with zero fitting: REPORTED, not claimed) and an exact gap 5-sqrt5.
Lensing corrections booked (three iterations, each caught by the
script's own gate before any verdict was read): (1) the first scans
counted the trivial family members q = 0, 2pi as flux closings;
(2) the gap metric measured consecutive eigvalsh differences = the
multiplets themselves (= M1) instead of the booked DISTINCT-LEVEL
closures; (3) the booked "A5/A4 coset 5-coloring" does NOT exist as
an equivariant pattern (face stabilizer is C5, C5 not a subgroup of
A4, so no A5-equivariant map from the 12 faces to the five A4 cosets
exists - a genuine structural fact the booking missed); the honest
substitute (proper 4-coloring gradient flux) ran instead.
THE MEASUREMENT: equivariant family gap 1.24 at w=1 (pure gauge as
booked); col4 structured flux min gap 1.24 at nontrivial points;
stage flux min 1.6e-3 (no closing); 1e5 random draws on the 19-dim
gauge-fixed flux torus: min gap 8.14e-3, ZERO closings <= 1e-9.
The structural wall booked BEFORE the run (a closed surface has no
Brillouin zone, hence no isolated linear doublers) is now MEASURED,
not asserted: the published fermion sector cannot live on the cell.

## V347 run note (post-run; verdict as booked)
SPHERE-ENGINE-CLOCK-VOID, wall 430 s. THE TRANSPLANT: five
swap-rotation bond stages on the face graph (cos b I + i sin b
X_pair, X the symmetric pair exchange; the first attempt used the
antisymmetric K, K^2 = -I gives hyperbolic cosh/sinh - caught by the
unitarity gate at 1.7 before any measurement was read) times the
UNTOUCHED published Cl(10,0) internal layer (32x32 read from
supercell), 384 dims, GPU path through cl10.engine routing.
T0a green: published cubic anchors re-derived inside this script (M
gap 0.0, R cluster-size signature [32,48], internal mirror 8.9e-16;
metrics copied VERBATIM from the V337 gate including its booked
gate-defect lesson: the +- mirror is the internal layer's pairing).
T0b green: unitarity 1.3e-15; backend invariance (circular-Hausdorff
multiset metric) -3.6e-15.
T0c SPECTRAL-LENS FINDING, booked as a repo-level fact: the sphere
step at the published angle G10 carries a TRUE 10-fold degenerate
phase cluster at theta = 0.2896744895 - certified by Riesz-projection
rank (whole plane 384.000, disk r=0.15 about the cluster 18.01,
integers as contour-projector traces must be). The general QR
eigensolver splits it into a cloud of radius ~0.3, and the ENGINE'S
OWN trace-identity verification of device eigvals ACCEPTS the wrong
answer (the two trace identities are blind to intra-cluster
splitting) - booked engine finding for the next library change (no
library change made now; a full re-gate costs 3.9 h and nothing
published uses this path yet). Also booked: phase-comparison
observables on spectra with masses of exactly-0.0/-0.0 phases MUST
compare multisets on the circle (sorted pairwise comparison is a bad
observable and faked failures up to 0.29). All physics here reads
the SCHUR diagonal with certificate: worst strict-triangular norm
over every diagonalized sphere step = 0.0e+00 exactly (normal
matrices; the diagonal IS the spectrum); the QR lens is quarantined
for this experiment.
M1 THE ELECTION MEASUREMENT: 24 rung functionals (n in 3..24,170,208)
x 1001 beta-points: rung zeros within the guard of g_obs or
g_cross_rung10: ZERO; smallest rung defect AT the anchors 0.005
(g_obs) and 0.009 (g_cross) - three orders away from a born clock;
booked LOOP crossings: none at the anchors. M2: closure order at
G10: NO-CLOSE-TO-240. The clock is not born on the cell, in the
full engine, in the published machinery.

## V348 run note (post-run; verdict as booked)
ELECTED-PARTIAL-1R-3F-3U. The ledger, anchor by anchor, against the
full transplant:
REPRODUCED (1): the NAMING 1/5 - faces/rotation group = 12/60 = 1/5
EXACT on the closed sphere (V345), charpoly EXACT (V346). The
combinatorial skeleton of the paper survives on the user's cell.
FAILS (3): the published clock pair g_obs/g_cross (V347 VOID: zero
rung zeros at the anchors, defects 0.005/0.009); the LOOP^rung
ladder contacts (no crossing within the guard); the isolated linear
doublers (V346 VOID: zero flux closings in 1e5 draws - the cell
carries no Brillouin zone, no Fermi velocities, no 12 M-channels).
PARTIAL (1): family number 3 - real as A5 kinematic multiplets of
the cell (two triplets), but the V91 count-3 FORCING is an engine
theorem with no port.
UNDEFINED = NOT REPRODUCED (1): masses and mixings - the Yukawa
machinery is defined on cubic sublattice sites; no port exists;
booked as a failure, never excused.
NOT-TRANSFERABLE (1): the tier-2 cubic-cell theorems (closed words,
rigidity wall, dictionary) speak about the old cell.
READING OF THE LEDGER, stated once and without rhetoric: the
12-face closed sphere reproduces the paper's COMBINATORICS (naming
1/5, exact golden spectrum, family 3 kinematic) and fails every
DYNAMICAL anchor (clock, ladder, doublers) and the entire matter
sector. The election claim, when the whole engine is actually
transplanted to the cell and judged by the published battery, is
ELECTED-PARTIAL-1R-3F-3U: the paper is NOT rewritten around the
cell, and V335's verdict - the denominator 5 NAMES, it does not
BORN - stands, now strengthened by the transplant itself. The
legacy was not lost; it is the battery that judged the candidate,
and it will judge the next one the same way.

## V349_icosahedral_walk (the user's "new geometry" in its strongest concrete form; booked then run)
Directive taken literally: BUILD the new geometry FOR the cell instead
of testing the cell against the old one. The one dynamical organ the
battery found missing (V348: no doublers, no clock) is matter with a
velocity; the minimal new geometry that can carry it while keeping
EXACTLY the 12 cells is a quantum walk: each face is a site with a
5-way internal coin (its five pentagon edges), the A5 rotation group
IS the symmetry geometry (the 12 faces are A5/A4-cosets, V339 T2;
the 12 published M-points and the 12 faces have the same orbit size
- the coincidence the user is pointing at). Construction (all
exact, no fitting): states = the 60 directed face-adjacencies;
shift S = traverse the bond; coin C(alpha) = direct sum over faces
of exp(i alpha (P_f + P_f^dagger)), P_f = the pentagon's OWN cyclic
edge order at face f (boundary cycle of V345's rotation-system walk:
proper rotations transport it, so the coin family is A5-covariant
BY CONSTRUCTION). U(alpha) = S C(alpha), 60x60 unitary.
T0 GATES (before any reading): T0a unitarity <= 1e-12 at sample
alphas AND spectrum on the unit circle (|lambda|-1 <= 1e-12); T0b
A5-equivariance CERTIFICATE: the generators p5, p2 of V339 T2 act
as permutations of the directed-edge basis (cyclic orders
transported), and [U(alpha), P_g] = 0 for every sampled alpha, worst
<= 1e-12; T0c the walk graph IS the face graph: 12 faces, 30 bonds,
pentagon boundary cycles closed (V345 counts).
M1 GAPLESS LOCUS: scan alpha in (0, 2pi) grid 1e-3: quasi-energy
zero set (min |phase| <= 1e-9, polished to 1e-12); multiplicity per
closing; compare to the published doubler count 12.
M2 PROTECTION TRIBUNAL (the question that decides "new geometry
matter"): at each gapless alpha*, break A5 with 1e3 random
perturbations (random bond phases AND +-5% coin-angle jitter at 10
levels of amplitude): a protected closing keeps |quasi-energy|
<= 1e-6 at every level (topological), a kinematic one opens a gap;
report the protection fraction.
M3 CLOCK BATTERY (unchanged house rule): every gapless alpha*
compared against the published anchors g_obs, g_cross_rung10 (V206
artifact, art()) and the fifth, BOTH as alpha* = g and alpha* =
2 pi g, guard 1e-9.
BRANCH WORDS: QW-GAPLESS-PROTECTED-CLOCK-CONTACT (an A5-protected
gapless multiplet AT a published anchor: elected-theory CANDIDATE,
booked consequence: the paper IS rewritten around the cell - never
an identity proof) / QW-GAPLESS-PROTECTED-CLOCK-VOID (protected
matter EXISTS on the cell but not at the desert: the new geometry
is real physics, still not THE elected theory) / QW-GAPLESS-
UNPROTECTED / QW-MASSIVE (no gapless locus at all: the cell cannot
carry massless matter even with an internal coin - the structural
wall stands, now at its most generous).
Cost booked: seconds (60x60).

## V349 run note (post-run; verdict as booked)
QW-GAPLESS-UNPROTECTED, wall 12 s. The user's "build the whole new
geometry FOR the cell" in its strongest concrete form: the 12-face
cell WITH internal matter - A5-covariant quantum walk (60 directed
face-bonds; coin = exp(i alpha (P+P^dag)) on each face's OWN pentagon
cyclic order; U = S.C). Gates all green BEFORE reading: T0a unitary
7.5e-15; T0b A5-EQUIVARIANCE CERTIFICATE (V339 T2 p5, p2 lifted to
the directed basis commute with U at 6.8e-16 - covariance proved,
not assumed); T0c V345 graph counts.
M1: the walk DOES have a gapless locus - 9 exact points in
alpha in (0,2pi), 3-5 zero modes each (the internal coin is what
V346's bare adjacency walk lacked). But max modes = 5, not 12: the
"12 faces = 12 doublers" coincidence fails even numerically.
M2 PROTECTION TRIBUNAL (the decisive measurement): at amplitude
1e-4 of random A5-breaking (bond phases + coin jitter), ALL 9
gapless points open gaps (<= 160 draws, worst-gap chain booked per
point): the cell's massless modes are KINEMATIC crossings, not
topological. The published doublers' defining property - protected
by the lattice translation + point-group symmetries (the Nielsen-
Ninomiya parity the paper's V187/V244 theorems enforce) - is
precisely what a 12-cell closed geometry cannot supply: it has no
translation symmetry to protect anything, so gaplessness needs
fine-tuning at every point. The tribunal measured the structural
wall at its most generous level.
M3: nearest gapless point 0.0211 from g_obs: no clock contact.
Verdict line: even the new geometry BUILT FOR the cell, with an
internal coin supplying the missing dynamical organ, produces mass-
less matter that is fragile and off-desert. The concrete "new
geometry" family - signal (V345), engine (V347), matter-walk
(V349) - is now measured at every constructible level.

## V350_backreaction_fixed_point (tier 1, hard; the 30 % class, booked exactly before the run)
Question: is the clock an OUTPUT of a self-consistent matter<->geometry
loop instead of an input angle? The one physical hypothesis the whole
hunt never tested (V349's walk failed for lack of backreaction). All
machinery is the PUBLISHED rung ladder (V206 verbatim: step(g) over the
45 Spin(10) planes, sector_hops(G, rho) with the V85 projectors of TL,
phases exp(i g e_lambda /2), F_n = cartan_norm ratio; sectors 0..2,
nmax 10), and rho is a first-class input of the published hop
mechanics already.
BACKREACTION (booked exactly, no free knob): excitation weights
n_lambda(g) = 2 - 2 cos(g e_lambda / 2) from the published hop
exponents e_lambda; link densities sigma_lambda = Tr_int(P_lambda)/32
unit-normalized; rho(g) = sum_lambda n_lambda(g) sigma_lambda / Z,
Z = sum n_lambda Tr sigma (rho Hermitian PSD trace 1 by construction).
MAPS (the loop is then a map g -> g'):
A: g' = bisection root of |F_10(g; rho_prev)| = LOOP on (0.20, 0.21),
   start g_cross published;
B: g' = bisection root of 3 (F_10/LOOP)^10 = COEF_OBS on
   (0.1985, 0.200), start g_obs published.
Iteration: g_k = mix(g_{k-1}, g'_k), mix in {1.0, 0.5, 0.2} (0.5/0.2
only if the previous mix fails to converge), bisection 27 steps
(width 1e-10), at most 40 iterations per map, converged when
|g_k - g_{k-1}| < 1e-10; last 50 iterates all distinct -> cycle ->
NO-FIXED-POINT.
T0 GATES before reading any fixed point: T0a rho = I/45 reproduces
the V206 published anchors g_obs and g_cross inside THIS script's
verbatim bisections to <= 1e-5 vs the V206 artifact (art()); T0b the
reference branch F_10 at REF_G matches the V113 published line to
<= 1e-6; T0c rho(g) Hermiticity/PSD/trace-1 <= 1e-12 at 5 sampled g.
BATTERY (house rule unchanged): every converged g* against g_obs,
g_cross_rung10, 1/5, guards 1e-9 (and contraction factor + distance
reported regardless).
BRANCH WORDS: FP-CONTACT (a fixed point ON a published anchor within
1e-9: the clock is self-consistency - booked consequence: the clock
chapter of the paper is rewritten around V350 and the whole battery
re-run) / FP-NEAR (converged, 1e-9 < d <= 1e-4: a real mechanism
candidate, booked: rerun with corrected rho, no paper claim until
contact) / FP-DRIFT (converged self-consistent clocks exist but off
desert) / NO-FIXED-POINT (cycles/divergence: self-consistency has no
clock) / MACHINERY-DRIFT.
Cost booked: 3-4k rung evaluations, wall <= 20 min.

## V351_fold_icosahedral_star (tier 2, exact structural; booked)
V336 certified the seam at GROUP level (2I/McKay/E8: Coxeter element of
order 30 exact, the 240 roots project on the Coxeter plane to 4 radii,
angle sets closed under +2pi/30 and +2pi/10 at 1e-25). V351 asks the
CELL-level question V336 did not: is the 12 of the user's cell a
PROJECTED ROOT CONFIGURATION? Recompute the V336 fold exactly (E8
simple roots over Fraction, Cartan Gram det 1, 240 roots enumerated,
Coxeter C^30 = I exact, adjugate null vector of C - lambda I at 50 dp,
projections). BOOKED CHECKS: T0 the V336 T3 census re-derived (4 radii,
multiplicities, angle closure <= 1e-25) and radii/multiplicity numbers
equal the V336 artifact; M1 per-radius angle census: does ANY radius
class carry a 12-point star whose angle multiset equals the icosahedron
vertex set projected on its Coxeter plane (the 6-antipodal-pair dodecagon
pattern {2pi k/12} up to rotation, OR the icosahedral 2-fold projection
star with the booked golden angle differences in Z[phi]) - tested as
exact angle algebra, not fitted.
BRANCH: FOLD-ICOSA-STAR (the cell's 12 is a root configuration of the
published engine: booked consequence: the naming claim of V335/V336
upgraded from group level to cell configuration; NO clock claim is
permitted by this script) / FOLD-NO-STAR (the fold is group-level only;
the cell's 12 is not born from the 240) / FOLD-DRIFT.
Cost booked: seconds.
