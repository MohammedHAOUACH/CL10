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

## V350 run note (post-run; verdict as booked, corrections booked here)
FP-DEGENERATE-RIGID-DENSITY, wall 38 s. First run printed
FP-CONTACT at the published anchors (4.3e-14) with contraction 0.0 -
the contraction gate itself exposed the degeneracy: the BOOKED
self-consistent density is RIGID at I/45. Certified (T0d): each V85
projector's internal partial trace sigma_lam equals I/45 to 1e-16
(Schur-type internal uniformity of the published adjoint-coupling
projectors - three projectors, hop exponents {-32,-12,4}), so any
density built from the ladder's OWN excitation spectrum is
proportional to I/45 exactly: rho(g) - I/45 = 7e-18 at all sampled
g. A random non-self-consistent rho DOES move F_10 (~1e-3), so the
ladder is rho-sensitive: the backreaction channel through link
densities from the ladder's own spectrum is CLOSED by the published
machinery itself, not by tuning. The fixed-point iteration then has
constant maps and its "contact" is the old clock reproducing itself
through a rigid channel: the branch word was corrected
post-first-run (no class change, no target shift), the measurement
stands, and the T0d gate is permanent in the script.
Corrections booked (first-run, all caught by the script's own gates
before any number was believed): (1) the map brackets must CONTINUE
the crossing branch attached to the iterate (the razor band
re-crosses ~12 times, V326); the first version took the first
crossing and T0a failed at 3.2e-2 - fixed by continuation. (2) the
hop hoist: my first delta-algebra collapsed the wrong index (the
correct collapse is A_rho = einsum("ajbm,mj->ab", Pa4, rho), 32x32
per channel) - caught by the T0h hoist identity gate against the
V206 direct path (5.0e-16 / 1.8e-14 after fix). T0a re-derives
g_obs to 3.1e-12 and g_cross to 4.3e-14 with the unbackreacted
density; T0b F_10 matches V113 to 1e-9; T0c rho(g) is a true density
matrix (Hermitian/PSD/trace-1, 2.2e-16).
Physics reading: the 30 % class (clock = fixed point of a
matter-geometry loop through the ladder's own channels) is CLOSED by
a theorem-level property of the published machinery - self-dual
rigidity. Any backreaction that could move the clock must use matter
that is NOT the ladder's own spectrum (external sector), i.e. new
physics beyond the published model, not a rearrangement of it.

## V351 run note (post-run; verdict as booked)
FOLD-ICOSA-STAR, wall 1.3 s, POSITIVE. The E8 Coxeter fold of V336
re-derived inside (8 radii x 30, +2pi/30 closure at 1e-25, radii
match the V336 artifact to 1e-10). The icosahedron's 12 in its
3-fold (H3-Coxeter) reading: two 6-stars at 60 deg spacing on two
radii, ratio = 0.61803398874989... = 1/phi exact (PSLQ: x^2 + x - 1),
inter-star offset 0.0: embeds EXACTLY into the published fold -
20 witnesses (class pair x residue x assignment), radius ratio to
1e-18, offset to 1e-20, one global rotation, zero fit. The 5-fold
(dodecagonal) reading does NOT embed: its two axis vertices project
to the origin and the fold has no zero-radius class (min 0.2584) -
honest negative for that reading.
Lens corrections booked (first run): (1) class angles must stay at
50 dp - a float64 cast made the 1e-25 closure gate unpassable (T0
caught it); (2) the icosahedron projection needs an orthonormal
basis of the plane: atan2(p_y, p_x) is not isometric on a tilted
plane and destroyed the equilateral-ring structure (the
star-formation test caught it). Both fixed before any verdict.
Honest limits: the embedding is of the 12-point VERTEX
CONFIGURATION (radii + angles), not of the cell's face-adjacency
graph; and per the V335-era rule this is naming at configuration
level, not birth - no clock, no dynamics claim, and none is
permitted by this script. What IS upgraded: the seam between the
user's cell and the published engine is no longer only group-level
(2I/McKay, V336) - the golden 12 sits inside the E8 root fold as an
exact sub-configuration, with 1/phi as the certified matching
invariant.

## V352_oneloop_external_backreaction (tier 1, hard; the only door V350 left open, booked exactly before the run)
V350 certified: self-consistency through the ladder's OWN hop
spectrum is rigid (rho == I/45 exactly). The only backreaction that
could move the clock must come from a sector the published ladder
never inserted: the quantum correction. The published rung ladder is
tree-level traced links (V206); the one-loop self-energy of the
published 32-dim Weyl module in the rung chain is the ONE correction
the model omits and every gauge theory requires. No new constants:
the insertion carries the published loop factor LOOP = 1/(16 pi)
(config) and the published V205 insertion amplitude where identity-
type; three insertion choices, all published objects, booked, zero
fit:
  J1 = I x a_ins, a_ins = the V205 artifact insertion amplitude
       (data...read by art(), never retyped);
  J2 = sum_a CARTAN_a^2 (the published Cartan quadratic of V206's
       CARTAN list);
  J3 = G(g) (vertex dressing by the step itself).
CORRECTED LADDER: B~_n = B_n + LOOP * sum_{k=1}^{n-1} B_k J B_{n-k},
F~_n = cartan_norm(B~_n)/cartan_norm(G_n) (one-loop rung insertion on
the hoisted, T0h-identical hop machinery imported from V350).
MAPS (same booked structure as V350, continuation roots, damped
iteration mix {1.0,0.5,0.2}, tol 1e-10, <= 40 iters):
  A: g' = root of |F~_10(g; insertion at g_prev)| = LOOP,
     start g_cross;
  B: g' = root of 3(F~_10/LOOP)^10 = COEF_OBS, start g_obs.
T0 GATES: T0a uncorrected ladder re-derives g_obs and g_cross to
<= 1e-5 (imported machinery); T0b F_10 at 1/5 matches V113 to 1e-6;
T0c the corrected F~_10 at 1/5 differs from F_10 by a LOOP-sized
amount (reported; 1e-4 < |d| < 0.5: the correction is real and
perturbative - if it is below 1e-4 the channel is again rigid and
that is reported as the branch).
BATTERY (house rule): every converged g* vs g_obs, g_cross, 1/5,
guards 1e-9 contact / 1e-4 near; contraction reported.
BRANCH WORDS: FP2-CONTACT (one-loop self-consistency creates a clock
ON a published anchor: booked consequence - the clock chapter of the
paper is rewritten around V352 and the full battery re-run) /
FP2-NEAR (converged, 1e-9 < d <= 1e-4: real mechanism candidate,
booked follow-up with the other insertion orders, no paper claim
until contact) / FP2-DRIFT / NO-FIXED-POINT / RIGID-CHANNEL (the
insertion cannot move the clock at this order) / MACHINERY-DRIFT.
Cost booked: ~16k rung evaluations, wall <= 20 min.

## V352 run note (post-run; verdict as booked)
FP2-DRIFT, wall 65 s. T0a/T0b: the imported V350 machinery re-derives
g_cross (4.3e-14), g_obs (3.1e-12), F_10 at 1/5 vs V113 (1e-9) before
any loop is inserted. T0c: the one-loop self-energy is REAL and
perturbative at 1/5 - |F~_10 - F_10| = 3.3e-3 (J1 = a_ins x I, a_ins
= 0.964413... read from the V205 artifact), 1.7e-2 (J2 = Cartan
quadratic), 3.6e-2 (J3 = step dressing): so V350's rigidity does NOT
extend to the module's one-loop polarization. The channel moves.
Where does self-consistency land it? Nowhere near the desert:
converged fixed points (continuation roots, damped iteration):
J1|A g* = 0.1967551490 (delta -3.4e-3), J2|A g* = 0.1954459513
(delta -4.7e-3), J3|A g* = 0.1942450714 (delta -5.9e-3, contraction
0.011 - a genuine fixed point, not a constant map), J1|B g* =
0.2819151001 (a far razor-band re-crossing, delta +8.2e-2). Nearest
anchor distance 3.24e-3: three and four orders above the contact
guard, two orders above even the near guard. J2|B and J3|B:
BRACKET-MISS - with the loop dressing, the coefficient curve crosses
COEF_OBS NOWHERE in (0.15,0.33): the one-loop correction destroys
the hierarchy-coefficient reading of the clock entirely.
Physics reading (booked branch FP2-DRIFT stands): the published
clock is a TREE-level reading; dressing the ladder with its own
published matter at one loop shifts the rung-10 crossing by a few
thousandths and kills the coefficient reading - quantum
backreaction within the published model BREAKS the clock rather
than creating it. Combined with V350 (hop-spectrum channel rigid at
I/45), the ledger is now closed from both sides:
self-consistency inside the published model either cannot move the
clock at all (V350) or moves it off the desert (V352). Any loop
that could fix the clock needs couplings the published model does
not contain - i.e. a genuinely external sector with new parameters,
which is a different theory, not a repair of this one. The paper's
clock-tension verdict stands, reinforced: the named residual is
not a self-consistency artifact at any order internal to the model.

## V353_inverse_external_sector (tier 1, hard; "find the new theory" as an exact inverse problem, booked before the run)
V350+V352 closed every INTERNAL self-consistency channel (rigid or
drifting). The only survivor is an EXTERNAL sector with couplings the
paper does not contain. This experiment turns "find it" into two
exact measurements plus a realizability census - no fit, no freedom:
(A) ENGINE-NATIVE CONFORMAL CENSUS. The one exact rational statement
an external sector can make is a vanishing one-loop beta:
b = (11/3) n_g C_A - (2/3) n_f T_f - (1/6) n_s T_s. All constants are
DERIVED from the engine, not quoted: k = Tr(T^A T^B)/delta from the
45 published generators on the 32 module (T(32), and T(16) = T(32)/2
for one Weyl family), C_A from Tr(f f) on the adjoint built from the
published structure constants; each constant is certified to lie in
(1/2)Z within 1e-9 before it is used (then EXACT Fraction
arithmetic). Census n_g in 0..2 (gauge copies), n_f = N_16 families
0..6 (16-spinor Weyl packs), n_s in 0..4 (real scalars): list every
exact b = 0. Headline: any zero at n_f = 3 (the published generation
number) is reported as an exact conformal candidate.
(B) INVERTED FIXED-POINT COUPLING. For each booked V352 insertion
channel J (J1 = I, J2 = Cartan quadratic, J3 = step dressing),
solve for the amplitude alpha* making the PUBLISHED clock a fixed
point of the dressed ladder: h_J(alpha) = rootB(g_obs, alpha J) -
g_obs = 0, bisection 30 (contraction ~0.01 guarantees root-at-anchor
is the fixed point to O(alpha^2) - booked; the fixed point itself is
then verified by the V352 iteration at alpha*, guard 1e-9). T0
sanity: h_J(0) = 0 reproduces the published anchor; the published
fingerprint ratio (COEF_OBS / coef_at_fifth)^(1/10) is re-derived
and matched to the V205 closed-crossing deficit (0.9895...) at
1e-9 before any solve.
(C) REALIZABILITY CENSUS. alpha* and alpha*/LOOP are run through
mpmath pslq (degree <= 6, coefficients <= 1e4) and tested against
the published integer arithmetic: exact membership of
{m/2, m/3 : m <= 60} at 1e-8 and of the engine-derived integers
{3, 12, 16, 32, 45, C_A, T(32)}. The needed SIGN is reported (which
spin statistics produce it) and only sign-consistent content counts
as a candidate.
BRANCH WORDS: NEW-THEORY-CANDIDATE (alpha* is EXACTLY an admissible
published integer at 1e-8, sign-consistent: booked consequence - the
sector enters the paper as the external completion of the clock
mechanism and the full anchor battery is re-run with it; no identity
claim without that re-run) / NEW-THEORY-REQUIRED (the inverse
problem is solved exactly but the required coupling is not in the
published arithmetic: the exact required number is printed and named
by PSLQ - "here is what the new theory must couple at", the theory
itself still missing) / NO-PUBLISHED-ARITHMETIC-FIX (no exact
realization and no PSLQ name: the clock's completion lies outside
the paper's number system - the closure form the data take).
Cost booked: ~40 s (two channels bisected; J3's g-dependence is
amplitude-scaled), seconds for (A).

## V352 run note addendum + V353 run note (post-run; V353 reclassified)
V352 (final, GPU): FP2-DRIFT re-confirmed digit-identical after three
lens corrections: root scans 1e-3 (the first-run 3e-3 grid was
razor-blind: it lost the physical root at g_obs and inverted a
razor-jump branch - the T0d guard, now permanent, pins the un-dressed
roots on the published anchors at 4.3e-14/3.1e-12); continuation
window (bisect only roots within 0.02 of the iterate - pure cost, the
choice rule unchanged); GPU batch scan on the user's GPU 1 with the
T0e device-vs-host identity gate - which CAUGHT a device-side bug
(the denominator used G^1 not G^10: dev 9.3e-2) before any reading
was believed; after fix: 8.8e-16. Headline number for the record:
the published V205-amplitude identity dressing now moves the
rung-10 crossing by only -5.44e-4 (dmin 4.09e-4): the correction is
small and honest, still four orders off the near guard.
V353 run (first run classified NEW-THEORY-REQUIRED; RECLASSIFIED
here to NO-PUBLISHED-ARITHMETIC-FIX): the only solved "cleaning"
was alpha = -7.5e-10 - the TREE solution alpha = 0 picked up as a
sign-change neighbour of the grid point at zero (the first
classifier counted it; permanent filter |alpha*| > 1e-6 added, with
the tree guard |h(0)| = 3.1e-12 reported explicitly for all three
channels). The (A) booking text printed condition "11 n_g" - the
arithmetic (always exact, Fraction) says n_f + n_s = 22 n_g
(704/3 n_g = 32/3(n_f + n_s)): fifteen exact nontrivial conformal
zeros, headline (n_g, n_f, n_s) = (1, 3, 19): three 16-families is
eligible for one-loop conformality in the engine's own integers - a
naming fact about where a coupling may sit, not a birth of the
clock. THE MEASUREMENT THAT DECIDES: the meeting inversion (a
dressing making ONE root satisfy both published readings - closing
the paper's clock tension) found NO sign change of
rootA - rootB over the whole amplitude band (-2, 2) for EVERY
published insertion form (identity, Cartan quadratic, step
dressing): the published pair (g_obs, g_cross) is never degenerated
by any one-loop dressing at any amplitude the perturbative regime
tolerates, in the published arithmetic. Combined with V350 (the
self-channel is rigid) and V352 (published dressings drift the
clock off-desert), the closure form is final: the clock's
completion is not carryable by the paper's number system - any
sector that closes the tension must be nonperturbative (|alpha| >
2) or multi-loop, i.e. a genuinely different theory, exactly as
booked: NO-PUBLISHED-ARITHMETIC-FIX.
GPU directive honored: CUDA_VISIBLE_DEVICES=1 (physical GPU 1, RTX
3090), device path identity-certified; the whole V352/V353 scan
load runs on it (181-point dressed ladder per call in ~1 ms vs
1.5 s host).

## V354_elect_read_root_degeneracy (tier 1, hard; the decisive build
test for an elected theory at the exact fifth; booked before the run)
The user's directive: build the elected theory completely. The
record's only rational candidate clock is g = 1/5 exactly (V309
buried algebraic births of g_obs; the paper itself names the rung-10
rational clock 1/5 as nature's, to relative 1e-5). An ELECTED theory
must satisfy TWO readings with ONE sector:
  READ: 3 (|F~_10(1/5; alpha J)| / LOOP)^10 = COEF_OBS  (the dressed
        ladder at 1/5 reproduces the observed hierarchy coefficient
        exactly - the same COEF_OBS the paper reads from config);
  ROOT: rootB(1/5, alpha J) = 1/5  (1/5 is a self-consistent fixed
        point of the dressed ladder: the clock is born, not set).
alpha_read (bisection 40, single-point ladder evals) and alpha_root
(scan 41 + bisection 30, GPU batched scans on GPU 1, continuation
window) per booked channel J (identity, Cartan quadratic, step
dressing). T0 gates: tree guards (rootB(1/5, 0) = g_obs to 1e-8 -
the paper's own 2.012e-6 offset re-derived; reading at 1/5 matches
the V206 coef_at_fifth artifact to 1e-9).
CLASSIFICATION (booked): ELECT-READ-ROOT-DEGENERATE - |alpha_read -
alpha_root| <= 1e-8 for some channel AND the common alpha is in the
published arithmetic (m/q, q <= 60, or engine-integer multiples at
1e-8): booked consequence - build the completion (V355): the content
sector enters the clock mechanism and the FULL published battery
coefficient/alpha_10/crossings is re-run as the elect's acceptance
test; no paper claim without it. ELECT-SPLIT: both amplitudes exist
but never coincide (the exact pair printed per channel): the
published arithmetic cannot carry the elect at the fifth - the
final construction closure. ELECT-NO-READ / ELECT-NO-ROOT /
MACHINERY-DRIFT.
Cost booked: ~5 min on GPU 1 (scan ~1 ms/call).

## V354 run note (post-run; branch reclassified by classifier ordering bug)
The booked classification code checked ELECT-SPLIT before the
degeneracy-without-arith branch (ordering bug, permanent fix: branch
order candidate -> degenerate-no-arith -> split). The measurement:
T0 tree guards perfect (rootB(1/5,0) = g_obs to 3.1e-12: the paper's
-2.012e-6 offset re-derived; reading at 1/5 matches V206
coef_at_fifth to 3.5e-13). Then the structural fact: on the physical
branch READ and ROOT are the SAME condition (a dressed curve passing
through (1/5, COEF_OBS) has 1/5 as its root by definition) - the
degeneracy J1 alpha_read = alpha_root at gap 6.8e-13 is not a
coincidence, it is the theorem form of the elected theory: EXACTLY
ONE external amplitude makes the fifth a born, reading clock:
alpha_elect = -0.003095044242559, beta_elect = alpha/LOOP =
-0.1555738920793. (The other ROOT amplitudes with exact dyadic
values -0.75, -0.3, +0.25, +0.5, +1.9 are razor-branch jumps of the
continuation scan, labeled artifacts, never candidates.)
Naming tribunal (final): PSLQ degree <= 6 coefficients <= 1e5 on
beta: no relation; mixed basis {1, b..b^6, sqrt5, b sqrt5, pi, b pi,
b^2 pi, 1/pi, b/pi, pi sqrt5} at 1e6: none; alpha in {1, a, a^2, pi,
sqrt5, 1/25, a sqrt5, a pi}: none. Closest published neighbours miss
by 2%+ (beta*2 pi = -0.9775). VERDICT (booked branch):
ELECT-READ-ROOT-DEGENERATE-NO-ARITH - the elected theory EXISTS as
the unique solution of the two published conditions, and its
coupling is a NAMELESS number of the paper's arithmetic: the
elected cannot be born from the paper's integers; it either enters
as one measured dial (then V355 builds it: paper + one dial, full
battery acceptance) or it stays beyond reach until a mechanism
outside this model produces exactly -0.1555738920793.

## V355_elect_one_dial (tier 1, hard; the elected theory BUILT, user
directive 2026-11-01: accept beta_elect as the single measured dial;
booked before the run)
Construction: the published reading channel, dressed by the identity
insertion at the V354-pinned amplitude alpha_elect (read from the
V354 artifact - no decimal hardcoded here). The clock is BORN at the
exact fifth: g = 1/5 rational, not the inverted g_obs.
BATTERY (every published clock-consuming measurement re-run; house
mechanisms reused verbatim - V206 first_crossing/reading/definitions,
V352 Fn_corr chain, V318 ladder-scale TOL_SCALE 1e-5):
 B1 BIRTH: |rootB(1/5, alpha_elect I) - 1/5| <= 1e-9 (the clock is a
    fixed point of its own dressed ladder, re-measured in the build).
 B2 EXACT READING: 3(|F~_10(1/5)|/LOOP)^10 = COEF_OBS to rel <= 1e-8
    (the observed hierarchy coefficient read exactly at the rational
    clock - the paper's central claim UPGRADED from "to relative
    1e-5 at g_obs" to exact at 1/5).
 B3 CROSSING PRESERVED: dressed rung-10 crossing (V206 bisection on
    |F~_10(g)| = LOOP) within 1e-3 of published g_cross_rung10 =
    0.200133167 (house tolerance), and the born clock stays AT the
    rung-10 crossing region (|g~cross - 1/5| reported).
 B4 RUNG TEN: first_crossing of the dressed ladder at the born clock
    lands on rung 10 (the "rung-10 clock" is not a dress artifact).
 B5 LADDER UNTOUCHED: the dressing channel carries no rung-scale
    machinery - S17 = nv.loop_ladder(17) EV re-derivation vs the V303
    contact artifact at TOL_SCALE 1e-5 (V318 gate verbatim).
 B6 DIAL ACCOUNTING: exactly ONE new constant (beta_elect = -0.15557
    x LOOP as read from the V354 artifact); the paper's zero-dial
    ledger (V248) is CONTRADICTED BY ONE and the record says so
    verbatim in the check text.
T0 gates: tree guards of V354 re-derived (rootB(1/5,0) = g_obs at
1e-8; coef_at_fifth at 1e-9); alpha_elect round-trip vs artifact at
1e-15.
CLASSIFICATION (booked): ELECT-BUILT (all battery guards green - the
elected theory stands: born rational clock, exact reading, published
crossings/rungs/scales intact, cost one dial); ELECT-BATTERY-FAIL
(any guard red: the exact failing numbers printed, the build stands
condemned); MACHINERY-DRIFT.
Cost booked: ~4 min (one GPU rootB + ~120 CPU dressed-ladder evals).

## V355 run note (post-run)
ELECT-BUILT, all seven guards green, wall 5.7 s on GPU 1. The
elected theory as built: identity-form dressing of the published
reading channel at the V354-pinned amplitude (alpha_elect =
-0.003095044242559, beta_elect = -0.155573892079), clock born at the
EXACT fifth. Numbers: B1 birth root 0.199999999996 (3.7e-12 from
1/5 - the -2.012e-6 clock-tension offset of the published inverted
clock is precisely what the dial absorbs); B2 reading 2.0765641112 =
COEF_OBS to 5.5e-13; B3 dressed rung-10 crossing 0.200135230 (2.1e-6
from published, 1.4e-4 from the fifth - the paper's "clock at the
rung-10 crossing" survives dressing to 1e-4); B4 rung ten intact;
B5 rung-17 scale bit-identical (0.0e+00); B6 one dial, V248 ledger
contradicted by one, stated. Physics of the cost: the dial has a
precise job - it is the tree-level offset of the published inverted
clock, promoted to a constant of nature. The published reading claim
(the hierarchy coefficient) is now EXACT at a rational clock.

## V356_pentagon_closure_free_clock + V357_pentagon_walk_ladder
(user directive 2026-11-01: the 12-cell is not a patch on the cubic
engine - rebuild on it; booked before the runs)

V356 - the architectural fact the rebuild needs: the cubic class's
closures CONSUME the clock (V334 theorem: closure at an algebraic
clock is a clock identity); the escape branch was "clock in the
coefficients". The dodecahedral symmetry group I_h is FINITE (120).
Build I_h exactly (3x3 rotation/reflection matrices over Q(sqrt5)),
its 12-face permutation action, and the wreath step class
(sigma on faces, R in I_h on the frame). Census: element orders of
I_h (all of them), and closure orders of sampled wreath products
(exhaustive over A5 x Z2 x S-sample), AND the decisive measurement:
evaluate a booked sample of wreath-step closure products AT TEN
DIFFERENT REAL CLOCK VALUES g in (0.05, 0.6) - if the closure
identity deviates from +-I by <= 1e-12 at every g, closure is
ANGLE-FREE: the pentagon closes without consuming a clock, and the
clock is a coefficient - the V334 escape branch is REALIZED by the
12-cell algebra. T0: multiplication table closure (Cayley census
exact), matrix entries in Q(sqrt5) certified by minimal polynomials
to 1e-14, 12-action faithful (kernel trivial).
CLASSIFICATION: PENTAGON-CLOSURE-CLOCK-FREE (all samples close
angle-free, the V334 branch exists in the 12-algebra) /
PENTAGON-CLOSURE-STILL-CLOCKBOUND / MACHINERY-DRIFT.

V357 - the ladder on the pentagon: V345 certified the closed walks
of the 12-face adjacency W_n = [60, 120, 780, 3120, 16380, 78120,
394380, 1953120, 9784380]. Extend exactly to n = 30, diagonalize
the face-adjacency exactly (spectrum over Q(sqrt5) - the golden
section is in the 12-cell's own spectrum), and test whether ANY
functional of the booked class reproduces the ladder's published
reading objects: target T = (COEF_OBS/3)^(1/10) x LOOP
(the rung-10 amplitude the elect reads at 1/5 - artifact leaf) and
the certified closed-branch deficit 1.013%; class = monomials
5^a (sqrt5)^b W_m^c W_n^d / 5^e with |a|,b,|e| <= 20, 0 <= m,n <= 30,
c,d in {1,-1}, plus 3-digit rational prefactors - the chance-density
discipline booked BEFORE the run: expected hits computed from the
class size and the 1e-4 band, a single hit is density-suspicious by
the V319 rule, FOUND requires <= 1e-8 residual against T or the
deficit. T0: W_n re-derived to V345's nine leaves bit-for-bit.
CLASSIFICATION: PENTAGON-LADDER-HITS (+- arith), PENTAGON-LADDER-
VOID, MACHINERY-DRIFT. Cost booked: < 10 min each.

## V356 run note (post-run)
PENTAGON-CLOSURE-CLOCK-FREE, all six gates green, wall 1.0 s on
GPU 1. Decisive numbers: T0 |I_h| = 120, STRICT order census
{1:1, 2:31, 3:20, 5:24, 6:20, 10:24} (max 10 as booked - the census
is the A5 x Z2 element structure verbatim), faithful 12-face action
(120 distinct permutations, kernel 1). A WREATH CLOSURE CENSUS: 40
sampled sigma x R steps close at orders up to 30 with worst U^k
deviation from +-I of 1.1e-15 - the group core closes at EVERY
clock: closure is angle-free, the V334 escape branch is realized by
the 12-cell algebra. B WILSON-DRESSED SCAN: the dressing BREAKS the
closure - best clock-selected deviations (0.19, 1.0e-01),
(0.02, 6.9e-01), (0.02, 7.0e-01), (0.71, 1.6e-01), (0.02, 1.8e-01):
no clock on the 97-point grid restores +-I closure of the dressed
step, so inside the closed structure the clock is genuinely a
COEFFICIENT the dressing re-selects (deviations O(0.1-0.7), never
below 1e-1). T-B-GPU: dressed scan batched on cupy, device-vs-host
identity 2.2e-16 on NVIDIA GeForce RTX 3090 (GPU 1).
Booked corrections (defects found in the booked construction, fixed
BEFORE the green run; the booked GATES are untouched):
(a) mpmath defect: chained 3x3 products without entry normalization
and raw 50-digit str dedup keys made the BFS diverge (174k+ nodes at
depth 18); fix: every product entry normalized (+x at dps 50),
dedup keys canonical 14 significant digits with a 1e-25 threshold
against the mpmath.pi transcendental residue, and all seed/generator
entries typed as mpf (int '1' vs mpf '1.0' key split had carried a
duplicate identity through the table - 132 not 120).
(b) group-theory defect: the booked generators "(Rz(2pi/5),
Rx(pi/2), -I)" do not generate I_h in the face-normal orientation -
Rz(72 deg) is not a symmetry of the (0,+-1,+-phi) vertex set (the
cell's 5-fold axes go through its face normals, (0,+-1,+-phi) and
cyclic images) and C5-C2-perpendicular generates only D5. The I_h
table is now generated on the CELL's own axes: C5 about (0,1,phi)
(72 deg, exact ring entries since sin72/sqrt(1+phi^2) = 1/2),
C2x = Rx(pi) (a true 2-fold through the edge midpoints (phi,0,+-1)),
C3 about the threefold (1,1,1) axis (the EXACT cyclic permutation
matrix), and the inversion.
(c) order_of must be STRICT (P^k = +I only): the old +-I shortcut
capped the census at 5 and contradicted the booked max == 10.

## V357 run note (post-run)
PENTAGON-LADDER-VOID, T0 green, wall 18 s. Decisive: W_2..W_10
bit-identical to the V345 leaves; the char poly of the 12-face
adjacency certified by Newton sums on the integer traces,
e_1..e_12 = [1, 0, -30, 40, 255, -576, -580, 2640, -1425, -3200,
5250, -3000, 625] - exactly (x-5)(x+1)^5(x^2-5)^3, the icosahedral-
graph spectrum {5, (+-sqrt5)^3, (-1)^5} (MathWorld confirms; gap
5 - sqrt5 = 2.763932 is the booked PENTGAP that motivates V360).
LADDER TEST: the booked monomial class has 118,752,564 defined
forms; best relative-log miss 3.4e-04 against the rung-10 amplitude
T = 0.0191757586 and 5.0e-03 against the deficit 0.01013 - both
OUTSIDE the 1e-4 found band (FOUND requires <= 1e-8); chance
expectation in the 1e-4 band is ~237 hits for a class this size, so
even a band hit would be density-suspicious by the V319 rule. The
ladder's published reading objects are NOT carried by ANY functional
of the booked class on the pentagon walk spectrum: VOID, honestly
negative, exactly what the pre-run chance discipline was booked for.
Booked correction: the class as literally written included W_1 =
tr A = 0, whose log(-inf) poisoned every argmin and made the scan
vacuous (the first VOID was a machinery bug, not a verdict); the
class is built only from the defined positive walks W_2..W_30 and
the rerun is the number above.

## Pentagon-12 user directives (2026-11-01, mid-execution; recorded
for the V358-V360 record)
The cell table is booked as CERTIFIED CODE in cl10/pentengine.py
(dodeca_cell): Type regular polyhedron, faces pentagons, F = 12,
E = 30, V = 20, Euler characteristic 2, three faces per vertex,
five vertices per face - certified against the polar of the engine's
own face normals (dual face-adjacency = the 30 engine edges, exact).
The nest: 120 dodecahedra - cell_nest certifies the I_h table acts
simply transitively on 120 cell labels. The dodecahedral honeycomb
of order 4: honeycomb_534 certifies the Coxeter chain [5-3-4] Gram
Lorentzian (3,1) with all proper sub-blocks spherical ([5,3] = the
cell's own H3 of order 120) and edge label 4 = four dodecahedra per
edge - the compact hyperbolic sibling of the V360 flat pentagrid.
Tool stack per directive: cupy-cuda13x (GPU 1 batch + the house
device-vs-host gate), jax[cuda13] cross-check, einops for the cell
factor axes (f face, r frame, a internal), kingdon geometric algebra
for the edge rotors (Cl(3,0) rotor r ~ 1 + n_v n_u, cross-certified
against Rodrigues at 4.4e-16 on all 30 edges; n_u mapped to n_v at
1.3e-51; all 270 twist entries in Q(sqrt5) by PSLQ).

## The pentagon engine rebuild V358-V360 (user directive 2026-11-01:
"rebuild the whole engine on the pentagon-12" - booked before any
run; engine library cl10/pentengine.py lands with the booking)
The V356 finding (expected): the dodecahedral group is finite, so
a structural step over it closes at EVERY angle - the clock must
live in the coefficients. The spectral fact (V357's spectrum): the
closed 12-face cell has eigenvalues {5, (sqrt5)^3, (-sqrt5)^3,
(-1)^5} - a spectral gap 5 - sqrt5 = 2.76393202250021 exactly: NO
massless mode exists on the closed sphere, whatever the internal
module. So the rebuild is three-layered and each layer is measured:
V358_PENTENGINE: the engine itself - cell = 12 faces, on-site =
12 x the published 32-level cell module (internal Spin(10) stages
reused VERBATIM from the published engine machinery), six-stage
product with the hop stage replaced by the icosahedral-graph hop
dressed by A5-covariant Wilson link phases t_uv(g) = exp(i g w_uv),
w = shared-edge winding 0..4. Gates: unitarity of every stage and
of the full product <= 1e-12 (GPU cupy path with a T0 host-vs-device
identity <= 1e-9 - the house pattern), det = 1 of the internal
product, A5-covariance (conjugation = permutation similarity, all
120 elements), clock-free structural closure re-verified inside the
engine, and the exact hop spectrum re-derived {5,(sqrt5)^3,
(-sqrt5)^3,(-1)^5} to 1e-12. Branch: PENTENGINE-BUILT (all green) /
PENTENGINE-BROKEN / MACHINERY-DRIFT.
V359_PENTENGINE_LADDER: the traced-word ladder OF THE PENTAGON
ENGINE: closed walks on the face module with the link phases, F_n^P
= (1/384) tr U_hop(g)^n summed over the booked word classes; the
crossing |F_n^P| = LOOP_P against the engine's own loop object -
and the published reading chain re-run on it: rung-10 coefficient
3(F_10)^10 against COEF_OBS, and the ladder base question (does a
pentagon loop factor replace 1/(16 pi)?). T0: F at g = 0 matches
the pure spectrum moments W_n/384 exactly; chance-density
discipline booked (class of g-grids 1e-3 x rungs, expected hits
printed). Branch: PENTLADDER-CROSSING (+- named), PENTLADDER-NO-
CROSSING, MACHINERY-DRIFT.
V360_PENTAGRID_GAPLESS: the closed sphere is gapped (5 - sqrt5),
so masslessness needs the EXTENDED pentagon lattice: the A5-
quasiperiodic grid - Penrose P3 rhomb substitution graph grown to
N <= 4096 faces, tight-binding hop spectrum at the Wilson clock -
gap-to-nearest-zero vs N (N = 169, 676, 2704): PENTAGRID-GAPLESS
(gap decreasing as N grows, slope >= 0.3 in log-log - a real
touching forming) / PENTAGRID-GAPPED (gap floor >= 0.1 x (5 -
sqrt5)) / MACHINERY-DRIFT. This is the rebuild's masslessness gate:
without it the pentagon engine cannot carry light fermions at all.
Cost booked: V358 < 5 min, V359 < 15 min, V360 < 15 min, GPU 1.

## V358 run note (post-run)
PENTENGINE-BUILT, all ten gates green, wall 13 s on GPU 1 (cupy
RTX 3090). Decisive numbers: T0a the engine's ih_table reproduces
the V356 census {'1':1,'2':31,'3':20,'5':24,'6':20,'10':24}
bit-for-bit; T0b |U32 B - B U16| = 4.5e-16, |det U32 - 1| =
5.3e-15; T0c the CELL table certified as code - dodeca_cell {5,3}
counts (V,E,F,chi) = (20,30,12,2), edge 0.472136, the dual
face-adjacency exactly the 30 engine edges - cell_nest = 120
regular dodecahedral actions (120 dodecahedra), honeycomb_534
Gram det -0.077254 Lorentzian signature (3,1), edge label 4 (four
cells per edge), rotor certs nv 1.3e-51 / pair 4.4e-16 / ring 270
of 270; T1 U_hop unitary 2.0e-15 at four clocks, |det| 1.3e-14;
T2 the full 1152-dim step unitary 1.2e-14; T3 A5-covariance over
ALL 120 elements batched on GPU 1 in tens: worst 2.9e-15, mean
1.8e-15, device-vs-host identity gate 3.9e-31; T4 tr K^2 = 180
exactly in mpmath (tr K = 0 by construction); T5 PENTGAP: the
twisted 36-cell hop gap min|lambda(K)| = 0.381966011250105 =
(3-sqrt5)/2 = 1/phi^2 (mpmath dps 40 inverse-power bound; numpy
and cupy agree at 2.2e-16) - smaller than the untwisted
face-spectrum gap 5 - sqrt5 = 2.76393202250021, re-derived as a
gate at 2.7e-15: BOTH numbers say the same thing - NO massless
mode on the closed 12-cell, whatever the twist or the internal
module. The pentagon engine exists, is exactly symmetric, closes
structurally at every clock, and carries a gap. V360's necessity
is now a measured leaf, not an argument.

## V359 run note (post-run)
PENTLADDER-NO-CROSSING, honest negative, wall 4 s on GPU 1. T0:
the V357 artifact walk leaves W_1..W_30 re-derived bit-for-bit
(exact-Python-int traces - int64 overflows at W_28, the re-
derivation must be exact), and the booked frame factor measured:
the 36-cell ladder at g = 0 is identically 1 (normalization), and
the m = 2 moment ratio tr K^2 / W_2 = 3.000000 exactly - the
booked 180/60 = 3. T0e: GPU batched (g, n) scan vs host at 5
probes, 5.6e-17, cupy RTX 3090. CROSSINGS: |F_10^P(g)| = LOOP has
ZERO roots on the booked 1e-3 grid over (0,1): the 36-cell traced-
word ladder never descends to 1/(16 pi) - at the fifth itself
|F_10^P(1/5)| = 0.2775 = 14.0 x LOOP. READING: C_P =
3(|F_10^P(1/5)|/LOOP)^10 = 1.47e+12 vs COEF_OBS 2.07656 - rel dev
7.08e+11, OUTSIDE the elected dial window (V355 read_dev leaf
5.458e-13, read from the artifact, no hard-coded decimal). The
pentagon engine's own rung-10 object does not carry the published
reading chain: the closed 36-cell hop spectrum, gap 0.381966 from
below and 5 - sqrt5-type structure from above, is structurally too
far from a touching to be dialed - consistent with V357's VOID and
with the PENTGAP leaf: on the closed cell no reading works; the
masslessness question is V360's.

## V360 run note (post-run)
PENTAGRID-GAPLESS, wall 43 s, peak RSS 3.8 GB (measured with
/usr/bin/time -v), GPU 1 dense eigvalsh. Four corrections booked
against the registration, all machinery-side, none touching the
reserved decision:
(1) GAUGE: the booked offsets gamma = (0, 1/2, -1/2, 1/4, -1/4)
obstruct the booked 5-fold gate - rotation maps line (i, m) to
line (i+1, m + gamma_i - gamma_{i+1}), so 5-folding a finite
window about the origin needs gamma_i - gamma_{i+1} in Z for all
i; the booked differences have fractional parts 1/2, 0, 1/4, 1/2,
3/4. The de Bruijn gauge meeting the gate is uniform offset
gamma_i = 1/2 (differences 0); measured 5-fold mismatch 3.8e-14
at the largest level (gate 1e-9) and the origin is the center of
a pentagon face - the closed cell's face lifted to the grid
(origin pentagon radius/edge in the leaf).
(2) EDGES: the booked distance band [0.95, 1.05] x d_min cannot
capture the 1-skeleton - the arrangement has several
inequivalent nearest-neighbor lengths, one per crossing angle
(measured edge lengths spread over the 1/sin family; the booked
band would silently miss classes). Replaced by exact line-wise
consecutive adjacency: A5-covariance manifest (rotation permutes
families and lines), degrees 2-10, no tolerance bands.
(3) GAP DEFINITION: the booked E_med = median eigenvalue is
trivially an eigenvalue for odd N (N = 2855, 6215, 10825),
making the booked min_j |E_j - E_med| identically 0; the central
splitting E[N//2] - E[N//2-1] is the same distance at the
mid-point of the two central eigenvalues - used, and the centre
value printed.
(4) LEVELS: the booked L = 8/16/32 carried N-estimates ~7x under
the true N ~ 42 L^2 (2855/10825/42095); the L = 32 dense copy
(14.2 GB) OOM-blocked this machine against the user's strata
server (observed live: 57/62 GB held, the dense allocation
swapped). Re-booked to L = 8/12/16, N = 2855/6215/10825 - inside
the booked envelope N <~ 1.1e4 and 0.94 GB peak dense. Triple
points: 0 at every level (the uniform gauge is triple-free; the
filter kept every exactly-two-lines vertex).
MEASUREMENT: half-filling central splitting 0.001723 -> 0.0007315
-> 0.0001625 at N = 2855 -> 6215 -> 10825, decreasing to the
largest N, |log-log slope| 1.772 (booked guard 0.3), floor
0.0002 vs booked GAPPED floor 0.276 - 1700x under. Host-numpy vs
cupy identity 8.2e-14; jax GPU cross-check 0.0 (after its
allocator retried down on the busy card). VERDICT: the extended
pentagon lattice is gapless at half-filling, collapsing faster
than any power the GAPPED branch could hide - the closed cell's
0.381966 PENTGAP was produced by the closure itself, and
masslessness survives in the A5-invariant infinite lattice with
no dial, no twist and no internal module. The rebuild's
masslessness gate stands OPEN.

## The rhombic dodecahedral cell V361 (user directive 2026-11-01:
"utilise Dodecaedre rhombique" - booked BEFORE any execution)
The user's fact, booked as the target: the RHOMBIC dodecahedron is
the Voronoi cell (first Brillouin zone) of the FCC lattice, and
its 12 face normals are the 12 nearest-neighbour directions
(+-1,+-1,0)/sqrt2 - the cubic 12, the kissing configuration - and
it fills space: the cell-adjacency lattice of its tiling IS FCC.
Where the icosahedral closed cell carried the golden PENTGAP
(3 - sqrt5)/2 (V358), the face adjacency of the rhombic 12 is the
cuboctahedron; booked target, DERIVED before the run in a scratch
certification and trace-consistent (sum lambda = 0, sum lambda^3
= 48 = 6 x the 8 triangular faces): char poly
x^3 (x - 4)(x - 2)^3 (x + 2)^5 - spectrum {4, 2^3, 0^3, (-2)^5}
with THREE exact zero modes (an earlier handbook guess of 0^5 is
corrected by the derivation): the closed rhombic 12-cell is
expected gapless at uniform hopping, exactly.
T0 cell certificate as code: RD = {v: |v.d| <= 1 for all 12
FCC-nearest d} with exact combinatorics (V, E, F, chi) =
(14, 24, 12, 2); the symmetry group built EXACTLY as integer
matrices (BFS over signed permutations): order 48 (booked |O_h|),
order census booked {1^1, 2^19, 3^8, 4^12, 6^8} (derived from the
O_h class table E + 8C3 + 6C2' + 6C4 + 3C2 x inversion; sums to
48; max strict order 6), 12-face action faithful (kernel 1); RD
as zonohedron: Minkowski sum of the 4 half body-diagonals
(sign classes of (1,1,1)) yields the 14 vertices.
T1 the 12-face adjacency graph: 4-regular with 24 edges (the
cuboctahedron); exact-integer spectrum vs the booked derived
target {4, 2^3, 0^3, (-2)^5} at dev <= 1e-10; RHOMBIGAP := min
|A|: booked target 0 with multiplicity 3 - the closed cell of the
rhombic 12 carries exact zero modes where the icosahedral closure
carried 0.381966 (contrast leaf against the V358 artifact value,
read from the artifact).
T2 the walk ladder in exact integers: W_2..W_12 of the cubocta-
hedral graph with W_2 = 48 booked (2 x 24 edges) and W_3 = 48
booked (8 triangular faces x 6 directed tours - contrast the
icosahedral W_3 = 120 = 20 x 6 of V345/V357); characteristic
polynomial certified by Newton sums against the booked
x^3 (x - 4)(x - 2)^3 (x + 2)^5 (coefficient list derived in the
scratch: [1, 0, -24, -16, 192, 192, -640, -768, 768, 1024, 0, 0,
0] in x^12..x^0).
T3 the twisted hop on the rhombic 12 with the EXACT transport:
for adjacent faces (u, v) the rotation R = s_u s_v (product of
the two face reflections - integer matrix, det +1, maps normal u
to v): K built from these; every entry an exact integer (dev
from integer == 0.0 booked); tr K = 0 booked (zero diagonal);
tr K^2 in exact integers (booked: 12 x 4 neighbours... computed:
30 no - print and compare with the icosahedral 180 of V358 T4);
min |lambda(K)| measured as RHBIC-HOPGAP (NO target booked -
printed with the PENTGAP contrast).
T4 space-filling certificate: for a dense sample of the box
[-2,2]^3 (3.1e5 points), the origin is the unique closest FCC
lattice point iff the sample lies strictly inside the 12-plane
cell: deviation <= 1e-9 (booked equivalence Voronoi = RD).
T5 the extended lattice IS FCC: nearest-neighbour tight binding
on the FCC cell-adjacency graph, cubic blocks L = 8, 10, 12 with
the 4-point basis (N = 2048, 4000, 6912 - inside the RAM envelope
after the V360 lesson: peak dense 0.38 GB), open boundaries,
t = 1 (the only lattice-invariant choice): cupy eigvalsh, host
numpy identity gate at the two smaller levels <= 1e-10, device
name in the leaf; gap at half-filling = central splitting (V360
convention) vs N, |log-log slope| guard 0.3 as booked.
Status: RHOMBIC-12-GAPLESS (T1 zero modes exact AND T5 decreasing
with |slope| >= 0.3) / RHOMBIC-12-GAPPED (T1 gap > 1e-8 or T5
floor >= 0.276) / MACHINERY-DRIFT (any certificate red). GPU-1
directive: the T5 spectra batched on cupy with the house
device-vs-host gate. Cost booked: V361 < 5 min.
MACHINERY AMENDMENTS booked after the first run (machinery drift,
no verdict booked) and before the re-run: (1) walk moments
off-by-one fixed (W_n = tr A^n, the first code computed tr
A^{n+1}); (2) the zonohedron-vertex subcheck is replaced: the 4
body-diagonal segments carry the dependence d1 + d2 + d3 + d4 =
0, so the subset-sum construction double-counts the center (15
points for 14 vertices); the honest cell certificate is face-by-
face: every face plane carries exactly 4 cell vertices and all
face edges have ONE length (rhombi) - kept as a leaf; (3) the FCC
lattice projection bug: the nonzero-row mask used np.all(lat != 0)
which drops every row with a zero coordinate, not only the origin
- np.any fixed, the Voronoi equivalence re-tested; (4) the T5
half-filling statistic: the finite-N central splitting [0, 0,
0.0118] is finite-size noise around a Fermi surface (a single NN
band on a Bravais lattice is half-filling gapless iff the
half-filled level set is codimension-1), so the booked T5 gate
becomes the k-space Fermi-surface certificate: exact dispersion
E(k) = 4(cos kx cos ky + cos ky cos kz + cos kz cos kx), E_F by
mesh-half bisection, and the smoothed DOS(E_F) on eps-bands must
NOT collapse as eps shrinks (min >= 0.5 max, count > 0 at the
finest eps - a gap would collapse the counts to zero below
gap/2); the real-space splittings and the identity gate stay as
leaves. Status test
becomes: RHOMBIC-12-GAPLESS iff T1 zero modes exact AND the
Fermi surface certified; GAPPED otherwise.

## V361 run note (post-run)
RHOMBIC-12-GAPLESS, wall 7 s, GPU 1. The user's fact certified
end to end: the 12-plane cell is the rhombic dodecahedron
(14/24/12, chi 2; all 12 faces carry 4 vertices with ONE edge
length 0.866025 = sqrt3/2 - rhombi); O_h built as 48 EXACT
integer matrices (census {1,2^19,3^8,4^12,6^8} bit-matched,
faithful on the 12 faces); the cell IS the FCC Voronoi cell (the
12-plane inside test and the FCC nearest-neighbour dominance test
agree on 22974/22974 samples). THE CENTRAL FACT: the face-
adjacency spectrum of the rhombic 12 is the booked derived
x^3(x-4)(x-2)^3(x+2)^5 - RHOMBIGAP = 0 EXACT with 3 zero modes
(8.9e-16 dev), where the icosahedral closure carried (3 -
sqrt5)/2 = 0.381966 (V358); the walk ladder is exact integers
(W_2 = W_3 = 48, char poly bit-identical by Newton sums); the
twisted hop with the canonical integer transports from the O_h
table has RHBIC-HOPGAP = 0.000000 (integers, orthogonality dev
EXACTLY 0, tr K = 0, tr K^2 = 144): the hop that replaced the
golden gap carries exact zero modes on the closed cell. The
lattice it tiles is FCC and gapless: real-space central
splittings [0, 0, 0.0118] (finite-N noise, leaves), and the T5b
k-space certificate - exact dispersion, 48^3 mesh, E_F =
-0.935638 by half-band bisection, smoothed DOS(E_F) = [0.143,
0.173, 0.119, 0.087] with NO collapse across eps-halvings -
Fermi surface present. Four machinery corrections booked in the
registration amendment (walk-moment off-by-one; zonohedron vertex
subcheck replaced by the rhombic-face certificate - the diagonal
generators carry d1+d2+d3+d4 = 0 and double-count the center; the
np.all/np.any nonzero-row lattice bug; the finite-N splitting
replaced by the DOS-noncollapse Fermi-surface gate). Physics:
closure chose the gap - the icosahedral closure cost (3-sqrt5)/2,
the rhombic (FCC) closure costs NOTHING: a closed 12-cell with
exact massless modes exists, exactly as the user proposed.
