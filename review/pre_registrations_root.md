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
