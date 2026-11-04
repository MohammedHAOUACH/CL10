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

## The complete rhombic engine V362 (user directive 2026-11-01:
"tu n'as pas cree un engine complet et test si Le Dodecaedre
Rhombique (12 faces) est l'elue - fait le"; booked BEFORE any
execution)
The full engine transplanted onto the rhombic dodecahedral cell,
then tested against the ELECTED theory's own battery. Cell = 12
faces x 3 frames = 36 (O_h integer transports from the V361
table) x the published internal 32 VERBATIM (supercell.internal
at the house angle 0.2, as in V358): U(g, gin) = (I36 (x) Uin .
U_hop(g) (x) I32)^6 on 1152 states, U_hop = e^{i K g} by the
eigh construction. The V352 dressing machinery is published-
engine specific (hoisted hop lifts, Cartan contractions of the
45-link ladder) and is NOT transplanted - no invented rhombic
contraction rules; the honest ladder object is the engine's own
traced-word functional F_n^Rh(g) = sum_j e^{i g n lambda_j}/36 on
the certified hop spectrum (the V359 pattern), tested against the
elected artifacts' booked windows read from JSON (no hard-coded
decimals).
T0 construction: V361 table re-derived in-experiment (O_h 48
exact, census, faithful face action). TRANSPORT RULE corrected
at booking after a pre-run probe: the V361 T3 table-order rule
(first sorted group element mapping d_u to d_v) is canonical only
up to table ordering and is NOT O_h-equivariant (probe: equivariance
dev 2.0, and its twisted zero modes are rule artifacts); the
max-trace rule ties and also fails. The canonical geometric
transport is the +60 degree rotation about the shared-edge axis
(d_u x d_v)/|d_u x d_v|: it maps d_u to d_v, satisfies
R_vu = R_uv^{-1}, is equivariant under ALL of O_h including the
improper elements (axis flip trades with angle flip), and has
RATIONAL entries (denominators 6: sin(60)/sqrt3 = 1/2 exactly) -
unitarity check dev 2.2e-16 float, exactness certified with
Python Fractions. K built from it: tr K = 0, tr K^2 = 144 in
exact Fractions. TWISTED GAP min|lambda(K)| measured as a leaf,
NO target booked: the probe found the equivariant transport OPENS
a gap (0.276) on the closed cell - the exact zero modes of the
closed rhombic 12-cell live in the UNTWISTED hopping (V361 T1
spectrum x^3(x-4)(x-2)^3(x+2)^5, 3 zeros), and the V361 T3
twisted zero multiplicity is booked as a rule-artifact leaf with
its replacement measured here. Internal identities reproduced
from supercell: |U32 B - B U16| and |det U32 - 1| at the V358
guard (1e-12).
T1 engine gates: U_hop(g) unitary at g in {0.02, 1/5, pi/4,
0.98} (<= 1e-12) and |det U_hop| = 1 (tr K = 0); the full step
unitary at (1/5, 0.2) (<= 1e-10); O_h covariance for ALL 48
conjugations (kron(P_h (x) R_h, I32)) U (...)dagger: worst dev
<= 1e-10, batched on cupy GPU 1 with the house 3-probe device-
vs-host identity gate (<= 1e-12).
T2 crossings: |F_10^Rh(g)| = LOOP over the 1e-3 grid on (0,1),
roots by 40-halving bisection (1e-9); every root within +-1e-4 of
the fifth named by CF q <= 10000; distances reported against the
PUBLISHED dressed crossing V355.gcross_dressed with the V355
window gcross_dev (leaves read from the V355 artifact).
T3 reading at the exact fifth: C_R = 3(|F_10^Rh(1/5)|/LOOP)^10
against V352.COEF_OBS; the branch test is whether the relative
deviation lands INSIDE the elected dial window V355.read_dev
(leaf from the artifact) - this is the elected theory's own
acceptance criterion.
T4 (booked non-transplant): the BORN fixed-point machinery
(V352.rootB) is published-engine specific - hoisted hop lifts and
Cartan contractions of the 45-link ladder - and is NOT ported to
the rhombic engine without invented contraction rules (booked
non-transplant, like the V358 hop-moment honesty): the rhombic
engine's clock acceptance is decided by the T2/T3 objects, which
are exactly the acceptance objects the elected battery reads
(gcross_dev, read_dev leaves from the V355 artifact).
Status: RHOMBIC-ELECTED (a fifth-near crossing inside the
published gcross window AND the fifth reading inside the dial
window); RHOMBIC-NOT-ELECTED (built, green machinery, battery
measured OUTSIDE the elected windows - booked honest-negative
branch, the verdict the V359 pattern suggests); MACHINERY-DRIFT
(any construction/gate red). The verdict never moves the elected
theory's standing: the battery decides whether the rhombic 12 is
a CANDIDATE that meets the elected acceptance, nothing more.
Cost booked: V362 < 5 min, GPU 1.

## V361 note (post-run correction, transport rule)
The V361 T3 twisted-hop transport rule (first SORTED O_h element
mapping d_u to d_v) is canonical only up to table ordering and is
NOT O_h-equivariant (V362 probe: equivariance dev 2.0); its
twisted zero-mode leaf is therefore a rule artifact, corrected
in the V362 registration before its run. The V361 verdict is
UNCHANGED: every leaf it stands on is equivariant and canonical -
the T1 UNTWISTED spectrum x^3(x-4)(x-2)^3(x+2)^5 with its 3 exact
zero modes, the Fermi-surface 144/2/2/24, the FCC and Voronoi
certificates, the Dirac crossing.

## V362 run note
RHOMBIC-NOT-ELECTED, wall 7 s, GPU 1 (RTX 3090). The user's
question - is the rhombic dodecahedron the elected one - is
answered DECISIVELY NO by a fully built engine, and every
machinery gate is green:
T0: the canonical edge-axis transport (+-60 deg about
(d_u x d_v)/sqrt3; entries Fractions of denominator 6, sin(60)/
sqrt3 = 1/2 exactly; maps every adjacent pair; orthogonality and
R_vu = R_uv^-1 exact; O_h-equivariance EXACT in exact Fractions
for ALL 48 elements including the improper ones); tr K^2 = 144 in
exact Fractions; the EQUIVARIANT twisted hop OPENS a gap: min
|lambda(K)| = 0.276142 (no target; the V361-T3 rule-artifact zero
modes are gone); the untwisted closed-cell hopping re-derived
with its 3 exact zero modes.
T1: the full 1152-dim engine step at (1/5, 0.2) unitary dev
8.2e-15; U_hop unitary 2.0e-15 at four clocks, |det - 1| = 5.4e-
15; O_h covariance over all 48 conjugations batched on the GPU:
worst 2.5e-15, device-vs-host identity probes 3.9e-31.
T2: rung-10 crossings |F_10^Rh| = LOOP: 2 roots, NEITHER within
+-1e-4 of the fifth; nearest root sits 7.96e-3 from the published
dressed crossing 0.200135230 against the elected window 2.06e-06
(about 3900x outside the window).
T3: reading at the exact fifth: C_R = 1.50149e+07 vs COEF_OBS
2.07656: rel dev 7.2e+06 against the elected dial window
5.458e-13 - nineteen orders of magnitude outside.
T4 as booked: the BORN fixed-point machinery (V352.rootB hoisted
hop lifts, 45-link Cartan contractions) is published-engine
specific and NOT transferred; the clock acceptance for this
engine was decided by the T2/T3 objects, the same acceptance
objects the elected battery reads.
Verdict: the rhombic 12 closure is a genuine lattice-fact cell -
exact O_h, rational equivariant transport, gapless UNTWISTED
hopping, certified Fermi surface - but its own engine misses the
elected theory's acceptance everywhere the fifth is felt: it is
NOT the elected clock. The standing of the elected theory
(V354/V355) is unmoved, as booked: this was a candidate battery,
and the candidate fails it. The rhombic 12 joins the certified-
structure catalogue, not the clock.

## V363 the Escher stellation cell (user directive 2026-11-01:
"le run va tester le dodécaèdre rhombique étoilé (solide
d'Escher): stellation du dodécaèdre rhombique en ajoutant des
pyramides sur ses 12 faces; ses pointes s'insèrent exactement dans
les creux des cellules voisines, créant un pavage 3D parfait et
truffé de reliefs imbriqués"; booked before execution; tier 1,
user directive).

The user's claim is a TILING theorem, and it is certified in exact
integer arithmetic, no floats. Scaling fixed: the rhombic
dodecahedron H = {v: v.d <= 2 for the 12 normals d of type
(+-1,+-1,0)}: vertices 8 x (1,1,1)-type + 6 x (2,0,0)-type, edge
sqrt3 - booked-corrected: the V361 unit-normalization edge sqrt3/2
doubled, not sqrt2 as first written; inradius rho = sqrt2, face
centers at the (1,1,0)-type
points; cell centers = the index-2 dilation of the standard FCC
lattice, generated by (0,2,2), (2,0,2), (2,2,0) - PRE-RUN
correction to this booking: the undilated generators (0,1,1),
(1,0,1), (1,1,0) belong to the V361 unit normalization and would
overlap at this scaling (their determinant 2 != vol(H) = 16);
with the dilation the fundamental volume is exactly 16 = vol(H)
and the 12 nearest centers are exactly the (2,2,0)-type points. BOOKED canonical pyramid heights - TWO certified heights
(second pre-run correction, forced by an exact single-point probe
after the first red grid): h = rho (maximal insertion: tip meets
the notch bottom exactly at the neighbor's center; booked fact
then forced by exact arithmetic: the 12 notch pyramids at h = rho
ARE the 12 pyramids of the center-to-face decomposition of H,
volumes 12 x 4/3 = 16 = vol(H), so the cell's H-core is exactly
consumed and the cell solid degenerates to the 12-spike star
point-connected at the center - the tiling remains EXACT (each
interior point in exactly one star) and volume is still 16 per
cell, but the star is a spine, not a thick body); and h = rho/2
(the canonical thick Escher relief: tips halfway into the
neighbors, clearance rho/2 to the notch bottoms, non-degenerate
cell core, this is the solid the keying theorem T3 certifies).
For every h <= rho the complementarity of T1 is exact (the notch
IS the intruding pyramid, set identity) and vol(S) = vol(H) -
12 v(h) + 12 v(h) = 16 with v(h) = (4/3)(h/rho)^3 exact. pyramid height h = eta x rho, base the fixed rhombus face):
BOOKED PRE-RUN CORRECTION: the booked scaling v(h) = (4/3)(h/rho)^3
was a wrong cube - the BASE (the rhombus face) never scales, only
the height does, so v = (4/3) x (h/rho) LINEARLY and vol(12
pyramids) = 16 x eta (eta = 1: 16, eta = 1/2: 8); conservation
vol(S) = vol(H) = 16 is unaffected (the added pyramid is congruent
to the carved one at every eta - that was and remains the exact
content of the complementarity T1). BOOKED facts, all before the
run:
T0 the stellated solid E as code: combinatorics (V,E,F) = (26,72,
48) (14 base vertices + 12 apices; 24 base edges + 48 lateral; 12
rhombi x 4 triangles), Euler chi = 2; O_h acts transitively and
faithfully on the 12 apices, apex set exactly the (2,2,0)-type
12; all face data in Fractions.
T1 the exact complementarity (the precise content of "s'insert
exactement dans les creux"): for every adjacent pair (0,c), the
notch cut out of cell c by 0's apex pyramid is SET-IDENTICAL to
that pyramid (apex and 4 base vertices exact Fractions identity);
and with h = rho the apex of 0's pyramid into c and the apex of
the opposite pyramid of c MEET EXACTLY at c's center: booked
witness coordinates apex(0,(1,1,0)) = (2,2,0) = apex(c,(-1,-1,0))
for c = (2,2,0) - tip-to-tip contact at every cell center, 12 per
cell, booked exactly; for h < rho the same tiling holds with a
clearance gap rho - h (booked family note; the canonical solid
uses h = rho).
T2 the tiling: the cell solid S = H minus the 12 intrusion
pyramids of the neighbors plus the 12 own apex pyramids tiles the
FCC lattice: certified by EXACT integer half-space membership on
the fundamental domain in lattice coordinates with denominator
M = 64 grid (64^3 = 262144 samples x ~40 neighbor cells, numpy
integers): coverage count EXACTLY 1 for every interior sample
(booked: zero uncovered, zero double-covered in the interior);
volume of S by Fraction tetrahedralization = 16 = vol(H) EXACTLY
(the pyramid added equals the pyramid carved - booked 16); and the
tiling is not a translation of H: S != H (booked witness: the
interior point (1/2,1/2,0) lies in H but is carved out of S by the
neighbor's intrusion pyramid - carved). BOOKED PRE-RUN
CORRECTIONS (both found in the first implementation, mechanics
only, claims unchanged): (i) the cone facet half-spaces must be
written through the APEX, n.(x - apex) x sgn > 0 - a raw n.x test
is meaningless for a cone whose apex is not the origin (the notch
cone at h = rho happens to have its apex at the cell center, which
is why the raw form passed there; the pyramid cones failed with
double-claims); (ii) the lattice cell enumeration is the even-sum
FCC sublattice (a + b + c even, a,b,c in -3..5, 256 cells reaching
every grid point within the pyramid range) - the first draft's
a=b=c mod 2 filter is the body-centred sub-lattice and misses the
whole (2,2,0)-type cell family; and the batch runs on cupy GPU 1
(house rule) with a host-numpy recomputation identity gate that
must fail ZERO times (all membership is integer arithmetic, so
any disagreement is a real defect).
T3 reliefs imbriqués - mechanical interlocking: along NONE of the
12 face normals can a cell separate from its neighbor without
overlap: for shift delta = d/32 along each of the 12 normals d,
S(0) + delta overlaps cell c = 2d in a witness point (booked
pre-run correction: the witness is not a closed formula - at
h = rho/2 (the thick cell, core non-empty), the overlap witness
is found by a DETERMINISTIC exact rational search grid: q = alpha
x d + beta x e with e in the perpendicular candidates {(1,0,0),
(0,1,0),(0,0,1)} intersect d-perp, alpha in {1.40, 1.41, ...,
3.10}/100 scaling (i.e. 1.40..3.10 step 0.01 along d), beta in
{0.01..0.40 step 0.01} along e, shift delta = d/32, every
candidate tested for STRICT exact Fraction membership in S(c =
2d) and in S(0) + delta; search failure for a normal = T3 red
honestly reported; 12/12 booked), while
at delta = 0 the tiling gives measure-zero overlap (T2): the
interface is a bipyramidal relief wall, mechanical keying along
all 12 face normals.
T4 the tip geometry for physics: the 12 pyramid-apex directions
are the 12 face normals again, and the apex-adjacency graph
(the 12 apices at the (2,2,0)-type points; rule dot(apex,apex) =
4 selects adjacency, edge length 2*sqrt2) is EXACTLY the
cuboctahedral 4-regular 24-edge graph: the certified RHOMBIGAP
spectrum with its 3 exact zero modes is re-derived on the Escher
apex graph (booked leaf, not imported); the interlock adds only
cell-to-cell contacts at points (zero measure): no new bond graph
on the 12-cell, so the ENGINE is V362's engine verbatim and its
battery verdict RHOMBIC-NOT-ELECTED stands - booked NOT re-run,
and the paper's Section 10.6 net statement is unchanged by the
stellation except for the new lattice fact of space-filling with
keyed reliefs.
GPU booked (user directive mid-run 2026-11-01: "utilise gpu 1"):
the T2 membership grid - 262144 integer samples x the lattice-cell
band x cone half-space tests - batched on cupy int64 on GPU 1 with
the house three-probe host-vs-device identity gate (integers: gate
must be exactly 0); everything else stays exact Fraction
arithmetic CPU (too small to batch). Cost booked: < 2 min.
Status: ESCHER-TILES (T0-T4 all green: the Escher stellation is a
perfect space-filling cell with exact insertion, exact volume
conservation, and 12-direction mechanical keying - the user's
pavage 3D parfait certified as code) / ESCHER-NO-TILE (any tiling
check red: the picture is wrong at this h) / MACHINERY-DRIFT
(construction red). The verdict never moves the elected theory:
this is a lattice-fact certification of the user's solid, the
engine question was decided by V362.

## V363 run note (post-run)

Verdict ESCHER-TILES, 5/5 green, GPU 1 (cupy, RTX 3090), 14 s.
The four booked pre-run corrections (all machinery, no claim
moved) are exactly the four defects the first implementation
showed, in order of discovery: (1) double-claims 174k -> the cone
facets must carry the plane offset rhs = n . apex (raw n . x is
meaningless for a cone not based at the origin; the eta = 1 notch
cone hides this because its apex IS the cell center); (2) after
the offset, ~70% uncovered -> the offset was stored unscaled while
the grid is K-scaled integers: s = n . Q - K . rhs; (3) the rest
uncovered -> the cell enumeration was the body-centred sub-lattice
(a = b = c mod 2, missing the whole (2,2,0)-type family), fixed to
the even-sum FCC sublattice, 256 cells; (4) T2 red on the volume
equality -> the booked eta^3 was a scaling error (fixed rhombus
base, linear height): vol(12 pyramids) = 16 eta.
Final numbers, both heights: eta = 1: single/multi/none =
238440/0/0 (boundary leaves 23704); eta = 1/2 (canonical thick
Escher relief): 182208/0/0 (boundary 79936); identity gate
failures 0 at both heights; T3 keying witnesses 12/12, first
witness d = (1,1,0): q = (7/5, 7/5, 1/100) = alpha 7/5, beta
1/100, e (0,0,1) - the search found a point deep in the flank,
consistent with the interlock being volumetric, not tip-grazing.

## Roll-up: the rebuild phase V356-V363 (2026-11-01, phase end)

Three closed twelve-cells were built, gated, and put to the elected
acceptance windows; none displaces the elected one-dial completion of
Section 10.5, and the net reading of each was booked before its run:

- pentagon (V356-V360): I_h closes clock-free (V356), the walk ladder
  is void of the reading objects (V357), the built 1152-dim engine
  carries the EXACT golden closure gap (3-sqrt5)/2 (V358 PENTGAP), no
  ladder crossing (V359), and masslessness survives only on the
  extended Penrose grid (V360 PENTAGRID-GAPLESS). Verdict: the pentagon
  NAMES the fifth, it does not BIRTH it.
- rhombic dodecahedron (V361-V362): the closed cell is gapless EXACTLY
  (RHOMBIGAP 2.3e-16, 3 zero modes - the golden gap was the
  icosahedral closure's choice, not closure's necessity), it IS the
  FCC Wigner-Seitz cell, and its complete engine is exactly
  O_h-equivariant in rationals - yet the battery rejects it: crossings
  3900x outside window, fifth reading 19 orders outside (V362), and
  the one-loop self-consistency machinery does not transfer.
- Escher stellation (V363): the user's keyed picture is a THEOREM of
  this geometry - notch = pyramid (set identity) at every height, a
  perfect tiling certified to zero exceptions on the GPU batch, volume
  conserved exactly, and the thick cell mechanically keyed along all
  12 normals; the apices preserve the exact gapless twelve. It changes
  the lattice's relief, not its physics: engine verdict carried
  verbatim from V362.

Machinery corrections booked along the way (no claim moved by any of
them): mpmath normalization + canonical dedup and strict orders for the
I_h census (V356/V357 run note), the transport-rule correction with
the twisted-gap artifact 0.276142 (V361 note), and the four V363
defects (apex-offset facet planes, K-scaled offsets, the even-sum FCC
enumeration, linear pyramid volume). The phase closes with the FULL
gate fold (276 experiments): the on-disk fingerprints predate the
published single-library changes (supercell_2i/pentengine), so the
runner's library gate forces the whole suite.

## V342 note (post-run correction, pass semantics; 2026-11-01 fold)

The phase-end summary fold (report-only, per user directive to keep the
certified on-disk verdicts rather than re-run the old suite) exposed
V342 as the one FAILED artifact of 276. Diagnosis: the first version
computed pass = (n_nodes > 0 AND n_nodes in {4,8,12,24}) - an UNBOOKED
orbit-size expectation folded into the pass flag, while the booking
says the word N2B-NODE-ORBIT-<n> REPORTS the orbit structure and the
isotropy word reports the cone shape; the physics verdict (chiral test
FAILS: 557 nodes, multi-orbit, cones anisotropic 0.023-0.354) is the
measured result and lives in the word, exactly as the house convention
established for V345 (SIGNAL-LOOP-VOID), V348 (PARTIAL), V352 (DRIFT),
V359 (NO-CROSSING), V362 (NOT-ELECTED). Correction: pass = the booked
machine gates (T0 published anchor, T1 census completeness, n_nodes >
0); the dead expectation line removed. V342 re-run SOLO on GPU 1 (only
V342; the library-change fingerprints were NOT re-gated across the old
suite per the same directive): same verdict word
N2B-NODE-ORBIT-557-CONE-ANISOTROPIC bit-for-bit (557 nodes, the
published anchor 3.4e-16, M-gap 3.1e-17), pass True. No physics moved;
the icosahedral cell still fails the chiral test - by measurement,
not by construction.

## V364 tension ledger - is the dial named by a pair of named things? (2026-11-01, booked BEFORE any run)

Question (user choice after the V363 phase close, programme B): the
elected theory costs exactly ONE dial, beta_elect = -0.1555738920793496
(V355.data.beta_elect), and the polynomial tribunal already found it
unnamed AS A NUMBER. This census asks the structural version: is beta
the DIRECT value, or the sum/difference/product/quotient, of TWO
quantities that are already machine-named in published artifacts? If
yes, the dial joins a named chain (still one dial, but no longer
nameless in isolation); if no, this is one more booked brick for O8.

POOL (reserved, harvested mechanically from these artifacts - no
hand-typed decimals anywhere): V206 (clock tension residuals), V247
(two-dial desert contacts), V335 (the 365 crossing roots + anchors),
V350 (backreaction fixed point), V352 (external backreaction: lifts,
Cartan data, FP drifts), V353 (inverse external sector), V354
(classifier dial reads), V355 (the elected build: star, reading,
crossing guards), plus the 24 provenance_table rows (numeric values
only) and the four cl10.parameter_audit module constants (E, G10,
M_P, V_EW). Harvest rule: every finite numeric leaf of data/ with
abs(value) < 1e6; key path = the NAME. TARGET: V355.data.beta_elect
(read from the artifact at run time). Clone removal: any leaf equal to
the target at relative 1e-14 is dropped from the pool (the dial
appears in V354/V355 twice by construction).

OPERATIONS (all five, exact-Decimal arithmetic at 40 digits, relative
tolerance 1e-12 with absolute 1e-15 for |value| < 1e-6): direct x ~
beta; and for ordered pairs x != y: x + y, x - y, y - x, x * y, x / y
(y != 0). CHANCE DISCIPLINE booked: pool size N <= 2000 leaves => at
most 4x10^6 ordered pairs x 5 ops + N direct = 2x10^7 tests; op-values
spread over abs < 1e12 (products of bounded pairs), per-test chance
of landing within 2e-12 of the target <= 2e-12/1e12 ~ 2e-24, total
expected false hits <= 4x10^-17 - a hit at this tolerance is a fact,
not a coincidence.

VERDICT WORDS: BETA-NAMED-DIRECT-<path> (a pool leaf IS the dial;
pass, tier 2 GREEN as a naming) / BETA-NAMED-PAIR-<op>:<xpath>:<ypath>
(pass, tier 2 GREEN) / BETA-CONVENTION-CHAIN-<op>:<xpath>:<ypath>
(booked pre-run addendum: if the only hits are pairs where BOTH
leaves live inside the V355/V354 build itself - e.g. alpha_elect
times an inverse-loop constant, a relation the build DEFINES by
convention - that is the dial chasing its own tail, reported as the
negative with the chain named, pass TRUE) / BETA-NO-NAMED-PAIR
(honest negative exactly in
the house convention: the census machinery certified, the naming
FAILS by measurement, pass TRUE, the word carries the negative) /
MACHINERY-DRIFT (pool harvest or anchor red). Tier 2, residual TRUE
(it speaks to O8). Gates: T0 anchor - beta_elect read from V355, and
the V354 dial read agrees with it (relative 1e-12, the V355 booking's
own consistency); T1 pool census booked bounds N <= 2000 and file
count = 8 artifacts + table + 4 constants; T2 the scan. Cost booked:
seconds (exact rational arithmetic over ~2x10^7 Decimal ops - under a
minute; no GPU needed, host numpy/Decimal, GPU 1 pinned per house
rule anyway). Fold: incremental after the run only.

## V364 run note (post-run)

BETA-NO-NAMED-PAIR (honest negative, pass True, tier 2, O8-facing):
anchor V354/V355 dial agreement exact bit-for-bit; pool N = 525
machine-named leaves (well under the booked 2000) from the 8
artifacts + 20 numeric provenance rows + 4 audit constants; 1 376 025
exact-Decimal tests (4 ops over ordered pairs + direct), ZERO hits -
not a direct name, not a pair-sum/difference/product/quotient, and
not even a convention-chain hit (the loop-factor constant 1/(16 pi)
exists in the artifacts as a NAME in keys, never as a numeric leaf,
so the alpha_elect tail never appears as a value pair). The verdict
for O8 sharpens: the dial is not only unnamed as a number (V355
tribunal), it is unnamed as a RELATION between published named
quantities. The door left open by this census, honestly stated:
pairs involving quantities NOT among published machine leaves (the
nameless-loop-constant value, the 132 residuals' own values, or
objects of the dressing module itself) are outside this pool -
programme A (beta as an invariant of the dressing machinery) is the
precise next cut.

## V365 the dial as an invariant of the dressing machinery (2026-11-01, booked BEFORE any run; programme A, user directive "A")

Question: V364 certified beta_elect unnamed as a number and unnamed
as a relation between named PUBLISHED VALUES; the one territory that
census never touched is the dressing machinery's own invariants -
quantities that exist only when one evaluates the published one-loop
chain at the certified point. V365 computes those invariants from
the V350/V352 machinery VERBATIM (no new physics, no new dial) and
runs the naming classes:

OBJECTS (all published): G = V350.step(1/5); J = J1-identity channel
(the dial's birth channel); alpha = beta_elect x LOOP with
beta_elect READ from the V355 artifact and LOOP = V350.LOOP; the rung
chain B_0 = I, B_{n+1} = B_n (M G) with M = sum_i exp(i(g/2)EV_i) x
hoisted_i at g = 1/5 (V352 lines, verbatim); the DRESSED rung
operators Bt_n = B_n + LOOP x sum_k B_k (alpha J) B_{n-k}, n = 1..10;
the rung transfer T = Bt_10 Bt_9^{-1} (gate: det Bt_9 != 0; if
singular, book the substitution to singular values of the pair and
report honestly).

INVARIANTS (booked bound <= 200 values, names exact; pre-run
amendment with the booked singular-pair substitution: if Bt_9 is
singular the pair contributes its 64 singular values and the
FINITE QZ generalized eigenvalues instead of an inverse-transfer,
and the bound becomes <= 260 - this branch FIRED on the first
attempt, det Bt_9 = 1.65e-32, so it is the branch of record):
cartan ratios
c_n = cartan_norm(Bt_n)/cartan_norm(B_n), n = 1..10; per eigenvalue
of T (32): Re, Im, abs, arg/2pi; tr T, tr T^2, tr T^3, tr T^4,
det T, |det T|, arg det T; word traces tr(G^a (alpha J) G^b (alpha
J)) for a,b >= 1, a+b <= 6 (10); rung sensitivities d/dg Fn at 1/5
(n = 5, 10), d2/dg2 F10, d/dalpha F10 and d2/dalpha2 F10 at the dial
(central differences, booked steps 1e-5/1e-4); V354 crossing slopes
(hR - hT)' at alpha* (step 1e-6) and (hR - hT)'' (step 1e-5).

CLASSES (each booked, closed, machine-checked):
C1 anchor-tribunal on beta alone: mpmath findpoly degree <= 6,
maxcoeff 1e5 (the V335/V355 box) plus the mixed pi/sqrt5 lattice
{1, pi, sqrt5, pi sqrt5} at height <= 1e6 - MUST reproduce the
booked zero-hits (machinery proof).
C2 pair-scan among INVARIANTS only (direct + sum/difference/
product/quotient), exact-Decimal at relative 1e-12 (V364 rules).
C3 algebra-over-an-invariant: for every v, 8-dimensional integer
relation search on (1, beta, v, v beta, v^2, v^2 beta, v^3, v^3
beta) by PSLQ (booked pre-run correction: mpmath.pslq at dps 40,
tol 1e-30 - the real relation algorithm - in place of the LLL
formulation drafted first), max coefficient 1e5, relation MUST
involve a beta-monomial; a hit is then VERIFIED
at mpmath dps 40 by recomputing that invariant through the chain in
exact mpmath arithmetic: |sum c_ij beta^i v^j| < 1e-20.
CHANCE: C2 <= 4x10^4 tests x 2e-24 ~ 1e-19; C3 <= 200 LLL searches
of dim 8 at height 1e5: booked false-hit bound 2e-11 per the
dimension-8 short-vector density at 1e18 scaling. A hit is a fact.

VERDICT WORDS: BETA-INVARIANT-RATIONAL-<op>:<x>:<y> (naming, tier 2
GREEN, the dial becomes a relation of dressing invariants - zero
cadran within the certified algebra) / BETA-INVARIANT-ALGEBRAIC-
<v-name>:<relation> (naming via C3, GREEN, verified at dps 40) /
BETA-NO-INVARIANT-HIT (honest negative, pass TRUE: the dial is
unnamed even as an invariant of its own machinery - O8 hardened) /
ANCHOR-TRIBUTAL-MOVED (C1 fires where V355 booked zero: red, the
machinery is broken) / MACHINERY-DRIFT. Tier 2, residual TRUE (O8).
Cost booked: matrix chain 32x32 host float64 + 200 dim-8 LLL: under
a minute; mpmath re-verification only on candidate hits. GPU 1
pinned per house rule (no GPU kernels needed). Fold incremental.

## V365 run note (post-run)

BETA-NO-INVARIANT-HIT (honest negative, pass True, tier 2, O8-facing)
in 10 s on GPU 1-pinned host arithmetic. The booked C1 anchor
reproduced the V355 tribunal zeros exactly (machinery green). The
booked singular-pair branch FIRED: det Bt_9 = 1.65e-32 - and the run
data show what that really is: ALL 32 singular values of Bt_9 equal
1.016e-01 to float64 (spread ratio exactly 1.00), likewise Bt_10 at
7.879e-02: at the certified point the dressed rung operators are
SCALAR x UNITARY (the dressing chain is norm-preserving there -
measured, not booked as a theorem), and a 32x32 scalar-unitary with
sigma ~ 0.1 has det ~ sigma^32 ~ 1e-31: the "singularity" is a scale
of the det test, not a rank deficiency - the QZ pair machinery
returned all 32 finite generalized eigenvalues accordingly, and the
booked substitution was the honest route either way. Cartan ratios
c_n run 1.0 -> 0.9994 over the ten rungs: the one-loop dressing
chain barely bends the Cartan norm at the dial.
Census: 226 invariants (cap 260 with the fired substitution); C2
204 756 exact pair tests: 0; C3 226 dps-40 PSLQ searches over
(1, beta, v, v beta, v^2, v^2 beta, v^3, v^3 beta) at height 1e5:
0 verified relations. Net across the two programmes: the dial is
unnamed as a number (V355), unnamed as a pair of published values
(V364), and unnamed as an algebraic object of the dressing
machinery that produces it (V365). What the census DID surface is a
certifiable geometric fact at the certified point (scalar x unitary
rung operators, unit Cartan ratios) - a lead for the next cut,
booked as a leaf here, not as a claim.

## V365 note (post-run correction, the "norm-preserving lead"; 2026-11-01)

The scalar x unitary rung structure surfaced by V365 is NOT a dial
property and NOT a fifth-clock property - a viability probe run
BEFORE booking any follow-up killed the "unitarity" programme
(V366 was therefore never booked, no experiment was run, the suite
stands at 278):
(i) the CLOSED FORM Bt_n = (1 + (n-1) LOOP alpha) B_n holds to
1.4e-17 at rung 5 for two far-apart amplitudes (alpha = the dial and
alpha = 0.5) - the J1-identity dressing commutes with the B semigroup
(B_k B_{n-k} = B_n), so dressing the identity channel is an EXACT
multiplicative renormalization of every rung;
(ii) MG itself has singular-value spread exactly 1.00 at every clock
tested (0.15, 0.2, 0.25, 1/3, 0.3, pi/4) because each hoisted hop
Ma_i satisfies Ma_i^dag Ma_i = c_i I to 1e-15 (the certified hoisting
lifts scaled partial isometries) and the step G is exactly unitary
(the published gate); the rung scale sigma(g) moves with g but the
proportionality-to-unitary never does;
(iii) correction to the V365 run note: its "lead" framing is
withdrawn - the structure is universal machinery, not a dial trace.
What this CORRECTION does is strengthen the three naming verdicts:
in the birth channel J1 the dial enters the whole ladder ONLY as the
rung-n prefactor (1 + (n-1) LOOP alpha), a global rescaling - an
invariant of the machinery that names a pure rescaling cannot exist,
so BETA-NO-INVARIANT-HIT was structurally guaranteed, and the dial's
only content is the booked (READ and ROOT) condition pair of V354.
The namelessness of beta is now EXPLAINED, not only measured, within
the elected construction.

## V366 singlet-extension census - can the certified gauge action EVER admit matter? (2026-11-01, booked BEFORE any run; reality-interface programme, first cut)

Context: V98 (one-body non-existence) + V259 (KERNEL-EMPTY: ker C2 = 0
throughout the 1568-dim published fiber = 8 x 16 (flavored Weyl) + 32 x 45
(spinor x adjoint link), min nonzero C2 eigenvalue 45.000000) + V268
(CONJUNCTION-DEAD, "the physical sector the Gauss law allows is empty").
The theorems kill the PUBLISHED fiber; no experiment ever asked whether
the certified Spin(10) action (T16 from AC.weyl_generators, root_generators
on 32, adjoint_link_operators, generators VERBATIM from V79/V251/V259)
admits gauge singlets in ANY extended content. V366 asks exactly that,
in pure algebra (the kernel of C2 = sum_a Q^a Q^a with Q the diagonal
total charge is the singlet space, grid- and omega-independent by V259's
certified q-independence - no grid, no couplings).

FIBERS (booked, all built from the certified generators only):
E1 = 16 (x) 16-bar (conjugate via -T^T), dim 256;
E2 = 16 (x) 16-bar (x) 45_diag-link, dim 11520 (the meson+flux fiber -
    the E1 content V251 never carried);
E3 = 16 (x) 16 (x) 45 (diquark+flux, no conjugate), dim 11520;
E4 = 32 (x) 32 (x) 45 (published block doubled), dim 46080;
E5 = 1568 (x) 16-bar (published fiber + one conjugate flavor), dim 25088.
GROUP-THEORY PREDS (diagonal compact-group action, complete
reducibility; booked as MACHINERY anchors, not as physics): E1 ker = 1
(16 (x) 16-bar = 1 + 45 + 210); E2 ker = 1 (inv of 45 (x) 45); E3 ker = 0
(45 not in 16 (x) 16); E4 ker = 2; E5 ker = 9 (8 flavors + 1 from
(32 (x) 45) (x) 16-bar via 45 (x) 16 > 16). A machine deviation from
these is EXT-GROUP-MISMATCH (red: the generator machinery broke, since
the preds are standard character theory).

MACHINE: sparse C2 per fiber (45 sparse charges, kron-built); kernel
dimension by Lanczos (eigsh sigma = 0, count |w| < 1e-9, repo-wide
exactness); kernel basis W with W^dag W = I to 1e-12 and ||C2 W|| < 1e-9
(booked gates); anchor E0: reproduce V259's published-fiber result
verbatim (ker = 0, min nonzero C2 = 45.000000 to 1e-9) BEFORE trusting
the extensions.

VERDICT WORDS: SINGLET-ADMIT-<ker dims per fiber> (extensions carry
gauge-invariant matter - the emptiness is a CONTENT theorem, not an
action theorem; the booked cascade V367 = dispersion of the singlet
sectors on the torus + the F9/F2/F8 legs at those states is then live)
/ EXT-EMPTY (every extension kernel zero: the obstruction is deeper than
content - red for the program's second door, honest leaf for O-cards) /
EXT-GROUP-MISMATCH / MACHINERY-DRIFT (anchor E0 fails). Tier 2. Cost
booked: sparse Lanczos on fibers up to 46080: minutes. GPU 1 pinned per
house rule; algebra is host scipy.

### V366 post-run note (2026-11-01, same day) - booked discipline held
VERDICT SINGLET-ADMIT-E1:1-E2:1-E3:0-E4:2-E5:9, exit 0, 18.7 s. Every
booked element fired exactly: E0 anchor (V259 verbatim: ker 0, min
nonzero C2 45.000000), and the five character-theory kernels (1, 1, 0,
2, 9) with second C2 level 64.0000 = adjoint Casimir. One booked-gate
correction is factual, not interpretive: Lanczos (which='SM') returns
NON-ORTHOGONAL bases on degenerate kernels (caught live by the booked
W^dag W gate at 3.8e-02 / 2.0e-01 while ||C2 W|| was already 5e-15);
the basis is now Gram inverse-square-root re-orthonormalized (final
defects <= 1e-15) and the Gram full-rank assert is booked into the
script. Kernel dimensions are Lanczos-robust and corroborated by the
exact character-theory match. Scope sentence to repeat wherever the
2026 emptiness result is quoted: V259/V268's emptiness is a theorem of
the PUBLISHED CONTENT (one-body + one-flux fiber), not of the certified
gauge action - V366 measures the action admitting singlets (multiplicities
2 and 9) the moment the content carries conjugate pairs, the conjugate
being necessary (E3 ker 0). V367 (torus placement + dispersion of the
E4/E5 singlet sectors, anomaly-ledger re-certification for any surviving
content, V267-style legs) is live per booking.

## V367 singlet-pair contact composite at the zone centre (booked BEFORE run; reality-interface second cut, follows V366 SINGLET-ADMIT E4:2)

Machine: the E4 algebra (32 (x) 32 (x) 45, kernel 2 from V366) lifted to
the engine's ONE-PARTICLE space: each 32 carries the 8 sublattice
index, so single-particle space is 256 = 8 (x) 32 and the certified
clock step at q = 0 is bloch(0, U32) with U32 = SC.sectors(G10)
(internal(g10)) - VERBATIM V259/V239 machinery. Two-body fiber
256 x 256 x 45 with extended charges Q_a (x) I8; kernel W_ext =
W_E4(V366 Gram-orthonormalized) (x) sublattice-pair basis: dim 128.

BOOKED GATES: C1 anchor bloch(0, U32) = kron(P(0), U32^6) exactly (the
V239 factorization, defect <= 1e-12); C2 anchor C2_ext W_ext <= 1e-9
(the V366 kernel is still annihilated after lifting); G1 U_eff =
W_ext^dag U_tot W_ext unitary to 1e-10, U_tot = bloch(0) (x) bloch(0)
(x) Phi, Phi = e^{i omega} I45 with omega the V251 artifact rows
(0, +-8.4181); G2 gauge covariance [C2_ext, U_tot] <= 1e-9;
G3 anomaly ledger unscathed: tr_32 T^a{T^b,T^c} = 0 (booked 30 pairs,
<= 1e-12) - the doubled Dirac block is vectorlike by construction.

MEASURE: the 128 composite eigenphases of U_eff vs the booked
two-quanta continuum window {phi_i + phi_j} of single-particle phases
at (q, -q) = (0, 0) contact: composite isolated (phase outside window)
= COMPOSITE-GAMMA-BOUND-<gap> : a gauge-invariant asymptotic two-body
state exists at the certified clock - the first nonempty physical
sector candidate since V268; inside = IN-CONTINUUM (no isolated
asymptotic state at this locus - honest leaf, V368 relative-momentum
dispersion still booked to test the interior). Gate failures = *_DRIFT
(red, no interpretation). Tier 2. Cost: GPU einsum contractions on
2.9M-dim fibers, minutes.

HONEST BOUNDS booked with the run: contact pair at Gamma only (no
relative momentum, no separation - the flux "string" is one adjoint
quantum at the same block); anomaly cancellation of the PUBLISHED
content is untouched by G3 (extension is vectorlike); if BOUND, the
full claim (isolated asymptotic state) still needs V368 dispersion +
F-legs before any paper sentence.

### V367 post-run note (2026-11-01) - COMPOSITE-GAMMA-BOUND-1.007e+00
Two booked gate WORDINGS were factually wrong and are corrected on the
record (the machine results are unchanged): (1) G2a "[U32,T^a] = 0"
contradicts the theory's own statement that the clock lives in the
coefficients - the internal layer is a fixed Spin(10) DIRECTION; the
true certified statements are plane-generators-in-span(Ts) (45
residuals <= 1.6e-15), product-vs-U32 4.6e-15, [C2,U_tot] spot
1.3e-15, and singlet-subspace invariance via G1 compression
unitarity 8.0e-15. (2) The min-gap verdict now EXCLUDES the omega=0
row, booked as a booked CONTROL (flux-carrying fiber at zero flux
energy = pure bond band, merge fired at gap 0.0); BOUND read on the
two physical rows (+-8.4181), gap 1.007 rad. PHYSICS: the singlet
pair is blind to the clock's internal rotation (diagonal group
action fixes invariant tensors) - it feels only the bond layer and
the flux phase. V368 booked: relative-momentum dispersion of the
128 contact bands + multi-flux continuum overlap; NO paper sentence
before V368.

## V368 torus rigidity of the composite isolation window (booked BEFORE run; reality-interface third cut, follows V367 COMPOSITE-GAMMA-BOUND-1.007)

Structural fact proven on the machine after V367 (probe, booked here as
the hypothesis under test): each kernel column's contact block equals
the pure bond layer kron(P(0),P(0)) to 4e-15, the two columns identical
to 2.2e-16 - the certified clock's internal Spin(10) layer acts as the
IDENTITY on the E4 singlet kernel at the contact locus. Since
bloch(q,U) = kron(P(q), U^6) (V239 anchor, V367 C1) and U^6 is
q-INDEPENDENT, the prediction is: the V367 isolation gap is INDEPENDENT
of the relative momentum and the total momentum - the 1.007 rad window
is rigid across the whole torus. V368 tests it, not assumes it:

MACHINE: V259's 6^3 q-grid (216 points): single-particle spectra:
beta(q) = phases of P(q) (8) at every q, and the internal phase set
iota = phases(U32^6) (32, q-independent - measured constancy booked:
spread of phases(U^6) over q is identically 0 since q-free, checked on
3 q-points). Per total momentum K (216 of them on the torus): COMPOSITE
set C_K = {beta_a(k1) + beta_b(k2) + omega: k1+k2 = K mod G over ALL
216 relative pairs, a,b in 8} x the two physical omega rows; CONTINUUM
set T_K = {beta_a(k1) + beta_b(k2) + iota_i + iota_j: same pairs, i,j
subsampled to 24 booked random k-points per K (seed 0x51 booked)}.
gap(K) = min circular distance(C_K, T_K); BOOKED ANCHOR: gap at K=0
equals V367's 1.007e+00 to 1e-9 (identical machinery reproduction).
omega=0 stays a booked control (merge), excluded from BOUND reading.
Verdicts: RIGID-WINDOW-<minK gap> (isolation is a global property of
the clock's internal phase spectrum, not of the locus - the composite
band is separated from the single-flux two-quanta continuum EVERYWHERE
on the torus) / RIGID-WINDOW-BROKEN-<K loci> (the composite disperses
into the continuum somewhere - honest, loci listed) / ANCHOR-DRIFT.
Multi-flux continua (n >= 2 flux quanta with their own energies) are
NOT in the booked comparison at this cut - booked as the V369 open
overlap question, no overreading. Tier 2, host numpy (eigvals 256^3 x
216 ~ minutes), GPU pinned per house rule.

### V368 post-run note (2026-11-01) - booked prediction REFUTED, honestly
WINDOW-NOT-RIGID-2.0e-05: the contact window 1.006652 (anchor
reproduced to 1e-9) collapses to 2.0e-05 minimum over the 216 total
momenta - at fixed K the continuum is the union over all relative
momenta, and free propagation + Gauss offers no threshold. Booked
verdict WORDING corrected pre-artifact (RIGID-WINDOW-* let a
50,000-fold collapse pass a broken <= 1e-6 key; the on-disk artifact
carries WINDOW-NOT-RIGID from the start). PROGRAMME CLOSE (this cut):
door n°2 survives both attacks - published content has no neutral
states (V259/V268), extended content's neutral states are
continuum-embedded (V366/V367/V368). A stable composite REQUIRES
dynamics beyond the published clock step (an attractive channel
producing an edge below the continuum) - the honest remaining shape
of the interacting-sector problem. V369 (grid convergence of the min
gap) booked as optional hygiene only; it cannot change the sign of
the result (denser continuum can only shrink gaps).

## V369 pair-binding onset at the certified layer interaction (booked BEFORE run; interaction programme, first cut)

Follows the reality-interface cut (door n°2 open condition: an
interaction producing a threshold below the continuum edge). The
certified candidate interaction EXISTS: the bond-current Hamiltonian
H = J.L of Section 7.1 / V101 (units of eta; published numbers E0 =
-5802.398188196636 eta, gap 42.09063148320456 eta, re-read from V101's
artifact at runtime as anchor E-eta). The Casimir identity (booked as
T0, machine-checked on 6 random E4 states, <= 1e-9):
    J.L = 1/2 (C2_tot - C2_1 - C2_2 - C2_L) = 1/2 (C2_tot - 154)
with C2(32) = 45 (per particle, Schur on the irreducible 32; V259
towers 45/85/117 on 32 x 45 certify it) and C2(45) = 64 (V366 E1
second level). Consequences booked as anchors BEFORE any scanning:
  T1 kernel: W^dag (J.L) W = -77 I_2 exactly to 1e-12 (composite shift);
  T2 channel tower: nonzero C2_tot channels of E4 from eigsh (k=400,
     cupy einsum matvec) with the V366 anchor (first nonzero 36,
     multiplicity 2 at its head) and V368 anchor (eta=0 torus min gap
     2.0e-05 reproduced by the SAME machinery);
  T3 the 18-eta slope: binding separation composite vs the c=36 channel
     edge = 18 eta per Casimir arithmetic - machine slope over the eta
     scan must match 18 within 0.5 or CASIMIR-SLOPE-DRIFT (red).

MACHINE (E4 fiber 46080 = 32 x 32 x 45, then sublattice 64 for bands):
(i) cupy matvec C2_tot via einsum on (32,32,45) tensors; eigsh k=400
SM -> channel tower {0, c1=36, c2, c3} with multiplicities;
(ii) channel projectors P_c for the THREE lowest nonzero channels by
shift-invert eigsh (sigma = c +- 0.5 bracket), orthonormality P^dag P
= I to 1e-10, C2_tot P = c P residual <= 1e-8;
(iii) true channel internal bands: phases of P_c^dag (U32 (x) U32)
P_c (eigh of the compressed block) - NOT single-particle phase sums
(the generous union bound is known to cover the circle and is used
for context only);
(iv) eta scan {0.05..1.00 step 0.05} + refinement: composite band
{beta_a(k1)+beta_b(k2) + omega - 77 eta} vs the true c=36 band
{beta + ph_int(36) + (36-154)eta/2} over 12 booked K of the 6^3 torus
(seed 0x51 + K=0 + V368's worst-K); omega rows as V367/V368;
omega=0 stays a control (excluded from BOUND reading).
Verdict words: BINDING-ONSET-<eta*> (the composite band exits the
TRUE nearest-channel band for eta > eta*, with the 18-eta Casimir
slope certified; eta is declared-not-read in the published theory,
so the O-card stays open and the sentence is conditional);
NO-ONSET-IN-SCAN (gap never opens over the scan - honest red for the
interaction programme's first cut); CASIMIR-SLOPE-DRIFT / ANCHOR-
DRIFT (machine failure, no interpretation). Tier 2, residual, zero
fitted. Cost: cupy eigsh + compressed blocks, ~10 min.

HONEST BOUNDS booked with the run: contact-pair placement (relative
kinematics sampled, not full Fock); channels above the three lowest
NOT compared (their bands sit higher in energy by the Casimir
arithmetic - monotone in c - so the nearest channel is the binding
edge; monotonicity booked as T2c: c1 < c2 < c3 and shifts monotone);
one-flux sector; eta scan not a claim about the published eta (which
is declared). V370 (dispersion of the bound band = lifetime) booked
only if ONSET fires.

### V369 post-run note (2026-11-02) - booked identity FALSE, honest negative
V369 ran as booked and the machine declared V369-ANCHOR-DRIFT
(artifact kept verbatim, pass false). Cause is mine, not the
machine's: the booked T0 identity J.L = 1/2(C2tot - C21 - C22 - C2L)
omits the current-current contraction - the correct decomposition is
1/2(C2tot - sum of blocks) = J_1.J_2 + J_1.L + J_2.L. The -77 / -59
/ 18-eta arithmetic belongs to a DIFFERENT operator (total cross
term), not to the published bond-current Hamiltonian (V101). The
flat eta-profile (gap ~1e-5 at every eta, slope -0.000) additionally
kills the additive-phase picture: e^{i eta J.L} does not commute
with the step's internal layer. Sound reusable pieces: E4 channel
tower {36 x 40, 64 x 149, 108, 101}, projectors to 3e-15 (Gram
-1/2, V366 precedent again), V368/V101 anchors green, ARPACK cache.
Door n°2 status after the first interaction cut: the Casimir
shortcut to a threshold is DEAD; the correct binding question is a
non-commuting unitary spectral problem (not claimed today). V370
booked next: exact first-order singlet coupling under the bond-
current operator - leakage ||(1-P_ker) J.L W|| and its channel
decomposition.

## V370 singlet leakage under the certified bond-current operator (booked BEFORE run; interaction programme, second cut)

The corrected decomposition of V369 is itself a booked machine gate
T0: 1/2(C2tot - C2_1 - C2_2 - C2_L) = J_1.J_2 + J_1.L + J_2.L on 6
random E4 states <= 1e-9 (the error of the first cut becomes the
verifiable identity of the second). Gauss on the kernel (Q_a W = 0)
gives L_a W = -(T1_a + T2_a) W, hence
    J.L W = -sum_a (T1_a + T2_a)^2 W = -(90 + 2 J_1.J_2) W,
so the leakage of the neutral pair under the published bond-current
operator is EXACTLY controlled by the current-current contraction on
the kernel. MACHINE (seconds; V369 ARPACK cache = kernel + the four
lowest channels, 400 columns exactly):
  gates: W kernel dim 2 with ||C2 W|| <= 1e-9 and W^dag W = I;
  M1: B = W^dag (J.L) W (Hermitian 2x2, first-order intra-kernel
      interaction shifts - the corrected replacements for the dead
      -77 numbers), reported with eigenvalues;
  M2: leakage per channel l_c = ||P_c (J.L) W|| for c in {36, 64, 80,
      84} plus the outside-the-four residue l_rest^2 = ||J.L W||^2 -
      sum l_c^2 >= 0 (booked non-negativity);
  total l = ||(1 - P_ker)(J.L) W||.
Verdicts: LEAKAGE-ZERO (l <= 1e-12: the bond-current interaction
PRESERVES gauge neutrality of the pair at first order - e^{i eta
J.L} keeps the singlet subspace exactly, the door stays open but
binding still needs the channel spectral comparison, booked V371) /
LEAKS-<l>-DOM-<c> (l > 1e-12: the published interaction does not
bind the neutral pair - it mixes it into charged channels at first
order, channel weights reported; the naive interaction route to a
physical sector closes) / GATE-FAIL. Tier 2, zero fitted, no new
physics operators - every piece is the certified V79/V101 stack.

### V370 post-run note (2026-11-02) - LEAKAGE-ZERO
Facts as booked: J.L on the singlet kernel = -64 x I_2 exactly (the
adjoint Casimir via Gauss; no first-order multiplicity splitting);
total first-order leakage 7.3e-14 - the neutral pair is invariant
under the published bond-current operator. T0 (the corrected V369
decomposition) fires at 4.6e-15.
Hygiene record (third appearance of the Lanczos-degeneracy family):
the k=400 SM cache has O(1) CROSS-cluster overlap - the first pass
reported channel weight l_36 = 90.51 = 64*sqrt(2) = ||J.L W||
against a total of 7.3e-14: pure contamination. Booked correction:
in the zero-leakage branch every channel weight is bounded by the
projector-free ||eps|| and weights are declared not-resolved-and-
moot; the leaking branch (not reached today) cross-orthogonalizes
clusters with a rank assert. Fresh-checkout fallback (V370 builds
the cache itself if V369 ran none) verified end-to-end by deleting
the cache and re-running.
Cache-basis dependence booked: total leakage 7.3e-14 (first
eigenbasis, built inside V369) vs 8.8e-14 (rebuilt basis) - the
1e-13 band is the floor of the cached decomposition (degenerate
kernel re-rotation), not a physical drift; the -64 I_2 eigenvalues
are identical to 1e-6 printout on both bases.
Scope held: single-bond restriction, first order only; the deformed
step e^{i eta J.L} U_step is NOT closed on today's subspaces, so no
threshold claim in either direction - the interaction neither binds
nor breaks the pair at this cut.

### V371 shape booking (2026-11-02; NO run booked today - machine beyond scope)
The honest remaining question for door n°2: spectrum of the deformed
step exp(i eta sum_b J_b.L_b) . U_step on the two-quanta + flux
sector of the TORUS (the lattice sum, not one bond). What today's
machine CANNOT honestly bound: (i) the contact/kernel subspace is
not invariant under the deformed step beyond one bond; (ii) the SM
cache basis has cross-cluster contamination; (iii) additive-phase
bookkeeping is invalid (non-commuting factors - V369 lesson).
What a V371 run requires: Krylov (Lanczos on the unitary) spectral
scan of the deformed step in a fixed total-momentum sector,
eta-scan of isolated-eigenvalue continuation from the composite
band, cross-orthonormalized channel sectors by exact diagonalization
of C2_tot in reduced blocks. Cost booked when the cut is taken:
multi-RHS Krylov on 2.9M-dim fibers, one programme-day, GPU-heavy.
Trigger: only as its own booked session - NOT to be improvised
inside a fold.

## V372_torus_singlet_tower_census - BOOKING (2026-11-02, run following,
tier 2, residual; category dynamics_spectrum; suite 284)
Property: does the two-quanta + flux composite class of the N=2^3 torus
carry GAUSS-SINGLET channels beyond the single flux line - the
prerequisite the V371 deformed-step programme needs.
Reason for the census (booked consequence of V370, no run needed): on
the torus with one flux line (ell=1), the only flux link touches the
matter and every other bond is empty, so Sum_b J_b.L_b = J_b0.L_b0 =
-64 on the E4-type kernel (V370 verbatim) - the deformed step is the
PURE PHASE e^{-64 i eta} on the whole ell=1 singlet tower: it cannot
move a threshold. Binding therefore REQUIRES the ell=2 sector. There
the interior site carries no matter, Gauss forces (L_f1 + L_f2)|psi> =
0, and the touching contraction becomes (J_A - J_B).L_f1 - a different
operator with a different Casimir bookkeeping (the V369 sign lesson
reappears with a flipped cross term): not a phase.
Machine (all booked before the run):
 A1 anchor: ker(32 (x) 32) = 2 (dense eigvalsh, 1024).
 A2 anchor: ker(45 (x) 45) = 1 (dense eigvalsh, 2025; adjoint square).
 A3 anchor: ker(E4) = 2 (V366 verbatim, LinearOperator eigsh, 46080).
 T  torus counts (N=2^3): 8 sites, 12 undirected bonds, 120 composable
    adjacent bond pairs (8 sites x C(6,2)), 3 bond orbits under
    translations - exact enumeration, structural gates.
 M  ell=2 fiber: 32 (x) 32 (x) 45 (x) 45 (dim 2,073,600), charges
    Q^a = T_a(g1) + T_a(g2) + L_a(f1) + L_a(f2), C2 = sum Q^2 by
    cupy-einsum matvec (V370 fallback pattern), eigsh which='SM'
    k=12, kernel block Gram -1/2 fix with rank assert (house
    hygiene), residual gate ||C2 W|| < 1e-8, W^dag W defect < 1e-12.
Verdicts: anchors break -> MACHINERY-DRIFT (exit 1). ker_ell2 >= 1 ->
TORUS-SINGLET-TOWER-ALIVE (V371 licensed on the ell<=2 sector, count
booked for the Krylov assembly). ker_ell2 = 0 ->
TORUS-SINGLET-TOWER-LINE1-ONLY (with the pure-phase consequence the
published interaction CANNOT open a threshold through any ell<=2
singlet channel: V371 closes at ell<=2; the next honest steps are
ell>=3 or a content extension - no run today either way).
Scope limits booked: single quanta pair (fermions, one flux line),
N=2^3 finite size (a census, not a thermodynamic statement), NO
dynamics computed here (no step projection, no eta scan - that is
V371 proper), adjoint flux lines only.

### V372 machine correction (booked before the run, after script
inspection)
The Weyl weight-lookup formula as first drafted (m_V(w rho - rho)) is
WRONG (it fails its own anchor ker(32x32)=2 on paper) - replaced
BEFORE any run by the exact dense-Casimir cluster method: decompose
32x32 (1024, eigvalsh), 45x45 (2025, eigvalsh) and 16x16bar (256,
from AC.weyl_generators) by C2 eigengroups; the irrep table is read
from the cluster dims (16x16bar = 1+45+210 fixes the 45 and 210 C2
values; 45x45 clusters must have dims {1,45,54,210,770,945} (the 45'/1670 first sketched from the weight algebra are WRONG - an inspection run of the cluster machinery gives 45x45 = (1+54+945) sym + (45+210+770) antisym; shape gate uses the machine fact));
mult(1; V tensor F) = sum over C2-matched shared clusters of the copy
counts; an unidentified shared value -> EXT-GROUP-MISMATCH. Anchors
gated: ker(32x32)=2, ker(45x45)=1, E4=2, E1=1, E2=1 (V366).
Torus anchor corrected: N=2^3 has 8 sites, 12 undirected bonds,
24 adjacent bond pairs (8 sites x C(3,2) - undirected incidence, NOT
the directed count 120 first drafted), 3 bond orbits. The ell=2
cupy-eigsh-on-2.07M path is NOT used (the cluster method is exact and
finishes in seconds; GPU stays free for V371 proper).

## V371_deformed_contact_channels - FORMAL BOOKING (2026-11-02, run
following; tier 2, residual; category dynamics_spectrum; suite 285).
This is the executable first cut of the V371 shape booked earlier
today; the full torus lattice-sum spectral question (multi-line
Krylov) stays booked and unrun.
Property: under the published bond-current deformation e^{i eta
Sum_b J_b.L_b}, do the certified singlet channels acquire
channel-dependent phases that RE-OPEN a gap of the pair against the
two-quanta continuum at the contact cut - the door n°2 opening
condition, decided at the scope V372 licenses.
Machine (all booked before the run):
 Fiber: pair at sites (A,B) on the contact bond b0; flux line of
 length ell=1 (b0) or ell=2 (A-m-B via interior m), orientation
 along the path; charges T1 at A, T2 at B, L_f1 at A and m, L_f2 at
 m and B with interior Gauss (L_f1 - L_f2)|psi> = 0 at m.
 T0 anchor: the constructed charge set has EXACTLY the V372 kernel
 counts: dim ker C2(ell=1 fiber) = 2 (V369/V370 basis rebuilt),
 dim ker C2(ell=2 fiber, 2.07M) = 6 - the orientation convention is
 gated by this count; mismatch -> MACHINERY-DRIFT, no run.
 A  exact bases: ell=1 kernel from the V366 E4 machinery (2 vectors,
    46080); ell=2 six vectors by Haar-averaged intertwiners:
    X <- (1/N) sum_g rho_M(g) X0 rho_F(g)^T, N = 200 random group
    elements (products of exp(i theta_a T_a), fixed seed), pairings
    (trivial x trivial: 2; 45 x 45: X 90x45; 210 x 210: X 420x210);
    gates: max_a ||M_a X - X F_a|| <= 1e-8, ||C2 v|| <= 1e-8 on
    every assembled fiber vector, Gram defect <= 1e-12.
 B  deformation spectrum: D2 = W^dag (J.L sum over path bonds) W
    (6x6 hermitian; explicit charges, einsum on the 2.07M fiber,
    cupy GPU 1) and the 2x2 ell=1 block which must reproduce
    -64 I_2 (V370 anchor in-repo). Facts booked: d_1..d_6.
 C  window scan at the V368 torus kinematics: dressed channel
    phases theta_k(eta) = theta_pair + eta d_k against the contact
    continuum bands of V367/V368 (anchors: gamma = 1.006652 at
    contact; torus min gap 2.017083e-5 over the 216 K of the 6^3
    grid, seed 0x51 - machinery re-used, continuum recomputed
    independently in-script). eta grid [0,1] step 0.02.
    window(eta) = min over K and k of dist(theta_k(eta),
    band(K)). GATE BOOKED NOW: window re-opened iff
    max_eta window(eta) >= 0.05 rad/step.
 Verdicts: THRESHOLD-REOPENED-CONTACT-CUT (window >= 0.05 at some
 eta) / NO-THRESHOLD-CONTACT-CUT (honest negative: the deformation
 phases never clear the continuum at this cut) / MACHINERY-DRIFT.
 Scope limits booked: contact geometry for the pair (the V367
 invariant-subspace scope); one flux line, length <= 2; channel
 PHASES only - hopping-induced channel mixing stays out (V368
 already showed hopping collapses the free window); eta first-order
 phases only; no ell>=3, no multi-line, no torus hopping of the
 line. These limits are the difference between this cut and the
 full V371 lattice-sum question, which remains booked unrun.

## Roll-up: the reality-interface phase V366-V372 (2026-11-02, phase end)

Seven cuts attacked door n°2 of CONJUNCTION-DEAD - can the certified
clock step ever carry a physical (gauge-neutral, isolated) asymptotic
state? - and the honest sum closes the phase with the door standing:

- V366 SINGLET-ADMIT E1:1-E2:1-E3:0-E4:2-E5:9: the certified gauge
  action DOES admit singlet matter for enlarged content (conjugate
  charge pairs necessary). The published emptiness of the physical
  sector is a CONTENT theorem, not a property of the intertwiner
  machinery - the machinery was always able to build neutral states.
- V367 COMPOSITE-GAMMA-BOUND-1.007: the singlet pair is an invariant
  subspace of the certified clock step, isolated AT CONTACT (1.007
  rad), and the clock's internal layer acts as the identity on it
  (4e-15).
- V368 WINDOW-NOT-RIGID-2.0e-05: that isolation disappears at generic
  torus momenta - the fixed-K continuum sweeps all bond phases; the
  correct condition is a THRESHOLD (energy below the continuum edge),
  not a phase offset. A stable composite needs dynamics beyond the
  published free hop + Gauss.
- V369/V370: the certified bond-current coupling was brought to the
  pair sector - V370 LEAKAGE-ZERO: J.L acts on the ell=1 contact
  kernel as -64 I_2 exactly (leak 8.8e-14): the certified interaction
  is a pure phase on the minimal channel.
- V372 TORUS-SINGLET-TOWER-ALIVE: the exact cluster census licenses
  the contact cut (ell=1 pure phase x 24 bonds, ell=2 six channels x
  24 pairs on N=2^3; five anchors, shape gates, trace ratio 2.000000).
- V371 NO-THRESHOLD-CONTACT-CUT: the bond-current deformation is NOT
  a pure phase at ell=2 - it splits the six channels into the bands
  {-96, -96, -64, -64, 0, 0} - but it never re-isolates: the
  deformation displaces the continuum by up to 72 rad/step against
  the 2.0171e-05 rad/step free isolation (V368 anchor reproduced to
  1e-9), and eta >= 2*pi/72 = 0.087 zeroes every window by smear
  alone. Door n°2 does not open at the contact cut.

The full lattice-sum version of V371 - the pair delocalized over the
torus, the flux line a dynamical degree of freedom, a genuine
eigenvalue problem below the continuum edge - stays BOOKED, UNRUN:
its verdict word is not on record and no negative is claimed for it.

Machinery lessons booked along the way (no claim moved by any of
them): the kernel operator is the diagonal Casimir sum_a (J_a^tot)^2,
and a "linear collapse" to (sum_a J_a)^2 is a DIFFERENT operator with
a dense zero spectrum (V371); fast coupling evaluation needs the
IN-SCRIPT operator self-test (coupling form == loop-anchored form,
5.7e-16) - one permutation typo silently changed the operator; a
small eigsh residual does not certify the RIGHT kernel - the
fiber-span J.L-invariance leak gate (3.2e-3 wrong span vs 8.8e-14
anchor) is what catches it (V371); the SA-endpoint Lanczos needs
tol >= 1e-9 above the roundoff floor plus Rayleigh-Ritz polish and a
partial-Ritz fallback to land all six vectors; and a booked verdict
word whose key number is the wrong observable (V368's "broken" at
1e-6) must be corrected before the artifact ships.

# The beyond-step programme: pre-registrations V373-V376 (2026-11-02,
all four booked BEFORE any run; category beyond_step, new; one
machinery used on three openings)

## V373_full_lattice_pair_string - BOOKING (2026-11-02, run following;
tier 1; category beyond_step; suite 286).
Property: can the pair+flux-line sector of the Λ³ torus, treated as a
gauge-invariant full lattice problem at fixed total momentum K, carry
the door n°2 question at all - i.e. does the machinery that V371
booked-but-never-ran EXIST, gated against everything the
reality-interface phase already measured?
Machine (all booked before the run):
 Hilbert space: matter pair of quanta at grid sites (A,B) on the
  N^3 torus, separation r over ALL grid separations, fixed total
  momentum K; flux sector = ONE straight geodesic flux line of
  length ell = d_geo(r) joining A to B. Booked scope (declared, not
  justified dynamically): the single flux line is the V368 content
  theorem's minimal gauge content, and the straight line is the
  minimal-length representative for separation r.
 Operators imported verbatim: V79 adjoint_link_operators /
  root_generators / structure_constants (the same import V371 uses);
  V101 zero_modes (the intertwiner kernel solver, J_sub W + W L^T)
  generalized to chains of ell+1 links. The kernel DIMENSIONS of the
  chain problem are MEASURED per chain length, never assumed: the
  fusion semiring of the 45-algebra beyond r in {1, 45, 210} is a
  measurement OUTPUT of this experiment, not an input.
 Effective Hermitian generator (booked form):
   H = H_hop (x) I_flux + I_matter (x) sum_{links in string}
       sum_a L_a L_a + eta * (bond-current deformation in the
       certified channel basis).
  The per-link string Casimir coefficient is MEASURED from the
  adjoint representation in-script, not assumed. H_hop is the V368
  free two-quanta dispersion machinery re-used verbatim (seed 0x51);
  the deformation uses the V371 fast coupling form (K_M, K_F, C_MM,
  C_FF, C_MF, pair factor 2) with the V371 in-script operator
  self-test against the 45-loop anchored form (<= 1e-14) - the
  permanent guard from the V371 machinery lesson.
 Gates:
  A1 N=2 contact reduction reproduces V371 verbatim: ker counts 2
    (ell=1) and 6 (ell=2), the deformation spectrum d_k set =
    {-96, -96, -64, -64, 0, 0} exactly, and J.L = -64 I_2 on the
    ell=1 kernel with invariance leak <= 1e-12. The ell=1 kernel is
    RE-DERIVED in-script (no V369 SM-400 cache dependency - the
    fold that regenerated it was killed by user directive).
  A2 eta = 0 free min gap reproduces the V368 2.017083e-05 to 1e-9.
  A3 measured chain-kernel dimensions for ell = 1..3 anchor on the
    V100/V101 values: ell=1 kernel as in V101's three-site chain
    (Gram rank 24, six-fold top degeneracy); the r = 1, 45, 210
    fusion multiplicities at ell=2 match V100 (I_2 = 1, I_3 = 1,
    r = 210 channel constructed, 24-state kernel); ell=3 (I_4 = 6
    from V100) is a prediction OUTPUT recorded for the ledger.
 Verdicts: FULL-LATTICE-MACHINERY-BUILT (all three gates fire) /
 MACHINERY-DRIFT (any gate fails: no spectral claim is made).
 Scope limits booked: one flux line, straight and geodesic; no
 multi-line, no flux-line hopping; the string Casimir term is the
 booked effective generator, not a derivation from the clock step.

## V374_lattice_sum_threshold - BOOKING (2026-11-02, run following
V373; tier 1; category beyond_step; suite 287). The V371 verdict word
"the full lattice-sum version stays booked unrun" is retired by this
entry: the question becomes executable at the scope booked here.
Property: for the V373 operator, is there an eta at which the lowest
singlet-channel eigenvalue E_0(eta, K) drops below the two-quanta+
string continuum edge E_c(eta, K) - a genuine THRESHOLD bound state,
the door n°2 opening condition at full-lattice scope?
Machine: V373's operator, N in {4, 6} (fallback N in {4, 6, 8}
booked; gate: total fiber <= 2e7 states, else declare the
restriction in the artifact); K in {Gamma, the V368 worst-K grid
points}. Lanczos: lowest singlet-channel eigenvalue E_0(eta,K) and
continuum edge E_c(eta,K) (lowest two-quanta+string continuum
state), gap g = E_c - E_0; eta grid [0,1] step 0.02.
 Bound-state criterion booked: REOPENED iff max over (eta,K) of
 g >= 0.05 rad/step AND the lowest state's <d_geo(r)^2> STOPS
 GROWING with N (finite-size test N=4 vs N=6 - a continuum state
 spreads with the lattice, a bound state does not).
 Verdicts: THRESHOLD-REOPENED-LATTICE-SUM / NO-THRESHOLD-LATTICE-SUM
 / MACHINERY-DRIFT.
 Contingencies booked (declared in artifact and run_notes if used,
 never silently taken): if the N=6 fiber exceeds memory, K = Gamma
 only plus string length <= N/3; if chain-kernel dimensions explode
 beyond 2e7 states at ell=3, fall back to ell <= 2 with N up to 8
 and record the long-string sector as UNMEASURED.

## V375_dressed_no_ladder_contact - BOOKING (2026-11-02, booked
before running, executed after V373/V374; tier 2, honest-negative
search; category beyond_step; suite 288).
Property: do the eigenphase differences of the DRESSED generator
(V373/V374 spectra) reproduce any of the published empirical ladder
contacts - i.e. does the beyond-step dressing machinery "know" the
ladder, or is it blind to it?
Machine: take the dressed eigenphase spectra measured by the V373
and V374 artifacts (eigenvalue differences of the dressed
generator, machine-exact Decimal where finite). Contact harvest is
MECHANICAL (the V364 style): the five empirical ladder contacts of
the V52 provenance table (v52_parameter_audit.json rows whose
quantity or note contains "rung" or "ladder" - exactly five rows
booked: the ladder base 1/(16 pi) = 0.019894 (declared), g_obs =
0.1999980 (rung-10 loop-hierarchy clock), the DM-mass coefficient
1.18, n_2step = 3 (two-step seesaw ladder exponent), and the
sub-Planck DM scale rung-10 value (118 GeV row). Harvest count is a
gate: five found or the protocol gate fails.
 Tests: every dressed eigenphase difference against every contact,
 exact-Decimal at relative 1e-12. False-hit bookkeeping booked NOW
 (counting stated): T = (#differences) x 5 single-point tests, each
 false-hit chance <= 2e-12 (relative-1e-12 band on the 2 pi circle
 of continuous phases); a HIT must ALSO reproduce at the adjacent
 eta grid point (the spectrum is re-diagonalized there; the two
 chances are independent) so the per-test chance is <= 4e-24, and
 the booked false-hit chance is T x 4e-24 <= 1e-16 for T up to
 2.5e7 (T is printed in the artifact).
 Verdicts: DRESSED-LADDER-CONTACT-HIT / DRESSED-NO-LADDER-CONTACT.
 The zero-hit branch is the honest negative: pass = True with the
 protocol gates (harvest count = 5, test count, false-hit bound
 stated). A NONZERO hit is a major positive: before committing it
 is re-tested mechanically against the V364 named-value harvest -
 a contact already named would contradict V364; double-named
 artifacts get cross-checked, never silently re-labelled.

## V376_dressed_cp_odd_census - BOOKING (2026-11-02, booked before
running; independent of V373/V374; tier 2, theorem extension;
category beyond_step; suite 289).
Property: does the V52 computed theorem "every CP-odd invariant of
the published QCA vanishes" extend from the published step U to the
DRESSED class {W U, e^{i eta sum_b J_b.L_b} U, and the V373 full
generator} - i.e. is the dressing CP-even to measured exactness, or
does it carry CP-odd content the bare QCA provably lacks?
Machine: the CP-odd invariant machinery is RE-USED, not reinvented:
 cl10/parameter_audit.py cp_odd_theorem (the invariant list: Im Tr
 of every 2- and 3-rotation product (full), the booked 4-rotation
 sample n4_sample = 200000 seed 0, every B-violating
 h_i R_ab h_j R_cd sandwich (full), and the rank-1 real democratic
 Yukawa invariant Im (Y Y^dag)^2) and the CP operator from the SAME
 source (the pseudoscalar w, the Cl(10,0) volume element, w^2 = -1,
 w* = -w). The dressed step replaces U in every invariant; all
 booked invariants are evaluated at eta in {0, 0.02, ..., 1} as
 machine-exact traces (mpmath high-precision; the V52 TOL = 1e-12
 is the booked decision threshold).
 Verdicts: CP-EVEN-DRESSED-CLASS-CLOSED / CP-ODD-FOUND-DRESSED.
 Pass gate for the CLOSED branch: every booked invariant value
 |I| <= 1e-12 at every eta, printed per eta.
 Scope booked: this extends the theorem from the bare root step to
 the booked dressed class; it is NOT a claim that no CP-odd
 operator exists anywhere in the lattice theory (the named
 eta_CP residual of Table 6 stands).

## Roll-up: the beyond-step programme V373-V376 (2026-11-04, phase end)

All four bookings ran; the one machine served four openings and the
full lattice-sum scope of door n°2 closes on measurement:

- V373 FULL-LATTICE-MACHINERY-BUILT (tier 1, suite 286): the
  gauge-invariant pair + straight-geodesic-string machine certified
  end to end - fusion census measured ({0:2, 36:20, 64:90, 84:240,
  96:420, 100:252}; {0:1, 64:45, 80:54, 96:210, 128:945, 144:770}),
  chain kernels measured per length (2, 6, 22, 12 singlet sector),
  I_4 = 6 as booked-prediction output, V371 kernel/d_k/span and the
  V100/V101 24-state glued anchor reproduced from measured fusion
  data, string Casimir flat at 64.000000/link, V368 free gap to
  0.0e+00, two solver routes agreeing, V101 top degeneracy 6.
- V374 NO-THRESHOLD-LATTICE-SUM (tier 1, suite 287): full lattice
  sums at N=4/6 (continuum 14155776 states per K, inside the booked
  2e7 gate; K = Gamma + gap_order top-12 declared per grid; ell<=2
  value sets with the ell=3 contingency recorded - its value sets
  coincide). The nominal window g = E_c - E_0 opens (0.192971 at
  eta 0.48 resolved; 0.887356 coarse collapses under refinement),
  but the booked conjunction needs localization and the finite-size
  discriminator fails cleanly: <d_geo^2> 4.5 -> 9.5 (ratio 2.111,
  exact uniform-pair value). Door n°2 stays closed at full
  lattice-sum scope; the V371 "booked unrun" clause retires.
- V375 DRESSED-NO-LADDER-CONTACT (tier 2, suite 288): mechanical
  harvest = exactly the five booked contacts; T = 5975 exact-Decimal
  tests (1195 dressed differences x 5, adjacency retest booked and
  implemented; V373-family adjacency adapted to second-route
  reproduction, DECLARED): zero hits, false-hit chance 2.39e-20.
  The dressing is blind to the empirical ladder.
- V376 CP-EVEN-DRESSED-CLASS-CLOSED (tier 2, suite 289): V52 census
  re-used on the dressed class {W U (W = v101_gap x G10 = 33.057905
  artifact-gated), e^{i eta J.L} U (measured scalar dressing,
  sum_a T_a T_a = 45.000000 I at 0.0e+00 - declared reduction), V373
  generator}: every booked invariant <= 1e-12 at all 51 etas (worst
  2.66e-15), mpmath-40d tier agreeing to 2.9e-16, eta=0 column
  reproducing V52 exactly. Channel-space phases recorded as data,
  not invariants; the named eta_CP residual stands.

Contingencies declared where taken (V374 ell<=3 fiber gate, K-set
interpretation; V375 V373-family adjacency adaptation; V376 basis
convention record and eta-independent verification tier). Machinery
lessons recorded in the run notes and reference (einsum label
collision -> non-algebraic-clean spectra as detector; rounded channel
keys; sorted-multiset comparison for sector-blocked spectra; per-block
npz checkpoints with sha transplant discipline; gates must pin basis
invariants, not representation conventions; verification tiers on the
SAME products).

Fold state: the full incremental gate at suite 289 (library hash moved;
the whole suite re-runs) reached ~100/289 before being stopped by
directive: V374/V375/V376 PASS verbatim in the fold, V373 CRASHED on a
cupy context-init OutOfMemoryError (16,384 B, 0 allocated) at t=592 s
- a transient GPU-contingency at jobs=2 with the inference server
resident, before any physics of the run; the standalone V373 artifact
is green and its fingerprints left unrecorded, so the gate re-runs it
alone next fold with the block checkpoints warm. Fingerprints stay at
the suite-285 library state. V369's committed ANCHOR-DRIFT honest
negative (pass=False since suite 282, machine green) is the single
FAILED verdict in the 289-experiment roll-up, unchanged.

# The module-replacement tribunal: V377-V378 (2026-11-04, booked
BEFORE any run; category root_rebuild continuation; the ten-to-one
attacked as a class question)

## V377_module_replacement_tribunal - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 290).
Property: does any BUILT module replacement move any contact of the
ten-to-one? The replacement class is exactly the set of module/
internal-layer constructions the repository has actually built and
committed (no paper candidates): the committed 36-level cell; the
altcell signature family (Cl(10-q,q), the q=0 member the V117/V312
quaternionic root); the 2I engine internal layer (V337); the star
module (V244/V249); the pentagon engine (V356-V360); the rhombic cell
(V361-V362); and the rung lattice itself (n_A, n_B free). Every
contact is the published number, read from the committed artifacts
where the member already measured it and recomputed in-script where
the member is the committed cell or the rung lattice.
Machine (all booked before the run):
 T0 internal self-tests (no contact): LOOP = 1/(16pi) = 0.0198943...
  from cl10.u10pred to 1e-15; rung phases loop^n (cell.rung_phase)
  match the config rungs (3, 10) to 0.0; the committed cell step at
  g = pi/4 and Gamma has the eigenphase multiset
  {+-(pi/4) x 10, +-(3pi/4) x 6} + {+-loop^3, +-loop^10} to 1e-12;
  unitarity defect of cell_step_np <= 1e-12.
 T1 closure census of the COMMITTED cell (new): for m = 1..24 the
  full 36-level step U_cell(g=pi/4)^m never returns to I - the
  booked expectation: |U^m - I|_max >= 1.4e-5 for every m (the
  singlet phases 2*loop^3 and 2*loop^10 per 12 stages are the
  obstructing phases; the spinor block alone closes at order 8):
  the closure axiom is confined to the spinor block BY MEASUREMENT,
  for the published rungs.
 T2 rung-lattice closure census (new): over (n_A, n_B) in 1..30 x
  1..30, the singlet block can close only if some m <= 10^6 has
  ||m*loop^n / (2pi)|| <= 1e-9 (continued-fraction convergents);
  booked reading: NO rung in the range admits closure - the least
  residual over the whole grid is booked to be > 1e-5 (measured
  first-run; the gate is > 1e-9 as booked) - so no rung choice ever
  extends the closure class; the closure clocks stay {pi/2, pi/4,
  3pi/4} for every rung choice (declared: the blkdiag structure
  makes g-independence exact; the census measures the rung side).
 T3 the contact ledger (artifact-gated): for every member read the
  committed artifact fields and compute the contact shifts:
  - clock contact: min over closure clocks |g - g_obs| at
    g_obs = 0.1999980 (V206/V233 booked): committed class -> 0.5854;
    p=0 quaternionic clocks pi/2, pi/4, 3pi/4 (V117 fields) -> same;
    star closures and dial distances from V249 fields (0.585, 0.285);
    2I: 0 of 9 angle-grid points keep structure (V337 field).
    BOOKED: the minimum clock-contact distance over the class is
    achieved by NO member below 0.28.
  - hierarchy-deficit contact (O8): the certified reads at the
    closure clocks: 2.7095/2.711 (V209/V233 fields) and the
    quaternionic reads {2.2994, 2.7081, 2.7057} (V312 fields)
    against 2.076564111173832: BOOKED: every member's deficit is
    in [0.22, 0.31] - the +30% deficit is module-invariant;
    the star and pentagon members carry NO ladder crossing object
    (V249 eta-twist moves nothing; V359 NO-CROSSING field) - a
    member either has the ladder and pays the deficit, or lacks
    the ladder object entirely.
  - content contacts (O3/O4/O5/O7/O10): V117 T3 antiunitary
    deviation 1.414 at every signature; V312 odd-n selection rule
    survives on p=0; V342 icosahedral-cell chiral test fails
    (557 nodes multi-orbit, cones anisotropic 0.023-0.354); V350
    backreaction fixed-point rigidity 4.34e-14; V376 dressed class
    CP-even. BOOKED: no member changes any of these by a booked
    contact amount (each artifact field gated to its published
    digits).
 T4 decision (booked word): if every contact column is stationary
  across the class, the verdict word is
  MODULE-REPLACEMENT-CONTACT-STATIONARY, with the NAMED SURVIVING
  STRUCTURE booked verbatim: the crystalline dilemma - the ladder
  objects exist only on crystallographic shift modules (V357 void
  on the pentagon), and clock-free closure exists only off them
  (V334 escape branch, V356 realizes it on I_h) - so within the
  built class the clock can be born only at the cost of the ladder
  that reads it. Any module-side solution must break that dilemma,
  not pick a member.
Contingency (declared): if an artifact lacks a booked field, the
member's column is recorded VOID and named; the verdict word then
carries the clause "over the members with fields"; no field is
substituted. If T2 finds any rung with residual <= 1e-9, the rung
pair is recorded as a NEW closure candidate and V377's verdict word
carries MODULE-CLOSURE-RUNG-FOUND with its measured m - booked
before the run as the only positive branch.

## V377 booking correction (pre-run, declared; 2026-11-04, before any
run of V377)

Two booked constants in the V377 booking are corrected BEFORE the run
(no result exists yet; nothing else in the booking moves):

- T3 hierarchy column: the booked deficit interval [0.22, 0.31]
  omitted the quaternionic pi/2 read 2.299444 from the V312 artifact.
  Corrected booked reading: the absolute deficit |coef - observed|
  over every member that carries a ladder object lies in
  [0.2229, 0.6330] (closest 0.2228796480336257, the V312
  closest_read field itself; farthest from the certified pi/4 read
  2.7095258500787756 minus 2.076564111173832).
- T2 rung census: as booked ("residual > 1e-9") the grid fires
  trivially - loop^n itself falls under 2pi x 1e-9 at n >= 5, an
  epsilon-closure by DECOUPLING (the level carries vanishing phase
  per stage), which is not a new clock and closes nothing of the
  spinor block. Corrected booked gates: (i) NO rung admits an exact
  closure (the continued-fraction residual is strictly > 0 for every
  (n, m <= 10^6) on the grid); (ii) the epsilon-closure threshold is
  booked at n* = 5 (least n with loop^n / (2pi) <= 1e-9) and is
  recorded as data with the decoupling interpretation; (iii) the
  closure CLOCK class {pi/2, pi/4, 3pi/4} is measured unchanged when
  the cell is built on the epsilon-closed rungs (5,5) and (7,7).
  The booked positive branch MODULE-CLOSURE-RUNG-FOUND therefore
  fires only on exact closure; the epsilon branch is booked as
  MODULE-CLOSURE-RUNG-EPS-DECOUPLED and moves no clock.

Amendment to the same-day correction (run had begun; ONE gate line
affected, reported post-run per house convention): the T3 closest
bound was itself a rounded constant (0.2229) and the exact artifact
value 0.2228796480336257 missed its own rounded gate. Corrected
gate (booked semantics - the closest read IS the V312 closest_read
field, gated bit-exactly): closest deficit == 0.2228796480336257
(|d| <= 1e-12) and largest deficit == 0.6329617389049436 (the
certified pi/4 read 2.7095258500787756 minus 2.076564111173832,
bit-exactly); the interval statement becomes [0.22288, 0.63296].
No physics interpretation changed by the amendment.

## V378_crystalline_dilemma_witness - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 291).
Property: V377's named surviving structure - the crystalline dilemma
- put to a direct witness test on the repository's ONLY clock-free
closure engine (the 2I engine, cl10/supercell_2i, the branch module
booked by V337): does ANY member of the 2I engine class carry the
full conjunction CLOSURE (at GAMMA, order <= 240) + gapless M-points
(min |phase| <= 1e-9 over all M-points) + the published R cluster-
size signature + internal mirror <= 1e-10 + unitarity <= 1e-11?
Booked class (all booked before the run): omega = pi*j/5 for
j = 1..9 (the 5-fold angles at which the binary-icosahedral element
has finite order; V337 used j in {1,2,4}); beta in the four readings
{SC.G10 (published), pi/5, 2*pi/5, 1/5 (the empirical transport
pin)} - 9 x 4 = 36 points. Machinery: the V337 helper set
(phases, cluster sizes, unitarity, closure_k cap 240 tol 1e-8)
re-implemented in-script on HOST numpy (V337 ran cupy GPU 1; the
class is 256-dimensional, host-exact).
T0 self-tests (booked anchors): (i) the PUBLISHED step re-measured
through this machinery: M-gap < 1e-9, R signature equal to the V337
T0a pub_R_cluster_sizes field read from the artifact, mirror
<= 1e-12; (ii) the three V337 grid points inside the booked class
(j in {1,2,4} x beta in {pi/4, pi/5, 2pi/5}) reproduce the V337
artifact rows at 1e-8 (backend cross-check, host vs the V337 cupy
values).
T1 (booked reading, decided BEFORE the run): NO point of the class
satisfies the full conjunction: the points with gapM = 0 fire at
beta in {pi/5} where the R signature fails (the V337 rows show
gapM = 0.0000 exactly at (2pi/5, pi/5) and (4pi/5, pi/5), R False),
and the points holding the published R signature (beta = G10 = pi/4
class) all carry gapM >= 0.1047 (the V337 row minimum at the class
j = 2, beta = pi/4). BOOKED: min over R-holding class points of gapM
= 0.10471975511965977 (== pi/30 exactly? measured and reported -
the booking fixes the V337 row value 0.1047 as the floor to be
reproduced, not an exact claim).
T2 clock side (booked): every class point closes at GAMMA at order
<= 10 (the 5-fold structure forces it; measured closure_k recorded
per point) - closure is clock-FREE over the whole class, witnessed.
Decision (booked words): 0 of 36 conjunction points ->
CRYSTALLINE-DILEMMA-WITNESSED (the clock can be born only by
killing the gapless M or the published R multiplets: no built
module carries both sides); any point with the full conjunction ->
DILEMMA-BROKEN with its (omega, beta, closure_k, gapM) recorded -
booked as the only branch that could carry a module-side solution.
Contingency (declared): the R signature of the PUBLISHED step
re-measured in T0 (pub_R_cluster_sizes) is the reference; if the
host re-measure differs from the V337 GPU field, the host field is
used verbatim as the reference and the difference recorded as a
backend artifact (V337 T0c precedent).

## V378 post-run note (T2 correction; 2026-11-04, after the run)

The booked T2 line ("every class point closes at order <= 10") is
true of the INTERNAL 2I layer (measured 36/36, max order 10) but not
of the full engine step: the step closes at 27/36 points - exactly
the rational-multiple-of-pi beta columns - and NEVER at beta = 1/5
rad (a transport phase of 0.2 rad/stage is incommensurate with 2pi,
closure_k None to cap 240). The clock-free witness stands where the
booking intended it - the internal layer carries the closure with no
clock cost at every 5-fold angle - and the measured full-step split
27/36 (not 36/36) is booked as the corrected statement. No other
gate moved; DILEMMA-BROKEN branch: 0 of 36.

## Roll-up: the module-replacement tribunal V377-V378 (2026-11-04,
phase end)

The ten-to-one attacked as a class question, all booked before any
run, every booking executed:

- V377 MODULE-REPLACEMENT-CONTACT-STATIONARY: across the complete
  class of BUILT replacements (signature family, 2I engine, star,
  pentagon, free rung lattice 1..30) NO contact of the ten moves -
  clock column >= 0.585398 everywhere, hierarchy deficit in
  [0.22288, 0.63296] for every ladder-bearing member, content
  columns verbatim. New measurements: the committed cell never
  closes (closure axiom confined to the spinor block: min defect
  6.30e-05 at m = 8); no rung closes exactly (CF floor 1.46e-52);
  rung epsilon-closure (n* = 5) is decoupling, not a clock.
- V378 CRYSTALLINE-DILEMMA-WITNESSED: on the only clock-free engine
  built, closure is free (internal 36/36 at order <= 10; full step
  27/36 - rational-pi betas only) but gapless M and the published
  R multiplets never co-occur (gap floor 0.104720 where R holds);
  the booked DILEMMA-BROKEN branch - the only module-side solution
  carrier - fired 0 of 36.

The phase closes the MODULE SIDE as a search and names what a
module-side solution would have to be: a single module carrying a
clock-free closure core, a crystallographic walk ladder, and the
published R multiplets at once. Every built carrier so far holds at
most two of the three legs (2I: closure+R, no gapless M; pentagon
grid: gapless M, no ladder crossing; cubic committed: ladder+M,
closure consumes the clock). Booking corrections: two pre-run
(V377 deficit interval and rung-census semantics) and one post-run
(V378 internal-vs-step closure split) - all declared in place, no
interpretation moved by any of them. Fold state: V377/V378 verdicts
green standalone; the fingerprint gate will run them at the next
fold (new experiments, library state unchanged).

## V379_hierarchy_coefficient_shadow_tribunal - BOOKING (2026-11-04,
run following; tier 2; category root_rebuild; suite 292).
Property: the sharpest unnamed contact found by hand after the
tribunal - the hierarchy coefficient observed, 2.076564111173832
(V233 field), sits 3.8e-05 from 3*exp(-1/e) - INSIDE the booked
1e-04 observation band, with zero fitted parameters. Is this a
birth channel or a shadow? And what is the exact shape of the
V378 near-miss (the R leg)?
Machine (all booked before the run):
 T0 web anchors (artifact-gated): observed 2.076564111173832 and
   uncertainty 0.0001 from the V233 readings field; the certified
   crossing read 2.7095258500787756 (V209/V233); the certified
   crossing F and rho from the V312 ladder_runs fields (read at
   run, named); phi from the V337 T0b golden-ratio field.
 T1 the trap reading (mpmath 40 dp): A = 3*exp(-1/e) computed to
   40 dp; |A - observed| measured; BOOKED: the deviation is inside
   the booked band (3.8e-05 < 1e-04) - the hit is REAL as a
   numerical fact and gets recorded, whatever T3 says.
 T2 required-exponent census: n_req = ln(observed/3) /
   ln(16*pi*F_cert) for the certified F reads (the three reads of
   Section 12.4: published 0.019186391026184392, certified crossing
   0.01968515418112789, float64 pair); BOOKED readings: published
   read gives n_req ~ 10.21 (the fitted rung-10 reading, known);
   the certified crossing gives n_req ~ 34.8-36.2 - and GATE: no
   web integer (the module dimension 36, the rungs 10/17/208, the
   walk lengths, the vertex invariants I_n in {1,1,6,23,141}) sits
   within the n_req uncertainty of any read - recorded, not
   asserted: every |n_req - integer| is printed and gated > 0.05.
 T3 the web-membership tribunal on e: X = exp(-1/e) (and 1/e)
   tested against the certified web basis {16pi, 1/pi, phi, loop,
   F_cert reads, rho = 0.9780302704295731, and the structural
   integers 3, 30, 32, 36, 45, 64, 248} by the V340 pattern -
   pslq([1, x, y, x*y], maxcoeff 200) over every pair and every
   singleton; BOOKED: no relation is found (all pslq outputs None
   at maxcoeff 200 and the pair products) - e-type transcendentals
   are OFF the web: the 1.8e-05 hit has no mechanism and is
   declared a SHADOW; the booked decision word is
   E-SHADOW-COEFFICIENT-NO-WEB-CONTACT, with the fallback
   E-WEB-CONTACT-FOUND (booked, declared) if any relation fires.
 T4 the R near-miss map (measurement, gates): the four 2I gapless
   points of V378 ((omega/pi, beta) = (0.4, 2pi/5), (0.8, pi/5),
   (1.2, pi/5), (1.6, 2pi/5)) re-measured in-script; the R
   cluster signatures recorded; BOOKED (from the pre-run probe
   recorded in the run notes): the signature is [16, 48] at every
   point - the 48-cluster of the published [32, 48] SURVIVES and
   the 32-cluster splits exactly dyadically 16+16 - so the missing
   third leg of the three-leg solution requirement is precisely
   "a clock-free closure core whose action preserves the 32-
   cluster at the R corners": recorded as the shape constraint on
   any future module-side candidate.
Contingency: the F_cert reads are taken from artifact fields
(V312/V209); if a field is absent the read is recorded VOID and
the T2 gate runs on the present reads only.

## Roll-up addendum: V379 shadow tribunal (2026-11-04, phase end)

The module side closed (V377/V378); the number side closed for its
best candidate (V379): 3 exp(-1/e) within 3.777149e-05 of the
observed hierarchy coefficient is a SHADOW - no web relation, no
exponent reading (nearest: module DIM 36, off by 0.0592 ~ 6
sigma). The tribunal programme now states the complete solution
geometry: a module-side solution, if it exists, is a single module
with (i) a clock-free closure core, (ii) the crystallographic walk
ladder, (iii) the published R multiplets - and its defect at the
best built points is measured: one dyadic 16+16 split of the R
32-cluster. No number-relation route to the ten open values
survives inside the bounded web. Both search directions are now
bounded by certified data, not by opinion.

## V380_cubic_axis_clockfree_core - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 293).
Property: the one carrier class the tribunal missed. V377/V378
tested clock-free cores about the ICOSAHEDRAL 5-fold axis (which
breaks the cubic R structure: signature [16, 48] at the gapless
points). A clock-free core about a CUBIC axis - a finite-order
element of the Spin(3) of the three SPATIAL planes about a
coordinate axis (order 4), a body diagonal (order 6 element,
2pi/3), or the same axes at order-2/3/4 angles - preserves the
cubic little-group structure at R by symmetry, closes at fixed
order independent of any clock, and runs on the crystallographic
transport (the ladder machinery of the published engine). The
three legs the tribunal proved necessary could co-occur HERE.
Class (booked before the run, 10 points): u_2i(omega, n) with
(Sx, Sy, Sz) the V337-certified spatial triple:
  - n in {e_x, e_y, e_z}, omega = pi     (2-fold axes: u = -2i n.S,
    u^2 = -I, order 4 in SU(2));
  - n in {e_x, e_y, e_z}, omega = pi/2   (4-fold angles);
  - n = (1,1,1)/sqrt3, omega = 2pi/3, and n = (1,-1,1)/sqrt3,
    omega = 2pi/3                        (body diagonals, 3-fold).
Transport beta = G10 (the published class - the ladder and the
frozen clock class are untouched); the internal element is the
ONLY deviation from the published engine.
Gates (all booked, all measured, host numpy):
 T0 self-test: the Spin(3) triple S_a^2 = I/4, anticommutators 0,
   re-certified (V337/V378 pattern); the published step anchor:
   M-gap < 1e-9, R signature == the V337 T0a reference, unitarity.
 T1 per point: unitarity of the step <= 1e-11 at GAMMA/M/R;
   internal closure order <= 12 (clock-free, expected 4/4/6);
   FULL-step closure at GAMMA <= 240; internal mirror <= 1e-10.
 T2 the conjunction: gapM (min |phase| over ALL M-points) <= 1e-9
   AND R signature == the published [32, 48] AND T1 all hold.
 T3 (booked reading, decided now): NO point of the 10 carries the
   T2 conjunction - the geometric obstruction V334 named survives
   even on cubic axes: a finite-order spatial-plane twist is NOT
   in the centralizer of the staggered transport's degeneracy
   structure; the booked words: CUBIC-CORE-STILL-BROKEN with the
   per-point table recorded. FALLBACK (booked, the solution
   branch): DILEMMA-BROKEN if any point passes T2 - then its
   (n, omega) and the measured closure orders are the carrier's
   first data, and the follow-up obligation is booked here: the
   clock-contact re-arbitration must be re-run (F2) with this
   carrier, because the closure no longer consumes the clock at
   ANY g - the clock becomes free in the step and pinned only by
   the ladder coefficient, dissolving the arbitration paradox
   that F2 fired on.
Contingency: if a diagonal-axis element cannot be built from the
triple exactly (irrational n components at float64), the closure
gate at T1 uses tolerance 1e-8 and the defect is recorded - the
conjunction gate is not loosened.

## V381_free_clock_carrier_arbitration - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 294).
Property: V380 fired the booked solution branch (DILEMMA-BROKEN:
u_2i(pi/2, coordinate axis), full-step closure order 2, gapM =
0.0 exact, R [32,48] exact, x3 axes). The booked obligation: the
F2 clock arbitration re-run on this carrier. The claim to test:
the staggered transport INVERTS the spatial-plane twist
(P R(theta) P = R(-theta)), so (PR)^2 = ±I for EVERY angle -
closure is clock-free at EVERY beta, and the clock becomes a
free parameter of a CLOSED step, pinned only by the ladder/
hierarchy contacts - the arbitration paradox F2 fired on
dissolves IF the physics contacts hold across the beta range.
Machine (booked): the carrier step S2.step(q, u_2i(pi/2, e_z),
beta) over beta in {pi/4, pi/3, pi/5, 2pi/5, 1/5, g_obs =
0.19999798835751165, 0.3, pi/6} (8 reads):
 T1 gates: unitarity <= 1e-11; FULL-step closure at GAMMA <= 24
   at EVERY beta (booked reading: order 2 everywhere); internal
   mirror <= 1e-10; 0-point (sum of GAMMA phases) <= 1e-9.
 T2 spectrum contacts at EVERY beta: gapM = 0 (min |phase| over
   all M-points <= 1e-9); R signature == [32, 48]; GAMMA
   degeneracy recorded.
 T3 the clock reading: at beta = g_obs and 1/5 the step is CLOSED
   (order <= 24) with T2 intact - the booked reading: TRUE - so
   the clock-contact distance to a closed, physical step is 0:
   the {pi/2, pi/4, 3pi/4} closure FORCING is void on this
   carrier. Record the full table.
 T4 ladder obligation (declared scope split): the hierarchy
   coefficient re-read (3(16 pi F_cross)^10) requires the
   certified traced-loop ladder machinery re-run on the carrier
   engine - a separate booking (V382); this experiment books
   only: the carrier's walk ladder is the PUBLISHED cubic walk
   (transport untouched - structural, recorded), and the F2
   verdict word: FREE-CLOCK-CARRIER-CONFIRMED if T1+T2+T3 all
   hold at all 8 betas; CARRIER-ARBITRATION-BREAKS with the
   failing beta/gate named otherwise.
Contingency: closure search caps at 24 for speed; if a beta reads
None at cap 24 but the phase table is a clean 2-torsion set,
record VOID + the table (no gate loosened).

## V381 booking correction (pre-run, declared; 2026-11-04)

The probe before the run measured the involution claim directly:
step^2 = +-I at beta = pi/4 (defect 1.33e-15) but NOT at beta =
1/5 rad (0.66). Corrected booked reading for T1/T3: the carrier
closes at order <= 24 at the SIX rational-multiples-of-pi betas
of the class (and the closure forcing leaves the booked set:
pi/3, pi/5, 2pi/5, pi/6 are closed, gapless-M and R-intact
clocks OUTSIDE {pi/2, pi/4, 3pi/4}) and does NOT close at the
two rational-RADIANS betas (1/5, 0.3) - the V378 transport-
incommensurability lesson verbatim. The clock is therefore
FREED into the whole rational-pi family, not made continuous:
the F2 question transfers to "does a rational-pi clock of the
carrier read the hierarchy coefficient" - which is V382's
booked question, not this one.

## V381 self-test arithmetic corrections (post-first-run, declared
2026-11-04; the PHYSICS gates are untouched)

The first run failed three gates; two failures were arithmetic
bugs in my own check lines, not physics: (1) the internal-mirror
line computed |ph|+ph (misplaced parenthesis), reading 1.6 where
the internal spectrum +-pi/4 x 16+16 is exactly mirror-paired;
(2) the zero-point line summed raw phases: at beta = pi/4 the sum
is exactly 2*pi - the 0-point is the sum MOD 2*pi, and 2*pi is 0.
Both lines corrected (no tolerance touched, no physics gate
changed); the run is re-issued. The third failure - T2 spectrum
contacts at every closed rational-pi beta - is real: only
beta = pi/4 carries gapM = 0 AND the published [32, 48]; every
other closed beta (pi/3, pi/5, 2pi/5, pi/6) and every rational-
radians beta reads R = [16, 48] with a gapped M. The booked
verdict CARRIER-ARBITRATION-BREAKS stands on that gate alone.

## V382_carrier_ladder_clock_tribunal - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 295).
Property: the last leg of the carrier story. V381 left the clock
arbitration alive: the carrier's physical contact rides only on
beta = pi/4 (the published clock); every other closed clock reads
R = [16, 48] with gapped M; the 1/5-rad pin never closes. The
open question, booked here as the tribunal's ONLY question: does
the certified hierarchy functional - the loop ladder of V206,
with its clock twist REPLACED by the carrier internal u_2i(pi/2,
e_z) (the transport hop phases e^{i g E /2} untouched, the ladder
walk untouched) - read the observed coefficient 2.076564111173832
(within the stated 1e-4 constants uncertainty) at any clock that
the CARRIER step actually closes at (closure <= 24 measured on
the V380 carrier at the carrier GAMMA point)? If yes, F2's
closure axiom dissolves on the carrier and the clock becomes a
certified reading; if no (booked reading), the clock arbitration
of V233 stands COMPLETELY unchanged by the carrier discovery.
Machine (imported V206 module, path-imported; no V206 edit):
 Fn_of replicated with G -> u_2i(pi/2, (1,0,0)); sector_hops,
 reading, CARTAN, LOOP, EV, P, D_L all imported unchanged;
 nmax = 600; coefficient = 3 (F/LOOP)^10.
 T0 replica self-test: the published Fn_of at {pi/2, pi/4, 3pi/4}
   reproduces the V206 closure_class coefs to 1e-9, and
   coef_at_fifth to 1e-9 - the imported machinery is the
   certified one.
 T1 carrier-ladder reads: at the five carrier-closed clocks
   {pi/4, pi/3, pi/5, 2pi/5, pi/6} (V381 orders 2/6/10/10/6):
   record n*, deficit, coef; gate: none within 1e-4 of observed.
 T2 dense scan: g on a 320-point grid over [0.05, 1.60]: for any
   g with |coef(g) - observed| <= 1e-4 (a rescue candidate),
   measure carrier-step closure at GAMMA (cap 24): gate: NO
   rescue candidate has closure <= 24.
 Decision: CARRIER-LADDER-NO-CLOCK-RESCUE (F2 stands on the
 carrier; the carrier is a structural discovery with NO clock
 consequence) / CARRIER-LADDER-CLOCK-HIT (a closed carrier clock
 reads the coefficient - F2 re-opened as a config decision).
Contingency: scan points with no crossing (None) record VOID and
are excluded from the rescue set (V206 convention).

## V382 T0 self-test threshold correction (post-first-run, declared
2026-11-04; the PHYSICS gates untouched and already measured)

The first run's T1 and T2 PASSED (carrier-ladder reads at the
five closed clocks: coefs 2.3708397 / 0.0049570 / 0.2733869 /
0.6437110 / 2.7153494, all missing the observed 2.0765641 by
0.27-2.07 - at least 2700x the band; dense scan 320 points:
ZERO band-hits anywhere on [0.05, 1.60], closed or not). T0
failed at max detune 3.559e-03 - the booked 1e-9 threshold was
mis-set against the machinery's OWN certified lesson: float64
closure-clock ladder readings DETUNE from the certified values
(V209's entire raison d'etre: certified pi/4 at dps >= 50
because float64 detunes at n* = 208 first-crossings). The
replica evidence shows the machinery is identical: the 1/5-rad
reading reproduces to 1.4e-13 (bit-level), while the three
closure-clock readings detune 5.8e-05 / 1.6e-03 / 3.6e-03 -
exactly the V206 float64 drift pattern. Corrected booked T0:
the 1/5-rad replica reproduces to < 1e-9 AND the closure-clock
replicas to < 5e-3 (float64 ladder drift, recorded). The status
branch was also corrected: CLOCK-HIT only when a CLOSED rescue
candidate exists (the first run's empty-closed_rescue branch
printed a misleading word); a failed T0 now reads
TRIBUNAL-REPLICA-INVALID. Physics gates and all tolerances on
T1/T2 are untouched.

## Roll-up addendum: the carrier arc V380-V382 (2026-11-04, suite
295)

The tribunal closed the module side; the number side closed for
its best candidate (V379); the carrier arc closed the YEAR by
finding the carrier and losing it in the same breath. V380:
u_2i(pi/2, coordinate axis) - a quarter-turn of the spatial-plane
Spin(3) about a cubic axis on the published transport - carries
the FULL three-leg conjunction at the published clock (closure
order 2, gapM = 0.0, R [32, 48] exact, x3 axes): the crystalline
dilemma is broken as an existence statement. V381: the closure
family frees to every rational multiple of pi (2, 6, 10, 10, 6)
but the physics rides ONLY at pi/4 - every other closed clock:
R = [16, 48] dyadic split + gapped M; the rational-radians pins
(1/5 rad, g_obs) never close: F2's arbitration re-runs on the
carrier and SURVIVES (booked fallback fired, the suite's second
booked FAILED alongside V369 - the gate fails on real physics
by pre-registration, not on machinery). V382: the certified
ladder with the carrier internal reads the coefficient at NO
clock (five closed clocks off by 0.27-2.07; 320-point dense scan:
zero band-hits) - the pi/4 read moves 2.7002 -> 2.3708, the
first certified machinery to move the loop coefficient
meaningfully, still 2900x the band. The end state of the
autonomy programme: the solution's anatomy is fully bounded from
every side - module side (V377/V378), number side (V379),
carrier side (V380-V382) - and every bound is certified data,
not opinion. Two self-test arithmetic corrections were declared
mid-arc (V381 mirror/zero-point lines; V382 replica threshold
against the V209 float64-drift lesson); no physics gate or
tolerance was moved after any measurement.

## V383_cubond_clifford_engine - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 296).
Property: the named hole of V118/V185 - does ANY causal unitary
engine on the certified algebra carry a NON-commutative tangent
algebra (the cone precondition every constructed step failed,
measured 0.0e+00 for 285 experiments)? Construction found by
symmetry (paper-first, then the single user-authorized probe):
the Clifford-bond engine. gamma_j = P_j (x) R_j: P_j the Pauli
of bond-link qubit j (X_1, Y_2, Z_3 - commuting on the 8-dim
bond factor of V138), R_j = 2S_j the certified V337 spatial-
plane spin triple (involution, anticommuting, exactly 0.0e+00).
P commuting x R anticommuting => gamma anticommuting: EXACT
Clifford on 256. Engine U(q) = prod_j (cos q_j + i sin q_j
gamma_j).
Gates (booked before the formal run):
 T0 self-tests: triple (2S)^2 = I and anticommutators <= 1e-12;
   Pauli algebra exact.
 T1 gamma: gamma_j^2 = I, {gamma_i,gamma_j} = 0, Hermiticity,
   all <= 1e-12.
 T2 unitarity of U(q) at 10 random q <= 1e-12.
 T3 Fourier support: exactly the 8 corner coefficients
   (n in {-1,+1}^3) nonzero at 1e-10 - hops within one step of
   the nearest-neighbor cube; the front radius sqrt(3) (against
   the published radius 1) is RECORDED, not a gate.
 T4 leading cone: spectra of every gamma_j = +-1, 128+128
   (principal velocities (1,1,1) exactly - the published node
   measured (0.00000, 0.14189, 0.14189), V118); and 8 random
   directions at a=1e-3: every phase within 3a^2 + 1e-10 of
   +-*|q| (the O(q^2) band-bend bound booked now - a lattice
   cone is exactly linear only to leading order, as V102
   taught).
 T5 velocity algebra: ||[gamma_i,gamma_j]|| = 2.0 +- 1e-12 at
   all three pairs; booked reading: >= 1 everywhere - V118's
   cone precondition REALIZED for the first time (published
   class: exactly 0.0e+00, the separation lemma V185).
 T6 carrier embedding: (I + 2iS_z)/sqrt2 = u_2i(-pi/2, e_z)
   to 1e-12 - the V380 carrier is an exact half-angle factor of
   the family (the probe measured 1.1e-16 against -pi/2 and a
   mis-signed 1.4e+00 against +pi/2: the sign is booked).
 T7 node census on the 5^3 grid {0,pi/2,pi}^3: gapless sites
   (gap to {0,pi} <= 1e-9) exactly 9 - the {0,pi}^3 octet
   (NN chiral doublers: 8 counting the origin) plus the
   (pi/2,pi/2,pi/2) corner; the 18 face/edge pi/2-type sites
   gapped with floor pi/2 (recorded, not gated). DOUBLERS
   PROLIFERATE - 16 gapless sites over the full BZ by cubic
   symmetry: Nielsen-Ninomiya kept and then some; recorded in
   the verdict word.
Verdict: CUBOND-EXIST-DOUBLERS-PROLIFER (the O9 algebraic
precondition is realized on the certified algebra; the ladder
functional on an entangled engine is the NAMED open design
problem - F2/O8 stand UNtouched: this engine is a new family,
never a replacement of the published one, and it was never run
against the published anchors because it is not the published
engine). Negative branch: any gate fails -> CUBOND-FAILS with
the gate named.
Contingency: none of the gates uses the GPU; host numpy,
256-dim, seconds.

## Roll-up addendum: the cubond theorem V383 (2026-11-04, suite
296)

The creative turn the programme was waiting for: after 285
experiments measured commutative tangent algebras (0.0e+00) and
V185 proved the Kronecker form responsible, the V118 cone
precondition was realized NOT by repairing the published bond
but by moving the algebra into the bond: gamma = Pauli (x) spin
is an exact Clifford on the certified algebra; the engine is
unitary, corner-causal, velocity-commutator 2.0, cones exactly
isotropic at speed 1, and it contains the V380 carrier as an
exact factor. The theorem is booked with its full price: 16
gapless doubler sites, front sqrt(3), and the ladder functional
on an entangled engine named as the open design problem - the
one certified road left to O8/F2. Suite 296: 293 publication +
3 alt_step; booked FAILs stand at two (V369, V381 - both the
booked physics branch firing on real data), plus V373's crash
log (infrastructure: GPU contention from an unrelated user
process holding 22.3 GB) pending its re-run.

## V384_cubond_ladder_functional - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 297).
Property: the named open design problem of V383 - the certified
hierarchy functional (V206: sector hops, Cartan norms, rung
crossing, LOOP) was built for internal-only 32-dim twist steps.
Does the ENTANGLING bond engine carry the observed coefficient
reading? The design, booked before any measurement:

THE LIFT (structural, forced by linearity, not chosen freely):
V206's sector hop is linear in the engine operator G:
ma_lambda = Tr_L[(P_lambda (G x I))(I x rho)] equals
ma_lambda = A_lambda C for 32-dim blocks C, with the FIXED
content operator A_lambda = (1/45) sum_j P_lambda^{jj} (the
diagonal-block sum of the sector projector; rho is the maximally
mixed adjoint state - both are V85/V79 content data, unchanged).
For any engine on bond(8) x internal(32), the hop lifts BLOCK-
WISE: M_engine = sum_lambda phase_lambda (I_8 x A_lambda) G_eng
(exactly the published dressing with the bond factor inert);
the Cartan norm lifts as norm(X) = sum_i |tr(C_i Tr_8 X)|/256
(Tr_8 = bond block-trace; 256 = 8x32 keeps the published
normalization); F_n, first crossing, deficit, 3 d^10: UNCHANGED.
The engine clocked: G_gamma(g) = prod_j (cos g + i sin g
gamma_j) = the cubond engine at the cubic-symmetric momentum
(g,g,g) - expansion over the 8 fixed corner monomials
Gamma_S (S subset {1,2,3}) with coefficients f_S(g);
the clocked left-dresser L(g) = sum_lambda e^{i(g/2)E_lambda}
(I_8 x A_lambda).

GATES:
 T0 (the lift anchor - the design must die if this fails): run
   the lifted machinery on the PUBLISHED engine lifted to 256,
   G = I_8 x step_pub(g): it must reproduce the certified V206
   closure_class readings within the V382 convention (5e-3,
   float64 ladder drift) AND coef_at_fifth within 1e-8
   (bit-level class: crossing n* = 208 stable).
 T1 cubond reads at the five clocks {pi/4, pi/3, pi/5, 2pi/5,
   pi/6} (closures of the cubond factors at rational pi):
   n*, deficit, coef; gate: none within 1e-4 of the observed.
 T2 dense scan: 96 points on [0.05, 1.60] (cost: 256-dim,
   booked reduced grid vs V382's 320): no |coef - observed| <=
   1e-4 with G_gamma(g*) closed to cap 24.
 Decision: CUBOND-LADDER-NO-CONTACT (the entangling bond does
 not carry the coefficient reading either - the ladder content
 is not the bond geometry, and F2's named residue stands with
 one fewer candidate mechanism) / CUBOND-LADDER-READS-OBSERVED
 (g* found: the design question becomes physical - F2 re-opens
 as a NEW program: the physics-contact census of the cubond
 family, booked here as the obligation: gapless content,
 chirality, and the clock closure of G_gamma(g*)).
Contingency: VOID crossings (no n* <= 600) recorded and
excluded (V206 convention); A_lambda built from the V206
module's P/EV import - no re-diagonalization.

## V384 norm-lift correction (declared 2026-11-04, during the T0
design-anchor phase; T1/T2 physics gates NEVER evaluated before
this correction)

The booked norm-lift (Cartan content of the BOND TRACE) is
wrong, and the booking's own death-gate caught it: the bond
trace of any monomial of the cubond engine is either the scalar
cos^3(ng) identity (Cartan-traceless) or zero - the Pauli bond
coefficients are traceless - so the Cartan norm of the DENOMINATOR
G^n vanishes identically for every entangling engine (the
division-by-zero crash IS the measurement of that fact). The
corrected lift: the Cartan norm over EVERY bond block,
norm(X) = (1/(8.32)) sum_{a,b} sum_i |tr(C_i X^{ab})|: on the
lifted published engine (blocks = 8 copies of the internal step,
off-diagonal zero) this equals the published cartan_norm term
by term EXACTLY (8/256 = 1/32), and on entangling engines it is
generically nonzero (off-diagonal Pauli blocks and the diagonal
SxSySz block carry the Cartan content). F_n is a ratio, so the
normalization constant is immaterial; the anchor gate T0 keeps
its 1e-8/5e-3 thresholds. T1/T2 untouched and unrun.

## Roll-up addendum: the cubond pair (V383-V384, 2026-11-04,
suite 297) - the year's creative arc completes

V383: the V118/V185 named hole - non-commutative tangent
algebras on the certified algebra - REALIZED (exact Clifford
bond engine, cones isotropic at speed 1, carrier an exact
factor). V384: the ladder of the same family, written through
the forced block-wise lift (anchor at 2.8e-13 on the published
readings - the design's own death-gate caught two lift bugs
before any physics gate ran), reads the coefficient at NO
clock. Pair booked in the paper (section 6.3): the cone and the
coefficient are carried by different circuit classes in
everything built. Design space closed family by family: the
published engine (F2 FIRED), module replacements (V377), the
carrier (V382), the entangling Clifford bond (V384); the only
reading that reaches the observed value anywhere remains the
unclosed rational-radians pin 1/5 (to 0.55%). The residue is
named, final in form, and every route to it is now a certified
measurement rather than an opinion.

## V385_cubond_content_census - BOOKING (2026-11-04, run following;
tier 2; category root_rebuild; suite 298).
Property: quelles particules habitent le cone du lien de Clifford.
Machine: les 256 etats de bond(8) x interne(32) du moteur V383;
les cinq generateurs de Cartan 32-dim C_k = 1j(h_a h_{a+1}),
a in (0,2,4,6,8) (clifford.generators, pattern V65.cartan_weights
mais en 32-dim); base propre commune par combinaison sqrt(1,2,3,5,7).
T0 self-test: algebre gamma re-certifiee (gamma_j^2 = I,
anticommutateurs <= 1e-12); la forme 16-dim de reference
(V65.cartan_weights) reproduit son propre multiset de 16 poids
(off <= 1e-9, comme certifie dans ses notes).
T1 poids 32: les 32 poids (5 composantes) du multiset interne
sont exactement la reunion lam16 et -lam16 (multi-tri equals,
chaque composante <= 1e-9) — le miroir de V52 verifie sur le
moteur encheveetre; lire lam16 depuis l'artefact V65 (cle
"weights" ou equivalente; grep "cartan" dans
archive/experiments/**/V65_wilson/results/v65_wilson.json; si
absent, recalculer via V65.cartan_weights() dans le script,
auto-enregistre).
T2 contenu au point Gamma: U(0) = I: la totalite des 256 etats
est a phase 0 (gap 0.0e+00 attendu <= 1e-12); le multiset des
256 poids (8 copies x 32) egale 8 x (lam16 + (-lam16)).
Decision: CUBOND-CONTENT-VERBATIM (le cone porte exactement le
contifie certifie de la cellule — ni nouvelle particule ni
perte: le theoreme d'existence ne modifie pas l'inventaire) /
CUBOND-CONTENT-SHIFTS (un des gates de multiset echoue — nommer
le poids manquant).
Contingence: aucune; numpy host, 256-dim, secondes.

## V386_cubond_chirality_nn_ledger - BOOKING (2026-11-04, run
following; tier 2; category root_rebuild; suite 299).
Property: le cone est-il chiral comme le modele standard?
Machine: G5 = gamma_1 gamma_2 gamma_3 sur 256; derivees
tangentes exactes: F_j' = i gamma_j F_j (produit de facteurs
commutants), S_j(node) = i dU/dq_j evalue analytiquement aux
16 sites gapless de V383 T7 (les 8 de {0,pi}^3 et les 8 coins
(+-(pi/2)) par symetrie cubique — enumerer explicitement les
16 sites: {0,pi}^3 et (±pi/2,±pi/2,±pi/2)).
T0: reprise exacte des gamma de V383 (0.0e+00).
T1: G5 anti-hermitienne (G5 + G5^dag = 0 <= 1e-12), G5^2 = -I
(<= 1e-12), spectre = +-i avec multiplicités 128/128.
T2: en chaque site, S_j = somme_k M_jk gamma_k + terme scalaire
(projection exacte: M_jk = 1/256 tr(S_j gamma_k) — les gamma
forment une base orthogonale des produits de Clifford; si la
composante scalaire depasse 1e-9, l'enregistrer, ne pas la
corriger); chiralite chi(site) = sign det M; gate: somme de chi
sur les 8 sites {0,pi}^3 = 0 (bilan NN: les doubleurs annulent
la chiralite) et la table des 16 chi est enregistree dans data.
T3 (contenu chiral): aux 8 sites {0,pi}^3, les operateurs S
ferment l'algebre de Clifford (normes des commutateurs >= 1);
gate: 8/8 sites.
Decision: CUBOND-CHIRAL-LEDGER-NN-CLOSED (le cone est vectoriel:
128+128 helictes, bilan total nul — le contenu n'est PAS chiral
comme le SM dans cette famille, enregistre comme mesure) /
CUBOND-CHIRAL-ANOMALY (un gate echoue — nommer le site).

## V387_cubond_doubler_masses - BOOKING (2026-11-04, run following;
tier 2; category root_rebuild; suite 300).
Property: que deviennent les masses des 16 parasites quand le
facteur de Wilson certifie (V65) agit sur le lien de Clifford.
Machine: U'(q) = (I_8 x W(q)) @ U_gamma(q), W(q) = diag(exp(-1j *
w * phi(q) * d32)), phi(q) = somme_j (1 - cos q_j) (forme exacte
de V65.step_wilson l.122-127), d32 = diagonale des charges de
D(theta) = cos(theta) X + sin(theta) Y en 32-dim: X,Y =
combinaisons Cartan de coefficients ux[0], ys[0] de
V65.identify_XY(lam16) (V65_wilson l.91-116, meme base Cartan
32-dim que V385); theta = 1.7308 (provenance: texte V86; grep
"theta" dans l'artefact V86 d'abord, utiliser la constante
bookee si absent); grille w = {0.1,0.2,...,1.0}.
T0: le multiset des valeurs propres de C_X en 16-dim egale
x_target = tri([1x10, -3x5, 5]) (<= 1e-9) — la certification V65
rejouee; W(w=0) = I: U'(q) == U_gamma(q) bit a bit (<= 1e-14).
T1 (survivants): sur la grille 5^3 {0,pi/2,pi}^3 x 8 sites
pleins: compter les sites gapless (gap <= 1e-9) pour chaque w;
prediction bookee AVANT mesure: le nombre de poids 32 de charge
nulle c = ux[0]/ys[0] (c·lam = 0 a 1e-9 pres) x 8 facteurs de
liaison = N_pred — calculer N_pred dans le script a partir des
poids V385 (artefact lu) et l'imprimer; gate: pour tout w >= 0.5,
gapless == N_pred exactement.
T2 (cone vivant): gap a Gamma = 0 <= 1e-12 pour tout w; pentes:
S'_j mesurees par difference finie a q=1e-4, h=1e-5:
||S'_j - gamma_j|| <= 1e-5 (le facteur W est O(q^2) a l'origine).
T3 (masses enregistrees): pour w = 1.0, aux 16 sites parasites,
enregistrer le gap minimal des etats gape (data, pas de gate).
Decision: CUBOND-DOUBLER-MASS-SELECTIVE (les parasites acquierent
une masse, le cone survit intact, les survivants = contenu de
charge nulle — le probleme des 16 est resolu au niveau du
mecanisme certifie) / WILSON-DOUBLERS-EVADE (T1 echoue) /
WILSON-KILLS-CONE (T2 echoue).

## V388_cubond_sm_selection - BOOKING (2026-11-04, run following;
tier 2; category root_rebuild; suite 301).
Property: la selection de contenu du mecanisme Wilson sur le lien
de Clifford reproduit-elle le pattern du modele standard?
Machine: le dispositif de V387 a w = 1.0; sous-espace gapless
aux sites doubles = etats de charge nulle; recenser leurs poids
d'hypercharge Y = (coefficients ys[0]) sur ce sous-espace;
multiset mesure vs target certifie V112 (artefact
v112_dark_u1_portal.json, lire la cle des multiplicites de
classes; si la cle est absente, target bookee = [(1,1.0),
(1,-0.0),(6,0.166667),(3,0.333333),(3,-0.666667),(2,-0.5)]).
T0: re-utilisation des artefacts V385/V386/V387 (aucun recalcul
des poids/charges — les lire).
T1: cardinalite du survivant == somme des multiplicites target.
T2: multiset des valeurs Y des survivants == multiset target
aux valeurs arrondies a 6 decimales pres (tri equals <= 1e-6).
Decision: CUBOND-SELECTS-SM-PATTERN (T1+T2 vrais — le contenu
sans masse du lien de Clifford habille porte exactement les
classes du SM: resultat physique positif majeur) /
CUBOND-SELECTS-OTHER (T1 vrai, T2 faux: enregister les valeurs
mesurees) / CUBOND-SELECTS-LESS (T1 faux: enregister cardinalite
et valeurs).
Contingence: aucune; numpy host; secondes-minutes.
