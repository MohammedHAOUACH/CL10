# V372 run notes - torus singlet tower census

Run: GPU 1 pinned, ~175 s, exit 0, tier 2, residual True.
VERDICT: **TORUS-SINGLET-TOWER-ALIVE** (ell=1: 2 channels x 12 bonds,
ell=2: 6 channels x 24 adjacent pairs -> 24 + 144 singlet states per
internal channel structure on the N=2^3 torus).

## Facts

1. Cluster decomposition (exact, dense eigvalsh of the pairwise C2):
   32x32 = {2 (trivial x2), 90 (45^2), 420 (210^2), 20 (10^2),
   240 (120+120bar), 252 (126+126bar)}; 45x45 = {1, 45, 54, 210,
   770, 945} (sym 1+54+945 = 1001, antisym 45+210+770 = 1025);
   16x16bar = {1, 45, 210}. All cluster SHAPES fire; the anchor
   multiplicities A1=2, A2=1, A3=2, E1=1, E2=1 all reproduce.
2. Shared C2 values between 32x32 and 45x45: exactly {0, 64 (the 45),
   96 (the 210)} with copy counts (2,1) each ->
   **ker(32x32 x 45x45) = 6** = 2 (trivial) + 2 (45) + 2 (210).
3. Casimir trace gate: tr C2(32) = 1440, tr C2(45) = 2880, ratio
   exactly 2 - the D5 weight model with unit u = 4 (C2(32) = 11.25u,
   C2(adjoint) = 16u, both matching the V366/V369 anchor 45.000000
   = 11.25 x 4 on the spinor and 64 = 16 x 4 on the adjoint).
4. Torus N=2^3: 8 sites, 12 bonds, 24 adjacent bond pairs, 3 bond
   orbits - the flux-line configuration counts V371 needs.
5. Pure-phase consequence (booked, from V370): the whole ell=1
   tower (24 states) carries Sum_b J_b.L_b = -64 - the published
   interaction cannot move a threshold there. The ell=2 channels
   (6 per pair, 144 states) are the first place the interaction is
   NOT a phase: interior Gauss makes the touching contraction
   (J_A - J_B).L_f1, and the 6 channels carry a non-scalar
   deformation matrix (measured by V371).

## Corrections booked before the run (see pre-registrations)

- Weyl weight-lookup formula dropped (fails its own anchor on
  paper) for the dense-Casimir cluster method.
- 45x45 decomposition: my weight sketch (45', 1670) WRONG; the
  machine says 770 + 945 (shape gate uses the machine fact).
- Torus adjacent pairs: 24 (undirected incidence), not 120
  (directed).
- T16 conjugate must use -T^T (V366 verbatim); -T^dagger computes
  16x16 (10+120+126), a different fiber.

## What it means for V371

The V371 Krylov/spectral question is licensed with a small exact
basis: 6 ell=2 intertwiners (Haar-averaged intertwiner construction,
residual-gated) + 2 ell=1 (V369/V370 verbatim) per geometry. The
deformed-step question becomes a 6x6 (plus 2x2) spectral problem
against the V367/V368 contact + torus continuum - no Krylov on
2.07M needed. Scope: contact + touching lines; torus hopping is
the V368-cited limitation and stays out.
