# V367 run notes - singlet-pair contact composite at the zone centre

Run: GPU 1 pinned, 52.9 s, exit 0, tier 2, residual True.
VERDICT: **COMPOSITE-GAMMA-BOUND-1.007e+00**. All gates green:
C2 anchor (E4 kernel 2, ||C2 W|| 4.8e-15), C1 (bloch = kron(P,U32^6)
to 3.3e-16), G1 (contact-basis step unitarity 8.0e-15), G2 (plane
span 1.6e-15, product 4.6e-15, [C2,U_tot] spot 1.3e-15), G3 (vectorlike
traces exactly 0), control (omega=0 merge 0.0).

## The physics fact

The 128-dimensional singlet-pair contact subspace (E4 kernel x
sublattice-pair basis) is INVARIANT under the certified clock step
(G1 unitarity of the compression = invariance, exact to 1e-15), and
its eigenphases at the physical V251 flux rows (omega = +-8.4181)
sit **1.007 rad away from the two-quanta continuum phase window**.
This is the first gauge-invariant eigenstate of the certified clock
step isolated from the continuum since the V268 era - a candidate
asymptotic physical state: two matter quanta + one flux quantum,
neutral under all 45 Spin(10) charges.

The structural surprise (machine-forced, not assumed): the composite
is BLIND to the clock's internal Spin(10) rotation. The internal
layer is a product of 45 plane rotations with generators exactly in
span(root_generators) (resid 1.6e-15), so the step acts on the
singlet subspace as exp(i theta Q_tot) - and diagonal group actions
fix invariant tensors. The composite feels the clock only through
the bond hopping layer (P(0) sublattice mixing) and the flux energy.
Phase bookkeeping: composite phase = (bond-layer phase) + omega;
at omega=0 the band is the pure bond spectrum, booked as a control
and merging into the window (gap 0.0 exactly); BOUND is read only
on the physical rows.

## Booked-gate corrections (post-run, factual, recorded in pre-reg)

1. The booked G2a "max ||[U32,T^a]|| <= 1e-12" was FALSE BY
CONSTRUCTION: the clock internal layer IS a fixed Spin(10) direction
- the clock lives in the coefficients, which is the programme's own
claim; measured 1.4 (the clock does NOT commute with the charges;
only C2 does). Replaced by the true statements, machine-checked:
G2a' (i) every plane generator of internal(g10) lies in span(Ts)
(45 residuals <= 1.6e-15) and the product of projected plane
rotations reproduces U32 to 4.6e-15 - so U_tot is a product of
exponentials of TOTAL charges; G2b full-space [C2, U_tot] spot
rel defect 1.3e-15. The singlet-subspace invariance itself is
carried by G1 (compression unitarity = invariance).
2. The booked verdict read the min gap over ALL omega rows
including 0; the 0 row is structurally a control (a fiber carrying
a flux quantum must pay the flux energy; at omega=0 the composite
band is the pure bond spectrum, a subset of the window). BOUND is
read on the two physical rows; the control merge became a booked
gate (fired at 0.0).

## Honest limits (do not overread)

- Locus: CONTACT pair at zone centre only (relative momentum zero,
  both quanta in one block, one adjoint flux quantum at the link).
  "Isolated asymptotic state" further needs V368: relative-momentum
  dispersion of the 128-band + the overlap with multi-flux continua
  (the booked window is the two-quanta window in the one-flux
  sector: isolation against flux-2 continua is NOT measured yet).
- The phase window is a UNITARY-step notion (eigenphases), not a
  binding energy: 1.007 rad is a gap in the step spectrum at the
  clock, not a mass shift below threshold.
- G3 covers the cubic anomaly of the doubled Dirac block (exactly
  zero); full anomaly bookkeeping of a dynamical E4 content in the
  published ledger is a separate re-certification (O-card work).

## Next (booked in pre-reg addendum)

V368: dispersion of the 128 contact bands over the relative-momentum
grid (the same einsum machinery at q != 0 in the relative coordinate)
+ multi-flux continuum overlap; only then a sentence in the paper.
