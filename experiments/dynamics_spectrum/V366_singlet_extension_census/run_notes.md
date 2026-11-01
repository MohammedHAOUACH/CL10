# V366 run notes - singlet-extension census (reality-interface first cut)

Run: GPU 1 pinned, 18.7 s, exit 0, tier 2, residual True.
VERDICT: **SINGLET-ADMIT-E1:1-E2:1-E3:0-E4:2-E5:9** - all booked
group-theory anchors fired exactly, E0 anchor reproduced V259 verbatim
(published-fiber ker C2 = 0, min nonzero eigenvalue 45.000000).

## What the machine established

The certified Spin(10) gauge action (generators VERBATIM from
V79/V251/V259: T16 = AC.weyl_generators, root_generators on 32,
adjoint_link_operators) is NOT the obstruction. The obstruction is
CONTENT: on the published V251 fiber (8 x 16 + 32 x 45) there is no
gauge singlet anywhere (V259 KERNEL-EMPTY, re-certified here as E0),
but with enlarged content singlets exist already at small fibers:

- E1 = 16 (x) 16-bar: ker 1 (the meson delta, 16 (x) 16-bar =
  1 + 45 + 210, second C2 level 64.0000 = adjoint Casimir).
- E2 = 16 (x) 16-bar (x) link-adjoint: ker exactly 1 - ONE
  gauge-invariant matter+flux state.
- E3 = 16 (x) 16 (x) link (no conjugate): ker 0 booked - the conjugate
  pair is necessary, not optional.
- E4 = 32 (x) 32 (x) 45 (published Dirac block DOUBLED): ker exactly 2.
  Two matter quanta + one flux unit make gauge singlets; one of each
  cannot.
- E5 = published 1568 (x) 16-bar: ker 9 = 8 (flavor singlets from the
  8 x 16 block against the conjugate) + 1 from (32 (x) 45) (x) 16-bar.

So V98/V259/V268's "the physical sector the Gauss law allows is empty"
is a theorem about the PUBLISHED one-body/one-flux fiber; it does not
extend to any content. Gauge-invariant composites start at two matter
quanta plus one flux unit, with explicit multiplicities (2 and 9).

## Machine hygiene (booked gates)

Lanczos with which='SM' on a degenerate kernel returns a non-orthogonal
basis (first attempt: E4 Gram defect 3.8e-02, E5 2.0e-01 while
||C2 W|| was already 5e-15): the kernel basis is re-orthonormalized by
Gram inverse square root (Gram min eig > 1e-8 booked assert); final
gates: ||C2 W|| <= 1.4e-14, W^dag W = I to 1e-15. Kernel DIMENSIONS are
Lanczos-robust (count |w| < 1e-9 with k = pred + 8 booked excess); the
dimension counts are corroborated independently by the exact match to
character theory on all five fibers.

## What it does NOT say (honest bounds)

- The fibers are single-block algebra (one link block, no lattice
  placement, no momentum, no dynamics): a singlet existing in the
  tensor algebra is not yet a physical asymptotic state of the theory.
- Enlarging content re-opens the anomaly ledger: the published
  cancellation is booked for the published content; E2/E4/E5 content
  must be re-checked against the 16/16-bar mirror anomaly bookkeeping
  before ANY physical claim (booked into V367).
- The V268 CONJUNCTION-DEAD legs (F2/F8/F9) are untouched by this
  experiment: V268 fired on theory-point conjunctions, not on kernels.

## Next (booked pre-V366 run: "V367 cascade live")

V367: place the E4/E5 singlet sectors on the torus (momentum, omega
rows as in V259), disperse them with the certified hop, and run the
V267-style ledger on those sectors - plus the anomaly-bookkeeping
re-certification for any content that survives dispersion.
