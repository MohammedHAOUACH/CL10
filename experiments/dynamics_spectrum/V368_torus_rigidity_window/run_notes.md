# V368 run notes - torus rigidity of the composite isolation window

Run: GPU 1 pinned, 10.4 s, exit 0, tier 2, residual True.
VERDICT: **WINDOW-NOT-RIGID-2.0e-05**. Anchor (contact gap 1.006652
vs V367 artifact, 1e-9) and omega=0 control (merge at 0.0) both fired.

## The booked prediction was REFUTED by the machine - cleanly

V367's contact isolation window does NOT extend rigidly over the
torus. At fixed total momentum K, the continuum is the union over
ALL relative momenta (216 of them on the 6^3 grid): the bond-pair
bands of continuum pairs at (k1', k2') != (k1, k2) land on top of
the shifted composite band. Min gap over the 216 total momenta:
2.0e-05 (max 4.1e-03, spread 4.0e-03) against 1.0066 at contact -
a 50,000-fold collapse. The 1.007 rad window was a CONTACT-LOCUS
feature of the phase comparison, not a property of the composite.

Why the rigidity argument broke: the factorization bloch =
kron(P(q), U^6) makes the internal-phase shift q-independent (that
part was measured true - V367 probe: kernel blocks = pure bond to
4e-15), but a bound state must dodge the continuum at the SAME
total momentum, and the fixed-K continuum sweeps all bond phases
with it. The correct bound-state condition is a threshold (energy
below the continuum edge), not a phase offset at fixed relative
kinematics - and the certified clock step is free propagation plus
the Gauss constraint, so no attractive channel exists to produce
one. The collapse is the generic mechanism, and the machine caught
it in 10 s.

## What the reality-interface programme now states (honest sum)

- V366: gauge-neutral matter states EXIST for enlarged content
  (kernels 1/1/0/2/9 exact; conjugate pairs necessary) - the
  published emptiness is a CONTENT theorem.
- V367: the singlet pair is an invariant subspace of the certified
  clock step and isolated AT CONTACT (1.007 rad), and the clock's
  internal layer acts as the identity on it (4e-15, probe recorded
  in the pre-reg addendum).
- V368: that isolation DISAPPEARS at generic torus momenta
  (2.0e-05): the composite disperses into the single-flux
  two-quanta continuum.
Conclusion: with the certified clock dynamics (free hop + Gauss),
the physical sector is empty in BOTH senses - the published content
has no gauge-neutral states at all (V259/V268), and the enlarged
content's neutral states are continuum-embedded (not asymptotic).
Door n°2 of CONJUNCTION-DEAD survives the first two attacks. A
stable composite would require DYNAMICS beyond the published clock
step (an interaction that produces a threshold below the continuum
edge) - that is the honest shape of the remaining open problem,
not a gap in this programme's machinery.

## Booked-gate correction on the record

The booked verdict word RIGID-WINDOW-BROKEN keyed "broken" to gap
<= 1e-6, which let a 50,000-fold window collapse pass under a
"window holds" name. Corrected before the artifact shipped (the
artifact on disk carries the corrected WINDOW-NOT-RIGID wording;
the first run's mislabeled status string was never committed).

## Limits

6^3 grid (V259 house grid); continuum subsampled at 24 booked
random relative points per K (seed 0x51) - finer grids can only
SHRINK the min gap (denser continuum), so 2.0e-05 is an upper
bound on the true minimum. One-flux sector only; multi-flux
continua would embed it further. None of this changes any
published number - all three experiments are residual tier 2.
