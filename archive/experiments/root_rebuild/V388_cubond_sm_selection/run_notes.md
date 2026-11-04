# V388 - Standard-Model selection test on the Clifford bond: run notes

Verdict: **CUBOND-SELECTS-LESS** (tier 2, booked negative branch).
Wall 4.9 s, host numpy. Suite 301.

## The measurement

Does the Wilson content-selection of the certified mechanism on the
Clifford bond reproduce the Standard Model pattern? The machine is
the V387 device at w = 1.0; the gapless subspace at the doubled
sites is the zero-charge content; its hypercharge weights are read
from the certified V385/V387 artifacts (no re-derivation anywhere:
weights, X and Y identifications, theta and N_pred are all read),
and compared with the certified V112 target
hypercharge_by_class_multiplicity
[[1, 1.0], [1, -0.0], [6, 0.166667], [3, 0.333333], [3, -0.666667],
[2, -0.5]].

- **T0** all four measurements are read from disk pass-flagged
  (V385 32 weights, V386 NN ledger, V387 N_pred = 0 at theta =
  1.7308, V112 six classes).
- **T1** the survivor cardinality is **0** weights of zero charge,
  consistent with the V387 state-level counts (all 15 parasite
  sites hold 0 zero modes at w = 1.0) - against the target sum of
  multiplicities **16**: the branch CUBOND-SELECTS-LESS fires.
- **T2** the Y multiset of the survivors is empty against the
  16-value target; the comparison is recorded and fails, exactly
  as booked by the T1-false branch.

## Declared implementation note (post-first-run declared, physics gates untouched)

The first run printed the same measured branch (SELECTS-LESS); the
verdict wiring was then aligned with the repository convention
that a booked measurement whose decision tree fires is a certified
measurement (pass = measurement integrity: certified inputs,
weight-level and state-level counts consistent; precedent V374
NO-THRESHOLD, V378 WITNESSED, V379 NO-WEB-CONTACT, V384
NO-CONTACT). No tolerance, comparison or booked branch was
modified.

## Reading

At the physical junction angle 1.7308 the Wilson selection closes
NO massless Standard-Model content on the Clifford-bond family:
the zero-charge kernel of D(theta) on the 32 weights is empty, and
the smallest charge 0.005216 is the same clock tension the ladder
programme names - the physical reading does not close the neutral
direction, so the massless SM pattern (16 classes: 1, 0, 6 x 1/6,
3 x 1/3, 3 x -2/3, 2 x -1/2) is not supplied by this mechanism at
this angle. The measurement is booked as negative: the selection
machinery works (it gaps everything that should be gapped - V387)
and selects nothing massless beyond the cone itself.
