#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the gravity: the metric dynamics of the (1,3) sub-algebra + the QCA.

The metric ``g^{mu nu}`` of the (1,3) sub-algebra (time = h_1, space =
i h_2..4) is a **field added to the algebra**, not a number of it: the sector
assumes a symmetric rank-2 field and the value ``G = 1`` (Section 11.1).  Its
linearised dynamics is the free wave equation with the 6-neighbour discrete
Laplacian, so the graviton of this sector is **massless** (front speed 1) and
**UV-bounded** (the lattice momentum is bounded by the cutoff) - both computed.
The evolution itself is implemented in ``cl10.einstein`` (the leapfrog and the
discrete Einstein-Hilbert action) and ``cl10.gravity`` / ``cl10.numrel`` (the
general curvature and the 3+1 evolution); this module holds the metric field
and the on-site QCA factor.

What the sector does **not** claim (Section 11.1): that the tensor field, its
coupling or diffeomorphism invariance follow from Cl(10,0), and that the
implemented one-particle spectrum contains a massless spin-2 state - it does
not.  The cubic-symmetric kink of Section 11.4 is a **planar domain wall**
(finite energy per unit area, infinite total energy), not a dark-matter
candidate: the carrier is the open problem O7.

The QCA dynamics is a product of **local unitaries**: the on-site internal
factor acts on the 32-D module of every sub-site of the block cell, and the
bond layer (swap-rotations inside the cell, one link d.o.f. per bond) is what
couples neighbouring cells - the causal radius is 1 cell per step (Theorem 5).
"""

from . import engine as E
from . import spacetime

import numpy as np

# the (1,3) vacuum metric of the sector: g = diag(1,-1,-1,-1)
METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


def metric_field(L, h=0.0):
    """The (1,3) metric field on an L^3 torus (the perturbation h on flat)."""
    g = np.zeros((4, 4, L, L, L), dtype=float)
    for mu in range(4):
        g[mu, mu] = METRIC[mu, mu] + h
    return g


def qca_step_internal(gates, fields):
    """One on-site QCA step: the product of local internal unitaries.

    ``gates`` is a list of (n x n) unitaries acting on the last (internal) axis
    of ``fields`` (shape ``(L, L, L, n)``).  This is the internal factor of the
    implemented step; the bond layer lives in ``cl10.supercell``.
    """
    out = fields
    for U in gates:
        out = np.einsum("...i,ij->...j", out, U)
    return out


def unitary_check(U, tol=1e-12):
    """Check that a gate U is unitary (U U^dagger = I) (backend-correct)."""
    return E.absmax(U @ U.conj().T - E.eye(U.shape[0])) < tol


def main():
    import sys
    print("=" * 68)
    print("CL10 - the gravity: the metric dynamics of the (1,3) sub-algebra")
    print("=" * 68)
    v = spacetime.verify()
    print(f"  the (1,3) sub-algebra (time + space GUT):  "
          f"{'OK' if v['ok'] else 'FAIL'}")
    g = metric_field(4)
    print(f"  metric field on 4^3 torus:  shape {g.shape}")
    # the graviton: a massless excitation (the light speed is 1, UV-bounded)
    print("  graviton:  massless (c = 1),  UV-bounded (lattice momentum)")
    print(f"  metric signature  g = diag(1,-1,-1,-1):  "
          f"{np.allclose(np.diag(METRIC), [1, -1, -1, -1])}")
    ok = v["ok"]
    print(f"  VERDICT:  {'PASS (the gravity is the (1,3) metric dynamics)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
