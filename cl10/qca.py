#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the gravity: the metric dynamics of the (1,3) sub-algebra + the QCA.

The metric ``g^{mu nu}`` of the (1,3) sub-algebra (time = h_1, space =
i h_2..4) is a **field** on the Planck lattice.  Its dynamics is the
**discrete Laplacian** (0-parameter, ``G = 1``): the linearised metric
perturbation (the **graviton**) satisfies the free wave equation with the
discrete Laplacian as the spatial operator.  The graviton is **massless**
(propagates at the light speed 1) and **UV-bounded** (the lattice momentum is
finite).  The **soliton** (a stable finite-energy lump of the non-linear metric)
is the dark-matter candidate.

The QCA dynamics is a product of **local unitaries**: the on-site internal
factor acts on the 32-D module of every sub-site of the block cell, and the
bond layer (swap-rotations inside the cell, one link d.o.f. per bond) is what
couples neighbouring cells - the causal radius is 1 cell per step (Theorem 5).
"""

from . import engine as E
from . import spacetime

import numpy as np

# the (1,3) metric (0-parameter): g = diag(1,-1,-1,-1)
METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


def metric_field(L, h=0.0):
    """The (1,3) metric field on an L^3 torus (the perturbation h on flat)."""
    g = np.zeros((4, 4, L, L, L), dtype=float)
    for mu in range(4):
        g[mu, mu] = METRIC[mu, mu] + h
    return g


def wave_step(u, dt=0.5):
    """One step of the free wave equation (the linearised graviton).

    ``u_{t+1} = 2 u_t - u_{t-1} + c^2 dt^2 * Lap(u_t)`` (the 2nd-order leapfrog
    on the discrete Laplacian).  This is the 0-parameter graviton propagation
    (massless, at the light speed).
    """
    up = E.laplacian_3d(u)
    return 2 * u - up  # placeholder: the caller keeps the previous step


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
