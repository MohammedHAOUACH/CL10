#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""The 2I engine (branch root-rebuild): the automaton's internal layer
carries a BINARY-ICOSAHEDRAL element instead of the 45-plane product.

The change booked by V337 (review/pre_registrations_root.md): the internal
unitary of the step becomes

    U_2I(omega) = cos(omega/2) I - 2 i sin(omega/2) (n1 Sx + n2 Sy + n3 Sz)

with Sx = T_34/2, Sy = -T_24/2, Sz = T_23/2 the Spin(3) triple built from
the THREE SPATIAL planes of the Clifford generators h_2, h_3, h_4
(T_ab = i h_a h_b: pairwise anticommuting, each S_a^2 = I/4 - certified in
the V337 T0b gate), and n the icosahedral 5-fold axis

    n = (0, 1, phi) / sqrt(1 + phi^2),   phi = (1 + sqrt5)/2,

the dodecahedron-vertex direction in the cubic frame. Because (n.S)^2 = I/4
the Rodrigues form is exact: U_2I is an order-10 element (U_2I^5 = -I), and
its adjoint on the spatial generator span is the SO(3) rotation of angle
omega about n - trace exactly 1 + 2 cos(2 pi/5) = phi. The digit 5 is
inside the engine's internal algebra.

The bond layer keeps the cubic stage algebra of the published automaton:
A5 has no transitive action on the 8 sub-sites (60/8 is not an integer), so
the cell cannot carry the icosahedral group - the change is internal-only,
and the step is built by the PUBLISHED machinery cl10.supercell.bloch,
which takes any internal unitary (the V239 bond_g10 dial decouples the hop
angle beta).
"""
import math

import numpy as np

from . import clifford
from . import engine as E
from . import supercell as SC

PHI = (1.0 + math.sqrt(5.0)) / 2.0
_NRM = math.sqrt(1.0 + PHI * PHI)
AXIS = (0.0, 1.0 / _NRM, PHI / _NRM)     # icosahedral 5-fold axis, cubic frame


def spin_triple():
    """(Sx, Sy, Sz) host complex128 32x32: the Spin(3) triple of the three
    SPATIAL planes: Sx = T_34/2, Sy = -T_24/2, Sz = T_23/2."""
    h = [E.asnumpy(m).astype(np.complex128) for m in clifford.generators()]
    def T(a, b):                          # 1-based generator labels
        return 1j * (h[a - 1] @ h[b - 1])
    return T(3, 4) / 2.0, -T(2, 4) / 2.0, T(2, 3) / 2.0


def u_2i(omega, axis=AXIS):
    """The order-family element: cos(w/2) I - 2 i sin(w/2) n.S (exact
    Rodrigues: (n.S)^2 = I/4)."""
    Sx, Sy, Sz = spin_triple()
    nS = axis[0] * Sx + axis[1] * Sy + axis[2] * Sz
    I = np.eye(32, dtype=np.complex128)
    return math.cos(omega / 2.0) * I - 2j * math.sin(omega / 2.0) * nS


def spatial_adjoint(U):
    """The 3x3 SO(3) matrix of the adjoint of U on (h_2, h_3, h_4):
    R_ab = Tr(h_a U h_b U^dag)/32 for the a,b in {2,3,4} frame."""
    h = [E.asnumpy(m).astype(np.complex128) for m in clifford.generators()]
    hs = [h[1], h[2], h[3]]
    Ud = U.conj().T
    R = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            R[a, b] = np.trace(hs[a] @ U @ hs[b] @ Ud).real / 32.0
    return R


def step(q, U, beta):
    """The one-step Bloch matrix with internal unitary U and hop angle beta,
    built by the published automaton machinery; returned on the host (the
    128/256-dim momentum matrices are measured faster there, house
    convention; GPU paths are engaged by V337's T0c gate and power scan)."""
    return np.asarray(E.asnumpy(SC.bloch(q, U, g10=beta,
                                         stages=SC.STAGES,
                                         bond_g10=beta)),
                      dtype=np.complex128)
