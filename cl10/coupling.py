#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 - the angle and the coupling: ONE source, TWO branches.

The step carries exactly one continuous number, the rotation angle ``g``.  The
publication keeps two readings of it apart (Sections 1.4, 7.1, 9, 15) and so
does the code:

``implemented``   ``g = 1/5``  (rotation-squared reading ``alpha_10 = g^2 = 1/25``)
    Provenance **convention**.  The rule ``g = 1/n_qubits`` (the five qubits of
    the 32-dimensional module) is not selected by the data: the ``alpha = 1/n^2``
    scan prefers ``n = 7`` (V56).  The phenomenology of Sections 12-14 is
    computed on this branch, and it is the default of the whole repository.

``closed``        ``g = pi/4`` (order 8)
    Provenance **algebraic**.  The closure axiom - the fundamental update must
    be an element of the root, ``e^{ig}`` a root of unity - admits exactly two
    angles over all rational ``p*pi/q`` (V131): ``pi/2`` and ``pi/4``.  At
    ``pi/4`` the 32 internal levels are exactly ``{+-pi/4 (10), +-3pi/4 (6)}``
    and the zone-centre step is an involution, ``M(0)^2 = I``.

``half``          ``g = pi/2`` (order 4, ``U = w``, the pseudoscalar)
    The second closed root.  The bond layer degenerates to a pure lattice
    translation and the gapless set is a line network (V134): no mixing, hence
    no mass generation.  Kept for reference, not used as a physical branch.

Select the branch with ``CL10_ANGLE=implemented|closed|half``; the default is
``implemented``, so every number already in the repository is unchanged.

Coupling conventions (Section 1.4, stated once and used everywhere):

    alpha_std(g) = g^2 / 4 pi    the standard normalisation, the one used
                                 wherever a coupling meets running
    alpha10(g)   = g^2           the rotation-squared reading, a factor-4pi
                                 convention, invoked only where named

At ``g = 1/5``: ``alpha_std = 0.00318``, ``alpha10 = 0.0400``.  At ``g = pi/4``:
``alpha_std = pi/64 = 0.0491``, ``alpha10 = pi^2/16 = 0.6169``.
"""

import math
import os

# the rotation-squared reading at the implemented angle (Section 1.4: a
# convention, not a number of the root)
ALPHA_10 = 1 / 25.0

G_IMPLEMENTED = 1 / 5.0
G_CLOSED = math.pi / 4
G_HALF = math.pi / 2

# the second coupling reading: g_0 = 1/2^n_space = 1/8 (Section 15: the reading
# of the 3 spatial directions — a convention, not a derivation)
G0 = 1 / 8.0

# name -> (angle, provenance, the experiment that fixes it)
BRANCHES = {
    "implemented": (G_IMPLEMENTED, "convention", "V56 (the scan prefers n = 7)"),
    "closed": (G_CLOSED, "algebraic (the closure axiom)", "V131 (order 8)"),
    "half": (G_HALF, "algebraic (the degenerate second root)", "V134 (pure translation)"),
}

DEFAULT_BRANCH = "implemented"


def branch():
    """The active branch, from ``CL10_ANGLE`` (default ``implemented``)."""
    name = os.environ.get("CL10_ANGLE", DEFAULT_BRANCH).strip().lower()
    if name in ("impl", "implemented", "1/5", "0.2"):
        name = "implemented"
    elif name in ("closed", "pi/4", "quantized"):
        name = "closed"
    elif name in ("half", "pi/2"):
        name = "half"
    if name not in BRANCHES:
        raise ValueError(f"CL10_ANGLE={name!r}: not one of {sorted(BRANCHES)}")
    return name


def angle(name=None):
    """The rotation angle of the given branch (default: the active one)."""
    return BRANCHES[name or branch()][0]


def provenance(name=None):
    """The provenance class of the angle on the given branch (Section 15)."""
    return BRANCHES[name or branch()][1]


def alpha_std(g=None):
    """The standard-normalised coupling ``alpha = g^2 / 4 pi``."""
    gg = angle() if g is None else g
    return gg ** 2 / (4.0 * math.pi)


def alpha10(g=None):
    """The rotation-squared reading ``alpha_10 = g^2`` (a named convention)."""
    gg = angle() if g is None else g
    return gg ** 2


# the angle of the active branch: the module-level default every builder uses
G10 = angle()
