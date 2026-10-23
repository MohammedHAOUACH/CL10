#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the root: the Clifford algebra Cl(10) of the 10 GUT directions.

The single object of the theory.  Ten mutually-anticommuting generators
``h_1 .. h_10`` (each ``h_a^2 = I``) generate a ``2**10 = 1024``-dimensional
complex algebra; the irreducible module is the **32-D** (complex) representation
with a **16-D Weyl** spinor (the matter) and its conjugate **16-bar**.

Construction:  the 10 generators are built from **5 qubits** with the
symplectic Pauli pattern — for ``k = 0..4``::

    g_{2k}   =  Y^k  (x) X  (x) I^{4-k}
    g_{2k+1} =  Y^k  (x) Z  (x) I^{4-k}

which gives exactly 10 Hermitian generators, each squaring to ``I`` and
**mutually anticommuting** (verified below).  This is the GUT ``Cl(10,0)``.
"""

from . import config as _config
from . import engine as E

import numpy as np
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml
# the structure of the root (config block `physics`; n_qubits is derived from
# n_internal by the loader): editing the yaml moves the module
NQ = int(_config.CFG.physics.n_qubits)      # 5 qubits: 2**5 = n_internal
NGEN = int(_config.CFG.physics.n_gen)       # the 10 Cl(10) generators
DIM = int(_config.CFG.physics.n_internal)   # the 32-D module

# the Pauli matrices (2x2, complex)
_X = E.Array_([[0, 1], [1, 0]])
_Y = E.Array_([[0, -1j], [1j, 0]])
_Z = E.Array_([[1, 0], [0, -1]])
_I = E.Array_([[1, 0], [0, 1]])


def _kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = E.kron(out, m)
    return out


def generators():
    """The 10 Cl(10) generators: a list of 10 (32x32) complex matrices.

    ``g_{2k} = Y^k (x) X (x) I^{4-k}``, ``g_{2k+1} = Y^k (x) Z (x) I^{4-k}``
    (5 qubits = 32-D).  Each squares to ``I``; all 10 anticommute.
    """
    out = []
    for k in range(NQ):
        prefix = [_Y] * k
        tail = [_I] * (NQ - 1 - k)
        out.append(_kron_n(prefix + [_X] + tail))
        out.append(_kron_n(prefix + [_Z] + tail))
    return out


def verify(tol=_TOL.strict):
    """0-parameter check of the Cl(10) relations (all 10 anticommute, square to I).

    Returns a dict of the residual errors (all must be ~0).
    """
    g = [E.asnumpy(m) for m in generators()]
    n = DIM
    Id = np.eye(n, dtype=complex)
    sq = max(float(np.abs(g[a] @ g[a] - Id).max()) for a in range(NGEN))
    ac = max(float(np.abs(g[a] @ g[b] + g[b] @ g[a]).max())
             for a in range(NGEN) for b in range(a + 1, NGEN))
    return {"n": len(g), "square_to_I": sq, "anticommute": ac,
            "ok": (sq < tol and ac < tol)}


def basis():
    """The 1024-D algebra basis: all 1024 ordered products of the 10 generators.

    (2**10 = 1024: the identity + 10 single + 45 double + 120 triple + ... +
    the volume element.)
    """
    import itertools
    g = [E.asnumpy(m) for m in generators()]
    I = np.eye(DIM, dtype=complex)
    out = [I]
    for r in range(1, NGEN + 1):
        for combo in itertools.combinations(range(NGEN), r):
            prod = I
            for a in combo:
                prod = prod @ g[a]
            out.append(prod)
    return out


def pseudoscalar():
    """The volume element (the product of all 10 generators) = the 16/16-bar split.

    ``w = h_1 h_2 ... h_10`` is odd and anticommutes with every single
    generator; it projects the 32-D module into the **16 Weyl** (+) and
    **16-bar** (-) halves.
    """
    g = [E.asnumpy(m) for m in generators()]
    w = np.eye(DIM, dtype=complex)
    for a in g:
        w = w @ a
    return w


def weyl_spinor():
    """The 16-D Weyl and 16-bar projectors from the pseudoscalar's eigenspaces.

    The pseudoscalar (the product of 10 generators) is **anti-Hermitian** (the
    10 generators anticommute), so it has **purely imaginary** eigenvalues
    ±i, each with multiplicity 16.  The (+i) eigenspace is the 16 Weyl, the
    (-i) the 16-bar.  We use ``eig`` (not ``eigh``, which assumes Hermitian).
    """
    w = pseudoscalar()
    vals, vecs = np.linalg.eig(w)
    plus = vecs[:, np.isclose(vals.imag, 1, atol=_TOL.loose)]
    minus = vecs[:, np.isclose(vals.imag, -1, atol=_TOL.loose)]
    return plus, minus


def main():
    import sys
    print("=" * 68)
    print("CL10 - the root: Cl(10) of the 10 GUT directions")
    print("=" * 68)
    v = verify()
    print(f"  10 generators (32-D):  square to I  {v['square_to_I']:.2e}")
    print(f"  all 10 anticommute            {v['anticommute']:.2e}")
    b = basis()
    print(f"  algebra basis:  {len(b)}  (2**10 = 1024)")
    w = pseudoscalar()
    print(f"  pseudoscalar^2  =  -I :  "
          f"{float(np.abs(w @ w + np.eye(32, dtype=complex)).max()):.2e}")
    p16, p16b = weyl_spinor()
    print(f"  Weyl 16 / 16-bar  dims:  {p16.shape[1]} / {p16b.shape[1]}")
    ok = (v["ok"] and len(b) == 1024 and p16.shape[1] == 16
          and p16b.shape[1] == 16)
    print(f"  VERDICT:  {'PASS (the Cl(10) is the root)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
