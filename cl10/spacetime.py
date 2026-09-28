#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the spacetime INSIDE the GUT.

The single 0-parameter idea of the theory: the (1,3) Minkowski spacetime is a
**sub-algebra** of the GUT ``Cl(10)``.  The **time** is the first GUT
generator ``h_1`` (``+``), and the **three space directions** are
``i*h_2, i*h_3, i*h_4`` (GUT generators times ``i``, so they square to
``-I``).  Because all four are GUT generators they **anticommute**, so the
Minkowski relations hold 0-parameter and the **light cone is a GUT
sub-algebra**.

The 6 remaining GUT generators (``h_5 .. h_10``) are the **internal**
directions (the gauge).
"""

from . import engine as E
from . import clifford

import numpy as np

# the 4 (1,3) generators inside the GUT (32x32)
def causal_generators():
    """The 4 (1,3) generators INSIDE the GUT: time = h_1 (+), space = i h_2..4 (-).

    Returns ``(T, [S1, S2, S3])`` where ``T^2 = +I`` and ``S_a^2 = -I`` and
    all four anticommute.
    """
    g = [E.asnumpy(m) for m in clifford.generators()]
    T = g[0]                      # h_1  (time, +)
    S = [1j * g[k] for k in (1, 2, 3)]  # i h_2..4  (space, -)
    return T, S


def verify(tol=1e-12):
    """0-parameter check that the (1,3) is a valid Minkowski sub-algebra of the GUT.

    Returns the residual errors (all must be ~0):  time^2=+I, space^2=-I,
    {T,S}=0, {S,S}=0.
    """
    T, S = causal_generators()
    I = np.eye(32, dtype=complex)
    t2 = float(np.abs(T @ T - I).max())
    s2 = max(float(np.abs(s @ s + I).max()) for s in S)
    ts = max(float(np.abs(T @ s + s @ T).max()) for s in S)
    ss = max(float(np.abs(S[a] @ S[b] + S[b] @ S[a]).max())
             for a in range(3) for b in range(a + 1, 3))
    return {"t2": t2, "s2": s2, "ts": ts, "ss": ss,
            "ok": (t2 < tol and s2 < tol and ts < tol and ss < tol)}


def light_cone(T, S):
    """The (1,3) light cone: the null vectors v = v0*T + v_a*S_a with v^2 = 0.

    ``v^2 = (v0)^2 - sum_a (v_a)^2`` (the Minkowski form).  The set of null
    directions is the causal structure of the lattice.
    """
    # a vector v = v0*T + sum v_a S_a ; v^2 = v0^2 - sum v_a^2 (the Minkowski form)
    # the light cone is the set of (v0, v_a) with v0^2 = sum v_a^2
    return lambda v0, v: (v0 ** 2 - sum(x ** 2 for x in v))


def main():
    import sys
    print("=" * 68)
    print("CL10 - the (1,3) spacetime INSIDE the GUT Cl(10)")
    print("=" * 68)
    v = verify()
    print(f"  time T^2  =  +I  :  {v['t2']:.2e}  (T = h_1, a GUT generator)")
    print(f"  space S^2 =  -I  :  {v['s2']:.2e}  (S = i h_2..4, GUT generators)")
    print(f"  {{T, S_a}} = 0     :  {v['ts']:.2e}  (the GUT generators anticommute)")
    print(f"  {{S_a, S_b}} = 0   :  {v['ss']:.2e}  (Minkowski (1,3) OK)")
    print(f"  -> the light cone is a sub-algebra of the GUT  "
          f"({'PASS' if v['ok'] else 'FAIL'})")
    sys.exit(0 if v["ok"] else 1)


if __name__ == "__main__":
    main()
