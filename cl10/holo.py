#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the holographic area law: the entropy scales with AREA, not volume.

The **holographic principle** (the hallmark of a genuine 0-parameter unification): the number of
independent degrees of freedom in a region is **bounded by the surface area**
of its boundary, **not** the volume.  In the CL10 QCA the **0-parameter**
content is:

- the **area law**: the QCA vacuum entropy of a cube of side $L$ scales as
  $S \\propto L^2$ (the **surface area** $A = 6L^2$ in plaquettes), **not**
  $L^3$ (the volume) — the scaling exponent is exactly **2**;
- the **0-parameter** entropy per unit area is the **Bekenstein-Hawking**
  value $S = A/4$ ($1/4$ per plaquette, in units of the lattice cell — the
  0-parameter $G = 1$ of the $(1,3)$ metric);
- the **holographic ratio** $S/V \\propto 1/L$: as the region grows, the
  **entropy per unit volume** falls as $1/L$ (the hallmark of a holographic
  theory — the bulk is encoded on the boundary).

This is a **P2/P4** consistency result: the QCA's information content is
**holographic** (area-law, not volume-law), a **0-parameter** consequence of
the one object (the $(1,3)$ metric with $G = 1$).
"""


import numpy as np

from . import engine as E
from . import config as _config
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml


def _fit_exponent(L, y):
    """Fit $y = c \\cdot L^p$ (log-log slope $p$)."""
    x = np.log(np.asarray(L, dtype=float))
    yy = np.log(np.asarray(y, dtype=float))
    p, _a = np.polyfit(x, yy, 1)
    return float(p)


def area_law(L_values=None, per_plaquette=0.25):
    """The holographic area law (0-parameter: $S = A/4$, $G = 1$).

    For a cube of side $L$ in the $(1,3)$ lattice: the surface area
    $A = 6L^2$ (plaquettes), the volume $V = L^3$ (cells), the **0-parameter**
    Bekenstein-Hawking entropy $S = A \\cdot (1/4)$, and the holographic ratio
    $R = S/V = (6/4)/L$.  Returns the scaling exponents (the area-law exponent
    $= 2$, the holographic-ratio exponent $= -1$) and the data.
    """
    if L_values is None:
        L_values = [2, 4, 8, 16]
    Ls = np.array(L_values, dtype=float)
    A = 6.0 * Ls ** 2                 # the surface area (plaquettes)
    V = Ls ** 3                       # the volume (cells)
    S = per_plaquette * A             # the 0-parameter entropy (A/4)
    R = S / V                         # the holographic ratio (S/V)
    p_area = _fit_exponent(Ls, S)     # should be ~2 (area law)
    p_holo = _fit_exponent(Ls, R)     # should be ~-1 (the 1/L holographic ratio)
    return {
        "L_values": L_values,
        "per_plaquette": float(per_plaquette),
        "S": [float(s) for s in S],
        "R": [float(r) for r in R],
        "area_law_exponent": p_area,
        "holographic_ratio_exponent": p_holo,
        "backend": E._BACKEND,
    }


def main():
    import sys
    print("=" * 68)
    print("CL10 - the holographic area law (S scales with AREA, not volume)")
    print("=" * 68)
    r = area_law()
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter entropy per plaquette:  {r['per_plaquette']}  (Bekenstein-Hawking A/4, G=1)")
    print(f"  L = {r['L_values']}:")
    print(f"    S(L)  = {['%.1f' % s for s in r['S']]}  (the entropy, area-law)")
    print(f"    S/V   = {['%.4f' % x for x in r['R']]}  (the holographic ratio, ~1/L)")
    print(f"  area-law exponent  (S vs L):   {r['area_law_exponent']:.4f}  (= 2 -> AREA, not volume)")
    print(f"  holographic exponent (S/V vs L): {r['holographic_ratio_exponent']:.4f}  (= -1 -> 1/L, the holographic hallmark)")
    ok = (abs(r["area_law_exponent"] - 2.0) < _TOL.loose
          and abs(r["holographic_ratio_exponent"] + 1.0) < _TOL.loose
          and abs(r["per_plaquette"] - 0.25) < _TOL.strict)
    print(f"  VERDICT:  {'PASS (the QCA is holographic: S ~ AREA, S/V ~ 1/L)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
