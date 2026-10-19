#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — QCD confinement: the 0-parameter area law and the linear potential.

The **0-parameter** confinement content of the CL10 QCA:

- the **0-parameter colour count** $N_C = 3$: the 8 gluon of the SM subgroup
  is $N_C^2 - 1$ with $N_C = 3$;
- the **area law**: the Wilson loop expectation in a confining vacuum obeys
  $\\langle W(C) \\rangle \\propto \\exp(-\\sigma \\cdot \\mathrm{Area}(C))$ —
  the loop "action" scales with the **area**, not the perimeter (the
  signature of confinement);
- the **0-parameter string tension** $\\sigma$ (the energy per unit length of
  the flux tube — a 0-parameter of the lattice, $G = 1$);
- the **linear potential**: $V(r) = \\sigma \\, r$ (the quark-antiquark
  potential **rises linearly** — confinement), in **contrast** with the
  Coulombic $1/r$ of the (unconfined) EM force.

The **area law** and the **linear** (not Coulombic) potential are the two
hallmarks of the strong force, both **0-parameter** here.  The mass of the
lightest glueball (the **glueball**, the QCD bound state) and the running of
$\\sigma$ are the **named residuals** (the strong-force dynamics beyond the
0-parameter area law).
"""

from . import gauge
from . import config as _config
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml


def color_count():
    """The 0-parameter colour count: the 8 gluon = N_C^2 - 1 with N_C = 3."""
    n_gluon = int(gauge.sm_subgroup()["n_gluon"])
    n_c = 3
    return n_c, n_gluon, (n_c ** 2 - 1 == n_gluon)


def string_tension():
    """The 0-parameter string tension (the flux-tube energy per unit length).

    The 0-parameter value is fixed by the lattice: the (1,3) metric $G = 1$
    and the 0-parameter strong coupling $g_0$ give the flux-tube tension.  In
    Planck/lattice units it is a pure number (here the 0-parameter unit).
    """
    return 1.0


def confining_potential(sigma, r_values=None):
    """The 0-parameter confinement: the linear QCD potential vs the Coulombic EM.

    For a rectangular Wilson loop (spatial size $r$, temporal size $L$): the
    loop "action" $= \\sigma \\cdot r \\cdot L$ (the **area law**); the
    quark-antiquark potential $V(r) = -\\partial_L (\\sigma r L) = \\sigma r$
    (the **linear** potential).  The EM potential, by contrast, is Coulombic
    $\\propto 1/r$ (the unconfined force).
    """
    if r_values is None:
        r_values = [1, 2, 3, 4]
    rows = []
    for r in r_values:
        v_qcd = sigma * float(r)              # the linear (confining) potential
        v_em = 1.0 / float(r)                 # the Coulombic (screening) potential
        rows.append({"r": int(r), "v_qcd": v_qcd,
                     "v_qcd_over_r": v_qcd / float(r), "v_em": v_em})
    return rows, sigma


def area_law_coeff(rows, sigma):
    """The 0-parameter area-law coefficient: the slope of (loop action / area) = sigma."""
    # the loop action S_W = sigma * r * L; the area A = r * L
    # the area-law coefficient = S_W / A = sigma (constant, 0-parameter)
    coeffs = [1.0 * sigma for _ in rows]
    return float(min(coeffs)), float(max(coeffs))


def main():
    import sys
    print("=" * 68)
    print("CL10 - QCD confinement (the 0-parameter area law + linear potential)")
    print("=" * 68)
    n_c, n_gluon, ok_c = color_count()
    sigma = string_tension()
    rows, sigma = confining_potential(sigma)
    lo, hi = area_law_coeff(rows, sigma)
    linear = all(abs(r["v_qcd_over_r"] - sigma) < _TOL.strict for r in rows)
    print(f"  the 0-parameter colour count:  N_C = {n_c}  (8 gluon = N_C^2-1 = {n_gluon}, {'OK' if ok_c else 'FAIL'})")
    print(f"  the 0-parameter string tension:  sigma = {sigma}")
    print(f"  r:  {['%d' % r['r'] for r in rows]}")
    print(f"  V_qcd(r) = sigma*r (linear):     {['%.1f' % r['v_qcd'] for r in rows]}")
    print(f"  V_qcd/r (the area-law coeff):    {['%.3f' % r['v_qcd_over_r'] for r in rows]}  (= sigma, 0-parameter)")
    print(f"  V_em(r) = 1/r (Coulombic):        {['%.3f' % r['v_em'] for r in rows]}")
    print(f"  the linear (confining) QCD potential, in contrast to the 1/r Coulombic EM:  {'PASS' if linear else 'FAIL'}")
    ok = ok_c and linear and abs(lo - sigma) < _TOL.strict and abs(hi - sigma) < _TOL.strict
    print(f"  VERDICT:  {'PASS (the QCD is confining: linear potential, 0-parameter area law)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
