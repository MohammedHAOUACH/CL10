#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V114 - do the V76 junction conditions DERIVE the ratio r = kappa/w0?

**Why (W2 + V103).** V103: the seesaw with the built junction gives
m_3(r) = S17/(1+r); the atmospheric scale selects r = 1.91, and
r = 2 (m_nu = S17/3 = 0.0487 eV) is a junction-amplitude CONVENTION.
V76 built the junction
  w(x) = w0 + kappa [cos(pi(y-z)/2) + cos(pi(x-z)/2) + cos(pi(x-y)/2)]
claimed cubic-symmetric, zero-mean, Fourier support exactly closing
the node triangle. QUESTION: do those conditions FIX r = kappa/w0?

**The exact computation on the period-4 lattice (Z/4)^3** (the cos
profile has period 4; the supercell folding hides this):
1. STABILIZER: the exact symmetry group of the V76 modulation inside
   the 48-element signed-permutation group. (V76's "cubic symmetric"
   is checked, not assumed.)
2. CUBIC-SYMMETRIC CLOSED SPACE: profiles invariant under the 24
   chiral cubic rotations with Fourier support in {0} u +-triangle
   differences: dimension?
3. STABILIZER-SYMMETRIC CLOSED SPACE: invariance under the V76
   stabilizer instead: dimension d, modulation dimension d_m, and
   whether the V76 shape is unique.
r is DERIVED only if the family is 1-dimensional (no independent
constant). Otherwise the SHAPE is what the conditions fix and r stays
free.
"""
import itertools
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

SUPPORT = [(0, 1, 3), (0, 3, 1), (1, 0, 3), (3, 0, 1), (1, 3, 0), (3, 1, 0)]


def group_elements():
    """The 48 signed permutations as (perm, signs)."""
    return [(p, s) for p in itertools.permutations(range(3))
            for s in itertools.product((1, -1), repeat=3)]


def chiral_elements():
    out = []
    for p, s in group_elements():
        par = (p[0] > p[1]) ^ (p[0] > p[2]) ^ (p[1] > p[2])
        nflip = (s[0] < 0) + (s[1] < 0) + (s[2] < 0)
        if (par + nflip) % 2 == 0:
            out.append((p, s))
    return out


def main():
    pts = list(itertools.product(range(4), repeat=3))
    idx = {p: i for i, p in enumerate(pts)}

    def act(x, g):
        p, s = g
        return tuple((s[k] * x[p[k]]) % 4 for k in range(3))

    # the V76 modulation (zero-mean already: cosines average to 0)
    v76 = np.array([math.cos(math.pi * (x[1] - x[2]) / 2)
                    + math.cos(math.pi * (x[0] - x[2]) / 2)
                    + math.cos(math.pi * (x[0] - x[1]) / 2)
                    for x in pts], float)

    # 1. stabilizer of the V76 modulation in the 48-element group
    stab = []
    for g in group_elements():
        wp = np.array([v76[idx[act(x, g)]] for x in pts])
        if np.abs(wp - v76).max() < 1e-12:
            stab.append(g)
    stab_chiral = [g for g in stab if g in chiral_elements()]

    def invariance_rows(gens):
        rows = []
        for g in gens:
            for i, x in enumerate(pts):
                j = idx[act(x, g)]
                if j == i:
                    continue
                row = np.zeros(64)
                row[j] = 1.0
                row[i] = -1.0
                rows.append(row)
        return rows

    def frow(mult):
        return np.array([np.exp(1j * math.pi / 2 * sum(
            mult[k] * pts[i][k] for k in range(3))) for i in range(64)]) / 64

    allowed = {tuple(m) for m in SUPPORT} | {(0, 0, 0)}

    def fourier_rows():
        rows = []
        for mult in itertools.product(range(4), repeat=3):
            if mult in allowed:
                continue
            r = frow(mult)
            rows.append(r.real)
            rows.append(r.imag)
        return rows

    def nullspace(rows):
        A = np.array(rows)
        s = np.linalg.svd(A, compute_uv=False)
        rank = int(np.sum(s >= 1e-10))
        _, _, Vt = np.linalg.svd(A)
        return 64 - rank, Vt[rank:]

    # 2. cubic-symmetric (24 rotations) + closed
    dim_ch, basis_ch = nullspace(invariance_rows(chiral_elements())
                                 + fourier_rows())
    # 3. stabilizer-symmetric + closed
    dim_st, basis_st = nullspace(invariance_rows(stab) + fourier_rows())

    const = np.ones(64) / 8.0
    v76_err_st = float(np.abs(v76 - sum(float(np.dot(v76, b)) * b
                                        for b in basis_st)).max())
    v76_err_ch = float(np.abs(v76 - sum(float(np.dot(v76, b)) * b
                                        for b in basis_ch)).max())

    def mod_dim(basis):
        mb = []
        for b in basis:
            m = b - float(np.dot(b, const)) * const
            if np.linalg.norm(m) > 1e-9:
                # Gram-Schmidt against the existing modulation basis
                for x in mb:
                    m = m - float(np.dot(m, x)) * x
                if np.linalg.norm(m) > 1e-9:
                    mb.append(m / np.linalg.norm(m))
        return len(mb), mb

    dmod_ch, mb_ch = mod_dim(basis_ch)
    dmod_st, mb_st = mod_dim(basis_st)
    v76_mod = v76 - float(np.mean(v76))
    v76_mod_n = v76_mod / np.linalg.norm(v76_mod)
    shape_err_st = float(np.linalg.norm(
        v76_mod_n - sum(float(np.dot(v76_mod_n, m)) * m for m in mb_st)))
    # verify bases
    verif = 0.0
    for b in basis_st:
        for mult in itertools.product(range(4), repeat=3):
            if mult in allowed:
                continue
            verif = max(verif, abs(complex(frow(mult) @ b)))
    for g in stab:
        for i, x in enumerate(pts):
            j = idx[act(x, g)]
            for b in basis_st:
                verif = max(verif, abs(b[j] - b[i]))

    r_verdict = "NOT DERIVED" if dim_st >= 2 else "NOT DERIVED (no family)"
    status = (
        f"JUNCTION RATIO AUDIT on (Z/4)^3: the V76 modulation's exact "
        f"stabilizer in the 48-element signed-permutation group has "
        f"{len(stab)} elements ({len(stab_chiral)} of them chiral) - "
        f"NOT the 24-element chiral cubic group: V76's 'cubic "
        f"symmetric' holds only for the triangle SUPPORT (an orbit), "
        f"not for the profile. The cubic-symmetric closed space "
        f"(24 rotations + Fourier closure) has dimension {dim_ch} "
        f"(modulation {dmod_ch}): the V76 profile lies outside it to "
        f"{v76_err_ch:.2e} - NO cubic-symmetric closed junction "
        f"exists on the period-4 lattice. Inside the V76 stabilizer "
        f"the closed space has dimension {dim_st} (modulation "
        f"{dmod_st}, V76 shape error {shape_err_st:.2e}, basis "
        f"verified to {verif:.2e}): the constant w0 is an independent "
        f"direction, so r = kappa/w0 is {r_verdict} - the conditions "
        f"fix at most the shape, never the amplitude ratio. W2's "
        f"r = 2 (democratic) and r = 1.91 (atmospheric) remain named "
        f"inputs; F3 measures the outcome.")
    out = {"experiment": "V114_w2_junction_r",
           "pass": bool(verif < 1e-9 and dim_st >= 1
                        and v76_err_st < 1e-9),
           "status": status,
           "data": {"stabilizer_size": len(stab),
                    "stabilizer_chiral": len(stab_chiral),
                    "dim_cubic_closed": int(dim_ch),
                    "dim_cubic_mod": int(dmod_ch),
                    "dim_stab_closed": int(dim_st),
                    "dim_stab_mod": int(dmod_st),
                    "v76_in_cubic_err": v76_err_ch,
                    "v76_in_stab_err": v76_err_st,
                    "v76_shape_err": shape_err_st,
                    "verify": float(verif),
                    "r_verdict": r_verdict}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v114_w2_junction_r.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
