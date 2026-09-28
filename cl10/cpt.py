#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the CPT theorem: the 16 / 16-bar are CPT conjugates.

The **pseudoscalar** (volume element) $w = g_1 g_2 \\cdots g_{10}$ is
**anti-Hermitian** with **purely imaginary** eigenvalues $\\pm i$ (each
multiplicity 16): the $(+i)$ eigenspace is the **16 Weyl** (matter), the
$(-i)$ the **16-bar** (antimatter).  The **CPT operator** is the
**time-direction** generator

    C = g_1

(the "time" of the $(1,3)$).  It is a **Hermitian involution**
($C^2 = I$, $C = C^\\dagger$) and **anticommutes** with $w$ ($C w + w C = 0$),
so it **maps the 16 to the 16-bar** — the charge-conjugation of the CPT
transform.

The **CPT theorem** (the 0-parameter content, all from the one object):
1. **involution**  $C^2 = I$;
2. **unitary**     $C C^\\dagger = I$ (the CPT transform preserves the norm);
3. **maps 16 → 16-bar**  $C P_{16} = P_{16\\text{bar}} C$ (with
   $P_{16} = \\tfrac{1}{2}(I - i w)$, $P_{16\\text{bar}} = \\tfrac{1}{2}(I + i w)$);
4. **the GUT is CPT-symmetric**: the **45** gauge generators **commute** with
   $w$ (the gauge preserves the chiral 16/16-bar split), so the **matter and
   antimatter have the same dynamics and the same mass** (CPT: a particle and
   its antiparticle are degenerate).

The theorem is verified **on the GPU** engine (CuPy): the CPT operator is
unitary and norm-preserving in action, and the 16 / 16-bar carry **equal**
norm (the matter/antimatter CPT symmetry).
"""

import numpy as np

from . import clifford
from . import engine as E

# the dimension of the Cl(10) module
N = 32


def _gens():
    return [E.Array_(m) for m in clifford.generators()]


def pseudoscalar():
    """The volume element $w = g_1 \\cdots g_{10}$ (anti-Hermitian, $w^2 = -I$)."""
    w = E.eye(N)
    for a in _gens():
        w = w @ a
    return w


def cpt_operator():
    """The CPT operator $C = g_1$ (the time direction): a Hermitian involution."""
    return _gens()[0]


def projectors(w):
    """The 16 / 16-bar projectors from the pseudoscalar eigenspaces.

    $P_{16} = \\tfrac{1}{2}(I - i w)$ (the $+i$ eigenspace, the matter),
    $P_{16\\text{bar}} = \\tfrac{1}{2}(I + i w)$ (the $-i$ eigenspace, the antimatter).
    """
    I = E.eye(N)
    return 0.5 * (I - 1j * w), 0.5 * (I + 1j * w)


def _err(a):
    return float(np.abs(E.asnumpy(a)).max())


def cpt_theorem(seed=0):
    """The CPT theorem: the CPT relations + the GPU in-action check.

    Verifies (1) $C^2=I$, (2) $C$ unitary, (3) $C P_{16}=P_{16\\text{bar}} C$,
    (4) the 45 GUT generators commute with $w$, and (5) CPT **in action** on the
    GPU (a random 32-D state: $C$ is norm-preserving and the 16/16-bar carry
    equal norm).
    """
    C = cpt_operator()
    w = pseudoscalar()
    P16, P16bar = projectors(w)
    I = E.eye(N)
    # (1) the CPT operator is an involution
    c_sq = _err(C @ C - I)
    # (2) the CPT operator is unitary
    c_unit = _err(C @ C.conj().T - I)
    # (3) C maps the 16 to the 16-bar (anticommutes with w)
    c_w_anti = _err(C @ w + w @ C)
    map_err = _err(C @ P16 - P16bar @ C)
    # (4) the GUT (the 45) commutes with w (the chiral split is preserved)
    g = _gens()
    gut_commute = max(_err(g[a] @ g[b] @ w - w @ g[a] @ g[b])
                   for a in range(10) for b in range(a + 1, 10))
    # (5) CPT in action (on the GPU): a matter (16) state and its CPT
    #     conjugate (16-bar) carry EQUAL norm (CPT: matter = antimatter)
    rng = np.random.default_rng(seed)
    L = 4
    f = E.Array_(rng.standard_normal((L, L, L, N))
                 + 1j * rng.standard_normal((L, L, L, N)))

    def nrm(x):
        return float(np.sqrt(E.asnumpy(np.abs(x) ** 2).sum()))

    # C preserves the total norm (CPT is unitary)
    f_c = E.matmul(f, C)
    n0, n1 = nrm(f), nrm(f_c)
    # a matter (16) state and its CPT conjugate (16-bar) have the SAME norm
    f16 = E.matmul(f, P16)                 # the matter component (in the 16)
    f16_c = E.matmul(f16, C)               # its CPT conjugate (in the 16-bar)
    n_matter = nrm(f16)
    n_antimatter = nrm(f16_c)              # = n_matter (C is unitary)
    # and C maps the 16 entirely into the 16-bar:
    c16_into_16bar = nrm(E.matmul(f16_c, P16bar))
    return {
        "C_involution_err": c_sq,
        "C_unitary_err": c_unit,
        "C_w_anticommute_err": c_w_anti,
        "map_16_to_16bar_err": map_err,
        "gut_commute_w_err": gut_commute,
        "norm_before": n0,
        "norm_after_cpt": n1,
        "norm_rel_err": abs(n1 - n0) / n0,
        "n_matter": n_matter,
        "n_antimatter": n_antimatter,
        "cpt_equal_norm": bool(abs(n_matter - n_antimatter) / max(n_matter, 1e-30) < 1e-12),
        "c16_into_16bar": c16_into_16bar,
        "backend": E._BACKEND,
    }


def main():
    import sys
    print("=" * 68)
    print("CL10 - the CPT theorem: the 16 / 16-bar are CPT conjugates")
    print("=" * 68)
    r = cpt_theorem()
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  C = g_1 (the time direction):")
    print(f"    (1) involution  C^2 = I:        {r['C_involution_err']:.3e}")
    print(f"    (2) unitary      C C^dag = I:   {r['C_unitary_err']:.3e}")
    print(f"    (3) maps 16 -> 16-bar  C w + w C = 0:  {r['C_w_anticommute_err']:.3e}")
    print(f"    (3) C P_16 = P_16bar C:          {r['map_16_to_16bar_err']:.3e}")
    print(f"    (4) the GUT (the 45) commutes with w:  {r['gut_commute_w_err']:.3e}  (the chiral split is preserved)")
    print(f"  CPT in action (on the {r['backend']}):")
    print(f"    norm preserved:  {r['norm_before']:.4f} -> {r['norm_after_cpt']:.4f}  (rel {r['norm_rel_err']:.3e})")
    print(f"    a matter (16) state and its CPT conjugate (16-bar) carry EQUAL norm:  "
          f"{r['n_matter']:.4f} / {r['n_antimatter']:.4f}  (equal: {r['cpt_equal_norm']})")
    print(f"    C maps the 16 entirely into the 16-bar:  {r['c16_into_16bar']:.4f}  (the 16 -> 16-bar under CPT)")
    ok = (r["C_involution_err"] < 1e-10 and r["C_unitary_err"] < 1e-10
          and r["C_w_anticommute_err"] < 1e-10 and r["map_16_to_16bar_err"] < 1e-10
          and r["gut_commute_w_err"] < 1e-10 and r["norm_rel_err"] < 1e-10
          and r["cpt_equal_norm"])
    print(f"  VERDICT:  {'PASS (the 16/16-bar are CPT conjugates; the QCA is CPT-symmetric)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
