#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the local QCA: the causal radius **implemented** (the 't Hooft
program).

**The 0-parameter construction.**  The **local** QCA step on the
$(L,L,L)\\times 32$ Planck lattice (the $(1,3)$ lattice, the 32-D Weyl
internal space per cell) is a **product of 2-cell (local) unitaries**:

1. **internal layer**: the 0-parameter 45 Spin(10) plane rotations
   ``[cos(g_10) I - i sin(g_10) T_ab]`` (V23) applied to **each** cell;

2. **bond layer** (the x, y, z directions **cycle** with the time step — one
   direction per step, so a signal gets **one cell of budget per step**;
   **staggered** — the 't Hooft staggering, the even / odd bonds
   alternating with the time step): the **0-parameter swap-rotation**

       U_bond = cos(g_10) I - i sin(g_10) P,

   where ``P`` is the **0-parameter swap** of the 32-D internal spaces of
   two **adjacent** cells (a permutation — **no** angle).  ``P`` is a
   Hermitian involution (``P = P^dagger``, ``P^2 = I``), so ``U_bond`` is
   **exactly unitary**.

The **only** number of the whole step is ``g_10`` (the 0-parameter GUT
coupling, a **number of the root** — V38).  The **causal radius 1** is
therefore a **property of the dynamics** (the step is a product of 2-cell
unitaries), **not** an axiom.

**The result (the light cone realized).**  A state localized on one cell
stays **inside** the **light cone** (the Manhattan ball of radius k after k
steps) — **causality is implemented** — and the **front** of the support
**reaches** the light cone **exactly** (the front speed is **c = 1**).
The light cone is a **structural** property (it does not depend on the
0-parameter coupling).
"""

import math

import numpy as np

from . import engine as E
from . import sm_qca
from . import toepred

N = 32                            # the internal (Weyl) space per cell
G10 = math.sqrt(toepred.ALPHA_10) # = 1/5, the 0-parameter GUT angle


def internal_unitary():
    """The 0-parameter internal (1-cell) unitary: the product of the 45
    Spin(10) plane rotations at the 0-parameter angle ``g_10`` (V23).
    A (32,32) unitary."""
    return sm_qca.gut_qca_step(sm_qca.plane_generators(), G10)


def internal_layer(fields, U):
    """Apply the 32x32 internal unitary ``U`` to **every** cell of
    ``fields`` (L,L,L,32) (the GPU matmul on the last axis)."""
    return np.einsum("...i,ij->...j", fields, U)


def _bond_rotate(a, b, theta):
    """The 0-parameter swap-rotation of two bond-end states (each
    ``(..., 32)``): ``U_bond = cos(theta) I - i sin(theta) P`` with ``P``
    the 0-parameter swap of the two cell states.  Returns the two rotated
    bond ends."""
    c = math.cos(theta)
    s = math.sin(theta)
    return (c * a - 1j * s * b), (-1j * s * a + c * b)


def bond_layer(fields, ax, even=True, theta=G10):
    """One **staggered** 0-parameter bond layer along the spatial axis
    ``ax`` (0 = x, 1 = y, 2 = z): the swap-rotation on the **even** bonds
    ``(0,1), (2,3), ...`` (the ``even=True``) or the **odd** bonds
    ``(1,2), (3,4), ...`` (the ``even=False``, the 't Hooft staggering).
    ``fields`` has shape (L,L,L,32) (a torus)."""
    m = np.moveaxis(fields, ax, 0)
    if not even:
        m = np.roll(m, -1, axis=0)
    L = m.shape[0]
    r = m.reshape(L // 2, 2, *m.shape[1:])
    a, b = r[:, 0], r[:, 1]
    a2, b2 = _bond_rotate(a, b, theta)
    r[:, 0] = a2
    r[:, 1] = b2
    out = r.reshape(*m.shape)
    if not even:
        out = np.roll(out, 1, axis=0)
    return np.moveaxis(out, 0, ax)


def local_step(fields, U_internal, t, theta=G10):
    """One **local** QCA **time step** (``t`` = 0, 1, 2, ...): the internal
    layer + **one** staggered bond layer, along the direction ``t mod 3``
    (the x, y, z directions **cycle** with the time step — so a signal
    gets **one** cell of budget **per** step, in **one** direction at a
    time).  The staggering alternates with ``t``: the even bonds at even
    ``t``, the odd bonds at odd ``t`` — the 't Hooft staggering."""
    ax = t % 3
    even = (t % 2 == 0)
    f = internal_layer(fields, U_internal)
    f = bond_layer(f, ax, even, theta)
    return f


def run(L=16, n_steps=8, x0=(8, 8, 8), theta=G10, bonds=True, seed=0,
        tol=1e-12):
    """Run the local QCA: a state **localized** on the single cell ``x0``
    (a random 32-D vector), for ``n_steps`` time steps.  ``bonds=False``
    runs the internal layer **only** (no 2-cell coupling — the control).
    Returns the per-step support (the cells with |amplitude| > ``tol``)
    and the data (the light-cone radius of the support at each step)."""
    rng = np.random.default_rng(seed)
    f = E.zeros((L, L, L, N))
    v = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    f[x0] = E.Array_(v / np.linalg.norm(v))
    U = internal_unitary()

    def _cell_norms(g):
        a = E.asnumpy(g)
        return np.sqrt(np.sum(np.abs(a) ** 2, axis=-1))

    steps = []
    for k in range(1, n_steps + 1):
        f = local_step(f, U, k - 1, theta) if bonds else internal_layer(f, U)
        nrm = _cell_norms(f)
        support = np.argwhere(nrm > tol)
        # the light-cone radius: the max Manhattan distance (the minimum
        # image on the torus) of the support from x0
        d = []
        for (x, y, z) in support:
            dd = 0
            for a0, a1 in ((x, x0[0]), (y, x0[1]), (z, x0[2])):
                delta = abs(a0 - a1)
                dd += min(delta, L - delta)
            d.append(dd)
        steps.append({
            "k": k,
            "n_support": int(len(support)),
            "max_radius": int(max(d)) if d else 0,
            "n_outside": int(sum(1 for r in d if r > k)),
        })
    return {
        "backend": E._BACKEND,
        "L": int(L),
        "n_steps": int(n_steps),
        "theta": float(theta),
        "g10": float(G10),
        "bonds": bool(bonds),
        "steps": steps,
        "final_norm": float(np.linalg.norm(E.asnumpy(f).ravel())),
    }


def verify(tol=1e-9):
    """0-parameter checks:
    1. the **local** step is **unitary** (the norm is preserved);
    2. the **light cone is realized**: a localized excitation stays **inside
       the light cone** (the Manhattan ball of radius k after k steps — the
       causal bound of the $(1,3)$) and its **front** reaches the light cone
       **exactly** (the front speed is **c = 1**);
    3. the light cone is a **structural** property: the **same** cone (the
       ball of radius k, the same front) is realized for **different**
       0-parameter angles — the cone is carried by the **2-cell operator
       structure**, not by the value of the angle;
    4. the bond layers are **non-trivial** (the support grows, the internal
       layer alone does not move the state).
    """
    L = 16
    n_steps = 8
    k_struct = 6
    x0 = (8, 8, 8)
    r_full = run(L, n_steps, x0, theta=G10, bonds=True, seed=0)
    r_int = run(L, n_steps, x0, theta=G10, bonds=False, seed=0)
    r_alt = run(L, k_struct, x0, theta=0.3, bonds=True, seed=0, tol=1e-9)

    # 1. the unitary (the norm is preserved, the initial norm = 1)
    unitary_err = abs(r_full["final_norm"] - 1.0)

    # 2. the light cone is realized (causality + the front, c = 1)
    causal_ok = True
    front_reached = []
    for s in r_full["steps"]:
        if s["n_outside"] != 0 or s["max_radius"] > s["k"]:
            causal_ok = False
        front_reached.append(s["max_radius"] == s["k"])
    causal_ok = causal_ok and all(front_reached)
    front_speed = (r_full["steps"][-1]["max_radius"] / r_full["steps"][-1]["k"])

    # 3. the structural (the same cone for a different 0-parameter angle)
    structural_ok = all(
        r_alt["steps"][i]["n_outside"] == 0
        and r_alt["steps"][i]["max_radius"] == r_full["steps"][i]["max_radius"]
        == r_alt["steps"][i]["k"]
        for i in range(k_struct)
    )

    # 4. the non-trivial (the bond layers move the state, the internal
    #    layer alone does not)
    nontrivial_ok = (r_full["steps"][-1]["n_support"] > 1
                     and r_int["steps"][-1]["n_support"] == 1)

    ok = (unitary_err < tol and causal_ok and structural_ok and nontrivial_ok)
    return {
        "backend": r_full["backend"],
        "unitary_err": float(unitary_err),
        "causal_ok": bool(causal_ok),
        "front_speed": float(front_speed),
        "structural_ok": bool(structural_ok),
        "nontrivial_ok": bool(nontrivial_ok),
        "support_growth": [s["n_support"] for s in r_full["steps"]],
        "internal_only_growth": [s["n_support"] for s in r_int["steps"]],
        "g10": r_full["g10"],
        "theta_tested": r_alt["theta"],
        "ok": bool(ok),
    }


def main():
    import sys
    r = verify()
    print("=" * 68)
    print("CL10 - the local QCA: the causal radius IMPLEMENTED")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter angle:  g_10 = {r['g10']}  (the 45 rotations + the swap-rotations)")
    print(f"  the unitary (the norm is preserved):  residual = {r['unitary_err']:.2e}")
    print(f"  the light cone REALIZED (the front speed):  c = {r['front_speed']:.3f}  (1 cell / step)")
    print(f"    the support after k steps (the ball of radius k):  {r['support_growth']}")
    print(f"    the internal layer alone (no 2-cell coupling):  {r['internal_only_growth']}  (the state does not move)")
    print(f"  the light cone is STRUCTURAL (independent of g_10, tested at {r['theta_tested']}):  {r['structural_ok']}")
    print(f"  VERDICT:  {'PASS (the causal radius 1 is a property of the local dynamics, not an axiom)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
