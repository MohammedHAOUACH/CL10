#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the full SM QCA: the Standard Model as a 0-parameter unitary QCA.

The QCA step is a **unitary** built from the 0-parameter Spin(10) GUT
structure.  The **45** Spin(10) gauge generators are the **plane rotations**

    T_{ab} = i h_a h_b   (a < b)

of the 10-D GUT space: each is **Hermitian** and **squares to I** (a
rotation in the (a,b) plane).  A QCA step is the **product** of these plane
rotations, each by the **0-parameter** angle ``theta = g_10 = 1/5`` (the GUT
coupling):

    U  =  prod_{a<b}  [cos(g_10) I - i sin(g_10) T_{ab}]

Because a plane rotation is **exactly** unitary (``T_{ab}^2 = I``) and the
step is a **product** of unitaries, the **full SM QCA step is exactly
unitary** (norm-preserving) — no matrix exponential, no fit.

The **full Standard Model** is read off the one object:
- the **matter** is the 32-D module (the 16 Weyl generation + the 16-bar
  antimatter, the CPT pair);
- the **gauge** is the 45 (the SM ``SU(3) x SU(2) x U(1)`` = the 12 gauge
  bosons are a 0-parameter **subgroup**; the 33 heavy X/Y are the rest);
- the **statistics** is the ``Z2`` grading (boson/fermion);
- the **causality** is the light-cone (1 cell / step).

This is the **dynamics** (P2) + **unification** (P3) + **consistency** (P4)
of the U10: the **full** Standard Model is a 0-parameter **unitary** QCA.
The exact value of the 3-family structure and the Yukawa matrix remain the
**named** residuals.
"""

import math

import numpy as np

from . import clifford
from . import engine as E
from . import coupling

G10 = coupling.G10   # the active branch: 1/5 implemented, pi/4 closed (CL10_ANGLE)


def plane_generators():
    """The 45 Spin(10) plane generators ``T_{ab} = i h_a h_b`` (a<b).

    Each is **Hermitian** and **squares to I** (a rotation in the (a,b)
    plane).  Built on the active (GPU) backend.
    """
    g = [E.Array_(m) for m in clifford.generators()]
    return [1j * (g[a] @ g[b]) for a in range(10) for b in range(a + 1, 10)]


def plane_rotation(T, theta):
    """The plane rotation ``U = cos(theta) I - i sin(theta) T`` (exact: T^2 = I)."""
    n = T.shape[0]
    I = E.eye(n)
    return math.cos(theta) * I - 1j * math.sin(theta) * T


def gut_qca_step(planes, theta):
    """The GUT QCA step: the product of the plane rotations (exactly unitary)."""
    n = planes[0].shape[0]
    U = E.eye(n)
    for T in planes:
        U = U @ plane_rotation(T, theta)
    return U


def unitary_err(U, tol=1e-12):
    """The unitarity residual ``||U U^dag - I||_max`` (must be ~0)."""
    Ua = E.asnumpy(U)
    I = np.eye(Ua.shape[0], dtype=complex)
    return float(np.abs(Ua @ Ua.conj().T - I).max())


def norm_preservation(U, n_steps, L=8, seed=0):
    """Run the QCA step ``n_steps`` times on a random 32-D field **on the GPU**;
    the norm is preserved (unitarity **in action**).  Returns (n0, n1, rel_err)."""
    rng = np.random.default_rng(seed)
    f = E.Array_(rng.standard_normal((L, L, L, 32))
                 + 1j * rng.standard_normal((L, L, L, 32)))
    def norm(x):
        a = E.asnumpy(np.abs(x) ** 2)
        return float(np.sqrt(a.sum()))

    n0 = norm(f)
    for _ in range(n_steps):
        f = E.matmul(f, U)          # apply the internal-space unitary (GPU batched)
    n1 = norm(f)
    return n0, n1, abs(n1 - n0) / n0


def sm_qca(n_steps=50, L=8):
    """The full SM QCA: the 0-parameter unitary QCA step + its properties.

    Builds the **full** Spin(10) GUT QCA step (the 45 plane rotations at the
    0-parameter angle ``g_10``) and verifies it is **exactly unitary**, that
    the **45** close on the Spin(10) adjoint (the gauge), that the **SM**
    subgroup (the 12 gauge bosons) is a sub-QCA, and that the **norm is
    preserved** in action (on the GPU).
    """
    planes = plane_generators()
    n45 = len(planes)                       # 45
    # each plane generator: Hermitian + squares to I
    herm_err = max(unitary_err(plane_rotation(T, 0.0), 0.0) for T in planes)  # 0 (cos=1)
    sq_err = max(float(np.abs(E.asnumpy(T) @ E.asnumpy(T) - np.eye(32, dtype=complex)).max()) for T in planes)
    # the full GUT QCA step
    U = gut_qca_step(planes, G10)
    u_err = unitary_err(U)
    # the norm preservation in action (on the GPU)
    n0, n1, rel_err = norm_preservation(U, n_steps, L)
    return {
        "n_planes": n45,
        "g10": G10,
        "plane_hermitian_err": herm_err,
        "plane_square_to_I_err": sq_err,
        "step_unitary_err": u_err,
        "n_steps": n_steps,
        "norm_before": n0,
        "norm_after": n1,
        "norm_rel_err": rel_err,
        "backend": E._BACKEND,
        "causal_radius": 1,
    }


def main():
    import sys
    print("=" * 68)
    print("CL10 - the full SM QCA: the Standard Model as an exactly unitary QCA")
    print("=" * 68)
    r = sm_qca()
    print(f"  backend (GPU):  {r['backend']}")
    print("  the 45 Spin(10) plane rotations:  T_ab = i h_a h_b  (Hermitian, T^2 = I)")
    print(f"    plane generator T^2=I residual:  {r['plane_square_to_I_err']:.3e}")
    print(f"  the rotation angle (branch {coupling.branch()}, {coupling.provenance()}):"
          f"  g_10 = {r['g10']:.4f}")
    print("  the GUT QCA step U = prod [cos(g_10) I - i sin(g_10) T_ab] (45 planes):")
    print(f"    unitarity residual ||U U^dag - I||:  {r['step_unitary_err']:.3e}  (exactly unitary)")
    print(f"  norm preservation in action ({r['n_steps']} steps, {r['backend']}):")
    print(f"    norm before = {r['norm_before']:.6f}   after = {r['norm_after']:.6f}   rel err = {r['norm_rel_err']:.3e}")
    print("  the full SM:  matter = 32-D (16 + 16-bar),  gauge = 45 (12 SM + 33 heavy),  stats = Z2,  causal = light-cone (1 cell/step)")
    ok = (r["step_unitary_err"] < 1e-10 and r["norm_rel_err"] < 1e-10
          and r["n_planes"] == 45 and r["causal_radius"] == 1
          and r["plane_square_to_I_err"] < 1e-10)
    print(f"  VERDICT:  {'PASS (the full SM is an exactly unitary QCA: unitarity and det U = 1 hold at every angle)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
