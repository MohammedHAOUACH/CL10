#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the QCA spectrum: the **one** 0-parameter object and its invariants.

The **single** object of the theory is the **GUT QCA step** ``U`` — the 32x32
product of the 45 ``Spin(10)`` plane rotations at the angle of
``cl10.coupling`` (the implemented branch ``g_10 = 1/5``; the angle itself is a
convention, the paper's O4).  The
**spectrum** of ``U`` — the **32 levels** of

    H = -i log(U),   levels in (-pi, pi]

— is the **one computed 0-parameter object** from which the dynamics is read
off.  This module computes that spectrum and its **0-parameter invariants**
(all verified on the array backend, no fit):

- **the 0-point = 0** (O1).  The 32 levels are **exactly +-paired** (16
  pairs of ``+E``/``-E``, to machine precision, **for every** ``g_10``), so
  the QCA vacuum 0-point ``sum(levels) = arg det(U) = 0``.  This is the
  **dynamical face** of the 512/512 ``Z2`` s-mechanism.
- **the SL(32) identity** (O1).  Every plane rotation
  ``R_ab = exp(-i g T_ab) = cos(g) I - i sin(g) T_ab`` is a **special**
  unitary (``det = 1``), because every ``T_ab = i h_a h_b`` is a
  **traceless** involution (``16``/``16`` ``+1``/``-1`` eigenvalues,
  ``tr = 0``) on the 32-D Weyl module.  The step is a **product** of
  special unitaries, so ``det(U) = 1`` **for every** ``g_10`` (and for
  **any** set of the 45 angles).  The 0-point = 0 is therefore a
  **group-theoretic identity of the root**, not a counting accident.
- **the CPT (16/16-bar) split** (O3/O5).  ``U`` **commutes with the
  pseudoscalar** (the 45 are even in the ``Z2`` grading), so it preserves
  the 16/16-bar (matter/antimatter) eigenspaces; the 0-point of **each**
  sector is 0.
- **the S-matrix** (O9): the step **is** the 0-parameter S-matrix (the
  't Hooft program) — ``S = U``.

The **nine residuals O1-O9** (the named exact values the 0-parameter
structure does not fix) are the **nine projections of this one object**
(see the publication, the 9-to-1 reduction).
"""

import numpy as np

from . import clifford
from . import engine as E
from . import sm_qca
from . import coupling

# the angle of the object (the only input): the branch of cl10.coupling
G10 = coupling.G10


def _asnp(x):
    """The numpy view of an engine array (32x32 is CPU-size)."""
    return E.asnumpy(x)


def levels(g10=None):
    """The 32-level QCA spectrum (the phases of ``U``), in ``(-pi, pi]``,
    sorted ascending.  ``U`` is the GUT step at ``g10`` (0-parameter
    ``1/5`` by default)."""
    g10 = G10 if g10 is None else g10
    U = _asnp(sm_qca.gut_qca_step(sm_qca.plane_generators(), g10))
    return np.sort(np.angle(np.linalg.eigvals(U)))


def cpt_block_phases(g10=None):
    """The CPT (16/16-bar) mirror of the GUT step, as two 16-element
    eigenphase lists (each sorted ascending).  ``U`` **commutes** with the
    pseudoscalar, so it is block-diagonal in the pseudoscalar (16/16-bar)
    basis: the **16 (matter)** block and the **16-bar (antimatter)** block.
    CPT is phase conjugation, so the 16-bar multiset is the **mirror** of
    the 16 (``omega_16bar = -omega_16``) — the 0-parameter CPT mirror.

    Note: each block carries a **mix** of + and - phases (the 16 has 9
    negative / 7 positive); the mirror is between the two **blocks**, not
    between the + and - halves of the full 32-level spectrum."""
    g10 = G10 if g10 is None else g10
    U = _asnp(sm_qca.gut_qca_step(sm_qca.plane_generators(), g10))
    w = _asnp(clifford.pseudoscalar())
    ev, V = np.linalg.eigh(1j * w)            # Hermitian, eigenvalues +-1
    idx16 = np.where(np.abs(ev + 1) < 1e-9)[0]    # iw = -1 <=> w = +i = the 16
    idx16b = np.where(np.abs(ev - 1) < 1e-9)[0]   # iw = +1 <=> w = -i = the 16-bar
    ph16 = np.sort(np.angle(np.linalg.eigvals(V[:, idx16].conj().T @ U @ V[:, idx16])))
    ph16b = np.sort(np.angle(np.linalg.eigvals(V[:, idx16b].conj().T @ U @ V[:, idx16b])))
    return {"ph16": ph16, "ph16b": ph16b}


def qca_spectrum(g10=None):
    """The one 0-parameter object (the 32-level QCA spectrum) and its
    0-parameter invariants: the 0-point, the SL(32) identity, the +/--;
    pairing, the CPT (16/16-bar) split, and the commutation with the
    pseudoscalar."""
    g10 = G10 if g10 is None else g10
    planes = sm_qca.plane_generators()
    U = _asnp(sm_qca.gut_qca_step(planes, g10))

    # the 32 levels (the phases of U)
    ph = np.sort(np.angle(np.linalg.eigvals(U)))
    n = len(ph)

    # the 0-point: the sum of the 32 levels (= arg det U)
    zero_point = float(ph.sum())

    # the SL(32) identity: det(U) = 1  (every rotation is a special unitary)
    detU = np.linalg.det(U)
    det_minus_1 = float(abs(detU - 1.0))
    det_arg = float(np.angle(detU))

    # the +/-- pairing: pair the sorted levels smallest <-> largest
    pair_resid = float(max(abs(ph[i] + ph[n - 1 - i]) for i in range(n // 2)))

    # the 45 generators are traceless on the 32-D module (16/16 +1/-1)
    h = [_asnp(g) for g in clifford.generators()]
    t45 = [1j * (h[a] @ h[b]) for a in range(10) for b in range(a + 1, 10)]
    max_tr = float(max(abs(np.trace(T)) for T in t45))

    # the CPT (16/16-bar) split: U commutes with the pseudoscalar, and the
    # 0-point of each sector is 0
    w = _asnp(clifford.pseudoscalar())
    iw = 1j * w  # Hermitian, eigenvalues +-1
    ev, V = np.linalg.eigh(iw)
    idx16 = np.where(np.abs(ev + 1) < 1e-9)[0]   # iw=-1 <=> w=+i = the 16 Weyl
    idx16b = np.where(np.abs(ev - 1) < 1e-9)[0]  # iw=+1 <=> w=-i = the 16-bar
    B16, B16b = V[:, idx16], V[:, idx16b]
    zp16 = float(np.sum(np.angle(np.linalg.eigvals(B16.conj().T @ U @ B16))))
    zp16b = float(np.sum(np.angle(np.linalg.eigvals(B16b.conj().T @ U @ B16b))))
    comm_w = float(np.abs(U @ w - w @ U).max())

    return {
        "g10": g10,
        "n_levels": n,
        "spectrum": ph,
        "zero_point": zero_point,
        "det_U_minus_1": det_minus_1,
        "det_U_arg": det_arg,
        "sl32": bool(det_minus_1 < 1e-10),
        "pairing_resid": pair_resid,
        "max_trace_T45": max_tr,
        "cpt_0point_16": zp16,
        "cpt_0point_16bar": zp16b,
        "commutes_with_pseudoscalar": comm_w,
        "max_level": float(abs(ph).max()),
    }


def verify(tol=1e-10):
    """0-parameter checks on the one object: the 0-point = 0 (the +/--
    pairing), the SL(32) identity (det U = 1, the 45 traceless), and the CPT
    (16/16-bar) split (U commutes with the pseudoscalar, each sector's
    0-point = 0)."""
    s = qca_spectrum()
    ok = (
        abs(s["zero_point"]) < tol
        and s["sl32"]
        and s["pairing_resid"] < tol
        and s["max_trace_T45"] < tol
        and s["commutes_with_pseudoscalar"] < 1e-9
        and abs(s["cpt_0point_16"]) < 1e-8
        and abs(s["cpt_0point_16bar"]) < 1e-8
    )
    return {"spectrum": s, "ok": bool(ok)}


def main():
    import sys
    print("=" * 68)
    print("CL10 - the QCA spectrum: the one 0-parameter object")
    print("=" * 68)
    r = verify()
    s = r["spectrum"]
    print(f"  the object : the 32-level spectrum of U (g_10 = {s['g10']})")
    print(f"  levels     : {s['n_levels']}  (max |level| = {s['max_level']:.4f})")
    print(f"  0-point    : {s['zero_point']:.3e}   (= 0: the 32 levels are +/---paired)")
    print(f"  SL(32)     : det(U) - 1 = {s['det_U_minus_1']:.3e}, "
          f"arg = {s['det_U_arg']:.3e}   (det U = 1 for every g_10)")
    print(f"  +/--pair   : max residual = {s['pairing_resid']:.3e}")
    print(f"  45 trace   : max |tr(T_ab)| = {s['max_trace_T45']:.3e}   (16/16 +-1)")
    print(f"  CPT split  : 16/16-bar, 0-point = {s['cpt_0point_16']:.2e} / "
          f"{s['cpt_0point_16bar']:.2e},  [U, w] = {s['commutes_with_pseudoscalar']:.2e}")
    print(f"  VERDICT    : {'PASS (the one object: 0-point = 0, det U = 1, CPT)' if r['ok'] else 'FAIL'}")
    sys.exit(0 if r["ok"] else 1)


if __name__ == "__main__":
    main()
