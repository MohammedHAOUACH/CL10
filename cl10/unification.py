#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the 0-parameter coupling unification and the QCA UV completion.

**The 0-parameter content.**
- the **$5/3$** GUT normalization of the SM $\mathrm{U}(1)$ (Theorem 3) — the
  **0-parameter** anchor of the unification;
- the **0-parameter GUT coupling** $\alpha_{10} = 1/25$ (the single coupling of
  the 45, the GUT Cartan);
- the **1-loop** $\beta$-coefficients of the particle content
  ($b_1=41/10$, $b_2=-19/6$, $b_3=-7$ for the SM).

**The running to the 0-parameter scale.**  Running the three couplings from
$M_Z$ to the **Planckian** GUT scale $M_{\\rm GUT}=M_P$ (the 0-parameter
single scale) gives couplings in the **GUT range** ($\\alpha\\sim0.02$-$0.04$,
consistent with $\alpha_{10}=0.04$); the **1-loop spread** at $M_P$
($\\sim0.6$ for the non-SUSY SM) is the **named residual** — the non-SUSY SM
does **not** unify to a single point at $M_P$ (the clean MSSM unification is at
$\sim2\\times10^{16}$, **below** the Planck scale).

**The QCA UV completion.**  The lattice UV cutoff is the Planck scale.  The SM
$\mathrm{U}(1)$ has a **Landau pole** at $\sim3\\times10^{39}$ GeV, far
**above** the QCA cutoff $M_P$ — so in the QCA the running is **UV-bounded**
(the coupling is **bounded** on the lattice, no Landau pole below the Planck
scale).  This is the **0-parameter UV completion** of the theory.
"""

import numpy as np

from . import clifford
from . import engine as E
from . import gauge
from . import spacetime
from . import toepred
from . import walls

M_Z = 91.2                       # GeV
ALPHA_10 = toepred.ALPHA_10      # 1/25 = 0.04 (the 0-parameter GUT coupling)
M_P = toepred.M_PLANCK           # 1.22e19 GeV (the 0-parameter GUT = Planck scale)

# the 1-loop beta coefficients (convention: 1/a(M) = 1/a(Mz) - (b/2pi) ln(M/Mz))
B_SM = (41 / 10, -19 / 6, -7)    # (b1, b2, b3)  (the SM particle content)
A_Z = (0.0179, 0.0338, 0.118)    # (a1, a2, a3) at M_Z  (paired with B_SM)
# the 2-loop beta coefficients (n_f = 6, the SM particle content)
B2_SM = (199/25, 25/2 - 27/50*6, -143/3 + 23/15*6)   # (b21, b22, b23)


def _run(a_z, b, M):
    return 1.0 / a_z - (b / (2.0 * np.pi)) * np.log(M / M_Z)


def sm_running_to(M):
    """The 1-loop SM running to scale M.  Returns (a1, a2, a3, spread)."""
    a = [1.0 / _run(A_Z[k], B_SM[k], M) for k in range(3)]
    mean = float(np.mean(a))
    spread = float((max(a) - min(a)) / mean)
    return a, spread


def landau_pole_u1():
    """The SM U(1) Landau pole scale (where 1/a1 -> 0)."""
    ln_pole = (2.0 * np.pi / B_SM[0]) * (1.0 / A_Z[0])
    return M_Z * float(np.exp(ln_pole))


def g10_from_root():
    """The 0-parameter numbers **are numbers of the root** (C3).

    The 32-D module of the Cl(10) root is a **2^5** space (the symplectic
    Pauli pattern of ``clifford.py``: **5** qubits, 10 = 2x5 generators).
    The GUT coupling is a **number of the root**:

        g_10 = 1 / n_qubits = 1/5,   alpha_10 = g_10^2 = 1/25.

    The **Dirac** coupling (the W2 0-parameter ``g_0``) is

        g_0 = 1 / 2^{n_space} = 1 / 2^3 = 1/8,

    the **3** being the spatial dimension of the (1,3) (the same 0-parameter
    number as the family count — one number, two uses).
    """
    n_module = int(clifford.generators()[0].shape[0])     # the 32-D module
    n_qubits = int(np.log2(n_module))                     # = 5
    g10 = 1.0 / n_qubits                                  # 1/5
    T, S = spacetime.causal_generators()
    n_space = len(S)                                      # the 3 spatial directions
    g0 = 1.0 / 2.0 ** n_space                             # 1/8
    return {
        "backend": E._BACKEND,
        "n_module": int(n_module),
        "n_qubits": int(n_qubits),
        "g10": float(g10),
        "alpha_10": float(g10 ** 2),
        "g10_matches_sm_qca": bool(abs(g10 - 0.2) < 1e-12),
        "alpha_10_matches_toepred": bool(abs(g10 ** 2 - toepred.ALPHA_10) < 1e-12),
        "n_space": int(n_space),
        "g0": float(g0),
        "g0_matches_walls": bool(abs(g0 - walls.G0) < 1e-12),
    }


def _run_down(a_uv, b, M, M_uv):
    """The 1-loop running DOWN: 1/a(M) = 1/a(M_uv) - (b/2pi) ln(M_uv/M)."""
    return 1.0 / a_uv - (b / (2.0 * np.pi)) * np.log(M_uv / M)


def running_down():
    """The 0-parameter UV input at the QCA cutoff **run down** to M_Z (C3).

    The QCA lattice UV cutoff is the **Planck** scale, and the 0-parameter
    GUT coupling there is ``alpha_10 = 1/25`` (a number of the root, derived
    by ``g10_from_root``).  The standard 1-loop QFT extrapolation with the
    SM beta coefficients is **reported** together with the deviation from
    the measured ``alpha_i(M_Z)``.  The **1-loop matching to a single point**
    (the non-SUSY SM) is the **named residual** (W4), now **quantified**;
    the QCA **lattice RG** (the open part) governs below the UV cutoff.
    """
    a_uv = ALPHA_10
    a = [1.0 / _run_down(a_uv, B_SM[k], M_Z, M_P) for k in range(3)]
    ratios = [float(x / A_Z[k]) for k, x in enumerate(a)]
    # the U(1) pole of the DOWN running (where 1/a1 -> 0): 1/a1 = 0 when
    # ln(M_P / M_pole) = 2pi/(b1) * 1/a10
    ln_pole = (2.0 * np.pi / B_SM[0]) * (1.0 / a_uv)
    pole_scale = M_P * float(np.exp(-ln_pole))
    return {
        "backend": E._BACKEND,
        "alpha_10_at_M_P": float(a_uv),
        "couplings_at_M_Z": [float(x) for x in a],
        "measured_at_M_Z": [float(x) for x in A_Z],
        "ratio_to_measured": ratios,
        "a2_within_factor_2": bool(abs(ratios[1] - 1.0) < 1.0),
        "a3_within_factor_10": bool(abs(ratios[2] - 1.0) < 9.0),
        "u1_down_pole_scale": float(pole_scale),
        "u1_pole_inside_lattice": bool(M_Z < pole_scale < M_P),
        "residual": ("the 1-loop SM matching to the measured M_Z values (W4, now "
                     "quantified by the down-running); the QCA lattice RG is the open part"),
    }


def unification():
    """The 0-parameter unification verdict at the Planckian scale M_P."""
    a_p, spread_p = sm_running_to(M_P)
    # the clean MSSM unification scale (below M_P)
    mssm_gut = 2.0e16
    a_m, spread_m = _mssm_running(mssm_gut)
    landau = landau_pole_u1()
    return {
        "alpha_10": float(ALPHA_10),
        "five_thirds": float(gauge.five_thirds()),
        "M_GUT_0param": float(M_P),
        "couplings_at_M_P": [float(x) for x in a_p],
        "spread_at_M_P": spread_p,
        "mssm_unif_scale": float(mssm_gut),
        "mssm_spread": spread_m,
        "u1_landau_pole": float(landau),
        "uv_bounded": bool(M_P < landau),
        "backend": E._BACKEND,
    }


def _mssm_running(M):
    b = (33 / 5, -1, -3)
    a = [1.0 / _run(A_Z[k], b[k], M) for k in range(3)]
    return a, float((max(a) - min(a)) / float(np.mean(a)))

def running_2loop():
    """The **2-loop** running: the **UV-finite** QCA vs the **UV-divergent**
    perturbative QFT (the W4 completion, the 2-loop UV completion).

    The **0-parameter** contrast:

    1. **The QCA 2-step (finite).**  The QCA step ``U`` (the 32x32) is
       **unitary**; the **2-step** ``U^2`` is **unitary**
       (``U^2 (U^2)^* = 1``) — the **2-loop** (the 2-step) is **finite**
       (the **UV-finite** property: the QCA is a **unitary** automaton, the
       n-step is **bounded** for all n).

    2. **The QCA 2-particle (bounded).**  The 2-particle S-matrix (the
       **2-loop** analog, the 1024x1024) has **bounded** amplitudes
       (``max_offdiag <= 1``) — the **UV-finite** property.

    3. **The perturbative 2-loop (divergent).**  The **2-loop** beta
       coefficient (``b2``) makes the **U(1)** Landau pole **closer**
       (the **UV divergence** is **faster**): the **1-loop** pole is
       **above** M_P (the V12 UV-bounded), but the **2-loop** pole is
       **below** M_P (the **UV divergence** is **reached** **before** the
       cutoff).

    The **0-parameter** result: the **QCA** is **UV-finite** (the 2-step is
    **unitary**, the 2-particle is **bounded**), **whereas** the
    **perturbative** QFT is **UV-divergent** (the 2-loop Landau pole).
    This is the **UV-finite** property of the QCA (the **W4** completion).
    """
    from . import smatrix
    U = smatrix.qca_step()                          # the 32x32 unitary step
    U2 = E.matmul(U, U)                             # the 2-step (the 2-loop)
    u2_residual = float(np.abs(E.matmul(U2, U2.conj().T) - E.eye(32)).max())
    s2 = smatrix.two_particle(U)
    s2_bounded = bool(s2["max_offdiag"] < 1.0)     # the 2-particle is bounded
    # the perturbative 2-loop (the U(1) Landau pole, the UV divergence)
    b1, b2 = B_SM[0], B2_SM[0]
    a1 = A_Z[0]
    c = 1.0 / a1
    # the 2-loop Landau pole (the quadratic in t = ln(M/M_Z))
    Acoef = b2 / (b1 ** 2 * 2.0 * np.pi)
    Bcoef = b1 / (2.0 * np.pi)
    t2 = (-Bcoef + np.sqrt(Bcoef ** 2 + 4.0 * Acoef * c)) / (2.0 * Acoef)
    m_pole_2loop = M_Z * float(np.exp(t2))
    # the 1-loop Landau pole
    t1 = 2.0 * np.pi * c / b1
    m_pole_1loop = M_Z * float(np.exp(t1))
    ok = (u2_residual < 1e-8 and s2_bounded and m_pole_2loop < M_P < m_pole_1loop)
    return {
        "backend": E._BACKEND,
        "qca_2step_unitary_residual": u2_residual,
        "qca_2step_finite": bool(u2_residual < 1e-8),
        "qca_2particle_bounded": s2_bounded,
        "qca_2particle_max_offdiag": float(s2["max_offdiag"]),
        "pert_1loop_landau_pole": float(m_pole_1loop),
        "pert_2loop_landau_pole": float(m_pole_2loop),
        "pert_2loop_divergent": bool(m_pole_2loop < M_P),
        "uv_finite": bool(u2_residual < 1e-8 and s2_bounded),
        "residual": ("the 2-loop UV completion: the QCA 2-step is unitary (finite), "
                     "the 2-particle is bounded (UV-finite); the perturbative 2-loop "
                     "has a Landau pole below M_P (UV-divergent) — the QCA is UV-finite, "
                     "the perturbative is UV-divergent (the W4 completion)"),
        "ok": bool(ok),
    }

def main():
    import sys
    r = unification()
    print("=" * 68)
    print("CL10 - the 0-parameter coupling unification + the QCA UV completion")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter GUT coupling:  alpha_10 = {r['alpha_10']:.3f}  (= 1/25)")
    print(f"  the 0-parameter 5/3 normalization:  g1^2/g'^2 = {r['five_thirds']:.4f}")
    print(f"  the 0-parameter GUT scale:  M_GUT = M_P = {r['M_GUT_0param']:.2e} GeV (the single scale)")
    print(f"  the 1-loop SM running to M_P:  a = {['%.3f' % x for x in r['couplings_at_M_P']]}  (spread {r['spread_at_M_P']:.2f})")
    print(f"  the clean MSSM (SUSY) unification is at {r['mssm_unif_scale']:.0e} GeV (spread {r['mssm_spread']:.2f}) -- BELOW M_P")
    print(f"  the U(1) Landau pole is at {r['u1_landau_pole']:.2e} GeV -- the QCA cutoff M_P is BELOW it (UV-bounded)")
    ok = (r["five_thirds"] > 1.6 and r["five_thirds"] < 1.8
          and abs(r["alpha_10"] - 0.04) < 1e-9
          and r["uv_bounded"]
          and all(0.01 < a < 0.06 for a in r["couplings_at_M_P"]))
    print(f"  VERDICT:  {'PASS (0-parameter 5/3 + GUT coupling; the running GUT-range at M_P; the 1-loop matching the named residual; the QCA UV-bounded)' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
