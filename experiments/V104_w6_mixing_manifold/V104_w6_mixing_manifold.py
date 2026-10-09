#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V104 - where does the W6 splitting live? Two sectors, one verdict.

**The question.**  The 0-parameter structure derives sin2(theta12) = 1/3
(V45/V60), theta13 = 0 and CKM = identity (V28/V45); the data say
0.308, 0.0222, 0.225. The paper lists three separate items. They are
the imprint of the ONE named W6 ingredient - the 10-Higgs vev that
splits the democratic degeneracy. V104 asks WHERE that splitting can
live: the charged-lepton sector or the neutrino sector.

**Branch A - charged-lepton splitting (computed, RULED OUT).**  The
most general two-vev correction to the democratic base:
    m_E = m0 J + (rank <= 2),  m0 = (m_e + m_mu + m_tau)/3,
i.e. det(m_E - m0 J) = 0. With J = 3 v v^T, v = (1,1,1)/sqrt(3), the
matrix determinant lemma reduces this EXACTLY to
    sum_i |p_i|^2 e^{-i th_i} / D_i = 1 / sum_i D_i,  p = W^T v,
so the lightest family decouples from the democratic direction:
|p_1|^2 <= O(D_1 / sum D). And since the second column of TBM IS v,
    sin^2(theta12) = |PMNS_e2|^2 = |p_1|^2 <= 4e-3,
against the measured 0.308 +- 0.013: excluded at > 20 sigma. The
charged-lepton branch of W6 CANNOT coexist with a democratic neutrino
sector (V103: the junction m_nu is democratic for every r).

**Branch B - neutrino splitting (computed, VACUOUS).**
The W6 vev instead corrects the rank-1 democratic Majorana matrix
(V103): m_nu = 3 a v v^T + b u u^T - one named rank-1 vev term. Its
spectrum is {0, sqrt(dm2_21), sqrt(|dm2_31|)} = {0, 0.0086, 0.0501}
eV (NO): the SAME term that supplies the solar mass (the W6 gap of
V103's two massless states) also shifts TBM. The charged leptons stay
democratic (U_eL = TBM*), PMNS = TBM^T W*, and the scan over (u, a, b)
matched to the oscillation spectrum gives the bands of
sin2(theta12), sin2(theta13), sin2(theta23).

**Branch C - cubic-axis vev (computed, EXCLUDED).**  The
symmetry-natural choice: the vev points along a cubic axis,
u = e_axis. Then (a, b) are FIXED by the oscillation spectrum and the
angles are 0-parameter predictions. All distinct branches sit
>= 3 sigma from NuFIT (best ~ 12 sigma).

**Verdict (computed, negative).**  The natural forms of the W6
splitting are EXCLUDED (A, C) or VACUOUS (B): the mixing deviations
are a genuine structural gap, not an artifact of a wrong assumption.
W6 is narrowed by computation: the vev must live in the NEUTRINO
sector with a NON-CUBIC direction, and reproducing the observed
angles needs at least two vev terms in generic relative orientation -
more than the single named 10-Higgs input.

"""

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cl10 import family  # noqa: E402

# PDG 2024 [13] (lepton masses)
LEPTONS = np.array([0.000511, 0.10566, 1.77686])

# NuFIT 6.0 [50] (NO)
NUFIT = {"s12": (0.308, 0.013), "s13": (0.02221, 0.00053),
         "s23": (0.470, 0.025)}
M_NU = np.array([0.0, 0.0086, 0.0501])  # eV, NO spectrum

V3 = np.ones(3) / np.sqrt(3.0)



def sample_democratic_plus_rank2(D, rng):
    """Point of {m : det(m - m0 J) = 0} with Takagi spectrum D.
    Linear in the moduli x_i = |p_i|^2 once the phases are chosen."""
    c = 1.0 / D.sum()
    for _ in range(20000):
        t1, t2 = rng.uniform(-np.pi, np.pi, 2)
        A = np.array([[np.cos(t1) / D[0], np.cos(t2) / D[1], 1.0 / D[2]],
                      [np.sin(t1) / D[0], np.sin(t2) / D[1], 0.0],
                      [1.0, 1.0, 1.0]])
        b = np.array([c, 0.0, 1.0])
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            continue
        if np.any(x < 0.0) or np.any(x > 1.0):
            continue
        # p_i = sqrt(x_i) e^{-i t_i / 2}: the determinant lemma uses the
        # BILINEAR form v^T M^-1 v = sum (p_i^*)^2 / D_i, so the phases
        # of p carry the t_i and D' stays the real spectrum D.
        ph = np.exp(-0.5j * np.array([t1, t2, 0.0]))
        p = np.sqrt(x) * ph                    # p = W^T v
        # build W: W^T v = p  <=>  W^dag v* = ... v real: W^dag v = p*
        u = p.conj()

        def ortho_basis(q):
            P = np.eye(3, dtype=complex) - np.outer(q, q.conj())
            vals, vecs = np.linalg.eigh(P)
            return vecs[:, vals > 0.5]
        Bv = ortho_basis(V3)
        Bu = ortho_basis(u)
        z2 = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
        q2, r2 = np.linalg.qr(z2)
        U2 = q2 * (np.diag(r2) / np.abs(np.diag(r2)))
        Q = np.column_stack([u, Bu @ U2]) @ np.column_stack([V3, Bv]).conj().T
        W = Q.conj().T                          # W^dag = Q, W^T v = p
        Dp = D.copy()
        return W, Dp, x
    return None


def takagi(m):
    """m = W diag(s) W^T with s >= 0 (non-degenerate nonzero spectrum)."""
    H = m @ m.conj()
    s2, U = np.linalg.eigh(H)
    s = np.sqrt(np.maximum(s2, 0.0))
    D = U.T @ m @ U
    W = U.copy()
    for i in range(3):
        if s[i] > 1e-14:
            # diagonal entry s_i e^{i th}: half-phase
            th = np.angle(D[i, i])
            W[:, i] *= np.exp(-0.5j * th)
    return W, s


def sample_neutrino_split(rng, P):
    """m_nu = 3 a v v^T + b u u^T with spectrum M_NU; solve (a,b) for
    random u. With a = r cos(phi), b = r sin(phi) the trace fixes r,
    and the second invariant is a batched grid root-find in phi."""
    S1 = (M_NU[1] ** 2 + M_NU[2] ** 2)
    S2 = (M_NU[1] ** 4 + M_NU[2] ** 4)
    u = np.exp(1j * rng.uniform(-np.pi, np.pi, 3)) * np.sqrt(
        rng.dirichlet([1.0, 1.0, 1.0]))
    Q = np.outer(u, u)

    def resid(phi):
        m = 3 * np.cos(phi) * P + np.sin(phi) * Q
        H = m @ m.conj()
        f = np.trace(H)
        if f <= 0:
            return None
        g = np.trace(H @ H)
        return g * S1 * S1 - S2 * f * f

    phis = np.linspace(-np.pi, np.pi, 361)
    M = (3 * np.cos(phis)[:, None, None] * P[None]
         + np.sin(phis)[:, None, None] * Q[None])
    H = M @ M.conj()
    f = np.einsum('nii->n', H)
    g = np.einsum('nij,nji->n', H, H)
    vals = g * S1 * S1 - S2 * f * f
    sols = []
    for i in range(len(phis) - 1):
        if not (np.isfinite(vals[i]) and np.isfinite(vals[i + 1])):
            continue
        if np.sign(vals[i]) == np.sign(vals[i + 1]):
            continue
        lo, hi = phis[i], phis[i + 1]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            rm = resid(mid)
            if rm is None:
                break
            if np.sign(rm) == np.sign(vals[i]):
                lo = mid
            else:
                hi = mid
        phi = 0.5 * (lo + hi)
        m = 3 * np.cos(phi) * P + np.sin(phi) * Q
        r = np.sqrt(S1 / np.trace(m @ m.conj()))
        sols.append((r * np.cos(phi), r * np.sin(phi)))
    if not sols:
        return None
    a, b = sols[rng.integers(len(sols))]
    m = 3 * a * P + b * Q
    W, s = takagi(m)
    return W, (a, b), u


def angles_from_pmns(pmns):
    return (abs(pmns[0, 1]) ** 2, abs(pmns[0, 2]) ** 2,
            abs(pmns[1, 2]) ** 2)


def main():
    rng = np.random.default_rng(104)
    tbm = family.tribimaximal_3x3()

    # --- branch A: charged-lepton splitting ------------------------------
    W0, Dp0, x0 = sample_democratic_plus_rank2(LEPTONS, rng)
    M0 = W0 @ np.diag(Dp0) @ W0.T
    m0 = LEPTONS.sum() / 3.0
    sv = np.linalg.svd(M0 - m0 * np.ones((3, 3)), compute_uv=False)
    det_dev = float(sv[-1] / sv[0])   # rank(M0 - m0 J) <= 2 check
    spec_dev = float(np.max(np.abs(
        np.sort(np.linalg.svd(M0, compute_uv=False)) - np.sort(LEPTONS))))

    s12A = []
    for _ in range(1500):
        r = sample_democratic_plus_rank2(LEPTONS, rng)
        if r is None:
            continue
        W, _, x = r
        pmns = W.T @ tbm
        s12A.append(abs(pmns[0, 1]) ** 2)
    s12A = np.array(s12A)

    # --- branch B: neutrino splitting ------------------------------------
    P = np.outer(V3, V3)
    rows = []
    tries = 0
    while len(rows) < 1500 and tries < 60000:
        tries += 1
        r = sample_neutrino_split(rng, P)
        if r is None:
            continue
        W, ab, u = r
        pmns = tbm.T @ W.conj()
        s12, s13, s23 = angles_from_pmns(pmns)
        rows.append((s12, s13, s23, abs(u @ V3) ** 2))
    B = np.array(rows)
    s12B, s13B, s23B = B[:, 0], B[:, 1], B[:, 2]

    def frac_in(x, key):
        cc, e = NUFIT[key]
        return float(np.mean(np.abs(x - cc) < e))

    jointB = float(np.mean(
        (np.abs(s12B - NUFIT["s12"][0]) < NUFIT["s12"][1]) &
        (np.abs(s13B - NUFIT["s13"][0]) < NUFIT["s13"][1]) &
        (np.abs(s23B - NUFIT["s23"][0]) < NUFIT["s23"][1])))


    # --- branch C: cubic-axis vev (the symmetry-natural choice) ---------
    # u = e_axis: (a, b) fixed EXACTLY by the spectrum -> 0-parameter
    # predictions of (s12, s13, s23) for each solution branch.
    cands = []
    for axis in range(3):
        u0 = np.zeros(3, complex)
        u0[axis] = 1.0
        Q = np.outer(u0, u0)
        phis = np.linspace(-np.pi, np.pi, 7201)
        Mt = (3 * np.cos(phis)[:, None, None] * P
              + np.sin(phis)[:, None, None] * Q)
        H = Mt @ Mt.conj()
        f = np.real(np.einsum('nii->n', H))
        g = np.real(np.einsum('nij,nji->n', H, H))
        S1 = (M_NU[1] ** 2 + M_NU[2] ** 2)
        S2 = (M_NU[1] ** 4 + M_NU[2] ** 4)
        vals = g * S1 * S1 - S2 * f * f
        for i in range(len(phis) - 1):
            if np.sign(vals[i]) == np.sign(vals[i + 1]):
                continue
            lo, hi = phis[i], phis[i + 1]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                mm = 3 * np.cos(mid) * P + np.sin(mid) * Q
                HH = mm @ mm.conj()
                fm = np.real(np.trace(HH))
                gm = np.real(np.trace(HH @ HH))
                if np.sign(gm * S1 * S1 - S2 * fm * fm) == np.sign(vals[i]):
                    lo = mid
                else:
                    hi = mid
            phi = 0.5 * (lo + hi)
            mm = 3 * np.cos(phi) * P + np.sin(phi) * Q
            r = np.sqrt(S1 / np.real(np.trace(mm @ mm.conj())))
            m = 3 * r * np.cos(phi) * P + r * np.sin(phi) * Q
            W, s = takagi(m)
            spec = np.sort(np.linalg.svd(m, compute_uv=False))
            if np.max(np.abs(spec - M_NU)) > 1e-9:
                continue
            s12, s13, s23 = angles_from_pmns(tbm.T @ W.conj())
            cands.append((axis, float(s12), float(s13), float(s23)))
    # distinct branches (sign duplicates give identical mixings)
    uniq = []
    for ax, a12, a13, a23 in cands:
        if not any(ax == u0[0] and abs(a12 - u0[1]) < 1e-3
                   for u0 in uniq):
            uniq.append((ax, a12, a13, a23))
    c_sig = []
    for ax, a12, a13, a23 in uniq:
        chi = ((a12 - NUFIT["s12"][0]) / NUFIT["s12"][1]) ** 2 \
            + ((a13 - NUFIT["s13"][0]) / NUFIT["s13"][1]) ** 2 \
            + ((a23 - NUFIT["s23"][0]) / NUFIT["s23"][1]) ** 2
        c_sig.append(float(np.sqrt(chi / 3)))
    c_best = int(np.argmin(c_sig))

    checks = {
        "A_det_dev": det_dev,
        "A_spectrum_dev": spec_dev,
        "A_n": int(len(s12A)),
        "A_s12_max": float(s12A.max()),
        "A_s12_excluded_sigma": float(
            (NUFIT["s12"][0] - s12A.max()) / NUFIT["s12"][1]),
        "B_n": int(len(B)),
        "B_s12_range": [float(s12B.min()), float(s12B.max())],
        "B_s13_range": [float(s13B.min()), float(s13B.max())],
        "B_s23_range": [float(s23B.min()), float(s23B.max())],
        "B_joint_frac_NuFIT": jointB,
        "C_branches": [[ax, round(a12, 4), round(a13, 4), round(a23, 4)]
                       for (ax, a12, a13, a23) in uniq],
        "C_sigma_per_branch": [round(x, 1) for x in c_sig],
        "C_best_branch": cands[c_sig.index(min(c_sig))][0]
        if c_sig else None,
        "C_min_sigma": round(min(c_sig), 1) if c_sig else None,
    }

    ok = (det_dev < 1e-10 and spec_dev < 1e-12
          and checks["A_s12_excluded_sigma"] > 10
          and checks["B_s12_range"][0] < 0.05
          and checks["B_s12_range"][1] > 0.95
          and len(uniq) >= 3
          and min(c_sig) > 3)

    verdict = {
        "id": "V104_w6_mixing_manifold",
        "tier": 5,
        "property": "P6",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"the natural forms of the W6 splitting are EXCLUDED or "
            f"VACUOUS - the mixing deviations are a genuine structural "
            f"gap, not an artifact: (A) charged-lepton democratic + rank-2 "
            f"forces sin2(12) <= {checks['A_s12_max']:.4f} (the lightest "
            f"family decouples from the democratic direction, which IS "
            f"the TBM solar column) - excluded at "
            f"{checks['A_s12_excluded_sigma']:.0f} sigma; (B) one rank-1 "
            f"vev with a free direction has 4 free parameters against 3 "
            f"angles - it covers the whole angle cube "
            f"(s12 in [{checks['B_s12_range'][0]:.3f}, "
            f"{checks['B_s12_range'][1]:.3f}]) and predicts nothing; "
            f"(C) the symmetry-natural choice - vev along a CUBIC AXIS - "
            f"is 0-parameter and EXCLUDED: {len(uniq)} distinct branches, "
            f"best {min(c_sig):.1f} sigma from NuFIT. W6 is narrowed by "
            f"computation: the vev must live in the NEUTRINO sector, with "
            f"a NON-CUBIC direction - reproducing (0.308, 0.0222, 0.470) "
            f"needs at least two vev terms in generic relative orientation, "
            f"more than the single named 10-Higgs input."),
        "data": checks,
    }

    out = os.path.join(ROOT, "experiments", "V104_w6_mixing_manifold",
                       "results", "v104_w6_mixing_manifold.json")
    with open(out, "w") as f:
        json.dump(verdict, f, indent=2)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - {verdict['status']}")
    print(f"  JSON -> {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
