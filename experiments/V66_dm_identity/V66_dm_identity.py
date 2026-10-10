#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V66 - the dark-matter identity scan (W7).

**The question.**  The ladder gives a DM MASS SCALE, m_DM = M_P (1/16pi)^10 =
118 GeV (V54), but no state of the implemented spectrum has been identified
as the carrier.  This experiment enumerates every candidate the formalism
actually supplies and tests the three DM properties - stability, neutrality,
mass - against the implemented step.

**The candidates.**
1. The neutrino (0.146 eV ladder): neutral, but hot and below the structure
   bound - excluded.
2. The desert vectorlike fragments (V65): colored or chiral - not DM.
3. The MIRROR doublers: the step preserves the 16 and the 16-bar separately
   (Theorem 4), so mirror number is an EXACT dynamical Z2 - the lightest
   16-bar doubler is EXACTLY stable.  The Wilson factor is internal, so the
   16-bar weights are the negatives of the 16 weights and the mirror doubler
   masses EQUAL the desert thresholds.  The lightest mirror state is the
   mirror Q-bar (colored, 1.3e8 GeV at the V65 root) - a stable colored
   relic, excluded; the neutral mirror N-bar sits at 2.0e10 GeV and is NOT
   the lightest.
4. The classical sector: the Planck-tension domain wall / Planck soliton -
   superheavy, overclosing.

**The verdict.**  No candidate is simultaneously stable, neutral and at
118 GeV.  The exact mirror Z2 is the only exact stability mechanism in the
theory, and with the V65 charge pattern it stabilizes a COLORED relic - the
mirror sector cannot be the DM as implemented.  The 118 GeV ladder state has
no identity: W7 stays open, with the tension quantified (118 GeV vs the
2.0e10 GeV spectral candidate, 8 orders).
"""

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import clifford
from cl10 import engine as E
from cl10 import supercell as SC
from cl10 import sm_qca
from cl10 import unification as UN

M_P = UN.M_P


def bar_blocks():
    """The 16-bar basis as the charge conjugate of the 16 basis: C = h_1
    anticommutes with w, so C maps the +i to the -i eigenspace, and state s
    of this basis is the conjugate of state s of the 16 basis (same order,
    opposite charges)."""
    B, U32, U16 = SC.sectors()
    C = E.asnumpy(clifford.generators()[0])
    Bb = C @ B
    return Bb, U32, Bb.conj().T @ U32 @ Bb


def bar_closure(n_stages=6, Lb=4):
    """16-bar closure of the IMPLEMENTED automaton step: leakage into the 16.
    The field lives on the block lattice, (Lb,Lb,Lb,8,32); the six stages of
    one step are applied one by one and the leakage is measured after each."""
    Bb, U32, _ = bar_blocks()
    B, _, _ = SC.sectors()
    Pb = np.kron(np.eye(SC.NSUB), Bb @ Bb.conj().T).reshape(
        SC.NSUB, 32, SC.NSUB, 32)
    rng = np.random.default_rng(1)
    v = rng.standard_normal(16) + 1j * rng.standard_normal(16)
    f = np.zeros((Lb, Lb, Lb, SC.NSUB, 32), dtype=complex)
    f[Lb // 2, Lb // 2, Lb // 2, 0] = Bb @ v
    leak = 0.0
    for even, j in SC.STAGES[:n_stages]:
        f = SC.internal_layer(f, U32)
        f = SC.bond_layer(f, even, j, SC.G10)
        comp = np.einsum("aibj,...bj->...ai", Pb, f)
        leak = max(leak, float(np.linalg.norm(f - comp) / np.linalg.norm(f)))
    return leak


def bar_weights_qr():
    """The 16-bar Cartan weights in the -i eigenspace basis (QR): the
    multiset must be the negatives of the 16 weights."""
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    w = E.asnumpy(clifford.pseudoscalar()).astype(np.complex128)
    ev, evec = np.linalg.eig(w)
    Bb = np.linalg.qr(evec[:, np.argsort(-np.imag(ev))[16:]])[0]
    Dk = [Bb.conj().T @ C @ Bb for C in cartan]
    M = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(M)
    lam = np.array([[float(np.real(np.vdot(V[:, s], Dk[k] @ V[:, s])))
                     for k in range(5)] for s in range(16)])
    return lam


def main():
    t0 = time.time()

    # ---- the exact mirror Z2 ------------------------------------------------
    leak = bar_closure()
    z2_exact = leak < 1e-12

    # ---- mirror weights = -weights ------------------------------------------
    # the 16 weights, computed directly
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    B, _, U16 = SC.sectors()
    Dk = [B.conj().T @ C @ B for C in cartan]
    Mm = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(Mm)
    lam16 = np.array([[float(np.real(np.vdot(V[:, s], Dk[k] @ V[:, s])))
                       for k in range(5)] for s in range(16)])
    lambar = bar_weights_qr()
    # the multiset of bar weights must equal the negatives of the 16 weights
    s16 = sorted(map(tuple, np.round(lam16, 9)))
    sbar = sorted(map(tuple, np.round(-lambar, 9)))
    conjugate = s16 == sbar

    # ---- mirror doubler masses = desert thresholds --------------------------
    # Wilson factor on the 16-bar one-step matrix with the V65 GG root charges
    v65 = json.load(open(os.path.join(ROOT, "experiments", "V65_wilson",
                                      "results", "v65_wilson.json")))
    gg = next(r for r in v65["data"]["part_B"].values()
              if r["su5"] == "Georgi-Glashow")
    root = next(r for r in gg["selective_roots"] if r["physical"])
    th, w = root["theta_rad"], root["w"]
    x_coef = np.array(gg["X_coef"])
    y_coef = np.array(v65["data"]["Y_coef"])
    # the physical 16-bar sector: the -i eigenspace (QR basis) with the
    # restricted Cartan weights lambar
    wmat = E.asnumpy(clifford.pseudoscalar()).astype(np.complex128)
    evq, evecq = np.linalg.eig(wmat)
    Bbq = np.linalg.qr(evecq[:, np.argsort(-np.imag(evq))[16:]])[0]
    U16b = Bbq.conj().T @ E.asnumpy(sm_qca.gut_qca_step(
        sm_qca.plane_generators(), SC.G10)).astype(np.complex128) @ Bbq
    # the 16-bar step matrix with CONJUGATE Wilson charges has the same phase
    # spectrum as the 16 one-step matrix (mirror Z2 + charge conjugation): the
    # mirror doubler masses EQUAL the desert thresholds
    def node_phases(U, ddiag):
        out = []
        for q in SC.M_POINTS:
            M = SC.bloch(q, U)
            ph = w * sum(1.0 - math.cos(qj) for qj in q)
            M = np.kron(np.eye(8, dtype=complex),
                        np.diag(np.exp(-1j * ph * ddiag))) @ M
            out.append(np.sort(np.angle(np.linalg.eigvals(M))))
        return np.concatenate(out)
    ph16 = node_phases(U16, np.cos(th) * (lam16 @ x_coef)
                       + np.sin(th) * (lam16 @ y_coef))
    phb = node_phases(U16b, np.cos(th) * (lambar @ x_coef)
                      + np.sin(th) * (lambar @ y_coef))
    spec_dev = float(np.abs(np.sort(ph16) - np.sort(-phb)).max())
    masses_equal = spec_dev < 1e-6  # relative to the node gap ~2e-4

    # ---- candidate table ----------------------------------------------------
    # mirror fragment masses from the V65 root (same thresholds)
    thr = root["thresholds_GeV"]
    # identify the lightest mirror fragment and its quantum numbers
    frs = gg["fragments"]
    Ms = {f["label"]: thr[f["label"]] for f in frs}
    lightest = min(frs, key=lambda f: Ms[f["label"]])
    def neutral(f):
        return f["Y"] == 0.0 and f["n"] == 1
    cand = [
        {"name": "nu (ladder)", "mass_GeV": 0.146, "neutral": True,
         "stable": "yes (lepton number)",
         "verdict": "EXCLUDED: hot, sum m_nu < 0.12 eV"},
        {"name": "desert fragments (V65)",
         "mass_GeV": min(Ms.values()), "neutral": False,
         "stable": "no (chiral, decays)",
         "verdict": "not DM: colored/chiral"},
        {"name": "mirror Q-bar doubler", "mass_GeV": Ms[lightest["label"]],
         "neutral": False, "stable": "EXACT (16-bar Z2)",
         "verdict": "EXCLUDED: stable COLORED relic"},
        {"name": "mirror N-bar doubler",
         "mass_GeV": Ms[next(k for k in Ms if k == "N")],
         "neutral": True, "stable": "EXACT (16-bar Z2)",
         "verdict": ("viable ONLY if heavier than every other mirror state "
                     "- NOT the case at the V65 root")},
        {"name": "domain wall / Planck soliton (classical)",
         "mass_GeV": M_P, "neutral": True, "stable": "classical",
         "verdict": "EXCLUDED: Planck tension, overcloses"},
    ]
    # thermal abundance of the 2e10 candidate: sigma v ~ m^2 / M_P^4
    mN = Ms["N"]
    sv = mN ** 2 / M_P ** 4
    omega_est = 0.12 * (3e-26 / sv)

    ok = z2_exact and conjugate and masses_equal
    verdict = {
        "experiment": "V66_dm_identity",
        "pass": ok,
        "status": ("W7 CANDIDATE SCAN COMPLETE, NO VIABLE CARRIER: the "
                   "mirror Z2 is exact (leak {leak:.1e}), the 16-bar weights "
                   "are the exact conjugates and the mirror doubler masses "
                   "equal the desert thresholds; the only exact-stability "
                   "candidate is the mirror sector, whose lightest state is "
                   "a COLORED Q-bar at {mQ:.1e} GeV (excluded) while the "
                   "neutral N-bar sits at {mN:.1e} GeV, 8 orders above the "
                   "118 GeV ladder scale, with naive thermal overclosure "
                   "x{ov:.0e}. The 118 GeV state has no identity in the "
                   "implemented spectrum.").format(
                       leak=leak, mQ=Ms[lightest["label"]], mN=mN,
                       ov=omega_est / 0.12),
        "data": {
            "mirror_Z2_leak": leak,
            "weights_conjugate": conjugate,
            "mirror_masses_equal": masses_equal,
            "candidates": cand,
            "sigma_v_est_GeV2": sv,
            "omega_overclosure_factor": omega_est / 0.12,
        },
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v66_dm_identity.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V66 - the dark-matter identity scan (W7)")
    print("=" * 68)
    print(f"  mirror 16-bar closure leak = {leak:.2e} -> exact Z2: {z2_exact}")
    print(f"  16-bar weights = -16 weights: {conjugate}")
    print(f"  mirror doubler spectrum = desert (dev {spec_dev:.2e}): "
          f"{masses_equal}")
    print("  candidates:")
    for c in cand:
        print(f"    {c['name']:42s} m = {c['mass_GeV']:.2e} GeV  "
              f"{c['verdict']}")
    print(f"  <sigma v> ~ m^2/M_P^4 = {sv:.1e} GeV^-2 "
          f"(need 3e-26): overclosure x{omega_est / 0.12:.0e}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - structural facts "
          f"verified; NO viable DM carrier at 118 GeV")
    print("  JSON -> " + os.path.join(outdir, "v66_dm_identity.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
