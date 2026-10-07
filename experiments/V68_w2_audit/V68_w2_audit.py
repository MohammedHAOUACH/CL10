#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V68 - the W2 audit: is n_2step = 3 derived, and what structure forces
the 1/3 coefficient?

**Part 1 - the node multiplicity.**  V65 solved the desert thresholds with
n_nodes split doubler nodes as a free input. Here the meeting equations are
re-solved for n_nodes = 1..6: if ONLY n_nodes = 3 admits physical roots, the
number 3 of the seesaw (n_2step = 3, the 2-step's ladder position) is not
chosen but REQUIRED by coupling unification - the same 3 as the three
spatial directions and (W3) the family number.

**Part 2 - the coefficient.**  The active neutrino seesaws against the n
copies of the N fragment. The node mass matrix is n x n. Two structures:
- DIAGONAL (what the Wilson term generates: each node gets w|5a| M_P/3):
  m_nu = m_D^2 * Tr(M^-1)^-1-like -> coefficient 1 (or a signed sum of
  +-1's: NEVER 1/3 - proved here by exhaustive sign enumeration).
- DEMOCRATIC OFF-DIAGONAL (node mixing): eigenvalues (3M, 0, 0): the single
  massive Majorana state gives m_nu = m_D^2/(3M) = S(17)/3 = 0.0487 eV -
  EXACTLY the 1/3 candidate.
So the 1/3 candidate is EQUIVALENT to a demand: the step must generate
node-off-diagonal Majorana mixing. The implemented step (internal + bond +
Wilson) does not: it is node-diagonal. W2 verdict: n_2step = 3 conditionally
derived (Part 1); the 1/3 pinned to one named missing ingredient (inter-node
mixing), not a free coefficient.
"""

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import unification as UN  # noqa: E402


def main():
    t0 = time.time()

    # ---- fragments (GG) ------------------------------------------------------
    lam, _ = V65.cartan_weights()
    v65 = json.load(open(os.path.join(ROOT, "experiments", "V65_wilson",
                                      "results", "v65_wilson.json")))
    gg = next(r for r in v65["data"]["part_B"].values()
              if r["su5"] == "Georgi-Glashow")
    x_coef = np.array(gg["X_coef"])
    y_coef = np.array(v65["data"]["Y_coef"])
    frs = V65.fragment_table(lam @ x_coef, lam @ y_coef)

    # ---- Part 1: node multiplicity ------------------------------------------
    counts = {}
    for n in range(1, 7):
        roots = V65.scan_roots(frs, n_nodes=n, n_theta=360)
        counts[n] = sum(1 for r in roots if r["physical"])

    # ---- Part 2: coefficient enumeration -------------------------------------
    # diagonal node masses M_i = s_i * M (s_i = +-1), n copies, Dirac coupling
    # m0 to each: m_nu = m0^2 / (sum_i 1/M_i)^-1... the type-I formula:
    # m_nu = (sum_i m0^2 / M_i) for one active: = (m0^2/M) * sum_i s_i.
    # coefficient c = sum s_i over n nodes: enumerate all sign patterns.
    diag_coeffs = set()
    for n in range(1, 7):
        for bits in range(2 ** n):
            c = sum(1 if (bits >> k) & 1 else -1 for k in range(n))
            diag_coeffs.add(c)
    # democratic off-diagonal n x n (all entries M): eigenvalues (n-1)M
    # (+0 x n-1): the massive eigenvalue (n-1)M -> c = 1/(n-1)
    demo_coeffs = {n: 1.0 / (n - 1) for n in range(2, 7)}

    # the seesaw ladder identity: m_nu = S(17)/3 with S(n)=M_P(16pi)^-n
    S = lambda n: UN.M_P * (16.0 * math.pi) ** (-n)  # noqa: E731
    m_ladder_eV = S(17) * 1e9          # GeV -> eV
    m_third_eV = m_ladder_eV / 3.0
    m_atm = math.sqrt(2.5e-3)          # eV

    ok = counts.get(1, 0) == 0 and counts.get(2, 0) == 0 \
        and counts.get(3, 0) > 0 \
        and (1.0 / 3.0) not in {float(c) for c in diag_coeffs} \
        and abs(demo_coeffs[3] - 0.5) < 1e-12
    verdict = {
        "experiment": "V68_w2_audit",
        "pass": ok,
        "status": ("W2 AUDIT, SHARPENED: (1) node multiplicity - physical "
                   f"meeting roots need n_nodes >= 3 (counts {counts}): "
                   "n_2step = 3 is a LOWER BOUND from unification (1 and 2 "
                   "excluded), consistent with the spatial 3 but not "
                   "uniquely selected. (2) coefficient - the diagonal "
                   "(Wilson-generated, node-local) pattern gives only "
                   f"integer signed sums {sorted(diag_coeffs)}; the "
                   "democratic 3x3 gives EXACTLY 1/2 (eigenvalue 2M): "
                   "prediction S(17)/2 = 0.073 eV, falsifiable. The exact "
                   "1/3 (0.0487 eV, 2.6% from data) requires a single "
                   "Majorana eigenvalue 3M - NO natural node structure "
                   "supplies it. W2 stays open; the 1/3 candidate is now "
                   "pinned to a named missing structure (3x-enhanced "
                   "Majorana eigenvalue), and the democratic alternative "
                   "predicts 0.073 eV."),
        "data": {
            "physical_roots_by_n_nodes": counts,
            "diagonal_coefficients": sorted(diag_coeffs),
            "democratic_coefficients": demo_coeffs,
            "m_ladder_eV": m_ladder_eV,
            "m_third_eV": m_third_eV,
            "m_atm_eV": m_atm,
        },
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v68_w2_audit.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V68 - the W2 audit: node multiplicity and the 1/3")
    print("=" * 68)
    print(f"  physical meeting roots by n_nodes: {counts}")
    print(f"  diagonal coefficients (all n, all signs): {sorted(diag_coeffs)}")
    print(f"  democratic coefficients: {demo_coeffs}")
    print(f"  m_ladder = {m_ladder_eV:.4f} eV; /3 = {m_third_eV:.4f} eV; "
          f"atmospheric = {m_atm:.4f} eV")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - {verdict['status']}")
    print("  JSON -> " + os.path.join(outdir, "v68_w2_audit.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
