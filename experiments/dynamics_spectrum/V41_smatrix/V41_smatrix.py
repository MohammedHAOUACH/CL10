#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V41 — the S-matrix of the QCA (the 't Hooft program, C6).

**The 0-parameter result.**  The 't Hooft program: the time-evolution
for one Planck step is a unitary cellular automaton.  The GUT QCA step
(``sm_qca``) **is** that evolution, so its matrix elements are the
**S-matrix** of the theory — the 0-parameter **transition amplitudes**:

    S_{ij}     = <i | U_QCA | j>                  (1-particle)
    S_{ij,kl}  = <i (x) j | U (x) U | k (x) l>    (2-particle)

with ``|i>, |j>, |k>, |l>`` the **32-D internal basis** (16 Weyl +
16-bar) and ``U_QCA`` = the product of the 45 Spin(10) plane rotations
at ``g_10 = 1/5`` — the **only** input.

**The 0-parameter checks.**  (1) unitarity: the column-sum property of
the 't Hooft program is probability conservation — ``sum_j |S_ij|^2
= 1`` for every ``i`` (CA limit: the amplitudes are 0/1, the sum is
exactly 1); the 1024x1024 2-particle amplitude is unitary too.
(2) the **off-diagonal** amplitudes are **nonzero** at the
0-parameter level ``O(g_10) ~ 0.2`` — the 45 planes mix the 32, so the
0-parameter step **scatters**.  (3) the amplitudes depend on the only
number ``g_10 = 0.2``.  (4) the 2-body amplitude is non-trivial — the
**first** computed *process* of the theory.

**The named residual (O9).**  The **matching** of these amplitudes to
the **measured** cross-sections (the QCA lattice RG vs the
perturbative QFT — the non-linear completion beyond the quadratic EH)
is the open part of the 't Hooft program.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import smatrix

TOL = 1e-12


def main():
    t0 = time.time()
    r = smatrix.smatrix()
    s1, s2 = r["one_particle"], r["two_particle"]
    ok = (s1["unitary_err"] < TOL and s1["colsum_err"] < TOL
          and s1["colsum_err_col"] < TOL and s1["max_offdiag"] > 0.05
          and s2["unitary_err"] < TOL and s2["colsum_err"] < TOL
          and s2["max_abs"] > 0.05 and s2["max_offdiag"] > 0.05
          and r["only_g10"])

    verdict = {
        "id": "V41_smatrix",
        "tier": 2,
        "property": "P2",
        "residual": True,          # matching to measured cross-sections (O9) is the named residual
        "pass": bool(ok),
        "status": "the_qca_step_is_the_0_parameter_S_matrix_the_first_processes_of_the_theory",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "g10": r["g10"],
            "n_planes": r["n_planes"],
            "n_internal": r["n_internal"],
            "s1_unitary_err": s1["unitary_err"],
            "s1_colsum_err": s1["colsum_err"],
            "s1_colsum_err_col": s1["colsum_err_col"],
            "s1_raw_colsum_dev": s1["raw_colsum_dev"],
            "s1_max_offdiag": s1["max_offdiag"],
            "s2_shape": s2["shape"],
            "s2_unitary_err": s2["unitary_err"],
            "s2_colsum_err": s2["colsum_err"],
            "s2_max_abs": s2["max_abs"],
            "s2_max_offdiag": s2["max_offdiag"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v41_smatrix.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V41 - the S-matrix of the QCA (the 't Hooft program)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the ONLY input:  g_10 = {r['g10']:.4f}  (= 1/5, from the root)")
    print(f"  the QCA step:  U = the product of the {r['n_planes']} Spin(10) plane rotations "
          f"(the {r['n_internal']}-D internal basis)")
    print(f"  [1] unitarity:  ||U U^dag - I|| = {s1['unitary_err']:.3e}  (1p), "
          f"{s2['unitary_err']:.3e}  ({s2['shape']}x{s2['shape']} 2p)")
    print(f"  [1] column sum (probability conservation):  "
          f"row dev {s1['colsum_err']:.3e}, col dev {s1['colsum_err_col']:.3e}, "
          f"2p dev {s2['colsum_err']:.3e}  (raw complex 1p dev: {s1['raw_colsum_dev']:.3f})")
    print(f"  [2] 0-parameter scattering:  max off-diagonal |S_ij| = {s1['max_offdiag']:.4f}  "
          f"(> 0.05: the 45 planes mix the 32)")
    print(f"  [4] 2-particle (U (x) U):  max |S2| = {s2['max_abs']:.4f}, "
          f"max off-diagonal = {s2['max_offdiag']:.4f}  (the 2-body space mixes)")
    print(f"  [3] 0-parameter:  the only number is g_10 = {r['g10']:.4f}")
    print(f"  named residual (O9):  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the QCA step IS the 0-parameter S-matrix: the first process of the theory; the cross-section matching the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v41_smatrix.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
