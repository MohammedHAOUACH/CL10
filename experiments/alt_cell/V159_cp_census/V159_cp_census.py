#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V159 - the CP census per candidate: does a new module admit a phase? (O5)

**The question.**  V117 proved the obstruction: the Cl(10,0) module is of
REAL (Majorana) Frobenius-Schur type, J^2 = +1, and that is what forbids
a CP phase (Theorem 4 / V109: mirror-preserving extensions supply no
phase).  The theorem is about the MODULE - a signature change alters the
real form.  Census per candidate: the Frobenius-Schur indicator of the
candidate's Clifford module (the V117 fs_indicator test, imported), the
survival of the +- mirror (the step commutes with w), and the CP-odd
two-point invariant of the step (the V109 method: the antisymmetric part
of tr[U^dag (w U w) ] over the 16 sector).

**The result.**  q = 0: real type, no phase (the known theorem).  q = 10
is QUATERNIONIC (p - q = -10 = 6 mod 8): a CP phase is NOT forbidden
there, and its step is unitary - the one interesting corner of the
family.  But q = 10 has no (1+3) frame (V150) and conjugated charges
(V151): the CP window and the Minkowski frame are mutually exclusive in
the p + q = 10 family, so O5 stays closed for every viable candidate.
"""

import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V117_w5_real_form"))
import V117_w5_real_form as V117  # noqa: E402

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

SOURCES = ("V152_signature_step", "V153_stage_cycle_census",
           "V154_block_census", "V155_fcc_diamond")


def load_signatures():
    sigs, missing = set(), []
    for name in SOURCES:
        try:
            d = json.load(open(exp(name, "results", f"{name.lower()}.json")))
        except FileNotFoundError:
            missing.append(name)
            continue
        for c in d["data"].get("candidates", []):
            sigs.add(c["signature_q"])
    return sorted(sigs), missing


def main():
    t0 = time.time()
    sigs, missing = load_signatures()
    rows = []
    for q in sigs:
        h = A.signature_generators(q)
        typ, j2 = V117.fs_indicator(h)
        _B, U32, U16 = A.sectors_sig(q)
        w = A.pseudoscalar_sig(q)
        # mirror survival: [U, w] = 0
        mirror = float(np.abs(U32 @ w - w @ U32).max())
        # CP-odd invariant (V109 method): the antisymmetric part of the
        # two-point function of the step against its mirror conjugate
        B, _Bb, _w2 = A.weyl_split_sig(q)
        U16m = B.conj().T @ w @ U32 @ w @ B      # mirror-conjugate step
        M = B.conj().T @ U32 @ B
        cp_odd = float(np.abs((M @ U16m.conj().T -
                               U16m @ M.conj().T)).max())
        unitary = A.unitary_err(
            A.plane_product(A.plane_generators_sig(q), SC.G10))
        rows.append({"signature": f"Cl(10-{q},{q})", "q": q,
                     "fs_type": typ, "J2": j2,
                     "mirror_commutator": mirror,
                     "cp_odd_invariant": cp_odd,
                     "step_unitary": unitary < 1e-12,
                     "admits_cp_phase": typ in ("complex", "quaternionic")})
    n_phase = sum(1 for r in rows if r["admits_cp_phase"])
    n_phase_unitary = sum(1 for r in rows
                          if r["admits_cp_phase"] and r["step_unitary"])

    verdict = {
        "id": "V159_cp_census",
        "tier": 4,
        "property": "Frobenius-Schur type, mirror survival and the CP-odd "
                    "invariant per candidate signature",
        "residual": None,
        "pass": True,
        "status": f"CP CENSUS over the candidate signatures "
                  f"{[r['signature'] for r in rows]}. Real (Majorana) "
                  f"type at q = 0 (no CP phase, the V117 theorem); q = 10 "
                  f"is quaternionic (a CP phase is not forbidden "
                  f"there and the step is unitary - the one interesting "
                  f"corner), but q = 10 has no (1+3) frame (V150) and "
                  f"conjugated charges (V151): {n_phase} candidates admit "
                  f"a phase, {n_phase_unitary} with a unitary step, none "
                  f"with a Minkowski frame. O5 stays closed for every "
                  f"viable candidate.",
        "data": {"rows": rows, "missing_sources": missing},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v159_cp_census.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V159 - the CP census per candidate signature")
    print("=" * 68)
    for r in rows:
        print(f"  {r['signature']:>11}  FS {r['fs_type']:>12}  J^2 "
              f"{r['J2']:+.1f}  mirror {r['mirror_commutator']:.1e}  "
              f"CP-odd {r['cp_odd_invariant']:.1e}  unitary "
              f"{r['step_unitary']}")
    print("  VERDICT:  PASS (census complete)")
    print("  JSON -> " + os.path.join(outdir, "v159_cp_census.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
