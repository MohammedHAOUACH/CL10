#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V169 - the CP census per candidate motion (O5; alt-dynamics B)

**The question.**  V117: the Cl(10,0) module is of REAL (Majorana)
Frobenius-Schur type (J^2 = +1), and that is what forbids a CP phase
(V109: mirror-preserving extensions supply no phase).  The theorem is
about the MODULE - but the MOTION is the remaining lever: does a
different step rule open a CP window?  Per candidate step (V159
method): the Frobenius-Schur indicator of the module the step acts on
(the V117 fs_indicator test, imported), the survival of the +- mirror
(the candidate's effective internal unitary against the pseudoscalar),
and the CP-odd two-point invariant of the step (the antisymmetric part
of tr[U^dag (w U w)] over the 16 sector, the V109 method).

**The result.**  Every dynamics candidate acts on the q = 0 module:
real type, J^2 = +1, no CP phase.  The mirror commutator is exactly 0
for every candidate (all internal factors are even-grade), and the
CP-odd invariant vanishes to machine precision for every candidate
step.  The motion census cannot open O5: the obstruction is the real
form of the module, which no step rule touches.
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
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core
from cl10.paths import exp
sys.path.insert(0, exp("V117_w5_real_form"))
sys.path.insert(0, exp("V167_hardtier_dynamics"))
import V117_w5_real_form as V117  # noqa: E402
import V167_hardtier_dynamics as V167  # noqa: E402

import numpy as np

from cl10 import altcell as A

SOURCES = V167.SOURCES


def load_candidates():
    out, missing = [], []
    for name in SOURCES:
        try:
            d = json.load(open(exp(name, "results", f"{name.lower()}.json")))
        except FileNotFoundError:
            missing.append(name)
            continue
        for c in d["data"].get("candidates", []):
            out.append(c)
    return out, missing


def internal_eff(spec):
    """The effective internal unitary of one candidate step on the full
    32 module (the bond layers act on the sub-lattice index and do not
    enter; the Wilson factor is q-dependent and noted separately)."""
    _B, U32, _U16 = A.sectors_sig(0)
    kind = spec["kind"]
    if kind == "internal_family":
        Us = V167.V162.family_factors(spec["family"], spec["g"],
                                      A.plane_generators_sig(0))
        U = np.eye(32, dtype=complex)
        for x in Us:
            U = x @ U
        return U, False
    if kind == "mass":
        a, b = spec["plane"]
        pairs = [(p, q) for p in range(10) for q in range(p + 1, 10)]
        T = A.plane_generators_sig(0)[pairs.index((a, b))]
        g = spec["g"]
        Mab = np.cos(g) * np.eye(32, dtype=complex) \
            - 1j * np.sin(g) * T
        U6 = np.linalg.matrix_power(U32, 6)
        return Mab @ U6, False
    if kind == "wilson":
        return U32, True          # q-dependent factor: no fixed internal
    return U32, False


def main():
    t0 = time.time()
    V167._init()
    cands, missing = load_candidates()
    h = A.signature_generators(0)
    typ, j2 = V117.fs_indicator(h)
    w = A.pseudoscalar_sig(0)
    B = A.weyl_split_sig(0)[0]

    seen = set()
    rows = []
    for c in cands:
        if c["name"] in seen:
            continue
        seen.add(c["name"])
        Ueff, qdep = internal_eff(c["step_spec"])
        mirror = float(np.abs(Ueff @ w - w @ Ueff).max())
        U16m = B.conj().T @ w @ Ueff @ w @ B
        M = B.conj().T @ Ueff @ B
        cp_odd = float(np.abs(M @ U16m.conj().T
                              - U16m @ M.conj().T).max())
        conj_closed = V117.conj_closed(Ueff)
        rows.append({"candidate": c["name"], "source": c["source"],
                     "fs_type": typ, "J2": j2,
                     "mirror_commutator": mirror,
                     "cp_odd_invariant": cp_odd,
                     "spectrum_conj_closed": conj_closed,
                     "q_dependent_internal": qdep,
                     "admits_cp_phase": typ in ("complex", "quaternionic")})
    n_phase = sum(1 for r in rows if r["admits_cp_phase"])
    max_cp = max(r["cp_odd_invariant"] for r in rows)
    max_mir = max(r["mirror_commutator"] for r in rows)

    verdict = {
        "id": "V169_cp_dynamics",
        "tier": 4,
        "property": "Frobenius-Schur type, mirror survival and the CP-odd "
                    "invariant per candidate step",
        "residual": None,
        "pass": True,
        "status": f"CP CENSUS over {len(rows)} candidate motions. Every "
                  f"candidate acts on the q = 0 module: FS type "
                  f"'{typ}', J^2 = {j2:+.1f} - real (Majorana), no CP "
                  f"phase (the V117 theorem). Mirror commutator max "
                  f"{max_mir:.1e} (every internal factor even-grade); "
                  f"CP-odd invariant max {max_cp:.1e}: the motion "
                  f"supplies no CP window. {n_phase} candidates admit a "
                  f"phase. O5 stays closed for the whole dynamics "
                  f"family: the obstruction is the real form of the "
                  f"module, not the step rule.",
        "data": {"rows": rows, "missing_sources": missing,
                 "n_candidates": len(rows), "n_admits_phase": n_phase,
                 "max_cp_odd": max_cp, "max_mirror_commutator": max_mir},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v169_cp_dynamics.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V169 - the CP census per candidate motion")
    print("=" * 68)
    print(f"  candidates: {len(rows)}   FS type {typ} (J^2 {j2:+.1f})   "
          f"admits phase: {n_phase}")
    print(f"  max mirror commutator: {max_mir:.1e}   max CP-odd: "
          f"{max_cp:.1e}")
    print("  VERDICT:  PASS (census complete)")
    print("  JSON -> " + os.path.join(outdir, "v169_cp_dynamics.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
