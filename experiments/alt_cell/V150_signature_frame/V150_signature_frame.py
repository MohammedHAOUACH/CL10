#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V150 - the (1+3) frame per signature: does Cl(10,0) cause the frame choice?

**The question (signature axis, 2).**  The publication reads the (1+3)
Minkowski subalgebra inside Cl(10,0) as time = h_1 (+), space = i h_2,3,4
(``cl10.spacetime.causal_generators``).  For every signature Cl(10-q,q)
enumerate the 8 frame classes - time square -/+ x 0..3 space generators
with square -I - and decide per class:
1. is the class REALIZABLE by an index choice in the signature?
2. ``w^2 = -prod(eps_a)`` (the 16/16-bar split: +-i for even q, +-1 for
   odd q; the eigenspace dimensions must be 16/16 for every signature);
3. the Minkowski-subalgebra residual (``altcell.frame_residual``): the
   rotations J_a = S_b S_c and boosts K_a = T S_a of the frame must close
   as so(1,3); the residual is 0 iff the frame is exactly (1+3).

**The result.**  The 16/16 split survives EVERY signature (w^2 =
-(-1)^q, dimensions 16/16 always).  The exact (1+3) frame (time = h_t
raw, space = i h_s, the ``causal_generators`` convention) exists iff the
signature has at least four generators of square +I - i.e. q <= 6: the
frame indices can be re-chosen, so the frame alone does NOT force q = 0.
For q >= 7 the best frame residual is 4 (no generator of the required
square).  Combined with V149 (only q = 0 and q = 10 keep the 45-plane
product unitary, and q = 10 has no (1+3) frame), the Euclidean signature
is the only one carrying BOTH the unitary step and the Minkowski frame.
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

import numpy as np

from cl10 import altcell as A


def frame_classes(q):
    """The 8 frame classes (time_sign, k_minus) and their realizability in
    Cl(10-q,q): time_sign = the square of the time generator, k_minus = the
    number of space generators with square -I."""
    eps = A.signature_squares(q)
    out = []
    for t_sign in (+1.0, -1.0):
        for k in range(4):
            # pick a time generator with square t_sign, then 3 space
            # generators with exactly k of square -1 among the remaining 9
            t_idx = next((i for i, e in enumerate(eps) if e == t_sign), None)
            rest = [e for i, e in enumerate(eps) if i != (t_idx if t_idx is not None else -1)]
            m_rest = rest.count(-1.0)
            realizable = t_idx is not None and k <= m_rest and (3 - k) <= (9 - m_rest)
            rep = None
            if realizable:
                space = ([i for i, e in enumerate(eps) if i != t_idx and e == -1.0][:k]
                         + [i for i, e in enumerate(eps) if i != t_idx and e == +1.0][:3 - k])
                rep = (t_idx, tuple(space))
            out.append({"time_sign": int(t_sign), "space_minus": k,
                        "realizable": bool(realizable),
                        "frame": rep})
    return out


def main():
    t0 = time.time()
    rows = []
    for q in range(11):
        B, Bb, w2 = A.weyl_split_sig(q)
        # verify the split numerically: w B = +- B (eigenvalue check)
        w = A.pseudoscalar_sig(q)
        ev_plus = float(np.abs(w @ B - B @ (B.conj().T @ w @ B)).max())
        dims = (B.shape[1], Bb.shape[1])
        row = {"signature": f"Cl(10-{q},{q})", "q": q,
               "w2": round(w2, 9), "weyl_dims": list(dims),
               "split_err": ev_plus,
               "canonical_residual": A.frame_residual(q, 0, (1, 2, 3)),
               "frame_classes": frame_classes(q)}
        rows.append(row)

    base = rows[0]
    ok_base = (base["w2"] == -1.0 and base["weyl_dims"] == [16, 16]
               and base["canonical_residual"] < 1e-10
               and base["split_err"] < 1e-10)
    all_split = all(r["weyl_dims"] == [16, 16] and r["split_err"] < 1e-10
                    for r in rows)
    # the exact (1+3) frame exists iff some realizable class has residual 0
    mink_ok = []
    for r in rows:
        q = r["q"]
        best = min((A.frame_residual(q, t, s)
                    for fc in r["frame_classes"] if fc["realizable"]
                    for (t, s) in [fc["frame"]]), default=None)
        r["best_frame_residual"] = best
        if best is not None and best < 1e-10:
            mink_ok.append(r["signature"])

    ok = (ok_base and all_split
          and mink_ok == [f"Cl(10-{q},{q})" for q in range(7)])

    verdict = {
        "id": "V150_signature_frame",
        "tier": 4,
        "property": "the (1+3) frame census: w^2, the 16/16 split and the "
                    "Minkowski-subalgebra residual per signature x frame class",
        "residual": None,
        "pass": bool(ok),
        "status": f"FRAME CENSUS. q=0 baseline reproduces the implemented "
                  f"selection (time h_1, space i h_2,3,4: residual "
                  f"{base['canonical_residual']:.1e}, w^2 = -1, 16/16). The "
                  f"16/16 split survives every signature (w^2 = -(-1)^q). "
                  f"Exact (1+3) frame realizable for: {mink_ok} (needs at "
                  f"least four +I generators); q >= 7 best residual 4. The "
                  f"frame axis alone does not force q = 0; with V149 "
                  f"(unitary product only q = 0, 10) the Euclidean "
                  f"signature is the only one with both.",
        "data": {"rows": rows, "minkowski_signatures": mink_ok},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v150_signature_frame.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V150 - the (1+3) frame per signature")
    print("=" * 68)
    for r in rows:
        nreal = sum(1 for fc in r["frame_classes"] if fc["realizable"])
        print(f"  {r['signature']:>11}  w^2 {r['w2']:+.0f}  dims "
              f"{r['weyl_dims']}  canonical residual "
              f"{r['canonical_residual']:.1f}  best {r['best_frame_residual']:.1f}"
              f"  realizable classes {nreal}/8")
    print(f"  exact (1+3) frame: {mink_ok}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v150_signature_frame.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
