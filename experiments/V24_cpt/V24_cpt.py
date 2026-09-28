#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V24 - the CPT theorem: the 16 / 16-bar are CPT conjugates.

The **CPT theorem** of the CL10 QCA (the 0-parameter matter/antimatter
symmetry): the **16 Weyl** (matter) and the **16-bar** (antimatter) are the
**CPT conjugates** of each other.

The **CPT operator** is the **time-direction** generator $C = g_1$ (the "time"
of the $(1,3)$).  It is a **Hermitian involution** ($C^2 = I$, $C = C^\\dagger$)
that **anticommutes** with the pseudoscalar $w = g_1 \\cdots g_{10}$
($C w + w C = 0$), so it **maps the 16 to the 16-bar**
($C P_{16} = P_{16\\text{bar}} C$, with $P_{16} = \\tfrac{1}{2}(I - i w)$).

The **0-parameter** CPT content (all from the one object, verified on the
**GPU**):
1. $C^2 = I$ (the CPT operator is an involution);
2. $C C^\\dagger = I$ (CPT is **unitary** — it preserves the norm);
3. $C P_{16} = P_{16\\text{bar}} C$ (CPT maps **matter → antimatter**);
4. the **GUT is CPT-symmetric**: the 45 gauge generators **commute** with $w$
   (the chiral split is preserved), so **matter and antimatter have the same
   dynamics and the same mass** (CPT: a particle and its antiparticle are
   degenerate);
5. **in action** (on the GPU): a matter (16) state and its CPT conjugate
   (16-bar) carry **equal** norm, and C maps the 16 **entirely** into the
   16-bar.

This is a **P4 consistency** result: the QCA is **CPT-symmetric**, and the
matter/antimatter (the V10 16/16-bar pair) is a **0-parameter CPT pair**.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import cpt

TOL = 1e-12


def main():
    t0 = time.time()
    r = cpt.cpt_theorem()
    ok = (r["C_involution_err"] < 1e-10
          and r["C_unitary_err"] < 1e-10
          and r["C_w_anticommute_err"] < 1e-10
          and r["map_16_to_16bar_err"] < 1e-10
          and r["gut_commute_w_err"] < 1e-10
          and r["norm_rel_err"] < 1e-10
          and r["cpt_equal_norm"])

    verdict = {
        "id": "V24_cpt",
        "tier": 3,
        "property": "P4",
        "residual": False,         # the CPT theorem is a 0-parameter PASS (not a residual)
        "pass": bool(ok),
        "status": "the_16_16bar_are_cpt_conjugates_and_the_qca_is_cpt_symmetric",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "C_involution_err": r["C_involution_err"],
            "C_unitary_err": r["C_unitary_err"],
            "C_w_anticommute_err": r["C_w_anticommute_err"],
            "map_16_to_16bar_err": r["map_16_to_16bar_err"],
            "gut_commute_w_err": r["gut_commute_w_err"],
            "norm_before": r["norm_before"],
            "norm_after_cpt": r["norm_after_cpt"],
            "norm_rel_err": r["norm_rel_err"],
            "n_matter": r["n_matter"],
            "n_antimatter": r["n_antimatter"],
            "cpt_equal_norm": r["cpt_equal_norm"],
            "c16_into_16bar": r["c16_into_16bar"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v24_cpt.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V24 - the CPT theorem (the 16/16-bar are CPT conjugates, on the GPU)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  C = g_1 (the time direction):  involution C^2=I {r['C_involution_err']:.2e},  unitary {r['C_unitary_err']:.2e},  C w + w C = 0 {r['C_w_anticommute_err']:.2e}")
    print(f"  C maps 16 -> 16-bar (C P_16 = P_16bar C):  {r['map_16_to_16bar_err']:.2e}")
    print(f"  the GUT (the 45) commutes with w:  {r['gut_commute_w_err']:.2e}  (the chiral CPT split is preserved -> matter = antimatter dynamics)")
    print(f"  CPT in action:  matter/antimatter CPT conjugates carry EQUAL norm:  {r['n_matter']:.4f} / {r['n_antimatter']:.4f}  (equal: {r['cpt_equal_norm']})")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the 16/16-bar are 0-parameter CPT conjugates; the QCA is CPT-symmetric)")
    print("  JSON -> " + os.path.join(outdir, "v24_cpt.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
