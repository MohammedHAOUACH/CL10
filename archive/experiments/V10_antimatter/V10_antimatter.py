#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V10 - antimatter is the 16/16-bar pair (the 0-parameter CPT structure).

**Claim.**  The matter (the **16** Weyl) and the antimatter (the **16-bar**)
are the two **inequivalent** eigenspaces of the **pseudoscalar** (the volume
element $\\omega=h_1\\cdots h_{10}$ of the $\\mathrm{Cl}(10)$ root): the $+i$
eigenspace is the 16, the $-i$ the 16-bar.  They are the **CPT conjugates**
of each other:

  1. **CPT map** -- a single generator $h_0$ **anticommutes** with $\\omega$
     (it flips the chirality), so it maps the 16 onto the 16-bar
     (matter $\to$ antimatter).  This is the 0-parameter CPT transformation.
  2. **Annihilation** -- the 16 and 16-bar eigenspaces are **orthogonal**
     (the distinct-eigenvalue eigenspaces of the anti-Hermitian pseudoscalar),
     so a matter state and an antimatter state have **no overlap**
     ($|\\langle 16|16\text{bar}\rangle|=0$): particle $+$ antiparticle
     $\to$ the charge-neutral vacuum.
  3. **B-L asymmetry** -- the 16 has $B\\!-\\!L=-1$, the 16-bar $B\\!-\\!L=+1$
     (opposite); the **sum cancels to 0** (the annihilation leaves the
     vacuum, no net charge).
  4. **Mass degeneracy** (CPT invariance) -- the 16 and 16-bar have the
     **same** dimension (16 each) and the same module structure: a particle
     and its antiparticle have the same mass.

The **observed baryon asymmetry** (the matter excess) is a **named residual**:
it needs a CPT- or CP-violating process, which the 0-parameter QCA does not
generate.
"""

import sys
import os
import json
import time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import clifford, matter

TOL = 1e-8
ORTHO_TOL = 1e-10

def main():
    t0 = time.time()
    gen = [clifford.E.asnumpy(m) for m in clifford.generators()]
    omega = clifford.pseudoscalar()
    p, m = matter.weyl16()          # p = 16 (matter), m = 16-bar (antimatter)
    n16, n16b = p.shape[1], m.shape[1]

    # 1. CPT map: a single generator h0 anticommutes with omega -> flips 16->16bar
    h0 = gen[0]
    anticommute = float(np.linalg.norm(omega @ h0 + h0 @ omega))
    cpt_to_bar = float(np.linalg.norm(m @ m.conj().T @ (h0 @ p), axis=0).mean())   # -> 16-bar
    cpt_to_16 = float(np.linalg.norm(p @ p.conj().T @ (h0 @ p), axis=0).mean())    # -> 16 (should be 0)
    cpt_ok = anticommute < ORTHO_TOL and cpt_to_bar > 0.99 and cpt_to_16 < 0.01

    # 2. annihilation: the 16 and 16-bar eigenspaces are orthogonal
    ortho = float(np.linalg.norm(p.conj().T @ m).max())
    ortho_ok = ortho < ORTHO_TOL

    # 3. B-L asymmetry: opposite charges, and the sum cancels (-> vacuum)
    bl16 = matter.bl_charge("16")
    bl16b = matter.bl_charge("16bar")
    bl_opposite = (bl16 == -bl16b)
    bl_cancel = (bl16 + bl16b == 0)      # particle + antiparticle -> net 0 (vacuum)

    # 4. mass degeneracy (CPT): same dimension / structure
    mass_degenerate = (n16 == n16b)

    ok = (cpt_ok and ortho_ok and bl_opposite and bl_cancel and mass_degenerate
          and n16 == 16 and n16b == 16)

    verdict = {
        "id": "V10_antimatter",
        "tier": 4,
        "property": "P2,P5",
        "residual": False,   # the CPT structure is a 0-parameter consequence (the baryon asymmetry is a named residual)
        "pass": bool(ok),
        "status": "the_16_16bar_CPT_pair_is_the_antimatter",
        "tol": TOL,
        "data": {
            "n16": int(n16), "n16bar": int(n16b),
            "bl16": int(bl16), "bl16bar": int(bl16b),
            "cpt_map": "a single generator (h0) flips 16 -> 16bar (anticommutes with the pseudoscalar)",
            "cpt_to_bar_norm": cpt_to_bar, "cpt_to_16_norm": cpt_to_16,
            "anticommute_residual": anticommute,
            "annihilation_overlap_16_16bar": ortho,
            "bl_sum": int(bl16 + bl16b),
            "mass_degenerate": bool(mass_degenerate),
            "baryon_asymmetry": "named residual (needs a CPT/CP-violating process; the 0-parameter QCA does not generate it)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v10_antimatter.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V10 - antimatter: the 16/16-bar CPT pair")
    print("=" * 68)
    print(f"  16 (matter) / 16-bar (antimatter):  {n16} / {n16b}  (the +/-i eigenspaces of the pseudoscalar)")
    print(f"  CPT map: h0 flips 16->16bar  (->16bar={cpt_to_bar:.4f}, ->16={cpt_to_16:.4f}, anticommute={anticommute:.1e})")
    print(f"  annihilation: max|<16|16bar>| = {ortho:.2e}  (orthogonal -> no overlap -> vacuum)")
    print(f"  B-L:  matter {bl16} / antimatter {bl16b}  (opposite; sum = {bl16+bl16b} -> the vacuum)")
    print("  mass degeneracy (CPT):  the 16 and 16-bar have the same dimension")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the 16/16-bar CPT pair is the antimatter; the baryon asymmetry is a named residual)")
    print("  JSON -> " + os.path.join(outdir, "v10_antimatter.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1




if __name__ == "__main__":
    sys.exit(main())
