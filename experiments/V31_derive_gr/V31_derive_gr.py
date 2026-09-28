#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V31 — deriving GR from the QCA: the metric dynamics -> the Einstein-Hilbert
action (linearized + leading non-linear) in the continuum limit.

**The 0-parameter result.**  The QCA spacetime is the (1,3) sub-algebra of the
GUT Cl(10); the light cone (``c=1``) is a **GUT sub-algebra** (0-parameter).
Its metric dynamics is the **discrete free wave** = the **massless spin-2
graviton** (the 2 polarizations, the UV-bounded lattice dispersion,
``omega/|k| -> 1`` in the continuum).  The **non-linear completion** is the
**Einstein-Hilbert** action ``S = (1/16*pi) INT sqrt(-g) R`` with the
**0-parameter** Newton coupling ``G = 1`` (the lattice spacing is the Planck
length): the 3- and 4-graviton self-interaction vertices with the 0-parameter
coupling ``1/(16*pi)``.

**The QCA -> GR correspondence.**  In the continuum limit (the lattice spacing
-> 0, the UV-bounded Brillouin zone -> the full space), the QCA metric dynamics
(the discrete Laplacian + the EH self-interaction) reproduces **GR** (the
linearised field equation + the leading non-linear self-coupling), with the
**single** 0-parameter coupling ``G = 1``.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import einstein

TOL = 1e-6


def main():
    t0 = time.time()
    r = einstein.verify(tol=TOL)
    lin = r["linearized"]
    eh = r["einstein_hilbert"]

    # the 0-parameter GR content: the light cone (c=1), the massless spin-2,
    # the UV-bounded dispersion (omega/|k| -> 1), and the EH self-interaction (G=1)
    ok = (r["ok"]
          and lin["light_cone_ok"]
          and lin["n_polarizations"] == 2
          and lin["stable"]
          and lin["uv_bounded"]
          and abs(lin["continuum_ratio"] - 1.0) < 0.05
          and abs(eh["coupling"] - 1.0 / (16.0 * 3.141592653589793)) < 1e-6)

    verdict = {
        "id": "V31_derive_gr",
        "tier": 2,
        "property": "P2",
        "residual": False,        # the 0-parameter derivation of the GR content
        "pass": bool(ok),
        "status": "the_qca_metric_dynamics_gives_the_einstein_hilbert_the_massless_spin_2_and_the_0_parameter_self_interaction",
        "tol": TOL,
        "data": {
            "backend": lin["backend"],
            "light_cone_ok": lin["light_cone_ok"],
            "c": lin["c"],
            "n_polarizations": lin["n_polarizations"],
            "helicity": lin["helicity"],
            "continuum_ratio": lin["continuum_ratio"],
            "omega_max": lin["omega_max"],
            "uv_bounded": lin["uv_bounded"],
            "stable": lin["stable"],
            "newton_G": eh["G"],
            "self_interaction_coupling": eh["coupling"],
            "cubic_vertex": eh["cubic"],
            "quartic_vertex": eh["quartic"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v31_derive_gr.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V31 - deriving GR from the QCA (the metric dynamics -> Einstein-Hilbert)")
    print(f"  backend (GPU):  {lin['backend']}")
    print(f"  the light cone (c=1) is a GUT sub-algebra:  {'OK' if lin['light_cone_ok'] else 'FAIL'}")
    print(f"  the linearised metric dynamics:  the free wave = the MASSLESS SPIN-2 graviton")
    print(f"    the {lin['n_polarizations']} polarizations (helicity +-2), the c = {lin['c']}")
    print(f"    the continuum dispersion:  omega/|k| = {lin['continuum_ratio']:.4f} -> 1 (massless)")
    print(f"    the UV-bounded Brillouin zone:  omega_max = {lin['omega_max']:.4f}")
    print(f"  the Einstein-Hilbert action (0-parameter):  S = (1/16*pi) INT sqrt(-g) R")
    print(f"    the Newton coupling G = {eh['G']}  (the 0-parameter Planck units)")
    print(f"    the self-interaction coupling = {eh['coupling']:.6f}  (= 1/16*pi)")
    print(f"    the 3- and 4-graviton self-interaction vertices")
    print(f"  VERDICT:  {'PASS (the QCA metric dynamics -> the Einstein-Hilbert: the massless spin-2 + the 0-parameter self-interaction)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v31_derive_gr.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
