#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V63 - the Brillouin-zone chirality test: is the 16 a **single** chiral
sector, or does the Nielsen-Ninomiya doubling survive the QCA dynamics (W9)?

**What had never been done.**  The publication claimed "no doublers" as an
argument (the module is on-site, so the Nielsen-Ninomiya theorem does not
apply *as stated*) and named the Brillouin-zone analysis as the open item.
This experiment **runs that analysis**: the exact 6-step Bloch Floquet of the
implemented local step (``cl10.local_qca``), validated against the real-space
Floquet on a 4^3 torus, scanned over the Brillouin zone.

**The four results.**
1. **Chirality is a dynamical symmetry.**  The bond layer acts on the *cell*
   index, so it commutes with the pseudoscalar: a state initialized in the 16
   stays in the 16 under the implemented step (leakage ~1e-15).
2. **The zone center is gapped.**  min|omega(q=0)| = 0.0365 per 6-step cycle:
   the would-be massless matter is NOT at q = 0.
3. **The gapless points are the doublers.**  12 gapless nodes at
   (±pi/2, ±pi/2, 0) and permutations - the classic doubler momenta - with
   linear dispersion |v| ~ 0.14.
4. **Nielsen-Ninomiya is satisfied**: the 12 node chiralities sum to 0.  The
   doubling is present and the theory has no mechanism (Wilson term, overlap,
   domain wall) that removes it.

**The named residual (W9).**  Removing the doublers - or showing that the
gapped zone-center sector, not the gapless doubler momenta, is the physical
fermion sector.
"""

import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import brillouin

TOL = 1e-9


def main():
    t0 = time.time()
    r = brillouin.verify()
    ok = r["ok"]

    verdict = {
        "id": "V63_bz_chirality",
        "tier": 2,
        "property": "the 16-sector Floquet spectrum over the Brillouin zone",
        "residual": True,          # the doublers are present: W9 is NOT closed
        "pass": bool(ok),
        "status": "chirality_is_a_dynamical_symmetry_but_the_doublers_are_present_"
                  "12_gapless_nodes_at_the_zone_half_momenta_chirality_sum_0",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "g10": r["g10"],
            "floquet_cycle": r["floquet_cycle"],
            "sector_leakage_16_to_16bar": r["sector_leakage_16_to_16bar"],
            "bloch_vs_torus_dev": r["bloch_vs_torus_dev"],
            "gap_q0_16": r["gap_q0_16"],
            "gap_q0_32": r["gap_q0_32"],
            "n_doubler_nodes": r["n_doubler_nodes"],
            "max_doubler_gap": r["max_doubler_gap"],
            "v_typical": r["v_typical"],
            "chirality_sum": r["chirality_sum"],
            "nodes": r["nodes"],
            "g10_scan": r["g10_scan"],
            "residual": "the doublers are not removed (no Wilson/overlap/domain-wall "
                        "term in the theory); the gap at q=0 is not protected in g_10",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v63_bz_chirality.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V63 - the Brillouin-zone chirality test (the W9 doubling check)")
    print("=" * 68)
    print(f"  backend:  {r['backend']}   g_10 = {r['g10']:.4f}   Floquet cycle = {r['floquet_cycle']} steps")
    print(f"  [1] 16-sector closure of the IMPLEMENTED step:  16->16bar leakage = "
          f"{r['sector_leakage_16_to_16bar']:.2e}  ({'chirality is a dynamical symmetry' if r['sector_leakage_16_to_16bar'] < TOL else 'BROKEN'})")
    print(f"  [2] the Bloch Floquet vs the real-space 4^3 torus:  max deviation = "
          f"{r['bloch_vs_torus_dev']:.2e}  ({'the construction is exact' if r['bloch_vs_torus_dev'] < 1e-9 else 'MISMATCH'})")
    print(f"  [3] the zone center:  min|omega(q=0)| = {r['gap_q0_16']:.5f} (16), "
          f"{r['gap_q0_32']:.5f} (32)  ->  {'GAPPED: the matter is not at q=0' if r['gap_q0_16'] > 1e-3 else 'gapless'}")
    print(f"  [4] the doubler momenta (±pi/2, ±pi/2, 0):  {r['n_doubler_nodes']} nodes, "
          f"max gap = {r['max_doubler_gap']:.2e}, linear velocity |v| = {r['v_typical']:.4f}  "
          f"({'GAPLESS doublers' if r['max_doubler_gap'] < 1e-3 else 'gapped'})")
    print(f"  [5] Nielsen-Ninomiya:  sum of the node chiralities = {r['chirality_sum']:.1f}  "
          f"({'left and right nodes cancel: the doubling survives' if abs(r['chirality_sum']) < 1e-9 else 'unbalanced'})")
    print("  [6] g_10 sensitivity of the gaps:  " +
          ", ".join(f"g10={s['g10']}: gap(q=0)={s['gap_q0']:.4f}/gap(doubler)={s['gap_doubler']:.4f}"
                    for s in r["g10_scan"]))
    print(f"  named residual (W9):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the test is exact and the four results hold; W9 remains OPEN: the doublers are present)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v63_bz_chirality.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
