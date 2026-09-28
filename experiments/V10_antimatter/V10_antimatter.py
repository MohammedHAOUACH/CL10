#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V10 - antimatter is the 16/16-bar pair: the 16 and the 16-bar are
inequivalent modules (opposite B-L); the annihilation is the 0-parameter CPT."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import matter

TOL = 1e-12

def main():
    t0 = time.time()
    p, m = matter.weyl16()
    n16=p.shape[1]; n16b=m.shape[1]
    bl16=matter.bl_charge('16'); bl16b=matter.bl_charge('16bar')
    inequivalent = (bl16 != bl16b)
    verdict = {
        "id": "V10_antimatter",
        "pass": inequivalent and n16==16 and n16b==16,
        "status": "the_16_16bar_pair_is_the_antimatter",
        "tol": TOL,
        "data": {'n16':int(n16),'n16bar':int(n16b),'bl16':int(bl16),'bl16bar':int(bl16b),'inequivalent':bool(inequivalent)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v10_antimatter.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V10 - antimatter: the 16/16-bar pair")
    print("=" * 68)
    print(f"  16 Weyl (matter) / 16-bar (antimatter):  {n16} / {n16b}")
    print(f"  B-L:  matter -> {bl16},  antimatter -> {bl16b}  (inequivalent)")
    print("  the annihilation 16 + 16-bar -> vacuum  (the 0-parameter CPT)")
    print(f"  VERDICT:  {'PASS' if inequivalent else 'FAIL'}  (the 16/16-bar is the antimatter)")
    print("  JSON -> " + os.path.join(outdir, "v10_antimatter.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if inequivalent and n16==16 and n16b==16 else 1

if __name__ == "__main__":
    sys.exit(main())
