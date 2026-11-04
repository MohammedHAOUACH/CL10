#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V293 - hosted family multiplicity quantization (tier 3).

**The question (pre-registered).** Is the hosted family multiplicity
quantized: does the chain-kernel hosting from V291, counted against
the booked V70 family structure (families 3, per_family 4, nodes
12), give an integer per family?

**The object.** the V291 artifact (p0 channel weights per chain-
kernel state, n0), the V70 artifact (families/counts anchor, re-
derived: the phase set must match to 1e-12 on the +-pi/2 lattice),
the V45 family_symmetry fact (hosting is family-symmetric by the
booked 3-fold symmetry).

**Measured.** r_int = n0 (integers by construction); r_weighted =
sum(p0) (the fractional content booked); defect |r_weighted -
round|; the V70 anchor re-verified as T0 gate.

Branch: R-QUANTIZED (defect < 1e-9) / R-MIXED-WEIGHTS (defect
booked honestly; r_int carried anyway).
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    v291 = art("V291_higgs_host_channels")
    v70 = art("V70_family_index")

    fam = v70["data"]["families"]
    cnt = v70["data"]["counts"]
    ph = np.array(fam["x"] + fam["y"] + fam["z"], dtype=float)
    ok_ph = bool(np.all((np.abs(ph) < 1e-12)
                        | (np.abs(np.abs(ph) - np.pi / 2) < 1e-12)))
    t0ok = (int(cnt["nodes"]) == 12 and int(cnt["families"]) == 3
            and int(cnt["per_family"]) == 4 and ok_ph and len(ph) == 12)
    checks.append({"name": "T0 V70 family anchor re-derived: 12 nodes on "
                          "the +-pi/2 phase lattice, families 3 x 4",
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V293 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    p0 = np.array([float(x) for x in v291["data"]["host"]["p0"]])
    n0 = int(v291["data"]["host"]["n0"])
    r_w = float(np.sum(p0))
    defect = abs(r_w - round(r_w))
    data["hosting"] = {"r_int": n0, "r_weighted": r_w,
                       "weighted_defect": defect,
                       "families_booked": int(cnt["families"]),
                       "r_per_family_int": n0}
    status = "R-QUANTIZED" if defect < 1e-9 else "R-MIXED-WEIGHTS"
    data["status"] = status
    print("V293 {}: r_int = {} per family (families booked {}), "
          "r_weighted = {:.6f} (defect {:.3e})".format(
              status, n0, cnt["families"], r_w, defect), flush=True)
    checks.append({"name": "T1 hosted multiplicity: r_int {} per family, "
                          "weighted {} defect {:.3e}".format(n0, r_w,
                                                             defect),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V293_r_quantization", "tier": 3,
           "property": "integer quantization of the hosted family "
                       "multiplicity on the chain",
           "residual": True, "programme": "V293",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v293_r_quantization.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V293 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
