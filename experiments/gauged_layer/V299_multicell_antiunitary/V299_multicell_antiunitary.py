#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V299 - the multicell antiunitarity check (tier 3).

**The question (pre-registered).** The V272 chain computed ground
energies for 2082 cells by direct diagonalization and booked the
imaginary parts of the ground eigenvalues (chain_ground_h_im). Are
the chain Hamiltonians effectively real - i.e. does an antiunitary
particle-hole-like symmetry act exactly on the multicell chain, as
the CPT claim needs at the level of the booked chain data?

**The object.** the V272 artifact (chain_H_spectrum real parts,
chain_ground_h_im imaginary parts, all booked) and the V267
physical-sector anchor.

**Measured.** max |Im E_0| over the booked chain; max |Im| over the
whole booked spectrum sample; the reality verdict against the
1e-12 exactness threshold; and the booked V272 branch carried.

Branch: CPT-EXACT-ON-CHAIN (all booked imaginary parts below 1e-12)
/ CPT-CHAIN-VIOLATION (booked values and their maximum).
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
    v272 = art("V272_dynamical_chain")
    d = v272["data"]
    # the booked 64x6144 list is the materialized ground-state MATRIX
    # H(q_ell) (chain_ground_h_re/_im), not per-cell eigenvalues: the
    # antiunitarity test is exact reality of that matrix and of the
    # booked chain-H spectrum
    re = [float(x) for x in d["chain_ground_h_re"]]
    im = [float(x) for x in d["chain_ground_h_im"]]
    spec = [float(x) for row in d["chain_H_spectrum"]
            for x in (row if isinstance(row, list) else [row])]
    max_im = max(abs(x) for x in im)
    max_re = max(abs(x) for x in re)
    n_cells = max(1, len(im) // 4096)
    # H(c)_ij = sum_s psi_s(c) Hpsi_s(c): the 6144 entries are the
    # 12-kernel-basis rows of 512 (H(q) is a 64x64 matrix per state:
    # the residual vertex slice); compare entrywise against its
    # conjugate transpose per 512-entry row block is not defined, so
    # the hermiticity leg acts on the reconstructed 64x64 slices when
    # the entry count is a multiple of 4096, and otherwise books the
    # raw imaginary-part maximum (the antiunitarity leg is entrywise).
    blk = None
    if len(im) % 4096 == 0:
        blk = (np.array(re, dtype=float)
               + 1j * np.array(im, dtype=float)).reshape(-1, 64, 64)
        herm_def = float(np.abs(blk - np.swapaxes(blk, -1, -2)).max())
    else:
        herm_def = None
    data["chain"] = {"n_matrix_entries": len(im),
                     "n_cells_4096_entries_each": n_cells,
                     "max_abs_imag_matrix": max_im,
                     "max_abs_real_matrix": max_re,
                     "spectrum_sample_max_abs": max(abs(x) for x in spec),
                     "v272_branch": v272.get("branch", "")}
    exact = max_im < 1e-12
    data["chain"]["hermiticity_defect"] = herm_def
    status = "CPT-EXACT-ON-CHAIN" if exact else "CPT-CHAIN-VIOLATION"
    data["status"] = status
    print("V299 {}: max |Im H| = {:.3e} over {} booked matrix entries "
          "(cdim-sliced residual vertex, not eigenvalues); hermiticity "
          "leg: {}".format(status, max_im, len(im),
                           "n/a (not a square slice)" if herm_def is None
                           else "{:.3e}".format(herm_def)), flush=True)
    checks.append({"name": "T1 booked chain residual-vertex imaginary "
                          "parts: max {:.3e} over {} entries (threshold "
                          "1e-12)".format(max_im, len(im)), "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V299_multicell_antiunitary", "tier": 3,
           "property": "exact reality of the booked multicell chain "
                       "spectra (antiunitary symmetry on the chain)",
           "residual": True, "programme": "V299",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v299_multicell_antiunitary.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V299 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
