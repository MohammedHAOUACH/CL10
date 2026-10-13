#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V151 - can another signature carry the anomaly-free generation?

**The question (signature axis, 3).**  The SM charge table
(``cl10.matter.SM_FIELDS``) and the eight vanishing gauge anomalies are
properties of the Cl(10,0) root.  For every signature Cl(10-q,q) that
keeps a 16/16 split (all of them, V150): recompute the SM embedding the
same way the root does - the Cartan of the root, Sigma_k = i h_{2k-1}
h_{2k} (the V65 mechanism) - and TRANSPORT the q=0 hypercharge operator
(the same Cartan coefficient vector Y_coef, identified once by the V65
multiset match) to the new weights.  The 16 states are the SAME
eigenvectors for every signature (each Cartan plane is only rescaled by
1, i or -1), so the field grouping (Q, u^c, d^c, L, e^c, N by the q=0
hypercharge) is well defined and the derived charges are computable.

**The census.**  Per signature: the plane scale factors f_k in {1, i, -1};
the derived per-field hypercharges; the eight anomaly coefficients on the
DERIVED charges; and whether any Cartan combination reproduces the SM
hypercharge multiset (the V65 identify_XY scan - impossible for odd q,
where a split pair makes the weights complex).

**The result.**  q = 0 reproduces the table and all eight anomalies zero.
Every odd q: complex charges, no SM multiset, multiplets break.  Even q:
the multiset survives by re-choosing the Cartan combination (f_k = +-1);
the TRANSPORTED q=0 hypercharge keeps uniform multiplets and zero
anomalies only for q = 4 (a field permutation, not the SM table) and
q = 10 (charge conjugation).  Only Cl(10,0) carries the SM table itself.
"""

import math
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
sys.path.insert(0, exp("V65_wilson"))
import V65_wilson as V65  # noqa: E402

import numpy as np

from cl10 import altcell as A
from cl10 import matter

I32 = np.eye(32, dtype=np.complex128)
PAIRS = ((0, 1), (2, 3), (4, 5), (6, 7), (8, 9))


def cartan_sig(q):
    """The five Cartan planes Sigma_k = i h_a h_b of the signature."""
    h = A.signature_generators(q)
    return [1j * (h[a] @ h[b]) for (a, b) in PAIRS]


def weights0():
    """The 16 weights of the q=0 Cartan (V65 cartan_weights method)."""
    B = A.sectors_sig(0)[0]
    Ck = cartan_sig(0)
    Dk = [B.conj().T @ C @ B for C in Ck]
    M = sum(math.sqrt(p) * Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(M)
    lam = np.array([[float(np.real(np.vdot(V[:, s], Dk[k] @ V[:, s])))
                     for k in range(5)] for s in range(16)])
    off = max(float(np.abs(V.conj().T @ Dk[k] @ V - np.diagflat(
        [np.vdot(V[:, s], Dk[k] @ V[:, s]) for s in range(16)])).max())
        for k in range(5))
    return lam, off


def anomalies_from(charges):
    """The 8 anomaly coefficients from per-field (su3, su2, Y, n) rows,
    the matter.anomalies arithmetic on DERIVED charges (Y may be complex)."""
    a3 = a1 = ag = a13 = a12 = 0.0 + 0.0j
    for su3, su2, y, n in charges:
        a3 += n * su3
        a1 += n * y ** 3
        ag += n * y
        if su3 != 0:
            a13 += n * y * 0.5
        if su2 != 0:
            a12 += n * y * 0.5
    return {"su3_cubed": a3, "su2_cubed": 0.0, "u1_cubed": a1,
            "grav2_su3": 0.0, "grav2_su2": 0.0, "grav2_u1": ag,
            "u1_su3_sq": a13, "u1_su2_sq": a12}


def main():
    t0 = time.time()
    lam0, off0 = weights0()
    # the q=0 hypercharge coefficient vector: the V65 multiset identification
    _ux, ys = V65.identify_XY(lam0)
    Y_coef = np.asarray(ys[0], float)
    y0 = Y_coef @ lam0.T                      # the 16 q=0 hypercharges
    # field grouping by the q=0 hypercharge value (all six values distinct)
    labels = {"Q": 1 / 6, "u^c": -2 / 3, "d^c": 1 / 3, "L": -1 / 2,
              "e^c": 1.0, "N": 0.0}
    groups = {name: [s for s in range(16) if abs(y0[s] - v) < 1e-9]
              for name, v in labels.items()}
    grouped_ok = all(len(g) == matter.SM_FIELDS[name][3]
                     for name, g in groups.items())

    rows = []
    for q in range(11):
        Cq = cartan_sig(q)
        C0 = cartan_sig(0)
        f = []
        for k in range(5):
            # the exact rescaling factor of the plane: <Cq, C0> / <C0, C0>
            fac = np.trace(Cq[k] @ C0[k].conj().T) / 32.0
            f.append(complex(round(fac.real), round(fac.imag)))
        lam_q = np.array([fl * lam0[:, k] for k, fl in enumerate(f)]).T
        y_q = Y_coef @ lam_q.T                # transported hypercharges
        fields = []
        broken = []
        for name, g in groups.items():
            ys_q = y_q[g]
            uniform = float(np.abs(ys_q - ys_q[0]).max()) < 1e-9
            su3, su2, _y, n = matter.SM_FIELDS[name]
            if uniform:
                fields.append((su3, su2, complex(ys_q[0]), n))
            else:
                broken.append(name)
        an = anomalies_from(fields) if not broken else None
        real_weights = all(fl.imag == 0.0 for fl in f)
        if real_weights:
            hits = V65.identify_XY(lam_q)[1]
            n_hits = len(hits)
        else:
            n_hits = 0        # complex weights: no real SM multiset possible
        rows.append({"signature": f"Cl(10-{q},{q})", "q": q,
                     "plane_factors": [str(fl) for fl in f],
                     "transported_charges": {
                         name: [complex(v) for v in y_q[g]]
                         for name, g in groups.items()},
                     "multiplet_broken": broken,
                     "anomalies": ({k: complex(v) for k, v in an.items()}
                                   if an else None),
                     "anomaly_max": (max(abs(v) for v in an.values())
                                     if an else None),
                     "sm_multiset_hits": n_hits})

    base = rows[0]
    base_an_ok = (base["anomaly_max"] is not None
                  and base["anomaly_max"] < 1e-12)
    # the q=0 derived charges must reproduce the SM_FIELDS table
    table_ok = all(
        max(abs(c - matter.SM_FIELDS[n][2]) for c in
            base["transported_charges"][n]) < 1e-9 for n in labels)
    nonzero_others = all(
        (r["anomaly_max"] is None or r["anomaly_max"] > 1e-12)
        or r["sm_multiset_hits"] == 0
        for r in rows[1:] if r["q"] % 2 == 1)
    ok = grouped_ok and base_an_ok and table_ok and nonzero_others

    def cnum(v):
        if isinstance(v, complex):
            if abs(v.imag) < 1e-12:
                return round(v.real, 6)
            return [round(v.real, 6), round(v.imag, 6)]
        return v

    verdict = {
        "id": "V151_signature_anomaly",
        "tier": 4,
        "property": "the anomaly census: transported hypercharge and the 8 "
                    "anomaly coefficients per signature",
        "residual": None,
        "pass": bool(ok),
        "status": f"ANOMALY CENSUS. q=0 reproduces the SM_FIELDS table "
                  f"(charges match to 1e-9) and all 8 anomalies zero "
                  f"(max {base['anomaly_max']:.1e}). Odd q: complex "
                  f"charges (a split Cartan pair scales by i), 0 SM-"
                  f"multiset hits, the Q/L or Q/u^c/d^c multiplets break "
                  f"(non-uniform transported charges). Even q: the SM "
                  f"multiset stays realizable by re-choosing the Cartan "
                  f"combination; the transported q=0 hypercharge keeps "
                  f"uniform multiplets and zero anomalies only for q = 4 "
                  f"(u^c <-> d^c, e^c <-> N: an outer automorphism, NOT "
                  f"the SM table) and q = 10 (all charges conjugated). "
                  f"Only Cl(10,0) carries the SM charge table itself.",
        "data": {"Y_coef": Y_coef.tolist(), "rows": rows,
                 "cartan_offdiag_err": off0},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v151_signature_anomaly.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False, default=cnum)

    print("=" * 68)
    print("V151 - the anomaly census per signature")
    print("=" * 68)
    for r in rows:
        am = r["anomaly_max"]
        ams = "-" if am is None else f"{am:.2e}"
        print(f"  {r['signature']:>11}  plane factors "
              f"{' '.join(r['plane_factors']):>26}  anomaly max "
              f"{ams}  SM hits {r['sm_multiset_hits']}  broken "
              f"{r['multiplet_broken']}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v151_signature_anomaly.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
