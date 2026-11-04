#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V288 - the Bekenstein truncation map (tier 2, hard, bookkeeping).

**The question (pre-registered).** Which flux truncation of the link
module carries the Bekenstein factor toward 1/4, or none - the
V99/V100 ln 45 vs 1/4 gap named as a number: what does the
truncation map of the three-site kernel do to it?

**The object (closed input set).** V101's kernel builder and ground
state machinery (imported by spec: build_kernel, the cached
results/kernel.npz with the glued Is and the channel-ordered
Ts = [r = 1 (2), r = 45 (2), r = 210 (2)], the Gram/K construction
and the generalized bond-current eigenproblem, the middle-cut
entropy pattern), the V99 Haar vertex-invariant region numbers
(ln 45, bekenstein_factor read live), V100's capacity S_capacity and
vertex_coefficient (read live). Truncations: rebuild the three-site
kernel with the link flux space restricted to the channel subsets
{1}, {1,45}, {1,45,210}, {45}, {210} - each subset slices the booked
Ts block (the channels are the exact fusion sectors of the builder,
so slicing re-runs the construction on the reduced channel space
exactly: states (a,t,b) with t in the subset, exact Gram, exact
bond-current H, generalized eigenproblem, ground-state middle-cut
entropy per the V101 pattern).

The factor booked per row: S_bond / (1/4) - the same units as the
V99 factor 4 ln 45 and the V100 capacity factor; the V100
vertex-corrected coefficient (1/3) ln I_6 is booked as the fixed
vertex reference next to each row (no per-subset I_6 is claimed).

T0: ln 45 = 3.8066624897703183 region from the V99 artifact and the
V101 ground state reproduced live (E0, degeneracy, gap, S_ground to
1e-6 relative); S_cap = 4.4998 and the factor numbers read live from
the V100/V99 artifacts at run time.

Acceptance: BEKENSTEIN-MAP booked (the map rows); BEKENSTEIN-CLOSED
only if some truncation lands within 2 % of 1/4 (S_bond in
[0.98, 1.02] x 0.25 in the booked nats units). No truncation is
chosen in advance; every row is booked.
"""
import importlib.util
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402
from scipy.linalg import eigh as geigh                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v101", os.path.join(exp("V101_bond_current_ground_state"),
                         "V101_bond_current_ground_state.py"))
V101 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V101)
V79 = V101.V79


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def kernel_for_subset(Is, Ts, t_idx):
    """Exact slice of the glued kernel to a channel subset."""
    return Is, Ts[t_idx]


def eigen_and_entropy(Is, Ts, Ls_np):
    """The V101 pattern: Gram N, bond-current H, generalized eigen-
    problem, max-mixed ground-state middle-cut entropy."""
    nI, nT = Is.shape[0], Ts.shape[0]
    states = [(a, t, b) for a in range(nI) for t in range(nT)
              for b in range(nI)]
    ns = len(states)
    P = np.einsum("xpa,ypa->xya", Is.conj(), Is)
    Q = np.einsum("tjab,ujab->tuab", Ts.conj(), Ts)
    K = np.einsum("xya,tuab->xytub", P, Q)
    Z = np.einsum("xki,yki->xyi", Is.conj(), Is)
    N = np.einsum("xytub,zwb->xtzyuw", K, Z).reshape(ns, ns)
    N = (N + N.conj().T) / 2
    evn = np.linalg.eigvalsh(N)
    rank = int(np.sum(evn > 1e-9 * evn.max()))

    Ts_np = V79.root_generators()
    I32 = np.eye(32, dtype=np.complex128)
    Js = [np.kron(Ts_np[A], I32) + np.kron(I32, Ts_np[A])
          for A in range(45)]
    UJ = np.zeros((45, nI, nI, 45, 45), dtype=np.complex128)
    for A in range(45):
        for y in range(nI):
            JI = Js[A] @ Is[y]
            for x in range(nI):
                UJ[A, x, y] = Is[x].conj().T @ JI
    V2 = np.einsum("tjab,ujcb->tuac", Ts.conj(), Ts)
    O2 = np.einsum("tjab,ujac->tubc", Ts.conj(), Ts)
    H = np.zeros((ns, ns), dtype=np.complex128)
    for s, (a, t, b) in enumerate(states):
        for s2, (a2, t2, b2) in enumerate(states):
            h = 0.0 + 0.0j
            if b == b2:
                for A in range(45):
                    h += np.einsum("xz,xy,yz->", UJ[A, a, a2],
                                   Ls_np[A], V2[t, t2])
            if a == a2:
                for A in range(45):
                    h += np.einsum("xz,xy,yz->", UJ[A, b, b2],
                                   Ls_np[A], O2[t, t2])
            H[s, s2] = h
    H = (H + H.conj().T) / 2
    ev, evec = geigh(H, N)
    ev = np.real(ev)
    E0 = float(ev[0])
    deg = int(np.sum(np.abs(ev - E0) < 1e-8 * max(1.0, abs(ev[-1]))))
    gap = float(ev[deg] - E0) if ns > deg else 0.0

    c = evec[:, :deg]
    Cw = (c @ c.conj().T) / deg
    rhoB = np.zeros((1024, 1024), dtype=np.complex128)
    for si, (a, t, b) in enumerate(states):
        for sj, (a2, t2, b2) in enumerate(states):
            w = Cw[si, sj]
            if abs(w) < 1e-18:
                continue
            kk = K[a, a2, t, t2].conj()
            rhoB += w * np.einsum("b,kb,lb->kl", kk,
                                  Is[b], Is[b2].conj())
    rhoB = (rhoB + rhoB.conj().T) / 2
    evr = np.linalg.eigvalsh(rhoB)
    evr = evr[evr > 1e-14]
    S = float(-np.sum(evr * np.log(evr)))
    return {"n_states": ns, "gram_rank": rank, "E0": E0, "deg": deg,
            "gap": gap, "S_bond": S,
            "spectrum": [float(x) for x in ev[:8]]}


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V288 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    ln45_art = float(art("V99_three_site_chain", "data.ln45"))
    factor_art = float(art("V99_three_site_chain",
                           "data.bekenstein_factor"))
    cap_art = float(art("V100_vertex_correction", "data.S_capacity"))
    vc_art = float(art("V100_vertex_correction", "data.vertex_coefficient"))
    a101 = art("V101_bond_current_ground_state", "data")

    Is, Ts = V101.build_kernel(t0)                # cached builder
    Ls_np = V79.adjoint_operators if hasattr(V79, "adjoint_operators") \
        else V79.adjoint_link_operators(
            V79.structure_constants(V79.root_generators())[0])
    log("kernel loaded: Is {} Ts {} (order r = 1, 45, 210)"
        .format(Is.shape, Ts.shape))

    # ---- T0: full kernel reproduces V101 ---------------------------------
    full = eigen_and_entropy(Is, Ts, Ls_np)
    dE = abs(full["E0"] - float(a101["E0"])) / max(1e-9, abs(
        float(a101["E0"])))
    t0ok = (dE < 1e-6 and abs(full["S_bond"] - float(a101["S_ground"]))
            < 1e-6 and full["deg"] == int(a101["degeneracy"])
            and abs(full["gap"] - float(a101["gap"])) < 1e-6
            and abs(ln45_art - math.log(45.0)) < 1e-12)
    data["t0"] = {"E0": full["E0"], "E0_art": a101["E0"],
                  "deg": full["deg"], "gap": full["gap"],
                  "S_bond": full["S_bond"],
                  "S_art": a101["S_ground"],
                  "ln45_art": ln45_art,
                  "bekenstein_factor_art": factor_art,
                  "S_capacity_art": cap_art,
                  "vertex_coefficient_art": vc_art,
                  "E0_rel_defect": dE}
    checks.append({"name": "T0 V101 reproduced: E0 {:.6f} (defect "
                          "{:.1e}), deg {}, gap {:.6f}, S {:.6f} vs "
                          "artefact {:.6f}; ln45 {:.13f}; references "
                          "factor {:.4f} S_cap {:.4f} vc {:.4f}"
                   .format(full["E0"], dE, full["deg"], full["gap"],
                           full["S_bond"], float(a101["S_ground"]),
                           ln45_art, factor_art, cap_art, vc_art),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V288 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the truncation map ----------------------------------------------
    subsets = {"1": (0, 2), "45": (2, 4), "210": (4, 6),
               "1+45": (0, 4), "1+45+210": (0, 6)}
    rows = []
    for name, (lo, hi) in subsets.items():
        r = eigen_and_entropy(Is, Ts[lo:hi], Ls_np)
        r["subset"] = name
        r["factor_to_quarter"] = r["S_bond"] / 0.25
        closed = 0.98 * 0.25 <= r["S_bond"] <= 1.02 * 0.25
        r["within_2pct_of_quarter"] = bool(closed)
        rows.append(r)
        log("subset {:9s} ns {:3d} rank {:3d} E0 {:+.4f} deg {} gap "
            "{:.4f} S {:.4f} factor {:.3f}"
            .format(name, r["n_states"], r["gram_rank"], r["E0"],
                    r["deg"], r["gap"], r["S_bond"],
                    r["factor_to_quarter"]))
    data["rows"] = rows
    closed_any = any(r["within_2pct_of_quarter"] for r in rows)
    status = "BEKENSTEIN-CLOSED" if closed_any else "BEKENSTEIN-MAP"
    data["status"] = status
    print("V288 {}: map booked over {} truncations; smallest "
          "factor {:.3f} at subset {}"
          .format(status, len(rows),
                  min(rows, key=lambda r: abs(r["factor_to_quarter"] -
                                              1.0))["factor_to_quarter"],
                  min(rows, key=lambda r: r["factor_to_quarter"])
                  ["subset"]), flush=True)
    checks.append({"name": "map rows: {} truncations booked; closed test "
                          "within 2 % of 1/4: {}".format(len(rows),
                                                         closed_any),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V288 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V288_bekenstein_truncation", "tier": 2,
           "property": "flux-truncation map of the three-site kernel "
                       "against the Bekenstein 1/4",
           "residual": True,
           "programme": "V288",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v288_bekenstein_truncation.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V288 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
