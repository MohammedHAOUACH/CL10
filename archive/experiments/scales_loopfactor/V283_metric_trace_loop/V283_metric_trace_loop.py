#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V283 - the metric trace loop (tier 2, hard, O8-new).

**The question (pre-registered).** Does the metric sector - the
directional quintet modes - traced in its ground state produce the
one-loop Cartan ratio near 1/(16 pi) in the closed two-hop traced
functional?

Construction (locked here, the pre-reg object made explicit):
 - the two-hop: V252's dressed contact resolvent operator at the
   booked floor row (implemented clock, star node 0,0,0, delta =
   1e-4, sign -, contact limit): D = G10^2 sum_s V diag(1/(z - E))
   V^dag sandwiched through the flux contact superoperator, on the
   128-dim cell; V252's exact sandwich/cnorm128/Ksuper blocks;
 - the quintet: V140's sym^2 machinery (SB basis, J6 rotation
   generators, M_dir at the committed clock ETA = G10) and V229's
   quintet vertex matrices Emat(gamma), gamma = 0..4; the quintet
   ground mode gamma* is the band of least integrated |phase| over
   V140's fib_dirs(100) at the committed clock (T0 booked row);
 - the insertion: the cell factor 8 is the binary sublattice triple
   (the supercell convention); the direction operators d_a := sigma_z
   on factor a (a = 0,1,2) give the only canonical 3 x Z_2 structure
   on 8, and S(gamma) := sum_ab Emat(gamma)_ab d_a d_b embeds the
   quintet vertex on the 128-cell (booked embedding choice);
 - the ladder: F_n := cnorm128(W^n) / cnorm128(D^n) with
   W = D + g10 S(gamma*), n <= 30, at N = 16 and N = 32 (the
   step-doubling stability row); the closest |16 pi F - 1| and its
   (n, N, direction) label.

T0: the booked V252 floor row re-derived live to 1e-6 relative, the
V140 quintet band multiplicities at the committed clock, the spin-2
insertion norms. Acceptance: T0 fail -> MACHINERY-DRIFT;
METRICLOOP-FOUND iff a stable (both N rows agreeing to 1e-6) n <= 30
has |16 pi F - 1| < 0.5; else METRICLOOP-NEGATIVE with the floor.
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")
from cl10.paths import exp                                # noqa: E402

sys.path.insert(0, exp("V79_gauged_bond_layer"))
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)
from cl10 import altcell as AC                            # noqa: E402
from cl10 import clifford                                 # noqa: E402
from cl10 import coupling                                 # noqa: E402
from cl10 import engine as E                              # noqa: E402
from cl10 import supercell as SC                          # noqa: E402

G10 = coupling.G10                                        # 1/5


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def load_mod(name):
    path = exp(name, "%s.py" % name)
    spec = importlib.util.spec_from_file_location(name.lower(), path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    checks, data = [], {}
    xp = E._E if E._BACKEND == "cupy" else np
    print("V283 - the metric trace loop [backend: {}]".format(
        E._BACKEND), flush=True)

    # ---- machinery (the V252 set, imported) ------------------------------
    Ts = root_generators()
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))        # 1440
    B, _U32, _U16_ = SC.sectors(G10)
    from cl10 import altcell as _AC
    _Bm, T16, _f16 = _AC.weyl_generators()
    T16d = xp.asarray(np.stack(T16))
    HNP = [E.asnumpy(m) for m in clifford.generators()]
    CARNP = [1j * (HNP[a] @ HNP[a + 1]) for a in (0, 2, 4, 6, 8)]
    CAR16 = [B.conj().T @ C @ B for C in CARNP]

    def cnorm128(M):
        M4 = M.reshape(8, 16, 8, 16)
        return float(sum(abs(np.einsum("ab,xbxa->", C, M4))
                         for C in CAR16))

    # Casimir tower reduced weights (V252 verbatim; contact limit only)
    tev, TU = np.linalg.eigh(TL)
    eigs = (-32.0, -12.0, 4.0)
    wts = {}
    for s in eigs:
        idx = np.abs(tev - s) < 1e-9
        P = TU[:, idx] @ TU[:, idx].conj().T
        r = np.einsum("aibi->ab", P.reshape(32, 45, 32, 45)) / 32.0
        wts[s] = np.real(np.diag(r))

    def Ksuper(d):
        return np.einsum("A,Aac,Abd->abcd", d,
                         E.asnumpy(T16d), E.asnumpy(T16d).conj())

    def sandwich(M4, K2):
        Mm = M4.transpose(0, 2, 4, 1, 3).reshape(-1, 256, 64)
        Q = K2.reshape(256, 256) @ Mm
        return (Q.reshape(-1, 16, 16, 8, 8).transpose(0, 3, 1, 4, 2)
                .reshape(-1, 128, 128))

    # ---- T0a: the V252 booked floor row re-derived live -------------------
    ft = art("V252_dynamical_self_energy_loop", "data.F_int")
    # booked floor row: implemented / star / (0,0,0), d1e-4, sign -, N16,
    # contact (argmin |16 pi F - 1| over the booked table, read live)
    best_row, best_dev, best_nodekey = None, None, None
    for nodekey, rows in ft.items():
        clk, cls, qs = nodekey.split("/")
        for rk, F in rows.items():
            dev = abs(16.0 * math.pi * F - 1.0)
            if best_dev is None or dev < best_dev:
                dlstr, Nstr, oname = rk.split("/")
                sign = -1.0 if Nstr.startswith("-") else 1.0
                best_row = (clk, cls, tuple(float(x) for x in
                                            qs.split(",")),
                            float(dlstr[1:]), sign,
                            int(Nstr.split("N")[1]), oname, F)
                best_dev = dev
                best_nodekey = nodekey
    clk, cls, kpt, dl, sign, Nrow, oname, a_F = best_row
    # node energy at the booked node, booked clock
    beta, g = G10, G10
    U16 = SC.sectors(g)[2]

    def bfn(q, b=beta, gg=g):
        return SC.bloch(q, U16, gg, bond_g10=b)
    Hk = AC.effective_hamiltonian(bfn, kpt)
    wv, Vv = np.linalg.eigh(Hk)
    En = float(wv[int(np.argmin(np.abs(wv)))])

    def dressed_operator(N):
        """(ct, D) = V252's tree norm and G10^2 contact sum at the
        booked (z, N) row, same chunking convention."""
        full = 32
        r = full // N
        shifts = np.array([(x * r, y * r, z * r)
                           for x in range(N) for y in range(N)
                           for z in range(N)]) * (2 * math.pi / full)
        z = En + sign * 1j * dl
        K2c = xp.asarray(Ksuper(np.ones(45)))
        acc = xp.zeros((128, 128), dtype=np.complex128)
        tree = None
        ch = 2048
        for c0 in range(0, len(shifts), ch):
            sl = slice(c0, min(c0 + ch, len(shifts)))
            blk = xp.asarray(np.stack(
                [AC.effective_hamiltonian(bfn, tuple(s))
                 for s in shifts[sl]]))
            Ec, Vc = xp.linalg.eigh(blk)
            del blk
            w = 1.0 / (z - Ec)
            M = (Vc * w[:, None, :]) @ Vc.conj().swapaxes(-2, -1)
            S = sandwich(M.reshape(-1, 8, 16, 8, 16), K2c)
            acc += xp.sum(S, axis=0)
            if tree is None:
                tree = E.asnumpy(S[0])
        return cnorm128(tree), E.asnumpy(acc) * G10 ** 2

    cta, Da = dressed_operator(Nrow)
    F16 = cnorm128(Da) / cta
    t0a_ok = abs(F16 - a_F) <= 1e-6 * abs(a_F)
    data["v252_floor_row"] = {"row": [clk, cls, list(kpt), dl, sign,
                                      Nrow, oname],
                              "F_artefact": a_F, "F_recomputed_N16": F16,
                              "abs_dev_16piF_minus_1": best_dev,
                              "En": En, "tree_norm_ct": cta}
    checks.append({"name": "T0a V252 booked floor row re-derived: "
                          "{} d{:.0e} sign {:+d} N{} contact: F {:.15e} "
                          "vs booked {:.15e} (rel {:.1e}), floor "
                          "|16 pi F - 1| = {:.4f}".format(
                              best_nodekey, dl, int(sign), Nrow, F16, a_F,
                              abs(F16 - a_F) / abs(a_F), best_dev),
                   "pass": bool(t0a_ok)})
    if not t0a_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V283 MACHINERY-DRIFT: the booked V252 floor row does not "
              "re-derive; nothing downstream interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0b: the V140 quintet at the committed clock --------------------
    V140 = load_mod("V140_dressed_gate_quintet")
    V229 = load_mod("V229_diagonal_bond_quintet")
    SB = V140.SB
    # quintet basis: the traceless part of the sym^2 span - the three
    # off-diagonal elements plus two normalized diagonal differences
    off = [M for M in SB if abs(np.trace(M)) < 1e-12]
    dd = [M for M in SB if abs(np.trace(M)) > 0.5]
    q_basis = [(dd[0] - dd[1]) / math.sqrt(2.0),
               (dd[1] - dd[2]) / math.sqrt(2.0)] + off
    bands_rows, mode_cost = [], [0.0] * 5
    dirs = V140.fib_dirs(100)
    Emats = [V229.Emat(gg) for gg in range(5)]
    for qq in dirs:
        Mq = V140.M_dir(qq)
        ph = np.sort(np.angle(np.linalg.eigvals(Mq)))
        bands_rows.append([float(p) for p in ph])
        for i in range(5):
            mode_cost[i] += float(abs(ph[i]))
    gstar = int(np.argmin(mode_cost))
    band_mult = {}
    for row in bands_rows:
        for p in row:
            key = round(p, 6)
            hit = next((kk for kk in band_mult
                        if abs(kk - key) < 1e-4), None)
            if hit is None:
                band_mult[key] = 1
            else:
                band_mult[hit] += 1
    ins_norm = float(max(np.linalg.norm(M) for M in Emats))
    data["quintet"] = {"mode_costs": [round(v / len(dirs), 6)
                                      for v in mode_cost],
                       "gamma_star": gstar,
                       "band_multiplicities": {str(k): v for k, v
                                               in sorted(band_mult.items())},
                       "insertion_norm_max": ins_norm,
                       "n_traceless_modes": len(q_basis)}
    checks.append({"name": "T0b quintet at committed clock: mode costs "
                          "{} -> gamma* = {}, insertion norms max "
                          "{:.4f}, 5 traceless sym^2 modes".format(
                              [round(v / len(dirs), 5) for v in
                               mode_cost], gstar, ins_norm),
                   "pass": bool(len(q_basis) >= 5 and ins_norm > 0.5)})

    # ---- T1: the traced two-hop ladder with the insertion -----------------
    # direction operators on the sublattice triple factor of the cell
    s_z = np.diag([1.0, -1.0])
    I2 = np.eye(2)
    d_ops = [np.kron(np.kron(s_z if a == 0 else I2,
                             s_z if a == 1 else I2),
                     s_z if a == 2 else I2) for a in range(3)]
    S_g = np.zeros((8, 8))
    for a in range(3):
        for b in range(3):
            S_g += Emats[gstar][a, b].real * (d_ops[a] @ d_ops[b])
    S128 = np.kron(S_g, np.eye(16))
    S128 = S128 / max(np.linalg.norm(S_g), 1e-30)   # unit insertion norm

    results = {}
    if Nrow == 16:
        ct16, D16 = cta, Da
        ct32, D32 = dressed_operator(32)
    else:
        ct32, D32 = cta, Da
        ct16, D16 = dressed_operator(16)
    for N, Dx, ct in ((16, D16, ct16), (32, D32, ct32)):
        Db = Dx
        W = Db + G10 * S128
        Xn = np.eye(128, dtype=np.complex128)
        Bn = np.eye(128, dtype=np.complex128)
        rows = {}
        for n in range(1, 31):
            Xn = Xn @ W
            Bn = Bn @ Db
            cb = cnorm128(Bn)
            cx = cnorm128(Xn)
            rows[n] = cx / cb if cb > 1e-300 else None
        results[N] = rows
    # stability row: F_n agreement between N16 and N32
    stab = max(abs(results[16][n] - results[32][n])
               for n in range(1, 31)
               if results[16][n] is not None
               and results[32][n] is not None)
    best = None
    for N, rows in results.items():
        for n, F in rows.items():
            if F is None or results[16][n] is None or results[32][n] is None:
                continue
            if abs(results[16][n] - results[32][n]) > 1e-6:
                continue
            dev = abs(16.0 * math.pi * F - 1.0)
            if best is None or dev < best[0]:
                best = (dev, n, N, F)
    data["ladder"] = {"N16": {str(n): (v if v is None else float(v))
                              for n, v in results[16].items()},
                      "N32": {str(n): (v if v is None else float(v))
                              for n, v in results[32].items()},
                      "step_doubling_disagreement": stab,
                      "gamma_star": gstar,
                      "insertion": "S(gamma*)/|S| on the sublattice-"
                                   "triple direction embedding"}
    found = best is not None and best[0] < 0.5
    status = "METRICLOOP-FOUND" if found else "METRICLOOP-NEGATIVE"
    data["status"] = status
    data["closest"] = ({"dev": best[0], "n": best[1], "N": best[2],
                        "F": best[3]} if best is not None else None)
    print("V283 {}: closest |16 pi F - 1| = {} over stable n <= 30 "
          "(gamma* {}, step-doubling spread {:.1e})".format(
              status,
              "none above floor" if best is None else
              "{:.4f} at n={} N={}".format(best[0], best[1], best[2]),
              gstar, stab), flush=True)

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V283 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V283_metric_trace_loop", "tier": 2,
           "property": "the V209-type Cartan-ratio functional on the "
                       "V252 dressed two-hop with the quintet "
                       "ground-state insertion at g10: does the metric "
                       "sector close the loop (tier 2, O8-new)",
           "residual": True,
           "programme": "V283",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v283_metric_trace_loop.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V283 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
