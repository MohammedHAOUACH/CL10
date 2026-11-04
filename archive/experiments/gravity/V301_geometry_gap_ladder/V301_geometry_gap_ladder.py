#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V301 - can the geometry ensemble's screening mass sit at the ladder
loop base? (tier 2, hard, O8-new.)
**The question (pre-registered).** Can the non-perturbative screening
mass of the geometry statistical ensemble - the only scale the
geometry sector demonstrably generates - sit at the ladder loop base:
is m_gap (lattice/Planck units, ell_P = 1 by Axiom 1) inside the
ladder band [1.5e-2, 3.0e-2] around 1/(16 pi) at any beta?
**The object (closed input set).** V300's machinery at L = 6 with
beta in {4, 8, 16, 32}, seed 301 (the extra weak-coupling point
chases the smallest mass); T0 = the V289 beta = 1 protocol-row
replay (seed-289 draw sequence, tol 1e-6) plus the booked area-law
slopes re-derived from the V289 artifact rows.
**Measured (pre-registered).** m_gap(beta) = the direction-1
line-correlator mass (log-linear on d = 1..3, defined iff all C > 0)
with jackknife error over replicas; the Wilson-loop mass m_W(beta) =
-ln<W(2)>/4 as the independent cross-read (the V289 W-accumulation
convention); direction 0 booked too; the band crossing (beta*) if
any; the minimum m_gap over the booked beta range.
**Branch.** GEOMGAP-IN-BAND iff some beta has m_gap in
[1.5e-2, 3.0e-2] (beta*, m_gap +/- jackknife, m_W cross-read booked);
GEOMGAP-OFF-BAND otherwise, the curve and the distance to the band
booked. A GEOMGAP-IN-BAND hit books a scale coincidence only: the
identification of the screening mass with the ladder base remains
the declared step O8 already states, and no reinterpretation is
added by this run.
"""
import glob
import importlib.util
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

_spec141 = importlib.util.spec_from_file_location(
    "v141", os.path.join(exp("V141_regular_link_quintet_coupling"),
                         "V141_regular_link_quintet_coupling.py"))
V141 = importlib.util.module_from_spec(_spec141)
_spec141.loader.exec_module(V141)

import cupy as cp                                         # noqa: E402

BETAS_T0 = [1.0]           # the L=4 seed-289 replay prefix
BETAS = [4.0, 8.0, 16.0, 32.0]
L4, NREP4 = 4, 32
L6, NREP6 = 6, 16
SWEEPS = 20000
EQUIL = 2000
SAMPLE_EVERY = 20
BASEPOINTS4 = 16
BASEPOINTS = 32
DVALS = (1, 2, 3)
SEED_REPLAY = 289
SEED = 301
BAND = (1.5e-2, 3.0e-2)   # the booked ladder band around 1/(16 pi)


def art(name):
    d = exp(name)
    want = name.lower() + ".json"
    hits = [h for h in glob.glob(os.path.join(d, "results", "*.json"))
            if os.path.basename(h).lower() == want]
    with open(hits[0]) as fh:
        return json.load(fh)


def area_fit(row):
    Rs = np.array([1.0, 4.0, 9.0])
    yv = np.log(np.array([row["W1"], row["W2"], row["W3"]]))
    A = np.vstack([Rs, np.ones(3)]).T
    return np.linalg.lstsq(A, yv, rcond=None)[0]


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V301 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    # ---- group / adjoint rep / tables (V289 construction, live) --------
    G = V141.build_group()
    ng = len(G)
    idx = {g: i for i, g in enumerate(G)}
    checks.append({"name": "T0a group order re-measured: {} == 1024"
                          .format(ng), "pass": ng == 1024})

    def _np(x):
        return x.get() if hasattr(x, "get") else np.asarray(x)
    from cl10 import clifford                              # noqa: E402
    from cl10 import sm_qca                                # noqa: E402
    Hc = [np.array(_np(m), dtype=np.complex128)
          for m in clifford.generators()]

    def rho(g):
        M = np.eye(32, dtype=np.complex128)
        for a in g[1]:
            M = M @ Hc[a]
        return g[0] * M

    Ts = [np.array(_np(T), dtype=np.complex128)
          for T in sm_qca.plane_generators()]
    normT = np.array([float(np.real(np.trace(Ta.conj().T @ Ta)))
                      for Ta in Ts])

    def adjoint_mat(gm):
        Uu = rho(gm)
        Ui = np.linalg.inv(Uu)
        Rd = np.zeros((45, 45))
        act = [Uu @ Ta @ Ui for Ta in Ts]
        resid = 0.0
        for A in range(45):
            for B in range(45):
                Rd[B, A] = float(np.real(
                    np.trace(Ts[B].conj().T @ act[A])) / normT[B])
            rec = sum(Rd[B, A] * Ts[B] for B in range(45))
            resid = max(resid, float(np.abs(act[A] - rec).max()))
        return Rd, resid

    order = [(1, ())]
    seen = {order[0]}
    parent = {}
    front = [order[0]]
    while front:
        nf = []
        for g in front:
            for s in V141.BIFLECTIONS:
                m = V141.mul(g, s)
                if m not in seen:
                    seen.add(m)
                    parent[m] = (g, s)
                    nf.append(m)
                    order.append(m)
        front = nf
    assert len(order) == ng
    Rm = np.zeros((ng, 45, 45))
    Rm[0] = np.eye(45)
    for g in order[1:]:
        pg, s = parent[g]
        Rs, _ = adjoint_mat(s)
        Rm[idx[g]] = Rm[idx[pg]] @ Rs
    chi = np.real(np.trace(Rm, axis1=1, axis2=2))
    multab = np.zeros((ng, ng), dtype=np.int32)
    for i, gi in enumerate(G):
        for j, gj in enumerate(G):
            multab[i, j] = idx[V141.mul(gi, gj)]
    invtab = np.array([idx[V141.invof(g)] for g in G], dtype=np.int32)
    props = V141.BIFLECTIONS + [V141.invof(s)
                                for s in V141.BIFLECTIONS]
    nprop = len(props)
    tab = np.zeros((ng, nprop), dtype=np.int32)
    for i, g in enumerate(G):
        for j, s in enumerate(props):
            tab[i, j] = multab[i, idx[s]]
    Rg = cp.asarray(Rm)
    Rgi = cp.asarray(Rm[invtab])
    tabc = cp.asarray(tab)
    multabc = cp.asarray(multab)
    invtabc = cp.asarray(invtab)
    chig = cp.asarray(chi)
    log("group/rep/tables built")

    def make_lattice(L):
        ns = L ** 3
        planes = [(0, 1), (0, 2), (1, 2)]
        plaqs = []
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    s0 = x + L * y + L * L * z
                    for (mu, nu) in planes:
                        dmu = (1, L, L * L)[mu]
                        dnu = (1, L, L * L)[nu]
                        xs = (s0 + dmu) % ns
                        xn = (s0 + dnu) % ns
                        plaqs.append((s0 * 3 + mu, xs * 3 + nu,
                                      xn * 3 + mu, s0 * 3 + nu))
        plaqs = np.array(plaqs, dtype=np.int64)
        nb = 3 * ns
        bond_pl = [[] for _ in range(nb)]
        for pi in range(len(plaqs)):
            for slot in range(4):
                bond_pl[int(plaqs[pi, slot])].append(pi)
        return ns, planes, plaqs, bond_pl

    def run_metropolis(L, NREP, betas, seed, corr):
        ns, planes, plaqs, bond_pl = make_lattice(L)
        nplaq = len(plaqs)
        nb = 3 * ns
        plaqs_c = cp.asarray(plaqs)

        def plaq_chain(uu, plist):
            b = plaqs_c[plist]
            return (Rg[uu[:, b[:, 0]]] @ Rg[uu[:, b[:, 1]]]) @ \
                   (Rgi[uu[:, b[:, 2]]] @ Rgi[uu[:, b[:, 3]]])

        def plaq_gidx(u, plist):
            b = plaqs[plist]
            g12 = multabc[u[:, b[:, 0]], u[:, b[:, 1]]]
            g34 = multabc[invtabc[u[:, b[:, 2]]], invtabc[u[:, b[:, 3]]]]
            return multabc[g12, g34]

        all_pl = np.arange(nplaq)
        reps_all = cp.arange(NREP)
        rng = np.random.default_rng(seed)
        out = []
        chk_max, nchk = 0.0, 0
        for beta in betas:
            u = cp.asarray(rng.integers(0, ng, (NREP, nb), dtype=np.int64))
            chiP = cp.real(cp.trace(plaq_chain(u, all_pl),
                                    axis1=2, axis2=3))
            gp = plaq_gidx(u, all_pl)
            chk_max = max(chk_max, float(cp.abs(chig[gp] - chiP).max()))
            nchk += NREP * nplaq
            acc = 0
            wrun, nsam = np.zeros((NREP, 3)), 0
            csum = np.zeros((NREP, 2, len(DVALS)))
            for sweep in range(SWEEPS):
                bond = int(rng.integers(nb))
                pidx = rng.integers(0, nprop, NREP)
                pidxc = cp.asarray(pidx)
                ubak = u[reps_all, bond].copy()
                plist = np.array([int(p) for p in bond_pl[bond]])
                oldc = chiP[:, plist]
                unew = tabc[u[reps_all, bond], pidxc]
                u[reps_all, bond] = unew
                newc = cp.real(cp.trace(plaq_chain(u, plist),
                                        axis1=2, axis2=3))
                accept = cp.asarray(rng.random(NREP)) < cp.exp(
                    beta * cp.sum(newc - oldc, axis=1) / 45.0)
                acc += int(cp.sum(accept))
                keep = cp.where(~accept)[0]
                if len(keep):
                    u[keep, bond] = ubak[keep]
                chiP[:, plist] = cp.where(accept[:, None], newc, oldc)
                if sweep >= EQUIL and sweep % SAMPLE_EVERY == 0:
                    cpl = rng.choice(nplaq, 4, replace=False)
                    gp = plaq_gidx(u, cpl)
                    cc = cp.real(cp.trace(plaq_chain(u, cpl),
                                          axis1=2, axis2=3))
                    chk_max = max(chk_max,
                                  float(cp.abs(chig[gp] - cc).max()))
                    nchk += 4 * NREP
                    sel = rng.integers(0, nplaq, BASEPOINTS4)
                    wr = np.zeros((NREP, 3))
                    for pi in sel:
                        mu, nu = planes[int(pi) % 3]
                        s0 = int(pi) // 3
                        dmu = (1, L, L * L)[mu]
                        dnu = (1, L, L * L)[nu]
                        for R in (1, 2, 3):
                            tot = cp.zeros(NREP)
                            cnt = 0
                            for sh in range(R):
                                for th in range(R):
                                    cur = (s0 + sh * dmu
                                           + th * dnu) % ns
                                    P = cp.eye(45)[None].repeat(NREP, 0)
                                    c = cur
                                    for _ in range(R):
                                        P = P @ Rg[u[:, c * 3 + mu]]
                                        c = (c + dmu) % ns
                                    for _ in range(R):
                                        P = P @ Rg[u[:, c * 3 + nu]]
                                        c = (c + dnu) % ns
                                    for _ in range(R):
                                        c = (c - dmu) % ns
                                        P = P @ Rgi[u[:, c * 3 + mu]]
                                    for _ in range(R):
                                        c = (c - dnu) % ns
                                        P = P @ Rgi[u[:, c * 3 + nu]]
                                    tot += cp.real(cp.trace(
                                        P, axis1=1, axis2=2))
                                    cnt += 1
                            wr[:, R - 1] += cp.asnumpy(
                                tot / (45.0 * cnt))
                            del P
                    wr /= BASEPOINTS4
                    wrun = wr     # V289 accumulation convention, verbatim
                    if corr:
                        bps = rng.integers(0, ns, BASEPOINTS)
                        for di, mu in enumerate((0, 1)):
                            dm = (1, L, L * L)[mu]
                            P = cp.eye(45)[None].repeat(
                                BASEPOINTS * NREP, 0)
                            s2r = ((bps[:, None] * NREP
                                    + np.arange(NREP)[None, :]).ravel())
                            for d in range(1, max(DVALS) + 1):
                                cols = np.repeat(
                                    (bps + (d - 1) * dm) % ns, NREP)
                                bondidx = cols * 3 + mu
                                P = P @ Rg[u[s2r, bondidx]]
                                if d in DVALS:
                                    di_ = DVALS.index(d)
                                    tr = cp.real(cp.trace(
                                        P, axis1=1, axis2=2)) / 45.0
                                    csum[:, di, di_] += cp.asnumpy(
                                        tr.reshape(BASEPOINTS, NREP)
                                        .sum(axis=0))
                    nsam += 1
            out.append({"beta": beta, "accept": acc / (SWEEPS * NREP),
                        "nsamples": nsam,
                        "W": wrun.mean(axis=0).tolist(),
                        "C": (csum / max(nsam, 1)).tolist()})
            log("beta {:5.1f}: accept {:.4f} nsam {}".format(
                beta, out[-1]["accept"], nsam))
        return out, chk_max, nchk

    # ---- T0a: V289 seed-289 replay, beta = 1 ---------------------------
    art289 = art("V289_geometry_mc")
    rows_art = {float(r["beta"]): r for r in art289["data"]["rows"]}
    t0rows, chk0, nchk0 = run_metropolis(L4, NREP4, BETAS_T0, SEED_REPLAY,
                                         corr=False)
    ra = rows_art[1.0]
    t0dev = max(abs(t0rows[0]["accept"] - float(ra["accept"])),
                abs(t0rows[0]["W"][0] - float(ra["W1"])))
    data["t0_replay"] = {"rows": t0rows,
                         "artifact": [ra["accept"], ra["W1"]],
                         "max_abs_dev": float(t0dev)}
    checks.append({"name": "T0a V289 seed-{} replay beta 1 accept/W1: max "
                          "dev {:.2e} (tol 1e-6)".format(SEED_REPLAY,
                                                         t0dev),
                   "pass": bool(t0dev < 1e-6)})

    # ---- T0b: area-law slope anchors, re-derived from artifact rows ----
    slo_art = {float(f["beta"]): float(f["area_slope"])
               for f in art289["data"]["area_fits"]}
    slo_red = {}
    for r in art289["data"]["rows"]:
        if all(r[k] > 0 for k in ("W1", "W2", "W3")):
            slo_red[float(r["beta"])] = float(area_fit(r)[0])
    dslo = max(abs(slo_red[b] - slo_art[b]) for b in slo_art)
    data["t0_area_fit"] = {"artifact": slo_art, "rederived": slo_red,
                           "max_abs_dev": float(dslo)}
    checks.append({"name": "T0b area-law slopes re-derived from the V289 "
                          "artifact rows: max dev {:.2e}".format(dslo),
                   "pass": bool(dslo < 1e-9)})
    if t0dev >= 1e-6 or dslo >= 1e-9:
        data["status"] = "MACHINERY-DRIFT"
        print("V301 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: L=6 gap measurement ---------------------------------------
    rows, chk1, nchk1 = run_metropolis(L6, NREP6, BETAS, SEED, corr=True)
    data["char_check"] = {"n_checked": int(nchk0 + nchk1),
                          "max_abs_defect": float(max(chk0, chk1))}

    def masses(C):
        C = np.asarray(C)
        bars = C.mean(axis=0)

        def fit(curv):
            ms = []
            for di in range(2):
                if np.any(curv[di] <= 0):
                    ms.append(None)
                else:
                    A = np.vstack([np.array(DVALS, dtype=float),
                                   np.ones(len(DVALS))]).T
                    sol = np.linalg.lstsq(A, np.log(curv[di]),
                                          rcond=None)[0]
                    ms.append(-float(sol[0]))
            return ms

        ms = fit(bars)
        merr = [None] * 2
        if all(m is not None for m in ms):
            jm = np.stack([fit(np.delete(C, rr, axis=0).mean(axis=0))
                           for rr in range(C.shape[0])])
            for di in range(2):
                v = jm[:, di].astype(float)
                merr[di] = float(np.sqrt(C.shape[0] / (C.shape[0] - 1)
                                         * ((v - v.mean()) ** 2).mean()))
        return bars, ms, merr

    vrows = []
    for r in rows:
        bars, ms, merr = masses(r["C"])
        mW = (float(-np.log(r["W"][1]) / 4.0) if r["W"][1] > 0
              else None)
        vrows.append({"beta": r["beta"], "accept": r["accept"],
                      "C": bars.tolist(), "m0": ms[0], "m0_err": merr[0],
                      "m_gap": ms[1], "m_gap_err": merr[1],
                      "W": r["W"], "m_W": mW,
                      "in_band": bool(ms[1] is not None
                                      and BAND[0] <= ms[1] <= BAND[1])})
        log("beta {:5.1f}: m0 {} m_gap {} m_W {}".format(
            r["beta"], ms[0], ms[1], mW))
    data["rows"] = vrows
    data["band"] = list(BAND)

    hits = [v for v in vrows if v["in_band"]]
    defined = [v for v in vrows if v["m_gap"] is not None]
    if hits:
        status = "GEOMGAP-IN-BAND"
        dmin = min(abs(v["m_gap"] - 1.0 / (16 * np.pi)) for v in hits)
    else:
        status = "GEOMGAP-OFF-BAND"
        dmin = (min(min(abs(v["m_gap"] - BAND[0]), abs(v["m_gap"] - BAND[1]))
                    for v in defined) if defined else None)
    data["status"] = status
    data["closest_band_dev"] = (None if dmin is None else float(dmin))
    print("V301 {}: m_gap(beta) {} band {}; m_W {} ; closest dev to band "
          "{}. Scale-coincidence note: even a hit books the coincidence "
          "only - the identification with the ladder base stays the "
          "declared O8 step".format(
              status,
              [(v["beta"], None if v["m_gap"] is None
                else round(v["m_gap"], 5)) for v in vrows],
              list(BAND),
              [(v["beta"], v["m_W"]) for v in vrows],
              data["closest_band_dev"]), flush=True)
    checks.append({"name": "adjoint character cross-check on {} plaquette "
                          "products: max defect {:.2e}".format(
                              nchk0 + nchk1, max(chk0, chk1)),
                   "pass": bool(max(chk0, chk1) < 1e-8)})
    fin = all(bool(np.all(np.isfinite([v["accept"]] + v["W"]
                                       + [x for row in v["C"]
                                          for x in row])))
              for v in vrows)
    checks.append({"name": "correlator table finite over {} betas".format(
        len(vrows)), "pass": bool(fin)})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok=ok)
    print("[V301 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V301_geometry_gap_ladder", "tier": 2,
           "property": "does the geometry ensemble screening mass land "
                       "in the ladder loop band",
           "residual": True, "programme": "V301",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v301_geometry_gap_ladder.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V301 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
