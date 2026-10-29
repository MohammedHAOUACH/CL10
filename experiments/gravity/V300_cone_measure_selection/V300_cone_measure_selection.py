#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V300 - does the geometry measure concentrate a light cone?
(tier 2, hard, F8/O9.)
**The question (pre-registered).** The paper states the finite measure
on the geometry group (V289, MC-TUNED) "does not measurably select"
the light cone. Does that measure in fact concentrate a light cone:
on the sampled ensemble, is the screening mass of adjoint line
correlators isotropic in lattice units - an effective cone speed
v_cone = xi_1/xi_0 within [0.95, 1.05] - or does the measure prefer
no cone at all (v_cone far from 1 at every beta with acceptance >
2 %)?
**The object (closed input set).** V289's machinery re-derived live
(V141 build_group order 1024, the 45-dim adjoint rep built
multiplicatively along the generation tree, the full
multiplication/inverse/proposal tables, the adjoint-character Wilson
action, single-bond biflection Metropolis): the L = 4 protocol is
replayed for the first two beta (1, 2) of the original seed-289 draw
sequence and their booked (accept, W1) rows reproduced to 1e-6 (T0;
the character cross-check's exactness re-measured); then a new
ensemble at L = 6 (216 sites, 648 bonds, 648 based plaquettes),
beta in {2, 4, 8, 16}, NREP = 16, SWEEPS = 20000 (EQUIL 2000,
SAMPLE_EVERY 20), BASEPOINTS = 32, seed 300.
**Measured (pre-registered).** T0 row replay; the line correlators
C_mu(d) = mean_s0 Re chi_45(u_s u_{s+mu} ...)/45 along straight
forward lines of length d = 1, 2, 3 in direction mu = 0 (the step
direction) and mu = 1, averaged over samples, jackknifed over
replicas; m_mu from the 3-point log-linear fit (defined iff all
three C > 0 - the V289 fit-stage contract applies); v_cone(beta) =
m_0/m_1 with its jackknife error.
**Branch.** CONE-MEASURE-SELECTED iff |v_cone - 1| < 0.05 at some
beta with accept > 2 % and both masses defined; CONE-MEASURE-
ANISOTROPIC if masses are defined but no beta lands in the band (the
whole v_cone(beta) curve booked, with the implemented generator
speed 1.570796 read from the V285 artifact as the comparison
number); CONE-MEASURE-NO-MASS if no beta admits defined masses (the
correlator table booked). Pass = the measurement is complete in
every branch.
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

BETAS_T0 = [1.0, 2.0]        # the L=4 seed-289 replay prefix
BETAS = [2.0, 4.0, 8.0, 16.0]
L4, NREP4 = 4, 32            # the V289 protocol, verbatim
L6, NREP6 = 6, 16
SWEEPS = 20000
EQUIL = 2000
SAMPLE_EVERY = 20
BASEPOINTS4 = 16             # V289 W-loop basepoints (replay)
BASEPOINTS = 32              # V300 correlator base sites
DVALS = (1, 2, 3)
SEED_REPLAY = 289
SEED = 300
BAND = 0.05


def art(name):
    d = exp(name)
    want = name.lower() + ".json"
    hits = [h for h in glob.glob(os.path.join(d, "results", "*.json"))
            if os.path.basename(h).lower() == want]
    with open(hits[0]) as fh:
        return json.load(fh)


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


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V300 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    # ---- the group, the 45-dim adjoint rep, the tables (V289 verbatim) -
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

    sys.path.insert(0, exp("V140_dressed_gate_quintet"))
    from V140_dressed_gate_quintet import structure_constants  # noqa: E402
    Ts = [np.array(_np(T), dtype=np.complex128)
          for T in sm_qca.plane_generators()]
    f, _k32 = structure_constants(Ts)
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

    def run_metropolis(L, NREP, betas, seed, wloops, corr, bp4):
        """The V289 sweep protocol; samples W loops (wloops, the V289
        accumulation verbatim when replay=True) and/or line correlators
        (corr). Returns per-beta dicts."""
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
                    if wloops:
                        sel = rng.integers(0, nplaq, bp4)
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
                                        P = cp.eye(45)[None].repeat(
                                            NREP, 0)
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
                        wr /= bp4
                        wrun = wr            # the V289 accumulation, verbatim
                    if corr:
                        bps = rng.integers(0, ns, BASEPOINTS)
                        for di, mu in enumerate((0, 1)):
                            dm = (1, L, L * L)[mu]
                            # line products for all (base, rep) at once,
                            # accumulating the three lengths
                            P = cp.eye(45)[None].repeat(
                                BASEPOINTS * NREP, 0)
                            s2r = ((bps[:, None] * NREP
                                    + np.arange(NREP)[None, :])
                                   .ravel())
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
            res = {"beta": beta, "accept": acc / (SWEEPS * NREP),
                   "nsamples": nsam}
            if wloops:
                # the V289 convention, verbatim: Wsam is the final
                # sample's per-replica windows; wbar = mean over replicas
                res["W"] = wrun.mean(axis=0).tolist()
            if corr:
                res["C"] = (csum / max(nsam, 1)).tolist()
            out.append(res)
            log("beta {:5.1f}: accept {:.4f} nsam {}".format(
                beta, res["accept"], nsam))
        return out, chk_max, nchk

    # ---- T0: replay of the V289 L=4 protocol, first two beta ----------
    art289 = art("V289_geometry_mc")
    rows_art = {float(r["beta"]): r for r in art289["data"]["rows"]}
    t0rows, chk0, nchk0 = run_metropolis(L4, NREP4, BETAS_T0, SEED_REPLAY,
                                         wloops=True, corr=False,
                                         bp4=BASEPOINTS4)
    devs = []
    for r in t0rows:
        ra = rows_art[r["beta"]]
        devs.append(abs(r["accept"] - float(ra["accept"])))
        devs.append(abs(r["W"][0] - float(ra["W1"])))
    t0dev = max(devs)
    data["t0_replay"] = {"rows": t0rows,
                         "artifact": [rows_art[b]["accept"]
                                      for b in BETAS_T0],
                         "max_abs_dev": float(t0dev)}
    checks.append({"name": "T0 V289 seed-{} replay of beta {} accept/W1 "
                          "rows: max dev {:.2e} (tol 1e-6)".format(
                              SEED_REPLAY, BETAS_T0, t0dev),
                   "pass": bool(t0dev < 1e-6)})
    if t0dev >= 1e-6:
        data["status"] = "MACHINERY-DRIFT"
        print("V300 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: the L=6 ensemble with line correlators --------------------
    rows, chk1, nchk1 = run_metropolis(L6, NREP6, BETAS, SEED,
                                       wloops=True, corr=True,
                                       bp4=BASEPOINTS4)
    data["char_check"] = {"n_checked": int(nchk0 + nchk1),
                          "max_abs_defect": float(max(chk0, chk1))}

    def mass_and_v(C):
        """C: (NREP, 2, 3) sample-mean correlators. Jackknife over
        replicas; log-linear 3-point fits; returns per-direction masses
        (None if any C <= 0) and v = m0/m1."""
        C = np.asarray(C)
        bars = C.mean(axis=0)                      # (2, 3)
        jackm = []
        for rr in range(C.shape[0]):
            jk = np.delete(C, rr, axis=0).mean(axis=0)
            jackm.append(jk)
        jackm = np.stack(jackm)

        def fit(curv):
            ms = []
            for di in range(2):
                if np.any(curv[di] <= 0):
                    ms.append(None)
                else:
                    A = np.vstack([np.array(DVALS, dtype=float),
                                   np.ones(len(DVALS))]).T
                    sol = np.linalg.lstsq(
                        A, np.log(curv[di]), rcond=None)[0]
                    ms.append(-float(sol[0]))
            return ms

        ms = fit(bars)
        out = {"C": bars.tolist(),
               "err": (np.sqrt(C.shape[0] / (C.shape[0] - 1)
                               * ((jackm - bars[None]) ** 2).mean(axis=0))
                       ).tolist()}
        jms = [fit(jk) for jk in jackm]
        vbar = None
        verr = 0.0
        if ms[0] is not None and ms[1] is not None:
            vbar = ms[0] / ms[1]
            vs = [j[0] / j[1] for j in jms
                  if j[0] is not None and j[1] is not None]
            if vs:
                vs = np.array(vs)
                verr = float(np.sqrt(C.shape[0] / (C.shape[0] - 1)
                                     * ((vs - np.mean(vs)) ** 2).mean()))
        out["m0"] = ms[0]
        out["m1"] = ms[1]
        out["v_cone"] = None if vbar is None else float(vbar)
        out["v_err"] = verr
        return out

    vrows = []
    for r in rows:
        mv = mass_and_v(r["C"])
        vrows.append({"beta": r["beta"], "accept": r["accept"],
                      "W": r["W"], "C": mv["C"], "C_err": mv["err"],
                      "m0": mv["m0"], "m1": mv["m1"],
                      "v_cone": mv["v_cone"], "v_err": mv["v_err"]})
        log("beta {:5.1f}: m0 {} m1 {} v {}".format(
            r["beta"], mv["m0"], mv["m1"], mv["v_cone"]))
    data["rows"] = vrows

    # ---- branch ---------------------------------------------------------
    cgw = 1.570796     # the booked implemented-generator speed, paper
                       # Section 11.6 (V208/V211), comparison number only
    defined = [v for v in vrows if v["v_cone"] is not None]
    inband = [v for v in defined if v["accept"] > 0.02
              and abs(v["v_cone"] - 1.0) < BAND]
    if inband:
        status = "CONE-MEASURE-SELECTED"
    elif defined:
        status = "CONE-MEASURE-ANISOTROPIC"
    else:
        status = "CONE-MEASURE-NO-MASS"
    data["status"] = status
    data["comparison_implemented_speed"] = cgw
    print("V300 {}: v_cone(beta) {} (band |v-1| < {}); masses defined "
          "at {} of {} betas; implemented-generator comparison {} (V285 "
          "booked)".format(
              status,
              [(v["beta"], None if v["v_cone"] is None
                else round(v["v_cone"], 4)) for v in vrows],
              BAND, len(defined), len(vrows), cgw), flush=True)
    checks.append({"name": "adjoint character cross-check on {} plaquette "
                          "products: max defect {:.2e}".format(
                              nchk0 + nchk1, max(chk0, chk1)),
                   "pass": bool(max(chk0, chk1) < 1e-8)})
    fin = all(bool(np.all(np.isfinite([v["accept"]] + v["W"]
                                       + [x for row in v["C"]
                                          for x in row])))
              for v in vrows)
    checks.append({"name": "correlator table finite over {} betas "
                          "(the mass/branch verdict is physics, carried "
                          "by the branch)".format(len(vrows)),
                   "pass": bool(fin)})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok=ok)
    print("[V300 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V300_cone_measure_selection", "tier": 2,
           "property": "does the V289 geometry measure concentrate a "
                       "light cone (line-correlator anisotropy)",
           "residual": True, "programme": "V300",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v300_cone_measure_selection.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V300 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
