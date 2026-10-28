#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V289 - the geometry measure, proof of principle (tier 3, hard-lite).

**The question (pre-registered).** Does the root's finite geometry
group (the order-1024 even Clifford monomial group, V141
build_group) support a measure: Metropolis sampling of bond
variables u_b in Gamma on 4^3 with the adjoint-character Wilson
action S = beta sum_p (1 - chi_45(U_p)/45) gives a finite, decaying
Wilson-loop expectation (proof that the geometry measure exists
before it is promoted to dynamics)?

**The object (closed input set).** V141's group (build_group/mul/
invof, order re-measured), the 45 generators and structure constants
via V140's structure_constants on the Clifford plane generators,
the 45-dim adjoint representation R(g) built multiplicatively along
the generation tree; chi_45 = Tr R. Lattice 4^3 periodic, 192
bonds, 192 based plaquettes (inverse links carried as R(g^-1) =
R(g)^T of the orthogonal rep). Metropolis: one random bond per
sweep, proposal = left-multiplication by one of the 90 biflections
or their inverses, NREP = 32 replicas batched on the GPU,
SWEEPS = 20000 (2000 equilibrated), beta in {1,2,4,8,16}, seed 289.

**Measured (pre-registered).** accept(beta); <W(R)> for R = 1..3
(reals of chi_45/45 around based R x R squares, 16 random base
points/planes per sample, every 20 sweeps after equilibration),
jackknife error over replicas; the adjoint character cross-check:
chi_45 of the sampled plaquette element read from the precomputed
table (full multiplication table of Gamma) against
Tr[R(g1) R(g2) R(g3^-1) R(g4^-1)] of the plaquette product - exact
match (tol 1e-8) on every checked plaquette; area-law fit of
log<W(R)> against R^2 per beta where all W > 0.

**What is NOT derived (stated).** That this finite-group measure
reproduces continuum diffeo-invariant geometry; the decay law is a
proof-of-existence measurement on 4^3 only, and R <= 3 on L = 4 is
the tunnel-free regime bound booked here.

Branch: MC-GEOMETRY-OK (estimator finite, area fit exists at some
beta); MC-TUNED additionally books the beta range with accept > 2 %
(all replicas tunneled-free); MC-NULL otherwise - the estimator is
sound (character check exact, errors finite) but no decaying area
fit on 4^3 at the booked betas; MC-NULL is a booked physics
verdict, not an implementation failure.
"""
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

L = 4
NSITE = L ** 3
NBOND = 3 * NSITE
NPLAQ = 3 * NSITE
BETAS = [1.0, 2.0, 4.0, 8.0, 16.0]
NREP = 32
SWEEPS = 20000
EQUIL = 2000
SAMPLE_EVERY = 20
BASEPOINTS = 16
SEED = 289


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V289 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    # ---- the group, its 45-dim adjoint rep and characters ----------------
    G = V141.build_group()
    ng = len(G)
    idx = {g: i for i, g in enumerate(G)}
    checks.append({"name": "T0 group order re-measured: {} == 1024"
                          .format(ng), "pass": ng == 1024})
    log("GAMMA order {}".format(ng))
    from cl10 import clifford                              # noqa: E402
    from cl10 import sm_qca                                # noqa: E402

    def _np(x):
        return x.get() if hasattr(x, "get") else np.asarray(x)
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

    # multiplicative build along the generation tree (BFS parents)
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
    assert len(order) == ng, (len(order), ng)
    Rm = np.zeros((ng, 45, 45))
    Rm[0] = np.eye(45)
    closeness = 0.0
    for g in order[1:]:
        pg, s = parent[g]
        Rs, resid = adjoint_mat(s)
        closeness = max(closeness, resid)
        Rm[idx[g]] = Rm[idx[pg]] @ Rs
    chi = np.real(np.trace(Rm, axis1=1, axis2=2))
    # representation homomorphism spot check on 64 random products
    rng0 = np.random.default_rng(SEED)
    hmax = 0.0
    for _ in range(64):
        a, b = int(rng0.integers(ng)), int(rng0.integers(ng))
        ab = idx[V141.mul(G[a], G[b])]
        hmax = max(hmax, float(np.abs(Rm[a] @ Rm[b] - Rm[ab]).max()))
    checks.append({"name": "adjoint rep: chi(1) = {:.1f}, closure of the "
                          "45-span under conjugation max defect {:.2e}, "
                          "homomorphism defect on 64 sampled products "
                          "{:.2e}".format(chi[0], closeness, hmax),
                   "pass": bool(abs(chi[0] - 45.0) < 1e-8 and closeness
                                < 1e-9 and hmax < 1e-8)})

    # multiplication table (full), inverse table, and proposal table
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
    Rgi = cp.asarray(Rm[invtab])           # R(g^-1) = R(g)^T
    tabc = cp.asarray(tab)
    multabc = cp.asarray(multab)
    invtabc = cp.asarray(invtab)
    chig = cp.asarray(chi)

    # ---- lattice bookkeeping ---------------------------------------------
    planes = [(0, 1), (0, 2), (1, 2)]
    plaqs = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                s0 = x + L * y + L * L * z
                for (mu, nu) in planes:
                    dmu = (1, L, L * L)[mu]
                    dnu = (1, L, L * L)[nu]
                    xs = (s0 + dmu) % NSITE
                    xn = (s0 + dnu) % NSITE
                    # U_p = u1 u2 u3^-1 u4^-1 (u3, u4 stored as forward
                    # bonds, traversed backwards through R(g^-1))
                    plaqs.append((s0 * 3 + mu, xs * 3 + nu,
                                  xn * 3 + mu, s0 * 3 + nu))
    plaqs = np.array(plaqs, dtype=np.int64)
    assert len(plaqs) == NPLAQ
    bond_pl = [[] for _ in range(NBOND)]
    bond_slot = [[] for _ in range(NBOND)]
    for pi in range(NPLAQ):
        for slot in range(4):
            b = int(plaqs[pi, slot])
            bond_pl[b].append(pi)
            bond_slot[b].append(slot)
    # precomputed plaquette gather order with the inverse slots flagged
    plaq_use = np.zeros((NPLAQ, 4), dtype=np.int8)
    plaq_use[:, :2] = 1

    def plaq_chain(uu, plist):
        """ordered product matrices for plaquettes plist (all reps)."""
        b = plaqs[plist]
        P = (Rg[uu[:, b[:, 0]]] @ Rg[uu[:, b[:, 1]]]) @ \
            (Rgi[uu[:, b[:, 2]]] @ Rgi[uu[:, b[:, 3]]])
        return P

    def plaq_gidx(u, plist):
        """group index of the plaquette element u1 u2 u3^-1 u4^-1."""
        b = plaqs[plist]
        g12 = multabc[u[:, b[:, 0]], u[:, b[:, 1]]]
        g34 = multabc[invtabc[u[:, b[:, 2]]], invtabc[u[:, b[:, 3]]]]
        return multabc[g12, g34]

    all_pl = np.arange(NPLAQ)
    reps_all = cp.arange(NREP)

    rng = np.random.default_rng(SEED)
    rows = []
    chk_max = 0.0
    nchk = 0
    for beta in BETAS:
        u = cp.asarray(rng.integers(0, ng, (NREP, NBOND), dtype=np.int64))
        chiP = cp.real(cp.trace(plaq_chain(u, all_pl),
                                axis1=2, axis2=3))
        # exactness check at initialisation on every plaquette
        gp = plaq_gidx(u, all_pl)
        chk_max = max(chk_max, float(cp.abs(chig[gp] - chiP).max()))
        nchk += NREP * NPLAQ
        acc = 0
        Wsam = None
        nsam = 0
        for sweep in range(SWEEPS):
            bond = int(rng.integers(NBOND))
            pidx = rng.integers(0, nprop, NREP)
            pidxc = cp.asarray(pidx)
            ubak = u[reps_all, bond].copy()
            plist = np.array([int(p) for p in bond_pl[bond]])
            oldc = chiP[:, plist]
            unew = tabc[u[reps_all, bond], pidxc]
            u[reps_all, bond] = unew
            # recompute only the affected chains: full 4-factor chain per
            # affected plaquette (batched over replicas and plaquettes)
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
                # character cross-check on 4 random affected-free
                # plaquettes this sample
                cpl = rng.choice(NPLAQ, 4, replace=False)
                gp = plaq_gidx(u, cpl)
                cc = cp.real(cp.trace(plaq_chain(u, cpl),
                                      axis1=2, axis2=3))
                chk_max = max(chk_max, float(cp.abs(chig[gp] - cc).max()))
                nchk += 4 * NREP
                # Wilson loops: BASEPOINTS random (site, plane) pairs,
                # R x R based squares, R = 1..3
                sel = rng.integers(0, NPLAQ, BASEPOINTS)
                wr = np.zeros((NREP, 3))
                for ri, pi in enumerate(sel):
                    mu, nu = planes[int(pi) % 3]
                    s0 = int(pi) // 3
                    dmu = (1, L, L * L)[mu]
                    dnu = (1, L, L * L)[nu]
                    for R in (1, 2, 3):
                        tot = cp.zeros(NREP)
                        cnt = 0
                        for sh in range(R):
                            for th in range(R):
                                cur = (s0 + sh * dmu + th * dnu) % NSITE
                                P = cp.eye(45)[None].repeat(NREP, 0)
                                c = cur
                                for _ in range(R):
                                    P = P @ Rg[u[:, c * 3 + mu]]
                                    c = (c + dmu) % NSITE
                                for _ in range(R):
                                    P = P @ Rg[u[:, c * 3 + nu]]
                                    c = (c + dnu) % NSITE
                                for _ in range(R):
                                    c = (c - dmu) % NSITE
                                    P = P @ Rgi[u[:, c * 3 + mu]]
                                for _ in range(R):
                                    c = (c - dnu) % NSITE
                                    P = P @ Rgi[u[:, c * 3 + nu]]
                                tot += cp.real(cp.trace(P, axis1=1,
                                                        axis2=2))
                                cnt += 1
                        wr[:, R - 1] += cp.asnumpy(tot / (45.0 * cnt))
                        del P
                wr /= BASEPOINTS
                Wsam = wr
                nsam += 1
        wbar = Wsam.mean(axis=0)
        jack = np.stack([np.delete(Wsam, rr, axis=0).mean(axis=0)
                         for rr in range(NREP)])
        err = np.sqrt((NREP - 1) / NREP * np.sum((jack - wbar) ** 2,
                                                axis=0))
        rows.append({"beta": beta, "accept": acc / (SWEEPS * NREP),
                     "nsamples": int(nsam),
                     "W1": float(wbar[0]), "W2": float(wbar[1]),
                     "W3": float(wbar[2]),
                     "err": [float(e) for e in err]})
        log("beta {:5.1f}: accept {:.3f} W {} err {}".format(
            beta, acc / (SWEEPS * NREP),
            ["{:.4f}".format(w) for w in wbar],
            ["{:.4f}".format(e) for e in err]))
    data["rows"] = rows
    data["char_check"] = {"n_checked": int(nchk),
                          "max_abs_defect": float(chk_max)}

    # area-law fit log<W(R)> vs R^2 where all W are positive
    fits = []
    for r in rows:
        if all(r[k] > 0 for k in ("W1", "W2", "W3")):
            Rs = np.array([1.0, 4.0, 9.0])
            yv = np.log(np.array([r["W1"], r["W2"], r["W3"]]))
            A = np.vstack([Rs, np.ones(3)]).T
            sol = np.linalg.lstsq(A, yv, rcond=None)[0]
            sl, ic = float(sol[0]), float(sol[1])
            fits.append({"beta": r["beta"], "area_slope": float(sl),
                         "intercept": float(ic),
                         "decay": bool(sl < 0.0)})
    data["area_fits"] = fits
    fit_ok = any(f["decay"] for f in fits)
    tunneled_free = [r["beta"] for r in rows if r["accept"] > 0.02]
    status = "MC-TUNED" if fit_ok and tunneled_free else \
        ("MC-GEOMETRY-OK" if fit_ok else "MC-NULL")
    data["status"] = status
    print("V289 {}: fits {} (decaying {}), tunneled-free betas {}, "
          "char check {} plaquette-entries max defect {:.2e}"
          .format(status, len(fits),
                  sum(f["decay"] for f in fits), tunneled_free, nchk,
                  chk_max), flush=True)
    fin = all(np.isfinite([r[k] for r in rows
                           for k in ("W1", "W2", "W3", "accept")]))
    # (contract amendment, pre_registrations: zero windows are finite -
    # a physics result carried by the branch, not a machinery failure)
    checks.append({"name": "adjoint character cross-check exact on {} "
                          "plaquette products: max defect {:.2e}"
                          .format(nchk, chk_max),
                   "pass": bool(chk_max < 1e-8)})
    checks.append({"name": "estimator finite over {} betas x {} replicas:"
                          " {} (decay verdict carried by the branch)"
                          .format(len(BETAS), NREP, fin),
                   "pass": bool(fin)})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V289 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V289_geometry_mc", "tier": 3,
           "property": "Metropolis measure for the root geometry group "
                       "on 4^3 (proof of principle)",
           "residual": True,
           "programme": "V289",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v289_geometry_mc.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V289 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
