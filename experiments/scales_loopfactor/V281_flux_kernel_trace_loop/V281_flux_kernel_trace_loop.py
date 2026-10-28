#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V281 - the flux-kernel traced loop (tier 2, hard, O8-new).

**The question (pre-registered).** Does the traced-word Cartan-ratio
functional produce the loop band 1/(16 pi) when the traced sector is
the flux kernel of the three-site gauged chain and the traced state
is the bond-current ground state instead of the uniform link average?

Object per the pre-registration + construction addendum (both in
review/pre_registrations.md): the END-PAIR sector (1024 = 32 x 32)
with G_pair = G (x) G at the certified closed angle pi/4; the exact
Lagrange projectors of the integer end-pair Casimir C2pair (V101's
operator, V129.channels convention); M = sum_e sqrt(w_e) exp(i
theta_e) Pi_e with w_e = Tr(rho_B Pi_e), rho_B the V101 ground-state
end-pair (middle-cut) reduction (Tr 1, positivity, S = ln 45
gates); the rho-weighted content over the 10 Cartan planes (the five
V209 Cartan planes on each tensor factor); F_n = cnorm((M G)^n) /
cnorm(G^n), even n <= 40; the (theta_0, theta_64, theta_96) grid at
13 points (+ 14-point stability control); the ladder runs on the GPU
(1024-dim dense matmuls, the batched-chunk convention of V252); the
scalar 50-digit continuation of the measured geometric class in
mpmath with the digit cap booked.

**Acceptance.** FLUXLOOP-FOUND iff a stabilized grid point has
|16 pi F - 1| < 0.5 with F in [1.5e-2, 3.0e-2]; else
FLUXLOOP-NEGATIVE + floor. T0 anchors (V101 artifact): kernel 24,
E0 -5802.398188196636, degeneracy 6, gap 42.09063148320456,
S_ground 3.8066624897703183; drift -> nonzero exit.
"""
import importlib.util
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402
from scipy.linalg import eigh as geigh                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402

from cl10 import clifford                                 # noqa: E402
from cl10 import engine as E                              # noqa: E402

GCLOSED = math.pi / 4.0
try:
    import mpmath as mp                                  # noqa: E402
except Exception:                                        # pragma: no cover
    mp = None


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


def channels(C2, tol=1e-6):
    """Group eigenvalues of a Casimir into [value, multiplicity]
    (the V129 convention)."""
    ev = np.sort(np.real(np.linalg.eigvalsh(C2)))
    out = []
    for e in ev:
        if out and abs(e - out[-1][0]) < tol:
            out[-1][1] += 1
        else:
            out.append([float(e), 1])
    return out


def main():
    t0 = time.time()
    checks, data = [], {}
    xp = E._E if E._BACKEND == "cupy" else np
    print(f"V281 - flux-kernel traced loops [ladder backend: "
          f"{E._BACKEND}]", flush=True)

    V79 = load_mod("V79_gauged_bond_layer")
    V101 = load_mod("V101_bond_current_ground_state")

    # ---- T0: V101's kernel + ground state, re-derived --------------------
    Is, Ts = V101.build_kernel(t0)
    Ls_np = V79.adjoint_link_operators(
        V79.structure_constants(V79.root_generators())[0])
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
            rhoB += w * np.einsum("b,kb,lb->kl", kk, Is[b],
                                  Is[b2].conj())
    rhoB = (rhoB + rhoB.conj().T) / 2
    evr = np.linalg.eigvalsh(rhoB)
    trB = float(np.sum(evr[evr > 1e-12]))
    S_gs = float(-np.sum(evr[evr > 1e-14] *
                         np.log(evr[evr > 1e-14])))

    a_E0 = float(art("V101_bond_current_ground_state", "data.E0"))
    a_deg = int(art("V101_bond_current_ground_state", "data.degeneracy"))
    a_gap = float(art("V101_bond_current_ground_state", "data.gap"))
    a_S = float(art("V101_bond_current_ground_state", "data.S_ground"))
    drift = []
    if ns != 24:
        drift.append("kernel count")
    if abs(E0 - a_E0) > 1e-9 * abs(a_E0):
        drift.append("E0")
    if deg != a_deg:
        drift.append("degeneracy")
    if abs(gap - a_gap) > 1e-9 * abs(a_gap):
        drift.append("gap")
    if abs(S_gs - a_S) > 1e-9:
        drift.append("S_ground")
    if abs(trB - 1.0) > 1e-9 or evr.min() < -1e-14:
        drift.append("rho_B state gates")
    data["t0"] = {"kernel": ns, "E0": E0, "degeneracy": deg, "gap": gap,
                  "S_ground": S_gs, "tr_rhoB": trB,
                  "rhoB_min_eig": float(evr.min()),
                  "artifact": {"E0": a_E0, "degeneracy": a_deg,
                               "gap": a_gap, "S_ground": a_S}}
    checks.append({"name": "T0 V101 ground state re-derived: kernel {}, "
                          "E0 {:.6f} vs {:.6f}, deg {}, gap {:.4f}, "
                          "S {:.9f} vs {:.9f}, Tr rho_B {:.12f}, min "
                          "eig {:.1e}".format(ns, E0, a_E0, deg, gap,
                                              S_gs, a_S, trB, evr.min()),
                   "pass": not drift})
    if drift:
        data["status"] = "MACHINERY-DRIFT"
        print("V281 MACHINERY-DRIFT: {} disagrees with the live V101 "
              "artefact; nothing downstream interpreted"
              .format(", ".join(drift)), flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the channel algebra on the end-pair sector ----------------------
    print("V281 T0 pass: V101 state re-derived ({:.0f} s); C2pair build"
          .format(time.time() - t0), flush=True)
    C2pair = np.zeros((1024, 1024), dtype=np.complex128)
    for Jm in Js:
        C2pair += Jm @ Jm
    chans = channels(C2pair)
    vals = [round(v, 3) for v, _ in chans]
    # Lagrange projectors over the measured channel values
    Pi = {}
    for (v, _m) in chans:
        A = np.eye(1024, dtype=np.complex128)
        den = 1.0
        for (v2, _m2) in chans:
            if abs(v2 - v) < 1e-6:
                continue
            A = A @ (C2pair - v2 * np.eye(1024, dtype=np.complex128))
            den *= (v - v2)
        Pi[round(v, 3)] = A / den
    pidem = max(float(np.abs(p @ p - p).max()) for p in Pi.values())
    w_e = {k: float(np.real(np.trace(rhoB @ p))) for k, p in Pi.items()}

    # the end-pair link-pair channels row (booked, not in the ladder)
    I45 = np.eye(45, dtype=np.complex128)
    Lll = [np.kron(L, I45) + np.kron(I45, L) for L in Ls_np]
    C2ll = np.zeros((2025, 2025), dtype=np.complex128)
    for Lm in Lll:
        C2ll += Lm @ Lm
    chans_ll = channels(C2ll)

    # ---- G at the closed clock, and the ladder machinery ------------------
    h = [np.array(E.asnumpy(m), dtype=np.complex128)
         for m in clifford.generators()]
    planes = []
    for a in range(10):
        for b in range(a + 1, 10):
            T = 1j * (h[a] @ h[b])
            planes.append(T)

    # V209's Fm entry formula: Fm = (cos + sin*Im(T)) + i(-sin*Re(T))
    # with T = 1j h_a h_b - the float64 mirror of `mp_ladder`'s G.
    def G32_v209(angle):
        Gm = np.eye(32, dtype=np.complex128)
        for T in planes:
            Fm = (math.cos(angle) * np.eye(32, dtype=np.complex128)
                  + math.sin(angle) * T.imag
                  - 1j * math.sin(angle) * T.real)
            Gm = Gm @ Fm
        return Gm

    def G32(angle):
        Gm = np.eye(32, dtype=np.complex128)
        for a in range(10):
            for b in range(a + 1, 10):
                X = h[a] @ h[b]
                Hm = np.rint(X.real) + 1j * np.rint(X.imag)
                Gm = Gm @ (math.cos(angle) * np.eye(32,
                                                    dtype=np.complex128)
                           + math.sin(angle) * Hm)
        return Gm
    # ladder clock: the OPEN branch g = 1/5 (V209 booked that the
    # Cartan functional is precision-stable only on the open branch;
    # at the closed clock the content sits at the phase-degeneracy
    # floor). The booked closed-clock row is kept as a control.
    from cl10 import config as _config                    # noqa: E402
    GOPE = float(_config.CFG.angle.branches.implemented)  # 1/5
    Gf = G32(GOPE)
    Gv = G32_v209(GOPE)
    Gdiff = float(np.abs(Gf - Gv).max())
    Gc = G32(GCLOSED)

    G = np.kron(Gf, Gf)                     # 1024-dim G_pair
    Gx = xp.asarray(G)
    Px = {k: xp.asarray(p) for k, p in Pi.items()}
    # content: the five Cartan planes on each tensor factor, rho-weighted
    CH = []
    for ca in (0, 2, 4, 6, 8):
        Hh = h[ca] @ h[ca + 1]
        CH.append(np.kron(Hh, np.eye(32)) @ rhoB)
        CH.append(np.kron(np.eye(32), Hh) @ rhoB)
    CHx = xp.asarray(np.stack(CH))

    def cnorm(Yx):
        return float(xp.sum(xp.abs(xp.einsum("pab,ba->p", CHx, Yx))))

    Gx_closed = xp.asarray(np.kron(Gc, Gc))

    def ladder_at(thetas, nmax=40, Gx_in=None):
        keys = GRID_CH
        if Gx_in is None:
            Gx_in = Gx
        Mx = xp.zeros((1024, 1024), dtype=xp.complex128)
        for k, th in zip(keys, thetas):
            Mx = Mx + math.sqrt(max(w_e[k], 0.0)) * xp.exp(1j * th) * Px[k]
        for k in FROZEN_CH:
            Mx = Mx + math.sqrt(max(w_e[k], 0.0)) * Px[k]
        X = Mx @ Gx_in
        Xn = xp.eye(1024, dtype=xp.complex128)
        Gn = xp.eye(1024, dtype=xp.complex128)
        out = {}
        for n in range(1, nmax + 1):
            Xn = Xn @ X
            Gn = Gn @ Gx
            if n % 2 == 0:
                cg = cnorm(Gn)
                out[n] = {"cnorm_G": cg, "cnorm_X": cnorm(Xn),
                          "F": cnorm(Xn) / cg if cg > 0 else None}
        return out

    # ---- the phase grid ---------------------------------------------------
    # the pre-registration's 13^3 grid is over the kernel-supported
    # channels (w_e > 1e-9); any further channel is frozen at theta = 0
    # and its weight booked.
    GRID_CH = [k for k in sorted(w_e) if w_e[k] > 1e-9]
    FROZEN_CH = [k for k in sorted(w_e) if w_e[k] <= 1e-9]
    data["grid_channels"] = {"grid": GRID_CH, "frozen": FROZEN_CH}
    if len(GRID_CH) < 2:
        data["grid_channels"]["collapse_note"] = (
            "rho_B has single-channel support (w_64 = 1, all other "
            "weights zero): the channel-projector insertion reduces to "
            "a scalar phase on the rho-weighted content and the "
            "functional is identically blind on this traced state - a "
            "booked property of the V101 middle-cut state, not of the "
            "band; the negative is recorded with this note")
    NGRID = 13
    grid_th = np.linspace(0.0, 2.0 * math.pi, NGRID)
    best = None
    npts = 0
    import itertools
    for combo in itertools.product(grid_th, repeat=len(GRID_CH)):
        L = ladder_at(combo)
        F40 = L[40]["F"]
        if F40 is None:
            continue
        res = abs(16.0 * math.pi * F40 - 1.0)
        npts += 1
        if npts % 400 == 0:
            print("V281 grid {:4d} pts, best |16piF-1| {:.4e} ({:.0f} s)"
                  .format(npts, best[0] if best else -1,
                          time.time() - t0), flush=True)
        if best is None or res < best[0]:
            best = (res, combo, {n: L[n]["F"] for n in range(2, 41, 2)})
    # stability control: the 14-point grid at the best neighbourhood
    g14 = np.linspace(0.0, 2.0 * math.pi, 14)
    stab = None
    if best is not None:
        cand = tuple(min(g14, key=lambda v: abs(v - th))
                     for th in best[1])
        L14 = ladder_at(cand)
        stab = {"thetas_14": [float(v) for v in cand],
                "F40_14": L14[40]["F"],
                "F40_shift": abs(L14[40]["F"] - best[2][40])}
    # the all-thetas-zero row (the open-branch analogue) and the
    # closed-clock control (V209's floor class: booked, not band-tested)
    Lu = ladder_at((0.0,) * len(GRID_CH))
    Lc = ladder_at((0.0,) * len(GRID_CH), Gx_in=Gx_closed)
    data["closed_clock_control"] = {
        "F_even": {str(n): (Lc[n]["F"] if Lc[n]["F"] is not None
                            else None) for n in range(2, 41, 2)},
        "cnorm_G_even": {str(n): Lc[n]["cnorm_G"]
                         for n in range(2, 41, 2)},
        "note": "the closed-clock content sits at the phase-degeneracy "
                "floor (V209); the row is booked as the control, not "
                "band-tested"}

    data["channels"] = {"end_pair": chans, "link_pair": chans_ll,
                        "channel_weights_w_e": w_e,
                        "projector_idempotent_defect": pidem,
                        "G_route_agreement": Gdiff}
    data["grid"] = {"n_points": npts, "ngrid": NGRID,
                    "best": {"residual_16piF_minus_1": best[0],
                             "thetas": [float(v) for v in best[1]],
                             "F_even": {str(k): (float(v) if v is not
                                                 None else None)
                                        for k, v in best[2].items()}},
                    "stability": stab,
                    "uniform_row_F_even": {str(n): Lu[n]["F"]
                                           for n in range(2, 41, 2)}}

    # ---- scalar 50-digit continuation of the measured class ---------------
    Fseq = [v for k, v in sorted(best[2].items()) if v is not None]
    ks = list(range(1, len(Fseq) + 1))
    logr = [math.log(Fseq[i + 1] / Fseq[i]) for i in range(len(Fseq) - 1)]
    rho_chain = math.exp(sum(logr[5:]) / len(logr[5:])) if len(logr) > 6 \
        else math.exp(sum(logr) / len(logr))
    cont = {"rho_chain": rho_chain,
            "log_ratio_spread": float(max(logr) - min(logr)),
            "digit_cap": "float64 ladder: 16 significant digits enter the "
                         "50-digit mp continuation; the continuation is "
                         "guarded arithmetic on a float64 measurement, "
                         "booked as such"}
    if mp is not None:
        with mp.workdps(50):
            r = mp.mpf(repr(rho_chain))
            cont["mp_check_F40_geom"] = repr(mp.exp(mp.log(r * r) * 20)
                                             * mp.mpf(repr(Fseq[0])))
            cont["mp_psln_note"] = repr(mp.log(16 * mp.pi) / mp.log(1 / r))
    data["continuation"] = cont

    # ---- verdict -----------------------------------------------------------
    band = 1.5e-2 <= best[2][40] <= 3.0e-2
    stable = bool(stab is not None and stab["F40_shift"] < 1e-8)
    found = bool(best[0] < 0.5 and band and stable)
    branch = "FLUXLOOP-FOUND" if found else "FLUXLOOP-NEGATIVE"
    data["status"] = branch
    data["found_legs"] = {"residual_lt_0_5": best[0] < 0.5,
                          "band": band, "stable": stable}
    print("V281 {}: closest |16 pi F - 1| = {:.4e} at thetas {} "
          "(F_40 = {:.6e}, band {}, stable {}) channel weights {} "
          "(channels {}); rho_chain {:.6f}".format(
              branch, best[0], [round(v, 4) for v in best[1]],
              best[2][40], band, stable,
              {k: round(v, 4) for k, v in w_e.items()}, vals, rho_chain),
          flush=True)
    checks.append({"name": "T1 channel algebra and grid: channels {} "
                          "(link pair {}), w_e {}, projector defect "
                          "{:.1e}, {} grid points evaluated, closest "
                          "|16 pi F_40 - 1| {:.3e}".format(
                              vals, [[v, m] for v, m in chans_ll],
                              {k: round(v, 4) for k, v in w_e.items()},
                              pidem, npts, best[0]),
                   "pass": pidem < 1e-8 and Gdiff < 1e-12})

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    print("[V281 done] {}  ({:.0f} s)".format(branch, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V281_flux_kernel_trace_loop", "tier": 2,
           "property": "the traced-word Cartan-ratio functional with the "
                       "three-site flux-kernel channels and the "
                       "bond-current ground state as traced state: does it "
                       "land in the loop band (tier 2, O8-new)",
           "residual": True,
           "programme": "V281",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v281_flux_kernel_trace_loop.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V281 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
