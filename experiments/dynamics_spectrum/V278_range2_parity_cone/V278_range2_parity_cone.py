#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V278 - range-2 circuits at the parity-mixed node (tier 2, hard, O10).

**The question (pre-registered).** Does the circuit class enlarged by
range-2 second-neighbour swap rotations reach the parity-mixed
A1g (+) T1u conical node of V191 inside the engine's own six-stage
product?

**The object (closed input set).** the V277 decoration family (the
V277 module imported verbatim: same triplet, same closed clock, same
18 decoration parameters) extended on the bond factor by two orbit
rotations per stage:
  axial: exp(i a_t G_{j_t}(q)) with (the pre-registration's sum_j
    generator entered stage-axis-wise, so every factor is exactly
    local: G_j(q) = sum_{e: bit_j(e)=0} [e^{+i q_j} |e (+) (1<<j)><e|
    + h.c.] is a disjoint perfect matching of the eight sub-sites, so
    the exponential is a 2x2 rotation per pair - Fourier degree 1,
    real-space reach 2 per stage, the same pair set and phase
    convention as the engine's inter-block swaps - the
    (plus/minus 2,0,0)-type orbit);
  diagonal: exp(i d_t D) with D = (1/2) sum_{i<j} sum_e
    [|e (+) (1<<i) (+) (1<<j)><e| + h.c.] - 12 pairs, q-independent,
    the (plus/minus 1, plus/minus 1, 0)-type intra-block double flips.
30 parameters in [0, 2 g10], g10 = 0.2 from cl10 config. The stage
factor is kron(S_t(pi/4) exp(i a_t G_{j_t}(q)) exp(i d_t D),
U16 W_t(q)); unitarity exact per factor.
Nelder-Mead (V186 tune polish), 150 starts (seed 278: start 1 = the
published decoration with a_t = d_t = 0; 149 uniform box draws),
maxfev 300, objective = the GAMMA gap (minimize).

**Measured.** T0: the V277 artefact's bare top-64 degree and uniform
triplet degree re-derived (1e-6); the published step's GAMMA
equivariance defect vs the live V192 artefact (1e-9); the 64+64
census. At the best point: the gap-subspace dimension m and size,
the O_h content multiplicities and A1g/T1u projector ranks (the V191
character machinery, W48 from `supercell.oh_subsite_matrices` - the
representation of the cubic group on THIS step's sub-sites), the
V186 tangent rho_comm and 179-Fibonacci speeds (v_max, A_top), the
degree of the gap multiplet, and the GAMMA equivariance defect of
the extended product.

**Acceptance.** CONE-RANGE2-FOUND iff at the best point: GAMMA gap
< 1e-10, the gap subspace (m = 4) has A1g rank >= 1 and T1u rank
>= 3 (parity-mixed pair), A_top < 0.1, rho_comm > 0.5, v_max <= 1.
Else CONE-RANGE2-ABSENT with the best (gap, A_top, rho_comm) triple
and the binding leg.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("CL10_BACKEND", "numpy")  # launcher may force cupy

import numpy as np                                        # noqa: E402
from scipy.linalg import expm                             # noqa: E402
from scipy.optimize import minimize                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
sys.path.insert(0, exp("V277_bond_decoration_cone_scan"))

import V277_bond_decoration_cone_scan as V277             # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import parpool                                  # noqa: E402
from cl10 import supercell as SC                          # noqa: E402

CFG = _config.CFG
G10 = float(CFG.angle.branches.implemented)
GCLOSED = math.pi / 4.0
NSTARTS, MAXFEV, SEED = 150, 300, 278
BOX = 2.0 * G10

# ---------------------------------------------------------------------------
# the range-2 bond generators (module level: inherited by the pool)
# ---------------------------------------------------------------------------
def G_ax(j):
    """G_j(q) at unit amplitude: the axial (plus/minus 2,0,0) matching,
    8x8, Fourier degree 1 (phase e^{+/- i q_j})."""
    def G(qj):
        M = np.zeros((8, 8), dtype=complex)
        bit = 1 << j
        for e in range(8):
            if (e & bit) == 0:
                o = e ^ bit
                M[o, e] = np.exp(1j * qj)
                M[e, o] = np.conj(M[o, e])
        return M
    return G


def D_diag():
    """The diagonal (plus/minus 1, plus/minus 1, 0) double-flip generator,
    12 pairs, q-independent."""
    M = np.zeros((8, 8), dtype=complex)
    for i in range(3):
        for j in range(i + 1, 3):
            for e in range(8):
                o = e ^ (1 << i) ^ (1 << j)
                M[o, e] += 0.5
                M[e, o] += 0.5
    return M


D_STATIC = D_diag()


def bloch30(q, p):
    """p = (dec18 flat, a0..a5, d0..d5): the decorated, range-2 extended
    six-stage step at the closed clock."""
    dec = [p[3 * t:3 * t + 3] for t in range(6)]
    aa = p[18:24]
    dd = p[24:30]
    M = np.eye(128, dtype=complex)
    for t, (even, j) in enumerate(SC.STAGES):
        S = SC.stage_matrix(even, j, q[j], GCLOSED, GCLOSED)
        B = S @ expm(1j * aa[t] * G_ax(j)(q[j])) @ expm(1j * dd[t]
                                                        * D_STATIC)
        W = V277.twist(q, tuple(dec[t]))
        M = np.kron(B, V277.U16 @ W) @ M
    return M


def gap_gamma(p):
    ev = np.linalg.eigvals(bloch30((0.0, 0.0, 0.0), p))
    return float(np.min(np.abs(np.angle(ev))))


def _start_worker(i):
    rng = np.random.default_rng(SEED + i)
    if i == 0:
        x0 = np.concatenate([np.full(18, G10), np.zeros(12)])
    else:
        x0 = np.concatenate([rng.uniform(0.0, BOX, 18),
                             rng.uniform(0.0, BOX, 12)])
    simplex = [x0] + [x0 + 0.2 * e for e in np.eye(30)]
    res = minimize(gap_gamma, x0, method="Nelder-Mead",
                   options={"maxfev": MAXFEV, "fatol": 1e-12,
                            "xatol": 1e-10,
                            "initial_simplex": np.array(simplex)})
    x = np.clip(res.x, 0.0, BOX)
    return {"start": i, "x": [float(v) for v in x], "gap": gap_gamma(x),
            "fev": int(res.nfev)}


# ---------------------------------------------------------------------------
# O_h content on the 128 space (V191 character machinery, W48 from
# supercell.oh_subsite_matrices - the cubic group on THIS step's sites)
# ---------------------------------------------------------------------------
_G48, _W48 = [], []
for (perm, signs), W in SC.oh_subsite_matrices():
    R3 = np.zeros((3, 3))
    for i, (p_, s) in enumerate(zip(perm, signs)):
        R3[i, p_] = s
    _G48.append(R3)
    _W48.append(np.kron(W, np.eye(16, dtype=complex)))


def class_of(R3):
    tr, det = round(np.trace(R3)), int(round(np.linalg.det(R3)))
    if det > 0:
        if tr == 3:
            return "E"
        if tr == 0:
            return "C3"
        if tr == 1:
            return "C4"
        w, vecs = np.linalg.eigh(R3)
        ax = vecs[:, int(np.argmax(w))]
        return "C2p" if sum(abs(x) > 1e-9 for x in ax) == 1 else "C2"
    if tr == -3:
        return "i"
    if tr == 0:
        return "S6"
    if tr == -1:
        return "S4"
    w, vecs = np.linalg.eigh(R3)
    nm = vecs[:, int(np.argmin(w))]
    return "sh" if sum(abs(x) > 1e-9 for x in nm) == 1 else "sd"


_CLS = [class_of(R) for R in _G48]
CHAR = {
    "A1g": dict(E=1, C3=1, C2=1, C2p=1, C4=1, i=1, S6=1, sd=1, sh=1, S4=1),
    "A2g": dict(E=1, C3=1, C2=-1, C2p=1, C4=-1, i=1, S6=1, sd=-1, sh=1, S4=-1),
    "Eg": dict(E=2, C3=-1, C2=0, C2p=2, C4=0, i=2, S6=-1, sd=0, sh=2, S4=0),
    "T1g": dict(E=3, C3=0, C2=-1, C2p=-1, C4=1, i=3, S6=0, sd=-1, sh=-1, S4=1),
    "T2g": dict(E=3, C3=0, C2=1, C2p=-1, C4=-1, i=3, S6=0, sd=1, sh=-1, S4=-1),
    "A1u": dict(E=1, C3=1, C2=1, C2p=1, C4=1, i=-1, S6=-1, sd=-1, sh=-1, S4=-1),
    "A2u": dict(E=1, C3=1, C2=-1, C2p=1, C4=-1, i=-1, S6=-1, sd=1, sh=-1, S4=1),
    "Eu": dict(E=2, C3=-1, C2=0, C2p=2, C4=0, i=-2, S6=1, sd=0, sh=-2, S4=0),
    "T1u": dict(E=3, C3=0, C2=-1, C2p=-1, C4=1, i=-3, S6=0, sd=1, sh=1, S4=-1),
    "T2u": dict(E=3, C3=0, C2=1, C2p=-1, C4=-1, i=-3, S6=0, sd=-1, sh=1, S4=1),
}


def content_of(subspace_Q, names=tuple(CHAR)):
    """O_h multiplicities of the subspace spanned by the columns of Q
    (128 x m): mult_n = (d_n/48) sum_g chi_n(g)^* Tr(Q^dag W(g) Q)."""
    out = {}
    for n in names:
        d = CHAR[n]["E"]
        s = 0.0 + 0.0j
        for k in range(48):
            WQ = _W48[k] @ subspace_Q
            tr = np.trace(subspace_Q.conj().T @ WQ)
            s += CHAR[n][_CLS[k]] * tr
        out[n] = s / 48.0 * d
    return {n: {"re": float(np.real(v)), "abs": float(abs(v))}
            for n, v in out.items()}


def projector_rank(name, Q):
    """rank of the isotypic projector P_n restricted to the subspace."""
    P = np.zeros((128, 128), dtype=complex)
    d = CHAR[name]["E"]
    for k in range(48):
        P += CHAR[name][_CLS[k]] * _W48[k]
    P = d * P / 48.0
    Psub = Q.conj().T @ P @ Q
    ev = np.linalg.eigvalsh(0.5 * (Psub + Psub.conj().T))
    return float(np.sum(ev > 0.5)), float(np.abs(P @ P - P).max())


def equivar_defect(p):
    M0 = bloch30((0.0, 0.0, 0.0), p)
    d = 0.0
    for W in _W48:
        d = max(d, float(np.abs(M0 @ W - W @ M0).max()))
    return d


def measure_point(p, m_target=4):
    out = {"gap_gamma": gap_gamma(p)}
    M0 = bloch30((0.0, 0.0, 0.0), p)
    ev, Vv = np.linalg.eig(M0)
    order = np.argsort(np.abs(np.angle(ev)))
    m_exact = int((np.abs(np.angle(ev[order])) < 1e-10).sum())
    sel = order[:m_target]
    Q, _ = np.linalg.qr(Vv[:, sel])
    out["m_exact_lt_1e_10"] = m_exact
    out["top4_gap"] = float(np.abs(np.angle(ev[order[m_target - 1]])))
    out["content"] = content_of(Q)
    ranks = {}
    for n in ("A1g", "T1u"):
        r, idem = projector_rank(n, Q)
        ranks[n] = {"rank": r, "projector_idempotent_defect": idem}
    out["projector_ranks"] = ranks
    # V186 tangent maps on the m_target-gap subspace
    h = 1e-5
    S = []
    for jj in range(3):
        qp, qm = [0.0] * 3, [0.0] * 3
        qp[jj], qm[jj] = h, -h
        Dd = (bloch30(tuple(qp), p) - bloch30(tuple(qm), p)) / (2 * h)
        S.append(1j * (Q.conj().T @ Dd @ Q))
    nrm = [float(np.linalg.norm(s)) for s in S]
    mean2 = sum(n * n for n in nrm) / 3.0
    comm = max(float(np.linalg.norm(S[a] @ S[b] - S[b] @ S[a]))
               for a, b in ((0, 1), (1, 2), (2, 0)))
    out["rho_comm"] = comm / mean2 if mean2 > 0 else 0.0
    idx = np.arange(179) + 0.5
    phi = np.arccos(1 - 2 * idx / 179)
    th = np.pi * (1 + 5 ** 0.5) * idx
    fdirs = np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th),
                      np.cos(phi)], axis=1)
    spd = []
    for kd in fdirs:
        Ax = kd[0] * S[0] + kd[1] * S[1] + kd[2] * S[2]
        spd.extend(sorted(abs(x) for x in np.linalg.eigvalsh(
            0.5 * (Ax + Ax.conj().T))))
    spd = np.array(spd)
    nz = spd[spd > 1e-6]
    if nz.size:
        out["v_min"] = float(nz.min())
        out["v_max"] = float(nz.max())
        out["A_top"] = float((nz.max() - nz.min()) / (nz.max() + nz.min()))
    # multiplet degree (top of the m_target gap levels vs radius)
    ws = []
    for r in V277.FINE_RAD:
        w = np.inf
        for d in V277.CUBIC_DIRS:
            ph = np.sort(np.abs(np.angle(np.linalg.eigvals(
                bloch30(tuple(r * x for x in d), p)))))
            w = min(w, float(ph[m_target - 1]))
        ws.append(w)
    out["degree_fine"] = V277.deg_fit(V277.FINE_RAD, ws)
    out["gaps_fine"] = [float(w) for w in ws]
    # unitarity at random q + equivariance defect
    rng = np.random.default_rng(278)
    ud = 0.0
    for _ in range(3):
        qq = tuple(float(x) for x in rng.uniform(-math.pi, math.pi, 3))
        Mm = bloch30(qq, p)
        ud = max(ud, float(np.abs(Mm @ Mm.conj().T
                                  - np.eye(128, dtype=complex)).max()))
    out["unitarity_defect"] = ud
    out["equivar_defect_gamma"] = equivar_defect(p)
    return out


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    t0 = time.time()
    checks, data = [], {}

    # ---- T0: anchors -----------------------------------------------------
    deg_art = art("V277_bond_decoration_cone_scan",
                  "data.machinery_128.bare_top64_degree")
    uni_art = art("V277_bond_decoration_cone_scan",
                  "data.machinery_128.uniform_triplet_start_degree")
    def bare_step(q):
        M = np.eye(128, dtype=complex)
        for t, (even, j) in enumerate(SC.STAGES):
            S = SC.stage_matrix(even, j, q[j], GCLOSED, GCLOSED)
            M = np.kron(S, V277.U16) @ M
        return M
    ws = [V277.top_m(bare_step, 64, r, V277.CUBIC_DIRS)
          for r in V277.FINE_RAD]
    bare_deg = V277.deg_fit(V277.FINE_RAD, ws)
    wsu = [V277.top_m(lambda q: V277.bloch_dec(q, [(G10, G10, G10)] * 6),
                      64, r, V277.CUBIC_DIRS) for r in V277.FINE_RAD]
    uni_deg = V277.deg_fit(V277.FINE_RAD, wsu)

    v192_art = float(art("V192_engine_two_angle_cone",
                         "data.equivariance.Gamma"))
    equi_pub = equivar_defect(np.zeros(30))
    t0_ok = bool(abs(bare_deg - deg_art) < 1e-6
                 and abs(uni_deg - uni_art) < 1e-6
                 and equi_pub < max(1e-9, 10 * v192_art))
    data["t0"] = {"bare_degree": bare_deg, "bare_degree_v277_art": deg_art,
                  "uniform_degree": uni_deg,
                  "uniform_degree_v277_art": uni_art,
                  "published_equivar_defect": equi_pub,
                  "v192_artifact_gamma_defect": v192_art}
    checks.append({"name": "T0 anchors: bare closed-clock top-64 degree "
                          "{:.6f} vs V277 artefact {:.6f}; uniform triplet "
                          "{:.6f} vs {:.6f}; published GAMMA equivariance "
                          "defect {:.2e} (V192 artefact {:.2e})"
                   .format(bare_deg, deg_art, uni_deg, uni_art, equi_pub,
                           v192_art),
                   "pass": t0_ok})
    if not t0_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V278 MACHINERY-DRIFT: anchor disagreement; nothing "
              "downstream interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- search -----------------------------------------------------------
    rows = parpool.pmap(_start_worker, list(range(NSTARTS)), chunk=4)
    rows.sort(key=lambda r: r["gap"])
    best = rows[0]
    data["search"] = {"starts": NSTARTS, "maxfev": MAXFEV, "seed": SEED,
                      "box": [0.0, BOX], "best": best, "top5": rows[:5],
                      "all_gaps": [r["gap"] for r in rows]}
    checks.append({"name": "T1 search: {} starts (seed {}), best GAMMA "
                          "gap {:.2e} at x = {}"
                   .format(NSTARTS, SEED, best["gap"],
                           [round(v, 4) for v in best["x"]]),
                   "pass": True})

    meas = measure_point(np.array(best["x"]))
    data["best_measure"] = meas
    legs = {
        "gap_lt_1e_10": meas["gap_gamma"] < 1e-10,
        "m4_parity_mixed": (meas["m_exact_lt_1e_10"] == 4
                            and meas["projector_ranks"]["A1g"]["rank"] >= 1
                            and meas["projector_ranks"]["T1u"]["rank"] >= 3),
        "A_lt_0_1": meas.get("A_top", 9.9) < 0.1,
        "rho_gt_0_5": meas["rho_comm"] > 0.5,
        "vmax_le_1": meas.get("v_max", 9.9) <= 1.0,
    }
    found = all(legs.values())
    binding = [k for k, v in legs.items() if not v]
    branch = "CONE-RANGE2-FOUND" if found else "CONE-RANGE2-ABSENT"
    data["status"] = branch
    data["acceptance_legs"] = legs
    data["binding_legs"] = binding
    print("V278 {}: gap {:.2e}, m_exact {}, A_top {}, rho_comm {:.4f}, "
          "v_max {}, content A1g {} T1u {}, legs failed {}"
          .format(branch, meas["gap_gamma"], meas["m_exact_lt_1e_10"],
                  meas.get("A_top"), meas["rho_comm"], meas.get("v_max"),
                  meas["projector_ranks"]["A1g"]["rank"],
                  meas["projector_ranks"]["T1u"]["rank"],
                  binding or ["none"]), flush=True)
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    print("[V278 done] {}  ({:.0f} s)".format(branch, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V278_range2_parity_cone", "tier": 2,
           "property": "does the 30-parameter range-2 (axial + diagonal "
                       "swap-rotation) decorated circuit class reach the "
                       "parity-mixed A1g(+)T1u conical node at GAMMA "
                       "(tier 2, O10)",
           "residual": True,
           "programme": "V278",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v278_range2_parity_cone.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V278 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
