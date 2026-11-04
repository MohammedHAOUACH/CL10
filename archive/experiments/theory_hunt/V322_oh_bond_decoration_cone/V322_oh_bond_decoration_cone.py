#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V322 - the O_h bond-decoration cone census (tier 2, hard, O10/O4,
theory-hunt; closes V191's stated open piece).

**The question (pre-registered).** Does the O_h-equivariant SITE-side
decoration class - exact orbit sums of the V191 range-1 hopping-orbit
census on the 2x2x2 cell module - applied to the published six-stage
step carry a conical GAMMA node (gap < 1e-10, fine degree >= 0.99,
v_max <= 1, rho_comm > 0.5; V277's acceptance convention), and does
any member reduce the step's O_h defect toward the stage-sum value
(V191: "the remaining work is a stage factorisation of that bond
layer")?

**Class (closed; census corrected BEFORE the first real run, booked
in the pre-registration).** Range-1 terms (i, j in 8 sites,
v in {-1,0,1}^3; 1728 terms) quotiented by O_h (48 site permutations
in V191's +-1 coordinate convention): 56 orbits, sizes 8 x 6,
24 x 30, 48 x 20. Members: D(q) = exp(i THETA sum_o w_o A_o(q)),
A_o = orbit Bloch sum (M + M^dag, V191 convention - orbits are closed
under dagger since -I is in O_h), w_o in {0,+-1}: identity + all
single-orbit members (112) + double-orbit members restricted to the
six size-8 orbits (60) = 173 members. THETA = config angle.alpha_10.
Decoration multiplies each stage's site block: M -> kron(D(q) @
S_t(q), U16) @ M; internal U16 and the clocks g = beta =
angle.branches.closed untouched. Provenance combinatorial +
config-declared: zero-fitted.

**Measurement amendment booked BEFORE the run.** V277 measured the
degree at m0 = 64 (the pinned 64+64 closed-clock spectrum, V192/V132)
because its W(0) = I inherited it; here D(GAMMA) != I in general, so
the gap multiplicity is measured per member: m* = #{|phase| < 1e-8}
at GAMMA. Both legs run and book: the pinned m0 = 64 convention and
the adaptive m*; the DECISION uses m* (honest multiplicity), the
pinned leg is booked for comparability.

**T0 (machinery identity).** SC.bond_product defect over V191's exact
MOM momenta reproduces the V191 ARTIFACT leaves product_o_h_defect
{M: 1.0, R: 0.7071067811865476, generic: 0.7044629490796265} to rel
1e-9 (the docstring's 1.07e-1 is stale narrative, the artifact leaf is
authority - booked); and the identity member (all w = 0) reproduces
V277's leaf bare_top64_degree 1.9985044777023986 to < 1e-4 and the
closed-clock GAMMA gap < 1e-10.

**Branch (pre-registered words).** CONE-OH-FOUND iff some member
satisfies ALL four legs (gap < 1e-10 AND degree_fine(m*) >= 0.99 AND
v_max <= 1 AND rho_comm > 0.5); CONE-OH-NEGATIVE otherwise. Both
book: top-10 by degree_fine(m*), the minimum-defect member (the
factorisation leg, booked not gated), and the selectivity table
(count members keeping gap < 1e-10 - which orbit sums gap the node).
"""
import itertools
import json
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402

from cl10 import config as _config                        # noqa: E402
from cl10 import supercell as SC                          # noqa: E402

CFG = _config.CFG
GCLOSED = float(CFG.angle.branches.closed)
THETA = float(CFG.angle.alpha_10)

_B16, _U32, U16C = SC.sectors(GCLOSED)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# ------------------------------------------------------------- the group --
G48 = []
for _perm in itertools.permutations(range(3)):
    for _sg in itertools.product((1, -1), repeat=3):
        _R = np.zeros((3, 3))
        for _i, (_p, _s) in enumerate(zip(_perm, _sg)):
            _R[_i, _p] = _s
        G48.append(_R)
NS = SC.bond_product((0.0, 0.0, 0.0)).shape[0]


def perm_of(R3):
    out = []
    for e in range(NS):
        x = np.array([2.0 * ((e >> k) & 1) - 1.0 for k in range(3)])
        y = R3 @ x
        out.append(int(round((y[0] + 1) / 2))
                   | (int(round((y[1] + 1) / 2)) << 1)
                   | (int(round((y[2] + 1) / 2)) << 2))
    return np.array(out)


PERM = [perm_of(R) for R in G48]
PW8 = []
for p in PERM:
    Wm = np.zeros((NS, NS))
    for i, j in enumerate(p):
        Wm[j, i] = 1.0
    PW8.append(Wm)


# --------------------------------------------------- the orbit census ----
def term_key(t):
    i, j, v = t
    best = None
    for R3, P in zip(G48, PERM):
        vv = tuple(int(round(x)) for x in R3 @ np.array(v))
        k = (int(P[i]), int(P[j]), vv)
        if best is None or k < best:
            best = k
    return best


def build_orbits():
    terms = [(i, j, v) for i in range(NS) for j in range(NS)
             for v in itertools.product((-1, 0, 1), repeat=3)]
    orbs = {}
    for t in terms:
        orbs.setdefault(term_key(t), []).append(t)
    return sorted(orbs.values(), key=lambda m: term_key(m[0]))


def bloch_orbit(mem, q):
    A = np.zeros((NS, NS), dtype=complex)
    for (i, j, v) in mem:
        A[i, j] += np.exp(1j * (q[0] * v[0] + q[1] * v[1] + q[2] * v[2]))
    return A + A.conj().T


def dec_site(q, w):
    """D(q) = exp(i THETA sum w_o A_o(q)), exact via the 8x8 eigh."""
    S8 = np.zeros((NS, NS), dtype=complex)
    for wo, orb in zip(w, ORBITS):
        if wo:
            S8 += wo * bloch_orbit(orb, q)
    ev, Vm = np.linalg.eigh(THETA * S8)
    return (Vm * np.exp(1j * ev)) @ Vm.conj().T


def step_dec(q, w):
    Dq = dec_site(q, w)
    M = np.eye(128, dtype=complex)
    for (even, j) in SC.STAGES:
        S = SC.stage_matrix(even, j, q[j], GCLOSED, GCLOSED)
        M = np.kron(Dq @ S, U16C) @ M
    return M


CUBIC_DIRS = []
for _d in ((1, 0, 0), (0, 1, 0), (0, 0, 1),
           (1, 1, 0), (1, 0, 1), (0, 1, 1),
           (1, -1, 0), (1, 0, -1), (0, 1, -1),
           (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)):
    _v = np.array(_d, float)
    CUBIC_DIRS.append(_v / np.linalg.norm(_v))
FINE_RAD = (0.3, 0.1, 0.03, 0.01, 0.003)


def phases(M):
    return np.sort(np.abs(np.arccos(np.clip(
        np.linalg.eigvalsh(0.5 * (M + M.conj().T)), -1.0, 1.0))))


def top_m(bf, m, r, dirs):
    w = np.inf
    for d in dirs:
        p = phases(bf(tuple(r * x for x in d)))
        w = min(w, float(p[m - 1]))
    return w


def deg_fit(radii, ws):
    with np.errstate(divide="ignore"):
        return float(np.polyfit(np.log(list(radii)), np.log(
            [float(w) for w in ws]), 1)[0])


def measure(w):
    bf = lambda q: step_dec(q, w)  # noqa: E731
    M0 = bf((0.0, 0.0, 0.0))
    ph0 = np.sort(np.abs(np.angle(np.linalg.eigvals(M0))))
    gap = float(ph0[0])
    mstar = int((ph0 < 1e-8).sum())
    out = {"gap_gamma": gap, "m_star": mstar}
    for tag, m0 in (("m64", 64), ("ms", max(mstar, 1))):
        ws = [top_m(bf, m0, r, CUBIC_DIRS) for r in FINE_RAD]
        out["degree_" + tag] = deg_fit(FINE_RAD, ws)
    # tangent maps on the m*-gap subspace at GAMMA (V186 route, V277 copy)
    ev, Vv = np.linalg.eig(M0)
    sel = np.argsort(np.abs(np.angle(ev)))[:max(mstar, 1)]
    Q, _ = np.linalg.qr(Vv[:, sel])
    h = 1e-5
    S = []
    for jj in range(3):
        qp, qm = [0.0] * 3, [0.0] * 3
        qp[jj], qm[jj] = h, -h
        Dd = (bf(tuple(qp)) - bf(tuple(qm))) / (2 * h)
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
        out["n_flat"] = int((spd <= 1e-6).sum())
    # O_h defect of the full 128-dim step at GAMMA + 8 fib |dq| = 0.1
    d8 = np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th),
                   np.cos(phi)], axis=1)[:8] * 0.1
    dmax = 0.0
    for R3, W8 in zip(G48, PW8):
        WP = np.kron(W8, np.eye(16))
        for dvec in d8:
            qq = tuple(float(x) for x in dvec)
            lhs = WP @ step_dec(qq, w) @ WP.T
            rhs = step_dec(tuple(R3 @ np.array(qq)), w)
            dmax = max(dmax, float(np.abs(lhs - rhs).max()))
    out["oh_defect"] = dmax
    return out


def main():
    t0 = time.time()
    checks, data = [], {}
    global ORBITS
    ORBITS = build_orbits()
    sizes = [len(o) for o in ORBITS]
    data["orbit_count"] = len(ORBITS)
    data["orbit_sizes_sorted"] = sorted(sizes)
    if len(ORBITS) != 56 or sorted(sizes) != [8] * 6 + [24] * 30 + [48] * 20:
        print("V322 MACHINERY-DRIFT: orbit census deviates from the "
              "pre-registration (sizes %s)" % sizes, flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0a: V191 engine defect on V191's own machinery and momenta ----
    MOM = {"Gamma": [(0.0, 0.0, 0.0)],
           "M": [(np.pi / 2, np.pi / 2, 0.0)],
           "R": [(np.pi / 2, np.pi / 2, np.pi / 2)],
           "generic": [(0.3, -0.2, 0.7), (0.9, 0.9, 0.1), (0.1, 0.4, -0.6)]}

    def defect(fn, qs):
        d = 0.0
        for R3, W in zip(G48, PW8):
            for q in qs:
                qq = np.array(q)
                d = max(d, float(np.abs(W @ fn(qq) @ W.T - fn(R3 @ qq)).max()))
        return d

    ed = {k: defect(SC.bond_product, v) for k, v in MOM.items()}
    leaf = art("V191_parity_dirac_cone", "data.engine.product_o_h_defect")
    t0a = {k: abs(ed[k] - float(leaf[k])) / max(float(leaf[k]), 1e-12)
           for k in ("M", "R", "generic")}
    t0a_ok = all(v < 1e-9 for v in t0a.values())
    data["t0_defect_rel"] = t0a
    data["t0_defect_mine"] = ed
    checks.append({"name": "T0a V191 engine-defect identity: SC.bond_product"
                   " defect over V191's MOM momenta vs the artifact leaves "
                   "{{M 1.0, R 0.7071, generic 0.70446}} rel {} (docstring "
                   "1.07e-1 is stale narrative; leaf authority booked)"
                   .format({k: "{:.1e}".format(v) for k, v in t0a.items()}),
                   "pass": bool(t0a_ok)})

    # ---- T0b: identity member reproduces V277 bare degree and the gap ---
    idrow = measure(tuple([0] * 10))
    deg_leaf = art("V277_bond_decoration_cone_scan",
                   "data.machinery_128.bare_top64_degree")
    t0b_ok = (abs(idrow["degree_m64"] - deg_leaf) < 1e-4
              and idrow["gap_gamma"] < 1e-12)
    data["t0_identity"] = {"degree_m64": idrow["degree_m64"],
                           "degree_leaf": deg_leaf,
                           "gap": idrow["gap_gamma"],
                           "oh_defect": idrow["oh_defect"]}
    checks.append({"name": "T0b identity member: degree_m64 {:.7f} vs V277 "
                   "leaf {:.7f} (diff {:.1e} < 1e-4), GAMMA gap {:.1e} < "
                   "1e-12, undecorated step defect {:.3e} (booked)".format(
                       idrow["degree_m64"], deg_leaf,
                       abs(idrow["degree_m64"] - deg_leaf),
                       idrow["gap_gamma"], idrow["oh_defect"]),
                   "pass": bool(t0b_ok)})
    if not (t0a_ok and t0b_ok):
        data["status"] = "MACHINERY-DRIFT"
        print("V322 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- the census: 173 members (amended closed class, booked) ---------
    members = []
    for o1 in range(len(ORBITS)):
        for s1 in (1, -1):
            w = [0] * len(ORBITS)
            w[o1] = s1
            members.append(tuple(w))
    idx8 = [i for i, o in enumerate(ORBITS) if len(o) == 8]
    if len(idx8) != 6:
        print("V322 BUG: size-8 orbit count %d != 6" % len(idx8), flush=True)
        _emit(data, checks, ok=False)
        return 1
    for a in range(6):
        for b in range(a + 1, 6):
            for s1 in (1, -1):
                for s2 in (1, -1):
                    w = [0] * len(ORBITS)
                    w[idx8[a]], w[idx8[b]] = s1, s2
                    members.append(tuple(w))
    members = [tuple([0] * len(ORBITS))] + members
    if len(members) != 173:
        print("V322 BUG: member count %d != 173" % len(members), flush=True)
        _emit(data, checks, ok=False)
        return 1
    rows = [idrow]
    for k, w in enumerate(members[1:], 1):
        r = measure(w)
        r["w"] = list(w)
        rows.append(r)
        if k % 25 == 0:
            print("  member %d/172 done (%.0f s)" % (k, time.time() - t0),
                  flush=True)
    data["n_members"] = len(rows)
    idrow["w"] = [0] * len(ORBITS)
    data["identity_row"] = idrow

    def oklegs(r):
        return (r["gap_gamma"] < 1e-10 and r["degree_ms"] == r["degree_ms"]
                and r["degree_ms"] >= 0.99 and r.get("v_max", 1e9) <= 1.0
                and r["rho_comm"] > 0.5)

    winners = [r for r in rows if oklegs(r)]
    gapless = [r for r in rows if r["gap_gamma"] < 1e-10]
    bestdeg = sorted(rows, key=lambda r: -(r["degree_ms"]
                                           if r["degree_ms"] == r["degree_ms"]
                                           else -1e9))[:10]
    mindef = min(rows, key=lambda r: r["oh_defect"])
    branch = "CONE-OH-FOUND" if winners else "CONE-OH-NEGATIVE"
    data["status"] = branch
    data["winners"] = [{k: v for k, v in r.items() if k != "w"}
                       | {"w": r["w"]} for r in winners]
    data["gapless_count"] = len(gapless)
    data["top10_degree"] = [{k: v for k, v in r.items() if k != "w"}
                            | {"w": r["w"]} for r in bestdeg]
    data["min_defect_member"] = {k: v for k, v in mindef.items()
                                 if k != "w"} | {"w": mindef["w"]}
    data["rows"] = [{k: (v if not isinstance(v, float) or v == v
                         else None) for k, v in r.items()}
                    for r in rows]
    checks.append({"name": "branch {}: {} members gapless (<1e-10), best "
                   "degree_ms {:.4f} (gap {:.1e}, v_max {}, rho_comm {}), "
                   "min O_h defect {:.3e} (V191 undecorated 128-step leg "
                   "booked {:.3e})".format(
                       branch, len(gapless),
                       bestdeg[0]["degree_ms"] if bestdeg[0]["degree_ms"]
                       == bestdeg[0]["degree_ms"] else float("nan"),
                       bestdeg[0]["gap_gamma"],
                       bestdeg[0].get("v_max", "n/a"),
                       bestdeg[0]["rho_comm"], mindef["oh_defect"],
                       idrow["oh_defect"]), "pass": True})
    data["wall_time"] = time.time() - t0
    _emit(data, checks, ok=all(c["pass"] for c in checks))
    print("[V322] {} gapless {}/172 min-def {:.2e} ({} s)".format(
        branch, len(gapless), mindef["oh_defect"], data["wall_time"]),
        flush=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V322_oh_bond_decoration_cone", "tier": 2,
           "property": "O_h-equivariant site-side bond-decoration cone "
                       "census on the six-stage step; cone acceptance legs "
                       "and the stage-factorisation equivariance table "
                       "(tier 2, hard, O10/O4)",
           "residual": True, "programme": "theory-hunt V322",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v322_oh_bond_decoration_cone.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V322 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
