#!/usr/bin/env python3
r"""V334 - the clock-closure rigidity theorem, machine-certified (tier 2, hard, O4/O8, theory-hunt).

The analytic proof that V333 booked as a conjecture, made real. The theorem
(pre-registration V334, booked before any run):

  (Clock-closure rigidity.) Let M(g) be any finite product of the V333
  stage factors at q = 0 (Gamma), any stage count, any ordering, any
  INTEGER clock multipliers. Each bond stage at Gamma is cos(v) I +
  sin(v) K with K entries in Z[i]; the internal factor U16C is frozen at
  the root clock pi/4 with algebraic entries; hence every entry of M(g)
  and M(g)^k is a Laurent polynomial in z = e^{ig} with algebraic
  coefficients. If M(g0)^k = +/- I at a NONZERO ALGEBRAIC clock g0 for
  some k >= 1, then M(g)^k = +/- I IDENTICALLY in g: every closure at a
  nonzero algebraic clock is a clock-identity closure, and no member's
  closure can be SOLVED by an algebraic clock. Proof: Hermite-Lindemann
  makes z0 = e^{i g0} transcendental over Qbar; the Laurent polynomials
  over Qbar vanish at z0, so they vanish identically. QED.

This experiment certifies the two machine ingredients:
A  entry-form gate: stages are exactly cos(v) I + sin(v) K (K v-
   independent, Gaussian-integer entries); U16C entries algebraic
   (2^-a (sqrt2)^e times Gaussian integers) - the Laurent/Qbar premise;
B  identity-closure census: the 4096 booked ray members sampled at the
   25 equispaced angles (Laurent degree <= 12 in z, DFT alias-free):
   identity steps (all non-constant Fourier coefficients <= 1e-11) -
   expected exactly the 64 mirror members (c_{j+3} = -c_j), closure
   order <= 16 of their constant value;
C  census consistency: the closure set at the rational clocks {1/5,
   1/3, 2/7} equals the identity set (theorem's forced inclusion,
   machine-confirmed) and matches V333's rational_ledger leaves; at the
   12 transcendental clocks {m LOOP/2, m odd 1..23} the closure set must
   contain the identity set (no upper bound is imposed there - any
   accidental closure is booked as a ledger).

Backend: numpy direct path (the V333 machinery by module import); no GPU
queue. Branches: CLOSURE-RIGID-CENSUS-MATCHED / CLOSURE-RIGID-CENSUS-
MISMATCH / MACHINERY-DRIFT (the experiment PASSes with either census
branch; only T0/form failures drop it).
"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
from cl10 import config as _config  # noqa: E402
from cl10 import supercell as SC  # noqa: E402

CFG = _config.CFG

_spec = importlib.util.spec_from_file_location(
    "v333mach",
    os.path.join(ROOT, "experiments", "theory_hunt",
                 "V333_clock_transcendence", "V333_clock_transcendence.py"))
V333 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V333)

build = V333.build
closure_min = V333.closure_min
NMEM = V333.NMEM
MULTS = V333.MULTS
CLOSE_EPS = V333.CLOSE_EPS
K_MAX = V333.K_MAX

NANG = 25                        # alias-free for Laurent degree <= 12
CONST_EPS = 1e-11                # identity-step Fourier gate
RATS = ("1/5", "1/3", "2/7")     # the booked rational consistency clocks
LOOP = float(CFG.physics.loop)   # 1/(16 pi), transcendental
TRANS = [("LOOP*%d/2" % m, m * LOOP / 2.0) for m in range(1, 24, 2)]


def art(name, key):
    """A sibling verdict JSON, read by experiment name (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def mirror_member(c):
    """The booked mirror class c_{j+3} = -c_j."""
    return all(c[3 + j] == -c[j] for j in range(3))


def gaussian_dev(x):
    """Distance of a complex number to the nearest Gaussian integer."""
    return math.hypot(x.real - round(x.real), x.imag - round(x.imag))


def algebraic_dev(x):
    """min over a <= 6, e <= 2 of the Gaussian-integer distance of
    2^-a (sqrt2)^e x: a membership test for Z[1/2, sqrt2, i]."""
    best = math.inf
    for a in range(7):
        for e in range(3):
            s = (2.0 ** (-a)) * ((math.sqrt(2.0)) ** e)
            best = min(best, gaussian_dev(x * s))
    return best


def main():
    t0wall = time.time()
    data = {"class": "the V333 booked ray class: 4096 members, per-stage "
                     "multipliers c_j in {-2,-1,1,2}, the six STAGES at "
                     "Gamma, internal factor U16C frozen at pi/4",
            "theorem": "clock-closure rigidity (Hermite-Lindemann): "
                       "closure at a nonzero algebraic clock forces a "
                       "clock-identity closure in this class",
            "backend": "numpy direct (V333 machinery by module import)"}
    checks = []

    # ---- gate A: the entry-form premise, exactly ---------------------
    worst_stage_dev = 0.0
    worst_kconst = 0.0
    for even, jj in SC.STAGES:
        Ks = []
        for v in (0.7, 1.3, 1.9):
            S = SC.stage_matrix(even, jj, 0.0, v, v)
            Kv = (S - math.cos(v) * np.eye(SC.NSUB, dtype=complex)) \
                / math.sin(v)
            Ks.append(Kv)
            worst_stage_dev = max(worst_stage_dev,
                                  max(gaussian_dev(x) for x in
                                      Kv.reshape(-1)))
        worst_kconst = max(worst_kconst,
                           float(np.abs(Ks[0] - Ks[1]).max()),
                           float(np.abs(Ks[0] - Ks[2]).max()))
    A1 = (worst_stage_dev <= 1e-13) and (worst_kconst <= 1e-13)
    worst_u16 = 0.0
    for x in V333.U16C.reshape(-1):
        worst_u16 = max(worst_u16, algebraic_dev(complex(x)))
    A2 = worst_u16 <= 1e-13
    data["entry_form"] = {"worst_gaussian_dev": worst_stage_dev,
                          "worst_stage_v_independence": worst_kconst,
                          "worst_U16C_algebraic_dev": worst_u16}
    checks.append({"name": "A stage entry form: K Gaussian (<= 1e-13) and "
                           "v-independent across v = 0.7/1.3/1.9",
                   "pass": bool(A1)})
    checks.append({"name": "A U16C entries in Z[1/2, sqrt2, i] (<= 1e-13)",
                   "pass": bool(A2)})
    if not (A1 and A2):
        data["status"] = "MACHINERY-DRIFT"
        print("V334 MACHINERY-DRIFT: entry-form premise fails", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0a: the V315 identity-closure leaves re-derived ------------
    booked = art("V315_closure_quantizer_census", "data.t0a_identity_closure")
    ident = (1, 1, 1, 1, 1, 1)
    cache = {}
    t0a = {}
    t0a_ok = True
    for label, g in (("pi*1/2", np.pi / 2), ("pi*1/4", np.pi / 4),
                     ("pi*3/4", 3 * np.pi / 4)):
        M = build(ident, g, cache)
        closed, k, dev, _, _ = closure_min(M)
        bk = booked[label]
        ddev = abs(dev - bk["rederived_dev"]) if dev is not None else 9e9
        ok = closed and k == bk["rederived_k"] and ddev <= 1e-13
        t0a_ok = t0a_ok and ok
        t0a[label] = {"k": k, "dev": dev, "v315_k": bk["rederived_k"],
                      "v315_dev": bk["rederived_dev"], "abs_diff_dev": ddev}
    checks.append({"name": "T0a V315 identity closures re-derived (k match, "
                           "dev |d| <= 1e-13)", "pass": bool(t0a_ok)})
    data["t0"] = {"t0a_identity_closures": t0a}
    if not t0a_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V334 MACHINERY-DRIFT: identity closures do not re-derive",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0c: unitarity of 16 seeded members at 3 clocks -------------
    unit = 0.0
    I128 = np.eye(128, dtype=complex)
    for i in range(16):
        c = V333.member_of((i * 7919 + 31) % NMEM)
        for g in (0.2, TRANS[3][1], math.pi / 4):
            M = build(c, g, {})
            unit = max(unit, float(np.abs(M.conj().T @ M - I128).max()))
    checks.append({"name": "T0c unitarity of 16 seeded members x 3 clocks "
                           "(< 1e-10)", "pass": unit < 1e-10})
    data["t0"]["unitarity_worst_16x3"] = unit
    if unit >= 1e-10:
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- B: the identity-closure census (25-angle alias-free DFT) ----
    fs = np.arange(NANG) - (NANG // 2)
    ms = np.arange(NANG)
    E = np.exp(-2j * np.pi * np.outer(fs, ms) / NANG) / NANG
    identity_idx = []
    identity_order = {}
    floor_nc = 0.0
    order_hist = {}
    for i in range(NMEM):
        c = V333.member_of(i)
        cache_i = {}
        Mm = np.empty((NANG, 128 * 128), dtype=complex)
        for m in range(NANG):
            Mm[m] = build(c, 2.0 * np.pi * m / NANG, cache_i).reshape(-1)
        F = E @ Mm
        const = F[NANG // 2].reshape(128, 128)
        nc = float(np.abs(np.concatenate([F[:NANG // 2],
                                          F[NANG // 2 + 1:]])).max())
        if nc > floor_nc and not mirror_member(c):
            floor_nc = nc
        if nc <= CONST_EPS:
            identity_idx.append(i)
            Uk = np.eye(128, dtype=complex)
            kk = 0
            for k in range(1, K_MAX + 1):
                Uk = Uk @ const
                dev = min(float(np.abs(Uk - I128).max()),
                          float(np.abs(Uk + I128).max()))
                if dev < CONST_EPS:
                    kk = k
                    break
            identity_order[str(i)] = kk
            order_hist[str(kk)] = order_hist.get(str(kk), 0) + 1
    ident_set = set(identity_idx)
    mirror_set = {i for i in range(NMEM)
                  if mirror_member(V333.member_of(i))}
    B_ok = ident_set == mirror_set and all(k >= 1 for k in
                                           identity_order.values())
    # order of the frozen internal factor A = U16C^6 (U16C is the 16-level
    # internal step; the full cell step is kron(8x8 bond, 16x16 internal))
    A6 = np.linalg.matrix_power(V333.U16C, 6)
    A6_order = 0
    Uk = np.eye(16, dtype=complex)
    for k in range(1, 33):
        Uk = Uk @ A6
        dev = min(float(np.abs(Uk - np.eye(16)).max()),
                  float(np.abs(Uk + np.eye(16)).max()))
        if dev < CONST_EPS:
            A6_order = k
            break
    data["identity_census"] = {
        "n_identity_steps": len(ident_set),
        "identity_set_equals_mirrors": bool(ident_set == mirror_set),
        "n_mirrors": len(mirror_set),
        "closure_order_histogram": order_hist,
        "A6_U16C_pow6_order": A6_order,
        "non_identity_fourier_floor": floor_nc,
        "constant_gate": CONST_EPS}
    checks.append({"name": "B identity steps == the 64 mirror members, all "
                           "with closure order <= 16", "pass": bool(B_ok)})

    # ---- C: closure sets at rational and transcendental clocks -------
    def closure_set_at(g):
        cache_c = {}
        closed = []
        for i in range(NMEM):
            c = V333.member_of(i)
            M = build(c, g, cache_c)
            okk, k, dev, _, _ = closure_min(M)
            if okk:
                closed.append((i, k, dev))
        return closed

    rat_rows = {}
    rat_ok = True
    for key in RATS:
        p, q = key.split("/")
        g = int(p) / int(q)
        closed = closure_set_at(g)
        cs = {i for i, _, _ in closed}
        leaf = art("V333_clock_transcendence",
                   "data.rational_ledger.%s" % key)
        same = cs == ident_set
        leaf_ok = (leaf["closed"] == len(closed)
                   and leaf["scalar_degenerate"] == len(closed)
                   and leaf["found"] == 0)
        rat_ok = rat_ok and same and leaf_ok
        rat_rows[key] = {"closed": len(closed),
                         "equals_identity_set": bool(same),
                         "ledger_leaf_closed": leaf["closed"],
                         "ledger_leaf_scalar": leaf["scalar_degenerate"],
                         "ledger_leaf_found": leaf["found"],
                         "extra": sorted(cs - ident_set)[:8],
                         "missing": sorted(ident_set - cs)[:8]}
    data["rational_consistency"] = rat_rows
    checks.append({"name": "C closure sets at 1/5, 1/3, 2/7 == the identity "
                           "set and match V333 ledger leaves (64/64/0)",
                   "pass": bool(rat_ok)})

    trans_rows = {}
    trans_extra_total = 0
    trans_contains_ok = True
    for name, g in TRANS:
        closed = closure_set_at(g)
        cs = {i for i, _, _ in closed}
        extra = sorted(cs - ident_set)
        missing = sorted(ident_set - cs)
        trans_extra_total += len(extra)
        trans_contains_ok = trans_contains_ok and (not missing)
        trans_rows[name] = {"g": g, "closed": len(cs),
                            "extra_n": len(extra), "extra": extra[:8],
                            "missing_n": len(missing)}
    data["transcendental_census"] = trans_rows
    checks.append({"name": "C transcendental clocks: identity set closes at "
                           "all 12 (containment; extras allowed and booked)",
                   "pass": bool(trans_contains_ok)})

    # ---- branch ------------------------------------------------------
    if ident_set == mirror_set and rat_ok and trans_contains_ok \
            and trans_extra_total == 0 and A1 and A2:
        status = "CLOSURE-RIGID-CENSUS-MATCHED"
    else:
        status = "CLOSURE-RIGID-CENSUS-MISMATCH"
    data["status"] = status
    data["transcendental_extras_total"] = trans_extra_total
    checks.append({"name": "census complete (4096 members: 25-angle DFT "
                           "+ 3 rational + 12 transcendental clocks)",
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0wall
    print("V334 %s: identity set %d (mirrors %d), closure orders %s, "
          "extras at transcendental clocks %d, floor %.2e, A6 order %d"
          % (status, len(ident_set), len(mirror_set), order_hist,
             trans_extra_total, floor_nc, A6_order), flush=True)
    _emit(data, checks, ok=ok)
    print("[V334 done] %s" % status, flush=True)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V334_clock_closure_rigidity", "tier": 2,
           "property": "the clock-closure rigidity theorem certified "
                       "on the class: entry form, identity census, "
                       "rational/transcendental closure sets "
                       "(tier 2, O4/O8)",
           "residual": True,
           "programme": "theory-hunt V334",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v334_clock_closure_rigidity.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V334] pass=%s branch=%s" % (out["pass"], out["branch"]),
          flush=True)


if __name__ == "__main__":
    sys.exit(main())
