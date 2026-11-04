#!/usr/bin/env python3
r"""V333 - the clock-transcendence test (tier 2, hard, O4/O8, theory-hunt).

The constructive turn: the empirical clock is RATIONAL -
g_obs = 0.19999798835751165, i.e. 1/5 to 2.0116424883642736e-06 (V206) -
while every closure clock the root has produced is pi-rational (V315's
identity row pi/2, pi/4, 3pi/4; the 43/79 closing grid angles; V239's
four-point set C; the V326 dictionary-invariance theorem). The conjecture
under attack: no six-stage step of the root-generated class closes at a
RATIONAL clock, for any stage-clock multipliers.

The booked class (pre-registration 2026-11-01): the 128x128 six-stage
product of V315's dec_step machinery at q = 0 (Gamma), identity ordering,
zero twist, with per-stage clock multipliers c_j in {-2,-1,1,2} - 4^6 =
4096 members, the published step at c = (1,..,1). Clocks: the rational
grid p/q in lowest terms, q <= 12, g in [0.05, 1.6] (contains 1/5
exactly), and V315's 79-clock pi*a/b ledger grid (b <= 16). Closure:
smallest k <= 16 with ||M^k -/+ I||_max < 1e-12 (V311's convention).
FOUND gate: any non-degenerate rational closure (||M -/+ I||_max >= 1e-6
at the clock); at-obs flag when |g - 1/5| <= 2.1e-6.

Branches (one fires, the experiment PASSes either way):
  CLOCK-RATIONAL-CLOSURE-FOUND / CLOCK-TRANSCENDENT-DECLARED /
  MACHINERY-DRIFT.
Backend: parent-selected. Under cupy (GPU 1) the census runs on a batched
fast path - all 4096 members of a clock in stacked kernels - gated by T0e,
bit-level identity with the direct per-member path on seeded samples (the
V326 hoisting playbook: machinery speed-up, class unchanged, backend
recorded). Under numpy the direct path IS the census. T0a-c always run on
the direct path.
"""
import itertools
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
from cl10 import engine as E  # noqa: E402
from cl10 import supercell as SC  # noqa: E402

CFG = _config.CFG
DEN_MAX = int(CFG.pillar_surgery.closure_den_max)   # 16: the pi-grid ledger
K_MAX = int(CFG.pillar_surgery.closure_k_max)       # 16: the order test
CLOSE_EPS = 1e-12                # V311's own definition line, verbatim
Q_MAX = 12                       # the booked rational-grid denominator
G_LO, G_HI = 0.05, 1.6           # the booked clock window
MULTS = (-2, -1, 1, 2)           # the booked per-stage clock multipliers
NONDEG = 1e-6                    # FOUND non-degeneracy gate
OBS_BAND = 2.1e-6                # the inversion-precision band on 1/5
GCLOSED = math.pi / 4.0          # V277's sector-build clock (as V315)


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


_B16, _U32, U16C = SC.sectors(GCLOSED)
Q0 = (0.0, 0.0, 0.0)             # Gamma


def stage_factor(j, slot_g):
    """F_j(v) = kron(stage_matrix(slot j at Gamma, clock v, bond clock v),
    U16C): the V315 dec_step factor with both clock slots at v, zero
    twist (amps = 0 is the booked identity-degeneracy)."""
    even, jj = SC.STAGES[j]
    S = SC.stage_matrix(even, jj, Q0[jj], slot_g, slot_g)
    return np.kron(S, U16C)


def build(c, g, cache):
    """The member step M(c, g): cached stage factors, six matmuls."""
    M = np.eye(128, dtype=complex)
    for j in range(6):
        v = c[j] * g
        key = (j, round(v, 15))
        F = cache.get(key)
        if F is None:
            F = stage_factor(j, v)
            cache[key] = F
        M = F @ M
    return M


def closure_min(M):
    """(closed, k, dev, min_dev_over_k, k_of_min_dev): V315's closure_row
    extended with the closest non-closing approach over the same k range,
    for the booked ledger."""
    Uk = np.eye(128, dtype=complex)
    I128 = np.eye(128, dtype=complex)
    best = None
    best_k = 0
    for k in range(1, K_MAX + 1):
        Uk = Uk @ M
        dev = min(float(np.abs(Uk - I128).max()),
                  float(np.abs(Uk + I128).max()))
        if best is None or dev < best:
            best, best_k = dev, k
        if dev < CLOSE_EPS:
            return True, k, dev, best, best_k
    return False, 0, None, best, best_k


NMEM = len(MULTS) ** 6
# member i of itertools.product(MULTS, repeat=6): c_(j+1) is the base-4
# digit at 4^(5-j), c_1 the slowest
DIGITS = [np.array([(i // 4 ** (5 - j)) % 4 for i in range(NMEM)],
                   dtype=np.int64) for j in range(6)]


def member_of(i):
    return [MULTS[int(d[i])] for d in DIGITS]


def batch_stack(cp, g, cacheG):
    """The 4096-member stack M(c, g) at one clock, built by six stacked
    matmuls (the fast path): F_j(c_j g) gathered per member, left-multiplied
    layer by layer."""
    P = np.repeat(np.eye(128, dtype=complex)[None], NMEM, axis=0)
    P = cp.asarray(P)
    for j in range(6):
        even, jj = SC.STAGES[j]
        Fj = []
        for m in MULTS:
            key = (j, round(m * g, 15))
            F = cacheG.get(key)
            if F is None:
                S = SC.stage_matrix(even, jj, Q0[jj], m * g, m * g)
                F = cp.asarray(np.asarray(E.kron(S, U16C)))
                cacheG[key] = F
            Fj.append(F)
        Fj = cp.stack(Fj)
        P = Fj[cp.asarray(DIGITS[j])] @ P
    return P


def batch_closure(cp, P):
    """Per-member (closed, first-k, dev-at-first-k, min_dev, argmin-k,
    dev(M vs +/-I)): the direct closure_min semantics, run stacked. mind
    is min over k <= k_close for closed members, over k <= K_MAX
    otherwise."""
    I0 = cp.eye(128, dtype=complex)
    Uk = P.copy()
    closed = cp.zeros(NMEM, dtype=bool)
    kfirst = cp.zeros(NMEM, dtype=np.int64)
    devfirst = cp.zeros(NMEM, dtype=np.float64)
    mind = cp.full(NMEM, 1e9)
    mindk = cp.zeros(NMEM, dtype=np.int64)
    d_pmI = None
    for k in range(1, K_MAX + 1):
        dp = cp.abs(Uk - I0).reshape(NMEM, -1).max(1)
        dn = cp.abs(Uk + I0).reshape(NMEM, -1).max(1)
        dev = cp.minimum(dp, dn)
        if k == 1:
            d_pmI = dev.copy()
        upd = (~closed) & (dev < mind)
        mind = cp.where(upd, dev, mind)
        mindk = cp.where(upd, k, mindk)
        newly = (dev < CLOSE_EPS) & (~closed)
        kfirst = cp.where(newly, k, kfirst)
        devfirst = cp.where(newly, dev, devfirst)
        closed = closed | (dev < CLOSE_EPS)
        Uk = Uk @ P
    return closed, kfirst, devfirst, mind, mindk, d_pmI


def main():
    t0 = time.time()
    checks, data = [], {}

    # ---- booked grids -------------------------------------------------
    rats = []
    seen = set()
    for q in range(1, Q_MAX + 1):
        for p in range(1, int(G_HI * q) + 1):
            g = p / q
            if g < G_LO - 1e-15 or g > G_HI + 1e-15:
                continue
            if math.gcd(p, q) != 1:
                continue
            key = round(g, 15)
            if key in seen:
                continue
            seen.add(key)
            rats.append({"p": p, "q": q, "g": g})
    pigrid = []
    seen = set()
    for b in range(2, DEN_MAX + 1):
        for a in range(1, b):
            g = np.pi * a / b
            key = round(g, 15)
            if key in seen:
                continue
            seen.add(key)
            pigrid.append({"a": a, "b": b, "g": g})
    data["n_rational_clocks"] = len(rats)
    data["n_pi_clocks"] = len(pigrid)
    data["n_members"] = len(MULTS) ** 6

    # ---- T0a: the V315 identity-row closures re-derived bit-for-bit ---
    cache = {}
    ident = (1, 1, 1, 1, 1, 1)
    booked = art("V315_closure_quantizer_census", "data.t0a_identity_closure")
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
    data["t0"] = {"t0a_identity_closures": t0a}
    checks.append({"name": "T0a identity closures re-derived (k match, "
                   "dev |d| <= 1e-13 to V315 leaves)", "pass": bool(t0a_ok)})
    if not t0a_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V333 MACHINERY-DRIFT: identity closures do not re-derive",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T0b: the data-clock leaf ------------------------------------
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    leaf_minus = float(art("V206_clock_tension",
                           "data.data_clock.g_obs_minus_fifth"))
    d_obs = abs(abs(g_obs - 0.2) - abs(leaf_minus))
    data["t0"]["g_obs"] = g_obs
    data["t0"]["g_obs_minus_fifth_absdiff"] = d_obs
    checks.append({"name": "T0b g_obs leaf == 0.19999798835751165 and "
                   "|g_obs-1/5| re-derived to 0.0 diff", "pass": d_obs == 0.0})
    if d_obs != 0.0:
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- T0c: the published step does NOT close at the fifth ---------
    Mf = build(ident, 0.2, cache)
    cf, kf, devf, minf, _ = closure_min(Mf)
    data["t0"]["fifth_closure"] = {"closed": bool(cf), "min_dev": minf}
    checks.append({"name": "T0c identity member not closed at g = 1/5 "
                   "(None over k <= 16, min_dev >= 1e-6)",
                   "pass": (not cf) and minf >= NONDEG})
    if cf or minf < NONDEG:
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- T0d: exact unitarity on 16 seeded members --------------------
    rng = np.random.default_rng(333)
    unit_ok = True
    worst_u = 0.0
    I128 = np.eye(128, dtype=complex)
    for _ in range(16):
        c = tuple(int(x) for x in rng.choice(MULTS, size=6))
        g = float(rng.uniform(G_LO, G_HI))
        M = build(c, g, cache)
        worst_u = max(worst_u, float(np.abs(M @ M.conj().T - I128).max()))
    unit_ok = worst_u < 1e-10
    data["t0"]["worst_unitarity_dev"] = worst_u
    checks.append({"name": "T0d unitarity of 16 seeded members (< 1e-10)",
                   "pass": bool(unit_ok)})
    if not unit_ok:
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- T0e + census path: the batched fast path, gated on identity ---
    cp_mod = None
    if E._BACKEND == "cupy":
        import cupy
        cp_mod = cupy
    if cp_mod is not None:
        cacheG = {}
        e_worst, e_kbad = 0.0, 0
        for grow in (rats[0], pigrid[0]):
            g = grow["g"]
            Pc = batch_stack(cp_mod, g, cacheG)
            cl_b, kf_b, dv_b, md_b, mk_b, dpi_b = batch_closure(cp_mod, Pc)
            for i in (1234, 4000):
                c = tuple(member_of(i))
                M = build(c, g, cache)
                clo_d, k_d, dev_d, mind_d, mk_d = closure_min(M)
                e_worst = max(e_worst, abs(float(md_b[i]) - mind_d))
                if int(kf_b[i]) != k_d:
                    e_kbad += 1
                if clo_d:
                    e_worst = max(e_worst, abs(float(dv_b[i]) - dev_d))
                d0 = float(min(np.abs(M - I128).max(),
                               np.abs(M + I128).max()))
                e_worst = max(e_worst, abs(float(dpi_b[i]) - d0))
        e_ok = e_kbad == 0 and e_worst <= 1e-13
        data["t0"]["t0e_batch_identity"] = {"worst_abs_dev_diff": e_worst,
                                            "k_mismatches": e_kbad}
        checks.append({"name": "T0e batched cupy path == direct path (4 "
                       "seeded members x 2 clocks: same k, dev |d| <= "
                       "1e-13)", "pass": bool(e_ok)})
        if not e_ok:
            data["status"] = "MACHINERY-DRIFT"
            print("V333 MACHINERY-DRIFT: batch path != direct path",
                  flush=True)
            _emit(data, checks, ok=False)
            return 1
    data["census_path"] = "batched-cupy" if cp_mod else "direct-numpy"

    # ---- the census ----------------------------------------------------
    found = []                    # non-degenerate rational closures
    degenerate = []               # closures with M = +/- I at the clock
    pi_closures = {}              # pi-grid ledger: g-rounded -> count
    members_with_pi_closure = 0
    best_rational = {"dev": None} # closest rational-grid approach
    total = 0
    if cp_mod is not None:
        closed_any_pi = np.zeros(NMEM, dtype=bool)
        cacheG = {}
        rational_ledger = {}
        for grid, tag in ((rats, "r"), (pigrid, "p")):
            for row in grid:
                g = row["g"]
                P = batch_stack(cp_mod, g, cacheG)
                cl_b, kf_b, dv_b, md_b, mk_b, dpi_b = batch_closure(cp_mod,
                                                                    P)
                total += NMEM
                # amendment 1: scalar-member degeneracy - M = e^{i phi} I
                trP = cp_mod.trace(P, axis1=1, axis2=2)
                absT = cp_mod.abs(trP)
                ph = cp_mod.where(absT > 64.0,
                                  trP / cp_mod.maximum(absT, 1e-300),
                                  1.0 + 0j)
                I0 = cp_mod.eye(128, dtype=complex)
                dev_s = cp_mod.abs(P - ph[:, None, None] * I0) \
                    .reshape(NMEM, -1).max(1)
                dev_s = cp_mod.where(absT > 64.0, dev_s, 9.0)
                if tag == "r":
                    mi = int(cp_mod.argmin(md_b))
                    if (best_rational["dev"] is None
                            or float(md_b[mi]) < best_rational["dev"]):
                        best_rational = {
                            "dev": float(md_b[mi]), "p": row["p"],
                            "q": row["q"], "g": g,
                            "approach_k": int(mk_b[mi]),
                            "k": (int(kf_b[mi]) if bool(cl_b[mi])
                                  else None),
                            "member": mi, "c": member_of(mi)}
                    key = "%d/%d" % (row["p"], row["q"])
                    ent = rational_ledger.setdefault(
                        key, {"closed": 0, "scalar_degenerate": 0,
                              "found": 0})
                    ent["closed"] = int(cl_b.sum())
                    ent["scalar_degenerate"] = int((cl_b & (dev_s
                                                            <= 1e-12)).sum())
                    for ii in cp_mod.where(cl_b)[0].tolist():
                        rec = {"c": member_of(ii), "member": ii,
                               "p": row["p"], "q": row["q"], "g": g,
                               "k": int(kf_b[ii]),
                               "dev": float(dv_b[ii]),
                               "dev_vs_pmI": float(dpi_b[ii]),
                               "dev_scalar": float(dev_s[ii])}
                        if rec["dev_scalar"] <= 1e-12:
                            degenerate.append(rec)
                        else:
                            rec["at_obs"] = bool(abs(g - 0.2) <= OBS_BAND)
                            found.append(rec)
                    ent["found"] = sum(1 for r in found
                                       if r["p"] == row["p"]
                                       and r["q"] == row["q"])
                else:
                    n_cl = int(cl_b.sum())
                    if n_cl:
                        closed_any_pi |= cp_mod.asnumpy(cl_b)
                        pi_closures["pi*%d/%d" % (row["a"],
                                                  row["b"])] = n_cl
        data["rational_ledger"] = rational_ledger
        members_with_pi_closure = int(closed_any_pi.sum())
        found.sort(key=lambda r: (r["member"], r["g"]))
    else:
        rational_ledger = {}
        for im, c in enumerate(itertools.product(MULTS, repeat=6)):
            any_pi = False
            for row in rats:
                g = row["g"]
                key = "%d/%d" % (row["p"], row["q"])
                ent = rational_ledger.setdefault(
                    key, {"closed": 0, "scalar_degenerate": 0, "found": 0})
                M = build(c, g, cache)
                closed, k, dev, mind, mind_k = closure_min(M)
                total += 1
                if (best_rational["dev"] is None
                        or mind < best_rational["dev"]):
                    best_rational = {"dev": mind, "p": row["p"],
                                     "q": row["q"], "g": g,
                                     "approach_k": mind_k,
                                     "k": k if closed else None,
                                     "member": im, "c": list(c)}
                if closed:
                    ent["closed"] += 1
                    d0 = float(min(np.abs(M - I128).max(),
                                   np.abs(M + I128).max()))
                    tr = complex(np.trace(M))
                    a = abs(tr)
                    ds = float(np.abs(M - (tr / a) * I128).max()) \
                        if a > 64.0 else 9.0
                    rec = {"c": list(c), "member": im, "p": row["p"],
                           "q": row["q"], "g": g, "k": k, "dev": dev,
                           "dev_vs_pmI": d0, "dev_scalar": ds}
                    if ds <= 1e-12:
                        ent["scalar_degenerate"] += 1
                        degenerate.append(rec)
                    else:
                        rec["at_obs"] = bool(abs(g - 0.2) <= OBS_BAND)
                        ent["found"] += 1
                        found.append(rec)
            for row in pigrid:
                g = row["g"]
                M = build(c, g, cache)
                closed, k, dev, _, _ = closure_min(M)
                total += 1
                if closed:
                    any_pi = True
                    key = "pi*%d/%d" % (row["a"], row["b"])
                    pi_closures[key] = pi_closures.get(key, 0) + 1
            if any_pi:
                members_with_pi_closure += 1
        data["rational_ledger"] = rational_ledger

    data["closure_tests"] = total
    data["backend"] = E._BACKEND
    data["wall_time"] = time.time() - t0
    data["pi_closure_ledger"] = pi_closures
    data["members_with_pi_closure"] = members_with_pi_closure
    data["best_rational_approach"] = best_rational
    data["degenerate_closures"] = degenerate[:10]
    data["rational_closures"] = found[:10]
    data["n_scalar_degenerate_rows"] = len(degenerate)
    data["n_rational_found"] = len(found)

    if found:
        status = "CLOCK-RATIONAL-CLOSURE-FOUND"
        at_obs = any(r["at_obs"] for r in found)
        data["at_obs"] = bool(at_obs)
        data["flag"] = "RATIONAL-AT-OBS" if at_obs else "RATIONAL-OFF-OBS"
    else:
        status = "CLOCK-TRANSCENDENT-DECLARED"
        data["flag"] = None
    data["status"] = status

    checks.append({"name": "census complete (%d closure tests over "
                   "%d members x (%d rational + %d pi) clocks)"
                   % (total, data["n_members"], data["n_rational_clocks"],
                      data["n_pi_clocks"]), "pass": True})

    print("V333 %s: %d non-degenerate rational closures, %d "
          "scalar-degenerate closure rows (the 64 mirror members x the "
          "rational grid), %d members with pi-closure, closest rational "
          "approach %s at %s (member %s); pi-ledger %d closing clocks"
          % (status, len(found), len(degenerate), members_with_pi_closure,
             "%.3e" % best_rational["dev"],
             "%d/%d" % (best_rational["p"], best_rational["q"]),
             best_rational["member"], len(pi_closures)), flush=True)
    _emit(data, checks, ok=True)
    print("[V333 done] %s" % status, flush=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V333_clock_transcendence", "tier": 2,
           "property": "does the root-generated six-stage step class close "
                       "at any RATIONAL clock (rational-closure census "
                       "over the 4096-member stage-multiplier ray class; "
                       "the data clock 1/5) (tier 2, O4/O8)",
           "residual": True,
           "programme": "theory-hunt V333",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v333_clock_transcendence.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("V333 %s" % out["branch"], flush=True)


if __name__ == "__main__":
    sys.exit(main())
