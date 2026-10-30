#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V325 - the channel triple against the certified crossing (tier 2,
hard, O8 + the C-site triple, theory-hunt).

**The question (pre-registered).** The certified root-clock crossing
(F_208 = 0.019692794505156912, certified-continuity anchor at
pi/4 + 1e-8, V321) reads the coefficient 2.710; the observed 2.0766
needs F to drop to F_req = LOOP (COEF_OBS/3)^0.1: a multiplier
mu_req = F_req/F_cert = 0.9737 on F. Does any booked 3-of-12 channel
triple (the V205/V316 W3 machinery AT THE ROOT CLOCK g = pi/4) supply
exactly that multiplier - the C-site triple used as a coefficient
supply AND a selection hit at site S3 in one question?

**Object (closed class, fixed pre-run).** f_s = |sum_j exp(i Re w_j)|
/ 128 at each of the 12 SC.M_POINTS nodes at g = pi/4
(angle.branches.closed), w = SC.phases(SC.bloch(q, U16, g10 = pi/4)) -
V316's W3 block verbatim with the clock replaced;
mu(triple) = (f_a + f_b + f_c) / (3 fbar), fbar the 12-factor mean
(V316's booked composition). All C(12,3) = 220 triples; combinatorial
provenance, nothing fitted, no target inside any member.

**Decision (targets enter only here).** CH-F208-CLOSED iff some triple
has |3 (mu F_cert / LOOP)^10 - COEF_OBS| <= coef_band_abs. Uniqueness
leg: the gap of the closing mu to the runner-up (the S3 selection
content). Density (V306 rule, razor by construction): expected hits
220 x 2 ln(1 + coef_band_abs/coef_cert) / class-log-range, booked;
any hit is supplier-grade iff expected < 0.5. CH-F208-STUBBORN
otherwise with the full mu table booked (top 20 by closeness plus the
class summary); the 30.4% gap is then bounded from the channel side.

**T0.** The 12 factors at g = 1/5 reproduce V316's booked f_nodes
leaves (data.w3_node_factors.factors) to < 1e-9; the plain rung-10
read at 1/5 reproduces V206's coef_at_fifth < 1e-9; the certified
F_208 anchor re-derives (mp_ladder at pi/4 + 1e-8, 90 dps, nmax 260,
V321's exact replay) to the V321 leaf < 1e-12.
"""
import itertools
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
from cl10.paths import exp                               # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))
sys.path.insert(0, exp("V209_root_clock_degeneracy"))

from V203_closed_branch_coefficient import hop_machinery, loop_ratios  # noqa
from cl10 import config as _config                        # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
G5 = float(CFG.angle.branches.implemented)                # 1/5, T0 clock
GP4 = float(CFG.angle.branches.closed)                    # pi/4, object
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
NMAX = 20


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def read_coef(Fv, LOOP):
    return 3.0 * (Fv / LOOP) ** 10


def node_factors(g):
    """V316's W3 block verbatim with the clock as argument."""
    _B, _U32, U16c = SC.sectors(g)
    out = {}
    for i, q in enumerate(SC.M_POINTS):
        w = SC.phases(SC.bloch(np.array(q, float), U16c, g10=g))
        out[str(i)] = float(abs(sum(np.exp(1j * float(np.real(x)))
                                    for x in w)) / 128.0)
    return out


def main():
    t0 = time.time()
    checks, data = [], {}
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)

    # ---- T0a: the 12 factors at 1/5 reproduce V316's leaves ------------
    f5 = node_factors(G5)
    f5_art = art("V316_weight_double_contact",
                 "data.w3_node_factors.factors")
    d5 = max(abs(f5[k] - float(f5_art[k]["factor"])) for k in f5)
    data["t0a_node_factors"] = {"abs_diff_max": d5}
    checks.append({"name": "T0a the 12 node factors at g = 1/5 reproduce "
                   "V316's booked f_nodes leaves (max |d| {:.1e} < 1e-9)"
                   .format(d5), "pass": bool(d5 < 1e-9)})

    # ---- T0b: plain rung-10 read at 1/5 = V206 coef_at_fifth -----------
    G, cartan, Ma, ev = hop_machinery(G5)
    F5, _rho = loop_ratios(G, cartan, Ma, [float(e) for e in ev], G5, NMAX)
    coef5_art = float(art("V206_clock_tension",
                          "data.data_clock.coef_at_fifth"))
    d5b = abs(read_coef(F5[10], LOOP) - coef5_art)
    data["t0b_coef_fifth"] = {"artifact": coef5_art, "abs_diff": d5b}
    checks.append({"name": "T0b plain rung-10 read at g = 1/5 reproduces "
                   "V206's booked coef_at_fifth {:.15f} (|d| {:.1e} < "
                   "1e-9)".format(coef5_art, d5b),
                   "pass": bool(d5b < 1e-9)})

    # ---- T0c: certified F_208 anchor, V321's exact replay --------------
    import mpmath as mp  # noqa: F401 - the expr strings are eval'd by mp_ladder
    from V209_root_clock_degeneracy import mp_ladder, skeleton
    N, blocks, dens, D, cart, planes, TL = skeleton()
    r_anchor = mp_ladder("mp.pi/4 + mp.mpf('1e-8')", 260, 90,
                         planes, blocks, dens, cart)
    cert_leaf = float(art("V321_crossing_resummation",
                          "data.certified_leaves.F_208"))
    F_cert = float(abs(r_anchor[208]["F"]))
    d5c = abs(F_cert - cert_leaf)
    data["t0c_anchor"] = {"leaf": cert_leaf, "rederived": F_cert,
                          "abs_diff": d5c}
    checks.append({"name": "T0c certified F_208 anchor at pi/4 + 1e-8 "
                   "re-derived to the V321 leaf (|d| {:.1e} < 1e-12)"
                   .format(d5c), "pass": bool(d5c < 1e-12)})

    # ---- the object: 12 factors at the root clock, 220 triples ---------
    f4 = node_factors(GP4)
    fvals = list(f4.values())
    fbar = sum(fvals) / len(fvals)
    spread = max(abs(v - fbar) for v in fvals) / fbar
    data["root_node_factors"] = {"clock": "pi/4 (angle.branches.closed)",
                                 "f_bar": fbar, "factors": f4,
                                 "spread_max_over_fbar": spread}
    coef_cert = read_coef(cert_leaf, LOOP)
    F_req = LOOP * (COEF_OBS / 3.0) ** 0.1
    mu_req = F_req / cert_leaf
    data["certified"] = {"F_cert": cert_leaf, "coef_cert": coef_cert,
                         "F_req": F_req, "mu_req": mu_req,
                         "coef_observed": COEF_OBS, "band": BAND_C}

    rows = []
    for T in itertools.combinations(range(12), 3):
        mu = sum(f4[str(i)] for i in T) / (3.0 * fbar)
        coef = read_coef(mu * cert_leaf, LOOP)
        res = abs(coef - COEF_OBS)
        rows.append({"triple": list(T), "mu": mu, "coef": coef,
                     "residual": res, "normalized": res / BAND_C,
                     "inside": bool(res <= BAND_C),
                     "mu_gap_to_req": abs(mu - mu_req)})
    rows.sort(key=lambda r: r["residual"])
    n_inside = sum(1 for r in rows if r["inside"])
    mus = [r["mu"] for r in rows]
    data["class"] = {"n_rows": len(rows), "n_inside": n_inside,
                     "mu_min": min(mus), "mu_max": max(mus),
                     "best": rows[0], "top20": rows[:20]}

    # ---- density (V306 rule, the pre-registered razor) ------------------
    logwidth = 2.0 * math.log(1.0 + BAND_C / coef_cert)
    logrange = (math.log(max(mus)) - math.log(min(mus))
                if max(mus) > min(mus) > 0 else 0.0)
    expected = (len(rows) * logwidth / logrange) if logrange > 0 else None
    data["density"] = {"band_log_width": logwidth,
                       "class_log_range": logrange,
                       "expected_hits": expected,
                       "supplier_grade_rule": "expected < 0.5"}

    # ---- verdict ---------------------------------------------------------
    if n_inside > 0:
        status = "CH-F208-CLOSED"
        gap = (rows[1]["mu"] - rows[0]["mu"]) if len(rows) > 1 else None
        data["selection_leg"] = {"n_inside": n_inside,
                                 "closing_triple": rows[0]["triple"],
                                 "gap_to_runner_up_mu": gap,
                                 "unique": bool(n_inside == 1)}
    else:
        status = "CH-F208-STUBBORN"
        data["selection_leg"] = {"n_inside": 0,
                                 "nearest_triple": rows[0]["triple"],
                                 "nearest_mu": rows[0]["mu"],
                                 "nearest_normalized": rows[0]["normalized"]}
    data["status"] = status
    data["wall_time"] = time.time() - t0

    checks.append({"name": "class census: exactly 220 triples at the root "
                   "clock ({} ran); {} inside band; nearest normalized "
                   "residual {:.4f}".format(len(rows), n_inside,
                                            rows[0]["normalized"]),
                   "pass": bool(len(rows) == 220)})

    ok = all(c["pass"] for c in checks)
    out = {"id": "V325_channel_triple_f208", "tier": 2,
           "property": "channel-triple multiplier class against the "
           "certified root-clock crossing (tier 2, hard, O8)",
           "residual": True, "programme": "theory-hunt V325",
           "branch": status, "data": data, "checks": checks,
           "status": status, "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v325_channel_triple_f208.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V325 done] {}".format(status))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
