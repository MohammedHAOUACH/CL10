#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V323 - the joint-module stabiliser against the rung-17 contact
(tier 2, hard, O8 + neutrino, theory-hunt).

**The question (pre-registered).** Does an element of the centre of
the joint-module stabiliser - a Casimir monomial built ONLY from the
booked joint spectrum {45, 85, 117} x {32, 288, 1120} (V79/V140
leaves), no target inside any member - give the rung-17 neutrino
contact |S17/3 x chi - m_atm| / m_atm inside the 5% ledger band, and
is the hit a DERIVATION (exact algebraic match to chi_req) or
supplied-by-density (V306/V319 rule)?

**Class (closed, 28 members, pre-fixed).** t = (45, 85, 117) sorted,
m = (32, 288, 1120): identity (1); t_i (3); t_i/t_j i != j (6);
t_i/(247/3) (3); t_i/109 with 109 the multiplicity-weighted mean
(3); t_i*t_j/t_k^2 distinct (3); sorted-gap set {dt32/dt21,
dt21/dt32, (dt21+dt32)/dt32}, dt21 = t2 - t1, dt32 = t3 - t2 (3);
(m_j t_j - m_i t_i)/(m_j t_j) ordered i != j (6). Exact rationals via
fractions.Fraction - members carry zero floating-point tuning.

**Reading (inherited).** leg' = S17/3 x chi (V316/V319 convention);
band = pillar_surgery.ledger_scale_rel_tol; chi_req =
m_atm/(S17/3) enters the DECISION ONLY.

**T0.** S17 re-derived (cl10.neutrino.loop_ladder(17) x EV_PER_GEV)
against the V303 leaf < 1e-9; the joint content {45, 85, 117} x
{32, 288, 1120} re-derived in-experiment with V79's own machinery
(charge C^A = T^A x I + I x L^A on 32 x 45, Q2 = sum C^A C^A,
eigvalsh clustering) against the V79 leaves; drift -> halt.

**Branch (pre-registered words).** S17-SEALED-EXACT (a member matches
chi_req to rel < 1e-9 - the -2.6962% becomes an algebraic property of
the joint module, success tier, re-provenance); S17-SEALED-IN-BAND
(band hits, none exact - graded by the chance-density count:
expected hits = N x band log width / class log range, supplier-grade
iff < 0.5); S17-STUBBORN otherwise. Full 28-row table booked in all
cases.
"""
import json
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))

import numpy as np                                        # noqa: E402

from cl10 import config as _config                        # noqa: E402
from cl10 import neutrino as nv                           # noqa: E402
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)

CFG = _config.CFG
BAND = float(CFG.pillar_surgery.ledger_scale_rel_tol)
M_ATM = float(CFG.neutrino.m_atm_ev)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def build_class():
    t = (45, 85, 117)
    m = (32, 288, 1120)
    out = []

    def add(tag, x):
        out.append((tag, Fraction(x)))

    add("identity", 1)
    for i in range(3):
        add("t%d" % (i + 1), t[i])
    for i in range(3):
        for j in range(3):
            if i != j:
                add("t%d/t%d" % (i + 1, j + 1), Fraction(t[i], t[j]))
    for i in range(3):
        add("t%d/(247/3)" % (i + 1), Fraction(3 * t[i], 247))
    for i in range(3):
        add("t%d/109" % (i + 1), Fraction(t[i], 109))
    for (i, j, k) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        add("t%d*t%d/t%d^2" % (i + 1, j + 1, k + 1),
            Fraction(t[i] * t[j], t[k] * t[k]))
    dt21, dt32 = t[1] - t[0], t[2] - t[1]
    add("dt32/dt21", Fraction(dt32, dt21))
    add("dt21/dt32", Fraction(dt21, dt32))
    add("(dt21+dt32)/dt32", Fraction(dt21 + dt32, dt32))
    mt = [m[i] * t[i] for i in range(3)]
    for i in range(3):
        for j in range(3):
            if i != j:
                add("(m%d t%d - m%d t%d)/(m%d t%d)" % (j + 1, j + 1,
                                                       i + 1, i + 1,
                                                       j + 1, j + 1),
                    Fraction(mt[j] - mt[i], mt[j]))
    return out


def main():
    t0 = time.time()
    checks, data = [], {}

    # ---- T0a: S17 -------------------------------------------------------
    S17 = float(nv.loop_ladder(17)) * float(nv.EV_PER_GEV)
    leaf = float(art("V303_junction_energy_r", "data.contact.S17_eV"))
    t0a = abs(S17 - leaf) / leaf
    checks.append({"name": "T0a S17 re-derived {:.7e} eV vs the V303 leaf "
                   "{:.7e} (rel {:.1e} < 1e-9)".format(S17, leaf, t0a),
                   "pass": bool(t0a < 1e-9)})
    data["t0_S17_rel"] = t0a

    # ---- T0b: joint content with V79's own machinery --------------------
    Ts = root_generators()
    import cl10.engine as E
    Ts = [np.array(E.asnumpy(T), dtype=complex) for T in Ts]
    f, _, _ = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    Ls = [np.array(E.asnumpy(L), dtype=complex) for L in Ls]
    I45 = np.eye(45, dtype=np.complex128)
    I32 = np.eye(32, dtype=np.complex128)
    Cs = [np.kron(T, I45) + np.kron(I32, L) for T, L in zip(Ts, Ls)]
    Q2 = sum(C @ C for C in Cs)
    ev = np.linalg.eigvalsh(Q2)
    clusters = []
    for e in ev:
        if clusters and abs(e - clusters[-1][0]) < 1e-6 * max(1.0, abs(e)):
            clusters[-1][1] += 1
        else:
            clusters.append([float(e), 1])
    levels = [round(c[0], 4) for c in clusters]
    mults = sorted(c[1] for c in clusters)
    lv_leaf = art("V79_gauged_bond_layer", "data.content.casimir_levels")
    mu_leaf = art("V79_gauged_bond_layer", "data.content.multiplicities")
    t0b = (levels == [round(float(x), 4) for x in lv_leaf]
           and mults == list(mu_leaf))
    checks.append({"name": "T0b joint content re-derived: levels {} mults "
                   "{} vs V79 leaves {} / {} - exact".format(
                       levels, mults, lv_leaf, mu_leaf),
                   "pass": bool(t0b)})
    if not (t0a < 1e-9 and t0b):
        data["status"] = "MACHINERY-DRIFT"
        print("V323 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- decision --------------------------------------------------------
    chi_req = Fraction(M_ATM) / (Fraction(repr(S17 / 3.0)))
    table = []
    exact_hits, band_hits = [], []
    for tag, chi in build_class():
        if chi <= 0:
            table.append({"tag": tag, "chi": str(chi), "status":
                          "NEGATIVE-CHI (booked, not acceptable)"})
            continue
        leg = abs(S17 / 3.0 * float(chi) - M_ATM) / M_ATM
        rel_exact = abs(float(chi) - float(chi_req)) / float(chi_req)
        st = "open"
        if rel_exact < 1e-9:
            st = "EXACT"
            exact_hits.append((tag, chi))
        elif leg <= BAND:
            st = "BAND-HIT"
            band_hits.append((tag, chi))
        table.append({"tag": tag, "chi": str(chi),
                      "chi_f": float(chi), "leg_rel": leg,
                      "normalized": leg / BAND, "status": st})
    data["chi_req"] = float(chi_req)
    data["S17_over_3_eV"] = S17 / 3.0
    data["m_atm_eV"] = M_ATM
    data["band"] = BAND
    data["table"] = table
    pos = [float(c) for _, c in build_class() if c > 0]
    import math
    log_range = math.log(max(pos) / min(pos))
    band_width = math.log(1.0 + BAND) - math.log(1.0 - BAND)
    N = len(build_class())
    exp_hits = N * band_width / log_range
    data["density"] = {"N": N, "log_range": log_range,
                       "band_log_width": band_width,
                       "expected_hits": exp_hits,
                       "supplier_grade": bool(exp_hits < 0.5)}
    if exact_hits:
        branch = "S17-SEALED-EXACT"
    elif band_hits:
        branch = "S17-SEALED-IN-BAND"
    else:
        branch = "S17-STUBBORN"
    data["status"] = branch
    data["exact_hits"] = [[t_, str(c)] for t_, c in exact_hits]
    data["band_hits"] = [[t_, str(c)] for t_, c in band_hits]
    data["wall_time"] = time.time() - t0
    checks.append({"name": "branch {}: chi_req {:.7f}; exact {} band-hits "
                   "{} {}; density N={} expected {:.2f} supplier-grade {}"
                   .format(branch, float(chi_req), len(exact_hits),
                           len(band_hits),
                           [t_ for t_, _ in band_hits], N, exp_hits,
                           exp_hits < 0.5), "pass": True})
    for r in table:
        if r.get("status") in ("EXACT", "BAND-HIT"):
            print("  {}: chi {} leg {:.4f} norm {:.3f}".format(
                r["tag"], r["chi"], r["leg_rel"], r["normalized"]),
                flush=True)
    _emit(data, checks, ok=all(c["pass"] for c in checks))
    print("[V323 done] {}".format(branch), flush=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V323_joint_stabiliser_s17", "tier": 2,
           "property": "joint-module centre monomials as rung-17 "
                       "multipliers against the atmospheric contact; "
                       "exact-derivation vs chance-density tiers "
                       "(tier 2, hard, O8 + neutrino)",
           "residual": True, "programme": "theory-hunt V323",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v323_joint_stabiliser_s17.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V323 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
