#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V376 - CP-odd census of the dressed class (booked BEFORE the run
in review/pre_registrations_root.md lines 3007-3032).

Question: does the V52 theorem "every CP-odd invariant of the
published QCA vanishes" extend from the bare step U to the DRESSED
class { W U, e^{i eta sum_b J_b.L_b} U, and the V373 full generator
} - is the dressing CP-even to measured exactness, or does it carry
CP-odd content the bare QCA provably lacks?

Machine (booked): the CP-odd invariant machinery is RE-USED, not
reinvented - cl10/parameter_audit.cp_odd_theorem's invariant list
(Im Tr of every 2- and 3-rotation product full, the 4-rotation
sample n4_sample = 200000 seed 0, every B-violating h_i R_ab h_j
R_cd sandwich full, and the rank-1 real democratic Yukawa
invariant) and the CP operator from the SAME source (pseudoscalar
w, the Cl(10,0) volume element, w^2 = -1, w* = -w). The dressed
step replaces U in every invariant; all evaluated at eta in {0,
0.02, ..., 1} with machine-exact traces verified by mpmath at 40
digits; decision threshold V52 TOL = 1e-12.

Class member identifications (gated, declared):
  W U: W is the published dressed hop of V251, W = v101_gap * G10
    (V251's omega rows ARE +-W; recomputed here from the V101 and
    coupling artifacts - gate). The step: the same 45 plane
    rotations at the dressed angle W (_gut_step form).
  e^{i eta sum_b J_b.L_b} U: on the certified internal 32-dim
    space the bond-current operator is the MEASURED scalar
    sum_a T_a T_a = C2 * I (C2 = 45, V251; re-measured here as a
    gate), so the dressing multiplies U by a global phase and the
    45 rotation factors are unchanged - declared, not assumed.
  V373 full generator: its internal-space action is the same
    measured scalar dressing (flux string terms act outside the
    internal invariant space; H_hop is a lattice momentum term) -
    declared reduction; the certified CHANNEL-space spectra
    (V371/V373 dk) are recorded as measured data with an explicit
    caveat: channel-space relative phases are not invariants of
    the booked internal class.
Verdicts: CP-EVEN-DRESSED-CLASS-CLOSED / CP-ODD-FOUND-DRESSED.
Scope booked: extends the theorem from the bare root step to the
booked dressed class; NOT a claim that no CP-odd operator exists
anywhere in the lattice theory (the named eta_CP residual of
Table 6 stands).
"""
import json
import math
import os
import sys
import time

import mpmath
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
from cl10 import coupling  # noqa: E402
from cl10 import parameter_audit as PA  # noqa: E402
from cl10.paths import exp  # noqa: E402

TWO_PI = 2.0 * np.pi
T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def mp_matrix(m):
    return mpmath.matrix([[mpmath.mpc(complex(x)) for x in row]
                          for row in m])


def census(fac, rng_seed=0, n4=200000):
    """the booked invariant list on a list of rotation factors
    (float64; the mpmath verification tier is separate)."""
    n = len(fac)
    maxJ2 = 0.0
    for a in range(n):
        for b in range(n):
            maxJ2 = max(maxJ2, abs(float(np.imag(
                np.trace(fac[a] @ fac[b])))))
    maxJ3 = 0.0
    for a in range(n):
        A2 = fac[a]
        for b in range(n):
            P2 = A2 @ fac[b]
            for c in range(n):
                maxJ3 = max(maxJ3, abs(float(np.imag(
                    np.trace(P2 @ fac[c])))))
    rng = np.random.default_rng(rng_seed)
    maxJ4 = 0.0
    for _ in range(n4):
        a, b, c, d = (int(x) for x in rng.integers(0, n, 4))
        maxJ4 = max(maxJ4, abs(float(np.imag(np.trace(
            fac[a] @ fac[b] @ fac[c] @ fac[d])))))
    return maxJ2, maxJ3, maxJ4


def sandwich_max(fac, h):
    """every B-violating h_i R_ab h_j R_cd sandwich, full."""
    n = len(fac)
    maxJb = 0.0
    for i in range(10):
        for k in range(n):
            Aik = h[i] @ fac[k]
            for j in range(10):
                Aj = Aik @ h[j]
                for l in range(n):
                    maxJb = max(maxJb, abs(float(np.imag(np.trace(
                        Aj @ fac[l])))))
    return maxJb


def yukawa_invariant():
    Y = np.ones((3, 3)) / 3.0
    YYd = Y @ Y.conj().T
    inv = YYd @ YYd
    return abs(float(np.imag(np.trace(inv))))


def mp_verify_samples(fac, mf, h, mh, rng, nsamp):
    """mpmath 40-digit verification tier: the SAME random 2/3/4-
    rotation products and sandwiches evaluated in float64 and in
    mpmath; returns (max float-vs-mp deviation, max |Im Tr| at
    mpmath precision)."""
    mpmath.mp.dps = 40

    def tr(M):
        return sum(M[i, i] for i in range(M.rows))

    dev = 0.0
    worst = 0.0
    n = len(fac)
    for _ in range(nsamp):
        pick = [int(x) for x in rng.integers(0, n, 4)]
        L = int(rng.integers(2, 5))
        Pf = fac[pick[0]].copy()
        Pm = mf[pick[0]].copy()
        for q in range(1, L):
            Pf = Pf @ fac[pick[q]]
            Pm = Pm * mf[pick[q]]
        fi = float(np.imag(np.trace(Pf)))
        mi = float(mpmath.im(tr(Pm)))
        dev = max(dev, abs(fi - mi))
        worst = max(worst, abs(mi))
        i, k, j, l = (int(rng.integers(0, 10)), int(rng.integers(0, n)),
                      int(rng.integers(0, 10)), int(rng.integers(0, n)))
        Sf = h[i] @ fac[k] @ h[j] @ fac[l]
        Sm = mh[i] * mf[k] * mh[j] * mf[l]
        fi = float(np.imag(np.trace(Sf)))
        mi = float(mpmath.im(tr(Sm)))
        dev = max(dev, abs(fi - mi))
        worst = max(worst, abs(mi))
    return dev, worst


def main():
    checks, data = [], {}

    # ------------------------------------------------------------
    # machinery re-use gates: the CP operator and the bare anchor
    # ------------------------------------------------------------
    U, R, planes = PA._gut_step()
    Rn = [PA._asnp(r) for r in R]
    h = [PA._asnp(g) for g in PA._gen()]
    w = h[0]
    for a in h[1:]:
        w = w @ a
    w2def = float(np.abs(w @ w + np.eye(32)).max())
    antidef = max(float(np.abs(w @ a + a @ w).max()) for a in h)
    # the booked w* = -w is the phase convention of the pre-rebuild
    # representation; the current published 32x32 Pauli source (the
    # SAME source the theorem machinery runs on, fold 285) carries a
    # REAL pseudoscalar, w* = +w. The sign convention is recorded,
    # the basis-invariant content (w^2 = -1, anticommutation with
    # all 10 Cl(10,0) generators) is gated at machine exactness; the
    # booked census invariants do not involve w at all.
    wstar = float(np.abs(w.conj() - w).max()) < 1e-12
    data["w_star_convention_measured"] = "+w (current Pauli basis)"\
        if wstar else "-w"
    a_w = w2def <= 1e-12 and antidef <= 1e-12
    checks.append(("CP operator from the same source: w^2 = -1 "
                   "(defect %.2e) and w anticommutes with all 10 "
                   "generators (defect %.2e), both <= 1e-12; w* "
                   "sign convention recorded as data (%s, the "
                   "booked -w is the pre-rebuild basis phase): %s"
                   % (w2def, antidef, data["w_star_convention_"
                                          "measured"], a_w), a_w))
    TOL = float(PA._TOL.strict)
    a_tol = TOL <= 1e-12
    checks.append(("V52 TOL = %.1e <= booked 1e-12 decision "
                   "threshold: %s" % (TOL, a_tol), a_tol))

    v52 = json_art("V52_parameter_audit")
    c52 = v52["data"]["cp_odd_theorem"]

    # bare eta = 0 column: same machine, must reproduce V52
    m2, m3, m4 = census(Rn, n4=int(c52["n4_sample"]))
    mb = sandwich_max(Rn, h)
    yk = yukawa_invariant()
    bare = {"2rot": m2, "3rot": m3, "4rot_sample": m4,
            "sandwich": mb, "yukawa": yk}
    a_bare = (
        abs(m2 - float(c52["max_cpodd_2rot"])) <= 1e-12
        and abs(m3 - float(c52["max_cpodd_3rot"])) <= 1e-12
        and abs(m4 - float(c52["max_cpodd_4rot_sample"])) <= 1e-12
        and abs(mb - float(c52["max_cpodd_b_sandwich"])) <= 1e-12
        and abs(yk - float(c52["cpodd_democratic_yukawa"])) <= 1e-12)
    checks.append(("eta = 0 column reproduces the published V52 "
                   "census exactly (max dev %.2e <= 1e-12): %s"
                   % (max(abs(bare["2rot"] - float(c52["max_cpodd_"
                                                   "2rot"])),
                          abs(bare["3rot"] - float(c52[
                              "max_cpodd_3rot"])),
                          abs(bare["4rot_sample"] - float(
                              c52["max_cpodd_4rot_sample"])),
                          abs(mb - float(c52["max_cpodd_b_"
                                             "sandwich"]))), a_bare),
                  a_bare))

    # ------------------------------------------------------------
    # class member 1: W U (W = v101_gap * G10, V251's dressed hop)
    # ------------------------------------------------------------
    v101 = json_art("V101_bond_current_ground_state")
    W = float(v101["data"]["gap"]) * float(coupling.G10)
    v251 = json_art("V251_gauged_two_particle_sector")
    om = float(v251["data"]["sensitivity"]["omega_rows"]["plus"])
    # historical gate: V251 ran under the pre-root-rebuild clock
    # (G10 = 0.2); under the current clock W recomputes - the
    # artifact omega rows gate the gap factor, not the product.
    gapfac = om / float(v101["data"]["gap"])
    data["W_dressed_hop"] = W
    data["W_historical_V251"] = om
    data["V251_omega_over_gap"] = gapfac
    a_W = abs(gapfac - 0.2) <= 1e-9
    checks.append(("W = v101_gap x G10 from artifacts; V251's "
                   "omega rows == gap x 0.2 (the then-clock G10, "
                   "recorded): current W = %.9f: %s" % (W, a_W),
                   a_W))
    cW, sW = math.cos(W), math.sin(W)
    I32 = np.eye(32, dtype=np.complex128)
    facW = [cW * I32 - 1j * sW * PA._asnp(T) for T in planes]
    unitW = 0.0
    UA = I32.copy()
    for r in facW:
        UA = UA @ r
    unitW = float(np.abs(UA @ UA.conj().T - I32).max())
    checks.append(("W U step unitary (defect %.2e <= 1e-12): %s"
                   % (unitW, unitW <= 1e-12), unitW <= 1e-12))
    w2, w3, w4 = census(facW, n4=200000)
    wb = sandwich_max(facW, h)
    member1 = {"2rot": w2, "3rot": w3, "4rot_sample": w4,
               "sandwich": wb, "yukawa": yk}

    # ------------------------------------------------------------
    # measured-scalar gate for members 2 and 3
    # ------------------------------------------------------------
    Csum = np.zeros((32, 32), dtype=np.complex128)
    for a in range(45):
        T = PA._asnp(planes[a])
        Csum += T @ T
    # planes are the so(10) plane generators; the current-current
    # contraction sum_a T_a T_a over the 45 roots equals the
    # quadratic Casimir on the spinor space - measure the value.
    c2 = complex(np.trace(Csum) / 32.0)
    scalar_def = float(np.abs(Csum - c2 * I32).max())
    a_sc = scalar_def <= 1e-12
    checks.append(("internal-space bond-current operator "
                   "sum_a T_a.T_a = C2 I with C2 = %.6f is a "
                   "MEASURED scalar (dev %.2e <= 1e-12): the "
                   "members 2/3 dressings are global phases on "
                   "this space, the 45 rotation factors unchanged: "
                   "%s" % (c2.real, scalar_def, a_sc), a_sc))
    data["internal_C2_measured"] = c2.real

    etas = [round(x, 3) for x in np.arange(0.0, 1.001, 0.02)]
    members = {"W U": member1,
               "exp(i eta J.L) U": bare,
               "V373 full generator (internal action)": bare}

    # mpmath 40-digit verification tier (booked "machine-exact
    # traces, mpmath high precision"): the SAME products evaluated
    # in both precisions. The factor sets are eta-INDEPENDENT by
    # the measured scalar gate, so the tier is computed once per
    # set and recorded for every eta - declared, not silent.
    rngV = np.random.default_rng(0x376)
    mh = [mp_matrix(x) for x in h]
    mRn = [mp_matrix(x) for x in Rn]
    mfW = [mp_matrix(x) for x in facW]
    devB, ampB = mp_verify_samples(Rn, mRn, h, mh, rngV, 150)
    devW, ampW = mp_verify_samples(facW, mfW, h, mh, rngV, 150)
    dev = max(devB, devW)
    per_eta_mp = {str(eta): {
        "exp(i eta J.L) U": ampB,
        "V373 full generator (internal action)": ampB,
        "W U": ampW} for eta in etas}
    data["mpmath_40d_per_eta_max"] = per_eta_mp
    data["mpmath_tier_note"] = ("factor sets measured eta-"
                                "independent (scalar dressing "
                                "gate); verification tier computed "
                                "once per factor set, recorded at "
                                "every booked eta")
    data["mpmath_float_bound"] = dev
    a_dev = dev <= 1e-14
    checks.append(("mpmath-40d verification vs float64 census: "
                   "deviation bound %.2e (<= 1e-14): %s"
                   % (dev, a_dev), a_dev))

    # the booked per-eta decision: invariant maxima + bound
    table = {}
    closed = True
    for eta in etas:
        row = {}
        for name, mem in members.items():
            mx = max(mem["2rot"], mem["3rot"], mem["4rot_sample"],
                     mem["sandwich"], mem["yukawa"],
                     per_eta_mp[str(eta)][name])
            row[name] = mx
            if mx + dev > 1e-12:
                closed = False
        table[str(eta)] = row
    data["per_eta_maxima"] = table

    # channel-space record (declared, NOT decisional)
    v371 = json_art("V371_deformed_contact_channels")
    v373 = json_art("V373_full_lattice_pair_string")
    dkch = ([float(x) for x in v371["data"]["D2_eigs"]]
            + [float(x) for x in v373["data"]["dk_factored"]
               ["chain_ell2"]])
    ch_im = max(abs(np.imag(np.sum(np.exp(1j * eta *
                           np.array(dkch)))))
                for eta in etas)
    data["channel_space_phase_record_max_abs_imTr"] = float(ch_im)

    okm = a_w and a_tol and a_bare and a_W and a_sc and a_dev
    if not okm:
        status = "MACHINERY-DRIFT (a machinery gate failed)"
        ok = False
    elif closed:
        status = ("CP-EVEN-DRESSED-CLASS-CLOSED (every booked "
                  "invariant |I| <= 1e-12 at every eta for all "
                  "three class members; the dressings are "
                  "measured scalars on the internal space and the "
                  "W-dressed step keeps the 0-parameter algebraic "
                  "cancellations; channel-space relative phases "
                  "recorded as data, not invariants of the booked "
                  "class; the named eta_CP residual stands)")
        ok = True
    else:
        status = ("CP-ODD-FOUND-DRESSED (a booked invariant "
                  "exceeds 1e-12 + mpmath bound at some eta; see "
                  "per_eta_maxima)")
        ok = True
    checks.append(("closed-branch gate: every member max across "
                   "all 51 etas <= 1e-12 (with verification "
                   "bound): %s" % closed, closed))

    worst_all = max(max(r.values()) for r in table.values())
    out = {"id": "V376_dressed_cp_odd_census",
           "property": ("does the CP-odd-vanishing theorem extend "
                        "from the bare step to the dressed class "
                        "{W U, e^{i eta J.L} U, V373 generator}"),
           "tier": 2,
           "checks": [[m, bool(b)] for m, b in checks],
           "data": data, "pass": ok, "status": status,
           "residual": float(worst_all)}
    with open(os.path.join(HERE, "results",
                           "v376_dressed_cp_odd_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for m, b in checks:
        print(("  ok  " if b else "  FAIL ") + m)
    print("[V376 done] %s" % status)
    log("wall %.1f s" % (time.time() - T0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
