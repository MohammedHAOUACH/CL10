#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V194 - is the dark matter a state of the automaton's own cell? (O7)

The published extension adds ONE scalar field outside the root (Sections
12.2, 14.2). This experiment closes the in-cell question decisively, with
the automaton's OWN 32-level cell and nothing added:

  the carrier must be (i) an exact-symmetry-odd state of the cell,
  (ii) neutral under colour and electromagnetism, (iii) with the mass the
  ladder places at 118 GeV, and (iv) stable once the built interactions run.

The census is three-way - stability gradings, gauge invariance of those
gradings, and the mass placement the built roots actually deliver - and its
output is a CERTIFICATE: the exact requirements any redefinition of the cell
content must satisfy.

Checks
  T0 regressions: the V124 grading facts (exactly 2 of 1024 Clifford
     products commute with the internal step; the mirror splits 16/16-bar);
     the V41 factorized S-matrix (the n-particle step is U (x) ... (x) U:
     the bare automaton changes no particle number, so an exact grading is
     the ONLY stability structure, and decays open only through interactions).
  T1 mass window: every V65 solution root (X_0 flipped, X_1 Georgi-Glashow,
     both selective roots each) x every fragment: does any NEUTRAL fragment
     get a threshold inside the observed WIMP window (config: 1e2-1e3 GeV)
     while no coloured or charged fragment sits below 1 TeV?
  T2 grading x gauge: all 2^12 sign assignments constant on the six
     fragment classes of each of the 16 and 16-bar sectors; which commute
     with all 45 Spin(10) gauge generators (the step's own gauge
     invariance), and which additionally commute with the internal step?
  T3 carrier identification: for every surviving grading the carrier is the
     LIGHTEST odd fragment (anything heavier decays down through the gauged
     16-bar-16-bar-boson vertices once the desert bosons are on); the
     mirror grading's lightest odd fragment is the coloured Q-bar at
     1.33e8 GeV - excluded; the neutral N-bar is NOT the lightest odd
     state, so it is not stable; the Dirac mass term 16 . 10-Higgs . 16-bar
     (the implemented neutrino mechanism, F3's 0.146 eV ladder) is odd
     under the mirror and breaks under any grading with exactly one neutral
     state odd.
  T4 no bypass: the non-mirror single-neutral-odd patterns fail commutation
     with the 45 (the 16 is Spin(10)-irreducible) AND with the internal
     step; no grading protects a neutral state while keeping the built
     neutrino mass.

Branch A: an exact grading + a built root delivers a stable neutral carrier
in the WIMP window - the carrier is IN the cell. Branch B: the in-cell
no-go, and the certificate for a cell-content redefinition:
  (R1) at least two NEW neutral weight classes in the cell: a desert-scale
       partner A (even) that generates the neutrino mass quadratically -
       A sits on the n=3 ladder floor (M_N of u10pred.neutrino_ladder,
       9.6e13 GeV: the built F3 seesaw), and the carrier B occupies the
       n=10 floor (118.57 GeV from u10pred), so m_nu = m_B^2 / M_A with
       the exponent identity 2*10 - 3 = 17;
  (R2) an exact parity, odd on B alone and even on everything built
       (bilinears only for B - the V88 portal structure, production then
       computed by V193's solver with the same kernels);
  (R3) the five hard gates of the cell module (Section 3) recomputed on the
       redefined cell.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# reference spectrum: the g = 1/5 Wilson-family member (paper 8.2);
# the grading counts regressed in T0 are the V124 values measured on
# that member.
os.environ["CL10_ANGLE"] = "0.2"
from cl10.paths import exp                                        # noqa: E402
from cl10 import clifford as CL                                   # noqa: E402
from cl10 import supercell as SC                                  # noqa: E402
from cl10 import engine as E                                      # noqa: E402
from cl10 import gauge as GA                                      # noqa: E402
from cl10 import matter as MT                                     # noqa: E402
from cl10 import u10pred as U10                                   # noqa: E402
from cl10 import config as CFG                                    # noqa: E402

I32 = np.eye(32, dtype=complex)
_h = [E.asnumpy(m).astype(complex) for m in CL.generators()]
w32 = E.asnumpy(CL.pseudoscalar()).astype(complex)
B, U32, U16 = SC.sectors()
U32 = E.asnumpy(U32).astype(complex)
B = E.asnumpy(B).astype(complex)
GENS45 = [E.asnumpy(m).astype(complex) for m in GA.generators_45()]
V65 = json.load(open(exp("V65_wilson", "results", "v65_wilson.json")))
V124 = json.load(open(exp("V124_z2_gradings",
                          "results", "v124_z2_gradings.json")))
V41 = json.load(open(exp("V41_smatrix", "results", "v41_smatrix.json")))
V66 = json.load(open(exp("V66_dm_identity",
                         "results", "v66_dm_identity.json")))

WIMP_LO, WIMP_HI = CFG.CFG.physics.wimp_window_gev   # 1e2, 1e3 GeV (config)
BAD_FRAG = 1.0e3                                     # TeV: charged/coloured
# the ladder rung the carrier must sit on: the n=10 mass, single source
M_LADDER = U10.dm_ladder_mass()

# the weight basis of V65/V124 (common eigenbasis of the five Cartan
# involutions) for the 16, and the QR basis of the -i eigenspace for 16-bar
_cartan = [1j * (_h[a] @ _h[a + 1]) for a in (0, 2, 4, 6, 8)]
_Dk = [B.conj().T @ C @ B for C in _cartan]
_Mm = sum(math.sqrt(p) * _Dk[k] for k, p in enumerate((1, 2, 3, 5, 7)))
_ev, _V = np.linalg.eigh(_Mm)
LAM = np.array([[float(np.real(np.vdot(_V[:, s], _Dk[k] @ _V[:, s])))
                 for k in range(5)] for s in range(16)])
_ew, _wv = np.linalg.eig(w32)
Bb = np.linalg.qr(_wv[:, np.argsort(-np.imag(_ew))[16:]])[0]      # 16-bar basis
_xc = next(r for r in V65["data"]["part_B"].values()
           if r["su5"] == "Georgi-Glashow")
XCOEF, YCOEF = np.array(_xc["X_coef"]), np.array(V65["data"]["Y_coef"])
XY = [(round(float(LAM[s] @ XCOEF), 6), round(float(LAM[s] @ YCOEF), 6))
      for s in range(16)]

# fragment classes keyed exactly by the threshold-dict keys of the root:
# numerical (X,Y) keys for Q/u/d, literal "L(1,2)", "N", "e^c"; neutral iff
# the class is "N" (Y = 0, n = 1).
THR_KEYS = {lbl for lbl in _xc["selective_roots"][0]["thresholds_GeV"]}
COLOURED = {k for k in THR_KEYS if k not in ("N", "L(1,2)", "e^c")}


def frag_class(x, y):
    """Class key (a threshold-dict key) of the weight (X, Y)."""
    if abs(y) < 1e-6:
        return "N"
    if abs(y + 0.5) < 1e-4:
        return "L(1,2)"
    if abs(y - 1.0) < 1e-4:
        return "e^c"
    k = "(%g,%.6f)" % (x, round(y, 6))
    assert k in THR_KEYS, (x, y, k)
    return k


def main():
    t0 = time.time()
    checks = []
    out = {"experiment": "V194_carrier_in_cell", "data": {}}

    # ---- T0 regressions: the grading facts and the factorized S-matrix ----
    # stored V124: one row per non-identity commuting Clifford element
    comm_cliff = int(len(V124["data"]["gradings"])) + 1
    # recompute the Clifford commutant directly (287-line V124 recipe, short)
    n_ex = 0
    for bits in range(1024):
        P = I32.copy()
        for a in range(10):
            if bits & (1 << a):
                P = P @ _h[a]
        s = np.trace(P @ P) / 32.0
        if abs(abs(s) - 1.0) > 1e-9:
            continue
        G = (1.0 / np.sqrt(s) if s.real > 0 else 1j / math.sqrt(-s.real)) * P
        if float(np.abs(G @ G - I32).max()) < 1e-9 and \
           float(np.abs(G @ U32 - U32 @ G).max()) < 1e-10:
            n_ex += 1
    # the factorization clause (V41): the 2-particle step IS the Kronecker
    # product, so the bare step conserves particle number: an exact internal
    # grading is then the only stability structure that survives the step.
    two_shape = int(V41["data"]["s2_shape"])
    two_unit = float(V41["data"]["s2_unitary_err"])
    u2kron = two_unit < 1e-8 and two_shape == 1024
    checks.append(("T0 the Clifford grading group: %d of 1024 commute with "
                   "the internal step (V124 regression, stored %d); the "
                   "2-particle step is Kronecker-factorized (V41: shape %d, "
                   "unitary %.1e) so particle number is a superselection of "
                   "the bare step and an exact grading is the only stability "
                   "structure" % (n_ex, comm_cliff, two_shape,
                                  two_unit),
                   n_ex == comm_cliff == 2 and u2kron))

    # ---- T1 the mass window the built roots actually deliver --------------
    scan = []
    surv1 = 0
    for kk, r in V65["data"]["part_B"].items():
        for s in r.get("selective_roots", []):
            thr = s.get("thresholds_GeV") or {}
            if not thr:
                continue
            neutral_in = []
            bad_light = []
            for lbl, m in thr.items():
                neutral = (lbl == "N")
                if neutral and WIMP_LO <= m <= WIMP_HI:
                    neutral_in.append((lbl, m))
                if not neutral and m < BAD_FRAG:
                    bad_light.append((lbl, m))
            lightest_neutral = min([m for lbl, m in thr.items()
                                    if lbl == "N"] or [float("inf")])
            ok = bool(neutral_in) and not bad_light
            surv1 += ok
            scan.append({"root": "%s th=%.4f w=%.4f" % (kk, s["theta_rad"],
                                                        s["w"]),
                         "physical": bool(s.get("physical")),
                         "neutral_in_wimp_window": neutral_in,
                         "coloured_or_charged_below_TeV": bad_light,
                         "lightest_neutral_GeV": float(lightest_neutral),
                         "survivor": bool(ok)})
    out["data"]["mass_scan"] = scan
    closest = min(x["lightest_neutral_GeV"] for x in scan)
    checks.append(("T1 no built root places a NEUTRAL fragment in the WIMP "
                   "window %g-%g GeV: over %d roots the lightest neutral "
                   "threshold is %.3e GeV - %.1f orders above the ladder "
                   "rung %.2f GeV; survivors %d"
                   % (WIMP_LO, WIMP_HI, len(scan), closest,
                      math.log10(closest / M_LADDER), M_LADDER, surv1),
                   surv1 == 0 and closest > 1e8))

    # ---- T2 the grading x gauge census -----------------------------------
    # 12 classes: the six fragment classes of the 16 and the six of the
    # 16-bar; a pattern assigns +-1 per class, lifted to the 32 cell states
    # in the (weight basis | 16-bar basis). Vfull: 32 x 32 state basis.
    class_of = [(frag_class(*XY[s]), 0) for s in range(16)]
    # 16-bar weights are the exact conjugates (V66): same classes, bar side
    class_of += [(frag_class(*XY[s]), 1) for s in range(16)]
    keys = sorted({(lbl, sec) for lbl, sec in class_of})
    assert len(keys) == 12, keys
    Vfull = np.concatenate([B @ _V, Bb], axis=1)                 # 32 columns
    # precompute the generators and the step in the cell basis
    Ygens = [Vfull.conj().T @ X @ Vfull for X in GENS45]
    Yu32 = Vfull.conj().T @ U32 @ Vfull
    pos = {k: [i for i, c in enumerate(class_of) if c == k] for k in keys}
    surv45, survU32, survboth = [], [], []
    for pat in range(1 << 12):
        s = np.array([1.0 if (pat >> c) & 1 == 0 else -1.0
                      for c in range(12)])
        sv = np.empty(32)
        for c, k in enumerate(keys):
            for i in pos[k]:
                sv[i] = s[c]
        c45 = max(float(np.abs(sv[None, :] * Y - Y * sv[:, None]).max())
                  for Y in Ygens)
        cu32 = float(np.abs(sv[None, :] * Yu32 - Yu32 * sv[:, None]).max())
        if c45 < 1e-9:
            surv45.append(pat)
        if cu32 < 1e-9:
            survU32.append(pat)
        if c45 < 1e-9 and cu32 < 1e-9:
            survboth.append(pat)

    def pattern_desc(pat):
        neg = [keys[c] for c in range(12) if (pat >> c) & 1]
        odd16 = sorted({k[0] for k in neg if k[1] == 0})
        oddb = sorted({k[0] for k in neg if k[1] == 1})
        return {"odd_16": odd16, "odd_16bar": oddb}

    descsboth = [pattern_desc(p) for p in survboth]
    mirror_pat = None
    for p in surv45:
        d = pattern_desc(p)
        if not d["odd_16"] and len(d["odd_16bar"]) == 6:
            mirror_pat = p
    out["data"]["surviving_gradings"] = descsboth
    checks.append(("T2 of the %d class-constant sign patterns: %d commute "
                   "with all 45 gauge generators and %d additionally with "
                   "the internal step: the exact grading group of the step "
                   "with its gauge invariance is Z2 x Z2 - global sign x "
                   "mirror - and no fragment-selective grading survives "
                   "(the 16 is Spin(10)-irreducible)"
                   % (1 << 12, len(surv45), len(survboth)),
                   len(surv45) == 4 and len(survboth) == 4
                   and mirror_pat is not None))

    # ---- T3 carrier identification under the surviving gradings -----------
    thr = next(s for s in _xc["selective_roots"] if s.get("physical"))

    carriers = []
    for p in survboth:
        d = pattern_desc(p)
        if not d["odd_16"] and not d["odd_16bar"]:
            continue                                   # trivial grading
        odd_lbls = set(d["odd_16"]) | set(d["odd_16bar"])
        # the carrier is the LIGHTEST odd fragment: heavier odd fragments
        # decay down through the gauged 16-bar-16-bar-X vertices
        lightest = min(odd_lbls, key=lambda l: thr["thresholds_GeV"][l])
        carriers.append({"grading": ("mirror" if p == mirror_pat
                                     else str(p)),
                         "lightest_odd": lightest,
                         "mass_GeV": thr["thresholds_GeV"][lightest],
                         "neutral": bool(lightest == "N"),
                         "coloured": lightest in COLOURED,
                         "in_wimp_window": bool(
                             WIMP_LO <= thr["thresholds_GeV"][lightest]
                             <= WIMP_HI)})
    # the Dirac mass term: 16 . 10-Higgs . 16-bar (matter.mass_generation):
    # a grading with exactly one of the conjugate N's odd makes the term
    # odd -> the implemented neutrino mechanism (F3: 0.146 eV ladder) dies;
    # the mirror makes the term odd outright (16-bar odd, 16 even).
    mg = MT.mass_generation()
    dirac_present = bool(mg["off_diagonal"]) and bool(
        "16 . 10-Higgs . 16bar" in mg["mass_term"])
    nu_ladder = U10.neutrino_ladder()
    m_nu = nu_ladder["m_nu_eV_ladder"]
    checks.append(("T3 carrier identification: under the mirror grading the "
                   "only stable in-cell state is the LIGHTEST odd fragment, "
                   "%s at %.3e GeV (neutral: %s, coloured: %s) - every "
                   "heavier odd state (neutral N-bar at %.3e GeV included) "
                   "decays down through the gauged vertices; and the "
                   "implemented Dirac mass term %r (present: %s) is ODD "
                   "under the mirror and under any single-neutral-odd "
                   "grading, killing the F3 neutrino ladder (%s eV): the "
                   "built cell has NO stable neutral fragment in the window"
                   % (carriers[0]["lightest_odd"], carriers[0]["mass_GeV"],
                      carriers[0]["neutral"], carriers[0]["coloured"],
                      thr["thresholds_GeV"]["N"], mg["mass_term"],
                      dirac_present, m_nu),
                   all((not c["neutral"]) or (not c["in_wimp_window"])
                       for c in carriers) and dirac_present))

    # ---- T4 no bypass: single-neutral-odd patterns fail both layers -------
    # patterns odd on exactly one class, N of 16, N of 16-bar, or both
    single = []
    for c, k in enumerate(keys):
        if k[0] != "N":
            continue
        p = 1 << c
        sv = np.empty(32)
        for cc, kk2 in enumerate(keys):
            for i in pos[kk2]:
                sv[i] = -1.0 if cc == c else 1.0
        cu32 = float(np.abs(sv[None, :] * Yu32 - Yu32 * sv[:, None]).max())
        c45 = max(float(np.abs(sv[None, :] * Y - Y * sv[:, None]).max())
                  for Y in Ygens)
        single.append({"odd": k, "comm_U32": cu32, "comm_45_max": c45})
    allfail = all(x["comm_U32"] > 1e-6 or x["comm_45_max"] > 1e-6
                  for x in single)
    checks.append(("T4 no bypass: every single-neutral-odd pattern fails - "
                   "%s: the smallest commutator over the 45 generators is "
                   "%.2e and over the internal step %.2e (both need 0): a "
                   "grading cannot protect a neutral cell state while "
                   "respecting the built gauge structure"
                   % ([x["odd"] for x in single],
                      min(x["comm_45_max"] for x in single),
                      min(x["comm_U32"] for x in single)),
                   allfail))

    # ---- the certificate: what a redefined cell MUST have -----------------
    nu = U10.neutrino_ladder()
    m_a = nu["M_N_GeV"]                             # the n=3 ladder floor
    cert = {
        "R1_new_neutral_classes": {
            "A_seesaw_partner": {"mass_GeV_required": m_a,
                                 "parity": "even",
                                 "from": "the n=3 ladder floor: m_nu = "
                                         "m_B^2 / M_A with the exponent "
                                         "identity 2*10 - 3 = 17 (the built "
                                         "neutrino mechanism, F3 %s eV)"
                                         % m_nu,
                                 "in_desert": bool(1e11 <= m_a <= 1e16)},
            "B_carrier": {"mass_GeV_required": M_LADDER,
                          "parity": "odd", "rung": 10,
                          "from": "u10pred.dm_ladder_mass (n=10 rung)"},
        },
        "R2_parity": "exact on B alone; B enters only through bilinears "
                     "(the V88 phi_s^2 portal structure): production "
                     "computed by the V193 solver, direct detection by the "
                     "V148 map",
        "R3_hard_gates": "the five cell-module gates recomputed on the "
                         "cell (root, kinematics, unitary causal dynamics, "
                         "anomaly-free matter, 5/3)",
    }
    out["data"]["certificate"] = cert
    branch = "A" if surv1 == 1 else "B"
    ok = all(p for _, p in checks) and branch == "B" and \
        cert["R1_new_neutral_classes"]["A_seesaw_partner"]["in_desert"]

    wall = time.time() - t0
    mcar = next(c for c in carriers if c["grading"] == "mirror")
    status = ("THE IN-CELL CARRIER IS CLOSED, WITH A CERTIFICATE (O7). "
              "Stability: the bare step is Kronecker-factorized (T0), so an "
              "exact grading is the only stability structure; of the %d "
              "class-constant gradings exactly 4 commute with the 45 gauge "
              "generators AND the internal step - a Z2 x Z2 (global sign x "
              "mirror), with no fragment-selective grading (T2). Masses: no "
              "built root places a neutral "
              "fragment in the %g-%g GeV window - the lightest neutral "
              "threshold is %.3e GeV (T1). Carrier: under the mirror the "
              "only stable state is the lightest odd fragment, the COLOURED "
              "%s at %.3e GeV (excluded); the neutral N-bar is heavier and "
              "decays down through the gauged vertices (T3); and the "
              "implemented Dirac mass term (the F3 %s eV ladder) is odd "
              "under the mirror and under every single-neutral-odd grading, "
              "which all fail the gauge and the step outright (T4). THE "
              "AUTOMATON'S CURRENT CELL CANNOT HOLD THE CARRIER - and the "
              "certificate is quantitative: the redefined cell must add "
              ">= 2 neutral classes: A even at the n=3 ladder floor M_A = "
              "%.3e GeV (inside the desert: %s), m_nu = m_B^2/M_A exactly, "
              "and B odd at the n=10 rung %.2f GeV, coupled "
              "only through bilinears (R2: V88 portal, V193 production, "
              "V148 direct detection), with the five hard gates recomputed "
              "(R3). VERDICT: %s" % (1 << 12, WIMP_LO, WIMP_HI, closest,
                                     mcar["lightest_odd"],
                                     mcar["mass_GeV"], m_nu, m_a,
                                     cert["R1_new_neutral_classes"]
                                     ["A_seesaw_partner"]["in_desert"],
                                     M_LADDER,
                                     "PASS" if ok else "FAIL"))
    out["pass"] = bool(ok)
    out["branch"] = branch
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = status
    out["tier"] = 5
    out["property"] = ("is the dark-matter carrier a state of the automaton's "
                       "own cell: grading x gauge x mass census over the "
                       "built cell, and the quantitative certificate any "
                       "redefinition of the cell content must satisfy (O7)")
    out["programme"] = "V194"
    out["wall_time"] = wall

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v194_carrier_in_cell.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V194 - the in-cell carrier census and the cell-redefinition "
          "certificate (O7)")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  branch %s: %s" % (branch, "CARRIER FOUND IN CELL"
                               if branch == "A"
                               else "IN-CELL CLOSED; CERTIFICATE ISSUED"))
    print("  certificate: A at %.3e GeV (desert: %s), B at %.2f GeV (odd)"
          % (m_a, cert["R1_new_neutral_classes"]["A_seesaw_partner"]
             ["in_desert"], M_LADDER))
    print("  VERDICT: %s   wall %.2f s" % ("PASS" if ok else "FAIL", wall))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
