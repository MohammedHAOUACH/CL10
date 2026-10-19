#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""V196 - the cell definition as a root-level fact (the 36-level cell).

**What this experiment pins.**  cl10/cell.py defines the automaton's cell:
the 32-D Dirac module S of Cl(10,0) plus the two neutral singlet sectors -
conjugate pairs at rungs A (3, the seesaw partner) and B (10, the dark
matter) - with the exact B-parity.  This is the axiom-level content of
Section 2 (Axiom 1), so the check sits in the algebra-root programme: the
cell is part of the root definition, and every other DM experiment (V194
census, V195 production) presupposes it.

Checks
  T0 the ladder identity between cell levels: m_A is the seesaw Majorana
     rung, m_B is the DM rung, and m_nu = m_B^2 / m_A reproduces the
     ladder neutrino mass exactly - the exponent identity 2*10-3 = 17
     becomes an identity between two levels of the cell.
  T1 the exact B-parity and the full-cell theorems: the 36-level step is
     unitary, G = diag(I32,+1,+1,-1,-1) is an exact Z2 commuting with the
     step and with all 45 plane generators extended by zero, the spectrum
     is symmetric under phase -> -phase, det U = 1 and the quasi-energy
     sum vanishes exactly (the singlet conjugate pairs are what keep the
     pairing, the determinant and the tree-level vacuum energy exact on
     the whole cell), and the spinor block is byte-identical to the
     published 32-level internal unitary (the singlets are additive:
     nothing in the geometry moves).
  T2 the rung is a phase: one automaton step = six stages, each applying
     the internal layer; the six-stage product multiplies the A pair by
     exp(+/- i*loop^3) and the B pair by exp(+/- i*loop^10) exactly - the
     ladder read as on-site phases - and the one-step Bloch matrix of the
     36-level cell factorises exactly as kron(P(q), U_cell^6): the
     singlets ride the same bond bands as the spinor module, and the
     16-sector spinor anchor (M-point gap 2.0029e-4) is unchanged.
  T3 the content of the cell: Tr Y = 0 and Tr Y^3 = 0 over the full
     36-level hypercharge multiset (the singlet levels contribute zero),
     and the 5/3 normalisation regression through the unchanged algebra.
  T4 the engine's own invariant list (cl10/supercell.verify) passes with
     the cell invariants included.

Verdict: the 36-level cell is a consistent, exactly symmetric object of
the engine, and the ladder rungs of A and B are literally the phases the
step applies to them.
"""

import json
import os
import time

import numpy as np

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from cl10 import cell as CE                                       # noqa: E402
from cl10 import supercell as SC                                  # noqa: E402
from cl10 import u10pred as U10                                   # noqa: E402
from cl10 import coupling as CPL                                  # noqa: E402
from cl10 import gauge as GA                                      # noqa: E402
from cl10 import config as CFG                                    # noqa: E402


def bond_product(q):
    """the 8x8 product of the six bond stages at block momentum q, in the
    order the step applies them (leftmost last)."""
    P = np.eye(8, dtype=complex)
    for even, j in SC.STAGES:
        P = SC.stage_matrix(even, j, q[j], SC.G10) @ P
    return P


def main():
    t0 = time.time()
    checks = []

    # ---- T0 the ladder identity between cell levels ---------------------
    m = CE.masses()
    nl = U10.neutrino_ladder()
    rel_nu = abs(m["m_nu_eV_identity"] / nl["m_nu_eV_ladder"] - 1.0)
    rung_ok = (CFG.CFG.cell.rung_A == 3 and CFG.CFG.cell.rung_B == 10)
    t0_ok = (m["m_A_GeV"] == nl["M_N_GeV"]
             and m["m_B_GeV"] == U10.dm_ladder_mass()
             and rel_nu < 1e-12 and rung_ok
             and CE.RUNG_A == 3 and CE.RUNG_B == 10)
    checks.append(("T0 the ladder identity between cell levels: m_A == the "
                   "seesaw Majorana rung M_N = %.4e GeV (exact), m_B == the "
                   "DM rung M_P*(1/16pi)^10 = %.4f GeV (exact), and "
                   "m_B^2/m_A = %.6f eV == the ladder neutrino mass (rel "
                   "%.1e): the exponent identity 2*10-3 = 17 holds between "
                   "two levels of the cell; rungs from the config block "
                   "cell: %s" % (m["m_A_GeV"], m["m_B_GeV"],
                                 m["m_nu_eV_identity"], rel_nu, rung_ok),
                   t0_ok))

    # ---- T1 the exact B-parity and the full-cell theorems ----------------
    Uc = CE.cell_step_np(CPL.G10)
    G = CE.parity_np()
    unit = float(np.abs(Uc @ Uc.conj().T - np.eye(CE.DIM)).max())
    g2 = float(np.abs(G @ G - np.eye(CE.DIM)).max())
    cU = float(np.abs(G @ Uc - Uc @ G).max())
    c45 = float(max(np.abs(G @ X - X @ G).max() for X in CE.gauge_ext_np()))
    from cl10 import engine as E
    U32 = E.asnumpy(SC.internal(SC.G10)).astype(complex)
    spinor_ident = float(np.abs(Uc[:32, :32] - U32).max())
    # the three structural theorems on the FULL cell (conjugate pairs):
    phc = np.sort(np.angle(np.linalg.eigvals(Uc)))
    pairing = float(np.abs(phc + phc[::-1]).max())
    detph = float(abs(np.angle(np.linalg.det(Uc))))
    vsum = float(abs(phc.sum()))
    t1_ok = (unit < 1e-12 and g2 < 1e-12 and cU < 1e-12 and c45 < 1e-12
             and spinor_ident == 0.0
             and pairing < 1e-12 and detph < 1e-12 and vsum < 1e-12)
    checks.append(("T1 the exact B-parity: the 36-level step is unitary to "
                   "%.1e; G = diag(I32,+1,+1,-1,-1) squares to I to %.1e, "
                   "commutes with the step to %.1e and with all 45 plane "
                   "generators extended by zero on the singlets to %.1e; "
                   "the full-cell theorems hold on all 36 levels: spectrum "
                   "symmetric under phase -> -phase (%.1e), det U = 1 "
                   "(phase %.1e), quasi-energy sum zero (%.1e) - the tree-"
                   "level vacuum energy vanishes exactly; the spinor block "
                   "is byte-identical to the "
                   "published internal unitary (diff %.1e): B cannot decay, "
                   "and the singlets move nothing"
                   % (unit, g2, cU, c45, pairing, detph, vsum,
                      spinor_ident), t1_ok))

    # ---- T2 the rung is a phase of the step ------------------------------
    Ust = CE.cell_step_np(CPL.G10)               # what one stage applies
    U6 = np.linalg.matrix_power(Ust, 6)          # the six-stage product
    phA = float(np.angle(U6[32, 32]))
    phAc = float(np.angle(U6[33, 33]))
    phB = float(np.angle(U6[34, 34]))
    phBc = float(np.angle(U6[35, 35]))
    dA = abs(phA - 6 * CE.phase_A())
    dAc = abs(phAc + 6 * CE.phase_A())
    dB = abs(phB - 6 * CE.phase_B())
    dBc = abs(phBc + 6 * CE.phase_B())
    dspin = float(np.abs(U6[:32, :32] - np.linalg.matrix_power(U32, 6)).max())
    qerr = {}
    for name, q in [("M", SC.M_POINTS[0]), ("GAMMA", SC.GAMMA),
                    ("R", SC.R_POINTS[0])]:
        Mc = SC.bloch(q, Ust)
        qerr[name] = float(np.abs(Mc - np.kron(bond_product(q), U6)).max())
    # the two singlet bands at M: the bond-band phases shifted by the rung
    pp = np.angle(np.linalg.eigvals(bond_product(SC.M_POINTS[0])))
    def band_shift(theta):
        return float(np.max(np.abs(np.angle(np.exp(1j * (pp + theta))) - pp)))
    shA = band_shift(6 * CE.phase_A())
    shB = band_shift(6 * CE.phase_B())
    gap_anchor = SC.gap(SC.M_POINTS[0], SC.sectors()[2])   # the 16 sector
    anchor_ok = abs(gap_anchor / 2.0029e-4 - 1.0) < 1e-5
    t2_ok = (dA < 1e-15 and dAc < 1e-15 and dB < 1e-15 and dBc < 1e-15
             and dspin < 1e-12
             and max(qerr.values()) < 1e-10
             and abs(shA - 6 * CE.phase_A()) < 1e-15
             and abs(shB - 6 * CE.phase_B()) < 1e-15 and anchor_ok)
    checks.append(("T2 the rung is a phase: the engine applies the internal "
                   "layer once per stage (six stages per step), so each "
                   "stage multiplies the A pair by exp(+/-i*%.6e) (loop^3 "
                   "= m_A/M_P = "
                   "%.6e) and the B pair by exp(+/-i*%.6e) (loop^10 = "
                   "m_B/M_P = %.6e) - "
                   "the mass-phase relation at the elementary time, each "
                   "level with its conjugate - and "
                   "the six-stage product carries phase %.6e for A+ (%.6e, "
                   "diff %.1e), %.6e for A- (diff %.1e), %.6e for B+ "
                   "(%.6e, "
                   "diff %.1e) and %.6e for B- (diff %.1e); the spinor block of the "
                   "six-stage product is the published six-stage spinor "
                   "(diff %.1e); the one-step Bloch matrix of the cell "
                   "factorises exactly as kron(P(q), U^6) at M/GAMMA/R "
                   "with errors %.1e/%.1e/%.1e - the singlets ride the "
                   "bond bands, shifted by their step phase (max shift "
                   "A %.6e, B %.6e at M) - and the spinor 16-sector "
                   "anchor is unchanged: M-point gap %.6e"
                   % (CE.phase_A(), CE.phase_A(), CE.phase_B(),
                      CE.phase_B(), phA, 6 * CE.phase_A(), dA, phAc, dAc,
                      phB, 6 * CE.phase_B(), dB, phBc, dBc, dspin,
                      qerr["M"], qerr["GAMMA"],
                      qerr["R"], shA, shB, gap_anchor), t2_ok))

    # ---- T3 the content of the cell --------------------------------------
    ys = []
    for y, n in CE.content():
        ys += [float(y)] * int(n)
    trY = sum(ys)
    trY3 = sum(y ** 3 for y in ys)
    nlev = len(ys)
    k53 = float(GA.five_thirds())
    t3_ok = (abs(trY) < 1e-12 and abs(trY3) < 1e-12 and nlev == 36
             and abs(k53 - 5.0 / 3.0) < 1e-12)
    checks.append(("T3 the content: %d hypercharge levels (the six fragment "
                   "classes, their conjugates, and the four singlet levels "
                   "A+, A-, B+, B- at Y = 0); Tr Y = "
                   "%.1e and Tr Y^3 = %.1e (both anomaly cancellations "
                   "unchanged by the singlets); the 5/3 normalisation of "
                   "the unchanged algebra: k = %.15f"
                   % (nlev, trY, trY3, k53), t3_ok))

    # ---- T4 the engine invariant list includes the cell ------------------
    vr = SC.verify()
    t4_ok = (bool(vr["ok"]) and vr["cell_unitary_err"] < 1e-12
             and vr["cell_parity_step_comm"] < 1e-12
             and vr["cell_parity_gauge_max"] < 1e-12)
    checks.append(("T4 cl10/supercell.verify() passes with the cell "
                   "invariants: ok = %s, unitarity %.1e, [G,U] %.1e, "
                   "[G,45] %.1e" % (vr["ok"], vr["cell_unitary_err"],
                                    vr["cell_parity_step_comm"],
                                    vr["cell_parity_gauge_max"]), t4_ok))

    ok = all(p for _, p in checks)
    branch = ("THE 36-LEVEL CELL IS AN ENGINE FACT: spinor module + the "
              "conjugate singlet pairs at rungs 3 and 10, exact B-parity, "
              "rungs as step phases, pairing/det/vacuum exact on the cell"
              if ok else "FAIL")

    out = {"experiment": "V196_cell_definition", "data": {
        "rungs": {"A": CE.RUNG_A, "B": CE.RUNG_B},
        "m_A_GeV": m["m_A_GeV"], "m_B_GeV": m["m_B_GeV"],
        "m_nu_identity_eV": m["m_nu_eV_identity"],
        "phase_A": CE.phase_A(), "phase_B": CE.phase_B(),
        "six_stage_phase_A": phA, "six_stage_phase_B": phB,
        "commutators": {"unitarity": unit, "G_squared": g2,
                        "G_U": cU, "G_45_max": c45,
                        "spinor_block_diff": spinor_ident},
        "full_cell_theorems": {"pairing": pairing, "det_phase": detph,
                               "vacuum_sum": vsum},
        "bloch_factorisation_err": qerr,
        "singlet_band_shift_at_M": {"A": shA, "B": shB},
        "spinor_M_gap": gap_anchor,
        "content": {"levels": nlev, "trY": trY, "trY3": trY3,
                    "five_thirds": k53},
    }}
    out["pass"] = bool(ok)
    out["branch"] = branch
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = ("the 36-level cell of cl10/cell.py verified at the "
                     "root: ladder identity m_nu = m_B^2/m_A between cell "
                     "levels (exact), B-parity exact against the step and "
                     "the 45 generators, the full-cell theorems exact - "
                     "spectrum +- symmetric, det U = 1, quasi-energy sum "
                     "zero - rung phases loop^3 and loop^10 "
                     "reproduced by the six-stage product on each "
                     "conjugate pair, exact "
                     "kron(P,U^6) transport, anomaly bookkeeping intact, "
                     "engine invariants green")
    out["tier"] = 2
    out["property"] = ("the automaton cell is the 32-D Cl(10,0) Dirac "
                       "module plus the two neutral singlet sectors - "
                       "conjugate pairs A (rung 3) and B (rung 10) - with "
                       "an exact B-parity, the full cell +- paired with "
                       "det = 1 and vanishing quasi-energy sum, and the "
                       "ladder "
                       "rung of a level is the phase the step applies to "
                       "it (Axiom 1)")
    out["programme"] = "V196"
    out["wall_time"] = round(time.time() - t0, 3)

    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v196_cell_definition.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V196 - the 36-level cell as an engine fact (Axiom 1 verified)")
    for n, p in checks:
        print(f"  [{'PASS' if p else 'FAIL'}] {n}")
    print(f"  {branch}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   wall {time.time()-t0:.2f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
