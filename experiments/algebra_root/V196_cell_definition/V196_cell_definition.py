#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""V196 - the cell definition as a root-level fact (the 34-level cell).

**What this experiment pins.**  cl10/cell.py defines the automaton's cell:
the 32-D Dirac module S of Cl(10,0) plus the two neutral singlet levels
A (rung 3, the seesaw partner) and B (rung 10, the dark matter), with the
exact B-parity.  This is the axiom-level content of Section 2 (Axiom 1), so
the check sits in the algebra-root programme: the cell is part of the root
definition, and every other DM experiment (V194 census, V195 production)
presupposes it.

Checks
  T0 the ladder identity between cell levels: m_A is the seesaw Majorana
     rung, m_B is the DM rung, and m_nu = m_B^2 / m_A reproduces the
     ladder neutrino mass exactly - the exponent identity 2*10-3 = 17
     becomes an identity between two levels of the cell.
  T1 the exact B-parity of the cell: the 34-level step is unitary,
     G = diag(I32,+1,-1) is an exact Z2 commuting with the step and with
     all 45 plane generators extended by zero, and the spinor block is
     byte-identical to the published 32-level internal unitary (the
     singlet levels are additive: nothing in the geometry moves).
  T2 the rung is a phase: one automaton step = six stages, each applying
     the internal layer; the six-stage product multiplies the singlet A
     by exp(i*loop^3) and B by exp(i*loop^10) exactly - the ladder read
     as on-site phases - and the one-step Bloch matrix of the 34-level
     cell factorises exactly as kron(P(q), U_cell^6): the singlets ride
     the same bond bands as the spinor module, and the 16-sector spinor
     anchor (M-point gap 2.0029e-4) is unchanged.
  T3 the content of the cell: Tr Y = 0 and Tr Y^3 = 0 over the full
     34-level hypercharge multiset (the singlets contribute zero), and
     the 5/3 normalisation regression through the unchanged algebra.
  T4 the engine's own invariant list (cl10/supercell.verify) passes with
     the cell invariants included.

Verdict: the 34-level cell is a consistent, exactly symmetric object of
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

    # ---- T1 the exact B-parity of the 34-level cell ---------------------
    U34 = CE.cell_step_np(CPL.G10)
    G = CE.parity_np()
    unit = float(np.abs(U34 @ U34.conj().T - np.eye(CE.DIM)).max())
    g2 = float(np.abs(G @ G - np.eye(CE.DIM)).max())
    cU = float(np.abs(G @ U34 - U34 @ G).max())
    c45 = float(max(np.abs(G @ X - X @ G).max() for X in CE.gauge_ext_np()))
    from cl10 import engine as E
    U32 = E.asnumpy(SC.internal(SC.G10)).astype(complex)
    spinor_ident = float(np.abs(U34[:32, :32] - U32).max())
    t1_ok = (unit < 1e-12 and g2 < 1e-12 and cU < 1e-12 and c45 < 1e-12
             and spinor_ident == 0.0)
    checks.append(("T1 the exact B-parity: the 34-level step is unitary to "
                   "%.1e; G = diag(I32,+1_A,-1_B) squares to I to %.1e, "
                   "commutes with the step to %.1e and with all 45 plane "
                   "generators extended by zero on the singlets to %.1e; "
                   "the spinor block is byte-identical to the published "
                   "internal unitary (diff %.1e): B cannot decay, and the "
                   "singlets move nothing"
                   % (unit, g2, cU, c45, spinor_ident), t1_ok))

    # ---- T2 the rung is a phase of the step ------------------------------
    Ust = CE.cell_step_np(CPL.G10)               # what one stage applies
    U6 = np.linalg.matrix_power(Ust, 6)          # the six-stage product
    phA = float(np.angle(U6[32, 32]))
    phB = float(np.angle(U6[33, 33]))
    dA = abs(phA - 6 * CE.phase_A())
    dB = abs(phB - 6 * CE.phase_B())
    dspin = float(np.abs(U6[:32, :32] - np.linalg.matrix_power(U32, 6)).max())
    qerr = {}
    for name, q in [("M", SC.M_POINTS[0]), ("GAMMA", SC.GAMMA),
                    ("R", SC.R_POINTS[0])]:
        M34 = SC.bloch(q, Ust)
        qerr[name] = float(np.abs(M34 - np.kron(bond_product(q), U6)).max())
    # the two singlet bands at M: the bond-band phases shifted by the rung
    pp = np.angle(np.linalg.eigvals(bond_product(SC.M_POINTS[0])))
    def band_shift(theta):
        return float(np.max(np.abs(np.angle(np.exp(1j * (pp + theta))) - pp)))
    shA = band_shift(6 * CE.phase_A())
    shB = band_shift(6 * CE.phase_B())
    gap_anchor = SC.gap(SC.M_POINTS[0], SC.sectors()[2])   # the 16 sector
    anchor_ok = abs(gap_anchor / 2.0029e-4 - 1.0) < 1e-5
    t2_ok = (dA < 1e-15 and dB < 1e-15 and dspin < 1e-12
             and max(qerr.values()) < 1e-10
             and abs(shA - 6 * CE.phase_A()) < 1e-15
             and abs(shB - 6 * CE.phase_B()) < 1e-15 and anchor_ok)
    checks.append(("T2 the rung is a phase: the engine applies the internal "
                   "layer once per stage (six stages per step), so each "
                   "stage multiplies A by exp(i*%.6e) (loop^3 = m_A/M_P = "
                   "%.6e) and B by exp(i*%.6e) (loop^10 = m_B/M_P = %.6e) - "
                   "the mass-phase relation at the elementary time - and "
                   "the six-stage product carries phase %.6e for A (6*loop^3 "
                   "= %.6e, diff %.1e) and %.6e for B (6*loop^10 = %.6e, "
                   "diff %.1e); the spinor block of the "
                   "six-stage product is the published six-stage spinor "
                   "(diff %.1e); the one-step Bloch matrix of the cell "
                   "factorises exactly as kron(P(q), U^6) at M/GAMMA/R "
                   "with errors %.1e/%.1e/%.1e - the singlets ride the "
                   "bond bands, shifted by their step phase (max shift "
                   "A %.6e, B %.6e at M) - and the spinor 16-sector "
                   "anchor is unchanged: M-point gap %.6e"
                   % (CE.phase_A(), CE.phase_A(), CE.phase_B(),
                      CE.phase_B(), phA, 6 * CE.phase_A(), dA, phB,
                      6 * CE.phase_B(), dB, dspin, qerr["M"], qerr["GAMMA"],
                      qerr["R"], shA, shB, gap_anchor), t2_ok))

    # ---- T3 the content of the cell --------------------------------------
    ys = []
    for y, n in CE.content():
        ys += [float(y)] * int(n)
    trY = sum(ys)
    trY3 = sum(y ** 3 for y in ys)
    nlev = len(ys)
    k53 = float(GA.five_thirds())
    t3_ok = (abs(trY) < 1e-12 and abs(trY3) < 1e-12 and nlev == 34
             and abs(k53 - 5.0 / 3.0) < 1e-12)
    checks.append(("T3 the content: %d hypercharge levels (the six fragment "
                   "classes, their conjugates, and A, B at Y = 0); Tr Y = "
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
    branch = ("THE 34-LEVEL CELL IS AN ENGINE FACT: spinor module + rung-3 "
              "and rung-10 singlets, exact B-parity, rungs as step phases"
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
        "bloch_factorisation_err": qerr,
        "singlet_band_shift_at_M": {"A": shA, "B": shB},
        "spinor_M_gap": gap_anchor,
        "content": {"levels": nlev, "trY": trY, "trY3": trY3,
                    "five_thirds": k53},
    }}
    out["pass"] = bool(ok)
    out["branch"] = branch
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = ("the 34-level cell of cl10/cell.py verified at the "
                     "root: ladder identity m_nu = m_B^2/m_A between cell "
                     "levels (exact), B-parity exact against the step and "
                     "the 45 generators, rung phases loop^3 and loop^10 "
                     "reproduced by the six-stage product, exact "
                     "kron(P,U^6) transport, anomaly bookkeeping intact, "
                     "engine invariants green")
    out["tier"] = 2
    out["property"] = ("the automaton cell is the 32-D Cl(10,0) Dirac "
                       "module plus the singlet levels A (rung 3) and B "
                       "(rung 10) with an exact B-parity, and the ladder "
                       "rung of a level is the phase the step applies to "
                       "it (Axiom 1)")
    out["programme"] = "V196"
    out["wall_time"] = round(time.time() - t0, 3)

    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v196_cell_definition.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V196 - the 34-level cell as an engine fact (Axiom 1 verified)")
    for n, p in checks:
        print(f"  [{'PASS' if p else 'FAIL'}] {n}")
    print(f"  {branch}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   wall {time.time()-t0:.2f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
