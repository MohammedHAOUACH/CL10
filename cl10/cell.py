"""The automaton cell: the Clifford-Dirac module plus the two singlet levels.

This module is the single definition of what lives on a cell (Section 2,
Axiom 1).  The cell's eight sub-sites each carry the module

    M = S  (+)  C|A>  (+)  C|B>        dim M = 34

where S is the 32-dimensional Dirac module of Cl(10,0) (Section 3, Appendix A)
and A, B are two neutral singlet levels: they transform trivially under the ten
Clifford generators - equivalently, every one of the 45 internal-plane
generators X_[a,b] acts on them as zero - and carry hypercharge 0.

A rung.  A cell level on ladder rung n has mass M_P*loop^n.  The elementary
time of the automaton is the STAGE (Axiom 1: one step = six stages = 6 t_P),
and the mass-phase relation says a level of mass m accumulates the phase
m*t: so the internal layer applies to a rung-n level the phase

    exp(i * loop^n) = exp(i * m / M_P)      per stage,

and the six-stage step applies it six times, exp(i * 6 loop^n), the phase
m * 6 t_P.  The rung assignment is the config block `cell` (rung_A = 3,
rung_B = 10): the ladder of Section 12.4, read as the on-site phases of the
two singlet levels - a level's mass is the phase it carries.

    A: rung 3, m_A = M_P*loop^3 = 2.4467e9 GeV   - the seesaw singlet
                                                        (Section 12.3)
    B: rung 10, m_B = M_P*loop^10 = 118.57 GeV   - the dark matter
                                                        (Section 12.2)

B-parity.  The internal Z_2 G = diag(I_32, +1, -1), odd on B alone, is exact:
the step is block-diagonal and G commutes with the 32-D spinor block (any of
its angles) and with all 45 plane generators extended by zero.  This is the
grading theorem of V194 read constructively: the exact Z_2's of the step that
commute with the 45 are products of the global sign, the mirror, and this
B-parity - and V194's search shows the spinor module itself contributes no
other fragment-selective one.  B being odd makes it the stable carrier; A being
even keeps the seesaw a boson-exchange process.

The spinor block of everything here (the 32-D U(g), the 45 generators, the
band structure) is exactly the engine of Sections 3-7; the two singlet levels
are additive and never touch the geometry (V194's "additive" mechanism, applied
as built).  Content and anomaly bookkeeping: content() gives the cell's
hypercharge multiset; the singlets contribute 0 to Tr Y and Tr Y^3.

Self-test: python -m cl10.cell
"""

import numpy as np

from . import clifford
from . import config as _config
from . import coupling as _coupling
from . import engine as _engine
from . import sm_qca as _sm
from . import u10pred as _u10

# the rung assignment (config block `cell`; Section 12.4 Table 4)
RUNG_A = int(_config.CFG.cell.rung_A)
RUNG_B = int(_config.CFG.cell.rung_B)
DIM = 34


def rung_phase(n):
    """the phase a rung-n level carries per stage: loop^n = (M_P*loop^n)/M_P"""
    return float(_u10.LOOP) ** int(n)


def phase_A():
    """B-parity-even singlet A: the phase loop^rung_A applied per stage"""
    return rung_phase(RUNG_A)


def phase_B():
    """B-parity-odd singlet B: the phase loop^rung_B applied per stage"""
    return rung_phase(RUNG_B)


def cell_step_np(g10=None):
    """the 34-D internal layer of one stage: blkdiag(U(g), e^{i*phi_A},
    e^{i*phi_B}), phi_X = loop^rung_X = m_X/M_P (the mass-phase relation at
    the elementary time t_P).  Host numpy; g10 defaults to the root angle
    (coupling.G10)."""
    g = _coupling.G10 if g10 is None else g10
    U32 = _engine.asnumpy(_sm.gut_qca_step(_sm.plane_generators(), g)).astype(complex)
    U = np.zeros((DIM, DIM), dtype=np.complex128)
    U[:32, :32] = U32
    U[32, 32] = np.exp(1j * phase_A())
    U[33, 33] = np.exp(1j * phase_B())
    return U


def parity_np():
    """the B-parity G = diag(I_32, +1, -1): exact by construction."""
    return np.diag([1.0] * 33 + [-1.0]).astype(np.complex128)


def gauge_ext_np():
    """the 45 internal-plane generators extended to the cell: they act on S as
    h_a h_b and as 0 on A, B - the singlets are singlets."""
    h = [_engine.asnumpy(g).astype(complex) for g in clifford.generators()]
    out = []
    for a in range(10):
        for b in range(a):
            X = np.zeros((DIM, DIM), dtype=np.complex128)
            X[:32, :32] = h[a] @ h[b]
            out.append(X)
    return out


def masses():
    """the two singlet masses (ladder facts, Section 12.4) and the seesaw
    identity m_nu = m_B^2/m_A they satisfy."""
    m_A = _u10.ladder_scale(RUNG_A)
    m_B = _u10.ladder_scale(RUNG_B)
    return {"m_A_GeV": m_A, "m_B_GeV": m_B,
            "m_nu_eV_identity": (m_B * m_B / m_A) * 1e9}


def content():
    """the cell's hypercharge content as (Y, multiplicity) classes: the six
    fragment classes of the 16 and their conjugates (the V65 Georgi-Glashow
    fragment table), plus the two singlets.  Sums: Tr Y = 0, Tr Y^3 = 0."""
    frags = [(1 / 6, 6), (-2 / 3, 3), (1 / 3, 3),
             (-1 / 2, 2), (1.0, 1), (0.0, 1)]
    out = []
    for y, n in frags:
        out.append((y, n))
        out.append((-y, n))
    out.append((0.0, 2))            # the singlet levels A and B
    return out


def verify():
    """the cell facts: unitarity, exact B-parity, the ladder identities between
    cell levels, and the anomaly bookkeeping of the content."""
    r = {}
    U = cell_step_np()
    G = parity_np()
    r["rung_A"] = RUNG_A
    r["rung_B"] = RUNG_B
    r["phase_A"] = phase_A()
    r["phase_B"] = phase_B()
    r["m_A_GeV"] = masses()["m_A_GeV"]
    r["m_B_GeV"] = masses()["m_B_GeV"]
    r["unitarity"] = float(np.abs(U @ U.conj().T - np.eye(DIM)).max())
    r["G_squared"] = float(np.abs(G @ G - np.eye(DIM)).max())
    r["parity_step_comm"] = float(np.abs(G @ U - U @ G).max())
    r["parity_gauge_max"] = float(max(np.abs(G @ X - X @ G).max()
                                      for X in gauge_ext_np()))
    # the ladder identities: the singlet masses ARE the neutrino-module rungs
    nl = _u10.neutrino_ladder()
    r["m_A_eq_M_N"] = (r["m_A_GeV"] == nl["M_N_GeV"])
    r["m_B_eq_dm_ladder"] = (r["m_B_GeV"] == _u10.dm_ladder_mass())
    r["m_nu_identity_rel"] = float(
        abs(masses()["m_nu_eV_identity"] / nl["m_nu_eV_ladder"] - 1.0))
    # content anomaly bookkeeping
    ys = [(y, n) for y, n in content()]
    r["trY"] = float(sum(y * n for y, n in ys))
    r["trY3"] = float(sum(y ** 3 * n for y, n in ys))
    r["ok"] = (r["unitarity"] < 1e-12 and r["G_squared"] < 1e-12
               and r["parity_step_comm"] < 1e-12 and r["parity_gauge_max"] < 1e-12
               and r["m_A_eq_M_N"] and r["m_B_eq_dm_ladder"]
               and r["m_nu_identity_rel"] < 1e-12
               and abs(r["trY"]) < 1e-12 and abs(r["trY3"]) < 1e-12)
    return r


def main():
    r = verify()
    print("CL10 cell (cl10/cell.py) - the 34-level cell module")
    print(f"  rungs: A at n = {r['rung_A']}, B at n = {r['rung_B']}")
    print(f"  m_A = {r['m_A_GeV']:.6e} GeV, m_B = {r['m_B_GeV']:.6e} GeV")
    print(f"  stage phases: A {r['phase_A']:.6e}, B {r['phase_B']:.6e}")
    print(f"  unitarity {r['unitarity']:.2e}, G^2-I {r['G_squared']:.2e}, "
          f"[G,U] {r['parity_step_comm']:.2e}, [G,45] max {r['parity_gauge_max']:.2e}")
    print(f"  ladder identities: m_A == M_N {r['m_A_eq_M_N']}, "
          f"m_B == dm rung {r['m_B_eq_dm_ladder']}, "
          f"m_nu identity rel {r['m_nu_identity_rel']:.2e}")
    print(f"  content: Tr Y = {r['trY']:.2e}, Tr Y^3 = {r['trY3']:.2e}")
    print("VERDICT:", "PASS" if r["ok"] else "FAIL")
    return r


if __name__ == "__main__":
    main()
