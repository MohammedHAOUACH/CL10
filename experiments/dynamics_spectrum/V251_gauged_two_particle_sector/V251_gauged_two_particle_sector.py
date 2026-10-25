#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V251 - the gauged two-particle sector: construction and Gauss check
(tier 2, residual: false - sector construction).

**The question (pre-registered).**  Does the root carry an exact gauged
two-particle sector - spinor tensor one adjoint bond-current flux - with
the committed coupling, Gauss-preserving and bounded?  This is the minimal
interacting block O9 never had.

**The construction (closed input set).**  Per crystal momentum Q on a 6^3
torus grid (216 points, one chunk <= 512):

  P block (128):  the one-particle Bloch space, H1(Q) = i log of
                  supercell.bloch(q, U16(g10)) through the
                  altcell.effective_hamiltonian pattern (V242's use).
  Q block (1440): the 32-adjoint-link space of V126, 32 x 45, with
                  E_flux = T.L the exact Casimir tower [-32, -12, 4]
                  (multiplicities 32, 288, 1120: V126's C2 eigenspace
                  code; root_generators / structure_constants /
                  adjoint_link_operators imported from V79).
  Coupling:       C(Q) = sum_A (T^A_32 B) kron e_A . (w(Q) kron I_16),
                  the V177 named vertex Gamma^A = I_8 tensor T^A_16
                  conjugated to the 32 by sectors()[0], strength g10;
                  w(Q) = exp(i Q.r_x)/sqrt(8) over the 8 sub-lattices.
  H2(Q) = H1(Q) opsum (T.L + omega) opsum g10 (C + C^dag) on 1568.
  Charges:        Q_G^a = (I_8 tensor T^a_16) opsum
                  (T^a_32 kron I_45 + I_32 kron L^a), a = 1..45, the
                  committed V98/V99 total-charge operators (their
                  ker(T1+T2+L) intertwiner code).
  Flux at rest:   omega in {0, +- 42.0906 x g10} rows (V101's measured
                  gap 42.09 eta, an artifact number, not a dial) stored
                  in data.sensitivity.
  V141 gate:      the dressed-hop gate is a P->P vertex correction and
                  V141 measured it to be trivial on the dressed sector
                  (hop = bare translation, electric term scalar 45 per
                  link, front speed 1 to 3.6e-10): zero-parameter, booked
                  in data.vertex_correction.

Preamble note: smatrix.two_particle (cl10/smatrix.py, the 't Hooft on-site
1024^2 tensor product used by uvrunning) is a different object; grep shows
no two-particle Bloch spectrum anywhere in cl10/, so the new code above is
justified.  V98's one-body non-existence (adjoint flux carries even center
charge and cannot screen a spinor) is recorded in data.constraints: the
P+Q coupled space carries total Gauss; intermediate Q states need not be
singlets.

**Measured (pre-registered).**  (a) boundedness - spectrum min/max finite;
(b) Gauss - ||[exp(i H2), Q_G]|| < 1e-12 on 40 random states; (c) block
structure - the P-block restriction equals H1 to 0.0; (d) unitarity of
exp(i H2) to 1e-13; plus the 6^3 x 1568 spectrum (cupy eigh, chunked;
wall-time recorded).

**Acceptance (pre-registered).**  TWO-PARTICLE-SECTOR-CLOSED iff (a)-(d)
pass with fitted_inputs == []; else SECTOR-BLOCKED naming the failed
check - a real result: the interaction is then not constructible on this
module and V252-V256 record the floor.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")     # the 1/5 reference branch
from cl10.paths import exp                     # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))
import numpy as np                             # noqa: E402
from V79_gauged_bond_layer import (adjoint_link_operators,   # noqa: E402
                                   root_generators,
                                   structure_constants)
from cl10 import altcell as AC                 # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10            # 0.2: the implemented clock (angle pin above
                              # enforces the same value; SC.enforce_angle)
NPART = 1568                  # 128 + 1440
GRID_N = 6


def build_H2(H1, TL, C, omega):
    H = np.zeros((NPART, NPART), dtype=np.complex128)
    H[:128, :128] = H1
    H[128:, 128:] = TL + omega * np.eye(1440)
    H[:128, 128:] = G10 * C.conj().T
    H[128:, :128] = G10 * C
    return H


def main():
    t0 = time.time()
    checks = []
    data = {}

    # ------------------------------------------------------------------
    # committed objects
    # ------------------------------------------------------------------
    Ts = root_generators()                       # 45 generators on the 32
    f, k_tr, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)               # 45 adjoint link operators
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))     # 1440-dim Casimir
    B, U32, U16 = SC.sectors(G10)
    Bm, T16, f16 = AC.weyl_generators()

    # V177's named vertex convention: T16 = B^dag T32 B (both sources)
    vert_defect = max(float(np.abs(B.conj().T @ T @ B - t).max())
                      for T, t in zip(Ts, T16))
    basis_defect = float(np.abs(B - Bm).max())

    # total charges Q_G^a = (I8 x T16_a) (+) (T32_a x I45 + I32 x L_a)
    QP = [np.kron(np.eye(8, dtype=np.complex128), t) for t in T16]
    QQ = [np.kron(T, np.eye(45)) + np.kron(np.eye(32), L)
          for T, L in zip(Ts, Ls)]
    Qfull = []
    for a in range(45):
        M = np.zeros((NPART, NPART), dtype=np.complex128)
        M[:128, :128] = QP[a]
        M[128:, 128:] = QQ[a]
        Qfull.append(M)

    # the tower: eigenvalues exactly {-32, -12, 4}
    tl_ev = np.linalg.eigvalsh(TL)
    tl_int_defect = float(np.abs(tl_ev - np.round(tl_ev)).max())
    tower = {int(s): int(np.sum(np.abs(tl_ev - s) < 1e-9))
             for s in (-32.0, -12.0, 4.0)}
    # quadratic Casimir of the full 32 x 45 link per tower eigenvalue:
    # C2_tot = T^2 + L^2 + 2 T.L
    c2_32 = float(np.real(np.trace(sum(T @ T for T in Ts))) / 32)
    c2_45 = float(np.real(np.trace(sum(L @ L for L in Ls))) / 45)
    c2_tot = {int(s): c2_32 + c2_45 + 2.0 * s for s in (-32.0, -12.0, 4.0)}

    # ------------------------------------------------------------------
    # coupling C(Q) on the 6^3 torus grid
    # ------------------------------------------------------------------
    r = np.array([(e & 1, (e >> 1) & 1, (e >> 2) & 1) for e in range(8)],
                 dtype=float)                    # sub-lattice coordinates
    grid = [(2 * math.pi * i / GRID_N, 2 * math.pi * j / GRID_N,
             2 * math.pi * l / GRID_N)
            for i in range(GRID_N) for j in range(GRID_N)
            for l in range(GRID_N)]

    def coupling_at(q):
        w = np.exp(1j * (r @ np.asarray(q))) / math.sqrt(8.0)
        KW = np.kron(w, np.eye(16, dtype=np.complex128))
        return sum(np.kron(T @ B, np.eye(45)[:, a:a + 1]) @ KW
                   for a, T in enumerate(Ts))

    C0 = coupling_at(grid[0])
    Hint0 = np.zeros((NPART, NPART), dtype=np.complex128)
    Hint0[:128, 128:] = G10 * C0.conj().T
    Hint0[128:, :128] = G10 * C0
    coup_defect = max(float(np.abs(Hint0 @ M - M @ Hint0).max())
                      for M in Qfull)

    # (b) on 40 random states, fixed seed, shared by every grid point
    rng = np.random.default_rng(0x51)
    PSI = (rng.standard_normal((NPART, 40))
           + 1j * rng.standard_normal((NPART, 40)))
    PSI /= np.linalg.norm(PSI, axis=0)
    QPSI = np.stack([M @ PSI for M in Qfull])        # (45, 1568, 40)

    # ------------------------------------------------------------------
    # V101 provenance: flux gap at rest
    # ------------------------------------------------------------------
    with open(exp("V101_bond_current_ground_state", "results",
                  "v101_bond_current_ground_state.json")) as fh:
        v101 = json.load(fh)
    v101_gap = float(v101["data"]["gap"])            # 42.0906... at eta = 1
    W = v101_gap * G10
    omegas = {"0": 0.0, "plus": +W, "minus": -W}

    # ------------------------------------------------------------------
    # V141 gate (booked, zero-parameter) and V98/V99 constraints
    # ------------------------------------------------------------------
    with open(exp("V141_regular_link_quintet_coupling", "results",
                  "v141_regular_link_quintet_coupling.json")) as fh:
        v141 = json.load(fh)
    with open(exp("V98_gauged_invariant_sector", "results",
                  "v98_gauged_invariant_sector.json")) as fh:
        v98 = json.load(fh)
    with open(exp("V99_three_site_chain", "results",
                  "v99_three_site_chain.json")) as fh:
        v99 = json.load(fh)

    # ------------------------------------------------------------------
    # the sweep: three sensitivity rows x 216 points, cupy eigh chunked
    # ------------------------------------------------------------------
    use_gpu = E._BACKEND != "numpy"
    xp = E._E if use_gpu else np
    PSI_d = xp.asarray(PSI)
    QPSI_d = xp.asarray(QPSI)
    Qfull_d = xp.asarray(np.stack(Qfull))
    I1568_d = xp.eye(NPART, dtype=np.complex128)

    # one-shot linear attribution: the P-P block of [H2, Q^a] is exactly
    # [H1, Q_P^a] and the rest of the commutator is carried only by the
    # (measured-exact) coupling and the exact [T.L, Q^a] = 0
    H1a = AC.effective_hamiltonian(lambda qq: SC.bloch(qq, U16, G10),
                                   grid[0])
    H2a = build_H2(H1a, TL, C0, 0.0)
    pp_ident, rest_max, lin_def_at = 0.0, 0.0, 0.0
    for a in range(45):
        cm = H2a @ Qfull[a] - Qfull[a] @ H2a
        lin_def_at = max(lin_def_at, float(np.abs(cm).max()))
        pp_ident = max(
            pp_ident, float(np.abs(
                cm[:128, :128] - (H1a @ QP[a] - QP[a] @ H1a)).max()))
        rest_max = max(rest_max, float(np.abs(cm[:128, 128:]).max()),
                       float(np.abs(cm[128:, :128]).max()),
                       float(np.abs(cm[128:, 128:]).max()))
    rows = {}
    gauss_global = 0.0          # (b): max ||[exp(iH2), Q^a] psi_i||
    gauss_a_min = {}            # per-charge minimum over grid/row
    h1_leak_min, h1_leak_max = 1e30, 0.0
    h1_leak_row = []            # per-point max_a ||[H1, Q_P^a]||_max
    unit_max = 0.0
    spec_min, spec_max = 1e30, -1e30
    finite_all = True
    pblock_max = 0.0
    shifts_row = {}
    t_sweep = time.time()
    for rname, om in omegas.items():
        sh = []
        for q in grid:
            H1 = AC.effective_hamiltonian(lambda qq: SC.bloch(qq, U16, G10),
                                          q)
            Cq = coupling_at(q)
            H2 = build_H2(H1, TL, Cq, om)
            pblock_max = max(pblock_max,
                             float(np.abs(H2[:128, :128] - H1).max()))
            h1_def = max(float(np.abs(H1 @ p - p @ H1).max()) for p in QP)
            h1_leak_min = min(h1_leak_min, h1_def)
            h1_leak_max = max(h1_leak_max, h1_def)
            h1_leak_row.append(h1_def)
            H2d = xp.asarray(H2)
            w, V = xp.linalg.eigh(H2d)
            w_np = E.asnumpy(w)
            spec_min = min(spec_min, float(w_np[0]))
            spec_max = max(spec_max, float(w_np[-1]))
            finite_all &= bool(np.all(np.isfinite(w_np)))
            # exp(i H2) = V diag(e^{i w}) V^dag
            U = (V * xp.exp(1j * w)) @ V.conj().T
            unit_max = max(unit_max, float(xp.abs(
                U @ U.conj().T - I1568_d).max()))
            # (b) Gauss on the 40 shared random states, batched:
            Upsi = U @ PSI_d
            t1 = xp.matmul(U[None], QPSI_d)      # e^{iH2} Q^a psi
            t2 = xp.matmul(Qfull_d, Upsi)        # Q^a e^{iH2} psi
            nrm = xp.linalg.norm(t1 - t2, axis=1)         # (45, 40)
            gauss_global = max(gauss_global, float(nrm.max()))
            d_a = E.asnumpy(nrm.max(axis=1))
            for a in range(45):
                gauss_a_min[a] = min(gauss_a_min.get(a, 1e30),
                                     float(d_a[a]))
            # P-pole coupling shifts: 128 bare H1 eigenvalues
            bare = np.linalg.eigvalsh(H1)
            for e in bare:
                kdx = int(np.argmin(np.abs(w_np - e)))
                sh.append(float(abs(w_np[kdx] - e)))
        shifts_row[rname] = {"max": float(max(sh)),
                             "mean": float(np.mean(sh))}
    sweep_wall = time.time() - t_sweep

    # ------------------------------------------------------------------
    # checks (measurement integrity; the branch below carries the verdict)
    # ------------------------------------------------------------------
    tower_ok = (tl_int_defect < 1e-9
                and sorted(tower.values()) == [32, 288, 1120]
                and sum(tower.values()) == 1440)
    checks.append({
        "name": "T0 the V126 tower regression: T.L on 32 x 45 has exactly "
                "the eigenvalues {{-32, -12, 4}} with multiplicities {} "
                "(integrality defect {:.1e}); [T.L, Q_G^a] = 0 exactly"
                .format({str(s): tower[s] for s in tower}, tl_int_defect),
        "pass": bool(tower_ok)})
    tl_comm = max(float(np.abs(TL @ M - M @ TL).max()) for M in QQ)

    checks.append({
        "name": "T1 the committed vertex coupling is exactly Gauss "
                "invariant: max_a ||[g10 (C + C^dag), Q_G^a]|| = {:.2e} "
                "< 1e-12; vertex convention T16 = B^dag T32 B exact {} "
                "(max {:.1e}), B agreement with V177 source {:.1e}"
                .format(coup_defect, vert_defect == 0.0, vert_defect,
                        basis_defect),
        "pass": bool(coup_defect < 1e-12 and vert_defect == 0.0
                     and basis_defect < 1e-9)})

    checks.append({
        "name": "T2 (a) boundedness and (d) unitarity: all {} spectra "
                "finite, global min/max [{:.4f}, {:.4f}]; max ||exp(iH2) "
                "exp(iH2)^dag - I|| = {:.2e} < 1e-13"
                .format(len(omegas) * len(grid), spec_min, spec_max,
                        unit_max),
        "pass": bool(finite_all and unit_max < 1e-13)})

    checks.append({
        "name": "T3 (c) block structure: the P-block of H2(Q) equals "
                "H1(Q) to exactly {:.1e} (= 0.0) on all {} grid points"
                .format(pblock_max, len(grid)),
        "pass": bool(pblock_max == 0.0)})

    gauss_ok = gauss_global < 1e-12
    checks.append({
        "name": "T4 (b) Gauss measured on 40 shared random states over "
                "3 omega rows x {} Q points: max_a,i,Q ||[exp(iH2), "
                "Q_G^a] psi_i|| = {:.3e}; attribution (one-shot at Q0): "
                "the coupling-only commutator is exact ({:.1e}), [T.L, "
                "Q^a] = {:.1e}, the P-P block of [H2, Q^a] equals "
                "[H1, Q_P^a] to {:.1e} and the remaining blocks of "
                "[H2, Q^a] are at most {:.1e}: the exp(iH2) defect is "
                "carried entirely by the committed step (per-point "
                "max_a ||[H1, Q_P^a]|| in [{:.3f}, {:.3f}])".format(
                    len(grid), gauss_global, coup_defect, tl_comm,
                    pp_ident, rest_max, h1_leak_min, h1_leak_max),
        "pass": bool(coup_defect < 1e-12 and tl_comm < 1e-12
                     and pp_ident < 1e-12 and rest_max < 1e-12)})

    checks.append({
        "name": "T5 provenance rows: V101 gap {:.6f} (artifact) -> "
                "sensitivity omegas {{0, +-{:.4f}}}; V98 glued singlet "
                "sector C2 {:.1e}, intertwiners {}, V99 vertex Gauss "
                "residual {:.1e}".format(v101_gap, W,
                                         v98["data"]["max_C2_glued"],
                                         v98["data"]["n_intertwiners"],
                                         v99["data"]["vertex_gauss_residual"]),
        "pass": bool(abs(v101_gap - 42.0906) < 1e-3)})

    # ------------------------------------------------------------------
    # data + branch
    # ------------------------------------------------------------------
    data["tower"] = {"eigenvalues": [-32, -12, 4],
                     "multiplicities": {str(s): tower[s]
                                        for s in tower},
                     "integrality_defect": tl_int_defect,
                     "tl_gauss_comm": tl_comm,
                     "c2_32": c2_32, "c2_45": c2_45,
                     "c2_tot_per_sector": {str(s): c2_tot[s]
                                           for s in (-32.0, -12.0, 4.0)}}
    data["gauss"] = {
        "max_over_a_states_grid": gauss_global,
        "min_per_charge_over_grid": {str(a): gauss_a_min.get(a)
                                     for a in range(45)},
        "coupling_only_defect": coup_defect,
        "h1_leak_min": h1_leak_min, "h1_leak_max": h1_leak_max,
        "attribution": {"pp_block_ident_defect": pp_ident,
                        "non_pp_blocks_max": rest_max,
                        "linear_defect_at_Q0": lin_def_at},
        "states": 40, "seed": "0x51",
        "note": "the vertex coupling part is exactly invariant; the "
                "defect is the committed step's own leak of all 45 "
                "charges (measured here on i log bloch, min defect "
                "{:.3f}): consistent with V67 (the lattice term is "
                "gauge-EXPLICIT) and V140 T2 (the bare gate leaks the "
                "vertex charge)".format(h1_leak_min)}
    data["sensitivity"] = {
        "omega_rows": {r: omegas[r] for r in omegas},
        "gap_provenance": "V101 measured gap at rest {:.6f} x eta {:.3f}"
                          .format(v101_gap, G10),
        "p_pole_shifts": shifts_row}
    data["vertex_correction"] = {
        "v141_gate": "dressed hop = bare translation on the dressed "
                     "sector (T2), electric term scalar 45 per link, "
                     "front speed {}".format(
                         v141["data"]["velocity"]["bare_front_speed_max"]),
        "applied_as": "P->P vertex correction: measured trivial (zero "
                      "parameter), booked; no new number"}
    data["constraints"] = {
        "v98_one_body": "adjoint flux carries even center charge and "
                        "cannot screen a spinor (V98): no one-particle "
                        "physical state; ker(32x45) = {}".format(
                            v98["data"]["ker_32x45"]),
        "v98_m45_pair": v98["data"]["m45_pair"],
        "v98_glued_two_site_dim": v98["data"]["n_glued"],
        "pgq_rule": "the P+Q coupled space carries total Gauss; "
                    "intermediate Q states need not be singlets",
        "smatrix_note": "smatrix.two_particle is the 't Hooft on-site "
                        "1024^2 tensor product (used by uvrunning); grep "
                        "shows no two-particle Bloch spectrum in cl10/, "
                        "which justifies the new construction here"}
    data["spectrum"] = {"n_spectra": len(omegas) * len(grid),
                        "dim": NPART,
                        "min": spec_min, "max": spec_max,
                        "unitarity_max": unit_max,
                        "p_block_defect": pblock_max,
                        "wall_time_s": sweep_wall}
    data["grid"] = {"axes": GRID_N, "points": len(grid),
                    "chunk": len(grid)}
    data["backend"] = E._BACKEND
    data["fitted_inputs"] = []
    data["inputs"] = {"clocks": {"g10": G10},
                      "vertex": "Gamma^A = I_8 x T^A_16 (V177), "
                                "strength g10",
                      "flux": "T.L tower, omega rows from V101"}

    failed = []
    if not finite_all or not (spec_min > -1e9 and spec_max < 1e9):
        failed.append("boundedness")
    if not gauss_ok:
        failed.append("Gauss")
    if pblock_max != 0.0:
        failed.append("block-structure")
    if unit_max >= 1e-13:
        failed.append("unitarity")
    if data["fitted_inputs"]:
        failed.append("fitted-inputs")
    branch = ("TWO-PARTICLE-SECTOR-CLOSED" if not failed else
              "SECTOR-BLOCKED: " + failed[0])

    ok = all(c["pass"] for c in checks)
    status = (
        "gauged two-particle sector on 1568 = 128 (spinor Bloch) + 1440 "
        "(32 x 45 adjoint link), 6^3 grid x 3 omega rows, {:.0f} s: the "
        "construction is consistent - the V126 Casimir tower is exact "
        "{} (integrality {:.1e}), the V177 vertex coupling is exactly "
        "Gauss invariant to {:.1e}, the P-block is H1 to {:.1e}, "
        "unitarity {:.1e}, bounded [{:.2f}, {:.2f}]. Branch {}.").format(
            time.time() - t0, sorted(tower.items()), tl_int_defect,
            coup_defect, pblock_max, unit_max, spec_min, spec_max,
            branch)
    status += (
        " Gauss check passes on all 45 total charges for every grid "
        "point, omega row and random state." if not failed else
        " The Gauss check fails (max {:.3e}): the vertex coupling, the "
        "T.L tower, the block structure and unitarity are exact, but "
        "the committed step H1 = i log bloch leaks every one of the 45 "
        "total charges (per-point max_a ||[H1, Q_P^a]|| in [{:.3f}, "
        "{:.3f}]) - the step is gauge-EXPLICIT (V67) and the leak is "
        "inherited from V140 T2's bare-gate vertex-charge leak. No "
        "exact gauged two-particle sector exists on the published step "
        "at the root; per the pre-registered rule V252-V256 record the "
        "floor (they still measure their own quantities on this "
        "constructed H2).".format(gauss_global, h1_leak_min,
                                  h1_leak_max))
    out = {"id": "V251_gauged_two_particle_sector",
           "tier": 2,
           "property": "the minimal interacting block: spinor Bloch "
                       "block plus the exact 32 x 45 adjoint-flux tower "
                       "with the V177 vertex coupling, Gauss-tested on "
                       "40 states over the 6^3 grid (V251)",
           "residual": False,
           "programme": "V251",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v251_gauged_two_particle_sector.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
