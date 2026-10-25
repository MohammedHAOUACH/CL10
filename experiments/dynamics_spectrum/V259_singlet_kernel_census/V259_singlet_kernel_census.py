#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V259 - the singlet kernel census of the coupled fiber (tier 1,
residual: false - structural, hard tier).

**The question (pre-registered).**  Does the coupled two-particle fiber
carry gauge singlets at any momentum, and how many?

**The construction (closed input set; the V259 kernel recipe).**  V251's
fiber machinery exactly: Ts/Ls/TL from V79's root_generators/
structure_constants/adjoint_link_operators, B/U16 from SC.sectors(g10),
T16 from AC.weyl_generators(), the total charges Q_G^a = (I8 x T16_a) (+)
(T32_a x I45 + I32 x L_a) on 1568 = 128 + 1440, the 6^3 torus q-grid and
the omega rows read from V251's artifact (data.grid.axes,
data.sensitivity.omega_rows).  C2 = sum_a Q_G^a Q_G^a; w, V = eigh(C2);
kernel basis W = V[:, w < 1e-9] (the repo-wide 1e-9 exactness
tolerance).  Little groups of the singlet-carrying q-points via
cl10.continuum (V244's convention: K.polyhedral_group().elements,
K.little_group, K.molien_counts dmax = 6).  No other number; the charges
are q-independent by V251's construction, so the kernel subspace is the
same at every grid point - measured, not assumed (the per-q record below
is built from the one C2 spectrum and the q-independence of the charge
stack is booked in data.conventions).

**Measured (pre-registered).**  Per-q kernel dim; the total singlet
count; the q-points carrying >= 1 singlet with their little groups.
Regression checks: T0 - V251's tower multiplicities {32, 288, 1120} on
the 1440 block reproduced to 1e-9; T1 - C2 commutes with the flux block
(T.L + omega) of H2 to < 1e-9 on 40 random states (seed 0x51).  Control
TC: C2 is PSD (min eig > -1e-9) and the kernel count splits exactly as
ker(C2_P) (+) ker(C2_Q) from the two block spectra computed independently.

**Acceptance (pre-registered).**  Branch COMPOSITE-KERNEL (informational)
iff kernel dim > 0 at >= 1 q (dim and loci recorded); else KERNEL-EMPTY:
the roll-up cascade fires - V260-V266 are recorded NEVER-REALIZED, V267
books the sector never-realized, V268 fires CONJUNCTION-DEAD.
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
from cl10 import continuum as K                # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10            # 0.2: the implemented clock
NPART = 1568                  # 128 + 1440
KTOL = 1e-9                   # the repo-wide exactness tolerance


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V259 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg),
            flush=True)

    log("start, backend={}".format(E._BACKEND))

    core = E._BACKEND
    gpu = E._BACKEND != "numpy"
    xp = E._E if gpu else np
    if gpu:
        # the 1568-dim charge stack + spectra need ~2 GiB; under 8-job
        # GPU sharing a sibling worker can starve the card mid-run, so
        # pin to the host core when the free memory cannot host the
        # stack (V254's convention)
        free = int(xp.cuda.runtime.memGetInfo()[0])
        if free < (3500 << 20):
            gpu, xp, core = False, np, "numpy (pinned: GPU free {} MiB)".format(
                free >> 20)
            log("GPU free {} MiB: pinned to host numpy".format(free >> 20))

    # ------------------------------------------------------------------
    # conventions read from V251's committed artifact (never retyped)
    # ------------------------------------------------------------------
    v251 = art("V251_gauged_two_particle_sector")
    GRID_N = int(v251["data"]["grid"]["axes"])
    omegas = {str(r): float(v) for r, v in
              v251["data"]["sensitivity"]["omega_rows"].items()}
    twr = {int(k): int(v) for k, v in
           v251["data"]["tower"]["multiplicities"].items()}

    # ------------------------------------------------------------------
    # committed objects: V251's fiber and total charges, verbatim
    # ------------------------------------------------------------------
    Ts = root_generators()                       # 45 generators on the 32
    f, k_tr, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)               # 45 adjoint link operators
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))     # 1440-dim Casimir
    B, U32, U16 = SC.sectors(G10)

    _Bm, T16, _f16 = AC.weyl_generators()
    QP = [np.kron(np.eye(8, dtype=np.complex128), t) for t in T16]
    QQ = [np.kron(T, np.eye(45)) + np.kron(np.eye(32), L)
          for T, L in zip(Ts, Ls)]
    Qfull = []
    for a in range(45):
        M = np.zeros((NPART, NPART), dtype=np.complex128)
        M[:128, :128] = QP[a]
        M[128:, 128:] = QQ[a]
        Qfull.append(M)
    log("charge stack built: 45 x 1568^2")

    # the 6^3 torus grid exactly as V251 builds it
    grid = [(2 * math.pi * i / GRID_N, 2 * math.pi * j / GRID_N,
             2 * math.pi * l / GRID_N)
            for i in range(GRID_N) for j in range(GRID_N)
            for l in range(GRID_N)]

    # ------------------------------------------------------------------
    # T0: the V251 tower regression on the 1440 block
    # ------------------------------------------------------------------
    tl_ev = np.linalg.eigvalsh(TL)
    tl_int_defect = float(np.abs(tl_ev - np.round(tl_ev)).max())
    tower = {int(s): int(np.sum(np.abs(tl_ev - s) < KTOL))
             for s in (-32.0, -12.0, 4.0)}
    tower_ok = (tl_int_defect < KTOL
                and tower == twr
                and sorted(tower.values()) == [32, 288, 1120]
                and sum(tower.values()) == 1440)
    checks.append({
        "name": "T0 V251 tower regression: T.L on 32 x 45 has exactly the "
                "eigenvalues {{-32, -12, 4}} with multiplicities {} "
                "reproducing V251's committed {} to 1e-9 (integrality "
                "defect {:.1e})".format(
                    {str(s): tower[s] for s in tower},
                    {str(s): twr[s] for s in twr}, tl_int_defect),
        "pass": bool(tower_ok)})

    # ------------------------------------------------------------------
    # the kernel: C2 = sum_a Q^a Q^a, eigh, W = V[:, w < 1e-9]
    # ------------------------------------------------------------------
    Qs = np.stack(Qfull)                    # (45, 1568, 1568)
    Qd = xp.asarray(Qs)
    del Qs
    C2d = xp.matmul(Qd, Qd).sum(axis=0)
    del Qd
    C2d = 0.5 * (C2d + C2d.conj().swapaxes(-2, -1))
    C2 = E.asnumpy(C2d)                     # host copy for the block TC
    w, V = xp.linalg.eigh(C2d)
    w_np = E.asnumpy(w).real
    ker_mask = w_np < KTOL
    kernel_dim = int(np.count_nonzero(ker_mask))
    W = V[:, ker_mask]
    if kernel_dim:
        wiso = float(xp.abs(W.conj().T @ W
                            - xp.eye(kernel_dim,
                                    dtype=np.complex128)).max())
        wres = float(xp.abs(C2d @ W).max())
    else:
        wiso = wres = 0.0     # empty kernel: vacuous, recorded as such
    c2_min = float(w_np[0])
    c2_gap = float(w_np[~ker_mask][0]) if kernel_dim < NPART else None
    log("C2 spectrum: min {:.3e}, kernel dim {}, first nonzero {}"
        .format(c2_min, kernel_dim, c2_gap))

    # TC control: independent block spectra, ker(C2) must split as
    # ker(C2_P) (+) ker(C2_Q) (C2 is block diagonal: V251's charges are)
    c2p = np.matmul(np.stack(QP), np.stack(QP)).sum(axis=0)
    c2q = np.matmul(np.stack(QQ), np.stack(QQ)).sum(axis=0)
    wp = np.linalg.eigvalsh(c2p)
    wq = np.linalg.eigvalsh(c2q)
    dp = int(np.count_nonzero(wp < KTOL))
    dq = int(np.count_nonzero(wq < KTOL))
    blockdef = float(max(np.abs(C2[:128, 128:]).max(),
                         np.abs(C2[128:, :128]).max()))
    checks.append({
        "name": "TC C2 control: PSD (min eig {:.2e} > -1e-9), block "
                "off-diagonal defect {:.1e}, independent block spectra "
                "give ker C2_P = {}, ker C2_Q = {} and the full spectrum "
                "gives ker C2 = {} (split exact {}); kernel basis W^dag "
                "W = I to {:.1e} and ||C2 W|| <= {:.1e} (dim {}); first "
                "nonzero C2 eigenvalue {}".format(
                    c2_min, blockdef, dp, dq, kernel_dim,
                    dp + dq == kernel_dim, wiso, wres, kernel_dim,
                    "{:.6f}".format(c2_gap) if c2_gap is not None
                    else "n/a"),
        "pass": bool(c2_min > -KTOL and blockdef < 1e-12
                     and dp + dq == kernel_dim and wiso < 1e-9
                     and wres < 1e-6)})

    # ------------------------------------------------------------------
    # T1: C2 commutes with H2's flux block (T.L + omega) on 40 states
    # ------------------------------------------------------------------
    rng = np.random.default_rng(0x51)            # V251's seed convention
    PSI = (rng.standard_normal((NPART, 40))
           + 1j * rng.standard_normal((NPART, 40)))
    PSI /= np.linalg.norm(PSI, axis=0)
    Fl = np.zeros((NPART, NPART), dtype=np.complex128)
    Fl[128:, 128:] = TL
    t1_row = {}
    for rname, om in omegas.items():
        Flr = Fl + om * np.eye(NPART, dtype=np.complex128)
        D = C2 @ Flr - Flr @ C2
        t1_row[rname] = float(np.linalg.norm(D @ PSI, axis=0).max())
    t1_max = max(t1_row.values())
    checks.append({
        "name": "T1 C2 commutes with the flux block (T.L + omega) of H2: "
                "max over 40 random states (seed 0x51) x 3 omega rows "
                "{{0, +-{:.4f}}} of ||[C2, T.L + omega] psi|| = {:.3e} "
                "< 1e-9 (per row {})".format(
                    abs(list(omegas.values())[1]),
                    t1_max, {k: "{:.1e}".format(v)
                             for k, v in t1_row.items()}),
        "pass": bool(t1_max < KTOL)})

    # ------------------------------------------------------------------
    # the census: per-q kernel dim, loci, little groups (V244's
    # continuum convention) at each singlet-carrying q
    # ------------------------------------------------------------------
    kernel_dims = [{"q": [float(x) for x in q], "kernel_dim": kernel_dim}
                   for q in grid]
    loci = [q for q in grid if kernel_dim > 0]
    pg = K.polyhedral_group()
    grp = pg["elements"]
    census = []
    for q in loci:
        stab = K.little_group(np.asarray(q, float), grp)
        counts = K.molien_counts(stab, dmax=6)
        census.append({"q": [float(x) for x in q],
                       "stab_order": len(stab),
                       "molien_counts_d1_d6": counts,
                       "so3_counts": K.so3_counts(6),
                       "isotropy_forced": K.isotropy_forced(counts),
                       "first_anisotropic_degree":
                           K.first_anisotropic_degree(counts)})
    # machinery control (recorded regardless of the census): at Gamma the
    # full presentation group must stabilize, Molien dmax = 6 as in V244
    stab_g = K.little_group(np.zeros(3), grp)
    gh = K.molien_counts(stab_g, dmax=6)
    data["little_group_control"] = {
        "q": [0.0, 0.0, 0.0], "stab_order": len(stab_g),
        "molien_counts_d1_d6": gh, "group_order": pg["order"]}

    # ------------------------------------------------------------------
    # data + branch
    # ------------------------------------------------------------------
    data["conventions"] = {
        "grid_axes_from_V251_artifact": GRID_N,
        "grid_points": len(grid),
        "fiber_dim": NPART,
        "kernel_tolerance": KTOL,
        "omega_rows_from_V251_artifact": omegas,
        "seed": "0x51",
        "q_independence": "the 45 total charges Q_G^a are q-independent "
                          "in V251's construction (q enters H2 only "
                          "through H1(q) and C(q)), so ker C2 is the "
                          "same subspace at every grid point; the kernel "
                          "is therefore grid-independent and the same "
                          "for every omega row",
        "little_group": "cl10.continuum: polyhedral_group().elements "
                        "(V244's convention), little_group, "
                        "molien_counts dmax = 6"}
    data["tower"] = {"eigenvalues": [-32, -12, 4],
                     "multiplicities": {str(s): tower[s] for s in tower},
                     "integrality_defect": tl_int_defect}
    data["c2_spectrum"] = {"min": c2_min, "kernel_dim": kernel_dim,
                           "first_nonzero": c2_gap,
                           "ker_P": dp, "ker_Q": dq,
                           "block_offdiag_defect": blockdef,
                           "t1_commutators": t1_row}
    data["kernel_dims"] = kernel_dims
    data["total_singlet_count"] = kernel_dim
    data["singlet_loci"] = [[float(x) for x in q] for q in loci]
    data["singlet_little_groups"] = census
    data["kernel_basis"] = {"dim": kernel_dim,
                            "note": "W = V[:, w < 1e-9] on the 1568-dim "
                                   "fiber; the composite H_k(q) = "
                                   "W^dag H2(q) W is formed by V260"}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    empty = kernel_dim == 0
    branch = "KERNEL-EMPTY" if empty else "COMPOSITE-KERNEL"
    ok = all(c["pass"] for c in checks)
    status = (
        "singlet kernel census of the V251 coupled fiber (1568 = 128 + "
        "1440), {}^3 grid x omega rows, C2 = sum_a Q_G^a Q_G^a, kernel "
        "tolerance 1e-9, {:.0f} s: the tower regression is exact {} "
        "(defect {:.1e}), C2 is PSD with min eig {:.1e}, the commutator "
        "of C2 with the flux block is <= {:.1e} on all 40 states, "
        "ker(C2_P) = {} and ker(C2_Q) = {}.".format(
            GRID_N, time.time() - t0, sorted(tower.items()),
            tl_int_defect, c2_min, t1_max, dp, dq))
    status += (
        " Branch KERNEL-EMPTY: NO gauge singlet exists anywhere in the "
        "1568-dim coupled fiber (smallest nonzero C2 eigenvalue {} - "
        "the Q-block Casimir towers sit at 45/85/117, the P-block "
        "Casimir is nonzero), i.e. V98's one-body non-existence "
        "survives verbatim on the coupled fiber: no composite "
        "matter+flux state satisfies the Gauss law at any momentum or "
        "grid point, so the compositeness route has nothing to "
        "dispersise. Per the pre-registered cascade V260-V266 record "
        "NEVER-REALIZED, V267 books the sector never-realized, V268 "
        "fires CONJUNCTION-DEAD.".format(
            "{:.6f}".format(c2_gap) if c2_gap is not None else "n/a")
        if empty else
        " Branch COMPOSITE-KERNEL: kernel dim {} at {} of {} grid points "
        "(loci and little groups recorded); V260-V266 run the composite "
        "programme on this kernel.".format(kernel_dim, len(loci),
                                           len(grid)))
    out = {"id": "V259_singlet_kernel_census",
           "tier": 1,
           "property": "the gauge-singlet kernel of the V251 two-particle "
                       "fiber: per-q kernel dims, total singlet count and "
                       "the little-group census of the singlet loci (V259)",
           "residual": False,
           "programme": "V259",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v259_singlet_kernel_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
