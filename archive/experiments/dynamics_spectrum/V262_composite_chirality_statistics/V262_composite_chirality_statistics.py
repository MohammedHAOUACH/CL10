#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V262 - the composite chirality/statistics census (tier 3,
residual: false - structural, hard tier).

**The question (pre-registered).**  What are the little-group and
statistics content of the composite gapless loci - could they carry the
electron?

**The construction (closed input set).**  At each V260 gapless locus:
V244's little-group/Molien census (cl10.continuum:
polyhedral_group().elements, little_group, molien_counts dmax = 6,
so3_counts, isotropy_forced, first_anisotropic_degree - V244's
l. 145-146/219-231 convention).  Statistics census (measured, not
assumed): the kernel basis W of the V259 recipe is decomposed under the
internal Spin(10) Cartan and the parity of the matter weight is recorded
per kernel state.  Cartan machinery note (the plan's convention rule):
the plan named V99/V100, whose committed code exposes only the C2
eigenspace channels (V99 l. 95-104); the Cartan-eigenspace machinery
actually committed in the repo is V65/V67's
(experiments/spacetime_matter/V67_breaking_chain.py l. 63-70):
h = clifford.generators() (10 gammas on the 16), Cartan planes
cartan = [1j*(h[a] @ h[a+1]) for a in (0, 2, 4, 6, 8)] conjugated by
B = SC.sectors()[0] (Dk = B^dag C B).  That source is used here; the
flux-block Cartan is built from the SAME 45 total charges QQ^a via the
V79 structure constants (generic-element centralizer, rank 5).  A
kernel state is odd (half-integer/fermionic) iff it carries P-block
(one-matter-quantum, spinorial centre) support - ||w_P|| > 0 - and
integer (bosonic) iff it is flux-supported with vanishing Cartan
weights; V98's centre obstruction says exactly which can occur.

**Cascade (pre-registered roll-up + V262 entry).**  V259 KERNEL-EMPTY
and V260 NEVER-REALIZED: no gapless locus exists, so the branch is
NO-LOCUS (the parity table is empty by construction).  The full census
below is the drafted non-cascade path.

**Acceptance (pre-registered).**  CHIRAL-COMPOSITE iff a gapless locus
has half-integer matter-parity content (impossible by V98 - if it
fires, the theory is rebuilt); BOSONIC-COMPOSITE-ONLY iff every kernel
state is even-parity (the centre obstruction re-derived at composite
level); else the NO-LOCUS cascade.
"""
import json
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
from cl10 import clifford                      # noqa: E402
from cl10 import continuum as K                # noqa: E402
from cl10 import coupling                      # noqa: E402
from cl10 import engine as E                   # noqa: E402
from cl10 import supercell as SC               # noqa: E402

G10 = coupling.G10
NPART = 1568
KTOL = 1e-9


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V262 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg),
            flush=True)

    log("start, backend={}".format(E._BACKEND))
    core = E._BACKEND
    gpu = E._BACKEND != "numpy"
    xp = E._E if gpu else np
    if gpu:
        free = int(xp.cuda.runtime.memGetInfo()[0])
        if free < (3500 << 20):
            gpu, xp, core = False, np, "numpy (pinned: GPU free {} MiB)".format(
                free >> 20)
            log("GPU free {} MiB: pinned to host numpy".format(free >> 20))

    v259 = art("V259_singlet_kernel_census")
    v260 = art("V260_composite_gap_ladder")
    loci = v260["data"]["gapless_loci"]

    # ------------------------------------------------------------------
    # the V259 recipe (shared construction)
    # ------------------------------------------------------------------
    Ts = root_generators()
    f, k_tr, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    _TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
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
    Qd = xp.asarray(np.stack(Qfull))
    C2d = xp.matmul(Qd, Qd).sum(axis=0)
    del Qd
    C2d = 0.5 * (C2d + C2d.conj().swapaxes(-2, -1))
    w, V = xp.linalg.eigh(C2d)
    km = E.asnumpy(w).real < KTOL
    kernel_dim = int(np.count_nonzero(km))
    W = V[:, km]
    log("kernel recomputed here: dim {}".format(kernel_dim))

    # ------------------------------------------------------------------
    # Cartan machinery, both blocks (V65/V67 convention; see docstring)
    # ------------------------------------------------------------------
    h = [E.asnumpy(m) for m in clifford.generators()]   # 10 gammas, 32
    assert all(x.shape == (32, 32) for x in h)
    cartan16 = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Dk = [B.conj().T @ C @ B for C in cartan16]            # Hermitian, 16
    Pk = [np.kron(np.eye(8, dtype=np.complex128), D)
          for D in Dk]                                    # on the 128

    # flux-block Cartan: centralizer of a generic element of the algebra
    rng = np.random.default_rng(0x62)
    c0 = rng.standard_normal(45)
    # ad_{Z} in the 45 basis via the committed structure constants:
    # adz[b, c] = sum_a c0_a f[a, b, c]; the Cartan directions are its
    # null space (centralizer of a generic element)
    adz = np.einsum("a,abc->bc", c0, f)
    u_ad, s_ad, vt_ad = np.linalg.svd(adz)
    cart_vecs = vt_ad[np.abs(s_ad) < 1e-9]   # null space (centralizer)
    rank = len(cart_vecs)                     # = 5 for Spin(10)
    CartanQ = [sum(v[a] * QQ[a] for a in range(45))
               for v in cart_vecs[:5]]
    cartan_comm = max(float(np.abs(A @ Bb - Bb @ A).max())
                      for i, A in enumerate(CartanQ)
                      for Bb in CartanQ[i + 1:])
    cartan_herm = max(float(np.abs(A - A.conj().T).max())
                      for A in CartanQ)
    dk_comm = max(float(np.abs(X @ Y - Y @ X).max())
                  for i, X in enumerate(Dk) for Y in Dk[i + 1:])
    dk_spec = sorted({round(float(s), 6) for X in Dk
                      for s in np.linalg.eigvalsh(X)})
    log("Cartan: flux rank {}, commutators {:.1e}/{:.1e}, hermiticity "
        "{:.1e}, P-block spectrum {}".format(
            rank, cartan_comm, dk_comm, cartan_herm, dk_spec))

    # ------------------------------------------------------------------
    # census (drafted path) + parity table
    # ------------------------------------------------------------------
    cascade = (v259["branch"] == "KERNEL-EMPTY" and kernel_dim == 0)
    pg = K.polyhedral_group()
    grp = pg["elements"]
    loci_census = []
    parity_rows = []
    Wn = E.asnumpy(W)
    if not cascade:
        for row in loci:
            q = np.asarray(row["q"], float)
            stab = K.little_group(q, grp)
            counts = K.molien_counts(stab, dmax=6)
            loci_census.append({
                "q": [float(x) for x in q], "stab_order": len(stab),
                "molien_counts_d1_d6": counts,
                "so3_counts": K.so3_counts(6),
                "isotropy_forced": K.isotropy_forced(counts),
                "first_anisotropic_degree":
                    K.first_anisotropic_degree(counts)})
        CartFull = []
        for X in Pk:
            M = np.zeros((NPART, NPART), dtype=np.complex128)
            M[:128, :128] = X
            CartFull.append(M)
        for X in CartanQ:
            M = np.zeros((NPART, NPART), dtype=np.complex128)
            M[128:, 128:] = X
            CartFull.append(M)
        for i in range(Wn.shape[1]):
            v = Wn[:, i]
            wP_n = float(np.linalg.norm(v[:128]))
            wts = [float(np.real(np.vdot(v, C @ v))) for C in CartFull]
            parity_rows.append({
                "state": i,
                "p_block_norm": wP_n,
                "cartan_weights": wts,
                "matter_parity": "odd" if wP_n > 1e-9 else "even",
                "weights_all_vanish": all(abs(x) < KTOL
                                          for x in wts)})
        if any(r["matter_parity"] == "odd" and
               not r["weights_all_vanish"] for r in parity_rows):
            branch = "CHIRAL-COMPOSITE"
        elif all(r["matter_parity"] == "even" for r in parity_rows):
            branch = "BOSONIC-COMPOSITE-ONLY"
        else:
            branch = "CHIRAL-COMPOSITE"
    else:
        branch = "NO-LOCUS"

    stab_g = K.little_group(np.zeros(3), grp)
    gh = K.molien_counts(stab_g, dmax=6)

    # ------------------------------------------------------------------
    # checks + data
    # ------------------------------------------------------------------
    checks.append({
        "name": "T0 cascade bookkeeping: V259 branch {} with recorded "
                "total_singlet_count {}, V260 branch {} with {} gapless "
                "loci, kernel recomputed here dim {}".format(
                    v259["branch"], v259["data"]["total_singlet_count"],
                    v260["branch"], len(loci), kernel_dim),
        "pass": bool(v259["branch"] == "KERNEL-EMPTY"
                     and int(v259["data"]["total_singlet_count"]) == 0
                     and kernel_dim == 0 and not loci)})
    checks.append({
        "name": "T1 Cartan machinery (runs on the cascade): the five "
                "flux-block Cartan combinations from the V79 structure-"
                "constant centralizer have rank {}, pairwise commutators "
                "{:.1e} < 1e-9, hermiticity {:.1e}; the V67-convention "
                "P-block Cartan Dk commute to {:.1e} with spectra {} "
                "(plane generators: eigenvalues +-1)".format(
                    rank, cartan_comm, cartan_herm, dk_comm, dk_spec),
        "pass": bool(rank == 5 and cartan_comm < 1e-9
                     and cartan_herm < 1e-9 and dk_comm < 1e-9
                     and dk_spec == [-1.0, 1.0])})
    checks.append({
        "name": "T2 little-group machinery control at Gamma: stab order "
                "{} (group order {}), molien d1-d6 {} (icosahedral "
                "invariants [0, 1, 0, 1, 0, 2])".format(
                    len(stab_g), pg["order"], gh),
        "pass": bool(len(stab_g) == 60 and gh == [0, 1, 0, 1, 0, 2])})
    checks.append({
        "name": "T3 statistics census state: cascade {} -> parity table "
                "has {} rows, loci census {} entries, branch {} matches "
                "the guard".format(cascade, len(parity_rows),
                                   len(loci_census), branch),
        "pass": bool((cascade and branch == "NO-LOCUS"
                      and not parity_rows and not loci_census)
                     or (not cascade and parity_rows))})

    data["conventions"] = {
        "cartan_source": "P block: V65/V67 convention, "
                         "experiments/spacetime_matter/"
                         "V67_breaking_chain.py l. 63-70 (clifford."
                         "generators, planes (0,2,4,6,8), conjugated by "
                         "SC.sectors()[0]); Q block: centralizer of a "
                         "generic element in the committed V79 "
                         "structure constants (rank 5); seed 0x62",
        "cartan_rank": rank,
        "parity": "odd (half-integer/fermionic) iff the kernel state "
                  "carries P-block (single-matter-quantum, spinorial "
                  "centre) support ||w_P|| > 0; even otherwise",
        "little_group": "V244 convention: polyhedral_group().elements, "
                        "molien_counts dmax = 6"}
    data["cartan_controls"] = {"rank": rank,
                               "cartanQ_comm": cartan_comm,
                               "cartanQ_hermiticity": cartan_herm,
                               "Dk_comm": dk_comm,
                               "Dk_spectrum": dk_spec}
    data["loci_census"] = loci_census
    data["parity_rows"] = parity_rows
    data["little_group_control"] = {"stab_order": len(stab_g),
                                    "molien_counts_d1_d6": gh}
    data["cascade"] = {
        "rule": "V259 KERNEL-EMPTY / no gapless loci -> NO-LOCUS "
                "(pre-registered acceptance third branch + roll-up)",
        "fired": bool(cascade),
        "v259_branch": v259["branch"], "v260_branch": v260["branch"]}
    data["backend"] = core
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0

    ok = all(c["pass"] for c in checks)
    status = (
        "composite chirality/statistics census ({:.0f} s): Cartan "
        "machinery verified live (flux Cartan rank {}, commutators "
        "{:.1e}; V67-convention P-Cartan spectra {}), little-group "
        "control at Gamma exact, V259 kernel dim recomputed {}, V260 "
        "gapless loci {}. Branch {}.").format(
            time.time() - t0, rank, cartan_comm, dk_spec, kernel_dim,
            len(loci), branch)
    status += (
        " There are no composite gapless loci and no kernel states at "
        "all, so the parity table is empty: the question whether the "
        "composites carry the electron (half-integer matter parity) has "
        "no instances on this construction - V98's centre obstruction "
        "is not even reached, because the singlet sector itself does "
        "not exist. NO-LOCUS by the pre-registered rule; the drafted "
        "per-locus Molien + Cartan-parity census ran not a line."
        if cascade else
        " {} kernel states, parities {}, loci census {} rows.".format(
            len(parity_rows),
            [r["matter_parity"] for r in parity_rows],
            len(loci_census)))
    out = {"id": "V262_composite_chirality_statistics",
           "tier": 3,
           "property": "the little-group Molien census and the "
                       "matter-parity (statistics) census of the "
                       "composite singlet kernel (V262)",
           "residual": False,
           "programme": "V262",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v262_composite_chirality_statistics.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
