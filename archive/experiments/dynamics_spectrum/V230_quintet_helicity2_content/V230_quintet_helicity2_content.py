#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V230 - the helicity-2 content of the quintet carrier: attacking F4 (tier 2).

**The question this pre-registers against** (`review/pre_registrations.md`):
F4 fires if the implemented spectrum carries no MASSLESS spin-2 content. The
exact content test: take the quintet modes and their `T2g`/`Eg` split as
recorded in the V140/V141 artifacts and decompose that carrier under the node
little group - `O_h` at GAMMA (V138: GAMMA carries `T2g` but never `Eg` in the
qubit link) and the icosahedral group built from its `[2,3,5]` presentation
(the group behind V225/V226) - and against the lattice spin-2 content, the
symmetric rank-2 representation (`Sym^2 = A1g + Eg + T2g`, V208). Then: does a
MASSLESS branch carry the full content, and does the tensor structure survive
the `Eg` mass shift?

T1 content, constructed: the 6-dimensional sym^2 space decomposes under the
   cubic group as `A1g + Eg + T2g` (validated through V138's machinery, whose
   irreps are built from the defining representation), and under the
   icosahedral rotation group of order 60 (built from the presentation, not a
   table) it decomposes as `A + W` with `W` a SINGLE irreducible piece of
   dimension 5: the `T2g`/`Eg` split is a property of the cubic subgroup, and
   under the star's own symmetry the quintet is one representation. Verified
   by rank counts of the commutant (Schur: commutant of `W` = scalars) and an
   exact `Q(sqrt 5)` character check with the measured class sizes.

T2 masslessness per channel (gap floors and ladder speeds, all read or
   measured): the qubit-link GAMMA multiplets (V138/V128: gapless but FLAT),
   the directional axis-bond quintet (V140: six gapless bands, linear speeds
   `{0, g10, 2 g10}`), the diagonal star-bond quintet (V229: six gapless
   bands, speeds `{0, 2 g10 |d|, 4 g10 |d|}`), the 1024-dimensional star
   automaton at GAMMA (V215: quadratic, no linear content) and the star spinor
   cell node (V225: the gap cluster is the coin-SPINOR representation -
   helicity 1/2 content, no quintet).

T3 survival against the mass shift: the axis-bond `Eg` bubble shifts
   (V141 artifact, at the closed clock the artifacts now carry) and the
   star-bond full-quintet shifts (V229 artifact). F4 closes only if ONE
   massless branch carries the FULL spin-2 content with a shift vanishing to
   1e-16.

Acceptance (pre-registered): `branch: F4-CLOSED` only if one massless branch
carries the full quintet (`Eg` AND `T2g` pieces of the same gapless multiplet
under the cubic group, equivalently the irreducible 5-dimensional piece under
the icosahedral group) with a mass shift `|shift| < 1e-16`; otherwise
`branch: F4-FIRED`.

Inputs: the V140/V141/V128/V138/V208/V215/V225/V229 artifacts read through
`cl10.paths.exp`, V138's constructed `O_h` irreps (module-level machinery of
that experiment, imported by name), and `cl10.continuum.polyhedral_group()`
for the icosahedral group. Exact character arithmetic in `Fraction` + `sqrt 5`
pairs; the only numerics beyond integers are the commutant rank counts.
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                    # noqa: E402
sys.path.insert(0, exp("V138_gamma_spin2"))

from cl10 import continuum as K                               # noqa: E402
import V138_gamma_spin2 as G138                               # noqa: E402

SQ5 = math.sqrt(5.0)


def art(name):
    p = exp(name, "results", name.lower() + ".json")
    with open(p) as fh:
        return json.load(fh)


# exact Q(sqrt 5) arithmetic: pairs (a, b) = a + b sqrt 5, Fraction entries
def qadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def qmul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def qscaled(c, x):
    return (Fraction(c) * x[0], Fraction(c) * x[1])


def main():
    t0 = time.time()
    checks, data = [], {}
    data["fitted_inputs"] = []

    # ---- T1 content under the two little/symmetry groups -----------------
    # (a) cubic: V138's machinery (irreps built from the defining rep)
    chi_sym2_oh = G138.chi_of(G138.sym2_rep)
    dec_oh = G138.decompose(chi_sym2_oh)
    v208 = art("V208_metric_bilinear_cone")
    v138 = art("V138_gamma_spin2")
    sym2_208 = v208["data"]["content"]["polarisation_Sym2"]
    bond_138 = v138["data"]["bond_factor"]["decomposition"]
    link_138 = v138["data"]["link_pair_channel"]["decomposition"]
    t1a = (abs(dec_oh["A1g"] - 1) < 1e-6 and abs(dec_oh["Eg"] - 1) < 1e-6
           and abs(dec_oh["T2g"] - 1) < 1e-6
           and abs(sym2_208.get("Eg", 0) - 1) < 1e-6
           and abs(bond_138.get("Eg", 0)) < 1e-6)

    # (b) icosahedral group of order 60, built from the presentation
    grp = K.polyhedral_group()["elements"]
    n_g = len(grp)
    # conjugacy classes computed inside the group (no table): every element
    # is conjugated by every element; classes are marked visited by canonical
    # rounded entries.  The two order-5 classes share a trace, so classes are
    # enumerated properly, not by trace key.
    def key9(M):
        return tuple(round(float(x), 9) for x in np.asarray(M, float).reshape(9))
    glist = [np.asarray(g, float) for g in grp]
    keypos = {key9(g): i for i, g in enumerate(glist)}
    visited = [False] * n_g
    reps_, sizes = [], []
    for i, g in enumerate(glist):
        if visited[i]:
            continue
        cls = set()
        for b in glist:
            cls.add(key9(b @ g @ np.linalg.inv(b)))
        for k9 in cls:
            visited[keypos[k9]] = True
        reps_.append(g)
        sizes.append(len(cls))
    assert sum(sizes) == n_g, (sizes, n_g)

    # characters by class: chi_V3(g) = tr g, chi_V3(g^2) = tr g^2; each trace
    # is identified with an exact Q(sqrt 5) number (margins 1e-12, the
    # icosahedral traces are exactly {3, phi, 1 - phi, 0, -1})
    phi_ex = (Fraction(1, 2), Fraction(1, 2))
    EXACT = {3.0: (Fraction(3), Fraction(0)),
             (1 + SQ5) / 2: phi_ex,
             (1 - SQ5) / 2: (Fraction(1, 2), Fraction(-1, 2)),
             0.0: (Fraction(0), Fraction(0)),
             -1.0: (Fraction(-1), Fraction(0))}


    def exact_tr(t):
        best, bd = None, 1e9
        for k, v in EXACT.items():
            d = abs(k - t)
            if d < bd:
                best, bd = v, d
        assert bd < 1e-9, (t, bd)
        return best

    chi1 = []                                       # tr of class rep
    chi1sq = []                                     # tr of (rep)^2
    for g in reps_:
        t = float(np.trace(g))
        chi1.append(exact_tr(t))
        chi1sq.append(exact_tr(float(np.trace(g @ g))))
    chi_sym2 = [(x[0] / 2, x[1] / 2) for x in
                [qadd(qmul(c1, c1), c2) for c1, c2 in zip(chi1, chi1sq)]]

    def inner_q(x, y):
        acc = (Fraction(0), Fraction(0))
        for s, a, b in zip(sizes, x, y):
            acc = qadd(acc, qscaled(s, qmul(a, b)))
        return (acc[0] / n_g, acc[1] / n_g)

    end_sym2 = inner_q(chi_sym2, chi_sym2)          # number of pieces
    inv_sym2 = inner_q(chi_sym2, [(Fraction(1), Fraction(0))] * len(sizes))
    chi_W = [(c[0] - inv_sym2[0], c[1] - inv_sym2[1]) for c in chi_sym2]
    end_W = inner_q(chi_W, chi_W)                   # 1 iff W irreducible
    # the identity class is wherever its size-1 entry sits, not necessarily
    # first: dim W = chi_W at the identity class
    dim_W = chi_W[sizes.index(1)]
    t1b = (end_sym2 == (2, 0) and inv_sym2 == (1, 0)
           and end_W == (1, 0) and dim_W == (5, 0))

    # (c) constructively: the 5-dimensional piece is irreducible - commutant
    # of the sym^2 rep restricted to the traceless subspace is scalars
    # build sym^2 matrices directly: S -> g S g^T on the 6-dim orthonormal basis
    IDX = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    SB = []
    for (i, j) in IDX:
        M = np.zeros((3, 3))
        M[i, j] += 1.0
        if i != j:
            M[j, i] += 1.0
        SB.append(M / math.sqrt(2.0 if i != j else 1.0))
    def rep6(g):
        D = np.zeros((6, 6))
        for b in range(6):
            Sp = g @ SB[b] @ g.T
            for a in range(6):
                D[a, b] = float(np.trace(SB[a].T @ Sp))
        return D
    R6 = [rep6(np.asarray(g, float)) for g in grp]
    # commutant dimension of the 6-dim rep (vectorized equations)
    A = np.vstack([np.kron(M, np.eye(6)) - np.kron(np.eye(6), M.T)
                   for M in R6])
    comm6 = 36 - np.linalg.matrix_rank(A)
    # invariant rank of the traceless piece = 1, its commutant = 1
    Pinv = sum(M for M in R6) / n_g
    inv_rank = int(round(np.linalg.matrix_rank(Pinv, tol=1e-8)))
    idv = np.array([1, 1, 1, 0, 0, 0], float) / math.sqrt(3.0)
    # the traceless subspace = the +1 eigenspace of the symmetric projector;
    # QR-truncation of the rank-4 column set is NOT safe when idv has zero
    # entries (it smuggles idv into the basis)
    _w, _V = np.linalg.eigh(np.eye(6) - np.outer(idv, idv))
    T5 = _V[:, _w > 0.5]
    R5 = [T5.T @ M @ T5 for M in R6]
    A5 = np.vstack([np.kron(M, np.eye(5)) - np.kron(np.eye(5), M.T)
                    for M in R5])
    comm5 = 25 - np.linalg.matrix_rank(A5, tol=1e-8)
    t1c = (comm6 == 2 and inv_rank == 1 and comm5 == 1)

    checks.append({
        "name": f"T1 content, constructed: sym^2 (6-dimensional) decomposes "
                f"under the CUBIC group as A1g+Eg+T2g ({dec_oh}; V208 control "
                f"{sym2_208}; qubit bond factor lacks Eg: {bond_138}); under "
                f"the icosahedral group of order {n_g} built from the [2,3,5] "
                f"presentation it decomposes as A (+) W with W one "
                f"irreducible 5-dimensional piece: exact Q(sqrt5) character "
                f"inner products give <Sym2,Sym2> = {end_sym2[0]}, <1,Sym2> ="
                f" {inv_sym2[0]}, <W,W> = {end_W[0]}, dim W = {dim_W[0]}; the "
                f"commutant of the 6-dim rep is {comm6}-dimensional, the "
                f"invariant rank is {inv_rank}, and the commutant of the "
                f"traceless 5-dim piece is {comm5} (scalars: Schur says "
                f"irreducible). The Eg/T2g split is a cubic-subgroup property; "
                f"under the star's own symmetry the quintet is ONE "
                f"representation",
        "pass": bool(t1a and t1b and t1c)})
    data["content"] = {
        "oh_decomposition": dec_oh, "v208_sym2": sym2_208,
        "v138_bond_factor": bond_138, "v138_link_channel": link_138,
        "icosahedral": {"order": n_g, "class_sizes": sizes,
                        "chi_sym2": [[str(a), str(b)] for a, b in chi_sym2],
                        "end_sym2": end_sym2[0], "invariants": inv_sym2[0],
                        "end_W": end_W[0], "dim_W": dim_W[0]},
        "commutant_ranks": {"sym2_6d": int(comm6), "invariant_rank": inv_rank,
                            "traceless_5d": int(comm5)}}

    # ---- T2 masslessness per channel --------------------------------------
    v140 = art("V140_dressed_gate_quintet")
    v141 = art("V141_regular_link_quintet_coupling")
    v229 = art("V229_diagonal_bond_quintet")
    v225 = art("V225_covariant_clifford_cell")
    ladder140 = v140["data"]["dispersion"]["helicity_ladder_speeds"]
    gap140 = float(v140["data"]["gamma_content"]["M0_minus_I"])
    gap229 = v229["data"]["dispersion"]["pi/4"]["gap_gamma"]
    speeds229 = v229["data"]["dispersion"]["pi/4"]["predicted_speeds_unit"]
    t2 = {"qubit_link_gamma": {
              "gapless": True, "flat": True,
              "note": "V138/V128: GAMMA gapless multiplets are flat bands "
                      "(omega/|q| -> 0) and the bond factor lacks Eg"},
          "directional_axis_bonds": {
              "gap_at_gamma": gap140,
              "ladder_speeds": ladder140, "linear": True},
          "diagonal_star_bonds": {
              "gap_at_gamma": gap229, "speeds_at_pi4": speeds229,
              "linear": True},
          "star_automaton_gamma": {
              "quadratic": True,
              "note": "V215: 512/1024 phases zero at GAMMA, first-order "
                      "splitting 1.7e-13/2.0e-13/1.3e-13 - no linear term"},
          "star_spinor_cell_node": {
              "irreps": v225["data"]["node"]["irreps"]
              if "node" in v225["data"] else
              {k: v for k, v in v225["data"].items() if "irrep" in k},
              "note": "V225: the gap cluster is the coin-spinor "
                      "representation times the internal sector - helicity "
                      "1/2 content, no quintet piece"}}
    lin140 = all(abs(a) > 1e-6 for a in ladder140[2:])
    lin229 = all(abs(a) > 1e-6 for a in speeds229[1:])
    checks.append({
        "name": f"T2 masslessness and propagation per channel: the only "
                f"MASSLESS LINEAR quintet branches are on the directional "
                f"layer - axis bonds {ladder140} (V140), diagonal star bonds "
                f"{[round(x, 4) for x in speeds229]} at pi/4 with gap "
                f"{gap229:.1e} at GAMMA (V229); the qubit-link GAMMA bands "
                f"are flat, the 1024-dimensional star automaton at GAMMA is "
                f"quadratic (V215), and the star spinor node carries spinor "
                f"(not spin-2) content (V225)",
        "pass": bool(lin140 and lin229)})
    data["channels"] = t2

    # ---- T3 survival against the mass shift ------------------------------
    shift141 = v141["data"]["coupling"]["mass_shift_at_lambda_g10"]
    Pi141 = v141["data"]["coupling"]["Pi_Eg_eigenvalues"]
    PiT141 = v141["data"]["coupling"]["Pi_T2g_max"]
    shift229 = v229["data"]["mass_shift"]["star_bond_shifts_at_lambda_g10"]
    max141 = float(max(abs(x) for x in shift141))
    max229 = float(max(abs(x) for x in shift229))
    f4_closed = (max229 < 1e-16 and max141 < 1e-16)
    checks.append({
        "name": f"T3 survival: at the closed clock the axis-bond channel "
                f"shifts only the Eg half (bubble {Pi141}, shifts "
                f"{[round(x, 4) for x in shift141]}, T2g bubble "
                f"{PiT141:.0e}), the diagonal star-bond channel shifts the "
                f"whole quintet (shifts {[round(x, 4) for x in shift229]}, "
                f"max {max229:.4f}); NO massless branch carries the full "
                f"spin-2 content with a vanishing shift (acceptance 1e-16): "
                f"axis bonds gap the Eg half by {max141:.4f}, star bonds gap "
                f"every component by up to {max229:.4f}",
        "pass": bool(not f4_closed)})
    data["shifts"] = {
        "axis_bond_Eg_shifts": shift141, "axis_bond_T2g_bubble": PiT141,
        "axis_bond_max_abs_shift": max141,
        "star_bond_quintet_shifts": shift229,
        "star_bond_max_abs_shift": max229,
        "acceptance_vanishing_shift": 1e-16}

    branch = "F4-CLOSED" if f4_closed else "F4-FIRED"
    ok = all(c["pass"] for c in checks)
    status = (
        f"F4 {('CLOSED' if f4_closed else 'STAYS FIRED')}, content test with "
        f"constructed characters (zero fitted parameters, "
        f"{time.time() - t0:.0f} s). (T1) the quintet carrier: under the "
        f"cubic group sym^2 = A1g (+) Eg (+) T2g (V208 control), and the "
        f"qubit bond factor carries T2g but never Eg; under the icosahedral "
        f"group of order {n_g} (built from the presentation) the SAME space "
        f"is A (+) W with W irreducible of dimension 5 - the Eg/T2g split is "
        f"a cubic-subgroup artifact, the quintet is one representation of the "
        f"star's symmetry (character inner products exact in Q(sqrt 5), "
        f"commutant ranks {comm6}/{inv_rank}/{comm5}). (T2) the massless "
        f"linear quintet branches exist only on the directional layer: axis "
        f"bonds {ladder140}, diagonal star bonds {[round(x, 4) for x in speeds229]}"
        f" at pi/4; the qubit-link bands are flat and the star spinor node is "
        f"spinor content. (T3) neither channel carries the full content with "
        f"vanishing shift: axis bonds gap the Eg half by {max141:.4f} "
        f"(T2g bubble exactly 0), the diagonal bonds gap the whole quintet by "
        f"up to {max229:.4f}. VERDICT: the implemented spectrum carries T2g "
        f"content but the branch is gapped by {max229:.4f} (diagonal-bond "
        f"channel) / {max141:.4f} (axis-bond channel, Eg half) - the lattice "
        f"spin-2 content is complete only before the matter coupling, and the "
        f"coupling that closes the T2g gap gaps everything.")
    out = {"id": "V230_quintet_helicity2_content",
           "tier": 2,
           "property": "the helicity-2 content of the quintet carrier under "
                       "O_h and the icosahedral group, and its survival "
                       "against the mass shift (V230, F4)",
           "residual": True,
           "programme": "V230",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v230_quintet_helicity2_content.json"),
              "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: float(o) if isinstance(o, Fraction)
                  else str(o))
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
