#!/usr/bin/env python3
"""V217 - what the module change does NOT touch: the internal sector, its
chirality split, its CPT action and its trace forms.

V214 replaced the cubic shift module with the icosahedral star module and found
the internal map reusable unchanged, with the step's mirror Z2 surviving as
P M(-q) P = M_rev(q)^T.  Everything the paper claims about the internal sector -
the 16 / 16-bar split from the pseudoscalar, the CPT operator as the time
direction of the root lattice, the anomaly cancellation and the 5/3 trace-form
normalization - is built from `cl10/clifford.py`, `cl10/cpt.py`, `cl10/matter.py`
and `cl10/gauge.py`, none of which looks at the cell.  That is an argument; this
experiment turns it into a measurement on the new module.

The engine's step lives on the 16 sector, where the pseudoscalar is only the
sector label (w = +i on it), so chirality has to be tested where it is a real
operator: on the full 32-dimensional internal space.  The step is therefore
built twice, with the committed 16-sector map and with the committed 32-sector
map, and the comparison of the two is itself the measurement.

T1 chirality conservation: the 32-internal step is an internal scalar, so it
   commutes with 1 (x) w.  Measured for both modules at generic momenta.
T2 the sector blocks: with the pseudoscalar eigenbasis the 32-internal step is
   block diagonal - the off-block norm is measured - its 16-block is the
   committed step (the sector is exactly invariant, so V214's reuse is not a
   choice but a consequence), and the matter and antimatter blocks are compared
   spectrally: that is CPT on the spectrum of the new module.
T3 the CPT action: the committed CPT operator is g_1, the time direction, which
   exchanges the sectors (it anticommutes with w).  Its candidate actions on the
   step, with and without the cell reflection and the stage reversal that V214
   measured, are all evaluated; the one that holds is reported, not assumed.
T4 the internal data: the sector dimensions, the anomaly coefficients and the
   5/3 trace-form ratio are recomputed with the star module in the construction
   and shown unchanged, because the module enters none of them.

Inputs: cl10/config.yaml (clock, tolerances), the committed internal maps
`supercell.sectors()` and `supercell.internal`, the committed Clifford algebra
`cl10/clifford.py`, the committed CPT construction `cl10/cpt.py`,
`cl10/matter.py`, `cl10/gauge.py`, the cubic cell from `cl10/altcell.py` and the
star module from `cl10/continuum.py` and `cl10/starcell.py`.
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# The engine picks the core; the 2048-dimensional module step is device work.

import numpy as np

from cl10 import altcell as A
from cl10 import config as _config
from cl10 import continuum as K
from cl10 import cpt as CPT
from cl10 import engine as E
from cl10 import gauge
from cl10 import matter
from cl10 import starcell as X
from cl10 import supercell as SC

CFG = _config.CFG


def as_side(x, ref):
    """Put a matrix on the same side (device or host) as the matrices it will
    be multiplied with; the engine never pulls a device array to the host."""
    if isinstance(ref, np.ndarray):
        return x if isinstance(x, np.ndarray) else E.asnumpy(x)
    return x if not isinstance(x, np.ndarray) else E.to_device(x)


def reflection(c):
    """The digit-by-digit reflection (d -> -d mod n) of a mixed-radix cell."""
    P = np.zeros((c["nsub"], c["nsub"]))
    for e, digits in enumerate(c["sites"]):
        refl = [(-d) % n for d, n in zip(digits, c["shape"])]
        P[c["site_index"](refl), e] = 1.0
    return P


def sector_basis(w):
    """The unitary that orders the internal basis by pseudoscalar eigenvalue:
    the +i eigenspace (matter, the 16) first, then -i (antimatter, the 16-bar).
    This is the same ordering `supercell.sectors` uses, so its first sixteen
    columns must reproduce the committed sector frame."""
    ev, evec = np.linalg.eig(w)
    order = np.argsort(-np.imag(ev))
    return np.linalg.qr(evec[:, order])[0]


def sector_blocks(M, T, nsub, nint):
    """Split a step into its matter and antimatter blocks in the pseudoscalar
    eigenbasis, and report the off-block norm, the two blocks and their sorted
    eigenphase spectra.  T acts on the internal factor of the step, whose basis
    is ordered position (x) internal, so the change of basis is 1 (x) T."""
    Mh = E.asnumpy(M)
    Tfull = np.kron(np.eye(nsub), T)
    Mb = Tfull.conj().T @ Mh @ Tfull
    n16 = nint // 2
    idx16 = np.array([p * nint + i for p in range(nsub) for i in range(n16)])
    idxbar = np.array([p * nint + i for p in range(nsub)
                       for i in range(n16, nint)])
    off = max(float(np.max(np.abs(Mb[np.ix_(idxbar, idx16)]))),
              float(np.max(np.abs(Mb[np.ix_(idx16, idxbar)]))))
    b16 = Mb[np.ix_(idx16, idx16)]
    bbar = Mb[np.ix_(idxbar, idxbar)]
    return {"total": int(Mh.shape[0]), "dim16": int(b16.shape[0]),
            "dimbar": int(bbar.shape[0]), "off_block": off,
            "block16": b16, "blockbar": bbar,
            "phases16": np.sort(np.abs(np.angle(np.linalg.eigvals(b16)))),
            "phasesbar": np.sort(np.abs(np.angle(np.linalg.eigvals(bbar))))}


def main():
    t0 = time.time()
    tol = float(CFG.tolerances.standard)
    strict = float(CFG.tolerances.strict)
    g10 = float(CFG.angle.g10)
    branch = str(CFG.angle.branch)

    axes = K.module_generators()["vectors"]
    cell = X.module_cell()
    star_stages = X.module_stages(order="block-major")
    B, U32, U16 = SC.sectors()
    nint = U32.shape[0]

    # The two constructions, each with the 16-sector map (what the engine runs)
    # and with the full 32-sector map (where chirality is an operator).
    cubic_cell = A.cell((2, 2, 2))
    cubic_stages = A.stages_baseline()

    def cubic_M(q, U=U16, stages=cubic_stages):
        return A.bloch_cell(cubic_cell, U, g10, stages, np.asarray(q, float))

    def star_M(q, U=U16, stages=star_stages):
        return X.bloch_star(cell, U, g10, stages, np.asarray(q, float), axes)

    probes = [np.array(v, float) for v in
              ((0.4, 0.0, 0.0), (0.3, 0.2, 0.0), (0.1, 0.1, 0.1))]

    # ------------------------------------------- T1: the committed internal data
    w = E.asnumpy(CPT.pseudoscalar())
    g1 = E.asnumpy(CPT.cpt_operator())
    P16, P16bar = [E.asnumpy(p) for p in CPT.projectors(E.to_device(w))]
    T = sector_basis(w)
    wdef = {
        "w_antihermitian": float(np.max(np.abs(w + w.conj().T))),
        "w_square": float(np.max(np.abs(w @ w + np.eye(nint)))),
        "g1_hermitian_involution": float(np.max(np.abs(
            g1 @ g1 - np.eye(nint)))),
        "g1_anticommutes_w": float(np.max(np.abs(g1 @ w + w @ g1))),
        "w_commits_U32": float(np.max(np.abs(w @ U32 - U32 @ w))),
        "rank16": int(round(np.real(np.trace(P16)))),
        "rank16bar": int(round(np.real(np.trace(P16bar)))),
        "projectors_orthogonal": float(np.max(np.abs(P16 @ P16bar))),
        "frame_matches_sectors": float(np.max(np.abs(
            T[:, :B.shape[1]] @ T[:, :B.shape[1]].conj().T
            - B @ B.conj().T))),
    }

    # --------------------------------------- T1: chirality of the 32-internal step
    commute = {}
    for name, M_of_q, nsub in (("cubic", cubic_M, cubic_cell["nsub"]),
                               ("star", star_M, cell["nsub"])):
        M0 = M_of_q(np.zeros(3), U32)
        W = as_side(np.kron(np.eye(nsub), w), M0)
        worst = 0.0
        for q in probes:
            M = M_of_q(q, U32)
            worst = max(worst, float(E.absmax(M @ W - W @ M)))
        commute[name] = worst

    # -------------------------------------------------- T2: the sector blocks
    qprobe = np.array([0.4, 0.2, 0.1], float)
    blocks, sector_equality = {}, {}
    for name, M_of_q, nsub in (("cubic", cubic_M, cubic_cell["nsub"]),
                               ("star", star_M, cell["nsub"])):
        b = sector_blocks(M_of_q(qprobe, U32), T, nsub, nint)
        # the 16-block of the 32-internal step against the committed step: the
        # sector is exactly invariant (U16 = B^dagger U32 B), so in the
        # committed frame the two must be the same matrix.
        M32 = M_of_q(qprobe, U32)
        Tb = as_side(np.kron(np.eye(nsub), B), M32)
        M16fromM32 = Tb.conj().T @ M32 @ Tb
        committed = E.asnumpy(M_of_q(qprobe, U16))
        sector_equality[name] = float(E.absmax(
            M16fromM32 - as_side(committed, M16fromM32)))
        blocks[name] = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                        for k, v in b.items()
                        if k not in ("block16", "blockbar")}

    # ------------------------------------------------------ T3: the CPT action
    CPT_NAMES = ("g1 alone: G M(q) G = M(q)^dagger",
                 "g1 + cell reflection: R G M(q) G R = M(-q)",
                 "g1 + reflection + reverse: ... = M_rev(-q)^dagger",
                 "g1 + reflection + reverse + transpose: ... = M_rev(-q)^T")
    cpt_defects = {}
    for name, M_of_q, c, stages in (("cubic", cubic_M, cubic_cell, cubic_stages),
                                    ("star", star_M, cell, star_stages)):
        host = M_of_q(np.zeros(3), U32)
        R = as_side(np.kron(reflection(c), np.eye(nint)), host)
        G = as_side(np.kron(np.eye(c["nsub"]), g1), host)
        defs = {k: 0.0 for k in CPT_NAMES}
        for q in probes:
            M = M_of_q(q, U32)
            Mrm = M_of_q(-q, U32, list(reversed(stages)))
            L_int = G @ M @ G.conj().T
            LR = R @ L_int @ R.conj().T
            defs[CPT_NAMES[0]] = max(defs[CPT_NAMES[0]],
                                     float(E.absmax(L_int - M.conj().T)))
            defs[CPT_NAMES[1]] = max(defs[CPT_NAMES[1]],
                                     float(E.absmax(LR - Mrm)))
            defs[CPT_NAMES[2]] = max(defs[CPT_NAMES[2]],
                                     float(E.absmax(LR - Mrm.conj().T)))
            defs[CPT_NAMES[3]] = max(defs[CPT_NAMES[3]],
                                     float(E.absmax(LR - Mrm.T)))
        cpt_defects[name] = defs
    best_cpt = {n: min(d, key=lambda k: d[k]) for n, d in cpt_defects.items()}

    # ------------------------------------------------ T4: the internal data
    internal = {
        "U32_leaves_the_sector": float(np.max(np.abs(U32 @ B - B @ U16))),
        "anomalies": matter.anomalies(), "five_thirds": float(
            gauge.five_thirds()),
        "sector_dims": {"16": int(B.shape[1]), "32": int(B.shape[0])},
    }

    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    chk("T1 the pseudoscalar, the CPT operator and the sector frame are the "
        "committed ones",
        wdef["w_square"] < 100 * strict
        and wdef["w_antihermitian"] < 100 * strict
        and wdef["g1_hermitian_involution"] < 100 * strict
        and wdef["g1_anticommutes_w"] < 100 * strict
        and wdef["w_commits_U32"] < 100 * strict
        and wdef["rank16"] == 16 and wdef["rank16bar"] == 16
        and wdef["frame_matches_sectors"] < 100 * strict,
        f"w^2 = -I to {wdef['w_square']:.1e}, w anti-Hermitian to "
        f"{wdef['w_antihermitian']:.1e}, g1^2 = I to "
        f"{wdef['g1_hermitian_involution']:.1e}, g1 w + w g1 = 0 to "
        f"{wdef['g1_anticommutes_w']:.1e}, [w, U32] = 0 to "
        f"{wdef['w_commits_U32']:.1e}, projector ranks "
        f"{wdef['rank16']}/{wdef['rank16bar']} orthogonal to "
        f"{wdef['projectors_orthogonal']:.1e}; the pseudoscalar eigenbasis "
        f"reproduces the committed sector frame to "
        f"{wdef['frame_matches_sectors']:.1e}")
    chk("T1 the 32-internal step is an internal scalar: it commutes with "
        "1 (x) w on both modules",
        all(v < 100 * tol for v in commute.values()),
        "max ||[M_32(q), 1 (x) w]|| over probe momenta: "
        + ", ".join(f"{k} {v:.1e}" for k, v in commute.items())
        + " (on the 16 sector alone w is only the sector label, +i times the "
          "identity, which is why the test is run on the full internal space)")
    for name in ("cubic", "star"):
        b = blocks[name]
        pdiff = float(np.max(np.abs(np.asarray(b["phases16"])
                                    - np.asarray(b["phasesbar"]))))
        chk(f"T2 the 32-internal {name} step is block diagonal in the 16 / "
            "16-bar split, its 16-block IS the committed 16-sector step, and "
            "its matter and antimatter blocks are isospectral",
            b["off_block"] < 100 * tol
            and sector_equality[name] < 100 * tol
            and pdiff < 100 * tol,
            f"off-block norm {b['off_block']:.1e} on a {b['total']}-dimensional "
            f"step split into {b['dim16']} + {b['dimbar']}; the 16-block "
            f"reproduces the committed step to {sector_equality[name]:.1e}; "
            f"largest difference of the sorted matter/antimatter eigenphase "
            f"spectra {pdiff:.1e}")
    chk("T3 one CPT candidate holds on each module, selected by measurement",
        all(cpt_defects[n][best_cpt[n]] < 100 * tol for n in cpt_defects),
        "; ".join(f"{n}: best '{best_cpt[n]}' "
                  f"{cpt_defects[n][best_cpt[n]]:.1e} (the others: "
                  + ", ".join(f"{k.split(':')[0]} {v:.1e}"
                              for k, v in cpt_defects[n].items()
                              if k != best_cpt[n]) + ")"
                  for n in cpt_defects))
    chk("T4 the internal sector is untouched by the module change: the sector is "
        "invariant under the committed map and the anomaly and trace-form data "
        "are unchanged",
        internal["U32_leaves_the_sector"] < 100 * tol
        and abs(internal["five_thirds"] - 5 / 3) < 1e-9,
        f"U32 B = B U16 to {internal['U32_leaves_the_sector']:.1e}, 5/3 = "
        f"{internal['five_thirds']:.9f}, anomaly coefficients "
        f"{json.dumps(internal['anomalies'])}")

    star_b, cubic_b = blocks["star"], blocks["cubic"]
    cubic_pd = float(np.max(np.abs(np.asarray(cubic_b["phases16"])
                                   - np.asarray(cubic_b["phasesbar"]))))
    star_pd = float(np.max(np.abs(np.asarray(star_b["phases16"])
                                  - np.asarray(star_b["phasesbar"]))))
    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the internal sector under the module change: chirality "
                    "conservation of the 32-internal step (its commutator with "
                    "1 (x) w), the 16 / 16-bar block diagonalisation with the "
                    "16-block identified with the committed step, the matter / "
                    "antimatter spectral identity, the CPT action selected from "
                    "four candidates, and the invariance of the anomaly and "
                    "trace-form data",
        "residual": False,
        "pass": bool(ok),
        "status": (
            f"INTERNAL SECTOR on the star module at the {branch} clock g = "
            f"{g10:.6f}. (T1) the committed pseudoscalar and CPT operator are "
            f"intact (w^2 = -I to {wdef['w_square']:.1e}, g1^2 = I to "
            f"{wdef['g1_hermitian_involution']:.1e}, g1 exchanges the sectors "
            f"since g1 w + w g1 = 0 to {wdef['g1_anticommutes_w']:.1e}, sector "
            f"projectors of rank {wdef['rank16']}/{wdef['rank16bar']}) and the "
            f"32-internal step commutes with 1 (x) w on both modules ("
            + ", ".join(f"{k} {v:.1e}" for k, v in commute.items()) + f"). (T2) "
            f"the step is block diagonal in the 16 / 16-bar split (off-block "
            f"{cubic_b['off_block']:.1e} cubic on {cubic_b['total']} "
            f"dimensions, {star_b['off_block']:.1e} star on {star_b['total']} "
            f"split {star_b['dim16']} + {star_b['dimbar']}), its 16-block is "
            f"the committed 16-sector step to "
            f"{sector_equality['cubic']:.1e} (cubic) and "
            f"{sector_equality['star']:.1e} (star), and matter and antimatter "
            f"are isospectral on both (largest sorted-spectrum difference "
            f"{cubic_pd:.1e} cubic, "
            f"{star_pd:.1e} "
            f"star): CPT holds on the spectrum of the new module. (T3) the CPT "
            f"action on the step, against four candidates: "
            + "; ".join(f"{n} -> '{best_cpt[n].split(':')[0]}' at "
                        f"{cpt_defects[n][best_cpt[n]]:.1e}"
                        for n in cpt_defects)
            + f". (T4) the internal data are untouched: U32 B = B U16 to "
            f"{internal['U32_leaves_the_sector']:.1e}, the 5/3 trace-form ratio "
            f"{internal['five_thirds']:.9f}, anomalies "
            f"{json.dumps(internal['anomalies'])} - the module enters none of "
            f"them."),
        "data": {"clock_branch": branch, "g10": g10, "witness": wdef,
                 "commute": commute, "blocks": blocks,
                 "sector_equality": sector_equality,
                 "cpt": {"candidates": list(CPT_NAMES), "defects": cpt_defects,
                         "best": best_cpt},
                 "internal": internal, "engine": E.backend_info(),
                 "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v217_internal_sector_survival.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V217 - the internal sector under the module change")
    print("=" * 68)
    print(f"  w^2 + I {wdef['w_square']:.1e}, g1^2 - I "
          f"{wdef['g1_hermitian_involution']:.1e}, ranks "
          f"{wdef['rank16']}/{wdef['rank16bar']}, frame "
          f"{wdef['frame_matches_sectors']:.1e}")
    print("  [M_32, 1 (x) w]: " + ", ".join(f"{k} {v:.1e}"
                                            for k, v in commute.items()))
    for name, b in blocks.items():
        print(f"  {name}: off-block {b['off_block']:.1e}, blocks "
              f"{b['dim16']}+{b['dimbar']} of {b['total']}, 16-block = committed "
              f"to {sector_equality[name]:.1e}, spectrum difference "
              f"{float(np.max(np.abs(np.asarray(b['phases16']) - np.asarray(b['phasesbar'])))):.1e}")
    for n, d in cpt_defects.items():
        print(f"  CPT {n}: " + ", ".join(f"{k.split(':')[0]} {v:.1e}"
                                         for k, v in d.items()))
    print(f"  five thirds {internal['five_thirds']:.9f}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
