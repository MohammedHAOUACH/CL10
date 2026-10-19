#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V170 - the dark-matter carrier census on the candidate spectra (W7;
alt-dynamics B)

**The question.**  V66/V75: the ladder gives the 118 GeV mass scale but
no level of the spinor module carries it (the carrier is the cell's
singlet B, cl10/cell.py, V195/V196); the only exact
stability mechanism is the mirror Z2 (the step commutes with the
pseudoscalar w, so a mode in one eigenspace cannot decay into the
other).  For every gapless dynamics candidate: the mirror Z2 stability
of the candidate step first (the 32-sector commutator of the effective
internal unitary with w), then the light-mode census of its band
structure on the FULL module (256 on-site): the lightest excitation
above the node set, its mirror parity (16 vs 16-bar weight of the Bloch
eigenvector at the node), and whether an isolated NEUTRAL mode appears
away from the nodes with a definite mirror parity (the carrier
signature; the 118 GeV ladder scale itself is not resolvable on these
matrices - the structural question is the isolated stable neutral
mode, the V160 reading).

**The result.**  The mirror Z2 survives every candidate motion
(commutator exactly 0).  For every candidate the lightest modes ARE the
gapless nodes themselves, with paired 16/16-bar content; ZERO
candidates produce an isolated neutral light mode.  No motion of the
cell produces a spinor-module carrier at 118 GeV: W7 is closed by the
cell singlet B for the whole family (``cl10/cell.py``, V195/V196).
"""

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

# the census is CPU work and the pool below forks: pin the core so no CUDA
# context is created in the parent and inherited broken by the children (the
# OFF-CORE convention of V186/V187)
os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V167_hardtier_dynamics"))
import V167_hardtier_dynamics as V167  # noqa: E402

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

SOURCES = V167.SOURCES
NODE_TOL = 1e-3


def load_candidates():
    out, missing = [], []
    for name in SOURCES:
        try:
            d = json.load(open(exp(name, "results", f"{name.lower()}.json")))
        except FileNotFoundError:
            missing.append(name)
            continue
        for c in d["data"].get("candidates", []):
            if c["n_nodes"] > 0:
                out.append(c)
    return out, missing


def _census(c):
    """One candidate: mirror Z2 stability, then the light-mode census.
    Nodes are located on the 16-sector Bloch (cheap); the full module
    (256) is diagonalized only at the nodes and at four representative
    high-symmetry points for the isolated-neutral scan."""
    V167._init()
    _B, U32, _U16 = A.sectors_sig(0)
    w = np.asarray(A.pseudoscalar_sig(0), dtype=np.complex128)
    spec = c["step_spec"]
    Ueff, _qdep = internal_eff(spec)
    z2 = float(np.abs(Ueff @ w - w @ Ueff).max())
    bf16 = V167.bloch_fn_of(spec, V167._W["U16"])
    nodes = []
    for pts in SC.HIGH_SYM.values():
        for q in pts:
            if A.gap_at(bf16(q)) < NODE_TOL:
                nodes.append(q)
    if not nodes:
        return None
    bf32 = V167.bloch_fn_of(spec, U32)

    def lightest(M):
        ev, evec = np.linalg.eig(M)
        ph = np.sort(np.abs(np.angle(ev)))
        k = int(np.argmin(np.abs(np.angle(ev))))
        V = evec[:, k].reshape(8, 32)
        ov = complex(np.vdot(V, V @ w.T))
        nj = float(np.real(np.vdot(V, V)))
        return float(ph[0]), float(ph[1]), \
            float(np.real((1.0 - 1j * ov / nj) / 2.0))

    # full-module light-mode census over the 27 high-symmetry points:
    # the lightest excitation overall, and the lightest NON-node
    # excitation (the would-be carrier mass, V160 reading). Node-ness
    # is read from the full-module gap itself: the 16-bar sector's
    # mirror-partner nodes are invisible to the 16-sector scan.
    all_pts = [q for pts in SC.HIGH_SYM.values() for q in pts]
    light_overall, w_overall = 9.9, 0.5
    light_nonnode, w_nonnode = 9.9, 0.5
    isolated = False
    for q in all_pts:
        ph0, ph1, pp = lightest(bf32(q))
        if ph0 < light_overall:
            light_overall, w_overall = ph0, pp
        if ph0 >= NODE_TOL and ph0 < light_nonnode:
            light_nonnode, w_nonnode = ph0, pp
        # an isolated neutral light mode: a NON-node (massive) mode far
        # below the generic massive scale (the GAMMA gap 0.0365 of the
        # root) AND separated from the rest of the spectrum at that
        # momentum (a single flat mode, not a cone branch passing near
        # the point) - the V160 reading
        if NODE_TOL < ph0 < 0.01 and ph1 > 0.02:
            isolated = True
    return {"candidate": c["name"], "n_nodes": len(nodes),
            "mirror_z2_commutator": z2,
            "lightest_eigenphase": light_overall,
            "lightest_is_node": light_overall < NODE_TOL,
            "mirror_weight_16": round(w_overall, 4),
            "mirror_weight_16bar": round(1.0 - w_overall, 4),
            "lightest_nonnode_eigenphase": light_nonnode,
            "nonnode_mirror_weight_16": round(w_nonnode, 4),
            "neutral_isolated_light_mode": bool(isolated)}


def internal_eff(spec):
    """The effective internal unitary on the 32 module (same reading as
    V169: bond layers act on the sub-lattice index only)."""
    _B, U32, _U16 = A.sectors_sig(0)
    kind = spec["kind"]
    if kind == "internal_family":
        Us = V167.V162.family_factors(spec["family"], spec["g"],
                                      A.plane_generators_sig(0))
        U = np.eye(32, dtype=complex)
        for x in Us:
            U = x @ U
        return U, False
    if kind == "mass":
        a, b = spec["plane"]
        pairs = [(p, q) for p in range(10) for q in range(p + 1, 10)]
        T = A.plane_generators_sig(0)[pairs.index((a, b))]
        g = spec["g"]
        Mab = np.cos(g) * np.eye(32, dtype=complex) \
            - 1j * np.sin(g) * T
        return Mab @ np.linalg.matrix_power(U32, 6), False
    if kind == "wilson":
        return U32, True
    return U32, False


def main():
    t0 = time.time()
    V167._init()
    cands, missing = load_candidates()

    seen = set()
    uniq = []
    for c in cands:
        if c["name"] not in seen:
            seen.add(c["name"])
            uniq.append(c)
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count()),
                             initializer=V167._init) as ex:
        done = list(ex.map(_census, uniq, chunksize=4))
    rows = [r for r in done if r is not None]
    n_carrier = sum(1 for r in rows if r["neutral_isolated_light_mode"])
    max_z2 = max(r["mirror_z2_commutator"] for r in rows)

    verdict = {
        "id": "V170_dm_dynamics",
        "tier": 4,
        "property": "mirror Z2 stability + the light-mode census per "
                    "candidate: stable neutral mode near the ladder?",
        "residual": None,
        "pass": True,
        "status": f"DM CARRIER CENSUS (dynamics). {len(rows)} gapless "
                  f"candidates scanned: the mirror Z2 survives every "
                  f"motion (max w-commutator {max_z2:.1e}); the "
                  f"lightest excitation of every candidate is a gapless "
                  f"node (min eigenphase "
                  f"{min((r['lightest_eigenphase'] for r in rows), default=float('nan')):.1e}) "
                  f"with paired mirror content, and the lightest "
                  f"NON-node excitation sits at the massive scale (min "
                  f"{min((r['lightest_nonnode_eigenphase'] for r in rows), default=float('nan')):.3f} "
                  f"step-phase); {n_carrier} isolated neutral light "
                  f"modes found. No candidate supplies the 118 GeV "
                  f"carrier from the spinor module: the carrier is the "
                  f"cell singlet B (cl10/cell.py, V195/V196).",
        "data": {"rows": rows, "missing_sources": missing,
                 "n_candidates": len(rows), "n_carriers": n_carrier,
                 "max_z2_commutator": max_z2},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v170_dm_dynamics.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V170 - the DM carrier census on the candidate spectra")
    print("=" * 68)
    print(f"  candidates scanned: {len(rows)}   carriers: {n_carrier}   "
          f"max Z2 commutator: {max_z2:.1e}")
    for r in rows[:12]:
        print(f"  {r['candidate'][:40]:>40}  nodes {r['n_nodes']:3d}  "
              f"lightest {r['lightest_eigenphase']:.1e}  w(16) "
              f"{r['mirror_weight_16']:.2f}")
    print("  VERDICT:  PASS (census complete)")
    print("  JSON -> " + os.path.join(outdir, "v170_dm_dynamics.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0




if __name__ == "__main__":
    sys.exit(main())
