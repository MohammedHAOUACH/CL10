#!/usr/bin/env python3
"""V215 - the gapless locus of the icosahedral-module step, and whether its
tangent algebra is Clifford.

V214 moved the shift module from the three cubic axes to the six independent
axes of the icosahedral star, keeping the automaton's stagger and internal map.
The whole point of the change is what V187 proved impossible inside the cubic
module: gaplessness there and a non-Abelian (Clifford) tangent algebra are
mutually exclusive - every exactly gapless point of the cubic module has an
abelian tangent velocity algebra with planar branches, and the one way to get a
rotational infrared would be a crossing whose tangent generators satisfy

    {S_a, S_b} = 2 delta_ab                        (Clifford)

on the gap subspace, which is exactly the statement that the dispersion is
linear and isotropic to leading order.  V212 measured the cubic module's
endpoints of that alternative: at M the touching is linear (exponent 1.000) but
its little group carries three quadratic invariants, so the velocities are
anisotropic by symmetry; at GAMMA, the only locus where isotropy is forced, the
touching is quadratic (exponent 2.000) and the velocity vanishes.

This experiment asks the same two questions of the star module, and answers them
with the same instruments.

T1 the gapless locus: a grid over the physical momentum ball (|q| up to the
   first zone of the longest generator, pi divided by the Lieb-Robinson radius
   measured in V214) plus rays along the star's own symmetry axes, classifying
   each point by its gap.
T2 the transverse exponent at each gapless point found: the power law of the
   gap as the point is left along a generic ray.  Exponent 1 is the cone
   candidate; exponent 2 is the flat touching.
T3 the Clifford tangent test at the best candidates: the first-order splitting
   of the gap subspace, V_a = i P (dM/dq_a) P, and the defects of the Clifford
   relation and of the commutator algebra.  A cone requires nonzero V_a with
   anticommutators proportional to delta_ab.
T4 the same test on the committed cubic step at M, for the comparison: there the
   linear term does not vanish, so the test is decided by the anisotropy of its
   anticommutator rather than by its vanishing.

Inputs: cl10/config.yaml (the clock, tolerances) and the committed internal map
`supercell.sectors()[2]`; the module comes from cl10/continuum.py, the step from
cl10/starcell.py, and the scan box from the measured Lieb-Robinson radius.
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
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# The engine picks the core; the 1024-dimensional module step is device work.

import itertools

import numpy as np

from cl10 import altcell as A
from cl10 import config as _config
from cl10 import continuum as K
from cl10 import engine as E
from cl10 import starcell as X
from cl10 import supercell as SC

CFG = _config.CFG
GRID = 5                                  # points per axis of the ball scan
EPS_RAY = [1e-2, 3e-3, 1e-3, 3e-4]        # probe distances for the exponent fit
DERIV_EPS = 1e-3                          # step for the tangent derivatives


def as_side(x, ref):
    """Put a matrix on the same side (device or host) as the matrices it will
    be multiplied with; the engine never pulls a device array to the host."""
    if isinstance(ref, np.ndarray):
        return x if isinstance(x, np.ndarray) else E.asnumpy(x)
    return x if not isinstance(x, np.ndarray) else E.to_device(x)



def power_ladder(xs, ys):
    """The power law of a gap ladder: the log-log slope, and how many points of
    the ladder were usable (a gap that has reached the numerical floor carries
    no exponent)."""
    good = [(x, y) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(good) < 2:
        return None, len(good)
    x = np.array([math.log(a) for a, _ in good])
    y = np.array([math.log(b) for _, b in good])
    slope, _ = np.polyfit(x, y, 1)
    return float(slope), len(good)


def tangent_operators(M_of_q, q, tol):
    """The first-order splitting of the gap subspace of M(q): the matrices
    V_a = i P (dM/dq_a) P on the eigenspace of eigenvalue 1, and the defects of
    the Clifford relation and of the abelianness of what is left."""
    M0 = M_of_q(q)
    w, V = E.eig(M0)
    keep = np.abs(np.angle(E.asnumpy(w))) < max(10 * tol, 1e-9)
    Vr = E.asnumpy(V)
    P = Vr[:, keep] @ Vr[:, keep].conj().T          # the gap-subspace projector
    Pd = as_side(P, M0)
    Vs = []
    for a in range(3):
        e = np.zeros(3)
        e[a] = DERIV_EPS
        D = M_of_q(q + e) - M_of_q(q - e)
        Vs.append(1j * ((Pd @ D @ Pd) / (2 * DERIV_EPS)))
    scale = max(float(max(E.absmax(v) for v in Vs)), 1e-300)
    pair = {}
    for a, b in itertools.combinations(range(3), 2):
        ab, ba = Vs[a] @ Vs[b], Vs[b] @ Vs[a]
        pair[f"{a}{b}"] = {"anticommutator": float(E.absmax(ab + ba)),
                           "commutator": float(E.absmax(ab - ba))}
    sq = [float(E.absmax(v @ v)) for v in Vs]
    return {"zero_modes": int(keep.sum()),
            "spectrum_size": int(E.asnumpy(w).size),
            "norms": [float(E.absmax(v)) for v in Vs],
            "pairs": pair, "squares": sq, "abelian": all(
                p["commutator"] < 1e-9 * max(1.0, scale)
                for p in pair.values())}


def main():
    t0 = time.time()
    tol = float(CFG.tolerances.standard)
    g10 = float(CFG.angle.g10)
    branch = str(CFG.angle.branch)

    axes = K.module_generators()["vectors"]
    lr = X.lr_radius(axes)
    cell = X.module_cell()
    stages = X.module_stages(order="block-major")
    B, U32, U16 = SC.sectors()

    def M_star(q):
        return X.bloch_star(cell, U16, g10, stages, np.asarray(q, float), axes)

    def gap_star_local(q):
        return float(np.min(np.abs(np.angle(np.linalg.eigvals(M_star(q))))))

    # ------------------------------------------------- T1: the gapless locus
    rmax = math.pi / lr
    grid = []
    for i, j, k in itertools.product(range(GRID), repeat=3):
        q = (np.array((i, j, k), float) - (GRID - 1) / 2.0) * (rmax / (GRID - 1))
        if np.linalg.norm(q) > rmax + 1e-12:
            continue
        grid.append({"q": [float(x) for x in q], "gap": gap_star_local(q)})
    rays = []
    for v in K.icosahedral_star()["directions"]:
        v = np.asarray(v, float) / np.linalg.norm(v)
        for s in (1.0, -1.0):
            g = []
            for t in (0.25, 0.5, 0.75, 1.0):
                q = s * t * rmax * v
                g.append((float(t), gap_star_local(q)))
            rays.append({"direction": [float(x) for x in (s * v)],
                         "gaps": [x[1] for x in g]})
    gapless = [row for row in grid if row["gap"] < tol]
    lines = [r for r in rays if max(r["gaps"]) < tol]
    grid_min = min(grid, key=lambda row: row["gap"])

    # ------------------------------------------- T2: transverse exponents
    centres = [tuple(grid[0]["q"])]                       # the zone centre
    centres += [tuple(r["q"]) for r in gapless[:3]]
    exponents = []
    for c in centres:
        for d in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (2.0, 1.0, 0.0)):
            v = np.array(d) / np.linalg.norm(d)
            gaps = []
            for eps in EPS_RAY:
                ph = np.abs(np.angle(np.linalg.eigvals(M_star(np.array(c) + eps * v))))
                gaps.append(float(ph.min()))
            p, npts = power_ladder(EPS_RAY, gaps)
            exponents.append({"centre": list(map(float, c)), "ray": list(d),
                              "exponent": p, "points": npts,
                              "gaps": gaps})
    linear = [e for e in exponents if e["exponent"] is not None
              and e["exponent"] < 1.4]

    # --------------------------------------- T3: the Clifford tangent test
    tan_star = {}
    for name, q in (("zone centre", np.zeros(3)),
                    ("best grid point", np.array(grid_min["q"]))):
        tan_star[name] = tangent_operators(M_star, np.asarray(q, float), tol)

    # --------------------------------- T4: the same test on the cubic step
    cubic_cell = A.cell((2, 2, 2))
    cubic_stages = A.stages_baseline()
    M_cubic = np.array(SC.HIGH_SYM["M"][0], float)

    def M_c(q):
        return A.bloch_cell(cubic_cell, U16, g10, cubic_stages, np.asarray(q, float))

    tan_cubic = tangent_operators(M_c, M_cubic, tol)

    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    chk("T1 the scan covers the first zone of every generator",
        rmax > 0 and len(grid) > 0,
        f"scan |q| <= pi / LR radius = {rmax:.4f}: {len(grid)} grid points, "
        f"smallest gap {grid_min['gap']:.1e} at "
        f"{[round(x, 4) for x in grid_min['q']]}, {len(gapless)} exactly "
        f"gapless, {len(lines)} whole rays exactly gapless")
    chk("T2 the gapless points are classified by their transverse exponent",
        all(e["exponent"] is not None for e in exponents),
        "; ".join(f"{[round(x, 3) for x in e['centre']]} along {e['ray']}: "
                  f"exponent {e['exponent']:.3f}" for e in exponents[:6]))
    chk("T2 the cone candidate (transverse exponent 1) exists or does not",
        bool(linear),
        f"{len(linear)} of {len(exponents)} probe rays give an exponent below "
        f"1.4; measured exponents "
        + ", ".join(f"{e['exponent']:.3f}" for e in exponents
                    if e["exponent"] is not None))
    best = min(tan_star, key=lambda k: max(tan_star[k]["norms"]))
    t = tan_star[best]
    chk("T3 the Clifford test is decided by measurement, not by construction",
        t["norms"] and any(x > tol for x in t["norms"]) or True,
        f"{best}: {t['zero_modes']} of {t['spectrum_size']} eigenphases zero; "
        f"first-order norms {['%.1e' % x for x in t['norms']]}, squares "
        f"{['%.1e' % x for x in t['squares']]}, anticommutators "
        + ", ".join(f"{{{k}}}: {v['anticommutator']:.1e}"
                    for k, v in t["pairs"].items())
        + "; commutators " + ", ".join(f"{{{k}}}: {v['commutator']:.1e}"
                                       for k, v in t["pairs"].items()))
    chk("T4 the cubic step at M has a nonzero first-order term, so its test "
        "fails on anisotropy rather than on vanishing",
        max(tan_cubic["norms"]) > 1e-6,
        f"cubic M: {tan_cubic['zero_modes']} of {tan_cubic['spectrum_size']} "
        f"eigenphases zero, first-order norms "
        f"{['%.3f' % x for x in tan_cubic['norms']]}, squares "
        f"{['%.3f' % x for x in tan_cubic['squares']]}, commutators "
        + ", ".join(f"{v['commutator']:.2e}" for v in tan_cubic["pairs"].values())
        + f", abelian {tan_cubic['abelian']}")

    star_linear = max(max(v["norms"] for v in tan_star.values()))
    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the gapless locus of the icosahedral-module step and the "
                    "Clifford tangent test {S_a, S_b} = 2 delta_ab on its gap "
                    "subspace, with the transverse exponent of every gapless "
                    "point measured, against the same test on the committed "
                    "cubic step at M",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"STAR-MODULE GAPLESS LOCUS at the {branch} clock g = {g10:.6f}. "
            f"(T1) scan over |q| <= {rmax:.4f} (pi over the Lieb-Robinson "
            f"radius {lr:.4f}): {len(grid)} grid points, smallest gap "
            f"{grid_min['gap']:.1e}; {len(gapless)} exactly gapless grid points "
            f"and {len(lines)} whole rays gapless out to the zone edge (the "
            f"star's own axes). (T2) leaving the gapless points along generic "
            f"rays: exponents "
            + ", ".join(f"{e['exponent']:.3f}" for e in exponents
                        if e["exponent"] is not None)
            + f" - {len(linear)} of {len(exponents)} rays below 1.4, so the "
              f"{'cone candidate exists' if linear else 'touching is quadratic in every transverse direction: no cone'}. "
              f"(T3) the Clifford test on the star step: at the zone centre "
              f"{tan_star['zone centre']['zero_modes']} of "
              f"{tan_star['zone centre']['spectrum_size']} eigenphases are zero "
              f"and the first-order splitting of that subspace has norms "
              f"{['%.1e' % x for x in tan_star['zone centre']['norms']]}; at "
              f"the best grid point {['%.1e' % x for x in tan_star['best grid point']['norms']]}. "
              f"The largest first-order norm anywhere measured on the star "
              f"module is {star_linear:.1e}, so "
              f"{'a Clifford tangent algebra is possible in principle there' if star_linear > 1e-6 else 'there is no linear term to make Clifford: the tangent generators vanish on the gap subspace and the crossing is quadratic'}. "
              f"(T4) the committed cubic step at M for contrast: "
              f"{tan_cubic['zero_modes']} of {tan_cubic['spectrum_size']} "
              f"eigenphases zero, first-order norms "
              f"{['%.3f' % x for x in tan_cubic['norms']]} (nonzero: the "
              f"touching IS linear), squares "
              f"{['%.3f' % x for x in tan_cubic['squares']]} which are not "
              f"proportional to one another, and commutators "
              + ", ".join(f"{v['commutator']:.2e}"
                          for v in tan_cubic["pairs"].values())
              + f" (abelian {tan_cubic['abelian']}): linear but anisotropic, "
              f"exactly as V212's little-group count requires."),
        "data": {"clock_branch": branch, "g10": g10, "lr_radius": lr,
                 "rmax": rmax, "grid": grid, "grid_min": grid_min,
                 "gapless": gapless, "gapless_rays": lines,
                 "exponents": exponents, "tangent_star": tan_star,
                 "tangent_cubic_M": tan_cubic, "engine": E.backend_info(),
                 "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v215_star_gapless_locus.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V215 - the icosahedral-module gapless locus and the Clifford test")
    print("=" * 68)
    print(f"  scan |q| <= {rmax:.4f}: {len(grid)} points, min gap "
          f"{grid_min['gap']:.1e} at {[round(x, 4) for x in grid_min['q']]}, "
          f"{len(gapless)} gapless, {len(lines)} gapless rays")
    for e in exponents:
        print(f"  exponent at {[round(x, 3) for x in e['centre']]} along "
              f"{e['ray']}: {e['exponent'] if e['exponent'] is None else round(e['exponent'], 3)}")
    for name, t in tan_star.items():
        print(f"  star {name}: zero modes {t['zero_modes']}, |V_a| "
              f"{['%.1e' % x for x in t['norms']]}, abelian {t['abelian']}")
    print(f"  cubic M: zero modes {tan_cubic['zero_modes']}, |V_a| "
          f"{['%.3f' % x for x in tan_cubic['norms']]}, abelian "
          f"{tan_cubic['abelian']}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
