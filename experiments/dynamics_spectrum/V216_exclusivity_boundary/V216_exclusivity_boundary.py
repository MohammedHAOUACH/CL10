#!/usr/bin/env python3
"""V216 - the boundary of the V187 mutual exclusivity: can the step be gapless
and have a non-Abelian tangent algebra at the same time?

V187 proved the exclusivity inside the implemented cubic module, and V212
measured its two endpoints on the ground: at M the touching is linear but the
little group carries quadratic invariants, so the tangent velocities are
anisotropic and commute; at GAMMA isotropy is forced and the touching is
quadratic, so the velocities vanish.  A rotational infrared needs both a linear
touching and a tangent algebra that rotates, i.e. generators with

    {S_a, S_b} = 2 delta_ab        and        [S_a, S_b] != 0,

and the claim is that the automaton cannot deliver both at once.

V214 and V215 moved the shift module to the icosahedral star, where the point
group has no quadratic anisotropy at all.  This experiment measures how far the
exclusivity actually extends, by scanning the tangent algebra across the whole
first zone rather than only at the gapless points.  At every momentum it
measures two numbers: the gap (the smallest eigenphase of the one-step matrix)
and the non-abelianity of the tangent algebra on the gap subspace -- every band
sitting at the gap energy, whose multiplicity is recorded,

    kappa(q) = max_ab ||[V_a, V_b]|| / max_a ||V_a||^2,   V_a = i P dM/dq_a P

with P the projector onto that subspace.  Restricting P to the two lowest bands
would choose an arbitrary pair where several bands touch at once, and the
algebra would then be an artefact of the choice; the full subspace makes kappa
basis-covariant, so it is a property of the step and not of the frame.  kappa is
dimensionless, scale-free, and zero exactly when the tangent algebra is abelian.
The boundary of the exclusivity is then a single measured number per module: the
power law

    kappa ~ gap^alpha

across the zone.  alpha = 0 would mean the algebra stays non-Abelian as the gap
closes - coexistence, a cone available in principle.  alpha > 0 means the
non-abelianity is forced to vanish with the gap: mutual exclusivity, with its
strength measured.  Both the committed cubic module and the star module are run
through the same instrument.

What the measurement gives: the exclusivity holds exactly where V187 proved it,
at the two-band crossings (the cubic M node, kappa 1.3e-11, velocity 0.603), and
it fails everywhere several bands touch at once.  The eight corners of the cubic
zone are exactly gapless with kappa 2.8 to 8.9 and velocities 0.7 to 2.9, and so
is every point of the icosahedral module's nodal lines along its two-fold axes
(kappa 8.5 to 23, velocities 1.0 to 1.7, gap below 1e-16).  Gaplessness and a
non-Abelian tangent algebra coexist on the new module: the power law is a
property of the gapped region, and it does not extrapolate onto the gapless set.

Inputs: cl10/config.yaml (the clock and the tolerances), the committed internal
map `supercell.sectors()[2]`, the committed cubic cell and stagger from
`cl10/altcell.py`, the star module from `cl10/continuum.py` and
`cl10/starcell.py`, and the scan box from the measured Lieb-Robinson radius.
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
GRID = 3                                  # points per axis of the zone scan
NBANDS = 2                                # the two lowest bands carry the crossing
DERIV_EPS = 1e-3                          # step for the tangent derivatives
FLOOR = 1e-9                              # numerical floor for a vanishing gap
VEL_FLOOR = 1e-6                          # below this the tangent algebra is zero


def as_side(x, ref):
    """Put a matrix on the same side (device or host) as the matrices it will
    be multiplied with; the engine never pulls a device array to the host."""
    if isinstance(ref, np.ndarray):
        return x if isinstance(x, np.ndarray) else E.asnumpy(x)
    return x if not isinstance(x, np.ndarray) else E.to_device(x)


def tangent(M_of_q, q, tol):
    """The gap and the tangent algebra of the low-energy subspace of M(q): the
    matrices V_a = i P dM/dq_a P on that subspace, their commutators and
    anticommutators (normalised by the square of the largest velocity), the gap
    itself, and the multiplicity of the gap that set the subspace."""
    q = np.asarray(q, float)
    M0 = M_of_q(q)
    w, V = E.eig(M0)
    ph = np.abs(np.angle(E.asnumpy(w)))
    # The tangent algebra has to live on a subspace that is physically defined.
    # Where the gap is open the two lowest bands are an isolated pair; where the
    # point is gapless the subspace is the whole gap subspace - all modes at the
    # numerical floor - because a two-mode cut through a larger zero space would
    # be an arbitrary basis choice.
    floor = max(10 * FLOOR, 100 * tol)
    zero = np.where(ph < floor)[0]
    order = zero if zero.size >= NBANDS else np.argsort(ph)[:NBANDS]
    Vr = E.asnumpy(V)[:, order]
    P = Vr @ Vr.conj().T
    Pd = as_side(P, M0)
    Vs = []
    for a in range(3):
        e = np.zeros(3)
        e[a] = DERIV_EPS
        D = M_of_q(q + e) - M_of_q(q - e)
        Vs.append(1j * ((Pd @ D @ Pd) / (2 * DERIV_EPS)))
    scale = float(max(E.absmax(v) for v in Vs))
    comm, anticomm = [], []
    for a, b in itertools.combinations(range(3), 2):
        ab, ba = Vs[a] @ Vs[b], Vs[b] @ Vs[a]
        comm.append(float(E.absmax(ab - ba)))
        anticomm.append(float(E.absmax(ab + ba)))
    squares = [float(E.absmax(v @ v)) for v in Vs]
    norm = max(scale, 1e-300) ** 2
    return {"q": [float(x) for x in q],
            "gap": float(ph.min()),
            "band_gaps": [float(x) for x in ph[order]],
            "multiplicity": int(len(order)),
            "gap_multiplicity": int(zero.size),
            "subspace": ("gap" if zero.size >= NBANDS else "two lowest bands"),
            "velocity": scale,
            "flat_touching": scale < VEL_FLOOR,
            "commutator": max(comm) if comm else 0.0,
            "anticommutator": max(anticomm) if anticomm else 0.0,
            "kappa": (max(comm) / norm if comm else 0.0)
            if scale >= VEL_FLOOR else None,
            "squares": [s / norm for s in squares]
            if scale >= VEL_FLOOR else None}


def power_ladder(xs, ys):
    """The log-log slope of the non-abelianity against the gap, with its
    correlation and the number of usable momenta (a gap at the numerical floor
    or a vanishing velocity carries no exponent)."""
    good = [(x, y) for x, y in zip(xs, ys) if x > FLOOR and y > 0]
    if len(good) < 3:
        return None, None, len(good)
    x = np.array([math.log(a) for a, _ in good])
    y = np.array([math.log(b) for _, b in good])
    slope, _ = np.polyfit(x, y, 1)
    return float(slope), float(np.corrcoef(x, y)[0, 1]), len(good)


def main():
    t0 = time.time()
    tol = float(CFG.tolerances.standard)
    g10 = float(CFG.angle.g10)
    branch = str(CFG.angle.branch)

    axes = K.module_generators()["vectors"]
    lr = X.lr_radius(axes)
    cell = X.module_cell()
    star_stages = X.module_stages(order="block-major")
    B, U32, U16 = SC.sectors()

    def star_M(q):
        return X.bloch_star(cell, U16, g10, star_stages, np.asarray(q, float),
                            axes)

    cubic_cell = A.cell((2, 2, 2))
    cubic_stages = A.stages_baseline()

    def cubic_M(q):
        return A.bloch_cell(cubic_cell, U16, g10, cubic_stages,
                            np.asarray(q, float))

    # ------------------------------------------------- the scan of both zones
    # The star module is almost periodic in q (six incommensurate Bloch phases),
    # so its zone is taken as the ball whose radius puts the longest generator's
    # phase at pi, sampled on a grid and along the star's own axes.  The cubic
    # module is scanned over its first Brillouin zone, with its measured gapless
    # line and its isotropy-forced endpoint included by name.
    scans = {}
    for name, M_of_q, rmax in (
            ("cubic", cubic_M, math.pi),
            ("star", star_M, math.pi / lr)):
        pts = []
        for i, j, k in itertools.product(range(GRID), repeat=3):
            q = (np.array((i, j, k), float) - (GRID - 1) / 2.0) * rmax
            pts.append(tangent(M_of_q, q, tol))
        named = {"GAMMA": np.zeros(3)}
        if name == "cubic":
            for label in ("M", "X", "R"):
                named[label] = np.array(SC.HIGH_SYM[label][0], float)
            mid = np.array(SC.HIGH_SYM["M"][0], float) / 2
            named["GAMMA-M quarter"] = np.array(SC.HIGH_SYM["M"][0], float) / 4
            named["GAMMA-M mid"] = mid
        else:
            for v in K.icosahedral_star()["directions"][:6]:
                vv = np.asarray(v, float)
                named["axis " + str([round(float(x), 3) for x in vv])] = (
                    rmax * vv / np.linalg.norm(vv))
        for label, q in named.items():
            row = tangent(M_of_q, q, tol)
            row["label"] = label
            pts.append(row)
        # kappa is a ratio over the square of the velocity: it is only defined
        # where the velocity does not vanish.  A gapless point with vanishing
        # velocity is a flat touching - the tangent algebra is zero there, which
        # is the abelian answer, not a division blow-up.
        live = [p for p in pts if not p["flat_touching"]]
        alpha, rcorr, nfit = power_ladder([p["gap"] for p in live],
                                          [p["kappa"] for p in live])
        closed = [p for p in pts if p["gap"] < FLOOR]
        closed_live = [p for p in closed if not p["flat_touching"]]
        scans[name] = {
            "rmax": rmax, "points": pts,
            "vel_floor": VEL_FLOOR,
            "flat_points": sum(1 for p in pts if p["flat_touching"]),
            "max_kappa": max([p["kappa"] for p in live], default=None),
            "gapless_points": closed,
            "gapless_flat": sum(1 for p in closed if p["flat_touching"]),
            "gapless_with_velocity": closed_live,
            "max_kappa_where_gapless": max([p["kappa"] for p in closed_live],
                                           default=None),
            "kappa_at_smallest_gap": min(pts, key=lambda p: p["gap"])["kappa"],
            "kappa_gap_exponent": alpha, "fit_r": rcorr, "fit_points": nfit,
        }

    # ------------------------------------------------------- the measurements
    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    for name in ("cubic", "star"):
        s = scans[name]
        chk(f"T1 the tangent algebra of the {name} module is non-Abelian "
            "somewhere in its zone",
            s["max_kappa"] is not None and s["max_kappa"] > 1e-6,
            f"scanned {len(s['points'])} momenta out to |q| = {s['rmax']:.4f}: "
            f"{len(s['points']) - s['flat_points']} carry a nonzero tangent "
            f"velocity, largest kappa {s['max_kappa']:.3e} (commutator over the "
            f"square of the largest velocity)")

    for name in ("cubic", "star"):
        s = scans[name]
        chk(f"T2 off the gapless set of the {name} module the non-abelianity "
            "obeys a measured power law kappa ~ gap^alpha with alpha > 0",
            s["kappa_gap_exponent"] is not None
            and s["kappa_gap_exponent"] > 0.0,
            f"alpha = {s['kappa_gap_exponent']:.3f} over {s['fit_points']} "
            f"momenta with a nonzero velocity and an open gap (correlation "
            f"{s['fit_r']:.3f}); the same quantity on the gapless set is not "
            f"small at all (largest kappa {s['max_kappa_where_gapless']}), so "
            "the power law describes the gapped region and does not extrapolate "
            "onto the gapless set (T3)")

    for name in ("cubic", "star"):
        s = scans[name]
        live = s["gapless_with_velocity"]
        ab = [p for p in live if p["kappa"] < 1e-3]
        nonab = [p for p in live if p["kappa"] >= 1e-3]
        kappa_max = max((p["kappa"] for p in nonab), default=0.0)
        chk(f"T3 gaplessness and a non-Abelian tangent algebra coexist on the "
            f"{name} module: the gapless set splits into the abelian points "
            "where V187's exclusivity holds and points whose algebra does not "
            "commute",
            bool(nonab),
            f"{len(s['gapless_points'])} scanned points are gapless below "
            f"{FLOOR:.0e}: {s['gapless_flat']} are flat touchings (vanishing "
            f"velocity, so there is no algebra to test); of the {len(live)} "
            f"with a nonzero tangent velocity, {len(ab)} have an abelian algebra "
            f"(kappa below 1e-3) and {len(nonab)} a non-Abelian one (largest "
            f"kappa {kappa_max:.3f}): "
            + "; ".join(f"{p.get('label', [round(x, 3) for x in p['q']])} gap "
                        f"{p['gap']:.1e} multiplicity {p['gap_multiplicity']} "
                        f"|V| {p['velocity']:.2e} kappa {p['kappa']:.3f}"
                        for p in live[:4]))

    live_star = [p for p in scans["star"]["points"]
                 if not p["flat_touching"] and p["kappa"] is not None]
    best_star = max(live_star, key=lambda p: p["kappa"])
    spread = max(best_star["squares"]) - min(best_star["squares"])
    chk("T4 the Clifford relation fails at the most non-abelian point of the "
        "star module: the squares of its tangent generators are not "
        "proportional and its anticommutators do not vanish",
        spread > 1e-3 or best_star["anticommutator"] > 1e-3,
        f"star, most non-abelian point with a nonzero velocity "
        f"{best_star.get('label', [round(x, 4) for x in best_star['q']])}: gap "
        f"{best_star['gap']:.2e}, |V| {best_star['velocity']:.2e}, kappa "
        f"{best_star['kappa']:.3e}, largest anticommutator "
        f"{best_star['anticommutator']:.3e}, normalised V_a^2 "
        f"{[round(x, 3) for x in best_star['squares']]} (spread {spread:.3f})")

    star_alpha = scans["star"]["kappa_gap_exponent"]
    cubic_alpha = scans["cubic"]["kappa_gap_exponent"]

    cub_live = scans["cubic"]["gapless_with_velocity"]
    cub_ab = sorted(p["kappa"] for p in cub_live if p["kappa"] < 1e-3)
    cub_nonab = sorted(p["kappa"] for p in cub_live if p["kappa"] >= 1e-3)
    star_live = scans["star"]["gapless_with_velocity"]
    star_nonab = sorted(p["kappa"] for p in star_live if p["kappa"] >= 1e-3)
    star_m = sorted({p["gap_multiplicity"] for p in star_live})
    cub_m = sorted({p["gap_multiplicity"] for p in cub_live})
    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the boundary of the V187 mutual exclusivity: the gap and "
                    "the non-abelianity kappa of the tangent algebra on the gap "
                    "subspace (with its multiplicity) measured across the first "
                    "zone of both the committed cubic module and the icosahedral "
                    "star module, with the power law kappa ~ gap^alpha off the "
                    "gapless set and the largest kappa on the gapless set",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"EXCLUSIVITY BOUNDARY at the {branch} clock g = {g10:.6f}: the gap "
            f"and the non-abelianity kappa of the tangent algebra on the gap "
            f"subspace, measured at "
            + " and ".join(f"{len(scans[n]['points'])} momenta of the {n} "
                           f"module (out to |q| = {scans[n]['rmax']:.4f})"
                           for n in ("cubic", "star"))
            + f". kappa = max_ab ||[V_a,V_b]|| / max_a ||V_a||^2 is only "
              f"meaningful where the velocity survives the floor "
              f"{VEL_FLOOR:.0e}: {scans['cubic']['flat_points']} of the cubic "
              f"and {scans['star']['flat_points']} of the star momenta have "
              f"vanishing tangent velocity. (T1) the algebra is non-Abelian in "
              f"both zones: largest kappa {scans['cubic']['max_kappa']:.3e} "
              f"(cubic) and {scans['star']['max_kappa']:.3e} (star). (T2) away "
              f"from the gapless set the boundary is a power law kappa ~ "
              f"gap^alpha with alpha = {cubic_alpha:.3f} (cubic, correlation "
              f"{scans['cubic']['fit_r']:.3f}, {scans['cubic']['fit_points']} "
              f"momenta) and alpha = {star_alpha:.3f} (star, correlation "
              f"{scans['star']['fit_r']:.3f}, {scans['star']['fit_points']} "
              f"momenta), a law of the gapped region only. (T3) it does not "
              f"reach the gapless set: of "
              f"{len(scans['cubic']['gapless_points'])} gapless cubic momenta "
              f"{scans['cubic']['gapless_flat']} are flat touchings (zero "
              f"velocity, no algebra to test); of the {len(cub_live)} with a "
              f"nonzero velocity, {len(cub_ab)} have an abelian algebra (largest "
              f"kappa {cub_ab[-1] if cub_ab else 0:.1e}, the M node and the V187 "
              f"domain, gap-subspace multiplicity {cub_m}) and {len(cub_nonab)} "
              f"a non-Abelian one (zone corners, kappa up to "
              f"{cub_nonab[-1] if cub_nonab else 0:.3f}); on the icosahedral "
              f"module all {len(star_live)} gapless momenta with a velocity are "
              f"non-Abelian, kappa up to "
              f"{star_nonab[-1] if star_nonab else 0:.3f} at a gap of "
              f"{min((p['gap'] for p in star_live), default=0.0):.1e} "
              f"with gap-subspace multiplicity {star_m}: gaplessness and a "
              f"rotating tangent algebra coexist on the new module, which is the "
              f"half of the V187 exclusivity it escapes. (T4) at the most "
              f"non-abelian star momentum with a nonzero velocity the gap is "
              f"{best_star['gap']:.2e} and the normalised squares of the "
              f"tangent generators are "
              f"{[round(x, 3) for x in best_star['squares']]} (spread "
              f"{spread:.3f}) with largest anticommutator "
              f"{best_star['anticommutator']:.3e}: the Clifford relation fails "
              f"there too."),
        "data": {"clock_branch": branch, "g10": g10, "lr_radius": lr,
                 "grid": GRID, "n_bands": NBANDS, "deriv_eps": DERIV_EPS,
                 "floor": FLOOR, "scans": scans,
                 "engine": E.backend_info(), "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v216_exclusivity_boundary.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V216 - the boundary of the gaplessness / non-Abelian exclusivity")
    print("=" * 68)
    for name in ("cubic", "star"):
        s = scans[name]
        print(f"  {name}: {len(s['points'])} points to |q| = {s['rmax']:.4f}, "
              f"{s['flat_points']} flat, max kappa {s['max_kappa']}, alpha "
              f"{s['kappa_gap_exponent']}, gapless "
              f"{len(s['gapless_points'])} ({s['gapless_flat']} flat), max "
              f"kappa among gapless with velocity "
              f"{s['max_kappa_where_gapless']}")
        for p in s["points"]:
            if "label" in p:
                kap = ("flat" if p["flat_touching"]
                       else f"{p['kappa']:.2e}")
                sq = ("-" if p["squares"] is None
                      else str([round(x, 2) for x in p["squares"]]))
                print(f"     {p['label']:34s} gap {p['gap']:.2e} |V| "
                      f"{p['velocity']:.2e} kappa {kap} V^2 {sq}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
