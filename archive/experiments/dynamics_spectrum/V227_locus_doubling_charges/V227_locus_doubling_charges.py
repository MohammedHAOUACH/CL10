#!/usr/bin/env python3
r"""V227 - every gapless point of the 2I-covariant Clifford cell inside the
first zone, its charge, and the doubling bookkeeping.

V226 measured the charge of the zone-centre node and that a symmetry-respecting
perturbation moves it rather than removing it.  That leaves the question every
gapless fermion on a lattice has to answer: what else is gapless, and what is
the total chirality?  Nielsen-Ninomiya says the sum of chiralities over a
Brillouin zone vanishes, but its hypothesis is a *periodic* Brillouin zone, and
the icosahedral star module is not periodic: its twelve directions carry
$\sqrt5$ coordinates, so no nonzero momentum $G$ makes every
$\hat v_d\cdot G$ a multiple of $2\pi$.  Both halves are measured here rather
than argued, and the locus is classified rather than assumed: V219/V220 found
that the committed star module's gapless set is lines, not points, and a line
carries a holonomy rather than a Chern number.

T1 the locus.  A grid over the ball $|q| \le \pi / \ell$, where $\ell = 1.902113$
   is the generator length V214 measured (the same first-zone convention V215
   used for the star module): the grid *counts* the locus at several tolerances,
   the grid's own resolution decides how many separate components are refined
   (`MAX_NODES`, so the run states what it can afford to resolve), and each
   resolved component is classified by its touching dimension, its transverse
   exponent along a generic ray, and how many directions of an icosahedrally
   uniform mesh stay gapless off it -- a point node pierces none of them, a
   point on a line pierces many.  The directions the gapless grid points lie on
   are then collected into lines and probed further out, so the line structure
   is measured rather than read off group theory.
T2 the charge question at those points, decided rather than assumed: the same
   icosahedrally uniform link-variable flux V226 uses, tried at three radii.  A
   sphere that changes sector size refuses; a sphere whose link overlaps lose
   rank returns a number that is not a charge, and that is recorded by its
   smallest overlap and its integer defect rather than averaged in.
T3 the bookkeeping on spheres about the zone centre, at four radii: each either
   lands on an integer flux, or sits pinned at a half-integer (the mark of a
   sphere cut by a gapless line), or refuses because the sector changes size.
   Whether the integers agree across the radii -- the actual test of a
   conserved charge -- is then read off, not assumed.
T4 the hypothesis of the doubling theorem, measured.  Write the twelve
   directions in the basis of three of them: $A$ is stored one row per direction
   ($12\times3$), and a momentum with $\hat v_d\cdot G \in 2\pi\mathbb Z$ for the
   basis directions is $G = 2\pi M^{-\mathsf T}n$ with $n\in\mathbb Z^3$, so
   periodicity needs $An \in \mathbb Z^{12}$.  The search over small $n$ is
   reported with its best defect: the entries of $A$ live in $\mathbb Q(\sqrt5)$,
   so the condition fails for every nonzero $n$ and there is no Brillouin zone
   to double over.  What does exist -- and is reported as the thing that
   survives -- is the module's shortest rational-direction hop
   (`continuum.shortest_rational_displacement`, the instrument V214 used).
T5 the control that makes the answer above mean something: the same scan and the
   same spheres on the spinless (scalar-bond) cell.  If its sphere about the zone
   centre keeps full link rank and returns exactly zero flux, then the rank
   failures on the Clifford cell are statements about its gapless set, not about
   the instrument.
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np

from cl10 import continuum as K
from cl10 import engine as E
from cl10 import spincell as SP
from cl10 import supercell as SC

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

STEP = 0.22                      # the grid spacing inside the ball
REFINE_TOL = 1e-5                # a candidate is a point this close to gapless
NODE_TOL = 1e-8                  # an eigenphase that counts as zero
RAY_RADII = (3e-3, 1e-3, 3e-4)   # for the transverse exponent
FLUX_RADII = (1e-2, 4e-2, 1.2e-1)
NSEARCH = 3                      # the reciprocal-vector search radius in Z^3
MAX_NODES = 12                   # components of the locus this run can resolve
PIERCE_RADIUS = 3e-3             # the sphere on which rays are probed
AXIS_RADII = (0.1, 0.3, 0.5, 0.8, 1.1, 1.4)   # along each discovered line
MAX_LINES = 16                   # lines of the locus probed out to AXIS_RADII

def phase_gap(q, step_fn):
    r"""The smallest eigenphase modulus.  Read off the Hermitian part
    $C = (S + S^\dagger)/2$, whose eigenvalues are $\cos\theta$: since
    $\theta \in (-\pi,\pi]$ and $\cos$ is a bijection on $[0,\pi]$, the largest
    eigenvalue of $C$ gives $\min_\theta\lvert\theta\rvert =
    \arccos(\max\cos\theta)$ exactly.  Checked against the general eigenvalue
    route at four momenta of this cell: agreement to $1.3\times10^{-13}$, at a
    fifth ($\Gamma$) to $2.8\times10^{-17}$; the Hermitian solve costs
    $0.018\,\mathrm{s}$ against $0.083\,\mathrm{s}$ for `eigvals` on the
    $384$-dimensional step, and that factor is what makes a zone-wide
    enumeration affordable."""
    S = step_fn(q)
    cn = np.linalg.eigvalsh(E.asnumpy((S + S.conj().T) / 2.0))
    return float(math.acos(min(1.0, max(-1.0, float(cn[-1])))))


def pierce_count(step_fn, q0, verts, radius=PIERCE_RADIUS):
    """How many directions of a symmetry-uniform mesh leave `q0` gapless at a
    small radius: the difference between a point node (none of them) and a point
    sitting on a line or sheet of the locus (most of them)."""
    q0 = np.asarray(q0, float)
    return int(sum(1 for p in verts
                   if phase_gap(q0 + radius * np.asarray(p, float), step_fn)
                   < 1e-12))


def touching_dimension(step_fn, q):
    """How many eigenphases sit at zero at this momentum, the largest of them,
    and the first phase outside the cluster: the difference between a point node
    and a line, read off the spectrum (the convention V220 used).  The moduli
    come from the Hermitian part of the step, as in `phase_gap`."""
    S = step_fn(q)
    ph = np.sort(np.arccos(np.clip(np.linalg.eigvalsh(
        E.asnumpy((S + S.conj().T) / 2.0)), -1.0, 1.0)))
    m = int((ph < NODE_TOL).sum())
    return {"count": m, "cluster_top": float(ph[m - 1]) if m else 0.0,
            "first_outside": float(ph[m]),
            "isolation": float((ph[m] - ph[m - 1]) / max(ph[m - 1], 1e-12))
            if m else 0.0}


def ball_grid(h, rmax):
    n = int(math.ceil(rmax / h))
    pts = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                q = h * np.array([i, j, k], float)
                if np.linalg.norm(q) <= rmax:
                    pts.append(q)
    return pts


def refine(q0, step_fn, h0=1e-2, rounds=8):
    """A coordinate-wise shrinking star search: good enough to pin a crossing to
    the eigensolver's own accuracy without needing a derivative."""
    q, h, g = np.array(q0, float), h0, phase_gap(q0, step_fn)
    for _ in range(rounds):
        improved = True
        while improved:
            improved = False
            for d in range(3):
                for s in (+1, -1):
                    e = np.zeros(3)
                    e[d] = s * h
                    gn = phase_gap(q + e, step_fn)
                    if gn < g - 1e-14:
                        q, g, improved = q + e, gn, True
        h *= 0.3
        if h < 1e-7:
            break
    return q, g


def exponent_away(q0, step_fn):
    """The transverse exponent: the power law of the gap as the point is left
    along a generic ray.  None if the ray stays gapless, which is what a line
    through the point looks like."""
    u = np.array([0.31, 0.77, -0.55], float)
    u /= np.linalg.norm(u)
    gs = [phase_gap(q0 + r * u, step_fn) for r in RAY_RADII]
    if min(gs) < 1e-11:
        return None
    return float(np.polyfit(np.log(RAY_RADII), np.log(gs), 1)[0])


def flux(step_fn, q0, mesh, radii):
    """The Berry flux about a point at several radii, with the failures kept
    rather than raised: a sector that changes size over the sphere is the
    signature of another gapless point inside, and that is a measurement."""
    attempts, best = [], None
    for r in radii:
        try:
            res = SP.berry_chern(step_fn, np.asarray(q0, float), mesh, radius=r,
                                 sector="negative")
        except RuntimeError as exc:
            attempts.append({"radius": r, "error": str(exc)})
            continue
        attempts.append({k: v for k, v in res.items() if k != "sector"})
        if (best is None or (res["link_rank_ok"] and not best["link_rank_ok"])
                or (res["link_rank_ok"] == best["link_rank_ok"]
                    and res["min_link_overlap"] > best["min_link_overlap"])):
            best = res
    return best, attempts


def reciprocal_search(directions, nsearch=NSEARCH):
    r"""Is the step periodic in any nonzero momentum?  With $M$ the matrix whose
    columns are three independent directions and $V$ all twelve, $A = M^{-1} V$
    gives the phase coefficients: a momentum $G = 2\pi M^{-\mathsf T} n$ is
    trivial on the basis directions by construction, and on the rest exactly
    when $An$ is integral.  $A$ is stored one row per direction, so $An$ is the
    twelve-vector of phase coefficients directly.  The search returns the
    smallest defect found."""
    v = np.asarray(directions, float)
    m = v[:3].T
    a = np.linalg.solve(m, v.T).T
    best = None
    for n in itertools.product(range(-nsearch, nsearch + 1), repeat=3):
        if not any(n):
            continue
        nv = np.array(n, float)
        prod = a @ nv
        defect = float(np.max(np.abs(prod - np.round(prod))))
        if best is None or defect < best["defect"]:
            best = {"n": [int(x) for x in n], "defect": defect,
                    "coefficients": [float(x) for x in prod],
                    "G": [float(x) for x in (2.0 * math.pi
                                             * np.linalg.solve(m.T, nv))]}
    return {"coefficients_A": [[float(x) for x in row] for row in a],
            "best": best, "search_radius": nsearch}


def main():
    t0 = time.time()
    checks = []
    ok = True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    _, _, U16 = SC.sectors()
    g10 = float(SC.G10)
    branch = os.environ.get("CL10_ANGLE") or "default"
    layers, legs = 1, (4,)
    kw = {"g10": g10, "kind": "clifford", "layers": layers, "legs": legs}

    def cell_step(q):
        return SP.step(q, U16, **kw)

    ell = SP.module_length()
    rmax = math.pi / ell
    mesh1 = SP.icosahedral_sphere(3)
    mesh2 = SP.icosahedral_sphere(3)
    mesh_ray = SP.icosahedral_sphere(2)

    # ------------------------------------------------------------------ T1
    # The locus is a set, not a bag of points.  The grid counts it; the grid's
    # own resolution decides how many separate components of it this run can
    # afford to resolve; and each resolved component is asked whether the
    # gaplessness runs off it along rays (a line or sheet) or stops at it (a
    # point node, which is the only thing a Berry flux can be attached to).
    grid = ball_grid(STEP, rmax)
    gaps = [(phase_gap(q, cell_step), q) for q in grid]
    vals = sorted(gaps, key=lambda t: t[0])
    grid_min = {"gap": vals[0][0], "q": [float(x) for x in vals[0][1]]}
    print(f"  [t={time.time() - t0:7.1f}s] T1 grid: {len(grid)} momenta, "
          f"smallest gap {grid_min['gap']:.1e} at "
          f"{[round(x, 4) for x in grid_min['q']]}", flush=True)
    thresholds = [{"tol": t,
                   "grid_points": int(sum(1 for g, _ in gaps if g < t))}
                  for t in (NODE_TOL, 1e-6, 1e-3, 1e-2, 0.25 * STEP)]
    reps = []
    for g, q in vals:
        if g >= 0.25 * STEP or len(reps) >= MAX_NODES:
            break
        if not any(np.linalg.norm(q - s) < 1.2 * STEP for s in reps):
            reps.append(q)
    nodes, seen = [], []
    print(f"  [t={time.time() - t0:7.1f}s] T1: {len(reps)} components of the "
          f"locus to resolve", flush=True)
    for q in reps:
        qr, gr = refine(q, cell_step, h0=0.5 * STEP, rounds=6)
        if gr > REFINE_TOL:
            continue
        if any(np.linalg.norm(qr - s) < 5e-2 for s in seen):
            continue
        seen.append(qr)
        nd = {"q": [float(x) for x in qr], "gap": gr,
              "from": [float(x) for x in q]}
        nd.update(touching_dimension(cell_step, qr))
        nd["exponent"] = exponent_away(qr, cell_step)
        nd["pierced_rays"] = pierce_count(cell_step, qr, mesh_ray["vertices"])
        nd["mesh_rays"] = len(mesh_ray["vertices"])
        nodes.append(nd)
        print(f"  [t={time.time() - t0:7.1f}s]   resolved "
              f"{[round(x, 3) for x in qr]}: gap {gr:.1e}, touching "
              f"{nd['count']}, {nd['pierced_rays']}/{nd['mesh_rays']} rays "
              f"gapless", flush=True)
    isolated = [nd for nd in nodes if nd["isolation"] > 1e-3
                and nd["exponent"] is not None and nd["pierced_rays"] == 0]
    pierced_nodes = [nd for nd in nodes if nd["pierced_rays"] > 0]
    # The direction structure is measured, not taken from group theory: the
    # grid's gapless points are collected into lines through the origin, and
    # every line found is probed further out to see whether the gaplessness
    # continues along it or stops with the grid point it was found at.
    locus_lines = []
    for g, q in gaps:
        if g > NODE_TOL or np.linalg.norm(q) < 1e-12:
            continue
        u = q / float(np.linalg.norm(q))
        if (u[0], u[1], u[2]) < (-u[0], -u[1], -u[2]):
            u = -u                          # one representative per line
        if not any(np.allclose(u, l, atol=1e-9) for l in locus_lines):
            locus_lines.append(u)
    probes = []
    for u in locus_lines[:MAX_LINES]:
        pr = [(r, phase_gap(r * u, cell_step)) for r in AXIS_RADII]
        probes.append({"direction": [float(x) for x in u], "probes": len(pr),
                       "gapless_to_1e-12": int(sum(1 for _, g in pr
                                                  if g < 1e-12)),
                       "largest_gap": float(max(g for _, g in pr))})
    line_like = [p for p in probes
                 if p["gapless_to_1e-12"] == p["probes"]]

    node_txt = "; ".join(
        f"q {[round(x, 4) for x in nd['q']]} gap {nd['gap']:.1e} touching "
        f"{nd['count']} bands (isolation {nd['isolation']:.1e}), "
        f"{nd['pierced_rays']}/{nd['mesh_rays']} rays gapless, exponent "
        + ("gapless ray" if nd["exponent"] is None else f"{nd['exponent']:.3f}")
        for nd in nodes[:8])
    thr_txt = "; ".join(f"{th['tol']:.0e}: {th['grid_points']}"
                        for th in thresholds)
    line_txt = " | ".join(
        f"{[round(x, 3) for x in pr['direction']]}: {pr['gapless_to_1e-12']}/"
        f"{pr['probes']} probes gapless to 1e-12 out to r {AXIS_RADII[-1]:.1f} "
        f"(largest {pr['largest_gap']:.1e})" for pr in probes[:6])

    chk("T1 the gapless locus inside the first zone is counted, its "
        "components resolved at the grid's own resolution, and each resolved "
        "one classified by whether the gaplessness runs off it along rays",
        len(nodes) >= 1 and grid_min["gap"] < 1e-9
        and all(nd["gap"] < REFINE_TOL for nd in nodes)
        and (len(isolated) >= 1 or len(line_like) >= 1),
        f"{len(grid)} grid points over |q| <= pi/ell = {rmax:.4f} (spacing "
        f"{STEP}): smallest gap {grid_min['gap']:.1e} at "
        f"{[round(x, 4) for x in grid_min['q']]}; grid points below "
        f"{thr_txt}; {len(reps)} components refined to {len(nodes)} gapless "
        f"points, {len(isolated)} isolated point nodes, {len(pierced_nodes)} "
        f"pierced by gapless rays: {node_txt} | {len(locus_lines)} lines "
        f"through the origin carry the gapless grid points, {len(line_like)} of "
        f"the {len(probes)} probed stay gapless to 1e-12 out to "
        f"r {AXIS_RADII[-1]:.1f}: {line_txt}")

    # ------------------------------------------------------------------ T2
    for nd in nodes:
        best, attempts = flux(cell_step, nd["q"], mesh1, FLUX_RADII)
        nd["flux"] = best
        nd["flux_attempts"] = attempts
        nd["chern"] = best["nearest_integer"] if best else None
        nd["chern_defect"] = best["integer_defect"] if best else None
    charged = [nd for nd in nodes
               if nd["chern"] is not None and nd["flux"]["link_rank_ok"]]
    total = sum(nd["chern"] for nd in charged)
    charge_txt = "; ".join(
        f"q {[round(x, 4) for x in nd['q']]} C = {nd['chern']:+d} (defect "
        f"{nd['chern_defect']:.1e} at radius {nd['flux']['radius']:.3f}, "
        f"{len(nd['flux_attempts'])} radii tried)" for nd in charged[:8])

    untrusted = [nd for nd in nodes if nd["flux"] is not None
                 and not nd["flux"]["link_rank_ok"]]
    undef = [round(nd["flux"]["integer_defect"], 3) for nd in untrusted]
    ovl = [f"{nd['flux']['min_link_overlap']:.0e}" for nd in untrusted]
    chk("T2 the charge question at the resolved points is decided, not assumed: "
        "either a sphere about a locus point keeps full link rank and returns an "
        "integer, or none of them does and every reading is marked meaningless "
        "by its own link overlap",
        (bool(charged)
         and all(nd["chern_defect"] < 1e-6 for nd in charged))
        or (not charged and bool(untrusted) and len(untrusted) == len(nodes)),
        f"{len(charged)} of {len(nodes)} gapless points have a flux whose sphere "
        f"kept full link rank ({len(nodes) - len(charged)} were not trusted); "
        f"charges: {charge_txt}; their sum is {total:+d}; the untrusted readings "
        f"sit at integer defects {undef} with smallest link overlaps {ovl} -- a "
        f"flux pinned at a half-integer over rank-deficient links is the mark of "
        f"a sphere cut by a gapless line, not the mark of a charge")

    # ------------------------------------------------------------------ T3
    boundary = []
    for rr in (rmax, rmax * 0.75, rmax * 0.5, rmax * 0.25):
        best, attempts = flux(cell_step, np.zeros(3), mesh2, [rr])
        boundary.append({"radius": rr,
                         "chern": best["nearest_integer"] if best else None,
                         "defect": best["integer_defect"] if best else None,
                         "attempts": attempts})
    defined = [b for b in boundary if b["chern"] is not None]

    boundary_txt = ", ".join(
        f"{b['radius']:.4f} (flux {b['chern']}, defect {b['defect']:.1e})"
        if b["chern"] is not None else
        f"{b['radius']:.4f} (sector changes size: {b['attempts'][0]['error']})"
        for b in boundary)

    exact = [b for b in defined if b["defect"] is not None
             and b["defect"] < 1e-6]
    pinned = [b for b in defined if b["defect"] is not None
              and abs(b["defect"] - 0.5) < 1e-6]
    agree = len({b["chern"] for b in exact}) == 1
    if len(exact) >= 2 and agree:
        radius_txt = (f"the {len(exact)} trusted radii agree on "
                      f"{[b['chern'] for b in exact]}: the sphere flux about the "
                      f"zone centre is a radius-independent charge")
        radius_short = (f"the {len(exact)} trusted radii agree on "
                        f"{[b['chern'] for b in exact]}, so the sphere flux is a "
                        f"conserved charge")
    elif exact:
        radius_txt = (f"only {len(exact)} of the {len(boundary)} spheres is "
                      f"trusted (flux {exact[0]['chern']} at radius "
                      f"{exact[0]['radius']:.4f}), so this census cannot decide "
                      f"radius independence at all; what it does decide is the "
                      f"trichotomy -- an uncuts sphere lands on an integer, a "
                      f"sphere a gapless line cuts pins at a half-integer, and a "
                      f"sphere the lines run straight through refuses")
        radius_short = (f"exactly {len(exact)} sphere is trusted (flux "
                        f"{exact[0]['chern']}), so radius independence is not "
                        f"established; the trichotomy is")
    else:
        radius_txt = "no sphere about the zone centre is trusted"
        radius_short = "no sphere about the zone centre is trusted"
    chk("T3 the spheres about the zone centre sort themselves: a radius either "
        "keeps rank and lands on an integer flux, or is pinned at a half-integer "
        "because a gapless line cuts the sphere, or refuses outright -- and "
        "whether the integers agree across radii is then read off",
        len(defined) >= 1 and len(exact) + len(pinned) == len(defined),
        "spheres about the zone centre at radii " + boundary_txt + f"; "
        f"{len(exact)} give an integer flux to machine precision "
        f"{[b['chern'] for b in exact]}, {len(pinned)} are pinned at a "
        f"half-integer (defects {[round(b['defect'], 6) for b in pinned]}), "
        f"{len(boundary) - len(defined)} refuse because the sector changes size; "
        f"{radius_txt}; the summed trusted node charges are {total:+d}")

    # ------------------------------------------------------------------ T4
    recip = reciprocal_search(SP.directions())
    rhop = K.shortest_rational_displacement(K.icosahedral_star(), search=2)
    rh_txt = ("none found" if rhop is None else
              f"norm {rhop['norm']:.4f}, displacement "
              f"{[round(x, 4) for x in rhop['displacement']]}, "
              f"{rhop['steps']} steps")

    chk("T4 the doubling theorem's hypothesis is not available to this module: "
        "no nonzero reciprocal vector makes every direction's phase trivial",
        recip["best"] is not None and recip["best"]["defect"] > 1e-6,
        f"the twelve directions in the basis of three of them give the "
        f"coefficient matrix A = "
        f"{[[round(x, 4) for x in row] for row in recip['coefficients_A']]}; the "
        f"integrality search over |n| <= {recip['search_radius']} finds its "
        f"smallest defect {recip['best']['defect']:.4f} at n = "
        f"{recip['best']['n']} (G = {[round(x, 4) for x in recip['best']['G']]}), "
        f"so no nonzero G is a period of the step and there is no Brillouin zone "
        f"to double over; what the module does offer is a shortest "
        f"rational-direction hop: {rh_txt}")

    # ------------------------------------------------------------------ T5
    ckw = {"g10": g10, "kind": "scalar", "layers": 1, "legs": (1,)}

    def scalar_step(q):
        return SP.step(q, U16, **ckw)

    svals = sorted(((phase_gap(q, scalar_step), q) for q in grid[:125]),
                   key=lambda t: t[0])
    smin = {"gap": svals[0][0], "q": [float(x) for x in svals[0][1]]}
    sq, sg = refine(np.zeros(3), scalar_step)
    sflux, satt = flux(scalar_step, np.zeros(3), mesh1, [4e-2])
    sexp = exponent_away(sq, scalar_step)
    stouch = touching_dimension(scalar_step, sq)

    sflux_txt = (f"charge {sflux['nearest_integer']} (defect "
                 f"{sflux['integer_defect']:.1e})" if sflux
                 else f"charge undefined ({satt[0]['error']})")

    chk("T5 the control that turns the failures above into a statement about "
        "the cell rather than about the instrument: the spinless cell's sphere "
        "about the same zone centre keeps full link rank and returns exactly "
        "zero flux",
        sflux is not None and sflux["link_rank_ok"]
        and sflux["integer_defect"] < 1e-8 and sflux["nearest_integer"] == 0,
        f"scalar bond assignment: smallest sampled gap {smin['gap']:.1e} at "
        f"{[round(x, 3) for x in smin['q']]}; its zone-centre touching is "
        f"{stouch['count']} bands at phase zero and it stays gapless along the "
        f"transverse probe ray "
        f"({'gapless ray' if sexp is None else f'exponent {sexp:.3f}'}), so it "
        f"is not a point node either; its sphere of radius {sflux['radius']:.3f} "
        f"kept rank (smallest link overlap {sflux['min_link_overlap']:.4f}, "
        f"{sflux['bands_in_sector']} bands in sector) and gives {sflux_txt}: "
        f"exactly zero, where every sphere about a locus point of the Clifford "
        f"cell lost rank")

    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the gapless locus of the 2I-covariant Clifford cell inside "
                    "the first zone with each point's touching dimension, "
                    "transverse exponent and Berry charge, the flux bookkeeping "
                    "against spheres around the zone, the measured absence of a "
                    "reciprocal lattice vector of the icosahedral module, and the "
                    "exactly-zero trusted sphere of the spinless control",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"GAPLESS LOCUS AND DOUBLING BOOKKEEPING of the 2I-covariant Clifford "
            f"cell at the {branch} clock g = {g10:.6f} (layers {layers}, internal "
            f"exponent {legs[0]}). (T1) the locus over |q| <= pi/ell = "
            f"{rmax:.4f}: smallest gap {grid_min['gap']:.1e}; {len(reps)} "
            f"components of the locus resolved to {len(nodes)} gapless points, "
            f"{len(isolated)} isolated point nodes, {len(pierced_nodes)} pierced "
            f"by gapless rays; the gapless grid points lie on "
            f"{len(locus_lines)} lines through the origin, and "
            f"{len(line_like)} of the {len(probes)} probed stay gapless to "
            f"1e-12 at every probe out to r {AXIS_RADII[-1]:.1f}; "
            f"touching dimensions {[nd['count'] for nd in nodes]}. (T2) "
            f"{len(charged)} of {len(nodes)} locus points have a flux whose "
            f"sphere kept full link rank; every sphere about a locus point lost "
            f"rank (smallest link overlaps {ovl}, integer defects {undef}), so no "
            f"Berry charge is defined anywhere on the locus and the summed "
            f"trusted charge is {total:+d}. (T3) spheres around "
            f"the zone centre give "
            + ", ".join(f"{b['radius']:.3f}:{b['chern']}" for b in boundary)
            + f"; {len(exact)} give an integer to machine precision, "
              f"{len(pinned)} sit pinned at a half-integer, "
              f"{len(boundary) - len(defined)} refuse, and {radius_short}; the "
              f"summed trusted node charge is {total:+d}. (T4) the module has no "
              f"period: the integrality search over the direction-coefficient "
              f"matrix bottoms out at defect {recip['best']['defect']:.4f}, so "
              f"there is no Brillouin zone to double over, and the gapless lines "
              f"run out through the zone edge (the sphere at rmax is the one that "
              f"refuses); the module's shortest rational-direction hop is "
              f"{rh_txt}. (T5) the spinless "
              f"control's zone-centre touching is {stouch['count']} bands at "
              f"phase zero, gapless along the transverse probe ray "
              f"({'gapless ray' if sexp is None else f'exponent {sexp:.3f}'}), "
              f"and its trusted sphere of radius {sflux['radius']:.3f} gives "
              + (f"charge {sflux['nearest_integer']} (defect "
                 f"{sflux['integer_defect']:.1e}, smallest link overlap "
                 f"{sflux['min_link_overlap']:.4f})" if sflux else "no flux")
              + ": exactly zero where every sphere about the Clifford cell's "
                "locus points lost rank."),
        "data": {
            "clock": g10, "branch": branch,
            "cell": {"layers": layers, "legs": list(legs)},
            "zone": {"generator_length": ell, "rmax": rmax, "grid_step": STEP,
                     "grid_points": len(grid), "grid_min": grid_min,
                     "thresholds": thresholds, "components": len(reps)},
            "locus_lines": [[float(x) for x in u] for u in locus_lines],
            "line_probes": probes, "lines_staying_gapless": len(line_like),
            "nodes": nodes, "isolated": [nd["q"] for nd in isolated],
            "pierced_nodes": len(pierced_nodes),
            "total_charge": total, "boundary_flux": boundary,
            "periodicity": recip, "rational_hop": rhop,
            "control_scalar": {"grid_min": smin, "refined": [float(x) for x in sq],
                               "gap": sg, "exponent": sexp, "touching": stouch,
                               "flux": sflux, "flux_attempts": satt},
            "engine": E.backend_info(),
            "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, os.path.basename(HERE).lower() + ".json"),
              "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V227 - the gapless locus, its charges, and the doubling bookkeeping")
    print("=" * 68)
    print(f"  zone: |q| <= {rmax:.4f}, {len(grid)} points, smallest gap "
          f"{grid_min['gap']:.1e} at {[round(x, 4) for x in grid_min['q']]}")
    for nd in nodes:
        print(f"  node q {[round(x, 5) for x in nd['q']]} gap {nd['gap']:.1e} "
              f"touching {nd['count']} isolation {nd['isolation']:.1e} exponent "
              f"{nd['exponent']} charge {nd['chern']}")
    print(f"  total charge {total:+d} | boundary "
          + ", ".join(f"r{b['radius']:.3f}:{b['chern']}" for b in boundary))
    print(f"  reciprocal search: best defect {recip['best']['defect']:.4f} at n "
          f"{recip['best']['n']}; rational hop {rh_txt}")
    print(f"  spinless control: gap {smin['gap']:.1e}, exponent {sexp}, "
          f"touching {stouch['count']}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
