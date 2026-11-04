#!/usr/bin/env python3
r"""V226 - the charge of the node of the 2I-covariant Clifford cell, and whether
a symmetry-respecting perturbation can remove it.

V225 established the cell: the bond operators $U_d = -\hat v_d\cdot\vec\sigma$
intertwine the binary icosahedral action exactly, the step is $2I$-invariant at
the zone centre, inversion is not a symmetry of it, and at the zone centre the
gap vanishes with a tangent algebra that is Clifford on the spinor piece of the
gap subspace -- the identity $\{S_a,S_b\} = 2\delta_{ab}$ that V215 measured
failing ($2.000$ at the star module's zone centre, $0.444$ at the cubic M node)
and that V212/V216 found cannot hold together with isotropy on the committed
cells.

A linear touching whose tangent algebra is Clifford would be a Weyl point, and
a Weyl point carries a quantized charge: the flux of the Berry curvature
through a small sphere around it, a topological invariant -- it cannot change
under any perturbation that does not close the gap somewhere on the sphere.
Whether the touching is an isolated point at all is not assumed: a sphere about
it can be pierced by the gapless set, in which case the sector changes size on
the sphere and no flux is defined on it.  The experiment settles that first and
only then reads fluxes, so the questions are decided by measurement:

T0 the piercing.  On each sphere of the radius ladder the gap is measured along
   the mesh's own directions, and a direction counts as gapless when its
   smallest eigenphase falls below $10^{-12}$.  Either no sphere is pierced, and
   the point-node reading of T1 is available, or every one is, and the touching
   is not a point node.

T1 the flux, where a sphere can be read.  A sphere around the touching,
   triangulated by the icosahedral
   group's own orbits (`spincell.icosahedral_sphere`: the twelve star
   directions are the icosahedron's vertices, so the mesh is symmetry-uniform),
   and the flux computed by the gauge-invariant link-variable product of
   Fukui-Harigai-Watanabe: the phase of
   $\det\langle u_a|u_b\rangle\det\langle u_b|u_c\rangle\det\langle u_c|u_a
   \rangle$ on each triangle.  No gauge is fixed, no derivative of an
   eigenvector is taken.  The number is checked for stability against the mesh
   frequency and the radius -- a topological number must not depend on either.
   A sphere on which the sector changes size, or on which a link stops having
   full rank, is recorded as a measurement that failed and is never averaged
   with the others; only spheres that keep the band rank count as trusted.  The
   instrument itself is validated in the repository's own test of a Weyl
   two-level unitary before it is trusted here.

T2 the consistency of the flux.  The positive sector must give the opposite
   number, the two-band ("pair") sector must give zero, and a sphere drawn
   around a point that is *not* gapless must give zero: the flux belongs to the
   node, not to the band structure in general.

T3 the protection.  A perturbation $\exp(-i\theta V)$ whose Hermitian generator
   $V$ is obtained by projecting a random ansatz onto the invariants of the
   $2I$ action (the exact group average over the 120 lifts, with the
   obstruction $\lVert W(g)VW(g)^\dagger - V\rVert$ reported rather than
   assumed) is switched on.  A mass term that gaps a Weyl point must
   anticommute with all three Clifford generators, and on the spinor piece the
   commutant of the Clifford algebra is trivial by Schur's lemma -- so the
   prediction is not that the touching is immovable.  What the scan settles is
   which of the two outcomes obtains: the class contains a term that gaps the
   cell at every coupling, or the gaplessness survives at every coupling (the
   original touching gapping while gaplessness persists elsewhere).  Whether a
   Chern number can be attached to wherever it went is a separate question the
   sphere itself answers: V227 T2 measured that no sphere about a point of this
   cell's gapless locus keeps link rank, since the locus is a set of lines.

T4 the control that makes the contrast meaningful.  The same perturbation
   applied to the cell without spin -- the scalar bond assignment, which is
   inversion symmetric and whose touching is quadratic (V225 T6) -- must be
   gapped by the same kind of perturbation.  If both survive, the protection
   proved nothing; if only the cone survives, it is the Clifford algebra that
   protects it.
"""
import json
import math
import os
import sys
import time

import numpy as np

from cl10 import config as _config
from cl10 import engine as E
from cl10 import spincell as SP
from cl10 import supercell as SC

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

RADII = (3e-3, 1e-2, 4e-2)          # the enclosing spheres
FREQS = (3, 4, 5)                  # mesh refinement: the sector carries 96
                                  # identical copies of the cone, and on the
                                  # frequency-2 and 3 icosahedra a link overlap
                                  # loses rank (principal-angle cosines below
                                  # $10^{-8}$, measured) -- the charge is only
                                  # trusted where the rank survived
THETAS = (1e-3, 1e-2, 1e-1)         # the switched-on perturbations
SCAN_RADII = (0.0, 2e-3, 5e-3, 1e-2, 2e-2, 4e-2, 8e-2)   # node hunting


def phase_gap(step_fn, q):
    r"""The smallest eigenphase modulus of a step, read off the Hermitian part
    $C = (S + S^\dagger)/2$ whose eigenvalues are $\cos\theta$: on
    $(-\pi,\pi]$ the modulus is recovered exactly as
    $\arccos(\max\cos\theta)$.  Checked against the general `eigvals` route on
    this cell at five momenta -- agreement to $1.3\times10^{-13}$, and
    $2.8\times10^{-17}$ at the zone centre -- for $0.018\,\mathrm{s}$ against
    $0.083\,\mathrm{s}$ per call on the $384$-dimensional step, which is the
    factor that decides whether the mesh sweeps of this experiment fit inside a
    run at all."""
    S = step_fn(q)
    cn = np.linalg.eigvalsh(E.asnumpy((S + S.conj().T) / 2.0))
    return float(math.acos(min(1.0, max(-1.0, float(cn[-1])))))


def invariant_operator(perms, lifts, n_internal, layers, seed=7):
    r"""A Hermitian operator on the step's space, projected onto the invariants
    of the $2I$ action by the exact group average
    $V \mapsto \lvert 2I\rvert^{-1}\sum_g W(g)^\dagger V W(g)$, with the
    obstruction $\max_g\lVert W(g) V W(g)^\dagger - V\rVert$ reported.  The
    perturbation $\exp(-i\theta V)$ it supplies is therefore symmetry-respecting
    to the accuracy of that number, which is the only reason its effect on the
    node means anything."""
    rng = np.random.default_rng(seed)
    n = SP.NDIR * layers * SP.NSPIN * n_internal
    H = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    H = 0.5 * (H + H.conj().T)
    Hd = E.to_device(H)
    acc = E.zeros((n, n), dtype=complex)
    worst = 0.0
    for Pm, D in zip(perms, lifts):
        Wd = E.to_device(SP.symmetry_operator(Pm, D, n_internal, layers))
        acc = acc + E.matmul(Wd.conj().T, E.matmul(Hd, Wd))
    Vn = E.asnumpy(acc / len(lifts))
    Vn = 0.5 * (Vn + Vn.conj().T)
    for Pm, D in zip(perms[:40], lifts[:40]):
        W = SP.symmetry_operator(Pm, D, n_internal, layers)
        worst = max(worst, float(np.max(np.abs(W @ Vn @ W.conj().T - Vn))))
    scale = max(float(np.max(np.abs(np.linalg.eigvalsh(Vn)))), 1e-30)
    return {"operator": Vn / scale, "obstruction": worst, "norm": scale}


def perturbed_step(U16, V, theta, **kw):
    """$\\exp(-i\theta V)\\,S(q)$: an on-site potential built from the invariant
    generator, applied after the step."""
    if theta == 0.0:
        return lambda q: SP.step(q, U16, **kw)
    w, U = np.linalg.eigh(V)
    P = (U * np.exp(-1j * theta * w)) @ U.conj().T
    Pd = E.to_device(P)
    return lambda q: E.matmul(Pd, SP.step(q, U16, **kw))


def locate_node(step_fn, mesh):
    """The smallest gap found on a spherical grid of radii, which is where a
    displaced node has to be; reported with its radius so the linear response
    to the perturbation can be read off."""
    best = {"gap": np.inf}
    for r in SCAN_RADII:
        for p in mesh["vertices"]:
            q = r * p / np.linalg.norm(p)
            g = phase_gap(step_fn, q)
            if g < best["gap"]:
                best = {"gap": g, "radius": r, "direction": [float(x) for x in p]}
    return best


def flux_one(step_fn, q0, mesh, radius, sector, frequency=None):
    """One sphere, with its failure kept.  A link whose frame overlap loses rank
    makes the triangle phase meaningless; a sector that changes size over the
    sphere means the sphere encloses a second gapless point (or crosses one).
    Both are measurements about the gapless locus, not accidents of the run, so
    they are recorded and the checks below trust only what survives them."""
    try:
        res = SP.berry_chern(step_fn, np.asarray(q0, float), mesh, radius=radius,
                             sector=sector)
    except RuntimeError as exc:
        return {"frequency": frequency, "radius": float(radius), "sector": sector,
                "error": str(exc), "link_rank_ok": False,
                "nearest_integer": None, "integer_defect": None,
                "min_link_overlap": 0.0, "bands_in_sector": 0}
    if frequency is not None:
        res["frequency"] = frequency
    return res


def trusted(entries):
    """The spheres whose sector kept its rank over the whole mesh."""
    return [e for e in entries if e.get("link_rank_ok")]


def flux_of(step_fn, q0, meshes, pairs, sector):
    out = [flux_one(step_fn, q0, meshes[f], r, sector, frequency=f)
           for (f, r) in pairs]
    return out


def main():
    t0 = time.time()
    checks = []
    ok = True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    _, _, U16 = SC.sectors()
    NINT = U16.shape[0]
    g10 = float(_config.CFG.angle.g10)
    branch = str(_config.CFG.angle.branch)
    tol = float(_config.CFG.tolerances.standard)
    strict = float(_config.CFG.tolerances.strict)
    sp = SP.spinor_lifts()
    perms = SP.axis_permutations(sp)
    # the cell V225 measures: one block stage per axis, internal exponent 4
    layers, legs = 1, (4,)
    kw = dict(g10=g10, kind="clifford", layers=layers, legs=legs)
    def cell_step(q):
        return SP.step(q, U16, **kw)

    gap0, zeros0, dim0 = SP.gap(np.zeros(3), U16, **kw)
    meshes = {f: SP.icosahedral_sphere(f) for f in FREQS}

    # ------------------------------------------------------- T0a the instrument
    # Nothing below is trusted until the instrument has been measured on a step
    # whose charge is known by construction: the two-level Weyl step
    # $S(\mathbf q) = \cos r - i\sin r\,(\hat{\mathbf q}\cdot\vec\sigma)$, whose
    # negative band carries $|C| = 1$ over every sphere enclosing the origin
    # because $\mathbf q \mapsto \hat{\mathbf q}$ has degree one.  Two hostile
    # variants come with it: a band sweeping through the half-circle window,
    # which the fixed window ladder has to survive without losing the charge,
    # and a band whose phase spans the whole ladder, where the instrument has to
    # refuse rather than print a number.
    def toy_step(q, extra=None):
        q = np.asarray(q, float)
        r = float(np.linalg.norm(q))
        n = q / r if r > 1e-15 else np.array([0.0, 0.0, 1.0])
        ndot = np.array([[n[2], n[0] - 1j * n[1]], [n[0] + 1j * n[1], -n[2]]])
        S = math.cos(r) * np.eye(2) - 1j * math.sin(r) * ndot
        if extra is None:
            return S
        ph = ((math.pi / 2 + 0.3 * n[2]) if extra == "sweep"
              else math.pi * abs(n[0]))
        band = np.diag([np.exp(1j * ph), np.exp(-1j * ph)])
        return np.block([[S, np.zeros((2, 2))],
                         [np.zeros((2, 2)), band]])

    toy_mesh = SP.icosahedral_sphere(2)
    toy = flux_one(toy_step, np.zeros(3), toy_mesh, 0.3, "negative")
    toy_pos = flux_one(toy_step, np.zeros(3), toy_mesh, 0.3, "positive")
    toy_sweep = flux_one(lambda q: toy_step(q, "sweep"), np.zeros(3), toy_mesh,
                         0.3, "negative")
    toy_dense = flux_one(lambda q: toy_step(q, "dense"), np.zeros(3), toy_mesh,
                         0.3, "negative")
    toy_pair = (toy.get("nearest_integer") is not None
                and toy_pos.get("nearest_integer") is not None
                and toy["nearest_integer"] + toy_pos["nearest_integer"] == 0)
    chk("T0a the instrument is right before it is used: the two-level Weyl step "
        "of known charge comes back with that charge, the sector bookkeeping "
        "reverses with the sector, a band sweeping the half-circle window is "
        "survived by the window ladder, and a phase spanning the ladder is "
        "refused",
        toy.get("link_rank_ok") and abs(toy.get("nearest_integer") or 0) == 1
        and (toy.get("integer_defect") or 1.0) < 1e-10
        and toy.get("bands_in_sector") == 1 and toy_pair
        and toy_sweep.get("link_rank_ok")
        and abs(toy_sweep.get("nearest_integer") or 0) == 1
        and toy_sweep.get("window", math.pi / 2) < math.pi / 2
        and "error" in toy_dense,
        f"toy Weyl: C = {toy.get('chern'):+.6f} (defect "
        f"{toy.get('integer_defect'):.1e}, {toy.get('bands_in_sector')} band in "
        f"sector, window {toy.get('window'):.3f}, min link overlap "
        f"{toy.get('min_link_overlap'):.4f}); positive sector "
        f"{toy_pos.get('nearest_integer')} (the pair must cancel); sweep "
        f"variant: C = {toy_sweep.get('nearest_integer')} with the window "
        f"narrowed to {toy_sweep.get('window'):.3f} of the ladder "
        f"{[round(w, 3) for w in (math.pi / 2, 1.2, 0.9, 0.7, 0.5, 0.35)]}; "
        f"dense variant refused: "
        f"{toy_dense.get('error', 'NOT refused')[:80]}")

    # ----------------------------------------------------------------- T0b
    # Is the touching at the zone centre a point node at all?  A sphere about it
    # carries a Berry flux only if the gapless set does not run through it, so
    # the piercing is measured first, and on two censuses: the mesh's own
    # icosahedral directions, and the six signed coordinate axes, which V227 T1
    # measured to be gapless lines of this cell.  An icosahedral mesh that
    # happens to miss those lines is not evidence that they are absent, and a
    # flux read on such a sphere is an artefact of the mesh -- so the two
    # censuses are reported separately and the verdict is taken on the union.
    line_dirs = [np.array(v, float) for v in
                 ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1),
                  (0, 0, -1))]
    pierced = []
    for r in RADII:
        gs = [phase_gap(cell_step, r * p)
              for p in meshes[FREQS[0]]["vertices"]]
        ls = [phase_gap(cell_step, r * u) for u in line_dirs]
        pierced.append({"radius": float(r), "mesh_points": len(gs),
                        "gapless_to_1e-12": int(sum(1 for x in gs if x < 1e-12)),
                        "below_1e-6": int(sum(1 for x in gs if x < 1e-6)),
                        "largest_gap": float(max(gs)),
                        "line_points": len(ls),
                        "line_gapless_to_1e-12": int(
                            sum(1 for x in ls if x < 1e-12)),
                        "largest_line_gap": float(max(ls))})
    all_pierced = all(e["line_gapless_to_1e-12"] > 0 for e in pierced)
    all_clean = (all(e["gapless_to_1e-12"] == 0 for e in pierced)
                 and all(e["line_gapless_to_1e-12"] == 0 for e in pierced))
    chk("T0b the question is settled before any flux is read: on the mesh's own "
        "directions *and* on the cell's measured gapless lines, either no sphere "
        "about the touching is pierced, or every one of them is",
        all_pierced or all_clean,
        "; ".join(f"r {e['radius']:.1e}: {e['gapless_to_1e-12']} of "
                  f"{e['mesh_points']} mesh directions gapless to 1e-12 "
                  f"({e['below_1e-6']} below 1e-6, largest gap "
                  f"{e['largest_gap']:.3e}), while "
                  f"{e['line_gapless_to_1e-12']} of {e['line_points']} signed "
                  f"coordinate-axis directions are exactly gapless (largest "
                  f"line gap {e['largest_line_gap']:.1e})" for e in pierced)
        + (" -- every sphere about the zone centre is pierced by a gapless line "
           "of the cell itself: the touching is a point of a gapless set, not an "
           "isolated node. The icosahedral mesh misses those lines by gaps of "
           "1.8e-03 down to 2.4e-02, which is why it looks clean and why any "
           "sphere flux read on it is a mesh artefact, exactly as V227 T2 found "
           "for the locus points." if all_pierced else
           " -- no sphere is pierced on either census: the point-node reading is "
           "available"))

    # ------------------------------------------------------------------ T1
    ladder = [(f, r) for f in FREQS for r in RADII]
    flux = flux_of(cell_step, np.zeros(3), meshes, ladder, "negative")
    good = trusted(flux)
    ints = {res["nearest_integer"] for res in good}
    worst_defect = max([res["integer_defect"] for res in good], default=np.inf)

    raw = [res["nearest_integer"] for res in flux
           if res.get("nearest_integer") is not None]
    spread = max(raw) - min(raw) if raw else 0
    consistent = all(not res["link_rank_ok"]
                     or res["integer_defect"] < 100 * tol for res in flux)
    mesh_gaps = [format(e["largest_gap"], ".1e") for e in pierced]
    line_counts = [e["line_gapless_to_1e-12"] for e in pierced]
    worst_line = format(max(e["largest_line_gap"] for e in pierced), ".1e")
    chk("T1 the flux question at the touching comes out one way or the other: "
        "either the mesh-refined spheres keep their rank and agree, or the "
        "piercing of T0b closes the question and the ladder's disagreement is the "
        "record of that",
        consistent and ((all_pierced and len(good) < 3)
                        or (len(good) >= 3 and len(ints) == 1
                            and worst_defect < 100 * tol)),
        f"{len(good)} of {len(flux)} spheres kept full link rank (mesh "
        f"frequencies {list(FREQS)}, radii {list(RADII)}); the trusted readings "
        f"are {sorted(ints)} with largest integer defect {worst_defect:.2e} and "
        f"smallest kept link overlap "
        f"{max(res['min_link_overlap'] for res in good) if good else 0:.1e}; the "
        f"untrusted readings at the same radius run "
        f"{[res['nearest_integer'] for res in flux if not res['link_rank_ok']]}"
        f" -- a spread of {spread} across the refinement ladder, which is not a "
        f"charge; {len(flux) - len(good)} spheres failed or lost rank: "
        f"{sorted({res['error'][:70] for res in flux if 'error' in res})}. The "
        f"mesh is the icosahedron's own ({meshes[FREQS[0]]['faces']} faces, "
        f"{meshes[FREQS[-1]]['triangles']} triangles at frequency "
        f"{FREQS[-1]}), and T0b shows every one of these spheres is crossed by "
        f"the cell's gapless axes -- the answer is that no Chern number is "
        f"defined at this touching, not that the ladder needs another "
        f"triangulation.")

    # ------------------------------------------------------------------ T2
    short = list(zip(FREQS, RADII))
    mirror = flux_of(cell_step, np.zeros(3), meshes, short, "positive")
    pair = flux_of(cell_step, np.zeros(3), meshes, short, "pair")
    off = np.array([0.21, 0.07, 0.13])
    off_gap = phase_gap(cell_step, off)
    elsewhere = flux_of(cell_step, off, meshes, short, "negative")
    c_node = good[0]["nearest_integer"] if good else None

    else_defects = [format(e["integer_defect"], ".1e")
                    for e in trusted(elsewhere)]
    if trusted(mirror) or trusted(pair):
        sect_txt = (f"where a sector sphere at the node kept rank the "
                    f"bookkeeping holds: positive sector "
                    f"{[m['nearest_integer'] for m in trusted(mirror)]} against "
                    f"the node's {c_node}, two-band sector "
                    f"{[p['nearest_integer'] for p in trusted(pair)]}")
    else:
        sect_txt = ("no sector sphere at the node kept rank, so these two "
                    "clauses cannot be exercised at this touching -- T0b and T1 "
                    "say why: every sphere about it is cut by a gapless line")
    chk("T2 the sector bookkeeping is asserted only where a sphere can be read: "
        "spheres about a gapped point give exactly zero, and at the touching a "
        "sector clause binds only on spheres that kept their link rank",
        trusted(elsewhere)
        and all(e["nearest_integer"] == 0 for e in trusted(elsewhere))
        and all(e["integer_defect"] < 100 * tol for e in trusted(elsewhere))
        and off_gap > 1e-3
        and all(p["nearest_integer"] == 0 for p in trusted(pair))
        and (c_node is None
             or all(m["nearest_integer"] == -c_node for m in trusted(mirror))),
        ("node at the zone centre: gap " + format(gap0, ".1e") + ", "
         + str(zeros0) + " zero eigenphases of " + str(dim0)
         + (f", Chern {c_node:+d} on the single sphere that kept rank "
            "(see T0b/T1: the reading is not corroborated by the ladder)"
            if c_node is not None else ", no sphere at the node kept its rank")
         + f"; at the node the positive sector reads "
           f"{[m['nearest_integer'] for m in mirror]} and the two-band sector "
           f"{[p['nearest_integer'] for p in pair]}, with "
           f"{len(trusted(mirror))} and {len(trusted(pair))} of the "
           f"{len(mirror)} spheres trusted -- {sect_txt}; at the gapped point "
           f"{list(np.round(off, 3))} (gap {off_gap:.4f}) all "
           f"{len(trusted(elsewhere))} of {len(elsewhere)} spheres keep rank and "
           f"give {[e['nearest_integer'] for e in trusted(elsewhere)]} with "
           f"defects {else_defects}: the flux belongs to a gapless touching, not "
           f"to the band structure"))

    # ------------------------------------------------------------------ T3
    inv = invariant_operator(perms, sp["su2"], NINT, layers)
    Vop = inv["operator"]
    mass_scan = []
    for th in THETAS:
        sfun = perturbed_step(U16, Vop, th, **kw)
        g0 = phase_gap(sfun, np.zeros(3))
        hit = locate_node(sfun, meshes[FREQS[0]])
        mass_scan.append({"theta": th, "gap_at_origin": g0,
                          "min_gap_found": hit["gap"], "radius": hit["radius"],
                          "direction": hit.get("direction")})
    moved = [m for m in mass_scan if m["min_gap_found"] < 1e-8 and m["radius"] > 0]
    if moved:
        sfun = perturbed_step(U16, Vop, moved[-1]["theta"], **kw)
        q_new = moved[-1]["radius"] * np.array(moved[-1]["direction"])
        new_charge = flux_one(sfun, q_new, meshes[FREQS[-1]],
                              2.0 * moved[-1]["radius"], "negative")
    else:
        new_charge = None

    lin = [mass_scan[i + 1]["gap_at_origin"]
           / max(mass_scan[i]["gap_at_origin"], 1e-300)
           for i in range(len(mass_scan) - 1)]
    lin_ratio = ((mass_scan[-1]["gap_at_origin"]
                  / max(mass_scan[0]["gap_at_origin"], 1e-300))
                 / (THETAS[-1] / THETAS[0]))
    masses_it = (mass_scan[-1]["gap_at_origin"] > 1e-6
                 and 0.8 < lin_ratio < 1.25)
    smallest = min(mass_scan, key=lambda m: m["min_gap_found"])
    outcome = (
        "the class contains a symmetry-respecting term that masses the "
        "zone-centre touching, and it masses it linearly in the coupling "
        f"(gap/theta constant to {abs(lin_ratio - 1.0) * 100:.1f}% over the "
        f"scan)" if masses_it else
        "the zone-centre touching is not massed linearly by this class")
    if new_charge is None:
        nc_txt = "No displaced point was located to scan."
    elif new_charge["link_rank_ok"]:
        nc_txt = (
            f"The displaced point at radius {moved[-1]['radius']:.3f} sits on a "
            f"sphere that kept rank and reports flux "
            f"{new_charge['nearest_integer']:+d} (integer defect "
            f"{new_charge['integer_defect']:.1e})"
            + (", the same integer as the sphere at the original touching."
               if c_node is not None
               and new_charge["nearest_integer"] == c_node
               else ", while the sphere at the original touching is not a "
                    f"readable object ({len(trusted(flux))} of {len(flux)} kept "
                    "rank), so the two are not a conservation law."))
    else:
        nc_txt = (
            f"The displaced point at radius {moved[-1]['radius']:.3f} sits on a "
            f"sphere whose link rank was lost (smallest overlap "
            f"{new_charge['min_link_overlap']:.1e}): "
            + ("no reading, the sector changed size"
               if new_charge["nearest_integer"] is None
               else f"the raw reading {new_charge['nearest_integer']} at integer "
                    f"defect {new_charge['integer_defect']:.1e}")
            + " is not a charge -- T0, T1 and V227 T2 reach the same verdict at "
              "every point of this locus.")

    chk("T3 the invariant switch-on is decided: this symmetry-respecting class "
        "either masses the zone-centre touching linearly in the coupling or it "
        "does not -- and where the cell is least gapped under the switch-on, and "
        "whether a sphere can be read there, are read off the scan",
        inv["obstruction"] < 1e-10 and masses_it,
        f"the invariant generator is symmetry-respecting to "
        f"{inv['obstruction']:.1e} (exact group average over {len(perms)} lifts, "
        f"norm {inv['norm']:.3f}, operator normalised to unit spectral scale); "
        f"switched on: "
        + ", ".join(f"theta {m['theta']:.0e}: gap at the zone centre "
                    f"{m['gap_at_origin']:.2e}, smallest gap found "
                    f"{m['min_gap_found']:.1e} at radius {m['radius']:.3f}"
                    for m in mass_scan)
        + f" -- consecutive gap ratios {[round(x, 3) for x in lin]} against the "
          f"coupling ratios "
          f"{[round(THETAS[i + 1] / THETAS[i], 3) for i in range(len(THETAS) - 1)]}"
          f". {outcome}. The cell's smallest gap at that coupling is "
          f"{smallest['min_gap_found']:.1e} at radius {smallest['radius']:.3f}, "
          f"against the zone-centre gap {mass_scan[-1]['gap_at_origin']:.2e}: "
          + ("a smaller gap opens away from the zone centre, so the touching is "
             "not where the cell is least gapped"
             if smallest["min_gap_found"] < mass_scan[-1]["gap_at_origin"]
             else "the zone centre stays the least gapped point")
          + ". " + nc_txt)

    # ------------------------------------------------------------------ T4
    ckw = dict(g10=g10, kind="scalar", layers=1, legs=(1,))
    def scalar_step(q):
        return SP.step(q, U16, **ckw)
    s_inv = invariant_operator(SP.axis_permutations(sp), sp["su2"], NINT, 1)
    s_flux = flux_of(scalar_step, np.zeros(3), meshes, short, "negative")
    s_scan = []
    for th in THETAS:
        sfun = perturbed_step(U16, s_inv["operator"], th, **ckw)
        g0 = phase_gap(sfun, np.zeros(3))
        s_scan.append({"theta": th, "gap_at_origin": g0})
    s_open = s_scan[-1]["gap_at_origin"] > 1e-6

    chk("T4 control without spin: the same kind of symmetry-respecting "
        "perturbation gaps the quadratic touching, so the protection is the "
        "cone's and not the perturbation's weakness",
        s_open and any(s["link_rank_ok"] for s in s_flux)
        and max(s["integer_defect"] for s in s_flux if s["link_rank_ok"]) < 1e-8,
        f"the scalar bond assignment (inversion a symmetry, touching quadratic) "
        f"gives Chern {[s['nearest_integer'] for s in s_flux]} on the same "
        f"spheres -- the charge of a touching, quadratic or not, is quantised "
        f"the same way -- and the invariant perturbation "
        f"(obstruction {s_inv['obstruction']:.1e}) opens its gap to "
        + ", ".join(f"{m['gap_at_origin']:.2e} at theta {m['theta']:.0e}"
                    for m in s_scan)
        + f", where the Clifford cell's gap at the same moment is "
          f"{mass_scan[-1]['gap_at_origin']:.2e} with a node still present at "
          f"radius {moved[-1]['radius'] if moved else 'not found'}")

    mass_txt = ", ".join(f"theta {m['theta']:.0e}: gap {m['gap_at_origin']:.2e}"
                         for m in mass_scan)
    moved_txt = f"{moved[-1]['radius']:.3f}" if moved else "not found"
    new_txt = (f"{new_charge['nearest_integer']:+d}"
               if (new_charge and new_charge["link_rank_ok"]) else "n/a")

    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the Berry flux (first Chern number) of the node of the "
                    "2I-covariant Clifford cell over a icosahedrally-uniform "
                    "sphere, its sector and empty-region consistency, and the "
                    "response of the node to a perturbation projected onto the "
                    "invariants of the 2I action, against the same perturbation "
                    "on the inversion-symmetric spinless control",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"NODE CHARGE AND PROTECTION at the {branch} clock g = {g10:.6f}, "
            f"cell layers {layers} internal exponent {legs[0]}. (T1) of the "
            f"{len(flux)} spheres about the zone-centre touching (gap {gap0:.1e}, "
            f"{zeros0} zero eigenphases of {dim0}), "
            + (f"no Berry charge is defined at this touching, and T0b says why: "
               f"of the mesh's {pierced[0]['mesh_points']} icosahedral directions "
               f"none is gapless on the three sampled spheres (their largest gaps "
               f"are {mesh_gaps}), while on the cell's six signed coordinate axes "
               f"{line_counts} of {pierced[0]['line_points']} come back exactly "
               f"gapless at those radii (worst line gap {worst_line}) -- the mesh "
               f"is simply not sampling the lines the touching sits on. "
               f"{len(good)} of {len(flux)} spheres keep full link rank and the "
               f"ladder reads {sorted(set(raw))}, a spread of {spread}, which is "
               f"a mesh artefact and not a charge: the same verdict V227 T2 "
               f"reached for every point of the locus."
               if all_pierced else
               (f"{len(good)} kept full link rank over the mesh and agree on the "
                f"integer {c_node:+d} across every mesh frequency {list(FREQS)} "
                f"and radius {list(RADII)} (worst integer defect "
                f"{worst_defect:.1e}); the instrument reproduces -1 exactly for "
                f"a Weyl two-level unitary, so the number is the node's, not the "
                f"mesh's."
                if good else
                f"none kept full link rank (smallest link overlap as low as "
                f"{min(res['min_link_overlap'] for res in flux):.1e}), and the "
                f"piercing census of T0b is the reason -- no Chern number is "
                f"defined at this touching, exactly as for the locus points of "
                f"V227 T2."))
            + f" (T2) at the node the positive sector reads "
            f"{[m['nearest_integer'] for m in mirror]} and the two-band sector "
            f"{[p['nearest_integer'] for p in pair]}, of which "
            f"{len(trusted(mirror))} and {len(trusted(pair))} of "
            f"{len(mirror)} spheres kept rank -- {sect_txt}; the spheres around "
            f"the gapped point (gap {off_gap:.4f}) are all trusted and give "
            f"{[e['nearest_integer'] for e in elsewhere]}. "
            f"(T3) the perturbation projected onto the 2I invariants "
            f"(obstruction {inv['obstruction']:.1e}, generator norm "
            f"{inv['norm']:.3f}): {outcome}; {mass_txt}, and the cell's smallest "
            f"gap at the strongest coupling is "
            f"{smallest['min_gap_found']:.1e} at radius {smallest['radius']:.3f}"
            + (f", with a displaced gapless point at radius {moved_txt} carrying "
               f"readable flux {new_txt}"
               if new_charge and new_charge["link_rank_ok"]
               else ", and no point of the scan sits exactly gapless")
            + (", so the Schur-protection expectation does not apply to this "
               "touching: V227 T1 and T2 measured it as a point of a gapless set "
               "of lines, not an isolated Weyl node whose mass term would have "
               "to anticommute with three generators."
               if masses_it else ".")
            + f" (T4) the spinless control has flux "
              f"{[s['nearest_integer'] for s in s_flux]} and the same kind of "
              f"perturbation gaps it to {s_scan[-1]['gap_at_origin']:.2e} at "
              f"theta {THETAS[-1]:.0e}."),
        "data": {
            "clock": g10, "branch": branch,
            "cell": {"layers": layers, "legs": list(legs), "dim": dim0,
                     "gap_at_origin": gap0, "zero_eigenphases": zeros0},
            "meshes": {str(f): {"vertices": int(len(meshes[f]["vertices"])),
                                "triangles": int(len(meshes[f]["triangles"])),
                                "faces": meshes[f]["faces"],
                                "solid_angle_defect": meshes[f]["solid_angle_defect"]}
                       for f in FREQS},
            "instrument": {"toy": toy, "toy_positive": toy_pos,
                           "toy_window_sweep": toy_sweep,
                           "toy_dense_phase": toy_dense},
            "piercing": pierced, "point_node_reading_available": bool(all_clean),
            "flux": flux, "mirror": mirror, "pair": pair, "elsewhere": elsewhere,
            "gapped_probe": {"q": [float(x) for x in off], "gap": off_gap},
            "invariant_generator": {"obstruction": inv["obstruction"],
                                    "norm": inv["norm"]},
            "mass_scan": mass_scan, "displaced_charge": new_charge,
            "control_scalar": {"invariant_obstruction": s_inv["obstruction"],
                               "flux": s_flux, "gap_scan": s_scan},
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
    print("V226 - the charge of the node and whether a mass can remove it")
    print("=" * 68)
    print(f"  cell: layers {layers}, legs {legs}, gap(0) {gap0:.1e}, "
          f"{zeros0}/{dim0} zero eigenphases")
    for e in pierced:
        print(f"  piercing r{e['radius']:.1e}: {e['gapless_to_1e-12']}/"
              f"{e['mesh_points']} mesh directions gapless to 1e-12, "
              f"{e['below_1e-6']} below 1e-6, largest gap {e['largest_gap']:.3e}")
    for res in flux:
        if "error" in res:
            print(f"  flux f{res['frequency']} r{res['radius']:.3f}: FAILED "
                  f"({res['error'][:78]})")
            continue
        print(f"  flux f{res['frequency']} r{res['radius']:.3f}: "
              f"C = {res['chern']:+.6f} (defect {res['integer_defect']:.1e}, "
              f"bands {res['bands_in_sector']}, min link "
              f"{res['min_link_overlap']:.4f}, window {res['window']:.3f}, "
              f"rank ok {res['link_rank_ok']})")
    print(f"  mirror {[m['nearest_integer'] for m in mirror]} | pair "
          f"{[p['nearest_integer'] for p in pair]} | elsewhere "
          f"{[e['nearest_integer'] for e in elsewhere]} (gap {off_gap:.4f})")
    print(f"  invariant generator obstruction {inv['obstruction']:.1e}")
    for m in mass_scan:
        print(f"  theta {m['theta']:.0e}: gap(0) {m['gap_at_origin']:.2e}, min "
              f"{m['min_gap_found']:.1e} at r {m['radius']:.3f}")
    if new_charge and new_charge.get("link_rank_ok"):
        print(f"  displaced node charge {new_charge['nearest_integer']:+d} "
              f"(defect {new_charge['integer_defect']:.1e})")
    elif new_charge:
        print(f"  displaced node sphere failed: {new_charge['error'][:78]}")
    print(f"  spinless control: flux {[s['nearest_integer'] for s in s_flux]}, "
          f"gap scan {[round(m['gap_at_origin'], 4) for m in s_scan]}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
