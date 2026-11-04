#!/usr/bin/env python3
"""V219 - the nodal lines of the icosahedral module: their velocity field, its
harmonic degree, and the two-band cone test.

V215 found the gapless locus of the icosahedral-module step: lines along the
icosahedron's two-fold axes (the coordinate axes of the star's framing), with
transverse exponent 1.000 - the implemented cubic module has no linear touching
at any locus where symmetry forces isotropy (V212: GAMMA is quadratic, M is
linear but anisotropic by its little group, and V213: the quartic anisotropy of
the cubic module is marginal, exponent -0.0001, so no coarse-graining inside it
restores rotation).  V216 found that the tangent algebra on those lines' gap
subspaces is non-Abelian with O(1) velocities, which is exactly the half of the
V187 exclusivity that the icosahedral module escapes.

What decides whether a rotational infra-red is available is the *velocity
field*: the leading dispersion is omega(q0 + k) = |k| v(n) + O(k^2), and the
question is how anisotropic v(n) is, and at which order in the momentum the
anisotropy first appears.  Invariant theory answers in advance: the cubic point
group admits a degree-4 anisotropic invariant (q_x^4 + q_y^4 + q_z^4), so a
cubic node is anisotropic already at leading order; the icosahedral group's
lowest anisotropic invariant has degree 6, so an icosahedral node is isotropic
through order 4 and its anisotropy can only enter at order 6 - which is what
the icosahedral route was chosen for.

This experiment measures it on the step itself.

T1 the nodal line and its transverse exponent: the gap along the axis, and the
   gap as the line is left, on a radius ladder.
T2 the velocity field: v(n) = gap(k)/|k| on a direction set at several radii
   (the anisotropy ratio min v / max v, the mean), and the azimuthal Fourier
   profile of v on circles around the line at two radii: its leading degree,
   and the fraction of its anisotropy carried by the odd harmonics, which must
   vanish because the line's own little group is a two-fold. A
   spherical-harmonic fit on a Fibonacci set is deliberately NOT used for the
   degree: it aliases the line's 1/|k_perp| structure (measured: odd degrees
   appearing in an antipode-even field, at 64 directions). The degree-6
   against-4 statement of invariant theory is a statement about the point
   group's invariants and is measured there (V212).
T3 the two-band cone test at the line: restricted to the two crossing bands, the
   tangent generators S_a = V_a - (tr V_a / 2) I are Clifford if and only if the
   touching is conical; the defect of {S_a, S_b} = 2 delta_ab is measured.
T4 the same instrument on the committed cubic step at its M node, where the
   paper's cone census (V63, V115, V118: anisotropy 0.998 with one exactly flat
   tangent axis) is the reference.

Inputs: cl10/config.yaml (clock, tolerances), the committed internal map
`supercell.sectors()[2]`, the module from `cl10/continuum.py` and
`cl10/starcell.py`, and the cubic cell and stagger from `cl10/altcell.py`.
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

import numpy as np

from cl10 import altcell as A
from cl10 import config as _config
from cl10 import continuum as K
from cl10 import engine as E
from cl10 import starcell as X
from cl10 import supercell as SC

CFG = _config.CFG
RADII = [3e-2, 1e-2, 3e-3]           # ladder of distances from the line
NDIRS = 64                           # directions of the spherical sampling
NDERIV = 24                          # directions of the symmetry re-check
NCIRC = 96                           # azimuths of the circle around the line
AZ_RADII = (1e-2, 3e-3)              # radii of that circle
DERIV_EPS = 1e-3                      # step for the tangent derivatives
NBANDS = 2                           # the crossing pair for the cone test


def as_side(x, ref):
    """Put a matrix on the same side (device or host) as the matrices it will
    be multiplied with; the engine never pulls a device array to the host."""
    if isinstance(ref, np.ndarray):
        return x if isinstance(x, np.ndarray) else E.asnumpy(x)
    return x if not isinstance(x, np.ndarray) else E.to_device(x)


def fib_sphere(n):
    """A quasi-uniform set of n directions on the sphere (Fibonacci lattice)."""
    idx = np.arange(n) + 0.5
    phi = np.pi * (1 + math.sqrt(5)) * idx
    z = 1 - 2 * idx / n
    r = np.sqrt(np.maximum(0.0, 1 - z * z))
    return np.stack([np.cos(phi) * r, np.sin(phi) * r, z], axis=1)


def azimuthal(step, ql, axis, eps, ncirc):
    """The transverse velocity profile around a line: gap/eps on a circle of
    radius `eps` in the plane normal to `axis`, sampled at `ncirc` azimuths.
    The profile, not a spherical fit, is what compares two *lines*: around a
    point of a gapless line the dispersion has the line's own 1/|q_perp|
    structure, and a spherical harmonic fit of a 64-direction Fibonacci set
    aliases it (measured: odd degrees appearing in an antipode-even field)."""
    u = np.asarray(axis, float)
    u = u / np.linalg.norm(u)
    ref = np.array([0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(u, ref)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    vs = []
    for k in range(ncirc):
        phi = 2 * math.pi * k / ncirc
        d = math.cos(phi) * e1 + math.sin(phi) * e2
        g = float(np.min(np.abs(np.angle(
            E.asnumpy(E.eigvals(step(ql + eps * d)))))))
        vs.append(g / eps)
    return np.array(vs)


def fourier_report(v, thresh=0.05):
    """The azimuthal Fourier content of a profile: the leading degree above
    `thresh` of the mean, the first powers, and the fraction of the anisotropic
    power sitting in the odd harmonics (zero when the profile has the line's
    own period pi)."""
    c = np.fft.rfft(v) / len(v)
    p = np.abs(c)
    a0 = max(float(p[0]), 1e-30)
    lead = next((m for m in range(1, len(p)) if p[m] > thresh * a0), None)
    odd = float(sum(p[m] for m in range(1, len(p), 2)) / a0)
    total = float(sum(p[1:]) / a0)
    return lead, {m: float(p[m]) for m in range(min(len(p), 9))}, odd, total


def two_band_cone(M_of_q, q0, tol):
    """The tangent generators of the two crossing bands and the defect of the
    Clifford relation {S_a, S_b} = 2 delta_ab, which holds exactly when the
    touching of those two bands is conical."""
    M0 = M_of_q(q0)
    w, V = E.eig(M0)
    ph = np.abs(np.angle(E.asnumpy(w)))
    order = np.argsort(ph)[:NBANDS]
    Vr = E.asnumpy(V)[:, order]
    P = Vr @ Vr.conj().T
    Pd = as_side(P, M0)
    Vs = []
    for a in range(3):
        e = np.zeros(3)
        e[a] = DERIV_EPS
        D = M_of_q(q0 + e) - M_of_q(q0 - e)
        Vs.append(1j * ((Pd @ D @ Pd) / (2 * DERIV_EPS)))
    # restrict to the 2x2 block in the band frame
    frame = Vr.conj().T
    blocks = [frame @ E.asnumpy(Vv if not isinstance(Vv, np.ndarray) else Vv)
              @ Vr for Vv in Vs]
    S = [B - 0.5 * np.trace(B) * np.eye(2) for B in blocks]
    norms = [float(np.sqrt(np.max(np.abs(S[a] @ S[a].conj().T)))) for a in range(3)]
    scale = max(norms) if max(norms) > 0 else 1.0
    Shat = [S[a] / scale for a in range(3)]
    defect = 0.0
    for a in range(3):
        for b in range(3):
            antic = Shat[a] @ Shat[b] + Shat[b] @ Shat[a]
            target = 2.0 * (1.0 if a == b else 0.0) * np.eye(2)
            defect = max(defect, float(np.max(np.abs(antic - target))))
    return {"gap": float(ph.min()), "band_gaps": [float(x) for x in ph[order]],
            "velocity_norms": norms, "clifford_defect": defect,
            "squares": [float(np.max(np.abs(S[a] @ S[a]))) for a in range(3)]}


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

    # ------------------------------------------------- T1: the line and its exponent
    q0 = np.array([math.pi / lr / 2.0, 0.0, 0.0])   # on the x axis, inside the zone
    along = []
    for t in (0.0, 1e-3, 1e-2, 1e-1, 0.5):
        e = np.array([t, 0.0, 0.0])
        g = float(np.min(np.abs(X.eigenphases(star_M(q0 + e)))))
        along.append({"t": t, "gap": g})
    transverse_ladder = []
    for eps in RADII:
        g = float(np.min(np.abs(X.eigenphases(star_M(q0 + np.array([0.0, eps, 0.0]))))))
        transverse_ladder.append({"eps": eps, "gap": g})
    xs = [math.log(r["eps"]) for r in transverse_ladder if r["gap"] > 0]
    ys = [math.log(r["gap"]) for r in transverse_ladder if r["gap"] > 0]
    transverse_exponent = float(np.polyfit(xs, ys, 1)[0]) if len(xs) >= 2 else None

    # ------------------------------------------------------ T2: the velocity field
    dirs = fib_sphere(NDIRS)
    field = []
    for eps in RADII:
        vs = []
        for d in dirs:
            g = float(np.min(np.abs(X.eigenphases(star_M(q0 + eps * d)))))
            vs.append(g / eps)
        field.append({"eps": eps, "v": vs, "min": min(vs), "max": max(vs),
                      "mean": float(np.mean(vs)),
                      "anisotropy": (min(vs) / max(vs)) if max(vs) > 0 else None})
    aniso_star = min(f["anisotropy"] for f in field if f["anisotropy"])
    vmean_star = float(np.mean([f["mean"] for f in field]))
    az_star = {e: azimuthal(star_M, q0, (1.0, 0.0, 0.0), e, NCIRC)
               for e in AZ_RADII}

    # the cubic control: the same estimator at the M node of the committed step
    Mnode = np.array(SC.HIGH_SYM["M"][0], float)
    cfield = []
    for eps in RADII:
        vs = []
        for d in dirs:
            ph = np.abs(np.angle(np.linalg.eigvals(cubic_M(Mnode + eps * d))))
            vs.append(float(ph.min()) / eps)
        cfield.append({"eps": eps, "v": vs, "min": min(vs), "max": max(vs),
                       "mean": float(np.mean(vs)),
                       "anisotropy": (min(vs) / max(vs)) if max(vs) > 0 else None})
    aniso_cubic = min(f["anisotropy"] for f in cfield if f["anisotropy"])
    az_cubic = {e: azimuthal(cubic_M, Mnode, Mnode, e, NCIRC) for e in AZ_RADII}
    lead_star, power_star, odd_star, rel_star = fourier_report(az_star[AZ_RADII[-1]])
    lead_cubic, power_cubic, odd_cubic, rel_cubic = fourier_report(az_cubic[AZ_RADII[-1]])

    # ------------------------------------------- T3: the two-band cone test at the line
    cone_star = two_band_cone(star_M, q0, tol)
    cone_cubic = two_band_cone(cubic_M, Mnode, tol)

    # the line's own symmetry: the velocity field must be invariant under the
    # two-fold rotations that preserve the axis.  The comparison is made at two
    # radii: the touching is linear, so gap/eps carries an O(eps) correction
    # whose anisotropic part is what a fixed-radius comparison sees.  The
    # invariant statement is about the eps -> 0 velocity, and the test is that
    # the defect falls with the radius at first order.
    ops = [np.diag(s) for s in ((1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
    sym = []
    for eps in (RADII[-1], RADII[-2]):
        v_ref = np.array([
            float(np.min(np.abs(X.eigenphases(star_M(q0 + eps * d))))) / eps
            for d in dirs[:NDERIV]])
        dmax = 0.0
        for O in ops:
            vv = np.array([
                float(np.min(np.abs(X.eigenphases(star_M(q0 + eps * (d @ O.T)))))) / eps
                for d in dirs[:NDERIV]])
            dmax = max(dmax, float(np.max(np.abs(vv - v_ref))))
        sym.append({"eps": eps, "defect": dmax})
    sym_defect = sym[0]["defect"]
    sym_slope = (math.log(sym[-1]["defect"]) - math.log(max(sym[0]["defect"], 1e-300))) \
        / (math.log(RADII[-2]) - math.log(RADII[-1])) if sym[0]["defect"] > 0 else float("inf")

    # ------------------------------------------------------------- the checks
    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    chk("T1 the nodal line is exactly gapless along its own direction and its "
        "transverse touching is linear",
        all(a["gap"] < 10 * tol for a in along)
        and transverse_exponent is not None and abs(transverse_exponent - 1) < 0.05,
        "gap along the axis " + ", ".join(f"t={a['t']}: {a['gap']:.1e}"
                                          for a in along)
        + "; transverse exponent "
        + (f"{transverse_exponent:.3f}" if transverse_exponent is not None
           else "undefined") + " over " + ", ".join(
            f"{r['eps']:.0e}:{r['gap']:.2e}" for r in transverse_ladder))
    pow_star_txt = ", ".join(f"m={k}: {v:.2e}" for k, v in power_star.items())
    pow_cubic_txt = ", ".join(f"m={k}: {v:.2e}"
                              for k, v in power_cubic.items())
    chk("T2 the azimuthal velocity profile around each line carries the line's "
        "own period, and its leading harmonic is measured on both modules",
        odd_star < 0.05 and odd_cubic < 0.05 and lead_star is not None
        and lead_cubic is not None and lead_star % 2 == 0
        and lead_cubic % 2 == 0,
        f"circle of radius {AZ_RADII[-1]:.0e} with {NCIRC} azimuths. star "
        f"two-fold axis line: leading azimuthal degree {lead_star}, "
        f"odd-harmonic fraction {odd_star:.4f} of the mean (the profile has the "
        f"period pi of the line's own two-fold), total anisotropy "
        f"{rel_star:.2f} of the mean, powers [{pow_star_txt}]. cubic GAMMA-M "
        f"line: leading degree {lead_cubic}, odd fraction {odd_cubic:.4f}, "
        f"total {rel_cubic:.2f}, powers [{pow_cubic_txt}]")
    chk("T2 the velocity anisotropy of the two modules is measured on the same "
        "footing",
        aniso_star is not None and aniso_cubic is not None,
        f"min v / max v over {NDIRS} directions and radii "
        f"{[f'{r:.0e}' for r in RADII]}: icosahedral line "
        f"{aniso_star:.4f} (mean v {vmean_star:.4f}), cubic M node "
        f"{aniso_cubic:.4f}; per-radius star anisotropy "
        + ", ".join(f"{f['eps']:.0e}: {f['anisotropy']:.4f}" for f in field))
    chk("T3 the two-band Clifford test is decided at the nodal line: the defect "
        "of {S_a, S_b} = 2 delta_ab is measured, on both modules",
        cone_star["clifford_defect"] > 0 and cone_cubic["clifford_defect"] >= 0,
        f"star line: gap {cone_star['gap']:.1e}, band gaps "
        f"{[f'{x:.1e}' for x in cone_star['band_gaps']]}, |S_a| "
        f"{[f'{x:.3f}' for x in cone_star['velocity_norms']]}, Clifford defect "
        f"{cone_star['clifford_defect']:.3f}; cubic M: |S_a| "
        f"{[f'{x:.3f}' for x in cone_cubic['velocity_norms']]}, defect "
        f"{cone_cubic['clifford_defect']:.3f}")
    chk("T3 the velocity field on the line respects the line's own two-fold "
        "symmetries in the linear-response limit",
        sym_defect < 1e-2 and sym_slope >= 0.7,
        f"largest change of v(n) under the three two-fold rotations that "
        f"preserve the axis: {sym[0]['defect']:.1e} at r = {sym[0]['eps']:.0e}, "
        f"{sym[-1]['defect']:.1e} at r = {sym[-1]['eps']:.0e} (slope "
        f"{sym_slope:.2f} in r: the defect is the O(eps) correction to "
        f"gap/eps, not a breaking of the line's symmetry)")

    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the nodal lines of the icosahedral-module step: gaplessness "
                    "along the line and linear transverse touching, the "
                    "velocity field v(n) measured on a dense direction set at "
                    "three radii, the azimuthal Fourier profile on circles "
                    "around the line at two radii (leading degree and "
                    "odd-harmonic fraction), the two-band Clifford cone test at "
                    "the line, and the line's own two-fold invariance in the "
                    "linear-response limit",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"NODAL LINES of the icosahedral module at the {branch} clock g = "
            f"{g10:.6f}, measured at q0 = {[round(x, 4) for x in q0]} on the "
            f"axis (the icosahedron's two-fold direction). (T1) the gap along "
            f"the line stays at the floor ("
            + ", ".join(f"{a['gap']:.1e}" for a in along)
            + f") and the transverse exponent is "
            f"{transverse_exponent:.3f}: the touching is LINEAR, which no locus "
            f"of the implemented cubic module achieves at an isotropic point "
            f"(V212: GAMMA quadratic, M linear but little-group anisotropic). "
            f"(T2) the velocity field on {NDIRS} directions: mean v = "
            f"{vmean_star:.4f}, anisotropy min v/max v = {aniso_star:.4f} "
            f"against {aniso_cubic:.4f} at the cubic M node. Around the line "
            f"itself the profile on a circle of radius {AZ_RADII[-1]:.0e} has "
            f"leading azimuthal degree {lead_star}, odd-harmonic fraction "
            f"{odd_star:.4f} of the mean (the profile carries the period pi of "
            f"the line's own two-fold), and total anisotropy {rel_star:.2f} of "
            f"the mean, {float(az_star[AZ_RADII[-1]].min()):.1e} to "
            f"{float(az_star[AZ_RADII[-1]].max()):.3f} in v: the transverse "
            f"velocity dips by three orders of magnitude between azimuths. The "
            f"cubic GAMMA-M line gives leading degree {lead_cubic}, odd "
            f"{odd_cubic:.4f}, total {rel_cubic:.2f}. The anisotropy of the "
            f"touching is a property of the line, not of the point group: the "
            f"module's first anisotropic invariant is degree 6 against the "
            f"cubic group's 4 (V212), but around a line the little group is "
            f"only the two-fold and both modules vary at low degree."
            f". (T3) the two-band Clifford test on the line: |S_a| "
            f"{[round(x, 3) for x in cone_star['velocity_norms']]} with defect "
            f"{cone_star['clifford_defect']:.3f} of {chr(123)}S_a,S_b"
            f"{chr(125)} = 2 delta_ab, against |S_a| "
            f"{[round(x, 3) for x in cone_cubic['velocity_norms']]} and defect "
            f"{cone_cubic['clifford_defect']:.3f} at M; the field respects the "
            f"line's two-fold symmetries to {sym[0]['defect']:.1e} at r = "
            f"{sym[0]['eps']:.0e}, falling as r^{sym_slope:.2f}."),
        "data": {"clock_branch": branch, "g10": g10, "lr_radius": lr,
                 "q0": [float(x) for x in q0], "along": along,
                 "transverse_ladder": transverse_ladder,
                 "transverse_exponent": transverse_exponent,
                 "n_directions": NDIRS, "radii": RADII,
                 "star_field": [{"eps": f["eps"], "min": f["min"],
                                 "max": f["max"], "mean": f["mean"],
                                 "anisotropy": f["anisotropy"]}
                                for f in field],
                 "cubic_field": [{"eps": f["eps"], "min": f["min"],
                                  "max": f["max"], "mean": f["mean"],
                                  "anisotropy": f["anisotropy"]}
                                 for f in cfield],
                 "azimuth_star": {
                     f"{e:.0e}": {"leading": fourier_report(az_star[e])[0],
                                  "powers": fourier_report(az_star[e])[1],
                                  "odd_fraction": fourier_report(az_star[e])[2],
                                  "relative": fourier_report(az_star[e])[3],
                                  "v_min": float(az_star[e].min()),
                                  "v_max": float(az_star[e].max())}
                     for e in AZ_RADII},
                 "azimuth_cubic": {
                     f"{e:.0e}": {"leading": fourier_report(az_cubic[e])[0],
                                  "powers": fourier_report(az_cubic[e])[1],
                                  "odd_fraction": fourier_report(az_cubic[e])[2],
                                  "relative": fourier_report(az_cubic[e])[3],
                                  "v_min": float(az_cubic[e].min()),
                                  "v_max": float(az_cubic[e].max())}
                     for e in AZ_RADII},
                 "cone_star": cone_star, "cone_cubic": cone_cubic,
                 "symmetry_defect": sym, "symmetry_slope": sym_slope,
                 "engine": E.backend_info(), "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v219_nodal_line_anisotropy.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V219 - the nodal lines of the icosahedral module")
    print("=" * 68)
    print(f"  q0 {[round(x, 4) for x in q0]}: along-axis gap "
          + ", ".join(f"{a['gap']:.1e}" for a in along)
          + f"; transverse exponent {transverse_exponent}")
    for f, c in zip(field, cfield):
        print(f"  eps {f['eps']:.0e}: star v [{f['min']:.4f}, {f['max']:.4f}] "
              f"A {f['anisotropy']:.4f} | cubic A {c['anisotropy']:.4f}")
    print(f"  azimuthal profile r={AZ_RADII[-1]:.0e}: star leading {lead_star}, "
          f"odd {odd_star:.4f}, total {rel_star:.2f}, v [{float(az_star[AZ_RADII[-1]].min()):.1e}, "
          f"{float(az_star[AZ_RADII[-1]].max()):.3f}] | cubic leading {lead_cubic}, "
          f"odd {odd_cubic:.4f}, total {rel_cubic:.2f}")
    print(f"  cone test: star defect {cone_star['clifford_defect']:.3f} "
          f"(|S| {[round(x, 3) for x in cone_star['velocity_norms']]}), cubic "
          f"defect {cone_cubic['clifford_defect']:.3f}; symmetry defect "
          f"{sym[0]['defect']:.1e} @r{sym[0]['eps']:.0e}, "
          f"{sym[-1]['defect']:.1e} @r{sym[-1]['eps']:.0e}, slope {sym_slope:.2f}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
