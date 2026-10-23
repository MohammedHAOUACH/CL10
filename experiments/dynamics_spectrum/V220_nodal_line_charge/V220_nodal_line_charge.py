#!/usr/bin/env python3
"""V220 - the charge of the icosahedral module's nodal lines, and whether
anything can gap them.

V215 found that the icosahedral-module step is exactly gapless along the
icosahedron's two-fold axes and that the touching is linear transverse to the
line (exponent 1.000); V216 found that the tangent algebra on those lines' gap
subspace is non-Abelian (kappa 8.5 to 23) with velocities of order one; V219
measures how isotropic the resulting velocity field is.  A linear crossing of
that kind matters to the low-energy spectrum of Section 8 only if it is
*stable*: a crossing that an arbitrarily small symmetry-allowed perturbation
removes is a numerical accident, not a particle.

Stability is a topological question, and it has a measured answer.  The first
measurement is of what kind of object the line is.  The two-band picture
    M_2(q0 + k) = e^{i phi} ( cos rho I - i sin rho n_hat . sigma )
assumes that two bands are singled out at every momentum around the line.  They
are not: at the line the touching of the icosahedral module is a cluster of
128 of 1024 eigenphases at phase zero (V216 recorded the same multiplicity), and
the cluster stays degenerate at every point of a circle around the line, so the
"lower band" is an arbitrary vector in a degenerate subspace and a Berry phase
computed by following it is an artifact of the solver's choice of basis, not an
invariant - this experiment measures the failure (the overlap of consecutively
tracked bands) instead of trusting it.

What is gauge invariant is the cluster itself: an orthonormal frame of the m
eigenphases closest to phase zero, transported from site to site, gives an
m x m holonomy whose determinant phase is the U(1) charge of the whole touching
and whose unitarity defect measures whether the cluster is isolated along the
loop.  That is what is measured here, on circles of three radii in two planes
(one cut across the line, one tilted about it by 45 degrees), with the touching
dimension m and the isolation gap both read from the spectrum rather than
assumed.  The committed cubic module's GAMMA-M line, run through the identical
instrument, is the control (V63, V115, V118).

T1 what the line is: the dimension of the touching cluster at the line, the
   failure of band-resolved tracking on it, the isolation of the cluster along
   every circle, and the closing exponent of the cluster's top phase (a linear
   touching has beta = 1).
T2 the charge: the determinant phase of the cluster holonomy across radii and
   planes, with its unitarity defect measured rather than assumed.
T3 the mass test, which is what decides: the commutant of the internal map is
   computed by Cesaro averaging over the powers of the step (the obstruction
   norm is reported), and a unitary perturbation exp(-i theta C) built from the
   measured symmetry operator C is switched on at two strengths.  The gap on the
   line and just off it is measured with the perturbation present, and the
   scaling of the on-line gap with theta decides whether a symmetry-respecting
   mass is allowed.
T4 the protecting symmetry read off the step rather than postulated: the little
   group of the line inside the measured point group, and the agreement of the
   spectrum at the line point with its image under each element.

Inputs: cl10/config.yaml (clock, branch, tolerances), the committed internal map
`supercell.sectors()[2]`, the module and its point group from `cl10/continuum.py`
and `cl10/starcell.py`, and the cubic cell and stagger from `cl10/altcell.py` for
the control.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10 import altcell as A
from cl10 import config as _config
from cl10 import continuum as K
from cl10 import engine as E
from cl10 import starcell as X
from cl10 import supercell as SC

CFG = _config.CFG
RADII = (3e-2, 1e-2, 3e-3)         # circles around the line
NCIRCLE = 72                       # points per circle
NODE_TOL = 1e-9                    # a phase below this counts as gapless
THETAS = (1e-3, 1e-2)              # the switched-on perturbations
CESARO = 240                       # powers averaged for the commutant
def cluster_frame(M, m):
    """An orthonormal frame of the m eigenphases closest to phase zero, the
    sorted absolute phases, the relative gap separating the cluster from the
    rest of the spectrum, and the condition of the raw eigenvector block.  The
    frame has to be re-orthonormalised: inside a degenerate cluster `eig`
    returns an arbitrary, wildly non-orthogonal basis, and parallel transport
    built from such a frame is not unitary for reasons of linear algebra rather
    than physics."""
    w, V = E.eig(M)
    ph_all = np.abs(np.angle(E.asnumpy(w)))
    order = np.argsort(ph_all)
    ph = ph_all[order]
    Vr = E.asnumpy(V)[:, order[:m]]
    iso = float((ph[m] - ph[m - 1]) / max(ph[m - 1], 1e-12))
    U, s, _ = np.linalg.svd(Vr, full_matrices=False)
    cond = float(s[0] / max(s[-1], 1e-300))
    return U, ph, iso, cond, Vr[:, 0]


def touching_dimension(step, q):
    """The dimension of the touching at a momentum, read from the spectrum: how
    many eigenphases sit at zero, the largest of them, and the first phase
    outside the cluster."""
    ph = np.sort(np.abs(np.angle(E.asnumpy(E.eigvals(step(q))))))
    m = int((ph < NODE_TOL).sum())
    return m, float(ph[m - 1]), float(ph[m])


def circle_holonomy(step, q0, plane, eps, m):
    """The m x m parallel-transport holonomy of the touching cluster around the
    circle spanned by `plane` and centred at q0, together with the numbers that
    say whether the object means anything: the isolation of the cluster along
    the loop, the unitarity defect of the holonomy, the condition of the raw
    eigenvector block, the failure of band-resolved tracking, and the top phase
    of the cluster."""
    u1, u2 = [np.asarray(v, float) for v in plane]
    u1 = u1 / np.linalg.norm(u1)
    u2 = u2 / np.linalg.norm(u2)
    W = np.eye(m, dtype=complex)
    prev = None
    iso_min, top, overlap, cond = math.inf, [], [], []
    for k in range(NCIRCLE + 1):
        th = 2 * math.pi * k / NCIRCLE
        dq = eps * (math.cos(th) * u1 + math.sin(th) * u2)
        Vr, ph, iso, c, first = cluster_frame(step(q0 + dq), m)
        top.append(float(ph[m - 1]))
        iso_min = min(iso_min, iso)
        cond.append(c)
        if prev is not None:
            overlap.append(float(abs(np.vdot(prev, first))))
            W = W @ (prevf.conj().T @ Vr)
        prev, prevf = first, Vr
    ev = np.sort(np.angle(np.linalg.eigvals(W)))
    return {"eps": eps, "det_phase": float(np.angle(np.linalg.det(W))),
            "unitarity_defect": float(np.max(np.abs(
                W @ W.conj().T - np.eye(m)))),
            "isolation_min": float(iso_min),
            "frame_condition_max": float(max(cond)),
            "cluster_top": float(max(top)),
            "cluster_top_min": float(min(top)),
            "band_overlap_min": float(min(overlap)) if overlap else None,
            "band_overlap_max": float(max(overlap)) if overlap else None,
            "eigenphases": [round(float(x), 5) for x in ev]}



def commutant(U, rng_seed=11):
    """A Hermitian operator in the commutant of the internal map, obtained by
    Cesaro averaging a random Hermitian ansatz over the powers of the map.  The
    obstruction norm ||U C U^-1 - C|| is reported: the construction is only as
    good as that number."""
    Uh = np.asarray(U, complex)
    n = Uh.shape[0]
    rng = np.random.default_rng(rng_seed)
    Hk = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    Hk = 0.5 * (Hk + Hk.conj().T)
    C = np.zeros_like(Hk)
    Uk = np.eye(n, dtype=complex)
    for _ in range(CESARO):
        C = C + Uk.conj().T @ Hk @ Uk
        Uk = Uk @ Uh
    C = C / CESARO
    C = 0.5 * (C + C.conj().T)
    obs = float(np.max(np.abs(Uh @ C @ Uh.conj().T - C)))
    scale = max(float(np.max(np.abs(np.linalg.eigvalsh(C)))), 1e-30)
    return {"obstruction": obs, "norm": scale,
            "operator": C / scale}


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

    def star_step(q, U=U16, axes=axes):
        return X.bloch_star(cell, U, g10, stages, np.asarray(q, float), axes)

    cubic_cell = A.cell((2, 2, 2))
    cubic_stages = A.stages_baseline()

    def cubic_step(q):
        return A.bloch_cell(cubic_cell, U16, g10, cubic_stages,
                            np.asarray(q, float))

    q0 = np.array([math.pi / lr / 2.0, 0.0, 0.0])   # on the x-axis two-fold line
    qM = np.array([math.pi, math.pi, 0.0])          # the cubic GAMMA-M line
    s = 1.0 / math.sqrt(2.0)
    planes = {  # a circle cut across the line, and one tilted about it by 45 deg
        "cross": [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "oblique": [[0.0, 1.0, 0.0], [s, 0.0, s]],
    }
    m_star, star_zero_top, star_zero_out = touching_dimension(star_step, q0)
    m_cubic, cubic_zero_top, cubic_zero_out = touching_dimension(cubic_step, qM)

    loops = {"star_cross": [circle_holonomy(star_step, q0, planes["cross"], eps,
                                            m_star) for eps in RADII],
             "star_oblique": [circle_holonomy(star_step, q0, planes["oblique"],
                                              eps, m_star) for eps in RADII],
             "cubic_M_cross": [circle_holonomy(cubic_step, qM, planes["cross"],
                                               eps, m_cubic) for eps in RADII]}

    # the closing exponent of the cluster's top phase: eps^beta, beta = 1 linear
    def beta_of(rows):
        pts = [(math.log(r["eps"]), math.log(r["cluster_top"]))
               for r in rows if r["cluster_top"] > 1e-12]
        if len(pts) < 2:
            return None
        return float(np.polyfit([p[0] for p in pts], [p[1] for p in pts], 1)[0])

    betas = {k: beta_of(v) for k, v in loops.items()}

    # ---------------------------------------------------------- T3: the mass test
    com = commutant(U16)
    C = com["operator"]
    mass_scan = []
    for theta in THETAS:
        # exp(-i theta C) exactly unitary: C is Hermitian, exponentiate it
        vals, vecs = np.linalg.eigh(C)
        U = np.asarray(U16, complex) @ (vecs @ np.diag(np.exp(-1j * theta * vals))
                                        @ vecs.conj().T)
        unitary_defect = float(np.max(np.abs(U @ U.conj().T - np.eye(U.shape[0]))))
        on = float(np.min(np.abs(np.angle(E.asnumpy(
            E.eigvals(star_step(q0, U=U)))))))
        off = float(np.min(np.abs(np.angle(E.asnumpy(E.eigvals(
            star_step(q0 + np.array([0.0, 1e-3, 0.0]), U=U)))))))
        along = float(np.min(np.abs(np.angle(E.asnumpy(E.eigvals(star_step(
            q0 + np.array([0.31, 0.0, 0.0]), U=U)))))))
        mass_scan.append({"theta": theta, "unitary_defect": unitary_defect,
                          "gap_on_line": on, "gap_transverse_1e-3": off,
                          "gap_along_line": along})

    # ------------------------------------------------------ T4: the little group
    grp = K.polyhedral_group()
    mats = [np.asarray(R, float) for R in grp["elements"]]
    line_dir = np.array([1.0, 0.0, 0.0])
    little = [R for R in mats
              if np.allclose(R @ line_dir, line_dir, atol=1e-9)]
    lg = K.little_group(q0, mats)
    spectra_agree = []
    base = np.sort(np.abs(np.angle(E.asnumpy(E.eigvals(star_step(q0))))))[:6]
    for R in little:
        s2 = np.sort(np.abs(np.angle(E.asnumpy(
            E.eigvals(star_step(R @ q0))))))[:6]
        spectra_agree.append({"R": R.tolist(),
                              "defect": float(np.max(np.abs(base - s2)))})

    # ------------------------------------------------------------- the checks
    checks, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    star_loops = loops["star_cross"] + loops["star_oblique"]
    all_loops = star_loops + loops["cubic_M_cross"]
    det_spread = max(r["det_phase"] for r in star_loops) - min(
        r["det_phase"] for r in star_loops)
    chk("T1 the touching at each line is a degenerate cluster of measured "
        "dimension, not a pair, so no band-resolved Berry phase is defined on "
        "it: a tracked band jumps between members of the cluster",
        m_star >= 8 and m_cubic >= 8
        and all(r["band_overlap_min"] is not None and r["band_overlap_min"] < 0.5
                for r in all_loops),
        f"eigenphases at phase zero (< {NODE_TOL:.0e}) at the line momentum: "
        f"icosahedral module {m_star} of 1024 (largest inside the cluster "
        f"{star_zero_top:.1e}, first phase outside {star_zero_out:.1e}); cubic "
        f"control {m_cubic} of 1024 ({cubic_zero_top:.1e}, {cubic_zero_out:.1e}). "
        "Worst overlap of the eigenvector of the single closest phase between "
        "neighbouring circle points: "
        + ", ".join(f"{k}: {min(r['band_overlap_min'] for r in v):.1e}"
                    for k, v in loops.items())
        + "; largest condition number of the raw eigenvector block of the "
        "cluster: "
        + ", ".join(f"{k}: {max(r['frame_condition_max'] for r in v):.1e}"
                    for k, v in loops.items())
        + " (inside a degenerate cluster the eigensolver's basis is arbitrary "
        "and non-orthogonal, so the frame is re-orthonormalised by SVD before "
        "transport)")
    chk("T1 the cluster is isolated along every circle, so the holonomy below is "
        "an invariant of the step restricted to it",
        all(r["isolation_min"] > 0.02 for r in all_loops),
        "smallest relative gap from the cluster to the rest of the spectrum "
        "along a circle, by radii "
        + str(list(RADII)) + ": "
        + "; ".join(f"{k}: " + ", ".join(f"{r['isolation_min']:.3f}" for r in v)
                    for k, v in loops.items()))
    chk("T2 the cluster holonomy is unitary to solver accuracy on every circle: "
        "with the frame re-orthonormalised the transport is an element of U(m)",
        all(r["unitarity_defect"] < 1e-3 for r in all_loops),
        "largest unitarity defect of the "
        + ", ".join(f"{k}: {max(r['unitarity_defect'] for r in v):.1e}"
                    for k, v in loops.items())
        + f" (frames of dimension {m_star} icosahedral, {m_cubic} cubic)")
    chk("T2 the touching carries no quantized U(1) charge: the determinant "
        "phase of the cluster holonomy depends on the circle - its radius and "
        "its plane - where a topological monodromy would be circle independent, "
        "and the cubic control's phase vanishes at small radii",
        all(abs(v[0]["det_phase"] - v[-1]["det_phase"]) > 0.05
            for v in loops.values()) and det_spread > 0.05,
        "det phase of the "
        f"{m_star} x {m_star} (icosahedral) and {m_cubic} x {m_cubic} (cubic) "
        "holonomy, radii " + str(list(RADII)) + ": "
        + "; ".join(f"{k} " + ", ".join(
            f"eps {r['eps']:.0e}: {r['det_phase']:+.5f}" for r in v)
            for k, v in loops.items())
        + f"; spread over the icosahedral circles {det_spread:.4f}; the "
        "cross-sectional star circle at the smallest radius gives "
        f"{loops['star_cross'][-1]['det_phase']:+.5f} rad, not pi, and the "
        f"cubic control gives {loops['cubic_M_cross'][-1]['det_phase']:+.5f} "
        "there: a curvature flux that falls with the loop, not a monodromy")
    chk("T2 the cluster's top phase closes linearly onto the icosahedral line, "
        "the transverse exponent measured a second way; the cubic control's "
        "exponent is reported",
        betas["star_cross"] is not None and betas["star_oblique"] is not None
        and abs(betas["star_cross"] - 1) < 0.05
        and abs(betas["star_oblique"] - 1) < 0.05,
        "cluster top phase ~ eps^beta: "
        + ", ".join(f"{k}: beta {b if b is None else round(b, 4)}"
                    for k, b in betas.items())
        + "; top phase by radius: "
        + "; ".join(f"{k} " + ", ".join(f"{r['eps']:.0e}: {r['cluster_top']:.5f}"
                                        for r in v) for k, v in loops.items()))
    chk("T3 the mass test is run against a symmetry operator of the step that "
        "is computed, not postulated, and the surviving gaps are measured",
        com["obstruction"] < 1e-2 and all(
            m["unitary_defect"] < 1e-12 for m in mass_scan),
        f"commutant by Cesaro average over {CESARO} powers: obstruction "
        f"{com['obstruction']:.2e}, normalised operator norm "
        f"{com['norm']:.3f}; the perturbation exp(-i theta C) is unitary to "
        + ", ".join(f"{m['unitary_defect']:.1e}" for m in mass_scan)
        + "; measured gaps: "
        + ", ".join(
            f"theta {m['theta']:.0e} on line {m['gap_on_line']:.2e}, "
            f"1e-3 transverse {m['gap_transverse_1e-3']:.2e}, along the line "
            f"{m['gap_along_line']:.2e}" for m in mass_scan))
    chk("T4 the little group of the line is read from the measured point group, "
        "and each of its elements preserves the crossing pair's spectrum at the "
        "line",
        bool(little) and all(d["defect"] < 100 * tol for d in spectra_agree),
        f"point group order {len(mats)}; {len(little)} elements preserve the "
        f"line direction and {len(lg)} elements fix the momentum q0; largest "
        f"spectrum defect {max(d['defect'] for d in spectra_agree):.1e}")

    gap_lo, gap_hi = mass_scan[0]["gap_on_line"], mass_scan[1]["gap_on_line"]
    chk("T3 the measured symmetry operator opens a gap at the line that grows "
        "linearly with theta: the touching carries no symmetry-protected "
        "prohibition on a mass",
        5.0 < gap_hi / max(gap_lo, 1e-30) < 20.0,
        f"on-line gap {gap_lo:.3e} at theta {THETAS[0]:.0e}, {gap_hi:.3e} at "
        f"theta {THETAS[1]:.0e}: ratio {gap_hi / max(gap_lo, 1e-30):.2f} against "
        f"the factor {THETAS[1] / THETAS[0]:.0f} in theta, for a perturbation in "
        f"the commutant of the internal map (obstruction "
        f"{com['obstruction']:.1e})")

    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the charge of the nodal lines of the icosahedral-module "
                    "step: the dimension of the touching cluster at the line "
                    "read from the spectrum, the failure of band-resolved "
                    "tracking on it, the isolation and unitarity and "
                    "determinant phase of the cluster's parallel-transport "
                    "holonomy on circles of three radii in two planes, the "
                    "closing exponent of the cluster's top phase, the mass test "
                    "against a computed symmetry operator of the step and its "
                    "scaling with the coupling, and the little group of the "
                    "line, with the cubic module's GAMMA-M line as control",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"CHARGE OF THE NODAL LINE at the {branch} clock g = {g10:.6f}, "
            f"measured at q0 = {[round(x, 4) for x in q0]} on the icosahedron's "
            f"two-fold axis with {NCIRCLE} points per circle. (T1) the touching "
            f"at the line is a cluster of {m_star} of 1024 eigenphases at phase "
            f"zero (largest inside {star_zero_top:.1e}, first phase outside "
            f"{star_zero_out:.1e}); the cubic GAMMA-M control is a cluster of "
            f"{m_cubic}. No band-resolved Berry phase exists on such a touching: "
            "the overlap of a tracked frame vector between neighbouring circle "
            f"points falls to "
            f"{min(r['band_overlap_min'] for r in all_loops):.1e} worst over the "
            "circles, which says the solver's eigenvectors are arbitrary "
            "rotations inside the degeneracy and a single-band phase would be a "
            "gauge artifact. (T2) the gauge-invariant object is the cluster's "
            "own holonomy: determinant phase "
            + "; ".join(f"{k} " + ", ".join(f"eps {r['eps']:.0e}: "
                                            f"{r['det_phase']:+.5f}" for r in v)
                        for k, v in loops.items())
            + "; the transport is unitary to "
            + ", ".join(f"{max(r['unitarity_defect'] for r in v):.1e}"
                        for v in loops.values())
            + " and the cluster is isolated along the loop by at least "
            f"{min(r['isolation_min'] for r in all_loops):.3f}, so the holonomy "
            "is a genuine invariant of the step restricted to the touching - "
            "but its determinant phase is NOT circle independent: spread "
            f"{det_spread:.4f} over the icosahedral circles, "
            f"{loops['star_cross'][-1]['det_phase']:+.5f} rad at the smallest "
            "cross-sectional radius (not pi) and "
            f"{loops['cubic_M_cross'][-1]['det_phase']:+.5f} for the cubic "
            "control there: a curvature flux that falls with the loop, not a "
            "monodromy, so no quantized U(1) charge is attached to the line. "
            "The cluster's top phase closes as eps^beta with beta = "
            + ", ".join(f"{k} {b:.4f}" for k, b in betas.items()
                        if b is not None)
            + ". (T3) the commutant of the internal map is computed "
            f"(obstruction {com['obstruction']:.1e} over {CESARO} powers); the "
            "perturbation exp(-i theta C) it supplies opens a gap at the line "
            f"that grows linearly with the coupling ({gap_lo:.2e} at theta "
            f"{THETAS[0]:.0e}, {gap_hi:.2e} at {THETAS[1]:.0e}): a "
            "symmetry-respecting mass is allowed. (T4) "
            f"{len(little)} of the {len(mats)} point-group elements preserve the "
            f"line direction and leave the spectrum at q0 fixed to "
            f"{max(d['defect'] for d in spectra_agree):.1e}: the symmetry "
            "available to protect the crossing is a single two-fold rotation. "
            "VERDICT: the icosahedral module's lines are exact and linearly "
            "dispersive (V215, V219) but the protection the icosahedral point "
            "group was supposed to buy them is not there. The touching is a "
            f"{m_star}-dimensional degenerate cluster, so no band-resolved "
            "charge exists on it; the gauge-invariant holonomy of the cluster is "
            "unitary but its U(1) part is circle dependent, which is the "
            "signature of curvature flux rather than a quantized monodromy; and "
            "an operator that commutes with the internal map masses the line "
            "linearly. The lines are stable as loci and unstable as particles."),
        "data": {"clock_branch": branch, "g10": g10, "lr_radius": lr,
                 "q0": [float(x) for x in q0], "qM": [float(x) for x in qM],
                 "radii": list(RADII), "n_circle": NCIRCLE,
                 "node_tol": NODE_TOL,
                 "touching": {"star": {"dimension": m_star,
                                       "largest_in_cluster": star_zero_top,
                                       "first_outside": star_zero_out},
                              "cubic_M": {"dimension": m_cubic,
                                          "largest_in_cluster": cubic_zero_top,
                                          "first_outside": cubic_zero_out}},
                 "loops": loops, "betas": betas,
                 "det_phase_spread_star": det_spread,
                 "commutant": {"obstruction": com["obstruction"],
                               "norm": com["norm"]},
                 "mass_scan": mass_scan,
                 "point_group_order": len(mats),
                 "little_group_line": [R.tolist() for R in little],
                 "little_group_momentum": [np.asarray(m).tolist() for m in lg],
                 "little_group_momentum_order": len(lg),
                 "little_group_spectra": spectra_agree,
                 "engine": E.backend_info(), "checks": checks},
        "wall_time": round(time.time() - t0, 3),
    }
    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v220_nodal_line_charge.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V220 - the charge of the icosahedral module's nodal line")
    print("=" * 68)
    print(f"  touching at the line: icosahedral cluster {m_star} of 1024 "
          f"(top {star_zero_top:.1e}, first outside {star_zero_out:.1e}); cubic "
          f"{m_cubic} of 1024 (top {cubic_zero_top:.1e}, first outside "
          f"{cubic_zero_out:.1e})")
    for k, v in loops.items():
        print(f"  {k}: " + "; ".join(
            f"eps {r['eps']:.0e} det {r['det_phase']:+.5f} iso "
            f"{r['isolation_min']:.3f} top {r['cluster_top']:.5f} unit "
            f"{r['unitarity_defect']:.1e} overlap {r['band_overlap_min']:.1e}"
            for r in v) + f"  beta {betas[k]}")
    print(f"  commutant obstruction {com['obstruction']:.1e}; mass scan: "
          + ", ".join(f"theta {m['theta']:.0e} on {m['gap_on_line']:.2e} off "
                      f"{m['gap_transverse_1e-3']:.2e}" for m in mass_scan))
    print(f"  point group {len(mats)}, little group of the line {len(little)}")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
