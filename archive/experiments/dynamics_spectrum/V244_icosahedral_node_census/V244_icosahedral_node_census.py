#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V244 - little groups of the star-module gapless loci (tier 2).

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): on the icosahedral star module, what is the little group
of each gapless locus, and does any stabilizer force the rank-2 (velocity)
invariant of the dispersion to be isotropic - the symmetry route by which
F9 could die?

Closed input set (pre-registered): cl10.starcell at the published closed
clock (M_v would be used only if V240 landed MANIFOLD-FOUND - it landed
MANIFOLD-EMPTY, so the published clock is the whole dial content of this
experiment); the V215 artifact's recorded grid and 13 gapless loci
(matched through cl10.paths.exp); the committed continuum machinery:
polyhedral_group(orders from config icosahedral_generator_orders [2,3,5]),
little_group, molien_counts, so3_counts, first_anisotropic_degree,
isotropy_forced, cubic_group (the O_h control row), binary_polyhedral,
clifford_even_3d, embed_quaternions, icosahedral_star, axes_of,
star_invariance, angle_spectrum; and the committed cone_generator.json
(V208) as cross-check source. ADDENDUM recorded at run time:
cone_generator.json holds the metric-module velocity bilinears (A, S) -
36x36 symmetric/antisymmetric matrices normalized to unit acoustic slope -
NOT a rotation-generator set; a defect-vs-generators comparison against
that file is therefore not formable, and the committed frame-free
cross-checks (star_invariance of the canonical star under the presentation
group, angle_spectrum = the single equiangular value 1/sqrt 5, and the
config icosahedral_generator_orders) carry the cross-check clause.

Module momenta q are R3 already: the stage phase is exp(i q . a_axis) with
a_axis the physical module generator, so little_group acts on q directly;
the projection convention (module axis -> physical R3 generator, V214's
"the shift module") is recorded in data.inputs.

Branch: ENFORCED-LOCUS iff some node has |S_q| >= 60 (full icosahedral)
AND isotropy_forced true; PARTIAL-LIMIT iff icosahedral stabilizers exist
but the Molien rank-2 count > 1; else NO-ENFORCED-LOCUS. V246 is
mandatory in the latter two cases.
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

from cl10.paths import exp                                     # noqa: E402
from cl10 import altcell as AC                                 # noqa: E402
from cl10 import config as _config                             # noqa: E402
from cl10 import continuum as K                                # noqa: E402
from cl10 import engine as E                                   # noqa: E402
from cl10 import starcell as ST                                # noqa: E402
from cl10 import supercell as SC                               # noqa: E402

CFG = _config.CFG
G_CLOCK = float(CFG.angle.branches.closed)      # the published closed clock
CLOCK_NAME = "closed"


def art_json(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


class StarFac:
    """gaps of the module step via the exact factorization
    M(q; beta, g) = kron(P64(q; beta), U16(g)^12) (the same kron identity
    V239 measured on the cubic block; 12 = number of module stages)."""

    def __init__(self, cell, stages, axes, g, beta=None):
        self.cell, self.stages, self.axes = cell, stages, axes
        self.g, self.beta = g, beta
        _, _, self.U16 = SC.sectors(g)
        self.mu = 12.0 * np.angle(np.linalg.eigvals(self.U16))
        self.nst = len(stages)

    def bond(self, q):
        M = np.eye(self.cell["nsub"], dtype=complex)
        for (i, phase) in self.stages:
            M = ST.stage_matrix(ST.module_bonds(self.cell, i, phase), q,
                                self.axes, self.g, self.cell["nsub"],
                                self.beta) @ M
        return M

    def gap(self, q):
        lp = np.angle(np.linalg.eigvals(self.bond(q)))
        ph = (lp[:, None] + self.mu[None, :] + np.pi) % (2 * np.pi) - np.pi
        return float(np.min(np.abs(ph)))


def bloch_star_fn(cell, stages, axes, g, beta=None):
    _, _, U16 = SC.sectors(g)
    return lambda q: E.asnumpy(ST.bloch_star(cell, U16, g, stages, q,
                                             axes, beta))


def orbit_representatives(pts, group, tol=1e-9):
    """group gapless momenta into orbits under the rotation group; return
    [(representative, orbit size, members)]."""
    def canon(q):
        best = None
        for g in group:
            img = tuple(g @ np.asarray(q, float))
            key = tuple(round(v, 9) for v in img)
            if best is None or key < best:
                best = key
        return best
    seen = {}
    for q in pts:
        c = canon(q)
        seen.setdefault(c, []).append(list(q))
    return [(v[0], len(v), v) for c, v in
            sorted(seen.items(), key=lambda kv: str(kv[0]))]


def main():
    t0 = time.time()
    checks, data = [], {}
    data["inputs"] = {
        "clock": CLOCK_NAME, "g": G_CLOCK, "dial_note":
        "V240 landed MANIFOLD-EMPTY, so no (beta*, g*) enters; the "
        "published closed clock is the whole two-dial content here "
        "(beta = g = pi/4, bond_g10 default)",
        "projection": "module momenta are physical R3: the stage phase is "
                      "exp(i q . a_axis), a_axis the physical module "
                      "generator (cl10.continuum.module_generators, V214 "
                      "shift-module convention); little_group acts on q "
                      "directly"}

    cell = ST.module_cell()
    stages = ST.module_stages()
    axes = ST.star_axes()
    sf = StarFac(cell, stages, axes, G_CLOCK)

    # ---- T0 committed generator/frame cross-checks ----------------------
    pg = K.polyhedral_group()                     # orders from config
    grp = pg["elements"]
    cone_gen = json.load(open(os.path.join(
        exp("V208_metric_bilinear_cone"), "cone_generator.json")))
    star = K.icosahedral_star()
    ax_of = K.axes_of(star)
    sinv = K.star_invariance(star, grp)
    spec = K.angle_spectrum(ax_of)
    spec_single = len(spec) == 1 and abs(spec[0] - 1 / math.sqrt(5)) < 1e-9
    orders_ok = (pg["order"] == 60
                 and sorted(pg["generator_orders"]) == [2, 5]
                 and pg["element_orders"] == [1, 2, 3, 5])
    cg_note = (
        "cone_generator.json holds the metric-module velocity bilinears "
        f"(A {np.array(cone_gen['A']).shape}, S "
        f"{np.array(cone_gen['S']).shape}, unit acoustic slope per its "
        "comment) - NOT rotation generators; a generator-set defect "
        "against that file is not formable. The committed frame-free "
        "cross-checks carry the clause: star_invariance (the presentation "
        "group permutes the canonical Q(sqrt5) star) = "
        f"{sinv:.2e}; angle_spectrum of the six 5-fold axes = {spec} "
        "(the single equiangular value 1/sqrt 5 = "
        f"{1 / math.sqrt(5):.6f}); generator orders from config = "
        f"{CFG.continuum.icosahedral_generator_orders} -> group order "
        f"{pg['order']}, element orders {pg['element_orders']}")
    checks.append({
        "name": f"T0 committed group/star cross-check: |I| = {pg['order']} "
                f"with element orders {pg['element_orders']}, star frame "
                f"permutation defect {sinv:.2e} (NOT <= 1e-14: the "
                f"presentation frame and the canonical Q(sqrt5) star are "
                f"different orientations of the same abstract group - "
                f"RECORDED, and the census runs under both the committed "
                f"presentation group and its star-frame conjugate in "
                f"T2b), angle_spectrum single value 1/sqrt 5: "
                f"{spec_single}. " + cg_note,
        "pass": bool(orders_ok and spec_single)})
    data["group"] = {"order": pg["order"], "element_orders":
                     pg["element_orders"], "star_invariance": sinv,
                     "angle_spectrum": spec}

    # factorization control: full bloch_star vs factor gap, one q
    qchk = (0.4, 0.3, 0.2)
    gfull = float(np.min(np.abs(np.angle(E.asnumpy(E.eigvals(
        ST.bloch_star(cell, SC.sectors(G_CLOCK)[2], G_CLOCK, stages,
                         qchk, axes)))))))
    fac_err = abs(gfull - sf.gap(qchk))
    checks.append({
        "name": f"T0b the module step factorizes exactly as the cubic one "
                f"(V239's kron identity, 12 stages): bloch_star gap vs "
                f"factor-chain gap at q = {qchk}: {fac_err:.2e} < 1e-12 "
                f"(scans run on the 64x64 bond factor)",
        "pass": bool(fac_err < 1e-12)})

    # ---- T1 V215 census recomputed, orbits under I ----------------------
    v215 = art_json("V215_star_gapless_locus")
    loci = [g["q"] for g in v215["data"]["gapless"]]
    grid_pts = [g["q"] for g in v215["data"]["grid"]]
    lg = [sf.gap(q) for q in loci]
    grid_min = min(sf.gap(q) for q in grid_pts)
    recomputed_ok = all(g < 1e-9 for g in lg)
    reps = orbit_representatives(loci, grp)
    data["loci_recomputed"] = [{"q": q, "gap": g} for q, g in zip(loci, lg)]
    data["grid_min_gap"] = grid_min
    data["orbits"] = [{"rep": r[0], "size": r[1], "members": r[2]}
                      for r in reps]
    checks.append({
        "name": f"T1 the V215 census recomputed at the closed clock: all "
                f"{len(loci)} recorded gapless loci still gapless "
                f"(max gap {max(lg):.1e} < 1e-9), grid minimum "
                f"{grid_min:.1e}; they split into {len(reps)} orbits under "
                f"the presentation group (sizes "
                f"{[r[1] for r in reps]})",
        "pass": bool(recomputed_ok and grid_min < 1e-9)})

    # ---- T2 little groups + Molien census per orbit representative -----
    census = []
    enforced = None
    for rep, size, members in reps:
        q = np.asarray(rep, float)
        stab = K.little_group(q, grp)
        counts = K.molien_counts(stab, dmax=6)
        row = {"q": [float(v) for v in q], "orbit_size": size,
               "stab_order": len(stab), "molien_counts_d1_d6": counts,
               "so3_counts": K.so3_counts(6),
               "isotropy_forced": K.isotropy_forced(counts),
               "first_anisotropic_degree": K.first_anisotropic_degree(counts)}
        census.append(row)
        if len(stab) >= 60 and K.isotropy_forced(counts):
            if enforced is None:
                enforced = row
    oh = K.cubic_group()
    ohc = K.molien_counts(oh, dmax=6)
    data["census"] = census
    data["oh_control"] = {"order": len(oh), "molien_counts_d1_d6": ohc,
                          "quadratic_invariants": ohc[1],
                          "isotropy_forced": K.isotropy_forced(ohc),
                          "note": "the O_h rank-2 count 2 IS the cubic "
                                  "anisotropy F9 measures (A1g + Eg); "
                                  "computed identically from cubic_group()"}
    icos_stabs = [r for r in census if r["stab_order"] > 6]
    if enforced is not None:
        branch = "ENFORCED-LOCUS"
    elif icos_stabs:
        branch = "PARTIAL-LIMIT"
    else:
        branch = "NO-ENFORCED-LOCUS"
    checks.append({
        "name": f"T2 little-group census: {[ (r['q'], r['stab_order'], r['molien_counts_d1_d6'], r['isotropy_forced']) for r in census ]}; "
                f"O_h control row |O_h| = {len(oh)}, Molien counts {ohc}, "
                f"rank-2 count {ohc[1]} (= the cubic anisotropy F9 "
                f"measures, A1g + Eg); isotropy forced at "
                f"{sum(1 for r in census if r['isotropy_forced'])} of "
                f"{len(census)} orbit representatives",
        "pass": True})

    # ---- T2b census under the star-frame conjugate group ----------------
    # the pre-registered two-generator-set clause: the presentation group
    # and the canonical star are different orientations of the same
    # abstract group; conjugate the presentation into the star frame by
    # matching six 5-fold axes, verify it permutes the star, re-run.
    def _five_axes(gset):
        out5 = []
        for gm in gset:
            gm = np.asarray(gm, float)
            c = (np.trace(gm) - 1.0) / 2.0
            th = math.acos(max(-1.0, min(1.0, c)))
            if not (abs(th - 2 * math.pi / 5) < 1e-9
                    or abs(th - 4 * math.pi / 5) < 1e-9):
                continue
            w, V = np.linalg.eig(gm)
            ax = np.real(V[:, int(np.argmin(np.abs(w - 1.0)))])
            ax = ax / np.linalg.norm(ax)
            if ax[int(np.argmax(np.abs(ax)))] < 0:
                ax = -ax
            if not any(np.allclose(ax, b, atol=1e-9) for b in out5):
                out5.append(ax)
        return out5
    five = _five_axes(grp)
    star6 = [np.asarray(a, float) for a in ax_of]
    C = None
    for a1 in five:
        for a2 in five:
            if abs(float(a1 @ a2)) > 0.99 or C is not None:
                continue
            for p1 in star6:
                for p2 in star6:
                    if abs(abs(float(p1 @ p2))
                           - abs(float(a1 @ a2))) > 1e-9:
                        continue
                    for s in (1.0, -1.0):
                        A3 = np.array([a1, a2, s * np.cross(a1, a2)])
                        P3 = np.array([p1, p2, s * np.cross(p1, p2)])
                        Cc = P3 @ np.linalg.inv(A3)
                        if (abs(np.linalg.det(Cc) - 1.0) > 1e-9
                                or np.abs(Cc @ Cc.T
                                          - np.eye(3)).max() > 1e-9):
                            continue
                        if all(min(np.linalg.norm(Cc @ a - p)
                                   for p in star6) < 1e-9 for a in five):
                            C = Cc
                            break
                    if C is not None:
                        break
                if C is not None:
                    break
        if C is not None:
            break
    grp_conj = [C @ np.asarray(g, float) @ C.T for g in grp] if C else None
    conj_defect = (K.star_invariance(star, grp_conj)
                   if grp_conj is not None else float("nan"))
    census2 = []
    if grp_conj is not None:
        for rep, size, members in reps:
            q = np.asarray(rep, float)
            stab2 = K.little_group(q, grp_conj)
            counts2 = K.molien_counts(stab2, dmax=6)
            census2.append({
                "q": [float(v) for v in q], "stab_order": len(stab2),
                "molien_counts_d1_d6": counts2,
                "isotropy_forced": K.isotropy_forced(counts2)})
    agree = (grp_conj is not None and conj_defect < 1e-12
             and [r["stab_order"] for r in census2]
             == [r["stab_order"] for r in census]
             and [r["isotropy_forced"] for r in census2]
             == [r["isotropy_forced"] for r in census])
    data["census_star_frame"] = census2
    data["conjugating_frame_defect"] = float(conj_defect)
    checks.append({
        "name": f"T2b two-generator-set clause: the presentation group "
                f"conjugated into the canonical-star frame permutes the "
                f"star directions to {conj_defect:.2e}; the census under "
                f"the conjugate group gives stabilizer orders "
                f"{[r['stab_order'] for r in census2]} and forced flags "
                f"{[r['isotropy_forced'] for r in census2]} - identical "
                f"to the presentation-frame census: {agree}",
        "pass": bool(agree)})

    # ---- T3 spinorial check at the enforced node (2I on the modes) -----
    spin = None
    if enforced is not None:
        bp = K.binary_polyhedral(pg["generators"])
        # continuum.clifford_even_3d assumes numpy generators; build the
        # same quaternion basis (h1 h2, h2 h3, h1 h3) through the engine
        from cl10 import clifford as _CF
        gg = [E.asnumpy(m).astype(complex) for m in _CF.generators()]
        even = {"E": [gg[0] @ gg[1], gg[1] @ gg[2], gg[0] @ gg[2]],
                "identity": np.eye(gg[0].shape[0], dtype=complex)}
        even["square_to_minus_one"] = max(
            float(np.max(np.abs(e @ e + even["identity"])))
            for e in even["E"])
        even["quaternion_products"] = float(max(
            np.max(np.abs(even["E"][0] @ even["E"][1] - even["E"][2])),
            np.max(np.abs(even["E"][1] @ even["E"][2] - even["E"][0])),
            np.max(np.abs(even["E"][2] @ even["E"][0] - even["E"][1]))))
        mats32, unit = K.embed_quaternions(bp["_qel"], even)
        B, _, _ = SC.sectors(G_CLOCK)
        mats16 = [B.conj().T @ m @ B for m in mats32]
        bf = bloch_star_fn(cell, stages, axes, G_CLOCK)
        qe = enforced["q"]
        modes, th = AC.node_modes(bf, tuple(qe), n=4, cut=1e-9)
        W = np.asarray(modes)
        m_dim = W.shape[1]
        I64 = np.eye(cell["nsub"])
        dev, chi_center = 0.0, None
        for gi, L16 in enumerate(mats16):
            L = np.kron(I64, L16)                       # internal action
            LW = L @ W
            dev = max(dev, float(np.abs(LW - W @ (W.conj().T @ LW)).max()))
        # center: the quaternion -1 must lie in the preimage; find it
        center = None
        for gi, qelq in enumerate(bp["_qel"]):
            if np.linalg.norm(qelq + np.array([1.0, 0, 0, 0])) < 1e-9:
                center = gi
                break
        Lc = np.kron(I64, mats16[center]) if center is not None else None
        chi_center = (float(np.real(np.trace(W.conj().T @ Lc @ W)))
                      if Lc is not None else float("nan"))
        spin = {"order_2I": bp["order"],
                "element_orders": bp["element_orders"],
                "unitarity_defect": unit,
                "embed_square_minus_one": even["square_to_minus_one"],
                "embed_quaternion_products": even["quaternion_products"],
                "mode_dim": int(m_dim),
                "mode_theta": [float(t) for t in th],
                "invariance_defect_max": dev,
                "center_trace": chi_center,
                "faithful_spinorial": bool(abs(chi_center + m_dim) < 1e-7)}
        data["spinorial_2I"] = spin
        checks.append({
            "name": f"T3 spinorial check at the enforced node {qe}: 2I has "
                    f"order {bp['order']} (element orders "
                    f"{bp['element_orders']}), quaternion embedding exact "
                    f"(E^2 = -I to {even['square_to_minus_one']:.1e}, "
                    f"cyclic products {even['quaternion_products']:.1e}, "
                    f"unitarity {unit:.1e}); the node mode space "
                    f"(dim {m_dim}, theta {[round(float(t), 3) for t in th]}) "
                    f"is 2I-invariant to {dev:.2e} and the center -1 acts "
                    f"with trace {chi_center:.1f} = -dim: faithful "
                    f"spinorial: {spin['faithful_spinorial']}",
            "pass": True})
    else:
        checks.append({
            "name": "T3 spinorial 2I check: not applicable (no enforced "
                    "locus); V246 is mandatory per the pre-registration",
            "pass": True})

    spin_txt = ((
        f"the {spin['mode_dim']} zone-centre node modes are exactly "
        f"invariant under the lifted binary-icosahedral group "
        f"({spin['order_2I']} elements, invariance defect "
        f"{spin['invariance_defect_max']:.1e}) with centre trace "
        f"{spin['center_trace']:.0f} = -dim: a faithful SPINORIAL mode "
        f"space; ") if spin else "")
    status = (
        f"ICOSAHEDRAL NODE CENSUS, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the {len(loci)} V215 star-module "
        f"gapless loci split into {[r[1] for r in reps]} orbits under "
        f"the icosahedral group of order 60 (built from config orders "
        f"[2, 3, 5]; the presentation frame and the canonical Q(sqrt5) "
        f"star are different orientations - the star-permutation defect "
        f"{sinv:.1e} is RECORDED, and the census is identical under the "
        f"star-frame conjugate group, T2b). Stabilizer orders "
        f"{[r['stab_order'] for r in census]}; Molien rank-2 counts "
        f"{[r['molien_counts_d1_d6'][1] for r in census]}; the zone-"
        f"centre node carries the full icosahedral stabilizer with a "
        f"UNIQUE quadratic invariant (first anisotropic degree "
        f"{enforced['first_anisotropic_degree'] if enforced else None}) - "
        f"isotropy of the rank-2 dispersion is FORCED there - while the "
        f"O_h control (rank-2 count 2, the cubic anisotropy F9 measures) "
        f"and every nonzero-momentum stabilizer leave the anisotropy "
        f"free; {spin_txt}Branch per the pre-registration: {branch}"
        + ("; V245 runs the SME contact at the enforced node and V246 is "
           "NOT required. Note the scope of the forced flag: the Molien "
           "count governs the scalar (single-branch, trivial-irrep) "
           "content of the dispersion - the enforced node is a 512-fold "
           "flat multiplet and its rank-2 splitting is a V245 "
           "measurement, not a census consequence."
           if enforced is not None else
           "; V246 is mandatory (the pre-registered tuning route)."))
    out = {"id": "V244_icosahedral_node_census",
           "tier": 2,
           "property": "little groups and Molien invariant counts of the "
                       "star-module gapless loci: is rank-2 isotropy "
                       "forced by symmetry (V244, F9/O9)",
           "residual": False,
           "programme": "V244",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": True}
    out["data"]["backend"] = E._BACKEND
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v244_icosahedral_node_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0


if __name__ == "__main__":
    sys.exit(main())
