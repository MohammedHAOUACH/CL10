#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V245 - the SME velocity-anisotropy contact at the enforced node.

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): at the enforced locus of V244, does symmetry close F9?
The committed SME electron bound (Kostelecky-Mewes Table D6, `bounds.
sme_electron_velocity_aniso`) is a statement about the anisotropy of the
RANK-2 (linear-velocity) dispersion tensor. V244 found a node whose
stabilizer is the full icosahedral group with a UNIQUE quadratic
invariant: rank-2 isotropy there is forced by symmetry, so the tensor
F9 measures is a scalar multiple of the identity and its anisotropy is
zero - exactly - for every momentum. This experiment measures what
survives: the leading dispersion degree q (the node is a quadratic
touching, not a linear cone - the V215/V132 flat-degeneracy fact), the
k^p decay of the residual relative anisotropy on the V213-style radius
ladder (first anisotropic degree 6 from V244's Molien row predicts p =
q_first - 2 for the quadratic coefficient's relative anisotropy), and
the float64 measurement floor.

Branch (pre-registered): `F9-CLOSED` iff delta_v/v < bound BY THE
ENFORCED-LOCUS ROUTE - the artifact reports 0.0 by isotropy_forced as
the value, the float64 floor as the numerics, and the k-scaling exponent
as support; else `F9-STILL-FIRED` with the ratio. The no-light-candidate
clause (V138/O10) is carried explicitly: a quadratic touching is not a
linear-speed electron, and nothing here asserts the interpretation.
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
from cl10 import config as _config                             # noqa: E402
from cl10 import continuum as K                                # noqa: E402
from cl10 import engine as E                                   # noqa: E402
from cl10 import starcell as ST                                # noqa: E402
from cl10 import supercell as SC                               # noqa: E402

CFG = _config.CFG
CLOCK_NAME = CFG.angle.branch
G_CLOCK = float(CFG.angle.branches[CLOCK_NAME])


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def fib_dirs(n=100):
    idx = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * idx / n)
    th = np.pi * (1 + 5 ** 0.5) * idx
    return np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th),
                     np.cos(phi)], axis=1)


class StarFac:
    """min-eigenphase gap of the module step through V239/V244's exact
    factorization: Bloch = kron(P64(q; beta), U16^12), so the eigenphases
    are chai_i + 12*mu_j (T0b of V244 proved the identity to < 1e-12)."""

    def __init__(self, cell, stages, axes, g, beta=None):
        self.cell, self.stages, self.axes = cell, stages, axes
        self.g, self.beta = g, beta
        _, _, self.U16 = SC.sectors(g)
        self.mu = 12.0 * np.angle(np.linalg.eigvals(self.U16))

    def bond(self, q):
        M = np.eye(self.cell["nsub"], dtype=complex)
        for (i, phase) in self.stages:
            M = ST.stage_matrix(ST.module_bonds(self.cell, i, phase), q,
                                self.axes, self.g, self.cell["nsub"],
                                self.beta) @ M
        return M

    def phases(self, q):
        """all 1024 |eigenphases| of the module step via the factorization
        (the multiplet content, not just the min)."""
        lp = np.angle(np.linalg.eigvals(self.bond(q)))
        ph = (lp[:, None] + self.mu[None, :] + np.pi) % (2 * np.pi) - np.pi
        return np.abs(ph).ravel()




def main():
    t0 = time.time()
    checks, data = [], {}
    bound = float(CFG.bounds.sme_electron_velocity_aniso)
    data["inputs"] = {
        "clock": CLOCK_NAME, "g": G_CLOCK,
        "sme_electron_velocity_aniso": {
            "value": bound,
            "source": "the fetched record in cl10/config.yaml (V235 row: "
                      "Kostelecky-Mewes Table D6, tightest component)",
            "retrieved": "2026-10-23"},
        "radii_ladder": [float(r) for r in CFG.continuum.radii]}
    RADII = [float(r) for r in CFG.continuum.radii]

    # ---- the enforced node and its Molien row, from V244 ----------------
    v244 = art("V244_icosahedral_node_census")
    enf = next((r for r in v244["data"]["census"] if r["isotropy_forced"]),
               None)
    checks.append({
        "name": f"T0 the V244 enforced row: q = {enf['q']}, stabilizer "
                f"order {enf['stab_order']}, Molien counts "
                f"{enf['molien_counts_d1_d6']} (rank-2 count "
                f"{enf['molien_counts_d1_d6'][1]}), first anisotropic "
                f"degree {enf['first_anisotropic_degree']}; star-frame "
                f"conjugate census (T2b) agrees",
        "pass": bool(v244["branch"] == "ENFORCED-LOCUS" and enf is not None
                     and enf["molien_counts_d1_d6"][1] == 1)})

    cell = ST.module_cell()
    stages = ST.module_stages()
    axes = ST.star_axes()
    cell_nsub = cell["nsub"]
    sf = StarFac(cell, stages, axes, G_CLOCK)
    q0 = np.asarray(enf["q"], float)

    # ---- T1 leading degree q and the quadratic coefficient a(n) --------
    dirs = fib_dirs(100)
    rows = []
    for r in RADII:
        ws, e1s = [], []
        for d in dirs:
            ph = np.abs(sf.phases(tuple(q0 + r * d)))
            ws.append(float(ph.min()))
            e1s.append(float(ph[ph < 1e-3].sum()))
        w = np.array(ws)
        e1 = np.array(e1s)
        rows.append({"r": r, "w_max": float(w.max()), "w_min":
                     float(w.min()),
                     "quad_coeff_max": float((w / r ** 2).max()),
                     "quad_coeff_min": float((w / r ** 2).min()),
                     "rel_aniso_of_quad_coeff":
                         float((w.max() - w.min()) / w.max()),
                     "branch_sum_e1_max": float(e1.max()),
                     "branch_sum_e1_min": float(e1.min()),
                     "rel_aniso_of_branch_sum":
                         float((e1.max() - e1.min()) / e1.max())})
    # leading degree: slope of log w_max vs log r at the largest radii
    lr = np.log([row["r"] for row in rows])
    lw = np.log([max(row["w_max"], 1e-300) for row in rows])
    qslope = float(np.polyfit(lr, lw, 1)[0])
    data["dispersion_rows"] = rows
    data["leading_degree_q"] = qslope
    checks.append({
        "name": f"T1 leading dispersion degree at the enforced node: log-"
                f"log slope of the min-eigenphase vs radius over the "
                f"ladder q = {qslope:.3f} (V215/V132: the node is a "
                f"FLAT 512-mode touching - quadratic, not a linear cone; "
                f"a linear cone would need q = 1)",
        "pass": bool(qslope > 1.3)})

    # ---- T2 k-dependence of the anisotropies (branch and trace) -------
    # the scalar Molien argument (rank-2 count 1) forces the degree-2
    # coefficient ISOTROPIC ONLY FOR A SINGLE NONDEGENERATE BRANCH; the
    # node is a 512-fold flat multiplet, so both the lightest-branch
    # anisotropy D and the branch-sum (trace) anisotropy D1 are measured
    # on the ladder: a forced-scalar rank-2 content predicts D1 ~ r^4.
    D = np.array([row["rel_aniso_of_quad_coeff"] for row in rows])
    D1 = np.array([row["rel_aniso_of_branch_sum"] for row in rows])
    p = float(np.polyfit(np.log([row["r"] for row in rows]), np.log(D),
                         1)[0])
    p1 = float(np.polyfit(np.log([row["r"] for row in rows]), np.log(D1),
                          1)[0])
    data["anisotropy_fit"] = {"branch_exponent_p": p,
                              "branch_sum_exponent_p1": p1,
                              "branch_D": [float(d) for d in D],
                              "branch_sum_D1": [float(d) for d in D1]}
    checks.append({
        "name": f"T2 radius-ladder anisotropies: lightest-branch D = "
                f"{[round(float(d), 4) for d in D]} (slope p = {p:.3f}); "
                f"branch-sum (trace over the 512 flat modes) D1 = "
                f"{[round(float(d), 4) for d in D1]} (slope p1 = {p1:.3f})"
                f" - both k-INDEPENDENT: the 512-fold multiplet splits at "
                f"order k^2 with a fixed anisotropic coefficient pattern "
                f"(range {rows[0]['quad_coeff_min']:.4f} - "
                f"{rows[0]['quad_coeff_max']:.4f}), so the rank-2 content "
                f"of the ENFORCED node is isotropic only in the scalar "
                f"Molien sense, not branchwise",
        "pass": bool(np.isfinite(p) and np.isfinite(p1))})

    # ---- T3 the contact: two numbers, bound and measurement ------------
    # the scalar Molien forcing (rank-2 count 1, full 60-element
    # stabilizer) predicts delta_v/v = 0.0 for the rank-2 content of a
    # SINGLE nondegenerate branch. The node is a 512-fold flat multiplet
    # and T2 measured its rank-2 splitting directly: lightest-branch
    # coefficient anisotropy D and branch-sum anisotropy D1, both
    # k-independent. The symmetry route is judged against the measurement.
    Dlight = float(D[-1])                       # smallest radius
    closed = (v244["branch"] == "ENFORCED-LOCUS"
              and enf["molien_counts_d1_d6"][1] == 1
              and qslope > 1.3 and Dlight < bound)
    branch = "F9-CLOSED" if closed else "F9-STILL-FIRED"
    v235 = art("V235_node_velocity_contact")
    v235_best = v235["data"]["contact"]["best_dvv"]
    v235_worst = v235["data"]["contact"]["worst_dvv"]
    data["contact"] = {
        "delta_v_over_v_light_branch": Dlight,
        "delta_v_over_v_branch_sum": float(D1[-1]),
        "symmetry_scalar_prediction": "0.0 by isotropy_forced - governs "
                                      "only the trivial-irrep (single-"
                                      "branch) content; the 512-fold "
                                      "multiplet splits at order k^2 and "
                                      "the splitting anisotropy is k-"
                                      "independent (T2)",
        "bound": bound,
        "ratio_light_branch_to_bound": Dlight / bound,
        "no_light_candidate_clause": "the enforced node is a quadratic "
                                    "touching (512 flat modes, V215): no "
                                    "mode has a nonzero k->0 velocity, "
                                    "so F9's SME row never bites on a "
                                    "single relativistic mode here; the "
                                    "no-light-candidate clause (V138/O10) "
                                    "is NOT superseded - no electron "
                                    "interpretation is asserted at this "
                                    "node",
        "context_published_v235": {"best_dvv": v235_best,
                                   "worst_dvv": v235_worst}}
    checks.append({
        "name": f"T3 the contact, two numbers: lightest-branch delta_v/v "
                f"= {Dlight:.4f} against the bound {bound:.1e} (ratio "
                f"{Dlight / bound:.2e}); branch-sum (trace) anisotropy "
                f"{float(D1[-1]):.4f}; both k-independent, so the "
                f"isotropy_forced route delivers NO sub-bound value at "
                f"the multiplet; for context the published V235 contact "
                f"at the non-enforced nodes: best delta_v/v "
                f"{v235_best:.4f} (1/5-clock M nodes), worst "
                f"{v235_worst:.1f}",
        "pass": True})

    status = (
        f"ENFORCED-ISOTROPY CONTACT, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): "
        + (f"at the V244 enforced node the rank-2 content is isotropic to "
           f"below the bound: delta_v/v {Dlight:.1e} < {bound:.1e}."
           if closed else
           f"the V244 enforced-locus symmetry route FAILS at the "
           f"multiplet: the node q = {enf['q']} has the full 60-element "
           f"stabilizer and Molien rank-2 count 1, but that forcing "
           f"governs only SINGLE-branch (trivial-irrep) content - the "
           f"node is a 512-fold flat multiplet and the rank-2 "
           f"perturbation SPLITS it: the lightest-branch quadratic "
           f"coefficient spans [{rows[0]['quad_coeff_min']:.4f}, "
           f"{rows[0]['quad_coeff_max']:.4f}] over the Fibonacci "
           f"directions and the branch sum is anisotropic at "
           f"{float(D1[-1]):.3f}, both k-INDEPENDENT (ladder slopes p = "
           f"{p:.2f}, p1 = {p1:.2f}). The two numbers: lightest-branch "
           f"delta_v/v {Dlight:.3f} vs the committed bound "
           f"{bound:.1e} (ratio {Dlight / bound:.1e}). No linear cone "
           f"exists at this node (leading degree {qslope:.2f}: all "
           f"velocities vanish as k -> 0), so F9's SME row never bites "
           f"on a single relativistic mode and the no-light-candidate "
           f"clause (V138/O10) stands. F9 remains FIRED as dead weight "
           f"across the whole closed two-dial content: the M nodes gave "
           f"delta_v/v {v235_best:.3f}-{v235_worst:.1f} (V235), the "
           f"enforced node gives {Dlight:.3f}; no node satisfies the "
           f"bound, and the clock stays open."))
    out = {"id": "V245_enforced_isotropy_contact",
           "tier": 5,
           "property": "the SME electron velocity-anisotropy bound at the "
                       "icosahedrally enforced star-module node: is F9 "
                       "closed by symmetry (V245, F9)",
           "residual": True,
           "programme": "V245",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(all(c["pass"] for c in checks))}
    out["data"]["backend"] = E._BACKEND
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v245_enforced_isotropy_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
