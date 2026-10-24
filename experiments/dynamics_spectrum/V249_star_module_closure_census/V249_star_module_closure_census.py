#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V249 - the closure census of the star module (tier 2, F2 clock side).

**The question.**  The supercell closure C (V239) quantized the clock to
pi-multiples and excluded both 1/5 (the fitted clock that reads the
decimals) and 1/2 (the V208 lightlike pin).  The star module's closure was
never censused - only the supercell's C was.  Is there a closed point of
the star step near either contact clock?

**The object.**  At q = 0 the star step with a decoupled bond angle beta and
a Cartan-twisted internal map is

    M_star(GAMMA; beta, g, eta) = kron(B(GAMMA; beta), (U16(g) W(eta))^12),

B(GAMMA; beta) = the product of the 12 phase-major stage matrices with
bond_g10 = beta at q = 0 (cl10/starcell.stage_matrix), U16(g) =
supercell.sectors(g)[2], and W(eta) = exp(i eta D), D = cos(theta) X +
sin(theta) Y the committed Wilson direction: the same five Cartan
generators 1j h_a h_{a+1}, a in {0,2,4,6,8}, that V67 uses with the V65
coefficients X_coef = (1,1,-1,-1,-1), Y_coef = (-1/4,-1/4,-1/6,-1/6,-1/6)
and the V65 root angle theta = 1.7307836435290485 (= config
desert.theta_root_rad), restricted to the 16 sector by the V67 projector
basis B exactly as V67 restricts the Cartans (`Dk = B^dag C B`).  No new
parameter: the twist reuses the committed direction; eta is a census dial.

**The factorization claim (measured, T1).**  `bloch_star` is a product of
kron(S_stage, U); by the mixed-product property it factorizes exactly as
M = kron(Pi_t S_t(q), U^12).  The 512-fold flat GAMMA node of V244 is a
consequence.  T1 measures that identity at 8 random q-points.

**The control (T4).**  V214 published the star step at the closed clock
g = pi/4 and recorded its GAMMA gap 1.04083e-17 and unitarity 1.55431e-14,
but never its closure order - no artifact records one.  This control
therefore reproduces V214's two numbers and measures the closure order of
the published step at (beta = g = pi/4, eta = 0) for the first time; the
census then carries that order.

**Closure machinery.**  M is unitary, so for M = kron(B, V12) the spectral
norm ||M^k - I||_2 = max_{j,l} |mu_j^k nu_l^k - 1| over the eigenvalues of
B and V12 = (U W)^12 - an exact evaluation, not an approximation.  T4
verifies it directly (1024-dim matrix_power + svd) at the control point.
The census scans k = 1..16 at tol 1e-12 on the 61x61 (beta, g) grid at the
eta ladder, refines every order-plateau boundary to 1e-10 by bisection,
and refines every closed grid point by 1-D descent (V239's machinery,
reimplemented here on the spectral form; no cl10/ core change).

**Acceptance (pre-registered).**  STAR-CLOCK-OPEN iff any closed point with
k <= 16 lies within 1e-6 of g = 0.2 (the clock that reads 2.0766/0.04) or
g = 0.5 (the V208 lightlike pin); else STAR-CLOCK-QUANTIZED with the
distance table.  Both wordings pre-registered; entries never edited.
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

from cl10 import config as _config                            # noqa: E402
from cl10 import clifford                                     # noqa: E402
from cl10 import engine as E                                  # noqa: E402
from cl10 import starcell as X                                # noqa: E402
from cl10 import supercell as SC                              # noqa: E402
from cl10.paths import exp                                    # noqa: E402

CFG = _config.CFG

CELL = X.module_cell()
STAGES = X.module_stages(order="phase-major")
AXES = np.asarray(X.star_axes(), float)
N_DIAL = 61
LO, HI = 0.02, math.pi / 2
CLOS_TOL = 1e-12
MAX_ORDER = 16
BRANCH_ANGLES = {  # the config `angle:` branch values: on no linspace grid
    "implemented": float(CFG.angle.branches.implemented),
    "closed": math.pi / 4.0,
    "half": math.pi / 2.0,
}
CLOCK_1_5 = 0.2          # the fitted clock that reads 2.0765641 / 0.04
CLOCK_1_2 = 0.5          # the V208 lightlike pin (g = 1/2, closure-dead on C)

# --------------------------------------------------------------------------
# the committed twist direction D (zero new parameters: V67's Cartans, V65's
# coefficients and root angle)
# --------------------------------------------------------------------------


def wilson_direction():
    """D restricted to the 16 sector: D = cos(th) X + sin(th) Y with
    X = sum_k x_coef[k] 1j h_{a_k} h_{a_k + 1}, Y likewise (V67:66-83), the
    Cartans conjugated to the 16 basis exactly as V67:66 does (`B^dag C B`).
    Returns (D16, th_root, x_coef, y_coef)."""
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    B16 = SC.sectors()[0]                       # the 16 basis of sectors()
    Dk = [B16.conj().T @ C @ B16 for C in cartan]
    v65 = json.load(open(exp("V65_wilson", "results", "v65_wilson.json")))
    gg = next(r for r in v65["data"]["part_B"].values()
              if r["su5"] == "Georgi-Glashow")
    x_coef = np.array(gg["X_coef"])
    y_coef = np.array(v65["data"]["Y_coef"])
    root = next(r for r in gg["selective_roots"] if r["physical"])
    th = float(root["theta_rad"])
    assert abs(th - float(CFG.get("desert.theta_root_rad"))) < 1e-12, \
        "V65 root angle disagrees with config desert.theta_root_rad"
    D = (math.cos(th) * sum(c * Mk for c, Mk in zip(x_coef, Dk))
         + math.sin(th) * sum(c * Mk for c, Mk in zip(y_coef, Dk)))
    return D, th, [float(c) for c in x_coef], [float(c) for c in y_coef]


def twist(D, eta):
    """W(eta) = exp(i eta D): D is Hermitian (a real Cartan combination), so
    the exponential is exact through its own eigen-decomposition."""
    ev, Vv = np.linalg.eigh(D)
    return (Vv * np.exp(1j * eta * ev)) @ Vv.conj().T


# --------------------------------------------------------------------------
# M_star and its spectral closure functional
# --------------------------------------------------------------------------


def bond_block(beta, q=(0.0, 0.0, 0.0)):
    """B(q; beta) = Pi_t S_t(q; beta): the 64-dim bond factor of the star
    step (the product of the 12 phase-major stage matrices)."""
    M = np.eye(CELL["nsub"], dtype=complex)
    for (i, phase) in STAGES:
        M = X.stage_matrix(X.module_bonds(CELL, i, phase), q, AXES,
                           0.0, CELL["nsub"], beta) @ M
    return M


class ClosureScan:
    """Spectral closure of M_star = kron(B, (U16 W)^12): the exact
    ||M^k - I||_2 for a unitary M, through the 64+16 eigenvalues."""

    def __init__(self):
        self._mu = {}      # beta -> eigenvalues of B(GAMMA; beta)
        self._nu = {}      # (g, eta) -> eigenvalues of (U16(g) W(eta))^12

    def mus(self, beta):
        if beta not in self._mu:
            self._mu[beta] = np.linalg.eigvals(bond_block(beta))
        return self._mu[beta]

    def nus(self, g, eta, D):
        key = (g, eta)
        if key not in self._nu:
            _B, _U32, U16 = SC.sectors(g)
            V = U16 @ twist(D, eta)
            self._nu[key] = np.linalg.eigvals(V) ** 12
        return self._nu[key]

    def err(self, beta, g, eta, k, D):
        mu = self.mus(beta) ** k
        nu = self.nus(g, eta, D) ** k
        return float(np.abs(mu[:, None] * nu[None, :] - 1.0).max())

    def label(self, beta, g, eta, D, tol=CLOS_TOL):
        for k in range(1, MAX_ORDER + 1):
            if self.err(beta, g, eta, k, D) < tol:
                return k
        return None

    def refine_boundary(self, f, lo, hi, tol=1e-10):
        flo, fhi = bool(f(lo)), bool(f(hi))
        if flo == fhi:
            return None
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if bool(f(mid)) == flo:
                lo = mid
            else:
                hi = mid
        return (lo, hi)

    def refine_point(self, beta0, g0, eta, k, D, tol=1e-10):
        """Minimize the k-th closure error along beta, then along g (the
        V239 descent, spectral form)."""
        b, gg, e = beta0, g0, self.err(beta0, g0, eta, k, D)
        if e < 1e-14:
            return b, gg, e
        for step in (1e-3, 1e-5, 1e-7, 1e-9, 1e-11):
            improved = True
            while improved:
                improved = False
                for db in (-step, step):
                    v = self.err(min(max(b + db, 0.0), math.pi), gg, eta, k, D)
                    if v < e:
                        b, e = b + db, v
                        improved = True
                for dg in (-step, step):
                    v = self.err(b, min(max(gg + dg, 0.0), math.pi), eta, k, D)
                    if v < e:
                        gg, e = gg + dg, v
                        improved = True
            if e < 1e-15:
                break
        return b, gg, e


# --------------------------------------------------------------------------


def direct_closure(beta, g, eta, D, k):
    """Direct 1024-dim check of the spectral claim at one point: build
    M via bloch_star with the twisted map, matrix-power it, report the
    entry-max and spectral norms of M^k - I."""
    _B, _U32, U16 = SC.sectors(g)
    U = U16 @ twist(D, eta)
    M = E.asnumpy(X.bloch_star(CELL, U, g, STAGES, (0.0, 0.0, 0.0), AXES,
                               bond_g10=beta))
    Mk = np.linalg.matrix_power(M, k)
    R = Mk - np.eye(1024, dtype=complex)
    return float(np.abs(R).max()), float(np.linalg.svd(R, compute_uv=False)[0])


def main():
    t0 = time.time()
    checks, data = [], {}
    rng = np.random.default_rng(249)
    D, th_root, x_coef, y_coef = wilson_direction()

    # ---- T1 the factorization identity kron(Pi S_t, U^12) = bloch_star ----
    fac = []
    for s in range(8):
        q = rng.uniform(-1.0, 1.0, 3)
        g = rng.choice([0.2, math.pi / 4, 0.37, 1.1, math.pi / 2])
        beta = rng.choice([0.2, math.pi / 4, 0.63, 1.21, math.pi / 2])
        _B, _U32, U16 = SC.sectors(float(g))
        M_fact = np.kron(bond_block(float(beta), tuple(q)), np.linalg.matrix_power(U16, 12))
        M_ref = E.asnumpy(X.bloch_star(CELL, U16, float(g), STAGES,
                                       tuple(float(x) for x in q), AXES,
                                       bond_g10=float(beta)))
        fac.append({"q": [float(x) for x in q], "beta": float(beta),
                    "g": float(g),
                    "max_abs": float(np.abs(M_fact - M_ref).max())})
    fac_max = max(s["max_abs"] for s in fac)
    checks.append({
        "name": f"T1 the star step factorizes exactly, "
                f"bloch_star(q; beta, g) = kron(Pi_t S_t(q; beta), U16(g)^12): "
                f"max element defect over 8 random q-points "
                f"{fac_max:.2e} < 1e-13 (the mixed-product identity; the "
                f"512-fold flat GAMMA node of V244 is its consequence)",
        "pass": bool(fac_max < 1e-13)})
    data["factorization"] = {"max_defect": fac_max, "samples": fac}

    # ---- T2 the closure census C_star on the (beta, g) grid x eta ladder --
    eta_values = [0.0, 0.05, 0.1, 0.2, th_root]
    scan = ClosureScan()
    dials = np.linspace(LO, HI, N_DIAL)
    labels, closed_grid = {}, []
    for eta in eta_values:
        for g in dials:
            for b in dials:
                k = scan.label(float(b), float(g), eta, D)
                labels[(eta, round(float(b), 6), round(float(g), 6))] = k
                if k is not None:
                    closed_grid.append((float(b), float(g), eta, k))
        # config branch values on both dials (the linspace misses them)
        for vb in BRANCH_ANGLES.values():
            for vg in BRANCH_ANGLES.values():
                k = scan.label(vb, vg, eta, D)
                if k is not None:
                    closed_grid.append((vb, vg, eta, k))
                    labels[(eta, round(vb, 6), round(vg, 6))] = k

    # plateau boundaries along each dial, bisection to 1e-10
    boundaries = []
    for eta in eta_values:
        for gi in range(N_DIAL):
            g = float(dials[gi])
            for bi in range(N_DIAL - 1):
                b1, b2 = float(dials[bi]), float(dials[bi + 1])
                if (labels[(eta, round(b1, 6), round(g, 6))]
                        != labels[(eta, round(b2, 6), round(g, 6))]):
                    br = scan.refine_boundary(
                        lambda x: scan.label(x, g, eta, D) is not None, b1, b2)
                    if br:
                        boundaries.append({"along": "beta", "eta": eta, "g": g,
                                           "beta_lo": br[0], "beta_hi": br[1]})
        for bi in range(N_DIAL):
            b = float(dials[bi])
            for gi in range(N_DIAL - 1):
                g1, g2 = float(dials[gi]), float(dials[gi + 1])
                if (labels[(eta, round(b, 6), round(g1, 6))]
                        != labels[(eta, round(b, 6), round(g2, 6))]):
                    br = scan.refine_boundary(
                        lambda y: scan.label(b, y, eta, D) is not None, g1, g2)
                    if br:
                        boundaries.append({"along": "g", "eta": eta, "beta": b,
                                           "g_lo": br[0], "g_hi": br[1]})

    # refine every closed grid point; C_star
    C_star, seen = [], []
    for (b, g, eta, k) in sorted(closed_grid,
                                 key=lambda t: (t[3], t[2], t[0], t[1])):
        rb, rg, e = scan.refine_point(b, g, eta, k, D)
        if any(abs(rb - s[0]) < 1e-6 and abs(rg - s[1]) < 1e-6
               and abs(eta - s[2]) < 1e-12 for s in seen):
            continue
        seen.append((rb, rg, eta))
        C_star.append({"beta": rb, "g": rg, "eta": eta, "order": int(k),
                       "closure_err": e})
    by_order = {k: sum(1 for v in labels.values() if v == k)
                for k in range(1, MAX_ORDER + 1)
                if any(v == k for v in labels.values())}
    n_none = sum(1 for v in labels.values() if v is None)
    checks.append({
        "name": f"T2 closure census of M_star(GAMMA; beta, g, eta) on the "
                f"{N_DIAL}x{N_DIAL} grid beta,g in [{LO:.2f}, pi/2] plus the "
                f"config branch values, at the eta ladder {eta_values} "
                f"(last entry = the V65/V67 root angle, the theta-coupled "
                f"choice; source recorded): grid points closing at some "
                f"k <= 16: {sum(by_order.values())} (by order {by_order}), "
                f"none at {n_none}; {len(boundaries)} order-plateau "
                f"boundaries refined to 1e-10 by bisection; C_star refines "
                f"to {len(C_star)} point(s), max refined closure error "
                f"{max((c['closure_err'] for c in C_star), default=0.0):.2e} "
                f"(spectral form, verified directly in T4)",
        "pass": bool(all(c["closure_err"] < CLOS_TOL for c in C_star))})
    data["closure"] = {
        "grid": [N_DIAL, float(LO), float(HI)],
        "eta_values": eta_values,
        "eta_source": "0, 0.05, 0.1, 0.2 are grid dials; the fifth is the "
                      "V65/V67 selective-root angle theta_rad = "
                      f"{th_root} (config desert.theta_root_rad) - the "
                      "theta-coupled choice through the committed direction",
        "grid_labels_none": n_none,
        "grid_labels_by_order": {str(k): v for k, v in by_order.items()},
        "boundary_flips": boundaries,
        "C_star": C_star,
    }

    # ---- T3 membership: distance to the two contact clocks ----------------
    def nearest(target):
        cand = [c for c in C_star if c["order"] <= MAX_ORDER]
        if not cand:
            return None
        c = min(cand, key=lambda c: abs(c["g"] - target))
        return {"g": c["g"], "beta": c["beta"], "eta": c["eta"],
                "order": c["order"], "distance": abs(c["g"] - target)}
    n15, n12 = nearest(CLOCK_1_5), nearest(CLOCK_1_2)
    dist_table = {}
    for eta in eta_values:
        cg = [c for c in C_star if c["eta"] == eta]
        dist_table[str(eta)] = {
            "closed_points": len(cg),
            "min_dist_to_0.2": min((abs(c["g"] - CLOCK_1_5) for c in cg),
                                   default=None),
            "min_dist_to_0.5": min((abs(c["g"] - CLOCK_1_2) for c in cg),
                                   default=None)}
    is_open = any(c["order"] <= MAX_ORDER and abs(c["g"] - CLOCK_1_5) < 1e-6
                  for c in C_star) or \
        any(c["order"] <= MAX_ORDER and abs(c["g"] - CLOCK_1_2) < 1e-6
            for c in C_star)
    def fmt(n):
        if n is None:
            return "none (no closed star point exists)"
        return (f"distance {n['distance']:.3e} at (beta, g, eta, k) = "
                f"({n['beta']:.6f}, {n['g']:.6f}, {n['eta']:.6f}, {n['order']})")
    checks.append({
        "name": f"T3 clock contacts: nearest closed star point to g = 0.2 "
                f"(reads 2.0766/0.04): {fmt(n15)}; nearest to g = 0.5 "
                f"(V208 lightlike pin): {fmt(n12)}; STAR-CLOCK-OPEN requires "
                f"some k <= 16 within 1e-6 of a pin",
        "pass": True})
    data["membership"] = {
        "clock_1_5": CLOCK_1_5, "clock_1_2": CLOCK_1_2,
        "nearest_1_5": n15, "nearest_1_2": n12,
        "per_eta": dist_table,
        "open": bool(is_open)}

    # ---- T4 control: the published star step at (pi/4, pi/4, eta = 0) -----
    v214 = json.load(open(exp("V214_star_module_step", "results",
                              "v214_star_module_step.json")))
    pm = v214["data"]["steps"]["phase-major"]
    k_pub = scan.label(math.pi / 4, math.pi / 4, 0.0, D)
    _B, _U32, U16 = SC.sectors(math.pi / 4)
    M_pub = E.asnumpy(X.bloch_star(CELL, U16, math.pi / 4, STAGES,
                                   (0.0, 0.0, 0.0), AXES,
                                   bond_g10=math.pi / 4))
    gap = float(np.min(np.abs(np.angle(np.linalg.eigvals(M_pub)))))
    unit = float(np.abs(M_pub @ M_pub.conj().T
                        - np.eye(1024, dtype=complex)).max())
    d_e, d_s = (direct_closure(math.pi / 4, math.pi / 4, 0.0, D, k_pub)
                if k_pub else (None, None))
    spec = (scan.err(math.pi / 4, math.pi / 4, 0.0, k_pub, D)
            if k_pub else None)
    ctrl_ok = (abs(gap - pm["gap_at_zero"]) < 1e-15
               and abs(unit - pm["unitarity"]) < 1e-13
               and np.abs(np.kron(bond_block(math.pi / 4),
                                  np.linalg.matrix_power(U16, 12))
                          - M_pub).max() < 1e-13)
    checks.append({
        "name": f"T4 control at (beta = g = pi/4, eta = 0): the kron form "
                f"reproduces the published bloch_star to "
                f"{float(np.abs(np.kron(bond_block(math.pi / 4), np.linalg.matrix_power(U16, 12)) - M_pub).max()):.2e}, "
                f"the GAMMA gap {gap:.3e} matches V214's {pm['gap_at_zero']:.3e} "
                f"(|diff| {abs(gap - pm['gap_at_zero']):.1e} < 1e-15), the "
                f"unitarity defect {unit:.3e} matches V214's "
                f"{pm['unitarity']:.3e}; the closure order of the published "
                f"star step (first measurement; V214 recorded no order) is "
                f"k = {k_pub}, spectral err {spec if spec is None else format(spec, '.3e')}, "
                f"direct 1024-dim check entry-max {d_e if d_e is None else format(d_e, '.3e')}, "
                f"spectral {d_s if d_s is None else format(d_s, '.3e')}",
        "pass": bool(ctrl_ok)})
    data["control"] = {
        "beta": math.pi / 4, "g": math.pi / 4, "eta": 0.0,
        "gap_at_gamma": gap, "v214_gap_at_zero": pm["gap_at_zero"],
        "unitarity": unit, "v214_unitarity": pm["unitarity"],
        "closure_order": k_pub, "spectral_err": spec,
        "direct_entry_max": d_e, "direct_spectral": d_s,
        "note": "V214 recorded the published star step at the closed clock "
                "but no closure order; T4 measures it here and reproduces "
                "V214's recorded numbers as the control"}

    # ---- the branch --------------------------------------------------------
    branch = "STAR-CLOCK-OPEN" if is_open else "STAR-CLOCK-QUANTIZED"
    ok = all(c["pass"] for c in checks)
    if is_open:
        status = (
            f"star module closure, STAR-CLOCK-OPEN ({time.time() - t0:.0f} s): "
            f"the icosahedral module's own closure admits the hierarchy pin / "
            f"the lightlike pin: the dead-zone diagnosis was supercell-specific. "
            f"C_star has {len(C_star)} closed point(s) within 1e-6 of "
            f"g = 0.2 (nearest {n15['distance']:.2e}) or g = 0.5 "
            f"(nearest {n12['distance']:.2e}); it is a theory point on C_star, "
            f"no fitting, and V255 uses it.")
    else:
        status = (
            f"star module closure, STAR-CLOCK-QUANTIZED "
            f"({time.time() - t0:.0f} s): no closed star clock reaches either "
            f"pin: the clock exclusion is a property of the module class, not "
            f"of one lattice. C_star: {len(C_star)} closed point(s) on the "
            f"{N_DIAL}x{N_DIAL} grid at the eta ladder, orders "
            f"{sorted(by_order)}; nearest closed point to g = 0.2 at distance "
            f"{None if n15 is None else format(n15['distance'], '.3e')}, to "
            f"g = 0.5 at distance "
            f"{None if n12 is None else format(n12['distance'], '.3e')}; "
            f"the factorization kron identity holds to {fac_max:.1e} and the "
            f"published step reproduces V214 exactly at the closed clock.")
    out = {"id": "V249_star_module_closure_census",
           "tier": 2,
           "property": "the star module's own closure census C_star(beta, g, "
                       "eta): the Cartan twist from the committed Wilson "
                       "direction, and the 1/5 and 1/2 clock contacts (V249)",
           "residual": True,
           "programme": "V249",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    out["data"]["backend"] = E._BACKEND
    out["data"]["fitted_inputs"] = []
    out["data"]["inputs"] = {
        "twist_direction": "D = cos(theta_root) X + sin(theta_root) Y, the "
                           "V67 Cartan combinations with the V65 "
                           f"X_coef {x_coef}, Y_coef {y_coef}, theta_root "
                           f"{th_root} (desert.theta_root_rad); restricted "
                           "to the 16 sector through sectors()[0] exactly "
                           "as V67:64-66 (no fallback needed)",
        "step": "cl10.starcell phase-major 12 stages, bond_g10 = beta",
        "internal": "supercell.sectors(g)[2]"}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v249_star_module_closure_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
