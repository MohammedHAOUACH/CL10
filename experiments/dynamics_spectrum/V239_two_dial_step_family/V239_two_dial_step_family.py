#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V239 - the two-dial step family and the closure set (tier 2).

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): the published step carries ONE angle - the bond (hop)
angle and the internal Spin(10) angle are the same dial, g10. The repair:
`M(q; beta, g) = prod_t kron(S_t(beta, q_j), U16(g))`, the
`altcell.bloch_stages` shape on the cubic block, with the bond angle beta
decoupled from the internal angle g. Is this family exactly unitary for
arbitrary dials, does it contain the published step as the exact diagonal
slice, and where does the closure axiom `M(GAMMA; beta, g)^k = I` hold?

1. **Unitarity** over the 41x41 dial grid beta, g in linspace(0.02, pi/2)
   at q in {GAMMA, X, M, (0.31, -0.77, 0.44)}: max ||M M^dag - I||_inf.
   The step is a product of two-level rotations and even internal layers,
   unitary at ANY angles - the measurement confirms it, bound 1e-13.
2. **The diagonal slice is exact**: bloch(q, U, g, bond_g10=g) equals
   bloch(q, U) BIT-IDENTICALLY (0.0, `np.array_equal`), and the published
   factorization holds for ANY (beta, g) - the stages and the internal
   layer commute across the kron, so
   `M(q; beta, g) = kron(P(q; beta), U(g)^6)` exactly, measured < 1e-14.
   The diagonal slice beta = g = g10 is the slice where ONE angle drives
   both hop and couplings: the published theory.
3. **The closure set C**: the smallest k in {2,4,8,16} with
   ||M(GAMMA; beta, g)^k - I|| < 1e-12 at every grid point (and at the
   config branch values, which the linspace grid misses), `none`
   otherwise; the order-plateau boundaries refined to 1e-10 by bisection
   along beta at fixed g and along g at fixed beta; the closure points
   themselves refined along both dials. C is stored as ordered
   [{beta, g, order}]. The published closed clock is the point
   (pi/4, pi/4) with order 2 (the artifact fact M(0)^2 = I of
   `supercell.verify`, which here reads P(GAMMA, pi/4)^2 = -I tensor
   U16(pi/4)^12 = -I).

Branch `FAMILY-OK` iff (1) < 1e-13, (2) is exact, and C contains
(pi/4, pi/4) with order 2; otherwise the EXPERIMENT fails (a broken
family is an engineering failure, not a physics verdict).
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
from cl10 import engine as E                                  # noqa: E402
from cl10 import supercell as SC                              # noqa: E402

CFG = _config.CFG
XP = E._E if E._BACKEND == "cupy" else np

N_DIAL = 41
LO, HI = 0.02, math.pi / 2
BRANCH_ANGLES = {  # the config `angle:` branch values: on no linspace grid
    "implemented": float(CFG.angle.branches.implemented),
    "closed": math.pi / 4.0,
    "half": math.pi / 2.0,
}
QS = [SC.GAMMA, SC.X_POINTS[0], SC.M_POINTS[0], (0.31, -0.77, 0.44)]
CLOS_TOL = 1e-12
ORDERS = (2, 4, 8, 16)


def unitarity_grid():
    """max ||M M^dag - I||_inf over the dial grid x QS, M from
    bloch(..., bond_g10=beta), batched on the active backend."""
    dials = np.linspace(LO, HI, N_DIAL)
    best, worst_at = 0.0, None
    for g in dials:
        U16 = SC.sectors(float(g))[2]
        for q in QS:
            Ms = XP.broadcast_to(XP.eye(128, dtype=complex),
                                 (N_DIAL, 128, 128)).copy()
            for even, j in SC.STAGES:
                L = XP.asarray(np.stack(
                    [np.kron(SC.stage_matrix(even, j, q[j], float(g),
                                             float(b)), U16)
                     for b in dials]))
                Ms = L @ Ms
            I128 = XP.eye(128, dtype=complex)
            err = float(np.abs(Ms @ XP.conj(XP.swapaxes(Ms, -1, -2))
                               - I128).max())
            if err > best:
                best = err
                worst_at = (float(g), q)
    return best, worst_at


def diag_controls(rng):
    """(a) bit-identity of bond_g10=g with the published bloch; (b) the
    ANY-(beta,g) factorization M = kron(P(q;beta), U^6) chained as the
    code chains it."""
    a_err = 0.0
    b_err = 0.0
    samples = []
    for _ in range(7):
        g = float(rng.uniform(LO, HI))
        q = tuple(float(v) for v in
                  rng.uniform(-math.pi, math.pi, 3))
        samples.append((g, q))
    samples += [(v, SC.GAMMA) for v in BRANCH_ANGLES.values()]
    for g, q in samples:
        _, _, U16 = SC.sectors(g)
        A = SC.bloch(q, U16, g, SC.STAGES, g)
        B = SC.bloch(q, U16, g)
        if not np.array_equal(A, B):            # EXACT equality required
            a_err = max(a_err, float(np.abs(A - B).max()))
        P = SC.bond_product(q, g, SC.STAGES, g)
        Uc = np.array(U16, dtype=complex)
        U6 = Uc.copy()
        for _ in range(5):
            U6 = U6 @ Uc
        b_err = max(b_err, float(np.abs(np.kron(P, U6) - A).max()))
        # off-diagonal pair too: the factorization holds for any (beta,g)
        beta = float(rng.uniform(LO, HI))
        Pb = SC.bond_product(q, g, SC.STAGES, beta)
        Mb = SC.bloch(q, U16, g, SC.STAGES, beta)
        b_err = max(b_err, float(np.abs(np.kron(Pb, U6) - Mb).max()))
    return a_err, b_err, samples


# --------------------------------------------------------------------------
# the closure machinery (imported by V240 - its label/refine primitives)
# --------------------------------------------------------------------------

def step_at_gamma(beta, g):
    _, _, U16 = SC.sectors(g)
    return SC.bloch(SC.GAMMA, U16, g, SC.STAGES, beta)


def closure_err(beta, g, k):
    M = step_at_gamma(beta, g)
    return float(np.abs(np.linalg.matrix_power(M, k)
                        - np.eye(128, dtype=complex)).max())


def closure_label(beta, g, orders=ORDERS, tol=CLOS_TOL):
    """The smallest k in `orders` closing M(GAMMA; beta, g) within tol, or
    None ('none')."""
    M = step_at_gamma(beta, g)
    Mk = np.eye(128, dtype=complex)
    for k in range(1, max(orders) + 1):
        Mk = Mk @ M
        if k in orders and np.abs(Mk - np.eye(128, dtype=complex)).max() < tol:
            return k
    return None


def refine_boundary(f, lo, hi, tol=1e-10):
    """Bisect the boolean flip of f along [lo, hi]; returns the bracket
    endpoint (1e-10 certified)."""
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


def refine_point(beta0, g0, k, tol=1e-10):
    """Minimize the k-th closure error along beta, then along g, by
    golden-free 1-D grid descent; returns the refined (beta, g, err)."""
    def err(b, gg):
        return closure_err(b, gg, k)
    b, gg, e = beta0, g0, err(beta0, g0)
    if e < 1e-14:                       # float64 round-off floor of the
        return b, gg, e                # 128-dim power: nothing below to find
    for step in (1e-3, 1e-5, 1e-7, 1e-9, 1e-11):
        improved = True
        while improved:
            improved = False
            for db in (-step, step):
                v = err(min(max(b + db, 0.0), math.pi), gg)
                if v < e:
                    b, e = b + db, v
                    improved = True
            for dg in (-step, step):
                v = err(b, min(max(gg + dg, 0.0), math.pi))
                if v < e:
                    gg, e = gg + dg, v
                    improved = True
        if e < 1e-16:
            break
    return b, gg, e


def main():
    t0 = time.time()
    checks, data = [], {}
    rng = np.random.default_rng(239)

    # ---- T1 unitarity over the dial grid -------------------------------
    uerr, worst = unitarity_grid()
    checks.append({
        "name": f"T1 the two-dial step is exactly unitary on the whole "
                f"{N_DIAL}x{N_DIAL} dial grid beta,g in [{LO:.2f}, pi/2] at "
                f"q = GAMMA, X, M and (0.31,-0.77,0.44): max ||M M^dag - I|| "
                f"= {uerr:.2e} < 1e-13 (worst at g = {worst[0]:.4f}, "
                f"q = {worst[1]})",
        "pass": bool(uerr < 1e-13)})
    data["unitarity"] = {"max_err": uerr, "worst_at": [worst[0], list(worst[1])],
                         "grid": [N_DIAL, float(LO), float(HI)],
                         "q_points": [list(q) for q in QS]}

    # ---- T2 the diagonal slice is the published step, exactly -----------
    a_err, b_err, samples = diag_controls(rng)
    checks.append({
        "name": f"T2 the diagonal slice bond_g10=g equals the published "
                f"bloch BIT-IDENTICALLY (max abs diff over 10 samples "
                f"{a_err}, required exactly 0.0 via np.array_equal), and "
                f"M(q; beta, g) = P(q; beta) kron U(g)^6 to {b_err:.2e} "
                f"< 1e-14 (the factorization holds for ANY pair, not only "
                f"the diagonal: the published slice beta = g = g10 is one "
                f"angle driving both hop and couplings)",
        "pass": bool(a_err == 0.0 and b_err < 1e-14)})
    data["diagonal_slice"] = {"bit_identity_max_abs": a_err,
                              "factorization_max_abs": b_err,
                              "samples": [[g, list(q)] for g, q in samples]}

    # ---- T3 the closure set C ------------------------------------------
    dials = np.linspace(LO, HI, N_DIAL)
    labels = {}
    closure_pts = []
    for g in dials:
        for b in dials:
            k = closure_label(float(b), float(g))
            labels[(round(float(b), 6), round(float(g), 6))] = k
            if k is not None:
                closure_pts.append((float(b), float(g), k))
    # the config branch values on both dials (the linspace misses them)
    extra_pts = []
    for vb in BRANCH_ANGLES.values():
        for vg in BRANCH_ANGLES.values():
            k = closure_label(vb, vg)
            if k is not None:
                extra_pts.append((vb, vg, k))
    # the published closed clock itself
    kp = closure_label(math.pi / 4, math.pi / 4)
    all_pts = closure_pts + extra_pts
    if kp is not None:
        all_pts.append((math.pi / 4, math.pi / 4, kp))

    n_none = sum(1 for v in labels.values() if v is None)
    by_order = {k: sum(1 for v in labels.values() if v == k)
                for k in ORDERS}

    # boundary refinement: sign flips of each order's closure predicate
    boundaries = []
    for gi in range(N_DIAL):
        for bi in range(N_DIAL - 1):
            b1, b2 = float(dials[bi]), float(dials[bi + 1])
            g = float(dials[gi])
            if labels[(round(b1, 6), round(g, 6))] != \
               labels[(round(b2, 6), round(g, 6))]:
                br = refine_boundary(
                    lambda x: closure_label(x, g) is not None, b1, b2)
                if br:
                    boundaries.append({"along": "beta", "g": g,
                                       "beta_lo": br[0], "beta_hi": br[1]})
    for bi in range(N_DIAL):
        for gi in range(N_DIAL - 1):
            g1, g2 = float(dials[gi]), float(dials[gi + 1])
            b = float(dials[bi])
            if labels[(round(b, 6), round(g1, 6))] != \
               labels[(round(b, 6), round(g2, 6))]:
                br = refine_boundary(
                    lambda y: closure_label(b, y) is not None, g1, g2)
                if br:
                    boundaries.append({"along": "g", "beta": b,
                                       "g_lo": br[0], "g_hi": br[1]})

    # refine each closure point found on the grid to its exact place
    C = []
    seen = []
    for (b, g, k) in sorted(all_pts, key=lambda t: (t[2], t[0], t[1])):
        rb, rg, e = refine_point(b, g, k)
        if any(abs(rb - s[0]) < 1e-6 and abs(rg - s[1]) < 1e-6 for s in seen):
            continue
        seen.append((rb, rg))
        C.append({"beta": rb, "g": rg, "order": int(k),
                  "closure_err": e})
    data["closure"] = {
        "grid_labels_none": n_none,
        "grid_labels_by_order": {str(k): by_order[k] for k in ORDERS},
        "boundary_flips": boundaries,
        "closure_set": C,
        "published_point": {"beta": math.pi / 4, "g": math.pi / 4,
                            "order": kp}}

    fam_ok = (kp == 2 and
              any(abs(c["beta"] - math.pi / 4) < 1e-6
                  and abs(c["g"] - math.pi / 4) < 1e-6
                  and c["order"] == 2 for c in C))
    checks.append({
        "name": f"T3 the closure axiom M(GAMMA; beta, g)^k = I measured on "
                f"the {N_DIAL}x{N_DIAL} grid and at the config branch "
                f"values: grid points closing at some k <= 16: "
                f"{sum(by_order.values())} (by order "
                f"{ {k: by_order[k] for k in ORDERS} }), none at "
                f"{n_none}; {len(boundaries)} order-plateau boundaries "
                f"refined to 1e-10 by bisection; the closure set refines "
                f"to {len(C)} point(s), the published closed clock "
                f"(pi/4, pi/4) closes at k = {kp} (the verify artifact "
                f"fact M(0)^2 = I reads P^2 = -I tensor U^12 = -I)",
        "pass": bool(fam_ok)})
    checks.append({
        "name": f"T4 C non-empty: {len(C)} refined closure point(s) with "
                f"max refined closure error "
                f"{max((c['closure_err'] for c in C), default=1.0):.2e}",
        "pass": bool(len(C) > 0)})

    ok = all(c["pass"] for c in checks)
    branch = "FAMILY-OK" if ok else "FAIL"
    status = (
        f"TWO-DIAL STEP FAMILY, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the step with the hop decoupled from "
        f"the couplings, M(q; beta, g) = kron(P(q; beta), U16(g)^6), is "
        f"exactly unitary for ANY angle pair (max defect {uerr:.1e}); the "
        f"published step is the diagonal slice beta = g = g10 and that "
        f"sentence is literal - the slice reproduces the published "
        f"matrices bit-for-bit (0.0). The closure axiom "
        f"M(GAMMA; beta, g)^k = I is no longer a single clock: it is the "
        f"measured set C with {len(C)} closed points (orders "
        f"{sorted({c['order'] for c in C})}), and it contains the "
        f"published closed clock (pi/4, pi/4) at order 2.")
    out = {"id": "V239_two_dial_step_family",
           "tier": 2,
           "property": "the two-dial step family: exact unitarity, the "
                       "published step as the exact diagonal slice, and "
                       "the closure set C(beta, g) (V239)",
           "residual": False,
           "programme": "V239",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    out["data"]["backend"] = E._BACKEND
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v239_two_dial_step_family.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
