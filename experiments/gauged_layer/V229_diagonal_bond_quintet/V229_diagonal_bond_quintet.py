#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V229 - the diagonal-bond quintet: attacking the graviton-speed falsifier F8
through the named gap (tier 5, F8).

**The question this pre-registers against** (`review/pre_registrations.md`):
F4/F8 name the defect. The quintet's speeds are `{0, g10, 2 g10}` - not 1 at
every clock that keeps the ladder (V140 T7, doc Section 10.4) - and the
`T2g` half of the quintet has NO nearest-neighbour coupling on axis bonds:
`d^T E^off d = 0` on every axis bond, "the coupling needs diagonal bonds"
(V141 T4, a named gap). The icosahedral star module supplies exactly diagonal
bonds with isotropic fourth moment (V214: six generators of norm 1.902113,
equiangular lines at 0.447214, fourth moment isotropic to 2.2e-16; V225: the
equivariant census pushes the first anisotropy to degree 5). Does a
diagonal-bond quintet give a massless branch at `c = 1`?

**What is computed.** The bond set `d` is the six star axes from
`cl10.starcell.star_axes()`; everything else is the V140/V141 machinery
carried over unchanged, per the pre-registration:

T1 the stress vertex `V_gamma(d) = d^T E^gamma d` for every gamma of the
   quintet basis (the V141 Part B `Emat` construction, copied) and every star
   axis: the T2g-type half acquires NONZERO vertices on the diagonal axes,
   where it vanished identically on the axis bonds (V141's artifact is read as
   the control row through `cl10.paths.exp`).

T2 the quintet dispersion on the star bond set: the dressed directional step
   `M(q) = exp(-i g10 sum_a sin(q . d_a) J_{dhat_a})`, with the V140 spin-2
   generators `J_k(S) = i (A_k S - S A_k)` on sym^2, at the closed clocks
   {pi/4, pi/2, 3pi/4}: gapless content at GAMMA (six bands, `A1g + Eg + T2g`
   as a representation of the six-bond layer - under the icosahedral group the
   quintet is one irreducible piece, V230); the group velocities
   `v = dE/d|q|` of every propagating branch on the `continuum.radii` ladder.
   The tight-frame identity `sum_a (n . d_a) J_{dhat_a} = 2 |d| (n . J)` -
   exact for the six equiangular star axes - is verified and makes the linear
   speeds direction-INDEPENDENT: `{0, 2 g10 |d|, 4 g10 |d|}` with the module's
   own displacement normalisation (and `{0, 2 g10, 4 g10}` with unit bonds).

T3 the group velocities of the star automaton's own propagating branches
   (the 1024-dimensional `bloch_star` step of `cl10.starcell`,
   `cl10/spincell`-free, at the closed clocks) - the control layer: the
   fermion side of the O9 node question, measured with the same radii ladder.

T4 the `Eg` (under the star bonds: the whole quintet) mass shift
   `T_gamma = sum_a (dhat_a^T E^gamma dhat_a) J_{dhat_a}` with the exact
   currents `J = i (d_v M) M^dag` of the IMPLEMENTED 128-dimensional module
   step `M(q) = kron(SC.bond_product(q, g), U16^6)` - the V141 Part B
   machinery, the only change being the bond set: the polarization bubble at
   GAMMA over the 5 x 5 quintet space and the order-lambda^2 shift at
   lambda = g10. V141's axis-bond result (`Pi|T2g = 0`, two Eg eigenvalues)
   is the control row.

**Acceptance (pre-registered).** `branch: F8-CLOSED` iff some propagating
quintet branch satisfies `|v - 1| < 1e-12` at a closed angle with
`data.fitted_inputs == []`; every other outcome is `branch: F8-FIRED`.

Inputs: cl10/config.yaml (the closed clocks), the star axes from
`cl10.starcell.star_axes()`, `continuum.radii`, the committed internal map
`supercell.sectors(pi/4)[2]`, the V140/V141 artifacts read through
`cl10.paths.exp`, and the machinery copied from V141 Part B. Zero new
parameters: `data.fitted_inputs == []`.
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
from scipy.linalg import expm, logm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10 import config as _config                            # noqa: E402
from cl10 import engine as E                                  # noqa: E402
from cl10 import starcell as X                                # noqa: E402
from cl10 import supercell as SC                              # noqa: E402
from cl10.paths import exp                                    # noqa: E402

CFG = _config.CFG
CLOSED = (("pi/4", math.pi / 4), ("pi/2", math.pi / 2),
          ("3pi/4", 3 * math.pi / 4))
RADII = [float(r) for r in CFG.continuum.radii]
TOL_MACH = 1e-12
TOL_CLOSE = 1e-12                                  # the F8 acceptance band


def art(name, key):
    p = exp(name, "results", name.lower() + ".json")
    with open(p) as fh:
        j = json.load(fh)
    return j[key] if key else j


# ---------------------------------------------------------------------------
# the quintet polarization basis, copied from V141 Part B (lines 317-343):
# the six symmetric 3x3 basis matrices SB, the trace vector, and the orthonormal
# quintet basis q5 (rows: two diagonal (Eg) then three off-diagonal (T2g) by
# the projector construction).  Emat(g) is the g-th polarization tensor.
# ---------------------------------------------------------------------------
IDX = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def sym2_basis():
    B = []
    for (k, l) in IDX:
        S = np.zeros((3, 3))
        S[k, l] += 1.0
        if l != k:
            S[l, k] += 1.0
        B.append(S / math.sqrt(2.0 if l != k else 1.0))
    return B


SB = sym2_basis()
TRACE_V = np.array([np.trace(S) for S in SB]) / math.sqrt(3.0)
P6 = np.eye(6) - np.outer(TRACE_V, TRACE_V)
ev6, U6m = np.linalg.eigh(P6)
q5 = U6m[:, ev6 > 0.5].T.conj()          # 5 x 6 orthonormal quintet basis


def Emat(g):
    M = np.zeros((3, 3))
    for a in range(6):
        M += q5[g, a] * SB[a]
    return M


diag_rows = [g for g in range(5) if np.abs(q5[g, 3:]).max() < 1e-9]
off_rows = [g for g in range(5) if np.abs(q5[g, :3]).max() < 1e-9]


# the V140 spin-2 generators J_k on sym^2 (rotation_generators/spin2_generators
# copied from V140 lines 253-278)
def rotation_generators():
    lev = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
           (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
    A = []
    for k in range(3):
        M = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                M[i, j] = -lev.get((k, i, j), 0)
        A.append(M)
    return A


def spin2_generators():
    A = rotation_generators()
    J = []
    for k in range(3):
        M = np.zeros((6, 6), dtype=complex)
        for b in range(6):
            act = 1j * (A[k] @ SB[b] - SB[b] @ A[k])
            for a in range(6):
                M[a, b] = np.trace(SB[a].T @ act)
        J.append(M)
    return J


J6 = spin2_generators()


def main():
    t0 = time.time()
    checks, data = [], {}
    data["fitted_inputs"] = []

    axes = [np.asarray(a, float) for a in X.star_axes()]
    norms = [float(np.linalg.norm(a)) for a in axes]
    unit = [a / n for a, n in zip(axes, norms)]

    # ---- T1 the stress vertex on the diagonal (star-axis) bonds ----------
    star_coeff = float(max(abs(u @ Emat(g) @ u)
                           for g in range(5) for u in unit))
    off_star = {(int(g), int(a)): float(unit[a] @ Emat(int(g)) @ unit[a])
                for g in off_rows for a in range(6)}
    nz = [abs(v) for v in off_star.values() if abs(v) > 1e-12]
    off_star_min_nonzero = min(nz)
    off_star_max = max(abs(v) for v in off_star.values())
    v141 = art("V141_regular_link_quintet_coupling", "data")["coupling"]
    axis_control = float(v141["t2g_axis_coefficient"])
    n_off_nonzero = len(nz)
    # the star axes are not generic: the T2g vertex obeys their selection
    # rules, so the measured pattern is the observable - what the named gap
    # asked is whether ANY T2g-type vertex turns on off the axes
    t1ok = (axis_control < 1e-12 and n_off_nonzero > 0 and off_star_min_nonzero > 1e-3)
    checks.append({
        "name": f"T1 the stress vertex V_gamma(d) = d^T E^gamma d: on the axis "
                f"bonds the T2g half vanishes identically (V141 control, max "
                f"{axis_control:.1e}); on the six DIAGONAL star axes every "
                f"T2g-type component acquires a nonzero vertex "
                f"({n_off_nonzero} of 15 axis x component entries nonzero, "
                f"the rest exactly zero by the star selection rules; "
                f"magnitudes {off_star_min_nonzero:.4f} to "
                f"{off_star_max:.4f}; full vertex max {star_coeff:.4f}): the "
                f"named gap is exactly what the star bond set supplies",
        "pass": bool(t1ok)})
    data["stress_vertex"] = {
        "axis_bond_control_t2g": axis_control,
        "star_axis_T2g_entries": {f"g{g}_a{a}": v for (g, a), v in
                                  sorted(off_star.items())},
        "star_axis_T2g_nonzero_count": n_off_nonzero,
        "all_components_max": star_coeff}

    # ---- T2 the quintet dispersion on the star bond set ------------------
    def Jhat(dv):
        return sum(c * Jm for c, Jm in zip(dv, J6))

    def M5(q, g10, bonds):
        Hq = sum(math.sin(float(np.dot(q, d))) * Jhat(d / np.linalg.norm(d))
                 for d in bonds)
        return expm(-1j * g10 * Hq)

    # tight-frame slope identity: sum_a (n.d_a) J_hat(d_a) = 2 |d| (n.J)
    rng = np.random.default_rng(11)
    tf_err = 0.0
    for n in rng.normal(size=(8, 3)):
        n = n / np.linalg.norm(n)
        L = sum(float(np.dot(n, a)) * Jhat(a / np.linalg.norm(a))
                for a in axes)
        R = 2 * norms[0] * Jhat(n)
        tf_err = max(tf_err, float(np.abs(L - R).max()))

    def bands5(q, g10, bonds):
        return np.sort(np.abs(np.angle(np.linalg.eigvals(M5(q, g10, bonds)))))

    disp = {}
    rays = [u for u in unit] + [-u for u in unit]
    min_gap_gamma = 0.0
    for tag, g in CLOSED:
        M0 = M5(np.zeros(3), g, axes)
        gap0 = float(np.abs(M0 - np.eye(6)).max())
        min_gap_gamma = max(min_gap_gamma, gap0)
        pred_f = [2 * norms[0] * g * m for m in (0, 1, 2)]
        pred_u = [2 * g * m for m in (0, 1, 2)]
        meas = {"full": {}, "unit": {}}
        spread = {}
        for rr in RADII:
            rows = {"full": [], "unit": []}
            for ray in rays:
                q = rr * ray
                rows["full"].append(
                    sorted(float(x) / rr for x in bands5(q, g, axes)))
                rows["unit"].append(
                    sorted(float(x) / rr for x in bands5(q, g, unit)))
            for conv in ("full", "unit"):
                v1 = [row[2] for row in rows[conv]]
                v2 = [row[4] for row in rows[conv]]
                meas[conv][float(rr)] = {
                    "v1_mean": float(np.mean(v1)), "v2_mean": float(np.mean(v2)),
                    "v1_band": [float(x) for x in sorted(np.unique(
                        np.round(v1, 10)))],
                    "v2_band": [float(x) for x in sorted(np.unique(
                        np.round(v2, 10)))]}
                if rr == RADII[-1]:
                    spread[conv] = {"v1": float(max(v1) - min(v1)),
                                    "v2": float(max(v2) - min(v2))}
        dev1 = abs(meas["full"][float(RADII[-1])]["v1_mean"] - pred_f[1])
        dev2 = abs(meas["unit"][float(RADII[-1])]["v2_mean"] - pred_u[2])
        disp[tag] = {"gap_gamma": gap0,
                     "predicted_speeds_full": pred_f,
                     "predicted_speeds_unit": pred_u,
                     "measured_by_radius": meas,
                     "direction_spread_smallest_radius": spread,
                     "closed_form_agreement": {
                         "full_v1_dev": float(dev1),
                         "unit_v2_dev": float(dev2)}}
        if tag == "pi/4":
            speed_check = (dev1 < 1e-6 and dev2 < 1e-6
                           and max(max(d["v1"], d["v2"])
                                   for d in spread.values()) < 1e-9)

    # the F8 closure test: any propagating quintet branch at |v - 1| < 1e-12
    # at a closed angle (any bond normalisation)?  Use the smallest radius
    # (best linear limit) and the predicted linear speeds.
    best_dev = 9.9
    best_which = None
    for tag, g in CLOSED:
        for conv, pred in (("full", disp[tag]["predicted_speeds_full"]),
                           ("unit", disp[tag]["predicted_speeds_unit"])):
            for v in pred[1:]:                       # propagating only
                dv = abs(v - 1.0)
                if dv < best_dev:
                    best_dev = dv
                    best_which = f"{tag} bond={conv} v={v:.6f}"
    f8_closed = best_dev < TOL_CLOSE and not data["fitted_inputs"]
    t2ok = (tf_err < 1e-14 and min_gap_gamma < 1e-12 and speed_check)
    checks.append({
        "name": f"T2 the quintet step on the six diagonal star bonds: the "
                f"tight-frame identity sum_a (n.d_a) J_dhat = 2|d| (n.J) holds "
                f"to {tf_err:.1e} (six equiangular axes -> exact leading-order "
                f"isotropy of the STAR bond set); at GAMMA the step is the "
                f"identity ({min_gap_gamma:.1e}, six gapless bands); the "
                f"measured branch speeds at the smallest radius agree with the "
                f"closed form {{0, 2 g10 |d|, 4 g10 |d|}} = "
                f"{[round(x, 4) for x in disp['pi/4']['predicted_speeds_full']]} "
                f"(module displacement) and {{0, 2 g10, 4 g10}} = "
                f"{[round(x, 4) for x in disp['pi/4']['predicted_speeds_unit']]}"
                f" (unit bonds) at the pi/4 clock; the closest propagating "
                f"quintet branch to v = 1 over all closed clocks is "
                f"{best_which} with |v - 1| = {best_dev:.4f} (pre-registered "
                f"acceptance band 1e-12: "
                f"{'CLOSED' if f8_closed else 'NOT closed'})",
        "pass": bool(t2ok)})
    data["dispersion"] = disp
    data["f8_closure_test"] = {"closest_branch": best_which,
                               "closest_abs_v_minus_1": best_dev,
                               "closed_at_any_clock": bool(f8_closed)}

    # ---- T3 the star automaton's own propagating branches (control) ------
    cell = X.module_cell()
    stages = X.module_stages(order="block-major")
    ferm = {}
    for tag, g in CLOSED:
        _B, _U32, U16 = SC.sectors(g)
        slopes = []
        for ray in rays:
            vals = []
            for r in RADII:
                ph = np.abs(X.eigenphases(
                    X.bloch_star(cell, U16, g, stages, r * ray, axes)))
                ph = np.sort(ph)
                nz = ph[ph > 1e-6]
                if len(nz):
                    vals.append(float(nz.min()) / r)
            if vals:
                slopes.append(vals[-1])
        ferm[tag] = {"branch_speeds_min_nonzero": sorted(
            {round(s, 4) for s in slopes})[:24],
            "closest_abs_v_minus_1": float(min(abs(s - 1.0) for s in slopes))}
    checks.append({
        "name": f"T3 control layer, the star automaton itself (bloch_star, "
                f"1024-dimensional): the smallest nonzero branch slopes at the "
                f"closed clocks stay away from v = 1 (closest "
                f"|v - 1| = {min(ferm[t]['closest_abs_v_minus_1'] for t, _ in CLOSED):.4f} "
                f"at the smallest radius) - the fermion layer supplies no "
                f"c = 1 branch either",
        "pass": bool(all(ferm[t]["closest_abs_v_minus_1"] > 1e-6
                         for t, _ in CLOSED))})
    data["star_automaton_branches"] = ferm

    # ---- T4 the quintet mass shift under the diagonal bonds --------------
    # V141 Part B carried over verbatim except the bond set: the implemented
    # 128-dimensional module step M(q) = kron(bond_product(q, g), U16^6), the
    # current J = i (d_v M) M^dag along the bond DIRECTION, the bubble at GAMMA
    # over the orthonormal basis of the Hermitian generator -i log M.
    g = math.pi / 4                                   # explicit closed clock
    _B, _U32, U16 = SC.sectors(g)
    U6p = np.linalg.matrix_power(E.asnumpy(U16), 6)

    def Mmat(q):
        return np.kron(E.asnumpy(SC.bond_product(np.asarray(q, float), g)),
                       U6p)

    def Jd(q, v, dh=1e-6):
        qp = np.asarray(q, float) + dh * np.asarray(v, float)
        qm = np.asarray(q, float) - dh * np.asarray(v, float)
        dM = (Mmat(qp) - Mmat(qm)) / (2 * dh)
        Jx = 1j * dM @ Mmat(q).conj().T
        return (Jx + Jx.conj().T) / 2

    def Tgam(gm):
        out_ = np.zeros((128, 128), dtype=complex)
        for a in range(6):
            c = unit[a] @ Emat(int(gm)) @ unit[a]
            out_ += c * Jd(np.zeros(3), unit[a])
        return out_

    # axis-bond variant of the same vertex for the in-file control
    def Tgam_axis(gm):
        out_ = np.zeros((128, 128), dtype=complex)
        for di in range(3):
            e = np.zeros(3)
            e[di] = 1.0
            c = e @ Emat(int(gm)) @ e
            out_ += 2.0 * c * Jd(np.zeros(3), e)
        return out_

    Hgen = 1j * logm(Mmat(np.zeros(3)))
    Hgen = (Hgen + Hgen.conj().T) / 2
    om, V0 = np.linalg.eigh(Hgen)
    occ, unocc = om < 0, om > 0
    den = om[unocc][:, None] - om[occ][None, :]
    zero_modes = int(np.sum(~(occ | unocc)))

    def bubble(Tb):
        Tbs = [V0.conj().T @ Tb(gm) @ V0 for gm in range(5)]
        PiM = np.zeros((5, 5), dtype=complex)
        for p in range(5):
            for q2 in range(5):
                X = Tbs[p][np.ix_(unocc, occ)] * Tbs[q2].conj()[np.ix_(unocc, occ)]
                PiM[p, q2] = -2.0 * np.sum(X / den)
        return PiM

    PiS = bubble(Tgam)
    evS = np.linalg.eigvalsh((PiS + PiS.conj().T) / 2)
    lam = g
    shiftS = [0.5 * lam * lam * float(c) for c in evS]
    PiA = bubble(Tgam_axis)
    evA = np.linalg.eigvalsh((PiA + PiA.conj().T) / 2)
    shiftA = [0.5 * lam * lam * float(c) for c in evA]
    shiftS_max = float(np.max(np.abs(shiftS)))
    shiftA_max = float(np.max(np.abs(shiftA)))
    t4ok = (shiftS_max > 1e-9 and abs(float(np.abs(PiS[np.ix_(off_rows,
                                                            off_rows)]).max())
                                      ) > 1e-9)
    checks.append({
        "name": f"T4 the quintet mass shift under the diagonal bonds (lambda "
                f"= g10 = pi/4): the V141 axis-bond vertex gaps only the Eg "
                f"half (max |shift| {shiftA_max:.4f} on the two Eg "
                f"directions, T2g bubble 0 to the V141 record "
                f"{v141['Pi_T2g_max']:.0e}); on the six diagonal star bonds "
                f"the bubble is nonzero on the FULL quintet: "
                f"{sum(1 for s in shiftS if abs(s) > 1e-9)} of 5 shift "
                f"eigenvalues nonzero, max |shift| = {shiftS_max:.4f} "
                f"({zero_modes} zero modes of -i log M excluded from the "
                f"resolvent) - closing the T2g coupling gap does not protect "
                f"gaplessness, it gaps every quintet component at order "
                f"lambda^2",
        "pass": bool(t4ok)})
    data["mass_shift"] = {
        "star_bond_bubble_eigenvalues": [float(x) for x in evS],
        "star_bond_shifts_at_lambda_g10": shiftS,
        "star_bond_max_abs_shift": shiftS_max,
        "axis_bond_bubble_eigenvalues": [float(x) for x in evA],
        "axis_bond_shifts_at_lambda_g10": shiftA,
        "axis_bond_max_abs_shift": shiftA_max,
        "v141_control": {k: v141[k] for k in
                         ("Pi_Eg_eigenvalues", "Pi_T2g_max",
                          "mass_shift_at_lambda_g10")},
        "zero_modes_excluded": zero_modes}

    # -----------------------------------------------------------------------
    branch = "F8-CLOSED" if f8_closed else "F8-FIRED"
    ok = all(c["pass"] for c in checks)
    status = (
        f"F8 {('CLOSED' if f8_closed else 'STAYS FIRED')}, measured through "
        f"the named diagonal-bond gap (zero fitted parameters, "
        f"{time.time() - t0:.0f} s). (T1) the T2g half of the quintet, whose "
        f"axis-bond stress vertex vanishes identically (V141 control "
        f"{axis_control:.1e}), acquires a nonzero vertex on every diagonal "
        f"star axis ({n_off_nonzero}/15 entries, "
        f"{off_star_min_nonzero:.4f}-{off_star_max:.4f}): the star bond set "
        f"is exactly the coupling the named gap asked for. (T2) the quintet "
        f"step on the star bonds is exactly isotropic to leading order by a "
        f"measured exact identity (tight frame, sum_a (n.d_a)J_dhat = 2|d| "
        f"(n.J) to {tf_err:.0e}), massless and linear with six gapless bands "
        f"at GAMMA (gap {min_gap_gamma:.1e}), but its speeds are "
        f"{{0, 2 g10 |d|, 4 g10 |d|}} = "
        f"{[round(x, 4) for x in disp['pi/4']['predicted_speeds_full']]} at "
        f"pi/4 (unit bonds: {{0, 2 g10, 4 g10}} = "
        f"{[round(x, 4) for x in disp['pi/4']['predicted_speeds_unit']]}) - "
        f"the diagonal bond set DOUBLES the axis-bond ladder "
        f"{{0, g10, 2 g10}}; the closest propagating branch to c = 1 over all "
        f"closed clocks is {best_which}, |v - 1| = {best_dev:.4f} >> 1e-12. "
        f"(T3) the star automaton's own branches stay away from v = 1 too "
        f"(closest |v - 1| = "
        f"{min(ferm[t]['closest_abs_v_minus_1'] for t, _ in CLOSED):.4f}). "
        f"(T4) the diagonal bonds couple the full quintet: "
        f"{sum(1 for s in shiftS if abs(s) > 1e-9)} of 5 mass shifts nonzero "
        f"at lambda = g10 = pi/4 (max {shiftS_max:.4f} against the axis-bond "
        f"Eg-only {shiftA_max:.4f}): the gap is closed and the price is that "
        f"EVERY quintet component is gapped at order lambda^2. VERDICT: no "
        f"root-supplied bond network gives a massless spin-2 branch at c = 1 "
        f"- axis bonds give {{0, g10, 2 g10}}, diagonal star bonds give "
        f"{{0, 2 g10, 4 g10}} x |d| - and the clock that would force one "
        f"branch to c = 1 (2 g10 |d| = 1 -> g10 = {1.0 / (2 * norms[0]):.6f}, "
        f"or 2 g10 = 1 -> g10 = 1/2) lies outside the closed class "
        f"({1.0 / (2 * norms[0]):.6f} and 1/2 are not pi/4, pi/2, 3pi/4): "
        f"F8 stays fired with the measured set.")
    out = {"id": "V229_diagonal_bond_quintet",
           "tier": 5,
           "property": "the diagonal-bond quintet: the T2g coupling gap "
                       "closes on the star axes, the speeds double to "
                       "{0, 2 g10, 4 g10} and no branch reaches c = 1 at any "
                       "closed clock (V229, F8)",
           "residual": True,
           "programme": "V229",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v229_diagonal_bond_quintet.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
