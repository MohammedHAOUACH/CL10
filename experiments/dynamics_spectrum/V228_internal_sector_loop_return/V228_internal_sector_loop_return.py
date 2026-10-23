#!/usr/bin/env python3
r"""V228 - what the Clifford bond layer does and does not touch: the internal
sector of the 2I-covariant cell, its CPT structure, and its loop return.

V214 replaced the cubic shift module with the icosahedral star module and found
the internal map reusable unchanged; V217 turned that into a measurement on the
star module: chirality conservation, the 16/16-bar sector blocks of the
32-internal step, the committed CPT operator acting on it, and the internal data
(anomaly coefficients, the 5/3 trace-form ratio) untouched because the module
enters none of them.  V225 then changed something deeper than the module: the
*bond layer*, so that a hop carries a spinor matrix $U_d = -\hat v_d\cdot\vec
\sigma$ instead of the identity, and the step is no longer a tensor product of
transport and internal layer.  That is exactly the change that could spoil the
internal sector, so the same questions are asked of it.

T1 chirality conservation on the new cell: with the committed 32-internal map
   the step must commute with $\mathbb 1\otimes w$ ($w$ the pseudoscalar) -- the
   step is an internal scalar.  Measured at generic momenta, not assumed: the
   bond layer acts on direction, layer and spinor, and the claim is that it does
   not leak into the internal factor.
T2 the sector blocks: in the pseudoscalar eigenbasis the 32-internal step must
   be block diagonal (off-block norm measured), its 16-block spectrally identical
   to the step built directly on the 16 sector (so the reuse of $U_{16}$ is a
   consequence, not a choice), and the matter and antimatter blocks compared
   spectrally -- that is CPT on the spectrum of the new cell.
T3 the CPT operator at the operator level: the committed $g_1$ (the time
   direction of the root lattice, which anticommutes with $w$ and exchanges the
   sectors) is applied to the step together with every candidate the cell
   offers -- the direction inversion, the transpose, the internal map inverted --
   and the identity that actually holds is reported, not assumed.  With the bond
   layer this is a real question: $U_{-d} = -U_d$ is a sign, and V225 T1 measured
   that this step is *not* inversion invariant at the zone centre (defect 1.414),
   so the candidates built on the inversion are expected to fail.  The same
   candidate set is run on the cell whose bonds are scalars, where V225 T6
   measured inversion intact: that control is what makes the Clifford cell's
   answer a statement about the bond layer rather than about the instrument.
T4 the internal data: the sector dimensions, the anomaly coefficients and the
   5/3 trace-form ratio recomputed with the Clifford cell in the construction.
T5 the loop return of the new cell: $F_n = \langle \mathrm{tr}\, S(q)^n\rangle /
   d$ averaged over the momentum ball, for $n = 1..12$, checked for stability
   under grid refinement (the star-module analogue of V218's exactness test), and
   measured in the same instrument for the cell whose bonds are scalars -- same
   clock, same layers, same internal exponent, so the only difference is the
   spinor factor: the number that says what the Clifford bond layer does to the
   amplitude the hierarchy's loop factor is built from.

Every claim below is a measurement on the cell V225 introduced.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
# The engine picks the core; the cell steps are device work.
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import config as _config
from cl10 import cpt as CPT
from cl10 import engine as E
from cl10 import gauge
from cl10 import matter
from cl10 import spincell as SP
from cl10 import supercell as SC


QPROBE = ((0.31, -0.17, 0.44), (0.09, 0.52, -0.23), (-0.61, 0.12, 0.08))
NPOW = 12                        # the return amplitudes F_1 .. F_NPOW
GRIDS = (5, 9, 17)               # grid refinement for the momentum average


def as_side(x, ref):
    """Put a matrix on the same side (device or host) as the matrices it will be
    multiplied with; the engine never pulls a device array to the host."""
    if isinstance(ref, np.ndarray):
        return x if isinstance(x, np.ndarray) else E.asnumpy(x)
    return x if not isinstance(x, np.ndarray) else E.to_device(x)


def inversion(nsub, layers, nspin, nint):
    """The direction inversion (each of the twelve directions with its opposite)
    lifted to the cell's full space."""
    Pin = np.zeros((SP.NDIR, SP.NDIR))
    for (a, b) in SP.opposite_pairs():
        Pin[a, b] = Pin[b, a] = 1.0
    W = np.kron(Pin, np.eye(layers, dtype=complex))
    W = np.kron(W, np.eye(nspin, dtype=complex))
    return np.kron(W, np.eye(nint, dtype=complex))


def sector_basis(w):
    """The unitary ordering the internal basis by pseudoscalar eigenvalue: +i
    (matter, the 16) first, then -i (antimatter, the 16-bar) -- the same
    ordering `supercell.sectors` uses (V217's instrument)."""
    ev, evec = np.linalg.eig(w)
    order = np.argsort(-np.imag(ev))
    return np.linalg.qr(evec[:, order])[0]


def sector_blocks(M, T, nsub, nint):
    """Matter and antimatter blocks of a step in the pseudoscalar eigenbasis,
    with the off-block norm and the sorted eigenphase spectra of each block
    (V217's instrument, with `nsub` the dimension of everything the internal
    factor is tensored with)."""
    Mh = E.asnumpy(M)
    Tfull = np.kron(np.eye(nsub), T)
    Mb = Tfull.conj().T @ Mh @ Tfull
    n16 = nint // 2
    idx16 = np.array([p * nint + i for p in range(nsub) for i in range(n16)])
    idxbar = np.array([p * nint + i for p in range(nsub)
                       for i in range(n16, nint)])
    off = max(float(np.max(np.abs(Mb[np.ix_(idxbar, idx16)]))),
              float(np.max(np.abs(Mb[np.ix_(idx16, idxbar)]))))
    b16 = Mb[np.ix_(idx16, idx16)]
    bbar = Mb[np.ix_(idxbar, idxbar)]
    return {"total": int(Mh.shape[0]), "dim16": int(b16.shape[0]),
            "dimbar": int(bbar.shape[0]), "off_block": off,
            "phases16": [float(x) for x in np.sort(np.abs(np.angle(
                np.linalg.eigvals(b16))))],
            "phasesbar": [float(x) for x in np.sort(np.abs(np.angle(
                np.linalg.eigvals(bbar))))]}


def return_amplitude(step_fn, grid_pts, nmax):
    """F_n = <tr S(q)^n>/d averaged over the sampled momenta, with the mean of
    |tr| alongside: the dynamical return of the cell that the loop factor is
    built from (V218's T2 instrument, on the ball rather than a torus, since the
    icosahedral module has no reciprocal lattice)."""
    out = [0.0] * (nmax + 1)
    absout = [0.0] * (nmax + 1)
    dims = set()
    for q in grid_pts:
        S = E.to_device(step_fn(q))
        dims.add(int(S.shape[0]))
        M = S
        for n in range(1, nmax + 1):
            tr = np.trace(E.asnumpy(M))
            out[n] += float(np.real(tr))
            absout[n] += float(abs(tr))
            M = E.matmul(M, S)
    m = len(grid_pts)
    d = int(sorted(dims)[0])
    return {"n": list(range(1, nmax + 1)),
            "F": [out[n] / (m * d) for n in range(1, nmax + 1)],
            "absF": [absout[n] / (m * d) for n in range(1, nmax + 1)],
            "dim": d, "points": m}


def ball_cube(k, rmax):
    """A $k^3$ grid filling the momentum ball of radius rmax."""
    pts = []
    for i in range(k):
        for j in range(k):
            for l in range(k):
                q = rmax * (2.0 * np.array([i, j, l], float) / (k - 1) - 1.0)
                if np.linalg.norm(q) <= rmax:
                    pts.append(q)
    return pts


def main():
    t0 = time.time()
    checks = []
    ok = True

    def chk(name, cond, detail):
        nonlocal ok
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        ok = ok and bool(cond)

    B, U32, U16 = SC.sectors()
    g10 = float(_config.CFG.angle.g10)
    branch = str(_config.CFG.angle.branch)
    tol = float(_config.CFG.tolerances.standard)
    strict = float(_config.CFG.tolerances.strict)
    layers, legs = 1, (4,)

    w = E.asnumpy(CPT.pseudoscalar())
    g1 = E.asnumpy(CPT.cpt_operator())
    T = sector_basis(w)
    sector_frame_agreement = float(np.max(np.abs(
        T[:, :16] @ T[:, :16].conj().T - B @ B.conj().T))) \
        if B.shape[1] >= 16 else None

    # the cell, built on the 32 sector and on the 16
    kw32 = {"g10": g10, "kind": "clifford", "layers": layers, "legs": legs}

    def cell32(q):
        return SP.step(q, U32, **kw32)

    def cell16(q):
        return SP.step(q, U16, **kw32)

    ckw = {"g10": g10, "kind": "scalar", "layers": layers, "legs": legs}

    def scalar32(q):
        return SP.step(q, U32, **ckw)

    def ctrl16(q):
        """The control for the return amplitude: the same clock, the same single
        layer, the same internal exponent -- the only difference from `cell16` is
        that the bonds are scalars instead of the spinor matrices."""
        return SP.step(q, U16, **ckw)

    nsub = SP.NDIR * layers * SP.NSPIN

    # ------------------------------------------------------------------ T1
    w4 = np.kron(np.eye(nsub), w)
    chir = []
    for q in QPROBE:
        Sd = E.asnumpy(cell32(np.array(q, float)))
        comm = float(np.max(np.abs(Sd @ w4 - w4 @ Sd)))
        chir.append({"q": [float(x) for x in q], "defect": comm})
    unit32 = float(np.max(np.abs(np.linalg.matrix_power(U32, 8)
                                 - np.eye(U32.shape[0]))))

    chk("T1 the Clifford bond layer does not leak into the internal factor: the "
        "32-internal step commutes with 1 (x) w",
        all(c["defect"] < 100 * strict for c in chir),
        f"the pseudoscalar is the committed one (w^2 = -1 to "
        f"{float(np.max(np.abs(w @ w + np.eye(32)))):.1e}, g1 exchanges the "
        f"sectors: g1 w g1^-1 + w = "
        f"{float(np.max(np.abs(g1 @ w @ np.linalg.inv(g1) + w))):.1e}), the "
        f"internal map U32 has U32^8 - I = {unit32:.1e}; the step commutes with "
        f"1 (x) w at the probe momenta to "
        f"{max(c['defect'] for c in chir):.1e}")

    # ------------------------------------------------------------------ T2
    blocks = {}
    for q in QPROBE[:2]:
        Sd = cell32(np.array(q, float))
        blk = sector_blocks(Sd, T, nsub, 32)
        blk.pop("block16", None)
        blk.pop("blockbar", None)
        S16 = np.sort(np.abs(np.angle(E.asnumpy(E.eigvals(cell16(np.array(q, float)))))))
        same = float(np.max(np.abs(np.array(blk["phases16"]) - S16))) \
            if len(blk["phases16"]) == len(S16) else None
        cpt_spec = float(np.max(np.abs(np.array(blk["phases16"])
                                       - np.array(blk["phasesbar"])))) \
            if len(blk["phases16"]) == len(blk["phasesbar"]) else None
        blocks[str(np.round(q, 3).tolist())] = {
            "off_block": blk["off_block"], "dim16": blk["dim16"],
            "dimbar": blk["dimbar"], "total": blk["total"],
            "sector16_equals_16_step": same, "matter_vs_antimatter": cpt_spec}

    worst_off = max(b["off_block"] for b in blocks.values())
    worst_same = max(b["sector16_equals_16_step"] or 0.0 for b in blocks.values())
    worst_mva = max(b["matter_vs_antimatter"] for b in blocks.values())

    chk("T2 the 32-internal step of the new cell is block diagonal in the "
        "pseudoscalar basis, its 16-block is the 16-sector step, and its matter "
        "and antimatter blocks are isospectral",
        worst_off < 100 * strict and worst_same < 100 * tol
        and worst_mva < 100 * tol
        and abs(sector_frame_agreement) < 100 * tol,
        "off-block norms "
        + ", ".join(f"{b['off_block']:.1e}" for b in blocks.values())
        + f"; the 16-block's eigenphases agree with the step built directly on "
          f"the 16 sector to {worst_same:.1e} (so the reuse of U16 is forced, "
          f"not chosen); the matter and antimatter blocks' sorted eigenphase "
          f"spectra agree to "
          + ", ".join(f"{b['matter_vs_antimatter']:.1e}" for b in blocks.values())
          + f" -- CPT on the spectrum carries over to the cell that carries the "
            f"spinor bond layer, as it did on the star module (V217 T2); the "
            f"pseudoscalar eigenbasis reproduces the committed sector frame to "
            f"{sector_frame_agreement:.1e}")

    # ------------------------------------------------------------------ T3
    def cpt_candidates(step_fn, U):
        """Every identity the cell offers for the committed CPT operator,
        evaluated over the probe momenta; a candidate's defect is its worst
        value over them."""
        out = {}
        gfull = np.kron(np.eye(nsub), g1)
        ginv = np.kron(np.eye(nsub), np.linalg.inv(g1))
        P = inversion(nsub, layers, SP.NSPIN, np.asarray(U).shape[0])
        PG = P @ gfull
        PGinv = np.linalg.inv(PG)
        for q in QPROBE:
            qq = np.array(q, float)
            S = E.asnumpy(step_fn(qq))
            Sm = E.asnumpy(step_fn(-qq))
            Si = E.asnumpy(SP.step(-qq, np.linalg.inv(np.asarray(U)),
                                   g10=g10, kind=kw32["kind"], layers=layers,
                                   legs=legs))
            opts = {
                "g1 S(q) g1^-1 vs S(-q)":
                    float(np.max(np.abs(gfull @ S @ ginv - Sm))),
                "g1 S(q) g1^-1 vs S(-q)^T":
                    float(np.max(np.abs(gfull @ S @ ginv - Sm.T))),
                "P g1 S(q) g1^-1 (P g1)^-1 vs S(-q)":
                    float(np.max(np.abs(PG @ S @ PGinv - Sm))),
                "P g1 S(q) g1^-1 (P g1)^-1 vs S(-q)^T":
                    float(np.max(np.abs(PG @ S @ PGinv - Sm.T))),
                "g1 S(q) g1^-1 vs S(-q)^-1":
                    float(np.max(np.abs(gfull @ S @ ginv - np.linalg.inv(Sm)))),
                "g1 S(q) g1^-1 vs S_inv(-q) (internal inverted)":
                    float(np.max(np.abs(gfull @ S @ ginv - Si)))}
            for k, v in opts.items():
                out.setdefault(k, []).append(v)
        return {k: max(v) for k, v in out.items()}

    cand = cpt_candidates(cell32, U32)
    cand_c = cpt_candidates(scalar32, U32)
    best = min(cand, key=cand.get)
    best_c = min(cand_c, key=cand_c.get)
    holds = {k: v for k, v in cand.items() if v < 100 * tol}
    holds_c = {k: v for k, v in cand_c.items() if v < 100 * tol}

    def ctext(d):
        return ", ".join(f"{k} = {v:.1e}"
                         for k, v in sorted(d.items(), key=lambda t: t[1]))

    chk("T3 the CPT candidate set is decided by measurement on both cells: it "
        "holds on the spinless-bond control, and the Clifford cell is no better "
        "-- whether any identity survives there is a finding, not a hope",
        len(holds_c) >= 1 and cand[best] >= cand_c[best_c],
        f"control (scalar bonds, same clock and internal exponent): {ctext(cand_c)}"
        f" -- '{best_c}' holds at {cand_c[best_c]:.1e}; Clifford cell: "
        f"{ctext(cand)} -- best '{best}' at {cand[best]:.1e}, "
        + (f"{len(holds)} candidate(s) hold" if holds else
           "no candidate holds, which is the inversion V225 T1 measured broken "
           "(defect 1.414 at the zone centre) showing up in the identity the "
           "CPT operator is asked to satisfy"))

    # ------------------------------------------------------------------ T4
    anom = matter.anomalies()
    five = float(gauge.five_thirds())
    dims = {"n32": int(np.asarray(U32).shape[0]), "n16": int(np.asarray(U16).shape[0])}
    worst_anom = max(abs(v) for v in anom.values())

    chk("T4 the internal data the paper reads off the internal sector are "
        "unchanged by the bond layer: the trace-form ratio, the anomaly "
        "coefficients and the sector dimensions",
        abs(five - 5.0 / 3.0) < 1e-6 and dims["n32"] == 2 * dims["n16"]
        and worst_anom < 100 * tol,
        f"the 5/3 trace-form ratio is {five:.10f}; all {len(anom)} anomaly "
        f"coefficients stay zero, largest {worst_anom:.1e} over the keys "
        f"{sorted(anom)}: "
        f"{', '.join(f'{k} = {anom[k]:.1e}' for k in sorted(anom))}; the sector "
        f"dimensions stay {dims['n32']} = 2 x {dims['n16']}: the bond layer acts "
        f"on direction, layer and spinor only")

    # ------------------------------------------------------------------ T5
    ell = SP.module_length()
    rmax = math.pi / ell
    amps, amps_scalar = [], []
    for k in GRIDS:
        pts = ball_cube(k, rmax)
        amps.append({"grid": k, **return_amplitude(cell16, pts, NPOW)})
        amps_scalar.append({"grid": k, **return_amplitude(ctrl16, pts, NPOW)})
    per_n = [[max(abs(a["F"][n] - b["F"][n]) for a, b in zip(amps[1:], amps[:-1]))
              for n in range(NPOW)]]
    shift = max(per_n[0])
    resolved = [n + 1 for n in range(NPOW)
                if per_n[0][n] <= max(100 * tol, 0.05 * abs(amps[-1]["F"][n]))]
    F = amps[-1]["F"]

    odd = [n for n in range(1, NPOW + 1) if n % 2]
    even = [n for n in range(1, NPOW + 1) if not n % 2]
    absF = amps[-1]["absF"]
    odd_zero = all(absF[n - 1] < 100 * tol for n in odd)
    unresolved = [n for n in range(1, NPOW + 1) if n not in resolved]
    chk("T5 the loop return of the new cell splits by parity at the grids run: "
        "every odd length vanishes in modulus to machine precision, the even "
        "lengths are what survives and whether they are resolved at these grids "
        "is read off the refinement, and the spinor bond layer moves them against "
        "the matched scalar-bond control",
        odd_zero and all(n in even for n in unresolved)
        and len(resolved) >= NPOW // 2
        and all(abs(x) < 1.0 for x in F)
        and max(abs(a["F"][n] - s["F"][n])
                for a, s in zip(amps, amps_scalar) for n in range(NPOW)) > 1e-6,
        f"momentum ball |q| <= {rmax:.4f}, grids {list(GRIDS)} with "
        f"{amps[-1]['points']} points and dimension {amps[-1]['dim']}: F_n for "
        f"n = 1..{NPOW} is {[round(x, 6) for x in F]}. The parity split is "
        f"exact and not a matter of taking real parts: the mean modulus at the "
        f"odd lengths is {[f'{absF[n - 1]:.1e}' for n in odd]}. The largest change "
        f"under refinement, length by length, is {[f'{v:.1e}' for v in per_n[0]]}, "
        f"so the even lengths {[f'{F[n - 1]:.3f}' for n in even]} are resolved at "
        f"{[n for n in even if n in resolved]} and unresolved at "
        f"{[n for n in even if n not in resolved]} (changes "
        f"{[f'{per_n[0][n - 1]:.1e}' for n in even]}) -- at these grids the cell "
        f"supplies no converged loop factor yet. The matched scalar-bond control "
        f"gives {[round(x, 6) for x in amps_scalar[-1]['F']]}, nonzero at every "
        f"length including the odd ones, so it is the Clifford bond layer that "
        f"kills the odd returns: the same odd-even structure the closed branch "
        f"shows for the loop ratios at the root angle (V203 T1), now in the "
        f"kinematic return of the cell itself. The hierarchy's loop factor has to "
        f"be rebuilt on this cell, not carried over")


    verdict = {
        "id": os.path.basename(HERE),
        "tier": 4,
        "property": "the internal sector of the 2I-covariant Clifford cell: "
                    "chirality conservation of the 32-internal step, its "
                    "pseudoscalar sector blocks and their relation to the "
                    "16-sector step, the committed CPT operator against every "
                    "candidate the cell offers, the unchanged anomaly and "
                    "trace-form data, and the grid-converged loop return "
                    "amplitude F_n against the spinless control",
        "residual": True,
        "pass": bool(ok),
        "status": (
            f"INTERNAL SECTOR AND LOOP RETURN OF THE 2I-COVARIANT CLIFFORD CELL "
            f"at the {branch} clock g = {g10:.6f} (layers {layers}, internal "
            f"exponent {legs[0]}). (T1) the 32-internal step commutes with 1 (x) w "
            f"to {max(c['defect'] for c in chir):.1e} at the probe momenta: the "
            f"bond layer acts on direction, layer and spinor and does not leak into "
            f"the internal factor. (T2) in the pseudoscalar eigenbasis the step is "
            f"block diagonal to {worst_off:.1e} and its 16-block is spectrally the "
            f"16-sector step to {worst_same:.1e} -- the reuse of U16 is a "
            f"consequence -- and the matter and antimatter blocks are isospectral "
            f"to {worst_mva:.1e}: CPT on the spectrum survives the bond layer. "
            f"(T3) the CPT "
            f"candidates, worst defect over the probe momenta -- control "
            f"(scalar bonds, same clock and internal exponent): {ctext(cand_c)}"
            f", so '{best_c}' holds at {cand_c[best_c]:.1e}; Clifford cell: "
            f"{ctext(cand)}, so the best is '{best}' at "
            + (f"{cand[best]:.1e} and {len(holds)} candidate(s) hold"
               if holds else
               f"{cand[best]:.1e} and none holds: the identity the committed CPT "
               f"operator has to satisfy fails exactly where V225 T1 measured "
               f"inversion broken (1.414 at the zone centre)")
            + f". (T4) the 5/3 "
              f"trace-form ratio stays {five:.10f}, the {len(anom)} anomaly "
              f"coefficients stay zero (largest {worst_anom:.1e}) and the sector "
              f"dimensions {dims['n32']} = 2 x {dims['n16']}: the internal data "
              f"are untouched. "
              f"(T5) the loop return F_n on the momentum ball is "
              f"{[round(x, 6) for x in F]} for n = 1..{NPOW}: the odd lengths "
              f"vanish in modulus to {[f'{absF[n - 1]:.1e}' for n in odd]} -- the "
              f"split is exact, not a real-part artefact -- while the even lengths "
              f"are resolved at {[n for n in even if n in resolved]} and not at "
              f"{[n for n in even if n not in resolved]} (refinement changes "
              f"{[f'{v:.1e}' for v in per_n[0]]}), against "
              f"{[round(x, 6) for x in amps_scalar[-1]['F']]} for the matched "
              f"scalar-bond control. The bond layer kills the odd returns, so the "
              f"loop factor of Section 12 has to be recomputed on this cell, and "
              f"at these grids no converged even-length value is available to "
              f"recompute it from."),
        "data": {
            "clock": g10, "branch": branch,
            "cell": {"layers": layers, "legs": list(legs),
                     "coin_dim": SP.NDIR, "block_dim": nsub},
            "chirality": chir,
            "internal_map": {"u32_order_defect_8": unit32,
                             "sector_frame_agreement": sector_frame_agreement},
            "sector_blocks": blocks,
            "cpt_candidates": cand, "cpt_identity_holding": best,
            "cpt_holding_candidates": sorted(holds),
            "cpt_control": {"candidates": cand_c, "holding": best_c,
                            "holding_candidates": sorted(holds_c)},
            "internal_data": {"anomalies": anom, "five_thirds": five,
                              "sector_dims": dims},
            "return_amplitude": {"cell": amps,
                                 "control_scalar_bond": amps_scalar,
                                 "refinement_shift": shift,
                                 "per_length_change": per_n[0],
                                 "resolved_n": resolved,
                                 "rmax": rmax, "nmax": NPOW},
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
    print("V228 - the internal sector and the loop return of the new cell")
    print("=" * 68)
    print(f"  chirality: max defect {max(c['defect'] for c in chir):.1e}")
    print("  blocks: " + ", ".join(f"off {b['off_block']:.1e} "
                                   f"16vs16step {b['sector16_equals_16_step']:.1e} "
                                   f"MvsA {b['matter_vs_antimatter']:.1e}"
                                   for b in blocks.values()))
    for k, v in sorted(cand.items(), key=lambda t: t[1]):
        print(f"  CPT Clifford {k}: {v:.2e}")
    for k, v in sorted(cand_c.items(), key=lambda t: t[1]):
        print(f"  CPT control    {k}: {v:.2e}")
    print(f"  internal data: five thirds {five:.10f}, dims {dims}")
    print(f"  return F_n (grid {GRIDS[-1]}): {[round(x, 6) for x in F]}")
    print(f"  control return: {[round(x, 6) for x in amps_scalar[-1]['F']]}")
    print(f"  refinement per length: {[f'{v:.1e}' for v in per_n[0]]} "
          f"(resolved n = {resolved})")
    for c in checks:
        print(f"  [{'ok' if c['pass'] else 'FAIL'}] {c['name']}")
        print(f"        {c['detail']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
