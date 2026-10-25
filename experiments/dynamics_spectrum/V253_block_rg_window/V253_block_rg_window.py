#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V253 - the block-RG window on the dressed sector (tier 2,
residual: true, O9/F9).

**The question (pre-registered).**  V213 measured the bare tangent
anisotropy marginal (exponent -0.0001); V177's static dressing did not
change it.  Does the dressed sector flow?

**The object (closed input set).**  V178's downfold ladder machinery
(altcell.lowdin_downfold, QSTAR = (pi/2, pi/2, 0), DQ = 1e-3, the 12
sphere_dirs(6) rays) applied to the V251 Sigma-corrected effective step,
cuts b = 1, 2, 4, 8, 16, i.e. cut = pi/(2b) - the energy-window ladder
the V178 artifact rows establish (its rows and A_full are the ladder
being extended).  H_dressed(q) = H1(q) + Re Sigma(q, z), Sigma being
V252's dynamical self-energy in its CONTACT limit on the P-only
dressed system (the roll-up clause fired: V251 landed SECTOR-BLOCKED,
the flux sector was never realized, so d_A = 1; N = 16 shift grid,
delta = 1e-2, z = E_node + i delta with E_node at QSTAR).  The static
Fock insertion of V178 - named by the same roll-up clause - is run as
a second control ladder (it reproduces V178's A_gauged rows and is
extended to b = 8, 16).  Both dressings rebuild the step as
exp(i H_dressed) so effective_hamiltonian recovers it exactly
(measured).  The pre-registration's "(n_low = 4)" node subspace is
MEASURED to have no exact cut window at M on the implemented clock - the
window ladder jumps 70, 39, 18, 10, 6 because a degenerate pair at
|theta| = 0.073191 enters together - so the ladder rows report their
measured subspaces and a fixed n_low = 4 row is added as an auxiliary
measurement (addendum, the V242 precedent; the ladder itself is
unchanged).  The config marginal band continuum.relevance_marginal_band
= 1e-2 is the comparator.

**Measured (pre-registered).**  v(b), A(b), m(b) per cut (bare and
dressed); the log-log exponents p_v, p_A against the band.

**Acceptance (pre-registered).**  Branch CONTINUUM-WINDOW iff
p_A < -1e-2 (anisotropy relevant: the coarse-grained cone rotates) -
report b*, the cut where A(b) < 8e-15 would extrapolate, honestly
labeled as an extrapolation with its exponent; else NO-WINDOW with
(p_v, p_A, the full b-table).
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")
from cl10.paths import exp                                     # noqa: E402

from cl10 import altcell as AC                                 # noqa: E402
from cl10 import config as _config                             # noqa: E402
from cl10 import coupling                                      # noqa: E402
from cl10 import supercell as SC                               # noqa: E402
from cl10 import engine as E                                  # noqa: E402

G10 = coupling.G10
QSTAR = (math.pi / 2, math.pi / 2, 0.0)
DQ = 1e-3
NSHIFT = 16
DELTA = 1e-2
BS = (1, 2, 4, 8, 16)


_HTARGET = {}


def _hs_chunk(rng):
    i0, i1 = rng
    H = _HTARGET["H"]
    bf = _HTARGET["bf"]
    sh = _HTARGET["shifts"]
    for i in range(i0, i1):
        H[i] = AC.effective_hamiltonian(bf, tuple(sh[i]))
def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def aniso(vs):
    return (max(vs) - min(vs)) / (max(vs) + min(vs))


def main():
    t0 = time.time()
    checks = []
    data = {}

    def log(msg):
        print("[V253 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg),
            flush=True)

    log("start, backend={}".format(E._BACKEND))
    gpu = E._BACKEND != "numpy"
    xp = E._E if gpu else np

    _B, _U32, U16 = SC.sectors(G10)
    _Bm, T16, _f16 = AC.weyl_generators()


    def bf(q):
        return SC.bloch(q, U16, G10)

    T16s = np.stack(T16)

    # shift set and its eigendecomposition (shared by every Sigma call)
    steps = np.arange(NSHIFT) * (2 * math.pi / NSHIFT)
    shifts = np.array([(x, y, z) for x in steps for y in steps
                       for z in steps])
    import multiprocessing as mproc
    from multiprocessing import shared_memory
    npts = NSHIFT ** 3
    shm = shared_memory.SharedMemory(
        create=True, size=npts * 128 * 128 * 16)
    Hs = np.ndarray((npts, 128, 128), dtype=np.complex128,
                    buffer=shm.buf)
    _HTARGET["H"] = Hs
    _HTARGET["shifts"] = shifts
    _HTARGET["bf"] = bf
    log("H-stack build: {} points on 48 fork procs".format(npts))
    with mproc.get_context("fork").Pool(48) as pool:
        pool.map(_hs_chunk, [(i0, min(i0 + 64, npts))
                             for i0 in range(0, npts, 64)])
    log("H-stack build done")
    Hsd = xp.asarray(Hs)
    # Hs and shm stay alive: an independent host recomputation of the
    # Sigma kernel follows after the cross-check
    # batched device eigh in 2048 blocks (the full-stack call's
    # LAPACK scratch does not fit the card)
    Eall = xp.empty((NSHIFT ** 3, 128), dtype=np.float64)
    Vall = xp.empty((NSHIFT ** 3, 128, 128), dtype=np.complex128)
    for c0 in range(0, NSHIFT ** 3, 2048):
        sl = slice(c0, c0 + 2048)
        Eall[sl], Vall[sl] = xp.linalg.eigh(Hsd[sl])
    del Hsd
    log("device batched eigh done")

    Hr = AC.effective_hamiltonian(bf, QSTAR)
    th_abs = np.abs(np.linalg.eigvalsh(Hr))
    E_node = float(th_abs[int(np.argmin(th_abs))])
    z_fix = E_node + 1j * DELTA
    d_fix = np.ones(45)     # contact limit: sector floor (V251 blocked)
    K2d = xp.asarray(np.einsum("A,Aac,Abd->abcd", d_fix, T16s,
                               T16s.conj()))

    def _sigma_build():
        acc = xp.zeros((128, 128), dtype=np.complex128)
        Km = K2d.reshape(256, 256)
        for c0 in range(0, len(shifts), 2048):
            Vc = Vall[c0:c0 + 2048]
            w = 1.0 / (z_fix - Eall[c0:c0 + 2048])
            M = (Vc * w[:, None, :]) @ Vc.conj().swapaxes(-2, -1)
            Mm = M.reshape(-1, 8, 16, 8, 16) \
                 .transpose(0, 2, 4, 1, 3).reshape(-1, 256, 64)
            Q = Km @ Mm
            S = (Q.reshape(-1, 16, 16, 8, 8).transpose(0, 3, 1, 4, 2)
                 .reshape(-1, 128, 128))
            acc += xp.sum(S, axis=0)
        return G10 ** 2 * E.asnumpy(acc)

    SIGMA = _sigma_build()   # the contact kernel is q-independent:
                             # one z, one kernel, computed once

    def sigma(q):
        return SIGMA
    log("Sigma contact kernel built")

    def bfd(q):
        """the dressed step: the unitary exp(i (H1 + Re Sigma)); the
        i-log recovery of the generator holds mod 2 pi (recorded)."""
        H = AC.effective_hamiltonian(bf, q)
        Hd = H + np.real(sigma(q))
        ev, U = np.linalg.eigh(Hd)
        return U @ np.diag(np.exp(1j * ev)) @ U.conj().T

    # dressing consistency: the device GEMM kernel against an
    # INDEPENDENT host recomputation (host eigh of the same H stack,
    # einsum-form sandwich - a different algebra path). The i-log
    # round trip is recorded, not gated: with Re Sigma outside the
    # first Floquet zone, i log exp(i X) = X holds only mod 2 pi on
    # the spectrum; the ladder acts on the unitary exp(i(H1 + Re S)).
    Hd_q = AC.effective_hamiltonian(bf, QSTAR) + np.real(SIGMA)
    rt = AC.effective_hamiltonian(bfd, QSTAR)
    roundtrip = float(np.abs(rt - Hd_q).max())
    sig_norm = float(np.max(np.abs(np.linalg.eigvalsh(
        0.5 * (SIGMA + SIGMA.conj().T)))))
    wrap = float(np.max(np.abs(np.linalg.eigvalsh(Hd_q))) / math.pi)
    npt = NSHIFT ** 3
    Eall_h = np.empty((npt, 128), dtype=np.float64)
    Vall_h = np.empty((npt, 128, 128), dtype=np.complex128)
    for i in range(npt):
        Eall_h[i], Vall_h[i] = np.linalg.eigh(Hs[i])
    K2h = np.einsum("A,Aac,Abd->abcd", d_fix, T16s, T16s.conj())
    acc_h = np.zeros((128, 128), dtype=np.complex128)
    for c0 in range(0, npt, 512):
        sl = slice(c0, c0 + 512)
        w = 1.0 / (z_fix - Eall_h[sl])
        M = np.einsum("kem,km,kfm->kef", Vall_h[sl], w,
                      Vall_h[sl].conj())
        acc_h += np.einsum(
            "abcd,kxcyd->kxayb", K2h,
            M.reshape(-1, 8, 16, 8, 16)).reshape(
                -1, 128, 128).sum(axis=0)
    sig_machine = float(np.abs(SIGMA - G10 ** 2 * acc_h).max()
                        / max(np.abs(acc_h).max(), 1e-30))
    del Eall_h, Vall_h, Hs
    shm.close()
    shm.unlink()
    log("Sigma machine cross-check: defect {:.2e}, |Re Sigma| max eig "
        "{:.2e}, Floquet-zone ratio {:.1f}".format(sig_machine,
                                                   sig_norm, wrap))

    dirs = AC.sphere_dirs(6)
    rows = []
    for b in BS:
        cut = math.pi / (2 * b)
        n_low = int(np.sum(th_abs < cut))
        edge = float(np.min(np.abs(th_abs - cut)))
        vs_bare, vs_dress = [], []
        for d in dirs:
            q = tuple(a + DQ * x for a, x in zip(QSTAR, d))
            _Hb, ev_b = AC.lowdin_downfold(bf, q, QSTAR, cut=cut)
            vs_bare.append(float(np.min(np.abs(ev_b))) / DQ)
            _Hd, ev_d = AC.lowdin_downfold(bfd, q, QSTAR, cut=cut)
            vs_dress.append(float(np.min(np.abs(ev_d))) / DQ)
        _Hd0, ev0 = AC.lowdin_downfold(bfd, QSTAR, QSTAR, cut=cut)
        _Hb0, evb0 = AC.lowdin_downfold(bf, QSTAR, QSTAR, cut=cut)
        rows.append({"b": b, "cut": cut, "n_low": n_low,
                     "window_margin": edge,
                     "A_bare": aniso(vs_bare),
                     "A_dressed": aniso(vs_dress),
                     "v_bare_max": max(vs_bare),
                     "v_dressed_max": max(vs_dress),
                     "m_bare": float(np.min(np.abs(evb0))),
                     "m_dressed": float(np.min(np.abs(ev0)))})

    # auxiliary fixed-dimension node subspace (the pre-reg n_low = 4)
    vs4b, vs4d = [], []
    for d in dirs:
        q = tuple(a + DQ * x for a, x in zip(QSTAR, d))
        _H, evb = AC.lowdin_downfold(bf, q, QSTAR, n_low=4)
        vs4b.append(float(np.min(np.abs(evb))) / DQ)
        _H, evd = AC.lowdin_downfold(bfd, q, QSTAR, n_low=4)
        vs4d.append(float(np.min(np.abs(evd))) / DQ)
    fixed4 = {"n_low": 4, "A_bare": aniso(vs4b),
              "A_dressed": aniso(vs4d), "v_bare_max": max(vs4b),
              "v_dressed_max": max(vs4d)}

    # control ladder: V178's static Fock insertion, extended
    sig = AC.fock_self_energy(bf, QSTAR, SC.G10, 16,
                              [(0.0, 0.0, 0.0)]
                              + [tuple(DQ * x for x in d) for d in dirs],
                              sub=(0.15, 0.6), nodes=SC.M_POINTS,
                              current_vertex=True)
    srows = []
    for b in BS:
        cut = math.pi / (2 * b)
        vs = []
        for i, d in enumerate(dirs):
            q = tuple(a + DQ * x for a, x in zip(QSTAR, d))
            Heff, ev = AC.lowdin_downfold(bf, q, QSTAR, cut=cut)
            s = float(np.real((sig[1 + i] - sig[0])[0, 0]))
            ev2 = np.linalg.eigvalsh(
                Heff + np.diag([s] + [0.0] * (Heff.shape[0] - 1)))
            vs.append(float(np.min(np.abs(ev2))) / DQ)
        srows.append({"b": b, "cut": cut, "A_static": aniso(vs)})
    log("ladder rows + static control done (dressed rows: {})"
        .format(len(rows)))

    # V178 ladder extension anchor: reproduce V178's bare rows b<=4
    v178 = art("V178_gauged_block_rg")
    v178_rows = v178["data"]["rows"]
    anchor = max(abs(r["A_bare"] - next(
        x for x in v178_rows if x["b"] == r["b"])["A_bare"])
        for r in rows[:3])
    anchor_g = max(abs(r["A_static"] - next(
        x for x in v178_rows if x["b"] == r["b"])["A_gauged"])
        for r in srows[:3])

    # log-log exponents over b = 1..16
    def slope(xs, ys):
        lx = np.log(np.array(xs, float))
        ly = np.log(np.maximum(np.array(ys, float), 1e-300))
        cx = lx - lx.mean()
        return float((cx * (ly - ly.mean())).sum() / (cx * cx).sum())

    p_v_bare = slope([r["b"] for r in rows], [r["v_bare_max"]
                                              for r in rows])
    p_v_dress = slope([r["b"] for r in rows], [r["v_dressed_max"]
                                               for r in rows])
    p_A_bare = slope([r["b"] for r in rows], [r["A_bare"] for r in rows])
    p_A_dress = slope([r["b"] for r in rows],
                      [r["A_dressed"] for r in rows])
    p_A_static = slope([r["b"] for r in srows],
                       [r["A_static"] for r in srows])

    band = float(_config.CFG.continuum.relevance_marginal_band)

    # b* extrapolation of the dressed anisotropy power law
    lp = (math.log(max(rows[-1]["A_dressed"], 1e-300))
          - p_A_dress * math.log(rows[-1]["b"]))
    _ex = (math.log(8e-15) - lp) / p_A_dress if p_A_dress < 0 else None
    bstar = (math.exp(_ex) if _ex is not None and _ex < 700.0
             else float("inf"))

    # ---------------- checks ---------------------------------------------
    checks.append({
        "name": "T0 V178 ladder extension anchors: the V253 bare rows "
                "reproduce V178's A_bare at b = 1, 2, 4 to {:.2e} and "
                "the static-insertion control rows reproduce V178's "
                "A_gauged to {:.2e} (V178 rows {}, A_full {}; ladder "
                "extended to b = 8, 16)".format(
                    anchor, anchor_g,
                    [round(r["A_bare"], 6) for r in rows[:3]],
                    round(v178["data"]["A_full"], 6)),
        "pass": bool(anchor < 1e-6 and anchor_g < 1e-6)})
    checks.append({
        "name": "T1 Sigma machine consistency: the device GEMM kernel "
                "matches an independent host recomputation (host eigh "
                "+ einsum-form sandwich on the same shift stack) to "
                "{:.2e}; i-log round trip {:.2e} recorded - max "
                "|Re Sigma| = {:.2e} puts the generator {:.1f}x "
                "outside the first Floquet zone (recovery mod 2 pi; "
                "the ladder acts on the unitary)".format(
                    sig_machine, roundtrip, sig_norm, wrap),
        "pass": bool(sig_machine < 1e-10)})
    checks.append({
        "name": "T2 ladder complete on all cuts (b = {}): measured "
                "window subspaces n_low {} - the pre-registered "
                "n_low = 4 has NO exact cut window at this node (the "
                "pair at |theta| = 0.073191 enters together; addendum "
                "recorded) - window margins {}; auxiliary fixed "
                "n_low = 4 row: A_bare {:.4f}, A_dressed {:.4f}"
                .format(list(BS), [r["n_low"] for r in rows],
                        ["{:.1e}".format(r["window_margin"])
                         for r in rows], fixed4["A_bare"],
                        fixed4["A_dressed"]),
        "pass": bool(len(rows) == 5 and all(
            r["window_margin"] > 1e-5 for r in rows))})
    checks.append({
        "name": "T3 exponents against the marginal band {}: p_v(bare) = "
                "{:.5f}, p_v(dressed) = {:.5f}, p_A(bare) = {:.5f}, "
                "p_A(dressed, primary) = {:.5f}, p_A(static control) "
                "= {:.5f}".format(band, p_v_bare, p_v_dress, p_A_bare,
                                  p_A_dress, p_A_static),
        "pass": all(math.isfinite(x) for x in (p_v_bare, p_v_dress,
                                               p_A_bare, p_A_dress,
                                               p_A_static))})

    # ---------------- branch + data --------------------------------------
    window = p_A_dress < -1e-2
    branch = ("CONTINUUM-WINDOW" if window
              else "NO-WINDOW") + ":sector-floor"
    data["rows"] = rows
    data["static_control_rows"] = srows
    data["fixed_nlow4_row"] = fixed4
    data["exponents"] = {"p_v_bare": p_v_bare, "p_v_dressed": p_v_dress,
                         "p_A_bare": p_A_bare, "p_A_dressed": p_A_dress,
                         "p_A_static": p_A_static}
    data["comparator_band"] = band
    data["bstar_extrapolation"] = {
        "value": bstar,
        "note": ("b* where A(b) would reach 8e-15 under the measured "
                 "power law: an EXTRAPOLATION of the fitted exponent, "
                 "not a measured cut" if math.isfinite(bstar) else
                 "no finite b* (dressed anisotropy not decreasing)")}
    data["sigma_conventions"] = {"grid_N": NSHIFT, "delta": DELTA,
                                 "d_mode": "contact d_A = 1 "
                                 "(sector floor: V251 SECTOR-BLOCKED)",
                                 "z": [E_node, DELTA],
                                 "origin": "V252's Sigma machine, P-only"}
    data["v178_anchor_dev"] = anchor
    data["v178_static_anchor_dev"] = anchor_g
    data["roundtrip_defect"] = roundtrip
    data["sigma_max_eig"] = sig_norm
    data["floquet_zone_ratio"] = wrap
    data["sigma_machine_defect"] = sig_machine
    data["backend"] = E._BACKEND
    data["fitted_inputs"] = []
    ok = all(c["pass"] for c in checks)
    log("checks: {}/{} pass".format(sum(bool(c["pass"])
                                        for c in checks), len(checks)))
    status = (
        "block-RG window on the dressed sector ({} s): V178's downfold "
        "ladder (A_bare reproduced to {:.1e}) extended to b = 8, 16 "
        "with V251/V252's dynamical Sigma inserted in the step "
        "(Sigma machine cross-checked to {:.1e}; generator wrap "
        "{:.1f}x recorded); measured window ladders n_low "
        "{} (the pre-reg n_low = 4 has no exact cut here - addendum). "
        "Branch {}: p_v = ({:.4f} bare, {:.4f} dressed), p_A = "
        "({:.4f} bare, {:.4f} dressed, {:.4f} static control) against "
        "the marginal band {}"
        "{}").format(
            time.time() - t0, anchor, sig_machine, wrap,
            [r["n_low"] for r in rows], branch, p_v_bare, p_v_dress,
            p_A_bare, p_A_dress, p_A_static, band,
            "; the coarse-grained cone ROTATES: b* extrapolating to "
            "8e-15 at b* = {:.3g} (labeled extrapolation with its "
            "exponent)".format(bstar) if window else
            ": the dressed sector does not flow - anisotropy is not "
            "relevant; full b-table recorded")
    out = {"id": "V253_block_rg_window",
           "tier": 2,
           "property": "V178's downfold ladder extended to b = 8, 16 on "
                       "the Sigma-dressed step: v(b), A(b), m(b) and "
                       "the log-log exponents against the marginal band"
                       " (V253)",
           "residual": True,
           "programme": "V253",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v253_block_rg_window.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
