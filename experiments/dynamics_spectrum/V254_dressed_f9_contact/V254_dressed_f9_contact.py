#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V254 - the dressed F9 contact at the deepest cut (tier 5,
residual: true, F9).

**The question (pre-registered).**  Does the dressed effective step
satisfy the SME electron velocity-anisotropy bound at the deepest
block-RG cut, measured (not extrapolated)?

**The object (closed input set).**  V253's artifact at b = 16: the same
conventions (P-only contact dressing of V252's machine, grid N = 16,
delta = 1e-2, z = E_node + i delta at the node) rebuild H_eff(q) =
i log of the b = 16 downfold of the dressed step; V235's exact
normalization dvv = (v_max - v_min)/v_max with altcell.cone_anisotropy
directions (12 rays, dq = 1e-3); the same node class V235 used (M,
implemented clock, all 12) AND the star loci (V244 orbits members);
the bound 8.0e-15 from config bounds.sme_electron_velocity_aniso.
T0 asserts the rebuild equals V253's recorded A(b = 16) row exactly.

**Measured (pre-registered).**  dvv at b = 16 per node; the A(b = 16)
value and V253's exponent as supporting numbers.

**Acceptance (pre-registered).**  Branch F9-DRESSED-CLOSED iff
dvv < 8.0e-15 measured (not extrapolated) with fitted_inputs == [];
else F9-STILL-FIRED with the ratio.
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
DQ = 1e-3


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


def main():
    t0 = time.time()
    checks = []
    data = {}
    gpu = E._BACKEND != "numpy"
    xp = E._E if gpu else np

    def log(msg):
        print("[V254 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg),
            flush=True)

    log("start, backend={}".format(E._BACKEND))

    core = E._BACKEND
    if gpu:
        # the b = 16 downfold stack needs ~2.7 GiB of its own; under 8-job
        # GPU sharing a sibling worker can starve the card mid-run, so pin
        # to the host core when the free memory cannot host the stack
        free = int(xp.cuda.runtime.memGetInfo()[0])
        if free < (3500 << 20):
            gpu, xp, core = False, np, "numpy (pinned: GPU free {} MiB)".format(
                free >> 20)
            log("GPU free {} MiB: pinned to host numpy".format(free >> 20))

    v253 = art("V253_block_rg_window")
    v252 = art("V252_dynamical_self_energy_loop")
    conv = v253["data"]["sigma_conventions"]
    NSHIFT = int(conv["grid_N"])
    DELTA = float(conv["delta"])
    assert "contact" in conv["d_mode"]          # sector-floor convention
    bound = float(_config.CFG.bounds.sme_electron_velocity_aniso)

    _B, _U32, U16 = SC.sectors(G10)
    _Bm, T16, _f16 = AC.weyl_generators()
    T16s = np.stack(T16)

    def bf(q):
        return SC.bloch(q, U16, G10)

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
    del Hs
    shm.close()
    shm.unlink()
    # batched device eigh in 2048 blocks
    Eall = xp.empty((NSHIFT ** 3, 128), dtype=np.float64)
    Vall = xp.empty((NSHIFT ** 3, 128, 128), dtype=np.complex128)
    for c0 in range(0, NSHIFT ** 3, 2048):
        sl = slice(c0, c0 + 2048)
        Eall[sl], Vall[sl] = xp.linalg.eigh(Hsd[sl])
    del Hsd
    log("device batched eigh done")
    K2d = xp.asarray(np.einsum("A,Aac,Abd->abcd", np.ones(45), T16s,
                               T16s.conj()))      # contact d_A = 1

    _SIGZ = {}

    def dressed(q, z):
        key = (z.real, z.imag)
        sig = _SIGZ.get(key)
        if sig is None:
            acc = xp.zeros((128, 128), dtype=np.complex128)
            Km = K2d.reshape(256, 256)
            for c0 in range(0, len(shifts), 2048):
                Vc = Vall[c0:c0 + 2048]
                w = 1.0 / (z - Eall[c0:c0 + 2048])
                M = (Vc * w[:, None, :]) @ Vc.conj().swapaxes(-2, -1)
                Mm = M.reshape(-1, 8, 16, 8, 16) \
                     .transpose(0, 2, 4, 1, 3).reshape(-1, 256, 64)
                Q = Km @ Mm
                S = (Q.reshape(-1, 16, 16, 8, 8)
                     .transpose(0, 3, 1, 4, 2).reshape(-1, 128, 128))
                acc += xp.sum(S, axis=0)
            sig = np.real(G10 ** 2 * E.asnumpy(acc))
            _SIGZ[key] = sig
        Hd = AC.effective_hamiltonian(bf, q) + sig
        ev, U = np.linalg.eigh(Hd)
        return U @ np.diag(np.exp(1j * ev)) @ U.conj().T

    def bfn_dressed(q, z):
        return dressed(q, z)

    def node_z(k):
        Hk = AC.effective_hamiltonian(bf, k)
        w = np.linalg.eigvalsh(Hk)
        En = float(np.min(np.abs(w)))   # V253's convention: |theta|
        return En + 1j * DELTA

    nodes = {"M": [tuple(q) for q in SC.M_POINTS]}
    v244 = art("V244_icosahedral_node_census")
    nodes["star"] = [tuple(m) for orb in v244["data"]["orbits"]
                     for m in orb["members"]]
    dirs = AC.sphere_dirs(6)
    cut = math.pi / (2 * 16)

    rows = []
    for cls, qs in nodes.items():
        for k in qs:
            z = node_z(k)
            vs = []
            for d in dirs:
                q = tuple(a + DQ * x for a, x in zip(k, d))
                _H, ev = AC.lowdin_downfold(
                    lambda qq: bfn_dressed(qq, z), q, k, cut=cut)
                vs.append(float(np.min(np.abs(ev))) / DQ
                          if len(ev) else 0.0)
            vmax, vmin = max(vs), min(vs)
            rows.append({"class": cls, "q": [float(x) for x in k],
                         "v_max": vmax, "v_min": vmin,
                         "dvv": (vmax - vmin) / vmax if vmax > 0
                                else 1.0,
                         "A_cone": (vmax - vmin) / (vmax + vmin)
                                   if vmax + vmin > 0 else 1.0})
            if len(rows) % 8 == 0:
                log("dressed cone rows: {}".format(len(rows)))

    # T0: the rebuild must equal V253's recorded b = 16 dressed row
    qstar = (math.pi / 2, math.pi / 2, 0.0)
    zq = node_z(qstar)
    vrec = []
    for d in dirs:
        q = tuple(a + DQ * x for a, x in zip(qstar, d))
        _H, ev = AC.lowdin_downfold(lambda qq: bfn_dressed(qq, zq),
                                    q, qstar, cut=cut)
        vrec.append(float(np.min(np.abs(ev))) / DQ
                    if len(ev) else 0.0)
    log("V253 anchor rebuild done")
    a_rec = (max(vrec) - min(vrec)) / (max(vrec) + min(vrec))
    anchor = abs(a_rec - v253["data"]["rows"][-1]["A_dressed"])

    best = min(rows, key=lambda r: r["dvv"])
    fired = best["dvv"] >= bound

    # ---------------- checks ---------------------------------------------
    checks.append({
        "name": "T0 rebuild control: the V254 b = 16 dressed effective "
                "step at QSTAR reproduces V253's recorded A(b = 16) "
                "dressed row to {:.2e} (same conventions: grid {}, "
                "delta {}, contact)".format(anchor, NSHIFT, DELTA),
        "pass": bool(anchor < 1e-9)})
    checks.append({
        "name": "T1 measured (not extrapolated): {} node classes x 12 "
                "rays per node, dq = {:.0e}, cut pi/32 - every dvv is a "
                "direct velocity spread of the b = 16 downfolded "
                "dressed step".format(len(rows), DQ),
        "pass": bool(len(rows) == 12 + len(nodes["star"])
                     and all(math.isfinite(r["dvv"]) for r in rows))})
    checks.append({
        "name": "T2 supporting numbers: V253 A(b = 16) dressed "
                "{:.6f}, p_A(dressed) {:.4f}, p_A(static control) "
                "{:.4f}; V252 loop rows {} (sector floor)".format(
                    v253["data"]["rows"][-1]["A_dressed"],
                    v253["data"]["exponents"]["p_A_dressed"],
                    v253["data"]["exponents"]["p_A_static"],
                    v252["branch"]),
        "pass": bool(v253["data"]["rows"][-1]["A_dressed"] > 0)})

    branch = "F9-STILL-FIRED" if fired else "F9-DRESSED-CLOSED"
    data["rows"] = rows
    data["best"] = best
    data["bound"] = bound
    data["ratio_best"] = best["dvv"] / bound
    data["support"] = {"A_b16_dressed": v253["data"]["rows"][-1]
                       ["A_dressed"],
                       "p_A_dressed": v253["data"]["exponents"]
                       ["p_A_dressed"]}
    data["anchor_dev"] = anchor
    data["conventions"] = conv
    data["backend"] = core
    data["fitted_inputs"] = []
    ok = all(c["pass"] for c in checks)
    log("checks: {}/{} pass".format(sum(bool(c["pass"])
                                        for c in checks), len(checks)))
    status = (
        "dressed F9 contact at the deepest cut ({:.0f} s): the b = 16 "
        "downfold of the dressed step, measured at 12 M nodes + {} "
        "star loci x 12 rays (rebuild reproduces V253's A(b = 16) to "
        "{:.1e}). Branch {}: best dvv = {:.4f} (node {}) vs the SME "
        "bound {:.1e} - ratio {:.3e}{}.".format(
            time.time() - t0, len(nodes["star"]), anchor, branch,
            best["dvv"], best["q"], bound, best["dvv"] / bound,
            ": the dressed step does not linearize the cone - F9 "
            "stays fired at the deepest measured cut" if fired else
            ": measured contact"))
    out = {"id": "V254_dressed_f9_contact",
           "tier": 5,
           "property": "the SME electron velocity-anisotropy contact "
                       "measured on the b = 16 dressed downfolded step "
                       "at every named node (V254)",
           "residual": True,
           "programme": "V254",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v254_dressed_f9_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
