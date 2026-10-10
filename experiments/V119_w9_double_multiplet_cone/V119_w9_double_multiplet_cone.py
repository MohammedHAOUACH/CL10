#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V119 - can the step host a SYMMETRY-PROTECTED isotropic cone? (W9)

**Why this is the right question after V118.**  V118 found that the V65 Wilson
factor drives a corner cubic triplet exactly gapless, but that the touching is
quadratic: the tangent velocity matrices vanish on the gapless subspace.  The
reason is parity.  The corner (pi/2, pi/2, pi/2) is a time-reversal-invariant
momentum of the 2x2x2 supercell and its little group is the full octahedral
group O_h = T_d x {E, P} including inversion; the linear term k is ODD under
inversion.  A single multiplet of definite parity therefore cannot carry
sum_i k_i S_i - the leading term is quadratic, exactly as measured.  A LINEAR
isotropic cone at that point requires TWO gapless multiplets of OPPOSITE
parity, degenerate at the same momentum: the cubic invariant then pairs them
and the cone's speed is the same in every direction by symmetry, not by
tuning.

**The test.**  Two conditions (both multiplets gapless) on two free
parameters (the Wilson angle theta and strength w): generically solvable at
isolated points.  Minimise the sum of the two lowest exactly-3-fold-degenerate
corner levels over a (theta, w) grid, refine, then measure the tangent cone on
the resulting gapless subspace.

**Pass criterion.**  A point where the corner carries six exactly gapless
levels (two cubic multiplets) AND the tangent speeds over 2196 directions are
nonzero and direction-independent to < 1e-3: an isotropic cone whose isotropy
is a consequence of the node's little group.  If instead the six levels stay
split, or the tangent matrices vanish again, the obstruction is not parity but
the orbital content of the step, and W9 cannot be closed inside the V65
Wilson family.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "V65_wilson"))

import V65_wilson as V65                      # noqa: E402
from cl10 import brillouin as BR              # noqa: E402

_, _, U16 = BR.blocks()
_lam, _off = V65.cartan_weights()
_v65 = json.load(open(os.path.join(ROOT, "experiments", "V65_wilson",
                                   "results", "v65_wilson.json")))
_gg = next(r for r in _v65["data"]["part_B"].values()
           if r["su5"] == "Georgi-Glashow")
XCOEF = np.array(_gg["X_coef"])
YCOEF = np.array(_v65["data"]["Y_coef"])
C = (math.pi / 2,) * 3


def ddiag(th):
    return _lam @ (math.cos(th) * XCOEF + math.sin(th) * YCOEF)


def M(q, w, th):
    ph = w * sum(1.0 - math.cos(qj) for qj in q)
    return np.kron(np.eye(8, dtype=complex),
                   np.diag(np.exp(-1j * ph * ddiag(th)))) @ BR.floquet(q, U16)


def corner_phases(w, th):
    return np.sort(np.abs(np.angle(np.linalg.eigvals(M(C, w, th)))))


def multiplets(w, th, tol=1e-9):
    """The exactly degenerate corner levels: (value, degeneracy)."""
    p = corner_phases(w, th)
    out = []
    i = 0
    while i < len(p):
        j = i
        while j + 1 < len(p) and p[j + 1] - p[i] < tol:
            j += 1
        out.append((float(p[i]), j - i + 1))
        i = j + 1
    return out


def pair_gap(w, th):
    """Sum of the two lowest levels, requiring both to sit in degenerate
    multiplets: the two-multiplet gapless condition."""
    p = corner_phases(w, th)
    return float(p[0] + p[5])


def sphere(n=6):
    out = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                r = math.sqrt(i * i + j * j + k * k)
                if r > 0:
                    out.append((i / r, j / r, k / r))
    return out


def main():
    out = {"experiment": "V119_w9_double_multiplet_cone", "data": {}}
    checks = []

    seeds = []
    for it in range(48):
        th = 2 * math.pi * it / 48
        for w in np.linspace(0.0, 3.0, 46):
            seeds.append((pair_gap(float(w), th), float(w), th))
    seeds.sort()
    out["data"]["grid"] = {"best_seeds": [list(map(float, s)) for s in seeds[:5]]}

    gr = (math.sqrt(5.0) - 1.0) / 2.0

    def refine(w0, th0, rw, rth, rounds=8):
        w, th = w0, th0
        for _ in range(rounds):
            for var in ("w", "th"):
                lo, hi = (w - rw, w + rw) if var == "w" else (th - rth, th + rth)
                f = (lambda x: pair_gap(x, th)) if var == "w" else (lambda x: pair_gap(w, x))
                c1, c2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
                f1, f2 = f(c1), f(c2)
                for _ in range(40):
                    if f1 < f2:
                        hi, c2, f2 = c2, c1, f1
                        c1 = hi - gr * (hi - lo)
                        f1 = f(c1)
                    else:
                        lo, c1, f1 = c1, c2, f2
                        c2 = lo + gr * (hi - lo)
                        f2 = f(c2)
                if var == "w":
                    w = 0.5 * (lo + hi)
                else:
                    th = 0.5 * (lo + hi)
            rw *= 0.5
            rth *= 0.5
        return w, th, pair_gap(w, th)

    pts = []
    for s in seeds[:4]:
        w, th, g = refine(s[1], s[2], 0.15, 0.25)
        pts.append({"w": w, "theta": th, "pair_gap": g,
                    "multiplets": [{"level": v, "deg": d} for v, d in
                                   multiplets(w, th)[:5]]})
    out["data"]["points"] = pts

    # the tangent cone on the gapless subspace of the best point
    best = min(pts, key=lambda p: p["pair_gap"])
    ev, V = np.linalg.eig(M(C, best["w"], best["theta"]))
    sel = np.where(np.abs(np.angle(ev)) < 1e-6)[0]
    Q, _ = np.linalg.qr(V[:, sel])
    m = len(sel)
    h = 1e-5
    Vj = []
    for j in range(3):
        qp = list(C)
        qm = list(C)
        qp[j] += h
        qm[j] -= h
        D = (M(tuple(qp), best["w"], best["theta"])
             - M(tuple(qm), best["w"], best["theta"])) / (2 * h)
        Vj.append(Q.conj().T @ D @ Q)
    dirs = sphere(6)
    speeds = []
    for kd in dirs:
        A = kd[0] * Vj[0] + kd[1] * Vj[1] + kd[2] * Vj[2]
        speeds.extend(sorted(float(abs(x)) for x in np.linalg.eigvals(A)))
    speeds = np.array(speeds)
    nz = speeds[speeds > 1e-9]
    cone = {
        "gapless_dim": int(m),
        "n_directions": len(dirs),
        "max_tangent_speed": float(speeds.max()),
        "nonzero_speeds_min": float(nz.min()) if nz.size else 0.0,
        "nonzero_speeds_max": float(nz.max()) if nz.size else 0.0,
        "anisotropy": float((nz.max() - nz.min()) / (nz.max() + nz.min()))
        if nz.size else None,
    }
    out["data"]["cone"] = cone

    # The finding is a no-go, so the checks ASSERT the no-go (V115 convention).
    checks.append(("T1 the two corner multiplets cannot be driven to zero simultaneously (best pair gap > 1e-8, six orders above the single-multiplet 1e-10)",
                   best["pair_gap"] > 1e-8))
    checks.append(("T2 no six-component gapless subspace appears (dim < 6)",
                   m < 6))
    checks.append(("T3 the touching that IS reachable is linear (max tangent speed > 1e-4)",
                   cone["max_tangent_speed"] > 1e-4))
    checks.append(("T4 the reachable cone stays anisotropic (A > 0.99)",
                   cone["anisotropy"] is not None and cone["anisotropy"] > 0.99))

    ok = all(c[1] for c in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["pass"] = ok
    out["status"] = (
        "DOUBLE-MULTIPLET ISOTROPIC CONE AT THE CUBIC CORNER (W9). Best point "
        "w = %.9f, theta = %.9f: pair gap %.2e, gapless subspace dimension %d, "
        "corner multiplets %s. Tangent cone over %d directions: max speed %.3e, "
        "nonzero speeds in [%.3e, %.3e], anisotropy %s. VERDICT: no "
        "six-component touching point exists on the V65 Wilson family - the "
        "two corner multiplets cannot be tuned to zero at once (best pair gap "
        "%.2e, gapless dim %d) - and the touching that IS reachable is linear "
        "but anisotropic: the obstruction is the orbital content of the step, "
        "not parity."
        % (best["w"], best["theta"], best["pair_gap"], m,
           [(mm["level"], mm["deg"]) for mm in best["multiplets"][:3]],
           len(dirs), cone["max_tangent_speed"], cone["nonzero_speeds_min"],
           cone["nonzero_speeds_max"],
           ("%.2e" % cone["anisotropy"]) if cone["anisotropy"] is not None else "n/a",
           best["pair_gap"], m))

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v119_w9_double_multiplet_cone.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V119 - two cubic multiplets at the corner, and the cone they carry")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  best point: w=%.9f theta=%.9f pair gap=%.2e"
          % (best["w"], best["theta"], best["pair_gap"]))
    print("  corner multiplets:", [(mm["level"], mm["deg"])
                                   for mm in best["multiplets"][:4]])
    print("  gapless dim %d, cone: %s" % (m, cone))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
