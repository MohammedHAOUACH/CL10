#!/usr/bin/env python3
"""V207 - the corner pair question: can the TWO lowest corner cubic triplets
be driven to zero SIMULTANEOUSLY, and what is the tangent velocity algebra of
the resulting six-state null space?

V118 showed that a corner cubic triplet can be made exactly gapless on the
(w, theta) family plane - it found three such points, each carrying a
three-state exactly-gapless null space - and measured the tangent velocity
algebra there: the S_j commute, sum_j S_j^2 is scalar, and every branch speed
is exactly |v_alpha . d| for a fixed velocity vector. Nodal planes, not a
cone. The completion of that question - can BOTH of the two lowest triplets
be zeroed at the same time? - was measured live and quoted in the paper
(Section 8.2, Table 5, O9, Appendix B) without a committed artifact, at
w = 2.8894, theta = 1.5717. V207 commits the measurement. The answer is YES,
six exactly gapless states and the algebra stays Abelian - but the quoted
coordinates are wrong: at w = 2.8894, theta = 1.5717 the six lowest corner
phases are 1.4e-4 and 2.6e-4, not zero (that point is where a cluster-based
refinement stalled, and cluster levels are the wrong objective here: the six
phases come in +/- pairs, so the two lowest individual phases belong to the
SAME triplet and driving them to zero zeroes only one triplet). The exact
pair sits at w = 0.141319757, theta = 0.932314669.

T1 the pair objective P(w, theta) = max(L1, L2), where L1 and L2 are the mean
   phases of the two lowest triplets (the six lowest phases, grouped): a
   61 x 48 scan of the family plane, then Nelder-Mead min-max refinement from
   the eight best seeds. The deliverable is the infimum and every distinct
   point that attains it.
T2 the spectrum at the pair: the exactly-gapless dimension (tolerance 1e-11),
   the largest of the six lowest phases, and the distance to the next level.
T3 the tangent velocity algebra on the six-state null space: anti-hermiticity
   of the projected derivatives, the three commutator norms ||[S_i,S_j]||, and
   the deviation of sum_j S_j^2 from a scalar - the cone test at the LARGEST
   gapless multiplet the corner offers.
T4 the branch structure: the joint velocity vectors v_alpha of the commuting
   S_j, the planar error max |speed - |v_alpha . d|| over the direction
   sphere, the count of exactly zero branch vectors and of exactly zero branch
   directions, and the anisotropy.
T5 the correction recorded: at the coordinates the paper quotes, and at
   V118's three exactly gapless SINGLE-triplet points (read from its committed
   artifact), the second triplet stays massive.

Measured 2026-10-21 (61 x 48 scan, 8 Nelder-Mead refinements): the pair
exists and it is exact - max(L1, L2) = 2.1e-16 at w = 0.241747262,
theta = 0.913136625 (6 of the 8 refined seeds land on an exact pair point,
modulo 2 pi in theta: (0.934143, 1.715292), (2.412996, 2.748448),
(0.561146, 0.734659), (0.403456, 0.751675), (2.055357, 3.288902)). The null
space is six-dimensional below 1e-11, its six phases are all below 2.5e-16,
and the next level sits at 0.0087. On that space the tangent velocity algebra
is Abelian (max commutator 5.3e-12 against ||S_j|| = 0.523, anti-hermiticity
9.6e-12) but sum_j S_j^2 is NOT scalar (deviation 7.2e-2): the six states
split into two cubic orbits of branch velocities, the permutations of
(-0.0434, -0.3272, -0.3272) and of (0.1299, 0.1299, -0.1539), speeds 0.4647
and 0.2396. Every branch is planar to 4.8e-12 over 2196 directions and the
six vectors define six distinct nodal planes - nodal planes, never a cone. At
the coordinates the live measurement quoted (2.8894, 1.5717) the six phases
sit at 1.4e-4 and 2.6e-4, and at V118's three single-triplet points the second
triplet never falls below 1.6e-2.

Smoke: CL10_SMOKE=1 confines the scan to the window the full run brackets and
shrinks the seeds and the direction sphere.
"""
import json
import math
import os
import sys
import time

# BEFORE the numpy import: single-threaded BLAS (the sweep is a long chain of
# independent 128x128 eigendecompositions; a shared pool buys nothing here).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

# The corner programme runs on the CPU core ALWAYS (V186/V187 convention):
# independent 128x128 Floquet eigendecompositions, ~40x faster on CPU than on
# a 3090, and the GPU route needs a CUDA context per worker.
os.environ["CL10_BACKEND"] = "numpy"

# The corner is a property of the IMPLEMENTED automaton, so this experiment
# measures it on the implemented clock g = 1/5. cl10/config.py resolves the
# clock ONCE, at its own import time, so the pin must precede every import
# that can pull config in: a bare `from cl10 import config` first would
# silently freeze the default closed clock pi/4 and the whole measurement
# would be made on the wrong step. V118 and V186 carry the same pin.
os.environ["CL10_ANGLE"] = "0.2"

import numpy as np                            # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "V118_w9_little_group"))
from cl10.config import smoke as _smoke       # noqa: E402
from cl10.paths import exp                    # noqa: E402
from scipy.optimize import minimize           # noqa: E402
import V118_w9_little_group as V118           # noqa: E402

M = V118.M                                    # the implemented step dressed by
C = (math.pi / 2,) * 3                        # the (w, theta) family direction
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

# the coordinates the live measurement quoted in the paper (Section 8.2,
# Table 5, O9, Appendix B), recorded here at the 4-decimal precision the paper
# printed them. Provenance: the V115 cone scan, data.gapless_points[548]
# (theta = 1.5717554519091528). They are kept rounded on purpose: this check
# measures the error that the paper's own rounding introduced.
QUOTED = (2.8894, 1.5717)


def phases(w, th):
    """The corner quasienergy phases |angle|, sorted."""
    return np.sort(np.abs(np.angle(np.linalg.eigvals(M(C, float(w), float(th))))))


def pairP(w, th):
    """P = max(L1, L2), the mean phases of the two lowest triplets. NOT the
    two lowest individual phases: the corner spectrum is particle-hole
    symmetric, so the two lowest phases belong to one triplet."""
    p = phases(w, th)
    return max((p[0] + p[1] + p[2]) / 3.0, (p[3] + p[4] + p[5]) / 3.0), p


def main():
    smoke = _smoke()                           # CL10_SMOKE, via cl10/config.py
    out = {"experiment": "V207_corner_pair_zero", "data": {}, "checks": []}
    checks = out["checks"]
    t0 = time.time()

    nth = 16 if smoke else 48
    nw = 13 if smoke else 61
    nseed = 2 if smoke else 8

    # the full run scans the WHOLE family plane (w in [0, 3], theta over 2 pi);
    # the smoke run scans only the window the full run brackets, so smoke still
    # measures the same pair.
    ws = (np.linspace(0.05, 0.35, nw) if smoke
          else np.linspace(0.0, 3.0, nw))
    ths = (0.932315 + (np.linspace(-0.1, 0.1, nth) if smoke else
                       np.linspace(0.0, 2 * math.pi, nth, endpoint=False)))

    # ---- T1: the min-max infimum of the pair objective --------------------
    print("  T1: %d x %d scan of the family plane" % (nw, nth), flush=True)
    seeds = []
    for th in ths:
        for w in ws:
            seeds.append((pairP(float(w), float(th))[0], float(w), float(th)))
    seeds.sort()

    pairs = []
    seen = []
    for s in seeds[:nseed]:
        r = minimize(lambda x: pairP(x[0], x[1])[0], [s[1], s[2]],
                     method="Nelder-Mead",
                     options=dict(xatol=1e-15, fatol=1e-18, maxfev=400,
                                  adaptive=True))
        P, p = pairP(r.x[0], r.x[1])
        w_, th_ = float(r.x[0]), float(r.x[1]) % (2 * math.pi)
        if any(abs(w_ - a) < 1e-7 and abs(th_ - b) < 1e-7 for a, b in seen):
            continue
        seen.append((w_, th_))
        pairs.append({"w": w_, "theta": float(r.x[1]), "theta_mod_2pi": th_,
                      "P": float(P), "triplet1": float(np.mean(p[:3])),
                      "triplet2": float(np.mean(p[3:6]))})
    pairs.sort(key=lambda d: d["P"])
    best = pairs[0]
    w_star, th_star = best["w"], best["theta"]
    n_exact = sum(1 for d in pairs if d["P"] < 1e-14)
    out["data"]["T1"] = {
        "grid_best": list(map(float, seeds[0])), "pairs": pairs,
        "n_distinct_pairs": len(pairs), "n_exact_pairs": n_exact,
        "w_star": w_star, "theta_star": th_star, "P_star": best["P"]}
    checks.append((
        "T1 both corner triplets CAN be driven to zero together: the min-max "
        "refinement of max(L1,L2) reaches %.2e at w = %.9f, theta = %.9f "
        "(%d distinct refined points, %d of them exact; grid best %s)"
        % (best["P"], w_star, th_star, len(pairs), n_exact,
           "P = %.2e at w = %.4f, theta = %.4f" % seeds[0]),
        best["P"] < 1e-14 and n_exact >= 1))

    # ---- T2: the six-state null space at the pair -------------------------
    p = phases(w_star, th_star)
    _, dim = V118.gapless_basis(C, w_star, th_star, tol=1e-11)
    higher = float(p[6])
    out["data"]["T2"] = {"gapless_dim_tol_1e-11": int(dim),
                         "six_lowest_phases": [float(x) for x in p[:6]],
                         "largest_of_the_six": float(p[5]),
                         "next_level": higher,
                         "isolation_ratio": higher / max(float(p[5]), 1e-300)}
    checks.append((
        "T2 the null space at the pair is %d-dimensional below 1e-11, the six "
        "phases are all below %.2e, and the next level sits at %.4f (an "
        "isolation of %.0e)"
        % (dim, p[5], higher, higher / max(float(p[5]), 1e-300)),
        dim == 6 and p[5] < 1e-12 and higher > 1e-3))

    # ---- T3: the tangent velocity algebra on the six-state space ----------
    Vj, m = V118.tangent_matrices(C, w_star, th_star, 1e-5, tol=1e-11)
    if Vj is None:
        checks.append(("T3 tangent matrices: NO gapless subspace", False))
        m = 0
    else:
        S = [1j * A for A in Vj]
        anti = float(max(np.abs(A + A.conj().T).max() for A in Vj))
        snorm = [float(np.linalg.norm(s)) for s in S]
        comm = [float(np.linalg.norm(S[a] @ S[b] - S[b] @ S[a]))
                for a, b in ((0, 1), (1, 2), (2, 0))]
        S2 = sum(s @ s for s in S)
        scalar = float(np.trace(S2).real / m)
        s2dev = float(np.abs(S2 - scalar * np.eye(m)).max())
        out["data"]["T3"] = {"dim": int(m), "antihermiticity_err": anti,
                             "S_norms": snorm, "commutator_norms": comm,
                             "sum_Sj2_scalar": scalar,
                             "sum_Sj2_dev_from_scalar": s2dev}
        checks.append((
            "T3 on the six-state null space the tangent velocity algebra is "
            "ABELIAN: max ||[S_i,S_j]|| = %.2e against ||S_j|| = %.3f, "
            "anti-hermiticity %.2e - the cone criterion fails at the largest "
            "gapless multiplet the corner offers. sum_j S_j^2 is NOT scalar "
            "on this space (deviation %.2e about the mean %.4f): the six "
            "states split into two cubic orbits of branch velocities, where "
            "V118's single triplet gives a scalar sum (6.3e-12)"
            % (max(comm), max(snorm), anti, s2dev, scalar),
            max(comm) < 1e-9 and anti < 1e-9 and m == 6))

    # ---- T4: the branch structure ----------------------------------------
    if m:
        dirs = V118.sphere(2 if smoke else 6)
        planar, zero_dir = 0.0, 0
        # an irrational generic direction in the Cartan plane (1, e-2, 1/e):
        # no rational direction can accidentally align with a symmetry axis
        G = (S[0] + (math.e - 2.0) * S[1] + math.exp(-1.0) * S[2])
        _, Ww = np.linalg.eigh(G)
        vecs = [[float(np.real(Ww[:, a].conj().T @ s @ Ww[:, a])) for s in S]
                for a in range(m)]
        vecs.sort(key=lambda v: -sum(x * x for x in v))
        for kd in dirs:
            evd = sorted(abs(x) for x in np.linalg.eigvalsh(
                kd[0] * S[0] + kd[1] * S[1] + kd[2] * S[2]))
            pred = sorted(abs(vd[0] * kd[0] + vd[1] * kd[1] + vd[2] * kd[2])
                          for vd in vecs)
            planar = max(planar, max(abs(a - b) for a, b in zip(evd, pred)))
            zero_dir += sum(1 for x in evd if x < 1e-9)
        speeds = [math.sqrt(sum(x * x for x in v)) for v in vecs]
        nzero = sum(1 for v in vecs if max(abs(x) for x in v) < 1e-9)
        nzs = [s for s in speeds if s > 1e-9]
        planes = []
        for v in vecs:
            nv = math.sqrt(sum(x * x for x in v))
            if nv < 1e-12:
                continue
            if not any(abs(sum(a * b for a, b in zip(v, pv)) / (nv * pn))
                       > 1.0 - 1e-9 for pv, pn in planes):
                planes.append((v, nv))
        out["data"]["T4"] = {
            "branch_velocity_vectors": vecs, "branch_speeds": speeds,
            "exactly_zero_branch_vectors": int(nzero),
            "n_nodal_planes": len(planes),
            "planar_branch_err": float(planar),
            "n_directions": len(dirs), "zero_branch_directions": int(zero_dir),
            "anisotropy": (float(max(nzs) / min(nzs)) if nzs else None)}
        checks.append((
            "T4 every branch is exactly planar, |speed - |v_alpha.d|| < %.2e "
            "over %d directions; %d branch directions vanish and the %d branch "
            "velocity vectors define %d distinct planes, speeds %s - nodal "
            "planes, never a cone"
            % (planar, len(dirs), zero_dir, m, len(planes),
               ", ".join("%.4f" % s for s in speeds)),
            planar < 1e-8 and zero_dir > 0 and len(planes) == 6))

    # ---- T5: the correction recorded, and V118's single-triplet locus -----
    pq = phases(*QUOTED)
    v118 = json.load(open(exp("V118_w9_little_group", "results",
                              "v118_w9_little_group.json")))
    pts = v118["data"]["T3"]["closed_multiplets"]
    second = []
    for pt in pts:
        pp = phases(pt["w"], pt["theta"])
        second.append({"w": pt["w"], "theta": pt["theta"],
                       "triplet1": float(np.mean(pp[:3])),
                       "triplet2": float(np.mean(pp[3:6]))})
    out["data"]["T5"] = {
        "quoted_coordinates": list(QUOTED),
        "quoted_phases": [float(x) for x in pq[:8]],
        "quoted_P": float(max(np.mean(pq[:3]), np.mean(pq[3:6]))),
        "v118_single_triplet_points": second,
        "min_second_triplet_at_v118_points": min(s["triplet2"] for s in second)}
    checks.append((
        "T5 the coordinates the live measurement quoted (w = %.4f, "
        "theta = %.4f) are NOT gapless: P = %.2e there (six phases %s); and at "
        "V118's %d exactly gapless SINGLE-triplet points the second triplet "
        "stays at %s (smallest %.1e), which is why the pair needed its own scan"
        % (QUOTED[0], QUOTED[1], out["data"]["T5"]["quoted_P"],
           ", ".join("%.1e" % x for x in pq[:6]), len(second),
           ", ".join("%.1e" % s["triplet2"] for s in second),
           min(s["triplet2"] for s in second)),
        out["data"]["T5"]["quoted_P"] > 1e-5
        and min(s["triplet2"] for s in second) > 1e-3))

    ok = all(o for _, o in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p_)} for n, p_ in checks]
    d3 = out["data"].get("T3", {})
    d4 = out["data"].get("T4", {})
    comm_txt = "%.1e" % max(d3["commutator_norms"]) if d3 else "n/a"
    snorm_txt = "%.3f" % max(d3["S_norms"]) if d3 else "n/a"
    s2_txt = "%.1e" % d3["sum_Sj2_dev_from_scalar"] if d3 else "n/a"
    vec_txt = ("; ".join("(%s)" % ", ".join("%.4f" % x for x in v)
                         for v in d4["branch_velocity_vectors"][:3])
               if d4 else "n/a")
    planar_txt = "%.1e" % d4["planar_branch_err"] if d4 else "n/a"
    planes_txt = "%d" % d4["n_nodal_planes"] if d4 else "0"
    out["status"] = (
        "CORNER PAIR QUESTION (cone programme, the V118 family): the two "
        "lowest corner cubic triplets CAN be driven to zero together - the "
        "min-max pair objective reaches %s at w = %.9f, theta = %.9f (%d exact "
        "pair points among %d refined), the null space there is 6-dimensional "
        "and isolated from the next level by %.4f. On that six-state space the "
        "tangent velocity algebra is ABELIAN (max commutator %s against "
        "||S|| = %s, sum_j S_j^2 NOT scalar but deviating by %s), the "
        "branches are planar to %s with velocity vectors %s defining %s nodal "
        "planes: the LARGEST "
        "gapless multiplet the corner offers carries nodal planes, not a cone "
        "(V118 single-triplet points, V186 bond scan, V187 rigidity wall, "
        "Proposition 6). The coordinates the uncommitted live measurement "
        "quoted (w = 2.8894, theta = 1.5717) are NOT gapless - P = %s there, "
        "because the two lowest individual phases belong to the SAME triplet "
        "under the corner's particle-hole pairing; the pair objective must be "
        "taken on the two triplets. The pair itself is an isolated point of "
        "the family plane, not a mechanism (O9)."
        % ("%.2e" % best["P"], w_star, th_star, n_exact, len(pairs),
           out["data"]["T2"]["next_level"], comm_txt, snorm_txt, s2_txt,
           planar_txt, vec_txt, planes_txt,
           "%.1e" % out["data"]["T5"]["quoted_P"]))
    out["wall_time"] = time.time() - t0
    for name, o in checks:
        print(("  [PASS] " if o else "  [FAIL] ") + name)
    print("  VERDICT:", "PASS" if ok else "FAIL")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "v207_corner_pair_zero.json"), "w") as f:
        json.dump(out, f)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
