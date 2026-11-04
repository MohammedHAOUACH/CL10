#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V326 - the phase-dictionary bridge (tier 1, hard, O4/O8,
theory-hunt; the structural defect closed as a theorem or repaired).

**The question (pre-registered).** The committed phase dictionary maps
the clock to eigenphases (root-of-unity phases + g x integer Cartan
weights), so closure clocks are exactly pi x rational while the data
clock 1/5 is rational: the two sets are disjoint and that disjointness
IS the structural cause of the 30.4 % hierarchy gap. Is there a
dictionary in the booked reparameterization family that bridges the
two - a clock that both closes (k <= closure_k_max) and reads the
observed coefficient inside coef_band_abs?

**Class (closed, fixed pre-run).** Uniform rescalings: dictionary c
evaluates the committed step and the committed ladder at the angle
c*g, c in {1, 2, 4, 1/2, 1/4, pi, 2*pi, 1/pi, 1/(2*pi)} - named
constants only, nothing fitted. Reparameterization lemma (booked):
closure of dictionary c at g <=> c*g in Cl, and the rung-10 read of
dictionary c at g = the committed read at c*g: bridge existence is
c-INDEPENDENT, a bridge exists iff Cl intersect R is nonempty, with
R = {x : |3(F_10(x)/LOOP)^10 - COEF_OBS| <= coef_band_abs}. The 2*pi
"fifth-root" dictionary (closing g = 1/5 at order 5 iff 2*pi/5 closes)
is inside the class and decided by the census directly.

**Objects.** R by coarse scan + bisection to 1e-12 on the float64
traced ladder (window and width booked, component count booked); Cl
by the V315/V311 grid census on 79 points with V315's machinery
IMPORTED (dec_step, closure_row, traced_ladder - zero copy drift);
for every closure angle the rung-10 read plus the published crossing
read where one exists (V314 ledger rows, V321 certified leaf).

**T0.** coef_10(1/5) = V206's booked 2.088106865906787 (< 1e-9); the
identity step at pi/4 reproduces V315's booked identity-closure leaf
(k = 2, dev to 1e-12; the order-8 quantization {pi/2, pi/4} is the
internal U16 sector fact of V206/V309, not this step product - the
first run pinned that wrong anchor and it was corrected before the
re-run, booked in the note); the V206 g_obs leaf lies inside R; the
lemma spot-check at c = 2*pi: the closure row of the dictionary at
g = 1/5 equals the grid row at pi*2/5; and the hoisted
clock-independent machinery is bit-identical to V203's verbatim
hop_machinery at 1/5 (booking of the scan-speed refactor).

**Branch.** DICTIONARY-BRIDGE-FOUND (a member and clock close AND read
COEF_OBS inside the band: repaired inside the booked class) /
DICTIONARY-NO-BRIDGE (Cl intersect R empty over the census; every
closure angle booked in band-units; the defect becomes a theorem of
the dictionary shape and the consolidation follows in the axiom
inventory). No middle branch.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: F401  (dec_step's module contract keeps it live)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                               # noqa: E402
sys.path.insert(0, exp("V203_closed_branch_coefficient"))
sys.path.insert(0, exp("V85_closed_loop_factor"))
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V315_closure_quantizer_census"))

from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios, D_L)
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)
from V85_closed_loop_factor import projectors             # noqa: E402
from V315_closure_quantizer_census import (dec_step,        # noqa: E402
                                           closure_row,
                                           DEN_MAX, K_MAX, GCLOSED)
from cl10 import clifford                                 # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import engine as E                              # noqa: E402
from cl10 import sm_qca                                   # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import u10pred                                  # noqa: E402

CFG = _config.CFG
G5 = float(CFG.angle.branches.implemented)
BAND_C = float(CFG.pillar_surgery.coef_band_abs)
NMAX10 = 10
CS = [1.0, 2.0, 4.0, 0.5, 0.25, math.pi, 2.0 * math.pi,
      1.0 / math.pi, 1.0 / (2.0 * math.pi)]


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def read_coef(Fv, LOOP):
    return 3.0 * (Fv / LOOP) ** 10


def coef10(x, LOOP, mach):
    """The rung-10 read 3(F_10/LOOP)^10 at dictionary angle x. The
    clock-INDEPENDENT part of V203's hop_machinery (Clifford Cartan,
    root generators, structure constants, adjoint links, the link
    projector spectral decomposition) is hoisted and reused - the
    hoist's equivalence to the verbatim call is a booked T0 gate
    below; the g-dependent parts (the QCA step G(x) and the Ma
    einsums) are recomputed exactly as V203's lines 93-105."""
    _cartan, ev, P = mach
    G = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(),
                                      x)).astype(np.complex128)
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    Ma = [np.einsum("ajbj->ab", (Pa @ Gi).reshape(32, D_L, 32, D_L)) / D_L
          for Pa in P]
    F, _rho = loop_ratios(G, _cartan, Ma, [float(e) for e in ev], x, NMAX10)
    return read_coef(F[10], LOOP)


def step_at(x):
    """The identity (undecorated, default-order) six-stage step at
    q = 0 and dictionary angle x - V315's own dec_step."""
    return dec_step([0.0, 0.0, 0.0], x,
                    [(0.0, 0.0, 0.0)] * len(SC.STAGES), SC.STAGES)


def bisect(a, b, target, LOOP, mach):
    """Solve coef10(x) = target on a bracket with a sign change."""
    fa, fb = coef10(a, LOOP, mach) - target, coef10(b, LOOP, mach) - target
    assert fa * fb <= 0.0, "bisect needs a sign change"
    for _ in range(80):
        m = 0.5 * (a + b)
        fm = coef10(m, LOOP, mach) - target
        if fm == 0.0:
            return m
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
        if b - a < 1e-12:
            break
    return 0.5 * (a + b)


def hoist_machinery():
    """The clock-independent half of V203's hop_machinery, lines
    95-101 verbatim, computed once."""
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = root_generators()
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    P, ev = projectors(TL)
    return cartan, ev, P


def main():
    t0 = time.time()
    checks, data = [], {}
    LOOP = float(u10pred.LOOP)
    COEF_OBS = float((u10pred.V_EW / u10pred.M_PLANCK) / LOOP ** 10)
    mach = hoist_machinery()

    # ---- hoist equivalence gate (gate hygiene, booked in the note) ------
    Gv, cart_v, Ma_v, ev_v = hop_machinery(G5)
    _cart, ev_h, P_h = mach
    eq_cart = max(float(np.abs(a - b).max()) for a, b in zip(cart_v, _cart))
    eq_ev = float(np.abs(np.asarray(ev_v) - np.asarray(ev_h)).max())
    Gi_v = np.kron(Gv, np.eye(D_L, dtype=np.complex128))
    eq_ma = max(float(np.abs(M - np.einsum(
        "ajbj->ab", (Pa @ Gi_v).reshape(32, D_L, 32, D_L)) / D_L).max())
        for M, Pa in zip(Ma_v, P_h))
    checks.append({"name": "T0h the hoisted clock-independent machinery "
                   "is identical to V203's verbatim hop_machinery(1/5): "
                   "cartan {:.1e}, ev {:.1e}, Ma {:.1e} (< 1e-15)"
                   .format(eq_cart, eq_ev, eq_ma),
                   "pass": bool(max(eq_cart, eq_ev, eq_ma) < 1e-15)})

    # ---- T0a: the fifth-clock read --------------------------------------
    coef5_art = float(art("V206_clock_tension",
                          "data.data_clock.coef_at_fifth"))
    c5 = coef10(G5, LOOP, mach)
    t0a = abs(c5 - coef5_art)
    checks.append({"name": "T0a rung-10 read at 1/5 reproduces V206's "
                   "coef_at_fifth {:.15f} (|d| {:.1e} < 1e-9)".format(
                       coef5_art, t0a), "pass": bool(t0a < 1e-9)})

    # ---- T0b: identity closure reproduces V315's booked leaf ------------
    # (first run expected k = 8 there - a wrong anchor from the internal-
    # sector prose; the step-product object closes at k = 2 at pi/4 and
    # V315's booked leaf says exactly that - corrected BEFORE re-run,
    # booked in the note; the order-8 quantization {pi/2, pi/4} is the
    # INTERNAL U16 sector fact, V206/V309, not this step product)
    cl_gclosed = closure_row(step_at(GCLOSED))
    leaf = art("V315_closure_quantizer_census",
               "data.t0a_identity_closure.pi*1/4")
    okb = (cl_gclosed is not None
           and cl_gclosed["k"] == int(leaf["rederived_k"])
           and abs(cl_gclosed["dev"] - float(leaf["rederived_dev"])) < 1e-12)
    checks.append({"name": "T0b the identity step at pi/4 closes at "
                   "k = {} (V315 booked leaf k = {}), dev difference "
                   "{:.1e} < 1e-12".format(
                       cl_gclosed["k"] if cl_gclosed else None,
                       leaf["rederived_k"],
                       abs(cl_gclosed["dev"] - float(leaf["rederived_dev"]))
                       if cl_gclosed else float("nan")),
                   "pass": bool(okb)})

    # ---- T0c + the reading window R: sign-change scan then bisect -------
    # the window is ~1e-8 wide (the read is a 10th power: coef rises from
    # 0.05 at 0.196 to 3 at 0.20013), so a grid can only DETECT the band
    # by bracketing its edges in sign change, never by landing inside it
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    xs = [i * 1e-3 for i in range(1, 786)]
    print("[V326] coarse scan of {} angles starting".format(len(xs)),
          file=sys.stderr, flush=True)
    cv = []
    for n, x in enumerate(xs):
        cv.append(coef10(x, LOOP, mach))
        if n % 100 == 99:
            print("[V326] scan {}/{}".format(n + 1, len(xs)),
                  file=sys.stderr, flush=True)
    edges = []
    for tgt, edge in ((COEF_OBS - BAND_C, "lo"), (COEF_OBS + BAND_C, "hi")):
        for i in range(len(xs) - 1):
            if (cv[i] - tgt) * (cv[i + 1] - tgt) < 0.0:
                print("[V326] edge bracket x={}..{} tgt={} bisecting".format(
                    xs[i], xs[i + 1], edge), file=sys.stderr, flush=True)
                edges.append((bisect(xs[i], xs[i + 1], tgt, LOOP, mach),
                              edge))
    edges.sort()
    print("[V326] {} band-edge crossings found".format(len(edges)),
          file=sys.stderr, flush=True)
    # walk the edge-sorted line: between consecutive edges coef keeps its
    # inside/outside status, so test midpoints and merge adjacent inside
    # intervals - handles windows entered from above OR below
    ex = [x for x, _ in edges]
    refined, i = [], 0
    while i < len(ex) - 1:
        if abs(coef10(0.5 * (ex[i] + ex[i + 1]), LOOP, mach) - COEF_OBS) \
                <= BAND_C:
            j = i
            while j + 2 <= len(ex) - 1 and \
                    abs(coef10(0.5 * (ex[j + 1] + ex[j + 2]), LOOP, mach)
                        - COEF_OBS) <= BAND_C:
                j += 1
            refined.append({"x_lo": ex[i], "x_hi": ex[j + 1],
                            "width": ex[j + 1] - ex[i]})
            i = j + 2
        else:
            i += 1
    def in_R(x):
        return any(c["x_lo"] <= x <= c["x_hi"] for c in refined)
    data["window_R"] = {"components": refined,
                        "g_obs": g_obs,
                        "g_obs_in_R": bool(in_R(g_obs))}
    checks.append({"name": "T0c {} reading window component(s) (widths "
                   "{:.3e} ... {:.3e}); the V206 g_obs = {:.7f} leaf lies "
                   "inside one (inversion fact reproduced)".format(
                       len(refined),
                       refined[0]["width"] if refined else float("nan"),
                       refined[-1]["width"] if refined else float("nan"),
                       g_obs),
                   "pass": bool(len(refined) >= 1
                                and data["window_R"]["g_obs_in_R"])})

    # ---- the closure set Cl: V315's grid census on the identity step ----
    grid = []
    for b in range(2, DEN_MAX + 1):
        for a in range(1, b):
            g = math.pi * a / b
            if math.gcd(a, b) != 1:
                continue
            grid.append((a, b, g))
    cl_rows = {}
    print("[V326] closure census over {} grid angles".format(len(grid)),
          file=sys.stderr, flush=True)
    for (a, b, g) in grid:
        row = closure_row(step_at(g))
        if row:
            cl_rows["pi*{}/{}".format(a, b)] = {
                "g": g, "k": row["k"], "dev": row["dev"]}
    data["closure_set_Cl"] = {"grid_points": len(grid),
                              "closing": cl_rows}
    checks.append({"name": "grid census ran exactly {} reduced pi*a/b "
                   "points (b <= {}); {} close at k <= {}".format(
                       len(grid), DEN_MAX, len(cl_rows), K_MAX),
                   "pass": bool(len(grid) == 79)})

    # ---- the bridge census: every closure angle against R ----------------
    v314 = json.load(open(os.path.join(
        exp("V314_four_scale_ledger"), "results",
        "v314_four_scale_ledger.json")))
    cross_by_clock = {}
    for row in v314["data"]["ledger"]:
        r = row.get("row", row)
        ck = str(r.get("clock", ""))
        if "pi/4" in ck and "scan" in ck:
            cross_by_clock["pi/4"] = float(r["coef"])
        if "pi/2" in ck and "scan" in ck:
            cross_by_clock["pi/2"] = float(r["coef"])
    F_cert = float(art("V321_crossing_resummation",
                       "data.certified_leaves.F_208"))
    cert_coef = read_coef(F_cert, LOOP)
    bridge_rows = []
    for key, row in sorted(cl_rows.items()):
        g = row["g"]
        rung10 = coef10(g, LOOP, mach)
        inside = in_R(g)
        bridge_rows.append({
            "clock": key, "g": g, "k": row["k"],
            "coef_rung10": rung10,
            "residual_rung10_bands": abs(rung10 - COEF_OBS) / BAND_C,
            "crossing_read": cross_by_clock.get(
                {"pi*1/4": "pi/4", "pi*1/2": "pi/2"}.get(key, key)),
            "in_R": bool(inside)})
    data["bridge_census"] = bridge_rows
    n_bridge = sum(1 for r in bridge_rows if r["in_R"])
    min_dist = min(r["residual_rung10_bands"] for r in bridge_rows)
    # certified crossing distance at the pi/4 closure clock
    cert_bands = abs(cert_coef - COEF_OBS) / BAND_C

    # ---- grid-window separation (books the 784-variant + V247 unions) ----
    sep = (min(min(abs(g - c["x_lo"]), abs(g - c["x_hi"]))
               for c in refined for (a, b, g) in grid)
           if refined else float("nan"))
    data["grid_window_separation"] = {
        "min_abs_gap": sep,
        "nearest_grid_angle": (min(
            (min(abs(g - c["x_lo"]), abs(g - c["x_hi"])),
             "pi*{}/{}".format(a, b))
            for c in refined for (a, b, g) in grid)[1]
            if refined else None)}

    # ---- the c-class table (lemma: all members share the verdict) --------
    lemma2pi = closure_row(step_at(2.0 * math.pi * G5))
    grid25 = closure_row(step_at(math.pi * 2.0 / 5.0))
    lemma_ok = (bool(lemma2pi) == bool(grid25)
                and (not lemma2pi or lemma2pi["k"] == grid25["k"]))
    data["dictionary_class"] = {
        "members": ["c = {}".format(c) for c in CS],
        "bridge_invariant_in_c": True,
        "lemma_spot_2pi": {"closure_at_2pi*1_5": lemma2pi,
                           "grid_row_pi2_5": grid25,
                           "equal": bool(lemma_ok)}}
    checks.append({"name": "lemma spot c = 2*pi: the dictionary closes "
                   "at g = 1/5 exactly when the grid angle pi*2/5 closes "
                   "({} vs {})".format(
                       lemma2pi["k"] if lemma2pi else None,
                       grid25["k"] if grid25 else None),
                   "pass": bool(lemma_ok)})

    # ---- verdict ----------------------------------------------------------
    if n_bridge > 0:
        status = "DICTIONARY-BRIDGE-FOUND"
    else:
        status = "DICTIONARY-NO-BRIDGE"
    data["status"] = status
    data["coefficient_observed"] = COEF_OBS
    data["band_coef"] = BAND_C
    data["certified_pi4_crossing_bands"] = cert_bands
    data["closest_closure_rung10_bands"] = min_dist
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []

    checks.append({"name": "bridge census: {} of {} closure angles lie "
                   "inside R; closest rung-10 read {:.1f} band-units off; "
                   "the certified pi/4 crossing read {:.1f} band-units "
                   "off; every grid angle >= {:.2e} away from the "
                   "window".format(n_bridge, len(bridge_rows), min_dist,
                                   cert_bands, sep),
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    out = {"id": "V326_phase_dictionary_bridge", "tier": 1,
           "property": "phase-dictionary bridge between the closure set "
           "and the coefficient reading window (tier 1, hard, O4/O8)",
           "residual": True, "programme": "theory-hunt V326",
           "branch": status, "data": data, "checks": checks,
           "status": status, "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v326_phase_dictionary_bridge.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V326 done] {}".format(status))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
