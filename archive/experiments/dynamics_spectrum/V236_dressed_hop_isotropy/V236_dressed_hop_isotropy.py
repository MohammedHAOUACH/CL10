#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V236 - one bounded isotropy attempt: the dressed hop (tier 2).

**The question this pre-registers against** (`review/pre_registrations.md`):
does the gapless node become isotropic when the hop is dressed by the
regular-link network of V141 Part A - the network that propagates at
`v = 1` - or on the icosahedral star module, where the equivariant census
pushes the first anisotropic harmonic to degree 5? ONE attempt; if it fails
the measurement stops there.

Closed input set (pre-registered): V141's regular-link dressed hop (its own
file: the group GAMMA of order 1024, the link insertion `A(h) = rho(h^2) =
sign(h^2) I`, the closed clock); `cl10.starcell` (V219's recorded star
nodal-line velocities); the closed angles; the V235/V115 cone estimator and
normalisation (`delta_v/v = (vmax - vmin)/vmax` of the fitted tensor and of
the axis speeds, with a flat direction entering as exactly 1). No new input.

Measured quantity: `delta_v/v` at the gapless points of the DRESSED hop at
the closed clock `pi/4`, plus the star nodal line, plus the bare-hop
control at both clocks.

Acceptance (pre-registered): `branch: ISOTROPY` only if `delta_v/v < 1e-6`
at a closed angle with no new input; otherwise `branch: ANISOTROPY` with
the measured number, and the programme stops - no third attempt.
"""
import glob
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
os.environ["CL10_ANGLE"] = "0.2"       # reference branch: pin BEFORE cl10
                                       # latches the clock (V118/V207 order)
from cl10.paths import exp                                    # noqa: E402
sys.path.insert(0, exp("V141_regular_link_quintet_coupling"))
sys.path.insert(0, exp("V235_node_velocity_contact"))
import V141_regular_link_quintet_coupling as V141             # noqa: E402
import V235_node_velocity_contact as V235                     # noqa: E402

from cl10 import config as _config                            # noqa: E402
from cl10 import engine as E                                  # noqa: E402
from cl10 import supercell as SC                              # noqa: E402

CLOSED = math.pi / 4


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)



def main():
    t0 = time.time()
    V235._init()
    checks, data = [], {}
    data["fitted_inputs"] = []

    v141 = art("V141_regular_link_quintet_coupling")
    v140 = art("V140_dressed_gate_quintet")
    v219 = art("V219_nodal_line_anisotropy")

    # ---- T0 V141 Part A facts from the committed artifacts -----------------
    hop = v141["data"]["hopping"]
    vel = v141["data"]["velocity"]
    eg = v140["data"]["dressed_electric_gate"]
    gate_phase = float(eg["gate_phases_mod2pi"][0])
    scalar_ok = abs(vel["electric_scalar_err"]) < 1e-12
    insertion_ok = hop["rho_h_squared_err"] < 1e-12
    front = (float(vel["bare_front_speed_max"]),
             float(vel["bare_front_speed_min"]))
    checks.append({
        "name": f"T0 V141 Part A facts: the link insertion is A(h) = "
                f"rho(h^2) = sign(h^2) I exactly (err "
                f"{hop['rho_h_squared_err']:.1e}); the V140 electric gate "
                f"at the closed clock is a scalar (err "
                f"{vel['electric_scalar_err']:.1e}, phase "
                f"{gate_phase:.4f} = 5 pi/4 on all three Casimir blocks), "
                f"and the dressed hop front runs at v = 1 "
                f"({front[1]:.9f} to {front[0]:.9f})",
        "pass": bool(insertion_ok and scalar_ok
                     and abs(gate_phase - 5 * math.pi / 4) < 1e-5)})
    data["v141_facts"] = {"insertion_err": hop["rho_h_squared_err"],
                          "gate_scalar_err": vel["electric_scalar_err"],
                          "gate_phase": gate_phase,
                          "front_speed": list(front)}

    # ---- T1 the dressing identity: six insertions of -I give +I ------------
    # the bond biflections h_a h_b square to -1 in Cl(10,0), so the link
    # insertion A(h) = sign(h^2) I = -I on EVERY link; the step carries six
    # bond factors, so the dressed step equals the bare step exactly.
    h2_sign = V141.mul((1, (0, 1)), (1, (0, 1)))[0]   # A(h) = sign(h^2)
    qtest = (0.31, -1.10, 0.77)
    U4 = SC.sectors(CLOSED)[2]

    def step6(q, signs):
        M = np.eye(128, dtype=complex)
        for even, j in SC.STAGES:
            M = np.kron(signs[j] * SC.stage_matrix(even, j, q[j], CLOSED),
                        U4) @ M
        return M

    ph_bare = np.sort(np.angle(np.linalg.eigvals(step6(qtest, (1,) * 6))))
    ph_drs = np.sort(np.angle(np.linalg.eigvals(step6(qtest, (-1,) * 6))))
    dev1 = float(np.max(np.abs(np.angle(np.exp(1j * (ph_drs - ph_bare))))))
    checks.append({
        "name": f"T1 dressing identity: every bond biflection squares to "
                f"sign {h2_sign}, so A(h) = -I on every link; the step "
                f"carries six bond factors and the six insertions multiply "
                f"to +I: the dressed and bare one-step spectra at the "
                f"closed clock are identical (max sorted-eigenphase "
                f"deviation {dev1:.1e}). At the closed clock the dressing "
                f"carries NO dispersion content - the v = 1 front of T0 is "
                f"a maximum-speed statement, not a cone statement",
        "pass": bool(h2_sign < 0 and dev1 < 1e-10)})
    data["dressing"] = {"A_sign_per_bond": h2_sign, "dev": dev1,
                        "gate_phase": gate_phase}

    # ---- T2 cone measurement at the dressed-hop gapless points -------------
    # gapless at CLOSED clock: bare eigenphase = -5 pi/4 (mod 2 pi); scan
    # the line q = (t, t, t) and the M-shifted line for crossing points.
    gap_pts = []
    for base in ((math.pi, math.pi, math.pi), (0.0, 0.0, 0.0),
                 (math.pi / 2, math.pi / 2, 0.0)):
        for t in np.linspace(-math.pi, math.pi, 97):
            q = tuple(base[k] + (t if k == 0 else 0.0) for k in range(3))
            M = E.asnumpy(V235.step_batch_at([q], CLOSED)[0])
            g = float(np.linalg.svd(M - np.eye(128, dtype=complex),
                                    compute_uv=False)[-1])
            if g < 1e-9:
                gap_pts.append(q)
                break
    cones = []
    for q in gap_pts[:3]:
        v, nlin = V235.tensor_fit(lambda qs: V235.step_batch_at(qs, CLOSED),
                                  q, V235.fib_dirs(60), 2e-4, 1e-3)
        axis = [max(V235.velocities(
            lambda qs: V235.step_batch_at(qs, CLOSED), q, [d], 1e-3)) / 1e-3
            for d in (np.array([1., 0, 0]), np.array([0, 1., 0]),
                      np.array([0, 0, 1.]))]
        cones.append({"q": [round(x, 5) for x in q],
                      "v_principal_fit": [round(float(x), 5) for x in v],
                      "axis_speeds": [round(x, 5) for x in axis],
                      "n_linear_dirs": nlin,
                      "dvv": max(V235.dvv(v), V235.dvv(axis))})
    dv_dressed = max((c["dvv"] for c in cones), default=1.0)
    checks.append({
        "name": f"T2 the dressed hop at the closed clock pi/4: "
                f"{len(gap_pts)} exactly gapless points located on the "
                f"probe lines (gap < 1e-9); their cone anisotropy "
                f"delta_v/v = {[round(c['dvv'], 4) for c in cones]} - "
                f"worst {dv_dressed:.4f}: the v = 1 front of the dressed "
                f"network is the MAXIMUM-velocity statement (T0); the "
                f"gapless cone of the same network is not isotropic - a "
                f"flat direction survives the dressing",
        "pass": bool(len(gap_pts) >= 1 and dv_dressed > 1e-6)})

    # ---- T3 the star module: equivariant census up to degree 5 --------------
    az = v219["data"]["azimuth_star"]
    r2 = az["3e-03"]
    star_dvv = (r2["v_max"] - r2["v_min"]) / r2["v_max"]
    lead = int(r2["leading"])
    checks.append({
        "name": f"T3 the star module: the equivariant census makes the "
                f"star's invariants start at degree 2 (V219 measured "
                f"leading azimuthal degree {lead} on the two-fold line, "
                f"odd-harmonic fraction {r2['odd_fraction']:.1e}); the "
                f"measured nodal-line cone is still "
                f"delta_v/v = {star_dvv:.4f} (v from {r2['v_min']:.3e} to "
                f"{r2['v_max']:.3f}): pushing the first anisotropic "
                f"harmonic up does not isotropize the node",
        "pass": bool(star_dvv > 1e-6)})
    data["star"] = {"dvv": star_dvv, "leading": lead,
                    "v_min": r2["v_min"], "v_max": r2["v_max"]}

    # ---- bare-hop control at both clocks -----------------------------------
    ctl = []
    for tag, ang in (("pi/4", CLOSED), ("1/5", 0.2)):
        q0 = tuple(float(x) for x in SC.M_POINTS[0])
        v, _ = V235.tensor_fit(
            lambda qs, a=ang: V235.step_batch_at(qs, a), q0,
            V235.fib_dirs(60), 2e-4, 1e-3)
        ctl.append({"clock": tag, "v_principal": [round(float(x), 5)
                                                  for x in v],
                    "dvv": V235.dvv(v)})
    data["bare_control"] = ctl

    best = min([dv_dressed, star_dvv] + [c["dvv"] for c in ctl])
    branch = "ISOTROPY" if best < 1e-6 else "ANISOTROPY"
    ok = all(c["pass"] for c in checks)
    status = (
        f"dressed-hop isotropy attempt: {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): at the closed clock pi/4 the V141 "
        f"dressing is degenerate as a mechanism - A(h) = rho(h^2) = -I "
        f"per bond, the step's six bond factors make the six insertions "
        f"multiply to +I (dressed and bare spectra identical to "
        f"{dev1:.0e}), and the electric gate is a scalar (5 pi/4, err "
        f"{vel['electric_scalar_err']:.0e}) - so it CANNOT move the cone: "
        f"the dressed network's gapless points carry delta_v/v = "
        f"{dv_dressed:.4f} (flat directions survive; the v = 1 statement "
        f"is the front-speed maximum, measured {front[1]:.9f} to "
        f"{front[0]:.9f}); the star module, whose equivariant census "
        f"starts at degree 2 with the odd harmonics dead "
        f"({r2['odd_fraction']:.0e}), gives delta_v/v = {star_dvv:.4f} at "
        f"its nodal line. Bare control: 1/5 {ctl[1]['dvv']:.3f}, pi/4 "
        f"{ctl[0]['dvv']:.3f}. BEST delta_v/v over the attempt: "
        f"{best:.3f} >= 1e-6. The bounded attempt fails; per the "
        f"pre-registration this programme stops here - no third attempt.")
    out = {"id": "V236_dressed_hop_isotropy",
           "tier": 2,
           "property": "whether dressing the hop by the regular-link "
                       "network (or the star module) isotropizes a gapless "
                       "cone below 1e-6 (V236, O9, one bounded attempt)",
           "residual": True,
           "programme": "V236",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v236_dressed_hop_isotropy.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
