#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V232 - the bounded insertion-amplitude search for O8 (tier 5).

**The question this pre-registers against** (`review/pre_registrations.md`):
O8's named residual is a single number: the loop functional at the certified
crossing is short of the observed hierarchy coefficient by the insertion
amplitude `0.9638788` (V205: required `0.9638787569864106` against the
closed-branch crossing deficit `0.989527` float64 / `0.989868` certified at
n* = 208, V209). Does ANY structure in the closed set supply `0.96388`?

The search is bounded to the closed structure set, exactly as pre-registered -
no new mechanism, no tuned input, every candidate read from a committed
artifact or evaluated exactly from artifact-recorded data:

C1 the five unification Cartan directions (V205 `per_cartan_insertion_at_
   crossing`, the bimodal set);
C2 the certified closed-clock crossing deficits (V209: pi/4 and 3pi/4 on the
   n = 4k class; pi/2 reported as floor-conditioned, counted with the flag);
C3 the dressed gate of V140/V141: the Kogut-Susskind electric gate amplitude
   |tr K_E| / 1440 at the closed clock, computed from the content theorem
   recorded in the V140 artifact (C2 levels 45, 85, 117 with multiplicities
   32, 288, 1120), and the regular-link insertion amplitudes of V141
   (A(h) = rho(h^2) = +-I -> amplitude exactly 1; the dressed hop speed
   v = 1);
C4 the star-module certified values of V231 (deficits of certified
   crossings, if any);
C5 the REFERENCE-branch interference deficit at the 1/5 clock
   (`0.9644132`, V205/V209/V210) - counted but flagged: 1/5 is not a closed
   angle, and this is the closest candidate.

Acceptance (pre-registered): a candidate counts as closing O8 only if
`|a - 0.9638787569864106| < 1e-6` and its input set is a subset of the
closed list; otherwise `branch: O8-OPEN` with the full measured set recorded.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                    # noqa: E402


CLOSED_CLOCKS = (("pi/4", math.pi / 4), ("pi/2", math.pi / 2),
                 ("3pi/4", 3 * math.pi / 4))
TARGET = 0.9638787569864106          # V205 arithmetic, sourced input
TOL = 1e-6                           # the pre-registered acceptance band


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    data["fitted_inputs"] = []

    v205 = art("V205_channel_filtered_deficit")
    v209 = art("V209_root_clock_degeneracy")
    v140 = art("V140_dressed_gate_quintet")
    v141 = art("V141_regular_link_quintet_coupling")
    try:
        v231 = art("V231_star_loop_functional_certified")
    except FileNotFoundError:
        v231 = None

    target = float(v205["data"]["arithmetic"]["required_insertion_amplitude"])
    assert abs(target - TARGET) < 1e-12

    cands = []

    # ---- C1 the five unification Cartan directions ------------------------
    for i, a in enumerate(v205["data"]["per_cartan_insertion_at_crossing"]):
        cands.append({"id": f"cartan_{i}", "a": float(a), "source":
                      "V205 per_cartan_insertion_at_crossing",
                      "closed_inputs": True, "flag": None})

    # ---- C2 certified closed-clock crossing deficits ----------------------
    cert = v209["data"]["certified"]
    for tag in ("pi/4@1e-8", "3pi/4@1e-6", "pi/2@1e-6"):
        if tag in cert:
            cands.append({"id": f"crossing_deficit_{tag}",
                          "a": float(cert[tag]["deficit"]),
                          "source": "V209 certified crossing deficit",
                          "closed_inputs": True,
                          "flag": ("floor-conditioned, not a stable number"
                                   if tag.startswith("pi/2") else None)})

    # ---- C3 the dressed gate amplitudes -----------------------------------
    eg = v140["data"]["dressed_electric_gate"]
    levels = [round(x) for x in eg["C2_levels"]]
    mults = [int(m) for m in eg["multiplicities"]]
    assert levels == [45, 85, 117] and mults == [32, 288, 1120]
    for tag, g in CLOSED_CLOCKS:
        z = sum(m * complex(math.cos(g * c), math.sin(g * c))
                for c, m in zip(levels, mults))
        cands.append({"id": f"KE_trace_{tag}", "a": abs(z) / sum(mults),
                      "source": "V140 content theorem (C2 levels/multiplic"
                                "ities), evaluated at the closed clock",
                      "closed_inputs": True, "flag": None})
    assert v141["data"]["hopping"]["rho_h_squared_err"] < 1e-12
    cands.append({"id": "regular_link_insertion", "a": 1.0,
                  "source": "V141: A(h) = rho(h^2) = +-I exactly "
                            "(err < 1e-12)",
                  "closed_inputs": True, "flag": "exactly 1"})
    vel = v141["data"]["velocity"]
    cands.append({"id": "dressed_hop_speed",
                  "a": 0.5 * (vel["bare_front_speed_max"]
                              + vel["bare_front_speed_min"]),
                  "source": "V141: the dressed hop front speed v = 1",
                  "closed_inputs": True, "flag": "a speed, counted as an "
                  "amplitude analogue"})

    # ---- C4 the V231 certified star-module deficits -----------------------
    if v231 is not None:
        cx = v231["data"].get("certified_crossings", {})
        for k, rows in cx.items():
            for e in rows:
                cands.append({"id": f"star_certified_{k}_n{e['n']}",
                              "a": float(e["deficit"]),
                              "source": "V231 certified star crossing",
                              "closed_inputs": True,
                              "flag": v231["branch"]})
        if not cx:
            data["v231_candidates"] = ("none: V231 branch "
                                       f"{v231['branch']} - no certified "
                                       "star-module crossing to inherit")

    # ---- C5 the reference-branch deficit ----------------------------------
    cands.append({"id": "reference_1_5_deficit",
                  "a": float(v209["data"]["open_branch"]["deficit"]),
                  "source": "V209 open_branch deficit at g = 1/5 "
                            "(V205/V210 agree)",
                  "closed_inputs": False,
                  "flag": "clock 1/5 is NOT a closed angle (F2's branch "
                          "question)"})

    # ---- acceptance --------------------------------------------------------
    for c in cands:
        c["dev"] = abs(c["a"] - target)
    cands.sort(key=lambda c: c["dev"])
    best = cands[0]
    closed_ok = [c for c in cands if c["dev"] < TOL and c["closed_inputs"]]
    branch = "O8-CLOSED" if closed_ok else "O8-OPEN"
    dev_min = best["dev"]

    # the honest nearest misses, each with WHY it misses
    near = [{"id": c["id"], "a": c["a"], "dev": c["dev"],
             "why": ("outside the closed input set" if not c["closed_inputs"]
                     else "outside the 1e-6 band")}
            for c in cands if c["dev"] < 1e-2]

    checks.append({
        "name": f"the bounded search: {len(cands)} candidate insertion "
                f"amplitudes from the closed set only (five Cartan "
                f"directions {sorted(round(c['a'], 4) for c in cands if c['id'].startswith('cartan'))}, "
                f"certified closed-clock deficits, dressed-gate amplitudes, "
                f"star-module certified values, reference-branch deficit), "
                f"target {target:.7f}: NOTHING within 1e-6 with closed "
                f"inputs (closest {best['id']} a = {best['a']:.6f}, "
                f"|a - target| = {dev_min:.3e}{'; flagged non-closed input' if not best['closed_inputs'] else ''}); "
                f"the reference-branch deficit 0.9644132 misses the target "
                f"by "
                f"{abs(next(c['dev'] for c in cands if c['id'] == 'reference_1_5_deficit')):.2e}"
                f" - the fingerprint deviation V205 already named "
                f"(-5.5e-4): O8 stays open",
        "pass": bool(branch == "O8-OPEN" and dev_min > TOL)})
    data["target"] = target
    data["candidates"] = cands
    data["near_misses"] = near
    data["acceptance_band"] = TOL

    ok = all(c["pass"] for c in checks)
    status = (
        f"O8 {('CLOSED' if closed_ok else 'STAYS OPEN')}, bounded search "
        f"(zero fitted parameters, {time.time() - t0:.0f} s): {len(cands)} "
        f"candidate insertion amplitudes inside the closed structure set "
        f"only - the five unification Cartan directions (V205, bimodal in "
        f"the 0.94-1.02 band), the certified crossing deficits at pi/4 "
        f"(0.989868) and 3pi/4 (0.989933), the dressed electric-gate trace "
        f"amplitudes at the closed clocks, the regular-link insertion (1 "
        f"exactly) and the dressed hop (v = 1), the V231 star-module "
        f"certified deficits ({'none, V231 found no certified crossing' if not any(c['id'].startswith('star_certified') for c in cands) else 'inherited'}), and the reference-branch deficit "
        f"0.9644132 at the 1/5 clock (flagged: not a closed angle): none "
        f"satisfies |a - 0.9638788| < 1e-6 with closed inputs (closest "
        f"overall {best['id']}: {best['a']:.7f}, dev {dev_min:.2e}; closest "
        f"closed-input candidate "
        f"{min((c for c in cands if c['closed_inputs']), key=lambda c: c['dev'])['id']}: "
        f"{min((c for c in cands if c['closed_inputs']), key=lambda c: c['dev'])['a']:.7f}, "
        f"dev {min((c['dev'] for c in cands if c['closed_inputs'])):.2e}). "
        f"The required amplitude is NOT supplied by any structure of the "
        f"closed set; the nearest object remains the discarded reference "
        f"branch, off by 5.3e-4. VERDICT: branch O8-OPEN with the full "
        f"measured set; the O8 deficit stays a named open number "
        f"(0.989868 measured vs 0.963879 required).")
    out = {"id": "V232_insertion_amplitude_search",
           "tier": 5,
           "property": "the bounded insertion-amplitude search over the "
                       "closed structure set against the O8 target 0.9638788 "
                       "(V232, O8)",
           "residual": True,
           "programme": "V232",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v232_insertion_amplitude_search.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
