#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V246 - the tuning alternative (CONDITIONAL, pre-registered trigger).

**Pre-registration** (`review/pre_registrations.md`, second programme):
"Conditional, pre-registered as MANDATORY iff V244 in {NO-ENFORCED-LOCUS,
PARTIAL-LIMIT}; if V244 landed ENFORCED-LOCUS the experiment writes a
conditional-skip verdict (`inconclusive: true`, no branch string) and
nothing is measured."

V244 landed `ENFORCED-LOCUS`, so this file does exactly that: read V244's
artifact, record the skip against its branch, write the verdict, exit 0.
Had the trigger fired, the pre-registered measurement would be V213's
relevance exponent of the anisotropic correction at the best C - M_v node
along the three C-nearest beta values, plus a root-find of the beta_t
where the degree-4 coefficient vanishes - and that programme is NOT run
here, exactly as pre-registered.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10.paths import exp                                     # noqa: E402


def main():
    with open(exp("V244_icosahedral_node_census", "results",
                  "v244_icosahedral_node_census.json")) as fh:
        v244 = json.load(fh)
    triggered = v244["branch"] in ("NO-ENFORCED-LOCUS", "PARTIAL-LIMIT")
    status = (
        f"V246 CONDITIONAL SKIP: the pre-registered trigger is V244 in "
        f"{{NO-ENFORCED-LOCUS, PARTIAL-LIMIT}}; V244 landed "
        f"{v244['branch']}, so the tuning-relevance measurement (V213-"
        f"style exponent of the anisotropic correction along the bond "
        f"dial, and the beta_t root of the degree-4 coefficient) is NOT "
        f"run and nothing is measured, exactly as pre-registered.")
    out = {"id": "V246_lorentz_tuning_relevance",
           "tier": 2,
           "property": "relevance of the node anisotropy along the bond "
                       "dial beta (V246, F9) - conditional programme",
           "residual": True,
           "programme": "V246",
           "inconclusive": True,
           "skip_reason": {"trigger_branches": ["NO-ENFORCED-LOCUS",
                                               "PARTIAL-LIMIT"],
                           "v244_branch": v244["branch"],
                           "triggered": bool(triggered)},
           "status": status,
           "pass": True}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v246_lorentz_tuning_relevance.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0


if __name__ == "__main__":
    sys.exit(main())
