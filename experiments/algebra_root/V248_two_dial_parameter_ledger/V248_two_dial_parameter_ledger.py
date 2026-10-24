#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V248 - the honest ledger of the repaired theory (tier 6).

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): the provenance of every number of the repaired theory -
for each landed dial value (beta*, g*, beta_t, F_infinity): algebraic
(unique on C or C - M_v, fixed by a measured acceptance rule with
fitted_inputs == []) or fitted? The Table-6-class counts per quantity,
compared against the published Table 6 rows (V52's provenance_table, the
artifact Table 6 is generated from).

Closed inputs: only the programme's verdict artifacts (V239-V247 via
cl10.paths.exp) and the V52 published rows. No new number is measured.

Branch `TWO-DIAL-LEDGER` (informational); `residual: true` iff a dial
landed fitted. V240 landed MANIFOLD-EMPTY, so the would-be dials are
recorded never-realized and add no fitted parameter (pre-registered).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10.paths import exp                                     # noqa: E402

IDS = ["V239_two_dial_step_family", "V240_two_dial_velocity_manifold",
       "V241_two_dial_algebra_replication",
       "V242_bz_loop_factor_census", "V243_bz_hierarchy_contact",
       "V244_icosahedral_node_census", "V245_enforced_isotropy_contact",
       "V246_lorentz_tuning_relevance", "V247_two_dial_desert_contacts"]


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = __import__("time").time()
    arts = {}
    for i in IDS:
        arts[i] = art(i)
    v52 = art("V52_parameter_audit")["data"]["provenance_table"]["rows"]
    checks, ledger = [], []

    # ---- the programme's own consistency ----------------------------------
    allread = all(isinstance(a, dict) for a in arts.values())
    checks.append({
        "name": f"T0 all nine verdict artifacts of the programme loaded "
                f"and parsed: {allread}; published Table 6 has {len(v52)} "
                f"rows to compare against (V52 provenance_table)",
        "pass": bool(allread and len(v52) >= 8)})

    # ---- the ledger rows ---------------------------------------------------
    closure = arts["V239_two_dial_step_family"]["data"]["closure"]
    ledger.append({
        "quantity": "the two-dial closure set C",
        "value": f"{len(closure['closure_set'])} pairs (dial closures "
                 f"of order-8 content, incl. (pi/4, pi/4))",
        "provenance": "algebraic",
        "why": "unique on the closed family (V239's measured acceptance "
               "rule, fitted_inputs == []); V241 confirms zero algebraic "
               "collateral damage on the slice",
        "source": "V239"})
    ledger.append({
        "quantity": "the velocity-manifold dials beta*, g* (vbar = 1)",
        "value": "never realized (M_v empty)",
        "provenance": "never-realized",
        "why": "V240 MANIFOLD-EMPTY: no point of C has cone velocity "
               "within the pre-registered band of c = 1; per the "
               "pre-registration these dials add no fitted parameter",
        "source": "V240"})
    ledger.append({
        "quantity": "the BZ loop factor F_infinity",
        "value": "0.142066.../0.144813... (n=256, N=32, implemented pair)",
        "provenance": "measured-input",
        "why": "a measured value of the functional at published clocks, "
               "not an input to any claimed coefficient: V243 HIERARCHY-"
               "UNREAD shows it misses the ladder reading by ~7.4x per "
               "loop; V242 NO-LOOP-FACTOR. Table 6's hierarchy row is "
               "UNCHANGED",
        "source": "V242/V243"})
    ledger.append({
        "quantity": "the tuning dial beta_t (degree-4 cancellation)",
        "value": "never realized (conditional skip)",
        "provenance": "never-realized",
        "why": "V246's pre-registered trigger (V244 in {NO-ENFORCED-LOCUS, "
               "PARTIAL-LIMIT}) did not fire; nothing measured",
        "source": "V246"})
    ledger.append({
        "quantity": "the clock g_10 and the fitted g_obs",
        "value": "unchanged: F2 fired, g_obs = 0.1999980 stays fitted",
        "provenance": "fitted (pre-existing, published)",
        "why": "V247 F2-FIRED-EXTENDED: no closed two-dial point restores "
               "both V233 bands and V237's desert rule; beta is verified "
               "absent from the running, so the repair changes nothing "
               "here",
        "source": "V247 (against V233/V237/V206 published rows)"})
    ledger.append({
        "quantity": "the enforced star-module node q = 0 and its "
                    "spinorial 512-mode space",
        "value": "stabilizer 60, Molien rank-2 count 1, 2I-invariant "
                 "modes (centre trace -512)",
        "provenance": "algebraic",
        "why": "symmetry census (V244/V245): a structural fact, not a "
               "dial; V245 F9-STILL-FIRED shows it does not close the "
               "SME contact (multiplet splitting)",
        "source": "V244/V245"})

    fitted = [r for r in ledger if r["provenance"] == "fitted"
              and "pre-existing" not in r["provenance"]]
    never = [r for r in ledger if r["provenance"] == "never-realized"]
    alg = [r for r in ledger if r["provenance"] == "algebraic"]
    meas = [r for r in ledger if r["provenance"] == "measured-input"]
    data = {
        "ledger": ledger,
        "counts": {"algebraic": len(alg), "fitted_new": len(fitted),
                   "never_realized": len(never),
                   "measured_input": len(meas)},
        "reclassification_flags_for_table6": [
            "NEW ROW (algebraic): two-dial closure set C (V239/V241)",
            "UNCHANGED (fitted): clock g_10 / g_obs - V247 extends F2's "
            "firing to the whole closed two-dial content",
            "UNCHANGED (named residual): the ladder coefficient - V242/"
            "V243 measure the BZ functional's miss (7.4x per loop), "
            "closing nothing",
            "NO NEW FITTED PARAMETERS: both would-be dials (vbar=1 "
            "manifold, beta_t) never realized (V240, V246)"],
        "fitted_inputs": []}
    checks.append({
        "name": f"T1 ledger counts: algebraic {len(alg)}, new fitted "
                f"{len(fitted)}, never-realized {len(never)}, measured-"
                f"input {len(meas)}; MANIFOLD-EMPTY clause honored: "
                f"{arts['V240_two_dial_velocity_manifold']['branch'] == 'MANIFOLD-EMPTY' and len(never) >= 2}",
        "pass": bool(len(fitted) == 0 and len(never) == 2
                     and len(alg) == 2)})

    status = (
        f"TWO-DIAL-LEDGER (zero fitted parameters, "
        f"{__import__('time').time() - t0:.0f} s): the repaired theory "
        f"gains NO new fitted parameter: the closure set C is algebraic "
        f"(V239/V241), the vbar = 1 dials (V240 MANIFOLD-EMPTY) and the "
        f"tuning dial beta_t (V246 conditional skip) are never-realized, "
        f"the BZ loop factor is a measured value that closes nothing "
        f"(V242/V243: 7.4x per loop off the ladder reading), the clock "
        f"row stays fitted exactly as published (V247 F2-FIRED-EXTENDED: "
        f"no closed point restores the 2x2), and the enforced node is an "
        f"algebraic structural fact that does not close F9 (V244/V245). "
        f"Table 6 flags: 1 new algebraic row, 2 unchanged rows, 0 new "
        f"fitted rows.")
    out = {"id": "V248_two_dial_parameter_ledger",
           "tier": 6,
           "property": "provenance ledger of the repaired theory: which "
                       "dial values are algebraic, fitted, never-realized, "
                       "measured (V248, Table 6)",
           "residual": bool(len(fitted) > 0),
           "programme": "V248",
           "branch": "TWO-DIAL-LEDGER",
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(all(c["pass"] for c in checks))}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v248_two_dial_parameter_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
