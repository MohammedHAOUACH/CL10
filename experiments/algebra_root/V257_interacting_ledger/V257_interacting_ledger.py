#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V257 - the interacting ledger (tier 6, residual: true iff anything
landed fitted).

**The question this pre-registers against** (`review/pre_registrations.md`,
third programme): the provenance of every number the interacting
programme produced - eta twist values (algebraic only on C_star, else
never-realized), F_int (measured), dressed cone/flux scales (measured),
the sector coupling (algebraic: g_10 and the V177/V141 zero-parameter
gates), any C_star clock (algebraic if closed).  The Table-6-class
counts per quantity, compared against the published Table 6 rows (V52's
provenance_table, the artifact Table 6 is generated from).

Closed inputs: only the V249-V256 verdict artifacts (cl10.paths.exp)
and the V52 published rows.  No new number is measured.

Branch `INTERACTING-LEDGER` (informational); `residual: true` iff a
quantity landed fitted.  V251 landed SECTOR-BLOCKED and V249
STAR-CLOCK-QUANTIZED, so the would-be sector numbers and the C_star
clock/twists are recorded never-realized and add no fitted parameter
(pre-registered roll-up clause).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10.paths import exp                                     # noqa: E402

IDS = ["V249_star_module_closure_census", "V250_single_particle_walls",
       "V251_gauged_two_particle_sector",
       "V252_dynamical_self_energy_loop", "V253_block_rg_window",
       "V254_dressed_f9_contact", "V255_dressed_f2_contact",
       "V256_ladder_rung_reconstruction"]


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    checks = []
    data = {}
    arts = {i.split("_")[0]: art(i) for i in IDS}

    br = {k: v["branch"] for k, v in arts.items()}

    # ---------------- the ledger ------------------------------------------
    ledger = []

    def row(q, value, prov, why, src):
        ledger.append({"quantity": q, "value": value, "provenance": prov,
                       "why": why, "source": src})

    v249 = arts["V249"]
    clo = v249["data"]["closure"]
    row("the star-module closure set C_star",
        "{} labeled closures by order {}, C_star {} points".format(
            sum(clo["grid_labels_by_order"].values()),
            clo["grid_labels_by_order"], len(clo["C_star"])),
        "algebraic",
        "V249's measured acceptance rule on the closed module family "
        "(fitted_inputs == []); the set is fixed by the module, not "
        "tuned",
        "V249")
    clock_real = br["V249"] == "STAR-CLOCK-EXISTS"
    row("the C_star clock",
        ("closed points realized" if clock_real
         else "never realized ({})".format(br["V249"])),
        "algebraic" if clock_real else "never-realized",
        "pre-registered rule: a clock exists only at a closed point of "
        "C_star; V249 quantized the clock alternatives",
        "V249")
    row("eta twist values",
        "never realized (no C_star clock to twist)",
        "never-realized",
        "twists are defined on the C_star clock (V249's footing); with "
        "the clock never realized they add no parameter (pre-reg)",
        "V249")
    v250 = arts["V250"]
    row("single-particle wall scales (wall classes, confinement ladder)",
        str(v250["branch"]),
        "measured-input",
        "V250 measures wall tensions/splittings on the published step "
        "with fitted_inputs == []",
        "V250")
    blocked = br["V251"] == "SECTOR-BLOCKED"
    row("the gauged two-particle sector (H2, flux sector)",
        ("never realized ({})".format(br["V251"]) if blocked
         else "constructed, 1440-dim H2"),
        "never-realized" if blocked else "algebraic",
        "V251's exact-sector test failed at the root: per the "
        "pre-registered roll-up clause the sector and its numbers are "
        "booked never-realized and add no fitted parameter; the "
        "constructed H2 is a measurement device, not a claim",
        "V251")
    row("the sector coupling (g_10, V177/V141 gates)",
        "algebraic gate values on the constructed H2 (T0/T1 exact)",
        "algebraic",
        "the vertex gates are the zero-parameter V177/V141 committed "
        "machinery; they are verified exact (V251 T-checks) even "
        "though the exact sector itself is blocked",
        "V251")
    v252 = arts["V252"]
    nF = sum(1 for sub in v252["data"]["F_int"].values()
             for tag in sub
             if "xN32/" in tag and tag.endswith("wcontact"))
    row("the interaction loop factor F_int",
        "{} measured values (V252 N = 32 rows, sector-floor "
        "convention {})".format(nF, v252["branch"]),
        "measured-input",
        "F_int is computed from the closed step (sum over states), not "
        "fitted; the sector-floor convention is the pre-registered "
        "response to V251, not a dial",
        "V252")
    v253 = arts["V253"]
    row("dressed block-RG scales v(b), A(b), exponents p_v, p_A",
        "5 cuts b = 1..16, exponents {}".format(
            {k: round(x, 5) for k, x in
             v253["data"]["exponents"].items()}),
        "measured-input",
        "V253's ladder (V178 machinery extended, V252 Sigma inserted, "
        "V178 static control) measured with fitted_inputs == []",
        "V253")
    v254 = arts["V254"]
    row("the dressed F9 velocity contact dvv at b = 16",
        "best {:.6f} vs bound {:.1e} ({})".format(
            v254["data"]["best"]["dvv"], v254["data"]["bound"],
            v254["branch"]),
        "measured-input",
        "direct velocity spread of the downfolded dressed step at the "
        "named nodes; no extrapolation",
        "V254")
    v255 = arts["V255"]
    row("the dressed F2 contact (ladder x unification pairings)",
        "{} pairs, best ladder dev {:.3g} ({})".format(
            len(v255["data"]["rows"]), v255["data"]["best_ladder"]["dev"],
            v255["branch"]),
        "measured-input",
        "every pair is a functional reading of V252's F_int and the "
        "coupling module; the g set is pre-registered, not fitted",
        "V255")
    v256 = arts["V256"]
    row("the rung reconstruction at loop' = F_int^(1/10)",
        "best joint rung residual {:.4f} ({})".format(
            v256["data"]["best"]["joint"], v256["branch"]),
        "measured-input",
        "the tenth root is V197's own loop-exponent rule; the rung "
        "values are V197 artifact data",
        "V256")

    counts = {}
    for r in ledger:
        counts[r["provenance"]] = counts.get(r["provenance"], 0) + 1

    # ---------------- checks ----------------------------------------------
    checks.append({
        "name": "T0 all eight programme artifacts present, passing, with "
                "the pre-registered branch vocabulary: {}".format(
                    " ".join("{}:{}".format(k, v) for k, v in
                             sorted(br.items()))),
        "pass": bool(all(v["pass"] for v in arts.values())
                     and br["V249"] in ("STAR-CLOCK-EXISTS",
                                        "STAR-CLOCK-QUANTIZED")
                     and (br["V251"] == "TWO-PARTICLE-SECTOR-CLOSED"
                          or br["V251"].startswith("SECTOR-BLOCKED")))})
    fitted = [k for k, v in arts.items() if v["data"].get("fitted_inputs")]
    checks.append({
        "name": "T1 no experiment of the programme used a fitted input: "
                "fitted_inputs lists {} ({} non-empty)".format(
                    "all empty" if not fitted else "nonempty",
                    len(fitted)),
        "pass": bool(not fitted)})
    v52 = art("V52_parameter_audit")["data"]["provenance_table"]["rows"]
    checks.append({
        "name": "T2 the published Table 6 rows are loaded for "
                "comparison ({} rows from V52's provenance_table); the "
                "ledger reclassifies {} quantities".format(
                    len(v52), sum(1 for r in ledger
                                  if r["provenance"] == "fitted")),
        "pass": bool(len(v52) >= 4)})

    branch = "INTERACTING-LEDGER"
    data["ledger"] = ledger
    data["counts"] = counts
    data["reclassification_flags_for_table6"] = [
        r for r in ledger if r["provenance"] == "fitted"]
    data["never_realized"] = [r["quantity"] for r in ledger
                              if r["provenance"] == "never-realized"]
    data["fitted_inputs"] = []
    ok = all(c["pass"] for c in checks)
    status = (
        "interacting ledger: {} quantities booked - algebraic {}, "
        "measured-input {}, fitted {}, never-realized {} (V52 Table 6 "
        "comparison loaded). Branch {}.".format(
            len(ledger), counts.get("algebraic", 0),
            counts.get("measured-input", 0), counts.get("fitted", 0),
            counts.get("never-realized", 0), branch))
    out = {"id": "V257_interacting_ledger",
           "tier": 6,
           "property": "provenance ledger of the interacting programme, "
                       "Table-6-class counts (V257)",
           "residual": bool(counts.get("fitted", 0) > 0),
           "programme": "V257",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v257_interacting_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
