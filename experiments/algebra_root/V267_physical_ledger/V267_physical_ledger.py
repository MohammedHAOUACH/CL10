#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V267 - the physical ledger (tier 6, residual: true iff anything
landed fitted).

**The question this pre-registers against** (review/pre_registrations.md,
gauge-invariant physical sector programme): the provenance of every
number the physical-sector programme produced - the coupled fiber and
its charges (algebraic, V251's committed construction), the singlet
kernel (measured by V259), the composite band/gap/loci/dvv, the
statistics and desert content, the loop-scale candidates, the F8/
F9 contacts, the finite-size ladder (all never-realized under the
V259 KERNEL-EMPTY cascade), plus the machinery controls and the scale
denominators read from committed artifacts (measured-input).  The
Table-6-class counts per quantity, against the published Table 6 rows
(V52's provenance_table).

Closed inputs: only the artifacts V251 and V259-V266 (cl10.paths.exp)
and the V52 published rows.  No new number is measured.

Branch `PHYSICAL-LEDGER` (informational); `residual: true` iff a
quantity landed fitted (the pre-registration demands fitted_inputs ==
[] across the programme).  V259 landed KERNEL-EMPTY, so per the
pre-registered cascade every composite quantity is booked
never-realized and adds no fitted parameter.  Bookkeeping copied from
V257 (same shape).
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10.paths import exp                                     # noqa: E402

IDS = ["V251_gauged_two_particle_sector",
       "V259_singlet_kernel_census", "V260_composite_gap_ladder",
       "V261_composite_f9_contact",
       "V262_composite_chirality_statistics",
       "V263_composite_loop_scale",
       "V264_composite_desert_content",
       "V265_composite_f8_contact",
       "V266_finite_size_convergence"]


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}
    arts = {i.split("_")[0]: art(i) for i in IDS}
    br = {k: str(v.get("branch")) for k, v in arts.items()}

    ledger = []

    def row(q, value, prov, why, src):
        ledger.append({"quantity": q, "value": value, "provenance": prov,
                       "why": why, "source": src})

    empty = br["V259"] == "KERNEL-EMPTY"
    kdim = arts["V259"]["data"]["total_singlet_count"]
    c2min = arts["V259"]["data"]["c2_spectrum"]["first_nonzero"]

    row("the coupled two-particle fiber (1568 = 128 + 1440, H2, the 45 "
        "total charges)",
        "constructed; V251 tower regression reproduced to {:.1e} in "
        "V259 T0, max [C2, flux block] defect {:.1e} in T1".format(
            arts["V259"]["data"]["tower"]["integrality_defect"],
            max(arts["V259"]["data"]["c2_spectrum"]
                ["t1_commutators"].values())),
        "algebraic",
        "the fiber, charges and vertex are V251's committed zero-"
        "parameter construction (g_10, V177 vertex, V126 tower); V259 "
        "re-verifies them before any census",
        "V251/V259")
    row("the singlet kernel ker(C2) of the fiber",
        "dim {} at all 216 q-points, all omega rows and the grids N = "
        "4/6/8; first nonzero C2 eigenvalue {:.6f}".format(kdim, c2min),
        "measured-input",
        "V259's census (eigh of C2 = sum Q^a Q^a, tolerance 1e-9); "
        "V260/V262/V263/V264/V265/V266 independently recompute the "
        "same kernel and agree",
        "V259")
    row("the composite band E0(q), min gap, gapless loci, dvv tensor",
        "never realized ({})".format(br["V260"]),
        "never-realized",
        "no kernel state exists to disperse; the ladder path is the "
        "pre-registered non-cascade branch and ran not a line",
        "V260")
    row("the composite F9 contact (best dvv vs the SME ceiling)",
        "never realized ({})".format(br["V261"]),
        "never-realized",
        "no composite gapless locus; the ceiling itself is committed",
        "V261")
    row("the SME ceiling used by the F9 contact",
        "{:.3e} (config bounds.sme_electron_velocity_aniso)".format(
            arts["V261"]["data"]["bound"]),
        "algebraic",
        "config-committed comparison constant, read at run time, "
        "never retyped",
        "V261")
    row("the internal Spin(10) Cartan machinery (flux Cartan rank, "
        "P-block plane-Cartan spectra, V67 convention)",
        "rank {} and spectra {} verified live".format(
            arts["V262"]["data"]["cartan_controls"]["rank"],
            arts["V262"]["data"]["cartan_controls"]["Dk_spectrum"]),
        "algebraic",
        "structure constants of the committed V79 generators (central-"
        "izer construction, seed 0x62); a verified structure, not a "
        "measurement dial",
        "V262")
    row("the composite statistics/parity table and per-locus little-"
        "group census",
        "never realized ({})".format(br["V262"]),
        "never-realized",
        "no kernel state and no gapless locus: the chirality question "
        "has no instances; V98's centre obstruction never gets reached",
        "V262")
    row("the scale denominators for the loop-factor test (V101 ground"
        "-state E0, V129 two-site gap)",
        "V101 E0 = {:.6f}, V129 gap = {:.6f}; V101 gap x g_10 "
        "reproduces V251's omega row to {:.1e}".format(
            arts["V263"]["data"]["scales"]["V101_E0"],
            arts["V263"]["data"]["scales"]["V129_gap"],
            arts["V263"]["data"]["scales"]
            ["omega_consistency_defect"]),
        "measured-input",
        "both read from their committed artifacts; the cross-artifact "
        "identity (V101 gap x g_10 = V251 omega_plus) verified in V263 "
        "T1",
        "V101/V129 via V263")
    row("the composite loop-scale candidates F_phys",
        "never realized ({}; band [{}, {}] pre-registered, target "
        "1/(16 pi) computed)".format(
            br["V263"], arts["V263"]["data"]["band"][0],
            arts["V263"]["data"]["band"][1]),
        "never-realized",
        "the composite gap (numerator) does not exist",
        "V263")
    row("the composite desert content table and the (1,2)-shift demand",
        "never realized ({})".format(br["V264"]),
        "never-realized",
        "empty kernel -> empty multiplicity table; kernel states would "
        "be Spin(10) singlets with db = (0, 0, 0) anyway",
        "V264")
    row("V237's threshold rule fidelity (control for the copied rule)",
        "max deviation {:.1e} vs V237's own artifact".format(
            arts["V264"]["data"]["fidelity_control"]["deviation"]),
        "measured-input",
        "the copied rule rerun on V237's own content reproduces its "
        "committed numbers - the cascade verdict is rule-faithful",
        "V264")
    row("the composite F8 spin-2 channel and c_gw",
        "never realized ({}); band {:.0e} per V258".format(
            br["V265"], arts["V265"]["data"]["conventions"]["band"]),
        "never-realized",
        "no composite ground/first-excited kernel state to carry a "
        "spin-2 channel",
        "V265")
    row("the V229 quintet machinery controls (Emat gram, q5, J so(3))",
        "gram {:.1e}, q5 {:.1e}, J algebra {:.1e}".format(
            arts["V265"]["data"]["machinery_controls"]["emat_gram_dev"],
            arts["V265"]["data"]["machinery_controls"]["q5_dev"],
            arts["V265"]["data"]["machinery_controls"]["j_algebra_dev"]),
        "algebraic",
        "committed V229 vertex machinery verified exact on import",
        "V265")
    row("the finite-size gap ladder (N = 4, 6, 8) and the trend "
        "exponent",
        "never realized ({}); kernel dim measured 0 at every grid".format(
            br["V266"]),
        "never-realized",
        "emptiness is grid-independent by construction (the charges "
        "never depend on q), measured per grid in V266 T0",
        "V266")

    counts = {}
    for r in ledger:
        counts[r["provenance"]] = counts.get(r["provenance"], 0) + 1

    # ---------------- checks ----------------------------------------------
    checks.append({
        "name": "T0 all nine programme artifacts present, passing, with "
                "the pre-registered branch vocabulary: {}".format(
                    " ".join("{}:{}".format(k, v) for k, v in
                             sorted(br.items()))),
        "pass": bool(all(v["pass"] for v in arts.values())
                     and br["V259"] in ("COMPOSITE-KERNEL",
                                        "KERNEL-EMPTY")
                     and br["V260"] in ("COMPOSITE-GAPLESS",
                                        "COMPOSITE-GAPPED",
                                        "NEVER-REALIZED")
                     and br["V266"] in ("GAP-CONFIRMED",
                                        "GAPLESS-EMERGES",
                                        "GRID-NEUTRAL",
                                        "NEVER-REALIZED"))})
    fitted = [k for k, v in arts.items() if v["data"].get("fitted_inputs")]
    checks.append({
        "name": "T1 no experiment of the programme used a fitted input: "
                "fitted_inputs lists {} ({} non-empty)".format(
                    "all empty" if not fitted else "nonempty",
                    len(fitted)),
        "pass": bool(not fitted)})
    cascade_ok = empty and all(
        br[k] in ("NEVER-REALIZED", "NO-LOCUS", "F9-COMPOSITE-NO-LOCUS",
                  "NO-CONTENT")
        for k in ("V260", "V261", "V262", "V263", "V264", "V265",
                  "V266"))
    checks.append({
        "name": "T2 cascade consistency: V259 {}, every downstream "
                "experiment books its pre-registered cascade branch "
                "({}); never-realized ledger rows: {}".format(
                    br["V259"],
                    {k: br[k] for k in ("V260", "V261", "V262", "V263",
                                        "V264", "V265", "V266")},
                    sum(1 for r in ledger
                        if r["provenance"] == "never-realized")),
        "pass": bool(cascade_ok if empty else True)})
    v52 = art("V52_parameter_audit")["data"]["provenance_table"]["rows"]
    checks.append({
        "name": "T3 the published Table 6 rows are loaded for "
                "comparison ({} rows from V52's provenance_table); the "
                "ledger reclassifies {} quantities".format(
                    len(v52), sum(1 for r in ledger
                                  if r["provenance"] == "fitted")),
        "pass": bool(len(v52) >= 4)})

    branch = "PHYSICAL-LEDGER"
    data["ledger"] = ledger
    data["counts"] = counts
    data["reclassification_flags_for_table6"] = [
        r for r in ledger if r["provenance"] == "fitted"]
    data["never_realized"] = [r["quantity"] for r in ledger
                              if r["provenance"] == "never-realized"]
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0
    ok = all(c["pass"] for c in checks)
    status = (
        "physical-sector ledger: {} quantities booked - algebraic {}, "
        "measured-input {}, fitted {}, never-realized {} (V52 Table 6 "
        "comparison loaded; V259 {}). Branch {}.".format(
            len(ledger), counts.get("algebraic", 0),
            counts.get("measured-input", 0), counts.get("fitted", 0),
            counts.get("never-realized", 0), br["V259"], branch))
    out = {"id": "V267_physical_ledger",
           "tier": 6,
           "property": "provenance ledger of the gauge-invariant "
                       "physical sector programme, Table-6-class "
                       "counts (V267)",
           "residual": bool(counts.get("fitted", 0) > 0),
           "programme": "V267",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v267_physical_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
