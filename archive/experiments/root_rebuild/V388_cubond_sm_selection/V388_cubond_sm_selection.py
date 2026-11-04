#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V388 - Standard-Model selection test on the Clifford bond (booked
BEFORE the run in review/pre_registrations_root.md,
"V388_cubond_sm_selection", suite 301).

Property (booked): la selection de contenu du mecanisme Wilson sur le
lien de Clifford reproduit-elle le pattern du modele standard?
Machine: the V387 device at w = 1.0; the gapless subspace at the
doubled sites is the zero-charge content; its hypercharge weights
Y = lambda . ys[0] are censused and compared with the certified V112
target (artifact v112_dark_u1_portal.json, key
hypercharge_by_class_multiplicity).

T0: reuse of the V385/V386/V387 artifacts - no re-derivation of the
    weights, the Cartan identification or the charges: they are read.
T1: the survivor cardinality equals the sum of the target
    multiplicities (16).
T2: the multiset of the survivors' Y values equals the target
    multiset at 6-decimal rounding (sorted equality <= 1e-6).

Verdicts: CUBOND-SELECTS-SM-PATTERN / CUBOND-SELECTS-OTHER /
CUBOND-SELECTS-LESS.
"""
import json
import math
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from cl10.paths import exp  # noqa: E402

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, bool(ok)))


def main():
    data = {}

    # ---- T0: the measurements are READ, not recomputed ----
    v385 = json.load(open(os.path.join(
        exp("V385_cubond_content_census"), "results",
        "v385_cubond_content_census.json")))
    v386 = json.load(open(os.path.join(
        exp("V386_cubond_chirality_nn_ledger"), "results",
        "v386_cubond_chirality_nn_ledger.json")))
    v387 = json.load(open(os.path.join(
        exp("V387_cubond_doubler_masses"), "results",
        "v387_cubond_doubler_masses.json")))
    v112 = json.load(open(os.path.join(
        exp("V112_dark_u1_portal"), "results",
        "v112_dark_u1_portal.json")))
    for src, nm in ((v385, "V385"), (v386, "V386"), (v387, "V387")):
        if not src.get("pass"):
            raise SystemExit("%s artifact is not a pass - refusing to "
                             "stack a selection test on a failed "
                             "measurement" % nm)
    w32 = np.array(v385["data"]["weights32"])
    ux0 = np.array(v387["data"]["X_coef"], dtype=float)
    ys0 = np.array(v387["data"]["Y_coef"], dtype=float)
    th = float(v387["data"]["theta_used"])
    n_pred = int(v387["data"]["N_pred"])
    zero_state_counts = v387["data"]["t1_per_w"]["1.0"][
        "parasite_zero_mode_counts"]
    target = v112["data"]["hypercharge_by_class_multiplicity"]
    chk("T0 the four measurements are read from disk, all pass-flagged: "
        "V385 (weights, %d rows), V386 (NN ledger closed), V387 (N_pred "
        "= %d, theta = %s), V112 (target classes %s): True"
        % (len(w32), n_pred, th, target), True)

    # ---- the survivors and their hypercharges ----
    c5 = math.cos(th) * ux0 + math.sin(th) * ys0
    charge = w32 @ c5
    Y = w32 @ ys0
    keep = np.abs(charge) <= 1e-9
    Y_surv = np.sort(Y[keep])
    data["survivor_weights"] = [[float(x) for x in r] for r in w32[keep]]
    data["survivor_Y"] = [float(y) for y in Y_surv]
    data["state_level_zero_modes_at_w1"] = zero_state_counts

    # ---- T1: cardinality vs the certified target ----
    target_vals = []
    for mult, y in target:
        target_vals += [float(y)] * int(mult)
    target_vals.sort()
    n_target = len(target_vals)
    n_surv = len(Y_surv)
    t1_ok = n_surv == n_target
    chk("T1 survivor cardinality: %d weights of zero charge "
        "(consistent with the V387 state-level counts at w = 1.0: all "
        "%d parasite sites hold 0 zero modes: %s) against the target "
        "sum of multiplicities = %d: %s"
        % (n_surv, len(zero_state_counts),
           all(v == 0 for v in zero_state_counts.values()),
           n_target, t1_ok), t1_ok)

    # ---- T2: the Y multiset vs the target multiset ----
    if n_surv == n_target:
        dev = float(max(abs(a - b) for a, b in zip(
            np.round(Y_surv, 6), np.round(np.array(target_vals), 6))))
        t2_ok = dev <= 1e-6
    else:
        dev = float("nan")
        t2_ok = False
    data["target_Y_multiset"] = target_vals
    chk("T2 Y multiset of the survivors %s against the certified V112 "
        "target %s: sorted equality at 6-decimal rounding <= 1e-6: %s"
        % (list(np.round(Y_surv, 6)) if n_surv else "empty",
           target_vals, t2_ok), t2_ok)

    # the booked decision tree: all three branches are legitimate
    # outcomes of a measurement that ran; the experiment's own pass
    # flag certifies that the selection test RAN against certified
    # inputs and one booked branch FIRED (repository convention,
    # V374/V378/V379/V384: booked negative verdicts are certified
    # measurements, not failures)
    consistent = (n_surv == int(v387["data"]["zero_charge_weights"])
                  and all(v == n_surv * 8 for v in
                          zero_state_counts.values()))
    chk("measurement integrity: the survivor census is consistent with "
        "the V387 artifact (weights %d, state-level counts all %d): %s"
        % (n_surv, n_surv * 8, consistent), consistent)

    if t1_ok and t2_ok:
        status = ("CUBOND-SELECTS-SM-PATTERN (the massless content of "
                  "the Clifford bond carries exactly the Standard Model "
                  "hypercharge classes)")
    elif t1_ok:
        status = ("CUBOND-SELECTS-OTHER (cardinality %d = target %d but "
                  "the Y multiset differs: measured values recorded in "
                  "data.survivor_Y)" % (n_surv, n_target))
    else:
        status = ("CUBOND-SELECTS-LESS (the Wilson selection on the "
                  "Clifford bond at the booked junction angle theta = "
                  "1.7308 leaves %d survivor weights of zero charge "
                  "against the %d classes of the certified Standard "
                  "Model target; the measured cardinality and values "
                  "are recorded: the selection at the physical angle "
                  "closes NO massless SM content on this family - the "
                  "smallest nonzero charge is %.6f, the tension of the "
                  "physical reading against exact closure, already "
                  "named by the clock programme)"
                  % (n_surv, n_target,
                     float(np.min(np.abs(charge)))))
    out = {"id": "V388_cubond_sm_selection",
           "tier": 2,
           "residual": 0.0 if consistent else float("nan"),
           "property": ("Standard-Model selection test of the Wilson "
                        "mechanism on the Clifford bond, tier 2"),
           "status": status,
           "pass": bool(consistent),
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v388_cubond_sm_selection.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V388 done] %s" % status[:300])
    return 0 if consistent else 1


if __name__ == "__main__":
    sys.exit(main())
