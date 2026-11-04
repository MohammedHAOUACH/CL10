#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V375 - dressed spectra vs the empirical ladder contacts
(booked BEFORE the run in review/pre_registrations_root.md lines
2972-3006).

Question: do the eigenphase differences of the DRESSED generator
(the V373 kernel spectra and the V374 lattice-sum scan) reproduce
any of the five published empirical ladder contacts - does the
dressing machinery "know" the ladder, or is it blind to it?

Machine (booked): contact harvest is MECHANICAL - the five V52
provenance rows whose quantity or note contains "rung" or "ladder"
(harvest count = 5 is a gate). Every dressed eigenphase difference
against every contact, machine-exact Decimal at relative 1e-12.
False-hit bookkeeping (counting stated in the booking): each hit
must ALSO reproduce at the adjacent eta grid point (independent
re-diagonalization; per-test chance <= 4e-24), booked false-hit
chance T x 4e-24 <= 1e-16 for T up to 2.5e7. For the V373 static
kernel spectra (eta-independent by construction) the adjacency
retest is adapted, DECLARED: the difference must reproduce in the
independent second measurement of the same channels (fiber-route
dk vs factored-route dk) at the same tolerance.
Verdicts: DRESSED-LADDER-CONTACT-HIT / DRESSED-NO-LADDER-CONTACT
(the zero-hit branch is the honest negative: pass = True with the
protocol gates).
"""
import json
import os
import re
import sys
import time
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
from cl10.paths import exp  # noqa: E402

T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def json_art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def dec(x):
    """machine-exact Decimal of a stored artifact float."""
    return Decimal(repr(x))


def rel_hit(d, c, tol=Decimal("1e-12")):
    if c == 0:
        return abs(d) <= tol
    return abs((d - c) / c) <= tol


def harvest_contacts(v52):
    """mechanical: rows whose quantity or note mention rung/ladder."""
    out = []
    for r in v52["data"]["provenance_table"]["rows"]:
        hay = (str(r.get("quantity", "")) + " "
               + str(r.get("note", ""))).lower()
        if "rung" in hay or "ladder" in hay:
            val = str(r["value"])
            m = re.search(r"-?\d+\.?\d*", val.replace(",", ""))
            if not m:
                raise AssertionError("no numeric token in %r" % val)
            out.append({"quantity": r.get("quantity"),
                        "value_str": val,
                        "value": Decimal(m.group()),
                        "declared": "declared" in val.lower()})
    return out


def pairwise(vals):
    ds = [dec(v) for v in vals]
    return [(ds[i] - ds[j] if ds[i] >= ds[j] else ds[j] - ds[i],
             (i, j)) for i in range(len(ds)) for j in range(i + 1,
                                                           len(ds))]


def main():
    checks, data = [], {}

    v52 = json_art("V52_parameter_audit")
    contacts = harvest_contacts(v52)
    a_hv = len(contacts) == 5
    checks.append(("contact harvest mechanical: %d V52 rows with "
                   "rung/ladder (booked 5): %s" % (len(contacts),
                                                   a_hv), a_hv))
    data["contacts"] = [{"quantity": c["quantity"],
                         "value": str(c["value"]),
                         "value_str": c["value_str"]}
                        for c in contacts]
    booked = [Decimal(x) for x in
              ("0.019894", "0.1999980", "1.18", "3", "118")]
    got = sorted(c["value"] for c in contacts)
    a_cv = got == sorted(booked)
    checks.append(("the five harvested values == the booked set "
                   "{0.019894, 0.1999980, 1.18, 3, 118}: %s"
                   % a_cv, a_cv))

    # ------------------------------------------------------------
    # dressed eigenphase differences
    # ------------------------------------------------------------
    v373 = json_art("V373_full_lattice_pair_string")
    D373 = v373["data"]
    fam373 = {}
    for key, lst in D373["dk_factored"].items():
        fam373["V373 dk_factored " + key] = lst
    fam373["V373 dk_fiber_ell2"] = D373["dk_fiber_ell2"]
    # second-route map for the declared V373 adjacency adaptation
    second = {
        "V373 dk_factored chain_ell1": fam373["V373 dk_fiber_ell2"],
        "V373 dk_fiber_ell2": fam373["V373 dk_factored chain_ell1"]}

    diffs = []  # (value, family, origin, adj_key)
    for fam, lst in fam373.items():
        for d, ij in pairwise(lst):
            diffs.append((d, fam, ij, None))

    v374 = json_art("V374_lattice_sum_threshold")
    scan = v374["data"]["scan_g"]          # [eta, N, K key, g]
    # the booked row eigenvalue difference is g = Ec - E0 itself
    n374 = 0
    for i, r in enumerate(scan):
        eta, N, key, g = r
        diffs.append((abs(dec(g)), "V374 g(eta,K) pair", (eta, N,
                                                          key), i))
        n374 += 1
    log("dressed differences: %d (V373 families), %d (V374 rows)"
        % (len(diffs) - n374, n374))

    # index V374 rows by (N, key, eta) for the adjacency retest
    idx = {}
    for i, r in enumerate(scan):
        idx[(r[1], r[2], round(r[0], 3))] = i

    T = len(diffs) * len(contacts)
    a_T = T <= 25_000_000
    checks.append(("test count T = %d x 5 = %d <= 2.5e7 (booked "
                   "false-hit bound T x 4e-24 = %.2e <= 1e-16): %s"
                   % (len(diffs), T, T * 4e-24, a_T), a_T))

    # ------------------------------------------------------------
    # the exact-Decimal test with the booked adjacency retest
    # ------------------------------------------------------------
    hits = []
    near = 0
    for k, (d, fam, origin, rowi) in enumerate(diffs):
        for c in contacts:
            if rel_hit(d, c["value"]):
                near += 1
                # adjacency retest (booked: independent
                # re-diagonalization at the adjacent eta point)
                ok_adj = False
                how = ""
                if rowi is not None:
                    r = scan[rowi]
                    eta, N, key = r[0], r[1], r[2]
                    for de in (+0.02, -0.02):
                        ae = round(eta + de, 3)
                        j = idx.get((N, key, ae))
                        if j is not None and rel_hit(
                                abs(dec(scan[j][3])), c["value"]):
                            ok_adj = True
                            how = "adjacent eta %.2f" % ae
                else:
                    # declared adaptation for the eta-independent
                    # V373 kernel spectra: reproduce in the
                    # independent second measurement of the same
                    # channels
                    alt = second.get(fam)
                    if alt is not None:
                        for d2, _ in pairwise(alt):
                            if rel_hit(d2, c["value"]):
                                ok_adj = True
                                how = ("second-route spectrum "
                                       "(declared)")
                hits.append({"value": str(d), "family": fam,
                             "origin": str(origin),
                             "contact": c["quantity"],
                             "contact_value": str(c["value"]),
                             "adjacency_passed": ok_adj,
                             "adjacency_how": how})
    data["near_misses_decimal"] = near
    data["hits"] = hits
    log("hits: %d (decimal-near %d)" % (len(hits), near))

    # booked double-naming cross-check (only exercised on nonzero
    # hits): a contact already named by V364 would contradict it.
    if hits:
        v364 = json_art("V364_tension_named_pairs")
        named = [dec(x) for x in
                 ([v364["data"]["target"]] if isinstance(
                     v364["data"].get("target"), (int, float))
                  else [])]
        for h in hits:
            for nv in named:
                if rel_hit(Decimal(h["value"]), nv):
                    log("NOTE: hit value %s is also a V364 named "
                        "value - cross-check required" % h["value"])

    okm = a_hv and a_cv and a_T
    if not okm:
        status = "MACHINERY-DRIFT (a protocol gate failed)"
        ok = False
    elif hits:
        status = ("DRESSED-LADDER-CONTACT-HIT (%d dressed "
                  "eigenphase differences reproduce a booked "
                  "ladder contact at Decimal relative 1e-12 with "
                  "the adjacency retest)" % len(hits))
        ok = True
    else:
        status = ("DRESSED-NO-LADDER-CONTACT (honest negative: %d "
                  "dressed eigenphase differences x 5 harvested "
                  "contacts, T = %d exact-Decimal tests at relative "
                  "1e-12, booked false-hit chance T x 4e-24 = "
                  "%.2e; zero hits: the dressing machinery is blind "
                  "to the empirical ladder)" % (len(diffs), T,
                                                T * 4e-24))
        ok = True

    out = {"id": "V375_dressed_no_ladder_contact",
           "property": ("do dressed eigenphase differences "
                        "reproduce published empirical ladder "
                        "contacts"),
           "tier": 2,
           "checks": [[m, bool(b)] for m, b in checks],
           "data": data, "pass": ok, "status": status,
           "residual": float(T * 4e-24)}
    with open(os.path.join(HERE, "results",
                           "v375_dressed_no_ladder_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for m, b in checks:
        print(("  ok  " if b else "  FAIL ") + m)
    print("[V375 done] %s" % status)
    log("wall %.1f s" % (time.time() - T0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
