#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the verification runner.

Runs every experiment under ``experiments/`` (in name order), reads each
JSON verdict, and writes a single summary to ``results/validation_summary.json``.

Usage::

    python3 run_validation.py          # run all experiments
    python3 run_validation.py V0 V1    # run a subset

Each experiment is a **0-parameter** check: it confirms a prediction of the
publication (Section 10) or reports an honest residual.  The verdict JSON has
``pass``, ``status``, the 0-parameter number, and ``wall_time``.
"""

import glob
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
EXPS = os.path.join(ROOT, "experiments")


def find_experiments():
    dirs = sorted(os.path.dirname(p) for p in glob.glob(os.path.join(EXPS, "*", "*.py")))
    return [d for d in dirs if d != EXPS]


def run_one(exp_dir):
    script = os.path.join(exp_dir, os.path.basename(exp_dir) + ".py")
    t0 = time.time()
    proc = subprocess.run([sys.executable, script], capture_output=True,
                          text=True, cwd=ROOT)
    time_ = time.time() - t0
    residir = os.path.join(exp_dir, "results")
    verdict = None
    if os.path.isdir(residir):
        for jf in glob.glob(os.path.join(residir, "*.json")):
            try:
                verdict = json.load(open(jf))
            except Exception:  # noqa: BLE001
                pass
    return {
        "id": os.path.basename(exp_dir),
        "dir": os.path.basename(exp_dir),
        "tier": verdict.get("tier") if verdict else None,
        "property": verdict.get("property") if verdict else None,
        "residual": bool(verdict.get("residual")) if verdict else False,
        "pass": bool(verdict["pass"]) if verdict else False,
        "status": verdict.get("status", "no_verdict") if verdict else "no_verdict",
        "returncode": proc.returncode,
        "wall_time": round(time_, 3),
        "verdict": verdict,
    }


def main():
    only = set(sys.argv[1:])
    exps = [d for d in find_experiments() if not only or os.path.basename(d) in only]
    print("=" * 68)
    print("CL10  —  0-parameter verification program  (the GPU engine)")
    print("=" * 68)
    import cl10.engine as E
    print("  backend:", E._BACKEND)
    results = []
    for d in exps:
        r = run_one(d)
        results.append(r)
        mark = "PASS" if r["pass"] else "FAIL"
        tier = f"T{r['tier']}" if r["tier"] is not None else "   "
        print(f"  [{mark}]  {r['id']:18s}  {tier}  {r['status']}  ({r['wall_time']} s)")
    n_pass = sum(1 for r in results if r["pass"])
    n = len(results)
    # the TOE verdict: the hard tiers (0-3) all PASS (V0-V7 + V11-V13).
    # The cosmology (tier 4) and the falsifiability (tier 5) are supporting.
    hard = [r for r in results if r["tier"] in (0, 1, 2, 3)]
    toe = bool(hard) and all(r["pass"] for r in hard)
    named = [r["id"] for r in results if r["residual"]]
    summary = {
        "n": n,
        "n_pass": n_pass,
        "pass": bool(n_pass == n),
        "toe": toe,
        "toe_criterion": "all of tiers 0-3 (V0-V7, V11-V13) PASS",
        "named_residuals": named,
        "experiments": results,
    }
    outdir = os.path.join(ROOT, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print(f"  ->  {n_pass}/{n}  PASS  ->  results/validation_summary.json")
    vstr = ("TOE (a consistent 0-parameter unification of matter, gauge, "
            "gravity)" if toe else "NOT a TOE (a hard tier failed)")
    print(f"  TOE VERDICT:  {vstr}")
    if named:
        print(f"  NAMED RESIDUALS:  {', '.join(named)}")
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main())
