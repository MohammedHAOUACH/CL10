#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the verification runner.

Runs every experiment under ``experiments/`` (in name order), reads each
JSON verdict, and writes a single summary to ``results/validation_summary.json``.

Usage::

    python3 run_validation.py              # run all experiments, in parallel
    python3 run_validation.py V0 V1        # run a subset
    python3 run_validation.py --jobs 4     # cap the worker pool (default: cores)
Each experiment is a **0-parameter** check: it confirms a prediction of the
publication (Section 10) or reports an honest residual.  The verdict JSON has
``pass``, ``status``, the 0-parameter number, and ``wall_time``.
"""

import argparse
import glob
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.abspath(__file__))
EXPS = os.path.join(ROOT, "experiments")


def find_experiments():
    dirs = sorted(os.path.dirname(p) for p in glob.glob(os.path.join(EXPS, "*", "*.py")))
    return [d for d in dirs if d != EXPS]


def read_verdict(exp_dir):
    residir = os.path.join(exp_dir, "results")
    if not os.path.isdir(residir):
        return None
    # un dossier peut contenir plusieurs JSON (verdict + meta de figures).
    # glob() n'est PAS trie : on choisit deterministement le JSON qui
    # EST un verdict (cle "pass" presente), en priorite celui nomme
    # d'apres l'experience ; les meta (fig*.meta.json) sont ignorees.
    jfs = sorted(glob.glob(os.path.join(residir, "*.json")))
    preferred = [j for j in jfs
                 if os.path.basename(j) == os.path.basename(exp_dir) + ".json"]
    for jf in preferred + [j for j in jfs if j not in preferred]:
        try:
            d = json.load(open(jf))
        except Exception:  # noqa: BLE001
            continue
        if isinstance(d, dict) and "pass" in d:
            return d
    return None


def report_one(exp_dir):
    # --report-only: le verdict deja ecrit sur disque, sans relancer le
    # script. Le resume global est un rollup des JSON: un run partiel ne
    # doit pas obliger a relancer les 138 experiences pour le reecrire.
    verdict = read_verdict(exp_dir)
    return {
        "id": os.path.basename(exp_dir),
        "dir": os.path.basename(exp_dir),
        "tier": verdict.get("tier") if verdict else None,
        "property": verdict.get("property") if verdict else None,
        "residual": bool(verdict.get("residual")) if verdict else False,
        "pass": bool(verdict["pass"]) if verdict else False,
        "status": verdict.get("status", "no_verdict") if verdict else "no_verdict",
        "returncode": None,
        "crashed": False,
        "wall_time": verdict.get("wall_time") if verdict else None,
        "verdict": verdict,
    }


def run_one(exp_dir):
    script = os.path.join(exp_dir, os.path.basename(exp_dir) + ".py")
    # one BLAS thread per subprocess: with --jobs 12 on 12 cores, the
    # default nested threading (12 jobs x 12 BLAS threads) thrashes and
    # made V110 8x slower inside the suite than standalone.
    env = os.environ.copy()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[v] = "1"
    t0 = time.time()
    proc = subprocess.run([sys.executable, script], capture_output=True,
                          text=True, cwd=ROOT, env=env)
    time_ = time.time() - t0
    verdict = read_verdict(exp_dir)
    # un crash du script invalide tout verdict JSON deja present sur disque :
    # sans ce garde-fou, une experience qui plante herite du PASS de sa
    # precedente execution (un resultat perime compte comme valide).
    crashed = proc.returncode != 0
    if crashed:
        tail = (proc.stderr or "").strip().splitlines()
        verdict = None
        status = f"CRASH (rc={proc.returncode}): {tail[-1] if tail else 'no stderr'}"
        # the one-line status cannot localise a backend-dependent crash: the
        # full traceback is kept next to the summary, per experiment.
        crashdir = os.path.join(ROOT, "results", "crashes")
        os.makedirs(crashdir, exist_ok=True)
        with open(os.path.join(crashdir, os.path.basename(exp_dir) + ".log"),
                  "w") as cf:
            cf.write(f"rc={proc.returncode}  backend={env.get('CL10_BACKEND', 'auto')}\n")
            cf.write(f"script={script}\n\n--- stderr ---\n{proc.stderr or ''}\n")
            cf.write(f"--- stdout ---\n{proc.stdout or ''}\n")
    else:
        status = verdict.get("status", "no_verdict") if verdict else "no_verdict"
    return {
        "id": os.path.basename(exp_dir),
        "dir": os.path.basename(exp_dir),
        "tier": None if crashed else (verdict.get("tier") if verdict else None),
        "property": None if crashed else (verdict.get("property") if verdict else None),
        "residual": False if crashed else (bool(verdict.get("residual")) if verdict else False),
        "pass": False if crashed else (bool(verdict["pass"]) if verdict else False),
        "status": status,
        "returncode": proc.returncode,
        "crashed": crashed,
        "wall_time": round(time_, 3),
        "verdict": verdict,
    }


def main():
    ap = argparse.ArgumentParser(description="the CL10 validation runner")
    ap.add_argument("--jobs", type=int, default=min(12, os.cpu_count() or 4),
                    help="parallel experiment workers (subprocesses)")
    ap.add_argument("only", nargs="*", help="experiment dir names to run")
    ap.add_argument("--report-only", action="store_true",
                    help="roll up the verdict JSONs already on disk; run nothing")
    args = ap.parse_args()
    only = set(args.only)
    exps = [d for d in find_experiments() if not only or os.path.basename(d) in only]
    print("=" * 68)
    print("CL10  —  0-parameter verification program  (the GPU engine)")
    print("=" * 68)
    import cl10.engine as E
    print(f"  backend: {E._BACKEND}   jobs: {args.jobs}")
    # experiments are independent subprocesses writing to their own
    # results/ dir: run them concurrently, REPORT EACH AS IT FINISHES so a
    # slow experiment is visible live; everything is also logged to
    # results/run_log.txt.
    t_start = time.time()
    logp = os.path.join(ROOT, "results", "run_log.txt")
    with open(logp, "w") as logf:
        def emit(line):
            print(line, flush=True)
            logf.write(line + "\n")
            logf.flush()
        done = []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pick = report_one if args.report_only else run_one
            futs = [pool.submit(pick, d) for d in exps]
            for fut in as_completed(futs):
                r = fut.result()
                done.append(r)
                mark = "PASS" if r["pass"] else "FAIL"
                tier = f"T{r['tier']}" if r["tier"] is not None else "   "
                el = time.time() - t_start
                emit(f" [{el:7.1f}s] [{mark}]  {r['id']:34s} {tier} "
                     f"{r['status']}  ({r['wall_time']} s)")
        results = sorted(done, key=lambda r: r["id"])
    n_pass = sum(1 for r in results if r["pass"])
    n = len(results)
    # the U10 verdict: the hard tiers (0-3) all PASS (V0-V7 + V11-V13).
    # The cosmology (tier 4) and the falsifiability (tier 5) are supporting.
    hard = [r for r in results if r["tier"] in (0, 1, 2, 3)]
    u10 = bool(hard) and all(r["pass"] for r in hard)
    named = [r["id"] for r in results if r["residual"]]
    # the named walls/residuals are REDUCED: each has a 0-parameter MECHANISM
    # that passes.  W2 is now completed on the 10-loop ladder (V57); W4 by the
    # UV-finite discrete RG (V58).
    reduced = {
        "V9_dark_energy": ("V17_w1_screening",
                           "the Z2 512/512 cancellation -> Lambda = 0"),
        "V15_seesaw": ("V57_neutrino",
                       "the 2-step seesaw on the 10-loop ladder (0-param M_N, m_D)"),
        "V16_families": ("V19_w3_families3",
                         "the 3 families = the 3 space"),
        "V47_2loop_uvfinite": ("V58_uvr",
                               "the UV-finite discrete RG (bounded tower U^(2^p))"),
    }
    reduced_ok = {w: any(r["id"] == mid and r["pass"] for r in results)
                  for w, (mid, _desc) in reduced.items()}
    summary = {
        "n": n,
        "n_pass": n_pass,
        "pass": bool(n_pass == n),
        "u10": u10,
        "u10_criterion": "all of tiers 0-3 (V0-V7, V11-V13) PASS",
        "named_residuals": named,
        "reduced_walls": [{"wall": w, "mechanism": f"{mid} ({desc})",
                            "mechanism_pass": reduced_ok[w]}
                           for w, (mid, desc) in reduced.items()],
        "experiments": results,
    }
    outdir = os.path.join(ROOT, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print(f"  ->  {n_pass}/{n}  PASS  ->  results/validation_summary.json")
    vstr = ("U10 (a 0-parameter unification of matter, gauge, and gravity)"
            if u10 else "NOT a 0-parameter unification (a hard tier failed)")
    print(f"  U10 VERDICT:  {vstr}")
    if named:
        print(f"  NAMED RESIDUALS (the exact values):  {', '.join(named)}")
    if reduced_ok and all(reduced_ok.values()):
        print(f"  THE {len(reduced)} NAMED RESIDUALS ARE REDUCED: each has a 0-parameter mechanism (PASS)")
        for w, (mid, desc) in reduced.items():
            print(f"    - {w:18s} reduced by {mid} ({desc})")
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main())
