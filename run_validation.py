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
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm

ROOT = os.path.dirname(os.path.abspath(__file__))
EXPS = os.path.join(ROOT, "experiments")


def find_experiments():
    dirs = sorted(os.path.dirname(p) for p in glob.glob(os.path.join(EXPS, "*", "*", "*.py")))
    return [d for d in dirs if d != EXPS]


def verdict_file(exp_dir):
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
            return jf
    return None


def read_verdict(exp_dir):
    jf = verdict_file(exp_dir)
    if jf is None:
        return None
    try:
        return json.load(open(jf))
    except Exception:  # noqa: BLE001
        return None


def stale_verdict(exp_dir):
    """True when the experiment crashed AFTER its verdict JSON was written. A
    crash leaves the previous verdict on disk, and a subset run folds the
    experiments it did not run from that disk copy: without this rule a crashed
    experiment is counted as the PASS of some earlier run."""
    jf = verdict_file(exp_dir)
    if jf is None:
        return False
    cl = os.path.join(ROOT, "results", "crashes",
                      os.path.basename(exp_dir) + ".log")
    try:
        return os.path.getmtime(cl) > os.path.getmtime(jf)
    except OSError:
        return False


def reported_backend(verdict):
    """The backend an experiment recorded in its verdict JSON, if it records one."""
    if not verdict:
        return None
    b = verdict.get("backend")
    if b:
        return b
    data = verdict.get("data")
    if isinstance(data, dict):
        return data.get("backend")
    return None


def report_one(exp_dir):
    # --report-only: le verdict deja ecrit sur disque, sans relancer le
    # script. Le resume global est un rollup des JSON: un run partiel ne
    # doit pas obliger a relancer les 112 experiences pour le reecrire.
    verdict = read_verdict(exp_dir)
    stale = bool(verdict) and stale_verdict(exp_dir)
    return {
        "id": os.path.basename(exp_dir),
        "dir": os.path.basename(exp_dir),
        "tier": verdict.get("tier") if verdict else None,
        "property": verdict.get("property") if verdict else None,
        "residual": bool(verdict.get("residual")) if verdict else False,
        "pass": bool(verdict["pass"]) if verdict and not stale else False,
        "status": (f"{verdict.get('status', 'no_verdict')} [STALE: this "
                   f"experiment crashed after the verdict was written]"
                   if verdict and stale else
                   (verdict.get("status", "no_verdict") if verdict else "no_verdict")),
        "returncode": None,
        "crashed": False,
        "stale": stale,
        "wall_time": verdict.get("wall_time") if verdict else None,
        "backend": reported_backend(verdict),
        "verdict": verdict,
    }


def _drain(stream, chunks, live):
    """Collect one pipe of a child's output.  tqdm draws its bars as ``\\r``
    frames with no newline, so the raw chunks are scanned for them: this is how
    the runner can show what a 45-minute scan is doing instead of a silent
    terminal (the cone scans draw their own bond/sweep bars)."""
    while True:
        try:
            # read1, not read: read(4096) blocks until the pipe holds 4096
            # bytes - for a slow scan that means EOF, and every frame would
            # arrive only after the experiment is over.
            c = stream.read1(4096)
        except (OSError, ValueError):        # the pipe closed under us
            break
        if not c:
            break
        chunks.append(c)
        if live is not None:
            for seg in c.split(b"\r"):
                s = seg.strip()
                if b"%|" in s:
                    live["progress"] = s.decode("utf-8", "replace")[-46:]


def run_one(exp_dir, live=None):
    script = os.path.join(exp_dir, os.path.basename(exp_dir) + ".py")
    # one BLAS thread per subprocess: with --jobs 12 on 12 cores, the
    # default nested threading (12 jobs x 12 BLAS threads) thrashes and
    # made V110 8x slower inside the suite than standalone.
    env = os.environ.copy()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[v] = "1"
    t0 = time.time()
    proc = subprocess.Popen([sys.executable, script], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, cwd=ROOT, env=env)
    out, err = [], []
    t_out = threading.Thread(target=_drain, args=(proc.stdout, out, None))
    t_err = threading.Thread(target=_drain, args=(proc.stderr, err, live))
    t_out.start()
    t_err.start()
    rc = proc.wait()
    t_out.join()
    t_err.join()
    stdout_s = b"".join(out).decode("utf-8", "replace")
    stderr_s = b"".join(err).decode("utf-8", "replace")
    time_ = time.time() - t0
    verdict = read_verdict(exp_dir)
    # un crash du script invalide tout verdict JSON deja present sur disque :
    # sans ce garde-fou, une experience qui plante herite du PASS de sa
    # precedente execution (un resultat perime compte comme valide).
    crashed = rc != 0
    if crashed:
        tail = stderr_s.strip().splitlines()
        verdict = None
        status = f"CRASH (rc={rc}): {tail[-1] if tail else 'no stderr'}"
        # the one-line status cannot localise a backend-dependent crash: the
        # full traceback is kept next to the summary, per experiment.
        crashdir = os.path.join(ROOT, "results", "crashes")
        os.makedirs(crashdir, exist_ok=True)
        with open(os.path.join(crashdir, os.path.basename(exp_dir) + ".log"),
                  "w") as cf:
            cf.write(f"rc={rc}  backend={env.get('CL10_BACKEND', 'auto')}\n")
            cf.write(f"script={script}\n\n--- stderr ---\n{stderr_s}\n")
            cf.write(f"--- stdout ---\n{stdout_s}\n")
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
        "returncode": rc,
        "crashed": crashed,
        "stale": False,
        "wall_time": round(time_, 3),
        "backend": reported_backend(verdict),
        "verdict": verdict,
    }


def main():
    ap = argparse.ArgumentParser(description="the CL10 validation runner")
    ap.add_argument("--jobs", type=int, default=min(12, os.cpu_count() or 4),
                    help="parallel experiment workers (subprocesses)")
    ap.add_argument("only", nargs="*", help="experiment dir names to run")
    ap.add_argument("--report-only", action="store_true",
                    help="roll up the verdict JSONs already on disk; run nothing")
    ap.add_argument("--no-bar", action="store_true",
                    help="plain one-line-per-experiment output, no progress bar")
    args = ap.parse_args()
    only = set(args.only)
    exps = [d for d in find_experiments() if not only or os.path.basename(d) in only]
    print("=" * 68)
    print("CL10  —  the verification program  (the GPU engine)")
    print("=" * 68)
    import cl10.engine as E
    print(f"  backend: {E._BACKEND}   jobs: {args.jobs}")
    # experiments are independent subprocesses writing to their own
    # results/ dir: run them concurrently, REPORT EACH AS IT FINISHES so a
    # slow experiment is visible live; everything is also logged to
    # results/run_log.txt.
    t_start = time.time()
    logp = os.path.join(ROOT, "results", "run_log.txt")
    # the live progress bar (stderr, so the verdict lines and the log stay
    # readable): 115 experiments from 0.03 s to 45 min, and a silent terminal
    # during the long tail reads as a hang.
    bar = None if args.no_bar or not sys.stderr.isatty() else tqdm(
        total=len(exps), desc="experiments", unit="exp", colour="green")
    with open(logp, "w") as logf:
        def emit(line):
            # through the bar when there is one: a bare print would tear it apart
            if bar:
                tqdm.write(line)
            else:
                print(line, flush=True)
            logf.write(line + "\n")
            logf.flush()
        done = []
        running = {}

        def tracked(d):
            """Run one experiment, registering it as in-flight so the bar can
            name it and show the progress bar it is drawing itself."""
            name = os.path.basename(d)
            slot = {"t0": time.time(), "progress": ""}
            running[name] = slot
            try:
                return report_one(d) if args.report_only else run_one(d, slot)
            finally:
                running.pop(name, None)
        def postfix():
            """Name what is still working, with the progress bar each child is
            drawing itself (the cone scans run for tens of minutes)."""
            if not bar:
                return
            now = time.time()
            live = sorted(running.items(), key=lambda kv: kv[1]["t0"])
            bar.set_postfix_str(
                "running: " + ", ".join(
                    f"{k} {int(now - s['t0'])}s"
                    + (f" [{s['progress']}]" if s["progress"] else "")
                    for k, s in live[:3]) if live else "idle")
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futs = [pool.submit(tracked, d) for d in exps]
            # the bar would otherwise move only when an experiment finishes: the
            # tail of the suite is three experiments over 55 minutes, and a line
            # frozen at 95% for 18 minutes reads as a hang.
            stop_tick = threading.Event()

            def tick():
                while not stop_tick.wait(2.0):
                    postfix()
            ticker = threading.Thread(target=tick, daemon=True)
            ticker.start()
            try:
                for fut in as_completed(futs):
                    r = fut.result()
                    done.append(r)
                    mark = "PASS" if r["pass"] else "FAIL"
                    tier = f"T{r['tier']}" if r["tier"] is not None else "   "
                    el = time.time() - t_start
                    b = r.get("backend")
                    btag = f" [{b}]" if b else ""
                    emit(f" [{el:7.1f}s] [{mark}]  {r['id']:34s} {tier} "
                         f"{r['status']}  ({r['wall_time']} s){btag}")
                    if bar:
                        bar.update(1)
                        postfix()
            finally:
                stop_tick.set()
                ticker.join(timeout=5)
            if bar:
                bar.close()
        if only:
            # a subset run must not shrink the rollup: the experiments that
            # were not run are folded in from the verdict JSON on disk.
            picked = {os.path.basename(d) for d in exps}
            done.extend(report_one(d) for d in find_experiments()
                        if os.path.basename(d) not in picked)
        results = sorted(done, key=lambda r: r["id"])
    n_pass = sum(1 for r in results if r["pass"])
    n = len(results)
    # the hard-tier verdict: tiers 0-3 all PASS - the one object, the
    # kinematics, the dynamics, the anomaly-free matter and the 5/3
    # normalisation.  Tier 4 (cosmology) and tier 5 (falsifiability) are
    # supporting.  The open problems O1-O9 of Section 16 are reported as the
    # named residuals below: they are listed, not hidden.
    hard = [r for r in results if r["tier"] in (0, 1, 2, 3)]
    hard_tier = bool(hard) and all(r["pass"] for r in hard)
    named = [r["id"] for r in results if r["residual"]]
    backends = [r.get("backend") for r in results]
    n_gpu = sum(1 for b in backends if b and b != "numpy")
    n_reported = sum(1 for b in backends if b)
    mismatch = [r["id"] for r in results
                if r["returncode"] is not None and r.get("backend")
                and r["backend"] != E._BACKEND]
    summary = {
        "n": n,
        "n_pass": n_pass,
        "pass": bool(n_pass == n),
        "backend_parent": E._BACKEND,
        "backend_requested": os.environ.get("CL10_BACKEND", "auto"),
        "n_backend_reported": n_reported,
        "n_backend_gpu": n_gpu,
        "hard_tier": hard_tier,
        "hard_tier_criterion": "all of tiers 0-3 PASS",
        "n_named_residuals": len(named),
        "named_residuals": named,
        "experiments": results,
    }
    outdir = os.path.join(ROOT, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print(f"  ->  {n_pass}/{n}  PASS  ->  results/validation_summary.json")
    vstr = ("hard tier MET (the root, the kinematics, the unitary causal "
            "dynamics, the anomaly-free matter, the 5/3)"
            if hard_tier else "hard tier FAILED (a tier 0-3 experiment failed)")
    print(f"  HARD-TIER VERDICT:  {vstr}")
    print(f"  BACKEND: parent {E._BACKEND} (requested {os.environ.get('CL10_BACKEND', 'auto')}); "
          f"{n_gpu}/{n_reported} experiments report a GPU backend")
    if mismatch:
        print(f"  OFF-CORE (an experiment pinned a different core than the parent, "
              f"by construction - the cone programme is CPU-only): {', '.join(mismatch)}")
    if named:
        print(f"  NAMED RESIDUALS ({len(named)}, the open problems O1-O9 of Section 16):")
        print(f"    {', '.join(named)}")
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main())
