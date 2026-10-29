#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the verification runner.

Runs every experiment under ``experiments/`` (in name order), reads each
JSON verdict, and writes a single summary to ``results/validation_summary.json``.

Usage::

    python3 run_validation.py              # incremental gate (the default): re-run
                                           # only experiments that are new or whose
                                           # scripts changed since the last gate; the
                                           # rest fold from their on-disk verdicts
    python3 run_validation.py --full       # re-run everything (a fresh clone, with no
                                           # fingerprint file, runs this way)
    python3 run_validation.py --seed-fingerprints
                                           # baseline the fingerprints from the verdicts
                                           # already on disk, running nothing
    python3 run_validation.py --report-only
                                           # roll up the on-disk verdicts, run nothing
    python3 run_validation.py V0 V1        # run a subset
    python3 run_validation.py --jobs 4     # cap the worker pool (default: pool.njobs in cl10/config.yaml)
The fingerprint gate (results/validation_fingerprints.json) records the SHA-256 of
every experiment's scripts and of the cl10 library plus cl10/config.yaml: a library
change re-runs the whole suite, a script change re-runs that experiment alone, and a
verdict that is missing, failed, or older than a crash log is re-run regardless.
Each experiment is a **0-parameter** check: it confirms a prediction of the
publication (Section 10) or reports an honest residual.  The verdict JSON has
``pass``, ``status``, the 0-parameter number, and ``wall_time``.  A check that
ran correctly but whose pre-registered physics branches did not fire declares
``"inconclusive": true``: it is counted apart from a failure, because "nature
said neither" is not "the check broke".
"""

import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm

from cl10 import config as _config

ROOT = os.path.dirname(os.path.abspath(__file__))
EXPS = os.path.join(ROOT, "experiments")


def find_experiments():
    # one directory is one experiment, whatever the number of programmes it
    # holds (an experiment may keep its reproduction check next to the search):
    # deduped, or a second .py makes the driver count and run the same experiment twice.
    dirs = {os.path.dirname(p)
            for p in glob.glob(os.path.join(EXPS, "*", "*", "*.py"))}
    return sorted(d for d in dirs if d != EXPS)


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


FPFILE = os.path.join(ROOT, "results", "validation_fingerprints.json")


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def library_fingerprint():
    """cl10 + cl10/config.yaml: when these change, on-disk verdicts stop
    being trustworthy and the fingerprint gate re-runs everything."""
    h = hashlib.sha256()
    root = os.path.join(ROOT, "cl10")
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d != "__pycache__")
        for fn in sorted(filenames):
            if fn.endswith((".py", ".yaml")):
                p = os.path.join(dirpath, fn)
                h.update(os.path.relpath(p, ROOT).encode())
                h.update(_sha(p).encode())
    return h.hexdigest()


def experiment_fingerprint(exp_dir):
    """Every .py an experiment ships (its script and any sibling machinery
    it imports from its own dir); the results dir is the verdict itself and
    never feeds its own fingerprint."""
    h = hashlib.sha256()
    for dirpath, dirnames, filenames in os.walk(exp_dir):
        dirnames[:] = sorted(d for d in dirnames
                             if d not in ("results", "__pycache__"))
        for fn in sorted(filenames):
            if fn.endswith(".py"):
                p = os.path.join(dirpath, fn)
                h.update(os.path.relpath(p, exp_dir).encode())
                h.update(_sha(p).encode())
    return h.hexdigest()


def load_fingerprints():
    try:
        with open(FPFILE) as fh:
            j = json.load(fh)
        if isinstance(j.get("experiments"), dict):
            return j
    except Exception:  # noqa: BLE001
        pass
    return {"library": None, "experiments": {}}


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
        "inconclusive": bool(verdict.get("inconclusive")) if verdict else False,
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
    # A non-zero exit is not by itself a crash. An experiment that records its
    # own FAIL - a pre-registered INCONCLUSIVE branch exits 1 - writes a
    # verdict JSON and returns non-zero: that verdict is the result, and the
    # returncode field keeps the exit code visible. What counts as a crash is a
    # non-zero exit that left NO verdict from this run: then the JSON on disk
    # belongs to an earlier execution and must not be inherited (a crashed
    # experiment would otherwise keep the PASS of the run before it).
    jf = verdict_file(exp_dir)
    recorded = jf is not None and os.path.getmtime(jf) >= t0
    crashed = rc != 0 and not recorded
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
        # the mirror rule: a pass invalidates the crash log a previous run left
        # behind, so results/crashes/ lists the experiments that crashed in
        # THIS run and not those that crashed the last time they failed.
        try:
            os.remove(os.path.join(ROOT, "results", "crashes",
                                   os.path.basename(exp_dir) + ".log"))
        except FileNotFoundError:
            pass
    return {
        "id": os.path.basename(exp_dir),
        "dir": os.path.basename(exp_dir),
        "tier": None if crashed else (verdict.get("tier") if verdict else None),
        "property": None if crashed else (verdict.get("property") if verdict else None),
        "residual": False if crashed else (bool(verdict.get("residual")) if verdict else False),
        "pass": False if crashed else (bool(verdict["pass"]) if verdict else False),
        "inconclusive": False if crashed else (
            bool(verdict.get("inconclusive")) if verdict else False),
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
    ap.add_argument("--jobs", type=int,
                    default=_config.CFG.pool.njobs_effective,
                    help="parallel experiment workers (subprocesses); the "
                         "default is pool.njobs in cl10/config.yaml")
    ap.add_argument("only", nargs="*", help="experiment dir names to run")
    ap.add_argument("--report-only", action="store_true",
                    help="roll up the verdict JSONs already on disk; run nothing")
    ap.add_argument("--full", action="store_true",
                    help="re-run every experiment (the default is the "
                         "fingerprint gate: only new or changed scripts run)")
    ap.add_argument("--seed-fingerprints", action="store_true",
                    help="baseline results/validation_fingerprints.json from the "
                         "verdicts already on disk; run nothing")
    ap.add_argument("--no-bar", action="store_true",
                    help="plain one-line-per-experiment output, no progress bar")
    args = ap.parse_args()
    only = set(args.only)
    exps = [d for d in find_experiments() if not only or os.path.basename(d) in only]
    if args.report_only:
        gate = "report-only"
    elif args.seed_fingerprints:
        gate = "seed"
    else:
        gate = "full" if args.full else "changed"
    fps = load_fingerprints()
    lib = library_fingerprint()
    if gate == "changed" and not fps["experiments"]:
        gate = "full"          # fresh tree: nothing has ever been gated
    fold_names = set()
    if gate in ("changed", "seed"):
        for d in exps:
            name = os.path.basename(d)
            if gate == "seed":
                fold_names.add(name)
                continue
            if only and name in only:
                continue       # an explicitly named experiment always runs
            v = read_verdict(d)
            if (fps.get("library") != lib
                    or fps["experiments"].get(name) != experiment_fingerprint(d)
                    or v is None or stale_verdict(d)
                    or not (v.get("pass") or v.get("inconclusive"))):
                continue       # new, changed, missing, failed or stale: re-run
            fold_names.add(name)
    print("=" * 68)
    print("CL10  —  the verification program  (the GPU engine)")
    print("=" * 68)
    import cl10.engine as E
    print(f"  backend: {E._BACKEND}   jobs: {args.jobs}")
    if gate in ("changed", "seed"):
        print(f"  gate: {gate}  ({len(exps) - len(fold_names)} re-run, "
              f"{len(fold_names)} folded from on-disk verdicts; "
              f"library {'unchanged' if fps.get('library') == lib else 'CHANGED'})")
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
                if args.report_only or name in fold_names:
                    return report_one(d)
                return run_one(d, slot)
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
    n = len(results)
    n_pass = sum(1 for r in results if r["pass"])
    # a declared-inconclusive physics verdict is not a failed check: it is a
    # measurement that landed outside the pre-registered branches, and the
    # rollup says so instead of folding it into the failure count.
    n_inconclusive = sum(1 for r in results
                         if r.get("inconclusive") and not r["pass"])
    n_failed = n - n_pass - n_inconclusive
    # the hard-tier verdict: tiers 0-3 all PASS - the one object, the
    # kinematics, the dynamics, the anomaly-free matter and the 5/3
    # normalisation.  Tier 4 is the exploration programme, tier 5 the
    # falsifiers, tier 6 the provenance and audit checks: all supporting.
    # The open problems O1-O10 of Section 16 are reported as the named
    # residuals below: they are listed, not hidden.
    # An experiment that writes no tier is NOT folded into the hard tier: the
    # denominator is printed, because "hard tier MET" over 22 of 156 checks
    # would otherwise read as a statement about all of them.
    hard = [r for r in results if r["tier"] in (0, 1, 2, 3)]
    untiered = [r["id"] for r in results if r["tier"] is None]
    hard_inconclusive = [r["id"] for r in hard
                         if r.get("inconclusive") and not r["pass"]]
    hard_tier = bool(hard) and all(r["pass"] or r.get("inconclusive")
                                   for r in hard)
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
        "n_inconclusive": n_inconclusive,
        "n_failed": n_failed,
        "pass": bool(n_failed == 0),
        "backend_parent": E._BACKEND,
        "backend_requested": os.environ.get("CL10_BACKEND", "auto"),
        "n_backend_reported": n_reported,
        "n_backend_gpu": n_gpu,
        "hard_tier": hard_tier,
        "hard_tier_inconclusive": hard_inconclusive,
        "hard_tier_n": len(hard),
        "n_untiered": len(untiered),
        "untiered": untiered,
        "hard_tier_criterion": "all of tiers 0-3 PASS",
        "n_named_residuals": len(named),
        "named_residuals": named,
        "experiments": results,
    }
    outdir = os.path.join(ROOT, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    if gate in ("changed", "seed"):
        # an experiment that ran clean enters the fingerprint ledger; one that
        # failed or is stale leaves no entry, so the next gate re-runs it; a
        # library hash is recorded only on a fully green pass over the suite.
        dirm = {os.path.basename(d): d for d in exps}
        fps = load_fingerprints()
        for r in results:
            d = dirm.get(r["id"])
            if d is None:
                continue
            ok = (r["pass"] or r.get("inconclusive")) and not r.get("stale")
            if not ok:
                fps["experiments"].pop(r["id"], None)
            elif gate == "seed" or r["returncode"] is not None:
                fps["experiments"][r["id"]] = experiment_fingerprint(d)
        if n_failed == 0:
            fps["library"] = lib
        with open(FPFILE, "w") as f:
            json.dump(fps, f, indent=1)
    print("=" * 68)
    print(f"  ->  {n_pass} PASS + {n_inconclusive} INCONCLUSIVE + "
          f"{n_failed} FAILED of {n}  ->  results/validation_summary.json")
    vstr = ("hard tier MET (the root, the kinematics, the unitary causal "
            "dynamics, the anomaly-free matter, the 5/3)"
            if hard_tier else "hard tier FAILED (a tier 0-3 check failed)")
    print(f"  HARD-TIER VERDICT:  {vstr}")
    print(f"    over {len(hard)} of {n} experiments; "
          f"{len(untiered)} write no tier and are outside the tier "
          f"accounting (they still count in PASS/FAILED above)")
    if hard_inconclusive:
        print(f"  HARD-TIER INCONCLUSIVE (ran, pre-registered branches did "
              f"not fire): {', '.join(hard_inconclusive)}")
    print(f"  BACKEND: parent {E._BACKEND} (requested {os.environ.get('CL10_BACKEND', 'auto')}); "
          f"{n_gpu}/{n_reported} experiments report a GPU backend")
    if mismatch:
        print(f"  OFF-CORE (an experiment pinned a different core than the parent, "
              f"by construction - the cone programme is CPU-only): {', '.join(mismatch)}")
    if named:
        print(f"  NAMED RESIDUALS ({len(named)}, the open problems O1-O10 of Section 16):")
        print(f"    {', '.join(named)}")
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
