#!/usr/bin/env bash
# run_experiments.sh - the whole CL10 experiment suite, in one run.
#
#   ./run_experiments.sh                    # every experiment under experiments/
#   ./run_experiments.sh V63_bz_chirality V102_lorentz_isotropy   # a subset (**exact** dir names)
#   JOBS=6 ./run_experiments.sh             # worker pool (default: pool.njobs in cl10/config.yaml)
#   FORCE=1 ./run_experiments.sh            # launch even if the card looks busy (the run then takes the CPU core)
#
# Each experiment is a subprocess.  This script checks the card first, then pins
# CL10_BACKEND=cupy for the whole run: a forced backend does NOT re-gate per
# child and does NOT fall back, so a worker that cannot get the card fails
# loudly (traceback in results/crashes/) instead of silently running the CPU
# path and reporting a GPU verdict it never earned.  Experiments that are CPU
# work by measurement (the cone programme, V186/V187) pin the CPU core
# themselves and are listed OFF-CORE in the summary.
set -uo pipefail
cd "$(dirname "$0")"

# the worker count and the free-memory gate come from cl10/config.yaml
_njobs() { .venv/bin/python -c "from cl10 import config; print(config.CFG.pool.njobs_effective)" 2>/dev/null || echo 8; }
_gate()  { .venv/bin/python -c "from cl10 import config; print(config.CFG.backend.gpu_min_free_mib)" 2>/dev/null || echo 2048; }
JOBS=${JOBS:-$(_njobs)}
FORCE=${FORCE:-0}
GATE=$(_gate)

mapfile -t GPUS < <(nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
                    --format=csv,noheader,nounits)
echo "== GPU =="
printf '%s\n' "${GPUS[@]}"
free0=$(awk -F', ' 'NR==1{print $3-$2}' <<<"${GPUS[0]}")
need=$(( GATE * JOBS ))              # the engine gate: this much free per worker
nexp=$(find experiments -mindepth 3 -maxdepth 3 -name '*.py' | xargs -n1 dirname | sort -u | wc -l)
echo "experiments: ${nexp}"

if [ "${free0}" -lt "${need}" ] && [ "${FORCE}" != 1 ]; then
    fit=$(( free0 / GATE ))
    echo
    echo "ABORT: ${free0} MiB free on GPU0 - cannot host ${JOBS} cupy workers"
    echo "(the engine gate needs ${GATE} MiB each; ${need} MiB required). The engine"
    echo "would silently fall back to numpy."
    if [ "${fit}" -ge 1 ]; then
        echo "Wait for the card, or lower JOBS (JOBS=${fit} fits now);"
    else
        echo "No worker fits on the card right now - wait for it to free up;"
    fi
    echo "set FORCE=1 to run anyway (numpy fallback per experiment)."
    exit 3
fi

# the backend follows the card: cupy only when it can host the workers, else the
# CPU core (both the engine gate and the direct-cupy experiments read CL10_BACKEND)
if [ "${free0}" -ge "${need}" ]; then
    export CL10_BACKEND=cupy
else
    export CL10_BACKEND=numpy
    echo "NOTE: ${free0} MiB free on GPU0 < ${need} MiB needed - running on the CPU core (CL10_BACKEND=numpy)."
fi

log="results/experiments_run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p results
echo
echo "== run_validation.py --jobs ${JOBS}   (log: ${log}) =="
# stdout (the per-experiment verdict lines) is tee'd into the log; stderr (the
# tqdm progress bar) stays on the terminal, so the bar redraws in place instead
# of filling the log with its frames.  --no-bar on the command line opts out.
.venv/bin/python run_validation.py --jobs "${JOBS}" "$@" > >(tee "${log}")
rc=$?

echo
echo "== summary (results/validation_summary.json) =="
.venv/bin/python - <<'EOF'
import json
d = json.load(open("results/validation_summary.json"))
print(f"{d['n_pass']}/{d['n']} PASS   backend-wide pass = {d['pass']}")
print("hard tier (0-3):", d["hard_tier"], "-", d["hard_tier_criterion"])
print("named residuals:", ", ".join(d["named_residuals"]))
fails = [e for e in d["experiments"] if not e.get("pass", False)]
if fails:
    print("FAILING:")
    for e in fails:
        print(f"  {e['id']:34s} rc={e.get('returncode')}  {str(e['status'])[:110]}")
slow = sorted(d["experiments"], key=lambda e: -e.get("wall_time", 0))[:8]
print("slowest:", ", ".join(f"{e['id']}({e['wall_time']}s)" for e in slow))
EOF
echo
echo "log: ${log}"
exit $rc
