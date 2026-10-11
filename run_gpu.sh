#!/usr/bin/env bash
# run_gpu.sh - the CL10 validation suite on the GPU (cupy backend).
#
#   ./run_gpu.sh                    # every experiment under experiments/
#   ./run_gpu.sh V63_bz_chirality V102_lorentz_isotropy   # a subset (**exact** dir names)
#   JOBS=6 ./run_gpu.sh             # worker pool (default 8)
#   FORCE=1 ./run_gpu.sh            # launch even if the card looks busy
#
# Each experiment is a subprocess; cl10/engine.py probes the card and falls
# back to numpy if fewer than 2 GiB are free, so a busy card silently runs
# the CPU path. This script checks the card first and says so out loud.
set -uo pipefail
cd "$(dirname "$0")"

export CL10_BACKEND=cupy
JOBS=${JOBS:-8}
FORCE=${FORCE:-0}

mapfile -t GPUS < <(nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
                    --format=csv,noheader,nounits)
echo "== GPU =="
printf '%s\n' "${GPUS[@]}"
free0=$(awk -F', ' 'NR==1{print $3-$2}' <<<"${GPUS[0]}")
need=$(( 2048 * JOBS ))            # the engine gate: 2 GiB free per worker
nexp=$(find experiments -mindepth 2 -maxdepth 2 -name '*.py' | xargs -n1 dirname | sort -u | wc -l)
echo "experiments: ${nexp}"

if [ "${free0}" -lt "${need}" ] && [ "${FORCE}" != 1 ]; then
    fit=$(( free0 / 2048 ))
    echo
    echo "ABORT: ${free0} MiB free on GPU0 - cannot host ${JOBS} cupy workers"
    echo "(the engine gate needs 2048 MiB each; ${need} MiB required). The engine"
    echo "would silently fall back to numpy."
    if [ "${fit}" -ge 1 ]; then
        echo "Wait for the card, or lower JOBS (JOBS=${fit} fits now);"
    else
        echo "No worker fits on the card right now - wait for it to free up;"
    fi
    echo "set FORCE=1 to run anyway (numpy fallback per experiment)."
    exit 3
fi

log="results/gpu_run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p results
echo
echo "== run_validation.py --jobs ${JOBS}   (log: ${log}) =="
.venv/bin/python run_validation.py --jobs "${JOBS}" "$@" 2>&1 | tee "${log}"
rc=${PIPESTATUS[0]}

echo
echo "== summary (results/validation_summary.json) =="
.venv/bin/python - <<'EOF'
import json
d = json.load(open("results/validation_summary.json"))
print(f"{d['n_pass']}/{d['n']} PASS   backend-wide pass = {d['pass']}")
print("U10 (tiers 0-3):", d["u10"], "-", d["u10_criterion"])
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
