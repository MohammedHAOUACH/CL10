#!/usr/bin/env bash
# run_V131.sh - the algebraic-closure experiment (V131) on the GPU (cupy backend).
#
#   ./run_V131.sh              # V131 only, JOBS=1 (single experiment, no pool)
#   ./run_V131.sh --full       # every experiment under experiments/
#   FORCE=1 ./run_V131.sh      # launch even if the card looks busy
#
# V131 closes the last open structural gap of the publication: the root
# quantizes its own step (the 0-parameter completion, Section 8.5).  It is a
# single-process experiment, so it needs only the engine's 2 GiB free-memory
# gate once, not once per worker.
set -uo pipefail
cd "$(dirname "$0")"

if [ "${1:-}" = "--full" ]; then
    shift
    exec ./run_gpu.sh "$@"
fi

export CL10_BACKEND=cupy
JOBS=${JOBS:-1}
FORCE=${FORCE:-0}

if command -v nvidia-smi >/dev/null 2>&1; then
    echo "== GPU =="
    nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
               --format=csv,noheader,nounits
    free0=$(nvidia-smi --query-gpu=memory.total,memory.used --format=csv,noheader,nounits | \
            awk -F', ' 'NR==1{print $1-$2}')
    if [ "${free0}" -lt 2048 ] && [ "${FORCE}" != 1 ]; then
        echo
        echo "ABORT: ${free0} MiB free on GPU0 - below the engine's 2 GiB gate."
        echo "The engine would silently fall back to numpy."
        echo "Wait for the card, or set FORCE=1 to run anyway (numpy path)."
        exit 3
    fi
fi

log="results/v131_gpu_$(date +%Y%m%d_%H%M%S).log"
mkdir -p results
echo "== V131_algebraic_closure  (log: ${log}) =="
.venv/bin/python -c "from cl10 import engine as E; print('engine backend:', E._BACKEND)"
.venv/bin/python experiments/V131_algebraic_closure/V131_algebraic_closure.py 2>&1 | tee "${log}"
rc=${PIPESTATUS[0]}

echo
if [ "${rc}" -ne 0 ]; then
    echo "== the run FAILED - the results JSON on disk is from an EARLIER run, not this one =="
else
    echo "== results (experiments/V131_algebraic_closure/results/v131_algebraic_closure.json) =="
    .venv/bin/python - <<'EOF'
import json
d = json.load(open("experiments/V131_algebraic_closure/results/v131_algebraic_closure.json"))
n = len(d["checks"])
p = sum(1 for c in d["checks"] if c["pass"])
print(f"{p}/{n} PASS")
for c in d["checks"]:
    if not c["pass"]:
        print("FAIL:", c["name"])
print("closed angles:", sorted(d["data"]["closure_scan"]["exactly_closed"].keys()),
      "(pi/2 order 4, pi/4 order 8)")
EOF
fi
echo
echo "log: ${log}"
exit $rc
