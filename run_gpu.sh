#!/usr/bin/env bash
# CL10 - full validation suite on the GPU (cupy backend).
# Frees nothing itself: make sure the GPU is idle, then run:
#   ./run_gpu.sh
# The log lands in gpu_run_<timestamp>.log; the verdict summary lands in
# results/validation_summary.json (read that back to me).
set -uo pipefail
cd "$(dirname "$0")"

export CL10_BACKEND=cupy
LOG="gpu_run_$(date +%Y%m%d_%H%M%S).log"
echo "backend: cupy | log: $LOG"

# sanity: cupy must see the GPU
.venv/bin/python -c "import cupy as cp; print('GPU:', cp.cuda.runtime.getDeviceProperties(0)['name'])" || {
  echo "cupy cannot see the GPU - aborting"; exit 1; }

# the 86 experiments (V0..V85), name order, fresh verdicts
.venv/bin/python run_validation.py 2>&1 | tee "$LOG"
rc=$?

echo
echo "exit code: $rc"
echo "summary:   $(pwd)/results/validation_summary.json"
tail -n 25 "$LOG"
exit $rc
