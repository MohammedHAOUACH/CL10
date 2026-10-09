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

# sanity: the suite needs the card nearly free (it peaked at ~23 GB);
# refuse to start while another process holds the memory
FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
if [ "${FREE:-0}" -lt 20000 ]; then
  echo "GPU has only ${FREE} MiB free - another process holds the card:"
  nvidia-smi --query-compute-apps=pid,used_memory --format=csv
  echo "free the GPU (or wait), then re-run ./run_gpu.sh"; exit 1
fi

# the 104 experiments (V0..V103), name order, fresh verdicts
# (V98/V99 read CL10_BACKEND too: their Lanczos matvecs and dense
#  diagonalizations run on the GPU)
.venv/bin/python run_validation.py 2>&1 | tee "$LOG"
rc=$?

echo
echo "exit code: $rc"
echo "summary:   $(pwd)/results/validation_summary.json"
tail -n 25 "$LOG"
exit $rc
