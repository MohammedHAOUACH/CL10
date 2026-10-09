#!/usr/bin/env bash
# run_gpu.sh - the full 117-experiment validation on the GPU (cupy backend).
# Usage: ./run_gpu.sh          (run when the GPU is free)
set -uo pipefail
cd "$(dirname "$0")"

export CL10_BACKEND=cupy

echo "== GPU state =="
nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu --format=csv,noheader
used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
total=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
if [ "${used:-0}" -gt $(( total * 3 / 10 )) ]; then
    echo "NOTE: GPU already ${used}/${total} MiB used - the engine probe"
    echo "falls back to numpy per experiment instead of OOM-crashing."
fi

echo
echo "== running run_validation.py (cupy, ${JOBS:-12} parallel jobs) =="
.venv/bin/python run_validation.py --jobs "${JOBS:-12}" 2>&1 | tee results/run_console.txt
rc=${PIPESTATUS[0]}

echo
echo "== summary (results/validation_summary.json) =="
.venv/bin/python - <<'EOF'
import json
d = json.load(open("results/validation_summary.json"))
print(f"{d['n_pass']}/{d['n']} PASS   backend-wide pass = {d['pass']}")
print("U10:", d["u10"])
print("named residuals:", d["named_residuals"])
fails = [e for e in d["experiments"] if not e.get("pass", False)]
if fails:
    print("FAILING:", ", ".join(e["id"] for e in fails))
EOF
exit $rc
