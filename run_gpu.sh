#!/usr/bin/env bash
# run_gpu.sh - the full 109-experiment validation on the GPU (cupy backend).
# Usage: ./run_gpu.sh          (run when the GPU is free)
set -uo pipefail
cd "$(dirname "$0")"

export CL10_BACKEND=cupy

echo "== GPU state =="
nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu --format=csv,noheader
used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
total=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
if [ "${used:-0}" -gt $(( total * 3 / 10 )) ]; then
    echo "WARNING: GPU already ${used}/${total} MiB used - another process"
    echo "may OOM cupy (this killed the last run). Ctrl-C now if so."
fi

echo
echo "== running run_validation.py (cupy) =="
.venv/bin/python run_validation.py
rc=$?

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
