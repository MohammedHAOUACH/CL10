#!/usr/bin/env bash
# run_gpu.sh - the alt_step experiments on the GPU.
# V184/V185 already have PASS JSONs; the remaining chain is V186
# (joint cell search, ~4 h GPU) then V187 (gate, ~2 min).
# Override the list:  ./run_gpu.sh V184_step_loss_landscape V185_...
# Logs: /tmp/alt_step_gpu/<exp>.log ; stops at the first failure.
set -u
cd ~/Desktop/code/CL10
source .venv/bin/activate

export CL10_BACKEND=numpy

# pick the GPU with the most free memory; abort if none has > 4 GiB
best=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits \
       | sort -t, -k2 -nr | head -1)
idx=$(echo "$best" | cut -d, -f1 | tr -d ' ')
free=$(echo "$best" | cut -d, -f2 | tr -d ' ')
if [ "$free" -lt 4096 ]; then
  echo "no GPU with >4 GiB free (best: idx $idx, $free MiB) - aborting"
  exit 1
fi
echo "using GPU $idx ($free MiB free)"
export CUDA_VISIBLE_DEVICES=$idx
export JAX_PLATFORM_NAME=cuda
nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader
python3 - <<'EOF' || { echo "GPU not visible to JAX - aborting"; exit 1; }
import jax
jax.config.update("jax_enable_x64", True)
print("JAX devices:", jax.devices())
assert jax.devices()[0].platform in ("gpu", "cuda"), "no CUDA device visible"
EOF

mkdir -p /tmp/alt_step_gpu
EXPS=${*:-V186_cell_step_joint_search V187_step_search_gate}
for V in $EXPS; do
  echo "==== $V ===="
  python3 "experiments/alt_step/$V/$V.py" 2>&1 \
    | tee "/tmp/alt_step_gpu/${V,,}.log"
  rc=${PIPESTATUS[0]}
  echo "---- $V exit $rc ----"
  [[ $rc -ne 0 ]] && break
done
echo "logs: /tmp/alt_step_gpu/"
