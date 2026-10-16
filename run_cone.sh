#!/usr/bin/env bash
# run_cone.sh - the cone programme (V185 separation lemma, V186 cone
# inverse problem) on the GPU.  Mirrors run_gpu.sh conventions.
#
#   ./run_cone.sh                # V185 then V186, GPU auto-picked
#   GPU=1 ./run_cone.sh          # pin a card
#   FORCE=1 ./run_cone.sh        # run even if the card looks busy
#   ONLY=V186 ./run_cone.sh      # one experiment
#
# The engine (cl10/engine.py) commits to cupy only when the card has
# >= 2 GiB free; otherwise it silently falls back to numpy.  This script
# checks the card first and says so out loud - free a card first (the
# LLM servers hold both 24 GiB cards most of the day).
set -uo pipefail
cd "$(dirname "$0")"

PY=.venv/bin/python
declare -A EXP=(
  [V185]="experiments/dynamics_spectrum/V185_separation_lemma/V185_separation_lemma.py"
  [V186]="experiments/dynamics_spectrum/V186_cone_inverse/V186_cone_inverse.py"
)
ORDER="V185 V186"
if [ -n "${ONLY:-}" ]; then ORDER="$ONLY"; fi

echo "== GPUs =="
nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu \
           --format=csv,noheader
if [ -z "${GPU:-}" ]; then
    GPU=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits \
          | sort -t, -k2 -n | tail -1 | cut -d, -f1 | tr -d ' ')
fi
free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits \
       -i "$GPU")
echo
echo "selected GPU $GPU: ${free} MiB free"
if [ "$free" -lt 2048 ] && [ "${FORCE:-0}" != 1 ]; then
    echo "ABORT: the engine needs >= 2048 MiB free on GPU $GPU (it would"
    echo "silently run the numpy path). Free a card, or FORCE=1 to run"
    echo "on CPU anyway."
    exit 3
fi
export CUDA_VISIBLE_DEVICES="$GPU"

logdir="results/cone_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$logdir"
rc_all=0
for tag in $ORDER; do
    exp="${EXP[$tag]:-}"
    if [ -z "$exp" ]; then echo "unknown experiment $tag"; rc_all=1; continue; fi
    echo
    echo "== $tag  ($(basename "$exp")) -> ${logdir}/${tag}.log =="
    /usr/bin/time -v "$PY" "$exp" > "${logdir}/${tag}.log" 2> "${logdir}/${tag}.time"
    rc=$?
    tail -1 "${logdir}/${tag}.time" | sed 's/^/    wall: /'
    grep -E "VERDICT|\[PASS\]|\[FAIL\]" "${logdir}/${tag}.log" | sed 's/^/    /'
    [ $rc -ne 0 ] && rc_all=$rc
done

echo
echo "logs: $logdir"
exit $rc_all
