#!/usr/bin/env bash
# run_all_gpus.sh - the whole verification tree, one command.
#
# Fans the 167 experiments across every GPU the host reports, one
# run_experiments.sh runner per card (CUDA_VISIBLE_DEVICES pins each runner to
# its card), and runs with the CPU core as well: the off-core suite (alt_*,
# V186-V192) executes on the host inside every runner.  When all halves are
# done it folds the on-disk verdicts into the complete summary with
# run_validation.py --report-only (no experiment re-runs) and prints the final
# verdict lines.
#
# usage:
#   ./run_all_gpus.sh              # all cards, 6 workers per card
#   JOBS=4 ./run_all_gpus.sh       # fewer workers per card
#   ./run_all_gpus.sh --dry-run    # show the split, launch nothing
set -u -o pipefail
cd "$(dirname "$0")"

JOBS_PER_CARD=${JOBS:-6}
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

# ---- the cards: one runner each, and the per-card memory gate --------------
mapfile -t GPUS < <(nvidia-smi --query-gpu=index,memory.used,memory.total \
                      --format=csv,noheader,nounits)
if [ "${#GPUS[@]}" -eq 0 ]; then
    echo "no GPU visible - use ./run_experiments.sh (CPU core) instead" >&2
    exit 1
fi
GATE=$(.venv/bin/python -c "from cl10 import config; print(config.CFG.backend.gpu_min_free_mib)" 2>/dev/null || echo 2048)
echo "== cards =="
NGPU=0
KEEP=()
for row in "${GPUS[@]}"; do
    IFS=', ' read -r idx used total <<<"$row"
    free=$(( total - used ))
    need=$(( GATE * JOBS_PER_CARD ))
    if [ "$free" -ge "$need" ]; then
        echo "  GPU$idx: ${free} MiB free  -> runner, ${JOBS_PER_CARD} workers"
        KEEP+=("$idx")
    else
        echo "  GPU$idx: ${free} MiB free  -> SKIPPED (needs ${need} for ${JOBS_PER_CARD} workers)"
    fi
    NGPU=$(( NGPU + 1 ))
done
if [ "${#KEEP[@]}" -eq 0 ]; then
    echo "ABORT: no card can host ${JOBS_PER_CARD} workers (${GATE} MiB each)." >&2
    echo "Free a card, lower JOBS (JOBS=$(( 1 )) fits a bare card), or use ./run_experiments.sh." >&2
    exit 1
fi

# ---- the split: interleaved so heavy scans spread over the cards ----------
mapfile -t ALL < <(find experiments -mindepth 3 -maxdepth 3 -name '*.py' \
                     | xargs -n1 dirname | sort -u | xargs -n1 basename)
echo "experiments: ${#ALL[@]}   cards used: ${#KEEP[@]}"
declare -a HALF
for c in "${!KEEP[@]}"; do HALF[$c]=""; done
for i in "${!ALL[@]}"; do
    c=$(( i % ${#KEEP[@]} ))
    HALF[$c]+="${ALL[i]} "
done

if [ "$DRY" -eq 1 ]; then
    for c in "${!KEEP[@]}"; do
        n=$(wc -w <<<"${HALF[$c]}")
        echo "  dry: GPU${KEEP[$c]} gets ${n} experiments"
    done
    echo "  dry: then ./run_experiments.sh --report-only"
    exit 0
fi

# ---- launch one runner per kept card --------------------------------------
TS=$(date +%Y%m%d_%H%M%S)
PIDS=(); LOGS=()
for c in "${!KEEP[@]}"; do
    log="results/run_all_gpus_gpu${KEEP[$c]}_${TS}.log"
    echo "  launching GPU${KEEP[$c]} runner -> ${log}"
    CUDA_VISIBLE_DEVICES="${KEEP[$c]}" JOBS="$JOBS_PER_CARD" \
        ./run_experiments.sh ${HALF[$c]} > "$log" 2>&1 &
    PIDS+=($!); LOGS+=("$log")
done

echo "running ${#PIDS[@]} halves - tail any log while waiting"
FAILHALF=0
for i in "${!PIDS[@]}"; do
    if ! wait "${PIDS[$i]}"; then
        echo "  half on GPU${KEEP[$i]} exited non-zero (see ${LOGS[$i]})"
        FAILHALF=1
    fi
done

# ---- fold every on-disk verdict into the complete summary ------------------
echo
echo "== final summary (report-only fold) =="
./run_experiments.sh --report-only | tail -14
[ "$FAILHALF" -eq 0 ] || echo "NOTE: at least one half exited non-zero - check the half logs above before trusting the fold."
