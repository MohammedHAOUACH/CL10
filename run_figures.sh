#!/usr/bin/env bash
# run_figures.sh - the CL10 publication figures on the GPU (cupy backend).
#
#   ./run_figures.sh                     # all 25 figures -> doc/figures/
#   ./run_figures.sh --only 22 23 24     # a subset of figure numbers
#   ./run_figures.sh --fast              # skip the young double-slit (fig 1)
#   ./run_figures.sh --cpu               # force the numpy backend on purpose
#   FORCE=1 ./run_figures.sh             # launch even if the card looks busy
#   NEED=4096 ./run_figures.sh           # raise the free-memory gate (MiB)
#
# make_figures pins the card from simulation/young/config.yaml (gpu: 1) and
# falls back to numpy below 2 GiB free. This script checks the pinned card
# first and says so out loud, so a busy card never silently produces CPU
# figures that disagree with the GPU validation run.
set -uo pipefail
cd "$(dirname "$0")"

GPU=$(awk '/^gpu:/{print $2}' simulation/young/config.yaml 2>/dev/null)
GPU=${GPU:-1}
FORCE=${FORCE:-0}

# the free-memory gate comes from cl10/config.yaml (the engine's own gate)
_gate() {
    .venv/bin/python -c "from cl10 import config; print(config.CFG.backend.gpu_min_free_mib)" \
        2>/dev/null || echo 2048
}
NEED=${NEED:-$(_gate)}        # the same gate the engine uses: backend.gpu_min_free_mib

echo "== GPU ${GPU} (the card make_figures pins) =="
nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
           --format=csv,noheader,nounits -i "${GPU}"
free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits -i "${GPU}")

if [ "${free}" -lt "${NEED}" ] && [ "${FORCE}" != 1 ]; then
    echo
    echo "ABORT: ${free} MiB free on GPU${GPU}; the gate needs ${NEED} MiB."
    echo "Below it make_figures silently switches to numpy, and the figures"
    echo "would no longer be the GPU products of the validation run."
    echo "Free the card, raise/lower NEED, set FORCE=1 to run anyway,"
    echo "or pass --cpu to choose numpy deliberately."
    exit 3
fi

log="results/figures_run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p results
echo
echo "== python -m make_figures $*   (log: ${log}) =="
# stdout (the per-figure log) is tee'd; stderr (the tqdm progress bar) stays on
# the terminal so the bar redraws in place instead of filling the log.
.venv/bin/python -m make_figures "$@" > >(tee "${log}")
rc=$?

echo
echo "== figures =="
ls -1 doc/figures/*.png | wc -l | tr -d ' ' | sed 's/^/png files: /'
ls -lt --time-style=+%H:%M:%S doc/figures/*.png | head -5
echo "log: ${log}"
exit ${rc}
