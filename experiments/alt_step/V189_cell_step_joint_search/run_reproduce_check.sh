#!/usr/bin/env bash
# V189 - reproduce the joint cell x step search's objective on the card.
#
#   ./run_reproduce_check.sh              # the (3,3,3) cell (M = 432), the one
#                                         # the carried verdict has not re-checked
#   ./run_reproduce_check.sh --all        # the four cells
#   WAIT=600 ./run_reproduce_check.sh     # wait up to 10 min for the card
#   FORCE=1 ./run_reproduce_check.sh      # run even if the card looks busy
#   CL10_JAX_CPU=1 ./run_reproduce_check.sh   # the CPU core (minutes, not seconds)
#
# The search itself (4 h of GPU) is not re-run: this re-evaluates the objective
# at the (theta*, q*) the verdict recorded, which is what decides whether the
# carried numbers are the numbers the current library produces.
set -uo pipefail
cd "$(dirname "$0")/../../.."   # the repo root: .venv and cl10 live there

NEED=${NEED:-$([ -x .venv/bin/python ] && .venv/bin/python -c \
    "from cl10 import config; print(config.CFG.backend.gpu_min_free_mib)" \
    2>/dev/null || echo 2048)}
WAIT=${WAIT:-0}
FORCE=${FORCE:-0}
# coarse gate: is any card roomy enough?  JAX orders its own devices, so the
# card it takes is not necessarily the one the engine pins for cupy; the
# programme repeats this check on the device it actually got and exits loudly
# there, which is the authoritative one.
gate() {
    nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits \
        | sort -rn | head -1
}

if [ "${CL10_JAX_CPU:-0}" != "1" ]; then
    echo "== cards =="
    nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
               --format=csv,noheader,nounits
    waited=0
    while :; do
        free=$(gate)
        if [ "${free}" -ge "${NEED}" ] || [ "${FORCE}" = "1" ]; then
            echo "roomiest card: ${free} MiB free (gate ${NEED} MiB) - launching"
            break
        fi
        if [ "${waited}" -ge "${WAIT}" ]; then
            echo
            echo "ABORT: ${free} MiB free on the roomiest card; the check needs"
            echo "${NEED} MiB for one JAX context. Free a card, raise NEED, set"
            echo "FORCE=1 to run anyway, WAIT=<s> to wait, or CL10_JAX_CPU=1 to"
            echo "check on the CPU core."
            exit 3
        fi
        echo "waiting for a card: ${free} MiB free, need ${NEED} (${waited}s / ${WAIT}s)"
        sleep 15
        waited=$((waited + 15))
    done
fi

log="experiments/alt_step/V189_cell_step_joint_search/results/reproduce_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "${log}")"
echo
echo "== V189_reproduce_objective.py $*   (log: ${log}) =="
.venv/bin/python experiments/alt_step/V189_cell_step_joint_search/V189_reproduce_objective.py "$@" 2>&1 | tee "${log}"
rc=${PIPESTATUS[0]}

echo
echo "== verdict =="
grep -E "^REPRODUCE:" "${log}" || echo "no verdict line (rc=${rc})"
exit "${rc}"
