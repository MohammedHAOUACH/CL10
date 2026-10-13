#!/usr/bin/env bash
# Run ONLY the new alt_interactions census (V171-V174) - the 134 already
# validated experiments are NOT re-run. Each script writes its own verdict
# JSON; the final --report-only just rolls up the existing JSONs (fast,
# no recomputation).
set -euo pipefail
cd "$(dirname "$0")"
source .venv/bin/activate
for v in V171_collision_vertex V172_self_energy_flow V173_block_rg V174_isotropy_gate; do
    echo "=== $v ==="
    python3 "experiments/alt_interactions/$v/$v.py"
done
python3 run_validation.py --report-only
