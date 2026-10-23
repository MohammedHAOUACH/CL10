#!/usr/bin/env bash
# run_altgpu.sh - GPU verification of the alt-interactions grid pass
# (numba CUDA Jacobi-Schur in cl10/altgpu.py). Run this when the GPU is
# free: it checks the batched CUDA eigensolver against numpy and re-runs
# the V172 Fock self-energy grid on the GPU (agreement < 1e-6 confirms
# the GPU route). Distinct from run_gpu.sh (the cupy validation suite).
set -euo pipefail
cd "$(dirname "$0")"
source .venv/bin/activate
python3 -m cl10.altgpu
