"""CL10 publication figures: 25 figures, all computed from the construction.

One command regenerates everything:

    .venv/bin/python -m make_figures            # all 25, GPU if cupy is there
    .venv/bin/python -m make_figures --only 7 12 19   # a subset
    .venv/bin/python -m make_figures --cpu      # force the numpy backend
    .venv/bin/python -m make_figures --fast     # skip the two slow GPU figures
                                              # (1 young, 2 hydrogen)

Modules group the figures by subject; each figure function is named
``fig<N>_<slug>`` and writes ``doc/figures/fig<N>_<slug>.png`` with the
number used in ``doc/final_publication.md``.  Fig. 1 (young double-slit) runs
the ``simulation/young`` QCA; Fig. 2 (hydrogen) is produced by
``experiments/phenomenology/V62_hydrogen_production``, invoked here as a
subprocess: it re-pins its own GPU and writes its PNG + meta JSON itself.

The GPU id is pinned from ``simulation/young/config.yaml`` *before* any CuPy
import (mirrors ``simulation/*/__main__.py``); ``CL10_BACKEND=numpy`` (or
``--cpu``) forces the CPU path.
"""
import os
import subprocess
import sys

from . import algebra, couplings, dynamics, residuals, scales
from .style import exp

from cl10 import config as _config      # the gate only: no backend import here

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (number, callable) in paper order
FIGURES = [
    (1, dynamics.fig1_young),
    (3, algebra.fig3_one_object),
    (4, couplings.fig4_uv_2loop),
    (5, algebra.fig5_brillouin_doublers),
    (6, couplings.fig6_wilson),
    (7, residuals.fig7_doubler_cones),
    (8, residuals.fig8_isotropy_search),
    (9, residuals.fig9_tangent_algebra),
    (10, residuals.fig10_gauged_layer),
    (11, residuals.fig11_vertex_ground_state),
    (12, dynamics.fig12_gravity_3d),
    (13, couplings.fig13_closed_loops),
    (14, residuals.fig14_loop_selection),
    (15, residuals.fig15_w7_portal),
    (16, residuals.fig16_w2_protection),
    (17, scales.fig17_ladder),
    (18, scales.fig18_residual_reduction),
    (19, couplings.fig19_couplings_mz),
    (20, couplings.fig20_desert),
    (21, couplings.fig21_equal_shift),
    (22, couplings.fig22_composite),
    (23, residuals.fig23_w4_cone),
    (24, residuals.fig24_w3_forcing),
    (25, residuals.fig25_residual_attacks),
]

YOUNG = 1
HYDROGEN = 2  # produced by experiments/phenomenology/V62_hydrogen_production


def setup_gpu(force_cpu=False):
    """Pin the GPU and pick the backend before anything imports CuPy."""
    if force_cpu:
        os.environ["CL10_BACKEND"] = "numpy"
        print("backend: numpy (--cpu)")
        return
    cfg = os.path.join(ROOT, "simulation", "young", "config.yaml")
    import yaml
    with open(cfg) as fh:
        gpu = yaml.safe_load(fh)["gpu"]      # the panel owns the card id
    os.environ.setdefault("CL10_BACKEND", "cupy")
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    try:
        import cupy  # noqa: F401
    except Exception:
        print("cupy unavailable -> backend: numpy")
        os.environ["CL10_BACKEND"] = "numpy"
        return
    print(f"backend: cupy | GPU pinned: {gpu}")
    # room check via nvidia-smi: creating a CUDA context can itself OOM
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.free",
             "--format=csv,noheader,nounits", "-i", str(gpu)],
            capture_output=True, text=True, timeout=10)
        free = int(out.stdout.strip())
        print(f"GPU {gpu} free: {free} MiB")
        gate = int(_config.CFG.backend.gpu_min_free_mib)
        if free < gate:
            print(f"GPU {gpu} has less than {gate} MiB free -> falling back to "
                  "numpy (free the card, raise backend.gpu_min_free_mib in "
                  "cl10/config.yaml, or edit gpu: in simulation/*/config.yaml)")
            os.environ["CL10_BACKEND"] = "numpy"
    except Exception:
        pass


def run_v62():
    """Fig. 2: the production hydrogen experiment writes its own PNG + meta."""
    script = exp("V62_hydrogen_production", "V62_hydrogen_production.py")
    rc = subprocess.call([sys.executable, script], cwd=ROOT)
    if rc != 0:
        raise SystemExit(f"V62 failed (rc={rc}): fig2_hydrogen.png NOT regenerated")


def main(only=None, cpu=False, fast=False):
    from .style import OUT
    from tqdm import tqdm
    setup_gpu(force_cpu=cpu)
    print("Generating CL10 publication figures ->", OUT)
    todo = [(num, fn) for num, fn in FIGURES
            if (only is None or num in only) and not (fast and num == YOUNG)]
    hydrogen = (only is None or HYDROGEN in only) and not fast
    # the bar goes to stderr (run_figures.sh keeps it on the terminal): the
    # figures take from a second to minutes, and Fig. 2 re-runs the whole V62
    # production experiment, so a silent terminal there reads as a hang.
    bar = tqdm(total=len(todo) + int(hydrogen), desc="figures", unit="fig",
               colour="green") if sys.stderr.isatty() else None
    n = 0
    for num, fn in todo:
        if bar:
            bar.set_postfix_str(f"fig {num}")
        fn()
        n += 1
        if bar:
            bar.update(1)
    if hydrogen:
        if bar:
            bar.set_postfix_str(f"fig {HYDROGEN} (the V62 hydrogen experiment)")
        run_v62()
        n += 1
        if bar:
            bar.update(1)
    if bar:
        bar.close()
    print(f"done: {n} figures -> {OUT}")
