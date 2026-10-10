"""CL10 publication figures: 39 figures, all computed from the construction.

One command regenerates everything:

    .venv/bin/python -m make_figures            # all 39, GPU if cupy is there
    .venv/bin/python -m make_figures --only 7 12 19   # a subset
    .venv/bin/python -m make_figures --cpu      # force the numpy backend
    .venv/bin/python -m make_figures --fast     # skip the two slow GPU figures
                                              # (18 young, 19 hydrogen/V62)

Modules group the figures by subject; each figure function is named
``fig<N>_<slug>`` and writes ``doc/figures/fig<N>_<slug>.png`` with the
number used in ``doc/publication.md``.  Fig. 18 (young double-slit) runs the
``simulation/young`` QCA; Fig. 19 (hydrogen) is produced by
``experiments/V62_hydrogen_production``, invoked here as a subprocess so it
re-pins its own GPU and writes its PNG + meta JSON itself.

The GPU id is pinned from ``simulation/young/config.yaml`` *before* any CuPy
import (mirrors ``simulation/*/__main__.py``); ``CL10_BACKEND=numpy`` (or
``--cpu``) forces the CPU path.
"""
import os
import subprocess
import sys

from . import algebra, couplings, dynamics, proton, residuals, scales

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (number, callable) in paper order
FIGURES = [
    (1, algebra.fig1_one_object),
    (2, scales.fig2_wimp_miracle),
    (3, algebra.fig3_512_512),
    (4, scales.fig4_qcd),
    (5, scales.fig5_baryon),
    (6, algebra.fig6_anomaly),
    (7, scales.fig7_ladder),
    (8, algebra.fig8_cpt_mirror),
    (9, scales.fig9_sm_head2head),
    (10, algebra.fig10_matter_content),
    (11, couplings.fig11_couplings_mz),
    (12, scales.fig12_ladder_sm),
    (13, scales.fig13_residual_reduction),
    (14, couplings.fig14_uv_2loop),
    (15, proton.fig15_proton_channels),
    (16, dynamics.fig16_qca_wave),
    (17, dynamics.fig17_gravity_3d),
    (18, dynamics.fig18_young),
    (20, couplings.fig20_desert),
    (21, couplings.fig21_wilson),
    (22, residuals.fig22_gauged_layer),
    (23, couplings.fig23_composite),
    (24, couplings.fig24_closed_loops),
    (25, couplings.fig25_equal_shift),
    (26, residuals.fig26_w4_closure),
    (27, residuals.fig27_w7_portal),
    (28, residuals.fig28_loop_selection),
    (29, residuals.fig29_w2_protection),
    (30, residuals.fig30_w3_forcing),
    (31, algebra.fig31_graphical_abstract),
    (32, residuals.fig32_w4_cone),
    (33, algebra.fig33_block_zone_nodes),
    (34, dynamics.fig34_bh_frames),
    (35, residuals.fig35_vertex_ground_state),
    (36, residuals.fig36_residual_attacks),
    (37, residuals.fig37_doubler_cones),
    (38, residuals.fig38_isotropy_search),
    (39, residuals.fig39_tangent_algebra),
]

YOUNG = 18
HYDROGEN = 19  # produced by experiments/V62_hydrogen_production


def setup_gpu(force_cpu=False):
    """Pin the GPU and pick the backend before anything imports CuPy."""
    if force_cpu:
        os.environ["CL10_BACKEND"] = "numpy"
        print("backend: numpy (--cpu)")
        return
    gpu = 1
    cfg = os.path.join(ROOT, "simulation", "young", "config.yaml")
    if os.path.exists(cfg):
        import yaml
        with open(cfg) as fh:
            gpu = yaml.safe_load(fh).get("gpu", 1)
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
        if free < 2048:
            print(f"GPU {gpu} has less than 2 GiB free -> falling back to numpy "
                  "(free the card or edit gpu: in simulation/*/config.yaml)")
            os.environ["CL10_BACKEND"] = "numpy"
    except Exception:
        pass


def run_v62():
    """Fig. 19: the production hydrogen experiment writes its own PNG + meta."""
    script = os.path.join(ROOT, "experiments", "V62_hydrogen_production",
                          "V62_hydrogen_production.py")
    rc = subprocess.call([sys.executable, script], cwd=ROOT)
    if rc != 0:
        raise SystemExit(f"V62 failed (rc={rc}): fig19_hydrogen.png NOT regenerated")


def main(only=None, cpu=False, fast=False):
    from .style import OUT
    setup_gpu(force_cpu=cpu)
    print("Generating CL10 publication figures ->", OUT)
    n = 0
    for num, fn in FIGURES:
        if only is not None and num not in only:
            continue
        if fast and num == YOUNG:
            continue
        fn()
        n += 1
    if (only is None or HYDROGEN in only) and not fast:
        run_v62()
        n += 1
    print(f"done: {n} figures -> {OUT}")
