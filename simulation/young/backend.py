#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The single computation core of the simulation.

Mirrors ``cl10/engine.py``: it picks **CuPy** (GPU) -> **NumPy** (CPU) at
import.  The **GPU id** is pinned by the runner *before* this module is
imported (``CUDA_VISIBLE_DEVICES`` from ``config.yaml``), so "use GPU 1" is
honoured here.  Everything the sim touches (the 3-D FDTD, the sponge, the
collapse) goes through this one module, so the whole thing runs on the GPU
out of the box.
"""

import os

import numpy as np

# ---------------------------------------------------------------------------
# backend selection:  CuPy  ->  NumPy   (the physics is identical on either)
# ---------------------------------------------------------------------------
_BACKEND = "numpy"
_E = np

try:  # the GPU backend (CuPy)
    import cupy as _cp  # noqa: F401
    _E = _cp
    _BACKEND = "cupy"
except Exception:  # noqa: BLE001
    _E = np
    _BACKEND = "numpy"

# the array type (cp.ndarray | np.ndarray)
Array = _E.ndarray


def backend_name():
    """Which backend is active: ``cupy`` (GPU) or ``numpy`` (CPU)."""
    return _BACKEND


def gpu_id():
    """The visible GPU id.  With ``CUDA_VISIBLE_DEVICES=N`` set, the process
    sees exactly one card and it is reported as device 0; we return the id
    that was pinned so the panel can show 'GPU N'."""
    return os.environ.get("CUDA_VISIBLE_DEVICES", "?")


def zeros(shape, dtype=None):
    """A zero array on the active backend (``shape`` a tuple)."""
    return _E.zeros(shape, dtype=dtype or _E.float32)


def ones(shape, dtype=None):
    return _E.ones(shape, dtype=dtype or _E.float32)


def arange(*args, **kwargs):
    return _E.arange(*args, **kwargs)


def exp(x):
    return _E.exp(x)


def cos(x):
    return _E.cos(x)


def clip(x, lo, hi):
    return _E.clip(x, lo, hi)


def minimum(a, b):
    return _E.minimum(a, b)


def meshgrid(*idx, indexing="ij"):
    return _E.meshgrid(*idx, indexing=indexing)


def zeros_like(x):
    return _E.zeros_like(x)


def asnumpy(x):
    """Move to the host (for PyVista / printing / JSON)."""
    if _BACKEND == "cupy":
        return _cp.asnumpy(x)
    return np.asarray(x)


def gpu_info():
    """A short one-line description of the active device (for the panel)."""
    if _BACKEND == "cupy":
        try:
            import cupy
            name = cupy.cuda.Device(0).name
            mem = cupy.cuda.Device(0).mem_info
            return f"GPU:{gpu_id()} {name} ({mem[0]/1e6:.0f} MiB free)"
        except Exception:  # noqa: BLE001
            return f"GPU:{gpu_id()} (cupy)"
    return "CPU (numpy)"
