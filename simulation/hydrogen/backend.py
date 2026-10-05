#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The backend selection of the hydrogen simulation.

Identical policy to ``simulation/young/backend.py``: **CuPy** (GPU) when
available, else **NumPy** (CPU).  The physics is identical on either; the
GPU id is pinned in ``__main__.py`` *before* this module is imported.
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
# the field precision (the simulation config may choose float64 for
# long-run robustness; the default float32 is the fast path)
FIELD = np.float32


def set_field_dtype(name):
    """Select the field precision: ``"float32"`` (default) or
    ``"float64"``.  Must be called **before** the QCA is built."""
    global FIELD
    FIELD = getattr(np, str(name).lower().replace("float", "float"))


def backend_name():
    """Which backend is active: ``cupy`` (GPU) or ``numpy`` (CPU)."""
    return _BACKEND


def gpu_id():
    """The visible GPU id.  With ``CUDA_VISIBLE_DEVICES=N`` set, the process
    sees it as device 0; this reports the *requested* id from the env."""
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


def sin(x):
    return _E.sin(x)


def clip(x, lo, hi):
    return _E.clip(x, lo, hi)


def log1p(x):
    return _E.log1p(x)


def where(cond, a, b):
    return _E.where(cond, a, b)
def minimum(a, b):
    return _E.minimum(a, b)


def maximum(a, b):
    return _E.maximum(a, b)


def zeros_like(x):
    return _E.zeros_like(x)


def full(shape, value, dtype=None):
    return _E.full(shape, value, dtype=dtype or _E.float32)


def asnumpy(x):
    """Move to the host (for PyVista / VisPy / printing / JSON)."""
    if isinstance(x, np.ndarray):
        return x
    return x.get() if hasattr(x, "get") else np.asarray(x)


def gpu_info():
    """A short one-line description of the active device (for the panel)."""
    if _BACKEND == "cupy":
        try:
            mem = _cp.cuda.Device().mem_info
            return (
                f"GPU cuda:{gpu_id()} "
                f"({mem[0] / 2**20:.0f}/{mem[1] / 2**20:.0f} MiB)"
            )
        except Exception:  # noqa: BLE001
            return f"GPU cuda:{gpu_id()}"
    return "CPU (numpy)"
