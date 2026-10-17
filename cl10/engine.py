#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the unique computation engine.

ONE core for the whole theory.  It picks the backend at import time:

    CuPy  (GPU, the preferred one)  ->  JAX (XLA / GPU)  ->  NumPy (CPU)

and exposes a single array type ``Array`` plus the primitives every module and
every experiment in the project uses.  All the physics code goes through this
one file, so the whole theory runs on the GPU (the RTX 3090) out of the box and
degrades to CPU cleanly.

Two rules the primitives enforce:

* **Stay on the device.**  Every array primitive (``roll``, ``einsum``,
  ``matmul``, ``kron``, ``eye``, ``zeros``) dispatches through the active
  backend.  Calling ``np.roll``/``np.einsum`` on a CuPy or JAX array silently
  copies it to the host, which is why the real-space automaton used to run on
  the CPU whatever the backend was.  Host transfers happen only where they are
  written down: ``asnumpy`` and the momentum-space eigensolvers (CuPy has no
  batched eigendecomposition, so the Bloch matrices are diagonalised on the
  host by design).
* **Double precision everywhere.**  JAX defaults to float32/complex64; the
  tolerances of this repository are 1e-12 and below, so the JAX branch enables
  x64 before the first array is built.

The declared stack (``ga`` in ``pyproject.toml``) is ``kingdon`` (GPU physics
utilities), ``einops`` (the tensor reshapes), ``jax[cuda13]`` and
``cupy-cuda13x`` (the GPU backends).  ``einops`` is used for the einsum-style
resolves when present; the rest fall back to NumPy.
"""

import os

import numpy as np

# ---------------------------------------------------------------------------
# backend selection:  CuPy  ->  JAX  ->  NumPy
# ---------------------------------------------------------------------------
_BACKEND = "numpy"
_E = np

# CL10_BACKEND=numpy forces the CPU core (the small algebra checks - 3x3, the
# mixings - do not need the GPU and must run even when the GPU is busy with
# another process).  Default: no effect, the GPU stays the preferred backend.
_FORCE_CPU = os.environ.get("CL10_BACKEND", "").strip().lower() in ("numpy", "cpu")

_CUPY_PRESENT = False
try:  # the preferred GPU backend (CuPy)
    if _FORCE_CPU:
        raise ImportError("CL10_BACKEND=numpy (the CPU override)")
    import cupy as _cp
    _CUPY_PRESENT = True
    # free-memory gate: a saturated card (another tenant at work) lets
    # tiny probe ops slip through but queues real work indefinitely -
    # require >= 2 GiB free before committing to the GPU.  The gate is
    # tunable via CL10_GPU_MIN_FREE_MIB for small-batch work on a card
    # that is full but idle (the CL10 batches are tens of MB).
    _free, _total = _cp.cuda.runtime.memGetInfo()
    _min_free = int(os.environ.get("CL10_GPU_MIN_FREE_MIB", "2048")) << 20
    if _free < _min_free:
        raise RuntimeError(f"GPU busy: {_free >> 20} MiB free")
    # probe: cuBLAS + cuSOLVER handles too - a full card must fall
    # back to numpy, not die later on CUBLAS_STATUS_ALLOC_FAILED.
    _probe = _cp.empty(1 << 22, dtype=_cp.float64)
    _cp.add(_probe[:1], _probe[1:2])
    del _probe
    _pa = _cp.eye(64, dtype=_cp.complex128)
    _cp.matmul(_pa, _pa)
    _cp.linalg.inv(_pa)
    del _pa
    _E = _cp
    _BACKEND = "cupy"
except Exception:  # noqa: BLE001
    # cupy importable but the probe failed => the card is CUDA-busy,
    # not absent: JAX would OOM the same way (and spams XLA retries).
    if not _FORCE_CPU and not _CUPY_PRESENT:
        try:  # the second GPU backend (JAX / XLA)
            os.environ.setdefault("XLA_PYTHON_CLIENT_MEM_FRACTION", "0.5")
            # the repository works at 1e-12 and below: JAX's float32/complex64
            # default would silently break every tolerance, so x64 is enabled
            # before the first array is built.
            os.environ.setdefault("JAX_ENABLE_X64", "1")
            import jax
            import jax.numpy as _jnp
            jax.config.update("jax_enable_x64", True)
            jax.config.update("jax_platform_name", "gpu")
            _pj = _jnp.eye(64, dtype=_jnp.complex128)
            (_pj @ _pj).block_until_ready()
            del _pj
            _E = _jnp
            _BACKEND = "jax"
        except Exception:  # noqa: BLE001
            _E = np
            _BACKEND = "numpy"
    else:
        _E = np
        _BACKEND = "numpy"

# the array type (cp.ndarray | jax.Array | np.ndarray)
Array = _E.ndarray if hasattr(_E, "ndarray") else _E.array

# optional: einops (the tensor reshapes) — used when importable
try:
    import einops
    _HAS_EINOPS = True
except Exception:  # noqa: BLE001
    einops = None
    _HAS_EINOPS = False

# optional: kingdon (the GPU physics utilities) — used when importable
try:
    import kingdon  # noqa: F401
    _HAS_KINGDON = True
except Exception:  # noqa: BLE001
    kingdon = None
    _HAS_KINGDON = False

# ---------------------------------------------------------------------------
# primitives (all backend-aware)
# ---------------------------------------------------------------------------

def Array_(*args, **kwargs):
    """Construct an array on the active backend."""
    return _E.array(*args, **kwargs)


def eye(n, dtype=None):
    if dtype is None:
        dtype = complex
    return _E.eye(n, dtype=dtype)


def zeros(shape, dtype=None):
    """A zero array (shape is a tuple or int, e.g. ``((L,L,L))`` or ``4``)."""
    if isinstance(shape, int):
        shape = (shape,)
    if dtype is None:
        dtype = complex
    return _E.zeros(shape, dtype=dtype)


def kron(a, b):
    return _E.kron(a, b)


def matmul(a, b):
    return _E.matmul(a, b)


def einsum(subscripts, *operands):
    """``einsum`` on the active backend (a NumPy call on a device array copies
    it to the host, so the QCA layers must go through here)."""
    return _E.einsum(subscripts, *operands)


def roll(a, shift, axis):
    """``roll`` on the active backend (the bond layer's neighbour shift)."""
    return _E.roll(a, shift, axis=axis)


def real(x):
    return _E.real(x) if _BACKEND != "numpy" else x.real


def imag(x):
    return _E.imag(x) if _BACKEND != "numpy" else x.imag


def absmax(x):
    return float(_E.max(_E.abs(x)).item())


def asnumpy(x):
    """Move to host (for printing / JSON)."""
    if _BACKEND == "cupy":
        return _cp.asnumpy(x)
    if _BACKEND == "jax":
        return np.asarray(x)
    return np.asarray(x)


def einops_rearrange(x, pattern_in, pattern_out):
    """einops rearrange, on the active backend when einops supports it."""
    if not _HAS_EINOPS:
        raise RuntimeError("einops not installed (pip install einops)")
    try:
        return einops.rearrange(x, pattern_in, pattern_out)
    except (TypeError, ValueError):
        # einops has no backend for this array type: the reshape is a pure
        # relabelling, so it is done on the host copy.
        return einops.rearrange(asnumpy(x), pattern_in, pattern_out)


# ---------------------------------------------------------------------------
# the QCA / lattice primitives (the shared dynamics)
# ---------------------------------------------------------------------------

def laplacian_3d(u):
    """The 6-neighbour discrete Laplacian on the **last 3** (spatial) axes of
    ``u`` (works for a ``(L,L,L)`` field or a ``(mu,nu,L,L,L)`` metric field)."""
    b = u.ndim - 3
    up = (roll(u, -1, b) + roll(u, 1, b)
          + roll(u, -1, b + 1) + roll(u, 1, b + 1)
          + roll(u, -1, b + 2) + roll(u, 1, b + 2))
    return (up - 6 * u)


def spectral_laplacian(k, L):
    """The spectral (momentum-space) Laplacian eigenvalue on a torus of size L."""
    return 4.0 * (np.sin(k[0] * np.pi / L) ** 2
                  + np.sin(k[1] * np.pi / L) ** 2
                  + np.sin(k[2] * np.pi / L) ** 2)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def backend_info():
    return {
        "backend": _BACKEND,
        "einops": _HAS_EINOPS,
        "kingdon": _HAS_KINGDON,
        "array": Array.__name__,
    }
