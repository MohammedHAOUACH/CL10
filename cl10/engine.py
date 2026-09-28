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

The declared stack (``ga`` in ``pyproject.toml``) is ``kingdon`` (GPU physics
utilities), ``einops`` (the tensor reshapes), ``jax[cuda13]`` and
``cupy-cuda13x`` (the GPU backends).  ``einops`` is used for the einsum-style
resolves when present; the rest fall back to NumPy.
"""

import numpy as np

# ---------------------------------------------------------------------------
# backend selection:  CuPy  ->  JAX  ->  NumPy
# ---------------------------------------------------------------------------
_BACKEND = "numpy"
_E = np

try:  # the preferred GPU backend (CuPy)
    import cupy as _cp
    _E = _cp
    _BACKEND = "cupy"
except Exception:  # noqa: BLE001
    try:  # the second GPU backend (JAX / XLA)
        import jax
        import jax.numpy as _jnp
        jax.config.update("jax_platform_name", "gpu")
        _E = _jnp
        _BACKEND = "jax"
    except Exception:  # noqa: BLE001
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


def zeros(*shape, dtype=None):
    if dtype is None:
        dtype = complex
    return _E.zeros(*shape, dtype=dtype)


def kron(a, b):
    return _E.kron(a, b)


def matmul(a, b):
    return _E.matmul(a, b)


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
    """einops rearrange (used when einops is available, else a NumPy fallback)."""
    if _HAS_EINOPS:
        return einops.rearrange(asnumpy(x), pattern_in, pattern_out)
    raise RuntimeError("einops not installed (pip install einops)")


# ---------------------------------------------------------------------------
# the QCA / lattice primitives (the shared dynamics)
# ---------------------------------------------------------------------------

def laplacian_3d(u):
    """The 6-neighbour discrete Laplacian on a 3D lattice (the metric operator)."""
    up = (np.roll(u, -1, axis=1) + np.roll(u, 1, axis=1)
          + np.roll(u, -1, axis=2) + np.roll(u, 1, axis=2)
          + np.roll(u, -1, axis=3) + np.roll(u, 1, axis=3))
    return (up - 6 * u)


def qca_step(fields, unitaries):
    """One QCA step: apply the local unitaries (a product of 2x2 or 4x4 gates).

    ``fields`` is a 4-D array (x, y, z, internal); each gate in ``unitaries`` is
    applied on the internal (last) axis.  The product of local unitaries is
    exactly the QCA dynamics of the theory.
    """
    out = fields
    for U in unitaries:
        n = U.shape[0]
        out = np.einsum("...ij,jk->...ik", out.reshape(*out.shape[:-1], -1, n)
                        if out.shape[-1] == n else out, U)
    return out


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
