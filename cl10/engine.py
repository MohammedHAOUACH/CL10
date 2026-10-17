#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the unique computation engine.

ONE core for the whole theory.  It picks the backend at import time and exposes
a single array type ``Array`` plus the primitives every module and every
experiment in the project uses, so the whole theory runs on the GPU out of the
box and degrades to CPU cleanly.

**The configuration is one file.**  ``cl10/config.yaml`` (loaded by
``cl10/config.py``) carries the backend order, the GPU free-memory gate, the
probe sizes, the routing thresholds below, the process-pool shape, the physics
inputs and the shared tolerances.  Nothing here hard-codes them; edit the YAML
and the whole engine follows.  The environment overrides the file
(``CL10_BACKEND``, ``CL10_GPU_MIN_FREE_MIB``, ``CL10_ANGLE``, ``CL10_NJOBS``,
``CL10_SMOKE``).

Three rules the primitives enforce:

* **GPU first, where the GPU wins.**  Every array primitive (``roll``,
  ``einsum``, ``matmul``, ``kron``, ``where``, the reductions ``summ`` /
  ``mean`` / ``sum_axis`` / ``field_norm``) works **where the data lives**: a
  device array is never pulled to the host by a host call (calling
  ``np.roll``/``np.einsum`` on a device array silently copies it, which is why
  the real-space automaton used to run on the CPU whatever the backend was),
  and a host array is never pushed to the device.  The linear-algebra helpers
  (``eigvals``, ``eigvalsh``, ``eig``, ``det``, ``inv``, ``norm``,
  ``matrix_rank``) route by ``routing.*``: above the configured size they run on
  the device, below it they go to the host because that is measured to be
  faster.  On this host (2x RTX 3090, 12 cores) a CuPy call on a 128x128 matrix
  costs ~50 ms against ~8 ms for NumPy, and CuPy has no batched
  eigendecomposition — so small matrices are host work and their parallelism
  comes from ``cl10/parpool.py``.  ``route_report()`` prints the policy.

* **Host transfers are written down.**  Only ``asnumpy`` and the routed
  eigensolvers move data to the host; the reductions (``field_norm``, ``summ``,
  ``mean``, ``absmax``) return a host scalar without moving the field.

* **Double precision everywhere.**  JAX defaults to float32/complex64; the
  tolerances of this repository are 1e-12 and below, so the JAX branch enables
  x64 before the first array is built.

The declared stack (``dependencies`` in ``pyproject.toml``) is ``kingdon`` (GPU
physics utilities), ``einops`` (the tensor reshapes), ``jax[cuda13]``,
``cupy-cuda13x`` (the GPU backends) and ``pyyaml`` (this configuration).
"""

import os

import numpy as np

from . import config as _config

_CFG = _config.CFG

# ---------------------------------------------------------------------------
# backend selection: the config sets the order, the probe decides the winner.
# An explicit force (CL10_BACKEND=cupy) skips the free-memory gate and never
# falls back: a GPU run that cannot get the card fails loudly instead of
# reporting a GPU verdict that was computed on the CPU.
# ---------------------------------------------------------------------------
_BACKEND = "numpy"
_E = np

_FORCE = (_CFG.get("backend.force") or "").strip().lower()


def _try_cupy():
    """CuPy.  In auto mode the card is gated on free memory and probed through
    cuBLAS + cuSOLVER.  A **forced** backend (``CL10_BACKEND=cupy``) skips the
    gate: the GPU runner has already reserved the aggregate budget, and a child
    that quietly fell back to NumPy would silently invalidate the GPU run."""
    import cupy as cp
    if not _FORCE:
        free, _total = cp.cuda.runtime.memGetInfo()
        need = int(_CFG.backend.gpu_min_free_mib) << 20
        if free < need:
            raise RuntimeError(f"GPU busy: {free >> 20} MiB free (< {need >> 20} MiB)")
    # a saturated card lets tiny probe ops through but queues real work
    # indefinitely, so the probe allocates and touches a real buffer
    probe = cp.empty(int(_CFG.backend.probe.elements), dtype=cp.float64)
    cp.add(probe[:1], probe[1:2])
    del probe
    n = int(_CFG.backend.probe.matrix_dim)
    mat = cp.eye(n, dtype=cp.complex128)
    cp.matmul(mat, mat)
    cp.linalg.inv(mat)
    del mat
    return cp, "cupy"


def _try_jax():
    """JAX / XLA, x64 before the first array is built."""
    os.environ.setdefault("XLA_PYTHON_CLIENT_MEM_FRACTION",
                          str(_CFG.backend.jax.mem_fraction))
    if _CFG.backend.jax.enable_x64:
        os.environ.setdefault("JAX_ENABLE_X64", "1")
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", bool(_CFG.backend.jax.enable_x64))
    jax.config.update("jax_platform_name", "gpu")
    n = int(_CFG.backend.probe.matrix_dim)
    mat = jnp.eye(n, dtype=jnp.complex128)
    (mat @ mat).block_until_ready()
    del mat
    return jnp, "jax"


def _try_numpy():
    return np, "numpy"


_TRIES = {"cupy": _try_cupy, "jax": _try_jax, "numpy": _try_numpy}

_gpu_unusable = False
for _name in list(_CFG.backend.priority):
    if _FORCE and _FORCE != _name and not (_FORCE == "cpu" and _name == "numpy"):
        continue                       # the config/env pinned another backend
    if _name == "jax" and _gpu_unusable:
        continue                       # the card is CUDA-busy, not absent: JAX would
                                       # OOM the same way (and spams XLA retries)
    try:
        _E, _BACKEND = _TRIES[_name]()
        break
    except Exception:
        if _FORCE == _name:
            raise                      # forced: a silent fallback would lie about the run
        if _name == "cupy":
            _gpu_unusable = True
        continue

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

def _device_work(operands):
    """True when a call must run on the device: the active backend is a GPU one
    **and** at least one operand is device-resident.  Host operands stay on the
    host — pulling them to the device and back is the accident this prevents."""
    if _BACKEND == "numpy":
        return False
    for a in operands:
        if isinstance(a, (int, float, complex, np.ndarray)):
            continue
        return True
    return False


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
    return _E.kron(a, b) if _device_work((a, b)) else np.kron(a, b)


def matmul(a, b):
    return _E.matmul(a, b) if _device_work((a, b)) else np.matmul(a, b)


def einsum(subscripts, *operands):
    """``einsum`` where the data lives: a NumPy call on a device array copies it
    to the host (which is why the QCA layers used to run on the CPU whatever the
    backend was), and a backend call on a host array copies it the other way."""
    return _E.einsum(subscripts, *operands) if _device_work(operands) \
        else np.einsum(subscripts, *operands)


def roll(a, shift, axis):
    """``roll`` where the data lives (the bond layer's neighbour shift)."""
    return _E.roll(a, shift, axis=axis) if _device_work((a,)) \
        else np.roll(a, shift, axis=axis)


def real(x):
    return x.real if _BACKEND == "numpy" or not _device_work((x,)) else _E.real(x)


def imag(x):
    return x.imag if _BACKEND == "numpy" or not _device_work((x,)) else _E.imag(x)


def absmax(x):
    return float(np.max(np.abs(x))) if not _device_work((x,)) \
        else float(_E.max(_E.abs(x)).item())


def asnumpy(x):
    """Move to host (for printing / JSON)."""
    if _BACKEND == "cupy":
        import cupy as cp
        return cp.asnumpy(x)
    return np.asarray(x)


def to_device(x):
    """Move a host array to the active backend."""
    if _BACKEND == "numpy":
        return np.asarray(x)
    return _E.asarray(np.asarray(x))


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
# the routing policy: GPU where the GPU wins, host where the host wins
# ---------------------------------------------------------------------------

def _last_dim(a):
    shape = getattr(a, "shape", ())
    return int(shape[-1]) if len(shape) >= 2 else 0


def _on_device(op, a):
    """Decide one linear-algebra call: the array must be device-resident, and
    the matrix must be at least the configured size for that class of op."""
    if not _device_work((a,)):
        return False
    n = _last_dim(a)
    if op in ("eig", "eigvals", "eigvalsh", "svd", "matrix_rank"):
        return n >= int(_CFG.routing.eig_gpu_min_dim)
    return n >= int(_CFG.routing.device_min_dim)


def eigvals(a):
    """The eigenvalues, on the device above ``routing.eig_gpu_min_dim``."""
    if _on_device("eigvals", a):
        return _E.linalg.eigvals(a)
    return np.linalg.eigvals(asnumpy(a))


def eigvalsh(a):
    """The eigenvalues of a Hermitian matrix (same routing as ``eigvals``)."""
    if _on_device("eigvalsh", a):
        return _E.linalg.eigvalsh(a)
    return np.linalg.eigvalsh(asnumpy(a))


def eig(a):
    """``eig`` (values + vectors); the vectors of a device matrix are returned
    on the host when the call is routed to the host."""
    if _on_device("eig", a):
        return _E.linalg.eig(a)
    return np.linalg.eig(asnumpy(a))


def eigh(a):
    if _on_device("eigvalsh", a):
        return _E.linalg.eigh(a)
    return np.linalg.eigh(asnumpy(a))


def det(a):
    if _on_device("det", a):
        return _E.linalg.det(a)
    return np.linalg.det(asnumpy(a))


def inv(a):
    if _on_device("inv", a):
        return _E.linalg.inv(a)
    return np.linalg.inv(asnumpy(a))


def norm(a, *args, **kwargs):
    """``linalg.norm``: reduced on the device when the array is on it, on the
    host otherwise (the same value either way)."""
    if _device_work((a,)):
        return float(_E.linalg.norm(a, *args, **kwargs).item())
    return float(np.linalg.norm(np.asarray(a), *args, **kwargs))


def field_norm(a):
    """The L2 (Frobenius) norm of a field, reduced where the field lives: a
    device field never comes back to the host for its own norm."""
    if _device_work((a,)):
        return float(_E.sqrt(_E.sum(_E.abs(a) ** 2)).item())
    return float(np.sqrt(np.sum(np.abs(np.asarray(a)) ** 2)))


def summ(a):
    """The total sum, reduced where the data lives (no host copy of a field)."""
    return float(_E.sum(a).item()) if _device_work((a,)) else float(np.sum(a))


def mean(a):
    """The mean, reduced where the data lives."""
    return float(_E.mean(a).item()) if _device_work((a,)) else float(np.mean(a))


def where(cond, a, b):
    """``where`` where the data lives (the mask follows the operands)."""
    if _device_work((cond, a, b)):
        return _E.where(cond, a, b)
    return np.where(cond, a, b)



def modulus(a):
    """The modulus, staying where the data lives."""
    return _E.abs(a) if _device_work((a,)) else np.abs(a)


def sum_axis(a, axis=None):
    """A partial sum along ``axis``, staying where the data lives."""
    return _E.sum(a, axis=axis) if _device_work((a,)) else np.sum(a, axis=axis)


def matrix_rank(a, tol=None):
    if _on_device("matrix_rank", a):
        return int(_E.linalg.matrix_rank(a, tol=tol))
    return int(np.linalg.matrix_rank(asnumpy(a), tol=tol))


def qr(a):
    if _on_device("qr", a):
        return _E.linalg.qr(a)
    return np.linalg.qr(asnumpy(a))


def route_report():
    """Where each class of work runs under the current configuration."""
    gpu = _BACKEND != "numpy"
    eig_min = int(_CFG.routing.eig_gpu_min_dim)
    dev_min = int(_CFG.routing.device_min_dim)
    return {
        "backend": _BACKEND,
        "lattice_fields_roll_einsum_matmul": _BACKEND,
        "square_algebra_min_dim": dev_min if gpu else 0,
        "eigendecomposition_min_dim": eig_min if gpu else 0,
        "eigendecomposition_below_threshold": "numpy (measured faster; no batched device eig)",
        "pool_njobs": int(_CFG.pool.njobs_effective),
    }


# ---------------------------------------------------------------------------
# the QCA / lattice primitives (the shared dynamics)
# ---------------------------------------------------------------------------

def laplacian_3d(u):
    """The 6-neighbour discrete Laplacian on the **last 3** (spatial) axes of
    ``u`` (works for a ``(L,L,L)`` field or a ``(mu,nu,L,L,L)`` metric field).
    Device-aware: it goes through ``roll``, so a device field never leaves the
    device."""
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
        "gpu_min_free_mib": int(_CFG.backend.gpu_min_free_mib),
        "config": _CFG.get("config_path"),
        **route_report(),
    }
