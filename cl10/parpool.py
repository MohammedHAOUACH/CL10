"""cl10/parpool.py - process-pool parallel map for the experiments.

The engine picks the DEVICE; this picks the WORK shape.  When the work is
"many independent small matrix ops" (per-bond grids of 128x128 eigvals)
a single device is the wrong unit of parallelism: cupy has no batched
eigvals, so the 'batched' path is a serial loop of GPU round-trips that
measured ~50 ms per call on this host, versus ~8 ms for a numpy call.
The engine's own answer for a saturated card is the numpy fallback; the
answer for a busy-but-useful CPU is a process pool (fork), OpenBLAS
inside each worker.

Usage:

    from cl10 import parpool
    results = parpool.pmap(my_module_level_fn, items,
                            progress=lambda d, n: print(...))

`fn` must be a MODULE-LEVEL function (it is pickled to the workers).
CL10_NJOBS overrides the worker count (default: min(8, cpus//2)).
"""
import multiprocessing as mp
import os


def njobs_default():
    n = int(os.environ.get("CL10_NJOBS", "0"))
    if n > 0:
        return n
    return max(1, min(8, (os.cpu_count() or 1) // 2))


# The worker function, held in a module global: a closure defined inside
# pmap() cannot be pickled, but with fork (Linux default) every child
# inherits the parent's globals, and the wrapper below is a module-level
# (picklable-by-reference) function.
_FN = None


def _init(fn):
    """Pool initializer: install the worker function (also covers
    forkserver/spawn, where the parent's globals are NOT inherited)."""
    global _FN
    _FN = fn


def _wrap(t):
    i, it = t
    return i, _FN(it)


def pmap(fn, items, njobs=None, chunk=None, progress=None):
    """Map `fn` over `items` with a process pool, results in input order.

    `fn` must be a module-level function.  Runs sequentially for njobs == 1
    or tiny item lists (smoke tests).  `progress(done, n)` fires after each
    completed item (throttled by the pool's chunking - pass chunk=1 to
    report per item).
    """
    global _FN
    _FN = fn
    njobs = njobs or njobs_default()
    n = len(items)
    if njobs <= 1 or n <= 2:
        out = []
        for i, it in enumerate(items):
            out.append(fn(it))
            if progress:
                progress(i + 1, n)
        return out

    # Python 3.14+ defaults to forkserver (globals not inherited by the
    # worker); force fork so the parent's modules are visible, and pass fn
    # through the initializer as belt-and-braces.
    try:
        ctx = mp.get_context("fork")
    except ValueError:  # pragma: no cover - platforms without fork
        ctx = mp
    pool = ctx.Pool(njobs, initializer=_init, initargs=(fn,))
    try:
        chunk = chunk or max(1, n // (4 * njobs))
        out = [None] * n
        done = 0
        for i, r in pool.imap_unordered(_wrap, list(enumerate(items)),
                                        chunksize=chunk):
            out[i] = r
            done += 1
            if progress:
                progress(done, n)
    finally:
        pool.close()
        pool.join()
    return out
