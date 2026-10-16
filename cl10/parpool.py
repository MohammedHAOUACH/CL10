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


def pmap(fn, items, njobs=None, chunk=None, progress=None):
    """Map `fn` over `items` with a process pool, results in input order.

    Runs sequentially for njobs == 1 or tiny item lists (smoke tests).
    `progress(done, n)` fires after each completed item (throttled by the
    pool's chunking - pass chunk=1 to report per item).
    """
    njobs = njobs or njobs_default()
    n = len(items)
    if njobs <= 1 or n <= 2:
        out = []
        for i, it in enumerate(items):
            out.append(fn(it))
            if progress:
                progress(i + 1, n)
        return out

    def wrap(t):
        i, it = t
        return i, fn(it)

    pool = mp.Pool(njobs)
    try:
        chunk = chunk or max(1, n // (4 * njobs))
        out = [None] * n
        done = 0
        for i, r in pool.imap_unordered(wrap, list(enumerate(items)),
                                        chunksize=chunk):
            out[i] = r
            done += 1
            if progress:
                progress(done, n)
    finally:
        pool.close()
        pool.join()
    return out
