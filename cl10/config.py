#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — the single configuration source.

One YAML file (`cl10/config.yaml`) carries every number the code used to hold as
a literal: the backend policy, the GPU gate, the routing thresholds that decide
where a task runs, the process-pool shape, the physics and Standard Model
inputs, the angle branches, the desert thresholds and the shared tolerance
bands.  Change it there and the whole engine picks it up — no file-by-file edit.

Precedence, highest first:

    1. the environment   CL10_BACKEND, CL10_GPU_MIN_FREE_MIB, CL10_ANGLE,
                         CL10_NJOBS, CL10_SMOKE
    2. cl10/config.yaml
    3. the built-in fallbacks below - which cover only how the engine runs and
       the root's own structure; the measured inputs exist in the file alone,
       so a missing one fails loudly instead of running on a stale number

Read it through the singleton `CFG`:

    from cl10 import config
    config.CFG.backend.gpu_min_free_mib      # 2048
    config.CFG.get("routing.eig_gpu_min_dim")
    config.CFG.angle.g10                     # the active branch value
    config.describe()                        # the effective values, flat

`config.reload()` re-reads the file (used by the tests and by any tool that
edits the configuration at runtime); `CFG` keeps its identity, so references held
by other modules stay valid.
"""

import copy
import functools
import math
import os
import pathlib

import yaml

PATH = pathlib.Path(__file__).with_name("config.yaml")

# the built-in fallbacks: how the engine runs (backend policy, the GPU gate, the
# routing thresholds, the pool shape, the tolerance bands) and the structure of
# the root object (the Cl(10) shape, the clock branches, the cell rungs). No
# measured input is duplicated here: a second copy of `m_planck_gev` or of a
# neutrino mass is exactly how a value drifts out of sync with the file that is
# supposed to own it, silently. With no config.yaml the engine still builds the
# algebra; the phenomenology modules fail on the missing key and say which one.
_DEFAULTS = {
    "backend": {
        "priority": ["cupy", "jax", "numpy"],
        "force": None,
        "gpu_min_free_mib": 2048,
        "probe": {"elements": 1 << 22, "matrix_dim": 64},
        "jax": {"mem_fraction": 0.5, "enable_x64": True},
    },
    "routing": {
        "device_min_dim": 64,
        "eig_gpu_min_dim": 512,
    },
    "pool": {"njobs": None, "njobs_cap": 8},
    "cell": {"rung_A": 3, "rung_B": 10},
    "physics": {
        "n_gen": 10,
        "n_space": 3,
        "n_subsites": 8,
        "n_internal": 32,
        "n_sector": 16,
        "loop_factor": None,
    },
    "angle": {
        "default_branch": "closed",
        "branches": {"implemented": 0.2, "closed": math.pi / 4, "half": math.pi / 2},
        "g0_from_n_space": True,
    },
    "tolerances": {"strict": 1e-12, "standard": 1e-9, "loose": 1e-6, "coarse": 1e-4},
}

# env variable -> (dotted config path, coercion)
_ENV_MAP = {
    "CL10_BACKEND": ("backend.force", str),
    "CL10_GPU_MIN_FREE_MIB": ("backend.gpu_min_free_mib", int),
    "CL10_NJOBS": ("pool.njobs", int),
    "CL10_ANGLE": ("angle.branch", str),
    "CL10_SMOKE": ("smoke", lambda v: str(v).strip() in ("1", "true", "yes")),
}


class Node:
    """A dict with attribute access (nested nodes stay nodes)."""

    __slots__ = ("_d",)

    def __init__(self, d):
        object.__setattr__(self, "_d", d)

    def __getattr__(self, key):
        try:
            value = self._d[key]
        except KeyError:
            raise AttributeError(f"no config key {key!r}") from None
        return Node(value) if isinstance(value, dict) else value

    def __getitem__(self, key):
        value = self._d[key]
        return Node(value) if isinstance(value, dict) else value

    def get(self, key, default=None):
        """`get("routing.eig_gpu_min_dim")` — dotted lookup, default if absent."""
        node = self
        for part in key.split("."):
            if not isinstance(node._d, dict) or part not in node._d:
                return default
            node = Node(node._d[part]) if isinstance(node._d[part], dict) else node._d[part]
        return node

    def keys(self):
        return self._d.keys()

    def as_dict(self):
        return copy.deepcopy(self._d)

    def __repr__(self):
        return f"Node({self._d!r})"


def _deep_merge(base, override):
    out = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def _set_dotted(tree, dotted, value):
    parts = dotted.split(".")
    node = tree
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    node[parts[-1]] = value


def _read_file():
    try:
        return yaml.safe_load(PATH.read_text()) or {}
    except FileNotFoundError:  # a checkout without the file still runs
        return {}


def _build():
    tree = _deep_merge(_DEFAULTS, _read_file())
    for var, (dotted, cast) in _ENV_MAP.items():
        raw = os.environ.get(var)
        if raw is None or str(raw).strip() == "":
            continue
        _set_dotted(tree, dotted, cast(str(raw).strip()))
    return tree


# the branch names the environment has always accepted; resolved once here so
# every reader (coupling, supercell, the experiments) sees one canonical name
_ANGLE_ALIASES = {
    "impl": "implemented", "implemented": "implemented",
    "1/5": "implemented", "0.2": "implemented",
    "closed": "closed", "pi/4": "closed", "quantized": "closed",
    "half": "half", "pi/2": "half",
}


def enforce_angle():
    """Fail loudly when the environment pins a clock this process did not
    have at config-import time.

    The clock is resolved exactly once, in ``_derive`` above. A script that
    sets ``CL10_ANGLE`` afterwards - typically a module imported after config,
    whose own preamble pins it - keeps measuring on the frozen clock and never
    finds out. Every module that binds the clock (``coupling``, ``supercell``)
    calls this at import, so the mismatch raises instead of shipping an
    artifact measured on a clock its caller did not ask for.
    """
    raw = os.environ.get("CL10_ANGLE")
    if raw is None or not str(raw).strip():
        return
    want = str(raw).strip().lower()
    want = _ANGLE_ALIASES.get(want, want)
    if want != CFG.angle.branch:
        raise RuntimeError(
            "CL10_ANGLE=%r pins the %r clock, but cl10.config froze the %r "
            "clock at its own import: set CL10_ANGLE before importing any cl10 "
            "module (the V118/V186/V187/V207 preambles do)."
            % (raw, want, CFG.angle.branch))


def _derive(tree):
    """The values that follow from the file instead of being written twice."""
    physics = tree["physics"]
    n_internal = int(physics["n_internal"])
    physics["n_qubits"] = n_internal.bit_length() - 1     # 32 -> 5
    angle = tree["angle"]
    branch = str(angle.get("branch") or angle["default_branch"]).strip().lower()
    branch = _ANGLE_ALIASES.get(branch, branch)
    branches = angle["branches"]
    if branch not in branches:
        raise ValueError(
            f"CL10_ANGLE / angle.default_branch = {branch!r}: "
            f"not one of {sorted(branches)}")
    angle["branch"] = branch
    angle["g10"] = float(branches[branch])
    n_space = int(tree["physics"]["n_space"])
    if angle.get("g0_from_n_space", True):
        angle["g0"] = 1.0 / (2.0 ** n_space)
    loop = tree["physics"].get("loop_factor")
    tree["physics"]["loop"] = float(loop) if loop else 1.0 / (16.0 * math.pi)
    phys = tree["physics"]
    if "v_over_mplanck_observed" not in phys:
        # the observed electroweak/Planck ratio is not a third input: it is the
        # quotient of the two measured constants above, so it is derived here
        # rather than written to the file at a lossier 5 digits
        phys["v_over_mplanck_observed"] = (float(phys["v_electroweak_gev"])
                                           / float(phys["m_planck_gev"]))
    if "hierarchy_coefficient_observed" not in phys:
        # the observed coefficient is not a fourth measured number either: it
        # is that quotient measured in units of the ladder factor, so the file
        # holds the two constants and the loader holds every quotient of them
        # (V52/V53/V108/V203 all compute this same value)
        phys["hierarchy_coefficient_observed"] = (
            phys["v_over_mplanck_observed"] / phys["loop"] ** 10)
    pool = tree["pool"]
    njobs = pool.get("njobs")
    pool["njobs_effective"] = int(njobs) if njobs else max(
        1, min(int(pool["njobs_cap"]), (os.cpu_count() or 1) // 2))
    tree["config_path"] = str(PATH)
    tree["config_present"] = PATH.exists()


CFG = Node({})


def reload():
    """Re-read config.yaml + the environment into the same CFG object."""
    tree = _build()
    _derive(tree)
    CFG._d.clear()
    CFG._d.update(tree)
    describe.cache_clear()
    return CFG


@functools.lru_cache(maxsize=1)
def describe():
    """The effective configuration, flat, for a report or a verdict file."""
    out = {}

    def walk(node, prefix=""):
        for key, value in node._d.items():
            name = f"{prefix}{key}"
            if isinstance(value, dict):
                walk(Node(value), name + ".")
            else:
                out[name] = value

    walk(CFG)
    return out


def smoke():
    """True when ``CL10_SMOKE`` asks for the shrunken experiment grids."""
    return bool(CFG.get("smoke", False))


reload()   # the effective configuration, ready at import
