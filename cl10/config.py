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
    3. the built-in defaults below (a missing key never crashes the engine)

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

# the built-in defaults: the same values as config.yaml (pool.njobs excepted -
# without the file the pool sizes itself), so a missing key or a missing file
# degrades to the documented behaviour instead of failing
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
    "physics": {
        "m_planck_gev": 1.220890e19,
        "proton_mass_gev": 0.93827208943,
        "v_electroweak_gev": 246.22,
        "omega_dm_h2": 0.120,
        "weak_scale_gev": 100.0,
        "hydrogen_rydberg_ev": 13.605693,
        "n_gen": 10,
        "n_space": 3,
        "n_subsites": 8,
        "n_internal": 32,
        "n_sector": 16,
        "loop_factor": None,
    },
    "units": {"gev_to_s": 6.582119569e-25, "s_to_yr": 3.15576e7,
              "gev2_to_cm2": 3.89379338e-28, "c_cm_s": 2.99792458e10},
    "sm": {
        "m_z_gev": 91.1876,
        "beta1": 41 / 10, "beta2": -19 / 6, "beta3": -7,
        "alpha1_mz": 0.0169, "alpha2_mz": 0.0338, "alpha3_mz": 0.118,
        "beta21": 199 / 25, "beta22": 25 / 2 - 27 / 50 * 6,
        "beta23": -143 / 3 + 23 / 15 * 6,
        "sphaleron": 28.0 / 79.0,
        "higgs_mass_gev": 125.25,
        "delta_b_obs": 6.12e-10,
        "m_unification_obs_gev": 2.0e16,
        "alpha_em_inv": 127.9,
    },
    "angle": {
        "default_branch": "implemented",
        "branches": {"implemented": 0.2, "closed": math.pi / 4, "half": math.pi / 2},
        "alpha_10": 0.04,
        "g0_from_n_space": True,
    },
    "desert": {"theta_root_rad": 1.7308, "m_x_gev": 1.8e16, "alpha_x": 0.04,
               "m_x_closed_gev": 1.235e17},
    "bounds": {"lz_sigma_si_cm2": 1.5e-47, "dd_sigma_si_cm2": 2.0e-47},
    "cosmology": {"g_star": 106.75, "hubble_coeff": 1.66,
                  "entropy_density_cm3": 2891.2, "rho_c_over_h2_gev_cm3": 1.0537e-5},
    "dm": {"f_n_nucleon": 3.0e-4, "sigma_fiducial_cm3s": 3.0e-26,
           "sigma_relic_target_cm3s": 2.2e-26},
    "neutrino": {"m_atm_ev": 0.0501, "m_sol_ev": 0.0086,
                 "sum_mnu_max_ev": 0.12},
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
