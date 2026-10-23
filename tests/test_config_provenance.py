"""Every config path the code reads is owned by cl10/config.yaml.

The publication's traceability claim is that each published number has one
named owner: it is declared in `cl10/config.yaml` (or measured by a named
experiment and read back from its verdict JSON), never restated as a decimal
somewhere else. A `CFG.a.b.c` lookup that resolves through neither the file nor
a derivation in `cl10/_derive` breaks that claim silently - the file does not
 say what the run used, and the fallback becomes the real input.

So each lookup must be one of:
  * a key in the YAML file, or
  * a value that `cl10/config.py::_derive` computes from declared keys, listed
    below with its formula. The formulas are re-checked here numerically, so a
    derived value cannot drift away from the inputs it claims to follow.
"""
import math
import os
import pathlib
import re

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
CFG_PATH = ROOT / "cl10" / "config.yaml"

# `CFG.a.b`, `CFG["a"]["b"]` and `CFG.get("a.b")` are the spellings in use.
_ATTR = re.compile(r"\bCFG\.([a-z_][A-Za-z0-9_]*(?:\.[A-Za-z0-9_]+)+)")
_GET = re.compile(r"\bCFG\.get\(\s*[\"']([A-Za-z0-9_.]+)[\"']")
_ITEM = re.compile(r"\bCFG((?:\[[\"'][A-Za-z0-9_]+[\"']\])+)")

SEARCH = ["cl10", "experiments", "make_figures", "simulation", "tests"]


def _raw():
    return yaml.safe_load(CFG_PATH.read_text()) or {}


def _sources():
    for part in SEARCH:
        for p in sorted((ROOT / part).rglob("*.py")):
            yield p


def _paths():
    out = {}
    for p in _sources():
        text = p.read_text()
        for m in _ATTR.finditer(text):
            out.setdefault(m.group(1), set()).add(p)
        for m in _GET.finditer(text):
            out.setdefault(m.group(1), set()).add(p)
        for m in _ITEM.finditer(text):
            dotted = ".".join(re.findall(r"[\"']([A-Za-z0-9_]+)[\"']", m.group(1)))
            out.setdefault(dotted, set()).add(p)
    return {k: v for k, v in out.items() if not k.startswith("_d")}


def _present(tree, dotted):
    node = tree
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return False
        node = node[part]
    return True


def _cfg():
    from cl10 import config
    config.reload()
    return config.CFG


def test_config_yaml_exists():
    assert CFG_PATH.is_file(), f"missing {CFG_PATH}"


def test_every_cfg_lookup_is_declared_in_the_yaml_file():
    raw = _raw()
    roots = set(raw)
    derived = set(_derived())
    missing = []
    for dotted, files in _paths().items():
        head = dotted.split(".")[0]
        if head not in roots:
            continue          # an unrelated local object that happens to be CFG
        if _present(raw, dotted) or dotted in derived:
            continue
        where = ", ".join(sorted(str(f.relative_to(ROOT)) for f in files)[:6])
        missing.append(f"CFG.{dotted}  ({where})")
    assert not missing, ("these config lookups exist in neither the file nor a "
                         "declared derivation, so nothing owns them:\n  "
                         + "\n  ".join(sorted(missing)))


def _derived():
    """The keys `_derive` computes, each with the formula it must obey."""
    raw = _raw()
    phys, angle, pool = raw["physics"], raw["angle"], raw["pool"]
    branch = str(angle.get("branch") or angle["default_branch"]).strip().lower()
    return {
        "physics.n_qubits": (float(int(phys["n_internal"]).bit_length() - 1),
                             "log2 of physics.n_internal"),
        "physics.loop": (float(phys["loop_factor"]) if phys.get("loop_factor")
                         else 1.0 / (16.0 * math.pi),
                         "physics.loop_factor"),
        "physics.v_over_mplanck_observed": (float(phys["v_electroweak_gev"])
                                            / float(phys["m_planck_gev"]),
                                            "v_electroweak_gev / m_planck_gev"),
        "physics.hierarchy_coefficient_observed": (
            (float(phys["v_electroweak_gev"]) / float(phys["m_planck_gev"]))
            / (float(phys["loop_factor"]) if phys.get("loop_factor")
               else 1.0 / (16.0 * math.pi)) ** 10,
            "v_over_mplanck_observed / loop**10"),
        "angle.branch": (branch, "angle.default_branch after alias resolution"),
        "angle.g10": (float(angle["branches"][branch]),
                      "angle.branches[angle.branch]"),
        "angle.g0": (1.0 / (2.0 ** int(phys["n_space"])), "1 / 2**n_space"),
        "pool.njobs_effective": (float(int(pool["njobs"]) if pool.get("njobs")
                                       else max(1, min(int(pool["njobs_cap"]),
                                                       (os.cpu_count() or 1) // 2))),
                                 "pool.njobs, else min(njobs_cap, cpus/2)"),
    }


def test_derived_values_follow_their_declared_formulas():
    cfg = _cfg()
    for dotted, (want, why) in _derived().items():
        node = cfg
        for part in dotted.split("."):
            node = getattr(node, part)
        if isinstance(want, str):
            assert str(node) == want, (f"CFG.{dotted} = {node!r} but its declared "
                                       f"derivation ({why}) gives {want!r}")
            continue
        got = float(node)
        assert abs(got - want) <= 1e-12 * max(1.0, abs(want)), (
            f"CFG.{dotted} = {got!r} but its declared derivation "
            f"({why}) gives {want!r}")


def test_derived_keys_are_not_secretly_written_in_the_file():
    # A key may be derived or declared, never both: two sources for one number
    # is exactly the defect this file exists to catch.
    raw = _raw()
    both = [d for d in _derived() if _present(raw, d)]
    assert not both, f"declared in config.yaml AND derived in _derive: {both}"


def test_cfg_lookup_count_is_meaningful():
    roots = set(_raw())
    hits = [d for d in _paths() if d.split(".")[0] in roots]
    assert len(hits) > 40, f"the CFG scan found only {len(hits)} config lookups"
