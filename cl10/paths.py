#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Where the experiments live.

The verification program is organised as ``experiments/<category>/<Vxx>/``.  An
experiment refers to a sibling (to import its module, to read its verdict JSON)
or to its own output by name, never by category: the category is resolved here,
so moving an experiment between categories moves a directory and nothing else.
"""
import glob
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPERIMENTS = os.path.join(ROOT, "experiments")


def exp(*parts):
    """Path inside ``experiments/<category>/<Vxx>/...``, resolved by name.

    ``exp("V65_wilson")`` is the experiment directory;
    ``exp("V65_wilson", "results", "v65_wilson.json")`` is a file inside it.
    Raises ``FileNotFoundError`` if the name matches no experiment or matches
    more than one.
    """
    hits = glob.glob(os.path.join(EXPERIMENTS, "*", *parts))
    if len(hits) != 1:
        raise FileNotFoundError(f"experiment path {parts!r}: {len(hits)} matches")
    return hits[0]
