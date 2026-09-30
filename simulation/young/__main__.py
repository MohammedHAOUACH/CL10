#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Entry point: ``python -m simulation.young``.

Pins **``CUDA_VISIBLE_DEVICES``** to the ``gpu`` id from ``config.yaml``
(default **1** -- the free card) *before* the CuPy import, so the whole
simulation runs on that GPU.  Then dispatches to the runner (headless video by
default, ``--interactive`` for the live window + parameter panel).
"""

import os
import sys

import yaml

_here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_here, "config.yaml")) as f:
    _cfg = yaml.safe_load(f)

# pin the GPU BEFORE the CuPy import (backend.py imports cupy)
os.environ["CUDA_VISIBLE_DEVICES"] = str(_cfg.get("gpu", 1))

from . import runner  # noqa: E402

runner.main(sys.argv[1:])
