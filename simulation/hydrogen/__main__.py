#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Entry point: ``python -m simulation.hydrogen``.

Pins **``CUDA_VISIBLE_DEVICES``** to the ``gpu`` id from ``config.yaml``
(default **1** -- the free card) *before* the CuPy import, so the whole
simulation runs on that GPU.  Then dispatches to the runner: the headless
video by default, or the live VisPy window with ``--interactive``.
"""

import os
import sys

import yaml

_here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_here, "config.yaml")) as f:
    _cfg = yaml.safe_load(f)

# pin the GPU BEFORE the CuPy import (backend.py imports cupy)
os.environ["CUDA_VISIBLE_DEVICES"] = str(_cfg.get("gpu", 1))

from . import backend as B  # noqa: E402

B.set_field_dtype(_cfg.get("field_dtype", "float32"))

from . import runner  # noqa: E402

runner.main(sys.argv[1:])
