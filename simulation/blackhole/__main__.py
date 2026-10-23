#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Entry point: ``python -m simulation.blackhole`` (headless video) or
``python -m simulation.blackhole --interactive`` (the live VisPy window
with the fire-particle keyboard shortcuts)."""

import sys

from . import runner

runner.main(sys.argv[1:])
