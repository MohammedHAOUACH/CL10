#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CL10 — Causal Clifford Theory: the Standard Model as a single Clifford root.

The single root object is the Clifford algebra ``Cl(10)`` of ten mutually
anticommuting generators.  The spacetime (1,3) is *inside* it (time = h_1,
space = i h_2..4); the matter is the 16-D Weyl spinor; the gauge is the 45
(Spin(10) adjoint); the statistics are the Z2 grading; gravity is the metric
dynamics of the (1,3) sub-algebra.  **0 free parameters in the algebra**;
the one 't Hooft angle of the dynamics carries a 0-parameter selection
rule (the *reading* is the audit's mapping convention).
"""

from . import engine, clifford, spacetime, matter, gauge, qca, u10pred  # noqa: F401

__all__ = ["engine", "clifford", "spacetime", "matter", "gauge", "qca", "u10pred"]
__version__ = "1.2.0"
