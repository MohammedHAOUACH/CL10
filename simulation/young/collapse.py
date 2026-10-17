#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The **collapse** -- the wave-packet reduction, made visible.

The screen must show the **reduction of the wave packet as the reality**:

1.  The **wave** (the superposition) is the evolving field ``psi``.
2.  A **detection** = a point ``(x, y)`` sampled from the on-screen
    probability ``psi**2`` (the Born rule, 0-parameter: the intensity *is*
    the distribution).
3.  The **collapse** = the wave packet **contracting** onto that point: the
    field is multiplied by a Gaussian that **shrinks** onto ``(x, y)`` over a
    few frames (the packet visibly reduces to the detected particle), then the
    particle **lands** on the screen and is added to the accumulation.
4.  Repeated over many events, the **detected points build the interference
    pattern** -- the *reality* (the measured pattern) emerges from the
    individual *reductions*.  This is the Young experiment as a **measurement
    process**, the 9-to-1 / Z2 projection of the CL10 root made visible.

The sampling and the contraction work on **host** numpy (the field is already
there for rendering); the heavy wave propagation stays on the GPU.
"""

import numpy as np


class Collapse:
    def __init__(self, cfg):
        self.show = bool(cfg["show_collapse"])
        self.frames = int(cfg["collapse_frames"])
        self.peak = float(cfg["collapse_peak"])
        self.n_detections = int(cfg["n_detections"])

        # the detected points (the accumulated "reality" on the screen)
        self.detections = []  # list of (x, y) grid indices
        # the current collapse event: (x, y, frame_index) or None
        self.event = None
        # the sigma range of the contraction (large -> small)
        self.sigma_start = float(cfg["source_width"])
        self.sigma_end = 1.5

    # -- a detection: sample a point from the on-screen psi**2 ---------------
    def sample(self, intensity):
        """Inverse-CDF sample of one point from the 2-D (Lx, Ly) intensity.
        Returns ``(x, y)`` grid indices, or ``None`` if the screen is empty."""
        I = np.asarray(intensity, dtype=np.float64)
        total = I.sum()
        if total <= 0.0:
            return None
        p = I.ravel() / total
        cdf = np.cumsum(p)
        u = np.random.random()
        flat = int(np.searchsorted(cdf, u, side="right"))
        flat = min(flat, p.size - 1)
        Ly = I.shape[1]  # I is (Lx, Ly); C-order ravel is x*Ly + y
        x = flat // Ly
        y = flat % Ly
        return int(x), int(y)

    def set_event(self, x, y):
        if self.show:
            self.event = (x, y, 0)

    def advance_event(self):
        if self.event is not None:
            x, y, f = self.event
            if f + 1 >= self.frames:
                self.event = None  # the contraction is over -> the particle lands
            else:
                self.event = (x, y, f + 1)

    def event_active(self):
        return self.event is not None

    # -- the contraction: the packet reducing onto the detected point ---------
    def contract(self, field):
        """The contracted field: ``field * (a shrinking Gaussian on x-y)``.
        ``field`` is a host (Lx, Ly, Lz) numpy array.  Returns a numpy array of
        the same shape (the reduced packet)."""
        if self.event is None:
            return field
        x0, y0, f = self.event
        Lx, Ly, Lz = field.shape
        # the sigma shrinks from sigma_start -> sigma_end over the frames
        frac = f / max(1, self.frames - 1)
        sigma = self.sigma_start * (self.sigma_end / self.sigma_start) ** frac
        # the peak grows as the packet shrinks (the "flash" of the reduction)
        pk = 1.0 + (self.peak - 1.0) * frac
        gx = np.arange(Lx, dtype=np.float32)
        gy = np.arange(Ly, dtype=np.float32)
        X, Y = np.meshgrid(gx, gy, indexing="ij")
        G = np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2.0 * sigma * sigma))
        return field * G[:, :, None] * pk

    def n_detected(self):
        return len(self.detections)

    def add(self, x, y):
        if self.n_detected() < self.n_detections:
            self.detections.append((x, y))
