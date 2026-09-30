#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The 3-D wave of the Young experiment.

A **real** scalar field ``psi`` on a big (Lx, Ly, Lz) box, propagated with the
**0-parameter discrete wave equation** (the same 6-neighbour Laplacian leapfrog
that CL10 uses for its massless graviton, ``cl10/qca.py::wave_step``):

    psi_{t+1} = 2 psi_t - psi_{t-1} + r^2 * (sum 6 neighbours - 6 psi_t)

with ``r = c dt/dx = 0.5`` (CFL in 3-D: r <= 1/sqrt(3)).  The **Gaussian
beam** is a driven boundary at ``z = 0``; a **sponge** absorbs the boundaries;
the **barrier** is a hard mask at ``z = barrier_z`` that passes only the two
**slits**.  The wave therefore propagates, diffracts through the two slits, and
interferes on the **detection screen** at ``z = screen_z``.

Everything (the Laplacian, the sponge, the barrier) is applied on the GPU
backend (``backend.py``).
"""

from . import backend as B


class WaveSim:
    def __init__(self, cfg):
        self.Lx, self.Ly, self.Lz = (int(v) for v in cfg["grid"])
        self.r = float(cfg["r"])
        self.wavelength = float(cfg["wavelength"])
        self.source_width = float(cfg["source_width"])
        self.source_amplitude = float(cfg["source_amplitude"])
        self.ramp_steps = int(cfg["ramp_steps"])
        self.barrier_z = int(cfg["barrier_z"])
        self.screen_z = int(cfg["screen_z"])
        self.slit_separation = float(cfg["slit_separation"])
        self.slit_width = float(cfg["slit_width"])
        self.slit_height = float(cfg["slit_height"])
        self.sponge_width = int(cfg["sponge_width"])

        self.psi = B.zeros((self.Lx, self.Ly, self.Lz))
        self.psi_prev = B.zeros((self.Lx, self.Ly, self.Lz))
        self.t = 0

        # the temporal frequency of the source, matched to the FDTD dispersion
        # so the driven boundary launches a clean propagating wave (no phase
        # pile-up).  phi = 2*pi/lambda is the spatial phase per cell.
        import math
        phi = 2.0 * math.pi / self.wavelength
        self.omega = math.asin(self.r * math.sin(phi))

        # the Gaussian source envelope (a 2-D mask on x-y at the source plane)
        gx = B.arange(self.Lx, dtype=B._E.float32)
        gy = B.arange(self.Ly, dtype=B._E.float32)
        X, Y = B.meshgrid(gx, gy)
        self.source_mask = B.exp(
            -((X - self.Lx / 2.0) ** 2 + (Y - self.Ly / 2.0) ** 2)
            / (2.0 * self.source_width ** 2)
        )

        self._build_sponge()
        self._build_barrier()

    # -- the absorbing boundary (a quadratic taper that kills reflections) ----
    def _build_sponge(self):
        """A quadratic absorbing taper: 1.0 in the bulk, -> 0 at every
        boundary (kills the FDTD reflections)."""
        s = float(self.sponge_width)
        X3, Y3, Z3 = B.meshgrid(
            B.arange(self.Lx, dtype=B._E.float32),
            B.arange(self.Ly, dtype=B._E.float32),
            B.arange(self.Lz, dtype=B._E.float32),
        )
        # distance (in cells) of each cell to its nearest boundary face
        d3 = B.minimum(
            B.minimum(
                B.minimum(X3, self.Lx - 1 - X3),
                B.minimum(Y3, self.Ly - 1 - Y3),
            ),
            B.minimum(Z3, self.Lz - 1 - Z3),
        )
        q = B.clip(1.0 - d3 / s, 0.0, 1.0)
        self.sponge = 1.0 - q * q  # = 1 in the bulk, -> 0 at the edge

    # -- the double-slit barrier mask (1 inside a slit, 0 elsewhere) ---------
    def _build_barrier(self):
        m = B.zeros((self.Lx, self.Ly), dtype=B._E.float32)
        xc = self.Lx / 2.0
        yc = self.Ly / 2.0
        for sign in (-1.0, +1.0):
            sc = xc + sign * self.slit_separation / 2.0
            x0 = max(0, int(sc - self.slit_width / 2.0))
            x1 = min(self.Lx, int(sc + self.slit_width / 2.0) + 1)
            y0 = max(0, int(yc - self.slit_height / 2.0))
            y1 = min(self.Ly, int(yc + self.slit_height / 2.0) + 1)
            m[x0:x1, y0:y1] = 1.0
        self.barrier = m
        # the slit rectangles (x0, x1, y0, y1) for the 3-D scene geometry
        self.slit_rects = []
        for sign in (-1.0, +1.0):
            sc = xc + sign * self.slit_separation / 2.0
            self.slit_rects.append(
                (
                    int(sc - self.slit_width / 2.0),
                    int(sc + self.slit_width / 2.0) + 1,
                    int(yc - self.slit_height / 2.0),
                    int(yc + self.slit_height / 2.0) + 1,
                )
            )

    # -- one FDTD step --------------------------------------------------------
    def step(self):
        p = self.psi
        pp = self.psi_prev
        r2 = self.r * self.r
        lap = B.zeros_like(p)
        lap[:-1, :, :] += p[1:, :, :]
        lap[1:, :, :] += p[:-1, :, :]
        lap[:, :-1, :] += p[:, 1:, :]
        lap[:, 1:, :] += p[:, :-1, :]
        lap[:, :, :-1] += p[:, :, 1:]
        lap[:, :, 1:] += p[:, :, :-1]
        lap -= 6.0 * p
        new = 2.0 * p - pp + r2 * lap
        # the sponge (kill reflections)
        new = new * self.sponge
        # the double-slit barrier (hard mask at the barrier plane)
        new[:, :, self.barrier_z] = new[:, :, self.barrier_z] * self.barrier
        # the driven Gaussian-beam source at z = 0 (ramped in to avoid a shock)
        ramp = 1.0 if self.t >= self.ramp_steps else (self.t / float(self.ramp_steps))
        self.psi_prev = p
        self.psi = new
        self.psi[:, :, 0] = (
            self.source_amplitude * ramp * self.source_mask * B.cos(self.omega * self.t)
        )
        self.t += 1
        return self.t

    # -- outputs to the host --------------------------------------------------
    def field(self):
        """The current field on the host (full resolution)."""
        return B.asnumpy(self.psi)

    def intensity(self):
        """The on-the-fly intensity ``psi**2`` on the host (full resolution)."""
        return B.asnumpy(self.psi * self.psi)

    def screen_intensity(self):
        """The 2-D (Lx, Ly) intensity on the detection screen (host numpy)."""
        return B.asnumpy(self.psi[:, :, self.screen_z] ** 2)
