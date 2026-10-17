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
import numpy as np


class WaveSim:
    def __init__(self, cfg):
        self.Lx, self.Ly, self.Lz = (int(v) for v in cfg["grid"])
        self.r = float(cfg["r"])
        self.wavelength = float(cfg["wavelength"])
        self.source_width = float(cfg["source_width"])
        self.source_amplitude = float(cfg["source_amplitude"])
        self.source_z = int(cfg["source_z"])
        self.ramp_steps = int(cfg["ramp_steps"])
        self.barrier_z = int(cfg["barrier_z"])
        self.screen_z = int(cfg["screen_z"])
        self.slit_separation = float(cfg["slit_separation"])
        self.slit_width = float(cfg["slit_width"])
        self.slit_height = float(cfg["slit_height"])
        self.sponge_width = int(cfg["sponge_width"])

        # the source model: "single" = one localized electron packet at a time
        # (the real Tonomura experiment); "continuous" = the old driven beam.
        self.source_mode = str(cfg["source_mode"])
        # the packet shape: "point" = a 3D Gaussian (a true single electron);
        # "sheet" = full in y (a vertical sheet, clean x-fringes like the ref).
        self.source_shape = str(cfg["source_shape"])
        self.source_y_width = float(cfg["source_y_width"])
        # the z-width of the emitted packet (a thin sliver of the beam).
        self.packet_z_width = float(cfg["packet_z_width"])

        self.psi = B.zeros((self.Lx, self.Ly, self.Lz))
        self.psi_prev = B.zeros((self.Lx, self.Ly, self.Lz))
        self.t = 0

        # the carrier: the spatial frequency is matched to the FDTD dispersion
        # (phi = 2*pi/lambda; omega = asin(r*sin(phi))) so the packet has a
        # clean central wavenumber (no phase pile-up).
        import math
        self.phi = 2.0 * math.pi / self.wavelength
        self.omega = math.asin(self.r * math.sin(self.phi))
        self._gx = B.arange(self.Lx, dtype=B._E.float32)
        self._gy = B.arange(self.Ly, dtype=B._E.float32)
        self._gz = B.arange(self.Lz, dtype=B._E.float32)
        # the z-part of the emitted packet: a Gaussian in z * the +z carrier.
        dz = self._gz - float(self.source_z)
        self._z_packet = (
            B.exp(-dz * dz / (2.0 * self.packet_z_width ** 2)) * B.cos(self.phi * dz)
        )
        # the "continuous" drive mask (a sheet: Gaussian in x, full in y).
        mx = B.exp(-((self._gx - self.Lx / 2.0) ** 2) / (2.0 * self.source_width ** 2))
        self.source_mask = mx[:, None] * B.ones((1, self.Ly), dtype=B._E.float32)
        self._t_since_emit = 0
        self._emit_count = 0
        # the total |psi|^2 of the last emitted packet (the wave is
        # conservative, so this only decays when the sponge absorbs it) --
        # used by absorbed() to decide when the next electron may fire.
        self._emit_energy = 0.0
        # the screen buffer: the on-screen |psi|^2 *integrated over time* while
        # the packet passes the detection screen.  A single instant of a
        # traveling FDTD wave is a small blob (the front), not the interference
        # pattern -- the fringes only emerge once the diffracted field has
        # crossed the slits and interfered across the screen.  Integrating the
        # intensity over the passage gives the time-averaged distribution,
        # which IS the Born-rule probability the real Tonomura screen records.
        self._screen_buf = np.zeros((self.Lx, self.Ly), dtype=np.float64)
        self._screen_armed = False  # start integrating once the front reaches it

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
        """One leapfrog step of the 0-parameter discrete wave equation.  In
        ``"continuous"`` mode the source plane is re-driven every step (the old
        beam); in ``"single"`` mode nothing is driven -- the packet that
        ``emit_packet`` launched simply propagates (and spreads)."""
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
        self.psi_prev = p
        self.psi = new
        if self.source_mode == "continuous":
            # the driven Gaussian-beam source (ramped in to avoid a shock)
            ramp = 1.0 if self.t >= self.ramp_steps else (self.t / float(self.ramp_steps))
            self.psi[:, :, self.source_z] = (
                self.source_amplitude * ramp * self.source_mask * B.cos(self.omega * self.t)
            )
        self._t_since_emit += 1
        self.t += 1
        return self.t

    # -- a single electron: launch one localized packet ------------------------
    def emit_packet(self):
        """Emit **one** electron: clear the (already ~0) field, then launch a
        localized 3-D wave packet (a Gaussian in x, y and z modulated by the +z
        carrier) at the source, with an **initial +z velocity** so it flies
        toward the slits and the screen.  This is the Tonomura single-electron
        source -- one packet in the apparatus at a time.  It is called **only
        after** the previous packet has been absorbed by the sponge (the field
        has decayed to ~0), so the clear here only removes numerical residue --
        there is **no mid-flight reset**: the electron rides all the way to the
        sponge between one emission and the next."""
        self.psi.fill(0)
        self.psi_prev.fill(0)
        # fresh screen buffer for the next electron (the previous one's
        # time-averaged pattern was already sampled before this fired)
        self._screen_buf.fill(0)
        self._screen_armed = False
        # the (x, y) envelope of the packet (a 3-D Gaussian, centered)
        if self.source_shape == "sheet":
            # full in y (a vertical sheet -> clean x-fringes, like the reference)
            ex = B.exp(-((self._gx - self.Lx / 2.0) ** 2) / (2.0 * self.source_width ** 2))
            xy = ex[:, None] * B.ones((1, self.Ly), dtype=B._E.float32)
        else:
            # a true single electron: a 3-D Gaussian (localized in x and y)
            ex = B.exp(-((self._gx - self.Lx / 2.0) ** 2) / (2.0 * self.source_width ** 2))
            ey = B.exp(-((self._gy - self.Ly / 2.0) ** 2) / (2.0 * self.source_y_width ** 2))
            xy = ex[:, None] * ey[None, :]
        # the full 3-D packet: the (x, y) envelope * the z-part (a carrier in z)
        self.psi[:] = self.source_amplitude * xy[:, :, None] * self._z_packet[None, None, :]
        # the initial +z velocity: the field at t-1 is the packet 1 cell behind
        # (so the leapfrog sees it moving toward +z, not back toward the source)
        self.psi_prev[:, :, :-1] = self.psi[:, :, 1:]
        self._t_since_emit = 0
        self._emit_count += 1
        # remember the packet's total energy (conservative until the sponge
        # takes it) -- absorbed() compares the live total against this.
        self._emit_energy = float((self.psi * self.psi).sum())
        return self

    # -- reset the field (the current electron is done) -----------------------
    def reset(self):
        """Clear the field (the current electron has been detected; ready to
        fire the next one)."""
        self.psi.fill(0)
        self.psi_prev.fill(0)
        self.t = 0
        self._t_since_emit = 0
        return self

    def at_screen(self):
        """Has the emitted packet's **front** reached the detection screen?  (the
        on-screen |psi|^2 is non-negligible.)  In ``"single"`` mode this triggers
        the one Born-rule *detection* of the electron (a point sampled from the
        on-screen intensity) -- the wave is **not** consumed here; it keeps
        propagating toward the sponge."""
        I = self.psi[:, :, self.screen_z] ** 2
        return float(I.max()) > 1e-6

    def absorbed(self, frac=0.02):
        """Has the packet been **absorbed by the sponge**?  The FDTD wave is
        **conservative**: a spreading packet keeps its total ``sum|psi|^2`` as
        it diffracts (only its *peak* falls), so the peak is a bad "gone"
        signal.  The total energy, by contrast, decays **only** when the
        sponge takes the wave at the boundary.  Once the live total has fallen
        below ``frac`` of the emitted energy, the sponge has absorbed the
        electron -- the *next* one may fire (no reset in between: the electron
        rode all the way to the sponge)."""
        if self._emit_energy <= 0.0:
            return False
        return float((self.psi * self.psi).sum()) < frac * self._emit_energy

    # -- the screen buffer: integrate the on-screen |psi|^2 over the passage --
    def arm_screen(self):
        """Start integrating the on-screen intensity (call when the front
        reaches the screen)."""
        self._screen_armed = True

    def accumulate_screen(self):
        """Add the current on-screen |psi|^2 to the running buffer (call every
        step while armed).  The buffer is the time-integrated intensity -- the
        true Born-rule probability distribution the Tonomura screen records."""
        if not self._screen_armed:
            return
        self._screen_buf += self.screen_intensity()

    def screen_pattern(self):
        """The time-integrated on-screen |psi|^2 (host numpy, Lx x Ly) -- the
        interference pattern to sample a Born-rule detection from."""
        return self._screen_buf

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
