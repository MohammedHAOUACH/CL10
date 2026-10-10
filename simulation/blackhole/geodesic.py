#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The geodesic engine of the 3-D black-hole simulation.

The physics is **the same as fig17** (``make_figures.fig17_gravity_3d``):
the Schwarzschild metric ``metric_schwarzschild(M, x)`` from
``cl10.gravity``, the Christoffel symbols by **central finite difference**
(``gravity._christoffel_at`` with the 0-parameter step ``H_FD``), and the
geodesic equation ``d^2 x^m/dl^2 = -Gamma^m_{ab} dx^a/dl dx^b/dl``.

The difference is the *integration*: fig17 hands the whole problem to
``solve_ivp`` (DOP853) and draws the finished curve.  Here the same ODE is
integrated **incrementally with a fixed-step RK4**, one chunk per rendered
frame, so each fired particle leaves a **trail that grows in real time** —
you watch the geodesic being constructed (the precessing rosette, the
plunge through the horizon, the deflected fly-by).

A particle is a timelike geodesic: the initial 4-velocity is built from a
LOCAL orthonormal velocity (v_r, v_t) at the start event, exactly as in
fig17 (``gamma / sqrt(1-2M/r0)`` for the t-component, ``gamma*v/r0`` for the
phi-component), so the constants of motion (E, L) emerge from the geometry,
they are never imposed.

Statuses: ``orbiting`` (still integrating), ``captured`` (r < 2M — it
crossed the horizon), ``escaped`` (r > escape_radius), ``stopped`` (trail
full or affine parameter exhausted).
"""

import numpy as np

from cl10 import gravity

H = gravity.H_FD          # the 0-parameter finite-difference step (fig17)


class Particle:
    """One fired test particle: its state u = (t, r, th, ph, dt, dr, dth, dph)
    and the growing 3-D trail (x, y, z) in well coordinates."""

    def __init__(self, M, preset, dlam, max_trail, depth=1.0):
        M = float(M)
        self.depth = float(depth)
        r0 = float(preset["r0"])
        vr, vt = float(preset["vr"]), float(preset["vt"])
        v2 = vr * vr + vt * vt
        if v2 >= 1.0:
            raise ValueError("local speed must stay < c")
        gamma = 1.0 / np.sqrt(1.0 - v2)
        f0 = 1.0 - 2.0 * M / r0
        ph0 = float(preset.get("phi0", 0.0))
        # the local orthonormal -> coordinate 4-velocity (fig17's recipe)
        self.u = np.array([0.0, r0, np.pi / 2, ph0,
                           gamma / np.sqrt(f0), vr * gamma,
                           0.0, vt * gamma / r0], dtype=float)
        self.M = M
        self.dlam = float(dlam)
        self.color = preset.get("color", "#ffd23f")
        self.name = preset.get("name", "particle")
        self.status = "orbiting"
        self.lam = 0.0
        self.trail = [self._xyz()]
        self._max_trail = int(max_trail)

    def _xyz(self):
        """The current position in the well coordinates (x, y, Phi)."""
        _, r, _, ph = self.u[:4]
        return (r * np.cos(ph), r * np.sin(ph),
                self.depth * -self.M / max(r, 2.0 * self.M))

    def _dudl(self, u):
        t, r, th, ph = u[:4]
        G = gravity._christoffel_at(
            lambda x: gravity.metric_schwarzschild(self.M, x),
            np.array([t, r, th, ph]), H)
        v = u[4:]
        dv = -np.einsum("mab,a,b->m", G, v, v)
        return np.concatenate([v, dv])

    def _stage_ok(self, u):
        """An RK4 stage is usable only while it stays outside the horizon:
        the Christoffels blow up at r <= 2M and a stage sampling there
        corrupts the whole step (the runaway seen as a fake 'escape')."""
        return u[1] > 2.0 * self.M * 1.001 and np.all(np.isfinite(u))

    def step(self, n):
        """Advance the geodesic by ``n`` fixed-step RK4 substeps, appending
        to the trail.  Returns the status after the chunk."""
        h = self.dlam
        for _ in range(n):
            if self.status != "orbiting":
                break
            u = self.u
            k1 = self._dudl(u)
            u2 = u + 0.5 * h * k1
            if not self._stage_ok(u2):
                self.trail.append(self._xyz())
                self.status = "captured"
                break
            k2 = self._dudl(u2)
            u3 = u + 0.5 * h * k2
            if not self._stage_ok(u3):
                self.trail.append(self._xyz())
                self.status = "captured"
                break
            k3 = self._dudl(u3)
            u4 = u + h * k3
            if not self._stage_ok(u4):
                self.trail.append(self._xyz())
                self.status = "captured"
                break
            k4 = self._dudl(u4)
            self.u = u + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            self.lam += h
            r = self.u[1]
            # the horizon: captured (the trail keeps its last visible point)
            if r <= 2.0 * self.M * 1.001:
                self.trail.append(self._xyz())
                self.status = "captured"
                break
            if r > self._escape:
                self.trail.append(self._xyz())
                self.status = "escaped"
                break
            self.trail.append(self._xyz())
            if len(self.trail) >= self._max_trail:
                self.status = "stopped"
                break
        return self.status

    _escape = 70.0


class BlackHole:
    """The simulation state: the hole (M = 1, the 0-parameter Planck-mass
    black hole) and the list of fired particles with their growing trails."""

    def __init__(self, cfg):
        self.M = float(cfg.get("M", 1.0))
        self.dlam = float(cfg.get("dlambda", 0.04))
        self.max_trail = int(cfg.get("max_trail", 6000))
        self._escape = float(cfg.get("escape_radius", 70.0))
        # the visual depth exaggeration of the well (rendering only: the
        # physics stays Phi = -M/r; the factor just makes the funnel read)
        self.depth = float(cfg.get("well_depth", 12.0))
        self.presets = cfg.get("presets", {})
        self.particles = []

    def fire(self, key):
        """Fire the preset particle ``key`` ('1'..'5', 'p').  Returns the
        Particle, or None for an unknown key."""
        preset = self.presets.get(str(key))
        if preset is None:
            return None
        p = Particle(self.M, preset, self.dlam, self.max_trail,
                     self.depth)
        p._escape = self._escape
        self.particles.append(p)
        return p

    def advance(self, steps_per_frame):
        """One rendered frame: every live particle integrates its chunk."""
        for p in self.particles:
            p.step(steps_per_frame)

    def counts(self):
        c = {"orbiting": 0, "captured": 0, "escaped": 0, "stopped": 0}
        for p in self.particles:
            c[p.status] += 1
        return c

    def horizon_ring(self, n=200):
        a = np.linspace(0.0, 2.0 * np.pi, n)
        z = np.full_like(a, self.depth * -self.M / (2.0 * self.M))
        return 2.0 * self.M * np.cos(a), 2.0 * self.M * np.sin(a), z

    def well(self, R, n=90):
        """The embedding well: X, Y, Phi = -M/r capped at the horizon
        (identical to fig17's left panel)."""
        g = np.linspace(-R, R, n)
        X, Y = np.meshgrid(g, g)
        rr = np.sqrt(X ** 2 + Y ** 2)
        return X, Y, self.depth * -self.M / np.maximum(rr, 2.0 * self.M)
