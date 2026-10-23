#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The **hydrogen atom as a CL10 QCA** -- a 0-parameter, honest model.

Two objects live on the Planck lattice:

* **the photon** -- the **free CL10 massless 3-D field**, in two parts
  (the same FDTD engine as ``simulation/young``):
  - ``F``     the **incoming** photon (fired at the atom; it passes and
    leaves the box, absorbed by the sponge);
  - ``F_out`` the **outgoing** photon (the **radiation** the atom emits;
    it is *not* coupled back into the atom -- it carries the emitted
    energy away, where the sponge absorbs it).

* **the atom** -- the **semiclassical Rabi / Fermi-golden-rule** model, the
  physics of an atom driven by a (classical) light field.  Its levels and
  wavefunctions are the **real hydrogen** solutions (0-parameter -- they
  are the known eigenstates of the hydrogen problem, not fitted):
  - the **energies** from the **Rydberg formula** ``E_n = -13.606/n**2``
    eV; the **transition frequencies** are the **gaps** between levels;
  - the **wavefunctions** are the **hydrogen eigenfunctions**
    ``1s, 2p, 3p, 3d`` (``R_nl(r) Y_lm``) evaluated on the 3-D lattice.

The **dynamics** (per lattice step):

* **absorption (the Fermi golden rule).**  The **local** photon amplitude
  ``F_local`` (the dipole approximation: the RMS of the field over the
  small atom, so the fast carrier oscillation averages to its envelope)
  drives the **dipole-allowed** transitions.  Each transition ``(i, j)``
  (``i`` the lower level) has a rate
  ``rate = A * F_local**2 * Lorentzian(omega_ph - (omega_j - omega_i))``,
  with a **Lorentzian** line shape (HWHM = ``gamma_line``, the natural
  linewidth set by the lifetime) and a **dipole matrix element** ``d_ij =
  <i|z|j>`` that encodes the **selection rules** (``Delta l = +/- 1``; the
  forbidden ones are ~0 and never fire).  Population flows ``i -> j``
  (absorption) and back (stimulated); for a ground-state atom it is pure
  absorption.  A photon **at** a gap (``delta ~ 0``) is absorbed (the atom
  is excited); one **between** the gaps (``|delta| >> gamma_line``) is not
  (it passes through).  This *is* the Rabi / Fermi-golden-rule selection
  rule, free of the Rabi cascade / power-broadening.

* **spontaneous emission (Einstein A).**  An excited level ``j`` decays at
  the rate ``gamma_sp``: its amplitude is damped and the **transition
  dipole** (which oscillates at the **gap** frequency) sources a photon
  into ``F_out`` via the Larmor law (the dipole *acceleration*).  The
  emitted wave is a decaying cosine at the gap -- the **natural line**.
  A **bare eigenstate** has a static (zero) transition dipole, so it
  **does not radiate** -- the ground state never self-decays.

* **the electron cloud** rendered in 3-D is the real part of the
  superposition ``psi = sum_k a_k e^{-i omega_k t} phi_k``: the **1s**
  is a single sphere, the **2p** a two-lobe dumbbell, the **3d** a four
  petal / doughnut shape -- the excitation **changes the cloud's shape**,
  the way it does in the real atom.

Nothing is fitted: the levels and shapes are the hydrogen physics, the
photon is the CL10 engine, and the only free knobs are the coupling
``A`` (the light-matter strength), ``h`` (the radiated source strength)
and ``gamma_sp`` (the radiative lifetime) -- the physical parameters of
the atom, not data fits.
"""

import math

import numpy as np

from . import backend as B
from cl10 import config as _config      # the shared measured inputs


# ---------------------------------------------------------------------------
# the hydrogen levels (0-parameter: the Rydberg formula)
# ---------------------------------------------------------------------------
# (label, n, l, m).  The energies are E_n = -RYDBERG_EV / n**2 eV.  Only the
# *gaps* matter for the transitions; we scale them into the lattice
# frequency band (a single unit-conversion, not a fit).
HYDROGEN_LEVELS = [
    ("1s", 1, 0, 0),
    ("2p", 2, 1, 0),
    ("3p", 3, 1, 0),
    ("3d", 3, 2, 0),
]
RYDBERG_EV = float(_config.CFG.physics.hydrogen_rydberg_ev)  # eV (cl10/config.yaml)
# the lattice scale of the spectrum (rad/step) is a unit choice: it lives in
# this simulation's config.yaml (lattice_gap_12, omega_base).


def _rydberg_gaps_eV():
    """The transition energies (eV) between every pair of levels, from the
    Rydberg formula."""
    n = [lv[1] for lv in HYDROGEN_LEVELS]
    E = [-RYDBERG_EV / nn ** 2 for nn in n]
    gaps = {}
    for i in range(len(n)):
        for j in range(len(n)):
            if i != j:
                gaps[(i, j)] = abs(E[j] - E[i])
    return gaps


def _hydrogen_phi(n, l, m, a0, r, costheta):
    """The (real, m = 0) hydrogen eigenfunction ``R_nl(r) Y_lm`` (a.u.),
    evaluated at radius ``r`` and ``cos(theta)``.  ``m`` must be 0."""
    # the radial part (standard, a0 = 1 in a.u.; we substitute r/a0)
    x = r / a0
    if (n, l) == (1, 0):
        R = 2.0 * np.exp(-x) / a0 ** 1.5
    elif (n, l) == (2, 1):
        R = x * np.exp(-x / 2.0) / (2.0 * math.sqrt(6.0) * a0 ** 1.5)
    elif (n, l) == (3, 1):
        R = x * (6.0 - x * x) * np.exp(-x / 3.0) / (81.0 * math.sqrt(6.0) * a0 ** 1.5)
    elif (n, l) == (3, 2):
        R = x * x * np.exp(-x / 3.0) / (81.0 * math.sqrt(30.0) * a0 ** 1.5)
    else:
        raise ValueError("no hydrogen radial for n=%d l=%d" % (n, l))
    # the spherical harmonic (m = 0 only)
    if l == 0:
        Y = 1.0 / math.sqrt(4.0 * math.pi)
    elif l == 1:
        Y = math.sqrt(3.0 / (4.0 * math.pi)) * costheta
    elif l == 2:
        Y = math.sqrt(5.0 / (16.0 * math.pi)) * (3.0 * costheta ** 2 - 1.0)
    else:
        raise ValueError("no spherical harmonic for l=%d" % l)
    return R * Y


# ---------------------------------------------------------------------------
# the atom + photon
# ---------------------------------------------------------------------------
class AtomQCA:
    """The hydrogen atom (the Rabi levels) coupled to the CL10 photon
    field (``F`` incoming, ``F_out`` emitted)."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.Lx, self.Ly, self.Lz = (int(v) for v in cfg["grid"])
        self.r = float(cfg["r"])
        self.a0 = float(cfg["bohr_radius"])
        # the light-matter coupling: A sets the absorption strength
        # (the Fermi-golden-rule rate scale) and h the emitted-photon
        # source strength.  These are the free knobs of the model (the
        # atom's light-matter strength and its radiative character).
        self.A = float(cfg["coupling_absorb"])
        self.h = float(cfg["coupling_emit"])
        # the Einstein-A decay rate (per step) of an excited level.
        self.gamma_sp = float(cfg["reaction_gamma_max"])
        # the natural linewidth (HWHM, rad/step) of the absorption line:
        # 0 ties it to the radiative lifetime (the energy-time uncertainty),
        # so a narrow line is a long-lived level.
        gl = float(cfg["line_width"])
        self.gamma_line = gl if gl > 0.0 else self.gamma_sp
        self.sponge_width = int(cfg["sponge_width"])
        self.photon_width = float(cfg["photon_width"])
        # the lattice scale of the atom's spectrum (a unit choice, not a fit)
        self.gap12 = float(cfg["lattice_gap_12"])
        self.omega_base = float(cfg["omega_base"])
        # the dipole-approximation widths, in cells
        self.emission_width = float(cfg["emission_width"])
        self.local_weight_width = float(cfg["local_weight_width"])

        self.cx, self.cy, self.cz = self.Lx // 2, self.Ly // 2, self.Lz // 2
        self.hx, self.hy, self.hz = (
            (self.Lx - 1) / 2.0, (self.Ly - 1) / 2.0, (self.Lz - 1) / 2.0
        )

        # -- the lattice geometry (the radius, the polar angle, the z arm) --
        ix = B.arange(self.Lx, dtype=B.FIELD)[:, None, None]
        iy = B.arange(self.Ly, dtype=B.FIELD)[None, :, None]
        iz = B.arange(self.Lz, dtype=B.FIELD)[None, None, :]
        rx = ix - self.hx
        ry = iy - self.hy
        rz = iz - self.hz
        self.R = (rx ** 2 + ry ** 2 + rz ** 2) ** 0.5
        # cos(theta) about the z axis (the dipole axis), bounded
        self.costheta = B.clip(rz / (self.R + 1e-6), -1.0, 1.0)
        self.z = rz  # the dipole arm (odd)
        self.dip = self.z  # the dipole pattern (the z coordinate)
        # the *emission* dipole pattern: the z arm localised to the atom
        # (the radiation is emitted from the atom's location; the raw z
        # grows without bound and would blow up the F_out source far away).
        ew = self.emission_width  # the emission source width (cells)
        self.dip_emit = self.z * B.exp(-self.R ** 2 / (2.0 * ew * ew))
        # -- the Coulomb well (rendering only: the potential the electron
        # sits in; a softened -1/r, the hydrogen 1/r shape).
        self.W = -3.0 / (self.R + 2.0)

        # -- the level frequencies (the Rydberg gaps, scaled to the lattice) --
        self.omega = []
        self.labels = [lv[0] for lv in HYDROGEN_LEVELS]
        n = [lv[1] for lv in HYDROGEN_LEVELS]
        E = [-RYDBERG_EV / nn ** 2 for nn in n]
        E1 = E[0]
        for k in range(len(n)):
            # the gap from the 1s base, scaled so (1s->2p) = self.gap12
            gap_eV = abs(E[k] - E1)
            # scale: the 1s->2p gap maps to self.gap12
            gap_lat = gap_eV * (self.gap12 / (abs(E[1] - E1)))
            self.omega.append(self.omega_base + gap_lat)
        # the (i, j) transition frequencies (rad/step)
        self.gaps = {}
        for i in range(len(self.omega)):
            for j in range(len(self.omega)):
                self.gaps[(i, j)] = abs(self.omega[j] - self.omega[i])
        # -- the energies in eV (for the HUD stick diagram) --
        self.level_E_eV = E  # the level energies E_n = -RYDBERG_EV/n**2
        self.gaps_eV = _rydberg_gaps_eV()  # (i, j) -> eV

        # -- the wavefunction shapes (the hydrogen eigenfunctions on the
        # lattice), normalised to unit L2 --
        self.phi = []
        for (label, nn, ll, mm) in HYDROGEN_LEVELS:
            p = _hydrogen_phi(nn, ll, mm, self.a0, self.R, self.costheta)
            pn = B.asnumpy(p)
            nrm = math.sqrt(float((pn * pn).sum()))
            if nrm > 0:
                p = p / nrm
            self.phi.append(p)

        # -- the dipole matrix elements <i|z|j> (the Rabi couplings; they
        # encode the dipole selection rules -- forbidden ones are ~0) --
        self.dipole_mat = {}
        for i in range(len(self.phi)):
            for j in range(len(self.phi)):
                v = B.asnumpy(self.phi[i] * self.dip * self.phi[j])
                self.dipole_mat[(i, j)] = float(v.sum())

        # -- the local photon-amplitude weight (dipole approximation) --
        flw = self.local_weight_width
        w = B.exp(-self.R ** 2 / (2.0 * flw * flw))
        self.F_local_weight = w / float(B.asnumpy(w.sum()))

        # -- the sponge (a 3-D quadratic absorbing taper) --
        # 1.0 in the bulk, -> 0 within ``sponge_width`` of every face: it
        # kills the FDTD reflections (the young _build_sponge pattern).
        sw = float(self.sponge_width)
        ix = B.arange(self.Lx, dtype=B.FIELD)[:, None, None]
        iy = B.arange(self.Ly, dtype=B.FIELD)[None, :, None]
        iz = B.arange(self.Lz, dtype=B.FIELD)[None, None, :]
        d3 = B.minimum(B.minimum(ix, self.Lx - 1.0 - ix),
                       B.minimum(iy, self.Ly - 1.0 - iy))
        d3 = B.minimum(d3, B.minimum(iz, self.Lz - 1.0 - iz))
        q = B.clip(1.0 - d3 / sw, 0.0, 1.0)
        self.sponge = 1.0 - q * q  # = 1 in the bulk, -> 0 at the edge

        # -- the photon fields (the CL10 FDTD 3-D fields) --
        shp = (self.Lx, self.Ly, self.Lz)
        self.F = B.zeros(shp, dtype=B.FIELD)
        self.F_prev = B.zeros(shp, dtype=B.FIELD)
        self.F_out = B.zeros(shp, dtype=B.FIELD)
        self.F_out_prev = B.zeros(shp, dtype=B.FIELD)

        # -- the atom state (the complex Rabi amplitudes, rotating frame) --
        self.n_levels = len(HYDROGEN_LEVELS)
        self.a = [1.0 + 0.0j] + [0.0 + 0.0j] * (self.n_levels - 1)
        self.omega_ph = 0.0   # the incoming photon's carrier frequency
        self.t = 0
        self.D = 0.0          # the transition dipole (the radiating observable)
        self._d_hist = []
        self._ext = 0.0
        self._ext_drive_fn = None

        self._clear_bookkeeping()

    # ------------------------------------------------------------------
    # state preparation
    # ------------------------------------------------------------------
    def _clear_bookkeeping(self):
        self.t = 0
        self.D = 0.0
        self._d_hist = []
        self.omega_ph = 0.0
        self._ext = 0.0
        self._ext_drive_fn = None

    def _clear_photons(self):
        shp = (self.Lx, self.Ly, self.Lz)
        self.F = B.zeros(shp, dtype=B.FIELD)
        self.F_prev = B.zeros(shp, dtype=B.FIELD)
        self.F_out = B.zeros(shp, dtype=B.FIELD)
        self.F_out_prev = B.zeros(shp, dtype=B.FIELD)

    def set_level(self, k):
        """Reset the atom to the **pure** k-th level (the k-th hydrogen
        state), with the ground (0) the default."""
        self.a = [0.0 + 0.0j] * self.n_levels
        if 0 <= k < self.n_levels:
            self.a[k] = 1.0 + 0.0j
        self._clear_photons()
        self._clear_bookkeeping()

    def set_ground(self):
        """The atom in its ground state (the 1s)."""
        self.set_level(0)

    # ------------------------------------------------------------------
    # one QCA step of the coupled system
    # ------------------------------------------------------------------
    def step(self):
        r2 = self.r * self.r
        E = B._E

        # -- the local photon amplitude (the dipole approximation) --
        # the *envelope* of the field at the atom: the RMS of the field
        # over the dipole weight, sqrt(sum(w * F^2)).  This is the carrier
        # amplitude (a plain average of F would cancel the fast oscillation
        # and vanish -- the field is a wave, its *intensity* F^2 carries
        # the envelope).  The Rabi / Fermi-golden drive scales with this.
        F_local = float(B._E.sqrt(E.sum(self.F_local_weight * self.F * self.F)))
        # the Rabi drive acts only while the photon is in the box; once it
        # has left we stop driving (and stop the detuning phase buildup).
        if float(E.max(B._E.abs(self.F))) < 1e-4:
            self.omega_ph = 0.0

        # -- absorption (the Fermi golden rule) --
        # the photon drives the dipole-allowed transitions at a rate
        # proportional to the field squared and a Lorentzian in the
        # detuning (the natural line shape).  A photon *at* a gap is
        # absorbed (the atom is excited); one *between* the gaps is not
        # (it passes through).  This is the Rabi / Fermi-golden-rule
        # selection rule, free of the Rabi cascade / power-broadening.
        self._absorb(F_local)

        # -- spontaneous emission (Einstein A) --
        # the excited levels that have a dipole channel to a lower level
        # decay into it at the rate gamma_sp; the emitted photon is a
        # dipole radiating at the *gap* frequency (a damped cosine -- the
        # natural line).  This is the vacuum-driven Einstein-A process:
        # a *bare* excited eigenstate decays (its transition dipole to the
        # lower level radiates), which a Larmor <psi|z|psi> source would
        # miss (the diagonal dipole of an eigenstate is exactly zero).
        S, self.D = self._emission_source()
        # -- the incoming photon (free; it passes the atom and leaves) --
        F_new = (
            2.0 * self.F - self.F_prev + r2 * self._lap(self.F)
        ) * self.sponge
        # -- the emitted photon (free + the dipole source; no re-coupling) --
        F_out_new = (
            2.0 * self.F_out - self.F_out_prev
            + r2 * self._lap(self.F_out) + S
        ) * self.sponge

        # -- commit --
        self.F_prev = self.F
        self.F = F_new
        self.F_out_prev = self.F_out
        self.F_out = F_out_new
        self.t += 1
        return self.t

    def _absorb(self, F_local):
        """The **Fermi golden rule** absorption: the photon drives the
        dipole-allowed transitions at a rate proportional to the field
        squared and a **Lorentzian** in the detuning (the natural line
        shape, HWHM = ``gamma_line``).  Population flows from the lower
        to the upper level (absorption) and back (stimulated emission);
        for a ground-state atom it is pure absorption.  A photon *at* a
        gap (delta ~ 0) is absorbed; one *between* the gaps (|delta| >>
        gamma_line) is not -- the resonance selection rule."""
        if self.omega_ph <= 0.0 or abs(F_local) < 1e-4:
            return
        for i in range(self.n_levels):
            for j in range(self.n_levels):
                if self.omega[j] <= self.omega[i]:
                    continue  # j must be the upper level
                d_ij = self.dipole_mat[(i, j)]
                if abs(d_ij) < 1e-6:
                    continue  # dipole-forbidden
                gap = self.omega[j] - self.omega[i]
                delta = self.omega_ph - gap
                # the Lorentzian line shape (HWHM = gamma_line)
                gl = self.gamma_line
                L = (0.5 * gl) / (delta * delta + 0.25 * gl * gl)
                rate = self.A * (F_local * F_local) * L
                P_i = abs(self.a[i]) ** 2
                P_j = abs(self.a[j]) ** 2
                # the net transfer (absorption i->j minus stimulated j->i)
                dP = rate * (P_i - P_j)
                if abs(dP) < 1e-14:
                    continue
                self._transfer_population(i, j, dP)

    def _transfer_population(self, i, j, dP):
        """Transfer population ``dP`` from level ``i`` into ``j``
        (``dP > 0``) or from ``j`` into ``i`` (``dP < 0``), keeping the
        total norm and the levels' phases."""
        if dP > 0.0:
            src, dst = i, j
        else:
            src, dst = j, i
            dP = -dP
        P_src = abs(self.a[src]) ** 2
        if P_src < dP:
            dP = P_src  # cannot transfer more than is present
        if dP < 1e-14:
            return
        P_src_new = P_src - dP
        P_dst = abs(self.a[dst]) ** 2
        P_dst_new = P_dst + dP
        if P_src > 1e-14:
            self.a[src] *= math.sqrt(max(0.0, P_src_new) / P_src)
        else:
            self.a[src] = 0.0 + 0.0j
        if P_dst > 1e-14:
            self.a[dst] *= math.sqrt(P_dst_new / P_dst)
        else:
            self.a[dst] = math.sqrt(P_dst_new) + 0.0j

    def _emission_source(self):
        """Spontaneous emission (Einstein A): for each excited level that
        has a dipole channel to a *lower* level, (a) source the emitted
        photon into ``F_out`` (a dipole radiating at the gap frequency,
        amplitude ~ sqrt(population) * dipole -- a damped cosine, the
        natural line) and (b) transfer that level's population into the
        lower level at the rate ``gamma_sp``.
        Returns ``(S, D)``: the emitted source field and the (coarse)
        radiating dipole for the HUD."""
        S = B.zeros_like(self.F)
        D = 0.0
        for j in range(1, self.n_levels):
            p_j = abs(self.a[j]) ** 2
            if p_j < 1e-6:
                continue
            # the strongest dipole channel to a *lower* level
            best_i, best_d = None, 0.0
            for i in range(self.n_levels):
                if self.omega[i] >= self.omega[j] - 1e-9:
                    continue
                d = self.dipole_mat[(i, j)]
                if abs(d) > abs(best_d):
                    best_i, best_d = i, d
            if best_i is None or abs(best_d) < 1e-6:
                continue  # no dipole channel: this state is (meta)stable
            omega_ji = self.omega[j] - self.omega[best_i]
            D += abs(best_d) * math.sqrt(p_j)
            # the emitted dipole: the pattern * a damped cosine at the gap
            amp = self.h * best_d * math.sqrt(p_j)
            S = S + self.dip_emit * amp * math.cos(omega_ji * self.t)
            # the decay: population j -> lower level, the norm is kept
            self._decay_population(j, best_i, self.gamma_sp)
        return S, D

    def _decay_population(self, j, i, gamma):
        """Transfer population from level ``j`` into level ``i`` at the
        rate ``gamma`` (per step), keeping the total norm and the phases
        (the Einstein-A relaxation)."""
        p_j = abs(self.a[j]) ** 2
        if p_j < 1e-12:
            return
        p_i = abs(self.a[i]) ** 2
        dP = gamma * p_j
        p_j_new = p_j - dP
        p_i_new = p_i + dP
        self.a[j] *= math.sqrt(max(0.0, p_j_new) / p_j)
        if p_i > 1e-12:
            self.a[i] *= math.sqrt(p_i_new / p_i)
        else:
            self.a[i] = math.sqrt(p_i_new)

    # ------------------------------------------------------------------
    # the photon: fire one one-directional packet from the boundary
    # ------------------------------------------------------------------
    def fire_photon(self, omega, direction="x", amp=1.0, width=None):
        """Fire one photon packet of (lattice) frequency ``omega`` toward
        the atom, from the boundary in ``direction``.  The carrier
        wavenumber ``k0`` is the inverse lattice dispersion
        ``omega = 2 asin( r sin(k0/2) )``; the packet is a Gaussian
        envelope times a travelling carrier, launched one-directionally.
        Returns ``(omega, k0)``."""
        omega = min(float(omega), 0.99 * self._omega_max())
        k0 = 2.0 * math.asin(min(0.999, math.sin(omega / 2.0) / self.r))
        if width is None:
            width = self.photon_width if self.photon_width > 0 else max(6.0, 3.0 / k0)

        axis = direction[0]
        plus = not direction.startswith("-")
        sign = 1.0 if plus else -1.0
        L = {"x": self.Lx, "y": self.Ly, "z": self.Lz}[axis]
        s0 = (self.sponge_width + 2) if plus else (L - self.sponge_width - 3)
        ix = B.arange(self.Lx, dtype=B.FIELD)[:, None, None]
        iy = B.arange(self.Ly, dtype=B.FIELD)[None, :, None]
        iz = B.arange(self.Lz, dtype=B.FIELD)[None, None, :]
        cx, cy, cz = self.cx, self.cy, self.cz
        if axis == "x":
            T = k0 * (ix - s0) * sign
            G = amp * B.exp(-((ix - s0) ** 2 + (iy - cy) ** 2
                              + (iz - cz) ** 2) / (2.0 * width ** 2))
        elif axis == "y":
            T = k0 * (iy - s0) * sign
            G = amp * B.exp(-((ix - cx) ** 2 + (iy - s0) ** 2
                              + (iz - cz) ** 2) / (2.0 * width ** 2))
        else:
            T = k0 * (iz - s0) * sign
            G = amp * B.exp(-((ix - cx) ** 2 + (iy - cy) ** 2
                              + (iz - cz) ** 2) / (2.0 * width ** 2))
        # one-directional launch (the standard two-time-level initialization)
        self.F = G * B.cos(T)
        self.F_prev = G * B.cos(T + omega)
        self.omega_ph = float(omega)
        self._last_photon_omega = float(omega)
        self._last_photon_k = float(k0)
        return float(omega), float(k0)

    def fire_transition(self, i, j, amp=1.0, direction="x", width=None):
        """Fire a photon tuned to the **gap** between level ``i`` (lower)
        and level ``j`` (upper) -- the resonant frequency of that
        transition.  The photon is absorbed (the atom is driven ``i -> j``)
        only when this detuning is small.  Returns ``(omega, k0)``."""
        omega = abs(self.omega[j] - self.omega[i])
        return self.fire_photon(omega, direction=direction, amp=amp, width=width)

    def _omega_max(self):
        # the band edge of the free field (Brillouin corner, 6 neighbours)
        return 2.0 * math.asin(min(1.0, self.r * math.sqrt(3.0)))

    def _lap(self, f):
        """The 6-neighbour discrete Laplacian (zero outside the box)."""
        out = B.zeros_like(f)
        out[1:, :, :] += f[:-1, :, :]
        out[:-1, :, :] += f[1:, :, :]
        out[:, 1:, :] += f[:, :-1, :]
        out[:, :-1, :] += f[:, 1:, :]
        out[:, :, 1:] += f[:, :, :-1]
        out[:, :, :-1] += f[:, :, 1:]
        out -= 6.0 * f
        return out

    # ------------------------------------------------------------------
    # observables
    # ------------------------------------------------------------------
    def populations(self):
        """The level populations (the |a_k|**2, the Rabi probabilities)."""
        return [float(abs(a) ** 2) for a in self.a]

    def excited_fraction(self):
        """The fraction of the atom's population in the excited levels
        (1 = fully out of the ground state)."""
        p = self.populations()
        return 1.0 - min(1.0, p[0])

    def wavefunction(self):
        """The **real** 3-D electron field for rendering:
        ``Re[ sum_k a_k e^{-i omega_k t} phi_k ]``."""
        psi = B.zeros((self.Lx, self.Ly, self.Lz), dtype=B.FIELD)
        for k in range(self.n_levels):
            a = self.a[k]
            if abs(a) < 1e-12:
                continue
            ph = -self.omega[k] * self.t
            real = (a * complex(math.cos(ph), math.sin(ph))).real
            psi = psi + real * self.phi[k]
        return psi

    def electron_density(self):
        """The 3-D electron probability density ``|psi|**2`` on the active
        backend (``psi`` the real superposition of the level
        wavefunctions).  The 1s is a sphere, the 2p a two-lobe dumbbell,
        the 3d a doughnut -- the shapes that deform under excitation."""
        psi = self.wavefunction()
        return psi * psi

    def photon_density(self):
        """The 3-D photon intensity ``|F|**2 + |F_out|**2`` on the active
        backend (the incoming packet + the emitted radiation)."""
        return self.F * self.F + self.F_out * self.F_out

    def to_host(self):
        """(psi, F, F_out) on the host (for rendering / dumping)."""
        return (
            B.asnumpy(self.wavefunction()),
            B.asnumpy(self.F),
            B.asnumpy(self.F_out),
        )

    # ------------------------------------------------------------------
    # diagnostics
    # ------------------------------------------------------------------
    def photon_energy(self):
        E = B._E
        return float(
            E.sum(self.F * self.F + (self.F - self.F_prev) ** 2)
            + E.sum(self.F_out * self.F_out
                    + (self.F_out - self.F_out_prev) ** 2)
        )

    def atom_energy(self):
        """The atom's internal energy (the population-weighted level
        energy, in the Rydberg units)."""
        p = self.populations()
        n = [lv[1] for lv in HYDROGEN_LEVELS]
        E = [-RYDBERG_EV / nn ** 2 for nn in n]
        return float(sum(p[k] * E[k] for k in range(len(p))))
