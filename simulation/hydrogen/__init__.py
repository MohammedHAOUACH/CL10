# -*- coding: utf-8 -*-
"""CL10 / simulation/hydrogen -- the 3-D hydrogen atom as a CL10 QCA.

The **atom** is the **semiclassical Rabi / Fermi-golden-rule** model: its
levels and wavefunctions are the **real hydrogen** eigenstates (the Rydberg
energies ``E_n = -13.606/n**2`` eV and the ``R_nl(r) Y_lm`` shapes, 0-parameter
-- the known solutions of the hydrogen problem, not fitted). The **photon**
is the free CL10 massless 3-D field (the same 6-neighbour leapfrog as
``simulation/young``), in two parts: the **incoming** packet (``F``) and the
**emitted** radiation (``F_out``).

The **excitation** is **not imposed**: a photon packet is fired at the atom
and the **dynamics decides** whether it is absorbed -- the (Fermi-golden-rule)
dipole drive transfers population **only** when the packet frequency matches
a transition (gap) of the atom. The excited atom **decays** by radiating a
photon (the Einstein-A process: the transition dipole sources ``F_out`` at
the gap frequency). No Schrödinger integrator is solved numerically: the
levels are the exact eigenstates and each step is a discrete lattice update.
"""

__all__ = ["qca", "backend", "runner", "scene", "vispy_h"]
