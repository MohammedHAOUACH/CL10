#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the publication figure of the 3-D hydrogen atom (Fig. 26).

**Why this script exists.**  The previous Fig. 26 was a **screenshot of the live
VisPy window**: it carried the application title, a `[PAUSE]` flag, the keyboard
shortcut legend and the live HUD telemetry, and it had no axes, no spatial scale
and no labelled colour bar.  That is a screen grab, not a figure.

This script **computes** the figure: it runs the `simulation.hydrogen` engine (the
real Rydberg levels coupled to the CL10 free massless 3-D field) and plots the
fields the engine actually produces — `|psi|^2` and `|F|^2`, `|F_out|^2` — at
instants chosen by the **physics**, with axes in units of the Bohr radius `a0`, a
labelled colour bar, and no overlay text.  Nothing is drawn by hand and nothing is
taken from a window.

The panels (each the **mid-plane cut** through the nucleus, the plane `x = Lx/2`,
so the horizontal axis is the dipole axis `z` and the vertical axis is `x`):

  a) the ground state `1s` — the reference sphere;
  b) the incoming CL10 photon packet `|F|^2` at the moment it is strongest at the
     atom;
  c) the cloud **after** the resonant absorption — the deformed `1s`+`2p`
     superposition (the dumbbell).  The instant is chosen where the **cut of the
     density itself** is largest after the arrival, not where the excited fraction
     is largest: `wavefunction()` returns the **real** part of the superposition,
     so at a node of `Re[a_k exp(-i omega_k t)]` the density momentarily collapses
     everywhere, and a panel taken there renders blank (this was the first bug of
     this figure);
  d) the emitted field `|F_out|^2` at the moment its cut is strongest (the
     Einstein-A radiation).

Each panel is **normalised to its own peak** (and the colour bar says so): the `2p`
peak is orders of magnitude below the `1s` peak, so a shared scale renders the
deformed cloud invisible (the second bug of this figure).

Run::

    python3 make_figure_hydrogen.py                    # -> doc/figures/fig26_hydrogen.png
    python3 make_figure_hydrogen.py --grid 160 --out /tmp/p.png    # prototype
    CL10_BACKEND=numpy python3 make_figure_hydrogen.py             # force CPU

The grid, the selected instants and the populations are written next to the PNG as
`fig26_hydrogen.meta.json`, so the caption can state exactly what was computed.
"""

import os
import sys
import json

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulation.hydrogen import runner, qca as QCA, backend as B

OUT_PNG = os.path.join(ROOT, "doc", "figures", "fig26_hydrogen.png")
OUT_META = os.path.join(ROOT, "doc", "figures", "fig26_hydrogen.meta.json")


def midplane(field_host):
    """The mid-plane cut through the nucleus: the plane `x = Lx/2`, plotted with
    the dipole axis `z` horizontal and `x` vertical."""
    sl = field_host[field_host.shape[0] // 2, :, :]   # (Ly, Lz)
    return np.ascontiguousarray(sl.T)                 # (Lz, Ly)


def axis_in_a0(arr2d, a0_cells):
    """The coordinate axis of a cut, in units of the Bohr radius."""
    n = arr2d.shape[0]
    return (np.arange(n) - (n - 1) / 2.0) / a0_cells


def render_panel(fig, ax, arr2d, cmap, label, a0_cells, bar_label):
    """Plot a cut normalised to its OWN peak, on a log scale, with real axes."""
    peak = float(arr2d.max())
    if peak <= 0.0:
        peak = 1.0
    # log2(1 + I/Ipeak): 0 at the background, 1 at this panel's own peak
    disp = np.log1p(np.maximum(arr2d, 0.0) / peak) / np.log(2.0)
    ax2 = axis_in_a0(arr2d, a0_cells)
    lo, hi = float(ax2[0]), float(ax2[-1])
    im = ax.imshow(disp, origin="lower", extent=(lo, hi, lo, hi), cmap=cmap,
                   aspect="equal", vmin=0.0, vmax=1.0)
    ax.set_xlabel(r"$z\;/\;a_0$", fontsize=9)
    ax.set_ylabel(r"$x\;/\;a_0$", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_title(label, fontsize=9.5, loc="left", pad=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.05)
    cb.set_label(bar_label, fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    return im


def main(argv=()):
    global OUT_PNG, OUT_META
    cfg = runner.load_config()
    # --out / --grid let the figure be prototyped on a small grid before the
    # production run, without touching the committed PNG.
    if "--out" in argv:
        OUT_PNG = os.path.abspath(argv[argv.index("--out") + 1])
        OUT_META = os.path.splitext(OUT_PNG)[0] + ".meta.json"
    if "--grid" in argv:
        n = int(argv[argv.index("--grid") + 1])
        cfg["grid"] = [n, n, n]
    B.set_field_dtype(cfg.get("field_dtype", "float32"))

    a0 = float(cfg["bohr_radius"])
    L = int(cfg["grid"][0])
    amp = float(cfg["photon_amp"])
    print(f"[fig26] backend={B.backend_name()}  grid={cfg['grid']}  a0={a0} cells",
          flush=True)

    sim = QCA.AtomQCA(cfg)
    sim.set_ground()
    gap12 = abs(sim.omega[1] - sim.omega[0])
    # the real flight distance: the packet is launched just inside the sponge
    launch = int(cfg["sponge_width"]) + 2
    flight = int(((L // 2) - launch) / float(cfg["r"]))
    n_steps = flight + 500
    print(f"[fig26] 1s->2p gap = {sim.gaps_eV[(0, 1)]:.3f} eV; the packet reaches the "
          f"atom near step {flight}; running {n_steps} steps", flush=True)

    # -- (a) the ground state, before anything is fired --
    cut_ground = midplane(B.asnumpy(sim.electron_density()))

    # -- fire the resonant photon and run ONCE, keeping only the three instants
    #    the physics selects (each snapshot is kept only while it is the record
    #    holder, so the memory stays bounded) --
    sim.fire_photon(gap12, direction="x", amp=amp)
    best = {"incoming": (-1.0, 0, None), "absorbed": (-1.0, 0, None),
            "emitted": (-1.0, 0, None)}
    t = 0
    while t < n_steps:
        sim.step(); t += 1
        f_at_atom = float(B.asnumpy(sim.F_local_weight * sim.F * sim.F).sum())
        exc = sim.excited_fraction()

        # (b) the packet at its strongest on the atom
        if f_at_atom > best["incoming"][0]:
            best["incoming"] = (f_at_atom, t, (midplane(B.asnumpy(sim.F * sim.F)),
                                               list(sim.populations()), exc))
        # (c) the cloud at its densest AFTER the arrival — never at a node of the
        #     real part of the superposition, which would render a blank panel
        if t > flight:
            cut = midplane(B.asnumpy(sim.electron_density()))
            if float(cut.max()) > best["absorbed"][0]:
                best["absorbed"] = (float(cut.max()), t, (cut,
                                                          list(sim.populations()), exc))
        # (d) the emitted field at its strongest cut
        if t > flight:
            cf = midplane(B.asnumpy(sim.F_out * sim.F_out))
            if float(cf.max()) > best["emitted"][0]:
                best["emitted"] = (float(cf.max()), t, (cf,
                                                        list(sim.populations()), exc))

        if t in (flight, flight + 100, flight + 250, flight + 400):
            print(f"[fig26]  t={t:5d}  excited={exc:.3f}  |F|^2_at_atom={f_at_atom:.4f}  "
                  f"max|F_out|^2={float(B.asnumpy(sim.F_out * sim.F_out).max()):.3e}",
                  flush=True)

    cut_in, pops_in, exc_in = best["incoming"][2]
    cut_abs, pops_abs, exc_abs = best["absorbed"][2]
    cut_out, pops_em, exc_em = best["emitted"][2]
    t_in, t_abs, t_em = best["incoming"][1], best["absorbed"][1], best["emitted"][1]
    print(f"[fig26] selected instants: incoming t={t_in} (|F|^2 at the atom "
          f"{best['incoming'][0]:.4f}), absorbed t={t_abs} (excited {exc_abs:.3f}), "
          f"emitted t={t_em} (max |F_out|^2 {best['emitted'][0]:.3e})", flush=True)

    # -- the figure --
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 9.8))
    a, b, c, d = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    psi_bar = r"$\log_2\!\left(1+|\psi|^2/|\psi|^2_{\max}\right)$"
    render_panel(fig, a, cut_ground, "viridis",
                 r"(a)  the ground state $1s$:  $|\psi|^2$", a0, psi_bar)
    render_panel(fig, b, cut_in, "inferno",
                 "(b)  the incoming CL10 photon:  $|F|^2$", a0,
                 r"$\log_2\!\left(1+|F|^2/|F|^2_{\max}\right)$")
    render_panel(fig, c, cut_abs, "viridis",
                 r"(c)  after the resonant absorption:  $|\psi|^2$  (the deformed cloud)",
                 a0, psi_bar)
    render_panel(fig, d, cut_out, "inferno",
                 r"(d)  the emitted field:  $|F_{\rm out}|^2$  (the Einstein-A radiation)",
                 a0, r"$\log_2\!\left(1+|F_{\rm out}|^2/|F_{\rm out}|^2_{\max}\right)$")

    fig.suptitle("the 3-D hydrogen atom driven by the CL10 photon\n"
                 "(mid-plane cut through the nucleus; each panel normalised to its own peak)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    plt.close(fig)
    print("[fig26] wrote", OUT_PNG, flush=True)

    meta = {
        "backend": B.backend_name(),
        "grid": list(cfg["grid"]),
        "bohr_radius_cells": a0,
        "cfl_r": float(cfg["r"]),
        "gap_1s2p_eV": float(sim.gaps_eV[(0, 1)]),
        "photon_omega_rad_per_step": float(gap12),
        "launch_cell": launch,
        "flight_steps": flight,
        "steps_run": n_steps,
        "instant_incoming_step": t_in,
        "instant_absorbed_step": t_abs,
        "instant_emitted_step": t_em,
        "field_at_atom_at_incoming": best["incoming"][0],
        "populations_at_incoming": {sim.labels[k]: pops_in[k] for k in range(len(pops_in))},
        "excited_fraction_at_incoming": exc_in,
        "populations_at_absorption": {sim.labels[k]: pops_abs[k] for k in range(len(pops_abs))},
        "excited_fraction_at_absorption": exc_abs,
        "populations_at_emission": {sim.labels[k]: pops_em[k] for k in range(len(pops_em))},
        "excited_fraction_at_emission": exc_em,
        "coupling_absorb": float(cfg["coupling_absorb"]),
        "coupling_emit": float(cfg["coupling_emit"]),
        "reaction_gamma_max": float(cfg["reaction_gamma_max"]),
        "photon_width_cells": float(cfg["photon_width"]),
        "photon_amp": amp,
        "note": ("computed by make_figure_hydrogen.py from the simulation.hydrogen engine; "
                 "no window screenshot, no overlay text; the atom is the real hydrogen "
                 "(an input of atomic physics), the photon is the CL10 free massless field"),
    }
    with open(OUT_META, "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print("[fig26] wrote", OUT_META, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
