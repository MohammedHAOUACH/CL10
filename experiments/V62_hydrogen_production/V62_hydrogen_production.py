#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — V62: the production 3-D hydrogen experiment on the GPU (regenerates Fig. 26).

**What this experiment is.**  V61 runs the *same model* on a reduced ``120^3``
grid as a fast check inside the validation suite.  This is the **production**
version: the full **296^3** simulation on the **GPU** (CuPy, the card pinned in
``simulation/hydrogen/config.yaml`` -- ``gpu: 1``, the free card), running the
**complete photon script** of the hydrogen experiment:

    Lyman-alpha (1s->2p)   Lyman-alpha again   1s->3p   2p->3d   off-resonant control

the same sequence the headless video plays (``simulation/hydrogen/runner.py``).
Every claim of Fig. 26 is measured at production resolution, and **the figure
itself is regenerated here** (the four mid-plane panels + the meta JSON), so the
published figure is a *product of the experiment*, like the other figures of the
suite, not a screenshot and not a top-level script.

**What is measured (the content of Fig. 26, on the production grid):**

1. **the baseline** — the atom alone (no photon) does not excite itself.
2. **the resonant absorption** — a photon fired exactly at a gap is absorbed:
   the excited fraction peaks in each resonant phase (Lyman-alpha twice,
   1s->3p, 2p->3d), and the electron cloud **deforms** — its density quadrupole
   goes from the spherical 1s (anisotropy ~ 0) to the two-lobe 2p dumbbell
   (the lobes along the dipole axis stretch the z direction, anisotropy > 0).
3. **the off-resonant control** — a photon fired between the gaps is mostly
   **not** absorbed: the excited fraction stays a small fraction of the
   resonant peak (not zero — the Fermi-golden-rule rate carries a Lorentzian
   tail; real physics, stated, not hidden).
4. **the Einstein-A radiation** — the excited atom radiates: the emitted field
   ``F_out`` is non-zero at the gap frequency, and the excited population
   decays back toward the 1s after the photon passes.
5. **the stability** — the atomic norm is conserved and the FDTD fields stay
   bounded (the UV-finite property, at the production CFL number).

**Honest provenance (what this is NOT).**  Same as V61: the **atom is an
input** — its levels and shapes are the **real hydrogen** eigenstates (Rydberg
energies, ``R_nl(r) Y_lm``), **not** derived from the Cl(10) root.  What the
root supplies is the **photon** — the free massless 3-D field, its
propagation, its resonance with the atom, and its emission.  So this
experiment **demonstrates the photon sector**; it is **not** a 0-parameter
derivation of hydrogen.  The three coupling knobs (``A``, ``h``, ``gamma``)
are physical parameters of the atom, **not** fitted to any measured number.

**Run** (from the repo root):

    python3 experiments/V62_hydrogen_production/V62_hydrogen_production.py

-> the verdict ``experiments/V62_hydrogen_production/results/v62_hydrogen_production.json``
+ the figure ``doc/figures/fig26_hydrogen.png`` and the run metadata
``experiments/V62_hydrogen_production/results/fig26_hydrogen.meta.json``.
"""

import os
import sys
import json
import time

import numpy as np
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

# Pin the GPU *before* the CuPy import (mirrors
# simulation/hydrogen/__main__.py): the configured card (``gpu``, default 1 --
# the free card) rather than the first visible one (which may be occupied by
# another process) or the CPU.  An explicit CUDA_VISIBLE_DEVICES set by the
# caller wins.
HYDRO_DIR = os.path.join(ROOT, "simulation", "hydrogen")
with open(os.path.join(HYDRO_DIR, "config.yaml")) as _f:
    _hcfg = yaml.safe_load(_f)
os.environ.setdefault("CUDA_VISIBLE_DEVICES", str(_hcfg.get("gpu", 1)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulation.hydrogen import runner, qca as QCA, backend as B

OUT_PNG = os.path.join(ROOT, "doc", "figures", "fig26_hydrogen.png")
OUT_META = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results",
                        "fig26_hydrogen.meta.json")


# ---------------------------------------------------------------------------
# the figure helpers (the mid-plane cut through the nucleus, the plane
# x = Lx/2, the dipole axis z horizontal)
# ---------------------------------------------------------------------------
def midplane(field_host):
    """The mid-plane cut through the nucleus (the plane ``x = Lx/2``), with
    the dipole axis ``z`` horizontal and ``x`` vertical."""
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


# ---------------------------------------------------------------------------
# the density quadrupole (the shape of the cloud, measured, not assumed)
# ---------------------------------------------------------------------------
def quadrupole_anisotropy(sim, rho):
    """The electron-density quadrupole anisotropy about the dipole axis ``z``:
    ``(Qzz - (Qxx + Qyy)/2) / tr(Q)``.  It is ~0 for the **spherical** 1s and
    > 0 for the **2p dumbbell** (the two lobes along z stretch the z direction).
    The shape of the cloud is a property of the lattice wavefunctions, so this
    is the quantitative statement of 'the excitation deforms the cloud'."""
    ix = B.arange(sim.Lx, dtype=B.FIELD)[:, None, None] - sim.hx
    iy = B.arange(sim.Ly, dtype=B.FIELD)[None, :, None] - sim.hy
    iz = B.arange(sim.Lz, dtype=B.FIELD)[None, None, :] - sim.hz
    qxx = float(B.asnumpy((rho * ix * ix).sum()))
    qyy = float(B.asnumpy((rho * iy * iy).sum()))
    qzz = float(B.asnumpy((rho * iz * iz).sum()))
    tr = qxx + qyy + qzz
    if tr <= 0.0:
        return 0.0, (qxx, qyy, qzz)
    return (qzz - 0.5 * (qxx + qyy)) / tr, (qxx, qyy, qzz)


# ---------------------------------------------------------------------------
# one production run: the full photon script, the per-phase records, and the
# physics-selected instants of the Lyman-alpha phase (for the figure)
# ---------------------------------------------------------------------------
def run_production(cfg):
    sim = QCA.AtomQCA(cfg)
    sim.set_ground()
    L = int(cfg["grid"][0])
    r = float(cfg["r"])
    launch = int(cfg["sponge_width"]) + 2
    flight = int(((L // 2) - launch) / r)
    n_steps = int(cfg["n_steps"])
    amp = float(cfg["photon_amp"])

    script = runner.build_script(cfg)
    print("[v62] backend=%s  gpu=%s  grid=%s  dtype=%s"
          % (B.backend_name(), B.gpu_info(), cfg["grid"],
             cfg.get("field_dtype", "float32")), flush=True)
    print("[v62] photon script (t, (i,j)): %s" % script, flush=True)
    print("[v62] 1s->2p gap = %.3f eV; the packet reaches the atom ~%d steps "
          "after a shot; running %d steps"
          % (sim.gaps_eV[(0, 1)], flight, n_steps), flush=True)

    # (a) the ground state, before anything is fired
    cut_ground = midplane(B.asnumpy(sim.electron_density()))
    anis_1s, moments_1s = quadrupole_anisotropy(sim, sim.electron_density())

    # the record holders for the fig26 panels (the first Lyman-alpha phase)
    best = {"incoming": (-1.0, 0, None), "absorbed": (0.0, 0, None),
            "emitted": (-1.0, 0, None)}

    phases = []
    cur = None
    baseline_max_excited = 0.0
    norm_min, norm_max = 1.0, 0.0
    field_max = 0.0
    t = 0
    fired = 0

    while t < n_steps:
        sim.step()
        t += 1

        # fire the next scripted photon (recording the outgoing state of the
        # previous phase right before the reset)
        if fired < len(script) and t >= script[fired][0]:
            if cur is not None:
                cur["t_end"] = t
                cur["final_excited"] = sim.excited_fraction()
                cur["final_populations"] = {
                    sim.labels[k]: v for k, v in enumerate(sim.populations())
                }
            i, j = script[fired][1]
            if i < 0:
                sim.set_ground()
                sim.fire_photon(0.13, direction="x", amp=amp)
                tag = "off-resonant control (0.13 rad/step, between the gaps)"
                omega = 0.13
            else:
                sim.set_level(i)
                omega = sim.fire_transition(i, j, amp=amp)[0]
                tag = "resonant %s->%s" % (sim.labels[i], sim.labels[j])
            cur = {
                "tag": tag,
                "upper": None,
                "omega_rad_per_step": float(omega),
                "t_fire": t,
                "peak_excited": 0.0,
                "t_peak": t,
                "peak_upper_pop": 0.0,
                "max_F_at_atom": 0.0,
                "max_F_out": 0.0,
            }
            if i >= 0:
                cur["upper"] = j
            phases.append(cur)
            fired += 1

        # the bookkeeping (sampled every 10 steps: the reductions are GPU
        # syncs, the fields are smooth on this timescale)
        if t % 10 == 0:
            nrm = float(sum(abs(a) ** 2 for a in sim.a))
            norm_min = min(norm_min, nrm)
            norm_max = max(norm_max, nrm)
            field_max = max(
                field_max,
                float(B.asnumpy((sim.F * sim.F).max())) ** 0.5,
                float(B.asnumpy((sim.F_out * sim.F_out).max())) ** 0.5)

        if cur is None:
            baseline_max_excited = max(
                baseline_max_excited, sim.excited_fraction())
            continue

        exc = sim.excited_fraction()
        f_at_atom = float(B.asnumpy((sim.F_local_weight * sim.F * sim.F).sum()))
        if exc > cur["peak_excited"]:
            cur["peak_excited"] = exc
            cur["t_peak"] = t
            if cur is phases[0]:
                # the dumbbell instant: the quadrupole of the *mixed*
                # cloud density at the maximum-excitation instant
                cur["quadrupole_at_peak"] = quadrupole_anisotropy(
                    sim, sim.electron_density())[0]
        if cur["upper"] is not None:
            cur["peak_upper_pop"] = max(
                cur["peak_upper_pop"], sim.populations()[cur["upper"]])
        cur["max_F_at_atom"] = max(cur["max_F_at_atom"], f_at_atom)
        cur["max_F_out"] = max(
            cur["max_F_out"], float(B.asnumpy((sim.F_out * sim.F_out).max())))

        # the first Lyman-alpha phase carries the fig26 panels: keep only the
        # record holders (the density cut peaks between the node instants of
        # the real part of the superposition)
        if cur is phases[0] and t > cur["t_fire"] + flight:
            if f_at_atom > best["incoming"][0]:
                best["incoming"] = (f_at_atom, t, (
                    midplane(B.asnumpy(sim.F * sim.F)),
                    list(sim.populations()), exc))
            if t % 2 == 0:
                rho = sim.electron_density()
                anis = quadrupole_anisotropy(sim, rho)[0]
                if anis > best["absorbed"][0]:
                    best["absorbed"] = (anis, t, (
                        midplane(B.asnumpy(rho)), list(sim.populations()), exc))
            if t % 4 == 0:
                cut_out = midplane(B.asnumpy(sim.F_out * sim.F_out))
                if float(cut_out.max()) > best["emitted"][0]:
                    best["emitted"] = (float(cut_out.max()), t, (
                        cut_out, list(sim.populations()), exc))

        if t % 300 == 0:
            print("[v62]  t=%5d  phase=%s  excited=%.3f  |F|^2@atom=%.4f"
                  % (t, cur["tag"].split()[0] + " " + cur["tag"].split()[1],
                     exc, f_at_atom), flush=True)

    if cur is not None:
        cur["t_end"] = n_steps
        cur["final_excited"] = sim.excited_fraction()
        cur["final_populations"] = {
            sim.labels[k]: v for k, v in enumerate(sim.populations())
        }

    # the 1s and 2p *pure* shapes (the deformation reference): the 1s is a
    # sphere (anisotropy ~ 0), the 2p a dumbbell along z (positive).
    anis_1s_pure, _ = quadrupole_anisotropy(sim, sim.phi[0] * sim.phi[0])
    anis_2p_pure, _ = quadrupole_anisotropy(sim, sim.phi[1] * sim.phi[1])

    return {
        "sim": sim,
        "flight": flight,
        "script": script,
        "n_steps": n_steps,
        "cut_ground": cut_ground,
        "best": best,
        "phases": phases,
        "baseline_max_excited": baseline_max_excited,
        "norm_min": norm_min,
        "norm_max": norm_max,
        "field_max": field_max,
        "anis_1s": anis_1s,
        "moments_1s": moments_1s,
        "anis_1s_pure": anis_1s_pure,
        "anis_2p_pure": anis_2p_pure,
        "lyman_quadrupole_at_peak": phases[0].get("quadrupole_at_peak", 0.0),
    }


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------
def checks(r):
    phases = r["phases"]
    resonant = [p for p in phases if p["tag"].startswith("resonant")]
    offres = [p for p in phases if p["tag"].startswith("off-resonant")]
    lyman = phases[0]

    per_phase = []
    for p in phases:
        peak = p["peak_excited"]
        final = p["final_excited"]
        per_phase.append({
            "tag": p["tag"],
            "omega_rad_per_step": p["omega_rad_per_step"],
            "t_fire": p["t_fire"],
            "t_peak": p["t_peak"],
            "t_end": p["t_end"],
            "peak_excited": peak,
            "peak_upper_pop": p["peak_upper_pop"],
            "final_excited": final,
            "max_F_at_atom": p["max_F_at_atom"],
            "max_F_out": p["max_F_out"],
            "final_populations": p["final_populations"],
            "decayed_after_peak": bool(final < peak - 0.02),
        })

    min_resonant_peak = min(p["peak_excited"] for p in resonant)
    offres_peak = offres[0]["peak_excited"]
    c1 = r["baseline_max_excited"] < 1e-6
    c2 = all(p["peak_upper_pop"] > 0.02 for p in resonant)
    c3 = lyman["peak_excited"] > 0.05
    c4 = offres_peak < 0.25 * min_resonant_peak and offres_peak < 0.05
    lyman_quad = r["lyman_quadrupole_at_peak"]
    c5 = (r["anis_2p_pure"] > 0.05
          and r["anis_1s_pure"] < 1e-3
          and lyman_quad > 0.05)
    c6 = (all(p["max_F_out"] > 1e-3 for p in resonant)
          and all(p["final_excited"] < p["peak_excited"] - 0.02 for p in resonant))
    c7 = (abs(r["norm_max"] - 1.0) < 1e-9 and r["norm_min"] > 1.0 - 1e-9
          and r["field_max"] < 1e6)

    return {
        "baseline_no_self_excitation": {
            "max_excited_fraction": r["baseline_max_excited"],
            "pass": bool(c1),
        },
        "resonant_absorption": {
            "phases": per_phase,
            "min_resonant_peak": min_resonant_peak,
            "lyman_alpha_peak": lyman["peak_excited"],
            "every_resonant_phase_absorbs": bool(c2),
            "lyman_alpha_excites_visibly": bool(c3),
        },
        "off_resonant_control": {
            "omega_rad_per_step": offres[0]["omega_rad_per_step"],
            "peak_excited": offres_peak,
            "suppression_vs_min_resonant": (
                min_resonant_peak / max(offres_peak, 1e-12)),
            "mostly_not_absorbed": bool(c4),
            "caveat": ("not zero: the Fermi-golden-rule rate carries a "
                       "Lorentzian tail (real physics, stated, not hidden)"),
        },
        "shape_deformation": {
            "anisotropy_1s_pure": r["anis_1s_pure"],
            "anisotropy_2p_pure": r["anis_2p_pure"],
            "anisotropy_1s_measured": r["anis_1s"],
            "lyman_alpha_peak_anisotropy": lyman_quad,
            "the_excitation_deforms_the_cloud": bool(c5),
        },
        "einstein_a_radiation": {
            "every_resonant_phase_radiates": bool(
                all(p["max_F_out"] > 1e-3 for p in resonant)),
            "the_excited_population_decays": bool(
                all(p["final_excited"] < p["peak_excited"] - 0.02 for p in resonant)),
            "max_F_out_per_phase": {p["tag"]: p["max_F_out"] for p in resonant},
            "pass": bool(c6),
        },
        "stability": {
            "norm_min": r["norm_min"],
            "norm_max": r["norm_max"],
            "field_max": r["field_max"],
            "pass": bool(c7),
        },
    }


# ---------------------------------------------------------------------------
# the figure (the four mid-plane panels of Fig. 26)
# ---------------------------------------------------------------------------
def make_figure(r, cfg):
    a0 = float(cfg["bohr_radius"])
    cut_in, pops_in, exc_in = r["best"]["incoming"][2]
    cut_abs, pops_abs, exc_abs = r["best"]["absorbed"][2]
    cut_out, pops_em, exc_em = r["best"]["emitted"][2]

    fig, axes = plt.subplots(2, 2, figsize=(9.8, 9.8))
    a, b, c, d = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    psi_bar = r"$\log_2\!\left(1+|\psi|^2/|\psi|^2_{\max}\right)$"
    render_panel(fig, a, r["cut_ground"], "viridis",
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

    meta = {
        "backend": B.backend_name(),
        "gpu": B.gpu_info(),
        "grid": list(cfg["grid"]),
        "bohr_radius_cells": a0,
        "cfl_r": float(cfg["r"]),
        "field_dtype": str(cfg.get("field_dtype", "float32")),
        "gap_1s2p_eV": float(r["sim"].gaps_eV[(0, 1)]),
        "photon_script": r["script"],
        "n_steps": r["n_steps"],
        "flight_steps": r["flight"],
        "instant_incoming_step": r["best"]["incoming"][1],
        "field_at_atom_at_incoming": r["best"]["incoming"][0],
        "instant_absorbed_step": r["best"]["absorbed"][1],
        "absorbed_quadrupole_anisotropy": r["best"]["absorbed"][0],
        "populations_at_absorption": {r["sim"].labels[k]: pops_abs[k]
                                      for k in range(len(pops_abs))},
        "excited_fraction_at_absorption": exc_abs,
        "instant_emitted_step": r["best"]["emitted"][1],
        "populations_at_emission": {r["sim"].labels[k]: pops_em[k]
                                    for k in range(len(pops_em))},
        "excited_fraction_at_emission": exc_em,
        "coupling_absorb": float(cfg["coupling_absorb"]),
        "coupling_emit": float(cfg["coupling_emit"]),
        "reaction_gamma_max": float(cfg["reaction_gamma_max"]),
        "photon_width_cells": float(cfg["photon_width"]),
        "photon_amp": float(cfg["photon_amp"]),
        "note": ("computed by the V62 production experiment (296^3, GPU) from the "
                 "simulation.hydrogen engine; no window screenshot, no overlay "
                 "text; the atom is the real hydrogen (an input of atomic "
                 "physics), the photon is the CL10 free massless field"),
    }
    os.makedirs(os.path.dirname(OUT_META), exist_ok=True)
    with open(OUT_META, "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print("[v62] wrote %s and %s" % (OUT_PNG, OUT_META), flush=True)
    return meta


# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    cfg = dict(runner.load_config())
    B.set_field_dtype(cfg.get("field_dtype", "float32"))

    r = run_production(cfg)
    ch = checks(r)
    meta = make_figure(r, cfg)

    ok = bool(
        ch["baseline_no_self_excitation"]["pass"]
        and ch["resonant_absorption"]["every_resonant_phase_absorbs"]
        and ch["resonant_absorption"]["lyman_alpha_excites_visibly"]
        and ch["off_resonant_control"]["mostly_not_absorbed"]
        and ch["shape_deformation"]["the_excitation_deforms_the_cloud"]
        and ch["einstein_a_radiation"]["pass"]
        and ch["stability"]["pass"]
    )

    verdict = {
        "id": "V62_hydrogen_production",
        "tier": 2,
        "property": "P2",
        "residual": True,
        "pass": ok,
        "status": ("the production 296^3 hydrogen experiment on the GPU: the resonant "
                   "photons (Lyman-alpha, 1s->3p, 2p->3d) are absorbed, the off-resonant "
                   "control is mostly not, the cloud deforms (1s sphere -> 2p dumbbell), "
                   "the excited atom radiates and decays, and the figure (Fig. 26) is "
                   "regenerated as a product of the experiment; the ATOM is an input "
                   "(real hydrogen), the PHOTON is CL10"),
        "data": {
            "backend": B.backend_name(),
            "gpu": B.gpu_info(),
            "grid": list(cfg["grid"]),
            "checks": ch,
            "figure": {
                "png": OUT_PNG,
                "meta": meta,
            },
            "provenance": {
                "atom_levels": ("the REAL hydrogen eigenstates (Rydberg) — an input of "
                                "atomic physics, NOT derived from the Cl(10) root"),
                "photon": ("the CL10 free massless 3-D field (the 6-neighbour leapfrog "
                           "engine, the same as the Young double-slit)"),
                "coupling_knobs": ("A, h, gamma: physical parameters of the atom, "
                                   "not fitted to any measured number"),
                "what_the_experiment_demonstrates": ("the photon sector of CL10 "
                                                     "propagates, resonates and is "
                                                     "emitted like the EM field; it "
                                                     "is NOT a 0-parameter derivation "
                                                     "of hydrogen"),
            },
            "ok": ok,
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "v62_hydrogen_production.json")
    with open(path, "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 72)
    print("V62 - the production 3-D hydrogen experiment on the GPU (Fig. 26)")
    print("=" * 72)
    print("  backend: %s   gpu: %s   grid: %s"
          % (verdict["data"]["backend"], verdict["data"]["gpu"],
             verdict["data"]["grid"]))
    ra = ch["resonant_absorption"]
    for p in ra["phases"]:
        print("  [%s]  t_fire=%d  peak excited=%.3f (t=%d)  "
              "upper-pop peak=%.3f  final=%.3f  F_out max=%.3e  decayed=%s"
              % (p["tag"], p["t_fire"], p["peak_excited"], p["t_peak"],
                 p["peak_upper_pop"], p["final_excited"], p["max_F_out"],
                 p["decayed_after_peak"]))
    print("  baseline (atom alone, no photon): max excited fraction "
          "%.2e  (pass: %s)"
          % (ch["baseline_no_self_excitation"]["max_excited_fraction"],
             ch["baseline_no_self_excitation"]["pass"]))
    orc = ch["off_resonant_control"]
    print("  off-resonant control: peak excited %.4f  (suppression "
          "%.1fx vs the weakest resonant; mostly-not-absorbed: %s)"
          % (orc["peak_excited"], orc["suppression_vs_min_resonant"],
             orc["mostly_not_absorbed"]))
    sd = ch["shape_deformation"]
    print("  shape deformation: anisotropy 1s=%.2e (pure)  2p=%.4f (pure)  "
          "Lyman-peak=%.4f  (dumbbell: %s)"
          % (sd["anisotropy_1s_pure"], sd["anisotropy_2p_pure"],
             sd["lyman_alpha_peak_anisotropy"],
             sd["the_excitation_deforms_the_cloud"]))
    st = ch["stability"]
    print("  stability: norm [%g, %g]  field max %.3e  (pass: %s)"
          % (st["norm_min"], st["norm_max"], st["field_max"], st["pass"]))
    print("  PROVENANCE: %s" % verdict["data"]["provenance"]["atom_levels"])
    print("              %s" % verdict["data"]["provenance"]["photon"])
    print("  VERDICT:  %s" % ("PASS" if ok else "FAIL"))
    print("  JSON -> " + path)
    print("  Wall time = %.2f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
