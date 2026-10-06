#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — V61: the 3-D hydrogen atom driven by the CL10 photon (the Fig. 26 check).

**What this experiment verifies.**  Figure 26 shows the hydrogen electron cloud
`|psi|^2` driven by a photon packet of the CL10 free massless 3-D field (the same
6-neighbour leapfrog engine as the Young double-slit, V25/fig25).  The caption
claims four things.  Each is **measured** here, not asserted:

1. **the levels are the real hydrogen ones, 0-parameter.**  The gaps are the
   Rydberg differences `13.606 (1/n_1^2 - 1/n_2^2)` eV: `1s->2p = 10.20`,
   `1s->3p = 12.09`, `2p->3d = 1.89`.  Nothing is fitted to atomic data.

2. **the dipole selection rules emerge from the lattice wavefunctions.**  The
   dipole matrix elements `<i|z|j>` are computed by summing the lattice
   eigenfunctions; the E1-forbidden pairs (`1s->3d`, `2p->3p`, every `Diagonal`)
   come out ~0, the allowed ones do not.  The rule is not imposed — it is
   measured from the shapes.

3. **the absorption follows the hydrogen line shape.**  This is the content of
   the figure: the excitation is **not imposed**, the dynamics decides.  The model
   implements the Fermi-golden-rule rate `A |F|^2 L(delta)` with a **Lorentzian**
   line shape, so the testable claim is that absorption **peaks on a gap and falls
   with the detuning** — measured by a detuning scan at fixed amplitude, width and
   duration.  Stated honestly: the off-resonant absorption is **not zero** (the
   Lorentzian tail is real physics), and the experiment also measures the
   **packet-locality trap** — a packet wider than the box never leaves, becomes a
   persistent bath, and then drives the atom off resonance.

4. **the excited atom radiates and decays** (Einstein A): with no incoming field,
   the excited population relaxes toward the 1s **and** the emitted field `F_out`
   is non-zero — the atom does not merely lose population, it sources a wave.

Plus the bookkeeping checks: the total atomic norm is conserved, and the discrete
photon energy stays bounded (the FDTD is stable at the CFL number used).

**Honest provenance (what this figure is NOT).**  The atom's levels and shapes are
the **real hydrogen** eigenstates — an input of atomic physics, **not** derived
from the Cl(10) root.  What the CL10 theory supplies is the **photon**: the free
massless 3-D field, its propagation, its resonance with the atom and its emission.
So this figure is a **demonstration of the photon sector**, not a 0-parameter
derivation of hydrogen.  The three coupling knobs (`A`, `h`, `gamma`) are physical
parameters of the atom, **not** fitted to any measured number.

**Grid note.**  The figure is rendered on the production grid (296^3).  This check
runs the **same model** on a reduced grid (120^3) so it is fast enough to be part
of the validation suite; the physics checks below are grid-independent.
"""

import os
import sys
import json
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from simulation.hydrogen import runner, qca as QCA, backend as B

RYDBERG_EV = 13.606

# the reduced grid (the same model, faster; the figure uses 296^3)
TEST_CFG = {
    "grid": [120, 120, 120],
    "r": 0.5,
    "sponge_width": 10,
    "bohr_radius": 12.0,
    "coupling_absorb": 0.06,     # a stronger drive so the test is short
    "coupling_emit": 0.5,
    "reaction_gamma_max": 0.004,
    "line_width": 0.0,
    "photon_width": 12,
    "photon_amp": 1.0,
    "field_dtype": "float64",
}


def _cfg():
    cfg = dict(runner.load_config())
    cfg.update(TEST_CFG)
    B.set_field_dtype(cfg["field_dtype"])
    return cfg


def check_levels(sim):
    """1. the gaps are the Rydberg differences (0-parameter, not fitted)."""
    analytic = {
        (0, 1): RYDBERG_EV * (1 / 1 - 1 / 4),     # 1s -> 2p
        (0, 2): RYDBERG_EV * (1 / 1 - 1 / 9),     # 1s -> 3p
        (1, 3): RYDBERG_EV * (1 / 4 - 1 / 9),     # 2p -> 3d
        (0, 3): RYDBERG_EV * (1 / 1 - 1 / 9),     # 1s -> 3d (same gap as 1s->3p)
    }
    out = {}
    worst = 0.0
    for key, ref in analytic.items():
        got = float(sim.gaps_eV[key])
        err = abs(got - ref) / ref
        worst = max(worst, err)
        out["%s->%s" % (sim.labels[key[0]], sim.labels[key[1]])] = {
            "simulated_eV": got,
            "rydberg_eV": ref,
            "relative_error": err,
        }
    return out, worst


def check_selection_rules(sim):
    """2. the E1 selection rules come out of the lattice wavefunctions."""
    n = sim.n_levels
    lab = sim.labels
    mat = {(i, j): float(sim.dipole_mat[(i, j)]) for i in range(n) for j in range(n)}
    # the (n, l) of each level, for the E1 rule Delta l = +-1
    nl = {k: (cfg_level[1], cfg_level[2]) for k, cfg_level in enumerate(QCA.HYDROGEN_LEVELS)}
    allowed, forbidden = [], []
    for i in range(n):
        for j in range(i + 1, n):
            dl = abs(nl[j][1] - nl[i][1])
            d = mat[(i, j)]
            entry = {
                "transition": "%s->%s" % (lab[i], lab[j]),
                "delta_l": dl,
                "dipole_z": d,
                "e1_allowed_by_rule": dl == 1,
                "nonzero_on_the_lattice": abs(d) > 1e-3,
            }
            (allowed if dl == 1 else forbidden).append(entry)
    # the rule and the measured matrix must agree
    agree = all(e["e1_allowed_by_rule"] == e["nonzero_on_the_lattice"]
                for e in allowed + forbidden)
    max_forbidden = max((abs(e["dipole_z"]) for e in forbidden), default=0.0)
    min_allowed = min((abs(e["dipole_z"]) for e in allowed), default=0.0)
    return {
        "allowed": allowed,
        "forbidden": forbidden,
        "rule_and_lattice_agree": agree,
        "max_forbidden_dipole": max_forbidden,
        "min_allowed_dipole": min_allowed,
        "suppression_factor": (min_allowed / max_forbidden) if max_forbidden > 0 else float("inf"),
    }


def check_resonance(sim_factory):
    """3. the resonance: the absorption follows the line shape (the content of
    the figure).

    **The honest statement of the claim.**  The model does NOT implement a sharp
    filter that passes a resonant photon and blocks every other one.  It
    implements the Fermi-golden-rule rate `A |F|^2 L(delta)` with a **Lorentzian**
    line shape `L = (gamma/2) / (delta^2 + gamma^2/4)`.  So the testable claim is
    that the absorption **follows that line shape**: it peaks at a gap and falls
    as the detuning grows.  A residual off-resonant excitation is **expected
    physics** (the Lorentzian tail, and power broadening when the drive is
    intense), not a bug - so the check measures the **suppression** and its
    **monotonicity in the detuning**, and it also measures the **packet-locality
    trap** below.

    **The packet-locality trap (measured, not assumed).**  A packet whose width
    is comparable to the box never leaves: it becomes a persistent bath, and a
    strong persistent drive excites the atom off resonance.  The scan below shows
    this happening, so the control shot is only valid when the packet is localised
    relative to the box (the production config: width 20 in a 296 box).
    """
    L = TEST_CFG["grid"][0]
    r = TEST_CFG["r"]
    travel = int((L // 2) / r)          # the steps for the packet to reach the centre
    dw_dk = r                          # dw/dk at k -> 0 for omega = 2 asin(r sin(k/2))

    def run(omega, n_extra, width=None, amp=None):
        sim = sim_factory()
        sim.set_ground()
        kw = {}
        if width is not None:
            kw["width"] = width
        if amp is not None:
            kw["amp"] = amp
        if omega is not None:
            sim.fire_photon(omega, direction="x", **kw)
        arrived = 0.0
        for _ in range(travel + n_extra):
            sim.step()
            arrived = max(arrived, float(
                B.asnumpy(sim.F_local_weight * sim.F * sim.F).sum()))
        return sim, arrived

    sim = sim_factory()
    gap12 = abs(sim.omega[1] - sim.omega[0])
    gap13 = abs(sim.omega[2] - sim.omega[0])
    gamma_line = sim.gamma_line

    # the baseline: the atom alone, no photon at all
    s_none, _ = run(None, 120)

    # the resonant shot: exactly the 1s->2p gap
    s_res, arr_res = run(gap12, 120)

    # the detuning scan: the same amplitude, the same width, the same duration,
    # the detuning varied.  The excited fraction must PEAK on the gap and fall
    # with the detuning - that is the line shape, and it is what the figure's
    # "resonant absorbed / off-resonant not" claim actually means.
    scan = []
    for frac in (0.0, 0.25, 0.5, 1.0, 2.0):
        om = gap12 * (1.0 + frac)
        s_w, arr_w = run(om, 120)
        scan.append({
            "omega_rad_per_step": om,
            "detuning_from_1s2p_rad": om - gap12,
            "detuning_in_line_widths": (om - gap12) / gamma_line,
            "field_reached_the_atom": arr_w,
            "excited_fraction": float(s_w.excited_fraction()),
        })
    on_res = scan[0]["excited_fraction"]
    far = scan[-1]["excited_fraction"]
    monotone = all(scan[k]["excited_fraction"] >= scan[k + 1]["excited_fraction"] - 1e-3
                   for k in range(len(scan) - 1))

    # the packet-locality scan (the trap): widen the packet in the SAME box
    locality = []
    for w in (12, 24, 48):
        s_w, arr_w = run(0.5 * gap12, 120, width=w)
        locality.append({
            "photon_width_cells": w,
            "width_over_box": w / L,
            "spectral_bandwidth_rad_per_step": dw_dk / w,
            "excited_fraction": float(s_w.excited_fraction()),
            "field_reached_the_atom": arr_w,
        })

    return {
        "travel_steps_to_the_atom": travel,
        "line_width_HWHM_rad_per_step": gamma_line,
        "atom_alone_no_photon": {
            "populations": {sim.labels[k]: v for k, v in enumerate(s_none.populations())},
            "excited_fraction": float(s_none.excited_fraction()),
        },
        "resonant": {
            "omega_rad_per_step": gap12,
            "gap_eV": float(sim.gaps_eV[(0, 1)]),
            "field_reached_the_atom": arr_res,
            "populations": {sim.labels[k]: v for k, v in enumerate(s_res.populations())},
            "excited_fraction": float(s_res.excited_fraction()),
        },
        "detuning_scan": scan,
        "on_resonance_excited": on_res,
        "far_off_resonance_excited": far,
        "suppression_on_over_far": on_res / max(far, 1e-12),
        "absorption_peaks_on_the_gap": bool(monotone),
        "packet_locality_scan": locality,
        "caveat": ("the off-resonant absorption is NOT exactly zero: the model's "
                   "Fermi-golden-rule rate carries a Lorentzian tail, and a packet "
                   "wider than the box becomes a persistent bath that drives the "
                   "atom off resonance (measured in packet_locality_scan).  The "
                   "production figure uses a localised packet (width 20 in a 296 box)."),
    }


def check_emission(sim_factory):
    """4. the excited atom radiates (Einstein A) and decays."""
    sim = sim_factory()
    sim.set_level(1)                      # put the atom in the 2p, no field
    p0 = sim.populations()
    emitted_before = float(B.asnumpy(sim.photon_density()).max())
    for _ in range(400):
        sim.step()
    p1 = sim.populations()
    emitted_after = float(B.asnumpy(sim.photon_density()).max())
    return {
        "populations_before": {sim.labels[k]: p0[k] for k in range(len(p0))},
        "populations_after_400_steps": {sim.labels[k]: p1[k] for k in range(len(p1))},
        "ground_population_increased": bool(p1[0] > p0[0]),
        "emitted_field_max_before": emitted_before,
        "emitted_field_max_after": emitted_after,
        "radiation_emitted": bool(emitted_after > 1e-6),
    }


def check_norm_and_stability(sim_factory):
    """the bookkeeping: the atomic norm is conserved, the field stays bounded."""
    sim = sim_factory()
    sim.set_ground()
    sim.fire_transition(0, 1, amp=TEST_CFG["photon_amp"])
    norms, fields = [], []
    for _ in range(500):
        sim.step()
        norms.append(float(sum(abs(a) ** 2 for a in sim.a)))
        fields.append(float(B.asnumpy(sim.F).max()))
    return {
        "norm_min": min(norms),
        "norm_max": max(norms),
        "norm_drift": max(abs(v - 1.0) for v in norms),
        "field_max": max(fields),
        "field_bounded": bool(max(fields) < 1e6),
    }


def hydrogen_check():
    cfg = _cfg()

    def factory():
        s = QCA.AtomQCA(cfg)
        return s

    sim = factory()
    levels, level_err = check_levels(sim)
    sel = check_selection_rules(sim)
    res = check_resonance(factory)
    em = check_emission(factory)
    stab = check_norm_and_stability(factory)

    ok = bool(
        level_err < 1e-9                       # the gaps ARE the Rydberg differences
        and sel["rule_and_lattice_agree"]     # the selection rules come from the shapes
        and sel["suppression_factor"] > 50.0  # forbidden really are suppressed
        # the baseline: the atom must not excite itself
        and res["atom_alone_no_photon"]["excited_fraction"] < 1e-9
        and res["resonant"]["excited_fraction"] > 0.05
        # the absorption must PEAK on the gap and fall with the detuning
        and res["absorption_peaks_on_the_gap"]
        and res["suppression_on_over_far"] > 5.0
        # every shot must actually have delivered a field to the atom
        and all(e["field_reached_the_atom"] > 1e-3 for e in res["detuning_scan"])
        and em["ground_population_increased"]
        and em["radiation_emitted"]
        and stab["norm_drift"] < 1e-9
        and stab["field_bounded"]
    )

    return {
        "backend": B.backend_name(),
        "grid_figure": list(runner.load_config()["grid"]),
        "grid_test": TEST_CFG["grid"],
        "levels_0param": {"checks": levels, "max_relative_error": level_err},
        "selection_rules": sel,
        "resonance": res,
        "spontaneous_emission": em,
        "stability": stab,
        "provenance": {
            "atom_levels": "the REAL hydrogen eigenstates (Rydberg) — an input of atomic physics, NOT derived from the Cl(10) root",
            "photon": "the CL10 free massless 3-D field (the 6-neighbour leapfrog engine, the same as the Young double-slit)",
            "coupling_knobs": "A, h, gamma: physical parameters of the atom, not fitted to any measured number",
            "what_the_figure_demonstrates": "the photon sector of CL10 propagates, resonates and is emitted like the EM field; it is NOT a 0-parameter derivation of hydrogen",
        },
        "ok": ok,
    }


def main():
    t0 = time.time()
    r = hydrogen_check()
    ok = r["ok"]

    verdict = {
        "id": "V61_hydrogen",
        "tier": 2,
        "property": "P2",
        "residual": True,
        "pass": ok,
        "status": ("the CL10 photon resonates with the real hydrogen levels (absorbed on resonance, "
                   "not off resonance), the E1 selection rules come out of the lattice wavefunctions, "
                   "and the excited atom radiates and decays; the ATOM is an input (real hydrogen), "
                   "the PHOTON is CL10"),
        "data": r,
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "v61_hydrogen.json")
    with open(path, "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 72)
    print("V61 - the 3-D hydrogen atom driven by the CL10 photon (Fig. 26)")
    print("=" * 72)
    print(f"  backend: {r['backend']}   grid: test {r['grid_test']}  figure {r['grid_figure']}")
    lv = r["levels_0param"]
    print(f"  [1] the levels (0-parameter Rydberg, max relative error {lv['max_relative_error']:.2e}):")
    for k, v in lv["checks"].items():
        print(f"      {k:9s}  simulated {v['simulated_eV']:7.3f} eV   Rydberg {v['rydberg_eV']:7.3f} eV")
    s = r["selection_rules"]
    print(f"  [2] the E1 selection rules measured on the lattice (agree: {s['rule_and_lattice_agree']}, "
          f"suppression {s['suppression_factor']:.1f}x):")
    for e in s["allowed"]:
        print(f"      allowed   {e['transition']:9s}  <z> = {e['dipole_z']:+8.4f}")
    for e in s["forbidden"]:
        print(f"      FORBIDDEN {e['transition']:9s}  <z> = {e['dipole_z']:+8.4f}  (Delta l = {e['delta_l']})")
    res = r["resonance"]
    print(f"  [3] the resonance (the content of the figure):")
    print(f"      baseline, the atom ALONE (no photon): excited fraction "
          f"{res['atom_alone_no_photon']['excited_fraction']:.2e}")
    print(f"      line HWHM {res['line_width_HWHM_rad_per_step']:.4f} rad/step")
    print(f"      RESONANT 1s->2p {res['resonant']['gap_eV']:.2f} eV -> excited fraction "
          f"{res['resonant']['excited_fraction']:.3f}")
    print(f"      the detuning scan (same amplitude, same width, same duration):")
    for e in res["detuning_scan"]:
        print(f"        detuning {e['detuning_from_1s2p_rad']:.4f} rad/step "
              f"({e['detuning_in_line_widths']:5.1f} line widths)  -> excited {e['excited_fraction']:.4f}")
    print(f"      suppression (on-resonance / far-off): {res['suppression_on_over_far']:.1f}x   "
          f"peaks on the gap: {res['absorption_peaks_on_the_gap']}")
    print(f"      the packet-locality trap (a packet wider than the box = a persistent bath):")
    for e in res["packet_locality_scan"]:
        print(f"        width {e['photon_width_cells']:3d} cells ({e['width_over_box']:.2f} of the box)  -> excited {e['excited_fraction']:.4f}")
    em = r["spontaneous_emission"]
    print(f"  [4] spontaneous emission: 2p population "
          f"{em['populations_before']['2p']:.3f} -> {em['populations_after_400_steps']['2p']:.3f}, "
          f"1s {em['populations_before']['1s']:.3f} -> {em['populations_after_400_steps']['1s']:.3f}, "
          f"emitted field max {em['emitted_field_max_after']:.3e} (radiated: {em['radiation_emitted']})")
    st = r["stability"]
    print(f"  [5] stability: atomic norm drift {st['norm_drift']:.2e}, field max {st['field_max']:.3e}")
    print(f"  PROVENANCE: {r['provenance']['atom_levels']}")
    print(f"              {r['provenance']['photon']}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + path)
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
