#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V317 - the selection conjunction (tier 3, F8/O9/O4, theory-hunt).

**The question (pre-registered).** Does ONE booked selector functional
uniquely select at >= 2 of the three unselected sites - S1: a member of
the 53-parameter lightcone family at the clock 1/5 (V313's artifact
members); S2: the six desert couplings kappa_s (the V105/V122/V145/V146
candidate families - cross-site consistency, not re-supply); S3: the
3-of-12 channel triple (the 12 chain-kernel channel indices in V272's
booked chain_H_spectrum order; the triples are the object)?

**The grid (12 cells, uniqueness definitions fixed pre-run).**
F1 locality T (V285 machinery rebuilt as V313 books; T_committed =
sum(profile(H)^2, r = 2..7)): S1 EVALUATED at the artifact members
(x_best and the committed member x = 0; V313's top5 books ONE non-penalty
member - booked, not invented later): unique iff T(member) < T_committed
strictly, interior by V313's own rule, one distinct minimum (V313's
n_starts_at_best leaf). S2, S3: N/A pre-run.
F2 bond-entropy weight G = 1/(4 ln 45) (V313's weight block): S1
weight-only, never a selector (V313's wording); S2 N/A; S3 EVALUATED:
G-weighted triple score, unique iff the gap to the runner-up > 1e-6
(booked: a uniform positive weight cannot reorder the triples - the
F2/F3 S3 cells share their gap).
F3 host energetics: S1 N/A; S2 EVALUATED via V247's v237_rule +
meet_solve at the V237-booked alpha_X against the four kappa_s families -
a family CLOSES iff its own booked desert-closure flag is true AND the
solved (zA, zB) e-folds are covered by the family's booked threshold
window where one exists; unique iff exactly one family closes. S3
EVALUATED: V302's ones_shift (sum of the chosen channels'
chain_H_spectrum entries), unique iff gap > 1e-6.
F4 relevance exponent (V213's p_A pattern): S1, S2 N/A; S3 EVALUATED:
mean over the triple's nodes of V213's power_ladder slope of
A_transverse over continuum.radii (tangent_field/lowest/velocity_spread
lifted verbatim, applied to the supercell step at clock 1/5 at the
triple's M nodes), unique iff gap > 1e-6.

**T0 (pre-registered).** V313's T_committed = 0.0011956422840809166
reproduced < 1e-9 from the same machinery (plus V313's Jacobian 107/53
gate); V247's booked (1,2)-shift row reproduced < 1e-9 from v237_rule at
alpha_std(pi/4). Failure: MACHINERY-DRIFT, nothing interpreted.

**Acceptance (pre-registered).** SELECTION-UNIFIED iff one functional
unique-selects at >= 3 sites; SELECTION-PARTIAL iff one functional at
exactly 2; SELECTION-DISPERSED iff every functional at <= 1. The full
12-cell table booked in all cases with the selecting member/triple/family
values when a cell is unique.
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402
sys.path.insert(0, exp("V285_most_local_lightlike_generator"))
sys.path.insert(0, exp("V247_two_dial_desert_contacts"))

import V285_most_local_lightlike_generator as V285        # noqa: E402
from V247_two_dial_desert_contacts import (v237_rule,     # noqa: E402
                                           meet_solve)
from cl10 import altcell as A                             # noqa: E402
from cl10 import config as _config                        # noqa: E402
from cl10 import coupling as _cpl                         # noqa: E402
from cl10 import supercell as SC                          # noqa: E402
from cl10 import unification as UN                        # noqa: E402

CFG = _config.CFG
G5 = float(CFG.angle.branches.implemented)
GAP_UNIQUE = 1e-6                        # pre-registered cell rule
TOL = float(CFG.tolerances.standard)
RADII = [float(r) for r in CFG.continuum.radii]
GRID = int(CFG.continuum.direction_grid)
G_ENTROPY = 1.0 / (4.0 * math.log(45.0))  # V313's weight block


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# ---- V213 machinery lifted verbatim (lowest / tangent_field /
# velocity_spread / power_ladder), applied to the supercell step at the
# clock - the S3 object's own step (the lift stated in the pre-reg) ----
_B, _U32, U16c = SC.sectors(G5)


def bloch(q):
    return SC.bloch(np.array(q, float), U16c, g10=G5)


def lowest(q):
    return float(np.min(np.abs(np.angle(np.linalg.eigvals(bloch(q))))))


def power_ladder(xs, ys):
    """log-log slope and prefactor over the usable (positive) points."""
    good = [(x, y) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(good) < 2:
        return None, None, len(good)
    x = np.array([math.log(a) for a, _ in good])
    y = np.array([math.log(b) for _, b in good])
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(math.exp(intercept)), len(good)


def velocity_spread(vs):
    """(v_max - v_min)/(v_max + v_min): the V115 anisotropy measure."""
    return (max(vs) - min(vs)) / (max(vs) + min(vs)) if len(vs) > 1 \
        else 0.0


def tangent_field(q0, r, dirs, tol):
    """The tangent dispersion at distance r from the node (V213
    verbatim)."""
    vs = [lowest(tuple(np.array(q0, float) + r * d)) / r for d in dirs]
    trans = [v for v in vs if v > tol]
    return {"r": r, "v_max": max(vs),
            "A_transverse": velocity_spread(trans),
            "E_mean": float(np.mean(vs)) * r}


def main():
    t0 = time.time()
    checks, data = [], {}

    def log(msg):
        print("[V317 t+{:5.0f}s] {}".format(time.time() - t0, msg),
              flush=True)
    log("start")

    # ================= F1 machinery: V313's V285 rebuild ================
    Rv, Dm = V285.module_reps()
    Bop, Kop = V285.velocity_ops()
    gen = json.load(open(os.path.join(
        exp("V208_metric_bilinear_cone"), "cone_generator.json")))
    H = [np.array(gen["A"][k]) + 1j * np.array(gen["S"][k])
         for k in range(3)]
    PROP = [i for i, R in enumerate(Rv) if np.linalg.det(R) > 0]

    def cov_project(Ts):
        outp = [np.zeros((36, 36), dtype=complex) for _ in range(3)]
        for i in PROP:
            D, R = Dm[i], Rv[i]
            DT = [D @ Tb @ D.conj().T for Tb in Ts]
            for a in range(3):
                for b in range(3):
                    outp[a] += R[a, b] * DT[b]
        return [o / 24. for o in outp]

    srng = np.random.default_rng(20261021)      # V285's basis seed

    def rherm():
        Aa = (srng.standard_normal((36, 36))
              + 1j * srng.standard_normal((36, 36)))
        return (Aa + Aa.conj().T) / 2.0
    rows, rank = [], 0
    for _ in range(4 * 160 + 60):
        v = V285.vec3(cov_project([rherm(), rherm(), rherm()]))
        r = int(np.linalg.matrix_rank(np.array(rows + [v]), tol=1e-9))
        if r > rank:
            rows.append(v)
            rank = r
        if rank == 160:
            break
    _, sv, Vt = np.linalg.svd(np.array(rows), full_matrices=False)
    Bmat = Vt[sv > 1e-9]

    def linearised(dH):
        d1, d2, d3 = dH
        return [V285.sym(H[0], d2) + V285.sym(d1, H[1]),
                V285.sym(H[0], d3) + V285.sym(d1, H[2]),
                V285.sym(H[1], d3) + V285.sym(d2, H[2]),
                V285.sym(H[0], d1) - V285.sym(H[1], d2),
                V285.sym(H[0], d1) - V285.sym(H[2], d3)]
    cols = []
    for v in Bmat:
        ms = linearised(V285.unvec(v))
        cols.append(np.concatenate([m.real.reshape(-1) for m in ms]
                                   + [m.imag.reshape(-1) for m in ms]))
    Jm = np.array(cols).T
    jac_rank = int(np.linalg.matrix_rank(Jm, tol=1e-9))
    family_dim = 160 - jac_rank

    profH = V285.profile(H)
    T_committed = sum(profH[r] ** 2 for r in range(2, 8))
    a_Tcommit = float(art("V285_most_local_lightlike_generator",
                          "data.t0.T_committed"))
    t0_f1 = abs(T_committed - a_Tcommit)
    data["t0_T_committed"] = {"rederived": T_committed,
                              "artifact": a_Tcommit,
                              "abs_diff": t0_f1,
                              "rank": jac_rank, "nullity": family_dim}
    checks.append({"name": "T0a the V285 machinery rebuilt: Jacobian rank "
                          "{} / nullity {} must be 107/53; T_committed "
                          "{:.6e} reproduces the booked {:.6e} (|d| "
                          "{:.1e} < 1e-9)".format(
                              jac_rank, family_dim, T_committed,
                              a_Tcommit, t0_f1),
                   "pass": bool(jac_rank == 107 and family_dim == 53
                                and t0_f1 < 1e-9)})
    wH = np.linalg.eigvalsh((H[0] + H[0].conj().T) / 2)
    b1 = int(np.sum(np.abs(wH - 1.0) < 1e-6))
    bm1 = int(np.sum(np.abs(wH + 1.0) < 1e-6))
    _, svJ, Vh = np.linalg.svd(Jm, full_matrices=False)
    NULLc = Vh[jac_rank:]
    V285._ST.update({"Bmat": Bmat, "NULL": NULLc,
                     "BOX": math.sqrt(V285.i3(H, H) / 3.0),
                     "BR": (b1, bm1)})
    vH = V285.vec3(H)
    cH = np.array([float(np.dot(vH, Bmat[j])) for j in range(160)])
    V285._ST["cH"] = cH
    J0 = np.zeros((len(V285._resid(cH)), 160))
    eps = 1e-6
    for j in range(160):
        cp = cH.copy()
        cp[j] += eps
        J0[:, j] = (V285._resid(cp) - V285._resid(cH)) / eps
    Q0, R0 = np.linalg.qr(J0)
    V285._ST.update({"J0": J0, "QR": (Q0, R0)})
    log("F1 machinery ready")

    # ---- F3 desert T0: the V247 booked (1,2)-shift row ----------------
    a_dz12 = float(art("V237_desert_threshold_matching",
                       "data.matching.dz12_efolds"))
    a_ax = float(art("V237_desert_threshold_matching",
                     "data.matching.alpha_X"))
    fires0, t0row = v237_rule(_cpl.alpha_std(math.pi / 4))
    t0_f3 = abs(t0row["dz12_efolds"] - a_dz12)
    data["t0_desert_row"] = {"dz12_rederived": t0row["dz12_efolds"],
                             "dz12_artifact": a_dz12,
                             "abs_diff": t0_f3, "fires": fires0,
                             "alpha_X_booked": a_ax}
    checks.append({"name": "T0b V247's booked (1,2)-shift row: v237_rule "
                          "at alpha_std(pi/4) reproduces dz12 {:.15g} "
                          "(|d| {:.1e} < 1e-9, host float64 rule; fires "
                          "= {})".format(a_dz12, t0_f3, fires0),
                   "pass": bool(t0_f3 < 1e-9)})
    if not (jac_rank == 107 and family_dim == 53 and t0_f1 < 1e-9
            and t0_f3 < 1e-9):
        data["status"] = "MACHINERY-DRIFT"
        print("V317 MACHINERY-DRIFT: the V285/V237 anchors do not "
              "reproduce; nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ==================== F1 x S1 ======================================
    x_best = [float(v) for v in art("V313_lightcone_dictionary_conjunction",
                                    "data.search.x_best")]
    n_at_best = int(art("V313_lightcone_dictionary_conjunction",
                        "data.search.n_starts_at_best"))
    top5 = [float(v) for v in art("V313_lightcone_dictionary_conjunction",
                                  "data.search.top5")]
    BOX = V285._ST["BOX"]

    def member_diag(x):
        xv = np.array(x, float)
        Tm = V285._obj(xv)
        g53 = np.zeros(53)
        hstep = 1e-5
        for a in range(53):
            xp = xv.copy()
            xm = xv.copy()
            xp[a] += hstep
            xm[a] -= hstep
            g53[a] = (V285._obj(xp) - V285._obj(xm)) / (2 * hstep)
        n_box = int(np.sum(np.abs(xv) > BOX - 1e-6))
        interior = (float(np.linalg.norm(g53))
                    < 1e-6 * max(1.0, T_committed) and n_box == 0)
        below = bool(Tm < T_committed)
        return {"T": float(Tm), "grad_norm": float(np.linalg.norm(g53)),
                "box_face_coords": n_box, "interior": interior,
                "below_T_committed": below}
    members = {"T_star member (V313 data.search.x_best)": x_best,
               "committed member (x = 0)": [0.0] * 53}
    s1_rows = {k: member_diag(x) for k, x in members.items()}
    single_distinct = (n_at_best == 1)
    s1_unique = any(r["below_T_committed"] and r["interior"]
                    for r in s1_rows.values()) and single_distinct
    data["F1S1"] = {"members": s1_rows,
                    "top5_leaf": top5,
                    "n_starts_at_best_leaf": n_at_best,
                    "single_distinct_minimum": single_distinct,
                    "unique": bool(s1_unique)}
    log("F1 x S1 evaluated")

    # ==================== F3 x S2 (kappa_s families) ===================
    mZ = float(UN.M_Z)
    xfix = math.log(float(CFG.desert.m_x_closed_gev) / mZ)
    sol = meet_solve(a_ax)                    # (x, zA, zB) e-folds
    zfolds = {"zA": float(sol[1]), "zB": float(sol[2])}
    fam_rows = {}

    def window_ok(zmin, zmax):
        return (zmin <= zfolds["zA"] <= zmax and
                zmin <= zfolds["zB"] <= zmax)

    # V105 cartan kappas: per-model desert_ok flags (point z0, no window)
    j105 = artj("V105_cartan_kappas")
    f105 = any(bool(r.get("desert_ok", False))
               for m in j105["data"]["models"].values()
               for r in m.values() if isinstance(r, dict))
    fam_rows["V105_cartan_kappas"] = {
        "flag_key": "data.models.*.<n>.desert_ok", "own_flag": f105,
        "window": None,
        "window_covered": None,
        "closes": bool(f105)}
    # V122 class-invariant kappas: free_kappa rows closed + [z_min, z_max]
    j122 = artj("V122_class_invariant_kappas")
    f122 = False
    w122 = None
    for emb, ev_ in j122["data"]["embeddings"].items():
        for row in ev_["free_kappa"]["rows"]:
            if row.get("closed") and row["z_min"] > 0 \
                    and row["z_max"] < xfix:
                if window_ok(row["z_min"], row["z_max"]):
                    f122 = True
                    w122 = [row["z_min"], row["z_max"], emb]
    fam_rows["V122_class_invariant_kappas"] = {
        "flag_key": "data.embeddings.<e>.free_kappa.rows "
                    "(closed, [z_min, z_max] in the desert e-fold window "
                    "(0, {:.3f}))".format(xfix),
        "own_flag": any(r.get("closed")
                        for e_ in j122["data"]["embeddings"].values()
                        for r in e_["free_kappa"]["rows"]),
        "window": w122,
        "window_covered": bool(f122), "closes": bool(f122)}
    # V145 superexchange kappas: meeting rows thresholds_in_desert
    j145 = artj("V145_superexchange_kappas")
    f145_flag = any(bool(r.get("thresholds_in_desert"))
                    for r in j145["data"]["meeting"])
    w145 = None
    f145 = False
    for r in j145["data"]["meeting"]:
        if bool(r.get("thresholds_in_desert")) and \
                window_ok(r["min_threshold"], r["max_threshold"]):
            f145 = True
            w145 = [r["min_threshold"], r["max_threshold"], r["n"]]
    fam_rows["V145_superexchange_kappas"] = {
        "flag_key": "data.meeting[*].thresholds_in_desert + "
                    "[min_threshold, max_threshold]",
        "own_flag": f145_flag, "window": w145,
        "window_covered": bool(f145), "closes": bool(f145)}
    # V146 two-number kappas: joint per_n in_desert flags
    j146 = artj("V146_two_number_kappas")
    f146 = any(bool(r.get("in_desert"))
               for r in j146["data"]["joint"].get("per_n", [])) \
        or j146["data"]["joint"].get("n_desert", 0) > 0
    fam_rows["V146_two_number_kappas"] = {
        "flag_key": "data.joint.per_n[*].in_desert / data.joint.n_desert",
        "own_flag": bool(j146["data"]["joint"].get("n_desert", 0) > 0),
        "window": None, "window_covered": None, "closes": bool(f146)}
    n_closed = sum(1 for r in fam_rows.values() if r["closes"])
    s2_unique = (n_closed == 1)
    data["F3S2"] = {"solved_z_folds": zfolds,
                    "desert_xfix_efolds": xfix,
                    "families": fam_rows, "n_closed": n_closed,
                    "unique": bool(s2_unique)}
    log("F3 x S2 evaluated")

    # ============ F2 x S3 and F3 x S3: chain-H energetics ==============
    E = [float(v) for v in art("V272_dynamical_chain",
                               "data.chain_H_spectrum")]
    assert len(E) == 12, "the chain-H spectrum is the booked 12"
    triples = list(itertools.combinations(range(12), 3))
    scores = sorted((sum(E[i] for i in T), T) for T in triples)
    gap = abs(scores[0][0] - scores[1][0])
    triple_best = list(scores[0][1])
    s3_energy_unique = (gap > GAP_UNIQUE)
    gapG = gap * G_ENTROPY          # uniform weight: same ordering
    data["F3S3"] = {"channel_energies": E, "best_triple": triple_best,
                    "best_score": scores[0][0], "runner_up_score":
                    scores[1][0], "gap": gap,
                    "unique": bool(s3_energy_unique)}
    data["F2S3"] = {"G": G_ENTROPY, "weight_only_note":
                    "a uniform positive weight cannot reorder the "
                    "triples (booked pre-run): G*score has the same "
                    "ordering and gap x G",
                    "best_triple": triple_best,
                    "best_score_weighted": G_ENTROPY * scores[0][0],
                    "gap_weighted": gapG,
                    "unique": bool(gapG > GAP_UNIQUE)}
    log("F2/F3 x S3 evaluated")

    # ==================== F4 x S3 (p_A of the triple) ==================
    dirs = A.sphere_dirs(GRID)
    p_nodes = {}
    for i, q in enumerate(SC.M_POINTS):
        ladder = [tangent_field(q, r, dirs, TOL) for r in RADII]
        pA, cA, nA = power_ladder([row["r"] for row in ladder],
                                  [row["A_transverse"] for row in ladder])
        p_nodes[str(i)] = {"p_A": pA, "usable_radii": nA}
    pv = [(p_nodes[str(i)]["p_A"] if p_nodes[str(i)]["p_A"] is not None
           else float("nan")) for i in range(12)]
    tscores = []
    for T in triples:
        vals = [pv[i] for i in T]
        if any(math.isnan(v) for v in vals):
            continue
        tscores.append((sum(vals) / 3.0, T))
    tscores.sort()
    p_gap = abs(tscores[0][0] - tscores[1][0])
    data["F4S3"] = {"node_p_A": p_nodes,
                    "best_triple": list(tscores[0][1]),
                    "best_p_A": tscores[0][0],
                    "runner_up_p_A": tscores[1][0],
                    "gap": p_gap,
                    "unique": bool(p_gap > GAP_UNIQUE),
                    "radii": RADII, "direction_grid": GRID,
                    "n_dirs": len(dirs)}
    log("F4 x S3 evaluated")

    # ==================== the 12-cell table and branch =================
    cells = {
        "F1": {"S1": ("evaluated", s1_unique),
               "S2": ("N/A (pre-run: T not applicable to a coupling "
                      "family)", False),
               "S3": ("N/A (pre-run: T has no channel argument)", False)},
        "F2": {"S1": ("weight-only (V313's booked wording; never a "
                      "selector)", False),
               "S2": ("N/A (pre-run)", False),
               "S3": ("evaluated", bool(gapG > GAP_UNIQUE))},
        "F3": {"S1": ("N/A (pre-run: host energetics has no lightcone "
                      "family argument)", False),
               "S2": ("evaluated", s2_unique),
               "S3": ("evaluated", s3_energy_unique)},
        "F4": {"S1": ("N/A (pre-run: p_A measures anisotropy at a fixed "
                      "member, not members)", False),
               "S2": ("N/A (pre-run)", False),
               "S3": ("evaluated", bool(p_gap > GAP_UNIQUE))}}
    per_fun = {f: sum(1 for s in ("S1", "S2", "S3")
                      if cells[f][s][1]) for f in cells}
    data["table"] = {f: {s: [cells[f][s][0], cells[f][s][1]]
                          for s in ("S1", "S2", "S3")}
                     for f in cells}
    data["sites_selected_per_functional"] = per_fun
    mx = max(per_fun.values())
    branch = ("SELECTION-UNIFIED" if mx >= 3 else
              "SELECTION-PARTIAL" if mx == 2 else "SELECTION-DISPERSED")
    data["status"] = branch
    checks.append({"name": "branch {}: 12-cell table (unique flags: {}); "
                          "sites per functional {}; F1 S1 {} / F3 S2 {} "
                          "families closed / F2+F3 S3 best triple {} gap "
                          "{:.2e} / F4 S3 best triple {} p_A gap "
                          "{:.2e}".format(
                              branch,
                              {f: [cells[f][s][1]
                                   for s in ("S1", "S2", "S3")]
                               for f in cells},
                              per_fun, "unique" if s1_unique else "not",
                              "unique" if s2_unique else "not",
                              triple_best, gap,
                              list(tscores[0][1]), p_gap),
                   "pass": True})
    print("V317 {}: sites per functional {}; F1S1 unique {}, F3S2 "
          "closed families {}, F2/F3 S3 triple {} (gap {:.2e}), F4 S3 "
          "triple {} (gap {:.2e}).".format(
              branch, per_fun, s1_unique,
              [k for k, v in fam_rows.items() if v["closes"]],
              triple_best, gap, list(tscores[0][1]), p_gap), flush=True)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, all(c["pass"] for c in checks))
    return 0 if all(c["pass"] for c in checks) else 1


def artj(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def _emit(data, checks, ok):
    out = {"id": "V317_selection_conjunction", "tier": 3,
           "property": "does one booked selector functional uniquely "
                       "select at >= 2 of the lightcone family, the "
                       "kappa_s families, and the 3-of-12 channel "
                       "triple (tier 3, F8/O9/O4)",
           "residual": True, "programme": "theory-hunt V317",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v317_selection_conjunction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V317 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
