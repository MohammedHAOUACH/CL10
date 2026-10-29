#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V311 - the quaternionic root step (tier 3, O5/O8, pillar-surgery).

**The question (pre-registered).** the real form Cl(10,0) carries no
CP-odd invariant (V117).  The quaternionic real forms of the SAME
complex algebra Cl_10(C) - p in pillar_surgery.quaternionic_p - do.
Does the automaton's step survive a signature swap?

**The object (closed input set).** generators by exact construction:
from cl10.clifford.generators() the 10 kron matrices Gamma_j;
Gamma'_j = Gamma_j for the first p generators (j < p, 0-indexed),
Gamma'_j = i Gamma_j for the rest - the real form Cl(p, 10 - p) of the
V117 table (anticommutation preserved, squares flip).  Plane generators
and step
copy the 4-line pattern of cl10/sm_qca.py (T'_ab = i Gamma'_a Gamma'_b,
plane_rotation U = cos(theta) I - i sin(theta) T', step = ordered
product) into this experiment as signature_step(p, g).  NOTE (booked
in the pre-registration): the Gamma_j -> i Gamma_j (j > p) construction
is a choice of real-form embedding inside the FIXED complex 32-dim
rep; the unitarity verdict is with respect to this embedding, so a
no-go is class-bounded and honest.

**Measured (pre-registered).** per p:
1. T0: the FS indicator of the module, fs_indicator imported from V117
   (line 83) - must read quaternionic with J^2 = -1, else MACHINERY-
   DRIFT.
2. the CP-odd test of the V117 T3 pattern: the intertwiner Q of the
   swapped generators and the deviation of the antiunitary image
   J U' J^-1 = Q U'^* Q^dagger... against U'^-1, computed on the step
   at the closed clock angle.branches.closed.
3. the unitarity residual ||U'U'^dagger - I||_max (the V208 unitary_err
   pattern, computed in-experiment, printed to full precision).
4. the closure census: for every g = pi*a/b with 0 < a < b <=
   pillar_surgery.closure_den_max, the smallest k <= closure_k_max
   with ||U'(g)^k -/+ I||_max < 1e-12; found closure clocks booked per
   signature.

**Acceptance (pre-registered).** STEP-UNITARY per signature (some p
keeps exact unitarity at tolerance.strict - its closure clocks booked)
/ STEP-HYPERBOLIC-DECLARED if all quaternionic signatures fail
unitarity: the signature surgery dies at the step level, a
theorem-grade booking (V312's skip is pre-decided off that word).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V117_real_form"))

from V117_real_form import fs_indicator  # noqa: E402
from cl10 import clifford  # noqa: E402
from cl10 import config as _config  # noqa: E402

CFG = _config.CFG
TOL = float(CFG.tolerances.strict)          # exactness band (1e-12)
G_CLOSED = float(CFG.angle.branches.closed)
PS = [int(p) for p in CFG.pillar_surgery.quaternionic_p]
DEN_MAX = int(CFG.pillar_surgery.closure_den_max)
K_MAX = int(CFG.pillar_surgery.closure_k_max)
CLOSE_EPS = 1e-12                           # closure vanishing band


def swapped_generators(p):
    """Gamma'_j = Gamma_j for the first p generators (0-indexed j < p),
    i Gamma_j for the rest: the real form Cl(p, 10-p), the V117
    quaternionic table."""
    h = [np.asarray(m, dtype=complex) for m in clifford.generators()]
    return [A if j < p else 1j * A for j, A in enumerate(h)]


def planes_of(hp):
    """The 45 plane generators T'_ab = i Gamma'_a Gamma'_b, a < b, in the
    sm_qca order."""
    return [1j * (hp[a] @ hp[b])
            for a in range(10) for b in range(a + 1, 10)]


def signature_step(p, g):
    """The copied 4-line step pattern at signature p and clock g: the
    ordered product of U_ab = cos(g) I - i sin(g) T'_ab."""
    import math
    n = 32
    In = np.eye(n, dtype=complex)
    U = In.copy()
    for T in planes_of(swapped_generators(p)):
        U = U @ (math.cos(g) * In - 1j * math.sin(g) * T)
    return U


def unitarity_err(U):
    return float(np.abs(U @ U.conj().T - np.eye(U.shape[0],
                                                dtype=complex)).max())


def intertwiner(hp, tol=1e-9):
    """The Q of the antilinear-commutant stack (V117 fs_indicator
    algebra), returned together with the type and J^2."""
    n = hp[0].shape[0]
    In = np.eye(n, dtype=complex)
    Mat = np.vstack([np.kron(A.T, In) - np.kron(In, np.conj(A))
                     for A in hp])
    _, s, vt = np.linalg.svd(Mat, full_matrices=False)
    null = vt[s < tol]
    if null.shape[0] == 0:
        return None, None
    C = null[0].reshape(n, n)
    if np.linalg.matrix_rank(C) < n:
        return None, None
    U, _, Vh = np.linalg.svd(C)
    return U @ Vh, C


def closure_census(p):
    """For every g = pi*a/b, 0 < a < b <= DEN_MAX: smallest k <= K_MAX
    with ||U(g)^k -/+ I||_max < CLOSE_EPS."""
    seen = {}
    found = {}
    worst_unit = 0.0
    for b in range(2, DEN_MAX + 1):
        for a in range(1, b):
            g = np.pi * a / b
            key = round(g, 15)
            if key in seen:
                continue
            seen[key] = (a, b)
            U = signature_step(p, g)
            unit = unitarity_err(U)
            worst_unit = max(worst_unit, unit)
            Uk = np.eye(32, dtype=complex)
            for k in range(1, K_MAX + 1):
                Uk = Uk @ U
                dev = min(float(np.abs(Uk - np.eye(32, dtype=complex))
                                .max()),
                          float(np.abs(Uk + np.eye(32, dtype=complex))
                                .max()))
                if dev < CLOSE_EPS:
                    found["pi*%d/%d" % (a, b)] = {"k": k, "dev": dev}
                    break
    return found, worst_unit, len(seen)


def main():
    checks, data = [], {}

    # ---- T0: every booked quaternionic root reads FS = -1 ---------------
    t0_ok = True
    types = {}
    for p in PS:
        hp = swapped_generators(p)
        typ, j2 = fs_indicator(hp)
        types[str(p)] = {"type": typ, "J2": (None if j2 is None
                                            else float(j2))}
        if typ != "quaternionic" or j2 is None or abs(j2 + 1) > 1e-6:
            t0_ok = False
    data["t0_fs"] = types
    checks.append({"name": "T0 FS indicator of every booked quaternionic "
                          "root p in {}: {}".format(
                              PS, "all quaternionic, J^2 = -1"
                              if t0_ok else types),
                   "pass": bool(t0_ok)})
    if not t0_ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V311 MACHINERY-DRIFT: a booked quaternionic root no longer "
              "reads FS = -1; nothing is interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- per signature: CP test, unitarity, closure census --------------
    signatures = {}
    any_unitary = False
    for p in PS:
        hp = swapped_generators(p)
        Q, _C = intertwiner(hp)
        U = signature_step(p, G_CLOSED)
        # V117 T3 pattern: antiunitary image of the ordered step against
        # the inverse
        cp_dev = (None if Q is None else
                  float(np.abs(Q.conj().T @ U.conj() @ Q
                               - np.linalg.inv(U)).max())
                  if np.linalg.matrix_rank(U) == 32 else "singular-step")
        unit = unitarity_err(U)
        squares = [float(np.real(np.trace(hp[j] @ hp[j]))) / 32.0
                   for j in range(10)]
        plane_sq = {}
        for T in planes_of(hp):
            key = round(float(np.real(np.trace(T @ T)) / 32.0), 6)
            plane_sq[key] = plane_sq.get(key, 0) + 1
        found, worst_unit, n_clocks = closure_census(p)
        is_unitary = unit < TOL
        any_unitary = any_unitary or is_unitary
        signatures[str(p)] = {
            "module": "Cl(%d,%d)" % (p, 10 - p),
            "square_signs": squares,
            "plane_T_sq_counts": plane_sq,
            "cp_antiunitary_dev": cp_dev,
            "unitarity_err_at_closed_clock": unit,
            "unitarity_worst_over_closure_grid": worst_unit,
            "closure_clocks": found,
            "n_closure_clocks_tested": n_clocks,
            "verdict": ("STEP-UNITARY" if is_unitary
                        else "STEP-HYPERBOLIC-DECLARED")}
        checks.append({"name": "p = {:d} (Cl({:d},{:d})): unitarity "
                              "residual at the closed clock {:.3e}, "
                              "worst over the {}-clock rational grid "
                              "{:.3e}; CP antiunitary deviation {}; "
                              "closure clocks: {}".format(
                                  p, p, 10 - p, unit, worst_unit,
                                  n_clocks, cp_dev,
                                  found if found else "none"),
                       "pass": True})
        print("  p={:d}: unitary_err {:.6e} (closed clock), worst {:.6e}, "
              "CP dev {}, closures: {}".format(
                  p, unit, worst_unit, cp_dev,
                  found if found else "none"), flush=True)
    data["signatures"] = signatures

    # ---- branch ----------------------------------------------------------
    branch = ("STEP-UNITARY" if any_unitary
              else "STEP-HYPERBOLIC-DECLARED")
    data["status"] = branch
    data["unitary_p"] = sorted(int(p) for p, s in signatures.items()
                               if s["verdict"] == "STEP-UNITARY")
    checks.append({"name": "branch {} over the quaternionic class p in "
                          "{} (unitarity band {:.0e} = "
                          "tolerances.strict)".format(branch, PS, TOL),
                   "pass": True})
    print("V311 {}: all p in {} - unitarity residuals {}; the "
          "quaternionic signature swap {} the step.".format(
              branch, PS,
              {p: "{:.3e}".format(s["unitarity_err_at_closed_clock"])
               for p, s in signatures.items()},
              "survives on p in {}".format(data["unitary_p"])
              if any_unitary else "does not survive (the step is "
                                  "non-unitary at every quaternionic "
                                  "signature)"), flush=True)
    _emit(data, checks, ok=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V311_quaternionic_root_step", "tier": 3,
           "property": "does the automaton step survive a swap to the "
                       "quaternionic real forms of Cl_10(C) "
                       "(tier 3, O5/O8)",
           "residual": True,
           "programme": "pillar-surgery V311",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v311_quaternionic_root_step.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V311 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
