#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 — the TBM regression test (the W6 mixing, honest).

**Pourquoi ce module existe.**  La version précédente de
`cl10.family.tribimaximal_3x3` contenait un **signe errone** dans `r23` : elle
produisait la matrice **bimaximale**, qui a les **mêmes modules carrés** que la
vraie TBM.  Tous les tests existants vérifiaient les **angles** (les modules) et
la **unitarité** — deux propriétés que la matrice fausse satisfaisait.  Le bug
était donc **invisible** : la matrice ne diagonalisait pas la base démocratique
`|1><1|`, et aucun test ne le vérifiait.

Ce module ferme ce trou.  Il vérifie, dans l'ordre :

1. **la diagonalisation** — `U^dagger |1><1| U` est diagonale (le test absent) ;
2. les **angles** TBM (1/3, 1/2, 0) ;
3. que la **matrice fausse** (le signe erroné restauré) **échoue** au test 1
   tout en **réussissant** les tests 2 et 3 — la preuve que le garde-fou est
   nécessaire, pas redondant ;
4. le **partage honnête** du secteur de mélange : `sin^2 12 = 1/3` est dérivé
   0-paramètre (le vecteur propre lourd est unique), `sin^2 23` et `sin^2 13`
   balayent `[0, 2/3]` sous la liberté O(2) du sous-espace 2-fold dégénéré ;
5. que le **CKM** calculé à partir de `H_u = H_d` démocratiques est l'**identité**.

Aucun de ces points n'est postulé : chacun est **calculé** ici.
"""

import os
import sys
import json
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from cl10 import family, yukawa_3x3

TOL = 1e-9


def _tbm_with_the_old_sign_error():
    """La matrice fausse : le signe erroné de `r23` restauré (bimaximale)."""
    s12, c12 = np.sqrt(1.0 / 3.0), np.sqrt(2.0 / 3.0)
    s23 = c23 = np.sqrt(0.5)
    r12 = np.array([[c12, s12, 0.0], [-s12, c12, 0.0], [0.0, 0.0, 1.0]])
    r23 = np.array([[1.0, 0.0, 0.0], [0.0, c23, s23], [0.0, -s23, c23]])   # the old sign
    return r23 @ r12


def _offdiag(U, J):
    D = U.conj().T @ J @ U
    return float(np.abs(D - np.diag(np.diag(D))).max())


def tbm_regression():
    """The full honest W6 mixing verdict, computed (not asserted)."""
    J = np.ones((3, 3), dtype=complex)

    good = family.tribimaximal_3x3()
    bad = _tbm_with_the_old_sign_error()

    # (1) diagonalization — the check the old suite lacked
    good_diag = _offdiag(good, J)
    bad_diag = _offdiag(bad, J)

    # (2)/(3) the angles and the unitarity: BOTH matrices pass these
    def angles(U):
        return (abs(U[0, 1]) ** 2, abs(U[1, 2]) ** 2, abs(U[0, 2]) ** 2)

    def unitary_err(U):
        return float(np.abs(U.conj().T @ U - np.eye(3)).max())

    good_angles, bad_angles = angles(good), angles(bad)
    target = (1.0 / 3.0, 0.5, 0.0)
    good_angles_ok = all(abs(a - t) < 1e-12 for a, t in zip(good_angles, target))
    bad_angles_ok = all(abs(a - t) < 1e-12 for a, t in zip(bad_angles, target))
    good_unitary_ok = unitary_err(good) < TOL
    bad_unitary_ok = unitary_err(bad) < TOL

    # (4) the honest split of the mixing sector
    honest = yukawa_3x3.mixing_from_democratic()

    # (5) the CKM from identical democratic matrices
    ckm = yukawa_3x3.quark_ckm()

    # the regression guard itself: the good matrix passes, the bad one fails
    guard_catches_the_bug = bool(good_diag < TOL and bad_diag > 1e-3)

    ok = bool(
        good_diag < TOL
        and good_angles_ok
        and good_unitary_ok
        and guard_catches_the_bug
        and bad_angles_ok and bad_unitary_ok      # the old tests were blind to it
        and honest["sin2_12_fixed"]
        and abs(honest["sin2_12_value"] - 1.0 / 3.0) < 1e-9
        and honest["sin2_23_range"][1] > 0.5      # genuinely not fixed
        and honest["sin2_13_range"][1] > 0.5
        and honest["tbm_diagonalizes_democratic"]
        and honest["tbm_in_the_O2_family"]
        and ckm["ckm_is_identity"]
    )

    return {
        "backend": "numpy",
        "check_1_diagonalization": {
            "good_TBM_max_offdiag": good_diag,
            "old_sign_error_max_offdiag": bad_diag,
            "good_diagonalizes": bool(good_diag < TOL),
            "old_sign_error_diagonalizes": bool(bad_diag < TOL),
        },
        "check_2_angles": {
            "good": [float(x) for x in good_angles],
            "old_sign_error": [float(x) for x in bad_angles],
            "target": [1/3, 1/2, 0.0],
            "good_ok": good_angles_ok,
            "old_sign_error_ok": bad_angles_ok,
        },
        "check_3_unitarity": {
            "good_err": unitary_err(good),
            "old_sign_error_err": unitary_err(bad),
            "good_ok": good_unitary_ok,
            "old_sign_error_ok": bad_unitary_ok,
        },
        "check_4_guard_is_necessary": {
            "guard_catches_the_bug": guard_catches_the_bug,
            "note": ("the old suite checked the ANGLES and the UNITARITY, both of "
                     "which the wrong-sign matrix satisfies; only the "
                     "DIAGONALIZATION test catches it"),
        },
        "check_5_honest_mixing_split": {
            "sin2_12_derived_0param": honest["sin2_12_fixed"],
            "sin2_12_value": honest["sin2_12_value"],
            "sin2_12_vs_observed": honest["sin2_12_vs_observed"],
            "sin2_23_O2_range": honest["sin2_23_range"],
            "sin2_13_O2_range": honest["sin2_13_range"],
            "tbm_is_in_the_O2_family": honest["tbm_in_the_O2_family"],
            "s3_transpositions": {k: v["reproduces_TBM"]
                                  for k, v in honest["s3_transpositions_equivalent"].items()},
        },
        "check_6_ckm": {
            "ckm_is_identity": ckm["ckm_is_identity"],
            "note": ("H_u = H_d both democratic -> the computed CKM is the identity; "
                     "a non-trivial CKM needs an up/down misalignment (the W6 residual)"),
        },
        "ok": ok,
    }


def main():
    t0 = time.time()
    r = tbm_regression()
    ok = r["ok"]

    verdict = {
        "id": "V60_tbm_regression",
        "tier": 2,
        "property": "P2",
        "residual": True,
        "pass": ok,
        "status": ("the_TBM_diagonalizes_the_democratic_base_(the_check_the_old_suite_lacked); "
                   "the_wrong_sign_matrix_is_caught; sin2_12=1/3_derived_0_parameter; "
                   "sin2_23_and_sin2_13_sweep_[0,2/3]_(O2_choice); the_CKM_from_identical "
                   "democratic_matrices_is_the_identity"),
        "tol": TOL,
        "data": r,
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v60_tbm_regression.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 72)
    print("V60 - the TBM regression test (the W6 mixing, honest)")
    print("=" * 72)
    c1 = r["check_1_diagonalization"]
    print("  [1] DIAGONALIZATION of |1><1| (the check the old suite lacked):")
    print(f"      the fixed TBM   : max off-diagonal = {c1['good_TBM_max_offdiag']:.2e}  -> {c1['good_diagonalizes']}")
    print(f"      the old sign bug: max off-diagonal = {c1['old_sign_error_max_offdiag']:.4f}  -> {c1['old_sign_error_diagonalizes']}")
    c2, c3 = r["check_2_angles"], r["check_3_unitarity"]
    print(f"  [2] the ANGLES (1/3, 1/2, 0): fixed TBM {c2['good_ok']}   old sign bug {c2['old_sign_error_ok']}")
    print(f"  [3] the UNITARITY          : fixed TBM {c3['good_ok']}   old sign bug {c3['old_sign_error_ok']}")
    print(f"  [4] -> the old tests were BLIND to the bug; the diagonalization guard catches it: "
          f"{r['check_4_guard_is_necessary']['guard_catches_the_bug']}")
    c5 = r["check_5_honest_mixing_split"]
    print("  [5] the honest mixing split:")
    print(f"      sin^2_12 = {c5['sin2_12_value']:.4f}  DERIVED 0-parameter "
          f"(observed {c5['sin2_12_vs_observed']['observed']}, ratio {c5['sin2_12_vs_observed']['ratio']:.3f})")
    print(f"      sin^2_23 sweeps {tuple(round(x,3) for x in c5['sin2_23_O2_range'])}  (O(2) freedom, NOT derived)")
    print(f"      sin^2_13 sweeps {tuple(round(x,3) for x in c5['sin2_13_O2_range'])}  (O(2) freedom, NOT derived)")
    print(f"      the S3 transpositions reproducing TBM: {c5['s3_transpositions']}")
    print(f"  [6] the CKM from H_u = H_d (both democratic) is the identity: "
          f"{r['check_6_ckm']['ckm_is_identity']}")
    print(f"  VERDICT:  {'PASS (the TBM diagonalizes the democratic base; the mixing split is honest; the guard catches the old bug)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v60_tbm_regression.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
