# Audit bibliographique — doc/publication.md (52 références)

Date : 2026-10-06. Méthode : chaque entrée vérifiée contre Crossref (DOI direct
ou query titre+auteur+volume+page), INSPIRE HEP et l'API arXiv. Un DOI n'est
retenu que si titre, auteurs, journal, volume, page et année correspondent.

## Corrections appliquées (erreurs réelles trouvées)

| # | Erreur initiale | Correction | Preuve |
|---|---|---|---|
| 4 | Arrighi/Nesme/Scherrer, "An introduction to QCA", Nat. Comput. 11, 546 (2012) | P. Arrighi, "An overview of quantum cellular automata", Nat. Comput. **18**, 885 (2019) | DOI 10.1007/s11047-019-09762-6 (Crossref : titre/vol/page exacts) |
| 5 | D'Ariano/Perinotti/Bahreyni/Ford, "Equivalence of QCA and circuits", JCP 390, 198 (2019) — n'existe pas sous cette forme | D'Ariano & Perinotti, Front. Phys. **12**, 120301 (2016) + Bisio et al., Found. Phys. **45**, 1203 (2015) | DOI 10.1007/s11467-016-0616-z ; INSPIRE |
| 7 | Fritzsch & Minkowski, "Unified interactions of strong and electroweak forces of quarks and leptons" | titre exact : "Unified interactions of leptons and hadrons", Ann. Phys. 93, 193 (1975) | DOI 10.1016/0003-4916(75)90211-0 |
| 19 | DOI erroné (91190-4 = erratum) | 10.1016/0370-2693(81)91026-1 = le no-go lui-même, PLB 105, 219 | Crossref |
| 23 | DOI erroné (BF01168405) | 10.1007/BF01647974, Manuscripta Math. 28, 235 (1979) | Crossref |
| 24 | DOI Elsevier inventé pour Sakharov | 10.1070/PU1991v034n05ABEH002497 (Sov. Phys. Usp. 34, 392, traduction anglaise) | INSPIRE (5917 citations) |
| 29 | 't Hooft NPB 256 727 : DOI absent | ajouté 10.1016/0550-3213(85)90418-3 | Crossref |
| 31 | "Leptogenesis with quasi-degenerate neutrinos", PRL 110, 191601 (Alloul et al.) attribué à tort comme référence unique | Blanchet & Di Bari, "The minimal scenario of leptogenesis", NJP 14, 125012 (2012) ; le papier quasi-dégénéré LHC est Alloul et al., EPJC 73, 2469 (2013) — remplacé par Pilaftsis & Underwood, "Resonant leptogenesis", NPB 692, 303 (2004) [40] | Crossref pour les deux |
| 33 | "M. Regge" | **T.** Regge, "General relativity without coordinates", Il Nuovo Cimento 19, 558 (1961) | DOI 10.1007/BF02733251 |
| 34 | Gundlach et al., "Constraint damping and the 3+1 formulation", PRD 68, 124025 (2003) — auteurs et référence inexacts | Gundlach, Martin-Garcia, Calabrese, Hinder, CQG **22**, 3767 (2005) | DOI 10.1088/0264-9381-22-17-025 |
| 36 | 't Hooft "Dimensional reduction of quantum gravity" (Plenum 1989) — introuvable | supprimé ; gardé Susskind, JMP 36, 6377 (1995) | DOI 10.1063/1.531249 |
| 41 | Claudson/Turner/Whisnant (BN+masse des neutrinos) — ne correspond pas à l'argument freeze-out cité | Kolb & Turner, *The Early Universe* (Addison-Wesley, 1990), Ch. 5 | traité classique |
| 46 | Hisano/Murayama/Yanagida NPB 374, 147 (1992) — référence mal citée | Hisano, Kobayashi, Muramatsu, Nagata, PLB **724**, 283 (2013) (facteurs A_L deux boucles) | DOI 10.1016/j.physletb.2013.06.030 |
| 48 | Kaplan "in Aspects of the Chiral Anomaly (World Scientific)" — chapitre inventé | Kaplan PLB 288, 342 (1992) + Shamir NPB 406, 90 (1993) + Narayanan-Neuberger PLB 417, 141 (1998) | Crossref |
| 49 | Machacek & Vaughn "II. Scalar model" — le papier II est le modèle scalaire pur ; le texte cite le calcul à 2 boucles général | I. Proof of the equivalence theorem, NPB 222, 83 (1983) + III. Scalar model with fermions, NPB 249, 70 (1985) | Crossref |
| 50 | de Gouvêa/Jazwinski/Rodejohann arXiv:2311.12798 "Future neutrino oscillations II" | remplacé par Esteban et al., **NuFIT 6.0**, PoS ICHEP2024, 446 (2025), arXiv:2410.05380 (cohérent avec [17]) | arXiv |

## Chiffres corrigés (liés aux refs)

- Borne Super-K p→μ⁺π⁰ : 7.0×10³³ yr → **1.6×10³⁴ yr** (Takenaka et al. 2020,
  PRD 102, 112011 [15]). Affecte : légende Fig. 22, make_figures.py (fig22).
- Borne p→e⁺K⁺ : 2.7×10³³ yr était invérifiable (aucune recherche SK dédiée
  e⁺K⁺). Fig. 22 utilise désormais la borne **p→ν̄K⁺ : 6.6×10³³ yr**
  (Abe et al. 2014, PRD 90, 072005 [52]), canal mesuré le plus proche, avec
  mention explicite dans la légende.
- Borne e⁺π⁰ : 2.4×10³⁴ yr confirmée ; make_figures.py fig20/fig21 utilisaient
  1.6×10³⁴ (ancienne borne) → corrigé à 2.4×10³⁴.
- Échelle GUT exclue : 5×10¹⁵ → **6×10¹⁵ GeV** (recalculée avec la borne
  2.4×10³⁴ : M = M_P·(τ_bound/τ₀)^¼ ≈ 6×10¹⁵). Affecte : §11.4, légendes
  Fig. 20/21, annotations fig20.
- KATRIN 0.45 eV (90% CL, 259 jours) : **confirmé** (Science 388, 180 (2025)).

## Références ajoutées (couverture des affirmations non sourcées)

[33] Regge (gravité discrète, §9.2) · [34] Gundlach Z4 (§9.3) ·
[35] Bekenstein/Hawking (§13.1) · [36] Susskind holographie (§13.1) ·
[37] Weinberg problème de la CC (§12.1) · [38] Morningstar-Peardon σ≈0.18 GeV²
(§13.2) · [39] Barbieri et al. + Buchmüller et al. (estimation leptogenèse,
§11.3) · [40] Pilaftsis-Underwood leptogenèse résonante (§11.3) ·
[41] Kolb-Turner (formule freeze-out, §12.2) · [42] Hall et al. freeze-in
(§12.2, phrase ajoutée) · [43] Fritzsch matrice démocratique (§11.5) ·
[44] Wolfenstein λ=0.225 (§11.5) · [45] Jungman et al. WIMP 100 GeV (§10.2) ·
[46] Hisano et al. A_L (§11.4) · [47] Tonomura fentes d'Young électron unique
(§7) · [48] Kaplan/Shamir/Neuberger fermions chiraux sur réseau (§8.4) ·
[49] Machacek-Vaughn (§11.1) · [50] NuFIT 6.0 (§11.5, Table 4) ·
[51] Nath-Fileviez Pérez (§11.4) · [52] SK ν̄K⁺ (Fig. 22).

Toutes les 52 entrées sont citées dans le corps du texte (aucune référence
orpheline). 46 DOIs vérifiés dans livrable/doi.json ; les 6 sans DOI sont des
livres ([9], [11], [18], [41]), un arXiv-only ([30]) et un compte-rendu PoS
avec arXiv ([50]).

## Fichiers touchés

- doc/publication.md : bibliographie réécrite (32→52), 23 citations insérées
  dans le texte, bornes/échelles corrigées.
- make_figures.py : bornes fig20/21/22 + annotations ; figs 20-22 régénérées.
- livrable/doi.json : 46 DOIs vérifiés.
- livrable/publication.{tex,pdf} : recompilé (55 p., 52 bibitems, 46 DOIs,
  aucun avertissement de citation non résolue).
