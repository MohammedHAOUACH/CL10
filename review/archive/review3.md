# Review3 — vérification indépendante des corrections de `review2_response.md`

**Date :** 6 octobre 2026. **Objet :** contrôle que les corrections annoncées dans
`review2_response.md` sont **réellement présentes** dans `doc/publication.md`
(2053 lignes), et re-exécution de `run_validation.py`.

**Méthode :** recherche littérale dans le fichier (pas lecture de la réponse),
re-exécution des modules, recalcul des nombres. Une affirmation d'application n'est
pas une vérification.

**Verdict en une phrase :** les **10 corrections sont bien appliquées** et le suite
**60/60 passe** (re-exécuté) ; il reste **une contradiction réelle non résolue**
(papier rank-1 vs code rank-3) et **deux points de transparence** (le `n` optimal des
données, et `28/79` étiqueté 0-paramètre).

---

## 1. Les 10 corrections : présentes, vérifiées ligne par ligne

| # | Correction annoncée | Présent ? | Emplacement vérifié |
|---|---|---|---|
| 1 | `α₁₀ = g₁₀² = 1/25` (pas `g₁₀²/4π`) | **oui** | L61, L1517 |
| 2 | base democratic = **rank-1** `|1><1|` | **oui** | L598, L1376, L1619, L1654 |
| 3 | convention unique **Q = T₃ + Y** | **oui** | L198, L280 |
| 4 | B–L = **generator**, singlet = **(1,−1)** | **oui** | L234, L246, L550, L558 |
| 5 | Landau pole **1.4×10³⁹ GeV** | **oui** | L1477, L1484 |
| 6 | `n_qubits = 5` (pas `n_gen = 5`) | **oui** | L1330, L1642, L1901 |
| 7 | « first computed **kinematic** object » | **oui** | L1347, L1840 |
| 8 | abstract : « strongest defensible sense » | **oui** | abstract |
| 9 | « one equation, two projections » | **oui** | L1597, L1899, L1904 |
| 10 | refs 11–13 Clifford GUT (Castro) | **oui** | bibliographie |

Note d'honnêteté : quatre de ces points avaient d'abord été signalés comme manquants
par ma recherche — c'étaient des **faux négatifs de regex** (unicode `T₃` vs `T_3`,
`10^{39}` vs `1.4×10³⁹`, markdown `*kinematic* object`). Les corrections sont bien là.

**Vérification du point 5 par exécution :** `cl10.unification.landau_pole_u1()` →
`1.385249013068412e+39` GeV. Le chiffre du papier est la valeur exécutée, pas une estimation.

**Vérification du point 8 :** aucune occurrence de « TOE » ni « Theory of Everything »
dans `doc/publication.md`. Le titre est « The Standard Model as a Single Clifford Root ».
Cohérent avec `results/validation_summary.json`, qui ne contient **plus** de clé `toe`
(alors que review2 citait `toe: true`). C'est le bon état.

**Vérification du suite :** `run_validation.py` → **60/60 PASS**. Re-exécuté, pas cru.

---

## 2. Le défaut réel restant : le papier dit rank-1, le code calcule rank-3

Le papier a été corrigé (rank-1, spectre `{3,0,0}`, une génération lourde). Le code ne l'a
pas été, et **deux modules s'opposent** :

- `cl10/yukawa.py` L34 : `Y = np.eye(n_gen)` → eigenvalues `(1,1,1)`, **rank 3**
  - L11 (docstring) : « the mass matrix is **non-singular** (rank 3) »
  - L42 : `"non_singular": bool(rank == n_gen)`
- `cl10/yukawa_3x3.py` L170 : `M = np.ones((n, n))` → `|1><1|`, **rank 1**, spectre `{3,0,0}`

Ce sont deux structures 0-paramètre **incompatibles** : `np.eye(3)` et `np.ones(3,3)` ne
peuvent pas être toutes deux « la » matrice de Yukawa 0-paramètre.

Conséquence concrète mesurée :

```
V28 as executed:            rank 3, non_singular True   -> PASS
V28 with democratic base:   rank 1, non_singular False  -> FAIL
```

Donc **V28 passe précisément parce qu'il utilise la structure que le papier a abandonnée.**
Le papier affirme une chose que le test vérifie par l'affirmation contraire.

Il y a de plus une tension physique interne au papier : L599/L1030 disent spectre `{3,0,0}`
(une génération massive, deux dégénérées à zéro), mais L1172 et le nom de V28 disent que la
structure donne « les masses des trois générations ». Un rank-1 ne donne pas trois masses.

**Correction proposée** (la plus simple et la plus honnête) :
1. dans `cl10/yukawa.py`, remplacer `Y = np.eye(n_gen)` par `Y = np.ones((n_gen, n_gen))` ;
2. remplacer le critère `non_singular` par le critère réel : `rank == 1` et spectre `{3,0,0}` ;
3. renommer V28 : « la structure de masse 0-paramètre est rank-1 (une génération lourde +
   deux dégénérées) » — et **ne plus** affirmer que les trois générations acquièrent une masse
   par cette structure, ce qui est le contenu du résidu nommé W6 ;
4. aligner la docstring L11.

Cela fait passer V28 de « passe parce que rank-3 » à « passe parce que rank-1 », ce qui est
ce que le papier dit.

---

## 3. Transparence : le `n` optimal des données est 6, pas 5

Le module `cl10/g10_selection.py` calcule lui-même :

```
best_n_1loop = 6        gap |1/36 - 0.0242| = 0.0036
n_qubits     = 5        gap |1/25 - 0.0242| = 0.0158
ratio_to_1loop = 1.6555
```

Le `n` le plus proche des données est **6**, soit **4.4× plus proche** que 5. Le papier dit
« right order, not exact » (L1517-1518) — c'est honnête — mais ne dit **jamais** que le scan
interne désigne 6. Le « nearest » de L1901 concerne un **autre** scan (`n = 1…20`, facteur
2.08 sur le rapport hiérarchie), pas le data-check de `α₁₀`.

Un referee verra `best_n_1loop = 6` dans le JSON et ne le trouvera pas dans le texte. Il faut
l'écrire explicitement : `n = 5` est un **nombre du root** (le module 32-D, `10 = 2×5`), pas
le meilleur ajustement ; le meilleur ajustement est 6, et c'est justement ce qui fait de la
lecture `1/n_qubits` une **convention de mapping** et non une prédiction. C'est cohérent avec
le cadrage déjà en place (L1350) — il manque juste la phrase.

---

## 4. Transparence : `28/79` n'est pas un nombre du root

`28/79` est le **facteur sphaleron du SM** (il dépend du contenu en degrés de liberté du SM,
donc du fait que le SM est dérivé). `cl10/eta_cp.py` L114-115 l'appelle « 0-parameter FACTORS ».

Le produit est correct et vérifié :

```
(28/79) × (g₀⁴/16π) × (1/16π)² = 6.813e-10   vs observé 6.1e-10   ratio 1.12
```

Mais le statut épistémique est différent des autres facteurs : `g₀ = 1/8` et `1/16π` sont des
nombres du root, `28/79` est un nombre **du secteur SM**. Défendable seulement si le contenu SM
est lui-même dérivé (V23/V26). Sinon c'est un input. Le papier devrait distinguer les deux,
comme il le fait déjà ailleurs — c'est une phrase, pas un changement de résultat.

---

## 5. Ce qui est bon et mérite d'être gardé tel quel

- La discipline de provenance : « 0-parameter in its **algebra** », « named residual »,
  « checks internes ≠ tests contre la nature », « falsifiable ». C'est la partie la plus
  solide du papier, et elle est appliquée de façon cohérente (337 occurrences de
  « 0-parameter », 28 de « named residual », 31 de « convention »).
- Le cadrage « one equation, two projections » est correct et protège contre la double compte.
- L'absence de claim TOE, y compris dans le JSON de validation.

---

## 6. Bilan

**Statut :** livrable sur la forme et sur l'honnêteté. Les 10 corrections de review2 sont
appliquées et le suite passe.

**Ce qui empêche encore d'appeler le papier « complet » :**

1. **la contradiction rank-1 / rank-3 entre `doc/publication.md` et `cl10/yukawa.py`** —
   c'est un défaut technique réel, V28 ne teste pas ce que le papier affirme ;
2. le data-check de `α₁₀` doit nommer `n = 6` comme l'optimum des données ;
3. `28/79` doit être étiqueté input du secteur SM, pas nombre du root.

Les points 2 et 3 sont des phrases. Le point 1 est une correction de code et de critère de
test — c'est le seul qui change un verdict.

---

## 7. Corrections appliquées (suite de l'audit, 6 octobre 2026)

Les trois points ci-dessus ont été **corrigés**, et l'audit a révélé au passage un
**quatrième défaut**, plus grave, que review2 n'avait pas vu.

### 7.1 Point 1 — rank-1 / rank-3 : corrigé

`cl10/yukawa.py` réécrit : `Y = np.ones((3,3))` = `|1><1|`, **rank 1**, spectre `{3,0,0}`.
`np.eye(3)` a disparu du projet. Le critère `non_singular` est remplacé par les critères
réels : `rank == 1`, `n_distinct_levels == 2`, `n_heavy == 1`, `n_light_2fold == 2`,
`democratic`, `singular`. V28 et V29 alignés.

V28 passe maintenant **parce que la base est rank-1**, ce que le papier affirme.

### 7.2 Point 2 — `n = 6` : corrigé

`cl10/g10_selection.py` expose `candidate_scan_gap`, `gap_root_n5`, `gap_best_n`,
`data_check_selects_n = False`. Le papier (§15), `doc/api.md`, `doc/reference.md` et V56
disent explicitement : l'optimum des données est **n = 6** (gap 0.15), **pas** le `n = 5`
du root (gap 0.66) ; la lecture `1/n_qubits` est une **convention de mapping**, pas un fit.

### 7.3 Point 3 — `28/79` : corrigé

Étiqueté partout comme **input du secteur SM** (le comptage des d.o.f. gauches), nombre du
root seulement **conditionnel à V23**. Le produit est renommé « produit de facteurs de
**trois provenances** » : root (`g₀ = 1/8`, `1/16π`), QFT (`ε_CP`), SM (`28/79`).
Valeur inchangée et revérifiée : `6.813e-10`, ratio **1.12**.

### 7.4 Point 4 — **découvert par l'audit** : la « TBM » du code n'était pas la TBM

En vérifiant la chaîne « base démocratique → TBM » (celle qui devait soutenir le point 1),
l'audit a trouvé que `cl10.family.tribimaximal_3x3()` contenait une **erreur de signe** dans
`r23` : elle produisait la matrice **bimaximale**, pas la TBM de Harrison-Perkins-Scott.

Le point grave n'est pas le signe, c'est que **le suite ne pouvait pas le voir**. Mesuré par
`experiments/V60_tbm_regression/` :

```
[1] DIAGONALISATION de |1><1|  (le test qui manquait)
    TBM corrigée  : max off-diagonal = 1.01e-16  -> True
    ancien signe  : max off-diagonal = 1.1547    -> False
[2] les ANGLES (1/3, 1/2, 0) : corrigée True    ancien signe True
[3] l'UNITARITE              : corrigée True    ancien signe True
[4] -> les anciens tests étaient AVEUGLES au bug ; seul le test de diagonalisation l'attrape
```

La matrice fausse avait les **bons modules** et était **unitaire** : elle passait V39, V45 et
V29. Elle ne diagonalisait tout simplement pas la base démocratique. **V60 est ajouté** comme
test de régression, et V39, V45 et `w6_verdict()` vérifient désormais la diagonalisation.

### 7.5 Point 4 (suite) — la conséquence physique : le partage honnête du mélange

La base démocratique `|1><1|` a un vecteur propre lourd **unique** `(1,1,1)/√3` et un
sous-espace léger **2-fold dégénéré**. Tout diagonalisateur est donc une rotation O(2) dans ce
sous-espace. Balayage sur 1441 rotations :

```
sin²θ12 = 1/3        FIXÉ 0-paramètre (le recouvrement avec le vecteur lourd unique)
sin²θ23  balaie [0, 2/3]      liberté O(2) — NON dérivé
sin²θ13  balaie [0, 2/3]      liberté O(2) — NON dérivé
les transpositions S3 qui reproduisent la TBM : (23) True, (12) False, (13) False
   or (12), (13), (23) sont trois symétries ÉQUIVALENTES de |1><1| : la base ne choisit pas (23)
le CKM calculé avec H_u = H_d démocratiques : l'IDENTITÉ, pas la TBM
```

**Théorème 13 réécrit** : seul `sin²θ12 = 1/3` est dérivé 0-paramètre. `θ23` et `θ13` sont un
**choix O(2)** (un input nommé). Le CKM non trivial exige un **misalignement** up/down que la
base symétrique ne fournit pas — c'est le résidu nommé W6. Le papier contient maintenant la
correction explicite de l'énoncé antérieur (« the CKM is the TBM », faux).

### 7.6 Ce que l'audit n'a PAS corrigé, parce que c'est une limite, pas un bug

- `sin²θ23` et `sin²θ13` ne sont pas prédits 0-paramètre. C'est un **théorème** de la
  construction (dégénérescence 2-fold), pas un défaut à réparer.
- Le CKM non trivial n'est pas obtenu sans brisure.
- `n = 5` n'est pas sélectionné par les données.

Ces trois faits sont maintenant **dits dans le papier**, avec les nombres qui les établissent.

### 7.7 Suite de validation finale

```
run_validation.py  ->  61/61 PASS   (exit 0, backend cupy, GPU)
```

`V60_tbm_regression` est le 61ᵉ. Compteurs mis à jour dans `doc/publication.md`,
`doc/reference.md` et `README.md` (60 → 61).

### 7.8 Verdict final

**Le papier est livrable.** Les trois points ouverts de la section 6 sont corrigés, le
quatrième défaut — celui que le suite rendait invisible — est corrigé **et** couvert par un
test qui aurait échoué avant. Les claims restants sont exacts dans leur provenance :
`sin²θ12 = 1/3` est dérivé, `θ23`/`θ13` sont un choix, le CKM est un résidu, `28/79` est un
input SM, `n = 5` est une convention de mapping.

La leçon méthodologique, à retenir : **vérifier les angles d'une matrice de mélange ne vérifie
pas la matrice.** Le test manquant était la diagonalisation de la base qu'elle est censée
diagonaliser.
