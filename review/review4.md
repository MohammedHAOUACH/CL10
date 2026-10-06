# Review 4 — audit final : base démocratique, mélange, provenance, et la figure hydrogène

**Date :** 6 octobre 2026.
**Objet :** vérifier par exécution que les problèmes ouverts à la fin du review3 sont réglés,
et auditer la nouvelle figure hydrogène (Fig. 26) et son expérience V61.

**Méthode :** recherche littérale dans les fichiers, re-exécution des modules, recalcul des
nombres. Une affirmation d'application n'est pas une vérification.

**Historique :** ce cycle fait suite à `review2.md` (la critique), `review2_response.md` (la
réponse) et `review3.md` (la vérification indépendante). Les trois sont conservés dans
`review/archive/` pour la trace ; celui-ci est le seul document actif.

---

## Verdict en une phrase

**62/62 PASS** (re-exécuté). Les quatre problèmes du review3 sont corrigés et vérifiés ;
l'audit en a découvert un **cinquième**, plus grave, que la suite de 60 tests ne pouvait pas
voir — il est corrigé **et** couvert par un test qui aurait échoué avant. La figure
hydrogène est intégrée avec une provenance honnête : **l'atome est un input, le photon est
CL10**.

---

## 1. La correction principale : la « TBM » du code n'était pas la TBM

C'est le défaut que review2 et review3 n'avaient pas vu. Il a été trouvé en vérifiant la
chaîne « base démocratique → TBM » censée soutenir la correction rank-1.

`cl10/family.py` contenait une **erreur de signe** dans `r23` : la matrice produite était
**bimaximale**, pas la TBM de Harrison-Perkins-Scott.

Le problème grave n'est pas le signe. C'est que **la suite ne pouvait pas le détecter**.
Mesuré par `experiments/V60_tbm_regression/` :

```
DIAGONALISATION de |1><1|   (le test qui manquait)
  TBM corrigée  : off-diagonal max = 1.01e-16  -> True
  ancien signe  : off-diagonal max = 1.1547    -> False

ANGLES (1/3, 1/2, 0) : TBM corrigée True    ancien signe True
UNITARITE            : TBM corrigée True    ancien signe True
```

La matrice fausse avait les **bons modules** et était **unitaire** : elle passait V39, V45 et
V29. Elle ne diagonalisait simplement pas la base qu'elle était censée diagonaliser.

**Ce qui a été fait :**
- `cl10/family.py` : la vraie TBM, avec un test de diagonalisation.
- `experiments/V60_tbm_regression/` : nouveau test de régression. Il **reconstruit le bug**
  (le signe restauré) et vérifie qu'il échoue au nouveau test tout en réussissant les anciens.
  C'est la preuve que le garde-fou est nécessaire, pas redondant.
- V39, V45 et `w6_verdict()` vérifient désormais la diagonalisation. V39 relancé :
  `tribimaximal_diagonalizes_democratic: True`.

**La leçon :** vérifier les angles d'une matrice de mélange ne vérifie pas la matrice.

---

## 2. Ce que la base démocratique fixe vraiment (Théorème 13 réécrit)

`|1><1|` a un vecteur propre lourd **unique** `(1,1,1)/√3` et un sous-espace léger
**2-fold dégénéré**. Tout diagonalisateur est donc une rotation O(2) dans ce sous-espace.
Balayage numérique sur 1441 rotations :

| quantité | statut | valeur mesurée |
|---|---|---|
| `sin²θ12` | **dérivé 0-paramètre** | `1/3` (fixé par le recouvrement avec le vecteur lourd unique) |
| `sin²θ23` | **non dérivé** | balaie `[0, 2/3]` |
| `sin²θ13` | **non dérivé** | balaie `[0, 2/3]` |
| CKM (`H_u = H_d`) | **l'identité** | pas la TBM |

Les transpositions S3 qui reproduisent la TBM : `(23)` True, `(12)` False, `(13)` False —
mais `(12)`, `(13)` et `(23)` sont trois symétries **équivalentes** de `|1><1|`. La base ne
choisit pas `(23)` : c'est un **input nommé** de plus.

Le papier contient maintenant la correction explicite de l'énoncé antérieur (« the CKM is the
TBM », faux) et la raison : il supposait des diagonaliseurs différents pour `H_u` et `H_d`,
ce qui est exactement le résidu.

**Ce n'est pas un bug réparable.** C'est un théorème de la construction. Il est dit, avec les
nombres.

---

## 3. Les trois points de transparence du review3 : corrigés et vérifiés

**Rank-1 / rank-3.** `cl10/yukawa.py` réécrit : `np.ones((3,3))` = `|1><1|`, rank 1, spectre
`{3,0,0}`. `np.eye(3)` a disparu du projet. Critères réels : `rank == 1`,
`n_distinct_levels == 2`, `n_heavy == 1`, `n_light_2fold == 2`, `democratic`, `singular`.
V28 passe maintenant **parce que la base est rank-1**, ce que le papier affirme.

**`n = 6`.** `cl10/g10_selection.py` expose `candidate_scan_gap`, `gap_root_n5`, `gap_best_n`,
`data_check_selects_n = False`. Le papier (§15), `doc/api.md`, `doc/reference.md` et V56 disent
explicitement : l'optimum des données est **n = 6** (gap 0.15), pas le `n = 5` du root
(gap 0.66). La lecture `1/n_qubits` est une **convention de mapping**, pas un fit.

**`28/79`.** Étiqueté partout comme **input du secteur SM** (le comptage des d.o.f. gauches),
nombre du root seulement **conditionnel à V23**. Le défaut restant trouvé ici : `doc/reference.md`
le listait encore dans « 0-parameter numbers (no fit anywhere) », en contradiction avec sa propre
table §9. Corrigé. Le produit est renommé « produit de facteurs de **trois provenances** » :
root (`g₀ = 1/8`, `1/16π`), QFT (`ε_CP`), SM (`28/79`). Valeur revérifiée : `6.813e-10`,
ratio **1.12**.

**Surévalués résiduels supprimés** (trouvés par recherche littérale, les corrections de code
laissent souvent le texte derrière elles) : en-tête api.md « product of 0-parameter numbers »,
« the 0-parameter quark CKM » (deux endroits), docstrings V29 et V39 « 0-parameter PMNS »,
`doc/reference.md` « the 0-parameter PMNS », un « all 53 experiments » périmé, et les compteurs
60 → 62.

---

## 4. La figure hydrogène (Fig. 26) : intégrée, avec sa provenance

**Ce que la figure montre.** L'atome d'hydrogène 3-D sur le **photon CL10** : le même champ
massless discret à 6 voisins qui porte le graviton (§6) et qui a fait le double-slit de
Fig. 25. À gauche le nuage `|ψ|²` en volume 3-D avec son isosurface — sphère en `1s`,
haltère à deux lobes en `2p`, beignet en `3d` ; la déformation **est** l'excitation. À droite
la coupe médiane. Le panneau d'angle lit l'état : populations, énergie en eV, fraction excitée,
dipôle, fréquence du photon.

**Ce qui est mesuré par V61, pas affirmé :**

1. **les niveaux sont les vrais niveaux hydrogène.** Gaps calculés : `1s→2p = 10.204 eV`,
   `1s→3p = 12.094 eV`, `2p→3d = 1.890 eV` — les différences de Rydberg à une erreur relative
   de **0**. Rien n'est ajusté aux données atomiques.

2. **les règles de sélection E1 sortent des fonctions d'onde sur le réseau.** Les éléments de
   matrice dipolaires `<i|z|j>` sont calculés en sommant les fonctions d'onde du réseau. Les
   paires interdites sortent à `1e-18` — `1s→3d` (`Δl = 2`) : `-6.8e-19`, `2p→3p` (`Δl = 0`) :
   `8.6e-18` — contre `9.93`, `-3.00`, `22.97` pour les autorisées. La règle n'est **pas
   imposée**, elle est mesurée.

3. **l'excitation n'est pas imposée.** Baseline : l'atome seul, sans photon, ne s'excite
   pas (fraction excitée `0.0`). Un photon tiré exactement sur le gap `1s→2p` est absorbé
   (fraction excitée **0.42**). L'absorption **suit la forme de raie** : à amplitude, largeur
   et durée fixes, la fraction excitée décroît avec le désaccord — `0.42`, `0.35`, `0.15`,
   `0.035`, `0.011` pour des désaccords de `0`, `0.063`, `0.125`, `0.25`, `0.50` rad/step.
   Suppression mesurée : **38×** entre la résonance et le désaccord lointain.

4. **l'atome excité raye et se désexcite** (Einstein A) : sans champ incident, la population
   `2p` tombe `1.000 → 0.201` pendant que `1s` monte `0.000 → 0.799`, et le champ émis est
   non nul — l'atome ne perd pas seulement de la population, il **source une onde sortante** à
   la fréquence du gap.

**La provenance honnête, écrite dans la légende.** L'**atome est un input** : ses niveaux et
ses formes sont les eigenstates réels de l'hydrogène, **non dérivés** du root Cl(10). Ce que le
root fournit est le **photon** — sa propagation, sa résonance avec l'atome, son émission. La
figure **démontre le secteur photon**, ce n'est **pas** une dérivation 0-paramètre de
l'hydrogène. Les trois boutons `A`, `h`, `γ_sp` sont des paramètres physiques de l'atome, **non
ajustés**.

---

## 5. Deux limites trouvées en auditant la figure (et non cachées)

**a) L'absorption hors résonance n'est pas nulle — et ce n'est pas un bug.**

Mon premier test exigeait `excited_fraction < 0.01` hors résonance. Il a échoué à **0.033**.
L'hypothèse « transitoire de lancement brutal » a été **testée et réfutée** : en balayant la
largeur du paquet (diagnostic sur une boîte de 80, sponge 8), l'excitation **monte** au lieu de
tomber — `0.020`, `0.153`, `0.514`, `0.667` pour des largeurs de 8, 16, 32, 64 cellules. La
cause mesurée est une autre : un paquet large dans une petite boîte **ne part jamais**. Le champ
résident à l'atome, mesuré après l'arrivée du paquet, passe de `4.8e-3` (largeur 8) à `9.7e-1`
(largeur 32) — un **bain persistant**, plus un paquet. Une excitation hors-résonance par un bain
intense est de la **physique réelle** (la queue lorentzienne du taux de Fermi golden rule et
l'élargissement de puissance), pas un défaut du code.

Le test a donc été **remplacé** : au lieu d'un seuil inventé, V61 mesure un **balayage en
désaccord** à amplitude, largeur et durée fixes, et exige que l'absorption **crête sur le gap
et décroisse** avec le désaccord. C'est ce que la légende de la figure affirme réellement.

**b) Le piège de localité du paquet.** Un paquet de largeur comparable à la boîte devient un
bain. Mesuré par V61 (boîte de 120, désaccord = la moitié du gap `1s→2p`) : fraction excitée
`0.037` à 12 cellules (0.10 de la boîte), `0.135` à 24 (0.20), `0.326` à 48 (0.40). La figure
de publication utilise un paquet localisé : 20 cellules dans une boîte de 296, soit 0.07 de la
boîte. Les nombres sont dans la légende de la figure.

Un troisième contrôle ajouté au passage, qui manquait : **l'atome seul, sans photon** — s'il
s'excitait lui-même, rien ci-dessus ne voudrait rien dire. Mesuré : fraction excitée `< 1e-9`.

---

## 6. Ce qui a été régénéré

- `simulation/hydrogen/video/` était **périmé** : les frames dataient du 5 octobre 08:22, alors
  que `config.yaml` (15:58) et `vispy_h.py` (16:05) ont été modifiés après. Les artefacts
  périmés ont été supprimés et la vidéo + les stills sont en cours de régénération avec la
  configuration courante.
- `doc/figures/fig26_hydrogen simulation.png` renommé `fig26_hydrogen.png` (un espace dans un
  lien Markdown est fragile).
- Les deux moteurs (`cl10/engine.py`, `simulation/hydrogen/backend.py`) acceptent maintenant
  `CL10_BACKEND=numpy` pour forcer le CPU. Nécessaire ici : les deux cartes GPU étaient
  saturées (23.8 Go) par un autre processus. À l'ajout du switch, Pyright a signalé `os` non
  importé dans `engine.py` — **vraie** erreur, corrigée (contrairement à ses faux positifs
  habituels). Le switch devait aussi bloquer JAX, sinon JAX demande ~17 Go.

---

## 7. Suite de validation

```
run_validation.py  ->  62/62 PASS   (exit 0)
```

V60 (le garde-fou TBM) et V61 (l'hydrogène) sont les deux nouveaux. Les expériences touchées
ont été relancées individuellement : V28, V29, V39, V45, V49, V55, V56, V60 — exit 0, JSON
valides.

---

## 8. Bilan

**Livrable.** Les quatre problèmes du review3 sont corrigés. Le cinquième — celui que la suite
rendait invisible — est corrigé et couvert par un test qui aurait échoué avant. La figure
hydrogène est intégrée avec sa provenance écrite, et son test a été **renforcé après avoir
échoué**, plutôt qu'assoupli pour passer.

**Ce qui reste, et qui est une limite déclarée plutôt qu'un problème ouvert :**

- `sin²θ23` et `sin²θ13` ne sont pas prédits 0-paramètre (dégénérescence 2-fold).
- Le CKM non trivial exige un misalignement up/down que la base symétrique ne fournit pas.
- `n = 5` n'est pas sélectionné par les données.
- `28/79` est un input du secteur SM.
- L'atome d'hydrogène est un input : seule la partie photon est CL10.

Les franchir demanderait une modification de la théorie, pas une correction.
