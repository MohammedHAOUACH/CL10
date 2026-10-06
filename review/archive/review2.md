# Review2 — avis complet sur `doc/publication.md` (CL10), état après les corrections

**Date :** 5 octobre 2026. **Objet :** `doc/publication.md` (1995 lignes, 25 figures), relu
intégralement, avec re-vérification indépendante des chiffres et re-exécution de
`run_validation.py` (60/60 PASS, `toe: true`).

**Verdict en une phrase :** le noyau algébrique est **correct, vérifié et reproductible**, et
la discipline d'honnêteté (audit des provenances, falsificateurs, « checks internes ≠ tests
contre la nature ») est la meilleure partie du papier ; mais la **substance phénoménologique
est un ordre de grandeur obtenu par deux choix (base 16π, exposant 10) plus trois coefficients
lus dans la donnée**, et il reste **quatre erreurs factuelles ou contradictions internes**
que la review1 n'a pas corrigées.

---

## 1. Ce qui est solide (et vérifié par moi)

- **L'algèbre est juste.** Cl(10) : 1024-D, module 32-D, pseudoscalar anti-hermitien
  (classe n ≡ 2 mod 8), split 16/16-bar, 16 = 10 ⊕ 5̄ ⊕ 1 de SU(5), 45 adjoint,
  argument de minimalité (Spin(8) n'a pas de 16 : ses irr. sont 8, 28, 35) — tout est
  standard et exact. Le 5/3 est le calcul usuel : Tr(Y²) = 5/6 sur le 5 de SU(5) ✓
  (recalculé : 0.8333). Les 8 coefficients d'anomalie s'annulent ✓ (recalculés).
  sin²θ_W = 3/8, m_W/m_Z = √(5/8) ✓ (prédiction SU(5) usuelle).
- **Les vrais résultats calculés** sont : l'identité SL(32) (det U = 1 pour tout g₁₀),
  le miroir spectral 16/16-bar (CPT comme conjugaison spectrale, indépendant de g₁₀),
  le théorème CP-odd (tous les invariants CP-impairs nuls en 0-paramètre), la tour
  U^(2^p) bornée. Ce sont des faits vérifiables, pas de la rhétorique.
- **La Section 15 est exemplaire.** Le passage « 0-parameter » → cinq classes de
  provenance (algebraic / free / convention / coincidence / named) est exactement ce
  qu'un reviewer demande. La phrase « les 60 checks vérifient la cohérence interne, pas
  la physique contre la nature » est la phrase la plus honnête du papier.
- **Reproductibilité :** `run_validation.py` ré-exécuté ici → 60/60, verdict identique.

## 2. Ce qui reste problématique (fond, pas forme)

### 2.1 L'admission d'honnêteté a déplacé le problème, ne l'a pas résolu
Le papier dit maintenant explicitement que 2.08, 1.18 et le ratio 1.12 sont des valeurs
**lues dans la donnée**. C'est bien. Mais alors la thèse phénoménologique se réduit à :
*choisir la base 16π, choisir l'exposant 10, obtenir des ordres de grandeur*. Un ordre de
grandeur n'est pas une prédiction — c'est une coïncidence à une décade. Le titre
« 0-parameter theory » n'est défendable que pour l'algèbre, et le papier le dit ;
l'abstract devrait le dire **d'abord**, avant les matches numériques.

### 2.2 « Une échelle, deux énoncés » est encore présenté comme une confirmation
Recalculé : `M_P(1/16π)^10 = 118.48 GeV` et `2.08 × 118.48 = 246.4 GeV = v`. La masse DM
est **par construction** v/2.08 : les « deux contraintes » sont la même équation écrite
deux fois. L'unicité de n = 10 n'existe que dans la fenêtre d'une décade (le critère
facteur-2 strict est vide, comme le dit maintenant le papier). La figure 18 est correcte
dans sa nouvelle formulation ; le texte de la Section 16 (« two constraints, one solution »)
garde un parfum de confirmation indépendante qui ne l'est pas.

### 2.3 Le WIMP miracle n'est pas un résultat de CL10
Ωh² ≈ 0.07 vs 0.12 est une propriété **générique** de tout relic thermique à ~100 GeV
(⟨σv⟩ ~ 1/m²). CL10 ne fait que *poser* m_DM = 118 GeV via le même facteur de boucle.
Le papier l'admet (« the WIMP miracle is a property of the freeze-out calculation »).
Donc le contenu DM du papier = la hiérarchie reformulée. Le candidat réel (le soliton) est
Planckien et sous-abondant de 10⁵² — le papier le dit. Le soliton est stable au sens de
l'ODE `u'' = u − u³` (stabilité classique), pas au sens quantique ; F7 est donc un
falsificateur faible.

### 2.4 Le « S-matrix » n'est pas une S-matrix
`S = U` sur un espace interne 32-D sans états asymptotiques, sans dépendance en énergie,
sans diffusion spatiale : c'est une matrice de rotation interne. La « cross-section »
σ₀ = P₂/M_P² est une **définition**, pas un calcul de processus. La borne de
Froissart-Martin est satisfaite trivialement (section constante et minuscule). Appeler cela
« the first computed process of the theory » surenchérit ; « the first computed *kinematic*
object » serait juste.

### 2.5 Le verdict TOE est auto-délivré
Les critères P0–P6 sont définis par les auteurs, les 60 checks sont pour l'essentiel des
vérifications d'identité (V8 : `√2 sech x` résout `u'' = u − u³` — identité ; V37 : l'action
est écrite à partir de la dynamique, l'Euler-Lagrange redonne la dynamique — tautologie
par construction). Le papier le reconnaît (§11 « what the 60 checks are and are not »).
Alors le verdict devrait s'appeler « cohérence interne complète », pas « TOE positive ».
La formulation « in the strongest defensible sense » est correcte mais l'abstract et la
conclusion utilisent le mot TOE sans cette protection partout.

### 2.6 Antériorité bibliographique manquante
Le 16 de Cl(10)/Spin(10) comme génération SM complète est connu depuis 1984
(Frampton, Nucl. Phys. B241 (1984) ; Grönblom & Hoebråten, Phys. Lett. B143 (1984) ;
Asherova et al.). La littérature « Clifford algebra GUT » est absente des références.
La nouveauté réelle du papier est le **QCA dynamique + l'audit de provenance**, pas le
contenu matière — un reviewer externe le signalera immédiatement.

## 3. Erreurs factuelles et contradictions internes (non corrigées par la review1)

1. **Théorème 6 — erreur arithmétique.** `α₁₀ = g₁₀²/4π` est faux : (1/5)²/4π = 0.00318,
   pas 0.04. Le code pose `ALPHA_10 = 1/25` indépendamment de g₁₀. Il faut choisir : soit
   α₁₀ = g₁₀² (sans 4π), soit admettre que 1/25 est une convention distincte. Le papier
   affirme les deux.
2. **Contradiction rank 3 / rank 1.** Section 10 : la matrice de masse 3×3 est
   « non-singular (rank 3): all three generations acquire mass ». Theorem 13 / W6 : la base
   démocratique a les valeurs propres {3,0,0}, donc **rank 1** — deux générations sans masse.
   Le code reflète la scission : `cl10/yukawa.py` utilise l'identité (rank 3),
   `cl10/yukawa_3x3.py` le tout-uns (rank 1). C'est la même contradiction dans deux modules.
   Il faut trancher : la structure 0-paramètre est rank 1 (démocratique), et « toutes les
   générations acquièrent une masse » doit être retiré ou requalifié en residual.
3. **Mélange de conventions d'hypercharge.** Proposition 1 annonce « Y in Q = T₃ + Y/2 »
   mais donne les valeurs de la convention Q = T₃ + Y (1/6, −2/3, …). Fig. 9 et le code
   (`cl10/matter.py`, `SM_FIELDS`) sont dans Q = T₃ + Y et annulent ✓. Theorem 3 est dans
   Q = T₃ + Y/2 ✓. Une seule convention doit être annoncée et tenue partout.
4. **Pôle de Landau U(1) : trois valeurs selon l'endroit.** Docstring `cl10/unification.py` :
   ~3×10³⁹ GeV ; code exécuté : 1.39×10³⁹ ; Fig. 2 : ~10³⁹. De plus la valeur dépend de la
   normalisation de Y (avec la normalisation GUT, le pôle 1-loop est ~10²⁵). Le point
   physique (le pôle est au-dessus de M_P en 1-loop, en dessous en 2-loop) est robuste ;
   les chiffres affichés ne le sont pas.
5. **Terminologie « n_gen » glissante.** `10 = n_gen` (générateurs de Cl(10)) partout, mais
   Fig. 18 dit « the naive n_gen = 5 satisfies neither » — 5 est le nombre de qubits, pas
   n_gen. Deux sens du même symbole dans le même papier.
6. **Fig. 23** : front à c = 0.5 (limite CFL du solveur 2-D) — le papier le signale lui-même,
   bien. Mais la figure illustre une équation 2-D, pas le QCA 3-D causal à rayon 1 ; à
   recadrer comme « illustration du caractère causal », pas comme le QCA.
7. **B–L mal défini.** « the 16 has B–L = −1 » n'est pas un nombre : B–L est un générateur de
   U(1) dont les valeurs varient sur les composantes du 16 (Q 1/3, u^c −1/3, L −1, e^c −1,
   N −1 ; somme sur les 16 composantes = −3). Ce qui est vrai est que le générateur **change
   de signe** entre le 16 et le 16-bar — donc « opposite B–L » est correct au niveau du module,
   mais la valeur unique −1 est une convention. Et « the 1 of the 16 is the sterile ν^c
   (B–L = 0) » est **faux** : ν^c a B–L = −1 (dans la décomposition SO(10) → SU(5)×U(1)_{B-L},
   le 1 porte −1). Le mécanisme d'annihilation de la Section 9 n'est pas affecté, mais ces deux
   phrases sont des erreurs factuelles.

## 4. Ce que je garderais si j'étais reviewer

- Sections 1–4 et 8 (algèbre, (1,3), 16, 45, 5/3, anomalies) : correctes, publiables en l'état.
- Section 15 (audit) + le théorème CP-odd + le miroir spectral + SL(32) : les contributions
  réelles du papier.
- Section 13 (Table 2) : la formulation « mechanism solved / exact value named » est la bonne
  discipline.
- À réécrire ou affaiblir : l'abstract (trop dense, le compte de paramètres y change trois
  fois), « the first computed process », « two independent constraints », le mot TOE dans la
  conclusion, la Section 5 (mesure/collapse : elle n'ajoute rien au-delà du postulat de
  projection, et le papier l'admet — elle pourrait être un paragraphe).

## 5. Note d'ensemble

C'est un papier **honnête pour un objet qui ne l'est presque jamais**. Sa vraie thèse est :
« une algèbre 0-paramètre produit la structure du SM, et tout ce qui est numérique est soit
un ordre de grandeur de boucle, soit un coefficient nommé lu dans la donnée ». Cette thèse est
défendable et vérifiée. Ce qui ne l'est pas : le mot « Theory of Everything », qui suppose que
les ordres de grandeur comptent comme des prédictions. Si le papier se renomme
« A 0-parameter Clifford-QCA framework for the SM structure, with an explicit parameter audit »,
il devient difficile à attaquer. Avec le mot TOE, il reste attaquable sur un point simple :
**trois nombres lus dans la donnée ne sont pas dérivés de la racine, et le papier le sait.**
