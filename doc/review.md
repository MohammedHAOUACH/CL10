# Revue critique de `doc/publication.md` (CL10)

**Date :** 5 octobre 2026. **Objet :** `doc/publication.md` (1859 lignes), avec vérification
des affirmations dans `cl10/` et `experiments/`.

**Méthode :** la suite `run_validation.py` a été exécutée intégralement (60 expériences,
~30 s, backend numpy) → 60/60 PASS, verdict `toe: true`. Les chiffres du texte ont été
recalculés indépendamment. Les références croisées (Fig. 1-25, Sections 1-16) ont été
vérifiées.

**Verdict en une phrase :** le noyau algébrique est correct et proprement vérifié ;
la couche phénoménologique (hiérarchie, matière noire, asymétrie baryonique, hiérarchie
des familles) repose sur **sept arguments circulaires** que le texte nomme « residuals
nommés » mais qui sont, dans le code, des coefficients calculés à partir des valeurs
observées.

---

## 1. Les arguments circulaires (classés par gravité)

### C1 — Le coefficient « 2.08 » est calculé depuis la valeur observée (le plus grave)

**Le texte (l. 405-411, 1056, 1216, 1485) :**
> « the pure 0-parameter prediction (1/16π)^10 = 9.7e-18 is a factor 2.08 below the
> observed 2.016e-17; the 2-digit coefficient 2.08 that closes the gap is a **named**
> residual, **not** a 0-parameter prediction »

**Le code (`cl10/four_limits.py:143`) :**
```python
coeff_hier = float(hier_obs / LOOP ** best["n"])     # ~2.08
coeff_dm   = float(best["dm_GeV"] / dm_obs)          # ~1.18
```
avec `hier_obs = V_EW / M_P` et `V_EW = 246.0` (`cl10/toepred.py:30`).

**Le cercle :** la hiérarchie est prédite comme `v/E_P = 2.08 (1/16π)^10`. Le nombre
`2.08` est obtenu en divisant la hiérarchie mesurée par le facteur de boucle. On insère
donc la valeur observée dans la formule, puis on déclare que le coefficient inséré est
un « residual nommé » et non un fit. La distinction entre « nommé » et « ajusté » est
rhétorique : le coefficient est *déduit de la donnée*. Le même schéma vaut pour 1.18
(matière noire) et pour le ratio 1.12 (asymétrie baryonique).

**Pourquoi c'est un problème pour la thèse centrale :** le titre annonce « 0-parameter ».
L'audit (Section 15) est louable — il classe honnêtement les nombres — mais il se
rattrape lui-même : la Section 16 déclare « 0 free parameters + 3 named », alors que ces
3 named sont exactement les 3 coefficients qui ont été calculés depuis les mesures. Le
compte final « 3 » n'est pas une réduction du nombre de paramètres ; c'est le nombre de
paramètres cachés dans le label « named ».

**Comment le rendre non-circulaire :** il faut *dériver* 2.08, 1.18 et le phase CP
à partir du spectre de `U(g_10)` (les 32 niveaux). Tant que ces trois nombres ne sortent
pas d'un calcul sur la racine, la théorie est 3-paramètres, et le dire explicitement
dans l'abstract serait plus fort que le défendre.

---

### C2 — « One factor, two results » : les deux contraintes ne sont pas indépendantes

**Le texte (l. 1451-1460, 1727-1732, Fig. 18) :**
> « the same loop count n = 10 is the unique n in 1..20 that simultaneously (i) puts the
> hierarchy at 2.08 (1/16π)^10 and (ii) puts the DM at the WIMP mass M_P(1/16π)^10 = 118
> GeV. **One 0-parameter number, two independent constraints, one solution** »

**Vérification numérique :**
- `M_P (1/16π)^10 = 118.48 GeV`
- `2.08 × 118.48 = 246.4 GeV = v` (le vev électrofaible, hard-codé à 246.0)

Les « deux contraintes indépendantes » sont la **même équation** : la masse de la matière
noire est, par construction, `v / 2.08`. Il n'y a pas deux contraintes satisfaites par un
même `n` ; il y a une contrainte (`(1/16π)^n ≈ v/M_P`) énoncée deux fois, une fois comme
ratio de masses, une fois comme masse absolue.

**Le cercle :** on présente la coïncidence hiérarchie ↔ WIMP comme une confirmation
indépendante de `n = 10`. Cette confirmation est analytiquement triviale : si la
hiérarchie est correcte, la masse WIMP l'est automatiquement, puisque `m_DM ≡ v/2.08`.

**Note technique :** le docstring de `unique_10loop` annonce « within a **factor 2** of
both », mais le code teste `hier_dev < 1.0` et `dm_dev < 1.0` en **décades** (facteur 10).
Avec un critère facteur 2 strict, l'ensemble est **vide** — aucune valeur de `n` ne
satisfait les deux. L'unicité de `n = 10` n'est sauvée que par la fenêtre large et par
l'ancrage `dm_obs = 100.0` GeV (valeur de référence arrondie, pas une mesure).

---

### C3 — La hiérarchie est expliquée par un choix de base et d'exposant qui la reproduit

**Le texte (l. 402-403, 1340-1350) :** l'exposant `-17` est « dérivé » comme
`10 × log10(16π) = 17.01`, où `10 = n_gen` (le nombre de générateurs de Cl(10)) et
`16π` est le facteur de boucle 1-loop.

**Le cercle :** `16π ≈ 50.27` et `50^10 ≈ 9.8e17`. N'importe quelle base de l'ordre de 50
à la puissance 10 donne 10^-17. Le calcul alternatif `(1/4π)^n = 2.016e-17` donne
`n = 15.19` — tout aussi « 0-parameter » (4π est aussi un facteur de QFT). Le choix de
la base `16π` **et** de l'exposant `10` est sélectionné *parce qu'il tombe juste*.

De plus, `n_gen = 10` n'est pas un nombre libre : Cl(10) a été choisi **parce que** son
module 32-D contient une génération SM complète. Le raisonnement est donc :
1. choisir Cl(10) pour reproduire la matière du SM ;
2. prendre le 10 de Cl(10) ;
3. obtenir la hiérarchie.

C'est une sélection post-hoc de la racine, présentée comme une conséquence de la racine.
Le texte est d'ailleurs honnête à l'intérieur (l. 1236-1248 : « the exponent -17 is **not**
a 0-parameter consequence of the one 32-level step ») mais l'abstract (l. 84-87) affirme
le contraire (« the hierarchy exponent -17 = the 10-loop factor »). L'abstract et le
corps se contredisent.

---

### C4 — Les 3 familles : le même nombre compté trois fois

**Le texte (l. 210-217) :**
> « The 3 observed families are a 0-parameter consequence of the one object in **three**
> equivalent ways: the 3 spatial directions of the (1,3); the 3 two-planes C(3,2) of the
> 3-space; and the rank 3 of the spatial metric. »

**Le cercle :** les trois « voies » sont trois reformulations combinatoires d'un seul
nombre (3), déjà introduit par le choix de la signature (1,3). Ce n'est pas trois
dérivations convergentes, c'est une seule hypothèse énoncée trois fois. Le texte le
reconnaît lui-même en nommant W3 (« the *why-3* — a derivation of why the family count
equals the spatial dimension ») : la dérivation est exactement ce qui manque.

**Point connexe (l. 149-152) :** « The choice of *which* four of the ten are the spacetime
is 0-parameter: any four mutually-anticommuting generators define a (1,3) sub-algebra,
and the GUT Cl(10) is the **minimal** unified root. »
- « any four define a (1,3) » prouve l'*existence* d'une (1,3), pas que *la* spacetime est
  celle-ci. La sélection est libre.
- « minimal » est faux sans preuve : Cl(8) (complexe) a un module irréductible 16-D et
  contient Cl(1,3). Cl(10) est minimal pour *Spin(10)*, pas pour la structure spacetime.

---

### C5 — Le programme de validation se valide lui-même

**Le texte (l. 631, 1143) :** « **all 60 experiments pass** … the **TOE verdict is positive** ».

**Le code (`run_validation.py`) :** le verdict est calculé comme
```
toe_criterion = "all of tiers 0-3 (V0-V7, V11-V13) PASS"
```
Les critères P0-P6 (Table 1) sont définis par les auteurs, et le verdict TOE est un
booléen dérivé de ces critères. Le critère est choisi pour coïncider avec les propriétés
que la construction a par définition.

**Nature des expériences :** la majorité sont des **identités mathématiques de la
construction**, pas des tests de la physique :
- V0 : `n_gen = 10`, `n_basis = 1024`, `n_weyl = 16` — théorèmes de l'algèbre de Clifford.
- V11 : annulation des anomalies du 16 de SU(5) — vrai pour **tout** GUT SU(5), pas
  spécifique à CL10.
- V6 : une équation d'onde discrète avec `c = 1` a `c = 1`.
- V8 : `u = √2 sech x` résout `u'' = u - u³` — identité analytique.
- V25 : `S = A/4` — c'est la formule de Bekenstein-Hawking, importée, pas dérivée du QCA.
- V37 : l'action discrète construite à partir de l'opérateur de la dynamique a pour
  Euler-Lagrange cette dynamique. Tautologie : l'action est écrite *à partir* de
  l'équation. C'est une identité correcte, pas une dérivation de la relativité générale.

**Ce qui est réellement testé contre la nature :** 5 comparaisons (durée de vie du
proton, `m_ν`, `ΔB`, `m_DM`, `v/E_P`) — et c'est exactement là que C1-C3 opèrent.

**Note de reproductibilité :** `results/validation_summary.json` sur disque est **stale**
(`n: 0, n_pass: 0, toe: false, experiments: []`) alors que le texte affirme 60/60. La suite
passe réellement (vérifié en exécution : 60/60 PASS, verdict TOE), mais l'artefact livré
contredit la revue. `results/` est dans `.gitignore` : un lecteur ne peut pas vérifier sans
exécuter le code.

---

### C6 — La matière noire : le mécanisme est emprunté, l'échelle est définie par la cible

**Le texte (l. 459, 1055, 1352-1360) :** `m_DM = M_P (1/16π)^10 = 118 GeV` → « the WIMP
miracle … gives the observed abundance (a 0-parameter cosmological result) ».

**Le cercle :** le « miracle WIMP » est un calcul standard du freeze-out du SM avec une
section efficace de portée faible. CL10 n'a pas de section efficace de portée faible : sa
seule échelle est `M_P`, et son candidat 0-parameter (le soliton) est de masse Planck —
le texte le reconnaît (gap `10^-52`, l. 846-849). Pour obtenir l'abondance observée, on
importe la formule WIMP standard et on applique la masse `118 GeV` déjà choisie pour
coïncider avec l'échelle électrofaible (C2). Le calcul de freeze-out n'est pas issu du
QCA.

Le texte est d'ailleurs contradictoire sur l'identité du candidat : Section 9 dit que le
DM **est** le soliton (l. 443-459), Section 12 (C7) dit que le soliton de masse Planck
**n'est pas** le DM observé (l. 846-849), Section 15 dit que le DM **est** le relic 10-loop
(l. 1300-1304). Trois réponses différentes à la même question.

---

### C7 — L'asymétrie baryonique : le théorème nie le mécanisme, la formule l'affirme

**Le texte (l. 1251-1261, 1462-1475) :** un théorème prouvé — tous les invariants CP-impairs
du QCA 0-parameter s'annulent, donc la dynamique satisfait **aucune** des trois conditions
de Sakharov et **ne peut pas** générer l'asymétrie baryonique. C'est le meilleur résultat
du papier.

**Puis (l. 801-811) :** `ΔB = (28/79) · (g_0^4/16π) · (1/16π)^2 = 6.8e-10`, ratio 1.12 avec
l'observé, « every number is 0-parameter ».

**Le cercle :** la formule est la leptogenèse standard, dont le mécanisme vient du SM, pas
du QCA. Ses ingrédients :
- `g_0 = 1/8 = 1/2^3` : l'audit lui-même classe `g_0` comme **convention** (l. 1205), pas
  comme une dérivée. Sa sensibilité est brutale : `g_0 = 1/4` → ratio 17.9 ; `g_0 = 1/16`
  → ratio 0.07. Le « match 12 % » dépend entièrement d'un choix non dérivé.
- `g_0^4/16π` comme asymétrie CP-impair 1-loop : ce n'est pas le résultat standard de la
  leptogenèse (l'asymétrie vient des diagrammes vertex/自 énergie avec les masses des `N`).
- `D = (1/16π)^2` « 2-loop dilution » : aucune dérivation, le choix du facteur est justifié
  par « loop counts 1 + 2 = 3 ».

Donc : le théorème dit « le QCA ne peut pas le faire », et le chiffre est produit par une
formule extérieure dont les paramètres sont choisis pour tomber juste. Ce n'est pas une
contradiction formelle (le texte nomme la phase CP comme residual), mais c'est le
meilleur exemple du pattern : *un mécanisme absent + une formule empruntée + un coefficient
nommé = une prédiction 0-parameter annoncée*.

---

## 2. Autres faiblesses (non circulaires mais importantes)

**I1 — Contradiction interne sur la constante cosmologique.** Section 9 : la
cancellation 512/512 donne `Λ_QCA = 0` (l. 467-468). Theorem 10 (l. 875-881) : le vide de
de Sitter donne `G_{μν} = -3H² g_{μν}`, donc `Λ ~ M_P^4` au Hubble de Planck. Le texte
ne réconcilie jamais les deux. Le « s-mechanism » (« the flat vacuum is a saddle point of
the energy, so the 0-point does not gravitate ») est affirmé, jamais démontré à partir de
l'action.

**I2 — Le graviton spin-2 est postulé, pas dérivé.** `cl10/einstein.py:84,113` compte
« 2 polarizations » à partir d'une équation d'onde **scalaire** sur le réseau. Une onde
scalaire discrète a une amplitude, pas deux hélicités. Le champ `h^{μν}` (perturbation
d'une métrique 4×4) est introduit à la main (l. 329-332). Rien dans Cl(10) n'impose que
ce champ se transforme comme un tenseur de masse nulle de spin 2. C'est aussi une
contradiction avec P0 (« no added fields ») : la métrique est un champ ajouté.

**I3 — Nielsen-Ninomiya (l. 851-858).** « The doubling theorem does not apply … the chiral
content is a property of the algebra, not of a lattice operator. » C'est une exemption par
catégorie, pas une démonstration. La chiralité est *assumée* dans le module interne 32-D ;
le théorème est écarté parce que le module n'est pas un opérateur de Dirac. Il faudrait
montrer que le spectre basse énergie du QCA local (V36) reproduit réellement 16 Weyl
chiraux sans doublement — ce n'est pas fait.

**I4 — La réduction « 9-to-1 » (l. 972-1010) est un renommage, pas une dérivation.**
« Les neuf residuals sont neuf projections d'un seul objet » signifie : « l'objet ne fixe
pas neuf choses ». Aucun mécanisme ne relie `Λ = 10^-122` à la matrice CKM via les 32
niveaux de `U(g_10)`. L'affirmation est rhétorique.

**I5 — `g_10 = 1/5` : le texte se contredit sur son propre statut.**
- Point 1 (l. 1224-1231) : `g_10` **n'est pas** fixé, « the theory is **1-parameter** in
  its dynamics ».
- L1 (l. 1439-1449) : `g_10 = 1/5` est un « root number » (`1/n_qubits`), donc « the theory
  has **no free parameters** ».
Le pont entre les deux est exactement la lecture `5 → 1/5` que l'audit qualifie lui-même
de « mapping choice, not a derivation » (l. 78, 1230). L'argument est : l'angle est libre ;
il n'est pas libre parce qu'on le lit comme `1/5` ; la lecture est une convention. Le
paramètre disparaît par changement d'étiquette.

**I6 — Falsifiabilité largement rhétorique (l. 608-626).**
- F3 : la hiérarchie est définie comme le ratio de deux nombres mesurés ; elle ne peut pas
  « entrer en conflit » avec la mesure.
- F5 : `τ_p ~ 4e47` yr, soit 13 ordres au-dessus de la borne. Aucune expérience
  foreseeable ne peut la toucher. De plus, `M_GUT = M_P` est déclaré **convention** (l. 1205) :
  avec `M_GUT = 10^16` GeV, `τ_p ~ 1.8e35` yr — la « prédiction » est un artefact de la
  convention.
- F6 : la borne est `ΔB ≪ 10^-11` ou `≫ 10^-8` ; la phase CP nommée peut absorber n'importe
  quelle valeur dans cet intervalle.
Aucun des falsificateurs n'est opérable.

**I7 — Références à vérifier.**
- Choptuik, PRL **70**, 1 (1993) → la pagination correcte est **70**, 9-12.
- « M. Bott, Periodicity in the theory of Clifford algebras, Comment. Math. Helv. **45**,
  197 (1970) » — à confirmer ; la périodicité de Bott est habituellement citée via
  *Adv. Math.* 1964 / Atiyah-Bott-Singer.
- Kay & Laflamme, CQG **24**, 5469 (2007) : le titre exact est « Causal set dynamics and
  quantum-gravitational field theory » (2007) — à vérifier.
- Fig. 7 renvoie à « Section 9.3 » qui n'existe pas (les sections ne vont pas en
  sous-sections).

**I8 — Erreurs de rédaction.**
- l. 831-832 : phrase coupée — « the 45 **mix** the / This is the first computed process… ».
- l. 1287-1288 : « This **refines** the CPT theorem of the CPT theorem ».
- l. 1565 : « the **Planck-mass** Planck-mass soliton ».
- l. 673 : « the the heavy-$N$ out-of-equilibrium ».
- l. 1728-1731 : Fig. 18 appelle la hiérarchie « the **1-loop** factor » et la DM « the
  **2-loop** factor » alors que partout ailleurs les deux sont le **10-loop**.
- l. 1793 : Fig. 23 annonce `c = 0.5` alors que la causalité du QCA est annoncée `c = 1`
  partout ailleurs (l. 112, 145, 337). Le dt = h/2 du CFL explique le 0.5 en unités
  réseau, mais le texte doit le dire.

---

## 3. Ce qui tient réellement (et mérite d'être mis en avant)

1. **L'algèbre est correcte et doublement vérifiée.** La construction 5-qubits
   symplectique-Pauli donne exactement 10 générateurs anticommutants, 1024-D, module 32-D,
   pseudoscalaire anti-Hermitien avec `w² = -I`, split 16/16-bar. Le cross-check symbolique
   (`kingdon`) contre le numérique est une vraie rigueur, rare dans ce type de papier.
2. **Le miroir spectral 16/16-bar (l. 1273-1292)** est un résultat propre : `[U, w] = 0`
   et `spec(U|_16bar) = conj(spec(U|_16))`, vérifié pour **toutes** les valeurs de `g_10`.
   C'est une propriété de l'algèbre, indépendante du paramètre. C'est le meilleur
   « 0-parameter » du papier.
3. **Le théorème CP-impair (l. 1251-1261)** est le résultat le plus intéressant : il prouve
   que le QCA 0-parameter ne peut pas produire l'asymétrie baryonique. C'est une limite
   réelle, démontrée numériquement sur un balayage complet (10×10×45×45).
4. **La vérification 2-route de la courbure 4-D** (l. 924-950) est solide : Riemann calculé
   par différences finies sur une métrique arbitraire, concordant avec l'analytique de Sitter
   à `O(h²)`, Schwarzschild retrouvé (`Ricci = 0`, `K = 48M²/r⁶`). C'est de la NR correcte.
5. **L'honnêteté de l'audit (Section 15)** est la meilleure partie du papier : classer les
   nombres en cinq statuts épistémiques est exactement ce que la plupart des « TOE »
   populaires refusent de faire. Le problème est que la Section 16 annule cette honnêteté
   en requalifiant les conventions en « 0-parameter ».

---

## 4. Recommandations concrètes

**À faire pour rendre la thèse défendable :**

1. **Réécrire l'abstract.** Annoncer : « 0 free parameters in the algebra; **3 named
   constants** (2.08, 1.18, δ_CP) whose values are not derived ». Le « 0-parameter » du
   titre est le point faible du papier ; l'honnêteté est un argument de vente plus fort.
2. **Supprimer C2.** « One factor, two results » doit devenir « one factor, one result,
   stated as a ratio and as a mass ». La prétendue indépendance des deux contraintes est
   fausse et attaquable en une ligne.
3. **Dériver ou supprimer les coefficients.** Soit 2.08 sort du spectre de `U(g_10)`, soit
   il faut l'appeler explicitement un *fit*. Le mot « named residual » ne change pas la
   nature épistémique d'un nombre calculé depuis la donnée.
4. **Réconcilier Λ = 0 et Λ ~ M_P⁴** (I1) ou supprimer l'une des deux affirmations.
5. **Dériver le spin-2** (I2) : montrer que le champ `h^{μν}` du réseau porte la
   représentation (2,0)⊕(0,2) du little group, sinon dire que la métrique est un champ
   ajouté (et abandonner P0 pour la gravité).
6. **Traiter Nielsen-Ninomiya sérieusement** (I3) : exhibiter le spectre basse énergie du
   QCA local (V36) et compter les modes chiraux.
7. **Rendre F1-F7 opérables** (I6) : au moins une falsification à portée expérimentale
   actuelle. Sinon, dire « falsifiable en principe » plutôt que « falsifiable ».
8. **Corriger les références et les coquilles** (I7, I8) — un reviewer sérieux commencera
   par là, et une référence fausse décrédibilise le reste.
9. **Commiter `results/validation_summary.json`** (I5) : la revue par les pairs exige que
   le verdict affiché soit vérifiable sans exécuter le code.

---

## 5. Avis global

Le document est **deux papiers distincts collés** :

- **Un bon papier de physique mathématique** (Sections 1-4, 12, 15) : une construction
  Cl(10) rigoureuse, doublement vérifiée, avec des résultats réels (16/16-bar, miroir
  spectral CPT, théorème CP-impair, NR 4-D 2-route). Cette partie est publiable telle
  quelle après nettoyage.

- **Un papier de phénoménologie circulaire** (Sections 7, 9, 13, 16) : la hiérarchie, la
  matière noire, l'asymétrie baryonique et les familles reposent sur des nombres choisis
  pour reproduire l'observation, puis rebaptisés « residuals nommés ». Le mécanisme
  d'échappatoire est identique dans les quatre cas : *un facteur de boucle dont la base et
  l'exposant sont sélectionnés par la cible, un coefficient calculé depuis la mesure, et
  l'affirmation que ce coefficient n'est pas un fit.*

La force du projet est la discipline de vérification (60 expériences, tolérances, 2 routes).
Sa faiblesse est que cette discipline s'applique à des identités internes, pas aux
prédictions qui comptent. Le « 60/60 PASS » ne valide pas la physique ; il valide que la
construction est cohérente avec elle-même.

**Verdict :** le noyau est solide et intéressant ; la thèse « Theory of Everything
0-parameter » n'est pas démontrée, et son argumentaire central est circulaire. Le papier
gagnerait énormément en crédibilité à se présenter comme « une construction algébrique
0-parameter de la matière, du gauge et de la gravité, avec trois constantes irréductibles
non dérivées » — c'est exactement ce que l'audit (Section 15) dit déjà, avant que la
Section 16 ne le nie.
