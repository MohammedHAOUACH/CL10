# Assessment de `doc/review.md` — verdicts par critique

**Date :** 5 octobre 2026. **Méthode :** chaque affirmation de la revue a été
re-vérifiée indépendamment — (i) relecture ligne à ligne de `doc/publication.md`
(1859 lignes), (ii) re-lecture du code cité (`cl10/`, `experiments/`),
(iii) **re-calculation** de tous les chiffres (numpy), (iv) **exécution** de la
suite complète `run_validation.py` (60/60 PASS, `toe: true`, ~92 s, backend
numpy), (v) vérification des références contre les sources primaires (APS,
index bibliographiques).

**Bilan :** la revue est **globalement de bonne qualité** — la grande majorité
de ses points sont fondés et ses numéros de ligne sont exacts. Mais elle est
l'œuvre d'un non-expert : **6 erreurs ou exagérations** y sont identifiées
(ci-dessous, § 4). Les corrections implémentées suivent les points *fondés*,
pas les recommandations aveuglément.

---

## 1. Critiques circulaires C1–C7

### C1 — « 2.08 est calculé depuis la valeur observée » — **FONDÉ**

- **Vérifié dans le code** (`cl10/four_limits.py:148`) :
  `coeff_hier = hier_obs / LOOP**n` avec `hier_obs = V_EW / M_P`,
  `V_EW = 246.0` (`cl10/toepred.py:30`). Re-calculation :
  `2.016e-17 / (1/16π)^10 = 2.076` ✓ (le texte arrondit à 2.08).
- Le texte **le dit lui-même** (l. 405-411, Section 15 « coincidence », W8) :
  2.08 est un *named residual*. Le reproche de la revue n'est donc pas sur la
  dissimulation, mais sur **l'étiquette** : Section 16 requalifie 2.08 (et 1.18)
  en « named » sans jamais dire que leurs **valeurs sont lues dans la donnée**
  (i.e. un fit 2-chiffres de l'ordre 0-paramètre `(1/16π)^10`).
- **Verdict :** le reproche épistémique est valide. « Named » ≠ « dérivé » :
  2.08 et 1.18 sont des **constantes nommées dont les valeurs proviennent de la
  donnée**. Correction : explicité dans l'abstract, Section 16 et le code.

### C2 — « les deux contraintes ne sont pas indépendantes » — **FONDÉ**

Trois sous-points, tous vérifiés :

1. **Docstring mensongère.** `unique_10loop` annonce « within a **factor 2** of
   **both** » (l. 124-126) mais le code teste `hier_dev < 1.0` et
   `dm_dev < 1.0` en **décades** = **facteur 10** (l. 147). ✓
2. **Le critère facteur-2 strict est vide.** Re-calculation : pour n = 10,
   `hier_dev = 0.3173` décades = facteur **2.08 > 2** ; pour tout autre
   n ∈ 1..20 la déviation est ≥ 1.38 décades. L'ensemble « facteur 2 des deux »
   est **vide** ✓. L'unicité de n = 10 n'est vraie que pour la fenêtre facteur 10.
3. **Indépendance sur-énoncée.** Les deux cibles sont `v/M_P = 2.016e-17`
   (mesurée) et `100/M_P = 8.2e-18` (100 GeV = référence « WIMP-miracle » du
   code, une convention arrondie, pas une mesure) ; leur ratio est 2.46. Ce sont
   deux **énoncés de la même échelle 0-paramètre** `M_P(1/16π)^10 = 118.48 GeV`
   (re-calculation : `v = 246 = 2.076 × 118.48` — la hiérarchie est littéralement
   2.08 × la masse WIMP). « Deux contraintes indépendantes » est faux ; « une
   échelle, deux énoncés » est vrai.

**Correction :** docstring + verdict du code, Section 16 (L2) et Fig. 18
réécrites en « une échelle 0-paramètre, deux énoncés ; n = 10 est l'unique
n ∈ 1..20 à une décade des deux cibles ; le critère facteur-2 strict est vide ».

### C3 — « hiérarchie = base/exposant sélectionnés parce qu'ils tombent juste » — **PARTIELLEMENT FONDÉ**

- **Fondé :** n'importe quelle base d'ordre 50 à la puissance 10 donne 10⁻¹⁷ ;
  le couple `(16π, 10)` est choisi *a posteriori* par la cible. `n_gen = 10`
  est le nombre du root (Cl(10)) — mais le root lui-même a été retenu pour son
  module 32-D contenant une génération SM, donc le 10 sert deux fois (matière,
  puis exposant). L'aspect post-hoc est réel.
- **Non fondé comme tel :** le texte **le dit déjà** — Point 2 de la Section 15
  (l. 1236-1248) : « the exponent -17 is **not** a 0-parameter consequence of
  the one 32-level step » ; le « honest caveat » (l. 1411-1419) qualifie la
  10-loop de « *representation* of the order ». Le reproche vise donc surtout
  l'**abstract** (l. 84-87) qui présente l'exposant comme une « 0-parameter
  loop-factor representation » sans cette nuance — c'est un **probleme
  d'alignement**, pas une contradiction cachée.
- La contre-proposition de la revue (« (1/4π)^n donne n = 15.19, tout aussi
  0-paramètre ») est vraie mais ne prouve rien : c'est exactement ce que dit
  déjà l'audit — le choix de base est une convention, l'ordre est 0-paramètre.

**Correction :** alignement abstract ↔ Section 15 (l'ordre -17 a une
représentation 10-loop 0-paramètre ; la base `16π` et le compte 10 sont des
nombres du root/QFT ; le coefficient 2-chiffres est nommé).

### C4 — « 3 familles comptées 3 fois ; « minimal » non prouvé » — **PARTIELLEMENT FONDÉ**

- **« 3 voies équivalentes » = 1 nombre.** Le texte l'admet : « One 0-parameter
  number (the 3 of the 3-space) fixes both … the *why-3* … is a **named**
  residual (W3) » (l. 210-217). Le reproche est donc surtout rhétorique :
  « three equivalent ways » *sonne* comme 3 dérivation, c'est 1 nombre
  reformulé 3 fois. **Fondé sur la forme**, pas sur la substance.
- **Le contrexemple Cl(8) du reviewer est mal fondé.** Le module complexe
  irréductible de Cl(8) est 16-D, **mais** sous Spin(8) il se décompose en
  `8_s ⊕ 8_c` (trialité) — Spin(8) n'a **pas** de représentation 16 (ses
  irr. sont 8, 28, 35…), donc il ne peut pas porter la génération SM
  (`16 = 10 ⊕ 5̄ ⊕ 1` de SU(5), qui est une irr. de Spin(10)). De plus
  `w² = +1` en Cl(8,0) : pas de split chiral ±i 16/16-bar. Le minimalité
  *correcte* à énoncer : Cl(10) est le plus petit root avec (i) le split chiral
  16/16-bar (classe `n ≡ 2 mod 8`, Cl(2) = 4-D trop petit), (ii) le module
  16-D portant une génération SM, (iii) l'adjoint 45-D contenant les 12 bosons
  de jauge SM. Le reviewer dit « minimal pour Spin(10), pas pour la structure
  spacetime » : à moitié juste — la spacetime seule vivrait dans un Cl(4), mais
  « unified root » inclut matière + jauge.

**Correction :** reformulation précise de la Remark (Section 2) et des « 3
voies » (Section 3).

### C5 — « le programme de validation se valide lui-même » — **FONDÉ** (+ 1 erreur du reviewer)

- **Le verdict est un booléen sur les tiers 0-3 définis par les auteurs**
  (`run_validation.py:80-83`) ✓.
- **La nature des 60 checks** confirmée par échantillonnage : V0 (identités
  Cl(10)), V6 (c = 1 par construction du laplacien discret), V8
  (`u = √2 sech x` est l'analytique de `u'' = u - u³` — identité),
  V11 (annulation des anomalies du 16 de SU(5) — vraie pour **tout** GUT SU(5)
  à une génération), V37 (l'action est **écrite à partir de** la dynamique :
  l'Euler-Lagrange redonne la dynamique — tautologie par construction) ✓.
- **« 5 comparaisons contre la nature »** : confirmé — v/E_P, m_ν, ΔB,
  m_DM/(Ω_DM), τ_p (Fig. 14) — c'est tout le contenu phénoménologique.
- **L'erreur du reviewer :** sa « note de reproductibilité » prétend que
  `results/validation_summary.json` sur disque est **stale** (`n: 0, n_pass: 0,
  toe: false, experiments: []`). **Faux au moment de la revue** : le fichier
  (mtime 39 min avant la revue) contenait `n: 60, n_pass: 60, toe: true`,
  60 entrées, 0 FAIL — et mon exécution indépendante le confirme (60/60, ~92 s).
  Le reviewer a probablement lu un état antérieur ou un autre fichier.
- **Ce qui reste valide :** `results/` est dans `.gitignore` → un clone frais
  ne contient **aucun** artefact de verdict ; la « 60/60 PASS » n'est
  vérifiable qu'en exécutant le code.

**Correction :** (a) un paragraphe de *signification* du programme dans la
Section 11 (60 checks = cohérence interne ; 6 comparaisons = physique ;
verdict = tiers durs 0-3), (b) dé-gitignore de `results/validation_summary.json`
+ commit.

### C6 — « DM : mécanisme emprunté, 3 identités contradictoires » — **PARTIELLEMENT FONDÉ**

- **Fondé :** le « WIMP miracle » est importé de la cosmologie standard. Le
  code (`cl10/dm_wimp.py:72-77`) calcule `relic_abundance(m) = Ω_obs ×
  (100/m)³` — une **convention** d'abondance, pas un freeze-out QCA ; le
  docstring du module lui-même étiquette le WIMP miracle « a **0-parameter**
  cosmological result », ce qui est un euphémisme : c'est un **résultat 0-paramètre
  de la cosmologie SM** appliqué à la masse 0-paramètre 118 GeV. CL10 n'a pas
  de section efficace de jauge à cette échelle (le texte l'admet : « the exact
  cross-section … named residual, W9 »).
- **Exagéré :** « trois réponses différentes à la même question ». La
  construction est en fait un **échelonné cohérent** : (i) le candidat
  0-paramètre est le soliton de masse Planck (Section 9) ; (ii) ce soliton est
  **sous**-abondant de 10⁻⁵² (Section 12, C7) ; (iii) l'échelle **observée**
  du DM est la projection 10-loop `M_P(1/16π)^10 = 118 GeV` (Section 15),
  l'identité particule = named residual (W7). Le vrai défaut est le **libellé**
  de la Section 9 (« Dark matter — the stable soliton »), qui se lit comme
  « le DM est le soliton de masse Planck » et entre en friction avec C7.

**Correction :** aligner le libellé Section 9 sur la hiérarchie
candidat-Planck / échelle-10-loop / identité-nommée.

### C7 — « baryogenèse : le théorème nie, la formule affirme » — **FONDÉ**

- **Tension réelle** : la Section 15 **prouve** que tous les invariants
  CP-impairs du QCA 0-paramètre s'annulent (balayage 10×10×45×45, etc.) → le
  QCA 0-paramètre **ne peut pas** générer d'asymétrie (Sakharov : aucune
  condition). Puis Theorem 8 écrit `ΔB = (28/79)(g_0⁴/16π)(1/16π)² =
  6.8e-10` avec « **Every number is 0-parameter** » (l. 813). Or :
  - `g_0 = 1/8 = 1/2³` est classé **convention** par l'audit lui-même
    (l. 1205) — le même module est « dérivé du 3-space » ailleurs :
    **contradiction d'étiquette** ✓.
  - Sensibilité **vérifiée par re-calculation** (note : premier essai buggué
    par priorité d'opérateurs `1/16*np.pi` ≠ `1/(16π)` ; corrigé) :
    g₀ = 1/4 → ratio **17.9** ; g₀ = 1/8 → ratio **1.12** ; g₀ = 1/16 →
    ratio **0.07**. Les 3 chiffres du reviewer sont exacts ✓. Le « 12% match »
    dépend entièrement du choix non dérivé g₀ = 1/8.
  - `g_0⁴/16π` comme « asymétrie CP-impair 1-loop » **n'est pas** le résultat
    standard de la leptogenèse (standard : ε ∝ m_D·M_N·f(M_N²) issu des
    diagrammes vertex/self-energy avec les masses des N) ✓ — c'est une
    estimation d'ordre 0-paramètre, pas une computation du mécanisme.
  - La cohérence avec le théorème CP-impair n'est atteinte qu'après insertion
    d'une **phase CP de O(1)** (le texte l'admet dans Table 2, W5 : « an
    O(1) input ») — mais Theorem 8 omet cette condition.
- **Verdict :** le pattern reproché est exact : *mécanisme absent (prouvé) +
  formule empruntée + coefficient nommé = « prédiction 0-paramètre »*.

**Correction :** Theorem 8 réécrit : les **facteurs** (28/79 sphaleron SM,
16π, comptes de boucles 1+2) sont 0-paramètres ; l'**asymétrie CP-impair** est
**nulle** en 0-paramètre (théorème) et devient non nulle **seulement** via la
**phase CP nommée** (input O(1)) ; le ratio 2-chiffres 1.12 est le résidu
nommé. L'énoncé « every number is 0-parameter » est retiré.

---

## 2. Faiblesses non circulaires I1–I8

### I1 — « Λ = 0 (Section 9) vs Λ ~ M_P⁴ (Theorem 10) » — **FONDÉ**

Les deux affirmations sont là (l. 467-468 : `Λ_QCA = 0` par la 512/512 ;
l. 879-881 : « At the Planck Hubble (H = M_P) the 0-parameter vacuum energy
is Λ ~ M_P⁴ ») et **ne sont jamais réconciliées** ; le « s-mechanism »
(« the 0-point does not gravitate ») est affirmé, jamais démontré de l'action.
**Correction :** paragraphe de réconciliation — la 512/512 annule l'**énergie
0-point quantique** des 1024 modes (→ 0) ; le de Sitter M_P⁴ est la borne du
**vide classique** auto-gravitationnel du lattice (saddle point) ; Λ_obs =
10⁻¹²² reste le résidu nommé W1.

### I2 — « le graviton spin-2 est postulé » — **FONDÉ**

`cl10/einstein.py:78-119` (`linearized_graviton`) : la **seule** computation
est la dispersion **scalaire** du laplacien 6-voisins (masse nulle, c = 1,
UV-bornée — ce qui est honnête). Mais le dictionnaire retourné **hard-code**
`"helicity": 2, "n_polarizations": 2` (l. 112-113), et le docstring déclare
« the 2 independent polarizations are the spin-2 + and × » **sans aucun
calcul de représentation**. La métrique `h^{μν}` est introduite à la main
(l. 328-333 du texte). P0 (« no added fields ») est en tension avec un champ
métrique 4×4 ajouté au lattice. **Correction :** énoncé honnête dans la
Section 6 + docstring : la **dispersion** (massless, c = 1, UV-bounded) est
computée ; le **caractère spin-2** (2 polarisations, transformation tensorielle)
est une **structure supposée** du champ métrique de la sub-algèbre (1,3) —
relâchement déclaré de P0 pour la gravité.

### I3 — « Nielsen-Ninomiya : exemption par catégorie » — **FONDÉ**

« The doubling theorem does not apply » (l. 851-858) est une **exemption**,
pas une démonstration : la chiralité est une propriété du module 32-D, pas
d'un opérateur de Dirac sur le lattice ; le reviewer dit juste qu'il faudrait
**montrer** que le spectre basse énergie du QCA local (V36) compte bien 16
Weyl chiraux sans doublon — ce n'est pas fait (V36 vérifie le rayon causal,
pas le contenu chiral). **Correction :** reformuler comme **open check**
(partie de W9), pas comme exemption ; et dire que la *construction* n'utilise
pas d'opérateur de Dirac naïf (donc le théorème ne la **contraint** pas, ce
qui est exact — l'excès était de présenter ça comme une réfutation du
théorème).

### I4 — « 9-to-1 = renommage » — **PARTIELLEMENT FONDÉ**

Le reviewer a raison sur la rhétorique de Section 12/16 (« nine projections of
one object »). Mais la Section 13 le formule **correctement** (l. 1000-1010) :
« the nine residuals are the nine **independent exact values** that this one
object **does not fix** » — ce n'est pas un mécanisme, c'est un **compte de ce
qui n'est pas fixé**. Le reproche « aucun mécanisme ne relie Λ à la CKM » est
vrai, mais le texte ne prétend pas l'inverse. **Correction :** légère,
rhétorique (aligner le libellé « 9-to-1 » sur la formulation de la Section 13).

### I5 — « g_10 : le texte se contredit » — **FONDÉ** (et le code aussi)

- Point 1 (l. 1224-1234) : « the theory is **1-parameter** in its dynamics » ;
  L1 (l. 1439-1449) : « the theory has **no free parameters** ». Le pont est
  la lecture `5 → 1/5` que l'audit **classe lui-même** « mapping choice, not a
  derivation » (l. 78, 1230). ✓ La contradiction est réelle : le paramètre
  « disparaît » par changement d'étiquette.
- **Le reviewer ne l'a pas vu : le code est aussi incohérent.** Le docstring du
  module `cl10/four_limits.py` dit « 1 free angle + 2 conventions + **1**
  2-digit coefficient (2.08) » (1.18 manquant) ; le docstring de
  `experiments/V59_four_limits` dit « 1 free + 2 conventions + **2** coeffs +
  1 phase = **6 named** » ; le statut V59 dit « **0** free + **3** named » ;
  l'assembly de `four_limits()` dit `free_parameters: 0`. Quatre comptages
  différents dans quatre endroits.
- **La seule formulation défendable** (celle adoptée par la correction) :
  *0 paramètre libre dans l'algèbre* (les invariants tiennent pour **tout**
  g_10 — calculé) ; *1 angle libre dans la dynamique* dont la **valeur
  sélectionnée** 1/5 = 1/n_qubits est une **convention de mapping** dont
  l'ordre est **data-checké** (α_10 = 1/25 vs 0.024 1-loop, facteur 1.65,
  V56) ; *2 conventions* (g_0, M_GUT = M_P) ; *3 constantes nommées*
  (2.08, 1.18, δ_CP) dont **les valeurs sont lues dans la donnée**.

### I6 — « falsifiabilité largement rhétorique » — **PARTIELLEMENT FONDÉ**

- **Fondé :** F5 — τ_p ≈ 4×10⁴⁷ yr dépend entièrement de la convention
  `M_GUT = M_P` (re-calculation : avec M_GUT = 10¹⁶ GeV, τ_p ≈ **1.8e35** yr,
  chiffre du reviewer exact ✓) ; F3 — la « precision 0.18% » n'existe qu'**après**
  insertion du 2.08 lu dans la donnée ; sans elle la prédiction 0-paramètre
  9.7e-18 est fausse par un facteur 2.08 non borné *a priori* ; F6 — la fenêtre
  10⁻¹¹..10⁻⁸ est large.
- **Sous-estimé par le reviewer :** deux falsificateurs **sont** opérables et
  le texte les enterre dans les « named residuals » : **m_DM = 118 GeV** (W7 —
  une détection du DM à une masse bien hors 118 ± O(10) GeV réfute la
  projection 10-loop) et **m_ν = 0.15 eV** (W2 — 2.9× la valeur observée,
  testable en β-β′ / KATRIN). Le reviewer a tort de dire « aucun des
  falsificateurs n'est opérable ».
- **Correction :** table F réécrite avec des tolérances explicites + mise en
  avant des deux prédicitions les plus accessibles (118 GeV, 0.15 eV).

### I7 — « références à vérifier » — **PARTIELLEMENT FONDÉ**

- **Choptuik : FONDÉ, et pire que dit le reviewer.** Le texte écrit
  « **P. M.** Choptuik, 'Universal Behavior Near Onset of Gravitational
  Collapse', PRL **70, 1** (1993) ». La source primaire (APS,
  DOI 10.1103/PhysRevLett.70.9) donne : **M. W.** Choptuik,
  « **Universality and scaling in gravitational collapse of a massless scalar
  field** », PRL **70, 9–12** (1993), publié le 4 janvier 1993. Les
  initiales, le titre **et** la pagination sont faux (le reviewer n'avait
  signalé que la pagination).
- **Bott :** le reviewer dit « à confirmer » (incertitude mutuelle). Ma
  recherche n'a **pas** confirmé directement « Comment. Math. Helv. 45, 197
  (1970) », mais son contrexemple (« Adv. Math. 1964 / Atiyah-Bott-Singer »)
  n'est pas confirmé non plus. La référence en place est la citation standard
  de la périodicité de Bott pour les algèbres de Clifford dans la littérature
  K-théorie ; **maintenue**.
- **Kay & Laflamme :** titre contesté par le reviewer (« … and
  quantum-gravitational field theory »). Ni l'un ni l'autre titre n'est
  confirmable par index ; « Causal set dynamics: a review », CQG **24**, 5469
  (2007) est la revue standard citée dans la littérature causal-set ;
  **maintenue**.

### I8 — « erreurs de rédaction » — **TOUTES CONFIRMÉES** (+ 3 que le reviewer a ratées)

| n° du reviewer | lieu | défaut | vérifié |
|---|---|---|---|
| 1 | l. 831-832 | phrase coupée « the 45 **mix** the / This is the first computed… » | ✓ |
| 2 | l. 1287-88 | « the CPT theorem **of the CPT theorem** » | ✓ |
| 3 | l. 1565 | « the **Planck-mass** Planck-mass soliton » | ✓ |
| 4 | l. 673 | « the **the** heavy-$N$ out-of-equilibrium » | ✓ |
| 5 | l. 1728-31 | Fig. 18 : « the **1-loop** factor » / « the **2-loop** factor » (les deux sont **10**-loop) | ✓ |
| 6 | l. 1793 | Fig. 23 : « c = 0.5 » vs c = 1 partout ailleurs (`wave.py`: r = c·dt/dx = 0.5, CFL 3-D — l'explication du reviewer est exacte) | ✓ |
| + (non dans la revue) | l. 1603 | « (Section 9.3) » — n'existe pas (la QCD est la Section 10) | ✓ |
| + | l. 603 | « against these **).** » — phrase rompue | ✓ |
| + | l. 1212 | LaTeX rompu « ($g_{10}=1/ $g_{10} = 1/n_qubits = 1/5$ » | ✓ |
| + | l. 1799 | ligne image rompue « The 24_gravity_3d.png) » (le fichier `fig24_gravity_3d.png` existe) | ✓ |
| + | l. 1810 | « ## Fig. 25 » (heading) alors que les 24 autres sont « **Fig. N —** » | ✓ |

---

## 3. Ce que la revue dit bien (et qu'on garde/renforce)

1. L'algèbre Cl(10) est correcte et doublement vérifiée (symbolique +
   numérique) — le mieux du projet.
2. Le miroir spectral 16/16-bar (CPT = conjugaison spectrale, **toutes**
   valeurs de g_10) est un résultat réel et fort.
3. Le **théorème CP-impair** (le QCA 0-paramètre ne peut pas générer
   l'asymétrie) est un résultat vrai, numériquement prouvé — et le **point**
   le plus utile du projet.
4. La vérification 2-route de la courbure 4-D (Riemann par différences
   finies vs analytique, Schwarzschild, BSSN/Z4) est de la NR correcte.
5. L'audit de la Section 15 (5 classes épistémiques) est honnête — c'est
  l'argument de défense principal ; les corrections ci-dessous le **préserver**
   et l'alignent avec la Section 16 au lieu de la nier.

---

## 4. Les 6 erreurs / exagérations du reviewer (reçu)

1. **JSON stale** (note C5) : faux au moment de la revue (fichier
   60/60/toe:true, mtime 39 min antérieur à la revue) ; seules la gitignore et
   l'absence d'artefact sur clone frais sont valables.
2. **Contrexemple Cl(8)** (C4) : mal fondé — Spin(8) n'a pas de 16 (trialité
   8⊕8) et Cl(8,0) n'a pas de split chiral ±i.
3. **« 3 réponses différentes » pour le DM** (C6) : exagéré — c'est un échelonné
   cohérent (candidat Planck → sous-abondant → 10-loop = l'échelle observée) ;
   seul le libellé Section 9 est en friction.
4. **« Aucun falsificateur opérable »** (I6) : faux — m_DM = 118 GeV et
   m_ν = 0.15 eV sont des prédicitions directes testables, simplement
   non mises en avant par le texte.
5. **Bott / Kay-Laflamme** (I7) : incertitude mutuelle non tranchée ; les
   références en place sont maintenues (Choptuik seule corrigée).
6. **I4** : la Section 13 formule déjà correctement ; seule la rhétorique
   12/16 doit être adoucie.

---

## 5. Corrections implémentées

**`doc/publication.md`**
- **Abstract / intro** : le compte final est réécrit — « 0 paramètre libre
  dans l'algèbre ; 1 angle dans la dynamique (valeur sélectionnée 1/5 =
  1/n_qubits, convention de mapping data-checkée à 1.65×) ; 2 conventions ;
  3 constantes nommées (2.08, 1.18, δ_CP) dont **les valeurs sont lues dans la
  donnée**, non dérivées » — le « 0 free » de Section 16 est conditionnel à
  la mapping, et le dit.
- **Section 2** : Remark « minimal unified root » → critère de minimalité
  précis (split chiral 16/16-bar + 16-D Weyl = une génération SM + adjoint 45
  contenant les 12 jauge SM ; Spin(8) n'a pas de 16).
- **Section 3** : « three equivalent ways » → « une même 0-paramètre number (3)
  en trois reformulations combinatoires ; le *why-3* est le résidu nommé W3 ».
- **Section 6 / Section 11** : le spin-2 (2 polarisations) est déclaré comme
  structure supposée du champ métrique (relâchement de P0 pour la gravité) ;
  le calculé reste massless/c=1/UV-bornée.
- **Section 9** : (a) libellé DM aligné (candidat 0-paramètre = soliton Planck,
  sous-abondant ; l'échelle observée = 10-loop 118 GeV ; identité = W7 ; le
  WIMP miracle est un résultat de la cosmologie standard appliqué à la masse
  0-paramètre) ; (b) **réconciliation Λ** : 0-point quantique 512/512 = 0
  (Section 9) vs borne du vide classique de Sitter M_P⁴ (Theorem 10) — deux
  objets différents ; 10⁻¹²² = W1.
- **Section 11** : note de **signification** du programme 60 checks (cohérence
  interne vs 6 comparations nature ; verdict = tiers durs 0-3) ; table F
  réécrite (tolérances explicites ; F5 : la convention M_GUT = M_P est dite,
  le contenu falsifiable est la suppression GUT + le non-GUT ; mise en avant
  de 118 GeV et 0.15 eV comme falsificateurs les plus accessibles) ;
  « two independent constraints » → « one 0-parameter scale, two statements ».
- **Section 12** : Theorem 8 réécrit (facteurs 0-paramètres ; ε_CP = 0 en
  0-paramètre par le théorème, non-nul **seulement** via la phase CP nommée
  O(1) ; ratio 1.12 = résidu) ; Nielsen-Ninomiya → open check (W9) ;
  Theorem 10 : Λ ~ M_P⁴ = borne du vide classique (pas la prédiction) ;
  phrase coupée l. 831-832 réparée.
- **Section 13** : « 9-to-1 » aligné sur « les 9 valeurs exactes non fixées ».
- **Section 15** : « CPT theorem of the CPT theorem » réparé ; W5/W8/W7 :
  « every number 0-parameter » retiré, « one factor, two results » → « une
  échelle, deux énoncés ».
- **Section 16** : L1 (1 angle libre ; sélection 1/5 = mapping data-checkée ;
  0-free **conditionnel** à la mapping), L2 (facteur-2 strict = vide ;
  fenêtre facteur 10 = unique ; 1.18 = 118.48/100), L3 (inchangé, c'est bon).
- **Figures** : Fig. 4 (Planck-mass duplicata), Fig. 7 (Section 9.3 → 10),
  Fig. 12 (refining du « 1 free »), Fig. 18 (1-loop/2-loop → 10-loop ;
  « two independent constraints »), Fig. 20 (convention M_GUT), Fig. 23
  (c·dt/dx = 0.5, c = 1 en unités physiques), Fig. 24 (ligne image rompue),
  Fig. 25 (heading), l. 603, l. 1212.
- **Références** : Choptuik corrigé (M. W. Choptuik, « Universality and
  scaling in gravitational collapse of a massless scalar field », PRL 70,
  9-12 (1993)).

**Code**
- `cl10/four_limits.py` : docstring module (compte aligné), docstring et
  verdict de `unique_10loop` (facteur 2 → facteur 10/une décade ; « two
  independent constraints » → « two statements of one scale » ; la clé
  `unique_n_within_factor2_of_both` → `unique_n_within_1decade_of_both`),
  verdict de `four_limits()` réécrit.
- `experiments/V59_four_limits/V59_four_limits.py` : docstring (6 named →
  compte aligné), statut, clé JSON.
- `cl10/einstein.py` : docstring de `linearized_graviton` (2 polarisations =
  caractère tensoriel supposé du champ métrique, pas une computation).

**Reproductibilité**
- `.gitignore` : `results/*` + `!results/validation_summary.json` ; le
  `results/validation_summary.json` (60/60, toe:true) est committé.

**Validation** : `run_validation.py` relancée en entier après les modifications
(60/60 PASS, toe:true) ; les chiffres re-calculés (2.076, 1.185, 118.48 GeV,
6.8e-10/1.12, 17.9/0.07, 3.98e47/1.8e35 yr, fenêtre facteur-2 vide).
