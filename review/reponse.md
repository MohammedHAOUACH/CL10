# Réponse au rapporteur GLM

**Manuscrit** : `doc/final_publication.md`
**Rapport visé** : `review/glmreview.md` (rejet avec possibilité de resoumission, 3/10)
**Méthode** : chaque affirmation du rapport a été reprise au manuscrit (numéro de ligne) ou au
dépôt (script, JSON de résultats, Crossref). Rien n'est accepté ou refusé sans preuve.

## 0. Verdicts

| point | verdict | statut |
|---|---|---|
| M1 échec du cœur empirique | **vrai comme description, faux en arithmétique** | rien à corriger au manuscrit |
| M2 critère des paramètres / axiome de fermeture | **partiellement vrai** | décision d'auteur (Q1, Q8) |
| M3 gravité sur-vendue | **partiellement faux** | reformulation du titre = décision d'auteur |
| M4 facteurs de boucle | **vrai sur les faits, faux sur l'accusation de masque** | rien à corriger |
| M5 section 15 manquante | **vrai** | **corrigé** |
| M6 doubleurs non résolus | **partiellement faux** (§8.1 avertit déjà) | rien à corriger |
| M7 couche jaugée déconnectée | **vrai sur les chiffres** ; la requalification est un choix éditorial | discuté |
| S1 légendes Fig. 1 / Fig. 2 | **vrai** | **corrigé** |
| S2 mot « CPT » | **partiellement vrai** | **corrigé** |
| S3 références | **vrai sur trois défauts** ; deux des références *suggérées* n'existent pas | **corrigé** |
| S4 chiffres en double emploi | **faux sur Table 4**, **déjà dit** pour $\tau_p$ et $\kappa_s$ | rien à corriger |
| S5 « nine-to-one » | **faux** (paragraphe, pas sous-section) | rien à corriger |
| S6 estimation baryonique | **vrai** ; suppression = choix éditorial | proposé |
| S7 franchise du secteur saveurs | **vrai** | **corrigé** |
| S8 figures 24/25 et lignes tronquées | **faux sur les figures et sur le PDF**, vrai sur la longueur des lignes | rien à corriger |
| S9 tolérances | **partiellement vrai** | **corrigé** |
| S10 titre | **vrai** | décision d'auteur, deux formulations proposées |

## 1. M1 — vrai comme description, faux en arithmétique

Tous les chiffres cités sont au manuscrit, et ce sont ceux que le manuscrit lui-même
avance : facteur $5.8$ et pôle de Landau à $2.2\times10^{9}$ GeV (abstract, l. 14 ;
§13.1 l. 562), abondance $0.798$ / $0.0246$ / $2.3\times10^{-5}$ / $7.6\times10^{-4}$
(Table 5 l. 683, §14.2 l. 731), $\sum m_\nu$ (l. 684), CKM identité (l. 696), tension
$5.4\times10^{40}$ (l. 14). Le manuscrit ne prétend pas le contraire : Table 1 marque le
critère des paramètres « **not met** » (l. 55).

Le décompte du rapport est faux. Table 5 compte **12** lignes à valeur mesurée, dont
**3** marquées « fails » ($\Omega_{\rm DM}$ l. 683, $\alpha_i(M_Z)$ l. 688, CKM l. 696),
**1** « excluded » (spectre du désert l. 693) et **1** « in tension » (neutrino l. 684) —
pas « 6 entrées, 4 fails ». La requalification en « cadre de construction » est ce que la
section 15 (audit de provenance) fait déjà, classe par classe.

## 2. M2 — l'axiome de fermeture

Le rapport a raison sur un point et se trompe sur deux.

- **Raison** : la fermeture est un axiome. Le manuscrit le nomme ainsi (« Call this the
  **closure axiom** », §9 l. 299) et n'en donne pas de dérivation. Sa motivation est
  énoncée (l'ontologie d'automate exige une mise à jour algébrique, pas analytique) ; elle
  n'est pas démontrée. C'est O-adjacent, et la réponse Q1 ci-dessous le dit sans détour.
- **Tort 1** : le manuscrit ne présente pas la fermeture comme un théorème qui « répare »
  le critère. Table 1 (l. 55) : « *not met* … The closure axiom removes the one free
  continuous parameter of the dynamics **at the price of moving that parameter's two
  empirical contacts** ». §9 verdict (l. 323) énumère ce que la théorie *perd*.
- **Tort 2** : « le papier contient deux théories sans critère interne pour choisir » est
  exact comme état de fait, mais le papier les sépare explicitement partout : l. 194
  (« Two readings of it coexist and are kept apart throughout »), Table 5 ligne
  $\alpha_{10}$ qui donne les deux valeurs et leurs deux statuts (l. 687), ligne $\tau_p$
  qui donne les trois branches et leurs statuts (l. 686). Le critère de choix n'existe pas
  — c'est un résidu, pas une dissimulation.

Chiffres vérifiés : $\pi/64 = 0.0491$, facteur $1.226$ (l. 687), $n^* = 208$
(l. 319, 323, 509, 757).

## 3. M3 — la gravité

**Faux sur la troisième exigence.** Le rapport demande que l'abstract mentionne l'absence
de spin-2 sans masse. L'abstract le dit, dans le paragraphe même qu'il cite : « *the
implemented step has no exactly massless one-particle state and therefore no graviton
candidate* » (l. 14). La mention demandée est déjà là.

**Vrai sur la formulation.** « *a lattice gravitational sector that reproduces general
relativity … to machine precision* » (l. 10) est vrai au sens où §11.3 valide
numériquement un secteur qui *implémente* la RG (Table 3, l. 412), mais la phrase se lit
comme une dérivation. Le correctif est une reformulation, pas un ajout de contenu ; il
dépend du choix de titre (S10).

**Ajout post-rapport (calcul nouveau, pas une réponse au rapporteur).** La question
« le graviton est-il ailleurs dans le spectre ? » a depuis reçu une réponse calculée, au
centre de zone où le petit groupe est $O_h$ tout entier et où un spin 2 de réseau est un
multiplet sans gap portant le quintet symétrique sans trace $E_g \oplus T_{2g}$ :
- le facteur de liaison de 8 dimensions, qui porte à lui seul l'action cubique, est
  $A_{1g} \oplus T_{2g} \oplus A_{2u} \oplus T_{1u}$ : il contient $T_{2g}$, jamais $E_g$ ;
- $\Gamma$ est exactement sans gap en $g_{10} = \pi/8, \pi/6, \pi/4, \pi/3$ (jamais en
  $1/5$, gapé de $3.65\times10^{-2}$) ; en chacun de ces angles son espace nul est
  covariant cubique et sans $E_g$ — $16\times(A_{1g}\oplus T_{2g})$, dimension 64, en
  $\pi/4$ — et sa dispersion est plate, $\omega/|\mathbf q|\to0$ ;
- l'espace symétrique de rang 2 des directions de liaison porte, lui,
  $A_{1g}\oplus E_g\oplus T_{2g}$ : le quintet existe sur le réseau, il est dans la couche
  de liaison directionnelle (§10), pas dans le lien qubit.

« Spin 2 est posé, non dérivé » devient donc un **énoncé calculé à deux obstructions**
(absence de $E_g$ à tout couplage ; bandes plates aux couplages sans gap) muni d'une voie
de sortie nommée. §8.3, §11.1 et le paragraphe « Where the construction is weakest » ont
été mis à jour en ce sens ; la revendication de titre (S10) reste une décision d'auteur.

**Le recensement « pas de mode exactement sans gap » était faux, et le rapporteur ne
l'a pas vu non plus.** V128 l'avait établi sur une grille $15^3$ (minimum de zone
$5.1\times10^{-4}$). V139 raffine à l'intérieur des trois plans nodaux : le pas de $1/5$
a une **orbite exactement sans gap de 24 points**, représentant
$q_0 = (-2.1890561, -1.4445961, 0)$, gap $4\times10^{-14}$, résidu au plancher de bruit
sous raffinement alors qu'un déplacement de $5\times10^{-2}$ dans le plan coûte
$8.8\times10^{-3}$ ; petit groupe d'ordre $2$ ; dispersion linéaire dans le plan,
$v\in[0.031, 0.219]$ (anisotropie $0.857$), gap linéaire hors du plan. Donc : des
excitations exactement sans masse existent à $g_{10}=1/5$ — l'abstract le niait, il est
corrigé — et le no-go spin-2 tient quand même : un petit groupe d'ordre 2 n'a que des
irreps de dimension 1, donc pas de paire d'hélicités. L'anisotropie du pas est mesurée
$0.857$ sur la branche exactement sans masse, et non $0.252$ (qui était celle des nœuds M
quasi sans gap). §8.3, Table 5 et Table 7 (O9) suivent.

## 4. M4 — les facteurs de boucle

Les faits du rapport sont exacts, vérifiés au manuscrit : $F_{10} = 0.019186$ mesuré
(l. 461, 475), $3F_{10}^{10} = 2.0881$ contre $2.0765$ (l. 479, Table 5 l. 681),
$55326$ mots sur $3^{10} = 59049$, $r = 0.2232$, $\mu = 2.981$ (l. 475), et l'extrapolation
au rang 17 qui manque l'exposant d'un facteur $0.066$ ($7.9\times10^{-31}$ contre
$1.2\times10^{-29}$, l. 475).

**L'accusation de masque est fausse.** Le manuscrit ne range pas ce nombre sous « fitted » :
il crée pour lui une classe, *conditional*, et nomme la condition dans la table même :
« the factor $3$ the one-split-node-per-family choice » (Table 6, l. 751), suivie de « *The
conditional class is neither a fitted number nor an algebraic one: its ingredients are all
measured, but the node-set choice that supplies the factor $3$ is an input* » (l. 754).
Le dilemme que pose le rapport (« soit ajusté, soit dérivé ») est exactement le tiers exclu
que le manuscrit nomme. La critique survivante est de fond, pas de comptabilité : le
mécanisme ne reproduit pas l'exposant 17 (facteur $0.066$), et c'est écrit.

## 5. M5 — la section 15 : vrai, corrigé

`git log -S "## 15." -- doc/final_publication.md` est vide : la section 15 n'a **jamais**
existé ; le trou de numérotation (14 → 16) était permanent. Le rapport a aussi raison sur
le chemin de l'outil : il est `review/tools/verify_refs.py`, pas `review/verify_refs.py`.

**Correction appliquée** (`review/tools/apply_ref_fixes.py`) : 16→15, 17→16, 18→17,
19→18, 20→19, avec les 9 renvois internes remappés (« Section 16 »×6 → « Section 15 »,
« Section 17 »×2 → « Section 16 », « Sections 17–18 » → « Sections 16–17 »).
Vérification dans le PDF rendu : en-têtes **XV. THE PARAMETER AUDIT / XVI. OPEN PROBLEMS /
XVII. FALSIFIERS / XVIII. DISCUSSION / XIX. CONCLUSION**, contigus, sans trou.

**`verify_refs.py` exécuté** : 40/53 conformes, 8 signalements. Triage : ce sont des
artefacts de rapprochement Crossref, pas des erreurs bibliographiques — PDG RPP
(*Phys. Rev. D* **110**, 03001 (2024)) et Planck 2018 (*Astron. Astrophys.* **641**, A6
(2020)) sont des volumes de revue que Crossref indexe sans le volume revendiqué ; Sakharov
1967 (*JETP Lett.* **5**, 24) est résolu par Crossref sur la traduction de 1998 ; Super-K
$\tau_p$ (PRD **90**, 072005 (2014)) et LZ (PRL **135**, 011802 (2025)) sont corrects, le
titre de la notice Crossref la plus proche diffère. L'outil doit apprendre les articles de
revue et les traductions ; aucun des 8 n'a été confirmé comme erreur.

## 6. M6 — les doubleurs

Le rapport affirme que §8.1 présente le remède sans signaler l'exclusion, avant que §13.1
ne l'exclue. §8.1 le signale, une phrase après la phrase citée : « *the meeting scale the
mechanism then predicts is in conflict with proton decay (Section 13.1)* » (l. 260), et
l. 262 nomme ce qui reste entrée (le choix de l'ensemble de nœuds). Le fond — ni
théoriquement (trois plans nodaux, pas de cône, §8.2) ni phénoménologiquement (le point de
rencontre est exclu) — est la conclusion du manuscrit, pas une thèse du rapporteur.

## 7. M7 — la couche jaugée

Chiffres exacts : tension $5.4\times10^{40}$ (l. 14), déficit de Bekenstein $4\ln 45 =
15.23$ (l. 691), coefficient corrigé de sommet $1.65$ nats/bond, facteur $6.6$ (l. 692),
« no element of the current construction produces $1/4$ » (l. 351). La recommandation de
présenter §10 comme résultat de construction indépendant est un choix de mise en avant :
le manuscrit classe déjà toutes ses lignes « computed » (jamais « phenomenology ») et
§10.2 énonce que $\Lambda_{\rm QCD}$ n'est pas produit (il faudrait $n\simeq11.6$,
non entier, l. 543).

## 8. Secondaires

**S1 — vrai, corrigé.** Le texte le disait (l. 160 : la réduction est « an input of the
interpretation (O9) » ; l. 166 : « The atom is an **input** »), les légendes non.
Fig. 1 (c) porte maintenant « the projection is an input of the interpretation (O9), not a
consequence of the unitary step » ; Fig. 2 « **supplied as an input** … The figure tests
the photon sector, not hydrogen ».

**S2 — partiellement vrai, corrigé.** « CPT » n'apparaît que deux fois, jamais dans
l'abstract ni dans une légende. §7.2 (l. 214) dit déjà que le mot surestime. Le §14.3
restait trop fort ; il dit maintenant : « Their degeneracy is that mirror and nothing more:
as Section 7.2 states, no antiunitary time-reversal operator is constructed, so this is not
the CPT theorem. »

**S3 — vrai sur trois défauts, corrigé ; deux des références suggérées n'existent pas.**

| défaut signalé | vérification | correction |
|---|---|---|
| réf. 2 ('t Hooft 1976, pseudoparticle) citée dans le groupe automates/causal sets | confirmé : jamais citée ailleurs que dans `[1–5]` | remplacée par K. G. Wilson, *Confinement of quarks*, **Phys. Rev. D 10**, 2445 (1974) — la phrase cite « lattice gauge theory » sans jamais la référencer |
| « Lieb–Robinson bound » utilisé deux fois sans référence | confirmé | E. H. Lieb & D. W. Robinson, *The finite group velocity of quantum spin systems*, **Commun. Math. Phys. 28**, 251 (1972), cité aux deux endroits |
| littérature Clifford-GUT absente | partiellement : la réf. Castro (ancienne 12) y était ; Atiyah–Bott–Shapiro et Penrose–Rindler n'étaient en effet pas des références d'unification | N. Furey, *Generations: Three Prints, in Colour*, **JHEP 10** (2014) 046, arXiv:1405.4601 ajoutée ; le groupe devient `[11–14]` |
| « H. Hisano » | confirmé | **Y.** Hisano (réf. 48) |

Le rapport ne doit pas être suivi sur ses propres suggestions : « Brower–Furman–Hughes »
et « Chen–Kaplan » ne résolvent à aucun article (Crossref ne renvoie Brower qu'avec
d'autres coauteurs ; aucun Chen–Kaplan correspondant). « Gamba–Marti–Smoladello » n'a pas
été retrouvé. Insérer ces références en l'état introduirait des citations fausses.

Les marqueurs de citation du corps ont été remappés à travers les deux insertions :
55 entrées contiguës, 55 citées, 0 orpheline, 0 marqueur hors plage ; les crochets
mathématiques (`[0,\pi/2]`, `[10^6,10^{19}]`, `[0.150,0.201]`, …) intacts.

**S4 — faux ou déjà dit.** Table 4 marque déjà la ligne électrofaible « **measured input
$v$** » (l. 539) : ce que demande le rapport y figure. Les trois $\tau_p$ sont déjà dans
une seule ligne de Table 5 avec leur statut par branche (l. 686). Pour les $\kappa_s$, le
texte dit textuellement « three meeting conditions for six unknowns — a solution family,
not a unique prediction » (l. 611) : les 3 degrés de liberté sont énoncés ; le compte
« 11 coefficients ajustés » de Table 6 compte des *entrées nommées*, ce qui est une autre
grandeur, et le manuscrit le dit.

**S5 — faux.** L'observation « nine-to-one » n'est pas une sous-section nommée : c'est un
paragraphe de §15 (l. 759) qui se disqualifie lui-même — « *That is a bookkeeping
observation … and it is not a derivation of any of their values. It should not be read as a
reduction of the theory's difficulty.* »

**S6 — vrai, décision d'editorial.** L'estimation $\Delta B = 6.8\times10^{-10}$ est dans
le corps ; Table 5 la qualifie déjà de « coincidental product » (l. 685). Supprimer du
corps et déplacer en annexe est défendable ; le manuscrit garde la version actuelle parce
que l'estimation est ce qui *mesure* l'écart entre les deux facteurs faux. À trancher.

**S7 — vrai, corrigé.** §13.5 se termine maintenant par : « **The flavour accounting is
therefore: one angle derived out of six, no mass derived, CKM predicted to the identity.** »

**S8 — faux sur les figures, faux sur le PDF, vrai sur la longueur.** Les figures 24 et 25
sont `fig24_w3_forcing` et `fig25_residual_attacks` : des graphiques, pas des tableaux
O1–O9. Le tableau O1–O9 est la **Table 7** de §16. Il est rendu sur deux pages réservées
(pages 35–36 du PDF, 555 pt et 528 pt de règle à règle dans la zone de texte), avec les 9
lignes complètes vérifiées ligne à ligne dans le PDF (les 71 lignes des 8 tableaux sont
toutes présentes). Les lignes longues existent (2231 à 2525 caractères aux lignes 311,
495, 772, 776) mais ne sont **pas tronquées** dans la source. La restructuration en
sous-sections reste une préférence éditoriale ; l'affirmation « le PDF et les revues ne
peuvent pas gérer ces lignes » est contredite par le build.

**S9 — partiellement vrai, corrigé.** §1.4 distingue maintenant les deux régimes : un
résidu $\le 10^{-15}$ sur une identité attendue exactement nulle est une **borne
supérieure**, pas une mesure ($\epsilon_{\rm mach}\simeq2\times10^{-16}$) ; une tolérance de
**recherche** est une résolution, énoncée là où elle est utilisée ($10^{-7}$ pour le scan
des 2726 points de §8.2, $10^{-6}$ pour le seuil d'isotropie), et jamais présentée comme
précision machine.

**S10 — vrai.** Le titre revendique « the Standard Model and gravity » alors que §11 pose
$g_{\mu\nu}$, $G=1$ et la normalisation $1/16\pi$. Décision d'auteur. Formulations exactes,
toutes deux fidèles au contenu :
- *Causal Clifford Theory: one anomaly-free Standard Model generation and a lattice
  gravitational sector from a single Clifford root* ;
- *Causal Clifford Theory: a Clifford-root quantum cellular automaton — Spin(10) content,
  the doubling problem, and the loop-factor ladder*.

## 9. Réponses aux huit questions

1. **Pourquoi un axiome ?** Parce que c'en est un, et le manuscrit l'appelle axiome. Aucun
   principe connu ne force la dynamique fondamentale à être algébrique ; la motivation
   énoncée est ontologique (un automate fondamental a un pas d'ordre fini). Le coût est
   chiffré (Table 1, §9 verdict). Ce qui départagerait les deux lectures n'existe pas
   encore : c'est le résidu, pas un défaut de rédaction.
2. **Nielsen–Ninomiya sur les lignes nodales.** Non refait à ce jour ; le manuscrit dit
   seulement que la convention par nœud est inapplicable et que le comptage naïf donne $-2$
   ou $-6$ selon la route. Un fait structurel est pourtant établi — vérifié par
   construction : les six droites $\{(t s_1, t s_2, 0)\}$, $t\in[-\pi/2,\pi/2]$, deux
   patronages de signes par plan de face, sont chacune **invariantes par inversion comme
   ensemble** ($t\mapsto-t$), et les 12 points M en sont des points — ce qui explique
   l'échec du comptage : la paire de Nielsen–Ninomiya d'un point n'existe pas, elle est à
   l'intérieur de chaque ligne. Le calcul qui répondrait à la
   question est spécifié : pour chaque ligne, sur un cercle de rayon $\rho$ dans le plan
   normal, évaluer le degré de l'application de Bloch restreinte au sous-espace sans gap
   (invariant $\mathbb{Z}$ ou $\mathbb{Z}_2$ de la ligne), puis vérifier que les lignes se
   ferment ou s'annihilent deux à deux. Point d'entrée du dépôt : `cl10.supercell`
   (`step_wilson`, `M_POINTS`) et le harnais de `experiments/.../V118_w9_little_group`.
   **Aucun résultat n'est revendiqué ici** : ce calcul n'a pas été fait.
3. **Que prédit la théorie avec les six $\kappa_s$ ?** Rien d'unique : le manuscrit le dit
   (l. 611). La falsifiabilité n'est pas perdue pour autant, elle se déplace : les
   $\kappa_s$ sont des séparations de seuils mesurables, et le **théorème à décalage égal**
   ($\sum db_1 = \sum db_2$ sur tout multiplet complet, §13.1) est une contrainte qui
   survit à l'extension et a déjà servi à exclure toutes les familles d'ansatz testées
   (résidu minimal $\ge 5.7$, l. 693). F2 porte sur la théorie implémentée ; l'extension est
   un objet nommé et testable.
4. **Le facteur 3.** Aucun mécanisme ne sélectionne « un nœud divisé par famille » ; c'est
   l'entrée nommée O3/O9, et la classe *conditional* de Table 6 existe pour que cela ne
   puisse pas être lu comme une prédiction.
5. **Masse 118 GeV sans porteur.** La prédiction de masse est celle de l'échelle
   $(1/16\pi)^{10}$ ; le porteur est une extension nommée (le scalaire pair/impair miroir),
   et le manuscrit dit des deux états : « no state of the implemented spectrum can be the
   dark matter » (abstract l. 14) et le dilemme portal LZ (§13.2, Fig. 15). F7 falsifie le
   couple (masse, portal) ; la prédiction n'est pas conditionnelle au spectre actuel, elle
   est conditionnelle à l'extension, ce qui est écrit.
6. **La section 15.** Elle n'a jamais existé ; trou de numérotation. Corrigé (renumérotation
   15–19, renvois remappés, PDF XV–XIX).
7. **Tension de corde et confinement QCD.** Oui, incompatible : $5.4\times10^{40}$ fois la
   valeur QCD, et $\Lambda_{\rm QCD}$ n'est pas produit (il faudrait un exposant
   $n\simeq11.6$ non entier, §10.2). « Constructible » et « compatible avec la
   phénoménologie » ne sont pas confondus dans le manuscrit : la couche est classée
   *computed*, jamais *phenomenology*.
8. **Laquelle des deux théories ?** Le manuscrit garde les deux et compare déjà leurs
   contacts dans une table unique (Table 5, lignes $\alpha_{10}$, $\tau_p$, recensement des
   doubleurs). Le choix n'est pas faisable avec les données actuelles : la branche fermée
   gagne le critère des paramètres (zéro nombre continu non fixé) et perd les deux contacts
   empiriques. La réponse honnête est celle du §9 verdict ; la reformulation en « une seule
   théorie » serait un faux.

## 10. Ce qui a été appliqué, et la preuve

| correctif | fichier |
|---|---|
| renumérotation 16→15 … 20→19 + 9 renvois | `doc/final_publication.md` |
| réf. 2 → Wilson ; Lieb–Robinson ajoutée et citée ×2 ; Furey ajoutée (groupe `[11–14]`) ; Y. Hisano | `doc/final_publication.md` |
| légendes Fig. 1 (projection = entrée, O9) et Fig. 2 (atome = entrée) | `doc/final_publication.md` |
| §14.3 : le miroir spectral n'est pas le théorème CPT | `doc/final_publication.md` |
| §1.4 : epsilon machine vs résolution de grille | `doc/final_publication.md` |
| §13.5 : phrase de synthèse du secteur saveurs | `doc/final_publication.md` |
| script de remap (numérotation + références + marqueurs de citation) | `review/tools/apply_ref_fixes.py` |

Build après correctifs : `bash publication/build.sh` → **39 pages**, `0` erreur,
`0` référence indéfinie, `0` `Overfull \vbox`. Références : 55 entrées contiguës, 55
citées, 0 orpheline. Sections rendues : XV, XVI, XVII, XVIII, XIX contiguës.

## 11. Ce que le rapport doit corriger de son côté

- figures 24/25 ne sont pas les tableaux O1–O9 (Table 7 l'est) ;
- « six entrées, quatre fails » dans Table 5 : 12 lignes à valeur mesurée, 3 « fails » ;
- « l'abstract ne dit pas l'absence de spin-2 sans masse » : il le dit (l. 14) ;
- « §8.1 ne signale pas l'exclusion » : il le fait (l. 260) ;
- « Table 4 ne marque pas $v$ comme entrée » : elle le marque (l. 539) ;
- « le PDF ne peut pas gérer les lignes longues » : la Table 7 est rendue sur deux pages
  réservées, toutes lignes présentes ;
- les références que le rapporteur suggère (« Brower–Furman–Hughes », « Chen–Kaplan »,
  « Gamba–Marti–Smoladello ») ne résolvent à aucun article vérifiable ;
- le chemin de l'outil est `review/tools/verify_refs.py`.

## 12. Ce qui reste ouvert après tout cela

Les objections de fond survivent à la vérification, et aucune n'est réglable par de la
rédaction : le désert SM ne fait pas se rencontrer les couplages (facteur $5.8$, pôle de
Landau à $2.2\times10^{9}$ GeV) ; la réparation exige un contenu que l'automate ne fournit
pas tel quel ; l'exposant 17 de l'échafaudage manque d'un facteur $0.066$ ; aucun état du
spectre implémenté n'est le porteur de matière noire ; aucune condition de Sakharov n'est
satisfaite ; le déficit de Bekenstein reste à un facteur $6.6$ ; la constante cosmologique,
le nombre de familles et la hiérarchie des masses ne suivent pas. Le manuscrit les énonce
comme ouverts (Table 6, Table 7) et le désaccord avec le rapporteur porte sur le fait de
savoir si les énoncer suffit à faire un article — pas sur leur existence.
