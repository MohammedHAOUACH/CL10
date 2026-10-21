# Risques de rejet et reste à faire

Analyse en posture de referee, menée le 2026-10-21 sur `doc/final_publication.md`
(54 pages, Table 1, 10 open problems, 147 expériences). Chaque preuve citée a été
relue dans le fichier ou dans l'artefact correspondant.

## 1. Risques rédhibitoires (sélectivité éditoriale)

| Rang | Faiblesse | Preuve dans le papier | Ce qu'écrit le referee | Correction |
|---|---|---|---|---|
| A | Aucun spectre bas-énergie réaliste : aucun fermion chiral léger, dispersif et isotrope n'est produit | §8.3 : à $\pi/4$ l'espace nul de $\Gamma$ est **plat** ($\omega/|\mathbf q|\to0$) et à $1/5$ $\Gamma$ est gappé de $3.65\times10^{-2}$ ; §9 : touchement **quadratique** sur les six droites, exposants $(1,1,2)$ ; §8.2/8.3 : à $1/5$ les seuls modes exactement sans masse ont un little group d'ordre 2 (pas de paire d'hélicité) et une anisotropie $0.857$ ; Table 4 : le secteur électrofaible est une **entrée mesurée** ; Table 5 : aucune paire $\pm1$/$\pm2$ nulle part | « Where is the electron ? Le dédoublement n'est pas supprimé, il est déplacé. Les nombres quantiques internes sont exacts, la réalisation spatio-temporelle des fermions observés n'est jamais établie. » | **O10 nommé** « the low-energy spectrum » + phrase dans le résumé |
| B | L'axiome de fermeture est en tension avec la seule mesure du coefficient qui départage les horloges | §12.1/O8 : l'empreinte fonctionnelle lit $g_{\rm obs}=0.199998$ — l'horloge non fermée $1/5$ à $10^{-5}$ — alors que la classe fermée s'épuise à $\{\pi/2,\pi/4\}$ et rate de $30.0\%$ ; F2 oppose $0.9895$ (branche fermée) à $0.96388$ (hiérarchie) ; « the closure axiom and the loop-hierarchy identification cannot both hold » | « L'axiome central est réfuté par la mesure du papier lui-même. » | Reformuler : l'axiome est une affirmation **structurelle** sur la dynamique sur-site, pas une hypothèse empirique sur l'angle choisi ; F2 est la mesure qui tranche |
| E | Aucune discussion de Bell / de contextuality, alors que le papier revendique un déterminisme exact de la mise à jour | le mot « Bell » (sensible à la casse), « Kochen », « contextual » sont absents ; la projection et la règle de Born sont déclarées entrées (§1.3, §6, §13.1) | « Un réseau unitaire localement déterministe qui prétend fonder la mécanique quantique doit dire d'où viennent les corrélations violant Bell. » | Paragraphe de périmètre dans §6 |
| K | Format : 50 pages, auteur unique sans affiliation, 9 problems ouverts, ~40 résidus nommés, Table 1 concluant « parameters: not met », vocabulaire idiosyncratique | Table 1 ; Table 6 (11 coefficients fittés, 2 conventions, 1 conditionnel) | « Un programme de recherche avec son tableau d'échecs, pas un résultat. » | Découpage en volets (voir §4) |

## 2. Objections de fond

- **Unification.** SU(2) hors de $5.8\times$, pôle SU(3) à $2.2\times10^{9}$ GeV — à
  l'intérieur de l'échelle du réseau — aucune des trois valeurs requises
  $(0.0299,0.0202,0.0191)$ atteinte ; la réparation garde deux seuils fittés,
  abandonne $M_{\rm GUT}=M_P$ et appelle une extension à six $\kappa_s$ fittés.
  Verdict du papier lui-même : « parameters: not met ». Le referee résume :
  un ajustement avec extensions, pas une prédiction.
- **L'échelle $(1/16\pi)^n$.** Exposant choisi par un scan $n=1..20$ ; coefficient
  en écart de $30.0\%$ ; QCD exclu ($n\simeq11.6$, non entier) ; règle des barreaux
  hors entier partout ailleurs ; $m_{\rm DM}=118.6$ GeV et le $\times1/3$
  neutronique sont deux coïncidences ; le second accord $\nu$
  ($0.0084/0.0086$) est explicitement concédé comme coïncidence.
- **Lorentz.** $A\ge0.997$ sur 2726 points, mur de rigidité $8.0\times10^{-3}$,
  « cônes : 0 ». Aucune confrontation aux bornes publiées de violation de Lorentz.
  Paragraphe manquant → ajouté (§3 ci-dessous, action 3).
- **Gravité.** Le contenu est maintenant calculé, pas seulement ajouté : §11.5
  (artefact `V208`) établit par caractères exacts qu'une amplitude de liaison
  unique ne peut pas être la métrique — $\mathrm{bond}(6) = A_{1g}\oplus E_g\oplus
  T_{1u}$ n'a pas le $T_{2g}$ du quintet — et que le hop, vecteur polaire, ne
  déplace pas la métrique au premier ordre (bloc nul à $2.7\times10^{-16}$ dans
  le secteur pair) ; le contenu est forcé, la dynamique non : $34$ auto-couplages
  autorisés, aucun choisi. §11.6 exhibe un cône exactement lightlike sur ce module
  (quintet sur le cône, aucune masse écrivable dans la classe) mais dans une
  famille de $53$ paramètres que la dynamique ne sélectionne pas, et la vitesse
  réellement prédite par la couche construite est $c_{\rm gw}/c = 2\pi/4 =
  1.5708$, exclue par GW170817 : nouveau falsificateur **F8**, déjà mordu, et
  l'horloge $g = 1/2$ qui réparerait tombe dans la zone morte de `V206`
  (échelle de boucles effondrée à $n = 3$). $G = 1$ reste une entrée ; le root
  ne fournit qu'un candidat, $1/(4\ln 45) = 6.567\times10^{-2}$ (§11.7).
  Table 1 : « one object: partial » inchangé.
- **Matière noire.** Le porteur est un secteur singulet ajouté nommément ;
  $\Omega h^2=0.120$ n'est obtenu qu'avec $\lambda_B=4.06\times10^{-14}$ fitté
  (les valeurs fournies par le root donnent $0.798$, $0.0246$, $2.3\times10^{-5}$,
  $7.6\times10^{-4}$).
- **Falsifiabilité surtout négative.** $\tau_p=2.8\times10^{39}$ yr hors d'atteinte,
  la version testable ($1.8\times10^{36}$ yr) étant déclarée convention abandonnée ;
  le barreau nu $\nu$ est déjà exclu ; F3 mord à courte échéance, et F8 mord déjà :
  la couche construite donne $c_{\rm gw}/c = 1.5708$ contre la borne
  $10^{-15}$ de GW170817 — le papier l'écrit au lieu de le masquer.
- **Auditabilité.** 147 expériences auto-produites, tolérances définies par
  l'auteur. Deux griefs de détail sont fermés : la libération est épinglée par
  le tag `v1.0` et le résultat du coin (§8.2, Table 5, O9) a son artefact
  commité, `V207`, avec le scan complet du plan $(w,\theta)$. Deux restent
  ouverts : le DOI archivé (Zenodo 10.5281/zenodo.23321228) et la release
  GitHub sont les instantanés `v1.0.2` / `v1.0.3`, donc une capture antérieure
  au tag `v1.0` que la note de disponibilité nomme ; et les tolérances de
  chaque verdict sont choisies par l'auteur, sans re-run indépendant fourni.

## 3. Actions exécutables dans le dépôt

| # | Action | Fichier(s) | Statut |
|---|---|---|---|
| 1 | Nommer **O10 « the low-energy spectrum »** : Table 7, résumé, §15/§16/§19, cascades de décompte (neuf → dix) | `doc/final_publication.md` | fait |
| 2 | Paragraphe Bell / contextuality / portée du déterminisme (§6) + entrées bibliographiques | `doc/final_publication.md`, `publication/doi.json` | fait |
| 3 | Paragraphe sur les bornes de violation de Lorentz, rattaché à O10 (§8.4) + entrée bibliographique | `doc/final_publication.md`, `publication/doi.json` | fait |
| 4 | Confrontation frontale aux programmes voisins (GUT de Clifford, QCA) : ce que CL10 ajoute et ce qu'il n'ajoute pas (§1.2) | `doc/final_publication.md` | fait |
| 5 | Recadrage de la tension axiome de fermeture / F2 (§9.2, §15, §19) | `doc/final_publication.md` | fait |
| 6 | Expérience commitée pour le zéro de paire du corner (scan $(w,\theta)$, espace nul de six états, commutateurs tangents, isolement, vitesses, planéité) : `V207` | `experiments/dynamics_spectrum/V207_corner_pair_zero/` | fait |
| 7 | Chaîne de décomptes 143 → 144 (127 → 128 publiées ; V207 écrit le verdict sans `tier`, donc 83 expériences hors tier), Table 9 / §8.3 / Table 5 / O9 rattachés à l'artefact V207 | `doc/final_publication.md`, `README.md`, `doc/reference.md`, `experiments/README.md` | fait |
| 8 | Épinglage de la libération : tag git annoté + note de disponibilité qui le nomme | `doc/final_publication.md`, git | fait (DOI Zenodo : action hors dépôt, voir §5) |
| 9 | Docstrings alignés sur leurs artefacts commités : `V187_rigidity_frontier.py` nomme le sommet $1.033$ (la liaison $\mathfrak{su}(2)$ d'indice $(5,8,9)$), le plancher $8.0\times10^{-3}$ et la médiane $-0.28$ — ceux de la re-exécution corrigée à l'horloge $1/5$, que le papier cite aussi ; `V186_cone_inverse.py` (le problème inverse du cône — il n'existe pas de script « rigidity_wall_scan ») ne chiffre son verdict que dans le status, et sa fuite radius-1, Wilson éteint, y est $8.7\times10^{-17}$ (le docstring disait $2\times10^{-16}$) | `experiments/dynamics_spectrum/V186_cone_inverse/V186_cone_inverse.py`, `experiments/dynamics_spectrum/V187_rigidity_frontier/V187_rigidity_frontier.py` | fait |
| 10 | Épinglage explicite de l'horloge, exigé par le tripwire : **25 scripts** posent `CL10_ANGLE = 0.2` avant tout import de `cl10` — les 20 que le tripwire faisait échouer, `V176` (seul de la chaîne V171–V179 à ne pas l'avoir : son artefact commité était passé en horloge fermée et contredisait son propre verdict documenté, $F(0)=5.000000$, poids $E_g+T_{2g}=5/6$), et `V174`, `V175`, `V179`, `V182`, qui lisaient le défaut en silence depuis l'intérieur des chaînes 1/5 — pour `V175` l'angle d'horloge *est* le couplage $\eta$ ($0.7854 \to 0.2$), donc tout son spectre de diffusion bouge ; le recensement `v70_regression` de `V91` est rendu conditionnel à la branche (exactement sans masse à $\pi/4$, les douze portant le même gap doubleur $2.003\times10^{-4}$ à $1/5$) ; re-exécution complète des 144 ($144/144$, 19 des 20 artefacts identiques au bit de physique près) et ré-ancrage des planchers machine (V196, V185) et des digits de `experiments/README.md` et des deux docs d'exploration (V172, V175, V177, V179) ; tag `v1.0` déplacé sur le commit final | `experiments/…` (25 scripts), `doc/final_publication.md`, `doc/paper1_core.md`, `experiments/README.md`, `doc/alt_gauged_exploration.md`, `doc/alt_interactions_exploration.md`, git | fait |
| 11 | Épinglage du PDF journal : `publication/publication.pdf` (14:29) était plus vieux que le maître `doc/final_publication.md` (16:11) — rebuild par `publication/build.sh` après la dernière passe de texte (52 pages ; le PDF porte $1.033$, $(5,8,9)$, le mur $8.0\times10^{-3}$, $2.699$, $0.019685$, O10, Kochen, la note de disponibilité nommant `v1.0`, et le décompte 144) | `publication/publication.pdf` | fait |
| 12 | Expérience commitée sur le secteur gravitationnel : contenu de la métrique (caractères exacts, le bilinéaire $\mathrm{bond}\otimes\mathrm{Sym}^2$ de $36$ états, commutant de dimension $34$), mur de parité hop (polaire, bloc pair nul) / dressing (axial, spectre $g\{0,\pm1,\pm2\}$), générateur de cône exactement lightlike commité `cone_generator.json` (régulier Clifford à $1.4\times10^{-13}$, sous-espace conique de dimension $8$, contenu $E\oplus T_1\oplus T_2$, poids du quintet $5.0000$, $4$ branches de pente $0$, $\mathrm{M}(0) = I$ exact), noyau quasi-local, dimension de la famille ($53$ dans un espace covariant de $160$, Jacobien de rang $107$), distance des triples implémentés ($0.0029$ / $0.0092$), zone morte de l'horloge $1/2$, et $G = 1/(4\ln 45)$ : `V208` ; sections §11.5–§11.7, falsificateur **F8**, ligne de Table 9, référence [60] | `experiments/gravity/V208_metric_bilinear_cone/`, `doc/final_publication.md` | fait |
| 13 | Chaîne de décomptes 145 → 147 (129 → 131 publiées ; `V209`, `V210` écrivent leur `tier`, donc $84$ expériences hors tier) | `README.md`, `doc/reference.md`, `experiments/README.md`, `doc/final_publication.md` (Annexe B, note de disponibilité) | fait |
| 14 | Ré-ancrage mécanique de tous les nombres du maître : balayage des nombres de `doc/final_publication.md` contre les 61 356 feuilles des verdicts (artefacts régénérés, suite rejouée : 147 PASS, palier dur atteint). Loi de composition mesurée là où la fonctionnelle n'est pas dégénérée — horloge ouverte $g=1/5$, dérive nulle entre 30 et 50 chiffres (`V210` : $F_{170}=1.3089\times10^{-21}$, $(F_{10})^{17}=6.4697\times10^{-30}$, rapport $2.02\times10^{8}$) ; les rapports de l'horloge racine à ces longueurs sont des planchers de la classe $4k{+}2$ ($0/0$ à l'angle de fermeture, `V209`) et sont retirés du texte. Pôles $\mathrm{U}(1)$ de la Fig. 4 calculés par la bibliothèque ($2.2\times10^{41}$ GeV à 1 boucle, $2.5\times10^{12}$ GeV à 2 boucles, `V38` ; docstring corrigée) ; lecture abandonnée $M_{\rm GUT}=M_P$ : $3.98\times10^{47}$ yr et échelle exclue $6.05\times10^{15}$ GeV (`V64`) ; rendement purement gravitationnel $\Gamma_{\rm grav}/H = 5.3\times10^{-53}$ (`V195`) ; plafond $\sum m_\nu < 0.0642$ eV (DESI DR2 + Planck PR4 + lentilles, réf. [61]) ; table des DOI corrigée (les entrées manuelles de `fetch_dois.py` étaient décalées de deux références — Castro pointait vers la PDG) et [61] résolu : $10.1103/\texttt{w9pk-xsk7}$ | `doc/final_publication.md`, `cl10/unification.py`, `experiments/couplings_unification/{V38,V64}`, `publication/{fetch_dois.py,doi.json}`, `publication/publication.pdf` (55 pages), `doc/final_publication.pdf` (67 pages) | fait |

## 4. Ce qui est fort et doit rester en tête d'affiche

Ce sont des négatifs structuraux prouvés, pas des ajustements : annulation des huit
anomalies, miroir spectral, $\det U=1$ vrais **pour toute valeur du couplage** ;
rayon causal exactement 1 mesuré ; le no-go de parité (Prop. 6) avec son contre-exemple
explicite de cône hors de la classe ($A=10^{-15}$, $\rho_{\rm comm}=1/\sqrt2$) et son
mur de rigidité mesuré sur 5292 points ; le théorème d'annulation des invariants CP-impairs
et son périmètre (la couche de liaison) ; le théorème equal-shift ; l'obstruction de centre
($\dim\ker=0$) et son élimination exacte par le lien régulier ; la loi d'aire exacte $\ln45$
avec son déficit de Bekenstein chiffré ; la règle de sélection impaire exacte, l'énumération
$55326/3^{10}$, et le fait que $n^*=208$ est un théorème du fonctionnel (polynôme cosinus,
63 classes de phase) et non un angle fitté.

## 5. Ce qui reste hors du dépôt (décisions auteur)

- **DOI et version épinglée publiquement.** Le DOI Zenodo est minté
  (concept 10.5281/zenodo.23321106, version 10.5281/zenodo.23321228) et la
  release GitHub existe, mais ils archivent les instantanés `v1.0.2` /
  `v1.0.3` : une nouvelle édition doit être déclenchée pour que le DOI pointe
  sur le commit que le tag `v1.0` de la note de disponibilité nomme.
- **Découpage en volets.** `doc/paper1_core.md`, `doc/paper2_gauge_gravity.md`,
  `doc/paper3_phenomenology.md` existent et sont régénérés depuis le maître par
  `review/tools/build_series.py` ; la décision de soumettre un volet plutôt que
  le tout (et le choix de la revue) reste un choix d'auteur. Le volet le plus
  défendable seul est le no-go de cône (Prop. 6 + mur de rigidité + recherche de
  famille d'étages) : il se publie aussi en matière condensée.
- **Affiliation / endorsement arXiv** pour un auteur indépendant.
