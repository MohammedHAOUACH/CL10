# Risques de rejet et reste à faire

Analyse en posture de referee, menée le 2026-10-21 sur `doc/final_publication.md`
(51 pages, Table 1, 10 open problems, 144 expériences). Chaque preuve citée a été
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
- **Lorentz.** $A\ge0.997$ sur 2726 points, mur de rigidité $4.8\times10^{-3}$,
  « cônes : 0 ». Aucune confrontation aux bornes publiées de violation de Lorentz.
  Paragraphe manquant → ajouté (§3 ci-dessous, action 3).
- **Gravité.** Champ symétrique ajouté à la main avec $G=1$ ; aucun graviton
  (spin 2 calculé négatif, §11.1, §10.4) ; quintet $E_g$ à $\lambda^2$ avec des
  vitesses $\pi/4,\pi/2\ne1$ ; tension de corde $5.4\times10^{40}\times$ QCD.
  Table 1 : « one object: partial ». Le titre et le résumé portent pourtant
  « gravitational sector » : sur-revendication perçue.
- **Matière noire.** Le porteur est un secteur singulet ajouté nommément ;
  $\Omega h^2=0.120$ n'est obtenu qu'avec $\lambda_B=4.06\times10^{-14}$ fitté
  (les valeurs fournies par le root donnent $0.798$, $0.0246$, $2.3\times10^{-5}$,
  $7.6\times10^{-4}$).
- **Falsifiabilité surtout négative.** $\tau_p=2.8\times10^{39}$ yr hors d'atteinte,
  la version testable ($1.8\times10^{36}$ yr) étant déclarée convention abandonnée ;
  le barreau nu $\nu$ est déjà exclu ; seul F3 mord à courte échéance.
- **Auditabilité.** 143 expériences auto-produites, tolérances définies par
  l'auteur, dépôt sans version épinglée (la note de disponibilité parle de
  « the archived release » sans tag, sans commit, sans DOI) ; et le résultat du
  coin (§8.2, Table 5, O9) — six états exactement sans masse à
  $w=2.8894,\ \theta=1.5717$ — provenait d'un calcul en direct sans artefact
 -commité, alors que la Table 9 l'attribue à la famille V118.

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
| 9 | Docstring de `V187_rigidity_frontier.py` et de `V186_rigidity_wall_scan.py` alignés sur leurs artefacts commités (le docstring citait 1.033 / $8.0\times10^{-3}$ / médiane $-0.28$, l'artefact 1.140 / $4.8\times10^{-3}$ / $-0.62$) | `experiments/dynamics_spectrum/V186…, V187…` | fait |

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

- **Dépôt Zenodo (DOI) et version épinglée publiquement.** Le tag `v1.0`
  est créé dans le dépôt ; le DOI Zenodo nécessite le compte auteur et ne peut
  être minté ici. La note de disponibilité reste donc au tag + au dépôt GitHub.
- **Découpage en volets.** `doc/paper1_core.md`, `doc/paper2_gauge_gravity.md`,
  `doc/paper3_phenomenology.md` existent et sont régénérés depuis le maître par
  `review/tools/build_series.py` ; la décision de soumettre un volet plutôt que
  le tout (et le choix de la revue) reste un choix d'auteur. Le volet le plus
  défendable seul est le no-go de cône (Prop. 6 + mur de rigidité + recherche de
  famille d'étages) : il se publie aussi en matière condensée.
- **Affiliation / endorsement arXiv** pour un auteur indépendant.
