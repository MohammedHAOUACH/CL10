# Risques de rejet et reste à faire

Analyse en posture de referee, menée le 2026-10-21 sur `doc/final_publication.md`
(73 pages, Table 1, dix problèmes ouverts, 162 expériences). Chaque preuve citée a
été relue dans le fichier ou dans l'artefact correspondant.

## 1. Risques rédhibitoires (sélectivité éditoriale)

| Rang | Faiblesse | Preuve dans le papier | Ce qu'écrit le referee | Correction |
|---|---|---|---|---|
| A | Aucun spectre bas-énergie réaliste : aucun fermion chiral léger, dispersif et isotrope n'est produit | §8.3 : à $\pi/4$ l'espace nul de $\Gamma$ est **plat** ($\omega/|\mathbf q|\to0$) et à $1/5$ $\Gamma$ est gappé de $3.65\times10^{-2}$ ; §9 : touchement **quadratique** sur les six droites, exposants $(1,1,2)$ ; §8.2/8.3 : à $1/5$ les seuls modes exactement sans masse ont un little group d'ordre 2 (pas de paire d'hélicité) et une anisotropie $0.857$ ; Table 4 : le secteur électrofaible est une **entrée mesurée** ; Table 5 : aucune paire $\pm1$/$\pm2$ nulle part | « Where is the electron ? Le dédoublement n'est pas supprimé, il est déplacé. Les nombres quantiques internes sont exacts, la réalisation spatio-temporelle des fermions observés n'est jamais établie. » | **O10 nommé** « the low-energy spectrum » + phrase dans le résumé |
| B | L'axiome de fermeture est en tension avec la seule mesure du coefficient qui départage les horloges | §12.1/O8 : l'empreinte fonctionnelle lit $g_{\rm obs}=0.199998$ — l'horloge non fermée $1/5$ à $10^{-5}$ — alors que la classe fermée lit $2.7002$ ($\pi/4$) et $2.7084$ ($3\pi/4$), à $30.0\%$ et $30.4\%$ de la valeur observée ; la lecture $\pi/2$ ($2.2994$, $10.7\%$) est un artefact de plancher sans limite à l'angle de fermeture ($n^*=102$ est une longueur $4k{+}2$, V209) et n'est pas un nombre ; F2 oppose $0.989868$ (branche fermée) à $0.96388$ (hiérarchie) ; « the closure axiom and the loop-hierarchy identification cannot both hold » | « L'axiome central est réfuté par la mesure du papier lui-même. » | Reformuler : l'axiome est une affirmation **structurelle** sur la dynamique sur-site, pas une hypothèse empirique sur l'angle choisi ; F2 est la mesure qui tranche |
| E | Aucune discussion de Bell / de contextuality, alors que le papier revendique un déterminisme exact de la mise à jour | le mot « Bell » (sensible à la casse), « Kochen », « contextual » sont absents ; la projection et la règle de Born sont déclarées entrées (§1.3, §6, §13.1) | « Un réseau unitaire localement déterministe qui prétend fonder la mécanique quantique doit dire d'où viennent les corrélations violant Bell. » | Paragraphe de périmètre dans §6 |
| K | Format : 73 pages, auteur unique sans affiliation, dix problèmes ouverts, ~40 résidus nommés, Table 1 concluant « parameters: not met », vocabulaire idiosyncratique | Table 1 ; Table 6 (11 coefficients fittés, 2 conventions, 1 conditionnel) | « Un programme de recherche avec son tableau d'échecs, pas un résultat. » | Découpage en volets (voir §4) |

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
  « cônes : 0 ». Le papier affronte les bornes publiées de violation de Lorentz
  (§8.4, « the experimental bound that would bite ») et montre qu'aucune borne
  actuelle ne mord sur le secteur sans masse.
- **Gravité.** Le contenu est calculé, pas seulement ajouté : §11.5
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
- **Auditabilité.** 162 expériences auto-produites, tolérances définies par
  l'auteur. Le résultat du coin (§8.2, Table 5, O9) a son artefact commité,
  `V207`, avec le scan complet du plan $(w,\theta)$, et chaque nombre du papier
  est adossé à une feuille d'artefact. Le grief qui reste ouvert : les
  tolérances de chaque verdict sont choisies par l'auteur, sans re-run
  indépendant fourni.

## 3. Ce qui est fort et doit rester en tête d'affiche

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

## 4. Ce qui reste hors du dépôt (décisions auteur)

- **DOI et version épinglée publiquement.** Le DOI Zenodo est minté
  (concept 10.5281/zenodo.23321106, version 10.5281/zenodo.23321228) et la
  release GitHub existe, mais elles archivent un instantané antérieur au
  commit que la note de disponibilité nomme : une nouvelle édition doit être
  déclenchée pour que le DOI pointe dessus.
- **Découpage en volets.** `doc/paper1_core.md`, `doc/paper2_gauge_gravity.md`,
  `doc/paper3_phenomenology.md` existent et sont régénérés depuis le maître par
  `review/tools/build_series.py` ; la décision de soumettre un volet plutôt que
  le tout (et le choix de la revue) reste un choix d'auteur. Le volet le plus
  défendable seul est le no-go de cône (Prop. 6 + mur de rigidité + recherche de
  famille d'étages) : il se publie aussi en matière condensée.
- **Affiliation / endorsement arXiv** pour un auteur indépendant.
