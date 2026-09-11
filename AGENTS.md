# Memoire du projet — SignalCancellation

## Contexte general
Depot GitHub : https://github.com/AndreAchim/FactorExclusion
Projet d'analyse factorielle (Signal Cancellation Analysis — SCA). Donnees typiques : matrices de correlations avec N autour de 175. Les variables R et sessions courantes sont dans `tf$R` et `tf$N`.

---

## Structure des fichiers

- Fonctions personnelles dans `R/` (projet) et certaines dans `C:/Users/Achim_A/Documents/`
- `corr_partielles.R` est maintenant dans `R/` (chargé par `devtools::load_all()`). La ligne correspondante dans `.Rprofile` a été commentée (août 2026).
- `prob_Rr.R` reste dans `Documents/` et est sourcé par `.Rprofile`
- `rss_article.R` dans `R/` est redondant depuis l'installation de `rticles`

---

## Dependances des fonctions principales

### `R/FE_methodes2.R` requiert :
`gen_data`, `corr_partielles`, `prob_Rr`, `partitionne_R`, `sat_sur_facteur`, `crit_R2`, `SC1`
(corr_partielles et prob_Rr sont dans Documents/, non dans R/)

### `R/scand3.R` requiert :
`init_SCA`, `as_paires_indicatrices`, `as_annule_paire`, `probNonDoublet`, `optim_paire_initiale`, `prodCorr`, `crit_R2`, `prob_R2`, `SC1`, `seq_k_dim`*, `test_k1`, `asPaireAvec`, `sat_sur_facteur`, `fareg`, `rmsd`, `optim_tuple`
(* `seq_k_dim.R` contient aussi `ajustePolarites` et `polarites_de_correlations`)

**Fonctions communes aux deux scripts :** `crit_R2`, `prob_R2`, `sat_sur_facteur`, `SC1`

---

## Simulations de taux de faux positifs — session du 7-8 aout 2026

### Chaine d'appel pour les tests de correlations residuelles
```
FE_tests -> FE_methodes2 -> partitionne_R -> prob_corr_residuelles
```

### `R/FE_tests.R` — version finale (17 aout 2026)
Fonction de simulation principale. Parametres : `FE` (matrice de saturations), `initial`, `saut`, `N` (negatif = graine fixe), `seuil`.

Structure de sortie de chaque element `out[[i]]` :
- `$sat` : vecteur de 2 saturations des ancres
- `$ns` : N (negatif = graine fixe ; `abs(ns)` = taille d'echantillon)
- `$stt` : resultat de `sim_stats(fe)`
- `$sigCR` : matrice logique nOK x 66 — detection a alpha=.05 pour chaque simulation valide
- `$OK_idx` : indices (1 a nSim) des simulations non-Heywood — indispensable pour l'appariement

Code cle :
```r
OK_idx <- which(!fe$Heywood)
sigCR  <- fe$pCR[OK_idx, ] < .05
out[[length(out)+1]] <- list(sat=sat, ns=ns, stt=sim_stats(fe), sigCR=sigCR, OK_idx=OK_idx)
```

**Corrections apportees :**
- `FE$pCR` -> `fe$pCR` (FE est une matrice, pas une liste)
- Ajout de `sigCR` et `OK_idx` pour l'analyse paired clogit

---

### Probleme statistique central : Var(atanh(xy)) != 1/(N-3)

`xy` = produit des saturations estimees via les ancres A et B. La transformation atanh stabilise la variance uniquement pour une correlation de Pearson ; `xy` etant une fonction de 5 correlations, `Var(atanh(xy))` ne vaut pas 1/(N-3).

**Solution implantee (8 aout 2026) — `R/prob_corr_residuelles.R`** :
- Methode delta complete avec formule de Pearson-Filon (Steiger, 1980)
- Signature : `prob_corr_residuelles(R, AB, N)`
- Retourne NA pour les paires a signes incoherents ou |xy| >= 1

**`R/sim_stats.R`** — correction NA : `mean(pr[,k]<.05, na.rm=TRUE)`

---

## Figure de biais des correlations residuelles — sessions des 15-16 aout 2026

### Contexte
`efg` est une liste de 16 elements produite par `FE_tests()` :
- Elements 1-15 : 15 conditions = croisement de 5 paires de saturations x 3 tailles d'echantillon (`$ns` negatif = graine fixe)
- Element 16 : `$CRatt` (66 valeurs attendues de population) et `$zro` (indices des 20 paires nulles)

### Pipeline atanh/tanh dans `extraire_cond(i)`
Les lignes 4-6 de `$stt$CR` stockent des deviations par rapport a `atanh(CRatt)`, transformees via tanh. Recuperation en atanh absolu :

```r
moy_atanh  <- atanh(d4) + atanh(att0)                           # mean(atanh(dt)) absolu
se_atanh   <- (atanh(d6) - atanh(d5)) / (2 * qt(0.975, df_t))  # SE en atanh
demi_ci999 <- qt(0.9995, df_t) * se_atanh                       # IC 99.9%
lower_ci   <- tanh(moy_atanh - demi_ci999)
sd_atanh   <- se_atanh * sqrt(df_t + 1)                         # +/- 1 ET
lower_sd   <- tanh(moy_atanh - sd_atanh)
```

Les paires a valeur attendue negative sont inversees (`* sign(att0)`) pour travailler sur les valeurs absolues, puis reordonnees par `|att|` croissant.

### Script de reproduction
**`R/fig_biais.R`** — a sourcer avec `efg` dans l'environnement :
```r
CRatt_all <- efg[[16]]$CRatt
zro       <- efg[[16]]$zro
source("R/fig_biais.R")   # produit fig15.pdf et fig15.png
```

### Description de la figure finale (`fig15.pdf` / `fig15.png`)
- **Format** : portrait 17 x 22 cm ; exports PDF vectoriel + PNG 300 dpi
- **15 panneaux** : 5 rangees (paires de saturations) x 3 colonnes (N = 125, 500, 2000)
- **Couches visuelles (de bas en haut)** :
  1. `grey60` — bande +/- 1 ET des simulations (variabilite individuelle)
  2. Ligne noire mince (`linewidth = 0.5`) — valeurs attendues de population
  3. `grey85` — barres IC 99.9% verticales, uniquement aux paires ou `att` sort de l'IC (biais significatif a alpha = 0.001)
  4. Ligne blanche mince (`linewidth = 0.25`) — moyenne des simulations, par-dessus tout
- **Lecture** : ligne blanche sur ligne noire = pas de biais ; ecart avec barre grey85 = biais important

---

## Analyse des taux de detection — sessions des 16-17 aout 2026

### Objectif
Documenter l'effet des saturations des ancres sur la puissance de detection des correlations residuelles, au-dela des effets connus de la taille d'effet et de N.

### Data frame `df_det`
Construit a partir de `efg` (15 conditions x 46 paires non-nulles = 690 lignes) :
- `n_det` / `n_tot` : nombre de detections sur nOK simulations valides (`sgCR[1, nonzro] * nOK`)
- `att` : |correlation residuelle attendue| (variable continue)
- `N_val` : taille d'echantillon (`abs(ns)`)
- `sat` : facteur a 5 niveaux — paires de saturations des ancres

### Modeles GLM quasibinomial (m0q a m3q)
- `m0q` : `log(att) * log(N_val)` — sans saturation
- `m3q` : `log(att)*sat + log(N_val)*sat + log(att):log(N_val)` — modele complet
- Surdispersion x31.6 (les 46 paires partagent les memes tirages -> non-independance)
- Inference correcte via test F quasibinomial : **F(12, 674) = 119.76, p < .001**

### Regression logistique conditionnelle (clogit) — 17 aout 2026

**Donnees :** `efg3` — version de `FE_tests` sauvegardant `sigCR` et `OK_idx`

Construction de `long_df` :
```r
groupes   <- list(N125=c(1,4,7,10,13), N500=c(2,5,8,11,14), N2000=c(3,6,9,12,15))
OK_commun <- Reduce(intersect, lapply(i_group, function(i) efg3[[i]]$OK_idx))
# strate_id = run_id + (pair_id - 1) * nRuns  (unique par run x paire x N)
```

Runs communs apres intersection : 930 (N=125), 964 (N=500), 1000 (N=2000)
`long_df` : 665 620 lignes x 7 colonnes, 133 124 strates

Modeles :
- `cl0` : `log(att) + log(N_val) + strata(strate_id)` — nul effectif (log(att) et log(N_val) sont constants au sein des strates, donc absorbes)
- `cl1` : `sat + strata(strate_id)`
- `cl2` : `sat * log(att) + sat * log(N_val) + strata(strate_id)`

**Resultats :**

| Comparaison | chi2 | df | p |
|---|---|---|---|
| Effet principal de sat (cl0 -> cl1) | 41 081 | 3 | < .001 |
| Interactions sat x log(att) et sat x log(N) (cl1 -> cl2) | 14 135 | 8 | < .001 |
| **Omnibus complet (cl0 -> cl2)** | **55 216** | **11** | **< .001** |

A rapporter dans l'article : **chi2(11) = 55 216, p < .001**
Les t individuels du summary() ne sont PAS la base correcte (valeurs extrapolees hors plage).

### Script de reproduction
**`R/fig_detection.R`** — a sourcer avec `df_det` et `m3q` dans l'environnement :
```r
source("R/fig_detection.R")   # produit fig_detection.pdf et fig_detection.png
```
- 3 panneaux empiles (N = 125, 500, 2000), courbes GLM ajustees par niveau de saturation
- Niveaux de gris + types de lignes, legende en anglais ("Anchor loadings"), ordre inverse
- Format 12 x 18 cm, PDF vectoriel + PNG 300 dpi
- Note : courbes = GLM (m3q) ; test formel = clogit (cl0-cl2)

---

## Resultats sur la correlation entre facteurs (session du 14 juin 2026)

**Contexte :** 10 variables, N = 175, deux facteurs (intelligence verbale ancre sur v1,v5 ; performance ancre sur v7,v9). Saturations estimees via `sat_sur_facteur`.

| Methode | r | Note |
|---|---|---|
| Correlation canonique (v1,v5 x v7,v9) | 0.358 | Sous-estimation — 4 variables seulement |
| `psych::fa` oblimin (toutes variables) | **0.479** | Estimation la plus defendable |
| Methode O/P (4 paires d'ancres) | 0.532 | Legere surestimation — ancres tres saturees |
| Ajustement joint sur 45 paires | -0.346 | Artefact — ne pas interpreter |

**Conclusion :** r ~ 0.48 (oblimin). `sat_sur_facteur` produit des saturations brutes -> artefact de double comptage si utilise conjointement.

---

## Environnement R Markdown / LaTeX (août 2026)

- `rmarkdown` v2.31 et `rticles` v0.27 installes (binaires)
- TeX Live 2025 present (`pdflatex` disponible)
- Template RSS disponible via File → New File → R Markdown → From Template → "Royal Statistical Society Article"
- Fichier manuscrit actif : `ManuscritExclFact.Rmd`

### Fonction `prob_corr_residuelles` — note methodologique
- Teste H0 : r(X,Y) = xy via la statistique D = atanh(r_XY) - atanh(xy)
- Var(D) calculee par methode delta : **g**' Σ **g**
  - Σ = covariance des 6 correlations brutes (formule Pearson-Filon, Steiger 1980) — sans facteur (1-r²)
  - **g** = gradient de D par rapport aux 6 correlations — contient les facteurs (1-r²) via derivees de atanh
- Signe positif du terme rAB dans **g** : deux negatifs s'annulent (signe moins dans D, derivee negative de xy/rAB)
- Distribution asymptotique normale (methode delta = approximation grands echantillons)

---

## Figures MQSHS2026

- `fig1` : nuage de points de 7 variables (A-G) dans un espace a deux facteurs, defini dans `R/MQSHS2026.R`
- `fig2` : `fig1` + segment plein de l'origine a la variable F (`geom_segment`)
- Sauvegarde avec `ggsave("MQSHS_2.tiff", plot = fig2, width = 5, height = 5)`

---

## Manuscrit ArticleRSS2026 — etat au 2 septembre 2026

### Fichier : `ArticleRSS2026/ArticleRSS2026.Rmd`
- Template `rticles::rss_article` (xelatex), classe `statsoc.cls`
- Double interligne via `\usepackage{setspace}` + `\doublespacing` dans `header-includes`
- Paragraphes en prose reformates a 80 caracteres (navigation editeur fonctionnelle)

### Structure de soumission (flottants apres les references)
Ordre : texte → `\clearpage` → References → `\clearpage` → Tableau 1 → `\clearpage` → Figure 1 → Figure 2 → `\clearpage` → Listing 1
- **Tableau 1** : chunk R avec `results='asis'`, `centering=FALSE` dans `kable()`, reference `\ref{tab:structure}` dans le texte
- **Figures** : blocs `\begin{figure}[H]` sans `\includegraphics` (pages de legendes + `\textit{Alt text:}` pour accessibilite). La ligne `\includegraphics` est en commentaire `%` pour reinsertion facile. Fichiers renommes `Figure1.pdf` et `Figure2.pdf`
- **Listing 1** : bloc `{=latex}` avec `\begin{lstlisting}` sans option `float`

### Contraintes statsoc.cls
- L'environnement `table` emballe le contenu dans `\hbox` (mode LR) : ne jamais mettre `\centering` ni `kable_styling()` dedans
- `endfloat` incompatible avec `\statsocwidth@` — ne pas utiliser
- `\floatplacement` et `\renewcommand{\topfraction}` retires

### Figures de simulation
- **`R/fig_biais.R`** → `Figure1.pdf` / `Figure1.png` (17x22 cm)
  - Ordre global des 46 paires : `|att|` croissant, puis `grand_moy` (moy. observee sur 15 conditions) en secondaire → variable `global_o2` precalculee avant `extraire_cond()`
- **`R/fig_detection.R`** → `Figure2.pdf` / `Figure2.png` (12x18 cm)
  - Axe x = `att` continu ; doublons = observations supplementaires dans le GLM, pas un probleme
- **`R/build_detection_data.R`** — construit `df_det` (690 lignes) et ajuste `m0q`, `m3q`
  - Utilise `colSums(cond$sigCR[, nonzro])` si `$sigCR` disponible, sinon `round(sgCR[1,] * nOK)`
  - Lit `CRatt_all` et `zro` depuis `efg[[16]]` au debut du script
- Les fichiers sont sauvegardes a la racine du projet

### Flux de reproduction complet
```r
AA <- FE_tests()                    # assign() met toujours efg global a jour
source("R/fig_biais.R")             # Figure1.pdf
source("R/build_detection_data.R")  # df_det, m0q, m3q
source("R/fig_detection.R")         # Figure2.pdf
```
Les trois `source()` sont dans `FE_tests.R` (commentes) — decommenter pour tout produire d'un coup.

### FE_tests.R — modifications septembre 2026
- Corrige `zro <- which(abs(CRatt)<1e-6)` en `abs(CRatt_all)`
- Ajoute `assign("efg", out, envir = .GlobalEnv)` avant les figures
- Le nom de la variable de retour n'importe pas (`AA <- FE_tests()` fonctionne)

### Fichiers de soumission separee (septembre 2026)

#### `ArticleRSS2026/ArticleRSS2026.Rmd` — version RSS Series B
- Moteur : `pdflatex` (XeLaTeX abandonne : `\mathrm` et `\textsc` invisibles en subscript)
- `header-includes` : `inputenc` (utf8), `fontenc` (T1), `setspace`/`doublespacing`, `\geometry{a4paper}`, `\renewcommand{\maketitle}{}`
- `includes: in_header: fix_abstract.tex` (dans la section `output`)
- **`fix_abstract.tex`** : definit `\rc` et `\Rc` + patche l'environnement `abstract` + redefinit `\keywords` pour sortie immediate

#### `ArticleRSS2026/FactorExclusion.Rmd` — version arXiv
- Moteur : `pdflatex`
- Memes encodages ; titre et abstract intacts (pas de `\renewcommand{\maketitle}{}`)
- `includes: in_header: fix_commands.tex`
- **`fix_commands.tex`** : definit seulement `\rc` et `\Rc`

#### Commandes `\rc` et `\Rc`
```latex
\newcommand{\rc}[1]{r_{\scriptscriptstyle #1}}
\newcommand{\Rc}[1]{R_{\scriptscriptstyle #1}}
```
- `\scriptscriptstyle` reduit les indices au niveau sub-subscript (~50% taille normale) -- c'est ce qui donne l'aspect petites capitales
- Sans `\mathrm`, les lettres sont en italique (defaut du mode math) : effet petites capitales italiques
- Definies dans des fichiers `.tex` passes via `--include-in-header` (verbatim) — **ne pas mettre dans `header-includes` YAML** : Pandoc traite les `\newcommand` du YAML comme des macros internes, causant une boucle infinie (PandocMacroLoop, code 91)
- Remplacement global `r_{` → `\rc{` : attention, remplace aussi dans les definitions elles-memes si elles sont dans le YAML

#### `RSS_table1/RSS_table1.Rmd` — tableau seul
- `output: pdf_document` avec `latex_engine: xelatex`, `header-includes: booktabs, threeparttable`
- `kable_styling(latex_options = "HOLD_position")` pour eviter que le tableau flotte hors page

#### Notes statsoc.cls complementaires
- `statsoc.cls` definit `\maketitle` ET appelle `\AtBeginDocument{\maketitle}` ; redefinir avec `\renewcommand{\maketitle}{}` dans le **preambule** (pas dans `\AtBeginDocument`) pour supprimer la page de titre
- `geometry` est deja charge par `statsoc.cls` : utiliser `\geometry{a4paper}` (pas `\usepackage[a4paper]{geometry}`)
- `\keywords{}` est place **apres** `\end{abstract}` par le template rticles : la redefinition dans `fix_abstract.tex` le fait sortir immediatement

### MiKTeX — conflit rpostback.exe
- `.Rprofile` place `C:/texlive/2025/bin/windows` en tete du PATH
- Evite l'erreur MiKTeX qui traite `rpostback.exe` comme un repertoire

### Dependances
- `expm` ajoute au `DESCRIPTION` (`Imports`) — requis par `sqrtm` dans `gen_data.R`, `prod_sim.R`, `prep_FP.R`, etc.
- Les `library()` dans les fichiers `R/` sont a remplacer par des appels qualifies (`expm::sqrtm`, `MASS::mvrnorm`) lors du prochain nettoyage du package
