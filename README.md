# FactorExclusion

**FactorExclusion** est un package R implementant la methode d'exclusion de facteur (FE) 
dans le cadre de l'analyse par annulation de signal (Signal Cancellation Analysis, SCA).

## Fonctionnalites principales

- Estimation des saturations factorielles a partir d'une paire de variables ancres (`sat_sur_facteur`)
- Calcul des correlations residuelles et de leurs p-valeurs via la methode delta de Pearson-Filon (`prob_corr_residuelles`, `partitionne_R`)
- Simulation des taux de faux positifs pour differentes combinaisons de saturations et de taille d'echantillon (`FE_tests`, `FE_methodes2`)
- Infrastructure d'annulation de signal par paires (`init_SCA`, `as_paires_indicatrices`, `as_annule_paire`)

## Installation

```r
devtools::install_github("AndreAchim/FactorExclusion")
```

## Utilisation typique

```r
library(FactorExclusion)

# Simulation complete (environ 15 conditions)
efg <- FE_tests()

# Partitionner une matrice de correlations avec les ancres 1 et 2
AS  <- pre_partition(R, N = 175, ote = c(1, 2))
res <- partitionne_R(AS, ote = c(1, 2))
```

## Scripts de figures

Les scripts de reproduction des figures de l'article sont dans `inst/scripts/` :

- `fig_biais.R` -> Figure 1 (biais des correlations residuelles)
- `build_detection_data.R` -> construit `df_det` et ajuste les modeles GLM
- `fig_detection.R` -> Figure 2 (taux de detection)

```r
# Reproduction complete
efg <- FE_tests()
source(system.file("scripts/fig_biais.R",           package = "FactorExclusion"))
source(system.file("scripts/build_detection_data.R", package = "FactorExclusion"))
source(system.file("scripts/fig_detection.R",        package = "FactorExclusion"))
```

## Auteur

- Andre Achim (UQAM)

## Licence

GPL (>= 3)
