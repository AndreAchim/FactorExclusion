# FactorExclusion

**FactorExclusion** est un package R implémentant la méthode d'exclusion de facteur (FE) 
dans le cadre de l'analyse par annulation de signal (Signal Cancellation Analysis, SCA).

## Fonctionnalités principales

- Estimation des saturations factorielles à partir d'une paire de variables ancres (`sat_sur_facteur`)
- Calcul des corrélations résiduelles et de leurs p-valeurs via la méthode delta de Pearson-Filon (`prob_corr_residuelles`, `partitionne_R`)
- Simulation des taux de faux positifs pour différentes combinaisons de saturations et de taille d'échantillon (`FE_tests`, `FE_methodes2`)
- Infrastructure d'annulation de signal par paires (`init_SCA`, `as_paires_indicatrices`, `as_annule_paire`)

## Installation

```r
# Via devtools (une fois le dépôt GitHub créé) :
devtools::install_github("AndreAchim/FactorExclusion")
```

## Utilisation typique

```r
library(FactorEclusion)

# Simulation complète (environ 15 conditions)
efg <- FE_tests()

# Partitionner une matrice de corrélations avec les ancres 1 et 2
AS  <- pre_partition(R, N = 175, ote = c(1, 2))
res <- partitionne_R(AS, ote = c(1, 2))
```

## Scripts de figures

Les scripts de reproduction des figures de l'article sont dans `inst/scripts/` :

- `fig_biais.R` → Figure 1 (biais des corrélations résiduelles)
- `build_detection_data.R` → construit `df_det` et ajuste les modèles GLM
- `fig_detection.R` → Figure 2 (taux de détection)

```r
# Exemple de reproduction complète
efg <- FE_tests()                                  # génère efg global
source(system.file("scripts/fig_biais.R",          package = "FactorExclusion"))
source(system.file("scripts/build_detection_data.R", package = "FactorExclusion"))
source(system.file("scripts/fig_detection.R",      package = "FactorExclusion"))
```

## Auteurs

- André Achim (UQAM)
- P.-O. Caron (TÉLUQ)

## Licence

CC BY 4.0
