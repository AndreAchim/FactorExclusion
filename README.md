# FactorExclusion

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22319448.svg)](https://doi.org/10.5281/zenodo.22319448)

**FactorExclusion** est un package R implémentant la méthode d'exclusion de facteur (FE) 
dans le cadre de l'analyse par annulation de signal (Signal Cancellation Analysis, SCA).

## Fonctionnalités principales

- Estimation des saturations factorielles à partir d'une paire de variables ancres (`sat_sur_facteur`)
- Calcul des corrélations résiduelles et de leurs p-valeurs via la méthode delta de Pearson-Filon (`prob_corr_residuelles`, `partitionne_R`)
- Simulation des taux de faux positifs pour différentes combinaisons de saturations et de taille d'échantillon (`FE_tests`, `FE_methodes2`)
- Infrastructure d'annulation de signal par paires (`init_SCA`, `as_paires_indicatrices`, `as_annule_paire`)

## Installation

```r
devtools::install_github("AndreAchim/FactorExclusion")
```

## Utilisation typique

```r
library(FactorExclusion)

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
# Reproduction complète
efg <- FE_tests()
source(system.file("scripts/fig_biais.R",           package = "FactorExclusion"))
source(system.file("scripts/build_detection_data.R", package = "FactorExclusion"))
source(system.file("scripts/fig_detection.R",        package = "FactorExclusion"))
```

## Citation

```bibtex
@software{achim_2026_factorexclusion,
  author  = {Achim, André},
  title   = {FactorExclusion},
  year    = {2026},
  doi     = {10.5281/zenodo.22319448},
  url     = {https://github.com/AndreAchim/FactorExclusion}
}
```

## Auteur

- André Achim (UQAM)

## Licence

GPL (>= 3)
