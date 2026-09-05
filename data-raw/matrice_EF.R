## Script de création de l'objet matrice_EF
## À exécuter depuis la racine du package FactorEclusion

matrice_EF <- matrix(
  c(.5, .7, .4, .5, .6, .3,  0, -.6, -.4, .5, .4,
     0,  0,  0, .6,  0, .4, .4,  .5,   0, .6, .2,
     0,  0,  0,  0, .4, .4, .6,   0,  .6,-.4, .5),
  ncol = 3
)

usethis::use_data(matrice_EF, overwrite = TRUE)
