pair_cancel <- function(p,A_V,probe,fun=max_abs,as_list=FALSE){
  # p est un scalaire
  # A_V est (n,2) deux colonnes dont la somme des carrés est égale à 1.0
  # probe est (n,) autant de colonnes qu'on veut de corréations avec le contraste
  # fun est la fonction à appliquer aux carrés des projections pour établir le critère
  # Si as_list=TRUE, la sortie est une liste plutôt que juste le critère
  # calcule A_V %*% c(p,-1), lui donne une somme de carrés de 1.0 et y projette probe
  # met ces projections au carré et en retourne la fonction fun dans $crit
  # le contraste normalisé est dans $contrast
  # et les projections (pas au carré) dans $ proj
  prd <- A_V %*% c(p,-1)
  d <- -1 / sqrt(t(prd) %*% prd)  # pour ramener le contraste à l'équivalent de c(-p,1)
  contrast <- as.vector(d) * prd
  proj <- t(probe) %*% contrast
  crit <- fun(proj * proj)
  if (as_list){
    return(list(crit=crit,contrast=contrast,proj=proj))
  }
  else
    return(crit)
}