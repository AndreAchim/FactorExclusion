crit_annule_paire <- function (p,A_V,probe,signe=NULL,as_list=FALSE){
# p est un scalaire sur une échelle logarithmique
# signe est celui de la corrélation entre les deux colonnes de A_V  
# A_V est (n,2) deux colonnes dont la somme des carrés est égale à 1.0
# probe est (n,) autant de colonnes qu'on veut de corréations avec le contraste
# fun est la fonction à appliquer aux carrés des projections pour établir le critère
# Si as_list=TRUE, la sortie est une liste plutôt que juste le critère
# calcule A_V %*% c(p,-1), lui donne une somme de carrés de 1.0 et y projette probe
# met ces projections au carré et en retourne la fonction fun dans $crit
# le contraste normalisé est dans $contrast
# et les projections (pas au carré) dans $ proj
  if (is.null(signe))
    signe <- sign(A_V[,1] %*% A_V[,2])
  # prd <- A_V %*% c(signe*exp(p),-1)
  # d <- -1 / sqrt(t(prd) %*% prd)  # pour ramener le contraste à l'équivalent de c(-p,1)
  # contrast <- as.vector(d) * prd
  contrast <- sc1(A_V %*% c(signe*exp(p),-1))
  proj <- t(probe) %*% contrast
  crit <- max_abs(proj)
#  crit <- t(proj) %*% proj
  if (as_list){
    return(list(crit=crit,contrast=contrast,proj=proj))
  }
  else
    return(crit)
}
