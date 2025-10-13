crit_annule_paire <- function (p,A_V,probe,signe=0,as_list=FALSE){
# p est un scalaire sur une échelle logarithmique lorsque abs(signe) est 1 
# signe est alors celui de la corrélation entre les deux colonnes de A_V
# Signe est 0 quand la corrélation n'est pas nettement positive ou négative 
# A_V est (n,2) deux colonnes dont la somme des carrés est égale à 1.0
# probe est (n,) autant de colonnes qu'on veut de corréations avec le contraste
# Si as_list=TRUE, la sortie est une liste plutôt que juste le critère
# calcule A_V %*% c(1,-p), lui donne une somme de carrés de 1.0 et y projette probe
# et en retourne le max(abs()) dans crit le critère à minimiser
# le contraste non-normalisé est dans $contrast
# les corrélations et les projections dans $correl et $proj
  if (is.null(signe)) error('le paramères signe st maintenant obligatoire dans crit_annule_paire')
  if (signe==0)
    contrast <- A_V %*% c(1,-p)    
  else
    contrast <- A_V %*% c(1,-signe*exp(p))
  corr <- t(probe) %*% sc1(contrast)
  crit <- max_abs(corr)
  if (as_list){
    return(list(crit=crit,contrast=contrast,proj=t(probe) %*% contrast,correl=corr))
  }
  else
    return(crit)
}
