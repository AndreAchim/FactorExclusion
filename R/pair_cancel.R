pair_cancel <- function(p,A_V,probe,as_list=FALSE){
  # p est un scalaire
  # A_V est (n,2) deux colonnes dont les sommes des carrés sont égale à 1.0
  # probe est (n,) autant de colonnes qu'on veut de corréations avec le contraste
  # Si as_list=TRUE, la sortie est une liste plutôt que juste le critère
  # calcule le contraste A_V %*% c(1,-p) et projette probe sur le contraste normalisé
  # produisant des corrélations
  # retourne dans $crit la plus grande absolue de ces corrélations
  # le contraste non-normalisé est retourné dans $contrast
  # les projections sur le contraste non-normalisé dans $proj
  # et les corrélations dans $corr
  contrast <- A_V %*% c(1,-p)
  if (!is.numeric(probe)) browser()
  corr <- t(probe) %*% sc1(contrast)
  crit <- max(abs(corr))
  if (as_list){
    browser()
    return(list(crit=crit,contrast=contrast,corr=corr,proj=t(probe) %*% contrast,contrast=contrast))
  }else
    return(crit)
}