exclusion <- function(AS,donnees){
  # donnees est une liste de rangs de variables, une cible suivie d'au moins deux variables d'exclusion
  # ou donnees est (ns,nv) avec toutes colonnes à une somme de carrés de 1.0
  #
  # D'abord vérifier si la cible corrèle avec les variables d'exclusion
  resultat <- list()
  donnees <- c(donnees[1],sort(donnees[-1]))
  nd <- length(donnees)
  Dat <- AS$GS[,donnees]
#  excl <- donnees[-1]
  ne <- nd-1
  li <- seq(-4,4,.25)
  cr <- li
  rg <- 2:nd
  for (k in 1:ne){
    Dt <- Dat[,c(1,rg[k],rg[-k])]
    for (j in 1:length(cr))
      cr[j] <- crit_exclusion(li[j],Dt)
    limites <- limites_balayage(li,cr)
    out <- optimize(crit_exclusion,limites,Dt)
    po <- out$minimum
    contrast <- Dt[,1]-po*Dt[,2]   ## k peut-il ici désigner plus d'une covariable?
    cor <- t(sc1(contrast)) %*% Dt[,3:ncol(Dt)]
    crit <- max(abs(cor))
    resultat[[k]] <- list(crit=crit,po=po,corr=cor,contrast=contrast)
  }
#  if (donnees[1]==7 && crit>.2) browser()
  return(resultat)
}



crit_exclusion <- function(p,Dat){
  # les colonnes de Dat  sont cible, exclusion, témoins
    contras <- sc1(Dat[,1]-p*Dat[,2])
    mxabs <- max(abs(t(contras) %*% Dat[,3:ncol(Dat)]))
  return(mxabs)
}
