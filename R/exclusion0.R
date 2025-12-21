exclusion0 <- function(AS,donnees){
  # donnees est une liste de rangs de variables, une cible suivie d'au moins deux variables d'exclusion
  # ou donnees est (ns,nv) avec toutes colonnes à une somme de carrés de 1.0
  #
  # D'abord vérifier si la cible corrèle avec les variables d'exclusion
  resultat <- list()
  donnees <- c(donnees[1],sort(donnees[-1]))
  ne <- length(donnees)-1
  Dat <- AS$GS[,donnees]
  excl <- donnees[-1]
  sat <-  AS$satur_var
  li <- seq(-4,4,.25)
  cr <- li
  for (k in 1:length(cr))
    cr[k] <- crit_exclusion(li[k],Dat,sat)
  limites <- limites_balayage(li,cr)
  out <- optimize(crit_exclusion,limites,Dat,sat)
  p <- out$minimum
  for (k in 1:ne){
    po <- p*prod(sat[-k])
    contrast <- Dat[,1]-po*Dat[,k+1]   ## k peut-il ici désigner plus d'une covariable?
    cor <- t(sc1(contrast)) %*% Dat[,-c(1,k+1)]
    crit <- max(abs(cor))
    resultat[[k]] <- list(crit=crit,po=po,corr=cor,contrast=contrast)
  }
  if (donnees[1]==7 && crit>.2) browser()
  return(resultat)
}



crit_exclusion <- function(p,Dat,satur){
  #Dat(ns,length(satur+1)); satur:saturations estimées des variables d'exclusion
  ne <- length(satur)
  mxabs <- 0
  Dt <- Dat[,-1]
  for (k in 1:ne){
    contras <- sc1(Dat[,1]-p*prod(satur[-k])*Dt[,k])
    mxabs <- max(mxabs,max(abs(t(contras) %*% Dt[,-k])))
  }
  return(mxabs)
}
