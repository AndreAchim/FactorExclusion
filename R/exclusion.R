exclusion <- function(AS,donnees){
  # donnees est une liste de rangs de variables, cible suivie d'au moins deux variables d'exclusion
  # ou donnees est (ns,nv) avec toutes colonnes à une somme de carrés de 1.0
  #
  # D'abord vérifier si la cible corrèle avec les variables d'exclusion
  resultat <- list()
  donnees <- c(donnees[1],sort(donnees[-1]))
  ne <- length(donnees)-1
  Dat <- AS$GS[,donnees]
  # co <- AS$R[donnees[1],donnees[-1]]
  # #    co <- t(Dat[,1]) %*% Dat[,-1]
  # pr <- prodCorr(min(abs(co)),AS$N)
  # p <- pr$p
  # if (p>.50){
  #   # si la cible ne corrèle pas avec les variables d'exclusion
  #   for (k in 1:ne){
  #     resultat[[k]] <- list(po=0,crit=max(abs(co)),contrasts=Dat[,1],corr=co)
  #   }
  # } else {
  #   # si corrélations
  excl <- donnees[-1]
  satu <- matrix(NA,nrow=ne,ncol=ne)
  for (j in 1:(ne-1))
    for (i in (j+1):ne){
      sat <- asSatPaire(AS,excl[sort(c(i,j))])
      satu[i,j] <- sat[1]
      satu[j,i] <- sat[2]
    }
  sat <-  rowSums(satu,na.rm=TRUE)/(ne-1)
  li <- seq(-4,4,.25)
  cr <- li
  for (k in 1:length(cr))
    cr[k] <- crit_exclusion(li[k],Dat,sat)
  limites <- limites_balayage(li,cr)
  out <- optimize(crit_exclusion,limites,Dat,sat)
  p <- out$minimum
#  browser()
  for (k in 1:ne){
    po <- p*prod(sat[-k])
    contrast <- Dat[,1]-po*Dat[,k+1]
    cor <- t(sc1(contrast)) %*% Dat[,1+excl[-k]]
    crit <- max(abs(cor))
    resultat[[k]] <- list(crit=crit,po=po,corr=cor,contrast=contrast)
  }
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
