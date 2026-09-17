exclusion <- function(AS,donnees,nom){
  # donnees est une liste de rangs de variables, une cible suivie d'au moins deux
  # variables d'exclusion
  nexcl <- length(AS$exclusions)
  noms <- dimnames(AS$GS)[[2]]
  donnees <- c(donnees[1],sort(donnees[-1])) # covariables par ordre croissant
  nd <- length(donnees)
  ne <- nd-1
  Dat <- AS$GS[,donnees]
  li <- seq(-4,4,.25)
  cr <- li
  rg <- 2:nd
  nm1 <- noms[donnees[1]]
  for (k in 1:ne){  # pour chaque covariable
    Dt <- Dat[,c(1,rg[k],rg[-k])] # la mettre en 2e position suivie des témoins
    nm <- paste0(nm1,noms[donnees[k+1]])
    for (j in 1:length(cr))
      cr[j] <- crit_exclusion(li[j],Dt)
    limites <- limites_balayage(li,cr)
    out <- optimize(crit_exclusion,limites,Dt)
    po <- out$minimum
    contrast <- as.matrix(Dt[,1]-po*Dt[,2])
    colnames(contrast) <- nm
    AS$GS <- cbind(AS$GS,contrast)
    cor <- t(sc1(contrast)) %*% sc1(Dt[,3:ncol(Dt)])
    crit <- max(abs(cor))
    nexcl <- nexcl+1
    AS$exclusions[[nexcl]] <- list(crit=crit,poids=po,correl=cor+0) # +0 pour ne garder que la valeur, pas le titre
    names(AS$exclusions)[[nexcl]] <- nm
  }
  return(AS)
}



crit_exclusion <- function(p,Dat){
  # les colonnes de Dat  sont cible, exclusion, témoins
  contras <- sc1(Dat[,1]-p*Dat[,2])
  mxabs <- max(abs(t(contras) %*% Dat[,3:ncol(Dat)]))
  return(mxabs)
}
