optim_tuple <- function(AS,melange){
  # melange contient les rangs de colonnes du tuple.
  TEMOINS <- as.matrix(AS$GS[,-melange]) # même si une seule colonne
  # if (length(AS$orphelines>0))    # on veut plutôt ne pas enlever les variables orphelines
  #   TEMOINS <- TEMOINS[,-AS$orphelines]
  Rinv <- solve(t(TEMOINS) %*% TEMOINS)
  VAR <- AS$GS[,melange]
  nm <- length(melange)
  out <- list()
  rg <- 1:nm
  poids <- matrix(NA,nrow=nm,ncol=nm)
  crit <- rep(NA,nm)
  basCrit <- 9e9
  meilPoids <- rep(NA,nm)
  for (k in 1:nm){
    par <- rep(0,nm-1)
    rg <- c(rg[-1],rg[1])
    VAR <- cbind(VAR[,-1],VAR[,1])  # mettre la variable k en dernier
    out <- optim(par,crit_R2,gr=NULL,VAR,TEMOINS,Rinv,method="Nelder-Mead")
    if (out$value<basCrit){
      basCrit <- out$value
      meilPoids[rg] <- c(out$par,-1)
    }
    crit[k] <- out$value
    cor <- abs(t(sc1(VAR %*% c(out$par,-1))) %*%  TEMOINS)
    po <- as.vector(c(out$par,-1))
    po <- po/as.numeric(sqrt(po %*% po))
    po[rg] <- po
    poids[,k] <- po
  }
  # mc <- which(crit==basCrit)[1]
  pR2 <- prob_R2(basCrit,ncol(TEMOINS),AS$N)
  meilPoids <- -meilPoids/abs(meilPoids[1])
  return(list(R2=basCrit,prob=pR2,meilPoids=meilPoids,poids=poids,crit=crit))
}

crit_tuple <- function(p,VAR,TEMOINS){
  max(abs(t(sc1(VAR %*% c(p,-1))) %*%  TEMOINS))
}