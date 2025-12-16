rapport_exclusion_rp <- function(SF,noms,pa=1){
  #  F <- SF$patron %*% SF$corrFct [,-1]
  # q05 <- rep(NA,6)
  # q01 <- rep(NA,6)
  Co <- rep(NA,1000)  # moy des corrélations
  Fi <- rep(NA,1000)  # moy des corrélations transformées
  Pr <- rep(NA,1000)  # moy des probabilités
  Lp <- rep(NA,1000)  # moy des probabilités transformées
  #  for (pa in 1:length(SF$paires)){
  Li <- SF$paires[[pa]]
  np <- ncol(Li)
  for (nom in noms) {
    print(nom)
    mt <- get(nom)[[pa]]
    N <- as.numeric(substring(nom,2))
    for (k in 1:np){
      # Co <- rep(NA,1000)  # moy des corrélations
      # Fi <- rep(NA,1000)  # moy des corrélations transformées
      # Pr <- rep(NA,1000)  # moy des probabilités
      # Lp <- rep(NA,1000)  # moy des probabilités transformées
      #      if (identical(nom,"L240") && k==8) browser()
      #      C <- mt[k,7:12,]
      #      sd <- apply(abs(C),2,signes_diff) ##### temporaire
      #      browser()
      #      nsd <- sum(sd)
      #      if (nsd>500) browser()
      # Co[sd] <- .5
      # Fi[sd] <- .5
      # Pr[sd] <- .5
      # Lp[sd] <- .5
      R <- colMeans(mt[k,7:12,])
      R <- R*sqrt((N-2)/(1-R*R))
      Co <- 2*pt(-abs(R),N-2)
      #      Co[!sd] <- pt(R,N-2)
      Z <- colMeans(atanh(mt[k,7:12,]))*sqrt(N-3)
      Fi <- 2*pnorm(-abs(Z))
      Pr <- colMeans(mt[k,13:18,])
      Lp <- exp(colMeans(log(mt[k,13:18,])))
      Zp <- 2*pnorm(colMeans(qnorm(mt[k,13:18,]/2)))
# 
#       sd <- apply(abs(mt[k,1:6,]),2,signes_diff) ##### temporaire
#       Co[sd] <- .5
#       Fi[sd] <- .5
#       Pr[sd] <- .5
#       Lp[sd] <- .5
#       
      r05 <- sum(Co<.05)
      z05 <- sum(Fi<.05)
      l05 <- sum(Lp<.05)
      p05 <- sum(Pr<.05)
      n05 <- sum(Zp<.05)
      r01 <- sum(Co<.01)
      z01 <- sum(Fi<.01)
      l01 <- sum(Lp<.01)
      p01 <- sum(Pr<.01)
      n01 <- sum(Zp<.01)
      C <- mt[k,1:6,]
      mCR <- rowMeans(C)
      C <- mt[k,7:12,]
      mr <- rowMeans(C)
      out <- sprintf('%d %d,CR %4.3f r %4.3f N<.05:%3d %3d %3d %3d %3d',Li[1,k],Li[2,k],mean(mCR),mean(mr),r05,z05,p05,l05,n05)
      if (pa==1)
        out <- paste(out,sprintf(' N<.01:%3d %3d %3d %3d %3d',r01,z01,p01,l01,n01))
      print(out,quote=FALSE)
    }
  }
}

signes_diff <- function(C) {
  any(sign(C[1])!=sign(C[-1]))
}