dist_z2 <- function(F,N=75,nrep=200,pF=NULL){
# inspection des distributions de z2 et en particulier du critère d'optimisation
# le z2 maximum
# nrep <- 200
# N <- 75
# F <- matrix(c(.6,.4,.8,.4,.7,0,0, 0,0,0,.7,.4,.5,.8),ncol=2)
# Quand F inclue des facteurs corrélés, pF peut être la matrice de patrons pour
# déterminer les paires exclusives au même facteur   
  nv <- nrow(F)
  nf <- ncol(F)
  if (is.null(pF)) pF <- F
  paires <- toutes_paires(pF)
  np <- nrow(paires)
  di <- c(nv-2,np,nrep)
  proj <- array(rep(NaN,prod(di)),di)
  z2 <- matrix(NaN,nrow=np,ncol=nrep)
  sqNm1 <- sqrt(N-1)
  browser()
  for (k in 1:nrep){
    dat <- scale(gen_data(F,N))/sqNm1
    rdat <- chol(t(dat) %*% dat)
    for (j in 1:np){
      a <- abs(paires[j,1])
      b <- paires[j,2]
      out <- annule_paire(rdat,a,b,N)
      if (length(out)>2){
# out <- list(crit=out$objective,po=exp(out$minimum),contrast=contraste,proj=t(contraste) %*% temoins)
        ctrst <- sc1(out$po * dat[,a] - dat[,b])
        mp <- mp+(mp>a)+(mp>b)
        prod <- ctrst * dat[,mp]
        tprod <- sign(prod) * sqrt(abs(prod))
        tt1 <- t.test(prod)
        tt2 <- t.test(tprod)
        
        
        proj[,j,k] <- out$proj * sqNm1
        z2[j,k] <- out$z2
      } else
        z2[j,k] <- out$po  # Inf pour r(ns) ou NaN pour minimum pas trouvé
#      aa <- 0
    }
  }
  out <- list(proj=proj,z2=z2,paires=paires)
  if (length(out)<1) browser()
  rapport_z2(out)
  return(out)
}