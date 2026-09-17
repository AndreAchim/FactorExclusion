optim_paire_initiale <- function(AS,ab){
  # le contraste est -a+p*b (le négatif de a-pb)
  autres <- setdiff(AS$pertinent,ab)
  col <- AS$GS[,ab]
  TEMOINS <- AS$GS[,autres]
  Rinv <- solve(t(TEMOINS) %*% TEMOINS)
  nt <- length(autres)
  # na <- ab[1]
  # nb <- ab[2]
  aa <- rep(NA,nt)
  bb <- rep(NA,nt)
  rAB <- AS$R[ab[1],ab[2]]
  for (k in 1:nt){
    if( rAB*prod(AS$R[ab,k]) <= 0){
      aa[k] <- 1
      bb[k] <- 0
    } else {
      aa[k] <- rAB * AS$R[ab[1],k] / AS$R[ab[2],k]
      bb[k] <- rAB * AS$R[ab[2],k] / AS$R[ab[1],k]
      if (is.nan(aa[k]) || is.nan(bb[k])) browser()
    }
  }
  aa <- sign(rAB) * sqrt(aa)
  bb <- sqrt(bb)
  pp <- median(bb / aa)
  if (is.nan(pp)) browser()
  # browser()
  out <- optim(pp,crit_R2,gr=NULL,col,TEMOINS,Rinv,method="BFGS")
  R2 <- out$value
  prob <- prob_R2(R2,nt,AS$N)
  po <- out$par
  contrast <- col %*% c(-1,po)
  corr <- t(sc1(contrast)) %*%  TEMOINS
  return(list(crit=R2,po=po,prob=prob,corr=corr,contrast=contrast))
}
