mon_Chol <- function(R){
  # R (nv,nv) fut produite par F %*% fctR % t(F)
  # retourne une matrice dont seulement les nf premières colonnes ne sont pas nulles
  nv <- nrow(S)
  mC <- matrix(0,nrow=nv,ncol=nv)
  mC[1,1] <- sqrt(R[1,1])
  # browser()
  for (k in 2:nv){
    for (j in 1:(k-1))
      mC[k,j] <- R[k,j] / mC[1,j]
    if (k<=nf)
      mC[k,k] <- sqrt(R[k,k] - sum(mC[k,1:(k-1)]^2))
  }

  #   Rinv <- solve(t(mC[,1:(k-1)]) %*% mC[,1:(k-1)])
  #   z <- t(y) %*% TEMOINS
  #   
  #   # if (k>2)
  #   #   for (j in 2:(k-1))
  #   #     if (mC[j,j]==0)
  #   #       mC[j,k] <- 0
  #   # else
  #     mC[j,k] <- (S[j,k] - (t(mC[1:(j-1),j]) %*% mC[1:(j-1),k])) / mC[j,j] 
  #   b <- j
  #   sc <- 1- t(mC[1:(k-1),k]) %*% mC[1:(k-1),k]
  #   if (sc < 1e-12) sc <- 0
  #   mC[k,k] <- sqrt(sc)
  #   a <- 0
  # }
  return(mC)
}
