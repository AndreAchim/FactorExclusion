mon_Chol <- function(S){
# si S n'a pas déjà des 1 dans sa diagonale, elle est transformée pour assurer cela
# retourne une matrice dont le triangle supérieur est rempli
  if (any(diag(S) != 1)){
    d <- diag(1/sqrt(diag(S)))
    S <- d %*% S %*% d 
  }
  nv <- nrow(S)
  mC <- matrix(0,nrow=nv,ncol=nv)
  mC[1,] <- S[1,]
  for (k in 2:nv){
    if (k>2)
      for (j in 2:(k-1))
        if (mC[j,j]==0)
          mC[j,k] <- 0
        else
          mC[j,k] <- (S[j,k] - (t(mC[1:(j-1),j]) %*% mC[1:(j-1),k])) / mC[j,j] 
    b <- j
    sc <- 1- t(mC[1:(k-1),k]) %*% mC[1:(k-1),k]
    if (sc < 1e-12) sc <- 0
    mC[k,k] <- sqrt(sc)
    a <- 0
  }
  return(mC)
}
  